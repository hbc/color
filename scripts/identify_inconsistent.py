#!/usr/bin/env python
"""Identify inconsistent calls between VCF replicates as starting point for analysis.
"""
import collections
import glob
import gzip
import pprint
import os
import re
import subprocess
import sys

import pybedtools
import toolz as tz
import numpy as np
from scipy.cluster import vq
import yaml

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.variation import bedutils, vcfutils

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    config["config"] = {}
    groups = organize_vcf_reps(glob.glob(tz.get_in(["inputs", "vcfs"], config)),
                               tz.get_in(["inputs", "namere"], config), config["remap"])
    groups = add_bams(glob.glob(tz.get_in(["inputs", "bams"], config)),
                      tz.get_in(["inputs", "namere"], config), groups, config["remap"])
    bed_file = bedutils.clean_file(tz.get_in(["inputs", "regions"], config), config) + ".gz"
    with utils.chdir(tz.get_in(["dirs", "work"], config)):
        groups = preprocess_vcfs(groups, bed_file, config["resources"])
        pprint.pprint(groups)
        incon = {}
        for name, fnames in groups.items():
            incon[name] = find_inconsistent(name, fnames["vcf"], bed_file)
        incon_check = []
        for name, info in sorted(incon.items(), key=lambda x: np.mean(x[1]["counts"]), reverse=True):
            print name, info["counts"]
            if np.mean(info["counts"]) > 100:
                incon_check.extend(investigate_high_counts(info["summary"], info["vcf_files"]))
        for to_check in incon_check:
            deconvolute_inconsistent(to_check, groups, bed_file)
        disc_bed = identify_shared_discordants(incon)
        filtered_bed = merge_filtered(incon)
        ann_bed = annotate_disc_bed(disc_bed, filtered_bed, config["annotations"])
        remain_disc = check_annotated_disc(ann_bed, config["annotations"])
        identify_discordant_reasons(remain_disc, incon)

        calculate_annotation_overlap(bed_file, filtered_bed, config["annotations"])
        #print ann_bed

# ## Identify causes of discordance with filters

def identify_discordant_reasons(discs, incon):
    """Interrogate discordants to identify potential causes.
    """
    all_reasons = collections.defaultdict(int)
    ctype_stats = collections.defaultdict(list)
    for (dchrom, dstart, dend, samples) in discs:
        print "---"
        #if len(samples) >= 5:
        if False:
            print dchrom, dstart, dend, samples
            for sample in samples:
                for chrom, start, calls in _get_isec_calls(dchrom, dstart, dend, sample, incon):
                    _check_problem_call([(chrom, start, calls, sample,
                                          incon[sample]["vcf_files"])])
        else:
            print dchrom, dstart, dend
            check_regions = []
            for sample in samples:
                ctypes = []
                for chrom, start, calls in _get_isec_calls(dchrom, dstart, dend, sample, incon):
                    call_counts = collections.defaultdict(int)
                    for c in list(calls):
                        call_counts[int(c)] += 1
                    if len(call_counts) > 1:
                        if call_counts[0] == 1:
                            ctype = "false negative"
                        elif call_counts[1] == 1:
                            ctype = "false positive"
                        else:
                            ctype = "mixed"
                        stats = _check_problem_call([(chrom, start, calls, sample,
                                                      incon[sample]["vcf_files"])])
                        ctype_stats[ctype].extend(stats)
                        ctypes.append((ctype, calls))
                        all_reasons[ctype] += 1
                print " -", sample, ctypes
    print dict(all_reasons)
    #pprint.pprint(dict(ctype_stats))

def _get_isec_calls(dchrom, dstart, dend, sample, incon):
    """Retrieve intersection calls for all samples in a specific region.
    """
    with open(tz.get_in([sample, "summary"], incon)) as in_handle:
        for line in in_handle:
            chrom, start, _, _, calls = line.strip().split()
            if chrom == dchrom and int(start) >= dstart and int(start) <= dend:
                yield chrom, start, calls

def _check_problem_call(call_info):
    """Identify potential issues with incorrect calls.
    """
    out = []
    for chrom, start, calls, sample, vcf_files in call_info:
        assert len(calls) == len(vcf_files)
        for cur_vcf in [x for c, x in zip(list(calls), vcf_files) if c == "1"]:
            call = _find_call(chrom, start, cur_vcf)
            parts = call.strip().split("\t")
            ftinfo = {key: val for (key, val) in zip(parts[-2].split(":"), parts[-1].split(":"))}
            depth = sum(int(x) for x in ftinfo["AD"].split(","))
            ref_pl = int(ftinfo["PL"].split(",")[0])
            qd = float(ref_pl) / float(depth)
            stats = {"depth": depth, "qd": "%.1f" % qd}
            info_want = ("GC=",)
            for info_item in parts[-3].split(";"):
                if info_item.startswith(info_want):
                    key, val = info_item.split("=")
                    stats[key] = val
            if depth > 15:
                print "  ", parts[:2] + parts[3:5] + parts[-2:], stats
            out.append(stats)
    return out

def _find_call(chrom, start, vcf_file):
    """Find a specific call at the given position in a bgzipped VCF file.
    """
    with gzip.open(vcf_file) as in_handle:
        for line in (l for l in in_handle if not l.startswith("#")):
            cur_chrom, cur_start = line.split("\t")[:2]
            if cur_chrom == chrom and cur_start == start:
                return line

# ## External annotations

def annotate_disc_bed(in_file, filtered_bed, annotations):
    """Annotate a BED file with potential causes from genomic regions.
    """
    out_file = "%s-annotate.bed" % utils.splitext_plus(in_file)[0]
    if not utils.file_exists(out_file):
        with file_transaction({}, out_file) as tx_out_file:
            names, files = [], []
            names.append("filtered")
            files.append(filtered_bed)
            for name, fname in annotations.items():
                names.append(name)
                files.append(fname)
            names_str = " ".join(names)
            files_str = " ".join(files)
            cmd = ("bedtools annotate -i {in_file} -files {files_str} "
                   "-names {names_str} > {tx_out_file}")
            do.run(cmd.format(**locals()), "Annotate discordant regions")
    return out_file

def check_annotated_disc(in_file, annotations):
    """Provide statistics on discordant variants removed by annotations.
    """
    ann_count = len(annotations) + 1  # all of the annotations plus filtered regions
    explained = collections.defaultdict(int)
    total = 0
    header = []
    remain = []
    with open(in_file) as in_handle:
        for line in in_handle:
            parts = line.rstrip().split("\t")
            lineid = parts[:-ann_count]
            if line.startswith("#"):
                anns = parts[-ann_count:]
                header = anns
            else:
                anns = [float(x) for x in parts[-ann_count:]]
                assert len(header) > 0
                gotit = False
                for i, val in enumerate(anns):
                    ann_name = header[i]
                    if val > 0:
                        explained[ann_name] += 1
                        gotit = True
                if not gotit:
                    chrom, start, end, samples = lineid
                    remain.append((chrom, int(start), int(end), samples.split(",")))
                total += 1
    print len(remain), total, dict(explained)
    return remain

def calculate_annotation_overlap(orig_bed, filtered_bed, annotations):
    """Calculate amount of original BED file falling in supplied annotations.
    """
    full_size = pybedtools.BedTool(gzip.open(orig_bed)).total_coverage()
    full_bed = pybedtools.BedTool(gzip.open(orig_bed))
    for name, fname in [("filtered", filtered_bed)] + annotations.items():
        remain_size = pybedtools.BedTool(gzip.open(orig_bed)).subtract(fname).total_coverage()
        print name, remain_size, full_size, "%.1f" % (float(remain_size) / full_size * 100.0)
        full_bed = full_bed.subtract(fname)
    remain_size = full_bed.total_coverage()
    print "Combined", remain_size, full_size, "%.1f" % (float(remain_size) / full_size * 100.0)

# ## Shared discordant variants

def identify_shared_discordants(incon):
    """Identify discordant variants shared in multiple samples.

    Looks for pervasive issues likely due to algorithmic/genome representation issues.
    Create a BED file of discordant calls and then merge these to identify regions.
    """
    dis_beds = [_isec_summary_to_bed(info["summary"], name, ext="-discordant")
                for name, info in incon.items()]
    work_dir = utils.safe_makedir(os.path.join(os.getcwd(), "cmpsum"))
    merge_disc_bed = _merge_discordant_beds(dis_beds, work_dir)
    return merge_disc_bed

def merge_filtered(incon):
    """Merge all positions filtered in the inputs to have a consistent callset.
    """
    work_dir = utils.safe_makedir(os.path.join(os.getcwd(), "cmpsum"))
    out_file = os.path.join(work_dir, "filtered-merged.bed")
    if not utils.file_exists(out_file):
        with file_transaction({}, out_file) as tx_out_file:
            bed_files_str = " ".join(info["excluded"] for info in incon.values())
            cmd = ("cat {bed_files_str} | sort -k1,1 -k2,2n | "
                   "bedtools merge -d 1 -c 4 -o distinct -i - > {tx_out_file}")
            do.run(cmd.format(**locals()), "Merge filtered regions")
    return out_file

def _merge_discordant_beds(bed_files, work_dir):
    out_file = os.path.join(work_dir, "discordant-merged.bed")
    if not utils.file_exists(out_file):
        with file_transaction({}, out_file) as tx_out_file:
            bed_files_str = " ".join(bed_files)
            cmd = ("cat {bed_files_str} | sort -k1,1 -k2,2n | "
                   "bedtools merge -d 100 -c 4 -o distinct -i - > {tx_out_file}")
            do.run(cmd.format(**locals()), "Merge discordant regions")
    return out_file

def _isec_summary_to_bed(isec_file, name, ext="", require_different=True):
    out_file = "%s%s.bed" % (utils.splitext_plus(isec_file)[0], ext)
    if not utils.file_exists(out_file):
        with file_transaction({}, out_file) as tx_out_file:
            with open(isec_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    for line in in_handle:
                        dline = _prep_discordant_line(line, name, require_different)
                        if dline:
                            out_handle.write("\t".join(str(x) for x in dline) + "\n")
    return out_file

def _prep_discordant_line(line, name, require_different=True):
    """Prepare information on a line from intersection if discordant.
    """
    chrom, start, ref, alt, cmpstr = line.strip().split()
    if not require_different or len(set(list(cmpstr))) > 1:
        start = int(start) - 1
        size = max(len(x) for x in [ref] + alt.split(","))
        return [chrom, start, start + size, name]

# ## Comparisons

def find_inconsistent(name, vcf_files, bed_file):
    """Find inconsistent calls in provided regions of interest for each group.
    """
    cmp_dir = utils.safe_makedir(os.path.join(os.getcwd(), "compare", name))
    isec_dir = os.path.join(cmp_dir, "isec")
    vcf_files_str = " ".join(vcf_files)
    target_count = len(vcf_files) - 1
    if not os.path.exists(isec_dir):
        with file_transaction({}, isec_dir) as tx_isec_dir:
            cmd = ("bcftools isec -f 'PASS,.' {vcf_files_str} -R {bed_file} "
                   "-n -{target_count} -p {tx_isec_dir} -O z")
            do.run(cmd.format(**locals()), "Intersection finding non-consistent calls")
    isec_summary = _calculate_summary(vcf_files, bed_file, cmp_dir)
    excluded = _calculate_excluded(vcf_files, bed_file, cmp_dir, name)
    inconsistent = []
    for i in range(len(vcf_files)):
        with gzip.open(os.path.join(isec_dir, "%04d.vcf.gz" % i)) as in_handle:
            inconsistent.append(sum(1 for l in in_handle if not l.startswith("#")))
    return {"counts": inconsistent,
            "vcf_files": vcf_files,
            "summary": isec_summary,
            "excluded": excluded}

def _calculate_summary(vcf_files, bed_file, cmp_dir):
    """Summarize all variants called in the VCF files as a bcftools isec output file.
    """
    file_list = os.path.join(cmp_dir, "input_files.txt")
    if not utils.file_exists(file_list):
        with file_transaction({}, file_list) as tx_out_file:
            with open(tx_out_file, "w") as out_handle:
                out_handle.write("\n".join(vcf_files))
    out_file = os.path.join(cmp_dir, "all_variants.txt")
    if not utils.file_exists(out_file):
        with file_transaction({}, out_file) as tx_out_file:
            vcf_files_str = " ".join(vcf_files)
            cmd = ("bcftools isec -f 'PASS,.' -n '+1' -o {tx_out_file} -R {bed_file} {vcf_files_str}")
            do.run(cmd.format(**locals()), "Variant comparison summary")
    return out_file

def _calculate_excluded(vcf_files, bed_file, cmp_dir, name):
    """Calculate a BED file of regions filtered in any input sample.
    """
    out_file = os.path.join(cmp_dir, "filtered_positions.txt")
    if not utils.file_exists(out_file):
        with file_transaction({}, out_file) as tx_out_file:
            vcf_files_str = " ".join(vcf_files)
            cmd = ("bcftools isec -f 'ColorCustom' -n '+1' -o {tx_out_file} -R {bed_file} {vcf_files_str}")
            do.run(cmd.format(**locals()), "Filtered summary")
    return _isec_summary_to_bed(out_file, name, require_different=False)

# ## Investigate comparison problems

def deconvolute_inconsistent(problem_vcf, groups, bed_file):
    """Identify correct sample classification for problematic VCF.
    """
    cmps = []
    for name, info in groups.items():
        for cmp_vcf in info["vcf"]:
            if cmp_vcf != problem_vcf:
                concordant = _count_concordant(problem_vcf, cmp_vcf, bed_file)
                cmps.append((concordant, cmp_vcf))
    cmps.sort(reverse=True)
    print problem_vcf
    for x in cmps[:5]:
        print " -", x

def _count_concordant(orig_vcf, cmp_vcf, bed_file):
    """Identify concordant calls between two VCF files.
    """
    cmd = ("bcftools gtcheck --GTs-only 1 --regions-file {bed_file} "
           "--genotypes {cmp_vcf} {orig_vcf} 2> /dev/null")
    bout = subprocess.check_output(cmd.format(**locals()), shell=True)
    parts = bout.split("\n")[-2].split()
    discordant = int(float(parts[1]))
    total = int(parts[3])
    return total - discordant

def investigate_high_counts(summary_file, vcf_files):
    counts = []
    with open(summary_file) as in_handle:
        for line in in_handle:
            calls = [int(x) for x in list(line.strip("\n").split("\t")[-1])]
            counts.append(calls)
    data = np.array(counts).transpose()
    centroids = vq.kmeans(data, 2)[0]
    codes = vq.vq(data, centroids)[0]
    # http://stackoverflow.com/questions/1518522/python-most-common-element-in-a-list
    common_code = max(set(codes), key=list(codes).count)
    problems = []
    for code, vcf_file in zip(codes, vcf_files):
        if code != common_code:
            print " -", os.path.basename(vcf_file)
            problems.append(vcf_file)
    return problems

# ## Pre-processing

def preprocess_vcfs(groups, bed_file, resources):
    """Pre-process VCFs to ensure they are bgzipped and tabix indexed.
    """
    out_dir = utils.safe_makedir(os.path.join(os.getcwd(), "prep"))
    filter_stats = []
    for name, fnames in groups.items():
        prep_vcfs = []
        for vcf_file in fnames["vcf"]:
            out_file = os.path.join(out_dir, os.path.basename(vcf_file))
            if not out_file.endswith(".vcf.gz"):
                out_file = out_file.replace("vcf.gz", ".vcf.gz")
            _prep_vcf(vcf_file, out_file)
            ann_file = _annotate_vcf(out_file, _find_bam_file(out_file,  fnames.get("bam", [])),
                                     bed_file, resources)
            filter_file, cur_stats = _filter_vcf(ann_file)
            filter_stats.append(cur_stats)
            prep_vcfs.append(filter_file)
        assert None not in prep_vcfs
        groups[name]["vcf"] = prep_vcfs
    _analyze_filtered_stats(filter_stats)
    return groups

def _analyze_filtered_stats(stats):
    pcts = [x["pct"] for x in stats]
    origs = [x["orig"] for x in stats]
    finals = [x["final"] for x in stats]
    removed = [x["orig"] - x["final"] for x in stats]
    print "Removed percent. range: %.1f%%-%.1f%%; median: %.1f%%" % (min(pcts), max(pcts), np.median(pcts))
    print "Removed counts. range: %s-%s; median: %s" % (min(removed), max(removed), np.median(removed))
    print "Original counts. range: %s-%s; median: %s" % (min(origs), max(origs), np.median(origs))
    print "Final counts. range: %s-%s; median: %s" % (min(finals), max(finals), np.median(finals))

def _prep_vcf(in_file, out_file):
    if not utils.file_exists(out_file):
        with file_transaction({}, out_file) as tx_out_file:
            cmd = ("gunzip -c {in_file} | bgzip -c > {tx_out_file}")
            do.run(cmd.format(**locals()), "Convert VCFs to bgzipped")
    vcfutils.bgzip_and_index(out_file, {})
    return out_file

def _filter_vcf(orig_file):
    """Filter VCF with bcftools, providing count summary of items removed.
    """
    expr = ('SUM(AD[*]) < 15 || '
            'PL[0] / SUM(AD[*]) <= 3.0 || '
            'GC < 27.0 || GC > 73.0')
    out_file = "%s-filter%s" % utils.splitext_plus(orig_file)
    if not utils.file_exists(out_file):
        with file_transaction({}, out_file) as tx_out_file:
            cmd = ("bcftools filter -O z -o {tx_out_file} "
                   "-e '{expr}' -s 'ColorCustom' {orig_file}")
            do.run(cmd.format(**locals()), "Hard filter VCF")
    vcfutils.bgzip_and_index(out_file, {})
    def count(f):
        with gzip.open(f) as h:
            return sum(1 for line in h if not line.startswith("#") and line.split("\t")[6] in  ["PASS", "."])
    removed_stats = {"orig": count(orig_file), "final": count(out_file)}
    removed_stats["pct"] = float(removed_stats["final"]) * 100.0 / removed_stats["orig"]
    return out_file, removed_stats

def _annotate_vcf(vcf_file, bam_file, bed_file, resources):
    """Provide GATK annotations missing from the original VCF.
    """
    assert bam_file, "BAM file not found for %s" % vcf_file
    out_file = "%s-annotated%s" % utils.splitext_plus(vcf_file)
    if not utils.file_exists(out_file):
        with file_transaction({}, out_file) as tx_out_file:
            gatk_jar = resources.get("gatk")
            ref_file = resources.get("ref_file")
            bed_file = _clean_bed_file(bed_file.replace(".bed.gz", ".bed"))
            cmd = ("java -jar {gatk_jar} -T VariantAnnotator -R {ref_file} "
                   "-L {bed_file} "
                   "-A GCContent --variant {vcf_file} --out {tx_out_file}")
            do.run(cmd.format(**locals()), "Annotate VCF")
    return out_file

def _clean_bed_file(orig_bed):
    out_bed = "%s-gatkclean%s" % utils.splitext_plus(orig_bed)
    if not utils.file_exists(out_bed):
        with file_transaction({}, out_bed) as tx_out_file:
            with open(orig_bed) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    for line in in_handle:
                        parts = line.split("\t")
                        if int(parts[2]) > int(parts[1]):
                            out_handle.write(line)
    return out_bed

def _find_bam_file(vcf_file, bam_files):
    search_file = os.path.basename(vcf_file).replace(".vcf.gz", ".bam")
    our_bam = [x for x in bam_files if os.path.basename(x) == search_file]
    if len(our_bam) > 0:
        assert len(our_bam) == 1
        return our_bam[0]
    else:
        return None

# ## Get files to work on

def add_bams(bam_files, name_re, groups, remap):
    """Add BAM information existing grouped VCFs.
    """
    pat = re.compile(name_re)
    for bam_file in bam_files:
        try:
            name = remap[utils.splitext_plus(os.path.basename(bam_file))[0]]
        except KeyError:
            name = pat.search(os.path.basename(bam_file)).group(0)
        if name in groups:
            cur_group = groups[name]
            try:
                cur_group["bam"].append(bam_file)
            except KeyError:
                cur_group["bam"] = [bam_file]
            groups[name] = cur_group
    return groups

def organize_vcf_reps(vcf_files, name_re, remap):
    """Retrieve VCFs analyzed as replicates.
    """
    pat = re.compile(name_re)
    by_name = collections.defaultdict(list)
    for vcf_file in vcf_files:
        try:
            name = remap[utils.splitext_plus(os.path.basename(vcf_file))[0]]
        except KeyError:
            name = pat.search(os.path.basename(vcf_file)).group(0)
        by_name[name].append(vcf_file)
    out = {}
    for name, vcfs in by_name.items():
        if len(vcfs) > 1:
            out[name] = {"vcf": vcfs}
    return out

if __name__ == "__main__":
    main(*sys.argv[1:])
