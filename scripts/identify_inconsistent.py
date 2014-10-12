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
        groups = preprocess_vcfs(groups)
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
        ann_bed = annotate_disc_bed(disc_bed, config["annotations"])
        check_annotated_disc(ann_bed, config["annotations"])
        calculate_annotation_overlap(bed_file, config["annotations"])
        print ann_bed

# ## External annotations

def annotate_disc_bed(in_file, annotations):
    """Annotate a BED file with potential causes from genomic regions.
    """
    out_file = "%s-annotate.bed" % utils.splitext_plus(in_file)[0]
    if not utils.file_exists(out_file):
        with file_transaction({}, out_file) as tx_out_file:
            names, files = [], []
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
    ann_count = len(annotations)
    explained = collections.defaultdict(int)
    remain = 0
    total = 0
    header = []
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
                    remain += 1
                    size = int(lineid[2]) - int(lineid[1])
                    print size, lineid
                total += 1
    print remain, total, dict(explained)

def calculate_annotation_overlap(orig_bed, annotations):
    """Calculate amount of original BED file falling in supplied annotations.
    """
    full_size = pybedtools.BedTool(gzip.open(orig_bed)).total_coverage()
    full_bed = pybedtools.BedTool(gzip.open(orig_bed))
    for name, fname in annotations.items():
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
    dis_beds = [_isec_summary_to_bed(info["summary"], name)
                for name, info in incon.items()]
    work_dir = utils.safe_makedir(os.path.join(os.getcwd(), "cmpsum"))
    merge_disc_bed = _merge_discordant_beds(dis_beds, work_dir)
    return merge_disc_bed

def _merge_discordant_beds(bed_files, work_dir):
    out_file = os.path.join(work_dir, "discordant-merged.bed")
    if not utils.file_exists(out_file):
        with file_transaction({}, out_file) as tx_out_file:
            bed_files_str = " ".join(bed_files)
            cmd = ("cat {bed_files_str} | sort -k1,1 -k2,2n | "
                   "bedtools merge -d 100 -c 4 -o distinct -i - > {tx_out_file}")
            do.run(cmd.format(**locals()), "Merge discordant regions")
    return out_file

def _isec_summary_to_bed(isec_file, name):
    out_file = "%s-discordant.bed" % utils.splitext_plus(isec_file)[0]
    if not utils.file_exists(out_file):
        with file_transaction({}, out_file) as tx_out_file:
            with open(isec_file) as in_handle:
                with open(tx_out_file, "w") as out_handle:
                    for line in in_handle:
                        dline = _prep_discordant_line(line, name)
                        if dline:
                            out_handle.write("\t".join(str(x) for x in dline) + "\n")
    return out_file

def _prep_discordant_line(line, name):
    """Prepare information on a line from intersection if discordant.
    """
    chrom, start, ref, alt, cmpstr = line.strip().split()
    if len(set(list(cmpstr))) > 1:
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
            cmd = ("bcftools isec {vcf_files_str} -R {bed_file} "
                   "-n -{target_count} -p {tx_isec_dir} -O z")
            do.run(cmd.format(**locals()), "Intersection finding non-consistent calls")
    isec_summary = _calculate_summary(vcf_files, bed_file, cmp_dir)
    inconsistent = []
    for i in range(len(vcf_files)):
        with gzip.open(os.path.join(isec_dir, "%04d.vcf.gz" % i)) as in_handle:
            inconsistent.append(sum(1 for l in in_handle if not l.startswith("#")))
    return {"counts": inconsistent,
            "vcf_files": vcf_files,
            "summary": isec_summary}

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
            cmd = ("bcftools isec -n '+1' -o {tx_out_file} -R {bed_file} {vcf_files_str}")
            do.run(cmd.format(**locals()), "Variant comparison summary")
    return out_file

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

def preprocess_vcfs(groups):
    """Pre-process VCFs to ensure they are bgzipped and tabix indexed.
    """
    out_dir = utils.safe_makedir(os.path.join(os.getcwd(), "prep"))
    for name, fnames in groups.items():
        prep_vcfs = []
        for vcf_file in fnames["vcf"]:
            out_file = os.path.join(out_dir, os.path.basename(vcf_file))
            if not out_file.endswith(".vcf.gz"):
                out_file = out_file.replace("vcf.gz", ".vcf.gz")
            _prep_vcf(vcf_file, out_file)
            prep_vcfs.append(out_file)
        groups[name]["vcf"] = prep_vcfs
    return groups

def _prep_vcf(in_file, out_file):
    if not utils.file_exists(out_file):
        with file_transaction({}, out_file) as tx_out_file:
            cmd = ("gunzip -c {in_file} | bgzip -c > {tx_out_file}")
            do.run(cmd.format(**locals()), "Convert VCFs to bgzipped")
    vcfutils.bgzip_and_index(out_file, {})
    return out_file

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
