#!/usr/bin/env python
"""Evaluate called CNVs organized by families.
"""
import collections
import glob
import json
import os
import re
import sys

import toolz as tz
import yaml

from bcbio import utils
from bcbio.provenance import do

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    sample_cnvs = read_cnvs(tz.get_in(["inputs", "cnvs"], config),
                            tz.get_in(["inputs", "namere"], config), config["remap"])
    all_calls = []
    passed_cnvs = []
    for sample, cnvs in sorted(sample_cnvs.iteritems()):
        cur_calls, cur_cnvs = summarize_cnvs(sample, cnvs)
        if cur_calls:
            all_calls.extend(cur_calls)
            for i, cnv in enumerate(cur_cnvs):
                cnv["iname"] = "%s-%s" % (cnv["name"], i)
                passed_cnvs.append(cnv)
    cnv_dir = utils.safe_makedir(os.path.join(os.getcwd(), "cnv_plots"))
    cnv_bed = combine_cnvs(all_calls, cnv_dir, tz.get_in(["resources", "ref_file"], config))
    background_bed = get_background_cnvs(tz.get_in(["inputs", "regions"], config), cnv_dir)
    for cnv in passed_cnvs:
        plot_cnv(cnv, cnv_bed, tz.get_in(["inputs", "bams"], config))
        plot_cnv(cnv, background_bed, tz.get_in(["inputs", "bams"], config), "-background")

def summarize_cnvs(sample, cnvs):
    by_replicate = collections.defaultdict(list)
    for cnv in cnvs:
        by_replicate[cnv["SAMPLE_ID"]].append(cnv)
    rep_names = sorted(by_replicate.keys())
    all_calls = []
    if len(rep_names) > 1:
        print "---", sample
        total_count = sum([len(xs) for xs in by_replicate.itervalues()])
        if total_count < 100:
            all_calls.extend(_examine_replicates(by_replicate, rep_names))
        else:
            for rep_name in rep_names:
                items = by_replicate[rep_name]
                print " ", rep_name, len(items)
    return all_calls, [xs[0] for xs in by_replicate.values()]

def _examine_replicates(cnvs_by_rep, rep_names):
    """Examine CNV calls organized by replicates.
    """
    calls = []
    for rep in rep_names:
        print " ", rep
        for cnv in cnvs_by_rep[rep]:
            if ((cnv["status"] == "filtered" and int(cnv["CopyNumber"]) > 2)
                  or float(cnv["mean.base.spacing"]) > 2.8):
                flag = " ColorFilter"
            else:
                flag = ""
            print "  %s %s %s %s %s%s" % (cnv["chr.id"], cnv["start"], cnv["end"], cnv["CopyNumber"],
                                          cnv["size"], flag)
            calls.append((cnv["chr.id"], int(cnv["start"]), int(cnv["end"])))
            if False:
                import pprint
                pprint.pprint(cnv)
    return calls

def plot_cnv(cnv, cnv_bed, bam_glob, ext=""):
    """Create a plot of coverage in all of the supplied CNV regions.
    """
    bam_files = [x for x in glob.glob(bam_glob) if x.find(cnv["SAMPLE_ID"]) >= 0]
    assert len(bam_files) == 1, (cnv["SAMPLE_ID"], bam_files)
    bam_file = bam_files[0]
    if not utils.file_exists(bam_file + ".bai"):
        do.run(["sambamba", "index", bam_file], "Index BAM file")
    cov_file = os.path.join(os.path.dirname(cnv_bed), "%s%s.bed" % (cnv["iname"], ext))
    cmd = ("sambamba view -f bam -L {cnv_bed} {bam_file} | "
           "bedtools coverage -abam - -b {cnv_bed} -d > {cov_file}")
    if not utils.file_exists(cov_file):
        do.run(cmd.format(**locals()), "Calculate coverage in regions")
        plot_script = os.path.join(os.path.dirname(__file__), "plot_cnv_coverage.py")
        do.run([sys.executable, plot_script, cov_file], "Plot coverage")

def combine_cnvs(calls, cnv_dir, ref_file):
    orig_file = os.path.join(cnv_dir, "called_cnvs.bed")
    out_file = "%s-merged%s" % os.path.splitext(orig_file)
    def by_chrom(calls):
        chrom, start, end = calls
        return (int(chrom.replace("chr", "")), int(start), int(end))
    with open(orig_file, "w") as out_handle:
        calls = sorted(calls, key=by_chrom)
        for chrom, start, end in calls:
            out_handle.write("%s\t%s\t%s\n" % (chrom.replace("chr", ""), start, end))
    cmd = ("bedtools slop -b 2000 -i {orig_file} -g {ref_file}.fai | bedtools merge -i - > {out_file}")
    do.run(cmd.format(**locals()), "Merge CNVs")
    return out_file

def get_background_cnvs(orig_bed, cnv_dir):
    """Subsample a set of background region coverage for comparison to CNVs.
    """
    out_file = os.path.join(cnv_dir, "background-sample.bed")
    if not utils.file_exists(out_file):
        cmd = "bedtools sample -n 15 -i {orig_bed} > {out_file}"
        do.run(cmd.format(**locals()), "Subset background")
    return out_file

def read_cnvs(in_file, name_re, remap):
    """Read CNV information from input file and group by sample.
    """
    pat = re.compile(name_re)
    by_sample = collections.defaultdict(list)
    with open(in_file) as in_handle:
        for line in in_handle:
            chrom, start, end, info = line.strip().split("\t")
            info = json.loads(info)
            try:
                name = remap[info["SAMPLE_ID"]]
            except KeyError:
                name = pat.search(info["SAMPLE_ID"]).group(0)
            info["name"] = name
            by_sample[name].append(info)
    return by_sample

if __name__ == "__main__":
    main(*sys.argv[1:])
