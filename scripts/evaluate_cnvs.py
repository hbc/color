#!/usr/bin/env python
"""Evaluate called CNVs organized by families.
"""
import collections
import json
import re
import sys

import toolz as tz
import yaml

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    sample_cnvs = read_cnvs(tz.get_in(["inputs", "cnvs"], config),
                            tz.get_in(["inputs", "namere"], config), config["remap"])
    for sample, cnvs in sample_cnvs.iteritems():
        summarize_cnvs(sample, cnvs)

def summarize_cnvs(sample, cnvs):
    by_replicate = collections.defaultdict(list)
    for cnv in cnvs:
        by_replicate[cnv["SAMPLE_ID"]].append(cnv)
    print "---", sample
    rep_names = sorted(by_replicate.keys())
    for rep_name in rep_names:
        items = by_replicate[rep_name]
        print " ", rep_name, len(items)

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
