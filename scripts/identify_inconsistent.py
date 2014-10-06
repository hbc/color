#!/usr/bin/env python
"""Identify inconsistent calls between VCF replicates as starting point for analysis.
"""
import collections
import glob
import pprint
import os
import re
import subprocess
import sys

import toolz as tz
import yaml

from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.provenance import do
from bcbio.variation import vcfutils

def main(config_file):
    with open(config_file) as in_handle:
        config = yaml.load(in_handle)
    groups = organize_vcf_reps(glob.glob(tz.get_in(["inputs", "vcfs"], config)),
                               tz.get_in(["inputs", "namere"], config))
    groups = add_bams(glob.glob(tz.get_in(["inputs", "bams"], config)),
                      tz.get_in(["inputs", "namere"], config), groups)
    with utils.chdir(tz.get_in(["dirs", "work"], config)):
        groups = preprocess_vcfs(groups)
    pprint.pprint(groups)

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

def add_bams(bam_files, name_re, groups):
    """Add BAM information existing grouped VCFs.
    """
    pat = re.compile(name_re)
    for bam_file in bam_files:
        name = pat.search(os.path.basename(bam_file)).group(0)
        if name in groups:
            cur_group = groups[name]
            try:
                cur_group["bam"].append(bam_file)
            except KeyError:
                cur_group["bam"] = [bam_file]
            groups[name] = cur_group
    return groups

def organize_vcf_reps(vcf_files, name_re):
    """Retrieve VCFs analyzed as replicates.
    """
    pat = re.compile(name_re)
    by_name = collections.defaultdict(list)
    for vcf_file in vcf_files:
        name = pat.search(os.path.basename(vcf_file)).group(0)
        by_name[name].append(vcf_file)
    out = {}
    for name, vcfs in by_name.items():
        if len(vcfs) > 1:
            out[name] = {"vcf": vcfs}
    return out

if __name__ == "__main__":
    main(*sys.argv[1:])
