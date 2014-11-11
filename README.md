Scripts for organizing and evaluating panel replicate variant calls to identify
and categorize non-consistent calls.

Running
-------
The analysis runs with a single command on the AWS machine with downloaded BAMs
and VCFs. It uses an isolated bcbio installation to provide a Python
installation with the libraries and third party tools:

   export PATH=/encrypted/permanent/bcbio/bin:$PATH
   git clone https://github.com/hbc/color.git
   cd /encrypted/permanent/work
   /encrypted/permanent/bcbio/anaconda/bin/python /path/to/color/scripts/identify_inconsistent.py /path/to/color/config/awsinputs.yaml

Beyond the inputs from Color (VCFs, BAMs and BED file of regions) we also use
two annotation files. They are located in `/encrypted/permanent/regions` and you
can fetch them from scratch with:

   wget -O - https://s3.amazonaws.com/gemini-annotations/hg19.rmsk.bed.gz | gunzip -c | sed 's/^chr//' | bgzip -c > GRCh37.rmsk.bed.gz
   tabix -f -p bed GRCh37.rmsk.bed.gz
   wget -O - https://github.com/lh3/varcmp/raw/master/scripts/LCR-hs37d5.bed.gz | gunzip -c | bgzip -c > LCR.bed.gz
   tabix -f -p bed LCR.bed.gz

Data files
----------

- `outputs/filter_counts.csv` -- Breakdown of the impact of filters on each of the
  samples and replicates. The categories are:

  - True positives (tp) -- Correct variants that are not filtered.
  - False positives (fp) -- Discordant variants that are not filtered. We
    took the most conservative metric for this, leaving in anything
    inconsistent across calls. This is conservative since we can get a
    handful of cases where you have, say, 3 samples a, b and c. b and c
    don't call a variant because of low depth, while it gets called in 'a'. In
    this case 'a' doesn't get filtered because it is fine in that sample,
    but will appear as a FP since the low depth in b and c cause it to get
    marked as discordant.
  - True negatives (tn) -- Variants the filters remove that are also
    discordant. The output splits these by the approach that filters them
    (depth, region or the full review criteria). The filters are
    cumulative, so each column reflects additional variants removed with
    the new filter.
  - False negatives (fn) -- Variants the filters remove that are not
    marked discordant. These are again split by approach.

- `outputs/orig_counts.csv` -- Breakdowns of true and false positives before
  applying any filtering. Columns are the same as in the `filter_counts.csv` file.

- `outputs/discordant_merged.bed` -- A combined BED file of all the discordant regions
  from all replicates, prior to filtering.

- `outputs/GRCh37-repeats.bed.gz` -- Merged BED file for repeat masked and low
  complexity regions used as inputs for filtering.

Filters
-------

The 3 filters used are:

- ColorRemoveDepth: Depth and Quality based filters that mark variants
  likely to be false positives.
- ColorRemoveRegions: Low complexity and repeat regions that mark
  variants likely to be false positives.
- ColorReview: More extensive low complexity and repeat region
  specifications that mark additional variants in regions that should
  get reviewed.

Add these filters to a VCF file with:

    echo '##INFO=<ID=RPT,Number=.,Type=String,Description="Repeat track overlaps">' > header-file.txt
    bcftools annotate -a /encrypted/permanent/regions/GRCh37-repeats.bed.gz \
             -h header-file.txt orig.vcf.gz -O z -o orig-repeats.vcf.gz
    tabix -p vcf orig-repeats.vcf.gz
    bcftools filter -O z -o orig-repeats-filter1.vcf.gz -m '+' -s 'ColorRemoveDepth' \
             -e 'SUM(AD[*]) < 15 || PL[0] / SUM(AD[*]) <= 3.0'
    tabix -p vcf orig-repeats-filter1.vcf.gz
    bcftools filter -O z -o orig-repeats-filter2.vcf.gz -m '+' -s 'ColorRemoveRegions' \
             -e 'RPT[*] = "lcr" && RPT[*] = "rmsk"'
    tabix -p vcf orig-repeats-filter2.vcf.gz
    bcftools filter -O z -o orig-repeats-filter3.vcf.gz -m '+' -s 'ColorReview' \
             -e 'SUM(AD[*]) < 15 || PL[0] / SUM(AD[*]) <= 3.0 || GC < 20.0 || GC > 77.0 || RPT[*] = "rmsk" || RPT[*] = "lcr"'
    tabix -p vcf orig-repeats-filter3.vcf.gz
