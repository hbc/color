"""Plot read coverage across regions calculated via bedtools coverage -d.

Usage:
  plot_cnv_coverage.py <bedtools coverage -d output file>
"""
import os
import sys

import matplotlib
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import pandas as pd

def main(in_file):
    out_file = "%s.pdf" % os.path.splitext(in_file)[0]
    with open(in_file) as in_handle:
        df = pd.read_table(in_handle, header=None, names=["chrom", "start", "end", "offset", "coverage"])
        df["position"] = df["start"] + df["offset"]
    with PdfPages(out_file) as pdf_out:
        sns.despine()
        for coords, regiondf in df.groupby(["chrom", "start", "end"]):
            region = "%s:%s-%s" % coords
            size = int(coords[-1]) - int(coords[-2])
            plot = regiondf.plot(x="position", y="coverage", kind="line", legend=False)
            plot.set_title("%s (%s bp)" % (region, size))
            plot.set_ylim(0)
            plot.get_xaxis().set_major_formatter(
                matplotlib.ticker.FuncFormatter(lambda x, p: '{:,}'.format(int(x))))
            plot.set_xlabel("")
            plot.set_ylabel("coverage")
            plot.get_xaxis().set_major_locator(matplotlib.ticker.MaxNLocator(6))
            pdf_out.savefig(plot.get_figure())

if __name__ == "__main__":
    main(sys.argv[1])
