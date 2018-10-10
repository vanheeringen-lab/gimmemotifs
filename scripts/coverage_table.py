#!/usr/bin/env python
import multiprocessing
import os
import sys
import argparse
from io import StringIO

import pysam
import numpy as np
import pandas as pd
from sklearn.preprocessing import scale

from gimmemotifs import __version__
from fluff.fluffio import load_heatmap_data

def make_table(peakfile, datafiles, window, log_transform=True, scale_table=True, rmdup=True, rmrepeats=True):
    for x in datafiles:
        if not os.path.isfile(x):
            print("ERROR: Data file '{0}' does not exist".format(x))
            sys.exit(1)
    for x in datafiles:
        if '.bam' in x and not os.path.isfile("{0}.bai".format(x)):
            print("Data file '{0}' does not have an index file. Creating an index file for {0}.".format(x))
            pysam.index(x)

    print("Loading data", file=sys.stderr)
    data = {}
    try:
        # Load data in parallel
        pool = multiprocessing.Pool(processes=12)
        jobs = []
        for datafile in datafiles:
            jobs.append(pool.apply_async(load_heatmap_data, args=(
            peakfile, datafile, 1, window // 2, window // 2, rmdup, False, rmrepeats,
            None, False, None)))
        for job in jobs:
            track, regions, profile, guard = job.get()
            data[os.path.splitext(track)[0]] = profile[:,0]
    except Exception as e:
        sys.stderr.write("Error loading data in parallel, trying serial\n")
        sys.stderr.write("Error: {}\n".format(e))
        for datafile in datafiles:
            track, regions, profile, guard = load_heatmap_data(peakfile, datafile, 1, window // 2,
                                                               window // 2, rmdup, False, rmrepeats,
                                                               None, False, None)
            data[os.path.splitext(track)[0]] = profile[:,0]
    
    # Create DataFrame with regions as index
    regions = ["{}:{}-{}".format(*region[:3]) for region in regions]
    df = pd.DataFrame(data, index=regions)

    if log_transform:
        print("Log transform", file=sys.stderr)
        df = np.log1p(df)
    if scale_table:
        print("Scale", file=sys.stderr)
        df[:] = scale(df, axis=0)

    return df


if __name__ == "__main__":
    description = """
    GimmeMotifs v{0} - coverage_table
    """.format(__version__)

    parser = argparse.ArgumentParser(
                                    description=description,
                                    )

    parser.add_argument("-p", "--peaks",
                         required=True,
                         dest="peakfile",
                         help="BED file containing peaks",
                         metavar="FILE",
                         default=None)
    parser.add_argument("-d", "--datafile",
                         required=True,
                         dest="datafiles",
                         help="data files (BAM, BED or bigWig format)",
                         metavar="FILE",
                         nargs='*')
    parser.add_argument("-w", "--window",
                         dest="window",
                         help="window size (default 200)",
                         metavar="WINDOW",
                        type=int,
                        default=200)
    parser.add_argument("-l", "--logtransform",
                         dest="log_transform",
                         help="Log transfrom",
                         default=False,
                         action="store_true")
    parser.add_argument("-s", "--scale",
                         dest="scale_table",
                         help="Scale per datafile",
                         default=False,
                         action="store_true")
    parser.add_argument("-D",
                         dest="rmdup",
                         help="keep duplicate reads (removed by default)",
                         default=True,
                         action="store_false")
    parser.add_argument("-R",
                         dest="rmrepeats",
                         help="keep reads with mapq 0 (removed by default) ",
                         action="store_false",
                         default=True)


    args = parser.parse_args()
    peakfile = args.peakfile
    datafiles = args.datafiles
    df = make_table(
            peakfile, 
            datafiles, 
            args.window,
            log_transform=args.log_transform,
            scale_table=args.scale_table,
            rmdup=args.rmdup,
            rmrepeats=args.rmrepeats
            )

output = StringIO()
df.to_csv(output, sep="\t", float_format="%0.5f")
output.seek(0)
print(output.read())
