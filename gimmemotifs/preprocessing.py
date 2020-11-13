# Copyright (c) 2009-2020 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.

""" Data preprocessing to create GimmeMotifs input. """
# Python imports
import os
import sys
import logging
import multiprocessing
from tempfile import NamedTemporaryFile

# External imports
import genomepy
import numpy as np
import pysam
from fluff.fluffio import load_heatmap_data
import pandas as pd
from pybedtools import BedTool
from sklearn.preprocessing import scale
import qnorm
from tqdm.auto import tqdm

# gimme imports
from gimmemotifs.utils import determine_file_type

logger = logging.getLogger("gimme.preprocessing")


def coverage_table(
    peakfile,
    datafiles,
    window,
    log_transform=True,
    normalization="none",
    top=0,
    topmethod="var",
    rmdup=True,
    rmrepeats=True,
    ncpus=12,
):
    for x in datafiles:
        if not os.path.isfile(x):
            print("ERROR: Data file '{0}' does not exist".format(x))
            sys.exit(1)
    for x in datafiles:
        if ".bam" in x and not os.path.isfile("{0}.bai".format(x)):
            print(
                "Data file '{0}' does not have an index file."
                " Creating an index file for {0}.".format(x)
            )
            pysam.index(x)

    logger.info("Loading data")
    data = {}
    try:
        # Load data in parallel
        pool = multiprocessing.Pool(processes=ncpus)
        jobs = []
        for datafile in datafiles:
            jobs.append(
                pool.apply_async(
                    load_heatmap_data,
                    args=(
                        peakfile,
                        datafile,
                        1,
                        window // 2,
                        window // 2,
                        rmdup,
                        False,
                        rmrepeats,
                        None,
                        False,
                        None,
                    ),
                )
            )
        for job in tqdm(jobs):
            track, regions, profile, guard = job.get()
            data[os.path.splitext(track)[0]] = profile[:, 0]
    except Exception as e:
        sys.stderr.write("Error loading data in parallel, trying serial\n")
        sys.stderr.write("Error: {}\n".format(e))
        for datafile in tqdm(datafiles):
            track, regions, profile, guard = load_heatmap_data(
                peakfile,
                datafile,
                1,
                window // 2,
                window // 2,
                rmdup,
                False,
                rmrepeats,
                None,
                False,
                None,
            )
            data[os.path.splitext(track)[0]] = profile[:, 0]

    # Create DataFrame with regions as index
    regions = ["{}:{}-{}".format(*region[:3]) for region in regions]
    df = pd.DataFrame(data, index=regions)

    if log_transform:
        logger.info("Log transform")
        df = np.log1p(df)
    if normalization == "scale":
        logger.info("Normalization by scaling")
        df[:] = scale(df, axis=0)
    if normalization == "quantile":
        logger.info("Normalization by quantile normalization")
        df = qnorm.quantile_normalize(df)
    else:
        logger.info("No normalization")

    if top > 0:
        if topmethod == "var":
            idx = df.var(1).sort_values().tail(top).index
        elif topmethod == "std":
            idx = df.std(1).sort_values().tail(top).index
        elif topmethod == "mean":
            idx = df.mean(1).sort_values().tail(top).index
        elif topmethod == "random":
            idx = df.sample(top).index
        else:
            raise ValueError(
                "unknown method {} for selecting regions".format(topmethod)
            )
        df = df.loc[idx]
    return df


def read_peak_file_to_df(fname):
    """
    Read a MACS2 summits.bed or narrowPeak file and return a DataFrame.

    Parameters
    ----------
    fname : str
        Filename.

    Returns
    -------
    pandas.DataFrame
        DataFrame with summits.
    """
    summit_header = ["chrom", "start", "end", "name", "value"]
    ftype = determine_file_type(fname)
    # Read MACS2 summit files
    if ftype == "narrowpeak":
        header = [
            "chrom",
            "start",
            "end",
            "name",
            "value",
            "strand",
            "signalValue",
            "pval",
            "qval",
            "peak",
        ]
        df = pd.read_table(fname, names=header, dtype={"chrom": "str"})
        df["chrom"] = df["chrom"].astype(str)

        # get the summit
        df["start"] = df["start"].astype(int) + df["peak"].astype(int)
        df["end"] = df["start"] + 1

        # qva
        df["value"] = df["qval"]
        df = df[summit_header]
    elif ftype == "bed":
        df = pd.read_table(fname, names=summit_header, dtype={"chrom": "str"})
        if ((df["end"] - df["start"]) != 1).sum() != 0:
            raise ValueError(f"{fname} does not contain summits.")
    else:
        raise ValueError(
            f"Can't determine file type of {fname}. "
            "Is the file a narrowPeak or summits.bed file?"
        )

    df["experiment"] = df.loc[0, "name"].replace("_peak_1", "")
    df["log_value"] = np.log1p(df["value"])
    df["log_value_scaled"] = scale(df[["log_value"]])
    return df


def combine_peaks(peaks, genome, window, scale_value):
    """
    Combine multiple MACS2 summit files and returns the summit
    with the maximum value.

    Parameters
    ----------
    peaks : list
        List with summit file names.

    genome : str
        Genome file name. Either a file with chromosome sizes or a genomepy
        genome name.

    window : int
        Window size. Summits will be extended to this size before merging.

    scale_value : bool
        Scale summit values before taking the maximum.

    Returns
    -------
    summits : pandas.DataFrame
        DataFrame with summits.
    """
    try:
        g = genomepy.Genome(genome)
        genome = g.sizes_file
    except Exception:
        pass

    dfs = [read_peak_file_to_df(fname) for fname in peaks]
    df_all = pd.concat(dfs)

    check_col = "log_value"
    if scale_value:
        check_col = "log_value_scaled"

    # store summit location + associated value in col4
    df_all["col4"] = (
        df_all["chrom"].astype(str)
        + ";"
        + df_all["start"].astype(str)
        + ";"
        + df_all["end"].astype(str)
        + ";"
        + df_all[check_col].astype(str)
    )

    tmp = NamedTemporaryFile(suffix=".all_peaks.bed", delete=False).name
    out = df_all[["chrom", "start", "end", "col4"]].sort_values(["chrom", "start"])
    out.to_csv(tmp, sep="\t", index=False, header=False)

    b = BedTool(tmp)
    all_summits = []
    # loop over merged peaks based on window size and collapse on col4 (summit + value)
    for f in b.slop(b=window // 2, g=genome).merge(c=4, o="collapse"):
        summits = [x.split(";") for x in f[3].split(",")]
        # only keep the highest summit
        all_summits.append(sorted(summits, key=lambda x: float(x[3]))[-1][:3])

    df = pd.DataFrame(all_summits, columns=["chrom", "start", "end"])
    return df
