# Copyright (c) 2009-2020 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.

""" Data preprocessing to create GimmeMotifs input. """
import logging
import multiprocessing as mp
import os
from tempfile import NamedTemporaryFile
from typing import Iterable

import genomepy
import numpy as np
import pandas as pd
import pysam
import qnorm
from fluff.track import Track  # noqa: biofluff
from pybedtools import BedTool
from sklearn.preprocessing import scale
from tqdm.auto import tqdm

from gimmemotifs import mytmpdir
from gimmemotifs.utils import determine_file_type

logger = logging.getLogger("gimme.preprocessing")


def coverage_table(
    peakfile,
    datafiles: Iterable,
    window=200,
    log_transform=True,
    normalization=None,
    top=-1,
    topmethod="var",
    rmdup=True,
    rmrepeats=True,
    ncpus=12,
):
    """
    Create a dataframe with peaks from the peakfile as index.
    For each datafile, add a column with reads in these peaks.
    Peaks are normalized to the given window size, and can be normalized,
    transformed and filtered as specified.
    """
    missing = [f for f in datafiles if not os.path.isfile(f)]
    if missing:
        print(f"Could not find {len(missing)} files: {','.join(missing)}")
        raise FileNotFoundError

    for f in datafiles:
        if ".bam" in f and not os.path.isfile(f"{f}.bai"):
            print(
                f"Data file '{f}' does not have an index file. "
                f"Creating an index file..."
            )
            pysam.index(x)  # noqa

    if normalization not in [None, "quantile", "scale"]:
        raise ValueError(f"Unknown method '{normalization}' for normalization")
    if topmethod not in ["var", "std", "mean", "random"]:
        raise ValueError(f"Unknown method '{topmethod}' for selecting regions")

    logger.info("Loading data")
    df = pd.DataFrame()
    peak_bed = peakfile2bedfile(peakfile, window=window)
    regions = _load_peaks(peak_bed)
    if ncpus == 1:
        for datafile in tqdm(datafiles):
            column = _load_reads(peak_bed, datafile, regions, rmdup, rmrepeats)
            df = df.merge(column, left_index=True, right_index=True, how="outer")

    else:
        pool = mp.Pool(processes=ncpus)
        try:
            jobs = []
            for datafile in datafiles:
                jobs.append(
                    pool.apply_async(
                        _load_reads, (peak_bed, datafile, regions, rmdup, rmrepeats)
                    )
                )
            pool.close()
            for job in tqdm(jobs):
                column = job.get()
                df = df.merge(column, left_index=True, right_index=True, how="outer")
            pool.join()
        except Exception as e:
            pool = None  # noqa: force garbage collection on orphaned workers
            raise e

    # normalize & transform
    if log_transform:
        logger.info("Log transform")
        df = np.log1p(df)
    if normalization is None:
        logger.info("No normalization")
    elif normalization == "scale":
        logger.info("Normalization by scaling")
        df[:] = scale(df, axis=0)
    elif normalization == "quantile":
        logger.info("Normalization by quantile normalization")
        # stay between 1-4 ncpus, after 4 highly diminishing returns
        df = qnorm.quantile_normalize(df, ncpus=sorted([1, ncpus, 4])[1])

    # select the top peaks, based on the specified method
    if top > 0:
        idx = None
        if topmethod == "var":
            idx = df.var(1).sort_values().tail(top).index
        elif topmethod == "std":
            idx = df.std(1).sort_values().tail(top).index
        elif topmethod == "mean":
            idx = df.mean(1).sort_values().tail(top).index
        elif topmethod == "random":
            idx = df.sample(top).index
        df = df.loc[idx]

    return df


def peakfile2bedfile(peakfile, bedfile=None, window=200):
    """
    Create a BED file with peak widths normalized to the window size
    """
    if bedfile:
        peak_bed = bedfile
    else:
        # deleted when gimme exits
        peak_bed = os.path.join(mytmpdir(), "peaks.bed")

    with open(peak_bed, "w+") as f:
        half_window = window // 2
        for line in open(peakfile):
            if line.startswith("#") or line[:5] == "track":
                continue
            vals = line.strip().split("\t")

            gene = ""
            score = 0
            strand = "+"
            if len(vals) > 3:
                gene = vals[3]
                if len(vals) > 4:
                    score = vals[4]
                    if len(vals) > 5:
                        strand = vals[5]

            center = (int(vals[2]) + int(vals[1])) // 2
            if strand == "+":
                start = center - half_window
                end = center + half_window
            else:
                start = center + half_window
                end = center - half_window

            if start >= 0:
                f.write(f"{vals[0]}\t{start}\t{end}\t{gene}\t{score}\t{strand}\n")

    return peak_bed


def _load_peaks(peak_bed):
    """read a bed file and return regions in 'chr:start-end' format"""
    regions = []
    with open(peak_bed) as f:
        for line in f:
            vals = line.split("\t")
            region = f"{vals[0]}:{vals[1]}-{vals[2]}"
            regions.append(region)
    return regions


def _load_reads(peak_bed, datafile, regions, rmdup=True, rmrepeats=True):
    """
    read a bam, bed or bw file and return the number of reads per peak
    """
    track = Track.load(
        datafile,
        rmdup=rmdup,
        rmrepeats=rmrepeats,
        fragmentsize=None,
    )
    result = track.binned_stats(peak_bed, nbins=1, split=True, rpkm=False)

    reads = [float(r[3]) for r in result]
    col = [os.path.splitext(os.path.basename(datafile))[0]]
    df = pd.DataFrame(reads, columns=col, index=regions)
    return df


def combine_peaks(peaks, genome, window=200, scale_value=False):
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
    except FileNotFoundError:
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

    df["log_value"] = np.log1p(df["value"])
    df["log_value_scaled"] = scale(df[["log_value"]])
    return df
