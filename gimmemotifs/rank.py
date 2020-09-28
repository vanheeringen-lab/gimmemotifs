#!/usr/bin/python
# Copyright (c) 2016-2019 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
"""Rank aggregation (includes wrapper for R RobustRankAgg)."""
from tempfile import NamedTemporaryFile
import subprocess as sp
import pandas as pd
import numpy as np
from scipy.stats import rankdata, norm

try:
    from scipy.special import factorial
except ImportError:
    from scipy.misc import factorial


def rankagg_R(df, method="stuart"):
    """Return aggregated ranks as implemented in the RobustRankAgg R package.

    This function is now deprecated.

    References:
        Kolde et al., 2012, DOI: 10.1093/bioinformatics/btr709
        Stuart et al., 2003,  DOI: 10.1126/science.1087447

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with values to be ranked and aggregated

    Returns
    -------
    pandas.DataFrame with aggregated ranks
    """
    tmpdf = NamedTemporaryFile()
    tmpscript = NamedTemporaryFile(mode="w")
    tmpranks = NamedTemporaryFile()

    df.to_csv(tmpdf.name, sep="\t", index=False)

    script = """
library(RobustRankAggreg);
a = read.table("{}", header=TRUE);
x = lapply(a, as.vector);
result = aggregateRanks(x, method="{}");
result$p.adjust = p.adjust(result$Score);
 write.table(result, file="{}", sep="\t", quote=FALSE, row.names=FALSE);
""".format(
        tmpdf.name, method, tmpranks.name
    )
    tmpscript.write(script)
    tmpscript.flush()

    p = sp.Popen(["Rscript", tmpscript.name], stdout=sp.PIPE, stderr=sp.PIPE)
    stderr, stdout = p.communicate()
    df = pd.read_table(tmpranks.name, index_col=0)
    return df["p.adjust"]


def sumStuart(v, r):
    k = len(v)
    l_k = np.arange(k)
    ones = (-1) ** l_k
    f = factorial(l_k + 1)
    p = r ** (l_k + 1)
    return np.dot(ones, v[::-1] * p / f)


def qStuart(r):
    N = (~r.isnull()).sum().sum()
    v = np.ones(N + 1)
    for k in range(N):
        v[k + 1] = sumStuart(v[: k + 1], r[N - k - 1])

    return factorial(N) * v[N]


def _rank_int(series, c=3.0 / 8, stochastic=True):
    # Based on code by Edward Mountjoy
    # See: https://github.com/edm1/rank-based-INT
    """Perform rank-based inverse normal transformation on pandas series.
    If stochastic is True ties are given rank randomly, otherwise ties will
    share the same value. NaN values are ignored.
    Args:
        param1 (pandas.Series):   Series of values to transform
        param2 (Optional[float]): Constand parameter (Bloms constant)
        param3 (Optional[bool]):  Whether to randomise rank of ties

    Returns:
        pandas.Series
    """

    # Check input
    assert isinstance(series, pd.Series)
    assert isinstance(c, float)
    assert isinstance(stochastic, bool)

    # Set seed
    np.random.seed(123)

    # Take original series indexes
    orig_idx = series.index

    # Drop NaNs
    series = series.loc[~pd.isnull(series)]

    # Get ranks
    if stochastic:
        # Shuffle by index
        series = series.loc[np.random.permutation(series.index)]
        # Get rank, ties are determined by their position in the series (hence
        # why we randomised the series)
        rank = rankdata(series, method="ordinal")
    else:
        # Get rank, ties are averaged
        rank = rankdata(series, method="average")

    # Convert numpy array back to series
    rank = pd.Series(rank, index=series.index)

    # Convert rank to normal distribution
    transformed = rank.apply(_rank_to_normal, c=c, n=len(rank))

    return transformed[orig_idx]


def _rank_to_normal(rank, c, n):
    # Standard quantile function
    x = (rank - c) / (n - 2 * c + 1)
    return norm.ppf(x)


def _rankagg_int(df):
    # Convert values to ranks
    df_int = df.apply(_rank_int)
    # Combine z-score using Stouffer's method
    df_int = (df_int.sum(1) / np.sqrt(df_int.shape[1])).to_frame()
    df_int.columns = ["z-score"]
    return df_int


def _rankagg_stuart(df):
    """
    Implementation is ported from the RobustRankAggreg R package

    References:
        Kolde et al., 2012, DOI: 10.1093/bioinformatics/btr709
        Stuart et al., 2003,  DOI: 10.1126/science.1087447
    """
    rmat = pd.DataFrame(index=df.iloc[:, 0])

    step = 1 / rmat.shape[0]
    for col in df.columns:
        rmat[col] = pd.DataFrame(
            {col: np.arange(step, 1 + step, step)}, index=df[col]
        ).loc[rmat.index]
    rmat = rmat.apply(sorted, 1, result_type="expand")
    p = rmat.apply(qStuart, 1)
    return pd.DataFrame({"score": p}, index=rmat.index)


def rankagg(df, method="int_stouffer", include_reverse=True, log_transform=True):
    """Return aggregated ranks.

    Stuart implementation is ported from the RobustRankAggreg R package

    References:
        Kolde et al., 2012, DOI: 10.1093/bioinformatics/btr709
        Stuart et al., 2003,  DOI: 10.1126/science.1087447

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame with values to be ranked and aggregated
    method : str, optional
        Either "int_stouffer" or "stuart". The "int_stouffer" method is based on combining z-scores
        from a inverse normal transform of ranks using Stouffer's method.

    Returns
    -------
    pandas.DataFrame with aggregated ranks
    """
    method = method.lower()
    if method not in ["stuart", "int_stouffer"]:
        raise ValueError("unknown method for rank aggregation")

    if method == "stuart":
        df_asc = pd.DataFrame()
        df_desc = pd.DataFrame()
        for col in df.columns:
            df_asc[col] = (
                df.sample(frac=1).sort_values(col, ascending=False).index.values
            )
            if include_reverse:
                df_desc[col] = (
                    df.sample(frac=1).sort_values(col, ascending=True).index.values
                )

        df_result = -np.log10(_rankagg_stuart(df_asc))
        if include_reverse:
            df_result += np.log10(_rankagg_stuart(df_desc))

        return df_result
    if method == "int_stouffer":
        return _rankagg_int(df)
