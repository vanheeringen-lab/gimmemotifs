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


def rankagg(df, method="stuart"):
    """Return aggregated ranks.

    Implementation is ported from the RobustRankAggreg R package

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
    rmat = pd.DataFrame(index=df.iloc[:, 0])

    step = 1 / rmat.shape[0]
    for col in df.columns:
        rmat[col] = pd.DataFrame(
            {col: np.arange(step, 1 + step, step)}, index=df[col]
        ).loc[rmat.index]
    rmat = rmat.apply(sorted, 1, result_type="expand")
    p = rmat.apply(qStuart, 1)
    return pd.DataFrame({"score": p}, index=rmat.index)
