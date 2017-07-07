#!/usr/bin/python
# Copyright (c) 2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
"""Rank aggregation (wrapper for R RobustRankAgg)."""
from tempfile import NamedTemporaryFile
import subprocess as sp
import pandas as pd

def rankagg(df, method="stuart"):
    """Return aggregated ranks as implemented in the RobustRankAgg R package.

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

    df.to_csv(tmpdf.name, sep="\t",index=False) 
 
    script = ''' 
library(RobustRankAggreg); 
a = read.table("{}", header=TRUE); 
x = lapply(a, as.vector); 
result = aggregateRanks(x, method="{}"); 
result$p.adjust = p.adjust(result$Score); 
 write.table(result, file="{}", sep="\t", quote=FALSE, row.names=FALSE); 
'''.format(tmpdf.name, method, tmpranks.name) 
    tmpscript.write(script) 
    tmpscript.flush() 
 
    p = sp.Popen(["Rscript", tmpscript.name], stdout=sp.PIPE, stderr=sp.PIPE)
    stderr, stdout = p.communicate()
    df = pd.read_table(tmpranks.name, index_col=0)
    return df["p.adjust"] 
