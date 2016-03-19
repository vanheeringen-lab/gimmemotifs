#!/usr/bin/python
# Copyright (c) 2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
import os
import sys

import pandas as pd
import numpy as np

from gimmemotifs.config import MotifConfig
from gimmemotifs.motif import read_motifs
from gimmemotifs.scanner import Scanner
from gimmemotifs.more import ( KSMoap, LassoMoap, LightningMoap, 
    MaraMoap, MoreMoap, ClassicMoap )

_cluster_methods = ["classic", "ks", "lightning"]
_value_methods = ["lasso", "mara", "more"]

def moap(inputfile, outfile, motiffile=None, pwmfile=None, genome=None, 
        cutoff=0.95, scoring="score", method="classic"):
    
    if scoring not in ['score', 'count']:
        raise ValueError("valid values are 'score' and 'count'")
    
    config = MotifConfig()

    m2f = None
    
    # read data
    df = pd.read_table(inputfile, index_col=0)
   

    if method in _cluster_methods:
        if df.shape[1] != 1:
            raise ValueError("1 column expected for {}".format(method))
    else:
        if np.dtype('object') in set(df.dtypes):
            raise ValueError(
                    "columns should all be numeric for {}".format(method))
        if method not in _value_methods:
            raise ValueError("method {} not valid".format(method))

    if motiffile is None:
        if genome is None:
            raise ValueError("need a genome")
        # check pwmfile
        if pwmfile is None:
            pwmfile = config.get_default_params().get("motif_db", None)
            if pwmfile is not None:
                pwmfile = os.path.join(config.get_motif_dir(), pwmfile)
        
        if pwmfile is None:
            raise ValueError("no pwmfile given and no default database specified")

        if not os.path.exists(pwmfile):
            raise ValueError("{} does not exist".format(pwmfile))

        try:
            motifs = read_motifs(open(pwmfile))
        except:
            sys.stderr.write("can't read motifs from {}".format(pwmfile))
            raise

        base = os.path.splitext(pwmfile)[0]
        map_file = base + ".motif2factors.txt"
        if os.path.exists(map_file):
            m2f = pd.read_table(map_file, index_col=0)

        # initialize scanner
        s = Scanner()
        sys.stderr.write(pwmfile + "\n")
        s.set_motifs(pwmfile)
        s.set_genome(genome)

        # scan for motifs
        sys.stderr.write("scanning for motifs\n")
        motif_names = [m.id for m in read_motifs(open(pwmfile))]
        scores = []
        if method == 'classic' or scoring == "count":
            for row in s.count(list(df.index), cutoff=cutoff):
                scores.append(row)
        else:
            for row in s.best_score(list(df.index)):
                scores.append(row)

        motifs = pd.DataFrame(scores, index=df.index, columns=motif_names)
    else:
        motifs = pd.DataFrame(motiffile)   

    
    clf = None
    if method == "ks":
        clf = KSMoap()
    if method == "lasso":
        clf = LassoMoap()
    if method == "lightning":
        clf = LightningMoap()
    if method == "mara":
        clf = MaraMoap()
    if method == "more":
        clf = MoreMoap()
    if method == "classic":
        clf = ClassicMoap()

    clf.fit(motifs, df)
    clf.act_.to_csv(outfile, sep="\t")

if __name__ == "__main__":
    infile = "/home/simon/git/gimmemotifs/test_rpkm_table.mm10.txt"
    #infile = "/home/simon/git/gimmemotifs/test_clustering.txt"
    moap(infile, "act.out.txt", genome="mm10", method="lasso", scoring="count")

