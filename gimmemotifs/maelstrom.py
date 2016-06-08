#!/usr/bin/python
# Copyright (c) 2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
import os
import subprocess as sp
import sys
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd

from gimmemotifs.background import RandomGenomicFasta
from gimmemotifs.config import MotifConfig
from gimmemotifs.moap import moap
from gimmemotifs.rank import rankagg

BG_LENGTH = 200
BG_NUMBER = 10000
FDR = 0.01

def moap_with_bg(input_table, genome, data_dir, method, scoring):
    # gimme_motifs config, to get defaults
    config = MotifConfig()
    
    threshold_file = None
    if scoring == "count":
        # Motif scanning threshold
        threshold_file = os.path.join(data_dir, "threshold.{}.txt".format(genome))
        if not os.path.exists(threshold_file):
        # Random sequences from genome
            index_dir = os.path.join(config.get_index_dir(), genome)
            bg_file = os.path.join(data_dir, "background.{}.fa".format(genome))
            if not os.path.exists(bg_file):
                m = RandomGenomicFasta(index_dir, BG_LENGTH, BG_NUMBER)
                m.writefasta(bg_file)
    
            pwmfile = config.get_default_params().get("motif_db")
            pwmfile = os.path.join(config.get_motif_dir(), pwmfile)
            
            cmd = "gimme threshold {} {} {} > {}".format(
                    pwmfile,
                    bg_file,
                    FDR,
                    threshold_file)
            sp.call(cmd, shell=True)
    
    outfile = os.path.join(data_dir,"activity.{}.{}.out.txt".format(
            method,
            scoring))

    moap(input_table, outfile=outfile, genome=genome, method=method,
            scoring=scoring, cutoff=threshold_file)

def moap_with_table(input_table, motif_table, data_dir, method, scoring):
    outfile = os.path.join(data_dir,"activity.{}.{}.out.txt".format(
            method,
            scoring))

    moap(input_table, outfile=outfile, method=method, scoring=scoring, 
            motiffile=motif_table)

def run_maelstrom(infile, genome, outdir, cluster=True, 
        score_table=None, count_table=None):

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    df = pd.read_table(infile, index_col=0)

    # Drop duplicate indices, doesn't work very well downstream
    df = df.loc[df.index.drop_duplicates(keep=False)]
    exps = []
    clusterfile = infile
    if df.shape[1] != 1:
        # More than one column
        exps += [
                ("mara", "count", infile),
                ("lasso", "score", infile),
                ]

        if cluster:
            clusterfile = os.path.join(outdir,
                    os.path.basename(infile) + ".cluster.txt")
            df = df.apply(scale, 0)
            names = df.columns
            df_changed = pd.DataFrame(index=df.index)
            df_changed["cluster"] = np.nan
            for name in names:
                df_changed.loc[(df[name] - df.loc[:,df.columns != name].max(1)) > 0.5, "cluster"] = name
            df_changed.dropna().to_csv(clusterfile, sep="\t")
    if df.shape[1] == 1 or cluster:
        exps += [
                ("rf", "score", clusterfile),
                ("classic", "count", clusterfile),
                ("ks", "score", clusterfile),
                ("lightning", "score", clusterfile),
                ]

    for method, scoring, fname in exps:
        try:
            sys.stderr.write("Running {} with {}\n".format(method,scoring))
            if scoring == "count" and count_table:
                moap_with_table(fname, count_table, outdir, method, scoring)
            elif scoring == "score" and score_table:
                moap_with_table(fname, score_table, outdir, method, scoring)
            else:
                moap_with_bg(fname, genome, outdir, method, scoring)
        
        
        except Exception as e:
            sys.stderr.write(
                    "Method {} with scoring {} failed\n{}\nSkipping\n".format(
                        method, scoring, e)
                    )
    
    dfs = {}
    for method, scoring,fname  in exps:
        t = "{}.{}".format(method,scoring)
        fname = os.path.join(outdir, "activity.{}.{}.out.txt".format(
                           method, scoring))
        dfs[t] = pd.read_table(fname, index_col=0, comment="#")
    
    
    df_p = pd.DataFrame(index=dfs.values()[0].index)
    names = dfs.values()[0].columns
    for e in names:
        df_tmp = pd.DataFrame()
        for method,scoring,fname in exps:
            k = "{}.{}".format(method, scoring)
            v = dfs[k]
            df_tmp[k] = v.sort_values(e, ascending=False).index.values
        
        df_p[e] = rankagg(df_tmp)
    df_p[names] = -np.log10(df_p[names])
    df_p.to_csv(os.path.join(outdir, "final.out.csv"), sep="\t")
    #df_p = df_p.join(m2f)
