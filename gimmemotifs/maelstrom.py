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
from sklearn.preprocessing import scale
from scipy.cluster import hierarchy
from scipy.spatial import distance

# Plotting
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')

from gimmemotifs.background import RandomGenomicFasta
from gimmemotifs.config import MotifConfig
from gimmemotifs.moap import moap
from gimmemotifs.rank import rankagg
from gimmemotifs.motif import read_motifs
from gimmemotifs.scanner import Scanner

BG_LENGTH = 200
BG_NUMBER = 10000
FDR = 0.01

def check_threshold(outdir, genome, scoring="count", pwmfile=None):
    # gimme_motifs config, to get defaults
    config = MotifConfig()
    
    threshold_file = None
    if scoring == "count":
        # Motif scanning threshold
        threshold_file = os.path.join(outdir, "threshold.{}.txt".format(genome))
        if not os.path.exists(threshold_file):
        # Random sequences from genome
            index_dir = os.path.join(config.get_index_dir(), genome)
            bg_file = os.path.join(outdir, "background.{}.fa".format(genome))
            if not os.path.exists(bg_file):
                m = RandomGenomicFasta(index_dir, BG_LENGTH, BG_NUMBER)
                m.writefasta(bg_file)
    
            if pwmfile is None:
                pwmfile = config.get_default_params().get("motif_db")
                pwmfile = os.path.join(config.get_motif_dir(), pwmfile)
            
            cmd = "gimme threshold {} {} {} > {}".format(
                    pwmfile,
                    bg_file,
                    FDR,
                    threshold_file)
            sp.call(cmd, shell=True)
        return threshold_file

def scan_to_table(input_table, genome, data_dir, scoring, pwmfile=None):
    threshold = check_threshold(data_dir, genome, scoring, pwmfile)
    
    config = MotifConfig()
    
    if pwmfile is None:
        pwmfile = config.get_default_params().get("motif_db", None)
        if pwmfile is not None:
            pwmfile = os.path.join(config.get_motif_dir(), pwmfile)

    if pwmfile is None:
        raise ValueError("no pwmfile given and no default database specified")

    df = pd.read_table(input_table, index_col=0)
    regions = list(df.index)
    s = Scanner()
    s.set_motifs(pwmfile)
    s.set_genome(genome)

    scores = []
    if scoring == "count":
        for row in s.count(regions, cutoff=threshold):
            scores.append(row)
    else:
        for row in s.best_score(regions):
            scores.append(row)
   
    motif_names = [m.id for m in read_motifs(open(pwmfile))]
    return pd.DataFrame(scores, index=df.index, columns=motif_names)

def moap_with_bg(input_table, genome, data_dir, method, scoring, pwmfile=None):
    threshold_file = check_threshold(data_dir, genome, scoring, pwmfile)
    
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

def safe_join(df1, df2):     
    tmp = df1.copy()
    tmp["_safe_count"] = range(df1.shape[0])     
    return tmp.join(df2).sort_values("_safe_count").drop( "_safe_count", 1)

def visualize_maelstrom(outdir, sig_cutoff=3, pwmfile=None):
    
    config = MotifConfig()
    if pwmfile is None:
        pwmfile = config.get_default_params().get("motif_db", None)
        pwmfile = os.path.join(config.get_motif_dir(), pwmfile)
    
    mapfile = pwmfile.replace(".pwm", ".motif2factors.txt")
    if os.path.exists(mapfile):
    
        m2f = pd.read_csv(mapfile, sep="\t", names=["motif","factors"], index_col=0) 
        m2f["factors"] = m2f["factors"].str[:50]
    else:
        motifs = [m.id for m in read_motifs(open(pwmfile))]
        m2f = pd.DataFrame({"factors": motifs}, index=motifs)

    sig_fname = os.path.join(outdir, "final.out.csv")
    df_sig = pd.read_table(sig_fname, index_col=0)
    f = np.any(df_sig >= sig_cutoff, 1)
    vis = df_sig[f]
    if vis.shape[0] == 0:
        sys.stderr.write("No motifs reach the threshold, skipping visualization.\n")
        return
    
    # cluster rows
    row_linkage = hierarchy.linkage(
        distance.pdist(vis, metric="euclidean"), 
        method='complete')
    idx = hierarchy.leaves_list(row_linkage)
    
    plt.figure(figsize=size)
    
    vis = safe_join(vis, m2f).set_index("factors")
    
    # size of figure
    size = [2 + vis.shape[1] * 0.4, 1.8 + vis.shape[0] * 0.3]
    
    cg = sns.heatmap(vis.iloc[idx], cmap="viridis", 
                        yticklabels=True, 
                       cbar_kws={"orientation":"horizontal"})
    _ = plt.setp(cg.yaxis.get_majorticklabels(), rotation=0)
    plt.title("Motif Relevance")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "motif.relevance.png"), dpi=300) 
   
    freq_fname = os.path.join(outdir, "motif.freq.txt")
    if os.path.exists(freq_fname):
        df_freq = pd.read_table(freq_fname, index_col=0) 
        df_freq = df_freq.T
        vis_freq = df_freq.loc[vis.iloc[idx].index]
        vis_freq = safe_join(vis_freq, m2f).set_index("factors")
        plt.figure(figsize=size)
        cg = sns.heatmap(vis_freq, cmap="viridis", 
                         yticklabels=True, vmin=0, vmax=0.2,
                           cbar_kws={"orientation":"horizontal"})
        #idx = cg.dendrogram_row.reordered_ind
        _ = plt.setp(cg.yaxis.get_majorticklabels(), rotation=0)
        plt.title("Motif Frequency")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "motif.frequency.png"), dpi=300) 
        
        plt.figure(figsize=size)
        
        bla = vis_freq.min(1)
        bla[bla < 0.01] = 0.01
        
        cg = sns.heatmap(np.log2(vis_freq.apply(lambda x: x / bla, 0)), 
                         yticklabels=True, vmin=-5, vmax=5,
                        cbar_kws={"orientation":"horizontal"})
        #idx = cg.dendrogram_row.reordered_ind
        _ = plt.setp(cg.yaxis.get_majorticklabels(), rotation=0)
        plt.title("Motif Enrichment")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "motif.enrichment.png"), dpi=300) 

def run_maelstrom(infile, genome, outdir, pwmfile=None, plot=True, cluster=True, 
        score_table=None, count_table=None):

    df = pd.read_table(infile, index_col=0)
    # Check for duplicates
    if df.index.duplicated(keep=False).any():
        raise ValueError("Input file contains duplicate regions! "
                         "Please remove them.")
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    # Create a file with the number of motif matches
    if not count_table:
        counts = scan_to_table(infile, genome, outdir, "count",
                pwmfile=pwmfile)
        count_table = os.path.join(outdir, "motif.count.txt.gz")
        counts.to_csv(count_table, sep="\t", compression="gzip")

    # Create a file with the score of the best motif match
    if not score_table:
        scores = scan_to_table(infile, genome, outdir, "score",
                pwmfile=pwmfile)
        score_table = os.path.join(outdir, "motif.score.txt.gz")
        scores.to_csv(score_table, sep="\t", float_format="%.3f", 
                compression="gzip")

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
                ("mwu", "score", clusterfile),
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
                moap_with_bg(fname, genome, outdir, method, scoring, pwmfile=pwmfile)
        
        
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
        try:
            dfs[t] = pd.read_table(fname, index_col=0, comment="#")
        except:
            sys.stderr.write("Activity file for {} not found!\n".format(t))
    
    df_p = pd.DataFrame(index=dfs.values()[0].index)
    names = dfs.values()[0].columns
    for e in names:
        df_tmp = pd.DataFrame()
        for method,scoring,fname in exps:
            k = "{}.{}".format(method, scoring)
            if k in dfs:
                v = dfs[k]
                df_tmp[k] = v.sort_values(e, ascending=False).index.values
        
        df_p[e] = rankagg(df_tmp)
    df_p[names] = -np.log10(df_p[names])
    df_p.to_csv(os.path.join(outdir, "final.out.csv"), sep="\t")
    #df_p = df_p.join(m2f)

    # Write motif frequency table
    
    if df.shape[1] == 1:
        mcount = df.join(pd.read_table(count_table, index_col=0))
        m_group = mcount.groupby(df.columns[0])
        freq = (m_group.sum() / m_group.count())
        freq.to_csv(os.path.join(outdir, "motif.freq.txt"), sep="\t")

    if plot:
        visualize_maelstrom(outdir, pwmfile=pwmfile)

