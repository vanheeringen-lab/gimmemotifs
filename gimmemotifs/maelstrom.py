# Copyright (c) 2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
import glob
import os
import subprocess as sp
import shutil
import sys
from tempfile import NamedTemporaryFile
import logging
from functools import partial

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
from gimmemotifs.moap import moap, Moap
from gimmemotifs.rank import rankagg
from gimmemotifs.motif import read_motifs
from gimmemotifs.scanner import Scanner
from gimmemotifs.report import maelstrom_html_report
from gimmemotifs.utils import join_max

from multiprocessing import Pool

BG_LENGTH = 200
BG_NUMBER = 10000
FPR = 0.01

logger = logging.getLogger("gimme.maelstrom")

def scan_to_table(input_table, genome, data_dir, scoring, pwmfile=None, ncpus=None):
    config = MotifConfig()
    
    if pwmfile is None:
        pwmfile = config.get_default_params().get("motif_db", None)
        if pwmfile is not None:
            pwmfile = os.path.join(config.get_motif_dir(), pwmfile)

    if pwmfile is None:
        raise ValueError("no pwmfile given and no default database specified")

    logger.info("reading table")
    if input_table.endswith("feather"):
        df = pd.read_feather(input_table)
        idx = df.iloc[:,0].values
    else:
        df = pd.read_table(input_table, index_col=0, comment="#")
        idx = df.index
    
    regions = list(idx)
    s = Scanner(ncpus=ncpus)
    s.set_motifs(pwmfile)
    s.set_genome(genome)
    nregions = len(regions)

    scores = []
    if scoring == "count":
        logger.info("setting threshold")
        s.set_threshold(fpr=FPR, genome=genome)
        logger.info("creating count table")
        for row in s.count(regions):
            scores.append(row)
        logger.info("done")
    else:
        s.set_threshold(threshold=0.0)
        logger.info("creating score table")
        for row in s.best_score(regions):
            scores.append(row)
        logger.info("done")
   
    motif_names = [m.id for m in read_motifs(open(pwmfile))]
    logger.info("creating dataframe")
    return pd.DataFrame(scores, index=idx, columns=motif_names)

def moap_with_bg(input_table, genome, data_dir, method, scoring, pwmfile=None, ncpus=None):
    outfile = os.path.join(data_dir,"activity.{}.{}.out.txt".format(
            method,
            scoring))

    moap(input_table, outfile=outfile, genome=genome, method=method,
            scoring=scoring, fpr=FPR, ncpus=ncpus)

def moap_with_table(input_table, motif_table, data_dir, method, scoring, ncpus=None):
    outfile = os.path.join(data_dir,"activity.{}.{}.out.txt".format(
            method,
            scoring))

    moap(input_table, outfile=outfile, method=method, scoring=scoring, 
            motiffile=motif_table, fpr=FPR, ncpus=ncpus)

def safe_join(df1, df2):     
    tmp = df1.copy()
    tmp["_safe_count"] = list(range(df1.shape[0]))
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
    df_sig = pd.read_table(sig_fname, index_col=0, comment="#")
    f = np.any(df_sig >= sig_cutoff, 1)
    vis = df_sig[f]
    if vis.shape[0] == 0:
        logging.info("No motifs reach the threshold, skipping visualization.\n")
        return
    
    # cluster rows
    row_linkage = hierarchy.linkage(
        distance.pdist(vis, metric="euclidean"), 
        method='complete')
    idx = hierarchy.leaves_list(row_linkage)
    
    plt.figure()
    
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
        df_freq = pd.read_table(freq_fname, index_col=0, comment="#") 
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

def _rank_agg_column(exps, dfs, e):
    tmp_dfs = [pd.DataFrame(), pd.DataFrame()]
        
    for i,sort_order in enumerate([False, True]):
        for method,scoring,fname in exps:
            k = "{}.{}".format(method, scoring)
            if k in dfs:
                v = dfs[k]
                tmp_dfs[i][k] = v.sort_values(e, ascending=sort_order).index.values
    return -np.log10(rankagg(tmp_dfs[0])) + np.log10(rankagg(tmp_dfs[1]))

def df_rank_aggregation(df, dfs, exps):
    df_p = pd.DataFrame(index=list(dfs.values())[0].index)
    df_negp = pd.DataFrame(index=list(dfs.values())[0].index)
    names = list(dfs.values())[0].columns
    pool = Pool(16)
    func = partial(_rank_agg_column, exps, dfs)
    ret = pool.map(func, names)
    pool.close()
    pool.join()
    
    for e, result in zip(names, ret):
        df_p[e] = result
        
    if df.shape[1] != 1:
        df_p = df_p[df.columns]

    return df_p

def run_maelstrom(infile, genome, outdir, pwmfile=None, plot=True, cluster=True, 
        score_table=None, count_table=None, methods=None, ncpus=None):
    logger.info("Starting maelstrom")
    if infile.endswith("feather"):
        df = pd.read_feather(infile)
        df = df.set_index(df.columns[0])
    else:
        df = pd.read_table(infile, index_col=0, comment="#")
    
    # Check for duplicates
    if df.index.duplicated(keep=False).any():
        raise ValueError("Input file contains duplicate regions! "
                         "Please remove them.")
    
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if methods is None:
        methods = Moap.list_predictors() 
    methods = [m.lower() for m in methods]

    shutil.copyfile(infile, os.path.join(outdir, "input.table.txt"))
    
    config = MotifConfig()
    # Default pwmfile
    if pwmfile is None:
        pwmfile = config.get_default_params().get("motif_db", None)
        if pwmfile is not None:
            pwmfile = os.path.join(config.get_motif_dir(), pwmfile)

    if pwmfile:
        shutil.copy2(pwmfile, outdir)
        mapfile = pwmfile.replace(".pwm", ".motif2factors.txt")
        if os.path.exists(mapfile):
            shutil.copy2(mapfile, outdir)
    
    # Create a file with the number of motif matches
    if not count_table:
        count_table = os.path.join(outdir, "motif.count.txt.gz")
        if not os.path.exists(count_table):
            logger.info("Motif scanning (counts)")
            counts = scan_to_table(infile, genome, outdir, "count",
                pwmfile=pwmfile, ncpus=ncpus)
            counts.to_csv(count_table, sep="\t", compression="gzip")
        else:
            logger.info("Counts, using: %s", count_table)

    # Create a file with the score of the best motif match
    if not score_table:
        score_table = os.path.join(outdir, "motif.score.txt.gz")
        if not os.path.exists(score_table):
            logger.info("Motif scanning (scores)")
            scores = scan_to_table(infile, genome, outdir, "score",
                pwmfile=pwmfile, ncpus=ncpus)
            scores.to_csv(score_table, sep="\t", float_format="%.3f", 
                compression="gzip")
        else:
            logger.info("Scores, using: %s", score_table)

    if cluster:
        cluster = False
        for method in methods:
            m = Moap.create(method, ncpus=ncpus)
            if m.ptype == "classification":
                cluster = True
                break
        if not cluster:
            logger.info("Skipping clustering, no classification methods")
    
    exps = []
    clusterfile = infile
    if df.shape[1] != 1:
        # More than one column
        for method in Moap.list_regression_predictors():
            if method in methods:
                m = Moap.create(method, ncpus=ncpus)
                exps.append([method, m.pref_table, infile])
                logger.debug("Adding %s", method)

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
        for method in Moap.list_classification_predictors():
            if method in methods:
                m = Moap.create(method, ncpus=ncpus)
                exps.append([method, m.pref_table, clusterfile])

    if len(exps) == 0:
        logger.error("No method to run.")
        sys.exit(1)

    for method, scoring, fname in exps:
        try:
            if scoring == "count" and count_table:
                moap_with_table(fname, count_table, outdir, method, scoring, ncpus=ncpus)
            elif scoring == "score" and score_table:
                moap_with_table(fname, score_table, outdir, method, scoring, ncpus=ncpus)
            else:
                moap_with_bg(fname, genome, outdir, method, scoring, pwmfile=pwmfile, ncpus=ncpus)
        
        except Exception as e:
            logger.warn("Method %s with scoring %s failed", method, scoring)
            logger.warn(e)
            logger.warn("Skipping")
            raise 
    dfs = {}
    for method, scoring,fname  in exps:
        t = "{}.{}".format(method,scoring)
        fname = os.path.join(outdir, "activity.{}.{}.out.txt".format(
                           method, scoring))
        try:
            dfs[t] = pd.read_table(fname, index_col=0, comment="#")
        except:
            logging.warn("Activity file for {} not found!\n".format(t))
   
    if len(methods) > 1:
        logger.info("Rank aggregation")
        df_p = df_rank_aggregation(df, dfs, exps)
        df_p.to_csv(os.path.join(outdir, "final.out.csv"), sep="\t")
    #df_p = df_p.join(m2f)

    # Write motif frequency table
    
    if df.shape[1] == 1:
        mcount = df.join(pd.read_table(count_table, index_col=0, comment="#"))
        m_group = mcount.groupby(df.columns[0])
        freq = m_group.sum() / m_group.count()
        freq.to_csv(os.path.join(outdir, "motif.freq.txt"), sep="\t")

    if plot and len(methods) > 1:
        logger.info("html report")
        maelstrom_html_report(
                outdir, 
                os.path.join(outdir, "final.out.csv"),
                pwmfile
                )
        logger.info(os.path.join(outdir, "gimme.maelstrom.report.html"))

class MaelstromResult():
    """Class for working with maelstrom output."""
    
    def __init__(self, outdir):
        """Initialize a MaelstromResult object from a maelstrom output 
        directory.
        
        Parameters
        ----------
        outdir : str
            Name of a maelstrom output directory.
        
        See Also
        --------
        maelstrom.run_maelstrom : Run a maelstrom analysis.
        """    
        if not os.path.exists(outdir):
            raise FileNotFoundError("No such directory: " + outdir)
            
        # Load motifs
        pwmfile = glob.glob(os.path.join(outdir, "*.pwm"))[0]
        with open(pwmfile) as fin:
            self.motifs = {m.id: m for m in read_motifs(fin)}
        
        self.activity = {}
        # Read individual activity files
        for fname in glob.glob(os.path.join(outdir, "activity*txt")):
            #print()
            _, name, mtype, _, _ = os.path.split(fname)[-1].split(".")
            self.activity["{}.{}".format(name, mtype)] = pd.read_table(
                fname, comment="#", index_col=0)
        
        # Read rank aggregation
        self.result = pd.read_table(
            os.path.join(outdir, "final.out.csv"), 
            comment="#", index_col=0
        )
        
        # Read motif results
        self.scores = pd.read_table(
            os.path.join(outdir, "motif.score.txt.gz"), 
            index_col=0
        )
        self.counts = pd.read_table(
            os.path.join(outdir, "motif.count.txt.gz"), 
            index_col=0
        )
        fname = os.path.join(outdir, "motif.freq.txt")
        if os.path.exists(fname):
            self.freq = pd.read_table(
                fname, 
                index_col=0
            )
        
        # Read original input file
        try:
            self.input = pd.read_table(
                os.path.join(outdir, "input.table.txt"), 
                index_col=0
            )
            if self.input.shape[1] == 1:
                self.input.columns = ["cluster"]
        except:
            pass
    
    def plot_heatmap(self, kind="final", min_freq=0.01, threshold=2, name=True, max_len=50, aspect=1, **kwargs):
        """Plot clustered heatmap of predicted motif activity.
        
        Parameters
        ----------
        kind : str, optional
            Which data type to use for plotting. Default is 'final', which will plot the 
            result of the rang aggregation. Other options are 'freq' for the motif frequencies,
            or any of the individual activities such as 'rf.score'.
            
        min_freq : float, optional
            Minimum frequency of motif occurrence.
            
        threshold : float, optional
            Minimum activity (absolute) of the rank aggregation result. 
        
        name : bool, optional
            Use factor names instead of motif names for plotting.
        
        max_len : int, optional
            Truncate the list of factors to this maximum length.
            
        aspect : int, optional
            Aspect ratio for tweaking the plot.
            
        kwargs : other keyword arguments
            All other keyword arguments are passed to sns.clustermap

        Returns
        -------
        cg : ClusterGrid
            A seaborn ClusterGrid instance.
        """
        
        filt = np.any(np.abs(self.result) >= threshold, 1) & np.any(np.abs(self.freq.T) >= min_freq, 1)
        
        idx = self.result[filt].index
        
        cmap = "RdBu_r" 
        if kind == "final":
            data = self.result
        elif kind == "freq":
            data = self.freq.T
            cmap = "Reds"
        elif kind in self.activity:
            data = self.activity[dtype]
            if kind in ["hypergeom.count", "mwu.score"]:
                cmap = "Reds"
        else:
            raise ValueError("Unknown dtype")
        
        #print(data.head())
        #plt.figure(
        m = data.loc[idx]
        if name:
            m["factors"] = [join_max(self.motifs[n].factors, max_len, ",", suffix=",(...)") for n in m.index]
            m = m.set_index("factors")
        h,w = m.shape
        cg = sns.clustermap(m, cmap=cmap, col_cluster=False, 
                            figsize=(2 + w * 0.5 * aspect, 0.5 * h), linewidths=1,
                           **kwargs)
        cg.ax_col_dendrogram.set_visible(False)
        plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0);
        return cg
        
        
    def plot_scores(self, motifs, name=True, max_len=50):
        """Create motif scores boxplot of different clusters.
        Motifs can be specified as either motif or factor names.
        The motif scores will be scaled and plotted as z-scores.
        
        Parameters
        ----------
        motifs : iterable or str
            List of motif or factor names.
        
        name : bool, optional
            Use factor names instead of motif names for plotting.
        
        max_len : int, optional
            Truncate the list of factors to this maximum length.
        
        Returns
        -------
        
        g : FacetGrid
            Returns the seaborn FacetGrid object with the plot.
        """
        if self.input.shape[1] != 1:
            raise ValueError("Can't make a categorical plot with real-valued data")
        
        if type("") == type(motifs):
            motifs = [motifs]
            
        plot_motifs = []
        for motif in motifs:
            if motif in self.motifs:
                plot_motifs.append(motif)
            else:
                for m in self.motifs.values():
                    if motif in m.factors:
                        plot_motifs.append(m.id)
        
        data = self.scores[plot_motifs].apply(scale, axis=0)
        if name:
            data = data.T
            data["factors"] = [join_max(self.motifs[n].factors, max_len, ",", suffix=",(...)") for n in plot_motifs]
            data = data.set_index("factors").T
        
        data = pd.melt(self.input.join(data), id_vars=["cluster"])
        data.columns = ["cluster", "motif", "z-score"]
        g = sns.factorplot(data=data, y="motif", x="z-score", hue="cluster", kind="box", aspect=2)
        return g
