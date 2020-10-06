# Copyright (c) 2016-2019 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
"""Class implementing the maelstrom method (Bruse & van Heeringen, 2018)

Examples
--------

run_maelstrom(input, "hg38", outdir)
mr = MaelstromResult(outdir)
"""
import glob
import os
import re
import shutil
import sys
import logging
from functools import partial

import numpy as np
import pandas as pd
from sklearn.preprocessing import scale
from scipy.stats import pearsonr
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.cluster import FeatureAgglomeration

# from scipy.spatial.distance import correlation

# Plotting
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns

sns.set_style("white")

from gimmemotifs.config import MotifConfig, DIRECT_NAME, INDIRECT_NAME
from gimmemotifs.moap import moap, Moap
from gimmemotifs.scanner import scan_regionfile_to_table
from gimmemotifs.rank import rankagg
from gimmemotifs.motif import read_motifs
from gimmemotifs.report import maelstrom_html_report
from gimmemotifs.utils import join_max, pfmfile_location

from multiprocessing import Pool

BG_LENGTH = 200
BG_NUMBER = 10000
FPR = 0.01

logger = logging.getLogger("gimme.maelstrom")


def moap_with_bg(
    input_table, genome, data_dir, method, scoring, pfmfile=None, ncpus=None
):
    outfile = os.path.join(data_dir, "activity.{}.{}.out.txt".format(method, scoring))

    moap(
        input_table,
        outfile=outfile,
        pfmfile=pfmfile,
        genome=genome,
        method=method,
        scoring=scoring,
        fpr=FPR,
        ncpus=ncpus,
    )


def moap_with_table(input_table, motif_table, data_dir, method, scoring, ncpus=None):
    outfile = os.path.join(data_dir, "activity.{}.{}.out.txt".format(method, scoring))

    moap(
        input_table,
        outfile=outfile,
        method=method,
        scoring=scoring,
        motiffile=motif_table,
        fpr=FPR,
        ncpus=ncpus,
    )


def safe_join(df1, df2):
    tmp = df1.copy()
    tmp["_safe_count"] = list(range(df1.shape[0]))
    return tmp.join(df2).sort_values("_safe_count").drop("_safe_count", 1)


def visualize_maelstrom(outdir, sig_cutoff=3, pfmfile=None):

    config = MotifConfig()
    if pfmfile is None:
        pfmfile = config.get_default_params().get("motif_db", None)
        pfmfile = os.path.join(config.get_motif_dir(), pfmfile)

    mapfile = pfmfile.replace(".pwm", ".motif2factors.txt")
    if os.path.exists(mapfile):

        m2f = pd.read_csv(
            mapfile, sep="\t", names=["motif", "factors"], index_col=0, comment="#"
        )
        m2f["factors"] = m2f["factors"].str[:50]
    else:
        motifs = [m.id for m in read_motifs(pfmfile)]
        m2f = pd.DataFrame({"factors": motifs}, index=motifs)

    sig_fname = os.path.join(outdir, "final.out.txt")
    df_sig = pd.read_table(sig_fname, index_col=0, comment="#")
    f = np.any(df_sig >= sig_cutoff, 1)
    vis = df_sig[f]
    if vis.shape[0] == 0:
        logger.info("No motifs reach the threshold, skipping visualization.\n")
        return

    # cluster rows
    row_linkage = hierarchy.linkage(pdist(vis, metric="euclidean"), method="complete")
    idx = hierarchy.leaves_list(row_linkage)

    plt.figure()

    vis = safe_join(vis, m2f).set_index("factors")

    # size of figure
    size = [2 + vis.shape[1] * 0.4, 1.8 + vis.shape[0] * 0.3]

    cg = sns.heatmap(
        vis.iloc[idx],
        cmap="viridis",
        yticklabels=True,
        cbar_kws={"orientation": "horizontal"},
    )
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
        cg = sns.heatmap(
            vis_freq,
            cmap="viridis",
            yticklabels=True,
            vmin=0,
            vmax=0.2,
            cbar_kws={"orientation": "horizontal"},
        )
        # idx = cg.dendrogram_row.reordered_ind
        _ = plt.setp(cg.yaxis.get_majorticklabels(), rotation=0)
        plt.title("Motif Frequency")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "motif.frequency.png"), dpi=300)

        plt.figure(figsize=size)

        bla = vis_freq.min(1)
        bla[bla < 0.01] = 0.01

        cg = sns.heatmap(
            np.log2(vis_freq.apply(lambda x: x / bla, 0)),
            yticklabels=True,
            vmin=-5,
            vmax=5,
            cbar_kws={"orientation": "horizontal"},
        )
        # idx = cg.dendrogram_row.reordered_ind
        _ = plt.setp(cg.yaxis.get_majorticklabels(), rotation=0)
        plt.title("Motif Enrichment")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, "motif.enrichment.png"), dpi=300)


def df_rank_aggregation(df, dfs, exps, method="int_stouffer"):
    df_p = pd.DataFrame(index=list(dfs.values())[0].index)
    names = list(dfs.values())[0].columns
    dfs = [
        pd.concat([v[col].rename(k, inplace=True) for k, v in dfs.items()], axis=1)
        for col in names
    ]
    pool = Pool(16)
    func = partial(rankagg, method=method)
    ret = pool.map(func, dfs)
    pool.close()
    pool.join()

    for name, result in zip(names, ret):
        df_p[name] = result

    if df.shape[1] != 1:
        df_p = df_p[df.columns]

    if method == "int_stouffer":
        df_p.columns = ["z-score " + c for c in df_p.columns]
    else:
        df_p.columns = ["activity " + c for c in df_p.columns]
    return df_p


def run_maelstrom(
    infile,
    genome,
    outdir,
    pfmfile=None,
    filter_redundant=True,
    filter_cutoff=0.8,
    plot=True,
    cluster=False,
    score_table=None,
    count_table=None,
    methods=None,
    ncpus=None,
    zscore=True,
    gc=True,
    center=False,
    aggregation="int_stouffer",
):
    """Run maelstrom on an input table.

    Parameters
    ----------
    infile : str
        Filename of input table. Can be either a text-separated tab file or a
        feather file.

    genome : str
        Genome name. Can be either the name of a FASTA-formatted file or a
        genomepy genome name.

    outdir : str
        Output directory for all results.

    pfmfile : str, optional
        Specify a PFM file for scanning.

    filter_redundant : bool, optional
        Create a non-redundant set of motifs based on correlation of motif scores in the input data.

    filter_cutoff : float, optional
        Cutoff to use for non-redundant motif selection. Default is 0.8.

    plot : bool, optional
        Create heatmaps.

    cluster : bool, optional
        If True and if the input table has more than one column, the data is
        clustered and the cluster activity methods are also run. Not
        well-tested.

    score_table : str, optional
        Filename of pre-calculated table with motif scores.

    count_table : str, optional
        Filename of pre-calculated table with motif counts.

    methods : list, optional
        Activity methods to use. By default are all used.

    ncpus : int, optional
        If defined this specifies the number of cores to use.

    zscore : bool, optional
        Use z-score normalized motif scores.

    gc : bool, optional
        Use GC% bins to normalize motif scores.

    center : bool, optional
        Mean-center the input table.

    aggregation: str, optional
        How to combine scores of the predictors. The default is "int_stouffer", for
        inverse normal transform followed by Stouffer's methods to combine z-scores.
        Alternatively, "stuart" performs rank aggregation and reports the -log10 of
        the rank aggregation p-value.
    """
    logger.info("Starting maelstrom")
    if infile.endswith("feather"):
        df = pd.read_feather(infile)
        df = df.set_index(df.columns[0])
    else:
        df = pd.read_table(infile, index_col=0, comment="#")

    # Check if the input is mean-centered
    if df.shape[1] > 1 and not np.allclose(df.mean(1), 0):
        if center:
            logger.info(
                "Input is not mean-centered, setting the mean of all rows to 0."
            )
            logger.info(
                "Use --nocenter if you know what you're doing and want to change this behavior."
            )
            logger.info(
                "Note that if you use count data (ChIP-seq, ATAC-seq) we recommend to "
                "first transform your data, for instance using log2(), and to normalize "
                "between samples. To create a table suitable for maelstrom you can use the "
                "coverage_table script included with GimmeMotifs."
            )
            df = df.sub(df.mean(axis=1), axis=0)
        else:
            logger.info("Input is not mean-centered, but --nocenter was specified.")
            logger.info(
                "Leaving the data as-is, but make sure this is what your really want."
            )

    # Check for duplicates
    if df.index.duplicated(keep=False).any():
        logger.warning("Input file contains duplicate regions!")
        logger.warning("These will be removed.")
        df = df.iloc[~df.index.duplicated(keep=False)]

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if methods is None:
        methods = Moap.list_predictors()
    methods = [m.lower() for m in methods]

    df.to_csv(os.path.join(outdir, "input.table.txt"), sep="\t")
    infile = os.path.join(outdir, "input.table.txt")

    # Copy the motif informatuon
    pfmfile = pfmfile_location(pfmfile)
    if pfmfile:
        shutil.copy2(pfmfile, outdir)
        mapfile = re.sub(".p[fw]m$", ".motif2factors.txt", pfmfile)
        if os.path.exists(mapfile):
            shutil.copy2(mapfile, outdir)

    # Create a file with the number of motif matches
    if count_table is None:
        count_table = os.path.join(outdir, "motif.count.txt.gz")
        if not os.path.exists(count_table):
            logger.info("motif scanning (counts)")
            counts = scan_regionfile_to_table(
                infile,
                genome,
                "count",
                pfmfile=pfmfile,
                ncpus=ncpus,
                zscore=zscore,
                gc=gc,
            )
            counts.to_csv(count_table, sep="\t", compression="gzip")
        else:
            logger.info("Counts, using: %s", count_table)

    # Create a file with the score of the best motif match
    if score_table is None:
        score_table = os.path.join(outdir, "motif.score.txt.gz")
        if not os.path.exists(score_table):
            logger.info("motif scanning (scores)")
            scores = scan_regionfile_to_table(
                infile,
                genome,
                "score",
                pfmfile=pfmfile,
                ncpus=ncpus,
                zscore=zscore,
                gc=gc,
            )
            scores.to_csv(
                score_table, sep="\t", float_format="%.3f", compression="gzip"
            )
        else:
            logger.info("Scores, using: %s", score_table)

    counts = pd.read_csv(count_table, index_col=0, comment="#", sep="\t")
    scores = pd.read_csv(score_table, index_col=0, comment="#", sep="\t")

    if filter_redundant:
        logger.info("Selecting non-redundant motifs")
        fa = FeatureAgglomeration(
            distance_threshold=filter_cutoff,
            n_clusters=None,
            affinity="correlation",
            linkage="complete",
            compute_full_tree=True,
        )
        fa.fit(scores)
        X_cluster = pd.DataFrame({"motif": scores.columns, "label": fa.labels_})
        X_cluster = X_cluster.join(scores.var().to_frame(name="var"), on="motif")
        selected_motifs = (
            X_cluster.sort_values("var")
            .drop_duplicates(subset=["label"], keep="last")["motif"]
            .values
        )

        nr_motif = (
            X_cluster.sort_values("var")
            .drop_duplicates(subset=["label"], keep="last")[["label", "motif"]]
            .set_index("label")
        )
        X_cluster = X_cluster.join(nr_motif, rsuffix="_nr", on="label")
        motif_map = X_cluster[["motif", "motif_nr"]].set_index("motif")

        scores = scores[selected_motifs]
        counts = counts[selected_motifs]
        score_table = os.path.join(outdir, "motif.nr.score.txt.gz")
        scores.to_csv(score_table, sep="\t", compression="gzip")
        count_table = os.path.join(outdir, "motif.nr.count.txt.gz")
        counts.to_csv(count_table, sep="\t", compression="gzip")

        m2f = pd.read_table(os.path.join(outdir, mapfile), comment="#")
        m2f = m2f.join(motif_map, on="Motif")
        m2f.loc[m2f["Motif"] != m2f["motif_nr"], "Curated"] = "N"
        m2f["Motif"] = m2f["motif_nr"]
        m2f = m2f.drop(columns=["motif_nr"])

        motifs = read_motifs(pfmfile)
        pfmfile = os.path.join(outdir, "nonredundant.motifs.pfm")
        with open(pfmfile, "w") as f:
            for motif in motifs:
                f.write(f"{motif.to_pfm()}\n")
        mapfile = pfmfile.replace(".pfm", ".motif2factors.txt")
        with open(mapfile, "w") as f:
            f.write(
                "# Note: this mapping is specifically created for this non-redundant set of motifs.\n"
            )
            f.write(
                "# It also includes factors for motifs that were similar, but this can be\n"
            )
            f.write("# specific to this analysis.\n")

        with open(mapfile, "a") as f:
            m2f.to_csv(f, index=False, sep="\t")
        logger.info(f"Selected {len(selected_motifs)} motifs")
        logger.info(f"Motifs: {pfmfile}")
        logger.info(f"Factor mappings: {mapfile}")

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
            clusterfile = os.path.join(
                outdir, os.path.basename(infile) + ".cluster.txt"
            )

            df[:] = scale(df, axis=0)
            names = df.columns
            df_changed = pd.DataFrame(index=df.index)
            df_changed["cluster"] = np.nan
            for name in names:
                df_changed.loc[
                    (df[name] - df.loc[:, df.columns != name].max(1)) > 0.5, "cluster"
                ] = name
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
            if scoring == "count":
                moap_with_table(
                    fname, count_table, outdir, method, scoring, ncpus=ncpus
                )
            elif scoring == "score":
                moap_with_table(
                    fname, score_table, outdir, method, scoring, ncpus=ncpus
                )

        except Exception as e:
            logger.warn("Method %s with scoring %s failed", method, scoring)
            logger.warn(e)
            logger.warn("Skipping")
            raise
    dfs = {}
    for method, scoring, fname in exps:
        t = "{}.{}".format(method, scoring)
        fname = os.path.join(outdir, "activity.{}.{}.out.txt".format(method, scoring))
        try:
            dfs[t] = pd.read_table(fname, index_col=0, comment="#")
        except FileNotFoundError:
            logger.warn("Activity file for {} not found!\n".format(t))

    if len(methods) > 1:
        logger.info("Rank aggregation")
        df_p = df_rank_aggregation(df, dfs, exps, method=aggregation)

        # Add percentage of input sequences with motif
        if df.shape[1] > 1:
            df_p["% with motif"] = counts[df_p.index].sum(0) / df.shape[0] * 100
        else:
            bla = counts.join(df).groupby(df.columns[0]).mean() * 100
            bla = bla.T
            bla = bla.rename(
                columns={col: f"{col} % with motif" for col in bla.columns}
            )
            df_p = df_p.join(bla)

        if df.shape[1] > 1:
            # Add correlation between motif score and signal
            logger.info("Correlation")
            for col in df.columns:
                df_p[f"corr {col}"] = 0
                for motif in df_p.index:
                    df_p.loc[motif, f"corr {col}"] = pearsonr(df[col], scores[motif])[0]

        df_p.to_csv(os.path.join(outdir, "final.out.txt"), sep="\t")
    # df_p = df_p.join(m2f)

    # Write motif frequency table

    if df.shape[1] == 1:
        mcount = df.join(counts)
        m_group = mcount.groupby(df.columns[0])
        freq = m_group.sum() / m_group.count()
        freq.to_csv(os.path.join(outdir, "motif.freq.txt"), sep="\t")

    if plot and len(methods) > 1:
        logger.info("html report")
        maelstrom_html_report(outdir, os.path.join(outdir, "final.out.txt"), pfmfile)
        logger.info(os.path.join(outdir, "gimme.maelstrom.report.html"))


def _get_factor_list(motif, indirect=False):
    factor_list = motif.factors[DIRECT_NAME]
    if indirect:
        factor_list += motif.factors[INDIRECT_NAME]

    return list(set([f.upper() for f in factor_list]))


class MaelstromResult:
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
        fnames = glob.glob(os.path.join(outdir, "nonredundant*.p[fw]m"))

        if len(fnames) == 0:
            fnames = glob.glob(os.path.join(outdir, "*.p[fw]m"))
        if len(fnames) > 0:
            pfmfile = fnames[0]
            with open(pfmfile) as fin:
                self.motifs = {m.id: m for m in read_motifs(fin)}

        self.activity = {}
        # Read individual activity files
        for fname in glob.glob(os.path.join(outdir, "activity*txt")):
            # print()
            _, name, mtype, _, _ = os.path.split(fname)[-1].split(".")
            self.activity["{}.{}".format(name, mtype)] = pd.read_table(
                fname, comment="#", index_col=0
            )

        # Read rank aggregation
        self.result = pd.read_table(
            os.path.join(outdir, "final.out.txt"), comment="#", index_col=0
        )
        self.correlation = self.result.loc[:, self.result.columns.str.contains("corr")]
        self.percent_match = self.result.loc[
            :, self.result.columns.str.contains("% with motif")
        ]
        self.result = self.result.loc[
            :,
            ~self.result.columns.str.contains("corr")
            & ~self.result.columns.str.contains("% with motif"),
        ]

        # Read motif results
        self.scores = pd.read_table(
            os.path.join(outdir, "motif.score.txt.gz"), index_col=0
        )
        self.counts = pd.read_table(
            os.path.join(outdir, "motif.count.txt.gz"), index_col=0
        )
        fname = os.path.join(outdir, "motif.freq.txt")
        if os.path.exists(fname):
            self.freq = pd.read_table(fname, index_col=0)

        # Read original input file
        try:
            self.input = pd.read_table(
                os.path.join(outdir, "input.table.txt"), index_col=0
            )
            if self.input.shape[1] == 1:
                self.input.columns = ["cluster"]
        except Exception:
            pass

    def plot_heatmap(
        self,
        kind="final",
        min_freq=0.01,
        threshold=2,
        name=True,
        indirect=True,
        figsize=None,
        max_number_factors=5,
        aspect=1,
        cmap="RdBu_r",
        **kwargs,
    ):
        """Plot clustered heatmap of predicted motif activity.

        Parameters
        ----------
        kind : str, optional
            Which data type to use for plotting. Default is 'final', which will
            plot the result of the rank aggregation. Other options are 'freq'
            for the motif frequencies, or any of the individual activities such
            as 'rf.score'.

        min_freq : float, optional
            Minimum frequency of motif occurrence.

        threshold : float, optional
            Minimum activity (absolute) of the rank aggregation result.

        name : bool, optional
            Use factor names instead of motif names for plotting.

        indirect : bool, optional
            Include indirect factors (computationally predicted or non-curated). Default is True.

        max_number_factors : int, optional
            Truncate the list of factors to this maximum size.

        figsize : tuple, optional
            Tuple of figure size (width, height).

        aspect : int, optional
            Aspect ratio for tweaking the plot.

        cmap : str, optional
            Color paletter to use, RdBu_r by default.

        kwargs : other keyword arguments
            All other keyword arguments are passed to sns.heatmap

        Returns
        -------
        cg : ClusterGrid
            A seaborn ClusterGrid instance.
        """

        filt = np.any(np.abs(self.result) >= threshold, 1)
        if hasattr(self, "freq"):
            filt = filt & np.any(np.abs(self.freq.T) >= min_freq, 1)
        else:
            filt = filt & (self.counts.sum() / self.counts.shape[0] > min_freq)

        idx = self.result.loc[filt].index

        if idx.shape[0] == 0:
            logger.warning("Empty matrix, try lowering the threshold")
            return

        if idx.shape[0] >= 100:
            logger.warning("The filtered matrix has more than 100 rows.")
            logger.warning(
                "It might be worthwhile to increase the threshold for visualization"
            )

        if kind == "final":
            data = self.result
        elif kind == "freq":
            if hasattr(self, "freq"):
                data = self.freq.T
                cmap = "Reds"
            else:
                raise ValueError(
                    "frequency plot only works with maelstrom output from clusters"
                )
        elif kind in self.activity:
            data = self.activity[kind]
            if kind in ["hypergeom.count", "mwu.score"]:
                cmap = "Reds"
        else:
            raise ValueError("Unknown dtype")

        m = data.loc[idx]

        if "vmax" in kwargs:
            vmax = kwargs.pop("vmax")
        else:
            vmax = max(abs(np.percentile(m, 1)), np.percentile(m, 99))

        if "vmin" in kwargs:
            vmin = kwargs.pop("vmin")
        else:
            vmin = -vmax

        if name:
            m["factors"] = [
                self.motifs[n].format_factors(
                    max_length=max_number_factors,
                    html=False,
                    include_indirect=indirect,
                    extra_str=",..",
                )
                for n in m.index
            ]
            m = m.set_index("factors")
        h, w = m.shape

        if figsize is None:
            figsize = (4 + m.shape[1] / 4, 1 + m.shape[0] / 3)

        fig = plt.figure(figsize=figsize)
        npixels = 30
        g = GridSpec(
            2, 1, height_ratios=(fig.get_figheight() * fig.dpi - npixels, npixels)
        )
        ax1 = fig.add_subplot(g[0, :])
        ax2 = fig.add_subplot(g[1, :])
        ax2.set_title("aggregated z-score")
        dm = pdist(m, metric="correlation")
        hc = linkage(dm, method="ward")
        leaves = dendrogram(hc, no_plot=True)["leaves"]
        cg = sns.heatmap(
            m.iloc[leaves],
            ax=ax1,
            cbar_ax=ax2,
            cbar_kws={"orientation": "horizontal"},
            cmap=cmap,
            linewidths=1,
            vmin=vmin,
            vmax=vmax,
            **kwargs,
        )
        plt.setp(cg.axes.xaxis.get_majorticklabels(), rotation=90)
        plt.tight_layout()
        # cg.ax_col_dendrogram.set_visible(False)
        # plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
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

        data = self.scores[plot_motifs]
        data[:] = data.scale(data, axix=0)
        if name:
            data = data.T
            data["factors"] = [
                join_max(self.motifs[n].factors, max_len, ",", suffix=",(...)")
                for n in plot_motifs
            ]
            data = data.set_index("factors").T

        data = pd.melt(self.input.join(data), id_vars=["cluster"])
        data.columns = ["cluster", "motif", "z-score"]
        g = sns.factorplot(
            data=data, y="motif", x="z-score", hue="cluster", kind="box", aspect=2
        )
        return g
