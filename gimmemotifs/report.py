# Copyright (c) 2009-2019 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
"""Reports (graphical and text) for motifs statistics."""
import os
import sys
from datetime import datetime
from multiprocessing import Pool
import re
import shutil
import logging

import jinja2
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

from gimmemotifs.comparison import MotifComparer
from gimmemotifs.fasta import Fasta
from gimmemotifs.motif import read_motifs
from gimmemotifs.config import MotifConfig, DIRECT_NAME, INDIRECT_NAME
from gimmemotifs.plot import roc_plot
from gimmemotifs.stats import calc_stats, add_star, write_stats
from gimmemotifs import __version__
from gimmemotifs.utils import motif_localization

logger = logging.getLogger("gimme.report")


def get_roc_values(motif, fg_file, bg_file, genome):
    """Calculate ROC AUC values for ROC plots."""
    try:
        stats = calc_stats(
            fg_file=fg_file,
            bg_file=bg_file,
            motifs=motif,
            genome=genome,
            stats=["roc_values"],
            ncpus=1,
        )
        (x, y) = list(stats.values())[0]["roc_values"]
        return None, x, y
    except Exception as e:
        print(motif)
        print(motif.id)
        print(str(e))
        raise


def create_roc_plots(pfmfile, fgfa, background, outdir, genome):
    """Make ROC plots for all motifs."""
    motifs = read_motifs(pfmfile, fmt="pwm", as_dict=True)
    ncpus = int(MotifConfig().get_default_params()["ncpus"])
    pool = Pool(processes=ncpus)
    jobs = {}
    for bg, fname in background.items():
        for m_id, m in motifs.items():

            k = "{}_{}".format(str(m), bg)
            jobs[k] = pool.apply_async(
                get_roc_values, (motifs[m_id], fgfa, fname, genome)
            )
    imgdir = os.path.join(outdir, "images")
    if not os.path.exists(imgdir):
        os.mkdir(imgdir)

    roc_img_file = os.path.join(outdir, "images", "{}_roc.{}.png")

    for motif in motifs.values():
        for bg in background:
            k = "{}_{}".format(str(motif), bg)
            error, x, y = jobs[k].get()
            if error:
                logger.error("Error in thread: %s", error)
                logger.error("Motif: %s", motif)
                sys.exit(1)
            roc_plot(roc_img_file.format(motif.id, bg), x, y)


def _create_text_report(inputfile, motifs, closest_match, stats, outdir):
    """Create text report of motifs with statistics and database match."""
    my_stats = {}
    for motif in motifs:
        match = closest_match[motif.id]
        my_stats[str(motif)] = {}
        for bg in list(stats.values())[0].keys():
            if str(motif) not in stats:
                logger.error("####")
                logger.error("{} not found".format(str(motif)))
                for s in sorted(stats.keys()):
                    logger.error(s)
                logger.error("####")
            else:
                my_stats[str(motif)][bg] = stats[str(motif)][bg].copy()
                my_stats[str(motif)][bg]["best_match"] = "_".join(
                    match[0].split("_")[:-1]
                )
                my_stats[str(motif)][bg]["best_match_pvalue"] = match[1][-1]

    header = ("# GimmeMotifs version {}\n" "# Inputfile: {}\n").format(
        __version__, inputfile
    )

    write_stats(my_stats, os.path.join(outdir, "stats.{}.txt"), header=header)


def _create_graphical_report(
    inputfile, pwm, background, closest_match, outdir, stats, best_id=None
):
    """Create main gimme_motifs output html report."""
    if best_id is None:
        best_id = {}

    logger.debug("Creating graphical report")

    class ReportMotif(object):
        """Placeholder for motif stats."""

        pass

    config = MotifConfig()

    imgdir = os.path.join(outdir, "images")
    if not os.path.exists(imgdir):
        os.mkdir(imgdir)

    motifs = read_motifs(pwm, fmt="pfm")

    roc_img_file = "%s_roc.%s"

    dbpwm = config.get_default_params()["motif_db"]
    pwmdir = config.get_motif_dir()

    dbmotifs = read_motifs(os.path.join(pwmdir, dbpwm), as_dict=True)

    report_motifs = []
    for motif in motifs:

        rm = ReportMotif()
        rm.id = motif.id
        rm.id_href = {"href": "#%s" % motif.id}
        rm.id_name = {"name": motif.id}
        rm.img = {"src": os.path.join("images", "%s.png" % motif.id)}
        motif.plot_logo(fname=os.path.join(outdir, "images/{}.png".format(motif.id)))

        # TODO: fix best ID
        rm.best = "Gimme"  # best_id[motif.id]

        rm.consensus = motif.to_consensus()
        rm.stars = int(
            np.mean([stats[str(motif)][bg].get("stars", 0) for bg in background]) + 0.5
        )

        rm.bg = {}
        for bg in background:
            rm.bg[bg] = {}
            this_stats = stats.get(str(motif), {}).get(bg)
            # TODO: fix these stats
            rm.bg[bg]["e"] = "%0.2f" % this_stats.get("enr_at_fpr", 1.0)
            rm.bg[bg]["p"] = "%0.2f" % this_stats.get("phyper_at_fpr", 1.0)
            rm.bg[bg]["auc"] = "%0.3f" % this_stats.get("roc_auc", 0.5)
            rm.bg[bg]["mncp"] = "%0.3f" % this_stats.get("mncp", 1.0)
            rm.bg[bg]["roc_img"] = {
                "src": "images/"
                + os.path.basename(roc_img_file % (motif.id, bg))
                + ".png"
            }
            rm.bg[bg][u"roc_img_link"] = {
                u"href": "images/"
                + os.path.basename(roc_img_file % (motif.id, bg))
                + ".png"
            }

        rm.histogram_img = {"data": "images/%s_histogram.svg" % motif.id}
        rm.histogram_link = {"href": "images/%s_histogram.svg" % motif.id}

        match_id = closest_match[motif.id][0]
        dbmotifs[match_id].plot_logo(
            fname=os.path.join(outdir, "images/{}.png".format(match_id))
        )

        rm.match_img = {"src": "images/{}.png".format(match_id)}
        rm.match_id = closest_match[motif.id][0]
        rm.match_pval = "%0.2e" % closest_match[motif.id][1][-1]

        report_motifs.append(rm)

    total_report = os.path.join(outdir, "gimme.denovo.html")

    star_img = os.path.join(config.get_template_dir(), "star.png")
    shutil.copyfile(star_img, os.path.join(outdir, "images", "star.png"))

    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader([config.get_template_dir()])
    )
    template = env.get_template("report_template.jinja.html")
    # TODO: title
    result = template.render(
        motifs=report_motifs,
        inputfile=inputfile,
        date=datetime.today().strftime("%d/%m/%Y"),
        version=__version__,
        bg_types=list(background.keys()),
    )

    with open(total_report, "wb") as f:
        f.write(result.encode("utf-8"))


def create_denovo_motif_report(
    inputfile, pfmfile, fgfa, background, locfa, outdir, params, stats=None
):
    """Create text and graphical (.html) motif reports."""
    logger.info("creating de novo reports")

    motifs = read_motifs(pfmfile, fmt="pwm")

    # ROC plots
    create_roc_plots(pfmfile, fgfa, background, outdir, params["genome"])

    # Closest match in database
    mc = MotifComparer()
    closest_match = mc.get_closest_match(motifs)

    if stats is None:
        stats = {}
        for bg, bgfa in background.items():
            for m, s in calc_stats(fg_file=fgfa, bg_file=bgfa, motifs=motifs).items():
                if m not in stats:
                    stats[m] = {}
                stats[m][bg] = s

    stats = add_star(stats)

    if not params:
        params = {}
    cutoff_fpr = params.get("cutoff_fpr", 0.9)
    lsize = np.median([len(seq) for seq in Fasta(locfa).seqs])

    # Location plots
    logger.debug("Creating localization plots")
    for motif in motifs:
        logger.debug("  {} {}".format(motif.id, motif))
        outfile = os.path.join(outdir, "images/{}_histogram.svg".format(motif.id))
        motif_localization(locfa, motif, lsize, outfile, cutoff=cutoff_fpr)

    # Create reports
    _create_text_report(inputfile, motifs, closest_match, stats, outdir)
    _create_graphical_report(
        inputfile, pfmfile, background, closest_match, outdir, stats
    )


def maelstrom_html_report(outdir, infile, pfmfile=None, threshold=2):
    df = pd.read_table(infile, index_col=0)
    df = df[np.any(abs(df) >= threshold, 1)]

    motifs = read_motifs(pfmfile)

    del df.index.name
    cols = df.columns

    motifs = read_motifs(pfmfile)
    idx = [motif.id for motif in motifs]
    direct = [
        ",".join(sorted(set([x.upper() for x in motif.factors[DIRECT_NAME]])))
        for motif in motifs
    ]
    indirect = [
        ",".join(sorted(set([x.upper() for x in motif.factors[INDIRECT_NAME]])))
        for motif in motifs
    ]
    m2f = pd.DataFrame({DIRECT_NAME: direct, INDIRECT_NAME: indirect}, index=idx)

    factor_cols = [DIRECT_NAME, INDIRECT_NAME]
    if True:
        for factor_col in factor_cols:
            f = m2f[factor_col].str.len() > 30
            m2f[factor_col] = (
                '<div title="'
                + m2f[factor_col]
                + '">'
                + m2f[factor_col].str.slice(0, 30)
            )
            m2f.loc[f, factor_col] += "(...)"
            m2f[factor_col] += "</div>"
        df = df.join(m2f)

    df["logo"] = [
        '<img src="logos/{}.png" height=40/>'.format(re.sub("[()/]", "_", x))
        for x in list(df.index)
    ]

    if not os.path.exists(outdir + "/logos"):
        os.makedirs(outdir + "/logos")
    for motif in motifs:
        if motif.id in df.index:
            motif.plot_logo(
                fname=outdir + "/logos/{}.png".format(re.sub("[()/]", "_", motif.id))
            )

    template_dir = MotifConfig().get_template_dir()
    js = open(
        os.path.join(template_dir, "sortable/sortable.min.js"), encoding="utf-8"
    ).read()
    css = open(
        os.path.join(template_dir, "sortable/sortable-theme-slick.css"),
        encoding="utf-8",
    ).read()
    df = df[factor_cols + ["logo"] + list(cols)]

    df_styled = df.style
    absmax = np.max((abs(df[cols].max().max()), abs(df[cols].min().min())))
    target = absmax * 1.75

    for col in cols:
        smin = df[col].min()
        smax = df[col].max()
        diff = smax - smin
        low = abs((-target - smin) / diff)
        high = (target - smax) / diff
        df_styled = df_styled.background_gradient(
            cmap="RdBu_r", low=low, high=high, subset=[col]
        )

    df_styled = df_styled.set_precision(3)
    df_styled = df_styled.set_table_attributes("data-sortable")
    df_styled = df_styled.render()
    df_styled = df_styled.replace(
        "data-sortable", 'class="sortable-theme-slick" data-sortable'
    )

    with open(outdir + "/gimme.maelstrom.report.html", "w", encoding="utf-8") as f:
        f.write("<head>\n")
        f.write("<style>{}</style>\n".format(css))
        f.write("</head>\n")
        f.write("<body>\n")
        f.write(df_styled)
        f.write("<script>{}</script>\n".format(js))
        f.write("</body>\n")


def roc_html_report(
    outdir,
    infile,
    pfmfile,
    outname="gimme.motifs.html",
    threshold=0.01,
    use_motifs=None,
    link_matches=False,
):
    df = pd.read_table(infile, index_col=0)
    del df.index.name
    df["corrected P-value"] = multipletests(df["P-value"], method="fdr_bh")[1]

    cols = [
        "logo",
        "# matches",
        "# matches background",
        "P-value",
        "log10 P-value",
        "corrected P-value",
        "ROC AUC",
        "PR AUC",
        "Enr. at 1% FPR",
        "Recall at 10% FDR",
    ]

    motifs = read_motifs(pfmfile)
    if use_motifs is not None:
        motifs = [m for m in motifs if m.id in use_motifs]

    idx = [motif.id for motif in motifs]
    df = df.loc[idx]
    direct = [",".join(motif.factors[DIRECT_NAME]) for motif in motifs]
    indirect = [",".join(motif.factors[INDIRECT_NAME]) for motif in motifs]
    m2f = pd.DataFrame({DIRECT_NAME: direct, INDIRECT_NAME: indirect}, index=idx)

    factor_cols = [DIRECT_NAME, INDIRECT_NAME]
    if True:
        for factor_col in factor_cols:
            f = m2f[factor_col].str.len() > 30
            m2f[factor_col] = (
                '<div title="'
                + m2f[factor_col]
                + '">'
                + m2f[factor_col].str.slice(0, 30)
            )
            m2f.loc[f, factor_col] += "(...)"
            m2f[factor_col] += "</div>"
        df = df.join(m2f)
        cols = factor_cols + cols

    df = df[df["corrected P-value"] <= threshold]

    if link_matches:
        df["# matches"] = (
            "<a href=motif_scan_results/"
            + df.index.to_series()
            + ".matches.bed>"
            + df["# matches"].astype(str)
            + "</a>"
        )

    df["logo"] = [
        '<img src="logos/{}.png" height=40/>'.format(re.sub(r"[^-_\w]+", "_", x))
        for x in list(df.index)
    ]

    df = df[cols]
    if not os.path.exists(outdir + "/logos"):
        os.makedirs(outdir + "/logos")
    for motif in motifs:
        if motif.id in df.index:
            motif.plot_logo(
                fname=outdir
                + "/logos/{}.png".format(re.sub(r"[^-_\w]+", "_", motif.id))
            )

    bar_cols = [
        "log10 P-value",
        "ROC AUC",
        "PR AUC",
        "MNCP",
        "Enr. at 1% FPR",
        "Recall at 10% FDR",
    ]
    template_dir = MotifConfig().get_template_dir()
    js = open(
        os.path.join(template_dir, "sortable/sortable.min.js"), encoding="utf-8"
    ).read()
    css = open(
        os.path.join(template_dir, "sortable/sortable-theme-slick.css"),
        encoding="utf-8",
    ).read()
    with open(os.path.join(outdir, outname), "w", encoding="utf-8") as f:
        f.write("<head>\n")
        f.write("<style>{}</style>\n".format(css))
        f.write("</head>\n")
        f.write("<body>\n")
        if df.shape[0] > 0:
            f.write(
                df.sort_values("ROC AUC", ascending=False)
                .style.bar(bar_cols)
                .set_precision(3)
                .set_table_attributes("data-sortable")
                .render()
                .replace("data-sortable", 'class="sortable-theme-slick" data-sortable')
            )
        else:
            f.write("No enriched motifs found.")
        f.write("<script>{}</script>\n".format(js))
        f.write("</body>\n")
