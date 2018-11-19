#!/usr/bin/python -W ignore
# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
"""Command line function 'roc'."""
from __future__ import print_function
from gimmemotifs.motif import read_motifs
from gimmemotifs.plot import roc_plot
from gimmemotifs.stats import calc_stats_iterator
from gimmemotifs.config import MotifConfig, DIRECT_NAME, INDIRECT_NAME
import numpy as np
import pandas as pd
import os
import re
import sys
from statsmodels.stats.multitest import multipletests

def html_report(outdir, infile, pwmfile, threshold=0.01):
    df = pd.read_table(infile, index_col=0)
    del df.index.name
    df["corrected P-value"] = multipletests(df["P-value"], method="fdr_bh")[1]
    
    cols = [
            "Logo",
            "# matches",
            "# matches background",
            "P-value",
            "log10 P-value",
            "corrected P-value",
            "ROC AUC",
            "PR AUC",
            "Enr. at 1% FPR",
            "Recall at 10% FDR"
    ]
   
    
    motifs = read_motifs(pwmfile)
    idx = [motif.id for motif in motifs]
    direct = [",".join(motif.factors[DIRECT_NAME]) for motif in motifs]
    indirect = [",".join(motif.factors[INDIRECT_NAME]) for motif in motifs]
    m2f = pd.DataFrame({DIRECT_NAME:direct, INDIRECT_NAME:indirect}, index=idx)

    factor_cols = [DIRECT_NAME, INDIRECT_NAME]
    if True:
        for factor_col in factor_cols:
            f = m2f[factor_col].str.len() > 30        
            m2f[factor_col] = '<div title="' + m2f[factor_col] + '">' + m2f[factor_col].str.slice(0,30) 
            m2f.loc[f, factor_col] += '(...)'
            m2f[factor_col] += '</div>'
        df = df.join(m2f)
        cols = factor_cols + cols
    
    df = df[df["corrected P-value"] <= threshold]
    
    df["Logo"] = ['<img src="logos/{}.png" height=40/>'.format(re.sub('[^-_\w]+', '_', x)) for x in list(df.index)]
    
    df = df[cols]
    if not os.path.exists(outdir + "/logos"):
        os.makedirs(outdir + "/logos")
    for motif in motifs:
        if motif.id in df.index:
            motif.to_img(outdir + "/logos/{}.png".format(re.sub('[^-_\w]+', '_', motif.id)), fmt="PNG")
    
    bar_cols = [
            "log10 P-value", "ROC AUC", "PR AUC", "MNCP",
            "Enr. at 1% FPR", "Recall at 10% FDR"
            ]
    template_dir = MotifConfig().get_template_dir()
    js = open(os.path.join(template_dir, "sortable/sortable.min.js"), encoding="utf-8").read()
    css = open(os.path.join(template_dir, "sortable/sortable-theme-slick.css"), encoding="utf-8").read()
    with open(outdir + "/gimme.roc.report.html", "w",  encoding="utf-8") as f: 
        f.write("<head>\n")
        f.write("<style>{}</style>\n".format(css))
        f.write("</head>\n")
        f.write("<body>\n")
        if df.shape[0] > 0: 
            f.write(df.sort_values("ROC AUC", ascending=False).style.bar(bar_cols).set_precision(3).set_table_attributes("data-sortable").render().replace("data-sortable", 'class="sortable-theme-slick" data-sortable'))
        else:
            f.write("No enriched motifs found.")
        f.write("<script>{}</script>\n".format(js))
        f.write("</body>\n")

def roc(args):
    """ Calculate ROC_AUC and other metrics and optionally plot ROC curve."""
    outputfile = args.outfile
    # Default extension for image
    if outputfile and not outputfile.endswith(".png"):
        outputfile += ".png"
    
    motifs = read_motifs(args.pwmfile, fmt="pwm")

    ids = []
    if args.ids:
        ids = args.ids.split(",")
    else:
        ids = [m.id for m in motifs]
    motifs = [m for m in motifs if (m.id in ids)]
    
    stats = [
            "phyper_at_fpr",
            "roc_auc", 
            "pr_auc", 
            "enr_at_fpr",
            "recall_at_fdr", 
            "roc_values",
            "matches_at_fpr",
            ]
    
    plot_x = []
    plot_y = []
    legend = []
    
    f_out = sys.stdout
    if args.outdir:
        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)
        f_out = open(args.outdir + "/gimme.roc.report.txt", "w")
    
    # Print the metrics
    f_out.write("Motif\t# matches\t# matches background\tP-value\tlog10 P-value\tROC AUC\tPR AUC\tEnr. at 1% FPR\tRecall at 10% FDR\n")
    
    
    for motif_stats in calc_stats_iterator(motifs, args.sample, args.background, 
            genome=args.genome, stats=stats, ncpus=args.ncpus):
    
        for motif in motifs:
            if str(motif) in motif_stats:
                if outputfile:
                    x, y = motif_stats[str(motif)]["roc_values"]
                    plot_x.append(x)
                    plot_y.append(y)
                    legend.append(motif.id)
                log_pvalue = np.inf
                if motif_stats[str(motif)]["phyper_at_fpr"] > 0:
                    log_pvalue = -np.log10(motif_stats[str(motif)]["phyper_at_fpr"])
                f_out.write("{}\t{:d}\t{:d}\t{:.2e}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.2f}\t{:0.4f}\n".format(
                      motif.id, 
                      motif_stats[str(motif)]["matches_at_fpr"][0], 
                      motif_stats[str(motif)]["matches_at_fpr"][1], 
                      motif_stats[str(motif)]["phyper_at_fpr"], 
                      log_pvalue, 
                      motif_stats[str(motif)]["roc_auc"], 
                      motif_stats[str(motif)]["pr_auc"], 
                      motif_stats[str(motif)]["enr_at_fpr"], 
                      motif_stats[str(motif)]["recall_at_fdr"],
                      ))
    f_out.close() 
    
    if args.outdir:
        html_report(
            args.outdir,
            args.outdir + "/gimme.roc.report.txt",
            args.pwmfile,
            0.01,
            )

    # Plot the ROC curve
    if outputfile:
        roc_plot(outputfile, plot_x, plot_y, ids=legend)
