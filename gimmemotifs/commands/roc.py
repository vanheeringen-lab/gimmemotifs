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
from gimmemotifs.stats import calc_stats
from gimmemotifs.config import MotifConfig
import numpy as np
import pandas as pd
import os
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
            "Enr. at 1% FPR",
            "Recall at 10% FDR"
    ]
    
    m2f = pwmfile.replace(".pwm", ".motif2factors.txt")
    if os.path.exists(m2f):
        sys.stderr.write("reading mapping\n")
        m2f = pd.read_table(m2f, index_col=0)
        m2f.columns = ["factors"]
        f = m2f["factors"].str.len() > 30        
        m2f["factors"] = '<div title="' + m2f["factors"] + '">' + m2f["factors"].str.slice(0,30) 
        m2f.loc[f, "factors"] += '(...)'
        m2f['factors'] += '</div>'
        df = df.join(m2f)
        cols = ["factors"] + cols
    
    df = df[df["corrected P-value"] <= threshold]
    
    df["Logo"] = ['<img src="logos/{}.png" height=40/>'.format(x) for x in list(df.index)]
    
    df = df[cols]
    if not os.path.exists(outdir + "/logos"):
        os.makedirs(outdir + "/logos")
    for motif in read_motifs(open(pwmfile)):
        if motif.id in df.index:
            motif.to_img(outdir + "/logos/{}.png".format(motif.id), fmt="PNG")
    
    bar_cols = [
            "log10 P-value", "ROC AUC", "MNCP",
            "Enr. at 1% FDR", "Max enr.", "Recall at 10% FDR"
            ]
    template_dir = MotifConfig().get_template_dir()
    js = open(os.path.join(template_dir, "sortable/sortable.min.js"), encoding="utf-8").read()
    css = open(os.path.join(template_dir, "sortable/sortable-theme-slick.css"), encoding="utf-8").read()
    with open(outdir + "/gimme.roc.report.html", "w",  encoding="utf-8") as f: 
        f.write("<head>\n")
        f.write("<style>{}</style>\n".format(css))
        f.write("</head>\n")
        f.write("<body>\n")
        f.write(df.sort_values("ROC AUC", ascending=False).style.bar(bar_cols).set_precision(3).set_table_attributes("data-sortable").render().replace("data-sortable", 'class="sortable-theme-slick" data-sortable'))
    
        f.write("<script>{}</script>\n".format(js))
        f.write("</body>\n")

def roc(args):
    """ Calculate ROC_AUC and other metrics and optionally plot ROC curve."""
    outputfile = args.outfile
    # Default extension for image
    if outputfile and not outputfile.endswith(".png"):
        outputfile += ".png"
    
    motifs = read_motifs(open(args.pwmfile), fmt="pwm")

    ids = []
    if args.ids:
        ids = args.ids.split(",")
    else:
        ids = [m.id for m in motifs]
    motifs = [m for m in motifs if (m.id in ids)]
    
    stats = [
            "phyper_at_fpr",
            "roc_auc", 
            "enr_at_fpr",
            "max_enrichment", 
            "recall_at_fdr", 
            "roc_values",
            "matches_at_fpr",
            ]
    
    motif_stats = calc_stats(motifs, args.sample, args.background, 
            genome=args.genome, stats=stats)

    plot_x = []
    plot_y = []
    legend = []
    
    f_out = sys.stdout
    if args.outdir:
        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)
        f_out = open(args.outdir + "/gimme.roc.report.txt", "w")
    
    # Print the metrics
    f_out.write("Motif\t# matches\t# matches background\tP-value\tlog10 P-value\tROC AUC\tEnr. at 1% FPR\tRecall at 10% FDR\n")
    for motif in motifs:
        if outputfile:
            x, y = motif_stats[str(motif)]["roc_values"]
            plot_x.append(x)
            plot_y.append(y)
            legend.append(motif.id)
        log_pvalue = np.inf
        if motif_stats[str(motif)]["phyper_at_fpr"] > 0:
            log_pvalue = -np.log10(motif_stats[str(motif)]["phyper_at_fpr"])
        f_out.write("{}\t{:d}\t{:d}\t{:.2e}\t{:.3f}\t{:.3f}\t{:.2f}\t{:0.4f}\n".format(
              motif.id, 
              motif_stats[str(motif)]["matches_at_fpr"][0], 
              motif_stats[str(motif)]["matches_at_fpr"][1], 
              motif_stats[str(motif)]["phyper_at_fpr"], 
              log_pvalue, 
              motif_stats[str(motif)]["roc_auc"], 
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
