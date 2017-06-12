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
            "roc_auc", 
            "mncp", 
            "enr_at_fpr",
            "max_enrichment", 
            "recall_at_fdr", 
            "roc_values"
            ]
    
    motif_stats = calc_stats(motifs, args.sample, args.background, 
            genome=args.genome, stats=stats)

    plot_x = []
    plot_y = []
    legend = []
    # Print the metrics
    print("Motif\tROC AUC\tMNCP\tEnr. at 5% FDR\tMax enr.\tRecall at 10% FDR")
    for motif in motifs:
        if outputfile:
            x, y = motif_stats[str(motif)]["roc_values"]
            plot_x.append(x)
            plot_y.append(y)
            legend.append(motif.id)
        print("{}\t{:.3f}\t{:.3f}\t{:.2f}\t{:0.2f}\t{:0.4f}".format(
              motif.id, 
              motif_stats[str(motif)]["roc_auc"], 
              motif_stats[str(motif)]["mncp"], 
              motif_stats[str(motif)]["enr_at_fpr"], 
              motif_stats[str(motif)]["max_enrichment"][0], 
              motif_stats[str(motif)]["recall_at_fdr"],
              ))
    
    # Plot the ROC curve
    if outputfile:
        roc_plot(outputfile, plot_x, plot_y, ids=legend)
