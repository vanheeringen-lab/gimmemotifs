#!/usr/bin/python -W ignore
# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

import sys
import os

from gimmemotifs.motif import read_motifs
from gimmemotifs.rocmetrics import ROC_values, ROC_AUC, MNCP, max_enrichment, enr_at_fdr, recall_at_fdr
from gimmemotifs.fasta import Fasta
from gimmemotifs.plot import roc_plot
from gimmemotifs.scanner import Scanner

def roc(args):
    """ Calculate ROC_AUC and other metrics and optionally plot ROC curve.
    """
    pwmfile = args.pwmfile
    fg_file = args.sample
    bg_file = args.background
    outputfile = args.outfile
    # Default extension for image
    if outputfile and   not outputfile.endswith(".png"):
        outputfile += ".png"
    
    motifs = read_motifs(open(pwmfile), fmt="pwm")

    s = Scanner()
    s.set_motifs(pwmfile)
    
    ids = []
    if args.ids:
        ids = args.ids.split(",")
    else:
        ids = [m.id for m in motifs]

    fg_total = dict([(m.id, []) for m in motifs])
    for scores in s.best_score(fg_file):
        for motif,score in zip(motifs, scores):
            fg_total[motif.id].append(score)
    
    bg_total = dict([(m.id, []) for m in motifs])
    for scores in s.best_score(bg_file):
        for motif,score in zip(motifs, scores):
            bg_total[motif.id].append(score)
   
    plot_x = []
    plot_y = []
    # Print the metrics
    print "Motif\tROC AUC\tMNCP\tEnr. at 5% FDR\tMax enr.\tRecall at 10% FDR"
    for motif_id in ids:
        fg_vals = fg_total[motif_id] 
        bg_vals = bg_total[motif_id]    
        (x, y) = ROC_values(fg_vals, bg_vals) 
        plot_x.append(x)
        plot_y.append(y)
        auc = ROC_AUC(fg_vals, bg_vals)
        mncp = MNCP(fg_vals, bg_vals)
        enr_fdr = enr_at_fdr(fg_vals, bg_vals)
        max_enr,score = max_enrichment(fg_vals, bg_vals)
        recall = recall_at_fdr(fg_vals, bg_vals, 0.1)
        print "%s\t%0.3f\t%03f\t%0.2f\t%0.2f\t%0.4f" % (
                motif_id, auc, mncp, enr_fdr, max_enr, recall)
    
    # Plot the ROC curve
    if outputfile:
        roc_plot(outputfile, plot_x, plot_y, ids=ids)
