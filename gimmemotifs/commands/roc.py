#!/usr/bin/python -W ignore
# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

import sys
import os

from gimmemotifs.motif import pwmfile_to_motifs
from gimmemotifs.rocmetrics import ROC_values, ROC_AUC, MNCP, max_enrichment, enr_at_fdr
from gimmemotifs.fasta import Fasta
from gimmemotifs.plot import roc_plot
from gimmemotifs.scan import scan

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
    
    motifs = dict([(x.id, x) for x in pwmfile_to_motifs(pwmfile)])

    ids = []
    if args.ids:
        ids = args.ids.split(",")
    else:
        ids = motifs.keys()

    fg_total = {}
    result = scan(fg_file, [motifs[x] for x in ids], 0.0, 1)    
    for key,m in result.items():
        fg_total[key.id.split("\t")[0]] = [matches[0][1] for matches in m.values()]
   
    bg_total = {}
    result = scan(bg_file, [motifs[x] for x in ids], 0.0, 1)    
    for key,m in result.items():
        bg_total[key.id.split("\t")[0]] = [matches[0][1] for matches in m.values()]
    
    plot_x = []
    plot_y = []
    # Print the metrics
    print "Motif\tROC AUC\tMNCP\tEnr. at 5% FDR\tMax enr."
    for id in ids:
        fg_vals = fg_total[id] 
        bg_vals = bg_total[id]    
        (x, y) = ROC_values(fg_vals, bg_vals) 
        plot_x.append(x)
        plot_y.append(y)
        auc = ROC_AUC(fg_vals, bg_vals)
        mncp = MNCP(fg_vals, bg_vals)
        enr_fdr = enr_at_fdr(fg_vals, bg_vals)
        max_enr,score = max_enrichment(fg_vals, bg_vals)
        print "%s\t%0.3f\t%03f\t%0.2f\t%0.2f" % (id, auc, mncp, enr_fdr, max_enr)
    
    # Plot the ROC curve
    if outputfile:
        roc_plot(outputfile, plot_x, plot_y, ids=ids)
