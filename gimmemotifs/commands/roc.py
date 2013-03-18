#!/usr/bin/python -W ignore
# Copyright (c) 2009-2010 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

import pp
import sys
import os

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from gimmemotifs.motif import pwmfile_to_motifs
from gimmemotifs.rocmetrics import ROC_values, ROC_AUC, MNCP, max_enrichment, enr_at_fdr
from gimmemotifs.fasta import Fasta
from gimmemotifs.plot import roc_plot

def get_scores(motif, file):
	from subprocess import Popen, PIPE
	from tempfile import NamedTemporaryFile

	pwm = NamedTemporaryFile()
	pwm.write(motif.to_pwm())
	pwm.flush()
	
	cmd = "pwmscan.py -i %s -p %s -c 0.0" % (file, pwm.name)
	out = Popen(cmd, shell=True, stdout=PIPE).stdout

	vals = []
	for line in out.readlines():
		vals.append(float(line.split("\t")[5]))
	return vals

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
    
    # Parallel processing
    n_cpu = 8
    job_server = pp.Server(n_cpu, secret="pumpkinrisotto")

    motifs = dict([(x.id, x) for x in pwmfile_to_motifs(pwmfile)])

    ids = []
    if args.ids:
    	ids = args.ids.split(",")
    else:
    	ids = motifs.keys()
	
    fg_jobs = {}
    bg_jobs = {}

    # Do the prediction
    for id in ids:
    	if motifs.has_key(id):
    		bg_jobs[id] = job_server.submit(get_scores, (motifs[id],bg_file,))
    		fg_jobs[id] = job_server.submit(get_scores, (motifs[id],fg_file,))
    	else:
    		print "Wrong id: %s" % id
    		sys.exit(1)

    plot_x = []
    plot_y = []
    # Print the metrics
    print "Motif\tROC AUC\tMNCP\tEnr. at 5% FDR\tMax enr."
    for id in ids:
    	fg_vals = fg_jobs[id]()	
    	bg_vals = bg_jobs[id]()	
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
