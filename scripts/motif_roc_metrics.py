#!/usr/bin/python -W ignore
# Copyright (c) 2009-2010 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

import pp
import sys
import os
from gimmemotifs.motif import pwmfile_to_motifs
from optparse import  OptionParser
from gimmemotifs.rocmetrics import * 
from gimmemotifs.fasta import Fasta

parser = OptionParser()
parser.add_option("-p", "--pwmfile", dest="pwmfile", help="File with pwms", metavar="FILE")
parser.add_option("-s", "--sample", dest="sample", help="Fasta formatted sample file", metavar="FILE") 
parser.add_option("-b", "--background", dest="background", help="Fasta formatted background file", metavar="FILE") 
parser.add_option("-i", "--ids", dest="ids", help="Comma-seperated list of motif ids to plot in ROC (default is all ids)", metavar="IDS") 

(options, args) = parser.parse_args()

if not options.pwmfile or not options.sample or not options.background:
	parser.print_help()
	exit()

if not os.path.exists(options.sample):
	print "File %s does not exist!" % options.sample
	exit()

if not os.path.exists(options.background):
	print "File %s does not exist!" % options.background
	exit()

def get_scores(motif, file):
	from gimmemotifs.fasta import Fasta
	result = motif.pwm_scan_score(Fasta(file), cutoff=0.0, nreport=1)
	vals = [sorted(x)[-1] for x in result.values()]
	return vals

job_server = pp.Server(secret="pumpkinrisotto")

pwmfile = options.pwmfile
fg_file = options.sample
bg_file = options.background
	
motifs = dict([(x.id, x) for x in pwmfile_to_motifs(pwmfile)])

ids = []
if options.ids:
	ids = options.ids.split(",")
else:
	ids = motifs.keys()
	
fg_jobs = {}
bg_jobs = {}

for id in ids:
	if motifs.has_key(id):
		bg_jobs[id] = job_server.submit(get_scores, (motifs[id],bg_file,))
		fg_jobs[id] = job_server.submit(get_scores, (motifs[id],fg_file,))
	else:
		print "Wrong id: %s" % id
		sys.exit()


print "Motif\tROC AUC\tMNCP\tMax f-measure\tSens @ max f-measure"

for id in ids:
	fg_vals = fg_jobs[id]()	
	bg_vals = bg_jobs[id]()	
	(x, y) = ROC_values(fg_vals, bg_vals) 
	auc = ROC_AUC(fg_vals, bg_vals)
	mncp = MNCP(fg_vals, bg_vals)
	max_f, y = max_fmeasure(x,y)
	print "%s\t%s\t%s\t%s\t%s" % (id,auc,mncp,max_f,y)
