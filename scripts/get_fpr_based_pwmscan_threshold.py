#!/usr/bin/python -W ignore
# Copyright (c) 2009-2011 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

import pp
import sys
import os
from gimmemotifs.motif import pwmfile_to_motifs
from gimmemotifs.c_metrics import pwmscan
from optparse import  OptionParser
from gimmemotifs.rocmetrics import * 
from gimmemotifs.fasta import Fasta
from scipy.stats import scoreatpercentile

parser = OptionParser()
parser.add_option("-p", "--pwmfile", dest="pwmfile", help="File with pwms", metavar="FILE")
parser.add_option("-i", "--inputfile", dest="inputfile", help="FASTA file with background sequences", metavar="FILE") 
parser.add_option("-f", "--fpr", dest="fpr", help="Desired fpr", type="float", metavar="FLOAT") 

(options, args) = parser.parse_args()

if not options.pwmfile or not options.inputfile or not options.fpr:
	parser.print_help()
	exit()

f = Fasta(options.inputfile)
motifs = pwmfile_to_motifs(options.pwmfile)

print "Motif\tScore\tCutoff"
for motif in motifs:
	pwm = motif.pwm
	scores = []
	min_score = motif.pwm_min_score()
	for name,seq in f.items():
		result = pwmscan(seq.upper(), pwm, min_score, 1, True)
		score = result[0][0]
		scores.append(score)
	opt_score = scoreatpercentile(scores, 100 - (100 * options.fpr))
	cutoff = (opt_score - min_score) / (motif.pwm_max_score() - min_score)
	print "%s\t%s\t%s" % (motif.id, opt_score , cutoff)
