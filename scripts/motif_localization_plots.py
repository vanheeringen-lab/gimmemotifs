#!/usr/bin/env python 
# Copyright (c) 2009-2010 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

from gimmemotifs.fasta import *
from gimmemotifs.motif import *
import pp
import sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p", "--pwmfile", dest="pwmfile", help="File with pwms", metavar="FILE")
parser.add_option("-f", "--fastafile", dest="fastafile", help="Fasta formatted file", metavar="FILE")
parser.add_option("-w", "--width", dest="width", help="Set width to W (default: determined from fastafile)", metavar="W")
parser.add_option("-i", "--ids", dest="ids", help="Comma-seperated list of motif ids to plot in ROC (default is all ids)", metavar="IDS")
parser.add_option("-c", "--cutoff", dest="cutoff", help="Cutoff to use (default 0.95)", type="float", default=0.95)

(options, args) = parser.parse_args()

def motif_localization(fastafile, motif, width, outfile, cutoff=0.9):
	NR_HIST_MATCHES = 100
	from gimmemotifs.utils import plot_histogram, ks_pvalue
	from gimmemotifs.fasta import Fasta
	from numpy import array

	matches = motif.pwm_scan(Fasta(fastafile), cutoff=cutoff, nreport=NR_HIST_MATCHES)
	if len(matches) > 0:
		ar = []
		for a in matches.values():
			ar += a
		matches = array(ar)
		p = ks_pvalue(matches, width - len(motif))
		plot_histogram(matches - width / 2 + len(motif) / 2, outfile, xrange=(-width / 2, width / 2), breaks=21, title="%s (p=%0.2e)" % (motif.id, p), xlabel="Position")
		return motif.id, p
	else:
		return motif.id, 1.0

if not options.fastafile and not options.pwmfile:
	parser.print_help()
	sys.exit()

fastafile = options.fastafile
pwmfile = options.pwmfile

lwidth = options.width
if not lwidth:
	f = Fasta(fastafile)
	lwidth = len(f.items()[0][1])
	f = None

job_server = pp.Server(secret="pumpkinrisotto")
jobs = []
motifs = pwmfile_to_motifs(pwmfile)
ids = [motif.id for motif in motifs]
if options.ids:
	ids = options.ids.split(",")

for motif in motifs:
	if motif.id in ids:
		outfile = os.path.join("%s_histogram" % motif.id)
		jobs.append(job_server.submit(motif_localization, (fastafile,motif,lwidth,outfile, options.cutoff), (),()))

for job in jobs:
	job()
