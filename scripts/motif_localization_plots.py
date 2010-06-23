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

if not options.fastafile and not options.pwmfile:
	parser.print_help()
	sys.exit()

def motif_localization(fastafile, motif, width, outfile, cutoff, bins=20):
	from tempfile import NamedTemporaryFile
	from subprocess import Popen
	from gimmemotifs.utils import make_gff_histogram, ks_pvalue

	temp = NamedTemporaryFile()
	temp.write(motif.to_pwm())
	temp.flush()

	tempgff = NamedTemporaryFile()
	
	cmd = "pwmscan.py -i %s -p %s -c %s > %s" % (fastafile, temp.name, cutoff, tempgff.name) 
	p = Popen(cmd, shell=True)
	p.communicate()

	make_gff_histogram(tempgff.name,outfile, width ,motif.id, bins)
	return motif.id, ks_pvalue(tempgff.name, width - len(motif))

fastafile = options.fastafile
pwmfile = options.pwmfile

lwidth = options.width
if not lwidth:
	f = Fasta(fastafile)
	lwidth = len(f.items()[0][1])
	f = None

job_server = pp.Server()
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
