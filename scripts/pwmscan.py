#!/usr/bin/python
# Copyright (c) 2009-2012 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

import sys
import re
from optparse import OptionParser,TitledHelpFormatter
import pp

from gimmemotifs.motif import *
from gimmemotifs.fasta import Fasta

VERSION = "1.2"
NREPORT = 1 
MAX_CPUS = 16
CHUNK = 1000

usage = "usage: %prog -i <FILE> [optional arguments]"
parser = OptionParser(version=VERSION, usage=usage, formatter=TitledHelpFormatter(max_help_position=40, short_first=1))
parser.add_option("-i", "--input", dest="inputfile", help="FASTA-formatted inputfile", metavar="FILE")
parser.add_option("-p", "--pwm", dest="pwmfile", help="PWM file with motifs", metavar="FILE")
parser.add_option("-n", "--nreport", dest="nreport", help="report the N best matches", metavar="N", default=NREPORT, type=int)
parser.add_option("-r", "--norc", dest="scan_rc", help="don't scan reverse complement (- strand)", default=True, action="store_false")
parser.add_option("-c", "--cutoff", dest="cutoff", help="motif score cutoff (fraction of maxscore, default 0.9)", metavar="", default=0.9, type=float)
parser.add_option("-b", "--bed", action="store_true", dest="bed", default=False, help="output bed format")

(options, args) = parser.parse_args()

if not (options.inputfile and options.pwmfile):
	parser.print_help()
	sys.exit(0)

inputfile = options.inputfile
nreport = options.nreport
cutoff = options.cutoff
bed = options.bed

def scan_fa_with_motif(fo, motif, cutoff, nreport, bed):
	return motif, motif.pwm_scan_all(fo, cutoff, nreport)

job_server = pp.Server(secret="beetrootsoup")
if job_server.get_ncpus() > MAX_CPUS:
	job_server.set_ncpus(MAX_CPUS)

jobs = []
motifs = pwmfile_to_motifs(options.pwmfile)
fa = Fasta(inputfile)
for motif in motifs:
	for i in range(0, len(fa), CHUNK):
		jobs.append(job_server.submit(scan_fa_with_motif, (fa[i:i + CHUNK], motif, cutoff, nreport, bed), (),()))

strandmap = {-1:"-",1:"+"}
	
for job in jobs:
	motif, result = job()
	for seq_id, matches in result.items():
		for (pos, score, strand) in matches:
			if bed:
				first = seq_id.split(" ")[0]	
				(chrom,loc) = first.split(":")
				if loc:
					(start, end) = map(int, loc.split("-"))
					print "%s\t%s\t%s\t%s" % (chrom, start + pos, start + pos + len(motif) , score)
				else:
					print "%s\t%s\t%s\t%s" % (seq_id, pos, pos +  len(motif), score)
			else:
				print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tmotif_name \"%s\" ; motif_instance \"%s\"" % (
					seq_id, 
					"pwmscan", 
					"misc_feature", 
					pos, pos + len(motif) , 
					score, 
					strandmap[strand], 
					".", 
					motif.id, 
					fa[seq_id][pos: pos + len(motif)]
				)

