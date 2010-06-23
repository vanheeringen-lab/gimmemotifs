#!/usr/bin/python
# Copyright (c) 2009-2010 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

import sys
from gimmemotifs.fasta import Fasta
from gimmemotifs.c_metrics import pwmscan
from gimmemotifs.utils import *
from gimmemotifs.motif import *
from optparse import OptionParser,TitledHelpFormatter
from math import log
import re

VERSION = "1.02"

nreport = 1 

usage = "usage: %prog -i <FILE> [optional arguments]"
parser = OptionParser(version=VERSION, usage=usage, formatter=TitledHelpFormatter(max_help_position=40, short_first=1))
parser.add_option("-i", "--input", dest="inputfile", help="FASTA-formatted inputfile", metavar="FILE")
parser.add_option("-p", "--pwm", dest="pwmfile", help="specify your own PWM file instead of MDmodule file", metavar="FILE")
parser.add_option("-n", "--nreport", dest="nreport", help="report the N best matches", metavar="N")
parser.add_option("-c", "--cutoff", dest="cutoff", help="motif score cutoff (fraction of maxscore, default 0.9)", metavar="", default=0.9)
parser.add_option("-b", "--bed", action="store_true", dest="bed", default=False, help="output bed format")

(options, args) = parser.parse_args()

if not (options.inputfile and (options.pwmfile or options.mdmodulefile)):
	parser.print_help()
	sys.exit(0)

inputfile = options.inputfile

if options.nreport:
	nreport = int(options.nreport)

cutoff = float(options.cutoff)

motifs = pwmfile_to_motifs(options.pwmfile)

bed = options.bed

f = Fasta(inputfile)
strandmap = {-1:"-",1:"+"}
for (id,seq) in f.items():
	for motif in motifs:
		pwm = motif.pwm
		c =  motif.pwm_min_score() + (motif.pwm_max_score() - motif.pwm_min_score()) * cutoff 
		result = pwmscan(seq.upper(), pwm, c, nreport)
		for (score, pos, strand) in result:
			if bed:
				first = id.split(" ")[0]	
				(chr,loc) = first.split(":")
				if loc:
					(start, end) = map(int, loc.split("-"))
					print "%s\t%s\t%s\t%s" % (chr, start + pos, start + pos + len(pwm) , score)
				else:
					print "%s\t%s\t%s\t%s" % (id, pos, pos +  len(pwm), score)
			else:
				print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tmotif_name \"%s\" ; motif_instance \"%s\"" % (
				id, 
				"pwmscan", 
				"misc_feature", 
				pos, pos + len(pwm) , 
				score, 
				strandmap[strand], 
				".", 
				motif.id, 
				seq[pos: pos + len(pwm)]
			)
