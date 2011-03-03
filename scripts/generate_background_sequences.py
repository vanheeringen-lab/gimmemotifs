#!/usr/bin/python
# Copyright (c) 2009-2011 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

import sys
import os
from optparse import OptionParser
from gimmemotifs.fasta import *
from gimmemotifs.background import *
from gimmemotifs.config import *
from gimmemotifs.utils import *

BG_TYPES = ["random", "matched_genomic"]

parser = OptionParser()
parser.add_option("-t", "--type", dest="bg_type", help="Type of background sequences to generate (%s)" % ",".join(BG_TYPES), metavar="TYPE")
parser.add_option("-i", "--inputfile", dest="inputfile", help="FASTA (random) or BED (matched_genomic) formatted inputfile", metavar="FILE")
parser.add_option("-o", "--outputfile", dest="outputfile", help="FASTA formatted outputfile", metavar="FILE")
parser.add_option("-n", "--number", dest="number", help="Number of sequence to generate compared to input (1 means same amount, default is 10)", metavar="NUMBER", default=10, type="int")
parser.add_option("-g", "--genome", dest="genome", help="If type is matched_genomic, this specified the organism for which matched genomic background sequences will be selected", metavar="GENOME", default="hg18")
parser.add_option("-m", "--markov_order", dest="markov_order", help="If type is random, this specifies the order of the Markov model (default 1)", metavar="N", default=1, type="int")

(options, args) = parser.parse_args()

inputfile = options.inputfile
out = options.outputfile
bg_type = options.bg_type

if not inputfile or not out or not options.bg_type:
	parser.print_help()
	sys.exit()

if not bg_type in BG_TYPES:
	print "The argument 'type' should be one of: %s" % (",".join(BG_TYPES))
	sys.exit()

if options.bg_type == "random":
	f = Fasta(inputfile)
	m = MarkovFasta(f, multiply=options.number, k=options.markov_order)
	m.writefasta(out)
elif options.bg_type == "matched_genomic":
	if not is_valid_bedfile(inputfile, 3):
		print "A BED file is required for background type 'matched_genomic'"
		sys.exit()
	
	config = MotifConfig()
	index_dir = os.path.join(config.get_index_dir(), options.genome)
	gene_file = os.path.join(config.get_gene_dir(), "%s.bed" % options.genome)
	if not os.path.exists(index_dir):
		print "Index for %s does not exist. Has the genome been indexed for use with GimmeMotifs?" % options.genome
		sys.exit()
	if not os.path.exists(gene_file):
		print "Can't find gene definition for %s (%s). See GimmeMotifs documentation on how to add gene files." % (options.genome, gene_file)
		sys.exit()
	l = median_bed_len(inputfile)	
	m = MatchedGenomicFasta(inputfile, gene_file, index_dir, length=l, multiply=options.number)
	m.writefasta(out)
