#!/usr/bin/env python 
# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

from gimmemotifs.genome_index import *
from gimmemotifs.config import *
import os
import sys
from optparse import OptionParser

VERSION = "0.2"
DEFAULT_GENOME = "hg18"

config = MotifConfig()
index_dir = config.get_index_dir()

parser = OptionParser(version="%prog " + VERSION)
parser.add_option("-i", "--inputfile", dest="inputfile", help="Inputfile in bed format", metavar="FILE")
parser.add_option("-f", "--fastafile", dest="fastafile", help="Fastafile output file", metavar="FILE")
parser.add_option("-s", "--start", dest="add_start", help="Extend start" , type="int", default=0)
parser.add_option("-e", "--end", dest="add_end", help="Extend end" , type="int", default=0)
parser.add_option("-g", "--genome", dest="genome", help="Genome (%s)" % DEFAULT_GENOME , default=DEFAULT_GENOME)
parser.add_option("-d", "--directional", dest="use_strand", help="Reverse complement the minus strand features", default=False, action="store_true")

(options, args) = parser.parse_args()

if not options.inputfile or not options.fastafile:
	parser.print_help()
	sys.exit()

if not os.path.exists(options.inputfile):
	print "File %s does not exist!" % options.inputfile
	sys.exit()

if not options.genome in available_genomes(index_dir):
	print "Sorry, %s is not supported!" % options.genome
	sys.exit()

track2fasta(os.path.join(index_dir, options.genome), options.inputfile, options.fastafile, options.add_start, options.add_end, use_strand=options.use_strand)
