#!/usr/bin/python
# Copyright (c) 2009-2010 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

import sys
from gimmemotifs.fasta import *
from gimmemotifs.background import *
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--inputfile", dest="inputfile", help="Fasta formatted inputfile", metavar="FILE")
parser.add_option("-o", "--outputfile", dest="outputfile", help="Fasta formatted outputfile", metavar="FILE")
parser.add_option("-n", "--number", dest="number", help="Number of sequence to generate compared to input (1 means same amount)", metavar="NUMBER", default=1, type="int")

(options, args) = parser.parse_args()

inputfile = options.inputfile
out = options.outputfile
if not inputfile or not out:
	parser.print_help()
	sys.exit(0)

f = Fasta(inputfile)
m = MarkovFasta(f)

m.writefasta(out)
