#!/usr/bin/env python
# Copyright (c) 2009-2010 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

from gimmemotifs.comparison import *
from gimmemotifs.motif import *
from optparse import  OptionParser
import sys

parser = OptionParser()
parser.add_option("-p", "--pwmfile", dest="pwmfile", help="File with pwms", metavar="FILE")
parser.add_option("-d", "--dbpwmfile", dest="dbpwmfile", help="Databse of pwms to match against", metavar="FILE")

(options, args) = parser.parse_args()
if not (options.pwmfile and options.dbpwmfile):
	parser.print_help()
	sys.exit()

sample = pwmfile_to_motifs(options.pwmfile)
db = pwmfile_to_motifs(options.dbpwmfile)

mc = MotifComparer()

result = mc.get_closest_match(sample, db, "partial", "wic", "mean")
for motif, match in result.items():
	print "%s\t%s\t%0.2f" % (motif, match[0], match[1][0])

#print "%s\t%s\t%0.2e" % (motif.id, result[0].id, 1 - result[1][0])
