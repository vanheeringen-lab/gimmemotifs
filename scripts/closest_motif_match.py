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
parser.add_option("-d", "--dbpwmfile", dest="dbpwmfile", help="Database of pwms to match against", metavar="FILE")

(options, args) = parser.parse_args()
if not (options.pwmfile and options.dbpwmfile):
	parser.print_help()
	sys.exit()

sample = dict([(m.id, m) for m in pwmfile_to_motifs(options.pwmfile)])
db = dict([(m.id, m) for m in pwmfile_to_motifs(options.dbpwmfile)])

mc = MotifComparer()

print "Motif\tMatch\tScore\tP-value"
result = mc.get_closest_match(sample.values(), db.values(), "partial", "wic", "mean")
for motif, match in result.items():
	pval, pos, orient = mc.compare_motifs(sample[motif], db[match[0]], "partial", "wic", "mean", pval=True)
	print "%s\t%s\t%0.2f\t%0.3e" % (motif, match[0], match[1][0], pval)

#print "%s\t%s\t%0.2e" % (motif.id, result[0].id, 1 - result[1][0])
