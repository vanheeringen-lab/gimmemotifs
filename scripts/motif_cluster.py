#!/usr/bin/env python
# Copyright (c) 2009-2010 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

from gimmemotifs.motif import *
from gimmemotifs.comparison import *
from gimmemotifs.cluster import *
import sys
import os

import kid
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-i", "--inputfile", dest="inputfile", help="Inputfile (PFM format)", metavar="FILE")
parser.add_option("-o", "--outdir", dest="outdir", help="Inputfile in bed format", metavar="FILE")
parser.add_option("-b", "--bg", dest="include_bg", help="bg (don't use this)", default=True, action="store_false")
parser.add_option("-s", "--singlestrand", dest="single", help="Don't compare reverse complements of motifs", default=False,action="store_true")
parser.add_option("-t", "--threshold", dest="threshold", help="Cluster threshold", default=0.95, type="float")

(options, args) = parser.parse_args()
if not options.inputfile or not options.outdir:
	parser.print_help()
	sys.exit()


revcomp = not options.single

outdir = os.path.abspath(options.outdir)
if not os.path.exists(outdir):
	os.mkdir(outdir)

trim_ic = 0.2
clusters = []
motifs = pwmfile_to_motifs(options.inputfile)
if len(motifs) == 1:
	clusters = [[motifs[0], motifs]]
else:
	tree = cluster_motifs(options.inputfile, "total", "wic", "mean", True, threshold=options.threshold, include_bg=options.include_bg)
	clusters = tree.getResult()

ids = []
mc = MotifComparer()

for cluster,members in clusters:
	cluster.trim(trim_ic)
	cluster.to_img(os.path.join(outdir,"%s.png" % cluster.id), format="PNG")
	ids.append([cluster.id, {"src":"%s.png" % cluster.id},[]])
	if len(members) > 1:
		scores = {}
		for motif in members:
			scores[motif] =  mc.compare_motifs(cluster, motif, "total", "wic", "mean", pval=True)	
		add_pos = sorted(scores.values(),cmp=lambda x,y: cmp(x[1], y[1]))[0][1]
		for motif in members:
			score, pos, strand = scores[motif]
			add = pos - add_pos
				
			if strand in [1,"+"]:
				pass
			else:
				#print "RC %s" % motif.id
				rc = motif.rc()
				rc.id = motif.id
				motif = rc
			#print "%s\t%s" % (motif.id, add)	
			motif.to_img(os.path.join(outdir, "%s.png" % motif.id.replace(" ", "_")), format="PNG", add_left=add)
	ids[-1][2] = [dict([("src", "%s.png" % motif.id.replace(" ", "_")), ("alt", motif.id.replace(" ", "_"))]) for motif in members]

kid.enable_import()
template_file = os.path.join(sys.prefix, "share/gimmemotifs/templates/cluster_template.kid")
template = kid.Template(file=template_file, motifs=ids)
f = open(os.path.join(outdir, "cluster_report.html"), "w")
f.write(template.serialize())
f.close()

f = open(os.path.join(outdir, "clustered_motifs.pwm"), "w")
if len(clusters) == 1 and len(clusters[0][1]) == 1:
	f.write("%s\n" % clusters[0][0].to_pwm())
else:
	for motif in tree.get_clustered_motifs():
		f.write("%s\n" % motif.to_pwm())
f.close()
