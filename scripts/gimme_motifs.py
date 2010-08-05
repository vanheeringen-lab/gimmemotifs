#!/usr/bin/env python 
# Copyright (c) 2009-2010 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
from optparse import OptionParser
from gimmemotifs.core import *
from gimmemotifs.config import *

VERSION = "0.60"
config = MotifConfig()
params = config.get_default_params()

parser = OptionParser(version="%prog " + VERSION)
parser.add_option("-i", "--inputfile", dest="inputfile", help="Inputfile in bed format", metavar="FILE")
parser.add_option("-k", "--keepintermediate", dest="keep_intermediate", help="Don't delete intermediate files", default=False, action="store_true")
parser.add_option("-n", "--name", dest="name", help="Give your analysis a name", metavar="NAME")
parser.add_option("-a", "--analysis",dest="analysis", help="Analysis type: small, medium, large, xl (%s)" % params["analysis"], metavar="ANALYSIS", default=params["analysis"])
parser.add_option("-g", "--genome", dest="genome", help="Genome version (%s)" % (params["genome"]), metavar="VERSION", default=params["genome"])
parser.add_option("-s", "--singlestrand", dest="single", help="Only predict motifs for single + strand (default is both)", action="store_true", default=False)
parser.add_option("-f", "--fraction", dest="fraction", help="Fraction of peaks to use for motif predicton (%s)" % params["fraction"], metavar="FRACTION", default=params["fraction"], type=float)
parser.add_option("-w", "--width", dest="width", help="Width to use for motif prediction (%s)" % params["width"], metavar="N", default=params["width"], type=int)
parser.add_option("-e", "--enrichment", dest="enrichment", help="Motif significance: enrichment cutoff (>%s)" % params["enrichment"], metavar="N", default=params["enrichment"], type=float)
parser.add_option("-p", "--pvalue", dest="pvalue", help="Motif significance: p-value cutoff (<%s)" % params["pvalue"], metavar="N", default=params["pvalue"], type=float)
parser.add_option("-b", "--background", dest="background", help="Background to determine significance genomic_matched,random (%s)" % params["background"], metavar="N", default=params["background"])
parser.add_option("-l", "--localization_width", dest="lwidth", help="Width to use for motif localization graphs (%s)" % params["lwidth"], metavar="N", default=params["lwidth"], type=int)
parser.add_option("-t", "--tools", dest="tools", help="Tools to use, any combination of %s (default %s)" % (params["available_tools"], params["tools"]), metavar="N", default=params["tools"])
parser.add_option("--max_time", dest="max_time", help="Time limit for motif prediction in hours (default: %s)" % str(params["max_time"]), metavar="HOURS", default=params["max_time"])
parser.add_option("-x", dest="weird_option", help="Do NOT use this!", default=False, action="store_true")

(options, args) = parser.parse_args()

if not options.inputfile:
	parser.print_help()
	sys.exit()

if not os.path.exists(options.inputfile):
	print "File %s does not exist!" % options.inputfile
	sys.exit()

background = [x.strip() for x in options.background.split(",")]
for bg in background:
	if not bg in ["genomic", "random", "genomic_matched"]:
		print "Invalid value for background argument"
		sys.exit()

if options.lwidth < options.width:
	sys.stderr.write("Warning: localization width is smaller than motif prediction width!")

available_tools = [t.strip() for t in params["available_tools"].split(",")]
tools = dict([(t,0) for t in available_tools])
for t in options.tools.split(","):
	tools[t.strip()] = 1

for tool in tools.keys():
	if tool not in available_tools:
		print "Sorry, motif prediction tool %s is not supported" % (tool)
		sys.exit()

if "matched_genomic" in background:
	gene_file = os.path.join(config.get_gene_dir(), options.genome)
	if not os.path.exists(gene_file):
		print "Sorry, genomic background for %s is not supported!" % options.genome
		sys.exit()


params = {
	"genome": options.genome,
	"width": options.width,
	"lwidth": options.lwidth,
	"tools": options.tools,
	"analysis": options.analysis,
	"pvalue": options.pvalue,
	"background": options.background,
	"enrichment": options.enrichment,
	"fraction": options.fraction,
	"use_strand": options.single,
	"keep_intermediate": options.keep_intermediate,
	"max_time": options.max_time,
	"weird_option": options.weird_option
}

gm = GimmeMotifs(options.name)
gm.run_full_analysis(options.inputfile, params)
