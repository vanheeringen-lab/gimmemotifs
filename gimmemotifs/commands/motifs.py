#!/usr/bin/env python 
# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

from __future__ import print_function
import sys
import os

from gimmemotifs.denovo import gimme_motifs
import gimmemotifs.config as cfg

def motifs(args):
    config = cfg.MotifConfig()
    params = config.get_default_params()

    if not os.path.exists(args.inputfile):
        print("File %s does not exist!" % args.inputfile)
        sys.exit(1)
    
    background = [x.strip() for x in args.background.split(",")]
    for bg in background:
        if not bg in (cfg.FA_VALID_BGS + cfg.BED_VALID_BGS):
            print("Invalid value for background argument")
            sys.exit(1)
        if "user" in bg and not args.user_background:
            print("Please specify a background file to use")
            sys.exit(1)
    
    if args.lwidth < args.width:
        sys.stderr.write("Warning: localization width is smaller than motif prediction width!")
    
    available_tools = [t.strip() for t in params["available_tools"].split(",")]
    tools = dict([(t,0) for t in available_tools])
    for t in args.tools.split(","):
        tools[t.strip()] = 1
    
    for tool in tools.keys():
        if tool not in available_tools:
            print("Sorry, motif prediction tool %s is not supported" % (tool))
            sys.exit(1)
    
    if "promoter" in background:
        gene_file = os.path.join(config.get_gene_dir(), args.genome)
        if not os.path.exists(gene_file):
            print("gene annotation for %s is missing, can't do background 'promoter'" % args.genome)
            sys.exit(1)
    
    params = {
        "genome": args.genome,
        "width": args.width,
        "lwidth": args.lwidth,
        "tools": args.tools,
        "analysis": args.analysis,
        "pvalue": args.pvalue,
        "background": args.background,
        "enrichment": args.enrichment,
        "fraction": args.fraction,
        "use_strand": args.single,
        "keep_intermediate": args.keep_intermediate,
        "max_time": args.max_time,
        "markov_model": args.markov_model,
        "user_background": args.user_background,
        "torque": args.torque,
    }
    
    gimme_motifs(args.inputfile, args.name, params)
