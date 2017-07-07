#!/usr/bin/python
# Copyright (c) 2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
"""Command line function 'maelstrom'"""
import os
from gimmemotifs.genome_index import check_genome
from gimmemotifs.maelstrom import run_maelstrom

def maelstrom(args):
    """Run the maelstrom method."""
    infile = args.inputfile
    genome = args.genome
    outdir = args.outdir
    pwmfile = args.pwmfile
    methods = args.methods

    if not os.path.exists(infile):
        raise ValueError("file {} does not exist".format(infile))

    if methods:
        methods = [x.strip() for x in methods.split(",")]

    # check if the genome exists and is indexed
    check_genome(genome)
                
    run_maelstrom(infile, genome, outdir, pwmfile, methods=methods)
