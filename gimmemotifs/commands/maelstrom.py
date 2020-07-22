#!/usr/bin/python
# Copyright (c) 2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
"""Command line function 'maelstrom'"""
import os
from gimmemotifs.maelstrom import run_maelstrom


def maelstrom(args):
    """Run the maelstrom method."""
    infile = args.inputfile
    genome = args.genome
    outdir = args.outdir
    pfmfile = args.pfmfile
    filter_redundant = args.filter_redundant
    filter_cutoff = args.filter_cutoff
    methods = args.methods
    ncpus = args.ncpus
    zscore = args.zscore
    center = args.center
    gc = args.gc
    aggregation = args.aggregation

    if not os.path.exists(infile):
        raise ValueError("file {} does not exist".format(infile))

    if methods:
        methods = [x.strip() for x in methods.split(",")]

    run_maelstrom(
        infile,
        genome,
        outdir,
        pfmfile,
        filter_redundant=filter_redundant,
        filter_cutoff=filter_cutoff,
        methods=methods,
        ncpus=ncpus,
        zscore=zscore,
        gc=gc,
        center=center,
        aggregation=aggregation,
    )
