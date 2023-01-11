#!/usr/bin/python
# Copyright (c) 2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
"""Command line function 'maelstrom'"""
import os

from numpy.random import RandomState

from gimmemotifs.maelstrom import run_maelstrom


def maelstrom(args):
    """Find differential motifs."""
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
    plot_all_motifs = args.plot_all_motifs
    plot_no_motifs = args.plot_no_motifs

    if not os.path.exists(infile):
        raise ValueError(f"file {infile} does not exist")

    if methods:
        methods = [x.strip() for x in methods.split(",")]

    random_state = None
    if args.seed is not None:
        random_state = RandomState(int(args.seed))

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
        plot_all_motifs=plot_all_motifs,
        plot_no_motifs=plot_no_motifs,
        random_state=random_state,
    )
