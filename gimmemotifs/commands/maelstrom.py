#!/usr/bin/python
"""Command line function 'maelstrom'"""
import os

from numpy.random import RandomState

from gimmemotifs.maelstrom import run_maelstrom


def maelstrom(args):
    """Find differential motifs."""
    infile = args.inputfile
    if not os.path.exists(infile):
        raise ValueError(f"file {infile} does not exist")

    methods = args.methods
    if methods:
        methods = [x.strip() for x in methods.split(",")]

    random_state = None
    if args.seed is not None:
        random_state = RandomState(int(args.seed))

    run_maelstrom(
        infile,
        args.genome,
        args.outdir,
        args.pfmfile,
        filter_redundant=args.filter_redundant,
        filter_cutoff=args.filter_cutoff,
        methods=methods,
        ncpus=args.ncpus,
        zscore=args.zscore,
        gc=args.gc,
        center=args.center,
        aggregation=args.aggregation,
        plot_all_motifs=args.plot_all_motifs,
        plot_no_motifs=args.plot_no_motifs,
        random_state=random_state,
        progress=not args.noprogress,
    )
