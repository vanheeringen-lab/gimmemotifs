#!/usr/bin/python
# Copyright (c) 2009-2019 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
"""Command line function 'scan'."""
import sys

from numpy.random import RandomState

from gimmemotifs.scanner import scan_to_file


def pfmscan(args):
    random_state = None
    if args.seed is not None:
        random_state = RandomState(int(args.seed))

    scan_to_file(
        args.inputfile,
        args.pfmfile,
        sys.stdout,
        nreport=args.nreport,
        fpr=args.fpr,
        cutoff=args.cutoff,
        bed=args.bed,
        scan_rc=args.scan_rc,
        table=args.table,
        score_table=args.score_table,
        bgfile=args.bgfile,
        genome=args.genome,
        ncpus=args.ncpus,
        zscore=args.zscore,
        gcnorm=args.gcnorm,
        random_state=random_state,
    )
