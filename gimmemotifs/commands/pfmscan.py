#!/usr/bin/python
# Copyright (c) 2009-2019 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
"""
Command line function 'scan'.
"""
import sys

from gimmemotifs.scanner import scan_to_file

MAX_CPUS = 16


def pfmscan(args):

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
        moods=args.moods,
        pvalue=args.pvalue,
        bgfile=args.bgfile,
        genome=args.genome,
        ncpus=args.ncpus,
        zscore=args.zscore,
        gcnorm=args.gcnorm,
    )
