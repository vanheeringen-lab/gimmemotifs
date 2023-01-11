# Copyright (c) 2009-2019 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
"""Command line function 'threshold'"""
import logging
import sys

from gimmemotifs.motif import read_motifs
from gimmemotifs.scanner import Scanner

logger = logging.getLogger("gimme.threshold")


def threshold(args):
    """Calculate motif score threshold for a given FPR."""
    if args.fpr < 0 or args.fpr > 1:
        logger.error("Please specify a FPR between 0 and 1")
        sys.exit(1)

    motifs = read_motifs(args.pfmfile)

    s = Scanner()
    s.set_motifs(args.pfmfile)
    s.set_background(fname=args.inputfile)
    s.set_threshold(args.fpr)

    print("Motif\tScore\tCutoff")
    for motif in motifs:
        min_score = motif.min_score
        max_score = motif.max_score
        opt_score = s.threshold[motif.id]
        if opt_score is None:
            opt_score = motif.max_score
        threshold = (opt_score - min_score) / (max_score - min_score)
        print(f"{motif.id}\t{opt_score}\t{threshold}")
