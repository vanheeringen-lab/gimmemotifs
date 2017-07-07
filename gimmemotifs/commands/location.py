# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
"""
Command line tool 'location'.

Creates a histogram of motif matches relative to sequence center.
"""
from gimmemotifs.fasta import Fasta
from gimmemotifs.motif import pwmfile_to_motifs
from gimmemotifs.utils import motif_localization
from gimmemotifs.mp import pool
import os

def location(args):
    """
    Creates histrogram of motif location.

    Parameters
    ----------
    args : argparse object
        Command line arguments.
    """
    fastafile = args.fastafile
    pwmfile = args.pwmfile

    lwidth = args.width
    if not lwidth:
        f = Fasta(fastafile)
        lwidth = len(f.items()[0][1])
        f = None

    jobs = []
    motifs = pwmfile_to_motifs(pwmfile)
    ids = [motif.id for motif in motifs]
    if args.ids:
        ids = args.ids.split(",")

    for motif in motifs:
        if motif.id in ids:
            outfile = os.path.join("%s_histogram" % motif.id)
            jobs.append(
                    pool.apply_async(
                        motif_localization, 
                        (fastafile,motif,lwidth,outfile, args.cutoff)
                        ))
    
    for job in jobs:
        job.get()
