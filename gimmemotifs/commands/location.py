# Copyright (c) 2009-2019 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
"""
Command line tool 'location'.

Creates a histogram of motif matches relative to sequence center.
"""
import os
from multiprocessing import Pool

from gimmemotifs.config import MotifConfig
from gimmemotifs.fasta import Fasta
from gimmemotifs.motif import read_motifs
from gimmemotifs.utils import motif_localization


def location(args):
    """
    Creates histrogram of motif location.

    Parameters
    ----------
    args : argparse object
        Command line arguments.
    """
    fastafile = args.fastafile
    pfmfile = args.pfmfile

    lsize = args.size
    if not lsize:
        f = Fasta(fastafile)
        lsize = len(f.items()[0][1])
        f = None

    jobs = []
    motifs = read_motifs(pfmfile)
    ids = [motif.id for motif in motifs]
    if args.ids:
        ids = args.ids.split(",")

    ncpus = int(MotifConfig().get_default_params()["ncpus"])
    pool = Pool(processes=ncpus, maxtasksperchild=1000)
    for motif in motifs:
        if motif.id in ids:
            outfile = os.path.join(f"{motif.id}_histogram")
            jobs.append(
                pool.apply_async(
                    motif_localization, (fastafile, motif, lsize, outfile, args.cutoff)
                )
            )
    pool.close()

    for job in jobs:
        job.get()
    pool.join()
