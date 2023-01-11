#!/usr/bin/python
# Copyright (c) 2013-2019 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
import logging
import os
import shutil
import sys
from tempfile import mkdtemp

import numpy as np
from genomepy import Genome

from gimmemotifs.fasta import Fasta
from gimmemotifs.motif import read_motifs
from gimmemotifs.plot import diff_plot
from gimmemotifs.scanner import Scanner

logger = logging.getLogger("gimme.diff")


def diff(args):

    infiles = args.inputfiles.split(",")
    bgfile = args.bgfile
    outfile = args.outputfile
    pfmfile = args.pfmfile
    cutoff = args.cutoff
    genome = args.genome
    minenr = float(args.minenr)
    minfreq = float(args.minfreq)

    tmpdir = mkdtemp()

    # Retrieve FASTA clusters from BED file
    if len(infiles) == 1 and infiles[0].endswith("bed"):
        if not args.genome:
            logger.error("Can't convert BED file without genome!")
            sys.exit(1)

        clusters = {}
        with open(infiles[0]) as f:
            for line in f:
                vals = line.strip().split("\t")
                clusters.setdefault(vals[4], []).append(vals[:3])

        infiles = []

        for cluster, regions in clusters.items():
            logger.info(f"Creating FASTA file for {cluster}")
            inbed = os.path.join(tmpdir, f"{cluster}.bed")
            outfa = os.path.join(tmpdir, f"{cluster}.fa")
            with open(inbed, "w") as f:
                for r in regions:
                    f.write(f"{r[0]}\t{r[1]}\t{r[2]}\n")
            Genome(genome).track2fasta(inbed, outfa)
            infiles.append(outfa)

    pfms = dict([(m.id, m) for m in read_motifs(pfmfile)])
    motifs = [m for m in pfms.keys()]
    names = [os.path.basename(os.path.splitext(fname)[0]) for fname in infiles]

    s = Scanner()
    s.set_motifs(pfmfile)
    s.set_threshold(threshold=cutoff)

    # Get background frequencies
    nbg = float(len(Fasta(bgfile).seqs))

    bgcounts = s.total_count(bgfile, nreport=1)
    bgfreq = [(c + 0.01) / nbg for c in bgcounts]

    # Get frequences in input files
    freq = {}
    counts = {}
    for fname in infiles:
        mcounts = s.total_count(fname, nreport=1)
        n = float(len(Fasta(fname).seqs))
        counts[fname] = mcounts
        freq[fname] = [(c + 0.01) / n for c in mcounts]

    freq = np.array([freq[fname] for fname in infiles]).transpose()
    counts = np.array([counts[fname] for fname in infiles]).transpose()

    # for row in freq:
    #    print freq

    diff_plot(
        motifs,
        pfms,
        names,
        freq,
        counts,
        bgfreq,
        bgcounts,
        outfile,
        minenr=minenr,
        minfreq=minfreq,
    )

    shutil.rmtree(tmpdir)
