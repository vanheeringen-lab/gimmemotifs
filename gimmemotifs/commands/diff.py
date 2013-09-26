#!/usr/bin/python
# Copyright (c) 2013 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
import sys
import os
import numpy as np
from gimmemotifs.scan import get_counts
from gimmemotifs.motif import pwmfile_to_motifs
from gimmemotifs.fasta import Fasta
from gimmemotifs.plot import diff_plot

def diff(args):

    infiles = args.inputfiles.split(",")
    bgfile = args.bgfile
    outfile = args.outputfile
    pwmfile = args.pwmfile
    cutoff = args.cutoff
    
    pwms = dict([(m.id, m) for m in pwmfile_to_motifs(pwmfile)])
    motifs = [m for m in pwms.keys()]
    names = [os.path.basename(os.path.splitext(f)[0]) for f in infiles]
    
    # Get background frequencies
    nbg = float(len(Fasta(bgfile).seqs))
    counts = get_counts(bgfile, pwms.values(), cutoff)
    bgfreq = [(counts[m] + 0.01) / nbg for m in motifs]
    
    # Get frequences in input files
    freq = {}
    for fname in infiles:
        counts = get_counts(fname, pwms.values(), cutoff)
        n = float(len(Fasta(fname).seqs))
        freq[fname] = [(counts[m] + 0.01) / n for m in motifs]
    
    freq = np.array([freq[fname] for fname in infiles]).transpose()
    
    for row in freq:
        print freq

    diff_plot(motifs, pwms, names, freq, bgfreq, outfile)
