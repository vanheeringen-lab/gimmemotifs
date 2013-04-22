# Copyright (c) 2009-2013 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

import sys
from gimmemotifs.fasta import Fasta
from gimmemotifs.motif import pwmfile_to_motifs
from gimmemotifs.c_metrics import pwmscan
from scipy.stats import scoreatpercentile

def threshold(args):
    if args.fdr < 0 or args.fdr > 1:
        print "Please specify a FDR between 0 and 1"
        sys.exit(1)

    f = Fasta(args.inputfile)
    motifs = pwmfile_to_motifs(args.pwmfile)

    print "Motif\tScore\tCutoff"
    for motif in motifs:
        pwm = motif.pwm
        scores = []
        min_score = motif.pwm_min_score()
        for name,seq in f.items():
            result = pwmscan(seq.upper(), pwm, min_score, 1, True)
            score = result[0][0]
            scores.append(score)
        opt_score = scoreatpercentile(scores, 100 - (100 * args.fdr))
        cutoff = (opt_score - min_score) / (motif.pwm_max_score() - min_score)
        print "%s\t%s\t%s" % (motif.id, opt_score , cutoff)
