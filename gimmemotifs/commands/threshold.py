# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

import sys

from scipy.stats import scoreatpercentile
import numpy as np

from gimmemotifs.fasta import Fasta
from gimmemotifs.motif import pwmfile_to_motifs
from gimmemotifs.scanner import Scanner

def threshold(args):
    if args.fdr < 0 or args.fdr > 1:
        print "Please specify a FDR between 0 and 1"
        sys.exit(1)

    motifs = pwmfile_to_motifs(args.pwmfile)
    
    s = Scanner()
    s.set_motifs(args.pwmfile)
    
    score_table = []
    for scores in s.best_score(args.inputfile):
        score_table.append(scores)

    print "Motif\tScore\tCutoff"
    for i,scores in enumerate(np.array(score_table).transpose()):
        motif = motifs[i]
        pwm = motif.pwm
        min_score = motif.pwm_min_score()
        if len(scores) > 0:
            opt_score = scoreatpercentile(scores, 100 - (100 * args.fdr))
            cutoff = (opt_score - min_score) / (
                    motif.pwm_max_score() - min_score)
            print "{0}\t{1}\t{2}".format(
                    motif.id, opt_score , cutoff)
        else:
            sys.stderr.write("Warning: no matches for {0}\n".format(motif.id))
            
