#!/usr/bin/env python
# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

from __future__ import print_function
from gimmemotifs.comparison import MotifComparer
from gimmemotifs.motif import pwmfile_to_motifs, Motif 
from gimmemotifs.plot import match_plot

def match(args):
    sample = dict([(m.id, m) for m in pwmfile_to_motifs(args.pwmfile)])
    db = dict([(m.id, m) for m in pwmfile_to_motifs(args.dbpwmfile)])

    mc = MotifComparer()
    result = mc.get_closest_match(sample.values(), db.values(), "partial", "wic", "mean")

    print("Motif\tMatch\tScore\tP-value")
    for motif, match in result.items():
        pval, pos, orient = mc.compare_motifs(sample[motif], db[match[0]], "partial", "wic", "mean", pval=True)
        print("%s\t%s\t%0.2f\t%0.3e" % (motif, match[0], match[1][0], pval))

    if args.img:
        plotdata = []
        for query, match in result.items():
            motif = sample[query]
            dbmotif = db[match[0]]
            pval, pos, orient = mc.compare_motifs(motif, dbmotif, "partial", "wic", "mean", pval=True)
            
            if orient == -1:
                tmp = dbmotif.id
                dbmotif = dbmotif.rc()
                dbmotif.id = tmp

            if pos < 0:
                tmp = motif.id
                motif = Motif([[0.25,0.25,0.25,0.25]] * -pos + motif.pwm)
                motif.id = tmp
            elif pos > 0:
                tmp = dbmotif.id
                dbmotif = Motif([[0.25,0.25,0.25,0.25]] * pos + dbmotif.pwm)
                dbmotif.id = tmp

            plotdata.append((motif, dbmotif, pval))
            match_plot(plotdata, args.img)
