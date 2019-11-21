#!/usr/bin/env python
# Copyright (c) 2009-2019 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.

from __future__ import print_function
from gimmemotifs.comparison import MotifComparer
from gimmemotifs.motif import read_motifs, Motif
from gimmemotifs.plot import match_plot


def match(args):
    sample = dict([(m.id, m) for m in read_motifs(args.pfmfile)])
    db = dict([(m.id, m) for m in read_motifs(args.dbpfmfile)])

    mc = MotifComparer()
    result = mc.get_best_matches(
        sample.values(), args.nmatches, db.values(), "partial", "seqcor", "mean"
    )

    plotdata = []
    print("Motif\tMatch\tScore\tP-value")
    for motif_name, matches in result.items():
        for match in matches:

            pval, pos, orient = mc.compare_motifs(
                sample[motif_name], db[match[0]], "partial", "seqcor", "mean", pval=True
            )
            print("%s\t%s\t%0.2f\t%0.3e" % (motif_name, match[0], match[1][0], pval))
            motif = sample[motif_name]
            dbmotif = db[match[0]]

            if args.img:
                if orient == -1:
                    tmp = dbmotif.id
                    dbmotif = dbmotif.rc()
                    dbmotif.id = tmp
                if pos < 0:
                    tmp = motif.id
                    motif = Motif([[0.25, 0.25, 0.25, 0.25]] * -pos + motif.pwm)
                    motif.id = tmp
                elif pos > 0:
                    tmp = dbmotif.id
                    dbmotif = Motif([[0.25, 0.25, 0.25, 0.25]] * pos + dbmotif.pwm)
                    dbmotif.id = tmp

                diff = len(motif) - len(dbmotif)
                if diff > 0:
                    dbmotif = Motif(dbmotif.pwm + [[0.25, 0.25, 0.25, 0.25]] * diff)
                else:
                    motif = Motif(motif.pwm + [[0.25, 0.25, 0.25, 0.25]] * -diff)

                plotdata.append((motif, dbmotif, pval))
    if args.img:
        match_plot(plotdata, args.img)
