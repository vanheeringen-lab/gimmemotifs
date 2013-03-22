# Copyright (c) 2009-2013 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

from gimmemotifs.fasta import Fasta
from gimmemotifs.motif import pwmfile_to_motifs
from gimmemotifs.utils import motif_localization
import pp
import sys

def location(args):
    fastafile = args.fastafile
    pwmfile = args.pwmfile

    lwidth = args.width
    if not lwidth:
        f = Fasta(fastafile)
        lwidth = len(f.items()[0][1])
        f = None

    job_server = pp.Server(secret="pumpkinrisotto")
    jobs = []
    motifs = pwmfile_to_motifs(pwmfile)
    ids = [motif.id for motif in motifs]
    if args.ids:
        ids = args.ids.split(",")

    for motif in motifs:
        if motif.id in ids:
            outfile = os.path.join("%s_histogram" % motif.id)
            jobs.append(job_server.submit(motif_localization, (fastafile,motif,lwidth,outfile, args.cutoff), (),()))
    
    for job in jobs:
        job()
