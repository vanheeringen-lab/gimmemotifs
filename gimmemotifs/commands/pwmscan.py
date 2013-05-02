#!/usr/bin/python
# Copyright (c) 2009-2012 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

import pp
import os
import sys

from gimmemotifs.motif import pwmfile_to_motifs
from gimmemotifs.fasta import Fasta

VERSION = "1.2"
NREPORT = 1 
MAX_CPUS = 16
CHUNK = 1000
DEFAULT_CUTOFF = 0.9

def scan_fa_with_motif(fo, motif, cutoff, nreport, bed):
    return motif, motif.pwm_scan_all(fo, cutoff, nreport)

def pwmscan(args):
    inputfile = args.inputfile
    nreport = args.nreport
    cutoff = args.cutoff
    bed = args.bed

    motifs = pwmfile_to_motifs(args.pwmfile)
    
    cutoffs = {}
    if os.path.isfile(cutoff):
        for i,line in enumerate(open(cutoff)):
            if line != "Motif\tScore\tCutoff\n":
                try:
                    motif,v,c = line.strip().split("\t")
                    c = float(c)
                    cutoffs[motif] = c
                except Exception as e:
                    sys.stderr.write("Error parsing cutoff file, line {0}: {1}\n".format(e, i + 1))
                    sys.exit(1)
    else:
        for motif in motifs:
            cutoffs[motif.id] = float(cutoff)

    for motif in motifs:
        if not cutoffs.has_key(motif.id):
            sys.stderr.write("No cutoff found for {0}, using default {1}\n".format(motif.id, DEFAULT_CUTOFF))
            cutoffs[motif.id] = DEFAULT_CUTOFF

    job_server = pp.Server(secret="beetrootsoup")
    if job_server.get_ncpus() > MAX_CPUS:
        job_server.set_ncpus(MAX_CPUS)

    jobs = []
    fa = Fasta(inputfile)
    for motif in motifs:
        for i in range(0, len(fa), CHUNK):
            jobs.append(job_server.submit(
                                          scan_fa_with_motif, 
                                          (fa[i:i + CHUNK], 
                                          motif, 
                                          cutoffs[motif.id], 
                                          nreport, 
                                          bed
                                         )
                                         , (),()))

    strandmap = {-1:"-",1:"+"}
    
    for job in jobs:
        motif, result = job()
        for seq_id, matches in result.items():
            for (pos, score, strand) in matches:
                if bed:
                    first = seq_id.split(" ")[0]    
                    (chrom,loc) = first.split(":")
                    if loc:
                        (start, end) = map(int, loc.split("-"))
                        print "%s\t%s\t%s\t%s\t%s\t%s" % (chrom, start + pos, start + pos + len(motif) , motif.id, score, strandmap[strand])
                    else:
                        print "%s\t%s\t%s\t%s\t%s\t%s" % (seq_id, pos, pos +  len(motif), motif.id, score, strandmap[strand])
                else:
                    print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tmotif_name \"%s\" ; motif_instance \"%s\"" % (
                        seq_id, 
                        "pwmscan", 
                        "misc_feature", 
                        pos, pos + len(motif) , 
                        score, 
                        strandmap[strand], 
                        ".", 
                        motif.id, 
                        fa[seq_id][pos: pos + len(motif)]
                    )

