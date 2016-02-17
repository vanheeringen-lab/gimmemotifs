#!/usr/bin/python
# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

import os
import sys

from gimmemotifs.fasta import Fasta
from gimmemotifs.config import MotifConfig
from gimmemotifs.utils import parse_cutoff
try: 
    from gimmemotifs.mp import pool
except:
    pass

CHUNK = 1000

def scan_fa_with_motif(fo, motif, cutoff, nreport, rc=True):
    #try:
    return motif, motif.pwm_scan_all(fo, cutoff, nreport, scan_rc=rc)
    #except:
    #    e = sys.exc_info()[0]
    #    msg = "Error calculating stats of {0}, error {1}".format(motif.id, e)
    #    sys.stderr.write("{0}\n".format(msg))

def scan_it(infile, motifs, cutoff, nreport=1, rc=True):
    # Get configuration defaults
    config = MotifConfig()
    # Cutoff for motif scanning, only used if a cutoff is not supplied
    default_cutoff = config.get_default_params()['scan_cutoff']
    # Number of CPUs to use
    ncpus =  config.get_default_params()['ncpus']
    
    cutoffs = parse_cutoff(motifs, cutoff, default_cutoff) 
    
    jobs = []
    fa = Fasta(infile)
    motifkey = dict([(m.id, m) for m in motifs])
    
    for motif in motifs:
        for i in range(0, len(fa), CHUNK):
            jobs.append(pool.apply_async(
                                          scan_fa_with_motif,
                                          (fa[i:i + CHUNK],
                                          motif,
                                          cutoffs[motif.id],
                                          nreport,
                                          rc,
                                          )))
    
        while len(jobs) > 10:
            job = jobs.pop(0) 
            motif, result = job.get()
            yield motifkey[motif.id], result

    for job in jobs:
        motif, result = job.get()
        yield motifkey[motif.id], result
 

def scan(infile, motifs, cutoff, nreport=1, it=False):
    # Get configuration defaults
    config = MotifConfig()
    # Cutoff for motif scanning, only used if a cutoff is not supplied
    default_cutoff = config.get_default_params()['scan_cutoff']
    # Number of CPUs to use
    ncpus =  config.get_default_params()['ncpus']
    
    cutoffs = parse_cutoff(motifs, cutoff, default_cutoff) 
    
    total_result = {}
    jobs = []
    fa = Fasta(infile)
    for motif in motifs:
        for i in range(0, len(fa), CHUNK):
            total_result[motif] = {}
            jobs.append(pool.apply_async(
                                          scan_fa_with_motif,
                                          (fa[i:i + CHUNK],
                                          motif,
                                          cutoffs[motif.id],
                                          nreport,
                                          )))
    motifkey = dict([(m.id, m) for m in motifs])
    for job in jobs:
        motif, result = job.get()
        
        total_result[motifkey[motif.id]].update(result)
   
    return total_result

def get_counts(fname, motifs, cutoff):
    counts = {}
    for m in motifs:
        counts.setdefault(m.id, 0)
    result = scan_it(fname, motifs, cutoff)
    for motif, r in result:
        ncount = len([x for x,y in r.items() if len(y) > 0])
        counts[motif.id] += ncount

    return counts

