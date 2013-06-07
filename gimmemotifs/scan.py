#!/usr/bin/python
# Copyright (c) 2009-2013 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

import pp
import os
import sys

from gimmemotifs.fasta import Fasta
from gimmemotifs.config import MotifConfig
from gimmemotifs.utils import parse_cutoff

CHUNK = 1000

def scan_fa_with_motif(fo, motif, cutoff, nreport):
    return motif, motif.pwm_scan_all(fo, cutoff, nreport)

def scan(infile, motifs, cutoff, nreport=1):
    # Get configuration defaults
    config = MotifConfig()
    # Cutoff for motif scanning, only used if a cutoff is not supplied
    default_cutoff = config.get_default_params()['scan_cutoff']
    # Number of CPUs to use
    ncpus =  config.get_default_params()['ncpus']
    
    cutoffs = parse_cutoff(motifs, cutoff, default_cutoff) 
    
    job_server = pp.Server(secret="beetrootsoup")
    if job_server.get_ncpus() > ncpus:
        job_server.set_ncpus(ncpus)
    
    total_result = {}
    jobs = []
    fa = Fasta(infile)
    for motif in motifs:
        for i in range(0, len(fa), CHUNK):
            total_result[motif] = {}
            jobs.append(job_server.submit(
                                          scan_fa_with_motif,
                                          (fa[i:i + CHUNK],
                                          motif,
                                          cutoffs[motif.id],
                                          nreport,
                                          ),
                                          (),()))
    motifkey = dict([(m.id, m) for m in motifs])
    for job in jobs:
        motif, result = job()
        total_result[motifkey[motif.id]].update(result)
    
    return total_result
