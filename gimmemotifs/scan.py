#!/usr/bin/python
# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

import os
import sys
from tempfile import mkdtemp

import MOODS.tools
import MOODS.parsers
import MOODS.scan
import numpy as np

from gimmemotifs.fasta import Fasta
from gimmemotifs.genome_index import rc
from gimmemotifs.config import MotifConfig
from gimmemotifs.utils import parse_cutoff
from gimmemotifs.motif import read_motifs
from gimmemotifs.utils import which
from gimmemotifs.c_metrics import pwmscan

try: 
    from gimmemotifs.mp import pool
except:
    pass

CHUNK = 1000

def scan_fasta_file_with_motifs(fastafile, motiffile, threshold, gfffile, scan_rc=True, nreport=1):
    error = None
    try:
        motifs = read_motifs(open(motiffile), fmt="pwm")
        fa = Fasta(fastafile)
        for motif in motifs:
            motif.pwm_scan_to_gff(fa, gfffile, nreport=nreport, cutoff=float(threshold), scan_rc=scan_rc, append=True)
    except Exception,e :
        error = e
    return error

def scan_fa_with_motif(fo, motif, cutoff, nreport, scan_rc=True):
    #try:
    return motif, motif.pwm_scan_all(fo, cutoff, nreport, scan_rc=scan_rc)
    #except:
    #    e = sys.exc_info()[0]
    #    msg = "Error calculating stats of {0}, error {1}".format(motif.id, e)
    #    sys.stderr.write("{0}\n".format(msg))

def scan_fa_with_motif_counts(fo, motifs, cutoffs, nreport, scan_rc=True):
    ret = []
    for name,seq in fo.items():
        counts = []
        for motif,cutoff in zip(motifs, cutoffs):
            result = pwmscan(seq.upper(), motif.pwm, cutoff, nreport, scan_rc)
            counts.append(len(result))
        ret.append((name, counts))
    return ret

def scan_fa_with_motif_score(fo, motifs, cutoffs, scan_rc=True):
    ret = []
    for name,seq in fo.items():
        scoress = []
        for motif,cutoff in zip(motifs, cutoffs):
            result = pwmscan(seq.upper(), motif.pwm, cutoff, 1, scan_rc)
            scores.append(result[0][0])
        ret.append((name, scores))
    return ret

def scan_fa_with_motif_moods(fo, motifs, matrices, bg, thresholds, nreport, scan_rc=True):
    
    scanner = MOODS.scan.Scanner(7)
    scanner.set_motifs(matrices, bg, thresholds)
    
    ret = []
    for name, seq in fo.items():
        l = len(seq)
        
        scan_seq = seq.upper()
        if scan_rc:
            scan_seq = "".join((scan_seq, "N"*50, rc(scan_seq)))
        results = scanner.scan_max_hits(scan_seq, nreport)
        for motif,result in zip(motifs, results):
            matches = []
            for match in result:
                strand = 1
                pos = match.pos
                if scan_rc:
                    if pos > l:
                        pos = l - (pos - l - 50) - len(motif) 
                        strand = -1
                matches.append((pos, match.score, strand))
            ret.append((motif, {name: matches}))

    return ret

def scan_fa_with_motif_moods_count(fo, motifs, matrices, bg, thresholds, nreport, scan_rc=True):
    scanner = MOODS.scan.Scanner(7)
    scanner.set_motifs(matrices, bg, thresholds)
  
    ret = []
    for name, seq in fo.items():
        l = len(seq)
        
        scan_seq = seq.upper()
        if scan_rc:
            scan_seq = "".join((scan_seq, "N"*50, rc(scan_seq)))
        results = scanner.counts_max_hits(scan_seq, nreport)
        ret.append((name, results))

    return ret

def scan_it(infile, motifs, cutoff, nreport=1, scan_rc=True):
    # Get configuration defaults
    config = MotifConfig()
    # Cutoff for motif scanning, only used if a cutoff is not supplied
    default_cutoff = config.get_default_params()['scan_cutoff']
    # Number of CPUs to use
    
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
                                          scan_rc,
                                          )))
    
        while len(jobs) > 10:
            job = jobs.pop(0) 
            motif, result = job.get()
            yield motifkey[motif.id], result

    for job in jobs:
        motif, result = job.get()
        yield motifkey[motif.id], result

def get_cutoffs(motifs, cutoff):
    # Get configuration defaults
    config = MotifConfig()
    # Cutoff for motif scanning, only used if a cutoff is not supplied
    default_cutoff = config.get_default_params()['scan_cutoff']
    
    cutoffs = parse_cutoff(motifs, cutoff, default_cutoff) 
    
    return [cutoffs[m.id] for m in motifs]

def scan_it_counts_table(infile, motifs, cutoff, nreport=1, scan_rc=True):
    cutoffs = get_cutoffs(motifs, cutoff)

    fa = Fasta(infile)
        
    jobs = []
    for i in range(0, len(fa), CHUNK):
        jobs.append(pool.apply_async(
                                          scan_fa_with_motif_counts,
                                          (fa[i:i + CHUNK],
                                          motifs,
                                          cutoffs,
                                          nreport,
                                          scan_rc,
                                          )))
    
    for job in jobs:
        for ret in job.get():
            yield ret

def scan_it_score_table(infile, motifs, cutoff, scan_rc=True):
    cutoffs = get_cutoffs(motifs, cutoff)

    fa = Fasta(infile)
        
    jobs = []
    for i in range(0, len(fa), CHUNK):
        jobs.append(pool.apply_async(
                                          scan_fa_with_motif_score,
                                          (fa[i:i + CHUNK],
                                          motifs,
                                          cutoffs,
                                          scan_rc,
                                          )))
    
    for job in jobs:
        ret = job.get()
        yield ret


def calc_threshold_moods(m, c):
    m_min = MOODS.tools.min_score(m)
    m_max = MOODS.tools.max_score(m)

    return m_min + (m_max - m_min) * c

def scan_it_moods(infile, motifs, cutoff, nreport=1, scan_rc=True, pvalue=None, count=False):
    tmpdir = mkdtemp()
    matrices = []
    pseudocount = 1e-3
    bgfile = "/home/simon/prj/atac_zebrafish/zf/random_w500.fa"
    bg = MOODS.tools.bg_from_sequence_dna("".join(Fasta(bgfile).seqs), 1)

    for motif in motifs:
        pfmname = os.path.join(tmpdir, "{}.pfm".format(motif.id))
        with open(pfmname, "w") as f:
            matrix = np.array(motif.pwm).transpose()
            for line in [" ".join([str(x) for x in row]) for row in matrix]:
                f.write("{}\n".format(line))

        matrices.append(MOODS.parsers.pfm_log_odds(pfmname, bg, pseudocount))
    
    thresholds = []
    if pvalue is not None:
        thresholds = [MOODS.tools.threshold_from_p(m, bg, float(pvalue)) for m in matrices]
        #sys.stderr.write("{}\n".format(thresholds))
    else:
        thresholds = [calc_threshold_moods(m, float(cutoff)) for m in matrices]
    
    scanner = MOODS.scan.Scanner(7)
    scanner.set_motifs(matrices, bg, thresholds)
    
    config = MotifConfig()
    ncpus =  int(config.get_default_params()['ncpus'])
    fa = Fasta(infile)
    chunk = 500
    if (len(fa) / chunk) < ncpus:
        chunk = len(fa) / (ncpus + 1)
    
    jobs = []
    func = scan_fa_with_motif_moods
    if count:
        func = scan_fa_with_motif_moods_count

    for i in range(0, len(fa), chunk):
        jobs.append(pool.apply_async(
                                          func,
                                          (fa[i:i + chunk],
                                          motifs,
                                          matrices,
                                          bg,
                                          thresholds,
                                          nreport,
                                          scan_rc,
                                          )))
    
    for job in jobs:
        for ret in job.get():
            yield ret

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

