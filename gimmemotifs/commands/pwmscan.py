#!/usr/bin/python
# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
"""
Command line function 'scan'.
"""
from __future__ import print_function
import os
import re

from gimmemotifs.motif import pwmfile_to_motifs
from gimmemotifs.utils import as_fasta 
from gimmemotifs.scanner import Scanner,scan_it_moods
from gimmemotifs.config import MotifConfig,GM_VERSION

MAX_CPUS = 16

def format_line(seq, seq_id, motif, score, pos, strand, bed=False, 
        seq_p=re.compile(r'([^\s:]+):(\d+)-(\d+)'), strandmap = {-1:"-",1:"+"}): 
    if bed:
        m = seq_p.search(seq_id)
        if m:
            chrom = m.group(1)
            start = int(m.group(2))
            return "{}\t{}\t{}\t{}\t{}\t{}".format(
                    chrom, 
                    start + pos, 
                    start + pos + len(motif) , 
                    motif.id, 
                    score, strandmap[strand])
        else:
            return "{}\t{}\t{}\t{}\t{}\t{}".format(
                    seq_id, 
                    pos, 
                    pos +  len(motif), 
                    motif.id, 
                    score, 
                    strandmap[strand])
    else:
        return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tmotif_name \"{}\" ; motif_instance \"{}\"".format(
            seq_id, 
            "pwmscan", 
            "misc_feature", 
            pos + 1,            # GFF is 1-based
            pos + len(motif), 
            score, 
            strandmap[strand], 
            ".", 
            motif.id, 
            seq[pos: pos + len(motif)]
        )

def scan_table(s, inputfile, fa, motifs, cutoff, bgfile, nreport, scan_rc, pvalue, moods):
    # header
    yield "\t{}".format("\t".join([m.id for m in motifs]))
    table = True
    if moods:
        result_it = scan_it_moods(inputfile, motifs, cutoff, bgfile,  nreport, scan_rc, pvalue, table)
        for seq_id, counts in result_it:
            yield "{}\t{}".format(seq_id, "\t".join([str(x) for x in counts]))
    else:
        # get iterator
        result_it = s.count(fa, nreport, scan_rc)
        # counts table
        for i, counts in enumerate(result_it):
            yield "{}\t{}".format(
                        fa.ids[i], 
                        "\t".join([str(x) for x in counts])
                        )
def scan_score_table(s, fa, motifs, scan_rc):
    
    s.set_threshold(threshold=0.0)
    # get iterator
    result_it = s.best_score(fa, scan_rc)
    # header
    yield "\t{}".format("\t".join([m.id for m in motifs]))
    # score table
    for i,scores in enumerate(result_it):
        yield "{}\t{}".format(
                    fa.ids[i], 
                    "\t".join(["{:4f}".format(x) for x in scores])
                    )

def scan_normal(s, inputfile, fa, motifs, cutoff, bgfile, nreport, scan_rc, pvalue, moods, bed):
    
    table = False
    if moods:
        result_it = scan_it_moods(inputfile, motifs, cutoff, bgfile, nreport, scan_rc, pvalue, table)
        for motif, d in result_it:
            for seq_id,matches in d.items():
                for pos,score,strand in matches:
                    yield format_line(fa[seq_id], seq_id, motif,
                            score, pos, strand, bed=bed)
    else:
        result_it = s.scan(fa, nreport, scan_rc)
        for i, result in enumerate(result_it):
            seq_id = fa.ids[i]
            seq = fa[seq_id]
            for motif, matches in zip(motifs, result):
                for (score, pos, strand) in matches:
                    yield format_line(seq, seq_id, motif, 
                               score, pos, strand, bed=bed)


def command_scan(inputfile, pwmfile, nreport=1, fpr=0.01, cutoff=None, 
        bed=False, scan_rc=True, table=False, score_table=False, moods=False, 
        pvalue=None, bgfile=None, genome=None):
    motifs = pwmfile_to_motifs(pwmfile)
    
    index_dir = None
    if genome is not None:
        index_dir = os.path.join(MotifConfig().get_index_dir(), genome) 
   
    fa = as_fasta(inputfile, index_dir)
    
    # initialize scanner
    s = Scanner()
    s.set_motifs(pwmfile)
    if not score_table:
        s.set_threshold(fpr=fpr, threshold=cutoff, 
            genome=genome, length=fa.median_length(), filename=bgfile)
    
    if table:
        it = scan_table(s, inputfile, fa, motifs, cutoff, bgfile, nreport, scan_rc, pvalue, moods)
    elif score_table:
        it = scan_score_table(s, fa, motifs, scan_rc) 
    else:
        it = scan_normal(s, inputfile, fa, motifs, cutoff, bgfile, nreport, scan_rc, pvalue, moods, bed)
    
    for row in it:
        yield row

def pwmscan(args):

    if args.fpr is None and args.cutoff is None:
        args.fpr = 0.01

    print("# GimmeMotifs version {}".format(GM_VERSION))
    print("# Input: {}".format(args.inputfile))
    print("# Motifs: {}".format(args.pwmfile))
    if args.fpr:
        if args.genome:
            print("# FPR: {} ({})".format(args.fpr, args.genome))
        elif args.bgfile:
            print("# FPR: {} ({})".format(args.fpr, args.bgfile))
    if args.cutoff:
        print("# Threshold: {}".format(args.cutoff))

    for line in command_scan(
            args.inputfile, 
            args.pwmfile, 
            nreport=args.nreport, 
            fpr=args.fpr,
            cutoff=args.cutoff, 
            bed=args.bed, 
            scan_rc=args.scan_rc, 
            table=args.table, 
            score_table=args.score_table,
            moods=args.moods,
            pvalue=args.pvalue,
            bgfile=args.bgfile,
            genome=args.genome,
            ):
        print(line)
