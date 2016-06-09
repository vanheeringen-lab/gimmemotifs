#!/usr/bin/python
# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

import os
import sys
import re

from gimmemotifs.motif import pwmfile_to_motifs
from gimmemotifs.utils import parse_cutoff, as_fasta 
from gimmemotifs.scanner import Scanner,scan_it_moods
from gimmemotifs.config import MotifConfig

MAX_CPUS = 16

def format_line(fa, seq_id, motif, score, pos, strand, bed=False):                

    strandmap = {-1:"-",1:"+"}
    p = re.compile(r'([^\s:]+):(\d+)-(\d+)')
    if bed:
        m = p.search(seq_id)
        if m:
            chrom = m.group(1)
            start = int(m.group(2))
            end = int(m.group(3))
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
            fa[seq_id][pos: pos + len(motif)]
        )

def command_scan(inputfile, pwmfile, nreport=1, cutoff=0.9, bed=False, 
        scan_rc=True, table=False, score_table=False, moods=False, 
        pvalue=None, bgfile=None, genome=None):
    motifs = pwmfile_to_motifs(pwmfile)
    
    index_dir = None
    if genome is not None:
        index_dir = os.path.join(MotifConfig().get_index_dir(), genome) 
    
    # initialize scanner
    s = Scanner()
    s.set_motifs(pwmfile)
    
    fa = as_fasta(inputfile, index_dir)
    
    if moods:
        result_it = scan_it_moods(inputfile, motifs, cutoff, bgfile, nreport, scan_rc, pvalue, table)
    else:
        result_it = s.scan(fa, nreport, scan_rc, cutoff)

    
    if table:
        # header
        yield "\t{}".format("\t".join([m.id for m in motifs]))
        
        if moods:
            result_it = scan_it_moods(inputfile, motifs, cutoff, bgfile,  nreport, scan_rc, pvalue, table)
            for seq_id, counts in result_it:
                yield "{}\t{}".format(seq_id, "\t".join([str(x) for x in counts]))
        else:
            # get iterator
            result_it = s.count(fa, nreport, scan_rc, cutoff)
            # counts table
            for i, counts in enumerate(result_it):
                yield "{}\t{}".format(
                        fa.ids[i], 
                        "\t".join([str(x) for x in counts])
                        )

    elif score_table:
        # get iterator
        result_it = s.best_score(fa, scan_rc)
        # header
        yield "\t{}".format("\t".join([m.id for m in motifs]))
        # score table
        for i,scores in enumerate(result_it):
            yield "{}\t{}".format(
                    fa.ids[i], 
                    "\t".join([str(x) for x in scores])
                    )

    else:
        if moods:
            for motif, d in result_it:
                for seq_id,matches in d.items():
                    for pos,score,strand in matches:
                        yield format_line(fa, seq_id, motif,
                                score, pos, strand, bed=bed)
        else:
            for i, result in enumerate(result_it):
                seq_id = fa.ids[i]
                for motif, matches in zip(motifs, result):
                    for (score, pos, strand) in matches:
                        yield format_line(fa, seq_id, motif, 
                                   score, pos, strand, bed=bed)

def pwmscan(args):
    inputfile = args.inputfile
    nreport = args.nreport
    cutoff = args.cutoff
    bed = args.bed
    scan_rc = args.scan_rc
    table = args.table
    score_table = args.score_table
    pwmfile = args.pwmfile
    bgfile = args.bgfile
    moods = args.moods

    for line in command_scan(
            inputfile, 
            pwmfile, 
            nreport, 
            cutoff, 
            bed, 
            scan_rc, 
            table, 
            score_table,
            moods,
            args.pvalue,
            bgfile,
            genome=args.genome,
            ):
        print line
