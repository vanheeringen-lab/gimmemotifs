#!/usr/bin/python
# Copyright (c) 2009-2013 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

import sys
import re

from gimmemotifs.fasta import Fasta
from gimmemotifs.motif import pwmfile_to_motifs
from gimmemotifs.utils import parse_cutoff 
from gimmemotifs.scan import scan, scan_it

VERSION = "1.2"
MAX_CPUS = 16

def pwmscan(args):
    inputfile = args.inputfile
    nreport = args.nreport
    cutoff = args.cutoff
    bed = args.bed

    motifs = pwmfile_to_motifs(args.pwmfile)
    result = scan_it(inputfile, motifs, cutoff, nreport)
   
    p = re.compile(r'([^\s:]+):(\d+)-(\d+)')
    fa = Fasta(inputfile)
    if args.table:
        table = {}
        for seq_id in fa.ids:
            table[seq_id] = {}

        for motif, result in result:
            for seq_id, matches in result.items():
                table[seq_id][motif] = len(matches)
        
        #mnames = [m.id for m in motifs]
        #print table
        print "\t{}".format("\t".join([m.id for m in motifs]))
        for seq_id in fa.ids:
            counts = [table[seq_id].get(m, 0) for m in motifs]
            print "{}\t{}".format(seq_id, "\t".join([str(x) for x in counts]))

    else:
        strandmap = {-1:"-",1:"+"}
        for motif, result in result:
            for seq_id, matches in result.items():
                for (pos, score, strand) in matches:
                    if bed:
                        m = p.search(seq_id)
                        if m:
                            chrom = m.group(1)
                            start = int(m.group(2))
                            end = int(m.group(3))
                            print "%s\t%s\t%s\t%s\t%s\t%s" % (chrom, start + pos, start + pos + len(motif) , motif.id, score, strandmap[strand])
                        else:
                            print "%s\t%s\t%s\t%s\t%s\t%s" % (seq_id, pos, pos +  len(motif), motif.id, score, strandmap[strand])
                    else:
                        print "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tmotif_name \"%s\" ; motif_instance \"%s\"" % (
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
    
