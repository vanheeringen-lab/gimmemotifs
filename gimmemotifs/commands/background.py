#!/usr/bin/python
# Copyright (c) 2009-2011 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

import sys
import os
from gimmemotifs.fasta import Fasta
import gimmemotifs.background as bg
from gimmemotifs.config import MotifConfig, BG_TYPES
from gimmemotifs.utils import *


def background(args):

    inputfile = args.inputfile
    out = args.outputfile
    bg_type = args.bg_type
    outformat = args.outformat.lower()

    if not bg_type in BG_TYPES:
        print "The argument 'type' should be one of: %s" % (",".join(BG_TYPES))
        sys.exit(1)

    if outformat == "bed" and bg_type == "random":
        print "Random background can only be generated in FASTA format!"
        sys.exit(1)
        
    # GimmeMotifs configuration for file and directory locations
    config = MotifConfig()
        
    # Genome index location for creation of FASTA files
    index_dir = os.path.join(config.get_index_dir(), args.genome)
    if bg_type in ["matched_genomic", "random_genomic", "random_promoter"] and outformat == "fasta":
        if not os.path.exists(index_dir):
            print "Index for %s does not exist. Has the genome been indexed for use with GimmeMotifs?" % args.genome
            sys.exit(1)
        
    # Gene definition
    gene_file = os.path.join(config.get_gene_dir(), "%s.bed" % args.genome)
    if bg_type in ["matched_genomic", "random_promoter"]:
        if not os.path.exists(gene_file):
            print "Can't find gene definition for %s (%s). See GimmeMotifs documentation on how to add gene files." % (args.genome, gene_file)
            sys.exit(1)
        
    l = args.length
    # Get median length
    if inputfile and os.path.exists(inputfile):
        if not bg_type in ["random"]:
            l = median_bed_len(inputfile)    
    
    if not l:
        print "Can't determine length. Specify either an inputfile or a sequence length"
        sys.exit(1)

    if bg_type == "random":
        f = Fasta(inputfile)
        m = bg.MarkovFasta(f, multiply=args.number, k=args.markov_order)
        m.writefasta(out)
    elif bg_type == "matched_genomic":
        if not is_valid_bedfile(inputfile, 3):
            print "A BED file is required for background type 'matched_genomic'"
            sys.exit()
        
        if outformat == "fasta":
            m = bg.MatchedGenomicFasta(inputfile, gene_file, index_dir, length=l, multiply=args.number)
            m.writefasta(out)
        else:
            bg.create_matched_genomic_bedfile(out, inputfile, gene_file, l, args.number, True)
    elif bg_type == "random_promoter":
        if outformat == "fasta":
            m = bg.PromoterFasta(gene_file, index_dir, length=l, n=args.number)
            m.writefasta(out)
        else:
            bg.create_promoter_bedfile(out, gene_file, l, args.number)
    elif bg_type == "random_genomic":
        if outformat == "fasta":
            m = bg.RandomGenomicFasta(index_dir, l, args.number)
            m.writefasta(out)
        else:
            bg.create_random_genomic_bedfile(out, index_dir, l, args.number)
        
