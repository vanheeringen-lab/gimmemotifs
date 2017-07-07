from __future__ import print_function
# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

import sys
import os
from gimmemotifs.fasta import Fasta
import gimmemotifs.background as bg
from gimmemotifs.genome_index import check_genome
from gimmemotifs.config import MotifConfig, BG_TYPES
from gimmemotifs.utils import number_of_seqs_in_file

def background(args):

    inputfile = args.inputfile
    out = args.outputfile
    bg_type = args.bg_type
    outformat = args.outformat.lower()
    length = args.length

    if bg_type not in BG_TYPES:
        print("The argument 'type' should be one of: %s" % (",".join(BG_TYPES)))
        sys.exit(1)

    if outformat == "bed" and bg_type == "random":
        print("Random background can only be generated in FASTA format!")
        sys.exit(1)
        
    if bg_type == "gc" and not inputfile:
        print("need a FASTA formatted input file for background gc")
        sys.exit(1)
    
    # GimmeMotifs configuration for file and directory locations
    config = MotifConfig()
        
    # Genome index location for creation of FASTA files
    index_dir = os.path.join(config.get_index_dir(), args.genome)
    if bg_type in ["gc", "genomic", "promoter"] and outformat == "fasta":
        check_genome(args.genome)

    # Gene definition
    gene_file = os.path.join(config.get_gene_dir(), "%s.bed" % args.genome)
    if bg_type in ["promoter"]:
        if not os.path.exists(gene_file):
            print("Can't find gene definition for %s (%s). See GimmeMotifs documentation on how to add gene files." % (args.genome, gene_file))
            sys.exit(1)
    
    # Number of sequences
    number = None
    if args.number:
        number = args.number
    elif inputfile:
        number = number_of_seqs_in_file(inputfile)
    else:
        sys.stderr.write("please provide either a number or an inputfile\n")
        sys.exit(1)
    
    if bg_type == "random":
        f = Fasta(inputfile)
        m = bg.MarkovFasta(f, n=number, k=args.markov_order)
        m.writefasta(out)
    elif bg_type == "gc":
        if outformat in ["fasta", "fa"]:
            m = bg.MatchedGcFasta(inputfile, args.genome, number=number)
            m.writefasta(out)
        else:
            bg.matched_gc_bedfile(out, inputfile, args.genome, number)
    elif bg_type == "promoter":
        if outformat in ["fasta", "fa"]:
            m = bg.PromoterFasta(gene_file, index_dir, length=length, n=number)
            m.writefasta(out)
        else:
            bg.create_promoter_bedfile(out, gene_file, length, number)
    elif bg_type == "genomic":
        if outformat in ["fasta", "fa"]:
            m = bg.RandomGenomicFasta(index_dir, length, number)
            m.writefasta(out)
        else:
            bg.create_random_genomic_bedfile(out, index_dir, length, number)
        
