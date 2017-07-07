#!/usr/bin/env python
# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

from __future__ import print_function
import sys
import os

from gimmemotifs.genome_index import GenomeIndex, create_bedtools_fa

def index(args):
    
    if not os.path.exists(args.indexdir):
        print("Index_dir %s does not exist!" % (args.indexdir))
        sys.exit(1)

    fasta_dir = args.fastadir
    index_dir = os.path.join(args.indexdir, args.indexname)

    g = GenomeIndex()
    g.create_index(fasta_dir, index_dir)

    create_bedtools_fa(index_dir, fasta_dir)
