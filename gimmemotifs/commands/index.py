#!/usr/bin/env python
# Copyright (c) 2009-2013 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

from gimmemotifs.genome_index import GenomeIndex
import sys

def index(args):
    if not os.path.exists(args.indexdir):
    	print "Index_dir %s does not exist!" % (args.indexdir)
    	sys.exit(1)

    fasta_dir = args.fastadir
    index_dir = os.path.join(args.indexdir, args.indexname)

    g = GenomeIndex()
    g = g.create_index(fasta_dir, index_dir)
