#!/usr/bin/env python
# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

import sys
import os
from tempfile import NamedTemporaryFile

import pybedtools

from gimmemotifs.genome_index import GenomeIndex
from gimmemotifs.shutils import find_by_ext
from gimmemotifs.config import FASTA_EXT



def index(args):
    
    if not os.path.exists(args.indexdir):
        print "Index_dir %s does not exist!" % (args.indexdir)
        sys.exit(1)

    fasta_dir = args.fastadir
    index_dir = os.path.join(args.indexdir, args.indexname)

    g = GenomeIndex()
    g.create_index(fasta_dir, index_dir)

    # Create genome FASTA file for use with bedtools
    with open(os.path.join(index_dir, "genome.fa"), 'w') as out:
        for f in find_by_ext(fasta_dir, FASTA_EXT):
            for line in open(f):
                out.write(line)

    test_chr = g.get_chromosomes()[0]
    tmp = NamedTemporaryFile()
    tmp.write("{}\t1\t2\n".format(test_chr))
    tmp.flush()
    
    b = pybedtools.BedTool(tmp.name)
    try:
        b.nucleotide_content(fi=os.path.join(index_dir, "genome.fa"))
    except pybedtools.helpers.BEDToolsError as e:
        if str(e).find("generating") == -1:
            raise

