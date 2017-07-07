#!/usr/bin/env python
# Copyright (c) 2014-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

from __future__ import print_function
import sys
import os
from gimmemotifs.genome_index import get_genome

def genome(args):
    
    if not os.path.exists(args.indexdir):
        print("Index_dir %s does not exist!" % (args.indexdir))
        sys.exit(1)

    if not os.path.exists(args.fastadir):
        print("FASTA dir %s does not exist!" % (args.fastadir))
        sys.exit(1)

    get_genome(args.genomebuild, args.fastadir, args.indexdir)
