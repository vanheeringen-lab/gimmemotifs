#!/usr/bin/env python
# Copyright (c) 2014 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

import sys
import os
import urllib
import urllib2
import re
import subprocess as sp
from tempfile import NamedTemporaryFile
from distutils.spawn import find_executable

from gimmemotifs.genome_index import GenomeIndex
from gimmemotifs.config import MotifConfig

UCSC_GENOME_URL = "http://hgdownload.soe.ucsc.edu/goldenPath/{0}/bigZips/chromFa.tar.gz"
UCSC_GENE_URL = "http://hgdownload.cse.ucsc.edu/goldenPath/{}/database/"
ANNOS = ["knownGene.txt.gz", "ensGene.txt.gz", "refGene.txt.gz"]

def genome(args):
    
    config = MotifConfig()
    
    if not os.path.exists(args.indexdir):
        print "Index_dir %s does not exist!" % (args.indexdir)
        sys.exit(1)

    if not os.path.exists(args.fastadir):
        print "FASTA dir %s does not exist!" % (args.fastadir)
        sys.exit(1)
    
    pred_bin = "genePredToBed"
    pred = find_executable(pred_bin)
    if not pred:
        sys.stderr.write("{} not found in path!\n".format(pred_bin))
        sys.exit(1)
    
    fa_split_bin = "faSplit"
    fa_split = find_executable(fa_split_bin)
    if not fa_split:
        sys.stderr.write("{} not found in path!\n".format(fa_split_bin))
        sys.exit(1)

    fastadir = args.fastadir
    genomebuild = args.genomebuild
    genome_dir = os.path.join(fastadir, genomebuild)
    index_dir = os.path.join(args.indexdir, args.genomebuild)

    # Check for rights to write to directory

    if not os.path.exists(genome_dir):
        try:
            os.mkdir(genome_dir)
        except:
            sys.stderr.write("Could not create genome dir {}\n".format(genome_dir))
            sys.exit(1)
    
    # Download gene file based on URL + genomebuild
    gene_file = os.path.join(config.get_gene_dir(), "%s.bed" % genomebuild)
    tmp = NamedTemporaryFile(delete=False, suffix=".gz")
    
    anno = []
    f = urllib2.urlopen(UCSC_GENE_URL.format(genomebuild))
    p = re.compile(r'\w+.Gene.txt.gz')
    for line in f.readlines():
        m = p.search(line)
        if m:
            anno.append(m.group(0))

    sys.stderr.write("Retrieving gene annotation for {}\n".format(genomebuild))
    url = ""
    for a in ANNOS:
        if a in anno:
            url = UCSC_GENE_URL.format(genomebuild) + a
            break
    if url:
        urllib.urlretrieve(
                url,
                tmp.name 
                )

        sp.call("zcat {} | cut -f2-11 | {} /dev/stdin {}".format(tmp.name, pred, gene_file), shell=True)

    else: 
        sys.stderr.write("No annotation found!")
  
    # download genome based on URL + genomebuild
    sys.stderr.write("Downloading {} genome\n".format(genomebuild))
    genome_fa = os.path.join(genome_dir, "chromfa.tar.gz".format(genomebuild))
    urllib.urlretrieve(
                UCSC_GENOME_URL.format(genomebuild),
                genome_fa
                )
    
    sys.stderr.write("Unpacking\n")
    cmd = "tar -C {0} -xvzf {1} && rm {1}".format(genome_dir, genome_fa)
    sp.call(cmd, shell=True)

    sys.stderr.write("Creating index\n")
    g = GenomeIndex()
    g = g.create_index(genome_dir, index_dir)
