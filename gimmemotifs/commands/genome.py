#!/usr/bin/env python
# Copyright (c) 2014-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
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
from glob import glob
from tempfile import NamedTemporaryFile
from distutils.spawn import find_executable
from gimmemotifs.genome_index import GenomeIndex
from gimmemotifs.config import MotifConfig
from gimmemotifs.fasta import Fasta

UCSC_GENOME_URL = "http://hgdownload.soe.ucsc.edu/goldenPath/{0}/bigZips/chromFa.tar.gz"
ALT_UCSC_GENOME_URL = "http://hgdownload.soe.ucsc.edu/goldenPath/{0}/bigZips/{0}.fa.gz"
UCSC_GENE_URL = "http://hgdownload.cse.ucsc.edu/goldenPath/{}/database/"
ANNOS = ["knownGene.txt.gz", "ensGene.txt.gz", "refGene.txt.gz"]


def check_genome_file(fname):
    if os.path.exists(fname):
        with open(fname) as f:
            for i in range(10):
                line = f.readline()
                if line.find("Not Found") != -1:
                    return False
            return True
    return False

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
    for genome_url in [UCSC_GENOME_URL, ALT_UCSC_GENOME_URL]:
        
        remote = genome_url.format(genomebuild)

        genome_fa = os.path.join(
                genome_dir, 
                os.path.split(remote)[-1]
                )

        sys.stderr.write("Trying to download {}\n".format(genome_url.format(genomebuild)))
        urllib.urlretrieve(
                genome_url.format(genomebuild),
                genome_fa
                )
        
        if not check_genome_file(genome_fa):    
            continue
        
        break

    if not check_genome_file(genome_fa):
        sys.stderr.write("Failed to download genome\n")
        sys.exit(1)

    sys.stderr.write("Unpacking\n")
    if genome_fa.endswith("tar.gz"):
        cmd = "tar -C {0} -xvzf {1} && rm {1}".format(genome_dir, genome_fa)
    else:
        cmd = "gunzip {0} && rm {0}".format(genome_fa)

    sp.call(cmd, shell=True, cwd=genome_dir)

    fa_files = glob("{}/*.fa".format(genome_dir))
    if len(fa_files) == 1:
        f = Fasta(fa_files[0])
        for n,s in f.items():
            with open("{}/{}.fa".format(n)) as f:
                f.write("{}\n{}\n".format(n,s))
    
        os.unlink(fa_files[0])

    sys.stderr.write("Creating index\n")
    g = GenomeIndex()
    g = g.create_index(genome_dir, index_dir)
