#!/usr/bin/python
# Copyright (c) 2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
import os
import sys

import pandas as pd

from gimmemotifs.config import MotifConfig
from gimmemotifs.motif import read_motifs
from gimmemotifs.scanner import Scanner

def command_moap(inputfile, motiffile=None, pwmfile=None, 
                        genome=None, method="classic"):
    
    config = MotifConfig()

    m2f = None
    
    # read data
    df = pd.read_table(inputfile, index_col=0)
    
    if motiffile is None:
        if genome is None:
            raise ValueError("need a genome")
        # check pwmfile
        if pwmfile is None:
            pwmfile = config.get_default_params().get("motif_db", None)
            if pwmfile is not None:
                pwmfile = os.path.join(config.get_motif_dir(), pwmfile)
        
        if pwmfile is None:
            raise ValueError("no pwmfile given and no default database specified")

        if not os.path.exists(pwmfile):
            raise ValueError("{} does not exist".format(pwmfile))

        try:
            motifs = read_motifs(open(pwmfile))
        except:
            sys.stderr.write("can't read motifs from {}".format(pwmfile))
            raise

        base = os.path.splitext(pwmfile)[0]
        map_file = base + ".motif2factors.txt"
        if os.path.exists(map_file):
            m2f = pd.read_table(map_file, index_col=0)

        # initialize scanner
        s = Scanner()
        sys.stderr.write(pwmfile + "\n")
        s.set_motifs(pwmfile)
        s.set_genome(genome)

        # scan for motifs
        sys.stderr.write("scanning for motifs\n")
        motif_names = [m.id for m in read_motifs(open(pwmfile))]
        scores = []
        for row in s.best_score(list(df.index)):
            scores.append(row)
        motifs = pd.DataFrame(scores, index=df.index, columns=motif_names)
    else:
        motifs = pd.DataFrame(motiffile)   

    
    clf = None
    if method == "ks":
        clf = KSMoap()
    if method == "ks":
        clf = KSMoap()

    clf.fit(motifs, df)
    print clf.act_

if __name__ == "__main__":
    infile = "/home/simon/git/gimmemotifs/test_rpkm_table.mm10.txt"
    command_moap(infile, genome="mm10", method="moap")

