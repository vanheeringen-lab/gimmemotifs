#!/usr/bin/env python
# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

from gimmemotifs.motif import pwmfile_to_motifs

def logo(args):
    inputfile = args.pwmfile
    
    motifs = pwmfile_to_motifs(inputfile)
    if args.ids:
        ids = args.ids.split(",")
        motifs = [m for m in motifs if m.id in ids]
    
    for motif in motifs:
        motif.to_img(motif.id, fmt="PNG")
