#!/usr/bin/env python 
# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

import sys
import os

from yaml import load

import gimmemotifs.config as cfg
from gimmemotifs.tools import *

def prediction(args):
    tool = args.tool
    infile =  args.infile
    outfile = args.outfile
    paramfile = args.paramfile

    try:
        g = dict([(k.lower(),v) for k,v in globals().items()])
        klass = g[tool.lower()]
    except:
        sys.stderr.write("Tool {0} not found!\n".format(tool))
        raise
    
    t = klass()

    if not t.is_installed():
        sys.stderr.write("Tool {0} not installed!\n".format(tool))
        
    if not t.is_configured():
        sys.stderr.write("Tool {0} not configured!\n".format(tool))
   
    params = {}
    if paramfile:
        params = load(open(paramfile))

    (motifs, stdout, stderr) = t.run(infile, ".", params)

    sys.stderr.write(stderr)
    sys.stdout.write(stdout)

    f = open(outfile, "w")
    for m in motifs:
        f.write("{0}\n".format(m.to_pfm()))
    f.close()
