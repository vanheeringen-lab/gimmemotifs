from __future__ import print_function

# Copyright (c) 2009-2019 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
from gimmemotifs.background import create_background_file


def background(args):
    create_background_file(
        outfile=args.outputfile,
        bg_type=args.bg_type,
        fmt=args.outformat,
        size=args.size,
        genome=args.genome,
        inputfile=args.inputfile,
        number=args.number,
    )
