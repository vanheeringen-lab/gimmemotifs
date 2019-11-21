#!/usr/bin/env python
# Copyright (c) 2009-2019 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
import os
import sys

from gimmemotifs.motif import read_motifs
from gimmemotifs.utils import pfmfile_location


def logo(args):
    if args.pfmfile is None and args.ids is None:
        name = os.path.splitext(os.path.split(pfmfile_location(None))[-1])[0]
        print(
            "Use the -i argument to specify which motif ids you want to use for logos."
        )
        print("If you really want to create logos for all of the motifs in the default")
        print("PFM file use the following command:")
        print(f"gimme logo -p {name}")
        sys.exit(1)
    inputfile = args.pfmfile

    motifs = read_motifs(inputfile)
    if args.ids:
        ids = args.ids.split(",")
        motifs = [m for m in motifs if m.id in ids]

    for motif in motifs:
        motif.plot_logo(
            fname="{}.png".format(motif.id), kind=args.kind, title=args.title
        )
