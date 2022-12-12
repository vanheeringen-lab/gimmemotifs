#!/usr/bin/env python
# Copyright (c) 2009-2019 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
import logging
import os
import sys

from gimmemotifs.motif import read_motifs
from gimmemotifs.utils import pfmfile_location

logger = logging.getLogger("gimme.logo")


def logo(args):
    if args.pfmfile is None and args.ids is None:
        name = os.path.splitext(os.path.split(pfmfile_location(None))[-1])[0]
        logger.error(
            "Use the -i argument to specify which motif ids you want to use for logos."
        )
        logger.error(
            "If you really want to create logos for all of the motifs in the default"
        )
        logger.error("PFM file use the following command:")
        logger.error(f"gimme logo -p {name}")
        sys.exit(1)
    inputfile = args.pfmfile

    motifs = read_motifs(inputfile)
    if args.ids:
        ids = args.ids.split(",")
        motifs = [m for m in motifs if m.id in ids]

    for motif in motifs:
        motif.plot_logo(fname=f"{motif.id}.png", kind=args.kind, title=args.title)
