# Copyright (c) 2009-2021 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
from gimmemotifs.orthologs import motif2factor_from_orthologs


def motif2factors(args):
    motif2factor_from_orthologs(
        new_reference=args.new_reference,
        extra_orthologs_references=args.new_reference,
        tmpdir=args.genomes_dir,
        outdir=args.outdir,
        strategy=args.strategy
    )
