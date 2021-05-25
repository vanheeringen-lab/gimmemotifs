# Copyright (c) 2009-2021 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
from gimmemotifs.orthologs import motif2factor_from_orthologs


def motif2factors(args):
    kwargs = {
        "new_reference": args.new_reference,
        "extra_orthologs_references": args.ortholog_references,
        "tmpdir": args.tmpdir,
        "outdir": args.outdir,
        "strategy": args.strategy,
        "database": args.database,
        "threads": args.threads,
    }
    kwargs = {k: v for k, v in kwargs.items() if v is not None}
    motif2factor_from_orthologs(**kwargs)
