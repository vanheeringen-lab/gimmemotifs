from gimmemotifs.orthologs import motif2factor_from_orthologs


def motif2factors(args):
    kwargs = {
        "new_reference": args.new_reference,
        "extra_orthologs_references": args.ortholog_references,
        "genomes_dir": args.genomes_dir,
        "tmpdir": args.tmpdir,
        "outdir": args.outdir,
        "strategy": args.strategy,
        "database": args.database,
        "threads": args.threads,
        "keep_intermediate": args.keep_intermediate,
    }
    kwargs = {k: v for k, v in kwargs.items() if v is not None}
    motif2factor_from_orthologs(**kwargs)
