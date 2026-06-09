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
