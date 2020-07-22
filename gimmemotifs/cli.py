#!/usr/bin/env python
# Copyright (c) 2013-2019 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
import os
import sys
import argparse
from gimmemotifs.config import MotifConfig, BG_TYPES, BED_VALID_BGS
from gimmemotifs import commands, __version__
from gimmemotifs.utils import check_genome


def cli(sys_args):
    config = MotifConfig()
    params = config.get_default_params()
    default_pfm_file = os.path.join(config.get_motif_dir(), params["motif_db"])
    default_pfm = params["motif_db"]

    description = """
    GimmeMotifs v{0}
    """.format(
        __version__
    )

    epilog = """
    commands:
        motifs      identify enriched motifs (known and/or de novo)
        scan        scan for known motifs
        maelstrom   find differential motifs
        match       find motif matches in database
        logo        create sequence logo(s)
        cluster     cluster similar motifs
        background  create a background file
        threshold   calculate motif scan threshold
        location    motif location histograms
        diff        compare motif frequency and enrichment
                    between fasta files

    type `gimme <command> -h` for more details
    """

    usage = "%(prog)s [-h] <subcommand> [options]"

    parser = argparse.ArgumentParser(
        usage=usage,
        description=description,
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    subparsers = parser.add_subparsers()  # title='subcommands', metavar="<command>")

    # gimme_motifs.py
    p = subparsers.add_parser("motifs")

    p.add_argument(
        "sample", help="FASTA, BED, narrowPeak or region file.", metavar="INPUT"
    )
    p.add_argument("outdir", metavar="OUTDIR", help="Output directory.")
    p.add_argument(
        "-b",
        "--background",
        help=(
            "Background type ({}) or a file with background "
            "sequences (FASTA, BED or regions)"
        ).format(",".join(BED_VALID_BGS)),
        metavar="BACKGROUND",
    )
    p.add_argument(
        "-g", dest="genome", help="Genome name or fasta file", metavar="GENOME"
    )
    p.add_argument(
        "--denovo",
        dest="known",
        help="Only use de novo motifs",
        default=True,
        action="store_false",
    )
    p.add_argument(
        "--known",
        dest="denovo",
        help="Only use known motifs",
        default=True,
        action="store_false",
    )
    p.add_argument(
        "--noreport",
        dest="report",
        help="Don't create a HTML report.",
        default=True,
        action="store_false",
    )
    p.add_argument(
        "--rawscore",
        dest="zscore",
        help="Don't z-score normalize motif scores",
        action="store_false",
        default=True,
    )
    p.add_argument(
        "--nogc",
        dest="gc",
        help="Don't use GC%% bins",
        action="store_false",
        default=True,
    )
    p.add_argument(
        "-N",
        "--nthreads",
        dest="ncpus",
        help="Number of threads (default %s)" % (params["ncpus"]),
        metavar="INT",
        type=int,
        default=int(params["ncpus"]),
    )

    # Specific arguments for known motifs
    known_grp = p.add_argument_group(title="optional arguments for known motifs")
    known_grp.add_argument(
        "-p",
        dest="pfmfile",
        help="PFM file with motifs." "(default: {0})".format(default_pfm),
        default=default_pfm_file,
        metavar="PFMFILE",
    )

    # Specific arguments for de novo motifs
    denovo_grp = p.add_argument_group(title="optional arguments for de novo motifs")
    denovo_grp.add_argument(
        "-t",
        "--tools",
        dest="tools",
        help="Tools to use, any combination of %s (default %s)"
        % (params["available_tools"], params["tools"]),
        metavar="N",
        default="MEME,Homer,BioProspector",
    )
    denovo_grp.add_argument(
        "-a",
        "--analysis",
        dest="analysis",
        help="Analysis type: small, medium, large, xl (xl)",
        metavar="ANALYSIS",
        default="xl",
    )
    denovo_grp.add_argument(
        "-k",
        "--keepintermediate",
        dest="keep_intermediate",
        help="Don't delete intermediate files",
        default=False,
        action="store_true",
    )
    denovo_grp.add_argument(
        "-S",
        "--singlestrand",
        dest="single",
        help="Only predict motifs for single + strand (default is both)",
        action="store_true",
        default=False,
    )
    denovo_grp.add_argument(
        "-f",
        "--fraction",
        dest="fraction",
        help="Fraction of peaks to use for motif predicton (%s)" % params["fraction"],
        metavar="FRACTION",
        default=params["fraction"],
        type=float,
    )
    denovo_grp.add_argument(
        "-s",
        "--size",
        dest="size",
        help=(
            "Region size to use for motif prediction ({}). "
            "Set to 0 to use the size of the input regions."
        ).format(params["size"]),
        metavar="N",
        default=params["size"],
        type=int,
    )

    p.set_defaults(func=commands.motifs)

    # pfmscan.py
    NREPORT = 1
    p = subparsers.add_parser("scan")
    p.add_argument(
        "inputfile", help="inputfile (FASTA, BED, regions)", metavar="INPUTFILE"
    )
    p.add_argument(
        "-g", "--genome", dest="genome", help="Genome", metavar="GENOME", default=None
    )
    p.add_argument(
        "-p",
        "--pfmfile",
        dest="pfmfile",
        help="PFM file with motifs " "(default: {0})".format(default_pfm),
        default=default_pfm_file,
        metavar="pfmfile",
    )
    p.add_argument(
        "-f",
        "--fpr",
        dest="fpr",
        help="FPR for motif scanning (default 0.01)",
        metavar="",
        default=None,
    )
    p.add_argument(
        "-B",
        "--bgfile",
        dest="bgfile",
        help="background file for threshold",
        metavar="",
        default=None,
    )
    p.add_argument(
        "-c",
        "--cutoff",
        dest="cutoff",
        help="motif score cutoff or file with cutoffs",
        metavar="",
        default=None,
    )
    p.add_argument(
        "-n",
        "--nreport",
        dest="nreport",
        help="report the N best matches",
        metavar="N",
        default=NREPORT,
        type=int,
    )
    p.add_argument(
        "-r",
        "--norc",
        dest="scan_rc",
        help="don't scan reverse complement (- strand)",
        default=True,
        action="store_false",
    )
    p.add_argument(
        "-b",
        "--bed",
        action="store_true",
        dest="bed",
        default=False,
        help="output bed format",
    )
    p.add_argument(
        "-t",
        "--table",
        dest="table",
        help="output counts in tabular format",
        action="store_true",
        default=False,
    )
    p.add_argument(
        "-T",
        "--score_table",
        dest="score_table",
        help="output maximum score in tabular format",
        action="store_true",
        default=False,
    )
    p.add_argument(
        "-z",
        "--zscore",
        dest="zscore",
        help="convert pfm logodds score to z-score",
        action="store_true",
        default=False,
    )
    p.add_argument(
        "--gc",
        dest="gcnorm",
        help="use GC frequency normalized z-score",
        action="store_true",
        default=False,
    )
    p.add_argument(
        "-N",
        "--nthreads",
        dest="ncpus",
        help="Number of threads (default %s)" % (params["ncpus"]),
        metavar="INT",
        type=int,
        default=int(params["ncpus"]),
    )
    p.add_argument(
        "-M",
        "--do_MOODS",
        dest="moods",
        help=argparse.SUPPRESS,
        # help="experimental: use MOODS for scanning",
        action="store_true",
        default=False,
    )
    p.add_argument(
        "-P",
        "--pvalue",
        dest="pvalue",
        help=argparse.SUPPRESS,
        # help="experimental: MOODS p-value cutoff",
        metavar="",
        type=float,
        default=None,
    )

    p.set_defaults(func=commands.pfmscan)

    p = subparsers.add_parser("maelstrom")
    p.add_argument(
        "inputfile", help="file with regions and clusters", metavar="INPUTFILE"
    )
    p.add_argument("genome", help="genome", metavar="GENOME")
    p.add_argument("outdir", help="output directory", metavar="DIR")
    p.add_argument(
        "-p",
        "--pfmfile",
        dest="pfmfile",
        help="PFM file with motifs " "(default: {0})".format(default_pfm),
        default=default_pfm_file,
        metavar="pfmfile",
    )
    p.add_argument(
        "--no-filter",
        dest="filter_redundant",
        help="Don't remove redundant motifs.",
        default=True,
        action="store_false",
    )
    p.add_argument(
        "-F",
        "--filter_cutoff",
        dest="filter_cutoff",
        help="Cutoff to select non-redundant motifs. Default is 0.8, increase this value to get fewer motifs.",
        default=0.8,
        type=float,
        metavar="FLOAT",
    )
    p.add_argument(
        "--nocenter",
        dest="center",
        help="Don't mean-center the rows by default",
        default=True,
        action="store_false",
    )
    p.add_argument(
        "-m",
        "--methods",
        dest="methods",
        help="Run with specific methods",
        default=None,
        metavar="NAMES",
    )
    p.add_argument(
        "-a",
        "--aggregation",
        dest="aggregation",
        help=(
            'How to combine motifs from individual methods. Default is "int_stouffer", '
            "for inverse normal transform of ranks, followed by Stouffer's method to combine "
            'z-scores. Alternatively, specify "stuart" for log-transformed rank aggregation '
            "p-values."
        ),
        default="int_stouffer",
        metavar="method",
    )
    p.add_argument(
        "-N",
        "--nthreads",
        dest="ncpus",
        help="Number of threads (default %s)" % (params["ncpus"]),
        metavar="INT",
        type=int,
        default=int(params["ncpus"]),
    )
    p.add_argument(
        "--rawscore",
        dest="zscore",
        help="Don't z-score normalize motif scores",
        action="store_false",
        default=True,
    )
    p.add_argument(
        "--nogc",
        dest="gc",
        help="Don't use GC%% bins",
        action="store_false",
        default=True,
    )

    p.set_defaults(func=commands.maelstrom)

    # closest_match.py
    p = subparsers.add_parser("match")
    p.add_argument("pfmfile", help="File with pfms", metavar="pfmfile")
    p.add_argument(
        "-d",
        dest="dbpfmfile",
        help="File with pfms to match against " "(default: {0})".format(default_pfm),
        default=default_pfm_file,
        metavar="DBFILE",
    )
    p.add_argument(
        "-n",
        dest="nmatches",
        help="Number of matches to return (default 1)",
        default=1,
        metavar="INT",
        type=int,
    )
    p.add_argument(
        "-o",
        dest="img",
        help="Output file with graphical report (png, svg, ps, pdf)",
        metavar="FILE",
    )
    p.set_defaults(func=commands.match)

    # pwm2logo.py
    p = subparsers.add_parser("logo")
    p.add_argument(
        "-p", "--pfmfile", help="PFM file with motifs", metavar="pfmfile", default=None
    )
    p.add_argument(
        "-i",
        "--ids",
        dest="ids",
        help="Comma-separated list of motif ids (default is all ids)",
        metavar="IDS",
    )
    p.add_argument(
        "-k",
        "--kind",
        dest="kind",
        help="Type of motif (information, frequency, energy or ensembl)",
        metavar="TYPE",
        default="information",
    )
    p.add_argument(
        "--notitle",
        dest="title",
        help="Don't include motif ID as title",
        default=True,
        action="store_false",
    )
    p.set_defaults(func=commands.logo)

    # motif_cluster.py
    p = subparsers.add_parser("cluster")
    p.add_argument("inputfile", help="Inputfile (PFM format)", metavar="INPUTFILE")
    p.add_argument("outdir", help="Name of output directory", metavar="OUTDIR")
    p.add_argument(
        "-s",
        dest="single",
        help="Don't compare reverse complements of motifs",
        default=False,
        action="store_true",
    )
    p.add_argument(
        "-t", dest="threshold", help="Cluster threshold", default=0.95, type=float
    )
    p.add_argument(
        "-N",
        "--nthreads",
        dest="ncpus",
        help="Number of threads (default %s)" % (params["ncpus"]),
        metavar="INT",
        type=int,
        default=int(params["ncpus"]),
    )

    p.set_defaults(func=commands.cluster)

    # generate_background_sequences.py
    p = subparsers.add_parser("background")
    p.add_argument("outputfile", help="outputfile", metavar="FILE")
    p.add_argument(
        "bg_type",
        help="type of background sequences to generate (%s)" % ",".join(BG_TYPES),
        metavar="TYPE",
    )
    p.add_argument(
        "-i", dest="inputfile", help="input sequences (BED or FASTA)", metavar="FILE"
    )
    p.add_argument(
        "-f",
        dest="outformat",
        help="output format (BED or FASTA",
        metavar="TYPE",
        default="fasta",
    )
    p.add_argument(
        "-s", dest="size", help="size of random sequences", metavar="INT", type=int
    )
    p.add_argument(
        "-n",
        dest="number",
        help="number of sequence to generate",
        metavar="NUMBER",
        default=10,
        type=int,
    )
    p.add_argument(
        "-g",
        dest="genome",
        help="genome version (not for type 'random')",
        metavar="GENOME",
    )
    p.add_argument(
        "-m",
        dest="markov_order",
        help="order of the Markov model (only for type 'random', default 1)",
        metavar="N",
        default=1,
        type=int,
    )
    p.set_defaults(func=commands.background)

    # get_fpr_based_pfmscan_threshold.py
    p = subparsers.add_parser("threshold")
    p.add_argument("pfmfile", help="File with pfms", metavar="pfmfile")
    p.add_argument(
        "inputfile", help="FASTA file with background sequences", metavar="FAFILE"
    )
    p.add_argument("fpr", help="Desired fpr", type=float, metavar="FPR")
    p.set_defaults(func=commands.threshold)

    # motif_localization_plots.py
    p = subparsers.add_parser("location")
    p.add_argument("pfmfile", help="File with pfms", metavar="pfmfile")
    p.add_argument("fastafile", help="Fasta formatted file", metavar="FAFILE")
    p.add_argument(
        "-s",
        dest="size",
        help="Set size to W (default: determined from fastafile)",
        metavar="INT",
        type=int,
    )
    p.add_argument(
        "-i",
        dest="ids",
        help="Comma-separated list of motif ids to plot (default is all ids)",
        metavar="IDS",
    )
    p.add_argument(
        "-c",
        dest="cutoff",
        help="Cutoff for motif scanning (default 0.95)",
        type=float,
        default=0.95,
    )
    p.set_defaults(func=commands.location)

    p = subparsers.add_parser("diff")
    p.add_argument(
        "inputfiles",
        help=(
            "FASTA-formatted inputfiles OR a BED file with an identifier in the 4th "
            "column, for instance a cluster number."
        ),
        metavar="FAFILES",
    )
    p.add_argument("bgfile", help="FASTA-formatted background file", metavar="BGFAFILE")
    p.add_argument("outputfile", help="outputfile (image)", metavar="PNGFILE")
    p.add_argument(
        "-p",
        "--pfmfile",
        dest="pfmfile",
        help="PFM file with motifs " "(default: {0})".format(default_pfm),
        default=default_pfm_file,
        metavar="pfmfile",
    )
    p.add_argument(
        "-c",
        "--cutoff",
        dest="cutoff",
        help="motif score cutoff or file with cutoffs (default 0.9)",
        metavar="",
        default=0.9,
    )
    p.add_argument(
        "-e",
        "--enrichment",
        dest="minenr",
        help="minimum enrichment in at least one of the datasets compared to background",
        metavar="MINENR",
        type=float,
        default=3,
    )
    p.add_argument(
        "-f",
        "--frequency",
        dest="minfreq",
        help="minimum frequency in at least one of the datasets",
        metavar="MINFREQ",
        type=float,
        default=0.01,
    )
    p.add_argument(
        "-g",
        "--genome",
        dest="genome",
        help=(
            "Genome; only necessary in combination with a BED file with clusters "
            "as inputfile."
        ),
        metavar="GENOME",
    )
    p.set_defaults(func=commands.diff)

    p.set_defaults(func=commands.logo)

    p = subparsers.add_parser("prediction")
    p.add_argument("tool", help="Specific motif prediction tool to run", metavar="NAME")
    p.add_argument("infile", help="Input FASTA file", metavar="FILE")
    p.add_argument("outfile", help="Output PFM file", metavar="FILE")
    p.add_argument(
        "-p",
        dest="paramfile",
        help="YAML file with paramaters",
        default=None,
        metavar="FILE",
    )
    p.set_defaults(func=commands.prediction)

    if len(sys_args) == 0:
        parser.print_help()
    elif sys_args[0] == "roc":
        print(
            "This command is deprecated. "
            "Use the following command for the same functionality:"
        )
        print()
        print("$ gimme motifs <inputfile> <outdir> --known")
        sys.exit(1)
    else:
        if len(sys_args) == 1:
            print(
                "\033[93mtype `gimme {} -h` for more details\033[0m\n".format(
                    sys_args[-1]
                )
            )
        args = parser.parse_args(sys_args)

        if hasattr(args, "genome"):
            if args.genome is not None:
                if not check_genome(args.genome):
                    print(
                        "Genome not found. Have you installed your genome with genomepy?"
                    )
                    print("See https://github.com/simonvh/genomepy for details.")
                    print("Alternatively, you can specify a FASTA file.")
                    exit(1)

        args.func(args)
