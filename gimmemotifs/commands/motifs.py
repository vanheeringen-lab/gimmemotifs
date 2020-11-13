#!/usr/bin/python -W ignore
# Copyright (c) 2009-2019 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
"""Command line function 'roc'."""
from __future__ import print_function
import os
import re
import sys
import shutil
import logging
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd

from gimmemotifs.background import create_background_file
from gimmemotifs.comparison import MotifComparer, select_nonredundant_motifs

from gimmemotifs.denovo import gimme_motifs
from gimmemotifs.motif import read_motifs
from gimmemotifs.stats import calc_stats_iterator
from gimmemotifs.report import roc_html_report
from gimmemotifs.scanner import scan_to_file
from gimmemotifs.utils import (
    determine_file_type,
    narrowpeak_to_bed,
    write_equalsize_bedfile,
)


logger = logging.getLogger("gimme.motifs")


def motifs(args):
    """ Calculate ROC_AUC and other metrics and optionally plot ROC curve."""
    if args.outdir is None:
        raise ValueError("an output directory is required!")
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    scan_dir = os.path.join(args.outdir, "motif_scan_results")
    if not os.path.exists(scan_dir):
        os.makedirs(scan_dir)

    file_type = determine_file_type(args.sample)
    outfile = os.path.join(args.outdir, f"input.w{args.size}.bed")
    sample = args.sample
    if file_type == "narrowpeak":
        narrowpeak_to_bed(args.sample, outfile, size=args.size)
        sample = outfile
    elif args.size and args.size > 0:
        if file_type == "fasta":
            logger.warn("size parameter will be ignored for FASTA input")
        elif file_type == "bed":
            write_equalsize_bedfile(args.sample, args.size, outfile)
            sample = outfile

    genome = args.genome
    if genome is None:
        args.zscore = False
        args.gc = False

    bgfile = None
    bg = args.background
    if bg is None:
        if genome is None:
            bg = "random"
        else:
            bg = "gc"

    if os.path.isfile(bg):
        bgfile = bg
        bg = "custom"
    else:
        # create background if not provided
        bgfile = os.path.join(args.outdir, "generated_background.{}.fa".format(bg))
        size = args.size
        if size <= 0:
            size = None
        if bg == "gc":
            logger.info("creating background (matched GC%)")
        else:
            logger.info("creating background (random)")

        create_background_file(
            bgfile,
            bg,
            fmt="fasta",
            genome=genome,
            inputfile=sample,
            size=size,
            number=10000,
        )

    pfmfile = args.pfmfile

    motifs = []
    if args.known:
        motifs = read_motifs(pfmfile, fmt="pfm")

    if args.denovo:
        gimme_motifs(
            sample,
            args.outdir,
            params={
                "tools": args.tools,
                "analysis": args.analysis,
                "background": bg,
                "custom_background": bgfile,
                "genome": args.genome,
                "size": args.size,
            },
        )
        denovo = read_motifs(os.path.join(args.outdir, "gimme.denovo.pfm"))
        mc = MotifComparer()
        result = mc.get_closest_match(denovo, dbmotifs=pfmfile, metric="seqcor")
        match_motifs = read_motifs(pfmfile, as_dict=True)
        new_map_file = os.path.join(args.outdir, "combined.motif2factors.txt")
        base = os.path.splitext(pfmfile)[0]
        map_file = base + ".motif2factors.txt"
        if os.path.exists(map_file):
            shutil.copyfile(map_file, new_map_file)

        motifs += denovo
        pfmfile = os.path.join(args.outdir, "combined.pfm")
        with open(pfmfile, "w") as f:
            for m in motifs:
                print(m.to_pwm(), file=f)

        with open(new_map_file, "a") as f:
            for m in denovo:
                print(
                    "{}\t{}\t{}\t{}".format(m.id, "de novo", "GimmeMotifs", "Y"), file=f
                )
                if result[m.id][0] in match_motifs:
                    for factor in match_motifs[result[m.id][0]].factors["direct"]:
                        print(
                            "{}\t{}\t{}\t{}".format(
                                m.id, factor, "inferred (GimmeMotifs)", "N"
                            ),
                            file=f,
                        )
    else:
        logger.info("skipping de novo")

    stats = [
        "phyper_at_fpr",
        "roc_auc",
        "pr_auc",
        "enr_at_fpr",
        "recall_at_fdr",
        "roc_values",
        "matches_at_fpr",
    ]

    f_out = sys.stdout
    if args.outdir:
        f_out = open(args.outdir + "/gimme.roc.report.txt", "w")

    # Print the metrics
    f_out.write(
        "Motif\t# matches\t% matches input\t# matches background\t%matches background\tP-value\tlog10 P-value\tROC AUC\tPR AUC\tEnr. at 1% FPR\tRecall at 10% FDR\n"
    )

    logger.info("creating motif scan tables")
    # ftype = determine_file_type(args.sample)
    # sample = args.sample
    # delete_sample = False
    # if ftype == "narrowpeak":
    #    f = NamedTemporaryFile(delete=False)
    #    logger.debug("Using {} as temporary BED file".format(f.name))
    #    narrowpeak_to_bed(args.sample, f.name, size=args.size)
    #    sample = f.name
    #    delete_sample = True

    # Create a table with the best score per motif for all motifs.
    # This has three reasons:
    # * Can be used to calculate statistics;
    # * Can be used to select a set of non-redundant motifs;
    # * These files are included in the output and can be used for further analyis.
    score_table = os.path.join(scan_dir, "input.motif.score.txt")
    bg_score_table = os.path.join(scan_dir, "background.motif.score.txt")
    for infile, outfile in [(sample, score_table), (bgfile, bg_score_table)]:
        scan_to_file(
            infile,
            pfmfile,
            filepath_or_buffer=outfile,
            score_table=True,
            genome=args.genome,
            zscore=True,
            gcnorm=True,
        )

    n_input = pd.read_csv(score_table, comment="#", sep="\t").shape[0]
    n_background = pd.read_csv(bg_score_table, comment="#", sep="\t").shape[0]

    logger.info("calculating stats")
    for motif_stats in calc_stats_iterator(
        motifs=pfmfile,
        fg_table=score_table,
        bg_table=bg_score_table,
        stats=stats,
        ncpus=args.ncpus,
    ):
        for motif in motifs:
            if str(motif) in motif_stats:
                log_pvalue = np.inf
                if motif_stats[str(motif)]["phyper_at_fpr"] > 0:
                    log_pvalue = -np.log10(motif_stats[str(motif)]["phyper_at_fpr"])
                f_out.write(
                    "{}\t{:d}\t{:.3f}\t{:d}\t{:.3f}\t{:.2e}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.2f}\t{:0.4f}\n".format(
                        motif.id,
                        motif_stats[str(motif)]["matches_at_fpr"][0],
                        motif_stats[str(motif)]["matches_at_fpr"][0] / n_input * 100,
                        motif_stats[str(motif)]["matches_at_fpr"][1],
                        motif_stats[str(motif)]["matches_at_fpr"][1]
                        / n_background
                        * 100,
                        motif_stats[str(motif)]["phyper_at_fpr"],
                        log_pvalue,
                        motif_stats[str(motif)]["roc_auc"],
                        motif_stats[str(motif)]["pr_auc"],
                        motif_stats[str(motif)]["enr_at_fpr"],
                        motif_stats[str(motif)]["recall_at_fdr"],
                    )
                )
    f_out.close()

    # Select a set of "non-redundant" motifs.
    # Using Recursive Feature Elimination, a set of motifs is selected that
    # best explains the peaks in comparison to the background sequences.
    nr_motifs = select_nonredundant_motifs(
        args.outdir + "/gimme.roc.report.txt",
        pfmfile,
        score_table,
        bg_score_table,
        tolerance=0.001,
    )

    # Provide BED files with motif scan results for the non-redundant motifs
    # At the moment this is not ideal, as scanning is now performed twice
    # for this set of non-redundant motifs.
    motif_dict = dict([(m.id, m) for m in motifs])
    for motif in nr_motifs:
        with NamedTemporaryFile(mode="w") as f:
            print(motif_dict[motif].to_pwm(), file=f)
            f.flush()
            safe_name = re.sub(r"[^a-zA-Z0-9\-]+", "_", motif)
            scan_to_file(
                sample,
                f.name,
                filepath_or_buffer=os.path.join(scan_dir, f"{safe_name}.matches.bed"),
                bed=True,
                fpr=0.01,
                genome=args.genome,
                zscore=True,
                gcnorm=True,
            )

    if args.report:
        logger.info("creating statistics report")
        if args.outdir:
            roc_html_report(
                args.outdir,
                args.outdir + "/gimme.roc.report.txt",
                pfmfile,
                threshold=0.01,
                outname="gimme.motifs.redundant.html",
                link_matches=False,
            )
            roc_html_report(
                args.outdir,
                args.outdir + "/gimme.roc.report.txt",
                pfmfile,
                threshold=0.01,
                use_motifs=nr_motifs,
                link_matches=True,
            )
            logger.info(
                f"gimme motifs final report: {os.path.join(args.outdir, 'gimme.motifs.html')}"
            )
