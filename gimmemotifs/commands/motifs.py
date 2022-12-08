#!/usr/bin/python -W ignore
# Copyright (c) 2009-2019 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
"""Command line function 'motifs'"""
import logging
import os
import re
import shutil
import sys
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
from tqdm.auto import tqdm

from gimmemotifs.background import create_background_file
from gimmemotifs.comparison import MotifComparer, select_nonredundant_motifs
from gimmemotifs.denovo import gimme_motifs
from gimmemotifs.motif import read_motifs
from gimmemotifs.report import roc_html_report
from gimmemotifs.scanner import scan_to_file
from gimmemotifs.stats import calc_stats_iterator
from gimmemotifs.utils import (
    determine_file_type,
    narrowpeak_to_bed,
    write_equalsize_bedfile,
)

logger = logging.getLogger("gimme.motifs")


def motifs(args):
    """Calculate ROC_AUC and other metrics and optionally plot ROC curve."""
    if args.outdir is None:
        raise ValueError("an output directory is required!")
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    scan_dir = os.path.join(args.outdir, "motif_scan_results")
    if not os.path.exists(scan_dir):
        os.makedirs(scan_dir)

    sample = args.sample
    size = int(args.size)
    base = os.path.basename(os.path.splitext(sample)[0])
    outfile = os.path.join(args.outdir, f"{base}_input_w{size}.bed")
    file_type = determine_file_type(sample)
    if file_type == "narrowpeak":
        if not os.path.isfile(outfile):
            narrowpeak_to_bed(sample, outfile, size=size)
        sample = outfile
    elif file_type == "fasta":
        logger.info(
            "Will use the sequences from the FASTA input. If size is specified, this will be ignored."
        )
    elif file_type == "bed":
        if not os.path.isfile(outfile):
            write_equalsize_bedfile(sample, size, outfile)
        sample = outfile

    genome = args.genome
    if genome is None:
        logger.warning("Genome is not specified!")
        logger.warning(
            "This means the z-score transformation and GC% normalization of the motif scores cannot be performed."
        )
        logger.warning("Will continue, but performance may be impacted.")
        args.zscore = False
        args.gc = False

    bgfile = args.background
    if bgfile and os.path.isfile(bgfile):
        bg = "custom"
    else:
        bg = "gc" if genome else "random"
        bgfile = os.path.join(args.outdir, f"{base}_background_{bg}.fa")
        if not os.path.isfile(bgfile):
            # create background if not provided
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
                number=10_000,
            )

    pfmfile = args.pfmfile

    motifs = []
    if args.known:
        motifs = read_motifs(pfmfile, fmt="pfm")

    if args.denovo:
        denovo = gimme_motifs(
            sample,
            args.outdir,
            params={
                "tools": args.tools,
                "analysis": args.analysis,
                "background": bg,
                "custom_background": bgfile,
                "genome": genome,
                "size": args.size,
                "fraction": args.fraction,
                "use_strand": not args.single,
            },
        )
        if len(denovo) == 0:
            logger.warning("No de novo motifs found!")
        else:
            logger.info(f"{len(denovo)} de novo motifs found!")

            mc = MotifComparer()
            result = mc.get_closest_match(denovo, dbmotifs=pfmfile, metric="seqcor")
            match_motifs = read_motifs(pfmfile, as_dict=True)

            new_map_file = os.path.join(args.outdir, "combined.motif2factors.txt")
            map_file = os.path.splitext(pfmfile)[0] + ".motif2factors.txt"
            if os.path.exists(map_file):
                shutil.copyfile(map_file, new_map_file)

            motifs += denovo

            # save known + de novo motifs to file
            pfmfile = os.path.join(args.outdir, "combined.pfm")
            with open(pfmfile, "w") as f:
                for m in motifs:
                    print(m.to_ppm(), file=f)

            with open(new_map_file, "a") as f:
                for m in denovo:
                    print(
                        "{}\t{}\t{}\t{}".format(m.id, "de novo", "GimmeMotifs", "Y"),
                        file=f,
                    )
                    if result[m.id][0] in match_motifs:
                        for factor in match_motifs[result[m.id][0]].factors["direct"]:
                            print(
                                "{}\t{}\t{}\t{}".format(
                                    m.id, factor, "inferred (GimmeMotifs)", "N"
                                ),
                                file=f,
                            )

    if len(motifs) == 0:
        logger.error("No motifs to report!")
        sys.exit()

    logger.info("creating motif scan tables")
    # Create a table with the best score per motif for all motifs.
    # This has three reasons:
    # * Can be used to calculate statistics;
    # * Can be used to select a set of non-redundant motifs;
    # * These files are included in the output and can be used for further analyis.
    score_table = os.path.join(scan_dir, "motif_score_input.txt")
    bg_score_table = os.path.join(scan_dir, "motif_score_background.txt")
    for infile, outfile in [(sample, score_table), (bgfile, bg_score_table)]:
        scan_to_file(
            infile,
            pfmfile,
            filepath_or_buffer=outfile,
            score_table=True,
            genome=genome,
            zscore=args.zscore,
            gc=args.gc,
            progress=None,
        )

    n_input = pd.read_csv(score_table, comment="#", sep="\t").shape[0]
    n_background = pd.read_csv(bg_score_table, comment="#", sep="\t").shape[0]

    logger.info("calculating stats")
    stats = [
        "phyper_at_fpr",
        "roc_auc",
        "pr_auc",
        "enr_at_fpr",
        "recall_at_fdr",
        "roc_values",
        "matches_at_fpr",
    ]
    it = calc_stats_iterator(
        motifs=pfmfile,
        fg_table=score_table,
        bg_table=bg_score_table,
        stats=stats,
        ncpus=args.ncpus,
    )

    # Print the metrics
    f_out = sys.stdout
    close = False
    if args.outdir:
        f_out = open(os.path.join(args.outdir, "gimme.roc.report.txt"), "w")
        close = True
    f_out.write(
        "Motif\t# matches\t% matches input\t# matches background\t%matches background\tP-value\tlog10 P-value\tROC AUC\tPR AUC\tEnr. at 1% FPR\tRecall at 10% FDR\n"
    )
    for motif_stats in it:
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
    if close:
        f_out.close()

    # Select a set of "non-redundant" motifs.
    # Using Recursive Feature Elimination, a set of motifs is selected that
    # best explains the peaks in comparison to the background sequences.
    nr_motifs = select_nonredundant_motifs(
        os.path.join(args.outdir, "gimme.roc.report.txt"),
        pfmfile,
        score_table,
        bg_score_table,
        tolerance=0.001,
    )

    # Provide BED files with motif scan results for the non-redundant motifs
    # At the moment this is not ideal, as scanning is now performed twice
    # for this set of non-redundant motifs.
    motif_dict = dict([(m.id, m) for m in motifs])
    for motif in tqdm(
        nr_motifs, unit=" motifs", desc="Creating BED files for non-redundant motifs"
    ):
        with NamedTemporaryFile(mode="w") as f:
            print(motif_dict[motif].to_ppm(), file=f)
            f.flush()
            safe_name = re.sub(r"[^a-zA-Z0-9\-]+", "_", motif)
            scan_to_file(
                sample,
                f.name,
                filepath_or_buffer=os.path.join(scan_dir, f"{safe_name}.matches.bed"),
                bed=True,
                fpr=0.01,
                genome=genome,
                zscore=args.zscore,
                gc=args.gc,
                bgfile=bgfile,
                progress=False,
            )

    if args.report:
        logger.info("creating statistics report")
        roc_html_report(
            args.outdir,
            os.path.join(args.outdir, "gimme.roc.report.txt"),
            pfmfile,
            threshold=0.01,
            outname="gimme.motifs.redundant.html",
            link_matches=False,
        )
        roc_html_report(
            args.outdir,
            os.path.join(args.outdir, "gimme.roc.report.txt"),
            pfmfile,
            threshold=0.01,
            use_motifs=nr_motifs,
            link_matches=True,
        )
        logger.info(
            f"gimme motifs final report: {os.path.join(args.outdir, 'gimme.motifs.html')}"
        )
