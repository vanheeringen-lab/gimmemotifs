#!/usr/bin/python -W ignore
# Copyright (c) 2009-2019 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
"""Command line function 'motifs'."""
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
from gimmemotifs.motif import gimme_motifs, read_motifs
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
    _motifs(
        sample=args.sample,
        analysis=args.analysis,
        tools=args.tools,
        outdir=args.outdir,
        size=args.size,
        fraction=args.fraction,
        single=args.single,
        zscore=args.zscore,
        gc=args.gc,
        genome=args.genome,
        background=args.background,
        pfmfile=args.pfmfile,
        known=args.known,
        denovo=args.denovo,
        report=args.report,
        ncpus=args.ncpus,
        seed=args.seed,
    )


def _motifs(
    sample,
    analysis,
    tools,
    outdir,
    size=200,
    fraction=0.2,
    single=False,
    zscore=True,
    gc=True,
    genome=None,
    background=None,
    pfmfile=None,
    known=True,
    denovo=True,
    report=True,
    ncpus=None,
    seed=None,
):
    """Calculate ROC_AUC and other metrics and optionally plot ROC curve."""
    if outdir is None:
        raise ValueError("an output directory is required!")
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    scan_dir = os.path.join(outdir, "motif_scan_results")
    if not os.path.exists(scan_dir):
        os.makedirs(scan_dir)

    random_state = None
    if seed:
        random_state = np.random.RandomState(int(seed))

    file_type = determine_file_type(sample)
    outfile = os.path.join(outdir, f"input.w{size}.bed")
    if file_type == "narrowpeak":
        narrowpeak_to_bed(sample, outfile, size=size)
        sample = outfile
    elif size and size > 0:
        if file_type == "fasta":
            logger.warning(
                "Will use the sequences from the FASTA input. "
                "If size is specified, this will be ignored."
            )
        elif file_type == "bed":
            write_equalsize_bedfile(sample, size, outfile)
            sample = outfile

    genome = genome
    if genome is None:
        logger.warning("Genome is not specified!")
        logger.warning(
            "This means the z-score transformation and GC% normalization "
            "of the motif scores cannot be performed."
        )
        logger.warning("Will continue, but performance may be impacted.")
        zscore = False
        gc = False

    if background is None:
        if genome is None:
            background = "random"
        else:
            background = "gc"

    if os.path.isfile(background):
        bgfile = background
        background = "custom"
    else:
        # create background if not provided
        bgfile = os.path.join(outdir, f"generated_background.{background}.fa")
        if size <= 0:
            size = None
        if background == "gc":
            logger.info("creating background (matched GC%)")
        else:
            logger.info("creating background (random)")

        create_background_file(
            bgfile,
            background,
            fmt="fasta",
            genome=genome,
            inputfile=sample,
            size=size,
            number=10000,
            random_state=random_state,
        )

    motifs = []
    if known:
        motifs = read_motifs(pfmfile, fmt="pfm")

    if denovo:
        denovo = gimme_motifs(
            sample,
            outdir,
            params={
                "tools": tools,
                "analysis": analysis,
                "background": background,
                "custom_background": bgfile,
                "genome": genome,
                "size": size,
                "fraction": fraction,
                "use_strand": not single,
            },
            random_state=random_state,
        )
    else:
        denovo = []

    if len(denovo) > 0:
        mc = MotifComparer()
        result = mc.get_closest_match(denovo, dbmotifs=pfmfile, metric="seqcor")
        match_motifs = read_motifs(pfmfile, as_dict=True)

        map_file = os.path.splitext(pfmfile)[0] + ".motif2factors.txt"
        new_map_file = os.path.join(outdir, "combined.motif2factors.txt")
        if os.path.exists(map_file):
            shutil.copyfile(map_file, new_map_file)

        motifs += denovo

        # save known + de novo motifs to file
        pfmfile = os.path.join(outdir, "combined.pfm")
        with open(pfmfile, "w") as f:
            for m in motifs:
                print(m.to_ppm(), file=f)

        with open(new_map_file, "a") as f:
            for m in denovo:
                print(f"{m.id}\tde novo\tGimmeMotifs\tY", file=f)
                found_motif = result[m.id][0]
                if found_motif not in match_motifs:
                    continue
                for factor in match_motifs[found_motif].factors["direct"]:
                    print(f"{m.id}\t{factor}\tinferred (GimmeMotifs)\tN", file=f)

    if len(motifs) == 0:
        logger.info("no motifs to report!")
        sys.exit()

    logger.info("creating motif scan tables")
    # Create a table with the best score per motif for all motifs.
    # This has three reasons:
    # * Can be used to calculate statistics;
    # * Can be used to select a set of non-redundant motifs;
    # * These files are included in the output and can be used for further analysis.
    score_table = os.path.join(scan_dir, "input.motif.score.txt")
    bg_score_table = os.path.join(scan_dir, "background.motif.score.txt")
    for infile, outfile in [(sample, score_table), (bgfile, bg_score_table)]:
        scan_to_file(
            infile,
            pfmfile,
            filepath_or_buffer=outfile,
            score_table=True,
            genome=genome,
            zscore=zscore,
            gcnorm=gc,
            random_state=random_state,
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
        ncpus=ncpus,
        random_state=random_state,
    )

    # Print the metrics
    roc_report = os.path.join(outdir, "gimme.roc.report.txt")
    with open(roc_report, "w") as f_out:
        header = [
            "Motif",
            "# matches",
            "% matches input",
            "# matches background",
            "% matches background",
            "P-value",
            "log10 P-value",
            "ROC AUC",
            "PR AUC",
            "Enr. at 1% FPR",
            "Recall at 10% FDR",
        ]
        f_out.write("\t".join(header) + "\n")

        for motif_stats in it:
            for motif in motifs:
                if str(motif) not in motif_stats:
                    continue
                stats = motif_stats[str(motif)]
                ph_fpr = stats["phyper_at_fpr"]
                log_pvalue = np.inf if ph_fpr <= 0 else -np.log10(ph_fpr)
                f_out.write(
                    "{}\t{:d}\t{:.3f}\t{:d}\t{:.3f}\t{:.2e}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.2f}\t{:0.4f}\n".format(
                        motif.id,
                        stats["matches_at_fpr"][0],
                        stats["matches_at_fpr"][0] / n_input * 100,
                        stats["matches_at_fpr"][1],
                        stats["matches_at_fpr"][1] / n_background * 100,
                        stats["phyper_at_fpr"],
                        log_pvalue,
                        stats["roc_auc"],
                        stats["pr_auc"],
                        stats["enr_at_fpr"],
                        stats["recall_at_fdr"],
                    )
                )

    # Select a set of "non-redundant" motifs.
    # Using Recursive Feature Elimination, a set of motifs is selected that
    # best explains the peaks in comparison to the background sequences.
    nr_motifs = select_nonredundant_motifs(
        roc_report,
        pfmfile,
        score_table,
        bg_score_table,
        tolerance=0.001,
        random_state=random_state,
    )

    # Provide BED files with motif scan results for the non-redundant motifs
    # At the moment this is not ideal, as scanning is now performed twice
    # for this set of non-redundant motifs.
    motif_dict = dict([(m.id, m) for m in motifs])
    logger.info("creating BED files with scan results")
    for motif in tqdm(nr_motifs):
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
                zscore=zscore,
                gcnorm=gc,
                bgfile=bgfile,
                random_state=random_state,
                progress=False,
            )

    if report:
        logger.info("creating statistics report")
        roc_html_report(
            outdir,
            roc_report,
            pfmfile,
            threshold=0.01,
            outname="gimme.motifs.redundant.html",
            link_matches=False,
        )
        roc_html_report(
            outdir,
            roc_report,
            pfmfile,
            threshold=0.01,
            use_motifs=nr_motifs,
            link_matches=True,
        )
        logger.info(
            f"gimme motifs final report: {os.path.join(outdir, 'gimme.motifs.html')}"
        )
