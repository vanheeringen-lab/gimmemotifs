# Copyright (c) 2009-2017 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
"""De novo motif prediction."""
import datetime
import os
import sys
import logging
import logging.handlers
import shutil
import numpy as np

from gimmemotifs.config import MotifConfig, BG_RANK, parse_denovo_params
from gimmemotifs import mytmpdir
from gimmemotifs.genome_index import track2fasta
from gimmemotifs.validation import check_denovo_input
from gimmemotifs.utils import ( divide_file, divide_fa_file, 
                                    write_equalwidth_bedfile )
from gimmemotifs.fasta import Fasta
from gimmemotifs.background import ( MarkovFasta, MatchedGcFasta,
                                    PromoterFasta, RandomGenomicFasta )
from gimmemotifs.stats import calc_stats, rank_motifs, write_stats
from gimmemotifs.report import create_denovo_motif_report
from gimmemotifs.motif import read_motifs

logger = logging.getLogger("gimme.denovo")

def prepare_denovo_input_bed(inputfile, params, outdir):
    """Prepare a BED file for de novo motif prediction.

    All regions to same size; split in test and validation set;
    converted to FASTA.

    Parameters
    ----------
    inputfile : str
        BED file with input regions.

    params : dict
        Dictionary with parameters.

    outdir : str
        Output directory to save files.
    """
    logger.info("preparing input (BED)")
    
    # Create BED file with regions of equal size
    width = int(params["width"])
    bedfile = os.path.join(outdir, "input.bed")
    write_equalwidth_bedfile(inputfile, width, bedfile)
    
    abs_max = int(params["abs_max"])
    fraction = float(params["fraction"])
    pred_bedfile = os.path.join(outdir, "prediction.bed")
    val_bedfile = os.path.join(outdir, "validation.bed")
    # Split input into prediction and validation set
    logger.debug(
                "Splitting %s into prediction set (%s) and validation set (%s)",
                bedfile, pred_bedfile, val_bedfile)
    divide_file(bedfile, pred_bedfile, val_bedfile, fraction, abs_max)

    config = MotifConfig()
    index_dir = os.path.join(config.get_index_dir(), params["genome"])
    
    for infile in [pred_bedfile, val_bedfile]:
        track2fasta(
            index_dir, 
            infile, 
            infile.replace(".bed", ".fa"), 
            )

    # Create file for location plots
    lwidth = int(params["lwidth"])
    extend = (lwidth - width) // 2
    
    track2fasta(
            index_dir, 
            val_bedfile, 
            os.path.join(outdir, "localization.fa"), 
            extend_up=extend, 
            extend_down=extend, 
            use_strand=params["use_strand"], 
            ignore_missing=True
            )

def prepare_denovo_input_fa(inputfile, params, outdir):
    """Create all the FASTA files for de novo motif prediction and validation.
    
    Parameters
    ----------
    """
    fraction = float(params["fraction"])
    abs_max = int(params["abs_max"])

    logger.info("preparing input (FASTA)")

    pred_fa = os.path.join(outdir, "prediction.fa")
    val_fa = os.path.join(outdir, "validation.fa")
    loc_fa = os.path.join(outdir, "localization.fa")

    # Split inputfile in prediction and validation set
    logger.debug(
        "Splitting %s into prediction set (%s) and validation set (%s)",
        inputfile, pred_fa, val_fa)

    divide_fa_file(inputfile, pred_fa, val_fa, fraction, abs_max)

    # File for location plots
    shutil.copy(val_fa, loc_fa)
    seqs = Fasta(loc_fa).seqs
    lwidth = len(seqs[0])
    all_same_width = not(False in [len(seq) == lwidth for seq in seqs])
    if not all_same_width:
        logger.warn(
            "PLEASE NOTE: FASTA file contains sequences of different lengths. "
            "Positional preference plots might be incorrect!")

def create_background(bg_type, fafile, outfile, genome="hg18", width=200, nr_times=10, user_background=None):
    """Create background of a specific type.

    Parameters
    ----------
    bg_type : str
        Name of background type.

    fafile : str
        Name of input FASTA file.

    outfile : str
        Name of output FASTA file.

    genome : str, optional
        Genome name.

    width : int, optional
        Size of regions.

    nr_times : int, optional
        Generate this times as many background sequences as compared to 
        input file.
    
    Returns
    -------
    nr_seqs  : int
        Number of sequences created.
    """
    width = int(width)
    config = MotifConfig()
    fg = Fasta(fafile)
    if bg_type in ["promoter", "genomic"]:
        index_dir = os.path.join(config.get_index_dir(), genome)

    if bg_type == "random":
        f = MarkovFasta(fg, k=1, n=nr_times * len(fg))
        logger.debug("Random background: %s", outfile)
    elif bg_type == "genomic":
        logger.debug("Creating genomic background")
        f = RandomGenomicFasta(index_dir, width, nr_times * len(fg))
    elif bg_type == "gc":
        logger.debug("Creating GC matched background")
        f = MatchedGcFasta(fafile, genome, nr_times * len(fg))
        logger.debug("GC matched background: %s", outfile)
    elif bg_type == "promoter":
        gene_file = os.path.join(config.get_gene_dir(), "%s.bed" % genome)

        logger.info(
                "Creating random promoter background (%s, using genes in %s)",
                genome, gene_file)
        f = PromoterFasta(gene_file, index_dir, width, nr_times * len(fg))
        logger.debug("Random promoter background: %s", outfile)
    elif bg_type == "user":
        bg_file = user_background
        if not bg_file:
            raise IOError(
                    "Background file not specified!")

        if not os.path.exists(bg_file):
            raise IOError(
                    "User-specified background file %s does not exist!",
                    bg_file)
        else:
            logger.info("Copying user-specified background file %s to %s.",
                    bg_file, outfile)
            f = Fasta(bg_file)
            l = np.median([len(seq) for seq in f.seqs])
            if l < (width * 0.95) or l > (width * 1.05):
                   logger.warn(
                    "The user-specified background file %s contains sequences with a "
                    "median length of %s, while GimmeMotifs predicts motifs in sequences "
                    "of length %s. This will influence the statistics! It is recommended "
                    "to use background sequences of the same length.", 
                    bg_file, l, width)
    
    f.writefasta(outfile)
    return len(f)

def create_backgrounds(outdir, background=None, genome="hg38", width=200, user_background=None):
    """Create different backgrounds for motif prediction and validation.

    Parameters
    ----------
    outdir : str
        Directory to save results.
    
    background : list, optional
        Background types to create, default is 'random'.

    genome : str, optional
        Genome name (for genomic and gc backgrounds).

    width : int, optional
        Size of background regions

    Returns
    -------
    bg_info : dict
        Keys: background name, values: file name.
    """
    if background is None:
        background = ["random"]
        nr_sequences = {}

    # Create background for motif prediction
    if "gc" in background:
        pred_bg = "gc"
    else:
        pred_bg = background[0]
    create_background(
                    pred_bg, 
                    os.path.join(outdir, "validation.fa"), 
                    os.path.join(outdir, "prediction.bg.fa"), 
                    genome=genome, 
                    width=width,
                    user_background=user_background)

    # Get background fasta files for statistics
    bg_info = {}
    nr_sequences = {}    
    for bg in background:
        fname = os.path.join(outdir, "bg.{}.fa".format(bg))
        nr_sequences[bg] = create_background(
                                        bg, 
                                        os.path.join(outdir, "validation.fa"), 
                                        fname, 
                                        genome=genome, 
                                        width=width,
                                        user_background=user_background)

        bg_info[bg] = fname
    return bg_info


def _is_significant(stats, metrics=None):
    """Filter significant motifs based on several statistics.
    
    Parameters
    ----------
    stats : dict
        Statistics disctionary object.

    metrics : sequence
        Metric with associated minimum values. The default is
        (("max_enrichment", 3), ("roc_auc", 0.55), ("enr_at_fpr", 0.55))
    
    Returns
    -------
    significant : bool
    """
    if metrics is None:
        metrics = (("max_enrichment", 3), ("roc_auc", 0.55), ("enr_at_fpr", 0.55))
    
    for stat_name, min_value in metrics:
        if stats.get(stat_name, 0) < min_value:
            return False

    return True
     
def filter_significant_motifs(fname, result, bg, metrics=None):
    """Filter significant motifs based on several statistics.

    Parameters
    ----------
    fname : str
        Filename of output file were significant motifs will be saved.

    result : PredictionResult instance
        Contains motifs and associated statistics.

    bg : str
        Name of background type to use.

    metrics : sequence
        Metric with associated minimum values. The default is
        (("max_enrichment", 3), ("roc_auc", 0.55), ("enr_at_f[r", 0.55))

    Returns
    -------
    motifs : list
        List of Motif instances.
    """
    sig_motifs = []
    with open(fname, "w") as f:
        for motif in result.motifs:
            stats = result.stats.get(
                    "%s_%s" % (motif.id, motif.to_consensus()), {}).get(bg, {}
                    ) 
            if _is_significant(stats, metrics):
                f.write("%s\n" % motif.to_pfm())
                sig_motifs.append(motif)
    
    logger.info("%s motifs are significant", len(sig_motifs))
    logger.debug("written to %s", fname)
    
    return sig_motifs

def best_motif_in_cluster(single_pwm, clus_pwm, clusters, fg_fa, background, stats=None, metrics=("roc_auc", "recall_at_fdr")):
    """Return the best motif per cluster for a clustering results.

    The motif can be either the average motif or one of the clustered motifs.

    Parameters
    ----------
    single_pwm : str
        Filename of motifs.

    clus_pwm : str
        Filename of motifs.

    clusters : 
        Motif clustering result.

    fg_fa : str
        Filename of FASTA file.

    background : dict
        Dictionary for background file names.

    stats : dict, optional
        If statistics are not supplied they will be computed.

    metrics : sequence, optional
        Metrics to use for motif evaluation. Default are "roc_auc" and 
        "recall_at_fdr".
    
    Returns
    -------
    motifs : list
        List of Motif instances.
    """
    # combine original and clustered motifs
    with open(single_pwm) as f1:
        with open(clus_pwm) as f2:
            motifs = read_motifs(f1) + read_motifs(f2)
    motifs = dict([(str(m), m) for m in motifs])

    # get the statistics for those motifs that were not yet checked
    clustered_motifs = []
    for clus,singles in clusters:
        for motif in set([clus] + singles):
            if str(motif) not in stats:
                clustered_motifs.append(motifs[str(motif)])
    
    new_stats = {}
    for bg, bg_fa in background.items():
        for m,s in calc_stats(clustered_motifs, fg_fa, bg_fa).items():
            if m not in new_stats:
                new_stats[m] = {}
            new_stats[m][bg] = s
    stats.update(new_stats)
    
    rank = rank_motifs(stats, metrics)

    # rank the motifs
    best_motifs = []
    for clus, singles in clusters:
        if len(singles) > 1:
            eval_motifs = singles
            if clus not in motifs:
                eval_motifs.append(clus)
            eval_motifs = [motifs[str(e)] for e in eval_motifs]
            best_motif = sorted(eval_motifs, key=lambda x: rank[str(x)])[-1]
            best_motifs.append(best_motif)
        else:
            best_motifs.append(clus)
        for bg in background:
            stats[str(best_motifs[-1])][bg]["num_cluster"] = len(singles)

    best_motifs = sorted(best_motifs, key=lambda x: rank[str(x)], reverse=True)
    return best_motifs
    
def rename_motifs(motifs, stats=None):
    """Rename motifs to GimmeMotifs_1..GimmeMotifs_N.
    
    If stats object is passed, stats will be copied."""
    final_motifs = []
    for i, motif in enumerate(motifs):
        old = str(motif)
        motif.id = "GimmeMotifs_{}".format(i + 1)
        final_motifs.append(motif)
        if stats:
            stats[str(motif)] = stats[old].copy()
    
    if stats:
        return final_motifs, stats
    else:
        return final_motifs

def gimme_motifs(inputfile, outdir, params=None, filter_significant=True, cluster=True, create_report=True):
    """De novo motif prediction based on an ensemble of different tools.

    Parameters
    ----------
    inputfile : str
        Filename of input. Can be either BED or FASTA.

    outdir : str
        Name of output directory.

    params : dict, optional
        Optional parameters.

    filter_significant : bool, optional
        Filter motifs for significance using the validation set.
    
    cluster : bool, optional
        Cluster similar predicted (and significant) motifs.

    create_report : bool, optional
        Create output reports (both .txt and .html).
 
    Returns
    -------
    motifs : list
        List of predicted motifs.     

    Examples
    --------

    >>> from gimmemotifs.denovo import gimme_motifs
    >>> gimme_motifs("input.fa", "motifs.out")
    """
    if outdir is None:
        outdir = "gimmemotifs_{}".format(datetime.date.today().strftime("%d_%m_%Y"))

    # Create output directories
    tmpdir = os.path.join(outdir, "intermediate")
    for d in [outdir, tmpdir]: 
        if not os.path.exists(d):
            os.mkdir(d) 

    # setup logfile
    logger = logging.getLogger("gimme")
    # Log to file
    logfile = os.path.join(outdir, "gimmemotifs.log")
    fh = logging.FileHandler(logfile, "w")
    fh.setLevel(logging.DEBUG)
    file_formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    fh.setFormatter(file_formatter)
    logger.addHandler(fh)
    logger = logging.getLogger("gimme.denovo.gimme_motifs")
    
    # Initialize parameters
    params = parse_denovo_params(params)
 
    # Check the input files
    input_type, background = check_denovo_input(inputfile, params)
   
    logger.info("starting full motif analysis")
    logger.debug("Using temporary directory %s", mytmpdir())
    
    
    # Create the necessary files for motif prediction and validation
    if input_type == "BED":
        prepare_denovo_input_bed(inputfile, params, tmpdir)
    elif input_type == "FASTA":
        prepare_denovo_input_fa(inputfile, params, tmpdir)
    else:
        logger.error("Unknown input type, shouldn't happen")
        sys.exit(1)

    # Create the background FASTA files
    background = create_backgrounds(tmpdir, background, params["genome"], params["width"], params.get("user_background", None))
    
    # Predict de novo motifs
    result = predict_motifs(
            os.path.join(tmpdir, "prediction.fa"),
            os.path.join(tmpdir, "prediction.bg.fa"),
            os.path.join(tmpdir, "all_motifs.pfm"),
            params=params,
            stats_fg=os.path.join(tmpdir, 'validation.fa'),
            stats_bg=background, 
            )

    if len(result.motifs) == 0:
        logger.info("finished")
        return []
    
    # Write statistics
    stats_file = os.path.join(tmpdir, "stats.{}.txt")
    write_stats(result.stats, stats_file)

    bg = sorted(background, key=lambda x: BG_RANK[x])[0]
    if filter_significant:
        motifs = filter_significant_motifs(
                os.path.join(tmpdir, "significant_motifs.pfm"),
                result, 
                bg)
        if len(motifs) == 0:
            logger.info("no significant motifs")
            return 

        pwmfile = os.path.join(tmpdir, "significant_motifs.pfm")
    else:
        logger.info("not filtering for significance")
        motifs = result.motifs
        pwmfile = os.path.join(tmpdir, "all_motifs.pfm")

    if cluster: 
        clusters = cluster_motifs_with_report(
                    pwmfile,
                    os.path.join(tmpdir, "clustered_motifs.pfm"),
                    outdir,
                    0.95,
                    title=inputfile)
        
        # Determine best motif in cluster
        best_motifs = best_motif_in_cluster(
                pwmfile,
                os.path.join(tmpdir, "clustered_motifs.pfm"),
                clusters, 
                os.path.join(tmpdir, 'validation.fa'), 
                background, 
                result.stats)
        
        final_motifs, stats = rename_motifs(best_motifs, result.stats)
    else:
        logger.info("not clustering")
        rank = rank_motifs(result.stats)
        sorted_motifs = sorted(motifs, key=lambda x: rank[str(x)], reverse=True)
        final_motifs, stats = rename_motifs(sorted_motifs, result.stats)

    with open(os.path.join(outdir, "motifs.pwm"), "w") as f:
        for m in final_motifs:
            f.write("{}\n".format(m.to_pwm()))
    
    if create_report:
        bg = dict([(b, os.path.join(tmpdir, "bg.{}.fa".format(b))) for b in background])

        create_denovo_motif_report(
                inputfile, 
                os.path.join(outdir, "motifs.pwm"), 
                os.path.join(tmpdir, "validation.fa"), 
                bg, 
                os.path.join(tmpdir, "localization.fa"), 
                outdir,
                params,
                stats,
                )
    
    with open(os.path.join(outdir, "params.txt"), "w") as f:
        for k,v in params.items():
            f.write("{}\t{}\n".format(k,v))
    
    if not(params.get("keep_intermediate")):
        logger.debug(
            "Deleting intermediate files. "
            "Please specifify the -k option if you want to keep these files.")
        shutil.rmtree(tmpdir)

    logger.info("finished")
    logger.info("output dir: %s", outdir) 
    if cluster:
        logger.info("report: %s", os.path.join(outdir, "motif_report.html"))

    return final_motifs

from gimmemotifs.cluster import cluster_motifs_with_report
from gimmemotifs.prediction import predict_motifs
