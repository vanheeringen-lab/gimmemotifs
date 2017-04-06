# Copyright (c) 2009-2017 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
from gimmemotifs.config import MotifConfig, BG_RANK, parse_denovo_params
from gimmemotifs import mytmpdir
from gimmemotifs.validation import check_denovo_input
from gimmemotifs.utils import divide_fa_file, motif_localization
from gimmemotifs.fasta import Fasta
from gimmemotifs.background import ( MarkovFasta, MatchedGcFasta,
                                    PromoterFasta, RandomGenomicFasta )
from gimmemotifs.stats import calc_stats
from gimmemotifs.report import create_denovo_motif_report
from gimmemotifs.motif import read_motifs
from tempfile import NamedTemporaryFile
from time import sleep

import os
import sys
import logging

logger = logging.getLogger()

def prepare_denovo_input_bed():
    self.prepare_input_bed(inputfile, params["genome"], params["width"], params["fraction"], params["abs_max"], params["use_strand"])


    # Create file for location plots
    index_dir = os.path.join(self.config.get_index_dir(), params["genome"])
    lwidth = int(params["lwidth"])
    width = int(params["width"])
    extend = (lwidth - width) / 2
    genome_index.track2fasta(index_dir, self.validation_bed, self.location_fa, extend_up=extend, extend_down=extend, use_strand=params["use_strand"], ignore_missing=True)

def write_stats(stats, fname):
    # Write stats output to file

    for bg in stats.values()[0].keys():
        f = open(fname.format(bg), "w")
        stat_keys = sorted(stats.values()[0].values()[0].keys())
        f.write("{}\t{}\n".format("Motif", "\t".join(stat_keys)))
    
        for motif in stats:
            m_stats = stats.get(str(motif), {}).get(bg)
            if m_stats:
                f.write("{}\t{}\n".format(
                    "_".join(motif.split("_")[:-1]), 
                    "\t".join([str(m_stats[k]) for k in stat_keys])
                    ))
            else:
                logger.warn("No stats for motif {0}, skipping this motif!".format(motif.id))
            #motifs.remove(motif)
        f.close()

    return
    self.motifs_with_stats = motifs

    f = open(self.ranks_file, "w")
    tools = dict((m.id.split("_")[2],1) for m in motifs).keys()
    f.write("Metric\tType\t%s\n" % ("\t".join(tools)))
    for stat in ["mncp", "roc_auc", "maxenr"]:
        best_motif = {}
        for motif in self.motifs_with_stats:
            d = stats.get("%s_%s" % (motif.id, motif.to_consensus()), {})
            val = d.get(stat, None)
            if val is None:
                continue
            name = motif.id.split("_")[2]
            if val > best_motif.setdefault(name, 0):
                best_motif[name] = val
        if len(best_motif.keys()) > 0:
            names = best_motif.keys()
            vals = [best_motif[name] for name in names]
            rank = rankdata(vals)
            ind = [names.index(x) for x in tools]

            f.write("%s\t%s\t%s\n" % (stat, "value", "\t".join([str(vals[i]) for i in ind])))
            f.write("%s\t%s\t%s\n" % (stat, "rank", "\t".join([str(rank[i]) for i in ind])))
    f.close()

def get_tempfile(name, tmpdir=mytmpdir):
    return "bla/{}".format(name)
    #
    #return os.path.join(tmpdir, "{}_{}".format(basename))


def prepare_denovo_input_fa(inputfile, params):
    """ Create all the bed- and fasta-files necessary for motif prediction and validation """
    width = int(params["width"])
    fraction = float(params["fraction"])
    abs_max = int(params["abs_max"])

    logger.info("preparing input (FASTA)")

    
    pred_fa = get_tempfile("prediction.fa")
    val_fa = get_tempfile("validation.fa")

    # Split inputfile in prediction and validation set
    logger.debug(
        "Splitting %s into prediction set (%s) and validation set (%s)",
        inputfile, pred_fa, val_fa)


    pred_num, val_num = divide_fa_file(inputfile, pred_fa, val_fa, fraction, abs_max)

#    # File for location plots
#    self.location_fa = self.validation_fa
#    fa = Fasta(self.location_fa)
#    seqs = fa.seqs
#    lwidth = len(seqs[0])
#    all_same_width = not(False in [len(seq) == lwidth for seq in seqs])
#    if not all_same_width:
#        self.logger.warn(
#                "PLEASE NOTE: FASTA file contains sequences of different lengths. "
#                "Positional preference plots will be incorrect!")

def create_background(bg_type, bedfile, fafile, outfile, genome="hg18", width=200, nr_times=10):
    fg = Fasta(fafile)
    if bg_type in ["promoter", "genomic"]:
        index_dir = os.path.join(self.config.get_index_dir(), genome)

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
        gene_file = os.path.join(self.config.get_gene_dir(), "%s.bed" % genome)

        logger.info(
                "Creating random promoter background (%s, using genes in %s)",
                genome, gene_file)
        f = PromoterFasta(gene_file, index_dir, width, nr_times * len(fg))
        logger.debug("Random promoter background: %s", outfile)
    elif bg_type == "user":
        bg_file = self.params["user_background"]
        if not os.path.exists(bg_file):
            raise IOError(
                    "User-specified background file %s does not exist!",
                    bg_file)
        else:
            logger.info("Copying user-specified background file %s to %s.",
                    bg_file, outfile)
            f = Fasta(bg_file)
            l = median([len(seq) for seq in fa.seqs])
            if l < width * 0.95 or l > width * 1.05:
                   logger.warn(
                    "The user-specified background file %s contains sequences with a "
                    "median length of %s, while GimmeMotifs predicts motifs in sequences "
                    "of length %s. This will influence the statistics! It is recommended "
                    "to use background sequences of the same length.", 
                    bg_file, l, width)

    f.writefasta(outfile)
    return len(f)

def create_backgrounds(background=None, genome="hg38", width=200):
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
                    get_tempfile("validation.bed"), 
                    get_tempfile("validation.fa"), 
                    get_tempfile("prediction.bg.fa"), 
                    genome=genome, 
                    width=width)

    # Get background fasta files for statistics
    bg_info = {}
    nr_sequences = {}    
    for bg in background:
        fname = get_tempfile("bg.{}.fa".format(bg))
        nr_sequences[bg] = create_background(
                                        bg, 
                                        get_tempfile("validation.bed"), 
                                        get_tempfile("validation.fa"), 
                                        fname, 
                                        genome=genome, 
                                        width=width)

        bg_info[bg] = fname
    return bg_info

def filter_significant_motifs(result, bg):
    # Determine significant motifs
    nsig = 0
    fname = get_tempfile("significant.pfm")
    f = open(fname, "w")
    sig_motifs = []
    for motif in result.motifs:
        stats = result.stats.get("%s_%s" % (motif.id, motif.to_consensus()), {}).get(bg, {}) 
        #if stats.get("maxenr", 0) >= 3 and stats.get("roc_auc", 0) >= 0.55 and stats.get('enr_fdr', 0) >= 2:
        if stats.get("roc_auc", 0) >= 0.55:
            f.write("%s\n" % motif.to_pfm())
            sig_motifs.append(motif)
            nsig += 1
    f.close()
    logger.info("%s motifs are significant", nsig)
    logger.debug("written to %s", fname)

    # ROC metrics of significant motifs
    #for bg in background:
    #    self._roc_metrics(self.significant_pfm, self.validation_fa, self.bg_file["fa"][bg], self.bg_file["roc"][bg])
    if nsig == 0:
        logger.info("no significant motifs found")
        return []
    
    return sig_motifs

def score_clustered_motifs(motifs, stats_fg, stats_bg):

    if isinstance(motifs, str):
        motifs = read_motifs(open(motifs))

    # Stars
    tmp = NamedTemporaryFile(dir=mytmpdir()).name
    p = PredictionResult(tmp, fg_file=stats_fg, background=stats_bg, do_counter=False)
    p.add_motifs(("clustering",  (motifs, "","")))
    p.wait_for_stats()
    
    while len(p.stats.keys()) < len(p.motifs):
        sleep(5)
    
    #print p.stats
    #print "num_cluster"
    #print num_cluster
    #for mid, num in num_cluster.items():
    #    for bg in p.stats[mid]:
    #    p.stats[mid]["numcluster"] = num

def best_motif_in_cluster(clusters, fg_fa, background, stats=None):
    pass

def gimme_motifs(inputfile, outdir, params=None, filter_significant=True, cluster=True, create_report=True):
    """ Full analysis: from bed-file to motifs (including clustering, ROC-curves, location plots and html report) """
    logger.info("starting full motif analysis")
    logger.debug("Using temporary directory {0}".format(mytmpdir()))

    # Initialize parameters
    params = parse_denovo_params(params)
 
    # Check the input files
    input_type, background = check_denovo_input(inputfile, params)
    
    if not os.path.exists(outdir):
        os.mkdir(outdir) 

    # Create the necessary files for motif prediction and validation
    if input_type == "BED":
        prepare_denovo_input_bed(inputfile, params)
    elif input_type == "FASTA":
        prepare_denovo_input_fa(inputfile, params)
    else:
        logger.error("Unknown input type, shouldn't happen")
        sys.exit(1)

    # Create the background FASTA files
    background = create_backgrounds(background, params["genome"], params["width"])
    
    # Predict de novo motifs
    result = predict_motifs(
            get_tempfile("prediction.fa"),
            get_tempfile("prediction.bg.fa"),
            get_tempfile("all_motifs.pfm"),
            params=params,
            stats_fg=get_tempfile('validation.fa'),
            stats_bg=background, 
            )

    # Write statistics
    stats_file = get_tempfile("stats.{}.txt")
    write_stats(result.stats, stats_file)

    bg = sorted(background, lambda x,y: cmp(BG_RANK[x], BG_RANK[y]))[0]
    if filter_significant:
        motifs = filter_significant_motifs(result, bg)
    else:
        motifs = result.motifs

    if cluster: 
        clusters = cluster_motifs_with_report(
                    get_tempfile("significant.pfm"),
                    get_tempfile("clustered.pfm"),
                    get_tempfile("."),
                    0.95)
        # Determine best motif in cluster
        #num_cluster, best_id = self._determine_best_motif_in_cluster(clusters, self.final_pwm, self.validation_fa, bg_file, self.imgdir)
        
        score_clustered_motifs(
                get_tempfile("clustered.pfm"),
                stats_fg=get_tempfile("validation.fa"), 
                stats_bg=background
                )
        
    if create_report:
        bg = dict([(b, get_tempfile("bg.{}.fa".format(b))) for b in background])

        create_denovo_motif_report(
                get_tempfile("clustered.pfm"), 
                get_tempfile("validation.fa"), 
                bg, 
                get_tempfile("validation.fa"), 
                "bla",
                params,
                )
    
    # print_params()
    # create_text_report()

    #if not(params["keep_intermediate"]):
    #    logger.debug(
    #            "Deleting intermediate files. "
    #            "Please specifify the -k option if you want to keep these files.")
    #    shutil.rmtree(self.tmpdir)

    #logger.info("finished")
    # TODO: fixme
    #logger.info("output dir: %s", "bla")#os.path.split(self.motif_report)[0])
    #logger.info("report: %s", os.path.split(self.motif_report)[-1])

    return self.motif_report

try:
    from gimmemotifs.mp import pool
except:
    pass
from gimmemotifs.comparison import MotifComparer
from gimmemotifs.cluster import cluster_motifs_with_report
from gimmemotifs.prediction import predict_motifs,PredictionResult
