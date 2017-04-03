# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
from gimmemotifs.config import MotifConfig, BG_RANK
from gimmemotifs import mytmpdir
from gimmemotifs.validation import check_denovo_input
from gimmemotifs.utils import divide_fa_file, motif_localization
from gimmemotifs.fasta import Fasta
from gimmemotifs.background import ( MarkovFasta, MatchedGcFasta,
                                    PromoterFasta, RandomGenomicFasta )
from gimmemotifs.plot import roc_plot

from gimmemotifs.rocmetrics import roc_values, roc_auc, mncp, max_fmeasure
from tempfile import NamedTemporaryFile
from time import sleep

import os
import sys
import logging

## for clustering
import jinja2
from datetime import datetime
from gimmemotifs.motif import read_motifs
from gimmemotifs.config import GM_VERSION
#
logger = logging.getLogger()

def get_roc_values(motif, fg_file, bg_file):
    try:
        fg_result = motif.pwm_scan_score(Fasta(fg_file), cutoff=0.0, nreport=1)
        fg_vals = [sorted(x)[-1] for x in fg_result.values()]

        bg_result = motif.pwm_scan_score(Fasta(bg_file), cutoff=0.0, nreport=1)
        bg_vals = [sorted(x)[-1] for x in bg_result.values()]

        (x, y) = roc_values(fg_vals, bg_vals)
        return None,x,y
    except Exception,e:
        error = e
        return error,[],[]

def parse_denovo_params(user_params=None):
    config = MotifConfig()
    
    if user_params is None:
        user_params = {}
    params = config.get_default_params()
    params.update(user_params)

    if params.get("torque"):
        from gimmemotifs.prediction_torque import pp_predict_motifs, PredictionResult
        logger.debug("Using torque")
    else:
        
        print "import"
        from gimmemotifs.prediction import pp_predict_motifs, PredictionResult
        logger.debug("Using multiprocessing")

    params["background"] = [x.strip() for x in params["background"].split(",")]

    logger.debug("Parameters:")
    for param, value in params.items():
        logger.debug("  %s: %s", param, value)

    # Maximum time?
    if params["max_time"]:
        try:
            params["max_time"] = float(params["max_time"])
        except Exception:
            logger.debug("Could not parse max_time value, setting to no limit")
            params["max_time"] = None

        if params["max_time"] > 0:
            logger.debug("Time limit for motif prediction: %0.2f hours" % max_time)
            params["max_time"] = 3600 * params["max_time"]
            logger.debug("Max_time in seconds %0.0f" % params["max_time"])
        else:
            logger.debug("Invalid time limit for motif prediction, setting to no limit")
            max_time = params["max_time"]
    else:
            logger.debug("No time limit for motif prediction")
    
    return params

def prepare_denovo_input_bed():
    self.prepare_input_bed(inputfile, params["genome"], params["width"], params["fraction"], params["abs_max"], params["use_strand"])


    # Create file for location plots
    index_dir = os.path.join(self.config.get_index_dir(), params["genome"])
    lwidth = int(params["lwidth"])
    width = int(params["width"])
    extend = (lwidth - width) / 2
    genome_index.track2fasta(index_dir, self.validation_bed, self.location_fa, extend_up=extend, extend_down=extend, use_strand=params["use_strand"], ignore_missing=True)

def predict_motifs(background, params):
    """ Predict motifs, input is a FASTA-file"""
    tools = dict(
            [
                (x.strip(), x in [y.strip() for y in  params["tools"].split(",")]) 
                    for x in params["available_tools"].split(",")
            ]
            )

    # Predict the motifs
    analysis = params["analysis"]
    logger.info("starting motif prediction (%s)", analysis)
    logger.info("tools: %s", 
            ", ".join([x for x in tools.keys() if tools[x]]))

    print background
    bg_name = sorted(background, lambda x,y: cmp(BG_RANK[x], BG_RANK[y]))[0]
    bg_file = get_tempfile("bg.{}.fa".format(bg_name))
    logger.debug("Using bg_file %s for significance" % bg_file)
    
    pfm_file = get_tempfile("all_motifs.pfm")

    result = pp_predict_motifs(
                    get_tempfile("prediction.fa"), 
                    pfm_file, 
                    analysis, 
                    params["genome"], 
                    params["use_strand"], 
                    get_tempfile("prediction.bg.fa"), 
                    tools, 
                    None, 
                    logger=logger, 
                    max_time=params["max_time"], 
                    stats_fg=get_tempfile("validation.fa"), 
                    stats_bg=background
                )

    motifs = result.motifs
    logger.info("predicted %s motifs", len(motifs))
    logger.debug("written to %s", pfm_file)

    if len(motifs) == 0:
        logger.info("no motifs found")
        sys.exit()
    
    return result

def write_stats():
    # Write stats output to file
    f = open(self.stats_file, "w")
    stat_keys = result.stats.values()[0].keys()
    f.write("%s\t%s\n" % ("Motif", "\t".join(stat_keys)))
    
    self.logger.debug(result.stats)
    
    for motif in motifs:
        stats = result.stats.get("%s_%s" % (motif.id, motif.to_consensus()), None)
        if stats:
            f.write("%s\t%s\n" % (motif.id, "\t".join([str(stats[k]) for k in stat_keys])))
        else:
            self.logger.warn("No stats for motif {0}, skipping this motif!".format(motif.id))
            #motifs.remove(motif)
    f.close()

    self.motifs_with_stats = motifs

    f = open(self.ranks_file, "w")
    tools = dict((m.id.split("_")[2],1) for m in motifs).keys()
    f.write("Metric\tType\t%s\n" % ("\t".join(tools)))
    for stat in ["mncp", "roc_auc", "maxenr"]:
        best_motif = {}
        for motif in self.motifs_with_stats:
            d = result.stats.get("%s_%s" % (motif.id, motif.to_consensus()), {})
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
        print stats
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

def create_roc_plots(pwmfile, fgfa, background, name, outdir): 
    motifs = dict([(m.id, m) for m in read_motifs(open(pwmfile), fmt="pwm")]) 
 
    jobs = {} 
    for bg,fname in background.items():
        for m_id, m in motifs.items(): 
            
            k = "{}_{}".format(str(m), bg)
            jobs[k] = pool.apply_async(
                                            get_roc_values, 
                                            (motifs[m_id], fgfa, fname,)
                                            )
 
    roc_img_file = os.path.join(outdir, "images", "{}_roc.{}.png") 
 
    for motif in motifs.values():
        for bg in background:
            k = "{}_{}".format(str(motif), bg)
            error, x, y = jobs[k].get() 
            if error: 
                logger.error("Error in thread: %s", error) 
                sys.exit(1) 
            roc_plot(roc_img_file.format(motif.id, bg), x, y) 

def _create_report(pwm, background, outdir, stats=None, best_id=None, closest_match=None):
    if stats is None:
        stats = {}
    if best_id is None:
        best_id = {}

    logger.debug("Creating graphical report")
    
    class ReportMotif:
        pass

    imgdir = os.path.join(outdir, "images")
    motifs = read_motifs(open(pwm), fmt="pwm")
    
    if closest_match:
        for m,match in self.closest_match.items():
            match[0].to_img(os.path.join(imgdir,"%s.png" % match[0].id), fmt="PNG")

    sort_key = sorted(background, lambda x,y: cmp(BG_RANK[x], BG_RANK[y]))[0]

    roc_img_file = "%s_roc.%s"

    # TODO: Implement different backgrounds
    sorted_motifs = sorted(motifs,
            cmp=lambda x,y: cmp(stats[str(y)], stats[str(x)])
            )

    report_motifs = []
    for motif in sorted_motifs:
        
        rm = ReportMotif()
        rm.id = motif.id
        rm.id_href = {"href": "#%s" % motif.id}
        rm.id_name = {"name": motif.id}
        rm.img = {"src":  os.path.join("images", "%s.png" % motif.id)}
        
        # TODO: fix best ID
        rm.best = "Gimme"#best_id[motif.id]

        rm.consensus = motif.to_consensus()
        rm.stars = stats[str(motif)].get("stars", 0)

        rm.bg = {}
        for bg in background:
            rm.bg[bg] = {}
            this_stats = stats.get(str(motif), {}).get(bg)
            print stats.get(str(motif), 0.0)
            # TODO: fix these stats
            rm.bg[bg]["e"] = "%0.2f" % this_stats.get("e", 1.0)
            rm.bg[bg]["p"] = "%0.2f" % this_stats.get("p", 1.0)
            rm.bg[bg]["auc"] = "%0.3f" % this_stats.get("roc_auc", 0.5)
            rm.bg[bg]["mncp"] = "%0.3f" % this_stats.get("mncp", 1.0)
            rm.bg[bg]["roc_img"] = {"src": "images/" + os.path.basename(roc_img_file % (motif.id, bg)) + ".png"}
            rm.bg[bg]["roc_img_link"] = {"href": "images/" + os.path.basename(roc_img_file % (motif.id, bg)) + ".png"}

        rm.histogram_img = {"data":"images/%s_histogram.svg" % motif.id}
        rm.histogram_link= {"href":"images/%s_histogram.svg" % motif.id}
        rm.match_img = ""#{"src":  "images/%s.png" % self.closest_match[motif.id][0].id}
        rm.match_id = ""#self.closest_match[motif.id][0].id
        rm.match_pval = "1.0"#"%0.2e" % self.closest_match[motif.id][1]

        report_motifs.append(rm)
    total_report = os.path.join(outdir, "motif_report.html")

    config = MotifConfig()
    env = jinja2.Environment(loader=jinja2.FileSystemLoader([config.get_template_dir()]))
    template = env.get_template("report_template.jinja.html")
    # TODO: title
    result = template.render(
                    expname="bla", 
                    motifs=report_motifs, 
                    inputfile="bla", 
                    date=datetime.today().strftime("%d/%m/%Y"), 
                    version=GM_VERSION)

    f = open(total_report, "w")
    f.write(result.encode('utf-8'))
    f.close()

def create_denovo_motif_report(pwmfile, fgfa, background, locfa, outdir, params):
    logger.info("creating report")

    # ROC plots
    create_roc_plots(pwmfile, fgfa, background, "bla",  'bla')
        
    # Location plots
    logger.debug("Creating localization plots")
    motifs = read_motifs(open(pwmfile), fmt="pwm")
    
    tmp = NamedTemporaryFile(dir=mytmpdir()).name
    p = PredictionResult(
            tmp, 
            logger=logger, 
            job_server=pool, 
            fg_file=fgfa, 
            background=background, 
            do_counter=False)
    p.add_motifs(("clustering",  (read_motifs(open(pwmfile)), "","")))
    while len(p.stats.keys()) < len(p.motifs):
        sleep(5)

    #for mid, num in num_cluster.items():
    #    p.stats[mid]["numcluster"] = num
  
    all_stats = {
            "mncp": [2, 5, 8],
            "roc_auc": [0.6, 0.75, 0.9],
            "maxenr": [10, 20, 30],
            "enr_fdr": [4, 8, 12],
            "fraction": [0.4, 0.6, 0.8],
            "ks_sig": [4, 7, 10],
            "numcluster": [3, 6, 9],
    }

 
    for k,v in p.stats.items():
        print k,v
    
    if not params:
        params = {}
    cutoff_fdr = params.get('cutoff_fdr', 0.9)
    lwidth = int(params.get('lwidth', 500))

    for motif in motifs:
        s = p.stats[str(motif)]
        outfile = os.path.join(outdir, "images/{}_histogram.svg".format(motif.id))
        motif_localization(locfa, motif, lwidth, outfile, cutoff=cutoff_fdr)

        #s["stars"] = int(mean([star(s[x], all_stats[x]) for x in all_stats.keys()]) + 0.5)
        #self.logger.debug("Motif %s: %s stars" % (m, s["stars"]))

    # Create report
    #self.print_params()
    #self._calc_report_values(pwmfile, background)
    _create_report(pwmfile, background, "bla", stats=p.stats)
    #self._create_text_report(pwmfile, background)
    
    logger.info("finished")
    # TODO: fixme
    logger.info("output dir: %s", "bla")#os.path.split(self.motif_report)[0])
    #logger.info("report: %s", os.path.split(self.motif_report)[-1])

def cluster_motifs_with_report(infile, outfile, outdir, threshold):
    # Cluster significant motifs
    
    motifs = read_motifs(open(infile), fmt="pwm")

    trim_ic = 0.2
    clusters = []
    if len(motifs) == 0:
        return []
    elif len(motifs) == 1:
        clusters = [[motifs[0], motifs]]
    else:
        logger.info("clustering significant motifs.")
        tree = cluster_motifs(
                infile,
                "total",
                "wic",
                "mean",
                True,
                threshold=float(threshold),
                include_bg=True,
                progress=False
                )
        clusters = tree.getResult()

    ids = []
    mc = MotifComparer()

    img_dir = os.path.join(outdir, "images")
    
    if not os.path.exists(img_dir):
        os.mkdir(img_dir)

    for cluster,members in clusters:
        cluster.trim(trim_ic)
        png = "images/{}.png".format(cluster.id)
        cluster.to_img(os.path.join(outdir, png), fmt="PNG")
        ids.append([cluster.id, {"src":png},[]])
        if len(members) > 1:
            scores = {}
            for motif in members:
                scores[motif] =  mc.compare_motifs(cluster, motif, "total", "wic", "mean", pval=True)
            add_pos = sorted(scores.values(),cmp=lambda x,y: cmp(x[1], y[1]))[0][1]
            for motif in members:
                score, pos, strand = scores[motif]
                add = pos - add_pos

                if strand in [1,"+"]:
                    pass
                else:
                   rc = motif.rc()
                   rc.id = motif.id
                   motif = rc
                #print "%s\t%s" % (motif.id, add)
                png = "images/{}.png".format(cluster.id.replace(" ", "_"))
                motif.to_img(os.path.join(outdir, png), fmt="PNG", add_left=add)
        ids[-1][2] = [dict([("src", png), ("alt", motif.id.replace(" ", "_"))]) for motif in members]

    expname = "test"
    inputfile = "input"
    config = MotifConfig()
    env = jinja2.Environment(loader=jinja2.FileSystemLoader([config.get_template_dir()]))
    template = env.get_template("cluster_template.jinja.html")
    result = template.render(
                expname=expname, 
                motifs=ids, 
                inputfile=inputfile, 
                date=datetime.today().strftime("%d/%m/%Y"), 
                version=GM_VERSION)

    cluster_report = os.path.join(outdir, "cluster_report.html")
    f = open(cluster_report, "w")
    f.write(result.encode('utf-8'))
    f.close()

    f = open(outfile, "w")
    if len(clusters) == 1 and len(clusters[0][1]) == 1:
        f.write("%s\n" % clusters[0][0].to_pwm())
    else:
        for motif in tree.get_clustered_motifs():
            f.write("%s\n" % motif.to_pwm())
    f.close()

    logger.debug("Clustering done. See the result in %s",
            cluster_report)
    return clusters


    # Determine best motif in cluster
    
    num_cluster, best_id = self._determine_best_motif_in_cluster(clusters, self.final_pwm, self.validation_fa, bg_file, self.imgdir)
 

def score_motifs():
    # Stars
    tmp = NamedTemporaryFile(dir=mytmpdir()).name
    p = PredictionResult(tmp, logger=self.logger, job_server=self.server, fg_file = self.validation_fa, bg_file = bg_file, do_counter=False)
    p.add_motifs(("clustering",  (read_motifs(open(self.final_pwm)), "","")))
    while len(p.stats.keys()) < len(p.motifs):
        sleep(5)

    #print "p.stats"
    #print p.stats
    #print "num_cluster"
    #print num_cluster
    for mid, num in num_cluster.items():
        p.stats[mid]["numcluster"] = num

    all_stats = {
        "mncp": [2, 5, 8],
        "roc_auc": [0.6, 0.75, 0.9],
        "maxenr": [10, 20, 30],
        "enr_fdr": [4, 8, 12],
        "fraction": [0.4, 0.6, 0.8],
        "ks_sig": [4, 7, 10],
        "numcluster": [3, 6, 9],
    }



def predict_denovo_motifs(inputfile, params=None, filter_significant=True, cluster=True, create_report=True):
    """ Full analysis: from bed-file to motifs (including clustering, ROC-curves, location plots and html report) """
    logger.info("starting full motif analysis")
    logger.debug("Using temporary directory {0}".format(mytmpdir()))

    # Initialize parameters
    params = parse_denovo_params(params)
 
    # Check the input files
    input_type, background = check_denovo_input(inputfile, params)
    
    
    if not os.path.exists("bla"):
        os.mkdir('bla') 

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
    
    print "HOEI", background
    # Predict de novo motifs
    result = predict_motifs(background, params)

    # Write statistics
    #write_stats()

    bg = sorted(background, lambda x,y: cmp(BG_RANK[x], BG_RANK[y]))[0]
    if filter_significant:
        motifs = filter_significant_motifs(result, bg)
    else:
        motifs = result.motifs

#    print motifs
    
    if cluster: 
        clusters = cluster_motifs_with_report(
                    get_tempfile("significant.pfm"),
                    get_tempfile("clustered.pfm"),
                    get_tempfile("."),
                    0.95)
  
    #print motifs
    
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
    
    sys.exit()
    score_motifs()
    
    if not(params["keep_intermediate"]):
        logger.debug(
                "Deleting intermediate files. "
                "Please specifify the -k option if you want to keep these files.")
        shutil.rmtree(self.tmpdir)

    logger.debug("Done")

    return self.motif_report
try:
    from gimmemotifs.mp import pool
except:
    pass
from gimmemotifs.comparison import MotifComparer
from gimmemotifs.cluster import cluster_motifs
from gimmemotifs.prediction import pp_predict_motifs,PredictionResult
