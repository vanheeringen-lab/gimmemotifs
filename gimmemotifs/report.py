# Copyright (c) 2009-2017 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
"""Reports (graphical and text) for motifs statistics."""
import os
import sys
from datetime import datetime
from multiprocessing import Pool
from tempfile import NamedTemporaryFile
import shutil
import logging

import jinja2
import numpy as np

from gimmemotifs.comparison import MotifComparer
from gimmemotifs.fasta import Fasta
from gimmemotifs.motif import read_motifs
from gimmemotifs.config import GM_VERSION,MotifConfig,BG_RANK
from gimmemotifs.plot import roc_plot
from gimmemotifs.rocmetrics import roc_values
from gimmemotifs.stats import calc_stats, add_star, write_stats
from gimmemotifs import mytmpdir
from gimmemotifs.utils import motif_localization

logger = logging.getLogger("gimme.report")

def get_roc_values(motif, fg_file, bg_file):
    """Calculate ROC AUC values for ROC plots."""
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

def create_roc_plots(pwmfile, fgfa, background, outdir):
    """Make ROC plots for all motifs."""
    motifs = dict([(m.id, m) for m in read_motifs(open(pwmfile), fmt="pwm")])

    ncpus = int(MotifConfig().get_default_params()['ncpus'])
    pool = Pool(processes=ncpus)
    jobs = {}
    for bg,fname in background.items():
        for m_id, m in motifs.items():

            k = "{}_{}".format(str(m), bg)
            jobs[k] = pool.apply_async(
                                            get_roc_values,
                                            (motifs[m_id], fgfa, fname,)
                                            )
    imgdir = os.path.join(outdir, "images")
    if not os.path.exists(imgdir):
        os.mkdir(imgdir)
    
    roc_img_file = os.path.join(outdir, "images", "{}_roc.{}.png")

    for motif in motifs.values():
        for bg in background:
            k = "{}_{}".format(str(motif), bg)
            error, x, y = jobs[k].get()
            if error:
                logger.error("Error in thread: %s", error)
                sys.exit(1)
            roc_plot(roc_img_file.format(motif.id, bg), x, y)

def _create_text_report(inputfile, motifs, closest_match, stats, outdir):
    """Create text report of motifs with statistics and database match."""
    my_stats = {}
    for motif in motifs:
        match = closest_match[motif.id]
        my_stats[str(motif)] = {}
        for bg in stats.values()[0].keys():
            my_stats[str(motif)][bg] = stats[str(motif)][bg].copy()
            my_stats[str(motif)][bg]["best_match"] = "_".join(match[0].split("_")[:-1])
            my_stats[str(motif)][bg]["best_match_pvalue"] = match[1][-1]
    
    header = "# GimmeMotifs version {}\n".format(GM_VERSION) + 
             "# Inputfile: {}\n".format(inputfile)

    write_stats(my_stats, os.path.join(outdir, "stats.{}.txt"), header=header)

def _create_graphical_report(inputfile, pwm, background, closest_match, outdir, stats, best_id=None):
    """Create main gimme_motifs output html report."""
    if best_id is None:
        best_id = {}

    logger.debug("Creating graphical report")
    
    class ReportMotif(object):
        """Placeholder for motif stats."""
        pass

    config = MotifConfig()
    
    imgdir = os.path.join(outdir, "images")
    if not os.path.exists(imgdir):
        os.mkdir(imgdir)
    motifs = read_motifs(open(pwm), fmt="pwm")
    
    roc_img_file = "%s_roc.%s"

    dbpwm = config.get_default_params()["motif_db"]
    pwmdir = config.get_motif_dir()
    dbmotifs = dict([(m.id, m) for m in read_motifs(open(os.path.join(pwmdir, dbpwm)))])

    report_motifs = []
    for motif in motifs:
        
        rm = ReportMotif()
        rm.id = motif.id
        rm.id_href = {"href": "#%s" % motif.id}
        rm.id_name = {"name": motif.id}
        rm.img = {"src":  os.path.join("images", "%s.png" % motif.id)}
        motif.to_img(os.path.join(outdir, "images/{}.png".format(motif.id)), fmt="PNG")
        
        # TODO: fix best ID
        rm.best = "Gimme"#best_id[motif.id]

        rm.consensus = motif.to_consensus()
        rm.stars = int(np.mean(
                [stats[str(motif)][bg].get("stars", 0) for bg in background]
                ) + 0.5)

        rm.bg = {}
        for bg in background:
            rm.bg[bg] = {}
            this_stats = stats.get(str(motif), {}).get(bg)
            # TODO: fix these stats
            rm.bg[bg]["e"] = "%0.2f" % this_stats.get("enr_at_fdr", 1.0)
            rm.bg[bg]["p"] = "%0.2f" % this_stats.get("phyper_at_fd", 1.0)
            rm.bg[bg]["auc"] = "%0.3f" % this_stats.get("roc_auc", 0.5)
            rm.bg[bg]["mncp"] = "%0.3f" % this_stats.get("mncp", 1.0)
            rm.bg[bg]["roc_img"] = {"src": "images/" + os.path.basename(roc_img_file % (motif.id, bg)) + ".png"}
            rm.bg[bg]["roc_img_link"] = {"href": "images/" + os.path.basename(roc_img_file % (motif.id, bg)) + ".png"}

        rm.histogram_img = {"data":"images/%s_histogram.svg" % motif.id}
        rm.histogram_link= {"href":"images/%s_histogram.svg" % motif.id}
        
        match_id = closest_match[motif.id][0]
        dbmotifs[match_id].to_img(os.path.join(outdir, "images/{}.png".format(match_id)), fmt="PNG")
    
        rm.match_img = {"src":  "images/{}.png".format(match_id)}
        rm.match_id = closest_match[motif.id][0]
        rm.match_pval = "%0.2e" % closest_match[motif.id][1][-1]

        report_motifs.append(rm)
    total_report = os.path.join(outdir, "motif_report.html")

    star_img = os.path.join(config.get_template_dir(), "star.png")
    shutil.copyfile(star_img, os.path.join(outdir, "images", "star.png"))

    env = jinja2.Environment(loader=jinja2.FileSystemLoader([config.get_template_dir()]))
    template = env.get_template("report_template.jinja.html")
    # TODO: title
    result = template.render(
                    motifs=report_motifs, 
                    inputfile=inputfile, 
                    date=datetime.today().strftime("%d/%m/%Y"), 
                    version=GM_VERSION)

    f = open(total_report, "w")
    f.write(result.encode('utf-8'))
    f.close()

def create_denovo_motif_report(inputfile, pwmfile, fgfa, background, locfa, outdir, params, stats=None):
    """Create text and graphical (.html) motif reports."""
    logger.info("creating reports")

    # ROC plots
    create_roc_plots(pwmfile, fgfa, background, outdir)
        
    # Location plots
    logger.debug("Creating localization plots")
    motifs = read_motifs(open(pwmfile), fmt="pwm")
    
    mc = MotifComparer()
    closest_match = mc.get_closest_match(motifs)
    
    if stats is None:
        stats = {}
        for bg, bgfa in background.items():
            for m, s in calc_stats(motifs, fgfa, bgfa).items():
                if m not in stats:
                    stats[m] = {}
                stats[m][bg] = s

    stats = add_star(stats)

    if not params:
        params = {}
    cutoff_fdr = params.get('cutoff_fdr', 0.9)
    lwidth = int(params.get('lwidth', 500))

    for motif in motifs:
        outfile = os.path.join(outdir, "images/{}_histogram.svg".format(motif.id))
        motif_localization(locfa, motif, lwidth, outfile, cutoff=cutoff_fdr)

    # Create reports
    _create_text_report(inputfile, motifs, closest_match, stats, outdir)
    _create_graphical_report(inputfile, pwmfile, background, closest_match, outdir, stats)
