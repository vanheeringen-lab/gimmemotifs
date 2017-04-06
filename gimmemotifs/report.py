# Copyright (c) 2009-2017 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
import os
import sys
from datetime import datetime
from multiprocessing import Pool
from tempfile import NamedTemporaryFile
import shutil

import jinja2
import numpy as np

from gimmemotifs.comparison import MotifComparer
from gimmemotifs.fasta import Fasta
from gimmemotifs.motif import read_motifs
from gimmemotifs.config import GM_VERSION,MotifConfig,BG_RANK
from gimmemotifs.plot import roc_plot
from gimmemotifs.rocmetrics import roc_values
from gimmemotifs.stats import calc_stats, add_star
from gimmemotifs import mytmpdir
from gimmemotifs.utils import motif_localization

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

def create_roc_plots(pwmfile, fgfa, background, name, outdir):
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

    roc_img_file = os.path.join(outdir, "images", "{}_roc.{}.png")

    for motif in motifs.values():
        for bg in background:
            k = "{}_{}".format(str(motif), bg)
            error, x, y = jobs[k].get()
            if error:
                #logger.error("Error in thread: %s", error)
                print("Error in thread: %s" % error)
                sys.exit(1)
            roc_plot(roc_img_file.format(motif.id, bg), x, y)

def _create_report(pwm, background, outdir, stats=None, best_id=None, closest_match=None):
    if stats is None:
        stats = {}
    if best_id is None:
        best_id = {}

    #logger.debug("Creating graphical report")
    
    class ReportMotif:
        pass

    config = MotifConfig()
    
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

    mc = MotifComparer()
    closest_match = mc.get_closest_match(motifs)

    dbpwm = config.get_default_params()["motif_db"]
    pwmdir = config.get_motif_dir()
    dbmotifs = dict([(m.id, m) for m in read_motifs(open(os.path.join(pwmdir, dbpwm)))])

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
        dbmotifs[match_id].to_img(os.path.join(outdir, "images/{}png".format(match_id)), fmt="PNG")
    
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
                    expname="bla", 
                    motifs=report_motifs, 
                    inputfile="bla", 
                    date=datetime.today().strftime("%d/%m/%Y"), 
                    version=GM_VERSION)

    f = open(total_report, "w")
    f.write(result.encode('utf-8'))
    f.close()

def create_denovo_motif_report(pwmfile, fgfa, background, locfa, outdir, params):
    #logger.info("creating report")

    # ROC plots
    create_roc_plots(pwmfile, fgfa, background, "bla",  'bla')
        
    # Location plots
    #logger.debug("Creating localization plots")
    motifs = read_motifs(open(pwmfile), fmt="pwm")
    
    tmp = NamedTemporaryFile(dir=mytmpdir()).name
    
    stats = {}
    for bg, bgfa in background.items():
        for m, s in calc_stats(motifs, fgfa, bgfa).items():
            if not m in stats:
                stats[m] = {}
            stats[m][bg] = s

    stats = add_star(stats)

    if not params:
        params = {}
    cutoff_fdr = params.get('cutoff_fdr', 0.9)
    lwidth = int(params.get('lwidth', 500))

    for motif in motifs:
        s = stats[str(motif)]
        outfile = os.path.join(outdir, "images/{}_histogram.svg".format(motif.id))
        motif_localization(locfa, motif, lwidth, outfile, cutoff=cutoff_fdr)

    # Create report
    _create_report(pwmfile, background, "bla", stats=stats)
