#!/usr/bin/env python
# Copyright (c) 2009-2016 Simon van Heeringen <s.vanheeringen@science.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

from distutils import sysconfig

from gimmemotifs.motif import pwmfile_to_motifs
from gimmemotifs.comparison import MotifComparer
from gimmemotifs.cluster import cluster_motifs
from gimmemotifs.config import MotifConfig
import sys
import os

import jinja2

def cluster(args):

    revcomp = not args.single

    outdir = os.path.abspath(args.outdir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    trim_ic = 0.2
    clusters = []
    motifs = pwmfile_to_motifs(args.inputfile)
    if len(motifs) == 1:
        clusters = [[motifs[0], motifs]]
    else:
        tree = cluster_motifs(args.inputfile, "total", "wic", "mean", True, threshold=args.threshold, include_bg=True)
        clusters = tree.getResult()
    
    ids = []
    mc = MotifComparer()

    sys.stderr.write("Creating images\n")
    for cluster,members in clusters:
        cluster.trim(trim_ic)
        cluster.to_img(os.path.join(outdir,"%s.png" % cluster.id), format="PNG")
        ids.append([cluster.id, {"src":"%s.png" % cluster.id},[]])
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
                    #print "RC %s" % motif.id
                    rc = motif.rc()
                    rc.id = motif.id
                    motif = rc
                #print "%s\t%s" % (motif.id, add)    
                motif.to_img(os.path.join(outdir, "%s.png" % motif.id.replace(" ", "_")), format="PNG", add_left=add)
        ids[-1][2] = [dict([("src", "%s.png" % motif.id.replace(" ", "_")), ("alt", motif.id.replace(" ", "_"))]) for motif in members]
    
    config = MotifConfig()
    env = jinja2.Environment(loader=jinja2.FileSystemLoader([config.get_template_dir()]))
    template = env.get_template("cluster_template.jinja.html")
    result = template.render(motifs=ids)

    with open(os.path.join(outdir, "cluster_report.html"), "w") as f:
        f.write(result.encode('utf-8'))

    f = open(os.path.join(outdir, "cluster_key.txt"), "w")
    for id in ids:
        f.write("%s\t%s\n" % (id[0], ",".join([x["alt"] for x in id[2]])))
    f.close()

    f = open(os.path.join(outdir, "clustered_motifs.pwm"), "w")
    if len(clusters) == 1 and len(clusters[0][1]) == 1:
        f.write("%s\n" % clusters[0][0].to_pwm())
    else:
        for motif in tree.get_clustered_motifs():
            f.write("%s\n" % motif.to_pwm())
    f.close()
