#!/usr/bin/env python
# Copyright (c) 2009-2019 Simon van Heeringen <s.vanheeringen@science.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
import logging
import os

import jinja2

from gimmemotifs.config import MotifConfig
from gimmemotifs.motif import read_motifs
from gimmemotifs.motif.cluster import cluster_motifs

from gimmemotifs.comparison import MotifComparer  # isort: skip

logger = logging.getLogger("gimme.cluster")


def _write_report(outdir, ids, tree, clusters):
    config = MotifConfig()
    env = jinja2.Environment(
        loader=jinja2.FileSystemLoader([config.get_template_dir()])
    )
    template = env.get_template("cluster_template.jinja.html")
    result = template.render(motifs=ids)

    with open(os.path.join(outdir, "gimme.clustered.html"), "w") as f:
        f.write(result)

    with open(os.path.join(outdir, "cluster_key.txt"), "w") as f:
        for motif_id in ids:
            f.write(f"{motif_id[0]}\t{','.join([x['alt'] for x in motif_id[2]])}\n")

    with open(os.path.join(outdir, "clustered_motifs.pfm"), "w") as f:
        if len(clusters) == 1 and len(clusters[0][1]) == 1:
            f.write(f"{clusters[0][0].to_ppm()}\n")
        else:
            for motif in tree.get_clustered_motifs():
                f.write(f"{motif.to_ppm()}\n")


def _create_images(outdir, clusters):
    ids = []
    mc = MotifComparer()
    trim_ic = 0.2

    logger.info("Creating images")
    for cluster, members in clusters:
        cluster.trim(trim_ic)
        cluster.plot_logo(fname=os.path.join(outdir, f"{cluster.id}.png"))
        ids.append([cluster.id, {"src": f"{cluster.id}.png"}, []])
        if len(members) > 1:
            scores = {}
            for motif in members:
                scores[motif] = mc.compare_motifs(
                    cluster, motif, "total", "wic", "mean", pval=True
                )
            add_pos = sorted(scores.values(), key=lambda x: x[1])[0][1]
            for motif in members:
                _, pos, strand = scores[motif]
                add = pos - add_pos

                if strand in [1, "+"]:
                    pass
                else:
                    rc = motif.rc()
                    rc.id = motif.id
                    motif = rc
                motif.plot_logo(
                    fname=os.path.join(outdir, f"{motif.id.replace(' ', '_')}.png"),
                    add_left=add,
                )
        ids[-1][2] = [
            dict(
                [
                    ("src", f"{m.id.replace(' ', '_')}.png"),
                    ("alt", m.id.replace(" ", "_")),
                ]
            )
            for m in members
        ]
    return ids


def cluster(args):

    outdir = os.path.abspath(args.outdir)
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    ncpus = args.ncpus

    clusters = []
    motifs = read_motifs(args.inputfile)
    if len(motifs) == 1:
        clusters = [[motifs[0], motifs]]
    else:
        tree = cluster_motifs(
            args.inputfile,
            "total",
            "wic",
            "mean",
            True,
            threshold=args.threshold,
            include_bg=True,
            ncpus=ncpus,
        )
        clusters = tree.getResult()

    ids = _create_images(outdir, clusters)
    _write_report(outdir, ids, tree, clusters)
