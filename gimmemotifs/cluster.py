# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
"""Module for motif clustering."""
import os
import sys
import logging

import jinja2
from datetime import datetime

# GimmeMotifs imports
from gimmemotifs.config import MotifConfig, GM_VERSION
from gimmemotifs.motif import read_motifs
from gimmemotifs.comparison import MotifComparer

logger = logging.getLogger("gimme.cluster")

class MotifTree(object):
    
    """Class MotifTree used by cluster_motifs"""
    
    def __init__(self, motif):
        self.motif = motif
        self.parent = None
        self.left = None
        self.right = None
        self.mergescore = None
        self.maxscore = 0
        self.frontier = False

    def setFrontier(self, arg, root):
        self.frontier = True 
        if self != root:
            self.parent.setFrontier(True, root)

    def checkMerge(self, root, threshold):
        if not self.frontier:# and self != root:
            if self.mergescore > threshold:
                if self.parent:
                    self.parent.checkMerge(root, threshold)
            else:
                self.setFrontier(True, root)
    
    def printFrontiers(self):
        if self.frontier:
            if self.left:
                self.left.printFrontiers()
                self.right.printFrontiers()
        else:
            self.motif.pwm = self.motif.pfm_to_pwm(self.motif.pwm)

    def get_clustered_motifs(self):
        if self.frontier:
            if self.left:
                return self.left.get_clustered_motifs() + self.right.get_clustered_motifs()
        else:
            return [self.motif]
    
    def getResult(self):
        if self.frontier:
            if self.left:
                return self.left.getResult() + self.right.getResult()
        else:
            return [[self.motif, self.recursive_motif()]]


    def recursive_name(self):
        if self.left:
            return self.left.recursive_name() + self.right.recursive_name()
        else:
            return [self.motif.id]

    def recursive_motif(self):
        if self.left:
            return self.left.recursive_motif() + self.right.recursive_motif()
        else:
            return [self.motif]

def cluster_motifs(motifs, match="total", metric="wic", combine="mean", pval=True, threshold=0.95, trim_edges=False, edge_ic_cutoff=0.2, include_bg=True, progress=True):
    """ 
    Clusters a set of sequence motifs. Required arg 'motifs' is a file containing
    positional frequency matrices or an array with motifs.

    Optional args:

    'match', 'metric' and 'combine' specify the method used to compare and score
    the motifs. By default the WIC score is used (metric='wic'), using the the
    score over the whole alignment (match='total'), with the total motif score
    calculated as the mean score of all positions (combine='mean').
    'match' can be either 'total' for the total alignment or 'subtotal' for the 
    maximum scoring subsequence of the alignment.
    'metric' can be any metric defined in MotifComparer, currently: 'pcc', 'ed',
    'distance', 'wic' or 'chisq' 
    'combine' determines how the total score is calculated from the score of 
    individual positions and can be either 'sum' or 'mean'
    
    'pval' can be True or False and determines if the score should be converted to 
    an empirical p-value

    'threshold' determines the score (or p-value) cutoff

    If 'trim_edges' is set to True, all motif edges with an IC below 
    'edge_ic_cutoff' will be removed before clustering

    When computing the average of two motifs 'include_bg' determines if, at a 
    position only present in one motif, the information in that motif should
    be kept, or if it should be averaged with background frequencies. Should
    probably be left set to True.

    """

    
    # First read pfm or pfm formatted motiffile
    if type([]) != type(motifs):
        with open(motifs) as f:
            motifs = read_motifs(f, fmt="pwm")
    
    mc = MotifComparer()

    # Trim edges with low information content
    if trim_edges:
        for motif in motifs:
            motif.trim(edge_ic_cutoff)
    
    # Make a MotifTree node for every motif
    nodes = [MotifTree(m) for m in motifs]
    
    # Determine all pairwise scores and maxscore per motif
    scores = {}
    motif_nodes = dict([(n.motif.id,n) for n in nodes])
    motifs = [n.motif for n in nodes]
    
    if progress:
        sys.stderr.write("Calculating initial scores\n")
    result = mc.get_all_scores(motifs, motifs, match, metric, combine, pval, parallel=True)
    
    for m1, other_motifs in result.items():
        for m2, score in other_motifs.items():
            if m1 == m2:
                if pval:
                    motif_nodes[m1].maxscore = 1 - score[0]
                else:
                    motif_nodes[m1].maxscore = score[0]
            else:
                if pval:
                    score = [1 - score[0]] + score[1:]
                scores[(motif_nodes[m1],motif_nodes[m2])] = score
               
    cluster_nodes = [node for node in nodes]
    ave_count = 1
    
    total = len(cluster_nodes)

    while len(cluster_nodes) > 1:
        l = sorted(scores.keys(), key=lambda x: scores[x][0])
        i = -1
        (n1, n2) = l[i]
        while n1 not in cluster_nodes or n2 not in cluster_nodes:
            i -= 1
            (n1,n2) = l[i]
        
        (score, pos, orientation) = scores[(n1,n2)]
        ave_motif = n1.motif.average_motifs(n2.motif, pos, orientation, include_bg=include_bg)
        
        ave_motif.trim(edge_ic_cutoff)
        ave_motif.id = "Average_%s" % ave_count
        ave_count += 1
        
        new_node = MotifTree(ave_motif)
        if pval:
            new_node.maxscore = 1 - mc.compare_motifs(new_node.motif, new_node.motif, match, metric, combine, pval)[0]
        else:
            new_node.maxscore = mc.compare_motifs(new_node.motif, new_node.motif, match, metric, combine, pval)[0]
            
        new_node.mergescore = score
        #print "%s + %s = %s with score %s" % (n1.motif.id, n2.motif.id, ave_motif.id, score)
        n1.parent = new_node
        n2.parent = new_node
        new_node.left = n1
        new_node.right = n2
        
        cmp_nodes = dict([(node.motif, node) for node in nodes if not node.parent])
        
        if progress:
            progress = (1 - len(cmp_nodes) / float(total)) * 100
            sys.stderr.write('\rClustering [{0}{1}] {2}%'.format(
                '#' * (int(progress) // 10), 
                " " * (10 - int(progress) // 10), 
                int(progress)))
        
        result = mc.get_all_scores(
                [new_node.motif], 
                list(cmp_nodes.keys()), 
                match, 
                metric, 
                combine, 
                pval, 
                parallel=True)
        
        for motif, n in cmp_nodes.items():
            x = result[new_node.motif.id][motif.id]
            if pval:
                x = [1 - x[0]] + x[1:]
            scores[(new_node, n)] = x
        
        nodes.append(new_node)

        cluster_nodes = [node for node in nodes if not node.parent]
     
    if progress:
        sys.stderr.write("\n") 
    root = nodes[-1]
    for node in [node for node in nodes if not node.left]:
        node.parent.checkMerge(root, threshold)
    
    return root

def cluster_motifs_with_report(infile, outfile, outdir, threshold, title=None):
    # Cluster significant motifs

    if title is None:
        title = infile

    with open(infile) as f:
        motifs = read_motifs(f, fmt="pwm")

    trim_ic = 0.2
    clusters = []
    if len(motifs) == 0:
        return []
    elif len(motifs) == 1:
        clusters = [[motifs[0], motifs]]
    else:
        logger.info("clustering %d motifs.", len(motifs))
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
            add_pos = sorted(scores.values(),key=lambda x: x[1])[0][1]
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
                png = "images/{}.png".format(motif.id.replace(" ", "_"))
                motif.to_img(os.path.join(outdir, png), fmt="PNG", add_left=add)
        ids[-1][2] = [dict([("src", "images/{}.png".format(motif.id.replace(" ", "_"))), ("alt", motif.id.replace(" ", "_"))]) for motif in members]

    config = MotifConfig()
    env = jinja2.Environment(loader=jinja2.FileSystemLoader([config.get_template_dir()]))
    template = env.get_template("cluster_template.jinja.html")
    result = template.render(
                motifs=ids,
                inputfile=title,
                date=datetime.today().strftime("%d/%m/%Y"),
                version=GM_VERSION)

    cluster_report = os.path.join(outdir, "cluster_report.html")
    with open(cluster_report, "wb") as f:
        f.write(result.encode('utf-8'))

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
