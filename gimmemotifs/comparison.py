# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

""" 
Module to compare DNA sequence motifs (positional frequency matrices)
"""

# Python imports
import sys
import os
from time import sleep
import random
import logging

# External imports
from scipy.stats import norm,entropy,chi2_contingency
from scipy.spatial import distance
from numpy import mean,std,array,sum

# GimmeMotifs imports
from gimmemotifs.motif import *
from gimmemotifs.config import *
from gimmemotifs.c_metrics import pwmscan,score
# pool import is at the bottom

# Try to import the fisim code, if it present
try:
    fisim_dir = "/usr/share/gimmemotifs/includes/fisim"
    if os.path.exists(fisim_dir):
        sys.path.append(fisim_dir)
        import fisim.Motif
except ImportError:
    pass

# Function that can be parallelized
def _get_all_scores(mc, motifs, dbmotifs, match, metric, combine, pval):
    try:
        scores = {}
        for m1 in motifs:
            scores[m1.id] = {}
            for m2 in dbmotifs:
                scores[m1.id][m2.id] = mc.compare_motifs(m1, m2, match, metric, combine, pval=pval)    
        return scores
    except Exception:
        logging.exception("_get_all_scores failed")


def akl(x,y):
    """ Calcylates the average Kullback-Leibler similarity
    See Mahony,  2007
    """
    return 10 - (entropy(x,y) + entropy(y,x)) / 2.0

def chisq(x,y):
    return chi2_contingency([x, y])[1]

def ssd(x,y):
    """ Sum of squared distances 
    """
    return sum([(a-b)**2 for a,b in zip(x,y)] )

def seqcor(m1,m2):
    l1 = len(m1)
    l2 = len(m2)

    l = max(l1, l2)

    # Create random sequence
    nucs = []
    L = 10 ** 4
    for i in range(L):
        nucs.append(random.choice(['A', 'C', 'T', 'G']))
    random_seq = "".join(nucs)

    # Scan random sequence
    result1 = pwmscan(random_seq.upper(), m1.pwm, m1.pwm_min_score(), len(random_seq), False, True)
    result2 = pwmscan(random_seq.upper(), m2.pwm, m2.pwm_min_score(), len(random_seq), False, True)

    # Return maximum correlation
    c = []
    for i in range(l1):
        c.append(1 - distance.correlation(result1[:L-l-i],result2[i:L-l]))
    for i in range(l2):
        c.append(1 - distance.correlation(result1[i:L-l],result2[:L-l-i]))
    return max(c)


class MotifComparer:
    def __init__(self):
        self.config = MotifConfig()
        self.metrics = ["pcc", "ed", "distance", "wic", "fisim"]
        self.combine = ["mean", "sum"]
        self._load_scores()
        # Create a parallel python job server, to use for fast motif comparison
        

    def _load_scores(self):
        self.scoredist = {}
        for metric in self.metrics:
            self.scoredist[metric] = {"total": {}, "subtotal": {}}
            for match in ["total", "subtotal"]:
                for combine in ["mean"]:
                    self.scoredist[metric]["%s_%s" % (match, combine)] = {}
                    score_file = os.path.join(self.config.get_score_dir(), "%s_%s_%s_score_dist.txt" % (match, metric, combine))
                    if os.path.exists(score_file):
                        for line in open(score_file):
                            l1, l2, m, sd = line.strip().split("\t")[:4]
                            self.scoredist[metric]["%s_%s" % (match, combine)].setdefault(int(l1), {})[int(l2)] = [float(m), float(sd)]
    
    def compare_motifs(self, m1, m2, match="total", metric="wic", combine="mean", pval=False):

        if metric == "fisim":
            return self.fisim(m1, m2)
        elif metric == "seqcor":
            return seqcor(m1, m2)
        elif match == "partial":
            if pval:
                return self.pvalue(m1, m2, "total", metric, combine, self.max_partial(m1.pwm, m2.pwm, metric, combine))
            elif metric in ["pcc", "ed", "distance", "wic", "chisq"]:
                return self.max_partial(m1.pwm, m2.pwm, metric, combine)
            else:
                return self.max_partial(m1.pfm, m2.pfm, metric, combine)

        elif match == "total":
            if pval:
                return self.pvalue(m1, m2, match, metric, combine, self.max_total(m1.pwm, m2.pwm, metric, combine))
            elif metric == "pcc":
                sys.stderr.write("Can't calculate PCC of columns with equal distribution!\n")
                return None
            elif metric in ["ed", "distance", "wic", "chisq"]:
                return self.max_total(m1.pwm, m2.pwm, metric, combine)
            else:
                return self.max_total(m1.pfm, m2.pfm, metric, combine)
                
        elif match == "subtotal":
            if metric in ["pcc", "ed", "distance", "wic", "chisq"]:
                return self.max_subtotal(m1.pwm, m2.pwm, metric, combine)
            else:
                return self.max_subtotal(m1.pfm, m2.pfm, metric, combine)
    

    def _check_length(self, l):
        # Set the length to a length represented in randomly generated JASPAR motifs 
        if l < 4:
            return 4
        if l == 13:
            return 14
        if l == 17:
            return 18
        if l == 19:
            return 20
        if l == 21:
            return 22
        if l > 22:
            return 30    
        return l    
    
    def pvalue(self, m1, m2, match, metric, combine, score):
        l1, l2 = len(m1.pwm), len(m2.pwm)
        
        l1 = self._check_length(l1)    
        l2 = self._check_length(l2)    
        
        m,s = self.scoredist[metric]["%s_%s" % (match, combine)][l1][l2]    
        
        try:
            [1 - norm.cdf(score[0], m, s), score[1], score[2]]
        except:
            print "HOEI: {0}".format(score)
            return 1
        return [1 - norm.cdf(score[0], m, s), score[1], score[2]]

    def fisim(self, m1, m2):
        try:
            m1 = m1.pfm
            m2 = m2.pfm
        except:
            pass
        
        fm1 = fisim.Motif.Motif(matrix=array(m1))
        fm2 = fisim.Motif.Motif(matrix=array(m2))
        score = fm1.fisim(fm2)
        if score[2]:
            return (score[0], score[1], 1)
        else:
            return (score[0], score[1], -1)

    def score_matrices(self, matrix1, matrix2, metric, combine):
        if metric in self.metrics and combine in self.combine:
            if metric == "fisim":
                s = self.fisim(matrix1, matrix2)[0]
            else:
               s = score(matrix1, matrix2, metric, combine)
            
            if s != s:
                return None
            else:
                return s
        
        else:
            if metric == "akl":
                func = akl
            elif metric == "chisq":
                func = chisq
            elif metric == "ssd":
                func = ssd
            else:
                try:
                    func = getattr(distance, metric)     
                except: 
                    raise Exception, "Unknown metric '{}'".format(metric)

            scores = []
            for pos1,pos2 in zip(matrix1,matrix2):
                scores.append(func(pos1, pos2))
            if combine == "mean":
                return mean(scores)
            elif combine == "sum":
                return sum(scores)
            else:
                raise "Unknown combine"

    def max_subtotal(self, matrix1, matrix2, metric, combine):
        scores = []
        min_overlap = 4 
        
        if len(matrix1) < min_overlap or len(matrix2) < min_overlap:
            return self.max_total(matrix1, matrix2, metric, combine)
    
        #return c_max_subtotal(matrix1, matrix2, metric, combine)

        for i in range(-(len(matrix2) - min_overlap), len(matrix1) - min_overlap + 1):
            p1,p2 = self.make_equal_length_truncate(matrix1, matrix2, i)
            s = self.score_matrices(p1, p2, metric, combine)
            if s:
                scores.append([s, i, 1])
    
        rev_matrix2 = [row[::-1] for row in matrix2[::-1]]
        for i in range(-(len(matrix2) - min_overlap), len(matrix1) - min_overlap + 1):
            p1,p2 = self.make_equal_length_truncate(matrix1, rev_matrix2, i)    
            s = self.score_matrices(p1, p2, metric, combine)
            if s:
                scores.append([s, i, -1])
        
        if not scores:
            return []
        return sorted(scores, key=lambda x: x[0])[-1]
    
    def max_partial(self, matrix1, matrix2, metric, combine):

        scores = []
    
        for i in range(-(len(matrix2) -1), len(matrix1)):
            p1,p2 = self.make_equal_length_truncate_second(matrix1, matrix2, i)    
            s = self.score_matrices(p1, p2, metric, combine)
            if s:
                scores.append([s, i, 1])
    
        rev_matrix2 = [row[::-1] for row in matrix2[::-1]]
        for i in range(-(len(matrix2) -1), len(matrix1)):
            p1,p2 = self.make_equal_length_truncate_second(matrix1, rev_matrix2, i)    
            s = self.score_matrices(p1, p2, metric, combine)
            if s:
                scores.append([s, i, -1])
        
        if not scores:
            return []
        return sorted(scores, key=lambda x: x[0])[-1]

    def max_total(self, matrix1, matrix2, metric, combine):
        scores = []
    
        for i in range(-(len(matrix2) -1), len(matrix1)):
            p1,p2 = self.make_equal_length(matrix1, matrix2, i)    
            s = self.score_matrices(p1, p2, metric, combine)
            if s:
                scores.append([s, i, 1])
    
        rev_matrix2 = [row[::-1] for row in matrix2[::-1]]
        for i in range(-(len(matrix2) -1), len(matrix1)):
            p1,p2 = self.make_equal_length(matrix1, rev_matrix2, i)    
            s = self.score_matrices(p1, p2, metric, combine)
            if s:
                scores.append([s, i, -1])
        
        if not scores:
            sys.stdout.write("No score {} {}".format(matrix1, matrix2))
            return []
        return sorted(scores, key=lambda x: x[0])[-1]
    
    def make_equal_length(self, pwm1, pwm2, pos, bg=None):
        if bg is None:
            bg = [0.25,0.25,0.25,0.25]
        
        p1 = pwm1[:]
        p2 = pwm2[:]
    
        if pos < 1:
            p1 = [bg for x in range(-pos)] + p1
        else:
            p2 = [bg for x in range(pos)] + p2
    
        diff = len(p1) - len(p2)
        if diff > 0:
            p2 += [bg for x in range(diff)]
        elif diff < 0:
            p1 += [bg for x in range(-diff)]
    
        return p1,p2
    
    def make_equal_length_truncate(self, pwm1, pwm2, pos):
        p1 = pwm1[:]
        p2 = pwm2[:]
    
        if pos < 0:
            p2 = p2[-pos:]
        elif pos > 0:
            p1 = p1[pos:]
        
        if len(p1) > len(p2):
            p1 = p1[:len(p2)]
        else:
            p2 = p2[:len(p1)]
        return p1, p2
    
    def make_equal_length_truncate_second(self, pwm1, pwm2, pos, bg=None):
        if bg is None:
            bg = [0.25,0.25,0.25,0.25]
        
        p1 = pwm1[:]
        p2 = pwm2[:]

        if pos < 0:
            p2 = p2[-pos:]
        else:
            p2 = [bg for x in range(pos)] + p2
            
        diff = len(p1) - len(p2)
        if diff > 0:
            p2 += [bg for x in range(diff)]
        elif diff < 0:
            p2 = p2[:len(p1)]
        return p1,p2

    def get_all_scores(self, motifs, dbmotifs, match, metric, combine, pval=False, parallel=True, trim=None):
    
            
        # trim motifs first, if specified
        if trim:
            for m in motifs:
                m.trim(trim)
            for m in dbmotifs:
                m.trim(trim)
        
        # hash of result scores
        scores ={}
        
        if parallel:    
            # Divide the job into big chunks, to keep parallel overhead to minimum
            # Number of chunks = number of processors available
            n_cpus = int(MotifConfig().get_default_params()["ncpus"])

            batch_len = len(dbmotifs) / n_cpus
            if batch_len <= 0:
                batch_len = 1
            jobs = []
            for i in range(0, len(dbmotifs), batch_len): 
                # submit jobs to the job server
                
                p = pool.apply_async(_get_all_scores, 
                    args=(self, motifs, dbmotifs[i: i + batch_len], match, metric, combine, pval))
                jobs.append(p)
            
            for job in jobs:
                # Get the job result
                result = job.get()
                # and update the result score
                for m1,v in result.items():
                    for m2, score in v.items():
                        if not scores.has_key(m1):
                            scores[m1] = {}

                        scores [m1][m2] = score
        
        else:
            # Do the whole thing at once if we don't want parallel
            scores = _get_all_scores(self, motifs, dbmotifs, match, metric, combine, pval)
        
        return scores

    def get_closest_match(self, motifs, dbmotifs, match, metric, combine, parallel=True):
        scores = self.get_all_scores(motifs, dbmotifs, match, metric, combine, parallel=parallel)
        for motif, matches in scores.items():
            scores[motif] = sorted(scores[motif].items(), cmp=lambda x,y: cmp(y[1][0], x[1][0]))[0]
        return scores

    def generate_score_dist(self, motifs, match, metric, combine):
        
        score_file = os.path.join(self.config.get_score_dir(), "%s_%s_%s_score_dist.txt" % (match, metric, combine))    
        f = open(score_file, "w")

        all_scores = {}
        for l in [len(motif) for motif in motifs]:
            all_scores[l] = {}

        sorted_motifs = {}
        for l in all_scores.keys():
            sorted_motifs[l] = [motif for motif in motifs if len(motif) == l]
        
        for l1 in all_scores.keys():
            for l2 in all_scores.keys():
                scores = self.get_all_scores(sorted_motifs[l1], sorted_motifs[l2], match, metric, combine)
                scores = scores.values()
                scores = [[y[0] for y in x.values() if y] for x in scores]
                scores = array(scores).ravel()
                f.write("%s\t%s\t%s\t%s\n" % (l1, l2, mean(scores), std(scores)))

        f.close()    

# import here is necessary as workaround
# see: http://stackoverflow.com/questions/18947876/using-python-multiprocessing-pool-in-the-terminal-and-in-code-modules-for-django
from gimmemotifs.mp import pool
