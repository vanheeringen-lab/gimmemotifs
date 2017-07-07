# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

""" 
Module to compare DNA sequence motifs (positional frequency matrices)
"""
from __future__ import print_function

# Python imports
import sys
import os
import random
import logging

# External imports
from scipy.stats import norm,entropy,chi2_contingency
from scipy.spatial import distance
import numpy as np

# GimmeMotifs imports
from gimmemotifs.config import MotifConfig
from gimmemotifs.c_metrics import pwmscan,score
from gimmemotifs.motif import parse_motifs
# pool import is at the bottom

try: 
    import copy_reg
    import types
    def _pickle_method(m):
        if m.im_self is None:
            return getattr, (m.im_class, m.im_func.func_name)
        else:
            return getattr, (m.im_self, m.im_func.func_name)

    copy_reg.pickle(types.MethodType, _pickle_method)
except:
    pass

# Create random sequence
nucs = []
L = 10 ** 4
for i in range(L):
    nucs.append(random.choice(['A', 'C', 'T', 'G']))
RANDOM_SEQ = "".join(nucs)

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
    return np.sum([(a-b)**2 for a,b in zip(x,y)] )

def seqcor(m1, m2, seq=None):
    l1 = len(m1)
    l2 = len(m2)

    l = max(l1, l2)

    if seq is None:
        seq = RANDOM_SEQ
    
    # Scan random sequence
    result1 = pwmscan(seq, m1.pwm, m1.pwm_min_score(), len(seq), False, True)
    result2 = pwmscan(seq, m2.pwm, m2.pwm_min_score(), len(seq), False, True)
    result1 = np.array(result1)
    result2 = np.array(result2)
    
    # Return maximum correlation
    c = []
    for i in range(l1):
        c.append(1 - distance.correlation(result1[:L-l-i],result2[i:L-l]))
    for i in range(l2):
        c.append(1 - distance.correlation(result1[i:L-l],result2[:L-l-i]))
    return max(c)


class MotifComparer(object):
    def __init__(self):
        self.config = MotifConfig()
        self.metrics = ["pcc", "ed", "distance", "wic"]
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
                        with open(score_file) as f:
                            for line in f:
                                l1, l2, m, sd = line.strip().split("\t")[:4]
                                self.scoredist[metric]["%s_%s" % (match, combine)].setdefault(int(l1), {})[int(l2)] = [float(m), float(sd)]
    
    def compare_motifs(self, m1, m2, match="total", metric="wic", combine="mean", pval=False):

        if metric == "seqcor":
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
        except Exception as e:
            print("Error with score: {}\n{}".format(score, e))
            return [1, np.nan, np.nan]
        return [1 - norm.cdf(score[0], m, s), score[1], score[2]]

    def score_matrices(self, matrix1, matrix2, metric, combine):
        if metric in self.metrics and combine in self.combine:
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
                    raise Exception("Unknown metric '{}'".format(metric))

            scores = []
            for pos1,pos2 in zip(matrix1,matrix2):
                scores.append(func(pos1, pos2))
            if combine == "mean":
                return np.mean(scores)
            elif combine == "sum":
                return np.sum(scores)
            else:
                raise ValueError("Unknown combine")

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
            p1 = [bg for _ in range(-pos)] + p1
        else:
            p2 = [bg for _ in range(pos)] + p2
    
        diff = len(p1) - len(p2)
        if diff > 0:
            p2 += [bg for _ in range(diff)]
        elif diff < 0:
            p1 += [bg for _ in range(-diff)]
    
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
            p2 = [bg for _ in range(pos)] + p2
            
        diff = len(p1) - len(p2)
        if diff > 0:
            p2 += [bg for _ in range(diff)]
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

            batch_len = len(dbmotifs) // n_cpus
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
                    for m2, s in v.items():
                        if m1 not in scores:
                            scores[m1] = {}
                        scores[m1][m2] = s
        
        else:
            # Do the whole thing at once if we don't want parallel
            scores = _get_all_scores(self, motifs, dbmotifs, match, metric, combine, pval)
        
        return scores

    def get_closest_match(self, motifs, dbmotifs=None, match="partial", metric="wic",combine="mean", parallel=True):
        """Return best match in database for motifs.

        Parameters
        ----------

        motifs : list or str
            Filename of motifs or list of motifs.

        dbmotifs ; list or str, optional
            Database motifs, default will be used if not specified.

        match : str, optional

        metric : str, optional

        combine : str, optional

        Returns
        -------

        closest_match : dict
        """

        if dbmotifs is None:
            pwm = self.config.get_default_params()["motif_db"]
            pwmdir = self.config.get_motif_dir()
            dbmotifs = os.path.join(pwmdir, pwm)
       
        motifs = parse_motifs(motifs)
        dbmotifs = parse_motifs(dbmotifs)

        dbmotif_lookup = dict([(m.id, m) for m in dbmotifs])

        scores = self.get_all_scores(motifs, dbmotifs, match, metric, combine, parallel=parallel)
        for motif in scores:
            scores[motif] = sorted(
                    scores[motif].items(), 
                    key=lambda x:x[1][0]
                    )[-1]
        
        for motif in motifs:
            dbmotif, score = scores[motif.id]
            pval, pos, orient = self.compare_motifs(
                    motif, dbmotif_lookup[dbmotif], match, metric, combine, True)
            
            scores[motif.id] = [dbmotif, (list(score) + [pval])]
        
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
                scores = [[y[0] for y in x.values() if y] for x in scores.values()]
                scores = np.array(scores).ravel()
                f.write("%s\t%s\t%s\t%s\n" % (l1, l2, np.mean(scores), np.std(scores)))

        f.close()    

# import here is necessary as workaround
# see: http://stackoverflow.com/questions/18947876/using-python-multiprocessing-pool-in-the-terminal-and-in-code-modules-for-django
try:
    from gimmemotifs.mp import pool
except:
    pass
