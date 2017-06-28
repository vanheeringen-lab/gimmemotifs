# Copyright (c) 2009-2017 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
""" Module to calculate motif scoring metrics.

Includes ROC AUC, MNCP, enrichment and others, which are calculated 
on the basis of motif scanning results.
"""

# External imports
from scipy.stats import stats, scoreatpercentile, kstest, fisher_exact
from sklearn.metrics import precision_recall_curve, roc_auc_score, roc_curve
import numpy as np

__all__ = [
    "recall_at_fdr",
    "fraction_fpr",
    "score_at_fpr",
    "enr_at_fpr",
    "max_enrichment",
    "phyper_at_fpr", 
    "mncp",
    "roc_auc",
    "roc_auc_xlim",
    "max_fmeasure",
    "ks_pvalue",
    "ks_significance",
]

def requires_scores(f):
    f.input_type = "score"
    return f

def requires_positions(f):
    f.input_type = "pos"
    return f

@requires_scores
def values_to_labels(fg_vals, bg_vals):
    """
    Convert two arrays of values to an array of labels and an array of scores.

    Parameters
    ----------
    fg_vals : array_like
        The list of values for the positive set.

    bg_vals : array_like
        The list of values for the negative set.

    Returns
    -------
    y_true : array
        Labels.
    y_score : array
        Values.
    """ 
    y_true = np.hstack((np.ones(len(fg_vals)), np.zeros(len(bg_vals))))
    y_score = np.hstack((fg_vals, bg_vals))
    
    return y_true, y_score

@requires_scores
def recall_at_fdr(fg_vals, bg_vals, fdr_cutoff=0.1):
    """
    Computes the recall at a specific FDR (default 10%).

    Parameters
    ----------
    fg_vals : array_like
        The list of values for the positive set.

    bg_vals : array_like
        The list of values for the negative set.
    
    fdr : float, optional
        The FDR (between 0.0 and 1.0).
    
    Returns
    -------
    recall : float
        The recall at the specified FDR.
    """
    if len(fg_vals) == 0:
        return 0.0
    
    y_true, y_score = values_to_labels(fg_vals, bg_vals)
    
    precision, recall, _ = precision_recall_curve(y_true, y_score)
    fdr = 1 - precision
    cutoff_index = next(i for i, x in enumerate(fdr) if x <= fdr_cutoff)
    return recall[cutoff_index]

@requires_scores
def matches_at_fpr(fg_vals, bg_vals, fpr=0.01):
    """
    Computes the hypergeometric p-value at a specific FPR (default 1%).

    Parameters
    ----------
    fg_vals : array_like
        The list of values for the positive set.

    bg_vals : array_like
        The list of values for the negative set.
    
    fpr : float, optional
        The FPR (between 0.0 and 1.0).
    
    Returns
    -------
    fraction : float
        The fraction positives at the specified FPR.
    """
    fg_vals = np.array(fg_vals)
    s = scoreatpercentile(bg_vals, 100 - fpr * 100)
    
    return [sum(fg_vals >= s), sum(bg_vals >= s)]

@requires_scores
def phyper_at_fpr(fg_vals, bg_vals, fpr=0.01):
    """
    Computes the hypergeometric p-value at a specific FPR (default 1%).

    Parameters
    ----------
    fg_vals : array_like
        The list of values for the positive set.

    bg_vals : array_like
        The list of values for the negative set.
    
    fpr : float, optional
        The FPR (between 0.0 and 1.0).
    
    Returns
    -------
    fraction : float
        The fraction positives at the specified FPR.
    """
    fg_vals = np.array(fg_vals)
    s = scoreatpercentile(bg_vals, 100 - fpr * 100)
    
    table = [
            [sum(fg_vals >= s), sum(bg_vals >= s)],
            [sum(fg_vals < s), sum(bg_vals < s)],
            ]
    
    return fisher_exact(table, alternative="greater")[1]

@requires_scores
def fraction_fpr(fg_vals, bg_vals, fpr=0.01):
    """
    Computes the fraction positives at a specific FPR (default 1%).

    Parameters
    ----------
    fg_vals : array_like
        The list of values for the positive set.

    bg_vals : array_like
        The list of values for the negative set.
    
    fpr : float, optional
        The FPR (between 0.0 and 1.0).
    
    Returns
    -------
    fraction : float
        The fraction positives at the specified FPR.
    """
    fg_vals = np.array(fg_vals)
    s = scoreatpercentile(bg_vals, 100 - 100 * fpr)
    return len(fg_vals[fg_vals >= s]) / float(len(fg_vals))

@requires_scores
def score_at_fpr(fg_vals, bg_vals, fpr=0.01):
    """
    Returns the motif score at a specific FPR (default 1%).

    Parameters
    ----------
    fg_vals : array_like
        The list of values for the positive set.

    bg_vals : array_like
        The list of values for the negative set.
    
    fpr : float, optional
        The FPR (between 0.0 and 1.0).
    
    Returns
    -------
    score : float
        The motif score at the specified FPR.
    """
    bg_vals = np.array(bg_vals)
    return scoreatpercentile(bg_vals, 100 - 100 * fpr)

@requires_scores
def enr_at_fpr(fg_vals, bg_vals, fpr=0.01):
    """
    Computes the enrichment at a specific FPR (default 1%).

    Parameters
    ----------
    fg_vals : array_like
        The list of values for the positive set.

    bg_vals : array_like
        The list of values for the negative set.
    
    fpr : float, optional
        The FPR (between 0.0 and 1.0).
    
    Returns
    -------
    enrichment : float
        The enrichment at the specified FPR.
    """
    pos = np.array(fg_vals)
    neg = np.array(bg_vals)
    s = scoreatpercentile(neg, 100 - fpr * 100)
    neg_matches = float(len(neg[neg >= s]))
    if neg_matches == 0:
        return float("inf")
    return len(pos[pos >= s]) / neg_matches * len(neg) / float(len(pos))

@requires_scores
def max_enrichment(fg_vals, bg_vals, minbg=2):
    """
    Computes the maximum enrichment.

    Parameters
    ----------
    fg_vals : array_like
        The list of values for the positive set.

    bg_vals : array_like
        The list of values for the negative set.
    
    minbg : int, optional
        Minimum number of matches in background. The default is 2.
    
    Returns
    -------
    enrichment : float
        Maximum enrichment.
    """
    scores = np.hstack((fg_vals, bg_vals))
    idx = np.argsort(scores)
    x = np.hstack((np.ones(len(fg_vals)), np.zeros(len(bg_vals))))
    xsort = x[idx]
    l_fg = len(fg_vals)
    l_bg = len(bg_vals)
    m = 0
    s = 0
    for i in range(len(scores), 0, -1):
        bgcount = float(len(xsort[i:][xsort[i:] == 0])) 
        if bgcount >= minbg:
            enr = (len(xsort[i:][xsort[i:] == 1]) / l_fg) / (bgcount / l_bg)
            if enr > m:
                m = enr
                s = scores[idx[i]]
    return m

@requires_scores
def mncp(fg_vals, bg_vals):
    """
    Computes the Mean Normalized Conditional Probability (MNCP).

    MNCP is described in Clarke & Granek, Bioinformatics, 2003.

    Parameters
    ----------
    fg_vals : array_like
        The list of values for the positive set.

    bg_vals : array_like
        The list of values for the negative set.
    
    Returns
    -------
    score : float
        MNCP score
    """
    fg_len = len(fg_vals)
    total_len = len(fg_vals) + len(bg_vals)

    if not isinstance(fg_vals, np.ndarray):
        fg_vals = np.array(fg_vals)
    if not isinstance(bg_vals, np.ndarray):
        bg_vals = np.array(bg_vals)
    
    fg_rank = stats.rankdata(fg_vals)
    total_rank = stats.rankdata(np.hstack((fg_vals, bg_vals)))

    slopes = []
    for i in range(len(fg_vals)):
        slope = ((fg_len - fg_rank[i] + 1) / fg_len ) / (
                (total_len - total_rank[i] + 1)/ total_len)
        slopes.append(slope)
    
    return np.mean(slopes)

@requires_scores
def roc_auc(fg_vals, bg_vals):
    """
    Computes the ROC Area Under Curve (ROC AUC)

    Parameters
    ----------
    fg_vals : array_like
        list of values for positive set

    bg_vals : array_like
        list of values for negative set
    
    Returns
    -------
    score : float
        ROC AUC score
    """
    # Create y_labels
    y_true, y_score = values_to_labels(fg_vals, bg_vals)
    
    return roc_auc_score(y_true, y_score)

@requires_scores
def roc_auc_xlim(x_bla, y_bla, xlim=0.1):
    """
    Computes the ROC Area Under Curve until a certain FPR value.

    Parameters
    ----------
    fg_vals : array_like
        list of values for positive set

    bg_vals : array_like
        list of values for negative set

    xlim : float, optional
        FPR value
    
    Returns
    -------
    score : float
        ROC AUC score
    """
    x = x_bla[:]
    y = y_bla[:]

    x.sort()
    y.sort()

    u = {}
    for i in x + y:
        u[i] = 1

    vals = sorted(u.keys())
    
    len_x = float(len(x))
    len_y = float(len(y))
    
    new_x = []
    new_y = []
    
    x_p = 0
    y_p = 0
    for val in vals[::-1]:
        while len(x) > 0 and x[-1] >= val:
            x.pop()
            x_p += 1
        while len(y) > 0 and y[-1] >= val:
            y.pop()
            y_p += 1
        new_y.append((len_x - x_p) / len_x)
        new_x.append((len_y - y_p) / len_y)
    
    #print new_x
    #print new_y
    new_x = 1 - np.array(new_x)
    new_y = 1 - np.array(new_y)
    #plot(new_x, new_y)
    #show()

    x = new_x
    y = new_y

    if len(x) != len(y):
        raise ValueError("Unequal!")

    if not xlim:
        xlim = 1.0

    auc = 0.0
    bla = zip(stats.rankdata(x), range(len(x)))
    bla = sorted(bla, key=lambda x: x[1])
    
    prev_x = x[bla[0][1]]
    prev_y = y[bla[0][1]]
    index = 1

    while index < len(bla) and x[bla[index][1]] <= xlim:

        _, i = bla[index]
        
        auc += y[i] * (x[i] - prev_x) - ((x[i] - prev_x) * (y[i] - prev_y) / 2.0)
        prev_x = x[i]
        prev_y = y[i]
        index += 1
    
    if index < len(bla):
        (rank, i) = bla[index]
        auc += prev_y * (xlim - prev_x) + ((y[i] - prev_y)/(x[i] - prev_x) * (xlim -prev_x) * (xlim - prev_x)/2)
 
    return auc

@requires_scores
def roc_values(fg_vals, bg_vals):
    """
    Return fpr (x) and tpr (y) of the ROC curve.

    Parameters
    ----------
    fg_vals : array_like
        The list of values for the positive set.

    bg_vals : array_like
        The list of values for the negative set.

    Returns
    -------
    fpr : array
        False positive rate.
    tpr : array
        True positive rate.
    """
    if len(fg_vals) == 0:
        return 0
    
    y_true, y_score = values_to_labels(fg_vals, bg_vals)
    
    fpr, tpr, _thresholds = roc_curve(y_true, y_score)
     
    return fpr, tpr

@requires_scores
def max_fmeasure(fg_vals, bg_vals):
    """
    Computes the maximum F-measure.

    Parameters
    ----------
    fg_vals : array_like
        The list of values for the positive set.

    bg_vals : array_like
        The list of values for the negative set.
    
    Returns
    -------
    f : float
        Maximum f-measure.
    """
    x, y = roc_values(fg_vals, bg_vals)
    x, y = x[1:], y[1:] # don't include origin
    
    p = y / (y + x)
    filt = np.logical_and((p * y) > 0, (p + y) > 0)
    p = p[filt]
    y = y[filt]
    
    f = (2 * p * y) / (p + y)
    if len(f) > 0:
        #return np.nanmax(f), np.nanmax(y[f == np.nanmax(f)])
        return np.nanmax(f)
    else:
        return None

@requires_positions
def ks_pvalue(fg_pos, bg_pos=None):
    """
    Computes the Kolmogorov-Smirnov p-value of position distribution.

    Parameters
    ----------
    fg_pos : array_like
        The list of values for the positive set.

    bg_pos : array_like, optional
        The list of values for the negative set.
    
    Returns
    -------
    p : float
        KS p-value.
    """
    if len(fg_pos) == 0:
        return 1.0
    a = np.array(fg_pos, dtype="float") / max(fg_pos)
    p = kstest(a, "uniform")[1]
    return p

@requires_positions
def ks_significance(fg_pos, bg_pos=None):
    """
    Computes the -log10 of Kolmogorov-Smirnov p-value of position distribution.

    Parameters
    ----------
    fg_pos : array_like
        The list of values for the positive set.

    bg_pos : array_like, optional
        The list of values for the negative set.
    
    Returns
    -------
    p : float
        -log10(KS p-value).
    """
    p = ks_pvalue(fg_pos, max(fg_pos))
    return -np.log10(p)
    

 
