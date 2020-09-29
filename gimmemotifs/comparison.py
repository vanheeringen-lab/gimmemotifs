# Copyright (c) 2009-2019 Simon van Heeringen <simon.vanheeringen@gmail.com>
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
import logging

# External imports
from scipy.stats import norm, entropy, chi2_contingency
from scipy.spatial import distance
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.feature_selection import RFE
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.metrics import average_precision_score, roc_auc_score
import pandas as pd


# GimmeMotifs imports
from gimmemotifs.config import MotifConfig
from gimmemotifs.c_metrics import pfmscan, score
from gimmemotifs.motif import parse_motifs, read_motifs
from gimmemotifs.utils import pfmfile_location

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
except Exception:
    pass

RCDB = (
    "AATGCCTGCCTCGCCCATATAAGCATCAAGGCATATTTATTACCGGCCGCATGGAACCTTTCCTCCCAATCGAACAAGTC"
    "AGGTGGCATCCGGGAACGACTGTTGACACCGACAAGAAGCCGACCCCGTAATGGGGAACCGAGTTTGGCTTAGCTGCTTC"
    "ATATGGCGTTCACTAGTAATATTCGGAGACAACTACGTCCGGTATCTAGGCGATGGGGTGAGCGTGGTATAAGTCCATTC"
    "AACCGTTATGACGGCTAAGGTCGCTACTCTCGTGAGCTATCTGTAGGACGAAGGGTTAGAGTCCCATGCAGAGAGGGGAA"
    "TCACTGAGAGCCACAATTCAGAGATATCACGCCGCGTACTGCAGCCTTTTCGGCTCCTCGGACTCGCTAATACGCGTCAG"
    "CTGGGAGCCCGCATATCTGAAAACCGCTAACTTCTCTCGAGCCCCTCAAGTACACTGGGAAAGTCCCCAATATGACCTGT"
    "CTCACCGATCTAATCCTACGCGATCCTATAGAGTATAGTGCTCCTTCATGCTGACTCCACAGTCACTCGCCGTCTCATGA"
    "CCGATGTTGGGGTACAAATCTCGTTCAATTCTAAACTTTTTAGTACCCATAGCATAGCGGGAATGGCCTTGTGTTCCGAC"
    "GCGAGGTACCGCCCCAAGTCGACCAAGTAGAATAAGTTAATGAGCCCTGTCATTAGACTATATTCCTGAAAGCATAAGCG"
    "AGCGATCAGCTATACATAATTAGTAGGTGATGTACAAGTTACGGTCTGACATGTTCAGCAATAGCAATGTTTGTGTGCAG"
    "AACTCACAGCTCTCGCGGCTGCTGCATCCACCTCTGGGTAGGAGGACACACGCCCGAAGACTTCAAATGGGCTGACCTCT"
    "CATAGGTAGTGCAGGGACTATGGTTACGAGTACTGTTTTTTCCTTCGGTAAAATAATACTAGCAACAGAAAATAGCCGGA"
    "TTGGTAGTTATCTTCGCACAATCCATCTTGCTGGACGGTTCGCGAGCTACTTCCGTTAATTGACTCTGGCCCGTAGTCCT"
    "CGATAGGCGTACGATAGACAAGCACAAATTTTTAATCACGGGTTGCGTCACCGGCACTCGGATGTAAAGGATGGTACGTT"
    "GGCGGCTTACTGATTCGCCAGTCTCAAGCATTTGGGGGTAATTGGTCACGTCTCGGAGGGCGACGACAAACGGTCCCCTT"
    "AAATCTATAATGTGATTGGCTGTGTTGGTGGCCAAAGTTCAACGGCTTGTGAGAGAGATGGGAATACAATCAGGTCGAAG"
    "CAACGAATACGTTACCTGCGAATGCTCAAGACAGCCGTAAGACGAACAGCTTTCGGTGCCCCACTAGCTGTAGCGGATGA"
    "GTGATCCTCGTTATAGCAGGCTCAAAACACACTACTATCTCGATTCGGTCCGGCGCGTTCTCCTGCGGGGGGCGATACTT"
    "AAGACCATAAAACATCTTCTTCCCCGATGCTAAACGCGAAATTACGAAGCGTTGGTCCGAAGCCTAAATTTCCTAGCACT"
    "TAGGAAGTGGCGGTTGTGAAGGCGAGCCAACTAGGGAGAATAGAACCCGCTCGGCGTTTGCTGCCAGAACGGCAGACCAA"
    "ATCGTATAGGCTCTTACCGCTTCCTCACACATTCGGCAGCCCCCAGGTTACACTTACGTACGGGATACCTCCATGGTGTG"
    "TGGGACGAGGTTTCTCTATGCCCTGATATTGACAGTCCTAAAGTGCGGGACTGATCTTGGTGACAGCGAGGGACCAGGCC"
    "AGCCAGCATATGTCTAACGCCCCGCCGACGATAACATTTTCGTGGTGAGTAGACCCTTTCGCCGCTAGCCCACGAGTGGT"
    "GCAAATAAAAAAGAAGTCCTTTCAACACTGAACTGGAGTTACATGGGTGCGGCTCACTATTTAATGGAATTGGCAGATTT"
    "CAAAATGGATAACGTCCTTATGTTTAGGTTTGGTATGAAACCCAGCGGCCCCCGAGCTTTTTCACAGTATAATCCCACTT"
    "AAAAGTGGATTACTTAGATGCGCCCGGGTCGCCAATGACGTGCCCTCGACACAATAACACGTTGAGGATTGCGAAACTTC"
    "CAGCTACCTGAGTCGGAGCGCGTAATCTCTAATAAGGCACTAATGTCGTATTTGGCGCCCAATGTCATCAACCCACTGCA"
    "CGCCTCGAGTCATACGTCTTCCAATCTTGACGTACATTGTCAGCGCGCAGTTCGTTCGACTGGTCGGATAGAGCACGGAA"
    "ATCGCACGGGGATCATGAATACTGAATGATGCGAACGGTATTACGTGGCAGGAAACGCACCATTGTACGGCAAACTTGAC"
    "TTTTGTTCTTCATTACCAGACGCACAGCAACTGGGGCCTCGTAGTGACGGGTATGCGCGGTGGGCATGCGTAGTAACTGC"
    "TGGGCCACGACCTCATAAGTAGTTCTAGGGTCAGAATGGTTGAGTCCACATCCATACCTTCTTAGGTCTGTTAGCTCTTC"
    "ACGTTTGAATAGGTTAAGCTGGTCTCCTATCTTACACCGTACTTACCACTTGTGCCATAATCATTTATATTATGCCAAAA"
    "CTTACATAGTATACCGACCGCCTACGTGAGTGGCTGAGAACCTACGGTTAACGCTATTCGCGGGGTTTGCGATGCACCGT"
    "CAAACCACTATCAAAGATTTTATAATAGCTACAAAACGGACTAGTTTTCGCTAGTCTCGACGAATGAGGCCAATAATTTG"
    "CGTGGCGATCTCTTTGTCCTTCCTGGTGTTGCACAGAACCATCTACGGCCGATAATCTTCCTTTGAACCGCGTCTTGCAC"
    "TGCTTAAGGATTTATCCCATTGATCTATTGTCGATTACCTCTTCCGAGCCGGCAGGTTTACTTCTAACTCTGTGTGTATT"
    "TACAGTTCCGTCTAGATCTGAGACCCATCACTAAGGCGTTGTCGCCGGCCCAACGGGTGCCAACAGCGCAGCTGTCCAAA"
    "CAAACCGGTAGCAAGCGATGTCCTATGTGGTTATTATCAACGACACCCTGTATAAACTCAGCTCCTAACAAGCTCCGGCA"
    "TTCTCATCTATGGCCCTCTGTCTAGGATACGCTAAGACAATCGTCAGATGGCGGGTAAACTGGCACTTGCCTTCCATAGG"
    "GCCCCTAGGACCCGCAAGGGCTGTATGAGCGCTATCATTGACGACTCTCAATAGTAAAAACGTAACGCGGATCGTTTCAA"
    "TGCGATCGAGTTCATTTGTCTGTAACGACCCTCTTTACGACGAGTTGTTGTACTGAGGAGTCTGCAGGCCGTTACTTGTC"
    "CCCGGCCTTACCCAATTGCTCACATCATGGCGCTTTAAATGGCTCCGTGAAATTGTAACACATGACTTACTCTTATGGAC"
    "CATTTCGGAATCTAGCGCCCCTGCGTGTACAGAATCGGATCAAGAAAGCGACGGCATGATTGCACGTGGTTCATCTGCGG"
    "TATGTACCTAGTGCGAGCAAGGTGGGAACAGGGAACTAAATAACGCATTCACACTTTGTATGTGACAAAATACCCCTTTT"
    "AGCTAGGGCTTTTACACGGTAGGCAAAGACGGGAGGTCCCAACCTGTTGTGCATTTCAGTCTTGGCAACATGCAACCGAC"
    "GTGGGGTTGGAATCGTGAACTATAAGGTACTCCTTGACAACATTAATTTTGGGAAGTCACCTCCGTTGCCTGTAGAAACG"
    "TTGCGCATGTAGACTGTATTAAGTTGAGAAGCGCTGCAATGGCATTGTTTTCCGATGGCCGCCGCAAAGTATTTCTTAAT"
    "CCGCCAGGTGAAGTGACCTATAAAGCTTCAGCTTCTGCCCAGTCCGCGCGCCCTTATCGATTTAGAGGGATAGGGACGGA"
    "CAGAAGTGTCTTTCCATCCAGACCCCTGTGTAGTCGTGCGCGAAGCTCTGTTTGATTACAGGAACAACACGGCGTCTCTC"
    "CGTCAGTGGTCGTCTGTGGTCTATGATCAGGAGTGGAGACTTTAAGCGGGGAGTCGAGTATCTTTTATCACAGGAGGCAA"
    "TTTACCTTGAGTATGTTAGATTGCTTTGACTGAGTACGTGTGACCGCGACGCTTACGCTGCGGATACTGCGGCACCTCGG"
    "CTAGAGCGCATTATTAAAGATCTCGCAGAGTGGGCGGAACCCAAAGCGTACCGTCGACATGCTTCCGGGTGGAACTCCTA"
    "GTCACAGATGTAGCACAGTGGAATAGTTTCACTTTCGTTTATTCCGCGATGATTCACTCCGATCATTATGATGGTGGAGG"
    "GAGTGTCACGGCACATCGACTAGATGTGTACGCAAGAATCCAGCGTTCGGATTCAATACATGTGAACACCTGCAACGCAG"
    "CCATCAGGCTGTTATAAGACTAAGATGCCTCAAAGTCGTACAGGTATATCATATAGAATCTGGGGAGGGGTGGTACCCGG"
    "TTCTTATTGATTTCGATGCGGTTAGTTTGAGTGTTTGCATTGCCCCGACGGGCAGGGCCGACTTAGTAAGAATATCCACA"
    "AGACCGGCTGAAACAGCATTGGTGCCTAGTACTTTCTAATTTATGAATCTTATCCGTGGGTGTGCGCTTAGACCTGGTTT"
    "AATTACACAACTTCAGTTCTCGCTCAGGCCCTTTAGATATTTGACCGTGAGAAAGAGCGGCTACACTATGAACGGATGGG"
    "CCTAATCAGACGGTCGTAGGGATCTACTGCCTGAAGATCCTTGTCATGCATTCCAGAGCGTAGACACATAACTGGTTGTC"
    "CATGTCCAGGCGTCCGTAGCGCAAAACCTTCAGAAATGAGATCCGTCGTCCACGATCAACATATACGTGCGGTCTAGCTT"
    "TGCAACTAAGCACTGTAGTGGACACTAGAAAACGAGCCTTGCCCAAGAGGCACGGCTCTGATCCATGAAATAGGAATCCT"
    "GCATTAGCAAAAAGGAGATGCATAACCCCTATATCCCCTAATTCCCTCATGGTTTTAAACTACATTACTGCTCTTGTAGA"
    "GACATCTCAAATTCTCTTAGAGCTCCCTGAGATTGTTGCGGAGCAGACGTCAGGACGGCCCGCCTTGGTTTCGAGCATCT"
    "GGTGGTCCTTGGCTATATACAACGTCACAATGAATGTATGGCTAGGCTACTAGACTCTTGAACAGATCCACGCCACATGG"
    "CTTTCCGGCCAGTAAAGTAGGGCGTATAACCGAACCCTATTAGTTACCGTATTGTTCATGTACTCTGACTAAAGCAGAGC"
    "AACCACCGAGATGACTACTGTACCGGGCTTATTTCGTAGCTTAGGCGCACTGATAGGTGGTGGGTCCGTGTCGGGAGTTT"
    "TTGACATAAACAACGATACCATGCCAGCGCCTGCGCCGTTTACGTCGCAGTGAAGCGAATTGCAAGCCCAAACGACGCAG"
    "TCCAACGAGTCGCTGCTAACGATGCCGAGTCTCTGTATCTCACTTATCTGCAAGTATGGTGATAGATAGAAGACCCGAAT"
    "TTTAGACGCTCGATGACCCAGTATGCTCTCCAGGTACGGATTTGCTTACCTACACGATGTATCATGTCTGCGTAAAACTG"
    "ACAGATTGATATGGAGCCTCTGCCTTTGGTGTCTGAACGTCTACCTTTGCTCCAAATTGATGTTCCAAAGAAAAGCACCC"
    "TTGTAAACAGGCCTCCAACAAGATATAATTGCGGTGCTGAGCCAGGGCGGCGCCTCCCTTGGGATCGAAGTTCGGCGCTC"
    "GTCACGACAGCAGGGTCCAGTCATTTTAACAGTGACATCAGCCTCCGCTCAACAAATAGAGGCGAAGTCGCGTGAGGCTA"
    "GTTAGTATCAGCAGATCGTAGAAGGGGGAGATAATGGCGAGACACGTGATGCAATACTCCGGGGGTGCACTCTTCTGGTT"
    "AAATTCCTTGCATAGGCATGTGCGGAAGGGACACGAGCACCAGTGGCCTATTGGATCCCGGTGTGGTAGAACAGTTAGGC"
    "AGCTCATGTTACTCATTGCGCCACCCAACTGCAAAAGAGTTTCCCGAACATTCTTTAGCCAGATCATCATCTTTCACCAG"
    "CGATTCAGCGAACACGAACGTAGCATGCACTTCGTTAAACGTCGTGTGGCAATATTTTTCGAACTAGAGGAATGTCTTCA"
    "ACTGAGCAGGACTCCGCGGGATTTTCATCGTAAGGTGCGCAGAAACTGTCTGGCGTCGATACAAGCGTGAATATGGTAAA"
    "GCCGGGACGCGCGAGTTAAAGGCGCGGGCGAACTCTCCCGGGGCTACCAGGGGAGCATGGCAGTCTGTCGAATTCATCCT"
    "AATATGTAGGGGGGGATTCGTGATCTCCGCCGGTTACCCTTAGAACGATTAATATCGCTATAGTACGAGCTCAATCAAGT"
    "CTGGGACATTTGCCCGGACTTGCTTGGCCGTGGAGCGGTCCACCGGAGGATGTCGCTCTCTTCGAATGGGTTACTGTCCC"
    "GTGTGCCCATGGCCAGGCTTGGAGTAGGCTGGAGCTAGCAGAATTACTCCCATAACAAAGCCTTCGAGGTCTTTATTGCG"
    "TTCCTCGCATGAGTTATTCGATATATGCAGTCGCATAGTTGATTCTCGGGCACTGGATCGGCTGGGGTCCCGAGTAAATC"
    "AGCGGTGATTTTTGGTTGGTTCGTATGGGGGAAATGTAATGCCCCTTGAAAAGACTCATATCCTGGGCTCTACGTATCTG"
    "GAGGAGGGTAGCGTGTTTAAGATCGATGGACGTGTATACGGGTCACCCGAGGGGGTCAAGCTTATCATCCCTGGGATGTT"
    "AACTGATGCCCGCTGCCCCCTGAACCCCCCACCCTACCTCACTGGCTTGAAGCCCCATGACGCTAGGTTGACGCCCAGAA"
    "AGGCCTACAATATACTTCATCACACCTCAATTGTCCGGAATGAAAGTTATATATTTCCACGTTCTTGGATTTCTATGTCC"
    "GCTGTGGCCGAAGTACCGAGCAGCACTCAGACTACCCCGAAACGGGGCCAACCGGATCCTGAGGTTCACCGTGGCTACGC"
    "CGGGCCTTCTCCATCAAAAGGTGAGGGCTAACCAATACCAGCCCGTGGTCATTCTAGCATTAAGCCTCGGTCGCACCGGT"
    "CCTGGCAAAATCGGCACGAGGGTGAATTTCTCGTCTATTAACACTTGGCGTAACAGGATGCTGCGATAAACCCGATTAGT"
    "CGGGTAGAGTGAGTCACGCATAAAGTTGCCCTAAGGAACTTAGCACGCACGCGGGTGTACCAGTACCAAGCTAATCTGAT"
    "AACCTCCTTACAGCACATTGAGTTGCGACAGTGTGCTCATCGGCGATTATTCAGGTTCGAGTGTGGAAGCCACCGTTGTG"
    "GGTATCAAGCGCGACACTGCCGCGGTAGATAAAGGTCTCGTATCGGGCCCACAGGCGACTATCCTTCTGATGAGGACGCT"
    "GTACTACTCGCAAGTTCTTTCGACGTCCCCCGTCCGCCTGACCATGTTTCTTCGTGCATACACCTTCCCTAAAACGCTCC"
    "TGGAAACACCGCCGTATGTCGGCAACCTTGGAAAATTGGACGCATCCCGATCGTCTCCCACGCGTTGACTATTACAACCC"
    "TCAGCGTGCCACTAAACCGATATCCGGAGAAAAACCCTGGATACACTCTCTAGGTGCAGCTTGCGCGTGTCAGGGTATAG"
    "ACTTAACAATTAGCCTTAACTATCGGACACCTAACGTTCCAGTAGGAACCAGCACGTTAAGAGGAGAGGCCGAGAGGATA"
    "ATTCGTACCCTAACCTTAATTCAAGTGATACGGAATTCTTCTCGAATCTCCAATGGAGTGCGTATGCCGCCTCAGAACAC"
    "TAACACCCCGGTCAGTAGCAGTCAGTCGGTCTCTTGCCGCTTGGGTCGAGCTGGCGCGATTGCCGGGTTAATCTACACCA"
    "CAAAGGGAAAAGTAAGGGATTGAAGGGCAGTTTGTAAGCTCAGAGGACCTACCATTAGGTATTCTGGACTAACGGGCGGG"
    "AGCGAAGGTCCGCATTTAGCGCGGCCTAGGTCGGGCTATCCAAAAGCGGAAACCGTGCACATATCGTGCTTCTAGAGTTG"
    "GTAACAACCTCTACTACAACTGTGTCATAGCCTGTCCTGTCGGACGATGGTCGCGATAGCACCTGTGCCTGGGTGATCAA"
    "ATCATACTACCGATTGACCTTGTTTCGTCATCCAACTTGGACCTGCTTGCCATCCCCACATTTACGGAGTTGACCAGAGA"
    "ACATAGCTTCCCGTTCCGGTCGATCACAAAAA"
)


logger = logging.getLogger("gimme.comparison")


# Function that can be parallelized
def _get_all_scores(mc, motifs, dbmotifs, match, metric, combine, pval):
    try:
        scores = {}
        for m1 in motifs:
            scores[m1.id] = {}
            for m2 in dbmotifs:
                scores[m1.id][m2.id] = mc.compare_motifs(
                    m1, m2, match, metric, combine, pval=pval
                )
        return scores
    except Exception:
        logging.exception("_get_all_scores failed")


def akl(p1, p2):
    """Calculates motif position similarity based on average Kullback-Leibler similarity.

    See Mahony, 2007.

    Parameters
    ----------
    p1 : list
        Motif position 1.

    p2 : list
        Motif position 2.

    Returns
    -------
    score : float
    """
    return 10 - (entropy(p1, p2) + entropy(p2, p1)) / 2.0


def chisq(p1, p2):
    """Calculates motif position similarity based on chi-square.

    Parameters
    ----------
    p1 : list
        Motif position 1.

    p2 : list
        Motif position 2.

    Returns
    -------
    score : float
    """
    return chi2_contingency([p1, p2])[1]


def ssd(p1, p2):
    """Calculates motif position similarity based on sum of squared distances.

    Parameters
    ----------
    p1 : list
        Motif position 1.

    p2 : list
        Motif position 2.

    Returns
    -------
    score : float
    """
    return 2 - np.sum([(a - b) ** 2 for a, b in zip(p1, p2)])


def seqcor(m1, m2, seq=None):
    """Calculates motif similarity based on Pearson correlation of scores.

    Based on Kielbasa (2015) and Grau (2015).
    Scores are calculated based on scanning a de Bruijn sequence of 7-mers.
    This sequence is taken from ShortCAKE (Orenstein & Shamir, 2015).
    Optionally another sequence can be given as an argument.

    Parameters
    ----------
    m1 : Motif instance
        Motif 1 to compare.

    m2 : Motif instance
        Motif 2 to compare.

    seq : str, optional
        Sequence to use for scanning instead of k=7 de Bruijn sequence.

    Returns
    -------
    score, position, strand
    """
    l1 = len(m1)
    l2 = len(m2)

    length = max(l1, l2)

    if seq is None:
        seq = RCDB

    L = len(seq)

    # Scan RC de Bruijn sequence
    result1 = pfmscan(seq, m1.pwm, m1.pwm_min_score(), len(seq), False, True)
    result2 = pfmscan(seq, m2.pwm, m2.pwm_min_score(), len(seq), False, True)

    # Reverse complement of motif 2
    result3 = pfmscan(seq, m2.rc().pwm, m2.rc().pwm_min_score(), len(seq), False, True)

    result1 = np.array(result1)
    result2 = np.array(result2)
    result3 = np.array(result3)

    # Return maximum correlation
    c = []
    for i in range(l1 - l1 // 3):
        c.append(
            [
                1
                - distance.correlation(
                    result1[: L - length - i], result2[i : L - length]
                ),
                i,
                1,
            ]
        )
        c.append(
            [
                1
                - distance.correlation(
                    result1[: L - length - i], result3[i : L - length]
                ),
                i,
                -1,
            ]
        )
    for i in range(l2 - l2 // 3):
        c.append(
            [
                1
                - distance.correlation(
                    result1[i : L - length], result2[: L - length - i]
                ),
                -i,
                1,
            ]
        )
        c.append(
            [
                1
                - distance.correlation(
                    result1[i : L - length], result3[: L - length - i]
                ),
                -i,
                -1,
            ]
        )

    return sorted(c, key=lambda x: x[0])[-1]


class MotifComparer(object):
    """Class for motif comparison.

    Compare two or more motifs using a variety of metrics. Probably the best
    metric to compare motifs is seqcor. The implementation of this metric
    is similar to the one used in Grau (2015), where motifs are scored
    according to the Pearson correlation of the scores along sequence. In this
    case a de Bruijn of k=7 is used.

    Valid metrics are:
    seqcor - Pearson correlation of motif scores along sequence.
    pcc - Pearson correlation coefficient of motif PFMs.
    ed - Euclidean distance-based similarity of motif PFMs.
    distance - Distance-based similarity of motif PFMs.
    wic - Weighted Information Content, see van Heeringen 2011.
    chisq - Chi-squared similarity of motif PFMs.
    akl - Similarity based on average Kullback-Leibler similarity, see Mahony, 2011.
    ssd - Sum of squared distances of motif PFMs.

    Examples
    --------
    mc = MotifComparer()

    # Compare two motifs
    score, pos, strand = mc.compare_motifs(m1, m2, metric="seqcor")

    # Compare a list of motifs to another list of motifs
    mc.get_all_scores(motifs, dbmotifs, match, metric, combine)

    # Get the best match for every motif in a list of reference motifs
    get_closest_match(motifs, dbmotifs=None)
    """

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
                    score_file = os.path.join(
                        self.config.get_score_dir(),
                        "%s_%s_%s_score_dist.txt" % (match, metric, combine),
                    )
                    if os.path.exists(score_file):
                        with open(score_file) as f:
                            for line in f:
                                l1, l2, m, sd = line.strip().split("\t")[:4]
                                self.scoredist[metric][
                                    "%s_%s" % (match, combine)
                                ].setdefault(int(l1), {})[int(l2)] = [
                                    float(m),
                                    float(sd),
                                ]

    def compare_motifs(
        self, m1, m2, match="total", metric="wic", combine="mean", pval=False
    ):
        """Compare two motifs.

        The similarity metric can be any of seqcor, pcc, ed, distance, wic,
        chisq, akl or ssd. If match is 'total' the similarity score is
        calculated for the whole match, including positions that are not
        present in both motifs. If match is partial or subtotal, only the
        matching psotiions are used to calculate the score. The score of
        individual position is combined using either the mean or the sum.

        Note that the match and combine parameters have no effect on the seqcor
        similarity metric.

        Parameters
        ----------
        m1 : Motif instance
            Motif instance 1.

        m2 : Motif instance
            Motif instance 2.

        match : str, optional
            Match can be "partial", "subtotal" or "total". Not all metrics use
            this.

        metric : str, optional
            Distance metric.

        combine : str, optional
            Combine positional scores using "mean" or "sum". Not all metrics
            use this.

        pval : bool, optional
            Calculate p-vale of match.

        Returns
        -------
        score, position, strand
        """
        if isinstance(metric, str):

            if metric == "seqcor":
                score, pos, orient = seqcor(m1, m2)
                if pval:
                    return [np.nan, pos, orient]
                else:
                    return [score, pos, orient]
            elif match == "partial":
                if pval:
                    return self.pvalue(
                        m1,
                        m2,
                        "total",
                        metric,
                        combine,
                        self.max_partial(m1.pwm, m2.pwm, metric, combine),
                    )
                elif metric in ["pcc", "ed", "distance", "wic", "chisq", "ssd"]:
                    return self.max_partial(m1.pwm, m2.pwm, metric, combine)
                else:
                    return self.max_partial(m1.pfm, m2.pfm, metric, combine)

            elif match == "total":
                if pval:
                    return self.pvalue(
                        m1,
                        m2,
                        match,
                        metric,
                        combine,
                        self.max_total(m1.pwm, m2.pwm, metric, combine),
                    )
                elif metric in ["pcc", "akl"]:
                    # Slightly randomize the weight matrix
                    return self.max_total(
                        m1.wiggle_pwm(), m2.wiggle_pwm(), metric, combine
                    )
                elif metric in ["ed", "distance", "wic", "chisq", "pcc", "ssd"]:
                    return self.max_total(m1.pwm, m2.pwm, metric, combine)
                else:
                    return self.max_total(m1.pfm, m2.pfm, metric, combine)

            elif match == "subtotal":
                if metric in ["pcc", "ed", "distance", "wic", "chisq", "ssd"]:
                    return self.max_subtotal(m1.pwm, m2.pwm, metric, combine)
                else:
                    return self.max_subtotal(m1.pfm, m2.pfm, metric, combine)
        else:
            return metric(m1, m2)

    def _check_length(self, length):
        # Set the length to a length represented in randomly generated JASPAR motifs
        if length < 4:
            return 4
        if length == 13:
            return 14
        if length == 17:
            return 18
        if length == 19:
            return 20
        if length == 21:
            return 22
        if length > 22:
            return 30
        return length

    def pvalue(self, m1, m2, match, metric, combine, score):
        l1, l2 = len(m1.pwm), len(m2.pwm)

        l1 = self._check_length(l1)
        l2 = self._check_length(l2)

        m, s = self.scoredist[metric]["%s_%s" % (match, combine)][l1][l2]

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
                except Exception:
                    raise Exception("Unknown metric '{}'".format(metric))

            scores = []
            for pos1, pos2 in zip(matrix1, matrix2):
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

        # return c_max_subtotal(matrix1, matrix2, metric, combine)

        for i in range(-(len(matrix2) - min_overlap), len(matrix1) - min_overlap + 1):
            p1, p2 = self.make_equal_length_truncate(matrix1, matrix2, i)
            s = self.score_matrices(p1, p2, metric, combine)
            if s:
                scores.append([s, i, 1])

        rev_matrix2 = [row[::-1] for row in matrix2[::-1]]
        for i in range(-(len(matrix2) - min_overlap), len(matrix1) - min_overlap + 1):
            p1, p2 = self.make_equal_length_truncate(matrix1, rev_matrix2, i)
            s = self.score_matrices(p1, p2, metric, combine)
            if s:
                scores.append([s, i, -1])

        if not scores:
            return []
        return sorted(scores, key=lambda x: x[0])[-1]

    def max_partial(self, matrix1, matrix2, metric, combine):

        scores = []

        for i in range(-(len(matrix2) - 1), len(matrix1)):
            p1, p2 = self.make_equal_length_truncate_second(matrix1, matrix2, i)
            s = self.score_matrices(p1, p2, metric, combine)
            if s:
                scores.append([s, i, 1])

        rev_matrix2 = [row[::-1] for row in matrix2[::-1]]
        for i in range(-(len(matrix2) - 1), len(matrix1)):
            p1, p2 = self.make_equal_length_truncate_second(matrix1, rev_matrix2, i)
            s = self.score_matrices(p1, p2, metric, combine)
            if s:
                scores.append([s, i, -1])

        if not scores:
            return []
        return sorted(scores, key=lambda x: x[0])[-1]

    def max_total(self, matrix1, matrix2, metric, combine):
        scores = []

        for i in range(-(len(matrix2) - 1), len(matrix1)):
            p1, p2 = self.make_equal_length(matrix1, matrix2, i)
            s = self.score_matrices(p1, p2, metric, combine)
            if s:
                scores.append([s, i, 1])

        rev_matrix2 = [row[::-1] for row in matrix2[::-1]]
        for i in range(-(len(matrix2) - 1), len(matrix1)):
            p1, p2 = self.make_equal_length(matrix1, rev_matrix2, i)
            s = self.score_matrices(p1, p2, metric, combine)
            if s:
                scores.append([s, i, -1])

        if not scores:
            sys.stdout.write("No score {} {}".format(matrix1, matrix2))
            return []
        return sorted(scores, key=lambda x: x[0])[-1]

    def make_equal_length(self, pwm1, pwm2, pos, bg=None):
        if bg is None:
            bg = [0.25, 0.25, 0.25, 0.25]

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

        return p1, p2

    def make_equal_length_truncate(self, pwm1, pwm2, pos):
        p1 = pwm1[:]
        p2 = pwm2[:]

        if pos < 0:
            p2 = p2[-pos:]
        elif pos > 0:
            p1 = p1[pos:]

        if len(p1) > len(p2):
            p1 = p1[: len(p2)]
        else:
            p2 = p2[: len(p1)]
        return p1, p2

    def make_equal_length_truncate_second(self, pwm1, pwm2, pos, bg=None):
        if bg is None:
            bg = [0.25, 0.25, 0.25, 0.25]

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
            p2 = p2[: len(p1)]
        return p1, p2

    def get_all_scores(
        self,
        motifs,
        dbmotifs,
        match,
        metric,
        combine,
        pval=False,
        parallel=True,
        trim=None,
        ncpus=None,
    ):
        """Pairwise comparison of a set of motifs compared to reference motifs.

        Parameters
        ----------
        motifs : list
            List of Motif instances.

        dbmotifs : list
            List of Motif instances.

        match : str
            Match can be "partial", "subtotal" or "total". Not all metrics use
            this.

        metric : str
            Distance metric.

        combine : str
            Combine positional scores using "mean" or "sum". Not all metrics
            use this.

        pval : bool , optional
            Calculate p-vale of match.

        parallel : bool , optional
            Use multiprocessing for parallel execution. True by default.

        trim : float or None
            If a float value is specified, motifs are trimmed used this IC
            cutoff before comparison.

        ncpus : int or None
            Specifies the number of cores to use for parallel execution.

        Returns
        -------
        scores : dict
            Dictionary with scores.
        """
        # trim motifs first, if specified
        if trim:
            for m in motifs:
                m.trim(trim)
            for m in dbmotifs:
                m.trim(trim)

        # hash of result scores
        scores = {}

        if parallel:
            # Divide the job into big chunks, to keep parallel overhead to minimum
            # Number of chunks = number of processors available
            if ncpus is None:
                ncpus = int(MotifConfig().get_default_params()["ncpus"])

            pool = Pool(processes=ncpus, maxtasksperchild=1000)

            batch_len = len(dbmotifs) // ncpus
            if batch_len <= 0:
                batch_len = 1
            jobs = []
            for i in range(0, len(dbmotifs), batch_len):
                # submit jobs to the job server

                p = pool.apply_async(
                    _get_all_scores,
                    args=(
                        self,
                        motifs,
                        dbmotifs[i : i + batch_len],
                        match,
                        metric,
                        combine,
                        pval,
                    ),
                )
                jobs.append(p)

            pool.close()
            for job in jobs:
                # Get the job result
                result = job.get()
                # and update the result score
                for m1, v in result.items():
                    for m2, s in v.items():
                        if m1 not in scores:
                            scores[m1] = {}
                        scores[m1][m2] = s

            pool.join()
        else:
            # Do the whole thing at once if we don't want parallel
            scores = _get_all_scores(
                self, motifs, dbmotifs, match, metric, combine, pval
            )

        return scores

    def get_best_matches(
        self,
        motifs,
        nmatches=1,
        dbmotifs=None,
        match="partial",
        metric="wic",
        combine="mean",
        parallel=True,
        ncpus=None,
    ):
        """Return best match in database for motifs.

        Parameters
        ----------
        motifs : list or str
            Filename of motifs or list of motifs.

        nmatches : int, optional
            Number of matches to return, default is 1.

        dbmotifs : list or str, optional
            Database motifs, default will be used if not specified.

        match : str, optional

        metric : str, optional

        combine : str, optional

        ncpus : int, optional
            Number of threads to use.

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

        scores = self.get_all_scores(
            motifs, dbmotifs, match, metric, combine, parallel=parallel, ncpus=ncpus
        )
        for motif in scores:
            scores[motif] = sorted(
                scores[motif].items(), key=lambda x: x[1][0], reverse=True
            )[:nmatches]

        ret_scores = {}
        for motif in motifs:
            ret_scores[motif.id] = []
            for dbmotif, match_score in scores[motif.id]:
                pval, pos, orient = self.compare_motifs(
                    motif, dbmotif_lookup[dbmotif], match, metric, combine, True
                )

                ret_scores[motif.id].append([dbmotif, (list(match_score) + [pval])])
        return ret_scores

    def get_closest_match(
        self,
        motifs,
        dbmotifs=None,
        match="partial",
        metric="wic",
        combine="mean",
        parallel=True,
        ncpus=None,
    ):
        """Return best match in database for motifs.

        Parameters
        ----------
        motifs : list or str
            Filename of motifs or list of motifs.

        dbmotifs : list or str, optional
            Database motifs, default will be used if not specified.

        match : str, optional

        metric : str, optional

        combine : str, optional

        ncpus : int, optional
            Number of threads to use.

        Returns
        -------
        closest_match : dict
        """
        scores = self.get_best_matches(
            motifs,
            nmatches=1,
            dbmotifs=dbmotifs,
            match=match,
            metric=metric,
            combine=combine,
            parallel=parallel,
            ncpus=ncpus,
        )

        return dict([k, v[0]] for k, v in scores.items())

    def generate_score_dist(self, motifs, match, metric, combine):

        score_file = os.path.join(
            self.config.get_score_dir(),
            "%s_%s_%s_score_dist.txt" % (match, metric, combine),
        )
        f = open(score_file, "w")

        all_scores = {}
        for motif_len in [len(motif) for motif in motifs]:
            all_scores[motif_len] = {}

        sorted_motifs = {}
        for motif_len in all_scores.keys():
            sorted_motifs[motif_len] = [
                motif for motif in motifs if len(motif) == motif_len
            ]

        for l1 in all_scores.keys():
            for l2 in all_scores.keys():
                scores = self.get_all_scores(
                    sorted_motifs[l1], sorted_motifs[l2], match, metric, combine
                )
                scores = [[y[0] for y in x.values() if y] for x in scores.values()]
                scores = np.array(scores).ravel()
                f.write("%s\t%s\t%s\t%s\n" % (l1, l2, np.mean(scores), np.std(scores)))

        f.close()


def select_nonredundant_motifs(
    roc_report, pfmfile, fg_table, bg_table, tolerance=0.001
):
    pfmfile = pfmfile_location(pfmfile)
    motifs = read_motifs(pfmfile)
    motif_dict = read_motifs(pfmfile, as_dict=True)

    mc = MotifComparer()

    df = pd.read_csv(roc_report, sep="\t", index_col=0)
    df = df[df["Enr. at 1% FPR"] >= 2]
    motifs = [m for m in motifs if m.id in df.index]

    cols = ["ROC AUC", "PR AUC", "Enr. at 1% FPR", "Recall at 10% FDR"]
    rank = df[cols].rank().mean(1).sort_values(ascending=False)

    redundant_motifs = []
    keep = []
    while df[~df.index.isin(redundant_motifs)].shape[0] > 0:
        motif = rank[~rank.index.isin(redundant_motifs)].head(1).index[0]
        keep.append(motif)

        result = mc.get_all_scores(
            [motif_dict[motif]],
            [m for m in motifs if m.id not in redundant_motifs],
            "partial",
            "seqcor",
            "mean",
        )
        result = result[motif]
        redundant_motifs += [m for m in result.keys() if result[m][0] >= 0.7]
    logger.debug(f"Selected {len(keep)} motifs for feature elimination")

    # Read motif scan results
    fg_table = pd.read_csv(fg_table, index_col=0, comment="#", sep="\t")
    bg_table = pd.read_csv(bg_table, index_col=0, comment="#", sep="\t")

    X = pd.concat((fg_table, bg_table), axis=0)
    y = np.hstack((np.ones(fg_table.shape[0]), np.zeros(bg_table.shape[0])))

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y,
        test_size=0.4,
        random_state=2,
        shuffle=True,
    )

    X_bla = X_train[keep]
    model = LogisticRegression(solver="liblinear", max_iter=500, penalty="l1")
    # = RandomForestClassifier(n_estimators=100)
    max_score = np.mean(
        cross_val_score(model, X_bla, y_train, cv=5, scoring="average_precision")
    )
    mean_scores = []
    step = 1

    logger.info("selecting non-redundant motifs")
    n_features = 1
    for i in range(1, X_bla.shape[1], step):
        rfe = RFE(model, i)
        fit = rfe.fit(X_bla, y_train)
        mean_score = np.mean(
            cross_val_score(
                model,
                X_bla.loc[:, fit.support_],
                y_train,
                cv=5,
                scoring="average_precision",
            )
        )
        if i > 1 and mean_score - mean_scores[-1] < (max_score * tolerance):
            n_features = i - 1
            break
        mean_scores.append(mean_score)

    rfe = RFE(model, n_features)
    fit = rfe.fit(X_bla, y_train)

    selected_features = X_bla.columns[fit.support_]
    model.fit(X_train.loc[:, selected_features], y_train)
    y_pred = model.predict_proba(X_test.loc[:, selected_features])[:, 1]
    pr_auc = average_precision_score(y_test, y_pred)
    roc_auc = roc_auc_score(y_test, y_pred)
    logger.info(
        f"selected {len(selected_features)} non-redundant motifs: ROC AUC {roc_auc:.3f}, PR AUC {pr_auc:.3f}"
    )
    return selected_features


# import here is necessary as workaround
# see: http://stackoverflow.com/questions/18947876/using-python-multiprocessing-pool-in-the-terminal-and-in-code-modules-for-django  # noqa: E501
try:
    from multiprocessing import Pool
except ImportError:
    pass
