# Copyright (c) 2009-2021 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
"""Comparison functions for Motif class"""
from math import log

import numpy as np
from scipy.stats import pearsonr

from gimmemotifs.utils import make_equal_length


def ic_pos(self, row1, row2=None):
    """Calculate the information content of one position.

    Returns
    -------
    score : float
        Information content.
    """
    if row2 is None:
        row2 = [0.25, 0.25, 0.25, 0.25]

    score = 0
    for a, b in zip(row1, row2):
        if a > 0:
            score += a * log(a / b) / log(2)
    return score


def other_ic_pos(self, row1, row2, bg=None):
    if bg is None:
        bg = [0.25, 0.25, 0.25, 0.25]

    score = 0
    score_a = 0
    score_b = 0

    for a, b, pbg in zip(row1, row2, bg):
        score += abs(a * log(a / pbg) / log(2) - b * log(b / pbg) / log(2))
        score_a += a * log(a / pbg) / log(2)
        score_b += b * log(b / pbg) / log(2)

    return (score_a + score_b) / 2 - score


def pcc(self, ppm1, ppm2, pos):
    """ """
    # xxCATGYT
    # GGCTTGYx
    # pos = -2
    ppm1, ppm2 = make_equal_length(ppm1, ppm2, pos, truncate="both")

    # Compute pearson correlation between aligned parts of the motif
    r = np.array(
        [
            pearsonr(x, y)[0]
            if not (np.all(x == 0.25)) or not (np.all(y == 0.25))
            else 0
            for x, y in zip(ppm1, ppm2)
        ]
    )
    r[np.isnan(r)] = 0

    return r.sum()


def ic(self, ppm1, ppm2, pos, bg=None, bg_factor=1):
    if bg is None:
        bg = [0.25, 0.25, 0.25, 0.25]

    # xxCATGYT
    # GGCTTGYx
    # pos = -2
    ppm1 = ppm1.copy()
    ppm2 = ppm2.copy()

    na = []
    if pos > 0:
        na = ppm1[:pos]
        ppm1 = ppm1[pos:]
    elif pos < 0:
        na = ppm2[:-pos]
        ppm2 = ppm2[-pos:]

    if len(ppm1) > len(ppm2):
        na += ppm1[len(ppm2) :]
        ppm1 = ppm1[: len(ppm2)]
    elif len(ppm2) > len(ppm1):
        na += ppm2[len(ppm1) :]
        ppm2 = ppm2[: len(ppm1)]

    # print "COMPARE"
    # print Motif(ppm1).to_consensus()
    # print Motif(ppm2).to_consensus()

    # Aligned parts of the motif
    score = 0
    for a, b in zip(ppm1, ppm2):
        score += (
            self.ic_pos(a) + self.ic_pos(b) - (self.ic_pos(a, b) + self.ic_pos(b, a))
        )

    # print "SCORE: %s" % score
    # Parts aligned to the background
    for x in na:
        score += (
            self.ic_pos(x) + self.ic_pos(bg) - (self.ic_pos(x, bg) + self.ic_pos(bg, x))
        ) * bg_factor

    #    print "SCORE WITH BG: %s" % score
    return score


def other_ic(self, ppm1, ppm2, pos, bg=None, bg_factor=0.8):
    if bg is None:
        bg = [0.25, 0.25, 0.25, 0.25]

    # xxCATGYT
    # GGCTTGYx
    # pos = -2
    ppm1 = ppm1.copy()
    ppm2 = ppm2.copy

    na = []
    if pos > 0:
        na = ppm1[:pos]
        ppm1 = ppm1[pos:]
    elif pos < 0:
        na = ppm2[:-pos]
        ppm2 = ppm2[-pos:]

    if len(ppm1) > len(ppm2):
        na += ppm1[len(ppm2) :]
        ppm1 = ppm1[: len(ppm2)]
    elif len(ppm2) > len(ppm1):
        na += ppm2[len(ppm1) :]
        ppm2 = ppm2[: len(ppm1)]

    # Aligned parts of the motif
    score = 0
    for a, b in zip(ppm1, ppm2):
        score += self.other_ic_pos(a, b)

    for x in na:
        score += self.other_ic_pos(x, bg) * bg_factor

    return score


def matrix_ic(self, ppm1, ppm2, bg=None):
    if bg is None:
        bg = [0.25, 0.25, 0.25, 0.25]

    # xxCATGYT
    # GGCTTGYx
    # pos = -2
    ppm1 = np.array(ppm1)
    ppm2 = np.array(ppm2)
    ppm2_rev = np.array([row[::-1] for row in ppm2[::-1]])
    bg = np.array(bg)

    a = ppm1 * np.log2(ppm1 / bg)
    b = ppm2 * np.log2(ppm2 / bg)

    b_rev = ppm2_rev * np.log2(ppm2_rev / bg)

    scores = []
    l1 = len(ppm1)
    l2 = len(ppm2)
    for pos in range(-(l2 - 1), l1):

        ppm1_start, ppm2_start = 0, 0
        ppm1_end, ppm2_end = l1, l2
        if pos > 0:
            ppm1_start = pos
            if l1 - pos > l2:
                ppm1_end = l2 + pos
            elif l1 - pos < l2:
                ppm2_end = l1 - pos
        elif pos < 0:
            ppm2_start = -pos
            if l2 + pos > l1:
                ppm2_end = l1
            elif l2 + pos < l1:
                ppm1_end = l2 + pos
        else:
            if l2 > l1:
                ppm2_end = l1
            elif l2 < l1:
                ppm1_end = l2

        score = np.sum(
            (np.sum(a[ppm1_start:ppm1_end], 1) + np.sum(b[ppm2_start:ppm2_end], 1)) / 2
            - np.sum(np.abs(a[ppm1_start:ppm1_end] - b[ppm2_start:ppm2_end]), 1)
        )
        scores.append([score, pos, 1])

        score = np.sum(
            (np.sum(a[ppm1_start:ppm1_end], 1) + np.sum(b_rev[ppm2_start:ppm2_end], 1))
            / 2
            - np.sum(np.abs(a[ppm1_start:ppm1_end] - b_rev[ppm2_start:ppm2_end]), 1)
        )
        scores.append([score, pos, -1])

    return max(scores, key=lambda x: x[0])


def max_ic(self, other, revcomp=True, bg_factor=0.8):
    ppm1 = self.ppm
    ppm2 = other.ppm

    scores = []

    for i in range(-(len(ppm2) - 1), len(ppm1)):
        scores.append((self.other_ic(ppm1, ppm2, i, bg_factor=bg_factor), i, 1))

    if revcomp:
        rev_ppm2 = [row[::-1] for row in ppm2[::-1]]
        for i in range(-(len(ppm2) - 1), len(ppm1)):
            scores.append(
                (self.other_ic(ppm1, rev_ppm2, i, bg_factor=bg_factor), i, -1)
            )

    return max(scores, key=lambda x: x[0])


def max_pcc(self, other, revcomp=True):
    ppm1 = self.ppm
    ppm2 = other.ppm

    scores = []

    for i in range(-(len(ppm2) - 1), len(ppm1)):
        scores.append((self.pcc(ppm1, ppm2, i), i, 1))

    if revcomp:
        rev_ppm2 = [row[::-1] for row in ppm2[::-1]]
        for i in range(-(len(ppm2) - 1), len(ppm1)):
            scores.append((self.pcc(ppm1, rev_ppm2, i), i, -1))

    # print scores
    return sorted(scores, key=lambda x: x[0])[-1]
