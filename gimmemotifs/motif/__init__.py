# Copyright (c) 2009-2019 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
"""Module contain core motif functionality"""
# Python imports
import os
import re
import sys
import random
from math import log, sqrt
from collections import Counter
from warnings import warn

from gimmemotifs.config import MotifConfig, DIRECT_NAME, INDIRECT_NAME
from gimmemotifs.utils import pfmfile_location

# External imports
try:
    import numpy as np
except ImportError:
    pass
import xxhash
import matplotlib.pyplot as plt
import logomaker as lm
import pandas as pd
import iteround


NUCS = "ACGT"


class Motif(object):

    """
    Representation of a transcription factor binding motif.

    Examples
    --------

    >>> motif = Motif([[0,1,0,0], [0.5,0,0,0.5], [0,0,1,0]])
    >>> print(motif.to_ppm())
    >
    0   1   0   0
    0.5 0   0   0.5
    0   0   1   0
    >>> print(motif.to_consensus())
    CwG

    """

    from ._scanning import (
        pwm_scan,
        pwm_scan_all,
        scan,
        scan_all,
        pwm_scan_score,
        pwm_scan_to_gff,
    )

    PSEUDO_PFM_COUNT = 1000  # JASPAR mean
    PSEUDO_PPM = 1e-6
    G = 0.25
    Z = 0.01

    # IUPAC table
    iupac = {
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
        "S": "CG",
        "R": "AG",
        "W": "AT",
        "Y": "CT",
        "K": "GT",
        "M": "AC",
        "B": "CGT",
        "H": "ACT",
        "D": "AGT",
        "V": "ACG",
        "N": "ACTG",
    }
    iupac_rev = {v: k for k, v in iupac.items()}
    iupac_ppm = {
        k: [int(nuc in v) / len(v) for nuc in "ACGT"] for k, v in iupac.items()
    }

    def __init__(self, pfm=None, places=4):

        self._places = places
        self.pfm = pfm

        self.factors = {DIRECT_NAME: [], INDIRECT_NAME: []}

        self.id = "unnamed_motif"
        self.config = MotifConfig()

    @property
    def pfm(self):
        return self._pfm

    @property
    def pwm(self):
        return self._ppm

    @property
    def ppm(self):
        return self._ppm

    @property
    def logodds(self):
        return self._logodds

    @pfm.setter
    def pfm(self, mtx):
        if mtx is not None and len(mtx) > 0:
            if np.sum(mtx[0]) > 2:
                self._pfm = [list(x) for x in mtx]
                self._ppm = self.pfm_to_ppm(mtx)
                self._ppm = [iteround.saferound(x, self._places) for x in self._ppm]
            else:
                self._ppm = [iteround.saferound(list(x), self._places) for x in mtx]
                self._pfm = [[n * self.PSEUDO_PFM_COUNT for n in col] for col in mtx]
        else:
            self._ppm = []
            self._pfm = []

        self._logodds = [
            [np.log(n / self.G + self.Z) for n in col] for col in self._ppm
        ]

        self._pfm = np.array(self._pfm)
        self._ppm = np.array(self._ppm)
        self._logodds = np.array(self._logodds)
        self._consensus = self.to_consensus(self.ppm)
        if len(self) > 0:
            self._max_score = self.logodds.max(1).sum()
            self._min_score = self.logodds.min(1).sum()
        else:
            self._max_score = 0
            self._min_score = 0

    @property
    def consensus(self):
        """Motif converted to consensus sequence.

        Returns
        -------
        str
            Consensus sequence.
        """
        if not hasattr(self, "_consensus"):
            self._consensus = self.to_consensus(self.ppm)

        return self._consensus

    def to_consensus(self, ppm=None, precision=4):
        """Convert position probability matrix to consensus sequence.

        Parameters
        ----------
        ppm : array_like, optional
            If not supplied, the ppm of the Motif object will be used.

        precision : int, optional
            Precision used for rounding.

        Returns
        -------
        str
            Consensus sequence.
        """
        if ppm is None:
            ppm = self.ppm

        if len(ppm) == 0:
            return ""

        consensus = ""
        for row in ppm:
            weights = sorted(zip(NUCS, row), key=lambda x: x[1])
            if (
                round(weights[-1][1], precision) >= 0.5
                and weights[-1][1] > 2 * weights[-2][1]
            ):
                consensus += weights[-1][0]
            elif (
                round(weights[-1][1], precision) + round(weights[-2][1], precision)
                >= 0.75
            ):
                consensus += self.iupac_rev[
                    "".join(sorted([weights[-1][0], weights[-2][0]]))
                ].lower()
            else:
                consensus += "n"

        return consensus

    @property
    def max_score(self):
        """Return the maximum logodds score.

        Returns
        -------
        score : float
            Maximum logodds score.
        """
        return self._max_score

    @property
    def min_score(self):
        """Return the minimum logodds score.

        Returns
        -------
        score : float
            Minimum logodds score.
        """
        return self._min_score

    def pwm_min_score(self):
        """Return the minimum PWM score.

        DEPRECATED: use min_score instead.

        Returns
        -------
        score : float
            Minimum PWM score.
        """
        warn(
            "The pwm_min_score() Function is deprecated and will be removed in future release."
            "Please use the min_score property instead.",
            DeprecationWarning,
        )
        return self.min_score

    def pwm_max_score(self):
        """Return the maximum PWM score.

        DEPRECATED: use max_score instead.

        Returns
        -------
        score : float
            Maximum PWM score.
        """
        warn(
            "The pwm_max_score() Function is deprecated and will be removed in future release."
            "Please use the max_score property instead.",
            DeprecationWarning,
        )

        return self.max_score

    def __getitem__(self, x):
        """
        Take slice of a motif and return as new Motif instance.

        Returns
        -------
        motif : Motif instance
            Slice of the motif.
        """
        return Motif(pfm=self.pfm[x])

    def __len__(self):
        """
        Return the motif length.

        Returns
        -------
        int
            Motif length.
        """
        return self.pfm.shape[0]

    def __repr__(self):
        return "{}_{}".format(self.id, self.to_consensus())

    @property
    def information_content(self):
        """Return the total information content of the motif.

        Returns
        -------
        float
            Motif information content.
        """
        return ((self.ppm * np.log2(self.ppm)).sum(1) + 2).sum()

    def score_kmer(self, kmer):
        """Calculate the log-odds score for a specific k-mer.

        Note: this is not necessarily the fastest way for scanning.

        Parameters
        ----------
        kmer : str
            String representing a kmer. Should be the same length as the motif.

        Returns
        -------
        score : float
            Log-odd score.
        """
        if len(kmer) != len(self):
            raise ValueError(
                f"Length of the k-mer should be the same as the motif length ({len(self)})"
            )

        score = self.logodds[np.arange(len(self)), [NUCS.index(n) for n in kmer]].sum()

        return score

    def pfm_to_ppm(self, pfm, pseudo=0.001):
        """Convert PFM with counts to a PFM with fractions (PPM).

        Parameters
        ----------
        pfm : array_like
            2-dimensional array_like with counts.
        pseudo : float
            Pseudocount used in conversion.

        Returns
        -------
        array_like
            2-dimensional array with probability count matrix
        """
        return [
            [(x + pseudo) / (float(np.sum(row)) + pseudo * 4) for x in row]
            for row in pfm
        ]

    def to_motevo(self):
        """Return motif formatted in MotEvo (TRANSFAC-like) format

        Returns
        -------
        m : str
            String of motif in MotEvo format.
        """
        m = "//\n"
        m += "NA {}\n".format(self.id)
        m += "P0\tA\tC\tG\tT\n"
        for i, row in enumerate(self.pfm):
            m += "{}\t{}\n".format(i, "\t".join([str(int(x)) for x in row]))
        m += "//"
        return m

    def to_transfac(self):
        """Return motif formatted in TRANSFAC format

        Returns
        -------
        m : str
            String of motif in TRANSFAC format.
        """
        m = "%s\t%s\t%s\n" % ("DE", self.id, "unknown")
        for i, (row, cons) in enumerate(zip(self.pfm, self.to_consensus())):
            m += "%i\t%s\t%s\n" % (i, "\t".join([str(int(x)) for x in row]), cons)
        m += "XX\n//"
        return m

    def to_meme(self):
        """Return motif formatted in MEME format

        Returns
        -------
        m : str
            String of motif in MEME format.
        """
        motif_id = self.id.replace(" ", "_")
        if motif_id == "":
            motif_id = "unnamed"
        m = "MOTIF %s\n" % motif_id
        m += "BL   MOTIF %s width=0 seqs=0\n" % motif_id
        m += "letter-probability matrix: alength= 4 w= %s nsites= %s E= 0\n" % (
            len(self),
            np.sum(self.pfm[0]),
        )
        m += "\n".join(["\t".join(["%s" % x for x in row]) for row in self.ppm])
        return m

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

    def pcc_pos(self, row1, row2):
        """Calculate the Pearson correlation coefficient of one position
        compared to another position.

        Returns
        -------
        score : float
            Pearson correlation coefficient.
        """
        mean1 = np.mean(row1)
        mean2 = np.mean(row2)

        a = 0
        x = 0
        y = 0
        for n1, n2 in zip(row1, row2):
            a += (n1 - mean1) * (n2 - mean2)
            x += (n1 - mean1) ** 2
            y += (n2 - mean2) ** 2

        if a == 0:
            return 0
        else:
            return a / sqrt(x * y)

    def rc(self):
        """Return the reverse complemented motif.

        Returns
        -------
        Motif instance
            New Motif instance with the reverse complement of the input motif.
        """
        m = Motif(pfm=self.pfm[::-1, ::-1])
        m.id = self.id + "_revcomp"
        return m

    def trim(self, edge_ic_cutoff=0.4):
        """Trim positions with an information content lower than the threshold.

        The default threshold is set to 0.4. The Motif will be changed in-place.

        Parameters
        ----------
        edge_ic_cutoff : float, optional
            Information content threshold. All motif positions at the flanks
            with an information content lower thab this will be removed.

        Returns
        -------
        m : Motif instance
        """
        left_idx = 0
        while left_idx < len(self) and self.ic_pos(self.ppm[left_idx]) < edge_ic_cutoff:
            left_idx += 1

        right_idx = len(self)
        while (
            right_idx > left_idx
            and self.ic_pos(self.ppm[right_idx - 1]) < edge_ic_cutoff
        ):
            right_idx -= 1

        self.pfm = self.pfm[left_idx:right_idx]
        return self

    def consensus_scan(self, fa):
        """Scan FASTA with the motif as a consensus sequence.

        Parameters
        ----------
        fa : Fasta object
            Fasta object to scan

        Returns
        -------
        matches : dict
            Dictionaru with matches.
        """
        regexp = "".join(
            ["[" + "".join(self.iupac[x.upper()]) + "]" for x in self.to_consensusv2()]
        )
        p = re.compile(regexp)
        matches = {}
        for name, seq in fa.items():
            matches[name] = []
            for match in p.finditer(seq):
                middle = (match.span()[1] + match.span()[0]) / 2
                matches[name].append(middle)
        return matches

  
    def plot_logo(
        self,
        kind="information",
        fname=None,
        title=True,
        ylabel=True,
        add_left=0,
        ax=None,
    ):
        """Plot motif logo

        Parameters
        ----------
        kind : str, optional
            Type of logo to plot, can be 'information', 'frequency', 'energy' or
            'ensembl'.
        fname : str, optional
            If fname is set, the plot will be saved with fname as filename.
        title : bool, optional
            Plot the motif id as the title.
        ylabel : bool, optional
            Plot the Y axis label.
        add_left : int, optional
            Add non-informative positions to the left (to align logo)
        """
        fig_height = 3
        fig_width = 0.45

        total = sum(self.pfm[0]) / 4
        pfm = [[total] * 4] * add_left + self.pfm
        matrix = pd.DataFrame(pfm, columns=["A", "C", "G", "T"])

        if kind == "ensembl":
            self.plot_ensembl_logo(fname=fname, title=title)
            return

        logo_params = {
            "information": {
                "df": lm.transform_matrix(
                    matrix, from_type="counts", to_type="information"
                ),
                "figsize": (fig_width * matrix.shape[0], fig_height),
                "show_spines": False,
                "vpad": 0.02,
            },
            "frequency": {
                "df": lm.transform_matrix(
                    matrix, from_type="counts", to_type="probability"
                ),
                "figsize": (fig_width * matrix.shape[0], fig_height),
                "show_spines": False,
                "vpad": 0.02,
                "font_name": "DejaVu Sans Mono",
            },
            "energy": {
                "df": lm.transform_matrix(
                    lm.transform_matrix(matrix, from_type="counts", to_type="weight"),
                    center_values=True,
                ),
                "figsize": (fig_width * matrix.shape[0], fig_height * 2),
                "fade_below": 0.7,
                "shade_below": 0.3,
                "flip_below": False,
                "show_spines": False,
            },
        }

        if ax is not None:
            logo_params[kind]["ax"] = ax
            del logo_params[kind]["figsize"]

        logo = lm.Logo(**logo_params[kind])

        if kind == "information":
            if ylabel:
                logo.ax.set_ylabel("Bits", fontsize=16)
            logo.ax.set_ylim(0, 2)
            logo.ax.set_yticks([0, 0.5, 1, 1.5, 2], minor=False)
        elif kind == "frequency":
            if ylabel:
                logo.ax.set_ylabel("Frequency", fontsize=16)
            logo.ax.set_ylim(0, 1)
        elif kind == "energy":
            if ylabel:
                logo.ax.set_ylabel(r"$\Delta \Delta G$/RT", labelpad=-1, fontsize=16)
        else:
            raise ValueError("Unknown motif visualization")

        logo.ax.set_xticks(range(matrix.shape[0]))
        logo.ax.set_xticklabels(range(1, matrix.shape[0] + 1))
        if title:
            logo.ax.set_title(self.id, fontsize=16)

        if fname:
            plt.savefig(fname, dpi=300)
            plt.close()
        else:
            return logo

    def plot_ensembl_logo(
        self, fname=None, ic=True, title=True, letters=True, height=2
    ):
        """Plot motif logo.

        This is an implementation of the logo presented here:
        http://www.ensembl.info/2018/10/15/new-ensembl-motif-features/

        Parameters
        ----------
        fname : str, optional
            If fname is set, the plot will be saved with fname as filename.
        ic : bool, optional
            Use the bit score. If this is set to False, the frequency
            will be used.
        title : bool, optional
            Plot the motif id as the title.
        letters : bool, optional
            Plot the nucleotides in the bars.
        height : float, optional
            Height of the plot.
        """
        width = 0.94

        ppm = self.ppm
        nucs = np.array(["A", "C", "G", "T"])

        neg_matrix = np.zeros((len(ppm), 4))

        pos_matrix = []
        nuc_ppm = []
        for row in ppm:
            if ic:
                ylabel = "bits"
                ic_row = []
                y_max = 2
                for p in row:
                    if p < 0.25:
                        ic_row.append(0)
                    else:
                        ic_row.append(p * np.log2((p) / 0.25))

            else:
                ic_row = row
                ylabel = "frequency"
                y_max = 1
            idx = np.argsort(ic_row)
            pos_matrix.append(np.array(ic_row)[idx])
            nuc_ppm.append(nucs[idx])

        colors = {
            "A": (0.308, 0.709, 0.280),
            "C": (0.145, 0.362, 0.6),
            "G": (0.969, 0.702, 0.172),
            "T": (0.841, 0.158, 0.224),
        }

        # x_max = np.max([np.sum(row) for row in pos_matrix])
        # x_min = -np.min([np.sum(row) for row in neg_matrix])
        # neg_matrix = neg_matrix / x_min * x_max

        plt.figure(figsize=(len(ppm) * 0.3, height))
        for (sign, matrix) in [(1, pos_matrix), (-1, neg_matrix)]:
            minbottom = np.zeros(len(matrix))
            alpha = 1
            if sign == -1:
                minbottom = np.array([np.sum(row) for row in matrix])
                alpha = 0.5

            # Print the bars
            for i in range(0, len(ppm[0])):

                pheight = [abs(r[i]) for r in matrix]
                bottom = minbottom + [sum(np.abs(r[:i])) for r in matrix]

                c = [colors[r[i]] for r in nuc_ppm]
                plt.bar(
                    range(1, len(ppm) + 1),
                    width=width,
                    height=pheight,
                    bottom=bottom,
                    color=c,
                    alpha=alpha,
                )

            if letters:
                # Print the letters
                for i in range(len(ppm)):
                    for n in range(4):
                        x = i + 1
                        y = matrix[i][n] / 2 + sum(matrix[i][:n])
                        nuc = nuc_ppm[i][n]
                        c = "white"
                        if abs(matrix[i][n]) * height >= 0.5:
                            plt.text(
                                x,
                                y,
                                nuc,
                                horizontalalignment="center",
                                verticalalignment="center",
                                fontsize=8 + 10 * (matrix[i][n] / y_max),
                                color=c,
                            )

        # Remove axis lines
        ax = plt.gca()
        for spine in ax.spines.values():
            spine.set_color("none")

        if title:
            plt.title(self.id)
        plt.xlim(0.47, len(self) + 0.5)
        plt.xticks(range(1, len(ppm) + 1))
        plt.ylabel(ylabel)

        if fname:
            plt.savefig(fname, dpi=300)
            plt.close()
        else:
            return ax

    def average_motifs(self, other, pos, orientation, include_bg=False):
        """Return the average of two motifs.

        Combine this motif with another motif and return the average as a new
        Motif object. The position and orientatien need to be supplied. The pos
        parameter is the position of the second motif relative to this motif.

        For example, take the following two motifs:
        Motif 1: CATGYT
        Motif 2: GGCTTGY

        With position -2, the motifs are averaged as follows:
        xxCATGYT
        GGCTTGYx

        Parameters
        ----------
        other : Motif object
            Other Motif object.
        pos : int
            Position of the second motif relative to this motif.
        orientation : int
            Orientation, should be 1 or -1. If the orientation is -1 then the
            reverse complement of the other motif is used for averaging.
        include_bg : bool , optional
            Extend both motifs with background frequencies (0.25) before
            averaging. False by default.

        Returns
        -------
        motif : motif object
            New Motif object containing average motif.
        """
        # xxCATGYT
        # GGCTTGYx
        # pos = -2
        pfm1 = self.pfm[:]
        pfm2 = other.pfm[:]

        if orientation < 0:
            pfm2 = [row[::-1] for row in pfm2[::-1]]

        pfm1_count = float(np.sum(pfm1[0]))
        pfm2_count = float(np.sum(pfm2[0]))

        if include_bg:
            if len(pfm1) > len(pfm2) + pos:
                pfm2 += [
                    [pfm2_count / 4.0 for x in range(4)]
                    for i in range(-(len(pfm1) - len(pfm2) - pos), 0)
                ]
            elif len(pfm2) + pos > len(pfm1):
                pfm1 += [
                    [pfm1_count / 4.0 for x in range(4)]
                    for i in range(-(len(pfm2) - len(pfm1) + pos), 0)
                ]

            if pos < 0:
                pfm1 = [
                    [pfm1_count / 4.0 for x in range(4)] for i in range(-pos)
                ] + pfm1
            elif pos > 0:
                pfm2 = [[pfm2_count / 4.0 for x in range(4)] for i in range(pos)] + pfm2

        else:
            if len(pfm1) > len(pfm2) + pos:
                pfm2 += [
                    [pfm1[i][x] / pfm1_count * (pfm2_count) for x in range(4)]
                    for i in range(-(len(pfm1) - len(pfm2) - pos), 0)
                ]
            elif len(pfm2) + pos > len(pfm1):
                pfm1 += [
                    [pfm2[i][x] / pfm2_count * (pfm1_count) for x in range(4)]
                    for i in range(-(len(pfm2) - len(pfm1) + pos), 0)
                ]

            if pos < 0:
                pfm1 = [
                    [pfm2[i][x] / pfm2_count * (pfm1_count) for x in range(4)]
                    for i in range(-pos)
                ] + pfm1
            elif pos > 0:
                pfm2 = [
                    [pfm1[i][x] / pfm1_count * (pfm2_count) for x in range(4)]
                    for i in range(pos)
                ] + pfm2

        pfm = [[a + b for a, b in zip(x, y)] for x, y in zip(pfm1, pfm2)]

        m = Motif(pfm)
        m.id = m.to_consensus()
        return m

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
        # xxCATGYT
        # GGCTTGYx
        # pos = -2
        ppm1 = ppm1[:]
        ppm2 = ppm2[:]

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
            score += self.pcc_pos(a, b)

        return score

    def ic(self, ppm1, ppm2, pos, bg=None, bg_factor=1):
        if bg is None:
            bg = [0.25, 0.25, 0.25, 0.25]

        # xxCATGYT
        # GGCTTGYx
        # pos = -2
        ppm1 = ppm1[:]
        ppm2 = ppm2[:]

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
                self.ic_pos(a)
                + self.ic_pos(b)
                - (self.ic_pos(a, b) + self.ic_pos(b, a))
            )

        # print "SCORE: %s" % score
        # Parts aligned to the background
        for x in na:
            score += (
                self.ic_pos(x)
                + self.ic_pos(bg)
                - (self.ic_pos(x, bg) + self.ic_pos(bg, x))
            ) * bg_factor

        #    print "SCORE WITH BG: %s" % score
        return score

    def other_ic(self, ppm1, ppm2, pos, bg=None, bg_factor=0.8):
        if bg is None:
            bg = [0.25, 0.25, 0.25, 0.25]

        # xxCATGYT
        # GGCTTGYx
        # pos = -2
        ppm1 = ppm1[:]
        ppm2 = ppm2[:]

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
                (np.sum(a[ppm1_start:ppm1_end], 1) + np.sum(b[ppm2_start:ppm2_end], 1))
                / 2
                - np.sum(np.abs(a[ppm1_start:ppm1_end] - b[ppm2_start:ppm2_end]), 1)
            )
            scores.append([score, pos, 1])

            score = np.sum(
                (
                    np.sum(a[ppm1_start:ppm1_end], 1)
                    + np.sum(b_rev[ppm2_start:ppm2_end], 1)
                )
                / 2
                - np.sum(np.abs(a[ppm1_start:ppm1_end] - b_rev[ppm2_start:ppm2_end]), 1)
            )
            scores.append([score, pos, -1])

        return sorted(scores, key=lambda x: x[0])[-1]

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

        return sorted(scores, key=lambda x: x[0])[-1]

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

    def _format_jaspar(self, version=1, header=True):
        rows = np.array(self.ppm).transpose()
        rows = [" ".join([str(x) for x in row]) for row in rows]
        if version == 2:
            rows = ["{} [{} ]".format(n, row) for n, row in zip(NUCS, rows)]

        str_out = "\n".join(rows)
        if header:
            str_out = "\n".join([self.id, str_out])

        return str_out

    def to_consensusv2(self):
        if self.consensus:
            return self.consensus
        else:
            consensus = ""
            for row in self.ppm:
                weights = sorted(zip(["A", "C", "G", "T"], row), key=lambda x: x[1])
                if weights[-1][1] >= 0.5:
                    if weights[-2][1] >= 0.25:
                        consensus += self.iupac_rev[
                            "".join(sorted([weights[-1][0], weights[-2][0]]))
                        ]
                    else:
                        consensus += weights[-1][0]
                elif weights[-1][1] + weights[-2][1] >= 0.75:
                    consensus += self.iupac_rev[
                        "".join(sorted([weights[-1][0], weights[-2][0]]))
                    ]
                elif weights[-1][1] + weights[-2][1] + weights[-3][1] >= 0.9:
                    consensus += self.iupac_rev[
                        "".join(
                            sorted([weights[-1][0], weights[-2][0], weights[-3][0]])
                        )
                    ]
                else:
                    consensus += "n"
            return consensus

    def to_pfm(self):
        if self.pfm:
            return ">%s\n%s" % (
                self.id,
                "\n".join(["\t".join(["%s" % x for x in row]) for row in self.pfm]),
            )
        else:
            pfm = [[n * self.PSEUDO_PFM_COUNT for n in col] for col in self.ppm]
            return ">%s\n%s" % (
                self.id,
                "\n".join(["\t".join(["%s" % x for x in row]) for row in pfm]),
            )

    def _ppm_to_str(self, precision=4):
        """Return string representation of ppm.

        Parameters
        ----------
        precision : int, optional, default 4
            Floating-point precision.

        Returns
        -------
        ppm_string : str
        """
        if not self.ppm:
            return ""

        fmt = "{{:.{:d}f}}".format(precision)
        return "\n".join(["\t".join([fmt.format(p) for p in row]) for row in self.ppm])

    def hash(self):
        """Return hash of motif.

        This is an unique identifier of a motif, regardless of the id.

        Returns:
        hash : str
        """
        return xxhash.xxh64(self._ppm_to_str(3)).hexdigest()

    def to_ppm(self, precision=4, extra_str=""):
        """Return ppm as string.

        Parameters
        ----------
        precision : int, optional, default 4
            Floating-point precision.

        extra_str |: str, optional
            Extra text to include with motif id line.

        Returns
        -------
        motif_str : str
            Motif formatted in ppm format.
        """
        motif_id = self.id

        if extra_str:
            motif_id += "_%s" % extra_str

        if not self.ppm:
            self.ppm = [self.iupac_ppm[char] for char in self.consensus.upper()]

        return ">%s\n%s" % (motif_id, self._ppm_to_str(precision))

    def to_img(self, fname, fmt="PNG", add_left=0, seqlogo=None, height=6):
        """Create a sequence logo using seqlogo.

        Create a sequence logo and save it to a file. Valid formats are: PNG,
        EPS, GIF and PDF.

        Parameters
        ----------
        fname : str
            Output filename.
        fmt : str , optional
            Output format (case-insensitive). Valid formats are PNG, EPS, GIF
            and PDF.
        add_left : int , optional
            Pad motif with empty positions on the left side.
        seqlogo : str
            Location of the seqlogo executable. By default the seqlogo version
            that is included with GimmeMotifs is used.
        height : float
            Height of the image
        """
        warn("Method to_img() is replaced by plot_logo()", DeprecationWarning)
        self.plot_logo(fname=fname, kind="information")

    def randomize(self):
        """Create a new motif with shuffled positions.

        Shuffle the positions of this motif and return a new Motif instance.

        Returns
        -------
        m : Motif instance
            Motif instance with shuffled positions.
        """
        random_pfm = [[c for c in row] for row in self.pfm]
        random.shuffle(random_pfm)
        m = Motif(pfm=random_pfm)
        m.id = "random"
        return m

    def randomize_dimer(self):
        length = len(self.pfm)
        random_pfm = []
        for _ in range(length / 2):
            pos = random.randint(0, length - 1)
            random_pfm += [[c for c in row] for row in self.pfm[pos : pos + 2]]
        motif = Motif(pfm=random_pfm)
        motif.id = "random"
        return motif

    def format_factors(
        self, max_length=5, html=False, include_indirect=True, extra_str=", (...)"
    ):
        if html:
            fmt_d = "<span style='color:black'>{}</span>"
            fmt_i = "<span style='color:#666666'>{}</span>"
        else:
            fmt_d = fmt_i = "{}"

        if hasattr(self, "factor_info"):
            fcount = Counter([x.upper() for x in self.factor_info.get("Factor", "")])
        else:
            fcount = Counter(self.factors[DIRECT_NAME] + self.factors[INDIRECT_NAME])

        direct = sorted(
            list(
                set(
                    [
                        x.upper() if x != "de novo" else x
                        for x in self.factors[DIRECT_NAME]
                    ]
                )
            ),
            key=lambda x: fcount[x],
            reverse=True,
        )

        indirect = []
        if include_indirect:
            indirect = sorted(
                list(
                    set(
                        [
                            x.upper()
                            for x in self.factors[INDIRECT_NAME]
                            if x.upper() not in direct
                        ]
                    )
                ),
                key=lambda x: fcount[x],
                reverse=True,
            )

        if len(direct) > max_length:
            show_factors = direct[:max_length]
        else:
            show_factors = direct[:]
            for f in sorted(indirect, key=lambda x: fcount[x], reverse=True):
                if f not in show_factors:
                    show_factors.append(f)
                if len(show_factors) >= max_length:
                    break

        if "de novo" in show_factors:
            show_factors = ["de novo"] + sorted(
                [f for f in show_factors if f != "de novo"],
                key=lambda x: fcount[x],
                reverse=True,
            )
        else:
            show_factors = sorted(show_factors, key=lambda x: fcount[x], reverse=True)

        factor_str = ",".join(
            [fmt_d.format(f) if f in direct else fmt_i.format(f) for f in show_factors]
        )

        if len(direct + indirect) > max_length:
            factor_str += extra_str

        if html:
            tooltip = ""
            if len(direct) > 0:
                tooltip += "direct: " + ",".join(sorted(direct))
            if len(indirect) > 0:
                if tooltip != "":
                    tooltip += "&#10;"
                tooltip += "predicted: " + ",".join(sorted(indirect))

            factor_str = '<div title="' + tooltip + '">' + factor_str + "</div>"

        return factor_str

    def sample(self, n_seqs):
        """Sample n_seqs random sequences from a motif. The sequences
        follow the distribution of the motif ppm.

        Parameters
        ----------
        n_seqs : int
            number of sequences to sample

        Returns
        -------
        sequences : List[str]
            A list of all the samples sequences
        """
        nucs = [
            random.choices(NUCS, weights=self.ppm[i], k=n_seqs)
            for i in range(len(self.ppm))
        ]
        seqs = ["".join(nuc) for nuc in zip(*nucs)]
        return seqs


def default_motifs():
    """Return list of Motif instances from default motif database."""
    config = MotifConfig()
    d = config.get_motif_dir()
    m = config.get_default_params()["motif_db"]

    if not d or not m:
        raise ValueError("default motif database not configured")

    fname = os.path.join(d, m)
    with open(fname) as f:
        motifs = read_motifs(f)

    return motifs


def motif_from_align(align):
    """Convert alignment to motif.

    Converts a list with sequences to a motif. Sequences should be the same
    length.

    Parameters
    ----------
    align : list
        List with sequences (A,C,G,T).

    Returns
    -------
    m : Motif instance
        Motif created from the aligned sequences.
    """
    width = len(align[0])
    nucs = {"A": 0, "C": 1, "G": 2, "T": 3}
    pfm = [[0 for _ in range(4)] for _ in range(width)]
    for row in align:
        for i in range(len(row)):
            pfm[i][nucs[row[i]]] += 1
    m = Motif(pfm)
    m.align = align[:]
    return m


def motif_from_consensus(cons, n=1200):
    """Convert consensus sequence to motif.

    Converts a consensus sequences using the nucleotide IUPAC alphabet to a
    motif.

    Parameters
    ----------
    cons : str
        Consensus sequence using the IUPAC alphabet.
    n : int , optional
        Count used to convert the sequence to a PFM.

    Returns
    -------
    m : Motif instance
        Motif created from the consensus.
    """
    width = len(cons)
    nucs = {"A": 0, "C": 1, "G": 2, "T": 3}
    pfm = [[0 for _ in range(4)] for _ in range(width)]
    m = Motif()
    for i, char in enumerate(cons):
        for nuc in m.iupac[char.upper()]:
            pfm[i][nucs[nuc]] = n / len(m.iupac[char.upper()])
    m = Motif(pfm)
    m.id = cons
    return m


def parse_motifs(motifs):
    """Parse motifs in a variety of formats to return a list of motifs.

    Parameters
    ----------

    motifs : list or str
        Filename of motif,  list of motifs or single Motif instance.

    Returns
    -------

    motifs : list
        List of Motif instances.
    """
    if isinstance(motifs, str):
        return read_motifs(motifs)
    elif isinstance(motifs, Motif):
        motifs = [motifs]
    else:
        if not isinstance(list(motifs)[0], Motif):
            raise ValueError("Not a list of motifs")

    return list(motifs)


def _add_factors_from_handle(motifs, handle):
    """Add factors to motifs.

    Reads the factor-motif association from a "motif2factors.txt" file.
    """
    if not (hasattr(handle, "name") and handle.name):
        return motifs

    base = os.path.splitext(handle.name)[0]
    map_file = base + ".motif2factors.txt"
    if not os.path.exists(map_file):
        return motifs

    m2f_direct = {}
    m2f_indirect = {}
    for line in open(map_file):
        if line.startswith("#"):
            continue
        try:
            motif, *factor_info = line.strip().split("\t")
            if len(factor_info) == 1:
                m2f_direct[motif] = factor_info[0].split(",")
            elif len(factor_info) == 3:
                if factor_info[2] == "Y":
                    m2f_direct[motif] = m2f_direct.get(motif, []) + [factor_info[0]]
                else:
                    m2f_indirect[motif] = m2f_indirect.get(motif, []) + [factor_info[0]]
        except Exception:
            pass

    m2f = pd.read_csv(map_file, sep="\t", comment="#", index_col=0)

    for motif in motifs:
        if motif.id in m2f.index:
            motif.factor_info = m2f.loc[motif.id]
        if motif.id in m2f_direct:
            motif.factors[DIRECT_NAME] = m2f_direct[motif.id]
        if motif.id in m2f_indirect:
            motif.factors[INDIRECT_NAME] = m2f_indirect[motif.id]

    for motif in motifs:
        for n in [DIRECT_NAME, INDIRECT_NAME]:
            motif.factors[n] = list(set(motif.factors[n]))

    return motifs


def _read_motifs_from_filehandle(handle, fmt):
    """
    Read motifs from a file-like object.

    Parameters
    ----------
    handle : file-like object
        Motifs.
    fmt : string, optional
        Motif format, can be 'pfm', 'transfac', 'xxmotif', 'jaspar' or 'align'.

    Returns
    -------
    motifs : list
        List of Motif instances.
    """
    if fmt.lower() == "pfm":
        motifs = _read_motifs_pfm(handle)
    elif fmt.lower() == "transfac":
        motifs = _read_motifs_transfac(handle)
    elif fmt.lower() == "xxmotif":
        motifs = _read_motifs_xxmotif(handle)
    elif fmt.lower() == "align":
        motifs = _read_motifs_align(handle)
    elif fmt.lower() == "jaspar":
        motifs = _read_motifs_jaspar(handle)
    elif fmt.lower() == "meme":
        motifs = _read_motifs_meme(handle)
    else:
        raise ValueError("Unknown format {}".format(fmt))

    # Remove everything after tab from motif identifiers
    for motif in motifs:
        motif.id = motif.id.split("\t")[0]

    motifs = _add_factors_from_handle(motifs, handle)

    return motifs


def read_motifs(infile=None, fmt="pfm", as_dict=False):
    """
    Read motifs from a file or stream or file-like object.

    Parameters
    ----------
    infile : string or file-like object, optional
        Motif database, filename of motif file or file-like object. If infile
        is not specified the default motifs as specified in the config file
        will be returned.

    fmt : string, optional
        Motif format, can be 'pfm', 'transfac', 'xxmotif', 'jaspar' or 'align'.

    as_dict : boolean, optional
        Return motifs as a dictionary with motif_id, motif pairs.

    Returns
    -------
    motifs : list
        List of Motif instances. If as_dict is set to True, motifs is a
        dictionary.
    """
    # Support old naming
    if fmt == "pwm":
        fmt = "pfm"

    if infile is None or isinstance(infile, str):
        infile = pfmfile_location(infile)
        with open(infile) as f:
            motifs = _read_motifs_from_filehandle(f, fmt)
    else:
        motifs = _read_motifs_from_filehandle(infile, fmt)

    if as_dict:
        motifs = {m.id: m for m in motifs}

    return motifs


def _read_motifs_pfm(handle):
    p = re.compile(
        r"(\d+(\.\d+)?(e-\d+)?)\s+(\d+(\.\d+)?(e-\d+)?)\s+(\d+(\.\d+)?(e-\d+)?)\s+(\d+(\.\d+)?(e-\d+)?)"  # noqa: E501
    )
    motifs = []
    pfm = []
    motif_id = ""
    seen_id = {}

    for n, line in enumerate(handle.readlines()):
        if line.startswith("#") or line.strip() == "":
            continue
        if line.startswith(">"):
            if pfm:
                motifs.append(Motif(pfm))
                motifs[-1].id = motif_id
                pfm = []
            motif_id = line.strip()[1:]
            seen_id[motif_id] = seen_id.get(motif_id, 0) + 1
            if seen_id.get(motif_id, 0) > 1:
                msg = "WARNING: multiple motifs with same id: {}\n".format(motif_id)
                sys.stderr.write(msg)
                motif_id += "_{}".format(seen_id[motif_id] - 1)

        else:
            m = p.search(line)
            if m:
                fractions = [float(m.group(x)) for x in (1, 4, 7, 10)]
                pfm.append(fractions)
            else:
                msg = "WARNING: can't parse line {}, ignoring:\n{}".format(n + 1, line)
                sys.stderr.write(msg)

    if len(pfm) > 0:
        motifs.append(Motif(pfm))
        motifs[-1].id = motif_id

    return motifs


def _read_motifs_jaspar(handle):
    p = re.compile(r"([ACGT])\s*\[?([^\]]+)\]?")
    motifs = []
    motif_id = ""
    pfm = {}
    for line in handle:
        line = line.strip()
        if len(line) == 0:
            continue

        if line.startswith(">"):
            motif_id = line[1:]
        if line[0] in NUCS:
            m = p.search(line)
            try:
                nuc = m.group(1)
                counts = re.split(r"\s+", m.group(2).strip())
                pfm[nuc] = [float(x) for x in counts]
                if nuc == "T":
                    motif = Motif(np.array([pfm[n] for n in NUCS]).transpose())
                    motif.id = motif_id
                    motifs.append(motif)
            except Exception:
                raise ValueError("Can't parse line\n" + line)

    if motif_id and motifs[-1].id != motif_id:
        motif = Motif(np.array([pfm[n] for n in NUCS]).transpose())
        motif.id = motif_id
        motifs.append(motif)

    return motifs


def _read_motifs_align(handle):
    motifs = []
    nucs = {"A": 0, "C": 1, "G": 2, "T": 3}
    motif_id = ""
    aligns = {}
    align = []
    for line in handle:
        if line.startswith(">"):
            if motif_id:
                aligns[motif_id] = align
            motif_id = line.strip()[1:]
            align = []
        else:
            align.append(line.strip())
    aligns[motif_id] = align

    for motif_id, align in aligns.items():

        width = len(align[0])
        pfm = [[0 for _ in range(4)] for _ in range(width)]
        for row in align:
            for i in range(len(row)):
                pfm[i][nucs[row[i]]] += 1
        m = Motif(pfm)
        m.align = align[:]
        m.pfm = pfm[:]
        m.id = motif_id
        motifs.append(m)
    return motifs


def _read_motifs_xxmotif(handle):
    motifs = []

    line = handle.readline()
    while line:
        while line and not line.startswith("Motif"):
            line = handle.readline()

        if line:
            mid = line.split(":")[0]
            freqs = []
            for _ in range(4):
                line = handle.readline()
                freqs.append([float(x) for x in line.strip().split("\t")])

            pfm = np.array(freqs).transpose()
            motif = Motif(pfm)
            motif.id = mid.replace(" ", "_")
            motifs.append(motif)

    return motifs


def _read_motifs_meme(handle):
    motifs = []
    line = handle.readline()
    while line:
        while line and not line.startswith("MOTIF"):
            line = handle.readline()

        if line:
            mid = line.strip().split(" ")[1]
            freqs = []
            while not line.startswith("letter"):
                line = handle.readline()
            while line:
                line = handle.readline()
                if not line.strip():
                    break
                row = [float(x) for x in re.split(r"\s", line.strip())]
                freqs.append(row)

            motif = Motif(freqs)
            motif.id = mid.replace(" ", "_")
            motifs.append(motif)

    return motifs


def _read_motifs_transfac(handle):
    p = re.compile(r"\d+\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s*\w?")
    p_id = re.compile(r"^(NA|ID|DE)\s+([^\s]+( [^\s]+)?)")
    motifs = []
    pfm = []
    motif_id = ""
    metadata = None
    for line in handle.readlines():
        m = p_id.search(line.strip())
        if m:
            motif_id = m.group(2)
            motif_id = re.sub(r"\s+", "_", motif_id)
        elif line.startswith("//"):
            if len(pfm) != 0:
                motifs.append(Motif(pfm))
                motifs[-1].id = motif_id
                if metadata is not None:
                    motifs[-1].metadata = metadata
            pfm = []
        elif p.search(line):
            m = p.search(line)
            pfm.append([float(x) for x in m.group(1, 2, 3, 4)])
        elif line.startswith("CC"):
            # CC tax_group:vertebrates; pubmed_ids:7592839; uniprot_ids:P30561,P53762; data_type:SELEX
            metadata = line.replace("CC ", "").split(";")
            metadata = dict([x.strip().split(":") for x in metadata])
    # If there's only one matrix, and the format is not complete
    if len(pfm) != 0:
        motifs.append(Motif(pfm))

    return motifs


def motifs_to_meme(motifs):
    m = "MEME version 3.0\n\nALPHABET= ACGT\n\nstrands: + -\n\n"
    m += "Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n"

    for motif in motifs:
        m += motif.to_meme() + "\n"
    return m