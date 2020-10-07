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
from gimmemotifs.c_metrics import pfmscan
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


class Motif(object):

    """
    Representation of a transcription factor binding motif.

    Examples
    --------

    >>> motif = Motif([[0,1,0,0], [0.5,0,0,0.5], [0,0,1,0]])
    >>> print(motif.to_pwm())
    >
    0   1   0   0
    0.5 0   0   0.5
    0   0   1   0
    >>> print(motif.to_consensus())
    CwG

    """

    PSEUDO_PFM_COUNT = 1000  # Jaspar mean
    PSEUDO_PWM = 1e-6
    G = 0.25
    Z = 0.01

    def __init__(self, pfm=None):
        if pfm is None:
            pfm = []

        if len(pfm) > 0:
            if np.sum(pfm[0]) > 2:
                self.pfm = [list(x) for x in pfm]
                self.pwm = self.pfm_to_pwm(pfm)
            else:
                self.pwm = [list(x) for x in pfm]
                self.pfm = [[n * self.PSEUDO_PFM_COUNT for n in col] for col in pfm]
            self.logodds = [
                [np.log(n / self.G + self.Z) for n in col] for col in self.pwm
            ]
        else:
            self.pwm = []
            self.pfm = []

        self.wiggled_pwm = None
        self.factors = {DIRECT_NAME: [], INDIRECT_NAME: []}
        self.seqs = []
        self.consensus = ""
        self.min_score = None
        self.max_score = None

        self.id = ""
        self.config = MotifConfig()

        self.nucs = "ACGT"

        self.iupac_rev = {
            "CG": "S",
            "AG": "R",
            "AT": "W",
            "CT": "Y",
            "GT": "K",
            "AC": "M",
            "CGT": "B",
            "ACT": "H",
            "AGT": "D",
            "ACG": "V",
        }

        self.iupac = {
            "A": ["A"],
            "C": ["C"],
            "G": ["G"],
            "T": ["T"],
            "S": ["C", "G"],
            "R": ["A", "G"],
            "W": ["A", "T"],
            "Y": ["C", "T"],
            "K": ["G", "T"],
            "M": ["A", "C"],
            "B": ["C", "G", "T"],
            "H": ["A", "C", "T"],
            "D": ["A", "G", "T"],
            "V": ["A", "C", "G"],
            "N": ["A", "C", "G", "T"],
        }
        self.iupac_pwm = {
            "A": [1, 0, 0, 0],
            "C": [0, 1, 0, 0],
            "G": [0, 0, 1, 0],
            "T": [0, 0, 0, 1],
            "S": [0, 0.5, 0.5, 0],
            "R": [0.5, 0, 0.5, 0],
            "W": [0.5, 0, 0, 0.5],
            "Y": [0, 0.5, 0, 0.5],
            "K": [0, 0, 0.5, 0.5],
            "M": [0.5, 0.5, 0, 0],
            "B": [0, 0.33, 0.33, 0.33],
            "H": [0.33, 0.33, 0, 0.33],
            "D": [0.33, 0, 0.33, 0.33],
            "V": [0.33, 0.33, 0.33, 0],
            "N": [0.25, 0.25, 0.25, 0.25],
        }

    def __getitem__(self, x):
        """
        Take slice of a motif and return as new Motif instance.

        Returns
        -------
        motif : Motif instance
            Slice of the motif.
        """
        m = Motif()
        if self.pwm:
            m.pwm = self.pwm[x]
        if self.pfm:
            m.pfm = self.pfm[x]
        if self.seqs:
            m.seqs = [seq[x] for seq in self.seqs]
        if self.consensus:
            m.consensus = self.consensus[x]
        return m

    def __len__(self):
        """
        Return the motif length.

        Returns
        -------
        len : int
            Motif length.
        """
        return len(self.to_consensus())

    def __repr__(self):
        return "{}_{}".format(self.id, self.to_consensus())

    def information_content(self):
        """Return the total information content of the motif.

        Return
        ------
        ic : float
            Motif information content.
        """
        ic = 0
        for row in self.pwm:
            ic += 2.0 + np.sum(
                [row[x] * log(row[x]) / log(2) for x in range(4) if row[x] > 0]
            )
        return ic

    def pwm_min_score(self):
        """Return the minimum PWM score.

        Returns
        -------
        score : float
            Minimum PWM score.
        """
        if self.min_score is None:
            score = 0
            for row in self.pwm:
                score += log(min(row) / 0.25 + 0.01)
            self.min_score = score

        return self.min_score

    def pwm_max_score(self):
        """Return the maximum PWM score.

        Returns
        -------
        score : float
            Maximum PWM score.
        """
        if self.max_score is None:
            score = 0
            for row in self.pwm:
                score += log(max(row) / 0.25 + 0.01)
            self.max_score = score

        return self.max_score

    def score_kmer(self, kmer):
        """Calculate the log-odds score for a specific k-mer.

        Parameters
        ----------
        kmer : str
            String representing a kmer. Should be the same length as the motif.

        Returns
        -------
        score : float
            Log-odd score.
        """
        if len(kmer) != len(self.pwm):
            raise Exception("incorrect k-mer length")

        score = 0.0
        d = {"A": 0, "C": 1, "G": 2, "T": 3}
        for nuc, row in zip(kmer.upper(), self.pwm):
            score += log(row[d[nuc]] / 0.25 + 0.01)

        return score

    def pfm_to_pwm(self, pfm, pseudo=0.001):
        """Convert PFM with counts to a PFM with fractions.

        Parameters
        ----------
        pfm : list
            2-dimensional list with counts.
        pseudo : float
            Pseudocount used in conversion.

        Returns
        -------
        pwm : list
            2-dimensional list with fractions.
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
        m += "\n".join(["\t".join(["%s" % x for x in row]) for row in self.pwm])
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
        m : Motif instance
            New Motif instance with the reverse complement of the input motif.
        """
        m = Motif()
        m.pfm = [row[::-1] for row in self.pfm[::-1]]
        m.pwm = [row[::-1] for row in self.pwm[::-1]]
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
        pwm = self.pwm[:]
        while len(pwm) > 0 and self.ic_pos(pwm[0]) < edge_ic_cutoff:
            pwm = pwm[1:]
            self.pwm = self.pwm[1:]
            self.pfm = self.pfm[1:]
        while len(pwm) > 0 and self.ic_pos(pwm[-1]) < edge_ic_cutoff:
            pwm = pwm[:-1]
            self.pwm = self.pwm[:-1]
            self.pfm = self.pfm[:-1]

        self.consensus = None
        self.min_score = None
        self.max_score = None
        self.wiggled_pwm = None

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

    def pwm_scan(self, fa, cutoff=0.9, nreport=50, scan_rc=True):
        """Scan sequences with this motif.

        Scan sequences from a FASTA object with this motif. Less efficient
        than using a Scanner object. By setting the cutoff to 0.0 and
        nreport to 1, the best match for every sequence will be returned.
        Only the position of the matches is returned.

        Parameters
        ----------
        fa : Fasta object
            Fasta object to scan.
        cutoff : float , optional
            Cutoff to use for motif scanning. This cutoff is not specifically
            optimized and the strictness will vary a lot with motif lengh.
        nreport : int , optional
            Maximum number of matches to report.
        scan_rc : bool , optional
            Scan the reverse complement. True by default.

        Returns
        -------
        matches : dict
            Dictionary with motif matches. Only the position of the matches is
            returned.
        """
        c = (
            self.pwm_min_score()
            + (self.pwm_max_score() - self.pwm_min_score()) * cutoff
        )
        pwm = self.pwm
        matches = {}
        for name, seq in fa.items():
            matches[name] = []
            result = pfmscan(seq.upper(), pwm, c, nreport, scan_rc)
            for _, pos, _ in result:
                matches[name].append(pos)
        return matches

    def pwm_scan_all(self, fa, cutoff=0.9, nreport=50, scan_rc=True):
        """Scan sequences with this motif.

        Scan sequences from a FASTA object with this motif. Less efficient
        than using a Scanner object. By setting the cutoff to 0.0 and
        nreport to 1, the best match for every sequence will be returned.
        The score, position and strand for every match is returned.

        Parameters
        ----------
        fa : Fasta object
            Fasta object to scan.
        cutoff : float , optional
            Cutoff to use for motif scanning. This cutoff is not specifically
            optimized and the strictness will vary a lot with motif lengh.
        nreport : int , optional
            Maximum number of matches to report.
        scan_rc : bool , optional
            Scan the reverse complement. True by default.

        Returns
        -------
        matches : dict
            Dictionary with motif matches. The score, position and strand for
            every match is returned.
        """
        c = (
            self.pwm_min_score()
            + (self.pwm_max_score() - self.pwm_min_score()) * cutoff
        )
        pwm = self.pwm
        matches = {}
        for name, seq in fa.items():
            matches[name] = []
            result = pfmscan(seq.upper(), pwm, c, nreport, scan_rc)
            for score, pos, strand in result:
                matches[name].append((pos, score, strand))
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

        pfm = self.pwm
        nucs = np.array(["A", "C", "G", "T"])

        neg_matrix = np.zeros((len(pfm), 4))

        pos_matrix = []
        nuc_pfm = []
        for row in pfm:
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
            nuc_pfm.append(nucs[idx])

        colors = {
            "A": (0.308, 0.709, 0.280),
            "C": (0.145, 0.362, 0.6),
            "G": (0.969, 0.702, 0.172),
            "T": (0.841, 0.158, 0.224),
        }

        # x_max = np.max([np.sum(row) for row in pos_matrix])
        # x_min = -np.min([np.sum(row) for row in neg_matrix])
        # neg_matrix = neg_matrix / x_min * x_max

        plt.figure(figsize=(len(pfm) * 0.3, height))
        for (sign, matrix) in [(1, pos_matrix), (-1, neg_matrix)]:
            minbottom = np.zeros(len(matrix))
            alpha = 1
            if sign == -1:
                minbottom = np.array([np.sum(row) for row in matrix])
                alpha = 0.5

            # Print the bars
            for i in range(0, len(pfm[0])):

                pheight = [abs(r[i]) for r in matrix]
                bottom = minbottom + [sum(np.abs(r[:i])) for r in matrix]

                c = [colors[r[i]] for r in nuc_pfm]
                plt.bar(
                    range(1, len(pfm) + 1),
                    width=width,
                    height=pheight,
                    bottom=bottom,
                    color=c,
                    alpha=alpha,
                )

            if letters:
                # Print the letters
                for i in range(len(pfm)):
                    for n in range(4):
                        x = i + 1
                        y = matrix[i][n] / 2 + sum(matrix[i][:n])
                        nuc = nuc_pfm[i][n]
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
        plt.xticks(range(1, len(pfm) + 1))
        plt.ylabel(ylabel)

        if fname:
            plt.savefig(fname, dpi=300)
            plt.close()
        else:
            return ax

    def pwm_scan_score(self, fa, cutoff=0, nreport=1, scan_rc=True):
        """Scan sequences with this motif.

        Scan sequences from a FASTA object with this motif. Less efficient
        than using a Scanner object. By setting the cutoff to 0.0 and
        nreport to 1, the best match for every sequence will be returned.
        Only the score of the matches is returned.

        Parameters
        ----------
        fa : Fasta object
            Fasta object to scan.
        cutoff : float , optional
            Cutoff to use for motif scanning. This cutoff is not specifically
            optimized and the strictness will vary a lot with motif lengh.
        nreport : int , optional
            Maximum number of matches to report.
        scan_rc : bool , optional
            Scan the reverse complement. True by default.

        Returns
        -------
        matches : dict
            Dictionary with motif matches. Only the score of the matches is
            returned.
        """
        c = (
            self.pwm_min_score()
            + (self.pwm_max_score() - self.pwm_min_score()) * cutoff
        )
        pwm = self.pwm
        matches = {}
        for name, seq in fa.items():
            matches[name] = []
            result = pfmscan(seq.upper(), pwm, c, nreport, scan_rc)
            for score, _, _ in result:
                matches[name].append(score)
        return matches

    def pwm_scan_to_gff(
        self, fa, gfffile, cutoff=0.9, nreport=50, scan_rc=True, append=False
    ):
        """Scan sequences with this motif and save to a GFF file.

        Scan sequences from a FASTA object with this motif. Less efficient
        than using a Scanner object. By setting the cutoff to 0.0 and
        nreport to 1, the best match for every sequence will be returned.
        The output is save to a file in GFF format.

        Parameters
        ----------
        fa : Fasta object
            Fasta object to scan.
        gfffile : str
            Filename of GFF output file.
        cutoff : float , optional
            Cutoff to use for motif scanning. This cutoff is not specifically
            optimized and the strictness will vary a lot with motif lengh.
        nreport : int , optional
            Maximum number of matches to report.
        scan_rc : bool , optional
            Scan the reverse complement. True by default.
        append : bool , optional
            Append to GFF file instead of overwriting it. False by default.
        """
        if append:
            out = open(gfffile, "a")
        else:
            out = open(gfffile, "w")

        c = (
            self.pwm_min_score()
            + (self.pwm_max_score() - self.pwm_min_score()) * cutoff
        )
        pwm = self.pwm

        strandmap = {-1: "-", "-1": "-", "-": "-", "1": "+", 1: "+", "+": "+"}
        gff_line = (
            "{}\tpfmscan\tmisc_feature\t{}\t{}\t{:.3f}\t{}\t.\t"
            'motif_name "{}" ; motif_instance "{}"\n'
        )
        for name, seq in fa.items():
            result = pfmscan(seq.upper(), pwm, c, nreport, scan_rc)
            for score, pos, strand in result:
                out.write(
                    gff_line.format(
                        name,
                        pos,
                        pos + len(pwm),
                        score,
                        strandmap[strand],
                        self.id,
                        seq[pos : pos + len(pwm)],
                    )
                )
        out.close()

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

    def pcc(self, pwm1, pwm2, pos):
        # xxCATGYT
        # GGCTTGYx
        # pos = -2
        pwm1 = pwm1[:]
        pwm2 = pwm2[:]

        na = []
        if pos > 0:
            na = pwm1[:pos]
            pwm1 = pwm1[pos:]
        elif pos < 0:
            na = pwm2[:-pos]
            pwm2 = pwm2[-pos:]

        if len(pwm1) > len(pwm2):
            na += pwm1[len(pwm2) :]
            pwm1 = pwm1[: len(pwm2)]
        elif len(pwm2) > len(pwm1):
            na += pwm2[len(pwm1) :]
            pwm2 = pwm2[: len(pwm1)]

        # Aligned parts of the motif
        score = 0
        for a, b in zip(pwm1, pwm2):
            score += self.pcc_pos(a, b)

        return score

    def ic(self, pwm1, pwm2, pos, bg=None, bg_factor=1):
        if bg is None:
            bg = [0.25, 0.25, 0.25, 0.25]

        # xxCATGYT
        # GGCTTGYx
        # pos = -2
        pwm1 = pwm1[:]
        pwm2 = pwm2[:]

        na = []
        if pos > 0:
            na = pwm1[:pos]
            pwm1 = pwm1[pos:]
        elif pos < 0:
            na = pwm2[:-pos]
            pwm2 = pwm2[-pos:]

        if len(pwm1) > len(pwm2):
            na += pwm1[len(pwm2) :]
            pwm1 = pwm1[: len(pwm2)]
        elif len(pwm2) > len(pwm1):
            na += pwm2[len(pwm1) :]
            pwm2 = pwm2[: len(pwm1)]

        # print "COMPARE"
        # print Motif(pwm1).to_consensus()
        # print Motif(pwm2).to_consensus()

        # Aligned parts of the motif
        score = 0
        for a, b in zip(pwm1, pwm2):
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

    def other_ic(self, pwm1, pwm2, pos, bg=None, bg_factor=0.8):
        if bg is None:
            bg = [0.25, 0.25, 0.25, 0.25]

        # xxCATGYT
        # GGCTTGYx
        # pos = -2
        pwm1 = pwm1[:]
        pwm2 = pwm2[:]

        na = []
        if pos > 0:
            na = pwm1[:pos]
            pwm1 = pwm1[pos:]
        elif pos < 0:
            na = pwm2[:-pos]
            pwm2 = pwm2[-pos:]

        if len(pwm1) > len(pwm2):
            na += pwm1[len(pwm2) :]
            pwm1 = pwm1[: len(pwm2)]
        elif len(pwm2) > len(pwm1):
            na += pwm2[len(pwm1) :]
            pwm2 = pwm2[: len(pwm1)]

        # Aligned parts of the motif
        score = 0
        for a, b in zip(pwm1, pwm2):
            score += self.other_ic_pos(a, b)

        for x in na:
            score += self.other_ic_pos(x, bg) * bg_factor

        return score

    def matrix_ic(self, pwm1, pwm2, bg=None):
        if bg is None:
            bg = [0.25, 0.25, 0.25, 0.25]

        # xxCATGYT
        # GGCTTGYx
        # pos = -2
        pwm1 = np.array(pwm1)
        pwm2 = np.array(pwm2)
        pwm2_rev = np.array([row[::-1] for row in pwm2[::-1]])
        bg = np.array(bg)

        a = pwm1 * np.log2(pwm1 / bg)
        b = pwm2 * np.log2(pwm2 / bg)

        b_rev = pwm2_rev * np.log2(pwm2_rev / bg)

        scores = []
        l1 = len(pwm1)
        l2 = len(pwm2)
        for pos in range(-(l2 - 1), l1):

            pwm1_start, pwm2_start = 0, 0
            pwm1_end, pwm2_end = l1, l2
            if pos > 0:
                pwm1_start = pos
                if l1 - pos > l2:
                    pwm1_end = l2 + pos
                elif l1 - pos < l2:
                    pwm2_end = l1 - pos
            elif pos < 0:
                pwm2_start = -pos
                if l2 + pos > l1:
                    pwm2_end = l1
                elif l2 + pos < l1:
                    pwm1_end = l2 + pos
            else:
                if l2 > l1:
                    pwm2_end = l1
                elif l2 < l1:
                    pwm1_end = l2

            score = np.sum(
                (np.sum(a[pwm1_start:pwm1_end], 1) + np.sum(b[pwm2_start:pwm2_end], 1))
                / 2
                - np.sum(np.abs(a[pwm1_start:pwm1_end] - b[pwm2_start:pwm2_end]), 1)
            )
            scores.append([score, pos, 1])

            score = np.sum(
                (
                    np.sum(a[pwm1_start:pwm1_end], 1)
                    + np.sum(b_rev[pwm2_start:pwm2_end], 1)
                )
                / 2
                - np.sum(np.abs(a[pwm1_start:pwm1_end] - b_rev[pwm2_start:pwm2_end]), 1)
            )
            scores.append([score, pos, -1])

        return sorted(scores, key=lambda x: x[0])[-1]

    def max_ic(self, other, revcomp=True, bg_factor=0.8):
        pwm1 = self.pwm
        pwm2 = other.pwm

        scores = []

        for i in range(-(len(pwm2) - 1), len(pwm1)):
            scores.append((self.other_ic(pwm1, pwm2, i, bg_factor=bg_factor), i, 1))

        if revcomp:
            rev_pwm2 = [row[::-1] for row in pwm2[::-1]]
            for i in range(-(len(pwm2) - 1), len(pwm1)):
                scores.append(
                    (self.other_ic(pwm1, rev_pwm2, i, bg_factor=bg_factor), i, -1)
                )

        return sorted(scores, key=lambda x: x[0])[-1]

    def max_pcc(self, other, revcomp=True):
        pwm1 = self.pwm
        pwm2 = other.pwm

        scores = []

        for i in range(-(len(pwm2) - 1), len(pwm1)):
            scores.append((self.pcc(pwm1, pwm2, i), i, 1))

        if revcomp:
            rev_pwm2 = [row[::-1] for row in pwm2[::-1]]
            for i in range(-(len(pwm2) - 1), len(pwm1)):
                scores.append((self.pcc(pwm1, rev_pwm2, i), i, -1))

        # print scores
        return sorted(scores, key=lambda x: x[0])[-1]

    def _format_jaspar(self, version=1, header=True):
        rows = np.array(self.pwm).transpose()
        rows = [" ".join([str(x) for x in row]) for row in rows]
        if version == 2:
            rows = ["{} [{} ]".format(n, row) for n, row in zip(self.nucs, rows)]

        str_out = "\n".join(rows)
        if header:
            str_out = "\n".join([self.id, str_out])

        return str_out

    def to_consensus(self, precision=4):
        if not self.consensus:
            consensus = ""
            for row in self.pwm:
                weights = sorted(zip(["A", "C", "G", "T"], row), key=lambda x: x[1])
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
            self.consensus = consensus

        return self.consensus

    def to_consensusv2(self):
        if self.consensus:
            return self.consensus
        else:
            consensus = ""
            for row in self.pwm:
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
            pfm = [[n * self.PSEUDO_PFM_COUNT for n in col] for col in self.pwm]
            return ">%s\n%s" % (
                self.id,
                "\n".join(["\t".join(["%s" % x for x in row]) for row in pfm]),
            )

    def _pwm_to_str(self, precision=4):
        """Return string representation of pwm.

        Parameters
        ----------
        precision : int, optional, default 4
            Floating-point precision.

        Returns
        -------
        pwm_string : str
        """
        if not self.pwm:
            return ""

        fmt = "{{:.{:d}f}}".format(precision)
        return "\n".join(["\t".join([fmt.format(p) for p in row]) for row in self.pwm])

    def hash(self):
        """Return hash of motif.

        This is an unique identifier of a motif, regardless of the id.

        Returns:
        hash : str
        """
        return xxhash.xxh64(self._pwm_to_str(3)).hexdigest()

    def to_pwm(self, precision=4, extra_str=""):
        """Return pwm as string.

        Parameters
        ----------
        precision : int, optional, default 4
            Floating-point precision.

        extra_str |: str, optional
            Extra text to include with motif id line.

        Returns
        -------
        motif_str : str
            Motif formatted in PWM format.
        """
        motif_id = self.id

        if extra_str:
            motif_id += "_%s" % extra_str

        if not self.pwm:
            self.pwm = [self.iupac_pwm[char] for char in self.consensus.upper()]

        return ">%s\n%s" % (motif_id, self._pwm_to_str(precision))

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

    def wiggle_pwm(self):
        if self.wiggled_pwm is None:
            self.wiggled_pwm = [
                np.array(row) + (np.random.random(4) / 1e6) for row in self.pwm
            ]
            self.wiggled_pwm = [list(row / np.sum(row)) for row in self.wiggled_pwm]

        return self.wiggled_pwm

    def format_factors(
        self, max_length=5, html=False, include_indirect=True, extra_str=", (...)"
    ):
        if html:
            fmt_d = "<span style='color:black'>{}</span>"
            fmt_i = "<span style='color:#666666'>{}</span>"
        else:
            fmt_d = fmt_i = "{}"

        if hasattr(self, "factor_info"):
            fcount = Counter([x.upper() for x in self.factor_info["Factor"]])
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
    pwm = {}
    for line in handle:
        line = line.strip()
        if len(line) == 0:
            continue

        if line.startswith(">"):
            motif_id = line[1:]
        if line[0] in "ACGT":
            m = p.search(line)
            try:
                nuc = m.group(1)
                counts = re.split(r"\s+", m.group(2).strip())
                pwm[nuc] = [float(x) for x in counts]
                if nuc == "T":
                    motif = Motif(np.array([pwm[n] for n in "ACGT"]).transpose())
                    motif.id = motif_id
                    motifs.append(motif)
            except Exception:
                raise ValueError("Can't parse line\n" + line)

    if motif_id and motifs[-1].id != motif_id:
        motif = Motif(np.array([pwm[n] for n in "ACGT"]).transpose())
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

            pwm = np.array(freqs).transpose()
            motif = Motif(pwm)
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
    pwm = []
    motif_id = ""
    metadata = None
    for line in handle.readlines():
        m = p_id.search(line.strip())
        if m:
            motif_id = m.group(2)
            motif_id = re.sub(r"\s+", "_", motif_id)
        elif line.startswith("//"):
            if len(pwm) != 0:
                motifs.append(Motif(pwm))
                motifs[-1].id = motif_id
                if metadata is not None:
                    motifs[-1].metadata = metadata
            pwm = []
        elif p.search(line):
            m = p.search(line)
            pwm.append([float(x) for x in m.group(1, 2, 3, 4)])
        elif line.startswith("CC"):
            # CC tax_group:vertebrates; pubmed_ids:7592839; uniprot_ids:P30561,P53762; data_type:SELEX
            metadata = line.replace("CC ", "").split(";")
            metadata = dict([x.strip().split(":") for x in metadata])
    # If there's only one matrix, and the format is not complete
    if len(pwm) != 0:
        motifs.append(Motif(pwm))

    return motifs


def motifs_to_meme(motifs):
    m = "MEME version 3.0\n\nALPHABET= ACGT\n\nstrands: + -\n\n"
    m += "Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n"

    for motif in motifs:
        m += motif.to_meme() + "\n"
    return m
