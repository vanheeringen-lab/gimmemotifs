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

    from ._comparison import ic, pcc, other_ic, matrix_ic, max_ic, max_pcc, ic_pos
    from ._plotting import plot_logo, plot_ensembl_logo, to_img
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

    def __init__(self, pfm=None, ppm=None, places=4):

        self._places = places

        if pfm is None:
            # PFM is not specified
            if ppm is not None and len(ppm) > 0:
                # PPM is specified
                self.ppm = ppm
            else:
                self.pfm = []
        else:
            if np.all(np.isclose(np.sum(pfm, 1), 1, atol=1e-3)) and ppm is None:
                # PFM is specified actually a PPM. We don't mind.
                self.ppm = pfm
            else:
                # PFM is really a PFM
                self.pfm = pfm
                if ppm is not None and len(ppm) > 0:
                    # And we got a PPM as well, so set it but don't update the PFM
                    self._set_ppm(ppm, update_pfm=False)

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

    def _set_ppm(self, mtx, update_pfm=True):
        if mtx is not None and len(mtx) > 0:
            if update_pfm:
                self._pfm = [[n * self.PSEUDO_PFM_COUNT for n in col] for col in mtx]
            self._ppm = [
                iteround.saferound([float(n) for n in x], self._places) for x in mtx
            ]
        else:
            self._ppm = []
            if update_pfm:
                self._pfm = []

        self._pfm = np.array(self._pfm)
        self._ppm = np.array(self._ppm)
        self._update_associated()

    def _update_associated(self):
        self._logodds = [
            [np.log(n / self.G + self.Z) for n in col] for col in self._ppm
        ]
        self._logodds = np.array(self._logodds)
        self._consensus = self.to_consensus(self.ppm)
        if len(self) > 0:
            self._max_score = self.logodds.max(1).sum()
            self._min_score = self.logodds.min(1).sum()
        else:
            self._max_score = 0
            self._min_score = 0

    def _set_pfm(self, mtx, update_ppm=True):
        if mtx is not None and len(mtx) > 0:
            self._pfm = [list(x) for x in mtx]
            self._set_ppm(self.pfm_to_ppm(mtx), update_pfm=False)

    @ppm.setter
    def ppm(self, mtx):
        self._set_ppm(mtx)

    @pfm.setter
    def pfm(self, mtx):
        self._set_pfm(mtx)

    # @pfm.setter
    # def pfm(self, mtx):
    #     if mtx is not None and len(mtx) > 0:
    #         if np.sum(mtx[0]) > 2:
    #             self._pfm = [list(x) for x in mtx]
    #             self._ppm = self.pfm_to_ppm(mtx)
    #             self._ppm = [iteround.saferound(x, self._places) for x in self._ppm]
    #         else:
    #             self._ppm = [iteround.saferound(list(x), self._places) for x in mtx]
    #             self._pfm = [[n * self.PSEUDO_PFM_COUNT for n in col] for col in mtx]
    #     else:
    #         self._ppm = []
    #         self._pfm = []

    #     self._logodds = [
    #         [np.log(n / self.G + self.Z) for n in col] for col in self._ppm
    #     ]

    #     self._pfm = np.array(self._pfm)
    #     self._ppm = np.array(self._ppm)
    #     self._logodds = np.array(self._logodds)
    #     self._consensus = self.to_consensus(self.ppm)
    #     if len(self) > 0:
    #         self._max_score = self.logodds.max(1).sum()
    #         self._min_score = self.logodds.min(1).sum()
    #     else:
    #         self._max_score = 0
    #         self._min_score = 0

    @property
    def logodds(self):
        return self._logodds

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

    @property
    def information_content(self):
        """Return the total information content of the motif.

        Returns
        -------
        float
            Motif information content.
        """
        # Ignore divide-by-zero errors in log2.
        # We only use the return from log2 if the input was positive,
        # so this error should not impact the calculation.
        with np.errstate(divide="ignore"):
            log_ppm = np.log2(self.ppm)
        return ((self.ppm * np.where(self.ppm > 0, log_ppm, 0)).sum(1) + 2).sum()

    @property
    def hash(self):
        """Return hash of motif.

        This is an unique identifier of a motif, regardless of the id.

        Returns
        -------
        str
            Hash of motif.
        """
        if not hasattr(self, "_hash") or self._hash is None:
            self._hash = xxhash.xxh64(self._ppm_to_str(3)).hexdigest()

        return self._hash

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

    def __lshift__(self, other):
        """Return the motif shifted left.

        Parameters
        ----------
        other : integer
            Number of positions to shift.

        Returns
        -------
        Motif instance
            Shifted Motif.
        """
        total = self.pfm[0].sum()
        bg = [[total / 4 for _ in range(4)]] * other
        m = Motif(
            pfm=np.vstack((self.pfm, bg)),
            ppm=np.vstack((self.ppm, [[0.25] * 4] * other)),
        )
        m.id = m.consensus
        return m

    def __rshift__(self, other):
        """Return the motif shifted right.

        Parameters
        ----------
        other : integer
            Number of positions to shift.

        Returns
        -------
        Motif instance
            Shifted Motif.
        """
        total = self.pfm[0].sum()
        bg = [[total / 4 for _ in range(4)]] * other
        m = Motif(
            pfm=np.vstack((bg, self.pfm)),
            ppm=np.vstack(([[0.25] * 4] * other, self.ppm)),
        )
        m.id = m.consensus
        return m

    def __invert__(self):
        """Return the reverse complemented motif.

        Returns
        -------
        Motif instance
            New Motif instance with the reverse complement of the input motif.
        """
        return self.rc()

    def __add__(self, other):
        """Return the average of two motifs.

        Combine this motif with another motif and return the average as a new
        Motif object. This method works on the pfm, which means that motifs with
        higher frequences will be weighed more heavily.

        Parameters
        ----------
        other : Motif object
            Other Motif object.

        Returns
        -------
        motif : motif object
            New Motif object containing average motif.
        """
        diff = len(self) - len(other)
        if diff > 0:
            new = Motif((other << diff).pfm + self.pfm)
        elif diff < 0:
            new = Motif((self << diff).pf + other.pfm)
        else:
            new = Motif(self.pfm + other.pfm)

        new.id = new.consensus
        return new

    def __mul__(self, other):
        """Return motif with pfm multiplied by an value.

        Parameters
        ----------
        other : int

        Returns
        -------
        motif : motif object
            New Motif object containing average motif.
        """
        return Motif(pfm=self.pfm * other, ppm=self.ppm)

    def __and__(self, other):
        """Return the average of two motifs.

        Combine this motif with another motif and return the average as a new
        Motif object. This method works on the ppm, which means that motifs will
        be weighed equally.

        Parameters
        ----------
        other : Motif object
            Other Motif object.

        Returns
        -------
        motif : motif object
            New Motif object containing average motif.
        """
        diff = len(self) - len(other)
        if diff > 0:
            new = Motif(ppm=((other << diff).ppm + self.ppm))
        elif diff < 0:
            new = Motif(ppm=((self << diff).ppm + other.ppm))
        else:
            new = Motif(ppm=(self.ppm + other.ppm))

        new.id = new.consensus
        return new

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

    def randomize(self):
        """Create a new motif with shuffled positions.

        Shuffle the positions of this motif and return a new Motif instance.

        Returns
        -------
        m : Motif instance
            Motif instance with shuffled positions.
        """
        warn(
            "The randomize() method is deprecated and will be removed in future release."
            "Please use the shuffle() method instead.",
            DeprecationWarning,
        )
        return self.shuffle()

    def shuffle(self):
        """Create a new motif with shuffled positions.

        Shuffle the positions of this motif and return a new Motif instance.

        Returns
        -------
        Motif instance
            Motif instance with shuffled positions.
        """
        m = Motif(pfm=np.random.permutation(self.pfm))
        m.id = f"{self.id}_shuffled"
        return m

    # To be checked, documented, tested and refactored after here

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

        # TODO: don't convert back to list, but make sure this works for arrays
        pfm1 = [x.tolist() for x in self.pfm]
        pfm2 = [x.tolist() for x in other.pfm]

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
        if len(self.pfm) > 0:
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
        if self.ppm is None or len(self.ppm) == 0:
            return ""

        fmt = "{{:.{:d}f}}".format(precision)
        return "\n".join(["\t".join([fmt.format(p) for p in row]) for row in self.ppm])

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

        if self.ppm is None or len(self.ppm) == 0:
            self.ppm = [self.iupac_ppm[char] for char in self.consensus.upper()]

        return ">%s\n%s" % (motif_id, self._ppm_to_str(precision))

    def to_pwm(self, precision=4, extra_str=""):
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
        warn(
            "The to_pwm() function is deprecated and will be removed in a future release."
            "Please use the to_ppm() function instead.",
            DeprecationWarning,
        )
        return self.to_ppm(precision=precision, extra_str=extra_str)

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

    def sample(self, n_seqs, rng=None):
        """Sample n_seqs random sequences from a motif. The sequences
        follow the distribution of the motif ppm.

        Parameters
        ----------
        n_seqs : int
            number of sequences to sample
        rng : np.random.Generator
            random number generator, optional

        Returns
        -------
        sequences : List[str]
            A list of all the samples sequences
        """
        if rng is None:
            rng = np.random.default_rng()
        cumsum = np.expand_dims(np.cumsum(self.ppm, axis=1), axis=1)
        unif = rng.uniform(0, 1, size=(self.ppm.shape[0], n_seqs, 1))
        # strings is aggregated into a (motif_len, n_seqs) array that will have
        # 0 for "A", 1 for "C", etc. (as in NUCS)
        strings = 4 - (unif < cumsum).sum(axis=2, dtype="i1")
        # the values are now replaced with the ASCII ordinals to be
        # reinterpreted as characters
        for i, n in enumerate(NUCS):
            strings[strings == i] = ord(n)
        # allocate a Python string for each of the results
        # TODO: will the following slicing be faster if the in-memory alignment
        # is correct (i.e., FORTRAN vs C)?
        return [strings[:, i].tobytes().decode("UTF-8") for i in range(n_seqs)]


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

    m = Motif(pfm=pfm)
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
    motifs : list or dict
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
            metadata = [x.strip().split(":") for x in metadata]
            metadata = dict([x for x in metadata if len(x) == 2])
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
