"""
Motif class
"""
import re
from collections import Counter
from warnings import warn

import iteround
import numpy as np
import xxhash

from gimmemotifs.config import DIRECT_NAME, INDIRECT_NAME, MotifConfig

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

    from ._comparison import (  # noqa
        ic,
        ic_pos,
        matrix_ic,
        max_ic,
        max_pcc,
        other_ic,
        pcc,
    )
    from ._plotting import plot_ensembl_logo, plot_logo  # noqa
    from ._scanning import pwm_scan_score, pwm_scan_to_gff, scan, scan_all  # noqa

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
        return f"{self.id}_{self.to_consensus()}"

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
        m += f"NA {self.id}\n"
        m += "P0\tA\tC\tG\tT\n"
        for i, row in enumerate(self.pfm):
            r = "\t".join([str(int(x)) for x in row])
            m += f"{i}\t{r}\n"
        m += "//"
        return m

    def to_transfac(self):
        """Return motif formatted in TRANSFAC format

        Returns
        -------
        m : str
            String of motif in TRANSFAC format.
        """
        m = f"DE\t{self.id}\tunknown\n"
        for i, (row, cons) in enumerate(zip(self.pfm, self.to_consensus())):
            row = "\t".join([str(int(x)) for x in row])
            m += f"{i}\t{row}\t{cons}\n"
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
        m = f"MOTIF {motif_id}\n"
        m += f"BL   MOTIF {motif_id} width=0 seqs=0\n"
        m += f"letter-probability matrix: alength= 4 w= {len(self)} nsites= {np.sum(self.pfm[0])} E= 0\n"
        m += "\n".join(["\t".join([str(x) for x in row]) for row in self.ppm])
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
            rows = [f"{n} [{row} ]" for n, row in zip(NUCS, rows)]

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
            pfm = self.pfm
        else:
            pfm = [[n * self.PSEUDO_PFM_COUNT for n in col] for col in self.ppm]
        rows = "\n".join(["\t".join([str(x) for x in row]) for row in pfm])
        return f">{self.id}\n{rows}"

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

        p = precision
        return "\n".join(["\t".join([f"{e:.{p}f}" for e in row]) for row in self.ppm])

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
            motif_id += f"_{extra_str}"

        if self.ppm is None or len(self.ppm) == 0:
            self.ppm = [self.iupac_ppm[char] for char in self.consensus.upper()]

        return f">{motif_id}\n{self._ppm_to_str(precision)}"

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
