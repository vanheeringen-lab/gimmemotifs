# Copyright (c) 2009-2021 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
"""Scanning functions for Motif class"""

from warnings import warn

from gimmemotifs.c_metrics import pfmscan, pwmscan


def pwm_scan(self, fa, cutoff=0.9, nreport=50, scan_rc=True):
    warn(
        "Method pwm_scan() is replaced by scan() and will be removed in the next release.",
        DeprecationWarning,
    )
    return self.scan(fa, cutoff=cutoff, nreport=nreport, scan_rc=scan_rc)


def pwm_scan_all(self, fa, cutoff=0.9, nreport=50, scan_rc=True):
    warn(
        "Method pwm_scan_all() is replaced by scan_all() and will be removed in the next release.",
        DeprecationWarning,
    )
    return self.scan_all(fa, cutoff=cutoff, nreport=nreport, scan_rc=scan_rc)


def scan(self, fa, cutoff=0.9, nreport=50, scan_rc=True):
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
    c = self.min_score + (self.max_score - self.min_score) * cutoff

    matches = {}
    for name, seq in fa.items():
        matches[name] = []
        result = pwmscan(seq.upper(), self.logodds.tolist(), c, nreport, scan_rc)
        for _, pos, _ in result:
            matches[name].append(pos)
    return matches


def scan_all(self, fa, cutoff=0.9, nreport=50, scan_rc=True):
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
    c = self.min_score + (self.max_score - self.min_score) * cutoff
    ppm = self.ppm
    matches = {}
    for name, seq in fa.items():
        matches[name] = []
        result = pfmscan(seq.upper(), ppm.tolist(), c, nreport, scan_rc)
        for score, pos, strand in result:
            matches[name].append((pos, score, strand))
    return matches


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
    c = self.min_score + (self.max_score - self.min_score) * cutoff
    ppm = self.ppm
    matches = {}
    for name, seq in fa.items():
        matches[name] = []
        result = pfmscan(seq.upper(), ppm.tolist(), c, nreport, scan_rc)
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

    c = self.min_score + (self.max_score - self.min_score) * cutoff
    ppm = self.ppm

    strandmap = {-1: "-", "-1": "-", "-": "-", "1": "+", 1: "+", "+": "+"}
    gff_line = (
        "{}\tpfmscan\tmisc_feature\t{}\t{}\t{:.3f}\t{}\t.\t"
        'motif_name "{}" ; motif_instance "{}"\n'
    )
    for name, seq in fa.items():
        result = pfmscan(seq.upper(), ppm.tolist(), c, nreport, scan_rc)
        for score, pos, strand in result:
            out.write(
                gff_line.format(
                    name,
                    pos,
                    pos + len(ppm),
                    score,
                    strandmap[strand],
                    self.id,
                    seq[pos : pos + len(ppm)],
                )
            )
    out.close()
