"""
Scanner class utility functions
"""
import os
import logging
import sys

from gimmemotifs.config import CACHE_DIR
from gimmemotifs.c_metrics import pwmscan  # noqa

logger = logging.getLogger("gimme.scanner")


def file_hash(fname):
    """very dirty, very quick method to identify a file."""
    name = os.path.splitext(fname)[0]
    name = os.path.basename(name)

    byte_size = os.path.getsize(fname)
    kilobyte_size = str(byte_size)[:-3]

    return hash(name + kilobyte_size)


def print_cluster_error_message():
    logger.error("Cache is corrupted.")
    logger.error(
        "This can happen when you try to run a GimmeMotifs tool in parallel on a cluster."
    )
    logger.error(f"To solve this, delete the GimmeMotifs cache directory: {CACHE_DIR}")
    logger.error("and then see here for a workaround:")
    logger.error(
        "https://gimmemotifs.readthedocs.io/en/master/faq.html#sqlite-error-when-running-on-a-cluster"
    )


def parse_threshold_values(motifs, cutoff):
    d = parse_cutoff(motifs, cutoff)
    threshold = {}
    for m in motifs:
        c = m.min_score + (m.max_score - m.min_score) * d[m.id]
        threshold[m.id] = c
    return threshold


def parse_cutoff(motifs, cutoff, default=0.95):
    """Provide either a file with one cutoff per motif* or a single cutoff.
    *gimme threshold can provide such a file

    Returns a dictionary with motif id as key and cutoff as value
    """
    cutoffs = {}
    if os.path.isfile(str(cutoff)):
        # cutoff is a table
        for i, line in enumerate(open(cutoff)):
            if line != "Motif\tScore\tCutoff\n":
                try:
                    split = line.strip().split("\t")
                    motif = split[0]  # assumption on position
                    cutoff = float(split[-1])  # assumption on position
                    cutoffs[motif] = cutoff
                except Exception:
                    logger.error(f"Error parsing cutoff file, line {line}")

        if len(cutoff) == 0:
            logger.error(f"Could not parse file {motifs}")
            sys.exit(1)

        for motif in motifs:
            if motif.id not in cutoffs:
                logger.error(
                    f"No cutoff found for {motif.id}, using default {default}\n"
                )
                cutoffs[motif.id] = default
    else:
        # cutoff is a value
        for motif in motifs:
            cutoffs[motif.id] = float(cutoff)

    return cutoffs


def scan_seq_mult(
    seqs, seq_gc_bins, motifs, nreport, scan_rc, motifs_meanstd=None, zscore=False
):
    """scan multiple sequences for motifs"""
    ret = []
    for seq, seq_gc_bin in zip(seqs, seq_gc_bins):
        result = scan_sequence(
            seq.upper(),
            seq_gc_bin,
            motifs,
            nreport,
            scan_rc,
            motifs_meanstd=motifs_meanstd,
            zscore=zscore,
        )
        ret.append(result)
    return ret


def scan_sequence(
    seq, seq_gc_bin, motifs, nreport, scan_rc, motifs_meanstd=None, zscore=False
):
    """scan single sequence for motifs"""
    ret = []
    for motif, cutoff in motifs:
        if cutoff is None:
            ret.append([])
        else:
            if zscore:
                m_mean, m_std = motifs_meanstd[seq_gc_bin][motif.id]
                result = pwmscan(
                    seq, motif.logodds.tolist(), motif.min_score, nreport, scan_rc
                )
                result = [[(row[0] - m_mean) / m_std, row[1], row[2]] for row in result]
                result = [row for row in result if row[0] >= cutoff]
            else:
                result = pwmscan(seq, motif.logodds.tolist(), cutoff, nreport, scan_rc)
            if cutoff <= motif.min_score and len(result) == 0:
                result = [[motif.min_score, 0, 1]] * nreport

            ret.append(result)

    return ret
