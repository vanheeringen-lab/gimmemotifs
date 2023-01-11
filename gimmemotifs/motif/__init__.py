"""
Motif core functions
"""
__all__ = [
    "gimme_motifs",
    "Motif",
    "read_motifs",
    "parse_motifs",
    "motif_from_consensus",
    "motif_from_align",
]

from gimmemotifs.motif.base import Motif
from gimmemotifs.motif.read import parse_motifs, read_motifs

from gimmemotifs.motif.denovo import gimme_motifs  # isort: skip (needs read_motif)


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


# def motifs_to_meme(motifs):
#     m = "MEME version 3.0\n\nALPHABET= ACGT\n\nstrands: + -\n\n"
#     m += "Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n"
#
#     for motif in motifs:
#         m += motif.to_meme() + "\n"
#     return m
