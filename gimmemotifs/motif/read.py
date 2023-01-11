"""
Read Motif class instances from file
"""
import logging
import os
import re

import numpy as np
import pandas as pd

from gimmemotifs.config import DIRECT_NAME, INDIRECT_NAME
from gimmemotifs.motif.base import NUCS, Motif
from gimmemotifs.utils import pfmfile_location

logger = logging.getLogger("gimme.motif")


def parse_motifs(motifs):
    """Parse motifs in a variety of formats to return a list of motifs.

    Parameters
    ----------
    motifs : list or str
        Filename of motif, list of motifs or single Motif instance.

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
        raise ValueError(f"Unknown format {fmt}")

    # Remove everything after tab from motif identifiers
    for motif in motifs:
        motif.id = motif.id.split("\t")[0]

    motifs = _add_factors_from_handle(motifs, handle)

    return motifs


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
    with open(map_file) as f:
        for line in f:
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
                        m2f_indirect[motif] = m2f_indirect.get(motif, []) + [
                            factor_info[0]
                        ]
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
                msg = f"multiple motifs with same id: {motif_id}"
                logger.warning(msg)
                motif_id += f"_{seen_id[motif_id] - 1}"

        else:
            m = p.search(line)
            if m:
                fractions = [float(m.group(x)) for x in (1, 4, 7, 10)]
                pfm.append(fractions)
            else:
                msg = f"can't parse line {n+1}, ignoring:\n{line}"
                logger.warning(msg)

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
