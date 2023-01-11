"""
Scanner core functions
"""
__all__ = [
    "scan_regionfile_to_table",
    "scan_to_file",
    "scan_to_best_match",
    "Scanner",
]

import logging
import os
import re
import sys

import numpy as np
import pandas as pd

from gimmemotifs import __version__
from gimmemotifs.config import MotifConfig
from gimmemotifs.motif import read_motifs
from gimmemotifs.scanner.base import FPR, Scanner
from gimmemotifs.utils import as_fasta

logger = logging.getLogger("gimme.scanner")


def _format_line(
    seq, seq_id, motif, score, pos, strand, bed=False, seq_p=None, strandmap=None
):
    if seq_p is None:
        seq_p = re.compile(r"([^\s:]+):(\d+)-(\d+)")
    if strandmap is None:
        strandmap = {-1: "-", 1: "+"}
    if bed:
        m = seq_p.search(seq_id)
        if m:
            chrom = m.group(1)
            start = int(m.group(2))
            return "{}\t{}\t{}\t{}\t{}\t{}".format(
                chrom,
                start + pos,
                start + pos + len(motif),
                motif.id,
                score,
                strandmap[strand],
            )
        else:
            return "{}\t{}\t{}\t{}\t{}\t{}".format(
                seq_id, pos, pos + len(motif), motif.id, score, strandmap[strand]
            )
    else:
        return '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tmotif_name "{}" ; motif_instance "{}"'.format(
            seq_id,
            "pfmscan",
            "misc_feature",
            pos + 1,  # GFF is 1-based
            pos + len(motif),
            score,
            strandmap[strand],
            ".",
            motif.id,
            seq[pos : pos + len(motif)],
        )


def scan_regionfile_to_table(
    input_table,
    genome,
    scoring,
    pfmfile=None,
    ncpus=None,
    zscore=True,
    gc=True,
    random_state=None,
    progress=None,
):
    """Scan regions in input table for motifs.
    Return a dataframe with the motif count/score per region.

    Parameters
    ----------
    input_table : str
        Filename of a table with regions as first column. Accepts a feather file.

    genome : str
        Genome name. Can be a FASTA file or a genomepy genome name.

    scoring : str
        "count" or "score".
        "count" returns the occurrence of each motif (with an FPR threshold of 0.01).
        "score" returns the best match score of each motif.

    pfmfile : str, optional
        Specify a PFM file for scanning (or a list of Motif instances).

    zscore : bool, optional
        Use z-score normalized motif scores. Only used if scoring="score".

    gc : bool, optional
        Equally distribute GC percentages in background sequences.

    ncpus : int, optional
        If defined this specifies the number of cores to use.

    random_state : numpy.random.RandomState object, optional
        make predictions deterministic (where possible).

    progress : bool or None, optional
        provide progress bars for long computations.

    Returns
    -------
    table : pandas.DataFrame
        DataFrame with motifs as column names and regions as index. Values
        are either motif counts or best motif match scores per region,
        depending on the 'scoring' parameter.
    """
    if pfmfile is None:
        config = MotifConfig()
        pfmfile = config.get_default_params().get("motif_db", None)
        if pfmfile is not None:
            pfmfile = os.path.join(config.get_motif_dir(), pfmfile)

    if pfmfile is None:
        raise ValueError("no pfmfile given and no default database specified")

    logger.info("reading table")
    if input_table.endswith("feather"):
        df = pd.read_feather(input_table)
        idx = df.iloc[:, 0].values
    else:
        df = pd.read_table(input_table, index_col=0, comment="#")
        idx = df.index

    regions = list(idx)
    if len(regions) >= 1000:
        random = np.random if random_state is None else random_state
        check_regions = random.choice(regions, size=1000, replace=False)
    else:
        check_regions = regions

    size = int(
        np.median([len(seq) for seq in as_fasta(check_regions, genome=genome).seqs])
    )
    s = Scanner(ncpus=ncpus, random_state=random_state, progress=progress)
    s.set_motifs(pfmfile)
    s.set_genome(genome)
    s.set_background(genome=genome, gc=gc, size=size)

    scores = []
    if scoring == "count":
        logger.info("setting threshold")
        s.set_threshold(fpr=FPR)
        logger.info("creating count table")
        for row in s.count(regions):
            scores.append(row)
        logger.info("done")
    else:
        s.set_threshold(threshold=0.0)
        msg = "creating score table"
        if zscore:
            msg += " (z-score"
            if gc:
                msg += ", GC%"
            msg += ")"
        else:
            msg += " (logodds)"
        logger.info(msg)
        for row in s.best_score(regions, zscore=zscore, gc=gc):
            scores.append(row)
        logger.info("done")

    motif_names = s.motif_ids

    logger.info("creating dataframe")
    dtype = "float16"
    if scoring == "count":
        dtype = int
    df = pd.DataFrame(scores, index=idx, columns=motif_names, dtype=dtype)

    return df


def scan_table(
    s,
    fa,
    motifs,
    nreport,
    scan_rc,
):
    # header
    yield "\t{}".format("\t".join([m.id for m in motifs]))
    # get iterator
    result_it = s.count(fa, nreport, scan_rc)
    # counts table
    for i, counts in enumerate(result_it):
        yield "{}\t{}".format(fa.ids[i], "\t".join([str(x) for x in counts]))


def scan_score_table(s, fa, motifs, scan_rc, zscore=False, gcnorm=False):

    s.set_threshold(threshold=0.0, gc=gcnorm)
    # get iterator
    result_it = s.best_score(fa, scan_rc, zscore=zscore, gc=gcnorm)
    # header
    yield "\t{}".format("\t".join([m.id for m in motifs]))
    # score table
    for i, scores in enumerate(result_it):
        yield "{}\t{}".format(fa.ids[i], "\t".join(["{:4f}".format(x) for x in scores]))


def scan_normal(
    s,
    fa,
    motifs,
    nreport,
    scan_rc,
    bed,
    zscore,
    gcnorm,
):
    result_it = s.scan(fa, nreport, scan_rc, zscore, gc=gcnorm)
    for i, result in enumerate(result_it):
        seq_id = fa.ids[i]
        seq = fa[seq_id]
        for motif, matches in zip(motifs, result):
            for (score, pos, strand) in matches:
                yield _format_line(seq, seq_id, motif, score, pos, strand, bed=bed)


def command_scan(
    inputfile,
    pfmfile,
    nreport=1,
    fpr=0.01,
    cutoff=None,
    bed=False,
    scan_rc=True,
    table=False,
    score_table=False,
    bgfile=None,
    genome=None,
    ncpus=None,
    zscore=False,
    gcnorm=False,
    random_state=None,
    progress=None,
):
    motifs = read_motifs(pfmfile)

    fa = as_fasta(inputfile, genome)

    # initialize scanner
    s = Scanner(ncpus=ncpus, random_state=random_state, progress=progress)
    s.set_motifs(pfmfile)

    if genome:
        s.set_genome(genome)

    if genome:
        s.set_background(
            genome=genome, fname=bgfile, size=fa.median_length(), gc=gcnorm
        )
    if bgfile:
        s.set_background(genome=genome, fname=bgfile, size=fa.median_length())

    if not score_table:
        s.set_threshold(fpr=fpr, threshold=cutoff)

    if table:
        it = scan_table(s, fa, motifs, nreport, scan_rc)
    elif score_table:
        it = scan_score_table(s, fa, motifs, scan_rc, zscore, gcnorm)
    else:
        it = scan_normal(
            s,
            fa,
            motifs,
            nreport,
            scan_rc,
            bed,
            zscore,
            gcnorm,
        )

    for row in it:
        yield row


def scan_to_file(
    inputfile,
    pfmfile,
    filepath_or_buffer=None,
    nreport=1,
    fpr=0.01,
    cutoff=None,
    bed=False,
    scan_rc=True,
    table=False,
    score_table=False,
    bgfile=None,
    genome=None,
    ncpus=None,
    zscore=True,
    gcnorm=True,
    random_state=None,
    progress=None,
):
    """Scan file for motifs.

    Parameters
    ----------
    inputfile : str
        path to FASTA, BED or regions file.

    pfmfile : str or list, optional
        Specify a PFM file for scanning (or a list of Motif instances).

    filepath_or_buffer : Any, optional
        where to write the output. If unspecified, writes to stdout.

    nreport : int , optional
        Maximum number of matches to report.

    fpr : float, optional
        Desired false positive rate, between 0.0 and 1.0.

    cutoff : float , optional
        Cutoff to use for motif scanning. This cutoff is not specifically
        optimized and the strictness will vary a lot with motif length.

    scan_rc : bool , optional
        Scan the reverse complement. default: True.

    table : bool, optional
        output motif counts in tabular format

    score_table : bool, optional
        output motif scores in tabular format

    bed : bool, optional
        outputs BED6 format, instead of GTF/GFF format (default).

    bgfile : str, optional
        FASTA file to use as background sequences. Required if no genome is given.

    genome : str, optional
        Genome name. Can be either the name of a FASTA-formatted file or a
        genomepy genome name. Required if no bgfile is given.

    zscore : bool, optional
        Use z-score normalized motif scores.

    gcnorm : bool, optional
        Equally distribute GC percentages in background sequences.

    ncpus : int, optional
        If defined this specifies the number of cores to use.

    random_state : numpy.random.RandomState object, optional
        make predictions deterministic (where possible).

    progress : bool or None, optional
        provide progress bars for long computations.
    """
    should_close = False
    if filepath_or_buffer is None:
        #  write to stdout
        fo = sys.stdout
    elif hasattr(filepath_or_buffer, "write"):
        # write to buffer (open file or stdout)
        fo = filepath_or_buffer
    else:
        # write to file
        file_name = os.path.expanduser(filepath_or_buffer)
        os.makedirs(os.path.dirname(file_name), exist_ok=True)
        fo = open(file_name, "w")
        should_close = True

    if fpr is None and cutoff is None:
        fpr = 0.01

    print(f"# GimmeMotifs version {__version__}", file=fo)
    print(f"# Input: {inputfile}", file=fo)
    print(f"# Motifs: {pfmfile}", file=fo)
    if fpr and not score_table:
        if genome is not None:
            print(f"# FPR: {fpr} ({genome})", file=fo)
        elif bgfile:
            print(f"# FPR: {fpr} ({bgfile})", file=fo)
    if cutoff is not None:
        print(f"# Threshold: {cutoff}", file=fo)
    if zscore:
        if gcnorm:
            print("# Scoring: GC frequency normalized z-score", file=fo)
        else:
            print("# Scoring: normalized z-score", file=fo)
    else:
        print("# Scoring: logodds score", file=fo)

    for line in command_scan(
        inputfile,
        pfmfile,
        nreport=nreport,
        fpr=fpr,
        cutoff=cutoff,
        bed=bed,
        scan_rc=scan_rc,
        table=table,
        score_table=score_table,
        bgfile=bgfile,
        genome=genome,
        ncpus=ncpus,
        zscore=zscore,
        gcnorm=gcnorm,
        random_state=random_state,
        progress=progress,
    ):
        print(line, file=fo)

    if should_close:
        try:
            fo.close()
        except Exception:
            pass


def scan_to_best_match(
    fname,
    motifs,
    ncpus=None,
    genome=None,
    score=False,
    zscore=False,
    gc=False,
    random_state=None,
    progress=None,
):
    """Scan a FASTA file for motifs.
    Return a dictionary with the best match per motif.

    Parameters
    ----------
    fname : str or Fasta
        Filename of a sequence file in FASTA format.

    motifs : str or list, optional
        Specify a PFM file for scanning (or a list of Motif instances).

    genome : str
        Genome name. Can be either the name of a FASTA-formatted file or a
        genomepy genome name.

    score : bool, optional
        return the best score instead of the best match

    zscore : bool, optional
        Use z-score normalized motif scores.

    gc : bool, optional
        Equally distribute GC percentages in background sequences.

    ncpus : int, optional
        If defined this specifies the number of cores to use.

    random_state : numpy.random.RandomState object, optional
        make predictions deterministic (where possible).

    progress : bool or None, optional
        provide progress bars for long computations.

    Returns
    -------
    result : dict
        Dictionary with motif as key and best score/match as values.
    """
    # Initialize scanner
    s = Scanner(ncpus=ncpus, random_state=random_state, progress=progress)
    s.set_motifs(motifs)
    s.set_threshold(threshold=0.0)
    if genome:
        s.set_genome(genome)

    if isinstance(motifs, str):
        motifs = read_motifs(motifs)

    logger.debug(f"scanning {fname}...")
    result = dict([(m.id, []) for m in motifs])
    if score:
        it = s.best_score(fname, zscore=zscore, gc=gc)
    else:
        it = s.best_match(fname, zscore=zscore, gc=gc)
    for scores in it:
        for motif, score in zip(motifs, scores):
            result[motif.id].append(score)

    # Close the pool and reclaim memory
    del s

    return result
