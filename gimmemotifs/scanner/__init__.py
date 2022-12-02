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
import numpy as np
import os
import pandas as pd
import re
import sys

from gimmemotifs import __version__
from gimmemotifs.motif import read_motifs
from gimmemotifs.scanner.base import Scanner, FPR
from gimmemotifs.utils import as_fasta, pfmfile_location

logger = logging.getLogger("gimme.scanner")


def scan_to_best_match(
    fasta,
    pfmfile=None,
    genome=None,
    score=False,
    zscore=False,
    gc=False,
    ncpus=None,
    random_state=None,
    progress=None,
):
    """Scan a FASTA file for motifs.
    Return a dictionary with the best match per motif.

    Parameters
    ----------
    fasta : str
        Filename of a sequence file in FASTA format.

    pfmfile : str or list, optional
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
    if not isinstance(pfmfile, list):
        pfmfile = pfmfile_location(pfmfile)

    # Initialize scanner
    s = Scanner(ncpus=ncpus, random_state=random_state, progress=progress)
    s.set_motifs(pfmfile)
    if genome:
        s.set_genome(genome)
    if genome and zscore:
        s.set_background(gc=gc)

    result = {m: [] for m in s.motif_ids}
    if score:
        it = s.best_score(fasta, zscore=zscore)
    else:
        it = s.best_match(fasta, zscore=zscore)
    for scores in it:
        for motif, score in zip(s.motif_ids, scores):
            result[motif].append(score)

    # Close the pool and reclaim memory
    del s

    return result


def scan_regionfile_to_table(
    input_table,
    genome,
    scoring,
    pfmfile=None,
    zscore=True,
    gc=True,
    ncpus=None,
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
    if not isinstance(pfmfile, list):
        pfmfile = pfmfile_location(pfmfile)

    if input_table.endswith("feather"):
        df = pd.read_feather(input_table)
        regions = list(df.iloc[:, 0].values)
    else:
        df = pd.read_table(input_table, index_col=0, comment="#")
        regions = list(df.index)

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
    s.set_background(gc=gc, size=size)

    # create an empty dataframe with
    # regions from the input table and motifs from the pfmfile
    df = pd.DataFrame(index=regions, columns=s.motif_ids, data=0.0, dtype="float16")

    logger.debug("setting threshold")
    if scoring == "count":
        s.set_threshold(fpr=FPR)
        logger.info("Creating count table")
        for idx, row in enumerate(s.count(regions)):
            df.iloc[idx] = row
        df = df.astype(int)
    else:
        scoring = "z-score" if zscore else "logodds"
        gc = ", GC% corrected" if gc else ""
        logger.info(f"Creating score table ({scoring}{gc})")
        for idx, row in enumerate(
            s.best_score(
                regions,
                zscore=zscore,
                dtype="float16",
            )
        ):
            df.iloc[idx] = row
    logger.info("Done")

    return df


def scan_to_file(
    inputfile,
    pfmfile=None,
    filepath_or_buffer=None,
    nreport=1,
    fpr=None,
    cutoff=None,
    scan_rc=True,
    count_table=False,
    score_table=False,
    bed=False,
    bgfile=None,
    genome=None,
    ncpus=None,
    zscore=True,
    gc=True,
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

    count_table : bool, optional
        output motif counts in tabular format

    score_table : bool, optional
        output motif scores in tabular format

    bed : bool, optional
        outputs BED6 format, instead of GTF/GFF format (default).

    bgfile : str
        FASTA file to use as background sequences. Required if no genome is given.

    genome : str
        Genome name. Can be either the name of a FASTA-formatted file or a
        genomepy genome name. Required if no bgfile is given.

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
    """
    if not isinstance(pfmfile, list):
        pfmfile = pfmfile_location(pfmfile)

    should_close = False
    if filepath_or_buffer is None:
        # write to stdout
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
        fpr = FPR

    print("# GimmeMotifs version {}".format(__version__), file=fo)
    print("# Input: {}".format(inputfile), file=fo)
    print("# Motifs: {}".format(pfmfile), file=fo)
    if fpr and not score_table:
        if genome is not None:
            print("# FPR: {} ({})".format(fpr, genome), file=fo)
        elif bgfile:
            print("# FPR: {} ({})".format(fpr, bgfile), file=fo)
    if cutoff is not None:
        print("# Threshold: {}".format(cutoff), file=fo)
    if zscore and not count_table:
        if gc:
            print("# Scoring: GC frequency normalized z-score", file=fo)
        else:
            print("# Scoring: normalized z-score", file=fo)
    else:
        print("# Scoring: logodds score", file=fo)

    motifs = read_motifs(pfmfile)
    fa = as_fasta(inputfile, genome)

    # initialize scanner
    s = Scanner(ncpus=ncpus, random_state=random_state, progress=progress)
    s.set_motifs(pfmfile)
    if genome or bgfile:
        s.set_genome(genome)
        s.set_background(fasta=bgfile, size=fa.median_length(), gc=gc)
    elif not score_table:
        logger.error("a genome or background file is required")
        sys.exit(1)

    if count_table:
        s.set_threshold(fpr=fpr, threshold=cutoff)
        it = _scan_count_table(s, fa, motifs, nreport, scan_rc)
    elif score_table:
        it = _scan_score_table(s, fa, motifs, scan_rc, zscore)
    else:
        s.set_threshold(fpr=fpr, threshold=cutoff)
        it = _scan_normal(s, fa, motifs, nreport, scan_rc, bed, zscore)

    for line in it:
        print(line, file=fo)

    if should_close:
        fo.close()


def _scan_count_table(s, fa, motifs, nreport, scan_rc):
    # get iterator
    result_it = s.count(fa, nreport, scan_rc)
    # header
    yield "\t{}".format("\t".join([m.id for m in motifs]))
    # counts table
    for i, counts in enumerate(result_it):
        yield "{}\t{}".format(fa.ids[i], "\t".join([str(x) for x in counts]))


def _scan_score_table(s, fa, motifs, scan_rc, zscore):
    # get iterator
    result_it = s.best_score(fa, scan_rc, zscore)
    # header
    yield "\t{}".format("\t".join([m.id for m in motifs]))
    # score table
    for i, scores in enumerate(result_it):
        yield "{}\t{}".format(fa.ids[i], "\t".join(["{:4f}".format(x) for x in scores]))


def _scan_normal(s, fa, motifs, nreport, scan_rc, bed, zscore):
    # get iterator
    result_it = s.scan(fa, nreport, scan_rc, zscore)
    # BED/GTF file
    for i, result in enumerate(result_it):
        seq_id = fa.ids[i]
        seq = fa[seq_id]
        for motif, matches in zip(motifs, result):
            for (score, pos, strand) in matches:
                yield _format_line(seq, seq_id, motif, score, pos, strand, bed)


def _format_line(seq, seq_id, motif, score, pos, strand, bed):
    loc = seq_id  # example seq_id: chr1:100-120
    start = pos  # relative to the seq_id start
    end = pos + len(motif)
    strand = {-1: "-", 1: "+"}[strand]

    if bed:
        # try to convert to chromosome + absolute coordinates
        m = re.search(r"([^\s:]+):(\d+)-(\d+)", seq_id)
        if m:
            loc = m.group(1)
            start += int(m.group(2))
            end += int(m.group(2))

        return "{}\t{}\t{}\t{}\t{}\t{}".format(
            loc,
            start,
            end,
            motif.id,
            score,
            strand,
        )

    return "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
        loc,
        "pfmscan",
        "misc_feature",
        start + 1,  # GFF is 1-based
        end,
        score,
        strand,
        ".",
        f'motif_name "{motif.id}" ; motif_instance "{seq[start:end]}"',
    )
