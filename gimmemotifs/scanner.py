import os
import re
import sys
import time
import math
from collections import Counter
from functools import partial
from tempfile import NamedTemporaryFile
import logging
from multiprocessing import Lock, Pool

from genomepy import Genome
from diskcache import Cache
import numpy as np
from sklearn.preprocessing import scale
import pandas as pd
import sqlite3
from tqdm.auto import tqdm

from gimmemotifs import __version__
from gimmemotifs.background import RandomGenomicFasta, gc_bin_bedfile
from gimmemotifs.config import MotifConfig, CACHE_DIR
from gimmemotifs.fasta import Fasta
from gimmemotifs.c_metrics import pwmscan  # noqa
from gimmemotifs.motif import read_motifs
from gimmemotifs.utils import parse_cutoff, as_fasta, file_checksum

logger = logging.getLogger("gimme.scanner")
config = MotifConfig()
FPR = 0.01
LOCK = Lock()


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
    progress=True,
):
    """Scan regions in input table with motifs.

    Parameters
    ----------
    input_table : str
        Filename of input table. Can be either a text-separated tab file or a
        feather file.

    genome : str
        Genome name. Can be either the name of a FASTA-formatted file or a
        genomepy genome name.

    scoring : str
        "count" or "score"

    pfmfile : str, optional
        Specify a PFM file for scanning (or a list of Motif instances).

    ncpus : int, optional
        If defined this specifies the number of cores to use.

    Returns
    -------
    table : pandas.DataFrame
        DataFrame with motif ids as column names and regions as index. Values
        are either counts or scores depending on the 'scoring' parameter.s
    """
    config = MotifConfig()

    if pfmfile is None:
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
        check_regions = np.random.choice(regions, size=1000, replace=False)
    else:
        check_regions = regions

    size = int(
        np.median([len(seq) for seq in as_fasta(check_regions, genome=genome).seqs])
    )
    s = Scanner(ncpus=ncpus)
    s.set_motifs(pfmfile)
    s.set_genome(genome)
    s.set_background(genome=genome, gc=gc, size=size)

    scores = []
    if scoring == "count":
        logger.info("setting threshold")
        s.set_threshold(fpr=FPR)
        logger.info("creating count table")
        for row in s.count(regions, progress=progress):
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
        for row in s.best_score(regions, zscore=zscore, gc=gc, progress=progress):
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
    progress=True,
):
    # header
    yield "\t{}".format("\t".join([m.id for m in motifs]))
    # get iterator
    result_it = s.count(fa, nreport, scan_rc, progress=progress)
    # counts table
    for i, counts in enumerate(result_it):
        yield "{}\t{}".format(fa.ids[i], "\t".join([str(x) for x in counts]))


def scan_score_table(s, fa, motifs, scan_rc, zscore=False, gcnorm=False, progress=True):

    s.set_threshold(threshold=0.0, gc=gcnorm)
    # get iterator
    result_it = s.best_score(fa, scan_rc, zscore=zscore, gc=gcnorm, progress=progress)
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
    progress=True,
):
    result_it = s.scan(fa, nreport, scan_rc, zscore, gc=gcnorm, progress=progress)
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
    progress=True,
):
    motifs = read_motifs(pfmfile)

    fa = as_fasta(inputfile, genome)

    # initialize scanner
    s = Scanner(ncpus=ncpus)
    s.set_motifs(pfmfile)

    if genome:
        s.set_genome(genome=genome)

    if genome:
        s.set_background(
            genome=genome, fname=bgfile, size=fa.median_length(), gc=gcnorm
        )
    if bgfile:
        s.set_background(genome=genome, fname=bgfile, size=fa.median_length())

    if not score_table:
        s.set_threshold(fpr=fpr, threshold=cutoff)

    if table:
        it = scan_table(
            s,
            fa,
            motifs,
            nreport,
            scan_rc,
            progress=progress,
        )
    elif score_table:
        it = scan_score_table(
            s, fa, motifs, scan_rc, zscore=zscore, gcnorm=gcnorm, progress=progress
        )
    else:
        it = scan_normal(
            s,
            fa,
            motifs,
            nreport,
            scan_rc,
            bed,
            zscore=zscore,
            gcnorm=gcnorm,
            progress=progress,
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
    progress=True,
):
    """Scan an inputfile with motifs."""
    should_close = False
    if filepath_or_buffer is None:
        fo = sys.stdout
    else:
        if hasattr(filepath_or_buffer, "write"):
            fo = filepath_or_buffer
        else:
            try:
                fo = open(os.path.expanduser(filepath_or_buffer), "w")
                should_close = True
            except Exception:
                logger.error(f"Could not open {filepath_or_buffer} for writing")
                sys.exit(1)

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
    progress=True,
):
    """Scan a FASTA file with motifs.

    Scan a FASTA file and return a dictionary with the best match per motif.

    Parameters
    ----------
    fname : str
        Filename of a sequence file in FASTA format.

    motifs : list
        List of motif instances.

    Returns
    -------
    result : dict
        Dictionary with motif scanning results.
    """
    # Initialize scanner
    s = Scanner(ncpus=ncpus)
    s.set_motifs(motifs)
    s.set_threshold(threshold=0.0)
    if genome:
        s.set_genome(genome)

    if isinstance(motifs, str):
        motifs = read_motifs(motifs)

    logger.debug(f"scanning {fname}...")
    result = dict([(m.id, []) for m in motifs])
    if score:
        it = s.best_score(fname, zscore=zscore, gc=gc, progress=progress)
    else:
        it = s.best_match(fname, zscore=zscore, gc=gc, progress=progress)
    for scores in it:
        for motif, score in zip(motifs, scores):
            result[motif.id].append(score)

    # Close the pool and reclaim memory
    del s

    return result


def parse_threshold_values(motif_file, cutoff):
    motifs = read_motifs(motif_file)
    d = parse_cutoff(motifs, cutoff)
    threshold = {}
    for m in motifs:
        c = m.min_score + (m.max_score - m.min_score) * d[m.id]
        threshold[m.id] = c
    return threshold


def scan_sequence(
    seq, seq_gc_bin, motifs, nreport, scan_rc, motifs_meanstd=None, zscore=False
):

    ret = []
    # scan for motifs
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


def scan_seq_mult(
    seqs, seq_gc_bins, motifs, nreport, scan_rc, motifs_meanstd=None, zscore=False
):
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


class Scanner(object):
    """
    scan sequences with motifs
    """

    def __init__(self, ncpus=None):
        self.config = MotifConfig()
        self._threshold = None
        self.genome = None
        self.background = None
        self.meanstd = {}
        self.gc_bins = [(0, 1)]
        self.fpr = None

        if ncpus is None:
            ncpus = MotifConfig().get_default_params()["ncpus"]
        self.ncpus = int(ncpus)

        if self.ncpus > 1:
            self.pool = Pool(processes=self.ncpus)

    def __del__(self):
        # Close the multiprocessing.Pool to release memory
        try:
            self.pool.close()
            self.pool.join()
        except (ValueError, AttributeError):
            pass

    def set_motifs(self, motifs):
        try:
            # Check if motifs is a list of Motif instances
            motifs[0].to_ppm()
            tmp = NamedTemporaryFile(mode="w", delete=False)
            for m in motifs:
                tmp.write(f"{m.to_ppm()}\n")
            tmp.close()
            motif_file = tmp.name
        except AttributeError:
            motif_file = motifs

        self.motifs = motif_file
        self.motif_ids = [m.id for m in read_motifs(motif_file)]

    def _meanstd_from_seqs(self, motifs, seqs):
        scan_motifs = [(m, m.min_score) for m in motifs]

        table = []
        for x in self._scan_sequences_with_motif(scan_motifs, seqs, 1, True):
            table.append([row[0][0] for row in x])

        for (motif, _), scores in zip(scan_motifs, np.array(table).transpose()):
            yield motif, np.mean(scores), np.std(scores)  # cutoff

    def _threshold_from_seqs(self, motifs, seqs, fpr):
        scan_motifs = [(m, m.min_score) for m in motifs]
        table = []
        seq_gc_bins = [self.get_seq_bin(seq) for seq in seqs]
        for gc_bin, result in zip(
            seq_gc_bins, self._scan_sequences_with_motif(scan_motifs, seqs, 1, True)
        ):
            table.append([gc_bin] + [row[0][0] for row in result])

        df = pd.DataFrame(table, columns=["gc_bin"] + [m.id for m in motifs])
        return df

    # the functions below are conceptually simpler (and faster?) but don't work due to
    # multiprocessing.Pool not being able to be used from another Pool.
    # Future improvement to fix this.
    # def _meanstd_from_seqs(self, motifs, seqs):
    #     # We need a new scanner instance, as this may be a subset of motifs
    #     s = Scanner()
    #     s.set_genome(self.genome)
    #     s.set_motifs(motifs)
    #     table = np.array(list(s.best_score(seqs, scan_rc=True, zscore=False, gc=False, progress=False)))

    #     for motif, scores in zip(motifs, table.transpose()):
    #         yield motif, np.mean(scores), np.std(scores)  # cutoff

    # def _threshold_from_seqs(self, motifs, seqs, fpr):
    #     # We need a new scanner instance, as this may be a subset of motifs
    #     s = Scanner()
    #     s.set_genome(self.genome)
    #     s.set_motifs(motifs)
    #     table = []
    #     seq_gc_bins = [self.get_seq_bin(seq) for seq in seqs]
    #     for gc_bin, row in zip(
    #         seq_gc_bins, list(s.best_score(seqs, scan_rc=True, zscore=False, gc=False, progress=False))
    #     ):
    #         print(row)
    #         table.append([gc_bin] + list(row)

    #     df = pd.DataFrame(table, columns=["gc_bin"] + [m.id for m in motifs])
    #     return df

    def set_meanstd(self, gc=False):
        if not self.background:
            self.set_background(gc=gc)

        self.meanstd = {}
        seqs = self.background.seqs
        if gc:
            seq_bins = [s.split(" ")[-1] for s in self.background.ids]
        else:
            seq_bins = ["0.00-1.00"] * len(seqs)
        if gc:
            bins = list(set(seq_bins))
        else:
            bins = ["0.00-1.00"]

        motifs = read_motifs(self.motifs)
        LOCK.acquire()

        try:
            with Cache(CACHE_DIR) as cache:
                scan_motifs = []
                for bin in bins:
                    if bin not in self.meanstd:
                        self.meanstd[bin] = {}
                    bin_seqs = [s for s, b in zip(seqs, seq_bins) if b == bin]

                    for motif in motifs:
                        k = "e{}|{}|{}".format(motif.hash, self.background_hash, bin)

                        results = cache.get(k)
                        if results is None:
                            scan_motifs.append(motif)
                        else:
                            self.meanstd[bin][motif.id] = results

                    if len(scan_motifs) > 0:
                        logger.debug("Determining mean and stddev for motifs.")
                        for motif, mean, std in self._meanstd_from_seqs(
                            scan_motifs, bin_seqs
                        ):
                            k = "e{}|{}|{}".format(
                                motif.hash, self.background_hash, bin
                            )
                            cache.set(k, [mean, std])
                            self.meanstd[bin][motif.id] = mean, std

                # Prevent std of 0
                # This should only happen in testing
                for motif in motifs:
                    stds = np.array(
                        [self.meanstd[gcbin][motif.id][1] for gcbin in bins]
                    )
                    idx = stds == 0
                    if True in idx:
                        std = np.mean(stds[~idx])
                        for gcbin in np.array(bins)[idx]:
                            k = "e{}|{}|{}".format(
                                motif.hash, self.background_hash, gcbin
                            )
                            mean = self.meanstd[gcbin][motif.id][0]
                            cache.set(k, [mean, std])
                            self.meanstd[gcbin][motif.id] = mean, std
        except sqlite3.DatabaseError:
            print_cluster_error_message()
            sys.exit(1)
        LOCK.release()

        for gc_bin in self.gc_bins:
            gc_bin = f"{gc_bin[0]:.2f}-{gc_bin[1]:.2f}"
            if gc_bin not in self.meanstd:
                valid_bins = []
                for b in self.gc_bins:
                    bstr = f"{b[0]:.2f}-{b[1]:.2f}"
                    if bstr in self.meanstd:
                        valid_bins.append(((b[0] + b[1]) / 2, bstr))

                v = float(gc_bin.split("-")[1])
                _, bstr = sorted(valid_bins, key=lambda x: abs(x[0] - v))[0]
                self.meanstd[gc_bin] = self.meanstd[bstr]

    def set_background(
        self, fname=None, genome=None, size=200, nseq=None, gc=False, gc_bins=None
    ):
        """Set the background to use for FPR and z-score calculations.

        Background can be specified either as a genome name or as the
        name of a FASTA file.

        Parameters
        ----------
        fname : str, optional
            Name of FASTA file to use as background.

        genome : str, optional
            Name of genome to use to retrieve random sequences.

        size : int, optional
            Size of genomic sequences to retrieve. The default
            is 200.

        nseq : int, optional
            Number of genomic sequences to retrieve.
        """
        if self.background:
            return

        size = int(size)

        if gc_bins is None:
            if gc:
                gc_bins = [(0.0, 0.2), (0.8, 1)]
                for b in np.arange(0.2, 0.799, 0.05):
                    gc_bins.append((b, b + 0.05))
            else:
                gc_bins = [(0, 1)]
        if nseq is None:
            nseq = max(10000, len(gc_bins) * 1000)

        if genome and fname:
            logger.debug("using genome for background")
            fname = None

        if fname:
            if not os.path.exists(fname):
                raise FileNotFoundError(f"Background file {fname} does not exist!")

            self.background = Fasta(fname)
            self.background_hash = file_checksum(fname)
            return

        if not genome:
            if self.genome:
                genome = self.genome
            else:
                raise ValueError("Need either genome or filename for background.")

        logger.debug(f"using background: genome {genome} with size {size}")
        LOCK.acquire()
        try:
            with Cache(CACHE_DIR) as cache:
                self.background_hash = "d{}:{}:{}:{}".format(
                    genome, int(size), gc, str(gc_bins)
                )
                c = cache.get(self.background_hash)
                if c:
                    fa, gc_bins = c
                else:
                    fa = None

                if not fa:
                    if gc:
                        with NamedTemporaryFile() as tmp:
                            logger.info(f"using {nseq} sequences")
                            gc_bin_bedfile(
                                tmp.name, genome, number=nseq, length=size, bins=gc_bins
                            )
                            fa = as_fasta(tmp.name, genome=genome)
                    else:
                        fa = RandomGenomicFasta(genome, size, nseq)
                    cache.set(self.background_hash, (fa, gc_bins))
        except sqlite3.DatabaseError:
            print_cluster_error_message()
            sys.exit(1)
        LOCK.release()

        self.background = fa
        if gc_bins:
            self.gc_bins = gc_bins

    @property
    def threshold(self):
        if self._threshold is None:
            self.set_threshold()
        return self._threshold

    def set_threshold(self, fpr=None, threshold=None, gc=False):
        """Set motif scanning threshold based on background sequences.

        Parameters
        ----------
        fpr : float, optional
            Desired FPR, between 0.0 and 1.0.

        threshold : float or str, optional
            Desired motif threshold, expressed as the fraction of the
            difference between minimum and maximum score of the PWM.
            Should either be a float between 0.0 and 1.0 or a filename
            with thresholds as created by 'gimme threshold'.

        """
        if threshold and fpr:
            raise ValueError("Need either fpr or threshold.")

        if threshold is None and fpr is None:
            if self.genome:
                fpr = 0.01
                logger.info(f"Using default FPR of {fpr}")
            else:
                threshold = 0.95
                logger.info(
                    f"Genome not specified, using default threshold of {threshold}."
                )
                logger.info("This is likely not ideal.")

        if fpr:
            fpr = float(fpr)
            if not (0.0 < fpr < 1.0):
                raise ValueError("Parameter fpr should be between 0 and 1")
            self.fpr = fpr

        if not self.motifs:
            raise ValueError("please run set_motifs() first")

        motifs = read_motifs(self.motifs)
        gc_bins = [f"{gc_bin[0]:.2f}-{gc_bin[1]:.2f}" for gc_bin in self.gc_bins]

        if threshold is not None:
            d = parse_threshold_values(self.motifs, threshold)
            self._threshold = pd.DataFrame(d, index=[0])
            self._threshold = self._threshold.join(
                pd.DataFrame(gc_bins, index=[0] * len(gc_bins), columns=["gc_bin"])
            )
            self._threshold = self._threshold.set_index("gc_bin")
            return

        if not self.background:
            try:
                self.set_background(gc=gc)
            except Exception:
                raise ValueError("please run set_background() first")

        seqs = self.background.seqs

        LOCK.acquire()
        try:
            with Cache(CACHE_DIR) as cache:
                scan_motifs = []
                self._threshold = None
                for motif in motifs:
                    k = "{}|{}|{:.4f}|{}".format(
                        motif.hash,
                        self.background_hash,
                        fpr,
                        ",".join(sorted(gc_bins)),
                    )
                    vals = cache.get(k)
                    if vals is None:
                        scan_motifs.append(motif)
                    else:
                        if self._threshold is None:
                            self._threshold = vals.to_frame()
                        else:
                            self._threshold[motif.id] = vals

                if len(scan_motifs) > 0:
                    logger.debug("determining FPR-based threshold")
                    df = self._threshold_from_seqs(scan_motifs, seqs, fpr).set_index(
                        "gc_bin"
                    )
                    if self._threshold is None:
                        self._threshold = df
                    else:
                        self._threshold = pd.concat((self._threshold, df), axis=1)
                    for motif in scan_motifs:
                        k = "{}|{}|{:.4f}|{}".format(
                            motif.hash,
                            self.background_hash,
                            fpr,
                            ",".join(sorted(gc_bins)),
                        )
                        cache.set(k, df[motif.id])
        except sqlite3.DatabaseError:
            print_cluster_error_message()
            sys.exit(1)
        LOCK.release()
        self.threshold_str = "{}_{}_{}_{}".format(
            fpr, threshold, self.background_hash, ",".join(sorted(gc_bins))
        )

    def set_genome(self, genome, genomes_dir=None):
        """
        Set the genome to converting regions to sequences

        Parameters
        ----------
        genome : str
            Path to the genome fasta or
            the genome name as found in the genomepy genomes directory.

        genomes_dir : str, optional
            Path to the genomepy genomes directory.
            Taken from the genomepy config if unspecified.
        """
        if not genome:
            return

        # raises error if checks fail
        genome = Genome(genome, genomes_dir, rebuild=False).filename

        self.genome = genome

    def count(self, seqs, nreport=100, scan_rc=True, progress=True):
        """
        count the number of matches above the cutoff
        returns an iterator of lists containing integer counts
        """
        for matches in self.scan(seqs, nreport, scan_rc, progress=progress):
            counts = [len(m) for m in matches]
            yield counts

    def total_count(self, seqs, nreport=100, scan_rc=True, progress=True):
        """
        count the number of matches above the cutoff
        returns an iterator of lists containing integer counts
        """

        count_table = [
            counts for counts in self.count(seqs, nreport, scan_rc, progress=progress)
        ]
        return np.sum(np.array(count_table), 0)

    def best_score(self, seqs, scan_rc=True, zscore=False, gc=False, progress=True):
        """
        give the score of the best match of each motif in each sequence
        returns an iterator of lists containing floats
        """
        self.set_threshold(threshold=0.0, gc=gc)
        for matches in self.scan(
            seqs, 1, scan_rc, zscore=zscore, gc=gc, progress=progress
        ):
            scores = np.array(
                [sorted(m, key=lambda x: x[0])[0][0] for m in matches if len(m) > 0]
            )
            yield scores

    def best_match(self, seqs, scan_rc=True, zscore=False, gc=False, progress=True):
        """
        give the best match of each motif in each sequence
        returns an iterator of nested lists containing tuples:
        (score, position, strand)
        """
        self.set_threshold(threshold=0.0)
        for matches in self.scan(
            seqs, 1, scan_rc, zscore=zscore, gc=gc, progress=progress
        ):
            yield [m[0] for m in matches]

    def get_seq_bin(self, seq):
        if len(str(seq)) == 0:
            gc = 0
        else:
            useq = seq.upper()
            gc = round((useq.count("G") + useq.count("C")) / len(useq), 2)
        if gc == 0:
            gc = 0.01
        for b_start, b_end in self.gc_bins:
            if gc > round(b_start, 2) and gc <= round(b_end, 2):
                return f"{b_start:.2f}-{b_end:.2f}"

        logger.error(f"Error determining seq: {seq}, bins: {str(self.gc_bins)}")
        raise ValueError()

    def get_motif_mean_std(self, gc_bin, motif):
        if gc_bin in self.meanstd:
            if motif not in self.meanstd[gc_bin]:
                raise ValueError("Motif mean and std not initialized")
        else:
            logger.warning(
                f"GC% {gc_bin} not present in genome, setting to closest GC% bin"
            )
            # Ideally this should not happen. If we get here, this means that
            # a sequence has a GC% that does not occur in the genome. This
            # can happen during tests with small genome, or when a sequence
            # is scanned that does not originate from the genome.
            # We use the closest GC% bin from the genome instead.
            valid_bins = []
            for b in self.gc_bins:
                bstr = f"{b[0]:.2f}-{b[1]:.2f}"
                if bstr in self.meanstd:
                    valid_bins.append(((b[0] + b[1]) / 2, bstr))

            v = float(gc_bin.split("-")[1])
            _, bstr = sorted(valid_bins, key=lambda x: abs(x[0] - v))[0]
            logger.warning(f"Using {bstr}")
            self.meanstd[gc_bin] = self.meanstd[bstr]
        return self.meanstd[gc_bin][motif]

    def scan(
        self, seqs, nreport=100, scan_rc=True, zscore=False, gc=False, progress=True
    ):
        """
        Scan a set of regions or sequences.
        """
        seqs = as_fasta(seqs, genome=self.genome)
        if zscore:
            if gc:
                if len(self.meanstd) <= 1:
                    self.set_meanstd(gc=gc)
            else:
                if len(self.meanstd) != 1:
                    self.set_meanstd(gc=gc)

        # progress bar
        pbar = tqdm(
            desc="scanning",
            unit=" sequences",
            total=len(seqs),
            disable=not progress,  # can be silenced
        )

        batch_size = 50000
        logger.debug("Scanning")
        for batch_idx in range(0, len(seqs), batch_size):
            it = self._scan_sequences(
                seqs.seqs[batch_idx : batch_idx + batch_size],
                nreport,
                scan_rc,
                zscore=zscore,
            )
            for result in it:
                yield result
                pbar.update(1)
        pbar.close()

    def get_gc_thresholds(self, seqs, motifs=None, zscore=False):
        MAX_SEQS = 20000

        # Simple case, only one threshold
        if np.all(self.threshold.nunique(axis=0) == 1):
            return self.threshold.iloc[0].to_dict()

        if len(seqs) > MAX_SEQS:
            seqs = np.random.choice(seqs, size=MAX_SEQS)

        if motifs is None:
            motifs = read_motifs(self.motifs)
        seq_gc_bins = [self.get_seq_bin(seq) for seq in seqs]

        gc_bin_count = Counter(seq_gc_bins)

        _threshold = self.threshold
        if zscore:
            grouped = _threshold.groupby(_threshold.index).apply(scale, axis=0)
            _threshold = pd.DataFrame(
                np.vstack(grouped.values),
                index=_threshold.index,
                columns=_threshold.columns,
            )

        nseqs = int(MAX_SEQS / np.sum(list(gc_bin_count.values())))
        t = {}
        maxt = pd.Series([m.max_score for m in motifs], index=_threshold.columns)
        # We do this in a loop as the DataFrame will get too big to fit in memory
        # when the difference between the number of sequences per gc_bin is very
        # high.
        _threshold = _threshold.reset_index()
        idx = np.hstack(
            [
                _threshold[_threshold[_threshold.columns[0]] == gc_bin]
                .sample(nseqs * count, replace=True, random_state=42)
                .index.values
                for gc_bin, count in gc_bin_count.items()
            ]
        )
        for motif in _threshold.columns[1:]:
            val = _threshold.loc[idx, motif].quantile(
                1 - self.fpr, interpolation="higher"
            )
            if val < maxt.loc[motif]:
                t[motif] = val
            else:
                t[motif] = None
        return t

    def _scan_sequences_with_motif(self, motifs, seqs, nreport, scan_rc):
        scan_func = partial(
            scan_seq_mult, motifs=motifs, nreport=nreport, scan_rc=scan_rc
        )
        for ret in self._scan_jobs(scan_func, seqs):
            yield ret[1]

    def _scan_sequences(self, seqs, nreport, scan_rc, zscore=False):
        thresholds = self.get_gc_thresholds(seqs, zscore=zscore)
        motifs = [(m, thresholds[m.id]) for m in read_motifs(self.motifs)]
        motifs_meanstd = None
        if zscore:
            motifs_meanstd = self.meanstd

        scan_func = partial(
            scan_seq_mult,
            motifs=motifs,
            nreport=nreport,
            scan_rc=scan_rc,
            motifs_meanstd=motifs_meanstd,
            zscore=zscore,
        )
        for _, ret in self._scan_jobs(scan_func, seqs):
            yield ret

    def _scan_jobs(self, scan_func, scan_seqs):
        if self.ncpus > 1:
            median_len = np.median([len(x) for x in scan_seqs])
            chunksize = 200000 // int(median_len)  # 1000 seqs for len 200
            # prepare for parallel processing
            k = 0
            max_queue_size = 2 * self.ncpus
            jobs = []

            # loop over each job/chunk, and keep adding them to the queue
            for i in range(math.ceil(len(scan_seqs) / chunksize)):
                batch_seqs = scan_seqs[i * chunksize : (i + 1) * chunksize]
                seq_gc_bins = [self.get_seq_bin(seq) for seq in batch_seqs]
                job = self.pool.apply_async(scan_func, (batch_seqs, seq_gc_bins))
                jobs.append(job)

                # if our queue is full, wait until oldest job finishes
                while (len(jobs) >= max_queue_size) and not jobs[0].ready():
                    time.sleep(0.05)

                # resolve oldest job if finished
                if jobs[0].ready():
                    for ret in jobs[0].get():
                        region = scan_seqs[k]
                        k += 1
                        yield region, ret
                    jobs = jobs[1:]

            # cleanup the last jobs that did not get resolved in the for loop
            while len(jobs) > 0:
                for ret in jobs[0].get():
                    region = scan_seqs[k]
                    k += 1
                    yield region, ret
                jobs = jobs[1:]
        else:
            # non-parallel job scanning
            batchsize = 1000
            for i in range((len(scan_seqs) - 1) // batchsize + 1):
                batch_seqs = scan_seqs[i * batchsize : (i + 1) * batchsize]
                seq_gc_bins = [self.get_seq_bin(seq) for seq in batch_seqs]

                for _j, ret in enumerate(scan_func(batch_seqs, seq_gc_bins)):
                    yield scan_seqs[i], ret
