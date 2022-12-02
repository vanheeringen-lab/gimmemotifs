import os
import re
import sys
import time
import math
from collections import Counter
from functools import partial
from tempfile import mkdtemp, NamedTemporaryFile
import logging
import multiprocessing as mp

# "hidden" features, in development
try:
    import MOODS.tools
    import MOODS.parsers
    import MOODS.scan
except ImportError:
    pass

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
from gimmemotifs.utils import parse_cutoff, as_fasta, file_hash, rc


try:
    import copy_reg
    import types

    def _pickle_method(m):
        if m.im_self is None:
            return getattr, (m.im_class, m.im_func.func_name)
        else:
            return getattr, (m.im_self, m.im_func.func_name)

    copy_reg.pickle(types.MethodType, _pickle_method)
except Exception:
    pass

# only used when using cache, should not be a requirement
try:
    from dogpile.cache import make_region
except ImportError:
    pass

logger = logging.getLogger("gimme.scanner")
config = MotifConfig()
FPR = 0.01
lock = mp.Lock()


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
    zscore=True,
    gc=True,
    ncpus=None,
    random_state=None,
    progress=None,
):
    """Scan regions in input table with motifs.

    Parameters
    ----------
    input_table : str
        Filename of input table. Can be either a tab-separated table or a
        feather file.

    genome : str
        Genome name. Can be either the name of a FASTA-formatted file or a
        genomepy genome name.

    scoring : str
        "count" or "score"

    pfmfile : str, optional
        Specify a PFM file for scanning (or a list of Motif instances).

    zscore : bool, optional
        Use z-score normalized motif scores.

    gc : bool, optional
        Use GC% bins for z-score.

    ncpus : int, optional
        If defined this specifies the number of cores to use.

    random_state : numpy.random.RandomState object, optional
        make predictions deterministic (where possible).

    progress : bool or None, optional
        provide progress bars for long computations.

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

    logger.debug("reading table")
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
        s.set_threshold(threshold=0.0)
        msg = "Creating score table"
        if zscore:
            msg += " (z-score"
            if gc:
                msg += ", GC%"
            msg += ")"
        else:
            msg += " (logodds)"
        logger.info(msg)
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


def _scan_table(
    s,
    inputfile,
    fa,
    motifs,
    cutoff,
    bgfile,
    nreport,
    scan_rc,
    pvalue,
    moods,
):
    # header
    yield "\t{}".format("\t".join([m.id for m in motifs]))
    table = True
    if moods:
        result_it = scan_it_moods(
            inputfile, motifs, cutoff, bgfile, nreport, scan_rc, pvalue, table
        )
        for seq_id, counts in result_it:
            yield "{}\t{}".format(seq_id, "\t".join([str(x) for x in counts]))
    else:
        # get iterator
        result_it = s.count(fa, nreport, scan_rc)
        # counts table
        for i, counts in enumerate(result_it):
            yield "{}\t{}".format(fa.ids[i], "\t".join([str(x) for x in counts]))


def _scan_score_table(s, fa, motifs, scan_rc, zscore=False, gcnorm=False):

    s.set_threshold(threshold=0.0, gc=gcnorm)
    # get iterator
    result_it = s.best_score(fa, scan_rc, zscore=zscore)
    # header
    yield "\t{}".format("\t".join([m.id for m in motifs]))
    # score table
    for i, scores in enumerate(result_it):
        yield "{}\t{}".format(fa.ids[i], "\t".join(["{:4f}".format(x) for x in scores]))


def _scan_normal(
    s,
    inputfile,
    fa,
    motifs,
    cutoff,
    bgfile,
    nreport,
    scan_rc,
    pvalue,
    moods,
    bed,
    zscore,
    gcnorm,
):

    table = False
    if moods:
        result_it = scan_it_moods(
            inputfile, motifs, cutoff, bgfile, nreport, scan_rc, pvalue, table
        )
        for motif, d in result_it:
            for seq_id, matches in d.items():
                for pos, score, strand in matches:
                    yield _format_line(
                        fa[seq_id], seq_id, motif, score, pos, strand, bed=bed
                    )
    else:
        result_it = s.scan(fa, nreport, scan_rc, zscore)
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
    moods=False,
    pvalue=None,
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
    s.set_genome(genome=genome)

    if genome or bgfile:
        s.set_background(fasta=bgfile, size=fa.median_length(), gc=gcnorm)

        if not score_table:
            s.set_threshold(fpr=fpr, threshold=cutoff)

    if table:
        it = _scan_table(
            s,
            inputfile,
            fa,
            motifs,
            cutoff,
            bgfile,
            nreport,
            scan_rc,
            pvalue,
            moods,
        )
    elif score_table:
        it = _scan_score_table(s, fa, motifs, scan_rc, zscore, gcnorm)
    else:
        it = _scan_normal(
            s,
            inputfile,
            fa,
            motifs,
            cutoff,
            bgfile,
            nreport,
            scan_rc,
            pvalue,
            moods,
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
    moods=False,
    pvalue=False,
    bgfile=None,
    genome=None,
    ncpus=None,
    zscore=True,
    gcnorm=True,
    random_state=None,
    progress=None,
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
        moods=moods,
        pvalue=pvalue,
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
    """Scan a FASTA file with motifs.

    Scan a FASTA file and return a dictionary with the best match per motif.

    Parameters
    ----------
    fname : str
        Filename of a sequence file in FASTA format.

    motifs : list or str
        List of motif instances or a pwm file.

    Returns
    -------
    result : dict
        Dictionary with motif scanning results.
    """
    # Initialize scanner
    s = Scanner(ncpus=ncpus, random_state=random_state, progress=progress)
    s.set_motifs(motifs)
    s.set_genome(genome)
    s.set_background(gc=gc)
    s.set_threshold(threshold=0.0)

    if isinstance(motifs, str):
        motifs = read_motifs(motifs)

    result = dict([(m.id, []) for m in motifs])
    if score:
        it = s.best_score(fname, zscore=zscore)
    else:
        it = s.best_match(fname, zscore=zscore)
    for scores in it:
        for motif, score in zip(motifs, scores):
            result[motif.id].append(score)

    # Close the pool and reclaim memory
    del s

    return result


def parse_threshold_values(motifs, cutoff):
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


def scan_fa_with_motif_moods(
    fo, motifs, matrices, bg, thresholds, nreport, scan_rc=True
):

    scanner = MOODS.scan.Scanner(7)
    scanner.set_motifs(matrices, bg, thresholds)

    ret = []
    for name, seq in fo.items():
        length = len(seq)

        scan_seq = seq.upper()
        if scan_rc:
            scan_seq = "".join((scan_seq, "N" * 50, rc(scan_seq)))
        results = scanner.scan_max_hits(scan_seq, nreport)
        for motif, result in zip(motifs, results):
            matches = []
            for match in result:
                strand = 1
                pos = match.pos
                if scan_rc:
                    if pos > length:
                        pos = length - (pos - length - 50) - len(motif)
                        strand = -1
                matches.append((pos, match.score, strand))
            ret.append((motif, {name: matches}))

    return ret


def scan_fa_with_motif_moods_count(
    fo, motifs, matrices, bg, thresholds, nreport, scan_rc=True
):
    scanner = MOODS.scan.Scanner(7)
    scanner.set_motifs(matrices, bg, thresholds)

    ret = []
    for name, seq in fo.items():
        scan_seq = seq.upper()
        if scan_rc:
            scan_seq = "".join((scan_seq, "N" * 50, rc(scan_seq)))
        results = scanner.counts_max_hits(scan_seq, nreport)
        ret.append((name, results))

    return ret


def calc_threshold_moods(m, c):
    m_min = MOODS.tools.min_score(m)
    m_max = MOODS.tools.max_score(m)

    return m_min + (m_max - m_min) * c


def scan_it_moods(
    infile, motifs, cutoff, bgfile, nreport=1, scan_rc=True, pvalue=None, count=False
):
    tmpdir = mkdtemp()
    matrices = []
    pseudocount = 1e-3
    # sys.stderr.write("bgfile: {}\n".format(bgfile))
    bg = MOODS.tools.bg_from_sequence_dna("".join(Fasta(bgfile).seqs), 1)

    for motif in motifs:
        pfmname = os.path.join(tmpdir, "{}.pfm".format(motif.id))
        with open(pfmname, "w") as f:
            matrix = np.array(motif.ppm).transpose()
            for line in [" ".join([str(x) for x in row]) for row in matrix]:
                f.write("{}\n".format(line))

        matrices.append(MOODS.parsers.pfm_log_odds(pfmname, bg, pseudocount))

    thresholds = []
    if pvalue is not None:
        thresholds = [
            MOODS.tools.threshold_from_p(m, bg, float(pvalue)) for m in matrices
        ]
        # sys.stderr.write("{}\n".format(thresholds))
    else:
        thresholds = [calc_threshold_moods(m, float(cutoff)) for m in matrices]

    scanner = MOODS.scan.Scanner(7)
    scanner.set_motifs(matrices, bg, thresholds)

    config = MotifConfig()
    ncpus = int(config.get_default_params()["ncpus"])
    fa = Fasta(infile)
    chunk = 500
    if (len(fa) / chunk) < ncpus:
        chunk = len(fa) / (ncpus + 1)

    jobs = []
    func = scan_fa_with_motif_moods
    if count:
        func = scan_fa_with_motif_moods_count

    pool = mp.Pool()
    for i in range(0, len(fa), chunk):
        jobs.append(
            pool.apply_async(
                func,
                (fa[i : i + chunk], motifs, matrices, bg, thresholds, nreport, scan_rc),
            )
        )

    for job in jobs:
        for ret in job.get():
            yield ret


class Scanner(object):
    """
    scan sequences with motifs
    """

    genome = None
    background = None
    background_hash = None
    threshold = None
    meanstd = {}
    gc_bins = []
    fpr = None
    motifs = None
    motif_ids = []
    _seed = None

    def __init__(self, ncpus=None, random_state=None, progress=None):
        self.config = MotifConfig()
        self.random_state = random_state
        self.progress = progress

        if ncpus is None:
            self.ncpus = int(MotifConfig().get_default_params()["ncpus"])
        else:
            self.ncpus = ncpus

        if self.ncpus > 1:
            self.pool = mp.Pool(processes=self.ncpus)

        self.use_cache = False
        if self.config.get_default_params().get("use_cache", False):
            self._init_cache()

    def __del__(self):
        # Close the pool because of memory leak
        try:
            if hasattr(self, "pool"):
                self.pool.close()
        except Exception:
            pass

    @property
    def random_state(self):
        return self._random_state

    @random_state.setter
    def random_state(self, random_state):
        self._random_state = random_state
        if random_state is not None:
            # can be used to set random.seed()
            self._seed = random_state.get_state()[1][0]

    @property
    def progress(self):
        """Whether to show progress bars for slow operations."""
        return self._progress

    @progress.setter
    def progress(self, progress):
        self._progress = progress
        if progress is None:
            self._disable = None
        else:
            self._disable = not progress

    def _init_cache(self):
        try:
            self.cache = make_region().configure(
                "dogpile.cache.pylibmc",
                expiration_time=3600,
                arguments={"url": ["127.0.0.1"], "binary": True},
            )
            self.use_cache = True
        except Exception as e:
            sys.stderr.write("failed to initialize cache\n")
            sys.stderr.write("{}\n".format(e))

    def set_motifs(self, motifs):
        try:
            # Check if motifs is a list of Motif instances
            motifs[0].to_ppm()
            tmp = NamedTemporaryFile(mode="w", delete=False)
            for m in motifs:
                tmp.write("{}\n".format(m.to_ppm()))
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

    def _threshold_from_seqs(self, motifs, seqs):
        scan_motifs = [(m, m.min_score) for m in motifs]
        seq_gc_bins = [self.get_seq_bin(seq) for seq in seqs]

        # progress bar
        pbar = tqdm(
            desc="Determining FPR-based threshold",
            unit=" sequences",
            total=len(seqs),
            disable=self._disable,  # can be silenced
        )

        table = []
        for gc_bin, result in zip(
            seq_gc_bins, self._scan_sequences_with_motif(scan_motifs, seqs, 1, True)
        ):
            table.append([gc_bin] + [row[0][0] for row in result])
            pbar.update(1)

        df = pd.DataFrame(table, columns=["gc_bin"] + [m.id for m in motifs])
        df.set_index("gc_bin", inplace=True)
        return df

    def set_meanstd(self):
        if not self.background:
            raise ValueError("please run set_background() first")

        if not self.motifs:
            raise ValueError("please run set_motifs() first")

        self.meanstd = {}
        motifs = read_motifs(self.motifs)

        lock.acquire()
        try:
            with Cache(CACHE_DIR) as cache:

                # for each bin, load the meanstd of each motif
                # or mark it as missing
                scan_bins = {}
                for gc_bin in self.gc_bins:
                    if gc_bin not in self.meanstd:
                        self.meanstd[gc_bin] = {}

                    for motif in motifs:
                        k = "e{}|{}|{}".format(motif.hash, self.background_hash, gc_bin)
                        results = cache.get(k)

                        if results is None:
                            # missing values
                            scan_bins.setdefault(gc_bin, []).append(motif)
                        else:
                            self.meanstd[gc_bin][motif.id] = results

                # generate missing values for each motif per gc_bin
                total_scans = sum(len(v) for v in scan_bins.values())
                if total_scans > 0:
                    pbar = tqdm(
                        desc="Determining mean stdev",
                        unit=" motifs",
                        total=total_scans,
                        disable=self._disable,  # can be silenced
                    )

                    for gc_bin, bin_motifs in scan_bins.items():
                        bin_seqs = [
                            s
                            for i, s in self.background.items()
                            if i.endswith(" " + gc_bin)
                        ]
                        for motif, mean, std in self._meanstd_from_seqs(
                            bin_motifs, bin_seqs
                        ):
                            k = "e{}|{}|{}".format(
                                motif.hash, self.background_hash, gc_bin
                            )
                            cache.set(k, [mean, std])
                            self.meanstd[gc_bin][motif.id] = mean, std
                            pbar.update(1)

                # Prevent std of 0
                # This should only happen in testing
                for motif in motifs:
                    stds = np.array(
                        [self.meanstd[gcbin][motif.id][1] for gcbin in self.gc_bins]
                    )
                    idx = stds == 0
                    if True in idx:
                        std = np.mean(stds[~idx])
                        for gcbin in np.array(self.gc_bins)[idx]:
                            k = "e{}|{}|{}".format(
                                motif.hash, self.background_hash, gcbin
                            )
                            mean = self.meanstd[gcbin][motif.id][0]
                            cache.set(k, [mean, std])
                            self.meanstd[gcbin][motif.id] = mean, std

        except sqlite3.DatabaseError:
            print_cluster_error_message()
            sys.exit(1)
        lock.release()

    def set_background(
        self,
        fasta=None,
        size=200,
        nseq=None,
        gc=False,
        gc_bins=None,
    ):
        """Set the background to use for FPR and z-score calculations.

        Background can be specified either as a genome name or as the
        name of a FASTA file.

        Parameters
        ----------
        fasta : str, optional
            Path to FASTA file to use as background sequences.
            If unspecified, background sequences are sampled from the genome.

        size : int, optional
            Size of genomic sequences to retrieve. The default
            is 200.

        nseq : int, optional
            Number of genomic sequences to retrieve.

        gc : bool, optional
            distribute background sequences over GC percentage bins.

        gc_bins : list, optional
            a list bins containing desired GC percentages
            (e.g. ["0.40-0.50", "0.50-0.50"]).

        Sets
        ----
        self.gc_bins
        self.background
        self.background_hash
        """
        if self.background:
            return

        if fasta:
            if not os.path.exists(fasta):
                raise FileNotFoundError(f"Background file {fasta} does not exist!")

            self.background = Fasta(fasta)
            self.background_hash = file_hash(fasta)
            logger.debug("using background fasta file")
            return

        if self.genome is None:
            raise ValueError(
                "please run set_genome() first, or provide a FASTA file as background."
            )

        size = int(size)
        if gc_bins is None:
            if gc:
                gc_bins = [
                    (0.0, 0.2),  # bigger
                    (0.2, 0.25),
                    (0.25, 0.3),
                    (0.3, 0.35),
                    (0.35, 0.4),
                    (0.4, 0.45),
                    (0.45, 0.5),
                    (0.5, 0.55),
                    (0.55, 0.6),
                    (0.6, 0.65),
                    (0.65, 0.7),
                    (0.7, 0.75),
                    (0.75, 0.8),
                    (0.8, 1.00),  # bigger
                ]
            else:
                gc_bins = [(0, 1)]
        self.gc_bins = ["{:.2f}-{:.2f}".format(*b) for b in gc_bins]
        if nseq is None:
            nseq = max(10_000, len(gc_bins) * 1000)

        lock.acquire()
        try:
            with Cache(CACHE_DIR) as cache:
                self.background_hash = "d{}:{}:{}:{}".format(
                    file_hash(self.genome), size, gc, str(gc_bins)
                )
                fa = cache.get(self.background_hash)
                if fa is not None:
                    logger.debug("using cached background sequences")

                else:
                    logger.info("generating background sequences")
                    if gc:
                        with NamedTemporaryFile() as tmp:
                            gc_bin_bedfile(
                                tmp.name,
                                self.genome,
                                nseq,
                                size,
                                gc_bins,
                                self.random_state,
                            )
                            fa = as_fasta(tmp.name, genome=self.genome)
                    else:
                        fa = RandomGenomicFasta(self.genome, size, nseq, self._seed)

                    cache.set(self.background_hash, fa)
        except sqlite3.DatabaseError:
            print_cluster_error_message()
            sys.exit(1)
        lock.release()

        self.background = fa

    def set_threshold(self, fpr=None, threshold=None):
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

        Sets
        ----
        self.threshold
        self.fpr (if fpr or self.genome is given)
        """
        if not self.background:
            raise ValueError("please run set_background() first")

        if not self.motifs:
            raise ValueError("please run set_motifs() first")

        if threshold and fpr:
            raise ValueError("Need either fpr or threshold.")

        if threshold is None and fpr is None:
            if self.genome:
                fpr = FPR
                logger.info(f"Using default FPR of {fpr}")
            else:
                threshold = 0.95
                logger.info(
                    f"Genome not specified, using default threshold of {threshold}."
                )
                logger.info("This is likely not ideal.")

        motifs = read_motifs(self.motifs)

        if threshold is not None:
            d = parse_threshold_values(motifs, threshold)
            self.threshold = pd.DataFrame(d, index=[0])
            self.threshold = self.threshold.join(
                pd.DataFrame(
                    self.gc_bins, index=[0] * len(self.gc_bins), columns=["gc_bin"]
                )
            )
            self.threshold.set_index("gc_bin", inplace=True)
            return

        if fpr is not None:
            fpr = float(fpr)
            if not (0.0 < fpr < 1.0):
                raise ValueError("Parameter fpr should be between 0 and 1")
            self.fpr = fpr

        lock.acquire()
        try:
            with Cache(CACHE_DIR) as cache:
                scan_motifs = []
                self.threshold = None
                for motif in motifs:
                    k = "{}|{}|{}|{}".format(
                        motif.hash,
                        self.background_hash,
                        fpr,
                        ",".join(sorted(self.gc_bins)),
                    )
                    vals = cache.get(k)
                    if vals is None:
                        scan_motifs.append(motif)
                    else:
                        if self.threshold is None:
                            self.threshold = vals.to_frame()
                        else:
                            self.threshold[motif.id] = vals

                if len(scan_motifs) > 0:
                    seqs = self.background.seqs
                    df = self._threshold_from_seqs(scan_motifs, seqs)
                    if self.threshold is None:
                        self.threshold = df
                    else:
                        self.threshold = pd.concat((self.threshold, df), axis=1)
                    for motif in scan_motifs:
                        k = "{}|{}|{}|{}".format(
                            motif.hash,
                            self.background_hash,
                            fpr,
                            ",".join(sorted(self.gc_bins)),
                        )
                        cache.set(k, df[motif.id])
        except sqlite3.DatabaseError:
            print_cluster_error_message()
            sys.exit(1)
        lock.release()

    def set_genome(self, genome, genomes_dir=None):
        """
        Set the genome to be used for:
            - converting regions to sequences
            - background for MOODS

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
        genome = Genome(genome, genomes_dir, build_index=False).filename

        self.genome = genome

    def count(self, seqs, nreport=100, scan_rc=True):
        """
        count the number of matches above the cutoff
        returns an iterator of lists containing integer counts
        """
        for matches in self.scan(seqs, nreport, scan_rc):
            counts = [len(m) for m in matches]
            yield counts

    def total_count(self, seqs, nreport=100, scan_rc=True):
        """
        count the number of matches above the cutoff
        returns an iterator of lists containing integer counts
        """

        count_table = [counts for counts in self.count(seqs, nreport, scan_rc)]
        return np.sum(np.array(count_table), 0)

    def best_score(
        self,
        seqs,
        scan_rc=True,
        zscore=False,
        dtype=None,
    ):
        """
        give the score of the best match of each motif in each sequence
        returns an iterator of lists containing floats
        """
        self.set_threshold(threshold=0.0)
        # use numpy's default dtype if None is given
        dtype = float if dtype is None else dtype
        for matches in self.scan(seqs, 1, scan_rc, zscore=zscore):
            scores = np.array(
                [sorted(m, key=lambda x: x[0])[0][0] for m in matches if len(m) > 0],
                dtype=dtype,
            )
            yield scores

    def best_match(self, seqs, scan_rc=True, zscore=False):
        """
        give the best match of each motif in each sequence
        returns an iterator of nested lists containing tuples:
        (score, position, strand)
        """
        self.set_threshold(threshold=0.0)
        for matches in self.scan(seqs, 1, scan_rc, zscore=zscore):
            yield [m[0] for m in matches]

    def get_seq_bin(self, seq):
        if not self.gc_bins:
            raise ValueError("please run set_background() first")

        if len(str(seq)) == 0:
            gc = 0.01
        else:
            useq = seq.upper()
            gc = round((useq.count("G") + useq.count("C")) / len(useq), 2)
        if gc == 0:
            gc = 0.01

        for gc_bin in self.gc_bins:
            b_start = float(gc_bin.split("-")[0])
            b_end = float(gc_bin.split("-")[1])
            if b_start < gc <= b_end:
                return gc_bin

        raise ValueError(f"Error determining seq: {seq}, bins: {self.gc_bins}")

    def get_motif_mean_std(self, gc_bin, motif):
        if not self.meanstd:
            raise ValueError("please run set_meanstd() first")

        if gc_bin not in self.meanstd:
            raise ValueError(
                f"gc_bin {gc_bin} not found. available bins: {list(self.meanstd.keys())}"
            )

        if motif not in self.meanstd[gc_bin]:
            raise ValueError(
                "Motif not found in this gc_bin. " "This should really not happen!"
            )

        return self.meanstd[gc_bin][motif]

    def scan(
        self,
        seqs,
        nreport=100,
        scan_rc=True,
        zscore=False,
    ):
        """
        Scan a set of regions or sequences.
        """
        seqs = as_fasta(seqs, genome=self.genome)
        if zscore and self.meanstd is None:
            self.set_meanstd()

        # progress bar
        pbar = tqdm(
            desc="Scanning",
            unit=" sequences",
            total=len(seqs),
            disable=self._disable,  # can be silenced
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
        if self.fpr is None:
            raise ValueError("please run set_threshold() first")

        # Simple case, only one threshold
        if np.all(self.threshold.nunique(axis=0) == 1):
            return self.threshold.iloc[0].to_dict()

        max_seqs = 20000
        if len(seqs) > max_seqs:
            random = np.random if self.random_state is None else self.random_state
            seqs = random.choice(seqs, size=max_seqs)

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

        nseqs = int(max_seqs / np.sum(list(gc_bin_count.values())))
        t = {}
        maxt = pd.Series([m.max_score for m in motifs], index=_threshold.columns)
        # We do this in a loop as the DataFrame will get too big to fit in memory
        # when the difference between the number of sequences per gc_bin is very
        # high.
        _threshold = _threshold.reset_index()
        idx = np.hstack(
            [
                _threshold[_threshold[_threshold.columns[0]] == gc_bin]
                .sample(nseqs * count, replace=True, random_state=self.random_state)
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
