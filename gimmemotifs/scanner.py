import os
import re
import sys
from collections import Counter
from tempfile import NamedTemporaryFile
import logging
from multiprocessing.managers import SharedMemoryManager
from multiprocessing.shared_memory import SharedMemory
from concurrent.futures import ProcessPoolExecutor, as_completed

from genomepy import Genome
from diskcache import Cache
import numpy as np
from sklearn.preprocessing import scale
import pandas as pd
import sqlite3
import multiprocessing as mp
from tqdm.auto import tqdm

from gimmemotifs import __version__
from gimmemotifs.background import RandomGenomicFasta, gc_bin_bedfile
from gimmemotifs.config import MotifConfig, CACHE_DIR
from gimmemotifs.fasta import Fasta
from gimmemotifs.c_metrics import pwmscan
from gimmemotifs.motif import read_motifs
from gimmemotifs.utils import parse_cutoff, as_fasta, file_checksum

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
    import xxhash
except ImportError:
    pass

from collections import namedtuple

ScanData = namedtuple(
    "ScanData", "motif_array motif_array_shape motif_array_dtype gc_bin_list seq_list"
)


logger = logging.getLogger("gimme.scanner")
config = MotifConfig()
FPR = 0.01
lock = mp.Lock()


def get_seq_bin(seq, bins):
    if len(str(seq)) == 0:
        gc = 0
    else:
        useq = seq.upper()
        gc = round((useq.count("G") + useq.count("C")) / len(useq), 2)
    if gc == 0:
        gc = 0.01
    for b_start, b_end in bins:
        if gc > round(b_start, 2) and gc <= round(b_end, 2):
            return "{:.2f}-{:.2f}".format(b_start, b_end)

    logger.error("Error determining seq: {}, bins: {}".format(seq, str(bins)))
    raise ValueError()


def motifs_to_shareable_array(smm, motifs, thresholds=None, zscore=False):

    max_len = max([len(m) for m in motifs])
    motif_a = np.zeros((len(motifs), max_len * 4 + 3), dtype=float)
    for i, m in enumerate(motifs):
        motif_a[i, 0] = len(m)
        if thresholds:
            motif_a[i, 2] = thresholds[i]
        else:
            if zscore:
                motif_a[i, 2] = -100
            else:
                motif_a[i, 2] = m.pwm_min_score()
        motif_a[i, 1] = m.pwm_min_score()
        motif_a[i, 3 : 3 + len(m) * 4] = np.array(m.logodds).flatten()

    shm = smm.SharedMemory(motif_a.nbytes)
    b = np.ndarray(motif_a.shape, dtype=motif_a.dtype, buffer=shm.buf)
    b[:] = motif_a[:]  # Copy the original data into shared memory

    return shm.name, motif_a.shape, motif_a.dtype


def pwms_from_scandata(scandata):
    shm = SharedMemory(scandata.motif_array)
    motifs = np.ndarray(
        shape=scandata.motif_array_shape,
        dtype=scandata.motif_array_dtype,
        buffer=shm.buf,
    )
    pwms = [
        motif[3 : 3 + int(motif[0]) * 4].reshape(int(motif[0]), 4).tolist()
        for motif in motifs
    ]
    return pwms


def get_scandata(smm, motifs, seqs, gc_bin_list, thresholds=None, zscore=False):
    seq_list = smm.ShareableList([seq.upper() for seq in seqs])
    motif_name, motif_shape, motif_dtype = motifs_to_shareable_array(
        smm,
        motifs,
        thresholds,
        zscore,
    )
    gc_bin_list = smm.ShareableList(gc_bin_list)
    scandata = ScanData(
        motif_array=motif_name,
        motif_array_shape=motif_shape,
        motif_array_dtype=motif_dtype,
        gc_bin_list=gc_bin_list,
        seq_list=seq_list,
    )
    return scandata


def scan_seqs(
    seqs,
    motifs,
    thresholds,
    gc_bin_list,
    nreport=1,
    scan_rc=True,
    motifs_meanstd=None,
    zscore=False,
):
    for seq in seqs:
        row = []
        seq = seq.upper()
        for i, (motif, cutoff) in enumerate(zip(motifs, thresholds)):
            pwm_min_score = motif.pwm_min_score()
            cutoff = cutoff
            if zscore:
                seq_gc_bin = get_seq_bin(seq, gc_bin_list)
                m_mean, m_std = motifs_meanstd[seq_gc_bin][i]
                result = pwmscan(
                    seq, motif.logodds, cutoff * m_std + m_mean, nreport, scan_rc
                )
                result = [[(row[0] - m_mean) / m_std, row[1], row[2]] for row in result]
                # result = [row for row in result if row[0] >= cutoff]
            else:
                logger.error(f"{seq} {motif.logodds} {cutoff} {nreport} {scan_rc}")
                result = pwmscan(seq, motif.logodds, cutoff, nreport, scan_rc)

            if len(result) == 0 and (cutoff is None or cutoff <= pwm_min_score):
                result = [[pwm_min_score, 0, 1]] * nreport
            row.append(result)

        yield row


def scan_seqs_worker(
    scandata, seq_ids, nreport=1, scan_rc=True, motifs_meanstd=None, zscore=False
):
    ret = []

    # Create logodds matrices from shared memory object
    shm = SharedMemory(scandata.motif_array)
    motifs = np.ndarray(
        shape=scandata.motif_array_shape,
        dtype=scandata.motif_array_dtype,
        buffer=shm.buf,
    )
    pwms = [
        motif[3 : 3 + int(motif[0]) * 4].reshape(int(motif[0]), 4).tolist()
        for motif in motifs
    ]
    gc_bin_list = [
        (scandata.gc_bin_list[i], scandata.gc_bin_list[i + 1])
        for i in range(0, len(scandata.gc_bin_list), 2)
    ]

    for idx in seq_ids:
        seq = scandata.seq_list[idx]
        row = []
        for i, pwm in enumerate(pwms):
            pwm_min_score = motifs[i][1]
            cutoff = motifs[i][2]

            if zscore:
                seq_gc_bin = get_seq_bin(seq, gc_bin_list)

                m_mean, m_std = motifs_meanstd[seq_gc_bin][i]
                result = pwmscan(seq, pwm, cutoff * m_std + m_mean, nreport, scan_rc)
                result = [[(row[0] - m_mean) / m_std, row[1], row[2]] for row in result]
                # result = [row for row in result if row[0] >= cutoff]
            else:
                result = pwmscan(seq, pwm, cutoff, nreport, scan_rc)

            if len(result) == 0 and (cutoff is None or cutoff <= pwm_min_score):
                result = [[pwm_min_score, 0, 1]] * nreport
            row.append(result)

        ret.append((idx, row))

    return ret


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
    input_table, genome, scoring, pfmfile=None, ncpus=None, zscore=True, gc=True
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
        Specify a PFM file for scanning.

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

    seqs = as_fasta(check_regions, genome=genome).seqs
    size = int(np.median([len(seq) for seq in seqs]))
    s = Scanner(ncpus=ncpus)
    s.set_motifs(pfmfile)
    s.set_genome(genome)
    s.set_background(genome=genome, gc=gc, size=size)

    scores = []
    if scoring == "count":
        logger.info("setting threshold")
        s.set_threshold(fpr=FPR)
        logger.info("creating count table")
        for row in s.count(seqs):
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
        for row in s.best_score(seqs, zscore=zscore, gc=gc):
            scores.append(row)
        logger.info("done")

    motif_names = [m.id for m in read_motifs(pfmfile)]

    logger.info("creating dataframe")
    dtype = "float16"
    if scoring == "count":
        dtype = int
    df = pd.DataFrame(scores, index=idx, columns=motif_names, dtype=dtype)

    return df


def scan_table(s, inputfile, fa, motifs, cutoff, bgfile, nreport, scan_rc, pvalue):
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
    result_it = s.best_score(fa, scan_rc=scan_rc, zscore=zscore, gc=gcnorm)
    # header
    yield "\t{}".format("\t".join([m.id for m in motifs]))
    # score table
    for i, scores in enumerate(result_it):
        yield "{}\t{}".format(fa.ids[i], "\t".join(["{:4f}".format(x) for x in scores]))


def scan_normal(
    s,
    inputfile,
    fa,
    motifs,
    cutoff,
    bgfile,
    nreport,
    scan_rc,
    pvalue,
    bed,
    zscore,
    gcnorm,
):

    result_it = s.scan(fa, nreport=nreport, scan_rc=scan_rc, zscore=zscore, gc=gcnorm)
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
    pvalue=None,
    bgfile=None,
    genome=None,
    ncpus=None,
    zscore=False,
    gcnorm=False,
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
            s, inputfile, fa, motifs, cutoff, bgfile, nreport, scan_rc, pvalue
        )
    elif score_table:
        it = scan_score_table(s, fa, motifs, scan_rc, zscore=zscore, gcnorm=gcnorm)
    else:
        it = scan_normal(
            s,
            inputfile,
            fa,
            motifs,
            cutoff,
            bgfile,
            nreport,
            scan_rc,
            pvalue,
            bed,
            zscore=zscore,
            gcnorm=gcnorm,
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
    pvalue=False,
    bgfile=None,
    genome=None,
    ncpus=None,
    zscore=True,
    gcnorm=True,
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
        pvalue=pvalue,
        bgfile=bgfile,
        genome=genome,
        ncpus=ncpus,
        zscore=zscore,
        gcnorm=gcnorm,
    ):
        print(line, file=fo)

    if should_close:
        try:
            fo.close()
        except Exception:
            pass


def scan_to_best_match(
    fname, motifs, ncpus=None, genome=None, score=False, zscore=False, gc=False
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

    if motifs is None or isinstance(motifs, str):
        motifs = read_motifs(motifs)

    logger.debug("scanning %s...", fname)
    result = dict([(m.id, []) for m in motifs])

    seqs = as_fasta(fname, genome=genome).seqs

    if score:
        it = s.best_score(seqs, zscore=zscore, gc=gc)
    else:
        it = s.best_match(seqs, zscore=zscore, gc=gc)
    for scores in it:
        for motif, score in zip(motifs, scores):
            result[motif.id].append(score)

    return result


def parse_threshold_values(motif_file, cutoff):
    motifs = read_motifs(motif_file)
    d = parse_cutoff(motifs, cutoff)
    threshold = {}
    for m in motifs:
        c = m.pwm_min_score() + (m.pwm_max_score() - m.pwm_min_score()) * d[m.id]
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
                    seq, motif.logodds, motif.pwm_min_score(), nreport, scan_rc
                )
                result = [[(row[0] - m_mean) / m_std, row[1], row[2]] for row in result]
                result = [row for row in result if row[0] >= cutoff]
            else:
                result = pwmscan(seq, motif.logodds, cutoff, nreport, scan_rc)
            if cutoff <= motif.pwm_min_score() and len(result) == 0:
                result = [[motif.pwm_min_score(), 0, 1]] * nreport

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

        if ncpus is None:
            self.ncpus = int(MotifConfig().get_default_params()["ncpus"])
        else:
            self.ncpus = ncpus

        if self.ncpus > 1:
            # try:
            #    ctx = mp.get_context('spawn')
            #    self.pool = ctx.Pool(processes=self.ncpus)
            # except AttributeError:
            self.pool = mp.Pool(processes=self.ncpus)

        self.use_cache = False
        if self.config.get_default_params().get("use_cache", False):
            self._init_cache()

    def __del__(self):
        # Close the pool because of memory leak
        if hasattr(self, "pool"):
            self.pool.close()

    def _init_cache(self):
        try:
            self.cache = make_region().configure(
                "dogpile.cache.pylibmc",
                expiration_time=3600,
                arguments={"url": ["127.0.0.1"], "binary": True}
                #    'dogpile.cache.dbm',
                #    expiration_time = 3600,
                #    arguments = {
                #        'filename': 'cache.dbm'
                #    }
            )
            self.use_cache = True
        except Exception as e:
            sys.stderr.write("failed to initialize cache\n")
            sys.stderr.write("{}\n".format(e))

    def set_motifs(self, motifs=None):
        if motifs is None:
            motifs = self.config.get_default_params().get("motif_db", None)

        try:
            # Check if motifs is a list of Motif instances
            motifs[0].to_pwm()
            tmp = NamedTemporaryFile(mode="w", delete=False)
            for m in motifs:
                tmp.write("{}\n".format(m.to_pwm()))
            tmp.close()
            motif_file = tmp.name
        except AttributeError:
            motif_file = motifs

        self.motifs = motif_file
        self.motif_ids = [m.id for m in read_motifs(motif_file)]
        self.checksum = {}
        if self.use_cache:
            chksum = xxhash.xxh64("\n".join(sorted(self.motif_ids))).digest()
            self.checksum[self.motifs] = chksum

    def _meanstd_from_seqs(self, motifs, seqs):
        table = []
        thresholds = [m.pwm_min_score() for m in motifs]
        for x in self.scan(
            seqs, motifs=motifs, thresholds=thresholds, nreport=1, scan_rc=True
        ):
            table.append([row[0][0] for row in x])

        for motif, scores in zip(motifs, np.array(table).transpose()):
            yield motif, np.mean(scores), np.std(scores)  # cutoff

    def _threshold_from_seqs(self, motifs, seqs, fpr):
        table = []
        seq_gc_bins = [self.get_seq_bin(seq) for seq in seqs]
        thresholds = [m.pwm_min_score() for m in motifs]
        for gc_bin, result in zip(
            seq_gc_bins,
            self.scan(
                seqs,
                motifs=motifs,
                thresholds=thresholds,
                nreport=1,
                scan_rc=True,
                progress=False,
            ),
        ):
            try:
                table.append([gc_bin] + [row[0][0] for row in result])
            except Exception:
                print(gc_bin)
                for row in result:
                    print(row)
                    row[0][0]

        df = pd.DataFrame(table, columns=["gc_bin"] + [m.id for m in motifs])
        return df

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
        lock.acquire()

        try:
            with Cache(CACHE_DIR) as cache:
                scan_motifs = []
                for bin in bins:
                    if bin not in self.meanstd:
                        self.meanstd[bin] = {}
                    bin_seqs = [s for s, b in zip(seqs, seq_bins) if b == bin]

                    for motif in motifs:
                        k = "e{}|{}|{}".format(motif.hash(), self.background_hash, bin)

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
                                motif.hash(), self.background_hash, bin
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
                                motif.hash(), self.background_hash, gcbin
                            )
                            mean = self.meanstd[gcbin][motif.id][0]
                            cache.set(k, [mean, std])
                            self.meanstd[gcbin][motif.id] = mean, std
        except sqlite3.DatabaseError:
            print_cluster_error_message()
            sys.exit(1)
        lock.release()

        for gc_bin in self.gc_bins:
            gc_bin = "{:.2f}-{:.2f}".format(*gc_bin)
            if gc_bin not in self.meanstd:
                valid_bins = []
                for b in self.gc_bins:
                    bstr = "{:.2f}-{:.2f}".format(b[0], b[1])
                    if bstr in self.meanstd:
                        valid_bins.append(((b[0] + b[1]) / 2, bstr))

                v = float(gc_bin.split("-")[1])
                _, bstr = sorted(valid_bins, key=lambda x: abs(x[0] - v))[0]
                logger.warn(f"Using {bstr}")
                self.meanstd[gc_bin] = self.meanstd[bstr]

        for gc_bin in self.meanstd:
            for i, motif_id in enumerate(self.motif_ids):
                self.meanstd[gc_bin][i] = self.meanstd[gc_bin][motif_id]

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
            logger.warn(
                "Genome and FASTA filename specified for background. Using custom FASTA file."
            )

        if fname:
            if not os.path.exists(fname):
                raise IOError("Background file {} does not exist!".format(fname))

            self.background = Fasta(fname)
            self.background_hash = file_checksum(fname)
            return

        if not genome:
            if self.genome:
                genome = self.genome
            else:
                raise ValueError("Need either genome or filename for background.")

        logger.debug("using background: genome {} with size {}".format(genome, size))
        lock.acquire()
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
                            logger.debug("using {} sequences".format(nseq))
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
        lock.release()

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

        if not self.motifs:
            raise ValueError("please run set_motifs() first")

        motifs = read_motifs(self.motifs)
        gc_bins = ["{:.2f}-{:.2f}".format(*gc_bin) for gc_bin in self.gc_bins]

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

        lock.acquire()
        try:
            with Cache(CACHE_DIR) as cache:
                scan_motifs = []
                self._threshold = None
                for motif in motifs:
                    k = "{}|{}|{:.4f}|{}".format(
                        motif.hash(),
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
                    logger.info("determining FPR-based threshold")
                    df = self._threshold_from_seqs(scan_motifs, seqs, fpr).set_index(
                        "gc_bin"
                    )
                    if self._threshold is None:
                        self._threshold = df
                    else:
                        self._threshold = pd.concat((self._threshold, df), axis=1)
                    for motif in scan_motifs:
                        k = "{}|{}|{:.4f}|{}".format(
                            motif.hash(),
                            self.background_hash,
                            fpr,
                            ",".join(sorted(gc_bins)),
                        )
                        cache.set(k, df[motif.id])
        except sqlite3.DatabaseError:
            print_cluster_error_message()
            sys.exit(1)
        lock.release()
        self.threshold_str = "{}_{}_{}_{}".format(
            fpr, threshold, self.background_hash, ",".join(sorted(gc_bins))
        )

    def set_genome(self, genome):
        """
        set the genome to be used for:
            - converting regions to sequences
        """
        if not genome:
            return

        # raises error if checks fail
        Genome(genome)

        self.genome = genome

    def count(self, seqs, nreport=100, scan_rc=True):
        """
        count the number of matches above the cutoff
        returns an iterator of lists containing integer counts
        """
        for matches in self.scan(seqs, nreport=nreport, scan_rc=scan_rc):
            counts = [len(m) for m in matches]
            yield counts

    def total_count(self, seqs, nreport=100, scan_rc=True):
        """
        count the number of matches above the cutoff
        returns an iterator of lists containing integer counts
        """

        count_table = [counts for counts in self.count(seqs, nreport, scan_rc)]
        return np.sum(np.array(count_table), 0)

    def best_score(self, seqs, scan_rc=True, zscore=False, gc=False):
        """
        give the score of the best match of each motif in each sequence
        returns an iterator of lists containing floats
        """
        self.set_threshold(threshold=0.0, gc=gc)
        for matches in self.scan(
            seqs, nreport=1, scan_rc=scan_rc, zscore=zscore, gc=gc
        ):
            scores = np.array(
                [sorted(m, key=lambda x: x[0])[0][0] for m in matches if len(m) > 0]
            )
            yield scores

    def best_match(self, seqs, scan_rc=True, zscore=False, gc=False):
        """
        give the best match of each motif in each sequence
        returns an iterator of nested lists containing tuples:
        (score, position, strand)
        """
        self.set_threshold(threshold=0.0)
        for matches in self.scan(
            seqs, nreport=1, scan_rc=scan_rc, zscore=zscore, gc=gc
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
                return "{:.2f}-{:.2f}".format(b_start, b_end)

        logger.error(
            "Error determining seq: {}, bins: {}".format(seq, str(self.gc_bins))
        )
        raise ValueError()

    def get_motif_mean_std(self, gc_bin, motif):
        if gc_bin in self.meanstd:
            if motif not in self.meanstd[gc_bin]:
                raise ValueError("Motif mean and std not initialized")
        else:
            logger.warn(
                "GC% {} not present in genome, setting to closest GC% bin".format(
                    gc_bin
                )
            )
            # Ideally this should not happen. If we get here, this means that
            # a sequence has a GC% that does not occur in the genome. This
            # can happen during tests with small genome, or when a sequence
            # is scanned that does not originate from the genome.
            # We use the closest GC% bin from the genome instead.
            valid_bins = []
            for b in self.gc_bins:
                bstr = "{:.2f}-{:.2f}".format(b[0], b[1])
                if bstr in self.meanstd:
                    valid_bins.append(((b[0] + b[1]) / 2, bstr))

            v = float(gc_bin.split("-")[1])
            _, bstr = sorted(valid_bins, key=lambda x: abs(x[0] - v))[0]
            logger.warn(f"Using {bstr}")
            self.meanstd[gc_bin] = self.meanstd[bstr]
        return self.meanstd[gc_bin][motif]

    def scan(
        self,
        seqs,
        thresholds=None,
        motifs=None,
        nreport=100,
        scan_rc=True,
        zscore=False,
        gc=False,
        progress=True,
    ):
        """
        Scan a set of regions or sequences.
        """
        if motifs is None:
            motifs = read_motifs(self.motifs)

        if isinstance(seqs, str):
            seqs = Fasta(seqs).seqs

        if isinstance(seqs, Fasta):
            seqs = seqs.seqs

        if zscore:
            if gc:
                if len(self.meanstd) <= 1:
                    self.set_meanstd(gc=gc)
            else:
                if len(self.meanstd) != 1:
                    self.set_meanstd(gc=gc)

        logger.debug("Scanning")
        if not thresholds:
            thresholds = self.get_gc_thresholds(seqs, motifs=motifs, zscore=zscore)
            thresholds = [thresholds.get(m.id, None) for m in motifs]

        flat_list = [float(item) for sublist in self.gc_bins for item in sublist]
        if self.ncpus == 1:
            for row in scan_seqs(
                seqs,
                motifs,
                thresholds,
                self.gc_bins,
                nreport=nreport,
                scan_rc=scan_rc,
                motifs_meanstd=self.meanstd,
                zscore=zscore,
            ):
                yield row
        else:
            with SharedMemoryManager() as smm:
                scandata = get_scandata(
                    smm, motifs, seqs, flat_list, thresholds, zscore
                )
                seq_ids = list(range(len(seqs)))
                batch = 200
                chunk = batch * self.ncpus
                if chunk > len(seqs):
                    chunk = len(seqs)

                if progress:
                    pbar = tqdm(total=len(seqs))
                with ProcessPoolExecutor(self.ncpus) as exe:
                    # We submit in chunks to keep memory use in check.
                    # If everything is submitted at once, memory explodes as the memory claimed
                    # by the futures is not released.
                    seq_idx = 0
                    
                    for j in range(0, len(seqs), chunk):
                        hold = {}
                        fs = [
                            exe.submit(
                                scan_seqs_worker,
                                scandata,
                                seq_ids[j + i : j + i + batch],
                                nreport=nreport,
                                scan_rc=scan_rc,
                                motifs_meanstd=self.meanstd,
                                zscore=zscore,
                            )
                            for i in range(0, chunk, batch)
                        ]
                        for future in as_completed(fs):
                            # Need to return the scanning results in order, but they may come back
                            # unordered.
                            rows = future.result()
                            hold[rows[0][0]] = rows
                            # for row in future.result():
                            #     if row[0] == seq_idx:
                            #         yield row[1]
                            #         seq_idx += 1
                            #     else:
                            #         hold[row[0]] = row[1]

                            # while seq_idx in hold:
                            #     yield hold[seq_idx]
                            #     seq_idx += 1

                            if progress:
                                pbar.update(batch)

                            del future
                        for idx in sorted(hold.keys()):
                            for i, row in hold[idx]:
                                yield row

                if progress:
                    pbar.close()

    def get_gc_thresholds(self, seqs, motifs=None, zscore=False):
        # Simple case, only one threshold
        if np.all(self.threshold.nunique(axis=0) == 1):
            return self.threshold.iloc[0].to_dict()

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

        nseqs = int(20000 / np.sum(list(gc_bin_count.values())))
        t = {}
        maxt = pd.Series([m.pwm_max_score() for m in motifs], index=_threshold.columns)
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
            val = _threshold.loc[idx, motif].quantile(0.99, interpolation="higher")
            if val < maxt.loc[motif]:
                t[motif] = val
            else:
                t[motif] = None
        return t
