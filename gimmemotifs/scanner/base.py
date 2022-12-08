"""
Scanner class
"""
import logging
import math
import os
import sqlite3
import time
from collections import Counter
from functools import partial
from multiprocessing import Lock, Pool
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
from diskcache import Cache
from genomepy import Genome
from sklearn.preprocessing import scale
from tqdm.auto import tqdm

from gimmemotifs.background import RandomGenomicFasta, gc_bin_bedfile
from gimmemotifs.config import CACHE_DIR, MotifConfig
from gimmemotifs.fasta import Fasta
from gimmemotifs.motif import read_motifs
from gimmemotifs.scanner.utils import (
    cluster_error,
    file_hash,
    parse_threshold_values,
    scan_seq_mult,
)
from gimmemotifs.utils import as_fasta

logger = logging.getLogger("gimme.scanner")
FPR = 0.01
THRESHOLD = 0.95
LOCK = Lock()


class Scanner(object):
    """
    scan sequences for motifs
    """

    genome = None
    background = None
    background_hash = None
    threshold = None
    meanstd = {}
    gc_bins = ["0.00-1.00"]
    fpr = None
    motifs = None
    motif_ids = []
    _seed = None

    def __init__(self, ncpus=None, random_state=None, progress=None):
        self.config = MotifConfig()
        self.random_state = random_state
        self.progress = progress

        if ncpus is None:
            ncpus = MotifConfig().get_default_params()["ncpus"]
        self.ncpus = int(ncpus)

        if self.ncpus > 1:
            self.pool = Pool(processes=self.ncpus)

    def __del__(self):
        # Close the Pool to release memory
        try:
            self.pool.close()
            self.pool.join()
        except (ValueError, AttributeError):
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

    def set_motifs(self, motifs):
        motif_file = motifs
        if hasattr(motifs[0], "to_ppm"):
            tmp = NamedTemporaryFile(mode="w", delete=False)
            for m in motifs:
                tmp.write("{}\n".format(m.to_ppm()))
            tmp.close()
            motif_file = tmp.name

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
        seq_gc_bins = [self.get_seq_gc_bin(seq) for seq in seqs]

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
        """mean and standard deviation of motif detection in background sequences
        per gc bin.

        Only required if zscore=True.
        """
        if not self.motifs:
            raise ValueError("please run set_motifs() first")

        if not self.background:
            raise ValueError("please run set_background() first")

        self.meanstd = {}
        motifs = read_motifs(self.motifs)

        LOCK.acquire()
        try:
            with Cache(CACHE_DIR) as cache:

                # for each bin, load the meanstd of each motif
                # or mark it as missing
                scan_gc_bins = {}
                for gc_bin in self.gc_bins:
                    if gc_bin not in self.meanstd:
                        self.meanstd[gc_bin] = {}

                    for motif in motifs:
                        k = "e{}|{}|{}".format(motif.hash, self.background_hash, gc_bin)
                        results = cache.get(k)

                        if results is None:
                            # missing values
                            scan_gc_bins.setdefault(gc_bin, []).append(motif)
                        else:
                            self.meanstd[gc_bin][motif.id] = results

                # generate missing values for each motif per gc_bin
                total_scans = sum(len(v) for v in scan_gc_bins.values())
                if total_scans > 0:
                    pbar = tqdm(
                        desc="Determining mean stdev",
                        unit=" motifs",
                        total=total_scans,
                        disable=self._disable,  # can be silenced
                    )

                    for gc_bin, gc_bin_motifs in scan_gc_bins.items():
                        gc_bin_seqs = [
                            seq
                            for seq in self.background.seqs
                            if self.get_seq_gc_bin(seq) == gc_bin
                        ]
                        if len(gc_bin_seqs) == 0:
                            # no background sequences with this GC%
                            pbar.update(len(gc_bin_motifs))
                            continue
                        for motif, mean, std in self._meanstd_from_seqs(
                            gc_bin_motifs, gc_bin_seqs
                        ):
                            k = "e{}|{}|{}".format(
                                motif.hash, self.background_hash, gc_bin
                            )
                            cache.set(k, [mean, std])
                            self.meanstd[gc_bin][motif.id] = mean, std
                            pbar.update(1)

                # fill gc_bins missing in the background with the
                # nearest gc_bin present in the background
                for gc_bin in self.gc_bins:
                    if len(self.meanstd[gc_bin]) == 0:
                        valid_bins = []
                        for bstr in self.gc_bins:
                            if len(self.meanstd.get(bstr, [])) > 0:
                                valid_bins.append(
                                    (sum(float(v) for v in bstr.split("-")) / 2, bstr)
                                )

                        v = float(gc_bin.split("-")[1])
                        _, bstr = sorted(valid_bins, key=lambda x: abs(x[0] - v))[0]
                        logger.debug(f"using mean stds of GC% {bstr} for GC% {gc_bin}")
                        self.meanstd[gc_bin] = self.meanstd[bstr]
                        for motif in motifs:
                            k = "e{}|{}|{}".format(
                                motif.hash, self.background_hash, gc_bin
                            )
                            cache.set(k, self.meanstd[gc_bin][motif.id])

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
            cluster_error()
        LOCK.release()

    def set_background(
        self,
        fasta=None,
        size=200,
        nseq=None,
        gc=False,
        gc_bins=None,
    ):
        """Set the background to use for FPR and z-score calculations.

        If fasta is given, this file is used as-is.
        Else, sequences are drawn from the genome.

        Parameters
        ----------
        fasta : str, optional
            Path to FASTA file to use as background sequences.
            If unspecified, background sequences are sampled from the genome.

        size : int, optional
            Size of sequences to retrieve. The default is 200.

        nseq : int, optional
            Total number of background sequences.

        gc : bool, optional
            Equally distribute GC percentages in background sequences.

        gc_bins : list, optional
            A list bins containing desired GC percentages
            e.g. [(0.0, 0.50), (0.50, 1.0)]

        Sets
        ----
        self.gc_bins
            a list of unique gc bins in string format
            (["0.00-1.00"] without gc correction)
        self.background
        self.background_hash
        """
        if fasta:
            if not os.path.exists(fasta):
                raise FileNotFoundError(f"Background file {fasta} does not exist!")
            self.background = Fasta(fasta)
            self.background_hash = file_hash(fasta)
            self.gc_bins = ["0.00-1.00"]
            if gc:
                logger.warning(
                    "No GC correction possible with a background fasta. Skipping."
                )
            return

        if self.genome is None:
            raise ValueError(
                "please run set_genome() first, or provide a FASTA file as background."
            )

        size = int(size)
        if gc_bins is None:
            if gc:
                gc_bins = [
                    (0.0, 0.2),  # wider
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
                    (0.8, 1.00),  # wider
                ]
            else:
                gc_bins = [(0, 1)]
        self.gc_bins = ["{:.2f}-{:.2f}".format(*b) for b in gc_bins]
        if nseq is None:
            nseq = max(10_000, len(gc_bins) * 1000)

        LOCK.acquire()
        try:
            with Cache(CACHE_DIR) as cache:
                self.background_hash = "d{}:{}:{}:{}".format(
                    file_hash(self.genome), size, gc, str(gc_bins)
                )
                fa = cache.get(self.background_hash)
                if fa is None:
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
            cluster_error()
        LOCK.release()

        self.background = fa

    def set_thresholds(self, fpr=None, threshold=None):
        """Set motif scanning threshold.

        If FPR is given, the threshold is based on background sequences.
        If a threshold file/value is given, no background sequences are used.

        Parameters
        ----------
        fpr : float, optional
            Desired false positive rate, between 0.0 and 1.0.

        threshold : float or str, optional
            Desired motif threshold, expressed as the fraction of the
            difference between minimum and maximum score of the PFM.
            Should either be a float between 0.0 and 1.0 or a filename
            with thresholds as created by 'gimme threshold'.

        Sets
        ----
        self.threshold
        self.fpr
        """
        if not self.motifs:
            raise ValueError("please run set_motifs() first")

        if threshold and fpr:
            raise ValueError("Need either fpr or threshold.")

        if threshold is None and fpr is None:
            if self.genome:
                fpr = FPR
                logger.info(f"Using default FPR of {fpr}")
            else:
                threshold = THRESHOLD
                logger.info(f"Using default threshold of {threshold}.")
                logger.info(
                    "This is likely not ideal, use set_genome() for improvements."
                )

        motifs = read_motifs(self.motifs)

        if threshold is not None:
            if len(self.gc_bins) > 1:
                self.gc_bins = ["0.00-1.00"]
                logger.warning(
                    "No GC correction possible with a threshold value. Skipping."
                )

            d = parse_threshold_values(motifs, threshold)
            self.threshold = pd.DataFrame(d, index=[0])
            self.threshold = self.threshold.join(
                pd.DataFrame(
                    self.gc_bins, index=[0] * len(self.gc_bins), columns=["gc_bin"]
                )
            )
            self.threshold.set_index("gc_bin", inplace=True)
            self.fpr = None
            return

        if not self.background:
            raise ValueError("please run set_background() first")

        if fpr is not None:
            fpr = float(fpr)
            if not (0.0 < fpr < 1.0):
                raise ValueError("Parameter fpr should be between 0 and 1")
            self.fpr = fpr

        LOCK.acquire()
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
            cluster_error()
        LOCK.release()

    def set_genome(self, genome=None, genomes_dir=None):
        """
        Set the genome to be used for:
            - converting regions to sequences
            - background for MOODS

        Parameters
        ----------
        genome : str or None
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

    def count(self, seqs, nreport=100, scan_rc=True):
        """
        count the number of matches above the threshold
        returns an iterator of lists containing integer counts
        """
        for matches in self.scan(seqs, nreport, scan_rc):
            counts = [len(m) for m in matches]
            yield counts

    def total_count(self, seqs, nreport=100, scan_rc=True):
        """
        count the number of matches above the threshold
        returns an array containing integer counts
        """
        count_table = [counts for counts in self.count(seqs, nreport, scan_rc)]
        return np.sum(np.array(count_table), 0)

    def best_score(self, seqs, scan_rc=True, zscore=False, dtype=None):
        """
        give the score of the best match of each motif in each sequence
        returns an iterator of lists containing floats
        """
        self.set_thresholds(threshold=0.0)
        # use numpy's default dtype if None is given
        dtype = float if dtype is None else dtype
        for matches in self.scan(seqs, 1, scan_rc, zscore):
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
        self.set_thresholds(threshold=0.0)
        for matches in self.scan(seqs, 1, scan_rc, zscore):
            yield [m[0] for m in matches]

    def get_seq_gc_bin(self, seq):
        if len(str(seq)) == 0:
            gc = 0.01
        else:
            useq = seq.upper()
            gc = round((useq.count("G") + useq.count("C")) / len(useq), 2)
        if gc == 0:
            gc = 0.01

        for gc_bin in self.gc_bins:
            b_start, b_end = [float(b) for b in gc_bin.split("-")]
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
        # Precalculate mean and stddev for z-score calculation
        if zscore and not all(gc_bin in self.meanstd for gc_bin in self.gc_bins):
            self.set_meanstd()
        seqs = as_fasta(seqs, genome=self.genome)

        # progress bar
        pbar = tqdm(
            desc="Scanning",
            unit=" sequences",
            total=len(seqs),
            disable=self._disable,  # can be silenced
        )

        batch_size = 50_000
        logger.debug("Scanning")
        for batch_idx in range(0, len(seqs), batch_size):
            it = self._scan_sequences(
                seqs.seqs[batch_idx : batch_idx + batch_size], nreport, scan_rc, zscore
            )
            for result in it:
                yield result
                pbar.update(1)
        pbar.close()

    def get_thresholds(self, seqs, motifs=None, zscore=False):
        if self.threshold is None:
            raise ValueError("please run set_thresholds() first")

        # Simple case: one unique threshold per motif
        # this is because set_thresholds(threshold) was used.
        if np.all(self.threshold.nunique(axis=0) == 1):
            return self.threshold.iloc[0].to_dict()

        if self.fpr is None:
            raise ValueError("please run set_thresholds() with fpr first")

        max_seqs = 20_000
        if len(seqs) > max_seqs:
            random = np.random if self.random_state is None else self.random_state
            seqs = random.choice(seqs, size=max_seqs, replace=True)

        if motifs is None:
            motifs = read_motifs(self.motifs)
        seq_gc_bins = [self.get_seq_gc_bin(seq) for seq in seqs]

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
            scan_seq_mult,
            motifs=motifs,
            nreport=nreport,
            scan_rc=scan_rc,
        )
        for ret in self._scan_jobs(scan_func, seqs):
            yield ret[1]

    def _scan_sequences(self, seqs, nreport, scan_rc, zscore):
        thresholds = self.get_thresholds(seqs, zscore=zscore)
        motifs = [(m, thresholds[m.id]) for m in read_motifs(self.motifs)]
        motifs_meanstd = None
        if zscore:
            if not self.meanstd:
                raise ValueError("please run set_meanstd() first")
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
            func = self._scan_jobs_multi_threaded
        else:
            func = self._scan_jobs_single_threaded

        for ret in func(scan_func, scan_seqs):
            yield ret

    def _scan_jobs_single_threaded(self, scan_func, scan_seqs):
        # non-parallel job scanning
        batchsize = 1000
        for i in range((len(scan_seqs) - 1) // batchsize + 1):
            batch_seqs = scan_seqs[i * batchsize : (i + 1) * batchsize]
            seq_gc_bins = [self.get_seq_gc_bin(seq) for seq in batch_seqs]

            for _j, ret in enumerate(scan_func(batch_seqs, seq_gc_bins)):
                yield scan_seqs[i], ret

    def _scan_jobs_multi_threaded(self, scan_func, scan_seqs):
        median_len = np.median([len(x) for x in scan_seqs])
        chunksize = 200_000 // int(median_len)  # 1000 seqs for len 200
        # prepare for parallel processing
        k = 0
        max_queue_size = 2 * self.ncpus
        jobs = []

        # loop over each job/chunk, and keep adding them to the queue
        for i in range(math.ceil(len(scan_seqs) / chunksize)):
            batch_seqs = scan_seqs[i * chunksize : (i + 1) * chunksize]
            seq_gc_bins = [self.get_seq_gc_bin(seq) for seq in batch_seqs]
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
