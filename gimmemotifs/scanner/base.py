"""
Scanner class
"""
import logging
import math
import os
import sqlite3
import sys
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
    file_hash,
    parse_threshold_values,
    print_cluster_error_message,
    scan_seq_mult,
)
from gimmemotifs.utils import as_fasta

logger = logging.getLogger("gimme.scanner")
FPR = 0.01
THRESHOLD = 0.95
LOCK = Lock()


class Scanner(object):
    """
    scan sequences with motifs
    """

    genome = None
    background = None
    background_hash = None
    meanstd = {}
    gc_bins = [(0, 1)]
    fpr = None
    motifs = None
    motif_ids = []
    _threshold = None
    _seed = None

    def __init__(self, ncpus=None, random_state=None, progress=None):
        self.config = MotifConfig()
        self.random_state = random_state
        self.progress = progress

        if ncpus is None:
            ncpus = self.config.get_default_params()["ncpus"]
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

    @property
    def random_state(self):
        return self._random_state

    @random_state.setter
    def random_state(self, random_state):
        self._random_state = random_state
        if random_state is not None:
            # can be used to set random.seed()
            self._seed = random_state.get_state()[1][0].item()

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

    def set_meanstd(self, gc=False):
        if not self.background:
            self.set_background(gc=gc)

        self.meanstd = {}
        seqs = self.background.seqs
        if gc:
            seq_bins = [s.split(" ")[-1] for s in self.background.ids]
            bins = sorted(set(seq_bins))
        else:
            seq_bins = ["0.00-1.00"] * len(seqs)
            bins = ["0.00-1.00"]

        motifs = read_motifs(self.motifs)
        LOCK.acquire()
        try:
            with Cache(CACHE_DIR) as cache:

                # for each bin, load the meanstd of each motif or mark it as missing
                scan_gc_bins = {}
                for bin in bins:
                    if bin not in self.meanstd:
                        self.meanstd[bin] = {}

                    for motif in motifs:
                        k = "e{}|{}|{}".format(motif.hash, self.background_hash, bin)
                        results = cache.get(k)
                        if results is None:
                            # missing values
                            scan_gc_bins.setdefault(bin, []).append(motif)
                        else:
                            self.meanstd[bin][motif.id] = results

                # generate missing values for each motif per gc_bin
                total_scans = sum(len(v) for v in scan_gc_bins.values())
                if total_scans > 0:
                    pbar = tqdm(
                        desc="Determining mean and stddev for motifs",
                        unit=" motifs",
                        total=total_scans,
                        disable=self._disable,  # can be silenced
                    )

                    for bin, scan_motifs in scan_gc_bins.items():
                        bin_seqs = [s for s, b in zip(seqs, seq_bins) if b == bin]
                        if len(bin_seqs) == 0:
                            # no background sequences with this GC%
                            pbar.update(len(scan_motifs))
                            continue
                        for motif, mean, std in self._meanstd_from_seqs(
                            scan_motifs, bin_seqs
                        ):
                            k = "e{}|{}|{}".format(
                                motif.hash, self.background_hash, bin
                            )
                            cache.set(k, [mean, std])
                            self.meanstd[bin][motif.id] = mean, std
                            pbar.update(1)

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

                if len(valid_bins) == 0:
                    logger.error(
                        f"No sequences found fo gc bin '{gc_bin}'. "
                        "This likely occurred because you set gc=True with a background file"
                    )
                    sys.exit()

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
                gc_bins = [(0.0, 0.2), (0.8, 1.0)]
                for b in np.arange(0.2, 0.799, 0.05):
                    gc_bins.append((round(b, 2), round(b + 0.05, 2)))
                gc_bins = sorted(gc_bins)
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
            self.background_hash = file_hash(fname)  # was file_checksum()
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
                                tmp.name, genome, nseq, size, gc_bins, self.random_state
                            )
                            fa = as_fasta(tmp.name, genome=genome)
                    else:
                        fa = RandomGenomicFasta(genome, size, nseq, self._seed)
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
                fpr = FPR
                logger.info(f"Using default FPR of {fpr}")
            else:
                threshold = THRESHOLD
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
                    df = self._threshold_from_seqs(scan_motifs, seqs)
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

    def set_genome(self, genome=None, genomes_dir=None):
        """
        Set the genome to converting regions to sequences

        Parameters
        ----------
        genome : str, optional
            Path to the genome fasta or
            the genome name as found in the genomepy genomes directory.

        genomes_dir : str, optional
            Path to the genomepy genomes directory.
            Taken from the genomepy config if unspecified.
        """
        if genome is not None:
            # raises error if checks fail
            genome = Genome(genome, genomes_dir, rebuild=False).filename

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

    def best_score(self, seqs, scan_rc=True, zscore=False, gc=False):
        """
        give the score of the best match of each motif in each sequence
        returns an iterator of lists containing floats
        """
        self.set_threshold(threshold=0.0, gc=gc)
        for matches in self.scan(seqs, 1, scan_rc, zscore=zscore, gc=gc):
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
        for matches in self.scan(seqs, 1, scan_rc, zscore=zscore, gc=gc):
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
            if round(b_start, 2) < gc <= round(b_end, 2):
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

    def scan(self, seqs, nreport=100, scan_rc=True, zscore=False, gc=False):
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
            desc="Scanning",
            unit=" sequences",
            total=len(seqs),
            disable=self._disable,  # can be silenced
        )

        batch_size = 50000
        for batch_idx in range(0, len(seqs), batch_size):
            it = self._scan_sequences(
                seqs.seqs[batch_idx : batch_idx + batch_size],
                nreport,
                scan_rc,
                zscore,
            )
            for result in it:
                yield result
                pbar.update(1)
        pbar.close()

    def get_gc_thresholds(self, seqs, motifs=None, zscore=False):
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
