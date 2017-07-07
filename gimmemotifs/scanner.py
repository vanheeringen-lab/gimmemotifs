import os
import re
import sys
import gc
from functools import partial
from tempfile import mkdtemp,NamedTemporaryFile
import logging
from multiprocessing import Pool
import six

# "hidden" features, in development
try:
    import MOODS.tools
    import MOODS.parsers
    import MOODS.scan
except ImportError:
    pass

from diskcache import Cache
import numpy as np
from scipy.stats import scoreatpercentile

from gimmemotifs.background import RandomGenomicFasta
from gimmemotifs.config import MotifConfig,CACHE_DIR
from gimmemotifs.fasta import Fasta
from gimmemotifs.genome_index import GenomeIndex
from gimmemotifs.c_metrics import pwmscan
from gimmemotifs.motif import read_motifs
from gimmemotifs.utils import parse_cutoff,as_fasta,file_checksum
from gimmemotifs.genome_index import rc,check_genome

try:
    import copy_reg
    import types
    def _pickle_method(m):
        if m.im_self is None:
            return getattr, (m.im_class, m.im_func.func_name)
        else:
            return getattr, (m.im_self, m.im_func.func_name)

    copy_reg.pickle(types.MethodType, _pickle_method)
except:
    pass

# only used when using cache, should not be a requirement
try:
    from dogpile.cache import make_region
    from dogpile.cache.api import NO_VALUE
    import xxhash
except ImportError:
    pass 

logger = logging.getLogger("gimme.scanner")
config = MotifConfig()

def scan_to_best_match(fname, motifs, ncpus=None, genome=None, score=False):
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
    
    if isinstance(motifs, six.string_types) and os.path.exists(motifs):
        with open(motifs) as f: 
            motifs = read_motifs(f)

    logger.debug("scanning %s...", fname)
    result = dict([(m.id, []) for m in motifs])
    if score:
        it = s.best_score(fname)
    else:
        it = s.best_match(fname)
    for scores in it:
        for motif,score in zip(motifs, scores):
            result[motif.id].append(score)

    return result

def parse_threshold_values(motif_file, cutoff):
    with open(motif_file) as f:
        motifs = read_motifs(f)
    d = parse_cutoff(motifs, cutoff)
    threshold = {}
    for m in motifs:
        c = m.pwm_min_score() + ( 
                    m.pwm_max_score() - m.pwm_min_score()
                ) * d[m.id]
        threshold[m.id] = c
    return threshold

def scan_sequence(seq, motifs, nreport, scan_rc):
    
    ret = []
    # scan for motifs
    for motif, cutoff in motifs:
        if cutoff is None:
            ret.append([])
        else:
            result = pwmscan(seq, motif.pwm, cutoff, nreport, scan_rc)
            if cutoff <= motif.pwm_min_score() and len(result) == 0:
                result = [[motif.pwm_min_score(), 0, 1]] * nreport
            ret.append(result)

    # return results
    return ret

def scan_region(region, genome_index, motifs, nreport, scan_rc):
    
    # retrieve sequence
    chrom,start,end = re.split(r'[:-]', region)
    seq = genome_index.get_sequence(chrom, int(start), int(end)).upper()
    
    return scan_sequence(seq, motifs, nreport, scan_rc)

def scan_seq_mult(seqs, motifs, nreport, scan_rc):
    ret = []
    for seq in seqs:
        result = scan_sequence(seq.upper(), motifs, nreport, scan_rc)
        ret.append(result)
    return ret

def scan_region_mult(regions, genome_index, motifs, nreport, scan_rc):
    ret = []
    for region in regions:
        result = scan_region(region, genome_index, motifs, nreport, scan_rc)
        ret.append(result)
    return ret


def scan_fa_with_motif_moods(fo, motifs, matrices, bg, thresholds, nreport, scan_rc=True):

    scanner = MOODS.scan.Scanner(7)
    scanner.set_motifs(matrices, bg, thresholds)

    ret = []
    for name, seq in fo.items():
        l = len(seq)

        scan_seq = seq.upper()
        if scan_rc:
            scan_seq = "".join((scan_seq, "N"*50, rc(scan_seq)))
        results = scanner.scan_max_hits(scan_seq, nreport)
        for motif,result in zip(motifs, results):
            matches = []
            for match in result:
                strand = 1
                pos = match.pos
                if scan_rc:
                    if pos > l:
                        pos = l - (pos - l - 50) - len(motif)
                        strand = -1
                matches.append((pos, match.score, strand))
            ret.append((motif, {name: matches}))

    return ret

def scan_fa_with_motif_moods_count(fo, motifs, matrices, bg, thresholds, nreport, scan_rc=True):
    scanner = MOODS.scan.Scanner(7)
    scanner.set_motifs(matrices, bg, thresholds)

    ret = []
    for name, seq in fo.items():
        l = len(seq)

        scan_seq = seq.upper()
        if scan_rc:
            scan_seq = "".join((scan_seq, "N"*50, rc(scan_seq)))
        results = scanner.counts_max_hits(scan_seq, nreport)
        ret.append((name, results))

    return ret

def calc_threshold_moods(m, c):
    m_min = MOODS.tools.min_score(m)
    m_max = MOODS.tools.max_score(m)

    return m_min + (m_max - m_min) * c

def scan_it_moods(infile, motifs, cutoff, bgfile, nreport=1, scan_rc=True, pvalue=None, count=False):
    tmpdir = mkdtemp()
    matrices = []
    pseudocount = 1e-3
    #sys.stderr.write("bgfile: {}\n".format(bgfile))
    bg = MOODS.tools.bg_from_sequence_dna("".join(Fasta(bgfile).seqs), 1)

    for motif in motifs:
        pfmname = os.path.join(tmpdir, "{}.pfm".format(motif.id))
        with open(pfmname, "w") as f:
            matrix = np.array(motif.pwm).transpose()
            for line in [" ".join([str(x) for x in row]) for row in matrix]:
                f.write("{}\n".format(line))

        matrices.append(MOODS.parsers.pfm_log_odds(pfmname, bg, pseudocount))

    thresholds = []
    if pvalue is not None:
        thresholds = [MOODS.tools.threshold_from_p(m, bg, float(pvalue)) for m in matrices]
        #sys.stderr.write("{}\n".format(thresholds))
    else:
        thresholds = [calc_threshold_moods(m, float(cutoff)) for m in matrices]

    scanner = MOODS.scan.Scanner(7)
    scanner.set_motifs(matrices, bg, thresholds)

    config = MotifConfig()
    ncpus =  int(config.get_default_params()['ncpus'])
    fa = Fasta(infile)
    chunk = 500
    if (len(fa) / chunk) < ncpus:
        chunk = len(fa) / (ncpus + 1)

    jobs = []
    func = scan_fa_with_motif_moods
    if count:
        func = scan_fa_with_motif_moods_count


    pool = Pool()
    for i in range(0, len(fa), chunk):
        jobs.append(pool.apply_async(
                                          func,
                                          (fa[i:i + chunk],
                                          motifs,
                                          matrices,
                                          bg,
                                          thresholds,
                                          nreport,
                                          scan_rc,
                                          )))

    for job in jobs:
        for ret in job.get():
            yield ret

class Scanner(object):
    """
    scan sequences with motifs
    """
    
    def __init__(self, ncpus=None):
        self.config = MotifConfig()
        self.threshold = None
        self.index_dir = None

        if ncpus is None:
            self.ncpus = int(MotifConfig().get_default_params()["ncpus"])
        else:
            self.ncpus = ncpus

        self.use_cache = False
        if self.config.get_default_params().get("use_cache", False):
            self._init_cache()
            
    def _init_cache(self):
        try:
            self.cache = make_region().configure(
                'dogpile.cache.pylibmc',
                expiration_time = 3600,
                arguments = {
                    'url':["127.0.0.1"],
                    'binary': True,
                    }
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

    def set_motifs(self, motifs):
        try:
            # Check if motifs is a list of Motif instances
            motifs[0].to_pwm()
            tmp = NamedTemporaryFile(mode="w", delete=False)
            for m in motifs:
                tmp.write("{}\n".format(m.to_pwm()))
            tmp.close()
            motif_file = tmp.name
        except AttributeError as e:
            motif_file = motifs

        self.motifs = motif_file
        with open(motif_file) as f:
            self.motif_ids = [m.id for m in read_motifs(f)]
        self.checksum = {}
        if self.use_cache:
            chksum = xxhash.xxh64("\n".join(sorted(self.motif_ids))).digest()
            self.checksum[self.motif_file] = chksum

    def _threshold_from_seqs(self, motifs, seqs, fpr):
        scan_motifs = [(m, m.pwm_min_score()) for m in motifs]
        
        table = []
        for x in self._scan_sequences_with_motif(scan_motifs, seqs, 1, True):
            table.append([row[0][0] for row in x])
                
        for (motif, _), scores in zip(scan_motifs, np.array(table).transpose()):
            min_score = motif.pwm_min_score()
            cutoff = 0
            if len(scores) > 0:
                opt_score = scoreatpercentile(scores, 100 - (100 * fpr))
                cutoff = (opt_score - min_score) / (motif.pwm_max_score() - min_score)
            yield motif, opt_score#cutoff


    def set_threshold(self, fpr=None, threshold=None, genome=None, 
                        length=200, filename=None):
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
        if threshold:
            if fpr:
                raise ValueError("Need either fpr or threshold.")
            if genome:
                sys.stderr.write(
                    "Parameter genome ignored when threshold is specified.\n"
                    "Did you want to use fpr?\n")
            if filename:
                sys.stderr.write(
                    "Parameter filename ignored when threshold is specified.\n"
                    "Did you want to use fpr?\n")

        if genome and filename:
            raise ValueError("Need either genome or filename.")
    
        if fpr:
            fpr = float(fpr)
            if not (0.0 < fpr < 1.0):
                raise ValueError("Parameter fpr should be between 0 and 1")
       
        if not self.motifs:
            raise ValueError("please run set_motifs() first")

        thresholds = {}
        with open(self.motifs) as f: 
            motifs = read_motifs(f)
        
        if threshold is not None:
            self.threshold = parse_threshold_values(self.motifs, threshold) 
            return
        
        if filename:
            if not os.path.exists(filename):
                raise IOError(
                        "File {} does not exist.".format(filename)
                        )
            
            bg_hash = file_checksum(filename)
            seqs = Fasta(filename).seqs
        elif genome:
            bg_hash = "{}\{}".format(genome, int(length))
        else:
            raise ValueError("Need genome or filename")

        with Cache(CACHE_DIR) as cache:
            scan_motifs = []
            for motif in motifs:
                k = "{}|{}|{:.4f}".format(motif.hash(), bg_hash, fpr)
           
                threshold = cache.get(k)
                if threshold is None:
                    scan_motifs.append(motif)
                else:
                    if np.isclose(threshold, motif.pwm_max_score()):
                        thresholds[motif.id] = None
                    else:
                        thresholds[motif.id] = threshold
                
            if len(scan_motifs) > 0:
                if genome:
                    check_genome(genome)    
                    sys.stderr.write("Determining threshold for fpr {} and length {} based on {}\n".format(fpr, int(length), genome))
                    index = os.path.join(config.get_index_dir(), genome)
                    fa = RandomGenomicFasta(index, length, 10000)
                    seqs = fa.seqs
                else: 
                    sys.stderr.write("Determining threshold for fpr {} based on {}\n".format(fpr, filename))
                for motif, threshold in self._threshold_from_seqs(scan_motifs, seqs, fpr):
                    k = "{}|{}|{:.4f}".format(motif.hash(), bg_hash, fpr)
                    cache.set(k, threshold)
                    if np.isclose(threshold, motif.pwm_max_score()):
                        thresholds[motif.id] = None
                    else:
                        thresholds[motif.id] = threshold
        
        self.threshold_str = "{}_{}_{}_{}_{}".format(fpr, threshold, genome,
                                        length, filename)
        self.threshold = thresholds

    def set_genome(self, genome):
        """
        set the genome to be used for:
            - converting regions to sequences
            - background for MOODS
        """
        if not genome:
            return
        
        if os.path.exists(genome) and os.path.isdir(genome):
            index_dir = genome
        else:
            index_dir = os.path.join(self.config.get_index_dir(), genome)
            if not os.path.exists(index_dir) or not os.path.isdir(index_dir):
                raise ValueError("index for {} does not exist".format(genome))
        self.index_dir = index_dir
    
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
    
    def best_score(self, seqs, scan_rc=True):
        """
        give the score of the best match of each motif in each sequence
        returns an iterator of lists containing floats
        """
        self.set_threshold(threshold=0.0)
        for matches in self.scan(seqs, 1, scan_rc):
            scores = [sorted(m, key=lambda x: x[0])[0][0] for m in matches if len(m) > 0]
            yield scores
 
    def best_match(self, seqs, scan_rc=True):
        """
        give the best match of each motif in each sequence
        returns an iterator of nested lists containing tuples:
        (score, position, strand)
        """
        self.set_threshold(threshold=0.0)
        for matches in self.scan(seqs, 1, scan_rc):
            top = [sorted(m, key=lambda x: x[0])[0] for m in matches]
            yield top
   
    def scan(self, seqs, nreport=100, scan_rc=True):
        """
        scan a set of regions / sequences
        """

        if not self.threshold:
            sys.stderr.write(
                "Using default threshold of 0.95. "
                "This is likely not optimal!\n"
                )
            self.set_threshold(threshold=0.95)

        seqs = as_fasta(seqs, index_dir=self.index_dir)
           
        it = self._scan_sequences(seqs.seqs, 
                    nreport, scan_rc)
       
        for result in it:
            yield result


    def _scan_regions(self, regions, nreport, scan_rc):
        index_dir = self.index_dir
        motif_file = self.motifs
        motif_digest = self.checksum.get(motif_file, None)

        # determine which regions are not in the cache 
        scan_regions = regions
        if self.use_cache:
            scan_regions = []
            for region in regions:
                key = str((region, index_dir, motif_digest, nreport, scan_rc))
                ret = self.cache.get(key)
                if ret == NO_VALUE:
                    scan_regions.append(region)
        
        # scan the regions that are not in the cache
        if len(scan_regions) > 0:
            
            genome_index = GenomeIndex(index_dir)
           
            with open(self.motifs) as f:
                motifs = [(m, self.threshold[m.id]) for m in read_motifs(f)]
            scan_func = partial(scan_region_mult,
                genome_index=genome_index,
                motifs=motifs,
                nreport=nreport,
                scan_rc=scan_rc)
    
            for region, ret in self._scan_jobs(scan_func, scan_regions):
                # return values or store values in cache
                if self.use_cache:
                    # store values in cache    
                    key = str((region, index_dir, motif_digest, nreport, scan_rc, self.threshold_str))
                    self.cache.set(key, ret)
                else:
                    #return values
                    yield ret
    
        if self.use_cache: 
            # return results from cache
            for region in regions:
                key = str((region, index_dir, motif_digest, nreport, scan_rc, self.threshold_str))
                ret = self.cache.get(key)
                if ret == NO_VALUE or ret is None:
                    raise Exception("cache is not big enough to hold all " 
                                    "results, try increasing the cache size "
                                    "or disable cache")
                yield ret
    
    def _scan_sequences_with_motif(self, motifs, seqs, nreport, scan_rc):
        scan_func = partial(scan_seq_mult,
            motifs=motifs,
            nreport=nreport,
            scan_rc=scan_rc)

        for ret in self._scan_jobs(scan_func, seqs):
            yield ret[1]

    def _scan_sequences(self, seqs, nreport, scan_rc):
        
        motif_file = self.motifs
        motif_digest = self.checksum.get(motif_file, None)
        
        scan_seqs = seqs
        if self.use_cache:
            # determine which sequences are not in the cache 
            hashes = dict([(s.upper(), xxhash.xxh64(s.upper()).digest()) for s in seqs])
            scan_seqs = []
        
            for seq,seq_hash in hashes.items():
                key = str((seq_hash, motif_digest, nreport, scan_rc, self.threshold_str))
                ret = self.cache.get(key)
                if ret == NO_VALUE or ret is None:
                    scan_seqs.append(seq.upper())
        
        # scan the sequences that are not in the cache
        if len(scan_seqs) > 0:
            with open(self.motifs) as f:
                motifs = [(m, self.threshold[m.id]) for m in read_motifs(f)]
            scan_func = partial(scan_seq_mult,
                motifs=motifs,
                nreport=nreport,
                scan_rc=scan_rc)
    
            for seq, ret in self._scan_jobs(scan_func, scan_seqs):
                if self.use_cache:
                    h = hashes[seq]
                    key = str((h, motif_digest, nreport, scan_rc, self.threshold_str))
                    self.cache.set(key, ret)
                else: 
                    yield ret
      
        if self.use_cache:
            # return results from cache
            for seq in seqs:
                key = str((hashes[seq.upper()], motif_digest, nreport, scan_rc, self.threshold_str))
                ret = self.cache.get(key)
                if ret == NO_VALUE or ret is None:
                    raise Exception("cache is not big enough to hold all " 
                                    "results, try increasing the cache size "
                                    "or disable cache")
                    
                yield ret
            
    def _scan_jobs(self, scan_func, scan_seqs):
        batchsize = 5000
        if self.ncpus > 1:
            for i in range((len(scan_seqs) - 1) // batchsize + 1):
                batch = scan_seqs[i * batchsize:( i+ 1) * batchsize]
                chunksize = len(batch) // self.ncpus + 1
                pool = Pool(processes=self.ncpus)
                jobs = []
                for j in range((len(batch) - 1) // chunksize + 1):
                    job = pool.apply_async(
                            scan_func, 
                            (batch[j * chunksize:(j + 1) * chunksize],)
                            )
                    jobs.append(job)
                
                for k,job in enumerate(jobs):
                    for ret in job.get():
                        region = batch[k]
                        yield region, ret
                pool.close()
                gc.collect()
        else:
            for i in range((len(scan_seqs) - 1) // batchsize + 1):
                for j,ret in enumerate(scan_func(scan_seqs[i * batchsize:( i+ 1) * batchsize])):
                    yield scan_seqs[i], ret


