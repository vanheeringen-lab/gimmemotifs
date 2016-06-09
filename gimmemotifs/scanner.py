import os
import re
import sys
from functools import partial
from tempfile import mkdtemp

# "hidden" features, in development
try:
    import MOODS.tools
    import MOODS.parsers
    import MOODS.scan
except ImportError:
    pass

import numpy as np

from gimmemotifs.config import MotifConfig
from gimmemotifs.fasta import Fasta
from gimmemotifs.genome_index import GenomeIndex
from gimmemotifs.c_metrics import pwmscan
from gimmemotifs.motif import read_motifs
from gimmemotifs.utils import parse_cutoff,get_seqs_type
from gimmemotifs.genome_index import rc

# only used when using cache, should not be a requirement
try:
    from dogpile.cache import make_region
    from dogpile.cache.api import NO_VALUE
    from cityhash import CityHash64
except ImportError:
    pass 

def load_motifs(motif_file, cutoff=0.95):
    motifs = read_motifs(open(motif_file))
    d = parse_cutoff(motifs, cutoff)
    cutoffs = []
    for m in motifs:
        c = m.pwm_min_score() + (m.pwm_max_score() - m.pwm_min_score()) * d[m.id]
        cutoffs.append(c)
    
    return zip(motifs, cutoffs)

def scan_sequence(seq, motifs, nreport, scan_rc):
    
    ret = []
    # scan for motifs
    for motif, cutoff in motifs:
        result = pwmscan(seq, motif.pwm, cutoff, nreport, scan_rc)
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
    
    def __init__(self):
        self.config = MotifConfig()

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
        self.motifs = motifs
        self.motif_ids = [m.id for m in read_motifs(open(motifs))]
        self.checksum = {}
        if self.use_cache:
            chksum = CityHash64("\n".join(sorted(self.motif_ids)))
            self.checksum[self.motifs] = chksum

    def set_threshold(self):
        """
        cutoff  # based on motif matrix
        p-value # moods only
        FDR     # based on bg file
        """
        
        pass

    def set_genome(self, genome):
        """
        set the genome to be used for:
            - converting regions to sequences
            - background for MOODS
        """
        index_dir = os.path.join(self.config.get_index_dir(), genome)
        if not os.path.exists(index_dir) or not os.path.isdir(index_dir):
            raise ValueError("index for {} does not exist".format(genome))
        self.index_dir = index_dir
    
    def count(self, seqs, nreport=100, scan_rc=True, cutoff=0.95):
        """
        count the number of matches above the cutoff
        returns an iterator of lists containing integer counts
        """
        for matches in self.scan(seqs, nreport, scan_rc, cutoff):
            counts = [len(m) for m in matches]
            yield counts
     
    def total_count(self, seqs, nreport=100, scan_rc=True, cutoff=0.95):
        """
        count the number of matches above the cutoff
        returns an iterator of lists containing integer counts
        """
        
        count_table = [counts for counts in self.count(seqs, nreport, scan_rc, cutoff)]
        return np.sum(np.array(count_table), 0)

    def best_score(self, seqs, scan_rc=True):
        """
        give the score of the best match of each motif in each sequence
        returns an iterator of lists containing floats
        """
        for matches in self.scan(seqs, 1, scan_rc, cutoff=0):
            scores = [sorted(m, lambda x,y: 
                                    cmp(y[0], x[0])
                                    )[0][0] for m in matches]
            yield scores
 
    def best_match(self, seqs, scan_rc=True):
        """
        give the best match of each motif in each sequence
        returns an iterator of nested lists containing tuples:
        (score, position, strand)
        """
        for matches in self.scan(seqs, 1, scan_rc, cutoff=0):
            top = [sorted(m, lambda x,y: 
                                    cmp(y[0], x[0])
                                    )[0] for m in matches]
            yield top
    
   
    def scan(self, seqs, nreport=100, scan_rc=True, cutoff=0.95):
        """
        scan a set of regions / sequences
        """

        # determine input type
        seqs_type = get_seqs_type(seqs)
        
        # Fasta object
        if seqs_type.startswith("fasta"):
            if seqs_type.endswith("file"):
                seqs = Fasta(seqs)
            
            it = self._scan_sequences(seqs.seqs, 
                    nreport, scan_rc, cutoff)
        # regions or BED
        else:
            if seqs_type == "regionfile":
                seqs = [l.strip() for l in open(seqs)]
            it = self._scan_regions(seqs, 
                    nreport, scan_rc, cutoff)
        
        for result in it:
            yield result


    def _scan_regions(self, regions, nreport, scan_rc, cutoff=0.95):
        index_dir = self.index_dir
        motif_file = self.motifs
        motif_digest = self.checksum.get(motif_file, None)

        # determine which regions are not in the cache 
        scan_regions = regions
        if self.use_cache:
            scan_regions = []
            for region in regions:
                key = str((region, index_dir, motif_digest, nreport, scan_rc, cutoff))
                ret = self.cache.get(key)
                if ret == NO_VALUE:
                    scan_regions.append(region)
        
        # scan the regions that are not in the cache
        if len(scan_regions) > 0:
            n = 12
            
            genome_index = GenomeIndex(index_dir)
           
            motifs = load_motifs(motif_file, cutoff)
     
            scan_func = partial(scan_region_mult,
                genome_index=genome_index,
                motifs=motifs,
                nreport=nreport,
                scan_rc=scan_rc)
    
            jobs = []
            chunksize = len(scan_regions) / n + 1
            for i in range(n):
                job = pool.apply_async(scan_func, (scan_regions[i * chunksize:( i+ 1) * chunksize],))
                jobs.append(job)
            
            # return values or store values in cache
            i = 0
            for job in jobs:
                for ret in job.get():
                    if self.use_cache:
                        # store values in cache    
                        region = scan_regions[i]
                        key = str((region, index_dir, motif_digest, nreport, scan_rc, cutoff))
                        self.cache.set(key, ret)
                    else:
                        #return values
                        yield ret
                    i += 1
    
        if self.use_cache: 
            # return results from cache
            for region in regions:
                key = str((region, index_dir, motif_digest, nreport, scan_rc, cutoff))
                ret = self.cache.get(key)
                if ret == NO_VALUE or ret is None:
                    raise Exception("cache is not big enough to hold all " 
                                    "results, try increasing the cache size "
                                    "or disable cache")
                yield ret
    
    def _scan_sequences(self, seqs, nreport, scan_rc, cutoff=0.95):
        
        motif_file = self.motifs
        motif_digest = self.checksum.get(motif_file, None)
        
        scan_seqs = seqs
        if self.use_cache:
            # determine which sequences are not in the cache 
            hashes = dict([(s.upper(), CityHash64(s.upper())) for s in seqs])
            scan_seqs = []
        
            for seq,seq_hash in hashes.items():
                key = str((seq_hash, motif_digest, nreport, scan_rc, cutoff))
                ret = self.cache.get(key)
                if ret == NO_VALUE or ret is None:
                    scan_seqs.append(seq.upper())
        
        # scan the sequences that are not in the cache
        if len(scan_seqs) > 0:
            n = 12
            motifs = load_motifs(motif_file, cutoff)
            scan_func = partial(scan_seq_mult,
                motifs=motifs,
                nreport=nreport,
                scan_rc=scan_rc)
    
            jobs = []
            chunksize = len(scan_seqs) / n + 1 
            for i in range(n):
                job = pool.apply_async(scan_func, (scan_seqs[i * chunksize:(i  + 1) * chunksize],))
                jobs.append(job)
            
            # store values in cache
            i = 0
            for job in jobs:
                for ret in job.get():
                    if self.use_cache:
                        h = hashes[scan_seqs[i]]
                        key = str((h, motif_digest, nreport, scan_rc, cutoff))
                        self.cache.set(key, ret)
                    else: 
                        yield ret
                    i += 1
        
        if self.use_cache:
            # return results from cache
            for seq in seqs:
                key = str((hashes[seq.upper()], motif_digest, nreport, scan_rc, cutoff))
                ret = self.cache.get(key)
                if ret == NO_VALUE or ret is None:
                    raise Exception("cache is not big enough to hold all " 
                                    "results, try increasing the cache size "
                                    "or disable cache")
                    
                yield ret

try: 
    from gimmemotifs.mp import pool
except ImportError:
    pass
