import os
import re
from multiprocessing import Pool
from functools import partial

from gimmemotifs.config import MotifConfig
from gimmemotifs.fasta import Fasta
from gimmemotifs.genome_index import GenomeIndex
from gimmemotifs.c_metrics import pwmscan
from gimmemotifs.motif import read_motifs
from gimmemotifs.utils import parse_cutoff

from dogpile.cache import make_region
from dogpile.cache.api import NO_VALUE

from cityhash import CityHash64

cache = make_region().configure(
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

def scan_regions(regions, index_dir, motif_file, nreport, scan_rc, cutoff=0.95, use_cache=False):
    
    # determine which regions are not in the cache 
    scan_regions = regions
    if use_cache:
        scan_regions = []
        for region in regions:
            key = str((region, index_dir, motif_file, nreport, scan_rc, cutoff))
            ret = cache.get(key)
            if ret == NO_VALUE:
                scan_regions.append(region)
    
    # scan the regions that are not in the cache
    if len(scan_regions) > 0:
        n = 12
        p = Pool(n)
        
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
            job = p.apply_async(scan_func, (scan_regions[i * chunksize:( i+ 1) * chunksize],))
            jobs.append(job)
        
        # store values in cache
        i = 0
        for job in jobs:
            for ret in job.get():
                if use_cache:
                    # store values in cache    
                    region = scan_regions[i]
                    key = str((region, index_dir, motif_file, nreport, scan_rc, cutoff))
                    cache.set(key, ret)
                else:
                    yield ret
                i += 1

    if use_cache: 
        # return results from cache
        for region in regions:
            key = str((region, index_dir, motif_file, nreport, scan_rc, cutoff))
            ret = cache.get(key)
            yield ret

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

def scan_sequences(seqs, motif_file, nreport, scan_rc, cutoff=0.95, use_cache=False):
    
    scan_seqs = seqs
    if use_cache:
        # determine which sequences are not in the cache 
        hashes = dict([(s.upper(), CityHash64(s.upper())) for s in seqs])
        scan_seqs = []
    
        for seq,seq_hash in hashes.items():
            key = str((seq_hash, motif_file, nreport, scan_rc, cutoff))
            ret = cache.get(key)
            if ret == NO_VALUE:
                scan_seqs.append(seq.upper())
    
    # scan the sequences that are not in the cache
    if len(scan_seqs) > 0:
        n = 12
        p = Pool(n)
        
        motifs = load_motifs(motif_file, cutoff)
        scan_func = partial(scan_seq_mult,
            motifs=motifs,
            nreport=nreport,
            scan_rc=scan_rc)

        jobs = []
        chunksize = len(scan_seqs) / n + 1 
        for i in range(n):
            job = p.apply_async(scan_func, (scan_seqs[i * chunksize:(i  + 1) * chunksize],))
            jobs.append(job)
        
        # store values in cache
        i = 0
        for job in jobs:
            for ret in job.get():
                if use_cache:
                    h = hashes[scan_seqs[i]]
                    key = str((h, motif_file, nreport, scan_rc, cutoff))
                    cache.set(key, ret)
                else: 
                    yield ret
                i += 1
    
    if use_cache:
        # return results from cache
        for seq in seqs:
            key = str((hashes[seq.upper()], motif_file, nreport, scan_rc, cutoff))
            ret = cache.get(key)
            yield ret

class Scanner(object):
    """
    scan sequences with motifs
    """
    
    def __init__(self):
        self.config = MotifConfig()

    def set_motifs(self, motifs):
        self.motifs = motifs

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
    
    def count(self, seqs, nreport=100, scan_rc=True):
        """
        count the number of matches above the cutoff
        returns an iterator of lists containing integer counts
        """
        for matches in self.scan(seqs, nreport, scan_rc):
            counts = [len(m) for m in matches]
            yield counts
    
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
    
    def _get_seqs_type(self, seqs):
        """
        automagically determine input type
        the following types are detected:
            - Fasta object
            - FASTA file
            - list of regions
            - region file
            - BED file
        """
        
        region_p = re.compile(r'^(.+):(\d+)-(\d+)$')
        if isinstance(seqs, Fasta):
            return "fasta"
        elif isinstance(seqs, list):
            if len(seqs) == 0:
                raise ValueError("empty list of sequences to scan")
            else:
                if region_p.search(seqs[0]):
                    return "regions"
                else:
                    raise ValueError("unknown region type")
        elif isinstance(seqs, str):
            if os.path.isfile(seqs):
                try:
                    f = Fasta(seqs)
                    return "fastafile"
                except:
                    pass
                try:
                    line = open(seqs).readline().strip()
                    if region_p.search(line):
                        return "regionfile"
                    else:
                        if line.split("\t"):
                            int(vals[1]), int(vals[2])
                            return "bedfile"
                except:
                    raise ValueError("unknown type for argument seqs")
            else:
                raise ValueError("no file found with name {}".format(seqs))
        else:
            raise ValueError("unknown type for argument seqs")
    
    def scan(self, seqs, nreport=100, scan_rc=True, cutoff=0.95):
        """
        scan a set of regions / sequences
        """

        # determine input type
        seqs_type = self._get_seqs_type(seqs)
        
        # Fasta object
        if seqs_type.startswith("fasta"):
            if seqs_type.endswith("file"):
                seqs = Fasta(seqs)
            
            it = scan_sequences(seqs.seqs, 
                    self.motifs, nreport, scan_rc, cutoff, use_cache=True)
        # regions or BED
        else:
            if seqs_type == "regionfile":
                seqs = [l.strip() for l in open(seqs)]
            it = scan_regions(seqs, self.index_dir,
                    self.motifs, nreport, scan_rc, cutoff, use_cache=True)
        
        for result in it:
            yield result
