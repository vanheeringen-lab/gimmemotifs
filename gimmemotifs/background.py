# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

""" 
Classes to generate background files in FASTA format. Two different methods are
included: MarkovFasta, which generates a background according to a 1st order 
Markov model and MatchedGenomicFasta, which generates a background with a 
similar genomic distribution as the input.

"""
from __future__ import division

# Python imports
import gzip
import os
import random
import re
import sys
from tempfile import NamedTemporaryFile
from random import choice
import logging

# External imports
import numpy as np
import pandas as pd
import pybedtools
from genomepy import Genome

# GimmeMotifs imports
from gimmemotifs import mytmpdir
from gimmemotifs.fasta import Fasta
from gimmemotifs.config import CACHE_DIR

logger = logging.getLogger("gimme.background")

def create_random_genomic_bedfile(out, genome, length, n):
    features = Genome(genome).get_random_sequences(n, length)

    # Write result to bedfile
    tmp = open(out, "w")
    for chrom,start,end in features:
        tmp.write("%s\t%d\t%d\n" % (chrom, start, end))
    tmp.flush()    

def create_promoter_bedfile(out, genefile, length, n):
    strand_map = {"+":True, "-":False, 1:True, -1:False, "1":True, "-1":False}

    features = []
    if genefile.endswith(".gz"):
        fin = gzip.open(genefile, "rt")
    else:
        fin = open(genefile)

    for line in fin:
        if not line.startswith("track") or line.startswith("#"):
            (chrom, start, end, _name, score, strand) = line[:-1].split("\t")[:6]
            start, end = int(start), int(end)
            strand= strand_map[strand]
            if strand:
                if start - length >= 0:
                    features.append([chrom, start - length, start, strand])
            else:
                features.append([chrom, end, end + length, strand])
    fin.close()
    
    if n < len(features):
        features = random.sample(features, n)
    else:
        sys.stdout.write("Too few promoters to generate %s random promoters! Just using all of them." % n)

    # Write result to temporary bedfile
    tmp = open(out, "w")
    for chrom, start, end, strand in sorted(features, key=lambda x: x[0]):
        tmp.write("%s\t%s\t%s\t0\t0\t%s\n" % (chrom, start, end, {True:"+",False:"-"}[strand]))
    tmp.flush()

class MarkovFasta(Fasta):
    """ 
    Generates a new Fasta object containing sequences using a 1st order Markov
    model, based on the input sequences. By default 10 times as many sequences 
    will be generated with the same length as the input sequences.

    Required arg 'fasta' is a Fasta object
    Optional arg 'length' can be used to generate sequences of a different length
    Optional arg 'n' specifies the number of sequences to generate
    Optional arg 'k' specifies the order of the Markov model, default is 1 for 1st
    order

    Returns a Fasta object
    
    Example:
        
        f = Fasta("input.fa")
        random_fa = MarkovFasta(f, multiply = 5)
        for id,seq in random_fa.items():
            print seq
    
    """
    
    def __init__(self, fasta, length=None, n=None, k=1, matrix_only=False):
        self.k = k

        # Initialize super Fasta object
        Fasta.__init__(self)
        
        # Initialize Markov transition matrix
        self._initialize_matrices(fasta.seqs, k=k)

        if matrix_only:
            return
        
        c = 0
        if not n:
            n = len(fasta)

        while len(self) < n:
            seq = choice(fasta.seqs)
            name = "random_Markov%s_%s" % (k,c)
            if length:
                random_seq = self._generate_sequence(length)
            else:
                random_seq = self._generate_sequence(len(seq))
            self.add(name, random_seq)    
            c += 1

    def _initialize_matrices(self, seqs, k=1, alphabet=None):
        if alphabet is None:
            alphabet = ['A','C','G','T']

        self.frequencies = {}
        kmercount = {}
        
        init = alphabet[:]
        for i in range(k - 1):
            new_init = []
            for x in init:
                for l in alphabet:
                    new_init.append(x + l)
            init = new_init[:]

        self.trans = dict([(word, dict([(l, 0.0) for l in alphabet])) for word in init])
        new_init = []
        for x in init:
            for l in alphabet:
                new_init.append(x + l)
        
        kmercount = dict([(word, 0) for word in new_init])
        lettercount =  dict([(word[:k], 0) for word in new_init])
        p = re.compile("^[%s]+$" % "".join(alphabet))
        total = 0
        for seq in seqs:
            seq = seq.upper()
            for i in range(len(seq) - k):
                if p.search(seq[i:i + k + 1]):
                    lettercount[seq[i:i + k]] += 1
                    kmercount[seq[i:i + k + 1]] += 1
                    total += 1
        
        for k,v in kmercount.items():
            self.trans[k[:-1]][k[-1]] = float(v)
        
        for k,v in self.trans.items():
            s = np.sum(np.array(list(v.values())))
            for x,y in v.items():
                v[x] = y / s
        
        self.init = {}
        total = float(np.sum(np.array(list(lettercount.values()))))
        for k,v in lettercount.items():
            self.init[k] = v / total
            
    def _generate_sequence(self, l):
        sequence = list(self._weighted_random(list(self.init.items())))
        for _ in range(l - self.k):
            sequence.append(
                    self._weighted_random(
                        list(self.trans["".join(sequence[-self.k:])].items())
                        ))
        return "".join(sequence)

    def _weighted_random(self, l):
        n = random.uniform(0,1)
        item = None
        for item, weight in l:
            if n < weight:
                break    
            else:
                n -= weight
        return item

def create_gc_bin_index(genome, fname):
    """Create index of GC content for a genome.

    Parameters
    ----------
    genome : str
        Genome name.
    fname : str
        Name of the index file.
    """
    logger.info("Creating index for genomic GC frequencies.")
    g = Genome(genome)
    fasta = g.filename
    sizes = g.props["sizes"]["sizes"]

    with NamedTemporaryFile() as tmp:
        pybedtools.BedTool().window_maker(g=sizes, w=100).nucleotide_content(fi=fasta).saveas(tmp.name)
        df = pd.read_csv(tmp.name, sep="\t", usecols=[0,1,2,4,9])

    df["w200"] = df.iloc[:,3].rolling(2, min_periods=2).mean()
    df["n200"] = df.iloc[:,4].rolling(2, min_periods=2).sum()
    df["w500"] = df.iloc[:,3].rolling(5, min_periods=5).mean()
    df["n500"] = df.iloc[:,4].rolling(5, min_periods=5).sum()
    cols = ["chrom", "start", "end", "w100", "n100", "w200", "n200", "w500", "n500"]
    df.columns = cols
    
    df.reset_index()[cols].to_feather(fname)

def gc_bin_bedfile(bedfile, genome, number, l=200, bins=None, random_state=None):
    """Create a BED file from different GC bins.
    
    Parameters
    ----------
    bedfile : str
        Name of the output BED file.
    genome : str
        Genome name.
    number : int
        Number of sequences to retrieve.
    l : int, optional
        Length of the sequences, default is 200.
    bins : list, optional
        GC frequency bins to use, for instance [(0,50),(50,100)]
    """
    if bins is None:
        bins = [(0.0, 0.2), (0.8, 1)]
        for b in np.arange(0.2, 0.799, 0.05):
            bins.append((b, b + 0.05))

    if number < len(bins):
        raise ValueError("Number of sequences requested < number of bins")

    fname = os.path.join(CACHE_DIR, "{}.gcfreq.feather".format(genome))
    if not os.path.exists(fname):
        create_gc_bin_index(genome, fname)
    
    df = pd.read_feather(fname)

    if l >= 100:
        col = "w{}".format(((l + 50)// 100) * 100)
    else:
        logger.warn("Can't create GC regions smaller than 100, using 100 instead")
        col = "w100"

    if not col in df.columns:
        df[col] = df.iloc[:,3].rolling(l//100, min_periods=l//100).mean()
        df[col.replace("w","n")] = df.iloc[:,3].rolling(l//100, min_periods=l//100).sum()

    df = df[df[col.replace("w","n")] < 0.1 * l]
    n = number // len(bins)
    
    with open(bedfile, "w") as f:
        pass

    with open(bedfile, 'a') as f:
        for b_start, b_end in bins:
            df_bin = df[(df[col] > b_start) & (df[col] <= b_end)]
            if df_bin.shape[0] > 0:
                df_bin = df_bin.sample(n, replace=True, random_state=random_state)
                df_bin["bin"] = "{:.2f}-{:.2f}".format(b_start, b_end)
                df_bin["end"] = df_bin["start"] + l
                df_bin[["chrom", "start", "end", "bin"]].to_csv(f, sep="\t", header=False, index=False)

def matched_gc_bedfile(bedfile, matchfile, genome, number):
    N_FRACTION = 0.1
    
    g = Genome(genome)
    genome_fa = g.filename
    try:
        fa = Fasta(matchfile)
        gc = [(seq.upper().count("C") + seq.upper().count("G")) / len(seq) for seq in fa.seqs]
        lengths = [len(seq) for seq in fa.seqs]
    except Exception:
        try:
            # pylint: disable=unexpected-keyword-arg
            fields = pd.read_csv(matchfile, comment="#", nrows=10, sep="\t").shape[1]
            bed = pybedtools.BedTool(matchfile)
            gc = [float(x[fields + 1]) for x in bed.nucleotide_content(fi=genome_fa)]
            lengths = [x.length for x in bed]
        except:
            sys.stderr.write("Please provide input file in BED or FASTA format\n")
            raise
    gc_hist,bins = np.histogram(gc, range=(0,1), bins=20)
    
    length = np.median(lengths)
    if np.std(lengths) > length * 0.05:
        sys.stderr.write("Sequences do not seem to be of equal length.\n")
        sys.stderr.write("GC% matched sequences of the median length ({}) will be created\n".format(length))

    if number:
        norm = number * gc_hist / (float(sum(gc_hist))) + 0.5
        inorm = norm.astype(np.int)

        s = np.argsort(norm - inorm)
        while sum(inorm) > number:
            if inorm[np.argmin(s)] > 0:
                inorm[np.argmin(s)] -= 1
            s[np.argmin(s)] = len(s)
        while sum(inorm) < number:
            inorm[np.argmax(s)] += 1
            s[np.argmax(s)] = 0
        gc_hist = inorm

    rnd = pybedtools.BedTool()
    out = open(bedfile, "w")
    #sys.stderr.write("Generating sequences\n")
    #sys.stderr.write("{}\n".format(number))
    
    # Create a file with chromosome sizes if it doesn't exist yet
    genome_size = genome_fa + ".sizes"
    del_size = False
    if not os.path.exists(genome_size):
        genome_size = NamedTemporaryFile().name
        del_size = True
        with open(genome_size, "w") as f:
            for seqname in g.keys():
                f.write("{}\t{}\n".format(seqname, len(g[seqname])))
    
    #sys.stderr.write("Done\n")
    features = []
    gc = []
    for bin_start, bin_end, count in zip(bins[:-1], bins[1:], gc_hist):
        #sys.stderr.write("CG {}-{}\n".format(bin_start, bin_end))
        if count > 0:
            rcount = 0
            c = 0            
            while count != rcount:
                c += 1
                if len(features) == 0 or c < 2:
                    # pylint: disable=unexpected-keyword-arg
                    n = number * (c ** 2 * 10)
                    if n < 10000:
                        n = 10000
                    r = rnd.random(l=length, n=n, g=genome_size).nucleotide_content(fi=genome_fa)
                    features += [f[:3] + [float(f[7])] for f in r if float(f[12]) <= length * N_FRACTION]
                    gc += [f[3] for f in features]
            
                for f in features:
                    if (f[3] >= bin_start and f[3] < bin_end):
                        out.write("{}\t{}\t{}\n".format(*f[:3]))
                        rcount += 1
                        if rcount >= count:
                            break
            
            if count != rcount:
                sys.stderr.write("not enough random sequences found for {} <= GC < {} ({} instead of {})\n".format(bin_start, bin_end, rcount, count))
    out.close()

    if del_size:
        os.unlink(genome_size)



class MatchedGcFasta(Fasta):
    """ 
    Generates a new Fasta object containing sequences randomly selected from
    the genome. These sequences are selected in such a way that the GC%
    distribution is similar to the GC% distribution og the input sequences.
   
    Required arg 'matchfile' is a BED or FASTA file
    
    Optional arg 'number' specifies the number of sequences to generate, 
    default is the number of input sequences/

    Returns a Fasta object
    
    """
    def __init__(self, matchfile, genome="hg19", number=None):
        # Create temporary files
        tmpbed = NamedTemporaryFile(dir=mytmpdir()).name
        tmpfasta = NamedTemporaryFile(dir=mytmpdir()).name
        
        # Create bed-file with coordinates of random sequences
        matched_gc_bedfile(tmpbed, matchfile, genome, number)
        
        # Convert track to fasta
        Genome(genome).track2fasta(tmpbed, fastafile=tmpfasta)

        # Initialize super Fasta object
        Fasta.__init__(self, tmpfasta)

        # Delete the temporary files
        os.remove(tmpbed)
        os.remove(tmpfasta)

class PromoterFasta(Fasta):
    """ 
    Generates a new Fasta object containing randomly selected promoters.
    A BED file of gene coordinates is used to extract sequences of a specified
    length upstream of the the TSS.
    
    Required arg 'genefile' is a file containing genes BED format (at least 6 
    columns including the strand information). 
    Required arg 'length' specifies the length 
    Required arg 'in' specifies the number of sequences to generate.

    Returns a Fasta object
    
    """
    def __init__(self, genefile, genome, length=None, n=None):
        length = int(length)

        # Create temporary files
        tmpbed = NamedTemporaryFile(dir=mytmpdir()).name
        tmpfasta = NamedTemporaryFile(dir=mytmpdir()).name
        
        # Create bed-file with coordinates of random sequences
        create_promoter_bedfile(tmpbed, genefile, length, n)
        
        # Convert track to fasta
        Genome(genome).track2fasta(tmpbed, fastafile=tmpfasta, stranded=True)

        # Initialize super Fasta object
        Fasta.__init__(self, tmpfasta)

        # Delete the temporary files
        os.remove(tmpbed)
        os.remove(tmpfasta)

class RandomGenomicFasta(Fasta):
    """ 
    Generates a new Fasta object containing randomly selected genomic regions.
    
    columns including the strand information). 
    Required arg 'length' specifies the length 
    Required arg 'in' specifies the number of sequences to generate.

    Returns a Fasta object
    
    """
    def __init__(self, genome, length=None, n=None):
        length = int(length)

        # Create temporary files
        tmpbed = NamedTemporaryFile(dir=mytmpdir()).name
        tmpfasta = NamedTemporaryFile(dir=mytmpdir()).name
        
        # Create bed-file with coordinates of random sequences
        create_random_genomic_bedfile(tmpbed, genome, length, n)
        
        # Convert track to fasta
        Genome(genome).track2fasta(tmpbed, fastafile=tmpfasta, stranded=True)

        # Initialize super Fasta object
        Fasta.__init__(self, tmpfasta)

        # Delete the temporary files
        os.remove(tmpbed)
        os.remove(tmpfasta)


