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
import os
import random
import re
import sys
from tempfile import NamedTemporaryFile
from random import choice

# External imports
import numpy as np
import pybedtools

# GimmeMotifs imports
from gimmemotifs import mytmpdir
from gimmemotifs.fasta import Fasta
from gimmemotifs.genome_index import track2fasta,get_random_sequences
from gimmemotifs.config import MotifConfig

def create_random_genomic_bedfile(out, index_dir, length, n):
    features = get_random_sequences(index_dir, n, length)

    # Write result to bedfile
    tmp = open(out, "w")
    for chrom,start,end in features:
        tmp.write("%s\t%d\t%d\n" % (chrom, start, end))
    tmp.flush()    

def create_promoter_bedfile(out, genefile, length, n):
    strand_map = {"+":True, "-":False, 1:True, -1:False, "1":True, "-1":False}

    features = []
    
    for line in open(genefile):
        if not line.startswith("track"):
            (chrom, start, end, _name, score, strand) = line[:-1].split("\t")[:6]
            start, end = int(start), int(end)
            strand= strand_map[strand]
            if strand:
                if start - length >= 0:
                    features.append([chrom, start - length, start, strand])
            else:
                features.append([chrom, end, end + length, strand])

    
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

def matched_gc_bedfile(bedfile, matchfile, genome, number):
    N_FRACTION = 0.1
    
    config = MotifConfig()
    index = os.path.join(config.get_index_dir(), genome)
    
    genome_size = os.path.join(index, "genome.size")
    genome_fa = os.path.join(index, "genome.fa")

    if not os.path.exists(genome_size) or not os.path.exists(genome_fa):
        raise RuntimeError("genome files not found, please re-index {} "  \
                "with a recent version of gimme index".format(genome))

    try:
        fa = Fasta(matchfile)
        gc = [(seq.upper().count("C") + seq.upper().count("G")) / len(seq) for seq in fa.seqs]
        lengths = [len(seq) for seq in fa.seqs]
    except Exception:
        try:
            # pylint: disable=unexpected-keyword-arg
            bed = pybedtools.BedTool(matchfile)
            gc = [float(x[4]) for x in bed.nucleotide_content(fi=genome_fa)]
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
    
    # pylint: disable=unexpected-keyword-arg
    r = rnd.random(l=length, n=number * 30, g=genome_size).nucleotide_content(fi=genome_fa)
    
    features = [f[:3] + [float(f[7])] for f in r if float(f[12]) <= length * N_FRACTION]
    gc = [f[3] for f in features]
    
    #sys.stderr.write("Done\n")
    for bin_start, bin_end, count in zip(bins[:-1], bins[1:], gc_hist):
        #sys.stderr.write("CG {}-{}\n".format(bin_start, bin_end))
        if count > 0:
            rcount = 0
            for f in features:
                if (f[3] >= bin_start and f[3] < bin_end):
                    out.write("{}\t{}\t{}\n".format(*f[:3]))
                    rcount += 1
                    if rcount >= count:
                        break

            if count != rcount:
                sys.stderr.write("not enough random sequences found for {} <= GC < {} ({} instead of {})\n".format(bin_start, bin_end, rcount, count))
    out.close()

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
        config = MotifConfig()
        index = os.path.join(config.get_index_dir(), genome)

        # Create temporary files
        tmpbed = NamedTemporaryFile(dir=mytmpdir()).name
        tmpfasta = NamedTemporaryFile(dir=mytmpdir()).name
        
        # Create bed-file with coordinates of random sequences
        matched_gc_bedfile(tmpbed, matchfile, genome, number)
        
        # Convert track to fasta
        track2fasta(index, tmpbed, tmpfasta)

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
    def __init__(self, genefile, index="/usr/share/gimmemotifs/genome_index/hg18", length=None, n=None):
        length = int(length)

        # Create temporary files
        tmpbed = NamedTemporaryFile(dir=mytmpdir()).name
        tmpfasta = NamedTemporaryFile(dir=mytmpdir()).name
        
        # Create bed-file with coordinates of random sequences
        create_promoter_bedfile(tmpbed, genefile, length, n)
        
        # Convert track to fasta
        track2fasta(index, tmpbed, tmpfasta, use_strand=True)

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
    def __init__(self, index="/usr/share/gimmemotifs/genome_index/hg18", length=None, n=None):
        length = int(length)

        # Create temporary files
        tmpbed = NamedTemporaryFile(dir=mytmpdir()).name
        tmpfasta = NamedTemporaryFile(dir=mytmpdir()).name
        
        # Create bed-file with coordinates of random sequences
        create_random_genomic_bedfile(tmpbed, index, length, n)
        
        # Convert track to fasta
        track2fasta(index, tmpbed, tmpfasta, use_strand=True)

        # Initialize super Fasta object
        Fasta.__init__(self, tmpfasta)

        # Delete the temporary files
        os.remove(tmpbed)
        os.remove(tmpfasta)


