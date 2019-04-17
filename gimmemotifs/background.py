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
from gimmemotifs.config import CACHE_DIR, BG_TYPES, MotifConfig
from gimmemotifs.utils import number_of_seqs_in_file

logger = logging.getLogger("gimme.background")

def create_background_file(outfile, bg_type, fmt='fasta', length=None,  genome=None, inputfile=None, number=10000):
    """
    Create a background file for motif analysis.

    Parameters
    ----------
    outfile : str
        Name of the output file.
    bg_type : str
        Type of background (gc, genomic, random or promoter). 
    fmt : str, optional
        Either 'fasta' or 'bed'.
    length : int, optional
        Length of the generated sequences, is determined from the inputfile if not 
        given.
    genome : str, optional
    inputfile : str, optional
    number : int, optional
    """
    fmt = fmt.lower()
    if fmt in ["fa", "fsa"]:
        fmt = "fasta"

    if bg_type not in BG_TYPES:
        print("The argument 'type' should be one of: %s" % (",".join(BG_TYPES)))
        sys.exit(1)

    if fmt == "bed" and bg_type == "random":
        print("Random background can only be generated in FASTA format!")
        sys.exit(1)
        
    if bg_type == "gc" and not inputfile:
        print("need a FASTA formatted input file for background gc")
        sys.exit(1)
    
    # GimmeMotifs configuration for file and directory locations
    config = MotifConfig()
        
    # Genome index location for creation of FASTA files
    if bg_type in ["gc", "genomic", "promoter"] and fmt == "fasta":
        if genome is None:
            print("Need a genome to create background file")
            sys.exit(1)
        Genome(genome)

    # Gene definition
    fname = Genome(genome).filename
    gene_file = fname.replace(".fa", ".annotation.bed.gz")
    if not gene_file:
        gene_file = os.path.join(config.get_gene_dir(), "{}.bed".format(genome))
    
    if bg_type in ["promoter"]:
        if not os.path.exists(gene_file):
            print("Could not find a gene file for genome {}".format(genome))
            print("Did you use the --annotation flag for genomepy?")
            print("Alternatively make sure there is a file called {}.bed in {}".format(genome, config.get_gene_dir()))
            sys.exit(1)

    # Number of sequences
    if number is None:
        if inputfile:
            number = number_of_seqs_in_file(inputfile)
            logger.info("Using {} of background sequences based on input file".format(number))
        else:
            number = 10000
            logger.info("Number of background sequences not specified, using 10,000 sequences")
    
    if bg_type == "random":
        f = Fasta(inputfile)
        m = MarkovFasta(f, n=number, k=1)
        m.writefasta(outfile)
    elif bg_type == "gc":
        if fmt == "fasta":
            m = MatchedGcFasta(inputfile, genome, number=number)
            m.writefasta(outfile)
        else:
            matched_gc_bedfile(outfile, inputfile, genome, number)
    elif bg_type == "promoter":
        if fmt == "fasta":
            m = PromoterFasta(gene_file, genome, length=length, n=number)
            m.writefasta(outfile)
        else:
            create_promoter_bedfile(outfile, gene_file, length, number)
    elif bg_type == "genomic":
        if fmt == "fasta":
            m = RandomGenomicFasta(genome, length, number)
            m.writefasta(outfile)
        else:
            create_random_genomic_bedfile(outfile, genome, length, number)

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

def create_gc_bin_index(genome, fname, min_bin_size=100):
    """Create index of GC content for a genome.

    Parameters
    ----------
    genome : str
        Genome name.
    fname : str
        Name of the index file.
    min_bin_size : int
        Minimum bin size (default 100). Warning: setting to a small value
        will result in a very large index file!
    """
    logger.info("Creating index for genomic GC frequencies.")
    g = Genome(genome)
    fasta = g.filename
    sizes = g.props["sizes"]["sizes"]

    with NamedTemporaryFile() as tmp:
        pybedtools.BedTool().window_maker(g=sizes, w=min_bin_size).nucleotide_content(fi=fasta).saveas(tmp.name)
        df = pd.read_csv(tmp.name, sep="\t", usecols=[0,1,2,4,9])

    cols = ["chrom", "start", "end", "w{}".format(min_bin_size), "n{}".format(min_bin_size)]
    for t in (2,5):
        df["w{}".format(min_bin_size * t)] = df.iloc[:,3].rolling(t, min_periods=t).mean()
        df["n{}".format(min_bin_size * t)] = df.iloc[:,4].rolling(t, min_periods=t).sum()
        cols += ["w{}".format(min_bin_size * t), "n{}".format(min_bin_size * t)]
    
    df.columns = cols
    df.reset_index()[cols].to_feather(fname)

def gc_bin_bedfile(bedfile, genome, number, l=200, bins=None, random_state=None, min_bin_size=100):
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

    fname = os.path.join(CACHE_DIR, "{}.gcfreq.{}.feather".format(os.path.basename(genome), min_bin_size))
    if not os.path.exists(fname):
        if not os.path.exists(CACHE_DIR):
            os.mkdirs(CACHE_DIR)
        create_gc_bin_index(genome, fname, min_bin_size=min_bin_size)
    
    df = pd.read_feather(fname)

    if l >= min_bin_size:
        col = "w{}".format(((l + min_bin_size // 2)// min_bin_size) * min_bin_size)
    else:
        logger.warn("For regions smaller than {}nt, GC% will not be exact".format(min_bin_size))
        col = "w{}".format(min_bin_size)

    if not col in df.columns:
        df[col] = df.iloc[:,3].rolling(l//min_bin_size, min_periods=l//min_bin_size).mean()
        df[col.replace("w","n")] = df.iloc[:,3].rolling(l//min_bin_size, min_periods=l//min_bin_size).sum()

    df = df[df[col.replace("w","n")] < 0.1 * l]
    n = number // len(bins)
    
    with open(bedfile, "w") as f:
        pass

    with open(bedfile, 'a') as f:
        for b_start, b_end in bins:
            #print(df.head())
            #print(col, b_start, b_end)
            df_bin = df[(df[col] > b_start) & (df[col] <= b_end)]
            #print(df_bin)
            if df_bin.shape[0] > 0:
                df_bin["start"] = df_bin["end"] - l 
                df_bin = df_bin[df_bin["start"] > 0]
                df_bin = df_bin.sample(n, replace=True, random_state=random_state)
                df_bin["bin"] = "{:.2f}-{:.2f}".format(b_start, b_end)
                df_bin[["chrom", "start", "end", "bin"]].to_csv(f, sep="\t", header=False, index=False)

def matched_gc_bedfile(bedfile, matchfile, genome, number, min_bin_size=100):
    """Create a BED file with GC% matched to input file.
    
    Parameters
    ----------
    bedfile : str
        Name of the output BED file.
    matchfile : str
        Name of input file (BED or FASTA format)
    genome : str
        Genome name.
    number : int
        Number of sequences to retrieve.
    """
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
            gc = np.array([float(x[fields + 1]) for x in bed.nucleotide_content(fi=genome_fa)])
            lengths = np.array([x.length for x in bed])
            gc = [round(x, 2) for x in gc]
        except:
            sys.stderr.write("Please provide input file in BED or FASTA format\n")
            raise
    
    # Get the median length of the sequences
    length = int(np.median(lengths))
    if np.std(lengths) > length * 0.05:
        sys.stderr.write("Sequences do not seem to be of equal length.\n")
        sys.stderr.write("GC% matched sequences of the median length ({}) will be created\n".format(length))

    bins = [(0.0, 0.2), (0.8, 1)]
    for b in np.arange(0.2, 0.799, 0.05):
        bins.append((b, b + 0.05))
   
    fraction = number / len(gc)
    gc = np.array(gc)
    #print("GC", gc)
    bin_count = []
    for b_start, b_end in bins:
        bin_count.append(int(np.sum((gc > round(b_start, 2)) & (gc <= round(b_end, 2))) * fraction))

    rest = number - sum(bin_count)
    for i in range(rest):
        bin_count[i] += 1
   
    nseqs = max(bin_count) * len(bins)
    
    with NamedTemporaryFile(delete=False) as tmp:
        gc_bin_bedfile(tmp.name, genome, nseqs, l=length, bins=bins, random_state=None, min_bin_size=min_bin_size)
        df = pd.read_csv(tmp.name, sep="\t", names=["chrom", "start", "end", "bin"])
        #print(tmp.name)
    with open(bedfile, "w") as f:
        pass
    with open(bedfile, 'a') as f:
        for (b_start, b_end), n in zip(bins, bin_count):
            #print(b_start, b_end, n)
            b = "{:.2f}-{:.2f}".format(b_start, b_end)
            df.loc[df["bin"] == b, ["chrom", "start", "end"]].sample(n).to_csv(f, sep="\t", header=False, index=False)

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


