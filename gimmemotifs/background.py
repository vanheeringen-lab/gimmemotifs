# Copyright (c) 2009-2010 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
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

# Python imports
import os
import random
import re
import sys
from tempfile import NamedTemporaryFile
# External imports
from numpy import array,sum
# GimmeMotifs imports
from gimmemotifs.fasta import Fasta
from gimmemotifs.genome_index import track2fasta,get_random_sequences

def create_matched_genomic_bedfile(out, bedfile, genefile, length, m, match_chromosome):
    strand_map = {"+":True, "-":False, 1:True, -1:False, "1":True, "-1":False}
    data = {}

    features = []
    
    for line in open(genefile):
        if not line.startswith("track"):
            (chr, start, end, name, score, strand) = line[:-1].split("\t")[:6]
            start, end = int(start), int(end)
            strand= strand_map[strand]
            data.setdefault(chr,{})[name] = {"chr":chr,"start":start, "end":end, "strand":strand,"up":0, "down":0}


    for chr in data.keys():
        prev_gene = None
        for gene,stuff in sorted(data[chr].items(),key=lambda x: data[chr][x[0]]["start"]):
            if prev_gene:
                if stuff["start"] < data[chr][prev_gene]["end"]:
                    if stuff["end"] > data[chr][prev_gene]["end"]:
                        if stuff["strand"]:    
                            data[chr][gene]["up"] = None 
                        else:
                            data[chr][gene]["down"] = None
                        if data[chr][prev_gene]["strand"]:
                            data[chr][prev_gene]["down"] = None
                        else:
                            data[chr][prev_gene]["up"] = None
                        prev_gene = gene
                    else: 
                        data[chr][gene]["up"] = None 
                        data[chr][gene]["down"] = None 
                else:
                    dist = stuff["start"] - data[chr][prev_gene]["end"]
                    if dist < 0:
                        print "ERROHR"
                        sys.exit(1)
                    if stuff["strand"]:    
                        data[chr][gene]["up"] = dist 
                    else:
                        data[chr][gene]["down"] = dist
                    if data[chr][prev_gene]["strand"]:
                        data[chr][prev_gene]["down"] = dist
                    else:
                        data[chr][prev_gene]["up"] = dist 
                    
                    prev_gene = gene
            else:
                prev_gene = gene
    targets = {}
    target_length = {}
    for line in open(bedfile):
        if not line.startswith("track"):
            vals = line.strip().split("\t")
            (chr, start, end) = vals[:3]
            if not targets.has_key(chr):
                targets[chr] = []
                target_length[chr] = []
            targets[chr].append((int(end) + int(start))/2)
            if length:
                target_length[chr].append(length)
            else:
                target_length[chr].append(int(end) - int(start))
            
    result = []
    data_all = []
    for chr in data.keys():
        for gene in data[chr].keys():
            data_all.append(data[chr][gene])
    
    for chr in targets.keys():
        if data.has_key(chr):
            vals = data[chr].values()
            locs = sorted([(x["start"], x["end"], x["strand"]) for x in data[chr].values()],key=lambda x: x[0])
            closest = [[None, None] for x in targets[chr]]
            j = 0
            for i,pos in enumerate(targets[chr]):
                while j >= len(locs) or (j > 0 and locs[j][1] > pos):
                    j -= 1
                while j < len(locs) and locs[j][1] < pos:
                    j += 1
                if j > 0:
                    j -= 1
                    if locs[j][2]:
                        if not closest[i][1] or pos - locs[j][1] < closest[i][1]:
                            closest[i][1] = pos - locs[j][1]
                    else:
                        if not closest[i][0] or pos - locs[j][1] < closest[i][0]:
                            closest[i][0] = pos - locs[j][1]
                    j += 1
                if j < len(locs):
                    if locs[j][0] <= pos and locs[j][1] >= pos:
                        closest[i] = [0,0]
                    else:
                        if locs[j][2]:
                            if not closest[i][0] or locs[j][0] - pos < closest[i][0]:
                                closest[i][0] = locs[j][0] - pos
                        else:
                            if not closest[i][1] or locs[j][0] - pos < closest[i][0]:
                                closest[i][1] = locs[j][0] - pos
                else: 
                    pass
            
            for i,c in enumerate(closest):
                genes = []
                if c[0] == 0 and c[1] == 0:
                    
                    genes = []
                    if match_chromosome:
                        genes = data[chr].values()
                    if len(genes) < m:
                        genes = data_all[:]

                    for x in range(m):
                        gene = random.choice(genes)
                        pos = random.randint(gene["start"], gene["end"])
                        features.append((chr, pos - target_length[chr][i]/2, pos + target_length[chr][i]/2))
                else:
                    if not c[1] or (c[0] and c[0] <= c[1]):
                        genes = []
                        if match_chromosome:
                            genes = [x for x in data[chr].values() if x["up"] > 2 * c[0]]
                        if len(genes) < m:
                            genes = [x for x in data_all if x["up"] and  x["up"] > 2 * c[0]]
                        if len(genes) >= m:
                            for x in range(m):
                                gene = random.choice(genes)
                                if gene["strand"]:
                                    pos = gene["start"] - c[0]
                                else:
                                    pos = gene["end"] + c[0]
                                
                                pos = random.randint(pos - int(0.1 * c[0]), pos + int(0.1 * c[0]))
                                features.append((gene["chr"], pos - target_length[chr][i]/2, pos + target_length[chr][i]/2))
                        else:
                            pass
                                
                    elif not c[0] or (c[1] and c[1] <= c[0]):
                        genes = []
                        if match_chromosome:
                            genes = [x for x in data[chr].values() if x["down"] > 2 * c[1]]
                        if len(genes) < m:
                            genes = [x for x in data_all if x["down"] and (x["down"] > 2 * c[1])]
                        if len(genes) >= m:
                            for x in range(m):
                                gene = random.choice(genes)
                                if gene["strand"]:
                                    pos = gene["end"] + c[1]
                                else:
                                    pos = gene["start"] - c[1]
                                pos = random.randint(pos - int(0.1 * c[1]), pos + int(0.1 * c[1]))
                                features.append((gene["chr"], pos - target_length[chr][i]/2, pos + target_length[chr][i]/2))
                        else:
                            pass
    
    # Write result to temporary bedfile
    tmp = open(out, "w")
    for chr, start, end in sorted(features, key=lambda x: x[0]):
        tmp.write("%s\t%s\t%s\n" % (chr, start, end))
    tmp.flush()

def create_random_genomic_bedfile(out, index_dir, length, n):
    features = get_random_sequences(index_dir, n, length)

    # Write result to bedfile
    tmp = open(out, "w")
    for chrom,start,end in features:
        tmp.write("%s\t%d\t%d\n" % (chrom, start, end))
    tmp.flush()    

def create_promoter_bedfile(out, genefile, length, n):
    strand_map = {"+":True, "-":False, 1:True, -1:False, "1":True, "-1":False}
    data = {}

    features = []
    
    for line in open(genefile):
        if not line.startswith("track"):
            (chr, start, end, name, score, strand) = line[:-1].split("\t")[:6]
            start, end = int(start), int(end)
            strand= strand_map[strand]
            if strand:
                if start - length >= 0:
                    features.append([chr, start - length, start, strand])
            else:
                features.append([chr, end, end + length, strand])

    
    if n < len(features):
        features = random.sample(features, n)
    else:
        sys.stdout.write("Too few promoters to generate %s random promoters! Just using all of them." % n)


    # Write result to temporary bedfile
    tmp = open(out, "w")
    for chr, start, end, strand in sorted(features, key=lambda x: x[0]):
        tmp.write("%s\t%s\t%s\t0\t0\t%s\n" % (chr, start, end, {True:"+",False:"-"}[strand]))
    tmp.flush()


class MarkovFasta(Fasta):
    """ 
    Generates a new Fasta object containing sequences using a 1st order Markov
    model, based on the input sequences. By default 10 times as many sequences 
    will be generated with the same length as the input sequences.

    Required arg 'fasta' is a Fasta object
    Optional arg 'length' can be used to generate sequences of a different length
    Optional arg 'multiply' specifies the number of sequences to generate, based
    Optional arg 'k' specifies the order of the Markov model, default is 1 for 1st
    order

    Returns a Fasta object
    
    Example:
        
        f = Fasta("input.fa")
        random_fa = MarkovFasta(f, multiply = 5)
        for id,seq in random_fa.items():
            print seq
    
    """
    
    def __init__(self, fasta, length=None, multiply=10, k=1, matrix_only=False):
        
        
        self.k = k

        # Initialize super Fasta object
        Fasta.__init__(self)
        
        # Initialize Markov transition matrix
        self._initialize_matrices(fasta.seqs, k=k)

        if matrix_only:
            return
        
        c = 0
        for seq in fasta.seqs:
            for i in range(multiply):
                id = "random_Markov%s_%s" % (k,c)
                if length:
                    random_seq = self._generate_sequence(length)
                else:
                    random_seq = self._generate_sequence(len(seq))
                self.add(id, random_seq)    
                c += 1

    def _initialize_matrices(self, seqs, k=1, skip_bad=True, alphabet=['A','C','G','T'], bad="n"):
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
            s = sum(array(v.values()))
            for x,y in v.items():
                v[x] = y / s
        
        self.init = {}
        total = float(sum(array(lettercount.values())))
        for k,v in lettercount.items():
            self.init[k] = v / total
            
    def _generate_sequence(self, l):
        sequence = list(self._weighted_random(self.init.items()))
        for i in range(l - self.k):
            sequence.append(self._weighted_random(self.trans["".join(sequence[-self.k:])].items()))
        return "".join(sequence)

    def _weighted_random(self, l):
        n = random.uniform(0,1)
        for item, weight in l:
            if n < weight:
                break    
            else:
                n -= weight
        return item

class MatchedGenomicFasta(Fasta):
    """ 
    Generates a new Fasta object containing sequences randomly selected from
    a set of genomic sequences. These sequences are selected in such a way
    that the general genomic distribution relative to the TSS of genes is
    similar as the input sequences. By default the random sequences will also
    have a similar distribution among chromosomes, unless the optional arg 
    'match_chromosome' is set to 'False'.
    
    Required arg 'bedfile' is a file containing genomic coordinates in BED format
    Required arg 'genefile' is a file containing genes BED format (at least 6 
    columns including the strand information). 

    Optional arg 'length' can be used to generate sequences of a different length
    Optional arg 'multiply' specifies the number of sequences to generate, based
    on the number of input sequences.

    Returns a Fasta object
    
    """
    def __init__(self, bedfile, genefile, index="/usr/share/gimmemotifs/genome_index/hg18", length=None, multiply=10, match_chromosome=True):
        self.match_chromosome = match_chromosome
        length = int(length)

        # Create temporary files
        tmpbed = NamedTemporaryFile().name
        tmpfasta = NamedTemporaryFile().name
        
        # Create bed-file with coordinates of random sequences
        create_matched_genomic_bedfile(tmpbed, bedfile, genefile, length, multiply, self.match_chromosome)
        
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
        tmpbed = NamedTemporaryFile().name
        tmpfasta = NamedTemporaryFile().name
        
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
        tmpbed = NamedTemporaryFile().name
        tmpfasta = NamedTemporaryFile().name
        
        # Create bed-file with coordinates of random sequences
        create_random_genomic_bedfile(tmpbed, index, length, n)
        
        # Convert track to fasta
        track2fasta(index, tmpbed, tmpfasta, use_strand=True)

        # Initialize super Fasta object
        Fasta.__init__(self, tmpfasta)

        # Delete the temporary files
        os.remove(tmpbed)
        os.remove(tmpfasta)


