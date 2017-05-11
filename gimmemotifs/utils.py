# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

""" Odds and ends that for which I didn't (yet) find another place """
from __future__ import print_function

# Python imports
import os
import re
import sys
import hashlib
import mmap
import random
import tempfile
from math import log
from subprocess import Popen
from tempfile import NamedTemporaryFile

# External imports
from scipy import special
import numpy as np
import pybedtools

# gimme imports
from gimmemotifs.fasta import Fasta
from gimmemotifs.plot import plot_histogram
from gimmemotifs.genome_index import track2fasta
from gimmemotifs.rocmetrics import ks_pvalue

lgam = special.gammaln

def phyper_single(k, good, bad, N):

    return np.exp(lgam(good+1) - lgam(good-k+1) - lgam(k+1) + lgam(bad+1) - lgam(bad-N+k+1) - lgam(N-k+1) - lgam(bad+good+1) + lgam(bad+good-N+1) + lgam(N+1))

def phyper(k, good, bad, N):
    """ Current hypergeometric implementation in scipy is broken, so here's the correct version """
    pvalues = [phyper_single(x, good, bad, N) for x in range(k + 1, N + 1)]
    return np.sum(pvalues)

def divide_file(fname, sample, rest, fraction, abs_max):
    with open(fname) as f:
        lines = f.readlines()
    #random.seed()
    random.shuffle(lines)

    x = int(fraction * len(lines))
    if x > abs_max:
        x = abs_max

    tmp = tempfile.NamedTemporaryFile(mode="w", delete=False)

    # Fraction as sample
    for line in lines[:x]:
        tmp.write(line)
    tmp.flush()

    # Make sure it is sorted for tools that use this information (MDmodule)
    stdout,stderr = Popen("sort -k4gr %s > %s" % (tmp.name, sample), shell=True).communicate()

    tmp.close()

    if stderr:
        print("Something went wrong.\nstdout: {}\nstderr; {}".format(
                stdout,
                stderr))
        sys.exit()

    # Rest
    f = open(rest, "w")
    for line in lines[x:]:
        f.write(line)
    f.close()

    #if os.path.exists(tmp.name):
    #    os.unlink(tmp.name)
    return x, len(lines[x:])    
    
def divide_fa_file(fname, sample, rest, fraction, abs_max):
    fa = Fasta(fname)
    ids = fa.ids[:]

    x = int(fraction * len(ids))
    if x > abs_max:
        x = abs_max

    sample_seqs = random.sample(ids, x)

    # Rest
    f_sample = open(sample, "w")
    f_rest = open(rest, "w")
    for name,seq in fa.items():
        if name in sample_seqs:
            f_sample.write(">%s\n%s\n" % (name, seq))
        else:
            f_rest.write(">%s\n%s\n" % (name, seq))
    f_sample.close()
    f_rest.close()
    
    return x, len(ids[x:])    

def write_equalwidth_bedfile(bedfile, width, outfile):
    """Read input from <bedfile>, set the width of all entries to <width> and 
    write the result to <outfile>.
    Input file needs to be in BED or WIG format."""

    BUFSIZE = 10000
    f = open(bedfile)
    out = open(outfile, "w")
    lines = f.readlines(BUFSIZE)
    line_count = 0
    while lines:
        for line in lines:
            line_count += 1
            if not line.startswith("#") and not line.startswith("track") and not line.startswith("browser"):
                vals = line.strip().split("\t")
                try:
                    start, end = int(vals[1]), int(vals[2])
                except ValueError:
                    print("Error on line %s while reading %s. Is the file in BED or WIG format?" % (line_count, bedfile))
                    sys.exit(1)

                start = (start + end) // 2 - (width // 2)
                # This shifts the center, but ensures the width is identical... maybe not ideal
                if start < 0:
                    start = 0
                end = start + width
                # Keep all the other information in the bedfile if it's there
                if len(vals) > 3:
                    out.write("%s\t%s\t%s\t%s\n" % (vals[0], start, end, "\t".join(vals[3:])))
                else:
                    out.write("%s\t%s\t%s\n" % (vals[0], start, end))
        lines = f.readlines(BUFSIZE)
    
    out.close()
    f.close()

class MotifMatch(object):
    def __init__(self, seq, name, instance, start, end, strand, score):
        self.sequence = seq
        self.motif_name = name
        self.motif_instance = instance
        self.start = start
        self.end = end
        self.strand = strand
        self.score = score

class MotifResult(object):
    def __init__(self):
        self.raw_output = ""
        self.datetime = ""
        self.command = ""
        self.fastafile = ""
        self.params = {}
        self.program = ""
        self.feature = ""
    
        self.sequences = {} 
        self.motifs = {}
    
        self.matches = {}

    def to_gff(self, gb_format=False):
        p = re.compile(r'([\w_]+):(\d+)-(\d+)')
        
        gff_output = ""    
        for seq, d in self.matches.items():
            for mms in d.values():
                for mm in mms:
                    print_seq = seq
                    (start, end) = (mm.start, mm.end)
                    if gb_format:
                        m = p.match(seq)
                        if m:
                            print_seq = m.group(1)
                            start = int(start) + int(m.group(2)) - 1
                            end = int(end) + int(m.group(2)) - 1
                    gff_output += "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %  (print_seq,
                        self.program,
                        self.feature,
                        start,
                        end,
                        mm.score,
                        mm.strand,
                        ".",
                        "motif_name \"%s\" ; motif_instance \"%s\"" % (mm.motif_name, mm.motif_instance))
        return gff_output[:-1]

    def seqn(self):
        return len(self.sequences.keys())

def parse_gff(gff_file):
    mr = MotifResult()
    total = 0
    
    BUFSIZE = 10000000
    with open(gff_file) as f:
        while 1:
            lines = f.readlines(BUFSIZE)
            if not lines:
                break
            for line in lines:
                vals = line.strip().split("\t")
                if len(vals) == 9:
                    (seq, program, feature, start, end, score, strand, bla, extra) = vals
            
                    (motif_name, motif_instance) = map(str.strip, extra.split(";"))
                    motif_name = motif_name.split(" ")[1][1:-1]
                    motif_instance = motif_instance.split(" ")[1][1:-1]
    
                    mr.sequences[seq] = 1
    
                    if not motif_name in mr.motifs:
                        mr.motifs[motif_name] = {}
                    if not seq in mr.motifs[motif_name]:
                        mr.motifs[motif_name][seq] = 0
                    mr.motifs[motif_name][seq] += 1
                else:
                    sys.stderr.write("Error parsing line in %s\n%s\n" % (gff_file, line))
            total += len(lines)
    return mr

def calc_motif_enrichment(sample, background, mtc=None, len_sample=None, len_back=None):
    """Calculate enrichment based on hypergeometric distribution"""
    
    INF = "Inf"


    if mtc not in [None, "Bonferroni", "Benjamini-Hochberg", "None"]:
        raise RuntimeError("Unknown correction: %s" % mtc)

    sig = {}
    p_value  = {}
    n_sample = {}
    n_back = {}
    
    if not(len_sample):
        len_sample = sample.seqn()
    if not(len_back):
        len_back = background.seqn()

    for motif in sample.motifs.keys():
        p = "NA"
        s = "NA"
        q = len(sample.motifs[motif])
        m = 0
        if(background.motifs.get(motif)):
            m = len(background.motifs[motif])
            n = len_back - m
            k = len_sample
            p = phyper(q - 1, m, n, k) 
            if p != 0:
                s = -(log(p)/log(10))
            else:
                s = INF
        else:
            s = INF
            p = 0.0

        sig[motif] = s
        p_value[motif] = p
        n_sample[motif] = q
        n_back[motif] = m
    
    if mtc == "Bonferroni":
        for motif in p_value.keys():
            if  p_value[motif] != "NA":
                p_value[motif] = p_value[motif] * len(p_value.keys())
                if p_value[motif] > 1:
                    p_value[motif] = 1
    elif mtc == "Benjamini-Hochberg":
        motifs = sorted(p_value.keys(), key=lambda x: -p_value[x])
        l = len(p_value)
        c = l
        for m in motifs:
            if  p_value[m] != "NA":
                p_value[m] = p_value[m] * l / c 
            c -= 1

    return (sig, p_value, n_sample, n_back)

def gff_enrichment(sample, background, numsample, numbackground, outfile):
    data_sample = parse_gff(sample)
    data_bg = parse_gff(background)
    (e,p,ns,nb) = calc_motif_enrichment(data_sample, data_bg, "Benjamini-Hochberg", numsample, numbackground)
    
    out = open(outfile, "w")
    out.write("Motif\tSig\tp-value\t# sample\t# background\tEnrichment\n")
    for m in e.keys():
        if nb[m] > 0:
            enrich = (ns[m] / float(numsample)) / (nb[m] / float(numbackground))
            out.write("%s\t%s\t%s\t%s\t%s\t%0.3f\n" % (m, e[m], p[m], ns[m], nb[m], enrich))
        else:
            out.write("%s\t%s\t%s\t%s\t%s\tInf\n" % (m, e[m], p[m], ns[m], nb[m]))
    out.close()

def is_valid_bedfile(bedfile, columns=6):
    f = open(bedfile)
    for i, line in enumerate(f.readlines()):
        if not (line.startswith("browser") or line.startswith("track")):
            vals = line.split("\t")
            
            # Gene file should be at least X columns
            if len(vals) < columns:    
                sys.stderr.write("Error in line %s: we need at least %s columns!\n" % (i, columns))
                return False

            # Check coordinates
            try:
                int(vals[1]), int(vals[2])
            except ValueError:
                sys.stderr.write("Error in line %s: coordinates in column 2 and 3 need to be integers!\n" % (i))
                return False
    
            if columns >= 6:
                # We need the strand
                if vals[5] not in ["+", "-"]:
                    sys.stderr.write("Error in line %s: column 6 (strand information) needs to be + or -" % (i))
                    return False
    
    f.close()
    return True

def median_bed_len(bedfile):
    f = open(bedfile)
    l = []
    for i, line in enumerate(f.readlines()):
        if not (line.startswith("browser") or line.startswith("track")):
            vals = line.split("\t")
            try:
                l.append(int(vals[2]) - int(vals[1]))
            except ValueError:
                sys.stderr.write("Error in line %s: coordinates in column 2 and 3 need to be integers!\n" % (i))
                sys.exit(1)
    f.close()
    return np.median(l)

def motif_localization(fastafile, motif, width, outfile, cutoff=0.9):
    NR_HIST_MATCHES = 100

    matches = motif.pwm_scan(Fasta(fastafile), cutoff=cutoff, nreport=NR_HIST_MATCHES)
    if len(matches) > 0:
        ar = []
        for a in matches.values():
            ar += a
        matches = np.array(ar)
        p = ks_pvalue(matches, width - len(motif))
        plot_histogram(
                matches - width / 2 + len(motif) / 2, 
                outfile, 
                xrange=(-width / 2, width / 2), 
                breaks=21, 
                title="%s (p=%0.2e)" % (motif.id, p), 
                xlabel="Position")
        return motif.id, p
    else:
        return motif.id, 1.0

def parse_cutoff(motifs, cutoff, default=0.9):
    """ Provide either a file with one cutoff per motif or a single cutoff
        returns a hash with motif id as key and cutoff as value
    """
    
    cutoffs = {}
    if os.path.isfile(str(cutoff)):
        for i,line in enumerate(open(cutoff)):
            if line != "Motif\tScore\tCutoff\n":
                try:
                    motif,_,c = line.strip().split("\t")
                    c = float(c)
                    cutoffs[motif] = c
                except Exception as e:
                    sys.stderr.write("Error parsing cutoff file, line {0}: {1}\n".format(e, i + 1))
                    sys.exit(1)
    else:
        for motif in motifs:
            cutoffs[motif.id] = float(cutoff)
    
    for motif in motifs:
        if not motif.id in cutoffs:
            sys.stderr.write("No cutoff found for {0}, using default {1}\n".format(motif.id, default))
            cutoffs[motif.id] = default
    return cutoffs
                
def _treesort(order, nodeorder, nodecounts, tree):
    # From the Pycluster library, Michiel de Hoon
        # Find the order of the nodes consistent with the hierarchical clustering
    # tree, taking into account the preferred order of nodes.
    nNodes = len(tree)
    nElements = nNodes + 1
    neworder = np.zeros(nElements)
    clusterids = np.arange(nElements)
    for i in range(nNodes):
        i1 = tree[i].left
        i2 = tree[i].right
        if i1 < 0:
            order1 = nodeorder[-i1-1]
            count1 = nodecounts[-i1-1]
        else:
            order1 = order[i1]
            count1 = 1
        if i2 < 0:
            order2 = nodeorder[-i2-1]
            count2 = nodecounts[-i2-1]
        else:
            order2 = order[i2]
            count2 = 1
        # If order1 and order2 are equal, their order is determined
        # by the order in which they were clustered
        if i1 < i2:
            if order1 < order2:
                increase = count1
            else:
                increase = count2
            for j in range(nElements):
                clusterid = clusterids[j]
                if clusterid == i1 and order1 >= order2:
                    neworder[j] += increase
                if clusterid == i2 and order1 < order2:
                    neworder[j] += increase
                if clusterid == i1 or clusterid == i2:
                    clusterids[j] = -i-1
        else:
            if order1 <= order2:
                increase = count1
            else:
                increase = count2
            for j in range(nElements):
                clusterid = clusterids[j]
                if clusterid == i1 and order1 > order2:
                    neworder[j] += increase
                if clusterid == i2 and order1 <= order2:
                    neworder[j] += increase
                if clusterid == i1 or clusterid == i2:
                    clusterids[j] = -i-1
    return np.argsort(neworder)

def number_of_seqs_in_file(fname):
    try:
        fa = Fasta(fname)
        return len(fa)
    except Exception:
        pass

    try:
        bed = pybedtools.BedTool(fname)
        return len([x for x in bed])
    except Exception:
        pass

    sys.stderr.write("unknown filetype {}\n".format(fname))
    sys.exit(1)

def get_seqs_type(seqs):
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
    elif isinstance(seqs, str) or isinstance(seqs, unicode):
        if os.path.isfile(seqs):
            try:
                Fasta(seqs)
                return "fastafile"
            except:
                pass
            try:
                with open(seqs) as f:
                    for line in f.readlines():
                        line = line.strip()
                        if not line.startswith("#"):
                            break
                
                if region_p.search(line):
                    return "regionfile"
                else:
                    vals = line.split("\t")
                    if len(vals) >= 3:
                        _, _ = int(vals[1]), int(vals[2])
                        return "bedfile"
                raise ValueError("unknown type")
            except:
                raise ValueError("unknown type")
        else:
            raise ValueError("no file found with name {}".format(seqs))
    else:
        raise ValueError("unknown type {}".format(type(seqs).__name__))

def as_fasta(seqs, index_dir=None):
    ftype = get_seqs_type(seqs)
    if ftype == "fasta":
        return seqs
    elif ftype == "fastafile":
        return Fasta(seqs)
    else:
        if index_dir is None:
            raise ValueError("need index_dir / genome to convert to FASTA")

        tmpfa = NamedTemporaryFile()
        
        if ftype == "bedfile":
            track2fasta(index_dir, seqs, tmpfa.name) 
        else:

            if ftype == "regionfile":
                with open(seqs) as f:
                    seqs = [l.strip() for l in f.readlines()]
            tmpbed = NamedTemporaryFile(mode="w")
            for seq in seqs:
                vals = re.split(r'[:-]', seq)
                tmpbed.write("{}\t{}\t{}\n".format(*vals))
            tmpbed.flush()
            track2fasta(index_dir, tmpbed.name, tmpfa.name) 
        return Fasta(tmpfa.name)

def file_checksum(fname):
    """Return md5 checksum of file.

    Note: only works for files < 4GB.

    Parameters
    ----------
    filename : str
        File used to calculate checksum.

    Returns
    -------
        checkum : str
    """
    size = os.path.getsize(fname)
    with open(fname, "r+") as f:
        checksum = hashlib.md5(mmap.mmap(f.fileno(), size)).hexdigest()
    return checksum
