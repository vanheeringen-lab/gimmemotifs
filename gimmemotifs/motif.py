# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
""" Core motif class """
# Python imports
import os
import re
import sys
import random
from math import log,sqrt
import subprocess as sp
from tempfile import NamedTemporaryFile
from warnings import warn
import six

from gimmemotifs import mytmpdir
from gimmemotifs.config import MotifConfig
from gimmemotifs.c_metrics import pwmscan

# External imports
try:
    import numpy as np
except ImportError:
    pass
import xxhash
import base64

class Motif(object):
    
    """
    Representation of a transcription factor binding motif.

    Examples
    --------

    >>> motif = Motif([[0,1,0,0], [0.5,0,0,0.5], [0,0,1,0]])
    >>> print(motif.to_pwm())
    >
    0   1   0   0
    0.5 0   0   0.5
    0   0   1   0
    >>> print(motif.to_consensus())
    CwG
    
    """
    
    PSEUDO_PFM_COUNT = 1000 # Jaspar mean
    PSEUDO_PWM = 1e-6

    def __init__(self, pfm=None):
        if pfm is None:
            pfm = []

        if len(pfm) > 0:
            if np.sum(pfm[0]) > 2:
                self.pfm = [list(x) for x in pfm]
                self.pwm = self.pfm_to_pwm(pfm)
            else:
                self.pwm = [list(x) for x in pfm]
                self.pfm = [[n * self.PSEUDO_PFM_COUNT for n in col] for col in pfm]
        else:
            self.pwm = []
            self.pfm = []
        
        self.factors = []
        self.seqs = []
        self.consensus = ""
        self.min_score = None
        self.max_score = None
        
        self.id = ""
        self.config = MotifConfig()
        self.seqlogo = self.config.get_seqlogo()

        self.nucs = "ACGT"
        
        self.iupac_rev = {
            'CG': 'S',
            'AG': 'R',
            'AT': 'W',
            'CT': 'Y',
            'GT': 'K',
            'AC': 'M',
            'CGT': 'B',
            'ACT': 'H',
            'AGT': 'D',
            'ACG': 'V'
        }

        self.iupac = {
            'A':['A'],
            'C':['C'],
            'G':['G'],
            'T':['T'],
            'S':['C','G'],
            'R':['A','G'],
            'W':['A','T'],
            'Y':['C','T'],
            'K':['G','T'],
            'M':['A','C'],
            'B':['C','G','T'],
            'H':['A','C','T'],
            'D':['A','G','T'],
            'V':['A','C','G'],
            'N':['A','C','G','T']
        }
        self.iupac_pwm = {
            'A':[1, 0, 0, 0],    
            'C':[0, 1, 0, 0],    
            'G':[0, 0, 1, 0],    
            'T':[0, 0, 0, 1],    
            'S':[0, 0.5, 0.5, 0],
            'R':[0.5, 0, 0.5, 0],
            'W':[0.5, 0, 0, 0.5],
            'Y':[0, 0.5, 0, 0.5],
            'K':[0, 0, 0.5, 0.5],
            'M':[0.5, 0.5, 0, 0],
            'B':[0, 0.33, 0.33, 0.33],
            'H':[0.33, 0.33, 0, 0.33],
            'D':[0.33, 0, 0.33, 0.33],
            'V':[0.33, 0.33, 0.33, 0],
            'N':[0.25, 0.25, 0.25, 0.25]
        }
    
    
    def __getitem__(self, x):
        """
        Take slice of a motif and return as new Motif instance.

        Returns
        -------
        motif : Motif instance
            Slice of the motif.
        """
        m = Motif()
        if self.pwm:
            m.pwm = self.pwm[x]
        if self.pfm:   
            m.pfm = self.pfm[x]
        if self.seqs:
            m.seqs = [seq[x] for seq in self.seqs]
        if self.consensus:
            m.consensus = self.consensus[x]
        return m
            
    def __len__(self):
        """
        Return the motif length.

        Returns
        -------
        len : int
            Motif length
        """
        return len(self.to_consensus())

    def __repr__(self):
        return "{}_{}".format(self.id, self.to_consensus())

    def information_content(self):
        ic = 0
        for row in self.pwm:
            ic += 2.0 + np.sum([row[x] * log(row[x])/log(2) for x in range(4) if row[x] > 0])
        return ic

    def pwm_min_score(self):
        """
        Return the minimum PWM score.

        Returns
        -------
        score : float
            Minimum PWM score.
        """
        if self.min_score is None:
            score = 0
            for row in self.pwm:
                score += log(min(row) / 0.25 + 0.01)
            self.min_score = score
        
        return self.min_score
   
    def pwm_max_score(self):
        """
        Return the maximum PWM score.

        Returns
        -------
        score : float
            Maximum PWM score.
        """
        if self.max_score is None:
            score = 0
            for row in self.pwm:
                score += log(max(row) / 0.25 + 0.01)
            self.max_score = score
        
        return self.max_score
    
    def score_kmer(self, kmer):
        if len(kmer) != len(self.pwm):
            raise Exception("incorrect k-mer length")
        
        score = 0.0
        d = {"A":0, "C":1, "G":2, "T":3}
        for nuc, row in zip(kmer.upper(), self.pwm):
            score += log(row[d[nuc]] / 0.25 + 0.01)

        return score

    def pfm_to_pwm(self, pfm, pseudo=0.001):
        return [[(x + pseudo)/(float(np.sum(row)) + pseudo * 4) for x in row] for row in pfm]

    def to_motevo(self):
        """
        Return motif formatted in MotEvo (TRANSFAC-like) format
        """
        m = "//\n"
        m += "NA {}\n".format(self.id)
        m += "P0\tA\tC\tG\tT\n"
        for i, row in enumerate(self.pfm):
            m += "{}\t{}\n".format(i, "\t".join([str(int(x)) for x in row]))
        m += "//"
        return m

    def to_transfac(self):
        m = "%s\t%s\t%s\n" % ("DE", self.id, "unknown")
        for i, (row, cons) in enumerate(zip(self.pfm, self.to_consensus())):
            m += "%i\t%s\t%s\n" % (i, "\t".join([str(int(x)) for x in row]), cons)
        m += "XX"
        return m

    def to_meme(self):
        motif_id = self.id.replace(" ", "_")
        m = "MOTIF %s\n" % motif_id
        m += "BL   MOTIF %s width=0 seqs=0\n"% motif_id
        m += "letter-probability matrix: alength= 4 w= %s nsites= %s E= 0\n" % (len(self), np.sum(self.pfm[0]))
        m +="\n".join(["\t".join(["%s" % x for x in row]) for row in self.pwm])
        return m

    def ic_pos(self, row1, row2=None):
        if row2 is None:
            row2 = [0.25,0.25,0.25,0.25]

        score = 0
        for a,b in zip(row1, row2):
            if a > 0:
                score += a * log(a / b) / log(2)
        return score

    def pcc_pos(self, row1, row2):
        mean1 = np.mean(row1)
        mean2 = np.mean(row2)

        a = 0
        x = 0
        y = 0
        for n1, n2 in zip(row1, row2):
            a += (n1 - mean1) * (n2 - mean2)
            x += (n1 - mean1) ** 2
            y += (n2 - mean2) ** 2
        
        if a == 0:
            return 0
        else:
            return a / sqrt(x * y)
    
    def rc(self):
        m = Motif()
        m.pfm = [row[::-1] for row in self.pfm[::-1]]
        m.pwm = [row[::-1] for row in self.pwm[::-1]]
        m.id = self.id + "_revcomp"
        return m

    def trim(self, edge_ic_cutoff=0.4):
        pwm = self.pwm[:]
        while len(pwm) > 0 and self.ic_pos(pwm[0]) < edge_ic_cutoff:
            pwm = pwm[1:]
            self.pwm = self.pwm[1:]
            self.pfm = self.pfm[1:]
        while len(pwm) > 0 and self.ic_pos(pwm[-1]) < edge_ic_cutoff:
            pwm = pwm[:-1]
            self.pwm = self.pwm[:-1]
            self.pfm = self.pfm[:-1]
        
        self.consensus = None 
        self.min_score = None
        self.max_score = None
        
        return self

    def consensus_scan(self, fa):
        regexp = "".join(["[" + "".join(self.iupac[x.upper()]) + "]" for x in self.to_consensusv2()])
        p = re.compile(regexp)
        matches = {}
        for name,seq in fa.items():
            matches[name] = [] 
            for match in p.finditer(seq):
                middle = (match.span()[1] + match.span()[0]) / 2
                matches[name].append(middle)
        return matches

    def pwm_scan(self, fa, cutoff=0.9, nreport=50, scan_rc=True):
        c = self.pwm_min_score() + (self.pwm_max_score() - self.pwm_min_score()) * cutoff        
        pwm = self.pwm
        matches = {}
        for name, seq in fa.items():
            matches[name] = [] 
            result = pwmscan(seq.upper(), pwm, c, nreport, scan_rc)
            for _,pos,_ in result:
                matches[name].append(pos)
        return matches
    
    def pwm_scan_all(self, fa, cutoff=0.9, nreport=50, scan_rc=True):
        c = self.pwm_min_score() + (self.pwm_max_score() - self.pwm_min_score()) * cutoff        
        pwm = self.pwm
        matches = {}
        for name, seq in fa.items():
            matches[name] = [] 
            result = pwmscan(seq.upper(), pwm, c, nreport, scan_rc)
            for score,pos,strand in result:
                matches[name].append((pos,score,strand))
        return matches

    def pwm_scan_score(self, fa, cutoff=0, nreport=1, scan_rc=True):
        c = self.pwm_min_score() + (self.pwm_max_score() - self.pwm_min_score()) * cutoff        
        pwm = self.pwm
        matches = {}
        for name, seq in fa.items():
            matches[name] = [] 
            result = pwmscan(seq.upper(), pwm, c, nreport, scan_rc)
            for score,_,_ in result:
                matches[name].append(score)
        return matches
            
    def pwm_scan_to_gff(self, fa, gfffile, cutoff=0.9, nreport=50, scan_rc=True, append=False):
        if append:
            out = open(gfffile, "a")
        else:    
            out = open(gfffile, "w")

        c = self.pwm_min_score() + (self.pwm_max_score() - self.pwm_min_score()) * cutoff        
        pwm = self.pwm

        strandmap = {-1:"-","-1":"-","-":"-","1":"+",1:"+","+":"+"}
        gff_line = ("{}\tpwmscan\tmisc_feature\t{}\t{}\t{:.3f}\t{}\t.\t"
                    "motif_name \"{}\" ; motif_instance \"{}\"\n")
        for name, seq in fa.items():
            result = pwmscan(seq.upper(), pwm, c, nreport, scan_rc)
            for score, pos, strand in result:
                out.write(gff_line.format( 
                    name, 
                    pos, 
                    pos + len(pwm), 
                    score, 
                    strandmap[strand], 
                    self.id, 
                    seq[pos:pos + len(pwm)]
                    ))
        out.close()

    def average_motifs(self, other, pos, orientation, include_bg=False):
        # xxCATGYT
        # GGCTTGYx
        # pos = -2
        pfm1 = self.pfm[:]
        pfm2 = other.pfm[:]

        if orientation < 0:
            pfm2 = [row[::-1] for row in pfm2[::-1]]
        
        pfm1_count = float(np.sum(pfm1[0]))
        pfm2_count = float(np.sum(pfm2[0]))
        
        if include_bg:
            if len(pfm1) > len(pfm2) + pos:
                pfm2 += [[pfm2_count / 4.0 for x in range(4)] for i in range(-(len(pfm1) - len(pfm2) - pos), 0)]
            elif len(pfm2) + pos > len(pfm1):
                pfm1 += [[pfm1_count / 4.0 for x in range(4)] for i in range(-(len(pfm2) - len(pfm1) + pos), 0)]
        
            if pos < 0:
                pfm1 = [[pfm1_count / 4.0 for x in range(4)] for i in range(-pos)] + pfm1
            elif pos > 0:
                pfm2 = [[pfm2_count / 4.0 for x in range(4)] for i in range(pos)] + pfm2
        
        else:
            if len(pfm1) > len(pfm2) + pos:
                pfm2 += [[pfm1[i][x] / pfm1_count * (pfm2_count) for x in range(4)] for i in range(-(len(pfm1) - len(pfm2) - pos), 0)]
            elif len(pfm2) + pos > len(pfm1):
                pfm1 += [[pfm2[i][x] / pfm2_count * (pfm1_count) for x in range(4)] for i in range(-(len(pfm2) - len(pfm1) + pos), 0)]
        
            if pos < 0:
                pfm1 = [[pfm2[i][x] / pfm2_count * (pfm1_count) for x in range(4)] for i in range(-pos)] + pfm1
            elif pos > 0:
                pfm2 = [[pfm1[i][x] / pfm1_count * (pfm2_count) for x in range(4)] for i in range(pos)] + pfm2

        pfm = [[a + b for a,b in zip(x,y)] for x,y in zip(pfm1, pfm2)]
        
        m = Motif(pfm)
        m.id = m.to_consensus()
        return m
    
    def other_ic_pos(self, row1, row2, bg=None):
        if bg is None:
            bg = [0.25,0.25,0.25,0.25]
        
        score = 0
        score_a = 0
        score_b = 0
        
        for a,b,pbg  in zip(row1, row2, bg):
            score += abs(a * log(a/pbg) / log(2) -  b * log(b/pbg) / log(2))
            score_a += a * log(a/pbg) / log(2)
            score_b += b * log(b/pbg) / log(2)

        return (score_a + score_b)/2 - score

    def pcc(self, pwm1, pwm2, pos):
        # xxCATGYT
        # GGCTTGYx
        # pos = -2
        pwm1 = pwm1[:]
        pwm2 = pwm2[:]

        na = []
        if pos > 0:
            na = pwm1[:pos]
            pwm1 = pwm1[pos:]
        elif pos < 0:
            na = pwm2[:-pos]
            pwm2 = pwm2[-pos:]
    
        if len(pwm1) > len(pwm2):
            na += pwm1[len(pwm2):]
            pwm1 = pwm1[:len(pwm2)]
        elif len(pwm2) > len(pwm1):
            na += pwm2[len(pwm1):]
            pwm2 = pwm2[:len(pwm1)]

        # Aligned parts of the motif
        score = 0
        for a,b in zip(pwm1, pwm2):
            score += self.pcc_pos(a, b)
        
        return score

    def ic(self, pwm1, pwm2, pos, bg=None, bg_factor=1):
        if bg is None:
            bg = [0.25,0.25,0.25,0.25]
        
        # xxCATGYT
        # GGCTTGYx
        # pos = -2
        pwm1 = pwm1[:]
        pwm2 = pwm2[:]

        na = []
        if pos > 0:
            na = pwm1[:pos]
            pwm1 = pwm1[pos:]
        elif pos < 0:
            na = pwm2[:-pos]
            pwm2 = pwm2[-pos:]
    
        if len(pwm1) > len(pwm2):
            na += pwm1[len(pwm2):]
            pwm1 = pwm1[:len(pwm2)]
        elif len(pwm2) > len(pwm1):
            na += pwm2[len(pwm1):]
            pwm2 = pwm2[:len(pwm1)]

        #print "COMPARE"    
        #print Motif(pwm1).to_consensus()
        #print Motif(pwm2).to_consensus()

        # Aligned parts of the motif
        score = 0
        for a,b in zip(pwm1, pwm2):
            score += self.ic_pos(a) +  self.ic_pos(b) - (self.ic_pos(a,b) + self.ic_pos(b,a))
        
        #print "SCORE: %s" % score
        # Parts aligned to the background
        for x in na:
            score += (self.ic_pos(x) + self.ic_pos(bg) - (self.ic_pos(x,bg) + self.ic_pos(bg,x)))  * bg_factor 
        
    #    print "SCORE WITH BG: %s" % score
        return score

    def other_ic(self, pwm1, pwm2, pos, bg=None, bg_factor=0.8):
        if bg is None:
            bg = [0.25,0.25,0.25,0.25]
        
        # xxCATGYT
        # GGCTTGYx
        # pos = -2
        pwm1 = pwm1[:]
        pwm2 = pwm2[:]

        na = []
        if pos > 0:
            na = pwm1[:pos]
            pwm1 = pwm1[pos:]
        elif pos < 0:
            na = pwm2[:-pos]
            pwm2 = pwm2[-pos:]
    
        if len(pwm1) > len(pwm2):
            na += pwm1[len(pwm2):]
            pwm1 = pwm1[:len(pwm2)]
        elif len(pwm2) > len(pwm1):
            na += pwm2[len(pwm1):]
            pwm2 = pwm2[:len(pwm1)]

        # Aligned parts of the motif
        score = 0
        for a,b in zip(pwm1, pwm2):
            score += self.other_ic_pos(a, b)
        
        for x in na:
            score += self.other_ic_pos(x, bg) * bg_factor 
        
        return score

    def matrix_ic(self, pwm1, pwm2, bg=None):
        if bg is None:
            bg = [0.25,0.25,0.25,0.25]
        
        # xxCATGYT
        # GGCTTGYx
        # pos = -2
        pwm1 = np.array(pwm1)
        pwm2 = np.array(pwm2)
        pwm2_rev = np.array([row[::-1] for row in pwm2[::-1]])
        bg = np.array(bg)
        
        a = pwm1 * np.log2(pwm1/bg)
        b = pwm2 * np.log2(pwm2/bg)
        
        b_rev = pwm2_rev * np.log2(pwm2_rev/bg)
        
        scores = []
        l1 = len(pwm1)
        l2 = len(pwm2)
        for pos in range(-(l2 -1), l1):
        
            pwm1_start, pwm2_start = 0,0
            pwm1_end, pwm2_end = l1,l2
            if pos > 0:
                pwm1_start = pos
                if l1 - pos > l2:
                    pwm1_end = l2 + pos
                elif l1 - pos < l2:
                    pwm2_end = l1 - pos
            elif pos < 0:
                pwm2_start = -pos
                if l2 + pos > l1:
                    pwm2_end = l1
                elif l2 + pos < l1:
                    pwm1_end = l2 + pos
            else:    
                if l2 > l1:
                    pwm2_end = l1
                elif l2 < l1:
                    pwm1_end = l2
            
            score = np.sum((np.sum(a[pwm1_start:pwm1_end],1) + 
                np.sum(b[pwm2_start:pwm2_end],1)) / 2 - 
                np.sum(
                    np.abs(a[pwm1_start:pwm1_end] - 
                        b[pwm2_start:pwm2_end]),1))
            scores.append([score, pos, 1])
            
            score = np.sum((np.sum(a[pwm1_start:pwm1_end],1) + 
                np.sum(b_rev[pwm2_start:pwm2_end],1)) / 2 - 
                np.sum(
                    np.abs(a[pwm1_start:pwm1_end] - 
                        b_rev[pwm2_start:pwm2_end]),1))
            scores.append([score, pos, -1])
        
        return sorted(scores, key=lambda x: x[0])[-1]


    def max_ic(self, other, revcomp=True, bg_factor=0.8):
        pwm1 = self.pwm
        pwm2 = other.pwm
        
        scores = []
        
        for i in range(-(len(pwm2) -1), len(pwm1)):
            scores.append((self.other_ic(pwm1, pwm2, i, bg_factor=bg_factor), i, 1))
        
        if revcomp:
            rev_pwm2 = [row[::-1] for row in pwm2[::-1]]
            for i in range(-(len(pwm2) -1), len(pwm1)):
                scores.append((self.other_ic(pwm1, rev_pwm2, i, bg_factor=bg_factor), i, -1))
    
        return sorted(scores, key=lambda x: x[0])[-1]

    def max_pcc(self, other, revcomp=True):
        pwm1 = self.pwm
        pwm2 = other.pwm
        
        scores = []
        
        for i in range(-(len(pwm2) -1), len(pwm1)):
            scores.append((self.pcc(pwm1, pwm2, i), i, 1))
        
        if revcomp:
            rev_pwm2 = [row[::-1] for row in pwm2[::-1]]
            for i in range(-(len(pwm2) -1), len(pwm1)):
                scores.append((self.pcc(pwm1, rev_pwm2, i), i, -1))
    
        #print scores
        return sorted(scores, key=lambda x: x[0])[-1]

    def _format_jaspar(self, version=1, header=True):
        rows = np.array(self.pwm).transpose()
        rows = [" ".join([str(x) for x in row]) for row in rows]
        if version == 2:
            rows = ["{} [{} ]".format(n,row) for n,row in zip(self.nucs, rows)]
        
        str_out = "\n".join(rows)
        if header:
            str_out = "\n".join([self.id, str_out])
        
        return str_out

    def to_consensus(self, precision=4):
        if not self.consensus:
            consensus = ""
            for row in self.pwm:
                weights = sorted(zip(["A","C","G","T"], row), key=lambda x: x[1])
                if round(weights[-1][1], precision) >= 0.5 and weights[-1][1] > 2 * weights[-2][1]:
                    consensus += weights[-1][0]
                elif round(weights[-1][1], precision) + round(weights[-2][1], precision) >= 0.75:
                    consensus +=  self.iupac_rev["".join(sorted([weights[-1][0], weights[-2][0]]))].lower()
                else:
                    consensus += "n"
            self.consensus = consensus
        
        return self.consensus

    def to_consensusv2(self):
        if self.consensus:
            return self.consensus
        else:
            consensus = ""
            for row in self.pwm:
                weights = sorted(zip(["A","C","G","T"], row), key=lambda x: x[1])
                if weights[-1][1] >= 0.5:
                    if weights[-2][1] >= 0.25:
                        consensus +=  self.iupac_rev["".join(sorted([weights[-1][0], weights[-2][0]]))]
                    else:        
                        consensus += weights[-1][0]
                elif weights[-1][1] + weights[-2][1] >= 0.75:
                    consensus +=  self.iupac_rev["".join(sorted([weights[-1][0], weights[-2][0]]))]
                elif weights[-1][1] + weights[-2][1] + weights[-3][1] >= 0.9:
                    consensus +=  self.iupac_rev["".join(sorted([weights[-1][0], weights[-2][0], weights[-3][0]]))]
                else:
                    consensus += "n"
            return consensus
            
    def to_pfm(self):
        if self.pfm:
            return ">%s\n%s" % (self.id, "\n".join(["\t".join(["%s" % x for x in row]) for row in self.pfm]))
        else:
            pfm = [[n * self.PSEUDO_PFM_COUNT for n in col] for col in self.pwm]
            return ">%s\n%s" % (self.id, "\n".join(["\t".join(["%s" % x for x in row]) for row in pfm]))

    def _pwm_to_str(self, precision=4):
        """Return string representation of pwm.

        Parameters
        ----------
        precision : int, optional, default 4
            Floating-point precision.

        Returns
        -------
        pwm_string : str
        """
        if not self.pwm:
            return ""
        
        fmt = "{{:.{:d}f}}".format(precision)
        return "\n".join(
                ["\t".join([fmt.format(p) for p in row])
                for row in self.pwm]
                )
   
    def hash(self):
        """Return hash of motif.

        Returns:
        hash : str
        """
        return xxhash.xxh64(self._pwm_to_str(3)).hexdigest()

    def to_pwm(self, precision=4, extra_str=""):
        """Return pwm as string.

        Parameters
        ----------
        precision : int, optional, default 4
            Floating-point precision.
        
        extra_str |: str, optional
            Extra text to include with motif id line.
        
        Returns
        -------
        motif_str : str
            Motif formatted in PWM format.
        """
        motif_id = self.id
        
        if extra_str:
            motif_id += "_%s" % extra_str

        if not self.pwm:
            self.pwm = [self.iupac_pwm[char]for char in self.consensus.upper()]

        return ">%s\n%s" % (
                motif_id, 
                self._pwm_to_str(precision)
            )

    def to_img(self, fname, fmt="EPS", add_left=0, seqlogo=None, height=6):
        """ Valid formats EPS, GIF, PDF, PNG """
        if not seqlogo:
            seqlogo = self.seqlogo
        if not seqlogo:
            raise ValueError("seqlogo not specified or configured")
        
        #TODO: split to_align function
        
        VALID_FORMATS = ["EPS", "GIF", "PDF", "PNG"]
        N = 1000
        fmt = fmt.upper()
        if not fmt in VALID_FORMATS:
            sys.stderr.write("Invalid motif format\n")
            return
        
        if fname[-4:].upper() == (".%s" % fmt):
            fname = fname[:-4]
        seqs = []
        if add_left == 0:
            seqs = ["" for i in range(N)]
        else:
            for nuc in ["A", "C", "T", "G"]:
                seqs += [nuc * add_left for i in range(N // 4)]

        for pos in range(len(self.pwm)):
            vals = [self.pwm[pos][0] * N]
            for i in range(1,4):
                vals.append(vals[i-1] + self.pwm[pos][i] * N)
            if vals[3] - N != 0:
                #print "Motif weights don't add up to 1! Error of %s%%" % ((vals[3] - n)/ n * 100)
                vals[3] = N
            for i in range(N):
                if i <= vals[0]:
                    seqs[i] += "A"
                elif i <= vals[1]:
                    seqs[i] += "C"
                elif i <= vals[2]:
                    seqs[i] += "G"
                elif i <= vals[3]:
                    seqs[i] += "T"
    
        f = NamedTemporaryFile(mode="w", dir=mytmpdir())
        for seq in seqs:
            f.write("%s\n" % seq)
        f.flush()
        makelogo = "{0} -f {1} -F {2} -c -a -h {3} -w {4} -o {5} -b -n -Y" 
        cmd = makelogo.format(
                              seqlogo, 
                              f.name, 
                              fmt, 
                              height,
                              len(self) + add_left, 
                              fname)
        sp.call(cmd, shell=True)
        
        # Delete tempfile
        #if os.path.exists(f.name):
        #    os.unlink(f.name)

    def randomize(self):
        random_pfm = [[c for c in row] for row in self.pfm]
        random.shuffle(random_pfm)
        m = Motif(pfm=random_pfm)
        m.id = "random"
        return m

    def randomize_dimer(self):
        l = len(self.pfm)
        random_pfm = []
        for _ in range(l / 2):
            pos = random.randint(0, l - 1)
            random_pfm += [[c for c in row] for row in self.pfm[pos:pos + 2]]
        m = Motif(pfm=random_pfm)
        m.id = "random"
        return m

def default_motifs():
    """Return list of Motif instances from default motif database."""
    config = MotifConfig()
    d = config.get_motif_dir()
    m = config.get_default_params()['motif_db']

    if not d or not m:
        raise ValueError("default motif database not configured")

    fname = os.path.join(d, m)
    with open(fname) as f:
        motifs = read_motifs(f)
    
    return motifs

def motif_from_align(align):
    width = len(align[0])
    nucs = {"A":0,"C":1,"G":2,"T":3}
    pfm =  [[0 for _ in range(4)] for _ in range(width)]
    for row in align:
        for i in range(len(row)):
            pfm[i][nucs[row[i]]] += 1
    m = Motif(pfm)
    m.align = align[:]
    return m

def motif_from_consensus(cons, n=12):
    width = len(cons)
    nucs = {"A":0,"C":1,"G":2,"T":3}
    pfm = [[0 for _ in range(4)] for _ in range(width)]
    m = Motif()
    for i,char in enumerate(cons):
        for nuc in m.iupac[char.upper()]:
            pfm[i][nucs[nuc]] = n / len(m.iupac[char.upper()])
    m = Motif(pfm)
    m.id = cons
    return m

def parse_motifs(motifs):
    """Parse motifs in a variety of formats to return a list of motifs.

    Parameters
    ----------

    motifs : list or str
        Filename of motif,  list of motifs or single Motif instance.

    Returns
    -------

    motifs : list
        List of Motif instances.
    """
    if isinstance(motifs, six.string_types):
        with open(motifs) as f:
            if motifs.endswith("pwm") or motifs.endswith("pfm"):
                motifs = read_motifs(f, fmt="pwm")
            elif motifs.endswith("transfac"):
                motifs = read_motifs(f, fmt="transfac")
            else: 
                motifs = read_motifs(f)
    elif isinstance(motifs, Motif):
        motifs = [motifs]
    else:
        if not isinstance(list(motifs)[0], Motif):
            raise ValueError("Not a list of motifs")
    
    return list(motifs)

def read_motifs(handle, fmt="pwm"):
    """ 
    Read motifs from a stream or file-like object.
    """

    if fmt.lower() == "pwm":
        motifs = _read_motifs_pwm(handle)
    if fmt.lower() == "transfac":
        motifs = _read_motifs_transfac(handle)
    if fmt.lower() == "xxmotif":
        motifs = _read_motifs_xxmotif(handle)
    if fmt.lower() == "align":
        motifs = _read_motifs_align(handle)
    if fmt.lower() == "jaspar":
        motifs = _read_motifs_jaspar(handle)
    
    if handle.name:
        base = os.path.splitext(handle.name)[0]
        map_file = base + ".motif2factors.txt"
        if os.path.exists(map_file):
            m2f = {}
            for line in open(map_file):
                try:
                    motif,factors = line.strip().split("\t")
                    m2f[motif] = factors.split(",")
                except:
                    pass
            for motif in motifs:
                if motif.id in m2f:
                    motif.factors = m2f[motif.id]
    return motifs

def _read_motifs_pwm(handle):
    p = re.compile(r'(\d+(\.\d+)?(e-\d+)?)\s+(\d+(\.\d+)?(e-\d+)?)\s+(\d+(\.\d+)?(e-\d+)?)\s+(\d+(\.\d+)?(e-\d+)?)')
    motifs = []
    pfm = []
    motif_id = ""
    seen_id = {}
    
    for n,line in enumerate(handle.readlines()):
        if line.startswith("#") or line.strip() == "":
            continue
        if line.startswith(">"):
            if pfm:
                motifs.append(Motif(pfm))
                motifs[-1].id = motif_id
                pfm = []
            motif_id = line.strip()[1:]
            seen_id[motif_id] = seen_id.get(motif_id, 0) + 1
            if seen_id.get(motif_id, 0) > 1:
                msg = "WARNING: multiple motifs with same id: {}\n".format(motif_id)
                sys.stderr.write(msg)
                motif_id += "_{}".format(seen_id[motif_id] - 1)
        
        else:
            m = p.search(line)
            if m:
                fractions =  [float(m.group(x)) for x in (1,4,7,10)]
                pfm.append(fractions)
            else:
                msg = "WARNING: can't parse line {}, ignoring:\n{}".format(n + 1, line)
                sys.stderr.write(msg)

    if len(pfm) > 0:
        motifs.append(Motif(pfm))
        motifs[-1].id = motif_id
            
    return motifs

def _read_motifs_jaspar(handle):
    p = re.compile("([ACGT])\s*\[?([^\]]+)\]?")
    motifs = []
    motif_id = ""
    pwm = {}
    for line in handle:
        line = line.strip()
        if len(line) == 0:
            continue
        
        if line.startswith(">"):
            motif_id = line[1:]
        if line[0] in "ACGT":
            m = p.search(line)
            nuc = m.group(1)
            counts = re.split(r'\s+', m.group(2).strip())
            pwm[nuc] = [float(x) for x in counts]
            if nuc == "T":
                motif = Motif(np.array([pwm[n] for n in "ACGT"]).transpose())
                motif.id = motif_id
                motifs.append(motif)
    if motif_id and motifs[-1].id != motif_id:
        motif = Motif(np.array([pwm[n] for n in "ACGT"]).transpose())
        motif.id = motif_id
        motifs.append(motif)
    
    return motifs
            

def _read_motifs_align(handle):
    motifs = []
    nucs = {"A":0,"C":1,"G":2,"T":3}
    motif_id = ""
    aligns = {}
    align = []
    for line in handle:
        if line.startswith(">"):
            if motif_id:
                aligns[motif_id] = align
            motif_id = line.strip()[1:]
            align = []
        else:
            align.append(line.strip())
    aligns[motif_id] = align

    for motif_id, align in aligns.items():

        width = len(align[0])
        pfm =  [[0 for _ in range(4)] for _ in range(width)]
        for row in align:
            for i in range(len(row)):
                pfm[i][nucs[row[i]]] += 1
        m = Motif(pfm)
        m.align = align[:]
        m.pfm = pfm[:]
        m.id = motif_id
        motifs.append(m)
    return motifs

def _read_motifs_xxmotif(handle):
    motifs = []
    
    line = handle.readline()
    while line:
        while line and not line.startswith("Motif"):
            line = handle.readline()
    
        if line:
            mid = line.split(":")[0]
            freqs = []
            for _ in range(4):
                line = handle.readline()
                freqs.append([float(x) for x in line.strip().split("\t")])

            pwm = np.array(freqs).transpose()
            motif = Motif(pwm)
            motif.id = mid.replace(" ", "_")
            motifs.append(motif)

    return motifs

def _read_motifs_transfac(handle):
    p = re.compile(r'\d+\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s*\w?')
    p_id = re.compile(r'^(NA|ID)\s+([^\s]+)')
    motifs = []
    pwm = []
    motif_id = ""
    for line in handle.readlines():
        m = p_id.search(line.strip())
        if m:
            motif_id = m.group(2)
        elif line.startswith("//"):
            motifs.append(Motif(pwm))
            motifs[-1].id = motif_id
            pwm = []
        elif p.search(line):
            m = p.search(line)
            pwm.append([float(x) for x in m.group(1,2,3,4)])
    
    # If there's only one matrix, and the format is not complete
    if len(pwm) != 0:
        motifs.append(Motif(pwm))

    return motifs

def motifs_to_meme(motifs):
    m = "MEME version 3.0\n\nALPHABET= ACGT\n\nstrands: + -\n\n"
    m += "Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n"

    for motif in motifs:
        m += motif.to_meme() + "\n"
    return m

def alignfile_to_motifs(fname):
    # this method should be deleted
    msg = "alignfile_to_motifs is deprecated, please use read_motifs"
    warn(msg, DeprecationWarning)

    return read_motifs(open(fname), fmt="align")    


def pwmfile_to_motifs(fname):
    # this method should be deleted
    msg = "pwmfile_to_motifs is deprecated, please use read_motifs"
    warn(msg, DeprecationWarning)
    
    return read_motifs(open(fname), fmt="pwm")    

def transfac_to_motifs(fname):
    # this method should be deleted
    msg = "transfac_to_motifs is deprecated, please use read_motifs"
    warn(msg, DeprecationWarning)
    
    return read_motifs(open(fname), fmt="transfac")    

def xxmotif_to_motifs(fname):
    # this method should be deleted
    msg = "xxmotif_to_motifs is deprecated, please use read_motifs"
    warn(msg, DeprecationWarning)
    
    return read_motifs(open(fname), fmt="xxmotif")    

if __name__ == "__main__":
    m = Motif()
