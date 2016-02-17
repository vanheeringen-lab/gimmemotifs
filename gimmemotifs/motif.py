# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

""" Core motif class """

# Python imports
import re
import os
import sys
import random
from math import log,sqrt,log
from tempfile import NamedTemporaryFile
from config import *
from subprocess import *

from gimmemotifs import mytmpdir

# External imports
try:
    from numpy import mean,sum
    import numpy as np
except:
    pass

class Motif:
    PSEUDO_PFM_COUNT = 1000 # Jaspar mean
    PSEUDO_PWM = 1e-6

    def __init__(self, pfm=[]):
        if len(pfm) > 0:
            if sum(pfm[0]) > 2:
                self.pfm = [list(x) for x in pfm]
                self.pwm = self.pfm_to_pwm(pfm)
            else:
                self.pwm = [list(x) for x in pfm]
                self.pfm = [[n * self.PSEUDO_PFM_COUNT for n in col] for col in pfm]
        else:
            self.pwm = []
            self.pfm = []
        self.seqs = []
        self.consensus = ""
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
        return len(self.to_consensus())

    def __repr__(self):
        return "%s\t%s" % (self.id, self.to_consensus())

    def information_content(self):
        ic = 0
        for row in self.pwm:
            ic += 2.0 + sum([row[x] * log(row[x])/log(2) for x in range(4) if row[x] > 0])
        return ic

    def pwm_min_score(self):
        score = 0
        for row in self.pwm:
            score += log(min(row) / 0.25 + 0.01)
        return score
   
    def score_kmer(self, kmer):
        if len(kmer) != len(self.pwm):
            raise Exception, "incorrect k-mer length"
        
        score = 0.0
        d = {"A":0, "C":1, "G":2, "T":3}
        for nuc, row in zip(kmer.upper(), self.pwm):
            score += log(row[d[nuc]] / 0.25 + 0.01)

        return score

    def pwm_max_score(self):
        score = 0
        for row in self.pwm:
            score += log(max(row) / 0.25 + 0.01)
        return score

    def pfm_to_pwm(self, pfm, pseudo=0.001):
        return [[(x + pseudo)/(float(sum(row)) + pseudo * 4) for x in row] for row in pfm]

    def to_transfac(self):
        m = "%s\t%s\t%s\n" % ("DE", self.id, "unknown")
        for i, (row, cons) in enumerate(zip(self.pfm, self.to_consensus())):
            m += "%i\t%s\t%s\n" % (i, "\t".join([str(int(x)) for x in row]), cons)
        m += "XX"
        return m

    def to_meme(self):
        id = self.id.replace(" ", "_")
        m = "MOTIF %s\n" % id
        m += "BL   MOTIF %s width=0 seqs=0\n" % id
        m += "letter-probability matrix: alength= 4 w= %s nsites= %s E= 0\n" % (len(self), sum(self.pfm[0]))
        m +="\n".join(["\t".join(["%s" % x for x in row]) for row in self.pwm])
        return m

    def ic_pos(self, row1, row2=[0.25,0.25,0.25,0.25]):
        score = 0
        for a,b in zip(row1, row2):
            if a > 0:
                score += a * log(a / b) / log(2)
        return score

    def pcc_pos(self, row1, row2):
        mean1 = mean(row1)
        mean2 = mean(row2)

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

    def consensus_scan(self, fa):
        regexp = "".join(["[" + "".join(self.iupac[x.upper()]) + "]" for x in self.to_consensusv2()])
        p = re.compile(regexp)
        matches = {}
        for id,seq in fa.items():
            matches[id] = [] 
            for match in p.finditer(seq):
                middle = (match.span()[1] + match.span()[0]) / 2
                matches[id].append(middle)
        return matches

    def pwm_scan(self, fa, cutoff=0.9, nreport=50, scan_rc=True):
        from gimmemotifs.c_metrics import pwmscan
        c = self.pwm_min_score() + (self.pwm_max_score() - self.pwm_min_score()) * cutoff        
        pwm = self.pwm
        strandmap = {"+":"+",1:"+","1":"+","-":"-",-1:"-","-1":"-"}
        matches = {}
        for id, seq in fa.items():
            matches[id] = [] 
            result = pwmscan(seq.upper(), pwm, c, nreport, scan_rc)
            for score,pos,strand in result:
                matches[id].append(pos)
        return matches
    
    def pwm_scan_all(self, fa, cutoff=0.9, nreport=50, scan_rc=True):
        from gimmemotifs.c_metrics import pwmscan
        c = self.pwm_min_score() + (self.pwm_max_score() - self.pwm_min_score()) * cutoff        
        pwm = self.pwm
        strandmap = {"+":"+",1:"+","1":"+","-":"-",-1:"-","-1":"-"}
        matches = {}
        for id, seq in fa.items():
            matches[id] = [] 
            #sys.stderr.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}".format(seq.upper(), pwm, c, nreport, scan_rc))
            result = pwmscan(seq.upper(), pwm, c, nreport, scan_rc)
            for score,pos,strand in result:
                matches[id].append((pos,score,strand))
        return matches

    def pwm_scan_score(self, fa, cutoff=0, nreport=1, scan_rc=True):
        from gimmemotifs.c_metrics import pwmscan
        c = self.pwm_min_score() + (self.pwm_max_score() - self.pwm_min_score()) * cutoff        
        pwm = self.pwm
        strandmap = {"+":"+",1:"+","1":"+","-":"-",-1:"-","-1":"-"}
        matches = {}
        for id, seq in fa.items():
            matches[id] = [] 
            result = pwmscan(seq.upper(), pwm, c, nreport, scan_rc)
            for score,pos,strand in result:
                matches[id].append(score)
        return matches
            
    def pwm_scan_to_gff(self, fa, gfffile, cutoff=0.9, nreport=50, scan_rc=True, append=False):
        #print "received", gfffile, cutoff, nreport, scan_rc, append
        from gimmemotifs.c_metrics import pwmscan
        if append:
            out = open(gfffile, "a")
        else:    
            out = open(gfffile, "w")

        c = self.pwm_min_score() + (self.pwm_max_score() - self.pwm_min_score()) * cutoff        
        pwm = self.pwm

        strandmap = {-1:"-","-1":"-","-":"-","1":"+",1:"+","+":"+"}
        for id, seq in fa.items():
            result = pwmscan(seq.upper(), pwm, c, nreport, scan_rc)
            for score, pos, strand in result:
                out.write("%s\tpwmscan\tmisc_feature\t%s\t%s\t%s\t%s\t.\tmotif_name \"%s\" ; motif_instance \"%s\"\n" % 
                    (id, pos, pos + len(pwm), score, strandmap[strand], self.id, seq[pos:pos + len(pwm)]))
        out.close()

    def average_motifs(self, other, pos, orientation, include_bg=False):
        # xxCATGYT
        # GGCTTGYx
        # pos = -2
        pfm1 = self.pfm[:]
        pfm2 = other.pfm[:]

        if orientation < 0:
            pfm2 = [row[::-1] for row in pfm2[::-1]]
        
        pfm1_count = 0.0
        for i in pfm1[0]:
            pfm1_count += i

        pfm2_count = 0.0
        for i in pfm2[0]:
            pfm2_count += i

        if orientation < 0:
            [row[::-1] for row in pfm2[::-1]]    

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
    
    def other_ic_pos(self, row1, row2, bg=[0.25,0.25,0.25,0.25]):
        score = 0
        score_a = 0
        score_b = 0
        max = 2
        for a,b,pbg  in zip(row1, row2, bg):
            #print  a  * log(a/pbg) /log(2),  b *  log(b /pbg) / log(2)
            #print a * log(a/pbg) / log(2),  b * log(b/pbg) / log(2), abs(a * log(a/pbg) / log(2) -  b * log(b/pbg) / log(2))
            score += abs(a * log(a/pbg) / log(2) -  b * log(b/pbg) / log(2))
            score_a += a * log(a/pbg) / log(2)
            score_b += b * log(b/pbg) / log(2)

        #print score_a, score_b, score
        #print
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

    def ic(self, pwm1, pwm2, pos, bg=[0.25,0.25,0.25,0.25], bg_factor=1):
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

    def other_ic(self, pwm1, pwm2, pos, bg=[0.25,0.25,0.25,0.25], bg_factor=0.8):
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

    def matrix_ic(self, pwm1, pwm2, bg=[0.25,0.25,0.25,0.25], bg_factor=0.8):
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
            
            score = np.sum((np.sum(a[pwm1_start:pwm1_end],1) + np.sum(b[pwm2_start:pwm2_end],1)) / 2 - np.sum(np.abs(a[pwm1_start:pwm1_end] - b[pwm2_start:pwm2_end]),1))
            scores.append([score, pos, 1])
            
            score = np.sum((np.sum(a[pwm1_start:pwm1_end],1) + np.sum(b_rev[pwm2_start:pwm2_end],1)) / 2 - np.sum(np.abs(a[pwm1_start:pwm1_end] - b_rev[pwm2_start:pwm2_end]),1))
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
            str_out = "\n".join(self.id, str_out)
        
        return str_out

    def to_consensus(self):
        if self.consensus:
            return self.consensus
        else:
            consensus = ""
            for row in self.pwm:
                weights = sorted(zip(["A","C","G","T"], row), key=lambda x: x[1])
                if weights[-1][1] >= 0.5 and weights[-1][1] > 2 * weights[-2][1]:
                    consensus += weights[-1][0]
                elif weights[-1][1] + weights[-2][1] >= 0.75:
                    consensus +=  self.iupac_rev["".join(sorted([weights[-1][0], weights[-2][0]]))].lower()
                else:
                    consensus += "n"
            return consensus

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

    def to_pwm(self, extra_str=""):
        id = self.id
        if extra_str:
            id += "_%s" % extra_str
        if self.pwm:
            return ">%s\n%s" % (id, "\n".join(["\t".join(["%s" % x for x in row]) for row in self.pwm]))
        
        pseudocount = 0.8
        nucs = ["A","C","G","T"]

        pwm = [self.iupac_pwm[char]for char in self.consensus.upper()]
        return ">%s\n%s" % (id, "\n".join(["\t".join(["%s" % x for x in row]) for row in pwm]))

    def to_img(self, file, format="EPS", add_left=0, seqlogo=None, height=6):
        """ Valid formats EPS, GIF, PDF, PNG """
        if not seqlogo:
            seqlogo = self.seqlogo
        if not seqlogo:
            raise ValueError, "seqlogo not specified or configured"
        
        #TODO: split to_align function
        
        VALID_FORMATS = ["EPS", "GIF", "PDF", "PNG"]
        N = 1000
        format = format.upper()
        if not format in VALID_FORMATS:
            sys.stderr.write("Invalid motif format\n")
            return
        
        if file[-4:].upper() == (".%s" % format):
            file = file[:-4]
        seqs = []
        if add_left == 0:
            seqs = ["" for i in range(N)]
        else:
            for nuc in ["A", "C", "T", "G"]:
                seqs += [nuc * add_left for i in range(N / 4)]

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
    
        f = NamedTemporaryFile(dir=mytmpdir())
        for seq in seqs:
            f.write("%s\n" % seq)
        f.flush()
        makelogo = "{0} -f {1} -F {2} -c -a -h {3} -w {4} -o {5} -b -n -Y" 
        cmd = makelogo.format(
                              seqlogo, 
                              f.name, 
                              format, 
                              height,
                              len(self) + add_left, 
                              file)
        call(cmd, shell=True)
        
        # Delete tempfile
        #if os.path.exists(f.name):
        #    os.unlink(f.name)

    def stats(self, fg_fa, bg_fa, logger=None):
        from gimmemotifs.rocmetrics import MNCP, ROC_AUC, max_enrichment, fraction_fdr, score_at_fdr, enr_at_fdr
        from gimmemotifs.fasta import Fasta
        from gimmemotifs.utils import ks_pvalue
        from numpy import array,std
        from math import log

        try:
            stats = {}
            fg_result = self.pwm_scan_all(fg_fa, cutoff=0.0, nreport=1, scan_rc=True)
            bg_result = self.pwm_scan_all(bg_fa, cutoff=0.0, nreport=1, scan_rc=True)
           
            pos = []
            for x in fg_result.values():
                if len(x):
                    pos.append(x[0][1])
                else:
                    pos.append(-100)

            neg = []
            for x in bg_result.values():
                if len(x):
                    neg.append(x[0][1])
                else:
                    neg.append(-100)
           
            stats["mncp"] = MNCP(pos, neg)
            stats["roc_auc"] = ROC_AUC(pos, neg)
            x,y = max_enrichment(pos, neg)
            stats["maxenr"] = x
            stats["scoreatmaxenr"] = y
            stats["fraction"] = fraction_fdr(pos, neg)
            stats["score_fdr"] = score_at_fdr(pos, neg)
            stats["enr_fdr"] = enr_at_fdr(pos, neg)
            stats["cutoff_fdr"] = (stats["score_fdr"] - self.pwm_min_score()) / (self.pwm_max_score() - self.pwm_min_score())

            pos = [x[0][0] for x in fg_result.values() if len(x)]
            p = ks_pvalue(pos, max(pos))
            stats["ks"] = p
            if p > 0:
                stats["ks_sig"] = -log(p)/log(10)
            else:
                stats["ks_sig"] = "Inf"
        

            return stats
        except Exception as e:
            raise
            #e = sys.exc_info()[0]
            msg = "Error calculating stats of {0}, error {1}".format(self.id, str(e))
            if logger:
                logger.error(msg)
            else:
                print msg


    def randomize(self):
        random_pfm = [[c for c in row] for row in self.pfm]
        random.shuffle(random_pfm)
        m = Motif(pfm=random_pfm)
        m.id = "random"
        return m

    def randomize_dimer(self):
        l = len(self.pfm)
        random_pfm = []
        for i in range(l / 2):
            pos = random.randint(0, l - 1)
            random_pfm += [[c for c in row] for row in self.pfm[pos:pos + 2]]
        m = Motif(pfm=random_pfm)
        m.id = "random"
        return m

def motif_from_align(align):
    width = len(align[0])
    nucs = {"A":0,"C":1,"G":2,"T":3}
    pfm =  [[0 for x in range(4)] for x in range(width)]
    for row in align:
        for i in range(len(row)):
            pfm[i][nucs[row[i]]] += 1
    m = Motif(pfm)
    m.align = align[:]
    return m

def motif_from_consensus(cons, n=12):
    width = len(cons)
    nucs = {"A":0,"C":1,"G":2,"T":3}
    pfm = [[0 for x in range(4)] for x in range(width)]
    m = Motif()
    for i,char in enumerate(cons):
        for nuc in m.iupac[char.upper()]:
            pfm[i][nucs[nuc]] = n / len(m.iupac[char.upper()])
    m = Motif(pfm)
    m.id = cons
    return m


def pwmfile_to_motifs(file):
    p = re.compile(r'(\d+(\.\d+)?(e-\d+)?)\s+(\d+(\.\d+)?(e-\d+)?)\s+(\d+(\.\d+)?(e-\d+)?)\s+(\d+(\.\d+)?(e-\d+)?)')
    motifs = []
    pfm = []
    id = ""
    for n,line in enumerate(open(file).readlines()):
        if line.startswith("#") or line.strip() == "":
            continue
        if line.startswith(">"):
            if pfm:
                motifs.append(Motif(pfm))
                motifs[-1].id = id
                pfm = []
            id = line.strip()[1:]
        else:
            m = p.search(line)
            if m:
                fractions =  [float(m.group(x)) for x in (1,4,7,10)]
                pfm.append(fractions)
            else:
                sys.stderr.write("WARNING: can't parse line %s, ignoring:\n%s" % (n + 1, line))

    if len(pfm) > 0:
        motifs.append(Motif(pfm))
        motifs[-1].id = id
            
    return motifs

def alignfile_to_motifs(file):
    motifs = []
    nucs = {"A":0,"C":1,"G":2,"T":3}
    pwm = []
    id = ""
    aligns = {}
    align = []
    for line in open(file):
        if line.startswith(">"):
            if id:
                aligns[id] = align
            id = line.strip()[1:]
            align = []
        else:
            align.append(line.strip())
    aligns[id] = align

    for id, align in aligns.items():

        width = len(align[0])
        pfm =  [[0 for x in range(4)] for x in range(width)]
        for row in align:
            for i in range(len(row)):
                pfm[i][nucs[row[i]]] += 1
        total = float(len(align))
        m = Motif(pfm)
        m.align = align[:]
        m.pfm = pfm[:]
        m.id = id
        motifs.append(m)
    return motifs

def xxmotif_to_motifs(fname):
    motifs = []
    
    f = open(fname)
    line = f.readline()
    while line:
        while line and not line.startswith("Motif"):
            line = f.readline()
    
        if line:
            mid = line.split(":")[0]
            freqs = []
            for i in range(4):
                line = f.readline()
                freqs.append([float(x) for x in line.strip().split("\t")])

            pwm = np.array(freqs).transpose()
            motif = Motif(pwm)
            motif.id = mid.replace(" ", "_")
            motifs.append(motif)

    return motifs

def transfac_to_motifs(file):
    p = re.compile(r'\d+\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s*\w?')
    motifs = []
    pwm = []
    id = ""
    for line in open(file).readlines():
        if line.startswith("ID"):
            if pwm:
                motifs.append(Motif(pwm))
                motifs[-1].id = id
                pwm = []
            try:
                id = line.strip().split(" ")[-1].split("\t")[-1]
            except:
                id = line.strip().split("\t")[1]
        elif p.search(line):
            m = p.search(line)
            pwm.append([float(x) for x in m.group(1,2,3,4)])
    
    motifs.append(Motif(pwm))
    motifs[-1].id = id
            
    return motifs

def motifs_to_meme(motifs):
    m = "MEME version 3.0\n\nALPHABET= ACGT\n\nstrands: + -\n\n"
    m += "Background letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n"

    for motif in motifs:
        m += motif.to_meme() + "\n"
    return m

if __name__ == "__main__":
    m = Motif()
