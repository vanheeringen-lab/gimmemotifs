# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

""" Module for simple motif prediction. Currently not in use! """

from gimmemotifs.fasta import *
import sys
from numpy import *
from string import maketrans
from tempfile import NamedTemporaryFile,mkdtemp
from subprocess import *
import os
from gimmemotifs.motif import *
from gimmemotifs.cluster import *


def nmer_predict(fastafile):
	from tempfile import NamedTemporaryFile,mkdtemp
	from gimmemotifs.fasta import Fasta
	from numpy import sum,histogram
	from subprocess import Popen,PIPE
	from gimmemotifs.motif import Motif,motif_from_align
	from gimmemotifs.cluster import cluster_motifs 
	from string import maketrans

	def rc(seq):
		t = maketrans("ATCG", "TAGC")
		return seq[::-1].translate(t)
	
	f = Fasta(fastafile)
	nmer = {}
	N = {6:4, 8:3,10:2,12:1}
	tmp = NamedTemporaryFile()
	abs_cutoff = len(f.items()) / 100.0 * 2 
	for check_n,cutoff in N.items():
		for id,seq in f.items():
			for i in range(len(seq) - check_n):
				n = seq[i: i + check_n]
				nmer.setdefault(n.upper(), []).append(i)

	for n,pos in nmer.items():
		if len(pos) > abs_cutoff:
			hist = histogram(pos, bins=9, range=(0,200))[0]	
			if sum(hist[3:6]) > sum(hist[0:3] * N[len(n)]) and  sum(hist[3:6]) > sum(hist[7:]) *  N[len(n)]:
				tmp.write(">%s\n" % n)
				for char in n:
					w = []
					for x in  ["A", "C", "G", "T"]:
						if x == char:
							w.append(len(pos))
						else:
							w.append(0)

					tmp.write("\t".join([str(x) for x in w]) + "\n")
	
	
	tmp.flush()
	tmpname = tmp.name
	
	tree = cluster_motifs(tmpname, "subtotal", "ed", "mean", False, threshold=-0.1, include_bg=False)	
	clusters = tree.getResult()

	def refine_by_scanning(motifs, fastafile):
		
		tmp_gff = NamedTemporaryFile()
		file_in = NamedTemporaryFile()
		for m in motifs:
			file_in.write("%s\n" % m.to_pfm())
		file_in.flush()
		
		cmd = "pwmscan.py -i %s -p %s -c 0.8 > %s" % (fastafile, file_in.name, tmp_gff.name)
		p = Popen(cmd, shell=True)
		stdout,stderr = p.communicate()

		aligns = {}
		for line in open(tmp_gff.name):	
			vals = line.strip().split("\t")
			motif,instance = [x.split(" ")[1].replace('"', "") for x in vals[8].split(" ; ")]
		
			if vals[6] == "+":
				aligns.setdefault(motif,[]).append(instance.upper())
			else:
				aligns.setdefault(motif,[]).append(rc(instance.upper()))

		tmp_out = NamedTemporaryFile()
		
		refined_motifs = []
		for id,align in aligns.items():
			if len(align) > 10:
				motif = motif_from_align(align)
				refined_motifs.append(motif)
		
		return refined_motifs
	
	motifs = refine_by_scanning([x[0] for x in clusters], fastafile)
	tmp4 = NamedTemporaryFile()
	for m in motifs:
		tmp4.write("%s\n" % m.to_pfm())
	tmp4.flush()


	motifs = []
	tree = cluster_motifs(tmp4.name, "total", "wic", "mean", True, threshold=0.95, include_bg=True)	
	clusters = tree.getResult()
	for i, (cluster,members) in enumerate(clusters):
		cluster.id = "Nmer_%s" % (i + 1)
		motifs.append(cluster)
	
	refined_motifs = refine_by_scanning(motifs, fastafile)
	for i,m in enumerate(refined_motifs):
		m.id = "WannaMotif_%s" % (i + 1)
	
	return refined_motifs, "", ""	

if __name__ == "__main__":

	motifs = nmer_predict(sys.argv[1])[0]
	for motif in motifs:
		print motif.to_pfm()
