# Copyright (c) 2009-2010 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

""" Odds and ends that for which I didn't (yet) find another place """

# Python imports
import os
import re
import sys
import random
import tempfile
from math import log
from string import strip
from subprocess import Popen, PIPE
# External imports
import numpy
from scipy import special
from gimmemotifs import tools
from gimmemotifs.fasta import *

lgam = special.gammaln

def phyper_single(k, good, bad, N):

	return numpy.exp(lgam(good+1) - lgam(good-k+1) - lgam(k+1) + lgam(bad+1) - lgam(bad-N+k+1) - lgam(N-k+1) - lgam(bad+good+1) + lgam(bad+good-N+1) + lgam(N+1))

def phyper(k, good, bad, N):
	""" Current hypergeometric implementation in scipy is broken, so here's the correct version """
	pvalues = [phyper_single(x, good, bad, N) for x in range(k + 1, N + 1)]
	return numpy.sum(pvalues)

def divide_file(file, sample, rest, fraction, abs_max):
	lines = open(file).readlines()
	random.seed()
	random.shuffle(lines)

	x = int(fraction * len(lines))
	if x > abs_max:
		x = abs_max

	tmp = tempfile.NamedTemporaryFile()

	# Fraction as sample
	for line in lines[:x]:
		tmp.write(line)
	tmp.flush()

	# Make sure it is sorted for tools that use this information (MDmodule)
	stdout,stderr = Popen("sort -k4gr %s > %s" % (tmp.name, sample), shell=True).communicate()

	tmp.close()

	if stderr:
		print "Something went wrong: %s" % stderr
		sys.exit()

	# Rest
	f = open(rest, "w")
	for line in lines[x:]:
		f.write(line)
	f.close()

	#if os.path.exists(tmp.name):
	#	os.unlink(tmp.name)
	return x, len(lines[x:])	
	
def divide_fa_file(file, sample, rest, fraction, abs_max):
	fa = Fasta(file)
	ids = fa.ids[:]

	x = int(fraction * len(ids))
	if x > abs_max:
		x = abs_max

	sample_seqs = random.sample(ids, x)

	# Rest
	f_sample = open(sample, "w")
	f_rest = open(rest, "w")
	for id,seq in fa.items():
		if id in sample_seqs:
			f_sample.write(">%s\n%s\n" % (id, seq))
		else:
			f_rest.write(">%s\n%s\n" % (id, seq))
	f_sample.close()
	f_rest.close()
	
	return x, len(ids[x:])	

def plot_histogram(values, outfile, xrange=None, breaks=10, title=None, xlabel=None, color=10):
	import matplotlib.pyplot as plt
	import matplotlib.cm as cm
	colors = [cm.Paired(256 / 11 * i) for i in range(11)]
	
	plt.clf()
	try:
		# matplotlib >= 0.99
		plt.hist(values, range=xrange, bins=breaks, color=colors[color], edgecolor="black")
	except:	
		plt.hist(values, range=xrange, bins=breaks)
	plt.xlim(xrange)

	if title:
		plt.title(title)
	
	plt.ylabel("Frequency")
	if xlabel:
		plt.xlabel(xlabel)
	if not outfile.endswith(".svg"):
		outfile += ".svg"
	plt.savefig(outfile, format="svg")

def make_gff_histogram(gfffile, outfile, l, title, breaks=21):
	import matplotlib.pyplot as plt
	data = []
	for line in open(gfffile):
		vals = line.strip().split("\t")
		data.append((int(vals[3]) + int(vals[4])) / 2)

	plt.hist(data, breaks)
	plt.title(title)
	plt.savefig(outfile)

def ks_pvalue(values, l):
	from scipy.stats import kstest
	from numpy import array
	if len(values) == 0:
		return 1.0
	a = array(values, dtype="float") / l
	return kstest(a, "uniform")[1]

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
				except:
					print "Error on line %s while reading %s. Is the file in BED or WIG format?" % (line_count, bedfile)
					sys.exit(1)

				start = (start + end) / 2 - (width / 2)
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

def which(file):
	if not os.environ.has_key("PATH") or not os.environ["PATH"]:
		path = os.defpath
	else:
		path = os.environ["PATH"]

	for p in [os.path.join(x, file) for x in path.split(os.pathsep)]:
		if os.access(p, os.X_OK) and not os.path.isdir(p):
			return p
	
	p = Popen("locate %s" % file, shell=True, stdout=PIPE, stderr=PIPE)
	(stdout, stderr) = p.communicate()
	if not stderr:
		for p in stdout.split("\n"):
			if os.path.basename(p) == file and os.access(p, os.X_OK) and not os.path.isdir(p):
				return p

def get_significant_motifs(motifs, fg_fasta, bg_fasta, e_cutoff=None, p_cutoff=None, save_result=None):
	pass
	
class MotifMatch:
	def __init__(self, seq, name, instance, start, end, strand, score):
		self.sequence = seq
		self.motif_name = name
		self.motif_instance = instance
		self.start = start
		self.end = end
		self.strand = strand
		self.score = score

class MotifResult:
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
		for seq, dict in self.matches.items():
			for motif, mms in dict.items():
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

def parse_gff(gff_file, lowmem=False):
	mr = MotifResult()
	total = 0
	f = open(gff_file)
	BUFSIZE = 10000000
	while 1:
		lines = f.readlines(BUFSIZE)
		if not lines:
			break
		for line in lines:
			(seq, program, feature, start, end, score, strand, bla, extra) = line[:-1].split("\t")
		
			(motif_name, motif_instance) = map(strip, extra.split(";"))
			motif_name = motif_name.split(" ")[1][1:-1]
			motif_instance = motif_instance.split(" ")[1][1:-1]

			mr.sequences[seq] = 1

			if not(mr.motifs.has_key(motif_name)):
				mr.motifs[motif_name] = {}
			if not(mr.motifs[motif_name].has_key(seq)):
				mr.motifs[motif_name][seq] = 0
			mr.motifs[motif_name][seq] += 1
				
		total += len(lines)
	return mr

def calc_motif_enrichment(sample, background, mtc=None, len_sample=None, len_back=None):
	"""Calculate enrichment based on hypergeometric distribution"""
	
	INF = "Inf"

	# Local imports to enable parellel Python calls
	from scipy.stats import hypergeom

	if mtc not in [None, "Bonferroni", "Benjamini-Hochberg", "None"]:
		raise RuntimeError, "Unknown correction: %s" % mtc

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
		motifs = p_value.keys()
		motifs.sort(cmp=lambda x,y: -cmp(p_value[x],p_value[y]))
		l = len(p_value)
		c = l
		for	m in motifs:
			if  p_value[motif] != "NA":
				p_value[m] = p_value[m] * l / c 
			c -= 1

	return (sig, p_value, n_sample, n_back)

def calc_enrichment(sample, background, len_sample, len_back, mtc=None):
	"""Calculate enrichment based on hypergeometric distribution"""
	
	INF = "Inf"

	# Local imports to enable parellel Python calls
	from scipy.stats import hypergeom

	if mtc not in [None, "Bonferroni", "Benjamini-Hochberg", "None"]:
		raise RuntimeError, "Unknown correction: %s" % mtc

	sig = {}
	p_value  = {}
	n_sample = {}
	n_back = {}
	
	for motif in sample.keys():
		p = "NA"
		s = "NA"
		q = sample[motif]
		m = 0
		if(background[motif]):
			m = background[motif]
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
		motifs = p_value.keys()
		motifs.sort(cmp=lambda x,y: -cmp(p_value[x],p_value[y]))
		l = len(p_value)
		c = l
		for	m in motifs:
			if  p_value[motif] != "NA":
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
			except:
				sys.stderr.write("Error in line %s: coordinates in column 2 and 3 need to be integers!\n" % (i))
				sys.exit(1)
	f.close()
	return numpy.median(l)


def locate_tool(tool, verbose=True):
	tool = re.sub(r'[^a-zA-Z]','',tool)
	m = eval("tools." + tool)()
	bin = which(m.cmd)
	if bin:
		print "Found %s in %s" % (m.name, bin)
		return bin
	else:
		print "Couldn't find %s" % m.name

