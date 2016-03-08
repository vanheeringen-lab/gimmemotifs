# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

""" Module to work with FASTA files """
import sys
import random
import re

class Fasta:

	def __init__(self, file=None, split_whitespace=False):
		""" Instantiate fasta object. Optional Fasta-formatted file as argument"""
		self.ids = []
		self.seqs = []
		p = re.compile(r'[^abcdefghiklmnpqrstuvwyzxABCDEFGHIKLMNPQRSTUVWXYZ]')
		if file:
			f = open(file, "r")
			c = f.read()
			f.close()
			if not(c.startswith(">")):
				raise IOError, "Not a valid FASTA file"
			
			for seq in c.split(">"):
				if len(seq) > 1:
					lines = seq.split("\n")
					seq_name = lines[0]
					if split_whitespace:
						seq_name = seq_name.split(" ")
					self.ids.append(seq_name)
					sequence = "".join(lines[1:])
					if p.match(sequence):
						raise IOError, "Not a valid FASTA file"
					self.seqs.append(sequence)
		
	def hardmask(self):
		""" Mask all lowercase nucleotides with N's """
		p = re.compile("a|c|g|t|n")
		for id in self.fasta_dict.keys():
			self.fasta_dict[id] = p.sub("N", self.fasta_dict[id])
		return self

	def get_random(self, n, l=None):
		""" Return n random sequences from this Fasta object """
		random_f = Fasta()
		if l:
			ids = self.ids[:]
			random.shuffle(ids)
			i = 0
			while (i < n) and (len(ids) > 0):
				id = ids.pop()
				if (len(self[id]) >= l):
					start = random.randint(0, len(self[id]) - l)
					random_f["random%s" % (i + 1)] = self[id][start:start+l]
					i += 1
			if len(random_f) != n:
				sys.stderr.write("Not enough sequences of required length")
				return
			else:
				return random_f

		else:
			choice = random.sample(self.ids, n)
			for i in range(n):
				random_f[choice[i]] = self[choice[i]]
		return random_f


	def __getitem__(self, idx):
		if isinstance(idx, slice):
			f = Fasta()
			f.ids = self.ids[idx][:]
			f.seqs = self.seqs[idx][:]
			return f
		elif idx in self.ids:
			return self.seqs[self.ids.index(idx)]
		else:
			return None

	def __repr__(self):
		return "%s sequences" % len(self.ids)

	def __len__(self):
		return len(self.ids)

	def __setitem__(self, key, value):
		if key in self.ids:
			self.seqs[self.ids.index(key)] = value
		else:
			self.ids.append(key)
			self.seqs.append(value)

	def __delitem__(self, key):
		i = self.ids.index(key)
		self.ids.pop(i)
		self.seqs.pop(i)
		
	def _format_seq(self, seq):
		return seq

	def add(self, id, seq):
		self.ids.append(id)
		self.seqs.append(seq)
	
	def has_key(self, key):
		if key in self.ids:
			return True
		else:
			return False

	def __str__(self):
		return "\n".join([">%s\n%s" % (id, self._format_seq(seq)) for id, seq in self.items()])

	def writefasta(self, file):
		""" Write sequences to FASTA formatted file"""
		f = open(file, "w")
		f.write(self.__str__())
		f.close()

	def items(self):
		return zip(self.ids, self.seqs)
