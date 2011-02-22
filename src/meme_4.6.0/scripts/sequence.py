#!@WHICHPYTHON@

import string, sys

#------------------ Alphabet -------------------

class Alphabet(object):
	"""Biological alphabet class. 
	This defines the set of symbols from which various objects can be built, e.g. sequences and motifs.
	The symbol set is immutable and accessed as a tuple.
	symstr: symbols in alphabet as either a tuple or string
	"""
	def __init__(self, symstr):
		"""Construct an alphabet from a string or tuple of characters.
		Lower case characters will be converted to upper case.
		Example:
		>>> alpha = sequence.Alphabet('ACGTttga')
		>>> alpha.getSymbols()
		will construct the DNA alphabet and output:
		('A', 'C', 'G', 'T')
		"""
		symlst = []
		for s in [str(sym).upper()[0] for sym in symstr]:
			if not s in symlst:
				symlst.append(s)
		self.symbols = tuple(symlst)
	
	def getSymbols(self):
		"""Retrieve a tuple with all symbols, immutable membership and order"""
		return self.symbols
	
	def isValidSymbol(self, sym):
		"""Check if the symbol is a member of alphabet"""
		return any([s==sym for s in self.symbols])
		
	def getIndex(self, sym):
		"""Retrieve the index of the symbol (immutable)"""
		for idx in range (len(self.symbols)):
			if self.symbols[idx] == sym:
				return idx
		raise RuntimeError("Symbol " + sym + " does not exist in alphabet")
	
	def isValidString(self, symstr):
		"""Check if the string contains only symbols that belong to the alphabet"""
		found = True
		for sym in symstr:
			if self.isValidSymbol(sym) == False:
				return False
		return True

	def getLen(self):
		"""Retrieve the number of symbols in (the length of) the alphabet"""
		return len(self.symbols)
	
# pre-defined alphabets that can be specified by their name
predefAlphabets = [
				   ("DNA"                , Alphabet('ACGT')),
				   ("RNA"                , Alphabet('ACGU')),
				   ("Extended DNA"       , Alphabet('ACGTYRN')),
				   ("Protein"            , Alphabet('ACDEFGHIKLMNPQRSTVWY')),
				   ("Extended Protein"   , Alphabet('ACDEFGHIKLMNPQRSTVWYX')),
				   ("TM Labels"          , Alphabet('MIO'))
				   ]

def getAlphabet(name):
	"""Retrieve a pre-defined alphabet by name.
	Currently, "Protein", "DNA", "RNA", "Extended DNA", "Extended Protein" and "TM Labels" are available.
	Example:
	>>> alpha = sequence.getAlphabet('Protein')
	>>> alpha.getSymbols()
	will retrieve the 20 amino acid alphabet and output the tuple:
	('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
	"""
	for (xname, xalpha) in predefAlphabets:
		if xname == name:
			return xalpha
	return None

#------------------ Sequence -------------------

class Sequence(object):
	"""Biological sequence class. Sequence data is immutable.
	
	data: the sequence data as a tuple or string
	alpha: the alphabet from which symbols are taken
	name: the sequence name, if any
	info: can contain additional sequence information apart from the name
	"""
	def __init__(self, sequence, alpha = None, name = "", seqinfo = ""):
		"""Create a sequence with sequence data. 
		Specifying the alphabet is optional, so is the name and info.
		Example:
		>>> myseq = sequence.Sequence('MVSAKKVPAIAMSFGVSF')
		will create a sequence with name "", and assign one of the predefined alphabets on basis of what symbols were used.
		>>> myseq.getAlphabet().getSymbols()
		will most likely output the standard protein alphabet:
		('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y')
		"""
		if type(sequence) is str:
			self.data = tuple(sequence.upper())
		elif type(sequence) is tuple:
			self.data = sequence
		elif type(sequence) is list:
			self.data = tuple([s.upper() for s in sequence])
		else:
			raise RuntimeError("Sequence data is not specified correctly: must be string or tuple")
		# Resolve choice of alphabet
		validAlphabet = False
		if alpha == None:                                   # Alphabet is not set, attempt to set it automatically...
			for (xname, xalpha) in predefAlphabets:         # Iterate through each predefined alphabet, in order
				if xalpha.isValidString( self.data ):        # This alphabet works, go with it
					self.alpha = alpha = xalpha
					validAlphabet = True
					break
		self.name = name
		self.info = seqinfo
		if validAlphabet == False:            # we were either unsuccessful above or the alphabet was specified so test it
			if type(alpha) is str:            # check if name is a predefined alphabet
				for (xname, xalpha) in predefAlphabets:   # Iterate through each predefined alphabet, check for name
					if (xname == alpha):
						alpha = xalpha
						break
			if type(alpha) is Alphabet:       # the alphabet is specified
				if alpha.isValidString(self.data) == False:
					raise RuntimeError("Invalid alphabet specified: "+"".join(alpha.getSymbols())+" is not compatible with sequence '"+"".join(self.data)+"'")
				else:
					self.alpha = alpha
			else:
				raise RuntimeError("Could not identify alphabet from sequence")
			
	#basic getters and setters for the class
	def getName(self):
		"""Get the name of the sequence"""
		return self.name
	def getInfo(self):
		"""Get additional info of the sequence (e.g. from the defline in a FASTA file)"""
		return self.info
	def getAlphabet(self):
		"""Retrieve the alphabet that is assigned to this sequence"""
		return self.alpha
	def setName(self, name):
		"""Change the name of the sequence"""
		self.name = name
	def setAlphabet(self, alpha):
		"""Set the alphabet, throws an exception if it is not compatible with the sequence data"""
		if type(alpha) is Alphabet:       
			if alpha.isValid( sequence ) == False:
				raise RuntimeError( "Invalid alphabet specified" )
	#sequence functions	
	def getSequence(self):
		"""Retrieve the sequence data (a tuple of symbols)"""
		return self.data
	def getString(self):
		"""Retrieve the sequence data as a string (copy of actual data)"""
		return "".join(self.data)
	def getLen(self):
		"""Get the length of the sequence (number of symbols)"""
		return len(self.data)
	def getSite(self, position, length = 1):
		"""Retrieve a site in the sequence of desired length. 
		Note that positions go from 0 to length-1, and that if the requested site 
		extends beyond those the method throws an exception.
		"""
		if position >= 0 and position <= self.getLen() - length:
			if length == 1:
				return self.data[position]
			else:
				return self.data[position:position+length]
		else:
			raise RuntimeError( "Attempt to access invalid position in sequence "+self.getName() )

	def nice(self):
		""" A short description of the sequence """
		print self.getName(), ":", self.getLen()
			
def readStrings(filename):
	""" Read one or more lines of text from a file--for example an alignment. 
	Return as a list of strings. 
	filename: name of file
	"""
	txtlist = []
	f = open(filename)
	for line in f.readlines():
		txtlist.extend(line.split())
	return txtlist
	
def readFASTA(filename, alpha = None):
	""" Read one or more sequences from a file in FASTA format.
	filename: name of file to load sequences from
	alpha: alphabet that is used (if left unspecified, an attempt is made to identify the alphabet for each individual sequence)
	"""
	seqlist = []
	seqname = None
	seqinfo = None
	seqdata = []
	fh = open(filename)
	thisline = fh.readline()
	while (thisline):
		if (thisline[0] == '>'): # new sequence
			if (seqname):        # take care of the data that is already in the buffer before processing the new sequence
				try:
					seqnew = Sequence(seqdata, alpha, seqname, seqinfo)
					seqlist.append(seqnew)
				except RuntimeError, e:
					print >> sys.stderr, "Warning: "+seqname+" is invalid (ignored): ", e
			seqinfo = thisline[1:-1]         # everything on the defline is "info"
			seqname = seqinfo.split()[0]     # up to first space
			seqdata = []
		else:  # pull out the sequence data
			cleanline = thisline.split()
			for line in cleanline:
				seqdata.extend(tuple(line.strip('*'))) # sometimes a line ends with an asterisk in FASTA files
		thisline = fh.readline()
		
	if (seqname):
		try:
			seqnew = Sequence(seqdata, alpha, seqname, seqinfo)
			seqlist.append(seqnew)
		except RuntimeError, e:
			print >> sys.stderr, "Warning: " + seqname + " is invalid (ignored): ", e
	else:
		raise RuntimeError("No sequences on FASTA format found in this file")
	fh.close()
	return seqlist

def _writeOneFASTA(sequence, filehandle):
	"""Write one sequence in FASTA format to an already open file"""
	filehandle.write(">" + sequence.getName()+"\n")
	data = sequence.getSequence()
	lines = ( sequence.getLen() - 1) / 60 + 1
	for i in range(lines):
		#note: python lets us get the last line (var length) free
		#lineofseq = data[i*60 : (i+1)*60] + "\n"
		lineofseq = "".join(data[i*60 : (i+1)*60]) + "\n"
		filehandle.write(lineofseq)
	
def writeFASTA(sequence, filename):
	"""Write a list (or a single) of sequences to a file in the FASTA format"""
	fh = open(filename, "w")
	if isinstance(sequence, Sequence):
		_writeOneFASTA(sequence, fh)
	else:
		for seq in sequence:
			if isinstance(seq, Sequence):
				_writeOneFASTA(seq, fh)
			else:
				print >> sys.stderr, "Warning: could not write " + seq.getName() + " (ignored)."
	fh.flush()	
	fh.close()

#------------------ Distrib -------------------

class Distrib(object):
	"""Class for storing a multinomial probability distribution over the symbols in an alphabet""" 
	def __init__(self, alpha, pseudo_count = 0.0):
		self.alpha = alpha
		self.tot = pseudo_count * self.alpha.getLen()
		self.cnt = [pseudo_count for _ in range( self.alpha.getLen() )]
		
	def count(self, syms = None ):
		"""Count an observation of a symbol"""
		if syms == None:
			syms = self.alpha.getSymbols()
		for sym in syms:
			idx = self.alpha.getIndex( sym )
			self.cnt[idx] += 1.0
			self.tot += 1

	def reset(self):
		"""Reset the distribution, that is, restart counting."""
		self.tot = 0
		self.cnt = [0.0 for _ in range( self.alpha.getLen() )]
		 
	def getFreq(self, sym = None):
		"""Determine the probability distribution from the current counts.
		The order in which probabilities are given follow the order of the symbols in the alphabet."""
		if self.tot > 0:
			if sym == None:
				freq = tuple([ y / self.tot for y in self.cnt ])
				return freq
			else:
				idx = self.alpha.getIndex( sym )
				return self.cnt[idx] / self.tot
		return None

	def pretty(self):
		"""Retrieve the probabilites for all symbols and return as a pretty table (a list of text strings)"""
		table = ["".join(["%4s " % s for s in self.alpha.getSymbols()])]
		table.append("".join(["%3.2f " % y for y in Distrib.getFreq(self)]))
		return table    

	def getSymbols(self):
		"""Get the symbols in the alphabet in the same order as probabilities are given."""
		return self.alpha.getSymbols()

	def getAlphabet(self):
		"""Get the alphabet over which the distribution is defined."""
		return self.alpha
	
#------------------ Motif (and subclasses) -------------------

class Motif(object):
	""" Sequence motif class--defining a pattern that can be searched in sequences.
	This class is not intended for direct use. Instead use and develop sub-classes (see below).
	"""
	def __init__(self, alpha):
		self.len = 0
		self.alpha = alpha
		
	def getLen(self):
		"""Get the length of the motif"""
		return self.len

	def getAlphabet(self):
		"""Get the alphabet that is used in the motif"""
		return self.alpha
	
	def isAlphabet(self, seqstr):
		"""Check if the sequence can be processed by this motif"""
		mystr = seqstr
		if type(seqstr) is Sequence:
			mystr = seqstr.getString()
		return self.getAlphabet().isValidString(mystr)
		
import re

class RegExp(Motif):
	"""A motif class that defines the pattern in terms of a regular expression"""
	def __init__(self, alpha, re_string):
		Motif.__init__(self, alpha)
		self.pattern = re.compile(re_string)
		
	def match(self, seq):
		"""Find matches to the motif in a specified sequence.
		The method is a generator, hence subsequent hits can be retrieved using next().
		The returned result is a tuple (position, match-sequence, score), where score is
		always 1.0 since a regular expression is either true or false (not returned).
		"""
		myseq = seq
		if not type(seq) is Sequence:
			myseq = Sequence(seq, self.alpha)
		mystr = myseq.getString()
		if not Motif.isAlphabet(self, mystr):
			raise RuntimeError("Motif alphabet is not valid for sequence " + myseq.getName())
		for m in re.finditer(self.pattern, mystr):
			yield (m.start(), m.group(), 1.0)
		
import math

class PWM(Motif):
	"""This motif subclass defines a pattern in terms of a position weight matrix.
	An alphabet must be provided. A pseudo-count to be added to each count is
	optional.  A uniform background distribution is used by default.
	"""
	def __init__(self, alpha):
		Motif.__init__(self, alpha)			# set alphabet of this multinomial distribution
		self.background = Distrib(alpha)		# the default background ...
		self.background.count(alpha.getSymbols())   	# ... is uniform
		
	def setFromAlignment(self, aligned, pseudo_count = 0.0):
		"""Set the probabilities in the PWM from an alignment.
		The alignment is a list of equal-length strings (see readStrings), OR
		a list of Sequence.
		"""
		self.cols = -1
		seqs = []
		# Below we create a list of Sequence from the alignment, 
		# while doing some error checking, and figure out the number of columns
		for s in aligned:
			# probably a text string, so we make a nameless sequence from it
			if not type(s) is Sequence: 
				s=Sequence(s, Motif.getAlphabet(self))
			else: 
			# it was a sequence, so we check that the alphabet in 
			# this motif will be able to process it
				if not Motif.isAlphabet(self, s):
					raise RuntimeError("Motif alphabet is not valid for sequence " + s.getName())
			if self.cols == -1:
				self.cols = s.getLen()
			elif self.cols != s.getLen():
				raise RuntimeError("Sequences in alignment are not of equal length")
			seqs.append(s)
		# The line below initializes the list of Distrib (one for each column of the alignment)
		self.counts = [Distrib(Motif.getAlphabet(self), pseudo_count) for _ in range(self.cols)]
		# Next, we do the counting, column by column
		for c in range( self.cols ):     # iterate through columns
			for s in seqs:               # iterate through rows
				# determine the index of the symbol we find at this position (row, column c)
				self.counts[c].count(s.getSite(c))   
			
	def setBackground(self, distrib):
		"""Set the background distribution"""
		if not distrib.getAlphabet() == Motif.getAlphabet(self):
			raise RuntimeError("Incompatible alphabets")
		self.background = distrib
		
	def getFreq(self):
		"""Get the probabilities for all positions in the PWM (a list of Distribs)"""
		return [y.getFreq() for y in self.counts]
    
	def pretty(self):
		"""Retrieve the probabilites for all positions in the PWM as a pretty table (a list of text strings)"""
		#table = ["".join(["%8s " % s for s in self.alpha.getSymbols()])]
		table = []
		for row in PWM.getFreq(self):
			table.append("".join(["%8.6f " % y for y in row]))
		return table    

	def logoddsPretty(self, bkg):
		"""Retrieve the (base-2) log-odds for all positions in the PWM as a pretty table (a list of text strings)"""
		table = []
		for row in PWM.getFreq(self):
			#table.append("".join(["%8.6f " % (math.log((row[i]+1e-6)/bkg[i])/math.log(2)) for i in range(len(row))]))
			table.append("".join(["%8.6f " % (math.log((row[i])/bkg[i])/math.log(2)) for i in range(len(row))]))
			#table.append("".join(["%8.6f " % row[i] for i in range(len(row))]))
		return table    


	def consensus_sequence(self):
		"""
		Get the consensus sequence corresponding to a PWM.
		Consensus sequence is the letter in each column
		with the highest probability.
		"""
		consensus = ""
		alphabet = Motif.getAlphabet(self).getSymbols()
		for pos in range(self.cols):
			best_letter = alphabet[0]
			best_p = self.counts[pos].getFreq(best_letter)
			for letter in alphabet[1:]:
				p = self.counts[pos].getFreq(letter)
				if p > best_p:
					best_p = p
					best_letter = letter
			consensus += best_letter
		return consensus
    

	def consensus(self):
		"""
		Get the consensus corresponding to a PWM.
		Consensus at each column of motif is a list of
		characters with non-zero probabilities.
		"""
		consensus = []
		for pos in range(self.cols):
			matches = [] 
			for letter in Motif.getAlphabet(self).getSymbols():
				p = self.counts[pos].getFreq(letter)
				if p > 0:
					matches += letter
			consensus.append(matches)
		return consensus

    
	def getScore(self, seq, start):
		"""Score this particular list of symbols using the PFM (background needs to be set separately)"""
		sum = 0.0
		seqdata = seq.getSequence()[start : start+self.cols]
		for pos in range(len(seqdata)):
			q = self.counts[pos].getFreq(seqdata[pos])
			if q == 0:
				q = 0.0001 # to avoid log(0) == -Infinity
			logodds = math.log(q / self.background.getFreq(seqdata[pos]))
			sum += logodds
		return sum

	def match(self, seq, _LOG0 = -10):
		"""Find matches to the motif in a specified sequence.
		The method is a generator, hence subsequent hits can be retrieved using next().
		The returned result is a tuple (position, match-sequence, score).
		The optional parameter _LOG0 specifies a lower bound on reported logodds scores.
		"""
		myseq = seq
		if not type(seq) is Sequence:
			myseq = Sequence(seq, self.alpha)
		if not Motif.isAlphabet(self, myseq):
			raise RuntimeError("Motif alphabet is not valid for sequence " + myseq.getName())
		for pos in range(myseq.getLen() - self.cols):
			score = PWM.getScore(self, myseq, pos)
			if score > _LOG0: 
				yield (pos, "".join(myseq.getSite(pos, self.cols)), score)

#------------------ Main method -------------------
# Executed if you run this file from the operating system prompt, e.g.
# > python sequence.py

if __name__=='__main__':
	alpha = getAlphabet('Extended DNA')
	#seqs = readFASTA('pos.fasta')
	seqs = []
	aln = readStrings('tmp0')
	#regexp = RegExp(alpha, '[AG]G.[DE]TT[AS].')
	pwm = PWM(alpha)
	pwm.setFromAlignment(aln)
	for row in pwm.pretty():
		print row
	for s in seqs:
		print s.getName(), s.getLen(), s.getAlphabet().getSymbols()
		for m in regexp.match( s ):
			print "pos: %d pat: %s %4.2f" % (m[0], m[1], m[2])
		for m in pwm.match( s ):
			print "pos: %d pat: %s %4.2f" % (m[0], m[1], m[2])
