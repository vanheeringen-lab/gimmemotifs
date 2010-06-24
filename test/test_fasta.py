from gimmemotifs.fasta import *
import tempfile
import unittest
import os

class TestFasta(unittest.TestCase):
	""" A test class for Fasta """

	def setUp(self):
		self.fasta_file = "test/data/fasta/test.fa"
		self.assert_(os.path.exists(self.fasta_file))
		self.assert_(Fasta(self.fasta_file))
		self.f = Fasta(self.fasta_file)
	
	def test1_index(self):
		""" Fasta as a dictionary """
		self.assertEquals(self.f["seq1"], "AAAA")
		self.assertEquals(self.f["seq2"], "ACGT")
		self.assertEquals(self.f["seq3"], "CCCCGGGG")
	
	def test2_items(self):
		""" Fasta.items() """
		self.assertEquals(len(self.f.items()), 3)
	
	def test3_writefasta(self):
		""" Write fasta-formatted file"""
		temp = tempfile.NamedTemporaryFile()
		tempname = temp.name
		self.f.writefasta(tempname)
		self.assertEquals(open(self.fasta_file).read().strip(),open(tempname).read().strip())

	
	def tearDown(self):
			pass

if __name__ == '__main__':
	unittest.main()
