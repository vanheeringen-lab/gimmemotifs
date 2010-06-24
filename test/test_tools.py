import unittest
import tempfile
import os
from gimmemotifs.tools import *

class TestMotifProgram(unittest.TestCase):
	""" A test class for motifprograms """

	def setUp(self):
		self.data_dir = "test/data/motifprogram"
		self.fa = os.path.join(self.data_dir, "test.fa")

	def test_meme(self):
		""" test meme """
		m = Meme()
		if m.is_installed():
			(motifs, stderr, stdout) =  m.run(self.fa, ".")
			#print "meme:", motifs
			self.assert_(len(motifs) > 0)
		else:
			sys.stderr.write("Skipping meme test\n")

	def test_MDmodule(self):
		""" test MDmodule """
		m = MDmodule()
		if m.is_installed():
			(motifs, stderr, stdout) =  m.run(self.fa, ".")
			#print motifs
			self.assert_(len(motifs) > 0)
		else:
			sys.stderr.write("Skipping MDmodule test\n")

	def test_Weeder(self):
		""" test Weeder """
		m = Weeder()
		if m.is_installed():
			(motifs, stderr, stdout) =  m.run(self.fa, ".")
			self.assert_(len(motifs) > 0)
		else:
			sys.stderr.write("Skipping Weeder test\n")

	def test_MotifSampler(self):
		""" test MotifSampler """
		m = MotifSampler()
		if m.is_installed():
			(motifs, stderr, stdout) =  m.run(self.fa, ".", {"background":os.path.join(self.data_dir, "test.bg")})
			#print motifs
			self.assert_(len(motifs) > 0)
		else:
			sys.stderr.write("Skipping MotifSampler test\n")

	#def test_gadem(self):
	#	""" test gadem """
	#	self.skipTest()
	#	m = Gadem()
	#	if m.is_installed():
	#		(motifs, stderr, stdout) =  m.run(self.fa, ".")
	#		#print motifs
	#		self.assert_(len(motifs) > 0)
	#	else:
	#		sys.stderr.write("Skipping gadem test\n")


	def tearDown(self):
		pass

if __name__ == '__main__':
	unittest.main()
