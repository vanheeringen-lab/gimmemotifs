import unittest
import tempfile
import os
from motiftools.motif import *
from motiftools.fasta import Fasta
from motiftools.utils import gff_enrichment 
from time import sleep

class TestMotifPwm(unittest.TestCase):
	""" A test class to test Motif pwmscan functionality and related things """

	def setUp(self):
		self.data_dir = "test/data/pwmscan"
		
		self.motif = pwmfile_to_motifs(os.path.join(self.data_dir, "TATA.pwm"))[0]
		self.prom = Fasta(os.path.join(self.data_dir, "promoters.fa"))
		self.prom_gff = os.path.join(self.data_dir, "promoters_result.gff")
		self.random = Fasta(os.path.join(self.data_dir, "random_sequences.fa"))
		self.random_gff = os.path.join(self.data_dir, "random_result.gff")
		self.enrichment = os.path.join(self.data_dir, "enrichment.txt")
		self.tmp = NamedTemporaryFile().name
	
	def test1_pwm_scan(self):
		""" Scan a FASTA file with PWM of motif """
		result = self.motif.pwm_scan(self.prom, nreport=1)

		# Every sequence should have a TATA match
		self.assertEquals(len(result.keys()), len(self.prom.items()))

	def test2_pwm_scan_to_gff(self):
		""" Scan a FASTA file with PWM of motif, and produce GFF """
		
		self.motif.pwm_scan_to_gff(self.prom, self.tmp)
		self.assertEquals(open(self.prom_gff).read(), open(self.tmp).read())

	def test3_gff_enrichment(self):
		""" Test gff_enrichment """
		self.motif.pwm_scan_to_gff(self.random, self.random_gff)
		gff_enrichment(self.prom_gff, self.random_gff, 316, 3160, self.tmp)
		self.assertEquals(open(self.enrichment).read(), open(self.tmp).read())

	def tearDown(self):
		pass

if __name__ == '__main__':
	unittest.main()
