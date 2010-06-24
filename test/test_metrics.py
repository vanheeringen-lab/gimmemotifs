import unittest
import tempfile
import os
from gimmemotifs.c_metrics import *

class TestMetrics(unittest.TestCase):
	""" A test class for column comparison metrics (metric.py) """

	def setUp(self):
		self.column_pairs = [
			([10,1,1,1], [1,10,1,1]),
			([4,3,2,1], [1,2,3,4]),
			([20,20,20,20], [2,2,2,2]),
			]
	
	def test_pcc(self):
		""" column metric: pearson correlation coefficient (pcc) """
		results = ['-3.333e-01', '-1.000e+00']
		for (col1, col2), result in zip(self.column_pairs[:-1], results):
			self.assertEquals(result, "%0.3e" % score([col1], [col2], "pcc", "mean"))
	
	def test_ed(self):
		""" column metric: euclidian distance (ed) """
		results = ['-1.273e+01', '-4.472e+00', '-3.600e+01']
		for (col1, col2), result in zip(self.column_pairs, results):
			self.assertEquals(result, "%0.3e" % score([col1], [col2], "ed", "mean"))

	def tearDown(self):
		pass

if __name__ == '__main__':
	unittest.main()
