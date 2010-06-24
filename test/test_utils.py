import unittest
import tempfile
import os
from gimmemotifs.utils import *

class TestUtils(unittest.TestCase):
	""" A test class to test utils functions """

	def setUp(self):
		pass

	def test1_phyper(self):
		""" Hypergeometric p-value """
		p = phyper(59, 500,500, 100)
		self.assertAlmostEquals(0.02238075, p)
		
		p = phyper(59, 5000, 5000, 100)
		self.assertAlmostEquals(0.02782685, p)

		p = phyper(59, 50000, 50000, 100)
		self.assertAlmostEquals(0.02838217, p)

	def tearDown(self):
		pass

if __name__ == '__main__':
	unittest.main()
