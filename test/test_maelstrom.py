import unittest
import tempfile
import os
from gimmemotifs.maelstrom import run_maelstrom

class TestMoap(unittest.TestCase):
    """ A test class to test maelstrom"""

    def setUp(self):
        self.outdir = "test/data/maelstrom"
        self.clusters = "test/data/moap/clusters.txt"

    def test1_maelstrom(self):
        """ Test Motif Activity by Ensemble Learning (maelstrom) """
        
        run_maelstrom(self.clusters, "mm10", self.outdir)
