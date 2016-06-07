import unittest
import tempfile
import os
from gimmemotifs.maelstrom import run_maelstrom
import pandas as pd
from glob import glob

class TestMoap(unittest.TestCase):
    """ A test class to test maelstrom"""

    def setUp(self):
        self.outdir = "test/data/maelstrom"
        self.clusters = "test/data/moap/clusters.txt"
        self.outfile = os.path.join(self.outdir, "final.out.csv")

    def test1_maelstrom(self):
        """ Test Motif Activity by Ensemble Learning (maelstrom) """
        
        run_maelstrom(self.clusters, "mm10", self.outdir)
        df = pd.read_table(self.outfile, index_col=0)
        self.assertEquals((623, 4), df.shape)

        for fname in glob(os.path.join(self.outdir, "changed*")):
            os.unlink(fname)
        os.unlink(self.outfile)

