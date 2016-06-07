import unittest
import tempfile
import os
from glob import glob
from shutil import rmtree
import pandas as pd

from gimmemotifs.maelstrom import run_maelstrom
from gimmemotifs.genome_index import get_genome,check_genome

class TestMoap(unittest.TestCase):
    """ A test class to test maelstrom"""

    def setUp(self):
        self.outdir = "test/data/maelstrom"
        self.clusters = "test/data/moap/clusters.txt"
        self.outfile = os.path.join(self.outdir, "final.out.csv")
        self.fadir = "test/data/genomes"

        try:
            check_genome("mm10")
        except Exception:
            get_genome("mm10", self.fadir)
    
    def test1_maelstrom(self):
        """ Test Motif Activity by Ensemble Learning (maelstrom) """
        
        run_maelstrom(self.clusters, "mm10", self.outdir)
        df = pd.read_table(self.outfile, index_col=0)
        self.assertEquals((623, 4), df.shape)

        for fname in glob(os.path.join(self.outdir, "changed*")):
            os.unlink(fname)
        os.unlink(self.outfile)
