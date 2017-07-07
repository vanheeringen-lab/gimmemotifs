import unittest
import tempfile
import os
from glob import glob
from shutil import rmtree
import pandas as pd

from gimmemotifs.config import MotifConfig
from gimmemotifs.maelstrom import run_maelstrom
from gimmemotifs.genome_index import get_genome,check_genome,GenomeIndex

class TestMoap(unittest.TestCase):
    """ A test class to test maelstrom"""

    def setUp(self):
        self.outdir = "test/data/maelstrom"
        self.clusters = "test/data/moap/clusters.txt"
        self.score_table = "test/data/moap/motifs.score.txt"
        self.count_table = "test/data/moap/motifs.count.txt"
        self.outfile = os.path.join(self.outdir, "final.out.csv")
    
    def test1_maelstrom(self):
        """ Test Motif Activity by Ensemble Learning (maelstrom) """
        
        run_maelstrom(self.clusters, "mm10", self.outdir,
                score_table=self.score_table, count_table=self.count_table,
                plot=False)
        df = pd.read_table(self.outfile, index_col=0, comment="#")
        print(df.shape)
        self.assertEquals((623, 4), df.shape)

        for fname in glob(os.path.join(self.outdir, "activity*")):
            os.unlink(fname)
        os.unlink(self.outfile)
        
        #run_maelstrom(self.clusters, "mm10", self.outdir)
        #df = pd.read_table(self.outfile, index_col=0, comment="#")
        #self.assertEquals((623, 4), df.shape)

        #for fname in glob(os.path.join(self.outdir, "activity*")):
        #    os.unlink(fname)
        #for fname in glob(os.path.join(self.outdir, "motifs*")):
        #    os.unlink(fname)
        #os.unlink(self.outfile)

if __name__ == '__main__':
    unittest.main()


