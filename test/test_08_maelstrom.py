import os
import unittest
from glob import glob

import pandas as pd

from gimmemotifs.maelstrom import run_maelstrom


class TestMaelstrom(unittest.TestCase):
    """A test class to test maelstrom"""

    def setUp(self):
        self.outdir = "test/data/maelstrom"
        self.clusters = "test/data/moap/4clusters.txt"
        self.score_table = "test/data/moap/4clusters.motifs.score.txt"
        self.count_table = "test/data/moap/4clusters.motifs.count.txt"
        self.outfile = os.path.join(self.outdir, "final.out.txt")

    def test1_maelstrom(self):
        """Test Motif Activity by Ensemble Learning (maelstrom)"""

        run_maelstrom(
            self.clusters,
            "mm10",
            self.outdir,
            filter_redundant=False,
            score_table=self.score_table,
            count_table=self.count_table,
            plot=False,
        )
        df = pd.read_table(self.outfile, index_col=0, comment="#")
        print(df.shape)

        self.assertEqual((623, 8), df.shape)

        # Filter redundant motifs
        run_maelstrom(
            self.clusters,
            "mm10",
            self.outdir,
            filter_redundant=True,
            score_table=self.score_table,
            count_table=self.count_table,
            plot=False,
        )
        df = pd.read_table(self.outfile, index_col=0, comment="#")
        print(df.shape)
        self.assertEqual((156, 8), df.shape)

        for fname in glob(os.path.join(self.outdir, "activity*")):
            os.unlink(fname)
        for fname in glob(os.path.join(self.outdir, "gimme.verte*")):
            os.unlink(fname)
        os.unlink(self.outfile)

        # run_maelstrom(self.clusters, "mm10", self.outdir)
        # df = pd.read_table(self.outfile, index_col=0, comment="#")
        # self.assertEquals((623, 4), df.shape)

        # for fname in glob(os.path.join(self.outdir, "activity*")):
        #    os.unlink(fname)
        # for fname in glob(os.path.join(self.outdir, "motifs*")):
        #    os.unlink(fname)
        # os.unlink(self.outfile)


if __name__ == "__main__":
    unittest.main()
