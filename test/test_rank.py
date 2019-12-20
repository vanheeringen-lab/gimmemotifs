import unittest
import tempfile
import os
import pandas as pd
from gimmemotifs.rank import rankagg


class TestRank(unittest.TestCase):
    """ A test class to test rank aggregation """

    def setUp(self):
        self.data_dir = "test/data/rank"
        self.fname = os.path.join(self.data_dir, "ranked.txt")
        self.rank_in = os.path.join(self.data_dir, "rank_input.txt")
        self.rank_out = os.path.join(self.data_dir, "rank_output.txt")

    def test1_rankagg(self):
        """ Test rank aggregation """
        df = pd.read_csv(self.fname, index_col=0, sep="\t")
        result = rankagg(df)
        self.assertEqual("AP2", result.sort_values("score").index[0])

    def test2_rankagg(self):
        """ Test Python implementation of rank aggregation """
        df = pd.read_csv(self.rank_in, index_col=0, sep="\t")
        result = rankagg(df)["score"].values
        ref = pd.read_csv(self.rank_out, index_col=0, sep="\t")["score"].values
        for v1, v2 in zip(ref, result):
            self.assertAlmostEqual(v1, v2)


if __name__ == "__main__":
    unittest.main()
