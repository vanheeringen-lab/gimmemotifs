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

    def test1_rankagg(self):
        """ Test rank aggregation """
        df = pd.read_table(self.fname, index_col=0)
        result = rankagg(df)
        self.assertEqual("AP2", result.sort_values().index[0])

if __name__ == '__main__':
    unittest.main()
        
