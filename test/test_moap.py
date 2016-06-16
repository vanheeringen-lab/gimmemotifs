import unittest
import tempfile
import os
from gimmemotifs.moap import moap

class TestMoap(unittest.TestCase):
    """ A test class to test Moap functionality"""

    def setUp(self):
        self.data_dir = "test/data/moap"
        self.clusters = os.path.join(self.data_dir, "clusters.txt")
        self.motifs_score = os.path.join(self.data_dir, "motifs.score.txt")
        self.motifs_count = os.path.join(self.data_dir, "motifs.count.txt")

    def test1_moap(self):
        """ Test motif activity prediction """
        
        for method in ["mwu", "rf", "ks", "lightning"]:
            df = moap(self.clusters,
                    method=method,
                    scoring="score",
                    motiffile=self.motifs_score,
                    )
            self.assertEquals((623, 4), df.shape)

        for method in ["classic"]:
            df = moap(self.clusters,
                    method=method,
                    scoring="count",
                    motiffile=self.motifs_count,
                    )
            self.assertEquals((623, 4), df.shape)
