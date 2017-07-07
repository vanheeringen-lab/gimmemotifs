import unittest
import tempfile
import os
from gimmemotifs.moap import moap

class TestMoap(unittest.TestCase):
    """ A test class to test Moap functionality"""

    def setUp(self):
        self.data_dir = "test/data/moap"
        self.clusters = os.path.join(self.data_dir, "clusters.txt")
        self.motifs_count = os.path.join(self.data_dir, "motifs.count.txt")
        self.motifs_score = os.path.join(self.data_dir, "motifs.score.txt")
        
        self.clusters2 = os.path.join(self.data_dir, "2clusters.txt")
        self.motifs_count2 = os.path.join(self.data_dir, "2clusters.motifs.count.txt")
        self.motifs_score2 = os.path.join(self.data_dir, "2clusters.motifs.score.txt")

    def test1_moap(self):
        """ Test motif activity prediction """
        
        for method in ["mwu", "rf", "lightningclassification"]:
            df = moap(self.clusters,
                    method=method,
                    scoring="score",
                    motiffile=self.motifs_score,
                    )
            self.assertEquals((623, 4), df.shape)

        for method in ["hypergeom"]:
            df = moap(self.clusters,
                    method=method,
                    scoring="count",
                    motiffile=self.motifs_count,
                    )
            self.assertEquals((623, 4), df.shape)

    def test2_moap(self):
        """ Test motif activity prediction for two clusters """
        
        for method in ["mwu", "rf", "lightningclassification"]:
            df = moap(self.clusters2,
                    method=method,
                    scoring="score",
                    motiffile=self.motifs_score2,
                    )
            self.assertEquals((623, 2), df.shape)

        for method in ["hypergeom"]:
            df = moap(self.clusters2,
                    method=method,
                    scoring="count",
                    motiffile=self.motifs_count2,
                    )
            self.assertEquals((623, 2), df.shape)

if __name__ == '__main__':
    unittest.main()

