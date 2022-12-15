import os
import unittest

import numpy as np

from gimmemotifs.maelstrom.moap import moap


class TestMoap(unittest.TestCase):
    """A test class to test Moap functionality"""

    def setUp(self):
        self.random_state = np.random.RandomState(1)
        self.data_dir = "test/data/moap"
        self.clusters4 = os.path.join(self.data_dir, "4clusters.txt")
        self.motifs_count4 = os.path.join(self.data_dir, "4clusters.motifs.count.txt")
        self.motifs_score4 = os.path.join(self.data_dir, "4clusters.motifs.score.txt")

        self.clusters2 = os.path.join(self.data_dir, "2clusters.txt")
        self.motifs_count2 = os.path.join(self.data_dir, "2clusters.motifs.count.txt")
        self.motifs_score2 = os.path.join(self.data_dir, "2clusters.motifs.score.txt")

    def test1_moap(self):
        """Test motif activity prediction"""

        # for method in ["bayesianridge", "mwu", "rf", "multitasklasso", "svr"]:
        for method in ["mwu", "rf"]:
            df = moap(
                self.clusters4,
                method=method,
                scoring="score",
                motiffile=self.motifs_score4,
            )
            self.assertEquals((623, 4), df.shape)

        # for method in ["hypergeom", "bayesianridge", "rf", "multitasklasso", "svr"]:
        for method in ["hypergeom"]:
            df = moap(
                self.clusters4,
                method=method,
                scoring="count",
                motiffile=self.motifs_count4,
            )
            self.assertEquals((623, 4), df.shape)

    def test2_moap(self):
        """Test motif activity prediction for two clusters"""

        # for method in ["bayesianridge", "mwu", "rf", "multitasklasso", "svr"]:
        for method in ["mwu", "rf"]:
            df = moap(
                self.clusters2,
                method=method,
                scoring="score",
                motiffile=self.motifs_score2,
            )
            self.assertEquals((623, 2), df.shape)

        # for method in ["hypergeom", "bayesianridge", "rf", "multitasklasso", "svr"]:
        for method in ["hypergeom"]:
            df = moap(
                self.clusters2,
                method=method,
                scoring="count",
                motiffile=self.motifs_count2,
            )
            self.assertEquals((623, 2), df.shape)


if __name__ == "__main__":
    unittest.main()
