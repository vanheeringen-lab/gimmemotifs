import unittest
import os
import numpy as np
from gimmemotifs.moap import moap


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

        expected_clusters_score = {
            "bayesianridge": 4,  # 1,
            # "xgboost",
            "mwu": 4,
            "rf": 4,
            "multitasklasso": 4,  # 1,
            "svr": 4,  # 1,
        }
        for method, nclust in expected_clusters_score.items():
            df = moap(
                self.clusters4,
                method=method,
                scoring="score",
                motiffile=self.motifs_score4,
                random_state=self.random_state,
            )
            self.assertEquals(
                df.shape[0],
                623,
                msg=f"{method} score, expected {623} sequences, but found {df.shape[0]}",
            )
            self.assertEquals(
                df.shape[1],
                4,
                msg=f"{method} score, expected {nclust} clusters but found {df.shape[1]}",
            )

        expected_clusters_count = {
            "hypergeom": 4,
            "bayesianridge": 4,  # 1,
            # "xgboost",
            "rf": 4,
            "multitasklasso": 4,  # 1,
            "svr": 4,  # 1,
        }
        for method, nclust in expected_clusters_count.items():
            df = moap(
                self.clusters4,
                method=method,
                scoring="count",
                motiffile=self.motifs_count4,
                random_state=self.random_state,
            )
            msg = f"{method} count, expected {nclust} clusters but found {df.shape[1]}"
            self.assertEquals(nclust, df.shape[1], msg=msg)

    def test2_moap(self):
        """Test motif activity prediction for two clusters"""

        expected_clusters_score = {
            "bayesianridge": 2,  # 1,
            # "xgboost",
            "mwu": 2,
            "rf": 2,
            "multitasklasso": 2,  # 1,
            "svr": 2,  # 1,
        }
        for method, nclust in expected_clusters_score.items():
            df = moap(
                self.clusters2,
                method=method,
                scoring="score",
                motiffile=self.motifs_score2,
                random_state=self.random_state,
            )
            msg = f"{method} score, expected {nclust} clusters but found {df.shape[1]}"
            self.assertEquals(nclust, df.shape[1], msg=msg)

        expected_clusters_count = {
            "hypergeom": 2,
            "bayesianridge": 2,  # 1,
            # "xgboost",
            "rf": 2,
            "multitasklasso": 2,  # 1,
            "svr": 2,  # 1,
        }
        for method, nclust in expected_clusters_count.items():
            df = moap(
                self.clusters2,
                method=method,
                scoring="count",
                motiffile=self.motifs_count2,
                random_state=self.random_state,
            )
            msg = f"{method} count expected {nclust} clusters but found {df.shape[1]}"
            self.assertEquals(nclust, df.shape[1], msg=msg)


if __name__ == "__main__":
    unittest.main()
