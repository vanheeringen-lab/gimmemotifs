import unittest
import tempfile
import os
from gimmemotifs.c_metrics import *
from gimmemotifs.rocmetrics import pr_auc


class TestMetrics(unittest.TestCase):
    """ A test class for column comparison metrics (metric.py) """

    def setUp(self):
        self.column_pairs = [
            ([10, 1, 1, 1], [1, 10, 1, 1]),
            ([4, 3, 2, 1], [1, 2, 3, 4]),
            ([20, 20, 20, 20], [2, 2, 2, 2]),
        ]

    def test_pcc(self):
        """ column metric: pearson correlation coefficient (pcc) """
        results = ["-3.333e-01", "-1.000e+00"]
        for (col1, col2), result in zip(self.column_pairs[:-1], results):
            self.assertEqual(result, "%0.3e" % score([col1], [col2], "pcc", "mean"))

    def test_ed(self):
        """ column metric: euclidian distance (ed) """
        results = ["-1.273e+01", "-4.472e+00", "-3.600e+01"]
        for (col1, col2), result in zip(self.column_pairs, results):
            self.assertEqual(result, "%0.3e" % score([col1], [col2], "ed", "mean"))

    def test_pr_auc(self):
        """ Test PR AUC"""
        fg_values = [4, 5, 3, 6, 5]
        bg_values = [1, 3, 5, 0, 1, 1, 3, 2]
        expect = 0.7850
        result = pr_auc(fg_values, bg_values)
        self.assertAlmostEqual(result, expect, 4)

    def tearDown(self):
        pass


if __name__ == "__main__":
    unittest.main()
