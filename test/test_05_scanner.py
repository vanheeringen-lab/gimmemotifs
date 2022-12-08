import unittest
import tempfile
import os
import numpy as np
from gimmemotifs.scanner import Scanner, scan_to_best_match
from gimmemotifs.fasta import Fasta
import pytest


class TestScanner(unittest.TestCase):
    """A test class to test scanner functionality"""

    def setUp(self):
        self.s = Scanner(
            ncpus=2,
            random_state=np.random.RandomState(1),
            progress=True,
        )
        self.data_dir = "test/data/scanner"

        self.motifs = os.path.join(self.data_dir, "motif.pwm")
        self.fa = os.path.join(self.data_dir, "test.fa")
        self.bed = os.path.join(self.data_dir, "test.bed")
        self.regions = os.path.join(self.data_dir, "test.txt")

        self.tmpdir = tempfile.mkdtemp()

    def test01_set_motifs(self):
        assert self.s.motifs is None
        assert len(self.s.motif_ids) == 0

        self.s.set_motifs(self.motifs)
        assert self.s.motifs == self.motifs
        assert len(self.s.motif_ids) > 0
        assert self.s.motif_ids[0] == "AP1"

    def test02_set_genome(self):
        assert self.s.genome is None

        self.s.set_genome(None)
        assert self.s.genome is None

        self.s.set_genome(self.fa)
        assert self.s.genome == os.path.abspath(self.fa)

    def test03_set_background(self):
        assert self.s.background is None
        with pytest.raises(ValueError):
            self.s.set_background()  # no genome, no bgfile

        # genome: all attributes are set
        self.s.set_genome(self.fa)
        self.s.set_background(size=2, nseq=10, gc=False)
        assert self.s.gc_bins == ["0.00-1.00"]  # gc=False
        assert self.s.background is not None
        assert self.s.background_hash is not None
        assert len(self.s.background.ids) == 10  # nseq
        assert len(self.s.background.seqs[0]) == 2  # size

        # caching works
        before = self.s.background_hash
        self.s.set_background(size=2, nseq=10, gc=False)
        after = self.s.background_hash
        assert before == after

        # bgfile: all attributes are set
        fname = "test/data/scan/scan_test_regions.fa"
        self.s.gc_bins = None
        self.s.background = None
        self.s.background_hash = None
        self.s.set_background(fasta=fname)
        assert self.s.gc_bins == ["0.00-1.00"]
        assert self.s.background is not None
        assert self.s.background_hash is not None
        assert "chr1:541337-541538" in self.s.background.ids

    def test04_set_meanstd(self):
        assert len(self.s.meanstd) == 0

        self.s.set_motifs("test/data/pwms/motifs.pwm")
        self.s.set_background("test/data/scan/scan_test_regions.fa")
        self.s.background["test_sequence 0.00-1.00"] = "ATTA"
        self.s.set_meanstd()
        assert len(self.s.meanstd) == 1
        assert len(self.s.meanstd["0.00-1.00"]) == 5  # number of motifs

    def test05_set_thresholds(self):
        assert self.s.threshold is None

        self.s.set_motifs("test/data/pwms/motifs.pwm")
        self.s.set_thresholds(threshold=0.0)
        assert self.s.threshold is not None
        assert self.s.fpr is None
        assert round(self.s.threshold.at["0.00-1.00", "M1500_1.01"], 4) == -17.4943

        with pytest.raises(ValueError):
            self.s.set_thresholds(fpr=0.02)  # fpr, but no background

        self.s.set_background("test/data/scan/scan_test_regions.fa")
        self.s.set_thresholds(fpr=0.02)
        assert self.s.threshold is not None
        assert self.s.fpr == 0.02
        assert self.s.threshold.shape == (58, 5)

    def test06_get_thresholds(self):
        f = Fasta(self.fa)

        self.s.set_motifs("test/data/pwms/motifs.pwm")
        self.s.set_background("test/data/scan/scan_test_regions.fa")
        self.s.set_thresholds(fpr=0.02)
        t = self.s.get_thresholds(f.seqs, zscore=False)
        assert round(t["M1500_1.01"], 4) == 6.7953
        tz = self.s.get_thresholds(f.seqs, zscore=True)
        assert round(tz["M1500_1.01"], 4) == 2.0026

    def test07_scan(self):
        """Scanner"""
        f = Fasta(self.fa)
        s = self.s
        s.set_motifs(self.motifs)
        for ncpus in [1, 2]:
            s.ncpus = ncpus

            s.set_thresholds(threshold=0.0)
            nmatches = [len(m[0]) for m in s.scan(f, 1, False)]
            self.assertEqual([1, 1, 1], nmatches)

            s.set_thresholds(threshold=0.99)
            nmatches = [len(m[0]) for m in s.scan(f.seqs, 1, False)]
            self.assertEqual([0, 1, 1], nmatches)

            s.set_thresholds(threshold=0.99)
            nmatches = [len(m[0]) for m in s.scan(f.seqs, 10, False)]
            self.assertEqual([0, 1, 2], nmatches)

            s.set_thresholds(threshold=0.99)
            nmatches = [len(m[0]) for m in s.scan(f.seqs, 10, True)]
            self.assertEqual([0, 2, 4], nmatches)

    def test08_total_count(self):
        f = Fasta(self.fa)
        self.s.set_motifs(self.motifs)
        self.s.set_thresholds(threshold=0.99)
        counts = self.s.total_count(f.seqs, 10, True)
        assert sum(counts) == 6

    def test09_scan_to_best_match(self):
        genome = os.path.join(self.data_dir, "genome.fa")
        for f in self.fa, self.bed, self.regions:
            result = scan_to_best_match(f, self.motifs, genome)
            scores = [-20.05276, 9.028887, 9.028887]
            self.assertIn("AP1", result)
            for score, match in zip(scores, result["AP1"]):
                self.assertAlmostEqual(score, match[0], 5)

    def test10_scan_to_best_score(self):
        result = scan_to_best_match(self.fa, self.motifs, score=True)
        scores = [-20.05276, 9.028887, 9.028887]
        self.assertIn("AP1", result)
        for score, match in zip(scores, result["AP1"]):
            self.assertAlmostEqual(score, match, 5)

    def tearDown(self):
        pass


if __name__ == "__main__":
    unittest.main()
