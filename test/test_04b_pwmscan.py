import unittest
from tempfile import NamedTemporaryFile
import os
from gimmemotifs.motif import read_motifs
from gimmemotifs.fasta import Fasta


class TestMotifPwm(unittest.TestCase):
    """A test class to test Motif pwmscan functionality and related things"""

    def setUp(self):
        self.data_dir = "test/data/pwmscan"

        with open(os.path.join(self.data_dir, "TATA.pwm")) as f:
            self.motif = read_motifs(f, fmt="pwm")[0]
        self.prom = Fasta(os.path.join(self.data_dir, "promoters.fa"))
        self.prom_gff = os.path.join(self.data_dir, "promoters_result.gff")
        self.random = Fasta(os.path.join(self.data_dir, "random_sequences.fa"))
        self.random_gff = os.path.join(self.data_dir, "random_result.gff")
        self.enrichment = os.path.join(self.data_dir, "enrichment.txt")
        self.tmp = NamedTemporaryFile().name

    def test1_scan(self):
        """Scan a FASTA file with PWM of motif"""
        result = self.motif.scan(self.prom, nreport=1)

        # Every sequence should have a TATA match
        self.assertEqual(len(result), len(self.prom))

    def test2_pwm_scan_to_gff(self):
        """Scan a FASTA file with PWM of motif, and produce GFF"""

        self.motif.pwm_scan_to_gff(self.prom, self.tmp)
        with open(self.tmp) as f:
            for line in f:
                vals = line.strip().split("\t")
                self.assertEqual(9, len(vals))
                self.assertTrue(int(vals[3]) > 0)
                self.assertTrue(int(vals[4]) > 0)
                self.assertTrue(float(vals[5]) > 5.25)
                self.assertTrue(float(vals[5]) < 9.06)
                self.assertIn(vals[6], ["+", "-"])

    def tearDown(self):
        pass


if __name__ == "__main__":
    unittest.main()
