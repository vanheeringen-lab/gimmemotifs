from __future__ import print_function
from io import StringIO
import unittest
import os
import numpy as np
from gimmemotifs.motif import Motif, read_motifs
from gimmemotifs.shutils import which


class TestMotif(unittest.TestCase):
    """ A test class for Motif """

    def setUp(self):
        self.data_dir = "test/data/motif"
        self.pwm = os.path.join(self.data_dir, "test.pwm")
        self.pwm2 = "test/data/pwms/motifs.pwm"
        self.jaspar = "test/data/pwms/test.jaspar"
        self.pfm = [
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
            [1, 1, 0, 0],
            [0, 1, 1, 0],
            [0, 0, 1, 1],
            [2, 0, 2, 0],
            [3, 0, 0, 3],
            [0, 1, 0, 1],
        ]

    def test1_motif_instance(self):
        """ Creation of Motif instance """
        m = Motif()

        self.assertTrue(type(m))

    def test2_motif_instance_pfm(self):
        """ Creation of Motif instance from pfm"""
        m = Motif(self.pfm)
        self.assertTrue(m)

    def test3_motif_length(self):
        """ Motif length """
        m = Motif(self.pfm)
        self.assertEqual(10, len(m))

    def test4_motif_consensus(self):
        """ Motif length """
        m = Motif(self.pfm)
        self.assertEqual("ACGTmskrwy", m.to_consensus())

    def test5_motif_to_img(self):
        """ Motif to img """
        seqlogo = which("seqlogo")
        if seqlogo:
            m = Motif(self.pfm)
            m.to_img("test/test.png", fmt="png", seqlogo=seqlogo)
            self.assertTrue(os.path.exists("test/test.png"))
            os.unlink("test/test.png")
        else:
            print("seqlogo not found, skipping.")

    def test6_pcc(self):
        pfm1 = [[5, 0, 0, 0], [0, 5, 0, 0], [0, 5, 0, 0], [0, 0, 0, 5]]
        pfm2 = [[5, 0, 0, 0], [0, 5, 0, 0], [0, 5, 0, 0], [0, 0, 0, 5]]

        m1 = Motif(pfm1)
        m2 = Motif(pfm2)

        self.assertEqual(4, m1.max_pcc(m2)[0])

    def test7__read_motifs_pwm(self):
        with open(self.pwm2) as f:
            motifs = read_motifs(f, fmt="pwm")

        motif_ids = [m.id for m in motifs]
        self.assertEqual(5, len(motif_ids))
        self.assertEqual(
            ["M1500_1.01", "M5659_1.01", "M5669_1.01", "M5715_1.01", "M5717_1.01"],
            motif_ids,
        )

    def test7__read_motifs_jaspar(self):
        with open(self.jaspar) as f:
            motifs = read_motifs(f, fmt="jaspar")

        my_motifs = ["MA0002.2", "MA0003.3", "MA0004.1", "MA0006.1"]

        my_lens = [6, 6, 11, 11]

        motif_ids = [m.id for m in motifs]
        self.assertEqual(4, len(motif_ids))
        self.assertEqual(my_motifs, motif_ids)
        self.assertEqual(my_lens, sorted([len(m) for m in motifs]))

    def test8_pwm_to_str(self):
        pwm = [[0.01, 0.01, 0.01, 0.97], [0.123, 0.456, 0.222, 0.199]]

        m = Motif(pwm)

        s2 = "0.01\t0.01\t0.01\t0.97\n0.12\t0.46\t0.22\t0.20"
        s3 = "0.010\t0.010\t0.010\t0.970\n0.123\t0.456\t0.222\t0.199"

        self.assertEqual(s2, m._pwm_to_str(precision=2))
        self.assertEqual(s3, m._pwm_to_str(precision=3))

    def test8_pwm_to_str_hash(self):
        pwm = [[0.01, 0.01, 0.01, 0.97], [0.123, 0.456, 0.222, 0.199]]
        m = Motif(pwm)
        h = "1f260320cac8c26a"
        self.assertEqual(h, m.hash())

        pwm = [
            [0.010000, 0.010000, 0.010000, 0.970000],
            [0.12300, 0.45600, 0.22200, 0.19900],
        ]
        m = Motif(pwm)
        self.assertEqual(h, m.hash())

    def test9_logodds_matrix(self):
        pwm = [[0.5, 0.4, 0.1, 0.0], [0.25, 0.25, 0.25, 0.25]]

        logodds = np.array(
            [
                [0.69813, 0.47623, -0.89160, -4.60517],
                [0.00995, 0.00995, 0.00995, 0.00995],
            ]
        )
        m = Motif(pwm)
        np.testing.assert_almost_equal(logodds, np.array(m.logodds), decimal=5)

    def test10_read_motifs(self):

        # Read motifs from file
        motifs = read_motifs(self.pwm2, fmt="pwm")
        self.assertEqual(5, len(motifs))

        # Read motifs from file as dictionary
        motifs = read_motifs(self.pwm2, fmt="pwm", as_dict=True)
        self.assertEqual(5, len(motifs))
        self.assertEqual(type({}), type(motifs))

    def test11_slice_motif(self):
        pfm = [
            [120, 0, 0, 0],
            [120, 0, 0, 0],
            [0, 60, 60, 0],
            [0, 0, 0, 120],
            [0, 0, 0, 120],
        ]

        m = Motif(pfm)
        m.to_consensus()

        # take slice
        m2 = m[1:-1]

        self.assertEqual("AST", m2.consensus.upper())
        self.assertEqual(pfm[1:-1], m2.pfm)

    def test_motif_export_import(self):
        pfm = [
            [120, 0, 0, 0],
            [120, 0, 0, 0],
            [0, 60, 60, 0],
            [0, 0, 0, 120],
            [0, 0, 0, 120],
        ]
        motif = Motif(pfm)
        motif.id = "test_motif"

        f = StringIO(motif.to_transfac())
        motif_from_file = read_motifs(f, fmt="transfac")[0]
        self.assertEqual("AASTT", motif_from_file.to_consensus().upper())
        self.assertEqual("test_motif", motif_from_file.id)

        f = StringIO(motif.to_meme())
        motif_from_file = read_motifs(f, fmt="meme")[0]
        self.assertEqual("AASTT", motif_from_file.to_consensus().upper())
        self.assertEqual("test_motif", motif_from_file.id)

        f = StringIO(motif.to_motevo())
        motif_from_file = read_motifs(f, fmt="transfac")[0]
        self.assertEqual("AASTT", motif_from_file.to_consensus().upper())
        self.assertEqual("test_motif", motif_from_file.id)

    def test_motif_from_alignment(self):
        align = "AACTT\n" "AAGTA\n" "AACTC\n" "AAGTG\n"
        f = StringIO(align)
        motif = read_motifs(f, fmt="align")[0]

        self.assertEqual("AASTN", motif.to_consensus().upper())

    def test_read_motifs_xxmotifs(self):
        fname = "test/data/motifprogram/xxmotif.pwm"
        motifs = read_motifs(fname, fmt="xxmotif")

        self.assertEqual(4, len(motifs))
        self.assertEqual(9, len(motifs[-1]))
        self.assertEqual("RGGCAWGYC", motifs[-1].to_consensus().upper())

    def tearDown(self):
        pass


if __name__ == "__main__":
    unittest.main()
