from __future__ import print_function
import unittest
import tempfile
import os
from gimmemotifs.motif import *
from gimmemotifs.shutils import which

class TestMotif(unittest.TestCase):
    """ A test class for Motif """

    def setUp(self):
        self.data_dir = "test/data/motif"
        self.pwm = os.path.join(self.data_dir, "test.pwm")
        self.pwm2 = "test/data/pwms/motifs.pwm"
        self.jaspar = "test/data/pwms/test.jaspar"
        self.pfm = [ 
            [1,0,0,0],
            [0,1,0,0],
            [0,0,1,0],
            [0,0,0,1],
            [1,1,0,0],
            [0,1,1,0],
            [0,0,1,1],
            [2,0,2,0],
            [3,0,0,3],
            [0,1,0,1]
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
        pfm1 = [[5,0,0,0],[0,5,0,0],[0,5,0,0],[0,0,0,5]]
        pfm2 = [[5,0,0,0],[0,5,0,0],[0,5,0,0],[0,0,0,5]]

        m1 = Motif(pfm1)
        m2 = Motif(pfm2)

        self.assertEqual(4, m1.max_pcc(m2)[0])
    
    def test7__read_motifs_pwm(self):
        with open(self.pwm2) as f:
            motifs = read_motifs(f, fmt="pwm")

        motif_ids = [m.id for m in motifs]
        self.assertEqual(5, len(motif_ids))
        self.assertEqual(["M1500_1.01","M5659_1.01","M5669_1.01","M5715_1.01", "M5717_1.01"], motif_ids)
     
    def test7__read_motifs_jaspar(self):
        with open(self.jaspar) as f:
            motifs = read_motifs(f, fmt="jaspar")

        my_motifs = [
                "MA0002.2\tRUNX1",
                "MA0003.3\tTFAP2A",
                "MA0004.1\tArnt",
                "MA0006.1\tAhr::Arnt"
        ]

        my_lens = [6,6,11,11]
        
        motif_ids = [m.id for m in motifs]
        self.assertEqual(4, len(motif_ids))
        self.assertEqual(my_motifs, motif_ids)
        self.assertEqual(my_lens, sorted([len(m) for m in motifs]))
    
    def test8_pwm_to_str(self):
        pwm = [
            [0.01, 0.01, 0.01, 0.97],
            [0.123, 0.456, 0.222, 0.199],
            ]
        
        m = Motif(pwm)

        s2 = "0.01\t0.01\t0.01\t0.97\n0.12\t0.46\t0.22\t0.20"
        s3 = "0.010\t0.010\t0.010\t0.970\n0.123\t0.456\t0.222\t0.199"

        self.assertEqual(s2, m._pwm_to_str(precision=2))
        self.assertEqual(s3, m._pwm_to_str(precision=3))
    
    def test8_pwm_to_str(self):
        pwm = [
            [0.01, 0.01, 0.01, 0.97],
            [0.123, 0.456, 0.222, 0.199],
            ]
        m = Motif(pwm)
        h = "1f260320cac8c26a"
        self.assertEqual(h, m.hash())
        
        pwm = [
            [0.010000, 0.010000, 0.010000, 0.970000],
            [0.12300, 0.45600, 0.22200, 0.19900],
            ]
        m = Motif(pwm)
        self.assertEqual(h, m.hash())

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
