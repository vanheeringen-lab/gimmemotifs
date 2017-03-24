import unittest
import tempfile
import os
from gimmemotifs.scanner import Scanner

class TestScanner(unittest.TestCase):
    """ A test class to test Motif pwmscan functionality and related things """

    def setUp(self):
        pass

    def testThreshold(self):
        s = Scanner()
        s.set_motifs("test/data/pwms/motifs.pwm")
        
        fname = "test/data/scan/scan_test_regions.fa"
        s.set_threshold(0.02, filename=fname)


    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
