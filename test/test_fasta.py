from gimmemotifs.fasta import *
import tempfile
import unittest
import os

class TestFasta(unittest.TestCase):
    """ A test class for Fasta """

    def setUp(self):
        self.fasta_file = "test/data/fasta/test.fa"
        self.assertTrue(os.path.exists(self.fasta_file))
        self.assertTrue(Fasta(self.fasta_file))
        self.f = Fasta(self.fasta_file)
    
    def test1_index(self):
        """ Fasta as a dictionary """
        self.assertEqual(self.f["seq1"], "AAAA")
        self.assertEqual(self.f["seq2"], "ACGT")
        self.assertEqual(self.f["seq3"], "CCCCGGGG")
    
    def test2_items(self):
        """ Fasta.items() """
        self.assertEqual(len(list(self.f.items())), 3)
    
    def test3_writefasta(self):
        """ Write fasta-formatted file"""
        temp = tempfile.NamedTemporaryFile()
        tempname = temp.name
        self.f.writefasta(tempname)
        with open(self.fasta_file) as f:
            with open(tempname) as f_ref:
                self.assertEqual(f.read().strip(), f_ref.read().strip())
    
    def tearDown(self):
            pass

if __name__ == '__main__':
    unittest.main()
