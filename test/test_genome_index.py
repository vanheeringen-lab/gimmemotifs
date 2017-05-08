import unittest
import tempfile
import os
import glob
from shutil import rmtree
from gimmemotifs.genome_index import *

class TestGenomeIndex(unittest.TestCase):
    """ A test class for GenomeIndex class """

    def setUp(self):
        self.g = GenomeIndex()
        self.fasta_dir = "test/data/genome_index"
        self.index_dir = tempfile.mkdtemp()
        self.temp = tempfile.NamedTemporaryFile()
        self.temp_file = self.temp.name
    
    def test_single_index(self):
        """ index a single fasta-file """
        # First fasta-file 
        fasta_file = os.path.join(self.fasta_dir, "chr1.fa")
        self.assertTrue(os.path.exists(fasta_file))
        self.g._make_index(fasta_file, self.temp_file)
        # Check if index exists
        self.assertTrue(os.path.exists(self.temp_file))
        # Read index file and check length
        with open(self.temp_file) as f:
            index_size = len(f.read())
        self.assertEqual(index_size, 7 * len(pack(self.g.pack_char, 0)))
        
        # Second fasta-file 
        fasta_file = os.path.join(self.fasta_dir, "chr3.fa")
        self.assertTrue(os.path.exists(fasta_file))
        self.g._make_index(fasta_file, self.temp_file)
        # Check if index exists
        self.assertTrue(os.path.exists(self.temp_file))
        # Read index file and check length
        with open(self.temp_file) as f:
            index_size = len(f.read())
        self.assertEqual(index_size, 4 * len(pack(self.g.pack_char, 0)))

    def test_create_index(self):
        """ index a directory of fasta-files """
        self.g.create_index(self.fasta_dir, self.index_dir)
        self.assertTrue(os.path.exists(os.path.join(self.index_dir, "chr1.index")))
        self.assertTrue(os.path.exists(os.path.join(self.index_dir, "chr2.index")))
        self.assertTrue(os.path.exists(os.path.join(self.index_dir, "chr3.index")))
        self.assertTrue(os.path.exists(os.path.join(self.index_dir, self.g.param_file)))

    def test_get_sequence(self):
        """ get_sequence should correctly retrieve sequences from the index """
        # Setup index
        self.g.create_index(self.fasta_dir, self.index_dir)
        # Check first line
        self.assertEqual(self.g.get_sequence("chr1", 0, 4), "AAAA")
        # Check second line
        self.assertEqual(self.g.get_sequence("chr1", 4, 8), "CCCC")
        # Check sequence spanning two lines
        self.assertEqual(self.g.get_sequence("chr1", 2, 6), "AACC")
        # Check complete sequence
        self.assertEqual(self.g.get_sequence("chr1", 0, 26), "AAAACCCCGGGGTTTTAAAACCCCGG")
        
    def test_get_sequence_bad_input(self):
        """ get_sequence should fail with invalid coordinates"""
        self.g.create_index(self.fasta_dir, self.index_dir)
        # Check coordinates greater than sequence end
        self.assertRaises(ValueError, self.g.get_sequence, "chr1", 0, 100)
        # Check coordinates greater than sequence end
        self.assertRaises(ValueError, self.g.get_sequence, "chr1", 100, 100)
        # Check negative start
        self.assertRaises(ValueError, self.g.get_sequence, "chr1", -100, 3)
        
    def test_track2fasta(self):    
        """ track2fasta should successfully convert bed-file to fasta-file """
        self.g.create_index(self.fasta_dir, self.index_dir)
        result = os.path.join(self.fasta_dir, "track2fasta_result")
        bedfile = os.path.join(self.fasta_dir, "test.bed")
        track2fasta(self.index_dir, bedfile, self.temp_file)
        with open(result) as f_new:
            with open(self.temp_file) as f_ref:
                self.assertEqual(f_new.read(), f_ref.read())
        ## Test track2fasta with extend
        #self.assertEqual(self.g.get_sequence("chr1", 2, 5, extend_up=2), "AAAACC")
        #self.assertEqual(self.g.get_sequence("chr1", 2, 5, extend_down=2), "AACCCC")
        #self.assertEqual(self.g.get_sequence("chr1", 2, 5, extend_up=2, extend_down=2), "AAAACCCC")
        
    def test_track2fasta_exons(self):
        """ track2fasta should convert bed12 to fasta"""
        from gimmemotifs.fasta import Fasta
        bedfile = os.path.join(self.fasta_dir, "genes.bed")
        fafile = os.path.join(self.fasta_dir, "genes.out")
        
        # Create index
        self.g.create_index(self.fasta_dir, self.index_dir)
        # Convert bed to fasta
        track2fasta(self.index_dir, bedfile, self.temp_file, use_strand=True)
        target = Fasta(fafile)
        test = Fasta(self.temp_file)
        for gene in test.ids:
            name = gene.split(" ")[-1]
            self.assertEqual(len(test[gene]), len(target[name]))
            self.assertEqual(test[gene].upper(), target[name].upper())
    
    def test_get_genome(self):
        """ test automatic install of genome """
        # pretty small genome
        genome = "sacCer3"
        fadir = tempfile.mkdtemp(prefix="gimme.")
        genome_dir = os.path.join(fadir, genome)
        index_dir = tempfile.mkdtemp(prefix="gimme.")
        get_genome(genome, fadir, index_dir)
        self.assertEqual(17, len(glob(os.path.join(genome_dir, "*.fa*"))))
        
        for d in fadir, index_dir:
            rmtree(d)
    
    def tearDown(self):
        for file in os.listdir(self.index_dir):
            os.remove(os.path.join(self.index_dir, file))
        os.rmdir(self.index_dir)


if __name__ == '__main__':
    unittest.main()
