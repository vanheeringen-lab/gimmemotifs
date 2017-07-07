import unittest
import tempfile
import os
from gimmemotifs.utils import *
from gimmemotifs.fasta import Fasta
from gimmemotifs.genome_index import GenomeIndex
from tempfile import mkdtemp
from shutil import rmtree

class TestUtils(unittest.TestCase):
    """ A test class to test utils functions """

    def setUp(self):
        self.genome_dir = "test/data/genome_index"
        self.datadir = "test/data/utils"
    
    def test1_phyper(self):
        """ Hypergeometric p-value """
        p = phyper(59, 500,500, 100)
        self.assertAlmostEqual(0.02238075, p)
        
        p = phyper(59, 5000, 5000, 100)
        self.assertAlmostEqual(0.02782685, p)

        p = phyper(59, 50000, 50000, 100)
        self.assertAlmostEqual(0.02838217, p)

    def test2_as_fasta(self):
        """ convert bed, regions, etc to Fasta """
        tmpdir = mkdtemp()

        g = GenomeIndex()
        g.create_index(self.genome_dir, tmpdir)

        fafile = os.path.join(self.datadir, "test.fa")
        fa = Fasta(fafile)
        bedfile = os.path.join(self.datadir, "test.bed")
        regionfile = os.path.join(self.datadir, "test.txt")
        with open(regionfile) as f:
            regions = [l.strip() for l in f]

        self.assertTrue(isinstance(as_fasta(fa), Fasta))
        self.assertTrue(isinstance(as_fasta(fafile), Fasta))

        self.assertTrue(isinstance(as_fasta(bedfile, tmpdir), Fasta))
        self.assertTrue(isinstance(as_fasta(regionfile, tmpdir), Fasta))
        self.assertTrue(isinstance(as_fasta(regions, tmpdir), Fasta))
        
        with self.assertRaises(ValueError):
            as_fasta(bedfile)
        
        rmtree(tmpdir)
    
    def test_checkum(self):
        fname = "test/data/fasta/test.fa"
        md5 = "a34798835d4110c34df45bbd8ed2f910"
        self.assertEqual(md5, file_checksum(fname))

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
