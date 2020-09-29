import unittest
import tempfile
import os
import glob
from gimmemotifs.utils import *
from gimmemotifs.fasta import Fasta
from genomepy import Genome
from tempfile import mkdtemp
from shutil import rmtree


class TestUtils(unittest.TestCase):
    """ A test class to test utils functions """

    def setUp(self):
        self.genomes_dir = "test/data/genome_index"
        self.datadir = "test/data/utils"

    def test1_phyper(self):
        """ Hypergeometric p-value """
        p = phyper(59, 500, 500, 100)
        self.assertAlmostEqual(0.02238075, p)

        p = phyper(59, 5000, 5000, 100)
        self.assertAlmostEqual(0.02782685, p)

        p = phyper(59, 50000, 50000, 100)
        self.assertAlmostEqual(0.02838217, p)

    def test2_as_fasta(self):
        """ convert bed, regions, etc to Fasta """
        tmpdir = mkdtemp()

        g = Genome("genome", genomes_dir=self.genomes_dir)

        fafile = os.path.join(self.datadir, "test.fa")
        fa = Fasta(fafile)
        bedfile = os.path.join(self.datadir, "test.bed")
        regionfile = os.path.join(self.datadir, "test.txt")
        with open(regionfile) as f:
            regions = [l.strip() for l in f]

        self.assertTrue(isinstance(as_fasta(fa), Fasta))
        self.assertTrue(isinstance(as_fasta(fafile), Fasta))

        self.assertTrue(isinstance(as_fasta(bedfile, g), Fasta))
        self.assertTrue(isinstance(as_fasta(regionfile, g), Fasta))
        self.assertTrue(isinstance(as_fasta(regions, g), Fasta))

        with self.assertRaises(ValueError):
            as_fasta(bedfile)

        rmtree(tmpdir)

    def test_checkum(self):
        fname = "test/data/fasta/test.fa"
        md5 = "a34798835d4110c34df45bbd8ed2f910"
        self.assertEqual(md5, file_checksum(fname))

    def test_join_max(self):
        size_in = range(25)
        size_out = [
            0,
            1,
            1,
            1,
            4,
            4,
            4,
            4,
            8,
            8,
            8,
            8,
            8,
            13,
            13,
            13,
            13,
            13,
            18,
            18,
            18,
            18,
            18,
            18,
            18,
        ]
        a = ["1", "22", "333", "4444", "5555"]

        for s_in, s_out in zip(size_in, size_out):
            self.assertEqual(len(join_max(a, s_in, ",")), s_out)

    def test_detect_filetype(self):
        for fname in glob.glob("test/data/filetype/*"):
            ftype = os.path.basename(fname).split(".")[0]
            self.assertEqual(ftype, determine_file_type(fname))

    def test_number_of_seqs_in_file(self):
        fnames = [
            "test/data/utils/test.fa",
            "test/data/utils/test.bed",
        ]

        for fname in fnames:
            self.assertEqual(3, number_of_seqs_in_file(fname))

    def test_median_bed_len(self):
        fname = "test/data/utils/test.narrowPeak"

        self.assertAlmostEqual(204.5, median_bed_len(fname))

    def test_check_genome(self):
        fname = "test/data/genome_index/test.bed"
        self.assertFalse(check_genome(fname))
        fname = "test/data/genome_index/genome/genome.fa"
        self.assertTrue(check_genome(fname))

    # narrowpeak_to_bed

    # get_jaspar_motif_info

    # divide_file

    # write_equalsize_bedfile

    # median_bed_len

    def tearDown(self):
        pass


if __name__ == "__main__":
    unittest.main()
