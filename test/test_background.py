import unittest
import tempfile
import os
from collections import Counter
from gimmemotifs.background import matched_gc_bedfile
from pybedtools import BedTool

class TestBackground(unittest.TestCase):
    """ A test class to test generation of background sequences"""

    def setUp(self):
        self.data_dir = "test/data/background"
        self.genome = "test/data/background/genome.fa"

    def test1_gc_matched_background(self):
        """ Generate GC matched background """
        tmp = tempfile.NamedTemporaryFile(delete=False)
    
        gc_bed1 = os.path.join( self.data_dir, "test_gc1.bed")
        gc_bed2 = os.path.join( self.data_dir, "test_gc2.bed")
        gc_bed3 = os.path.join( self.data_dir, "test_gc3.bed")
        
        for bed in [gc_bed1, gc_bed2, gc_bed3]:
            matched_gc_bedfile(tmp.name, bed, self.genome, 10)
            b = BedTool(tmp.name)
            gc = [f[4] for f in b.nucleotide_content(fi=self.genome)]
            gc = [str(round(float(x), 1)) for x in gc]
            c = Counter(gc)
            
            self.assertEquals(5, c["0.2"])
            self.assertEquals(5, c["0.4"])

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
