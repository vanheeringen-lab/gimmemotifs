import platform
import unittest
import tempfile
import os
from gimmemotifs.tools import __tools__, get_tool
from gimmemotifs.comparison import MotifComparer
from gimmemotifs.motif import motif_from_consensus

class TestMotifProgram(unittest.TestCase):
    """ A test class for motifprograms """

    def setUp(self):
        self.data_dir = "test/data/motifprogram"
        self.fa = os.path.join("test/data/denovo/input.fa")
        self.bg_fa = os.path.join("test/data/denovo/random.fa")

    def isInstalled(self, bin):
        return os.path.exists(bin) and os.access(bin, os.X_OK)
    
    def ap1_included(self, motifs):
        ap1 = motif_from_consensus("TGASTCA")
        mc = MotifComparer()
        for motif in motifs:
            match = mc.get_closest_match(ap1, motif)
            if match["TGASTCA"][1][3] < 1e-6:
                return True
        return False

    def test_tools(self):
        """Test motif prediction tools."""
        params = {
                "background": self.bg_fa,
                "organism": "hg38",
                "width": 7,
                }
        
        for tool_name in sorted(__tools__):
            if tool_name in [
                    "gadem",  # sometimes crashes on invalid pointer
                    "jaspar", 
                    "xxmotif", # takes too long 
                    ]:
                continue
           
            if platform.system() == "Darwin":
                # No support for osx
                if tool_name in [
                    "amd",
                    "hms",
                    "improbizer",
                    "motifsampler",
                    "posmo",
                    "meme", # until the bioconda package is fixed
                    "memew", # until the bioconda package is fixed
                    ]:
                    continue

            t = get_tool(tool_name)
            print("Testing {}...".format(t))
            
            (motifs, stderr, stdout) =  t.run(self.fa, params)
            self.assertTrue(self.ap1_included(motifs))

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
