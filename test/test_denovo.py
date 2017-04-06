import unittest
import tempfile
import os
from gimmemotifs.denovo import gimme_motifs
from gimmemotifs.motif import read_motifs, motif_from_consensus
from gimmemotifs.comparison import MotifComparer
import shutil
from tempfile import mkdtemp
from time import sleep

class TestScanner(unittest.TestCase):
    """ A test class to test scanner funcitonalitu """

    def setUp(self):
        self.outdir = mkdtemp()

    def test1_denovo(self):
        """ de novo motif prediction """
       
        
        gimme_motifs("test/data/denovo/input.fa", self.outdir,
            params={
                "tools":"Homer,MDmodule,BioProspector", 
                "fraction":0.5,
                "background":"random"
                },
            filter_significant=True,
            cluster=True)
       
        fnames = ["motifs.pwm", "motif_report.html", "cluster_report.html",
                    "params.txt", "stats.random.txt"]

        # Check if all output files are there
        for fname in fnames:
            self.assertTrue(os.path.exists(os.path.join(self.outdir, fname)))   
  
        # Check if correct motif is predicted
        predicted_motifs = read_motifs(open(os.path.join(self.outdir, "motifs.pwm")))
        ap1 = motif_from_consensus("TGASTCA")

        mc = MotifComparer()
        ap1_predicted = False
        for motif in predicted_motifs:
            match = mc.get_closest_match(ap, p)
            if match["TGASTCA"][3] < 1e-6:
                ap1_predicted = True
                break

        self.assertTrue(ap1_predicted)


    def tearDown(self):
        # remove output
        shutil.rmtree(self.outdir) 

if __name__ == '__main__':
    unittest.main()

