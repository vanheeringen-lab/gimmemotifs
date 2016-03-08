import unittest
import tempfile
import os
from gimmemotifs.cluster import cluster_motifs

class TestMotifPwm(unittest.TestCase):
    """ A test class to test motif clustering """

    def setUp(self):
        self.data_dir = "test/data/pwms"
        
        self.pwm = os.path.join(self.data_dir, "motifs.pwm")
    
    def test1_cluster_motifs(self):
        """ cluster a pwm file with motifs """
        # Run clustering
        tree = cluster_motifs(self.pwm, 
                                "total",
                                "wic",
                                "mean",
                                True,
                                threshold=0.95,
                                include_bg=True)
        
        clusters = tree.getResult()

        self.assertEquals(2, len(clusters))
        self.assertEquals([3,2], [len(c[1]) for c 
            in sorted(clusters, cmp=lambda x,y: cmp(len(x), len(y)))])

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
