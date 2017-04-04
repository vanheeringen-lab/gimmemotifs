import unittest
import tempfile
import os
from gimmemotifs.stats import *
from time import sleep

class TestStats(unittest.TestCase):
    """ A test class to test motif stats """

    def setUp(self):
        self.data_dir = "test/data/stats"
        
        self.motifs = os.path.join(self.data_dir, "motifs.pwm")
        self.fg_fa = os.path.join(self.data_dir, "p73.fa")
        self.bg_fa = os.path.join(self.data_dir, "random.w200.fa")
        self.stat_functions = [
            "recall_at_fdr",
            "fraction_fdr",
            "score_at_fdr",
            "enr_at_fdr",
            "max_enrichment",
            "mncp",
            "roc_auc",
            "roc_auc_xlim",
            "max_fmeasure",
            "ks_pvalue",
            "ks_significance",  
			]

    def test1_stats(self):
        """ Calculate motif statistics """
        stats = calc_stats(self.motifs, self.fg_fa, self.bg_fa)
        
        for f in self.stat_functions:
			self.assertIn(f, stats.values()[0])
        
        # Two motifs
        self.assertEquals(2, len(stats))

        m1 = "T-box_M1713_1.01" # not enriched
        m2 = "p53_Average_8"    # highly enriched

        self.assertLess(stats[m1]["roc_auc"] , 0.9)
        self.assertGreater(stats[m2]["roc_auc"] , 0.5)

        self.assertEquals(stats[m1]["recall_at_fdr"] , 0.0)
        self.assertGreater(stats[m2]["recall_at_fdr"] , 0.8)
    
        self.assertGreater(stats[m1]["ks_pvalue"] , 0.01)
        self.assertLess(stats[m2]["ks_pvalue"] , 0.001)
        
        self.assertGreater(stats[m1]["phyper_at_fdr"] , 0.1)
        self.assertLess(stats[m2]["phyper_at_fdr"] , 1e-16)
        
        # Only calculate specific statistic
        stats = calc_stats(self.motifs, self.fg_fa, self.bg_fa, ["roc_auc"])
        
        self.assertEquals(1, len(stats.values()[0]))
        
        self.assertLess(stats[m1]["roc_auc"] , 0.9)
        self.assertGreater(stats[m2]["roc_auc"] , 0.5)
    
    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
