import unittest
import tempfile
import os
from gimmemotifs.stats import calc_stats
from gimmemotifs.motif import read_motifs
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
            "fraction_fpr",
            "score_at_fpr",
            "enr_at_fpr",
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
        for ncpus in [1,2]:
            stats = calc_stats(self.motifs, self.fg_fa, self.bg_fa, ncpus=ncpus)
            
            for f in self.stat_functions:
                self.assertIn(f, list(stats.values())[0])
            
            # Two motifs
            self.assertEqual(2, len(stats))
    
            m1 = "T-box_M1713_1.01_CTAGGTGTGAA" # not enriched
            m2 = "p53_Average_8_CATGyCnGGrCATGy"    # highly enriched
    
            self.assertLess(stats[m1]["roc_auc"] , 0.9)
            self.assertGreater(stats[m2]["roc_auc"] , 0.5)
    
            self.assertEqual(stats[m1]["recall_at_fdr"] , 0.0)
            self.assertGreater(stats[m2]["recall_at_fdr"] , 0.8)
        
            self.assertGreater(stats[m1]["ks_pvalue"] , 0.01)
            self.assertLess(stats[m2]["ks_pvalue"] , 0.001)
            
            self.assertGreater(stats[m1]["phyper_at_fpr"] , 0.1)
            self.assertLess(stats[m2]["phyper_at_fpr"] , 1e-16)
            
            # Only calculate specific statistic
            stats = calc_stats(self.motifs, self.fg_fa, self.bg_fa, stats=["roc_auc"])
            
            self.assertEqual(1, len(list(stats.values())[0]))
            
            self.assertLess(stats[m1]["roc_auc"] , 0.9)
            self.assertGreater(stats[m2]["roc_auc"] , 0.5)
     
    def test2_stats_single_motif(self):
        """ Calculate motif statistics """
        
        m_id = "p53_Average_8_CATGyCnGGrCATGy"
        
        with open(self.motifs) as f:
            motifs = read_motifs(f)
        motif = [m for m in motifs if str(m) == m_id][0]
        
        stats = calc_stats(motif, self.fg_fa, self.bg_fa, stats=["roc_auc"])
        self.assertGreater(stats[m_id]["roc_auc"] , 0.9)
    
    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
