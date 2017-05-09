import unittest
import tempfile
import os
from gimmemotifs.prediction import PredictionResult
from gimmemotifs.motif import read_motifs
from gimmemotifs import rocmetrics
from time import sleep

class TestStats(unittest.TestCase):
    """ A test class to test motif stats """

    def setUp(self):
        self.data_dir = "test/data/stats"
        
        self.motifs = os.path.join(self.data_dir, "motifs.pwm")
        self.fg_fa = os.path.join(self.data_dir, "p73.fa")
        self.bg_fa = os.path.join(self.data_dir, "random.w200.fa")

    def test1_prediction_result(self):
        """ Calculates statistics of motifs """
        tmp = tempfile.NamedTemporaryFile().name
        
        p = PredictionResult(tmp, fg_file=self.fg_fa, background={"random":self.bg_fa})

        with open(self.motifs) as f:
            motifs = read_motifs(f)
        p.add_motifs((0, (motifs, "", "")))
        p.wait_for_stats()        
        self.assertEqual(2, len(p.stats))
        #for stat in rocmetrics.__all__:
        #    self.assertIn(stat, p.stats.values()[0]["random"])

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
