import unittest
import tempfile
import os
from gimmemotifs.comparison import MotifComparer
from time import sleep

class TestComparison(unittest.TestCase):
    """ A test class to test comparison funcitonality """

    def setUp(self):
        pass        

    def test1_closest_match(self):
        """ Closest match """
        mc = MotifComparer()

        pwm = "test/data/pwmscan/TATA.pwm"
        ret = mc.get_closest_match(pwm)
        
        self.assertIn('TATA-box', ret)
        
        match = ret['TATA-box']
        self.assertEqual('TBP_Average_1', match[0])
        
        scores = match[1]
        self.assertAlmostEqual(-0.3276, scores[0], 4)
        self.assertEqual(-1, scores[1])
        self.assertEqual(1, scores[2])
        self.assertAlmostEqual(4.134e-7, scores[3])

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
