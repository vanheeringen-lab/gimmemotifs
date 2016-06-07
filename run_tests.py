import os
import re
import sys
import unittest
from glob import glob

def get_tests():
    start_dir = os.path.join(os.path.dirname(__file__), "test")
    return unittest.TestLoader().discover(start_dir, pattern="test_*.py")

#libdirs = glob('build/lib.*')
#if len(libdirs) > 0:
#    p = [os.path.abspath(os.path.join(
#        os.path.dirname(sys.argv[0]), libdirs[0])
#        ), "test"]
#    sys.path = p + sys.path

suite = unittest.TestSuite()

for test in get_tests():
	suite.addTest(test)

runner = unittest.TextTestRunner(descriptions=1, verbosity=2)
ret = not runner.run(suite).wasSuccessful()
sys.exit(ret)
