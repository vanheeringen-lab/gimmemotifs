import os
import re
import sys
import unittest
from glob import glob

libdir = glob('build/lib.*')[0]
p = [os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]), libdir)), "test"]
sys.path = p + sys.path

p = re.compile(r'^test_\w+\.py$')
tests = [file[:-3] for file in os.listdir("test") if p.match(file)]
suite = unittest.TestSuite()

for test in tests:
	suite.addTest(unittest.defaultTestLoader.loadTestsFromName(test))

runner = unittest.TextTestRunner(descriptions=1, verbosity=2)
result = runner.run(suite)
