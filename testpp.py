import pp
import subprocess
import time
import sys

def run_long():
	bla = subprocess.Popen(["moan", "p63_prediction.fa","p63_bg_genomic_matched.fa"], stdout=subprocess.PIPE).communicate()[0]
	return "long"

def run_middle():
	p = subprocess.Popen(["p73scan.py","-i","p63_bg_genomic_matched.fa",">","bla"], stdout=subprocess.PIPE).communicate()
	return "middle"

def run_short():
	bla = subprocess.Popen(["ls"], stdout=subprocess.PIPE).communicate()[0]
	return "short"

done = []
def callback(arg):
	done.append(arg)
	

test = pp.Server()
jobs = []
jobs.append( test.submit(run_long, (), (), ("subprocess",), callback) )
jobs.append( test.submit(run_short, (), (), ("subprocess",), callback) )
jobs.append( test.submit(run_short, (), (), ("subprocess",), callback ) )
jobs.append( test.submit(run_short, (), (), ("subprocess",), callback ) )
jobs.append( test.submit(run_middle, (), (), ("subprocess",), callback ) )

#print Popen(["p73scan.py","-i", "p63_bg_genomic_matched.fa"], stdout=PIPE).communicate()
try:
	print done
	time.sleep(5)
	print done
	test.destroy()
	time.sleep(5)
	print done
	time.sleep(50)
	print done
except KeyboardInterrupt, e:
	print "Done so far", done
