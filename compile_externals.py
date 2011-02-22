import sys
import os
from subprocess import *
from platform import platform
from distutils import log

def compile_simple(name):
	path = "src/%s" % name
	
	if not os.path.exists(path):
		return
	
	try:
		Popen("gcc", stdout=PIPE, stderr=PIPE).communicate()
	except:
		return
	
	Popen(["gcc","-lm","-o%s" % name, "%s.c" % name], cwd=path, stdout=PIPE).communicate()
	if os.path.exists(os.path.join(path, name)):
		return True

def compile_configmake(name, binary):
	path = "src/%s" % name

	if not os.path.exists(path):
		return

	Popen(["./configure"], cwd=path, stdout=PIPE, stderr=PIPE).communicate()
	Popen(["make"], cwd=path, stdout=PIPE, stderr=PIPE).communicate()

	if os.path.exists(os.path.join(path, binary)):
		return True

def compile_moan():
	# Fix the Makefile 
	from re import compile
	p = compile(r'^CC=gcc.*$')
	makefile = "src/MoAn/Makefile"	
	if not os.path.exists(makefile):
		return 

	lines = open(makefile).readlines()
	
	f = open(makefile, "w")	
	for line in lines:
		if p.match(line):
			if platform().find("fedora") != -1:
				f.write("CC=gcc34\n")
			elif platform().find("gentoo") != -1:
				f.write("CC=gcc-4.1.2\n")
			elif platform().find("Ubuntu") != -1 or platform().find("debian"):
				f.write("CC=gcc-4.1\n")
			else:
				f.write(line)
		else:
			f.write(line)
	f.close()
	
	Popen(["make"], cwd="src/MoAn/", stdout=PIPE, stderr=PIPE).communicate()
	if os.path.exists("src/MoAn/moan"):
		return True

def compile_perl(path):
	if not os.path.exists(path):
		return
	Popen(["perl", "Makefile.PL"], cwd=path, stdout=PIPE).communicate()
	Popen(["make"], cwd=path, stdout=PIPE).communicate()
	return True

def print_result(result):
	if not result:
		log.info("... failed")
	else:
		log.info("... ok")
	
def compile_all():
	
	sys.stderr.write("compiling BioProspector")
	result = compile_simple("BioProspector")
	print_result(result)

	sys.stderr.write("compiling MDmodule")
	result = compile_simple("MDmodule")
	print_result(result)
	
	sys.stderr.write("compiling MEME")
	result = compile_configmake("meme_4.6.0", "src/meme.bin")
	print_result(result)
	
	sys.stderr.write("compiling GADEM")
	result = compile_configmake("GADEM_v1.3", "src/gadem")
	print_result(result)
	
	sys.stderr.write("compiling MoAn")
	result = compile_moan()
	print_result(result)
	
	sys.stderr.write("compiling trawler dependencies")
	result = compile_perl("src/Algorithm-Cluster-1.49")
	print_result(result)

	return
