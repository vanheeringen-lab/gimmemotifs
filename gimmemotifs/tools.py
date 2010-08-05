# Copyright (c) 2009-2010 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

""" Interface module for all motif programs """

# Python imports
import re
import os
import sys
import logging
from math import log,sqrt
from subprocess import *
from tempfile import NamedTemporaryFile,mkdtemp
import shutil
from string import maketrans
# External imports
import pp
# GimmeMotifs imports
from gimmemotifs.motif import * 

class MotifProgram:
	from gimmemotifs.config import MotifConfig
	config = MotifConfig()

	def __init__(self):
		pass	
	
	def bin(self):
		return self.config.bin(self.name)

	def dir(self):
		return self.config.dir(self.name)

	def is_configured(self):
		return self.config.is_configured(self.name)
	
	def is_installed(self):
		return self.is_configured() and os.access(self.bin(), os.X_OK)

	def run(self, fastafile, savedir, params={}):
		if not self.is_configured():
			raise ValueError, "%s is not configured" % self.name

		if not self.is_installed():
			raise ValueError, "%s is not installed or not correctly configured" % self.name

		try:
			return self._run_program(self.bin(), fastafile, savedir, params)
		except KeyboardInterrupt:
			return ([], "Killed", "Killed")

class BioProspector(MotifProgram):
	def __init__(self):
		self.name = "BioProspector"
		self.cmd = "BioProspector"

	def _run_program(self, bin, fastafile, savedir="", params={}):
		import os, tempfile, shutil
		from subprocess import Popen, PIPE
		
		default_params = {"single":False, "background":None, "analysis":"medium", "number":5, "width":10}
		default_params.update(params)
		
		prospector = bin
		
		fastafile = os.path.abspath(fastafile)
		
		# Background file is essential!
		if not default_params["background"]:
			print "Background file needed!"
			sys.exit()
		
		bgfile = os.path.abspath(default_params["background"])
		tmpdir = tempfile.mkdtemp()
		outfile = os.path.join(tmpdir, "bioprospector.out")	
		
		stdout = ""
		stderr = ""
		
		strand = 2
		if default_params["single"]:
			strand = 1

		cmd = "%s -i %s -W %s -d %s -b %s -r %s -o %s" % (
			prospector,
			fastafile,
			default_params["width"],
			strand,
			default_params["background"],
			default_params["number"],
			outfile)

		p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
		out,err = p.communicate()
		stdout += out
		stderr += err
		
		motifs = []
		
		if os.path.exists(outfile):
			f = open(outfile)
			motifs = self.parse(f)
			f.close()
		
		# remove temporary files
		if os.path.exists(tmpdir):
			shutil.rmtree(tmpdir)
		
		return motifs, stdout, stderr

	def parse(self, fo):
		motifs = []
		
		p = re.compile(r'^\d+\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)')
		align = []
		pwm = []
		id = ""
		for line in fo.readlines():
			if line.startswith("Motif #"):
				if pwm:
					m = Motif(pwm)
					m.id = "BioProspector_w%s_%s" % (len(m), id)
					motifs.append(m)
				id =  line.split("#")[1].split(":")[0]
				pwm = []
			else:
				m = p.search(line)
				if m:
					pwm.append([float(m.group(x))/100.0 for x in range(1,5)])

		if pwm:
			m = Motif(pwm)
			m.id = "BioProspector_w%s_%s" % (len(m), id)
			motifs.append(m)
		return motifs


class MoAn(MotifProgram):
	def __init__(self):
		self.name = "MoAn"
		self.cmd = "moan"

	def _run_program(self, bin, fastafile, savedir="", params={}):
		import os, tempfile, shutil
		from subprocess import Popen, PIPE
		
		default_params = {"single":False, "background":None, "analysis":"medium"}
		default_params.update(params)
		
		moan = bin
		
		fastafile = os.path.abspath(fastafile)
		
		# Background file is essential!
		if not default_params["background"]:
			print "Background file needed!"
			sys.exit()
		
		bgfile = os.path.abspath(default_params["background"])
		tmpdir = tempfile.mkdtemp()
		outfile = os.path.join(tmpdir, "moan.out")	
		
		current_path = os.getcwd()
		os.chdir(tmpdir)
		
		stdout = ""
		stderr = ""
	
		ranges  = {"xs":"5,5","small":"5,8", "medium":"5,10","large":"5,14", "xl":"5,20"}
		r = ranges[default_params["analysis"]]

		cmd = "%s -D -R %s %s %s > %s" % (moan, r, fastafile, bgfile, outfile)
		if default_params["single"]:
			cmd = "%s -R %s %s %s > %s" % (moan, r, fastafile, bgfile, outfile)
		p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
		out,err = p.communicate()
		stdout += out
		stderr += err
		
		os.chdir(current_path)
		motifs = []
		if os.path.exists(outfile):
			f = open(outfile)
			motifs = self.parse(f)
			f.close()
		
		print outfile
		# remove temporary files
		#if os.path.exists(tmpdir):
		#	shutil.rmtree(tmpdir)
		
		return motifs, stdout, stderr

	def parse(self, fo):
		motifs = []
		
		#160:  112  CACGTGC      7.25   chr14:32308489-32308689
		p = re.compile(r'\d+:\s+\d+\s+(\w+)\s+\d+\.\d+')
		align = []
		for line in fo.readlines():
			m = p.search(line)
			if m:
				align.append(m.group(1))

		if align:
			motifs = [motif_from_align(align)]
			motifs[0].id = "MoAn"	
		return motifs


class Improbizer(MotifProgram):
	def __init__(self):
		self.name = "Improbizer"
		self.cmd = "ameme"

	def _run_program(self, bin, fastafile, savedir="", params={}):
		import os, tempfile, shutil
		from subprocess import Popen, PIPE
		
		default_params = {"background":None, "number":10}
		default_params.update(params)
		
		ameme = bin
		
		fastafile = os.path.abspath(fastafile)
		
		# Not strictly necessary, but recommended
		if not default_params["background"]:
			print "Background file needed!"
			sys.exit()
		
		bgfile = os.path.abspath(default_params["background"])
		tmpdir = tempfile.mkdtemp()
		outfile = os.path.join(tmpdir, "improbizer.out.html")	
		
		current_path = os.getcwd()
		os.chdir(tmpdir)
		
		stdout = ""
		stderr = ""
		cmd = "%s good=%s bad=%s numMotifs=%s > %s" % (ameme, fastafile, bgfile, default_params["number"], outfile)
		p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
		out,err = p.communicate()
		stdout += out
		stderr += err
		
		motifs = []
		if os.path.exists(outfile):
			f = open(outfile)
			motifs = self.parse(f)
			f.close()
		
		os.chdir(current_path)
		# remove temporary files
		if os.path.exists(tmpdir):
			shutil.rmtree(tmpdir)
		
		return motifs, stdout, stderr

	def parse(self, fo):
		motifs = []
		p = re.compile(r'\d+\s+@\s+\d+\.\d+\s+sd\s+\d+\.\d+\s+(\w+)$')

		line = fo.readline()
		while line and line.find("Color") == -1:
			m = p.search(line)
			if m:
				pwm_data = {}
				for i in range(4):
					vals = [x.strip() for x in fo.readline().strip().split(" ") if x]
					pwm_data[vals[0].upper()] = vals[1:]
				pwm = []
				for i in range(len(pwm_data["A"])):
					pwm.append([float(pwm_data[x][i]) for x in ["A","C","G","T"]])
				motifs.append(Motif(pwm))
				motifs[-1].id = ">Improbizer_%s" % m.group(1)
			line = fo.readline()
		
		return motifs

class Trawler(MotifProgram):
	def __init__(self):
		self.name = "trawler"
		self.cmd = "trawler.pl"

	def _run_program(self, bin, fastafile, savedir="", params={}):
		import os, tempfile, shutil
		from subprocess import Popen, PIPE
		
		default_params = {"single":False, "background":None}
		default_params.update(params)
		
		trawler = bin
		
		fastafile = os.path.abspath(fastafile)
		if not default_params["background"]:
			print "Background file needed!"
			sys.exit()
		bgfile = os.path.abspath(default_params["background"])
		savedir = os.path.abspath(savedir)
		
		tmpdir = tempfile.mkdtemp()

		#savedir = "/tmp/trawler/"

		tmp = tempfile.NamedTemporaryFile()
		shutil.copy(fastafile, tmp.name)
		fastafile = tmp.name
	
		current_path = os.getcwd()
		os.chdir(self.dir())
		
		stdout = ""
		stderr = ""
		strand = "double"
		if default_params["single"]:
			strand = "single"
		cmd = "%s -sample %s -background %s -directory %s -strand %s" % (trawler, fastafile, bgfile, tmpdir, strand)
		p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
		out,err = p.communicate()
		stdout += out
		stderr += err
		
		os.chdir(current_path)
		motifs = []
		out_name = [dir for dir in os.listdir(tmpdir) if dir.startswith("tmp")][-1]
		out_file = os.path.join(tmpdir, out_name, "result", "%s.pwm" % out_name)
		if os.path.exists(out_file):
			motifs = pwmfile_to_motifs(os.path.join(tmpdir, out_name, "result", "%s.pwm" % out_name))
		
		# remove temporary files
		if os.path.exists(tmp.name):
			os.unlink(tmp.name)
		if os.path.exists(tmpdir):
			shutil.rmtree(tmpdir)
		
		return motifs, stdout, stderr

	def parse(self, fo):
		return []


class Weeder(MotifProgram):
	def __init__(self):
		self.name = "Weeder"
		self.cmd = "weederTFBS.out"

	def _run_program(self, bin,fastafile, savedir="", params={}):
		import os, tempfile, shutil
		from subprocess import Popen, PIPE
		
		default_params = {"analysis":"small", "organism":"HS", "single":False, "parallel":True}
		default_params.update(params)
		
		weeder = bin
		adviser = weeder.replace("weederTFBS", "adviser")
	
		
		dir = bin.replace("weederTFBS.out", "")
		if self.is_configured():
			dir = self.dir()

		freq_files = os.path.join(dir, "FreqFiles")
		if not os.path.exists(freq_files):
			raise ValueError, "Can't find FreqFiles directory for Weeder"
				

		fastafile = os.path.abspath(fastafile)
		savedir = os.path.abspath(savedir)

		tmp = tempfile.NamedTemporaryFile()
		name = tmp.name
		tmp.close()
		shutil.copy(fastafile, name)
		fastafile = name
	
		current_path = os.getcwd()
		os.chdir(dir)
		
		coms = ((8,2),(6,1))

		strand = "-S"
		if default_params["single"]:
			strand = ""
			
		if default_params["analysis"] == "xl":
			 coms = ((12,4),(10,3),(8,2),(6,1))
		elif default_params["analysis"] == "large":
			 coms = ((10,3),(8,2),(6,1))
		elif default_params["analysis"] == "medium":
			 coms = ((10,3),(8,2),(6,1))
		
		# TODO: test organism
		organism = default_params["organism"]
		stdout = ""
		stderr = ""
		
		def run_weeder_subset(weeder, fastafile, w, e, organism, strand):
			from subprocess import Popen,PIPE
			cmd = "%s -f %s -W %s -e %s -R 50 -O %s %s" % (weeder, fastafile, w, e, organism, strand)
			p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
			out,err = p.communicate()
			return out, err
			
		if default_params["parallel"]:
			job_server = pp.Server(secret="pumpkinrisotto")
			jobs = []
			for (w,e) in coms:
				jobs.append(job_server.submit(run_weeder_subset, (weeder, fastafile, w, e, organism, strand,), (), ()))
			for job in jobs:
				out,err = job()
				stdout += out
				stderr += err
		else:

			for (w,e) in coms:
				out,err = run_weeder_subset(weeder, fastafile, w, e, organism, strand)
				stdout += out
				stderr += err
	
		cmd = "%s %s" % (adviser, fastafile)
		#print cmd
		p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
		out,err = p.communicate()
		stdout += out
		stderr += err
		
		os.chdir(current_path)

		motifs = []
		if os.path.exists(fastafile + ".wee"):
			f = open(fastafile + ".wee")
			motifs = self.parse(f)
			f.close()
		
		for ext in [".wee", ".html", ".mix", ""]:
			if os.path.exists(fastafile + ext):
				os.unlink(fastafile + ext)

		return motifs, stdout, stderr

	def parse(self, fo):
		motifs = []
		nucs = {"A":0,"C":1,"G":2,"T":3}

		p = re.compile(r'(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)')
		pa = re.compile(r'(\s*\d+\s+.\s+([ACGT]+)\s+.+\))')
		pwm = []
		align = []
		c = 1
		for line in fo.readlines():
			m = p.search(line)
			if m:
				pwm.append([int(m.group(x)) for x in [2,3,4,5]])
			
			m = pa.search(line)
			if m:
				align.append(m.group(2))

			elif line.startswith("===="):
				motifs.append(Motif())
				#total = float(pwm[0][0] + pwm[0][1] + pwm[0][2] + pwm[0][3])
				#motifs[-1].pwm = [[x / total for x in row] for row in pwm]
				motifs[-1].id = "Weeder_%s" % c
				motifs[-1].align= align[:]
				
				width = len(align[0])
				pfm =  [[0 for x in range(4)] for x in range(width)]
				for row in align:
					for i in range(len(row)):
						pfm[i][nucs[row[i]]] += 1
				total = float(len(align)) 
				pwm = [[((x)/total) for x in row] for row in pfm]
				motifs[-1].pwm = pwm[:]
				motifs[-1].pfm = pfm[:]
				
				align = []
				c += 1
				pwm = []
					
		return motifs

class MotifSampler(MotifProgram):
	def __init__(self):
		self.name = "MotifSampler"
		self.cmd = "MotifSampler"

	def _run_program(self, bin, fastafile, savedir,params={}):
		import os, tempfile
		from subprocess import Popen,PIPE
		
		default_params = {"width":10, "background":"", "single":False, "number":10}
		default_params.update(params)
		
		background = default_params['background']
		width = default_params['width']
		number = default_params['number']

		if not background:
			raise Error, "No background specified for %s" % self.name

		fastafile = os.path.abspath(fastafile)
		savedir = os.path.abspath(savedir)
	
		tmp = tempfile.NamedTemporaryFile()
		pwmfile = tmp.name

		tmp2  = tempfile.NamedTemporaryFile()
		outfile = tmp2.name
	
		strand = 1
		if default_params["single"]:
			strand = 0

		# TODO: test organism
		cmd = "%s -f %s -b %s -m %s -w %s -n %s -o %s -s %s > /dev/null 2>&1" % (bin, fastafile, background, pwmfile, width, number, outfile, strand)
		#print cmd
		#p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
		#stdout, stderr = p.communicate()
		
		stdout,stderr = "",""
		p = Popen(cmd, shell=True)
		p.wait()

		motifs = []
		#if os.path.exists(pwmfile):
		#	motifs = self.parse(open(pwmfile))
		if os.path.exists(outfile):
			motifs = self.parse_out(open(outfile))
		
		# remove temporary files
		tmp.close()
		tmp2.close()
		
		return motifs, stdout, stderr

	def parse(self, fo):
		motifs = []

		pwm = []
		info = {}
		for line in fo.readlines():
			if line.startswith("#"):
				vals =  line.strip()[1:].split(" = ")
				if len(vals) > 1:
					info[vals[0]] = vals[1]
			elif len(line) > 1:
				pwm.append([float(x) for x in line.strip().split("\t")])
			else:
				motifs.append(Motif())
				motifs[-1].consensus = info["Consensus"]
				motifs[-1].width = info["W"]
				motifs[-1].id = info["ID"]
				motifs[-1].pwm = pwm[:]
				pwm = []
			
		return motifs


	def parse_out(self, fo):
		motifs = []
		nucs = {"A":0,"C":1,"G":2,"T":3}
		pseudo = 0.0 # Should be 1/sqrt(# of seqs)
		aligns = {}
		for line in fo.readlines():
			if line.startswith("#"):
				pass
			elif len(line) > 1:
				vals = line.strip().split("\t")
				id, site = [x.strip().split(" ")[1].replace('"',"") for x in vals[8].split(";") if x]
				#if vals[6] == "+":
				if site.upper().find("N") == -1:
					aligns.setdefault(id, []).append(site)
				#else:
				#	print site, rc(site)
				#	aligns.setdefault(id, []).append(rc(site))
						
		for id, align in aligns.items():
			#print id, len(align)

			width = len(align[0])
			pfm =  [[0 for x in range(4)] for x in range(width)]
			for row in align:
				for i in range(len(row)):
					pfm[i][nucs[row[i]]] += 1
			total = float(len(align))
			pwm = [[(x + pseudo/4)/total+(pseudo) for x in row] for row in pfm]
			m = Motif()
			m.align = align[:]
			m.pwm = pwm[:]
			m.pfm = pfm[:]
			m.id = id
			motifs.append(m)	
		return motifs

class MDmodule(MotifProgram):
	def __init__(self):
		self.name = "MDmodule"
		self.cmd = "MDmodule"

	def _run_program(self, bin, fastafile, savedir, params={}):
		from subprocess import Popen, PIPE
		import os, tempfile, shutil
		
		default_params = {"width":10, "number":10}
		default_params.update(params)
		
		fastafile = os.path.abspath(fastafile)
		savedir = os.path.abspath(savedir)
	
		tmpdir = tempfile.mkdtemp()
		new_file = os.path.join(tmpdir, "mdmodule_in.fa")
		shutil.copy(fastafile, new_file)
		
		fastafile = new_file
		pwmfile = fastafile + ".out"
	
		width = default_params['width']
		number = default_params['width']
	
		current_path = os.getcwd()
		os.chdir(tmpdir)	
		cmd = "%s -i %s -a 1 -o %s -w %s -t 10 -r %s" % (bin, fastafile, pwmfile, width, number)
		p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
		stdout,stderr = p.communicate()
		
		#stdout,stderr = "",""
		#p = Popen(cmd, shell=True)
		#p.wait()
			
		motifs = []
		if os.path.exists(pwmfile):
			motifs = self.parse(open(pwmfile))
		
		os.chdir(current_path)
		
		# remove temporary files
		shutil.rmtree(tmpdir)
		
		
		return motifs, stdout, stderr

	def parse(self, fo):
		motifs = []
		nucs = {"A":0,"C":1,"G":2,"T":3}
		p = re.compile(r'(\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)')
		pf = re.compile(r'>.+\s+[bf]\d+\s+(\w+)')

		pwm = []
		pfm = []
		align = []
		id = ""
		for line in fo.readlines():
			if line.startswith("Motif"):
				if id:
					motifs.append(Motif())
					motifs[-1].id = id
					motifs[-1].pwm = pwm
					motifs[-1].pfm = pfm
					motifs[-1].align = align
					pwm = []
					pfm = []
					align = []
				id = line.split("\t")[0]
			else: 
				m = p.search(line)
				if m:
					pwm.append([float(m.group(x))/100 for x in [2,3,4,5]])
				m = pf.search(line)
				if m:
					if not pfm:
						pfm = [[0 for x in range(4)] for x in range(len(m.group(1)))]
					for i in range(len(m.group(1))):
						pfm[i][nucs[m.group(1)[i]]] += 1
							
					align.append(m.group(1))
		
		if pwm:
			motifs.append(Motif())
			motifs[-1].id = id
			motifs[-1].pwm = pwm
			motifs[-1].pfm = pfm
			motifs[-1].align = align

		return motifs

class Gadem(MotifProgram):
	def __init__(self):
		self.name = "gadem"
		self.cmd = "gadem"

	def _run_program(self, bin, fastafile, savedir, params={}):
		from subprocess import Popen, PIPE
		import os, tempfile, shutil
		
		default_params = {}
		default_params.update(params)
		
		fastafile = os.path.abspath(fastafile)
		savedir = os.path.abspath(savedir)
	
		tmpdir = tempfile.mkdtemp()
		new_file = os.path.join(tmpdir, "gadem_in.fa")
		shutil.copy(fastafile, new_file)
		
		fastafile = new_file
		pwmfile = fastafile + ".pwm"
		outfile = fastafile + ".out"
	
		current_path = os.getcwd()
		os.chdir(tmpdir)	
		cmd = "%s -fseq %s -fpwm %s -fout %s" % (bin, fastafile, pwmfile, outfile)
		p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE) 
		stdout, stderr = p.communicate()
			
		motifs = []
		if os.path.exists(pwmfile):
			motifs = self.parse(open(pwmfile))
		
		os.chdir(current_path)
		# remove temporary files
		shutil.rmtree(tmpdir)
		
		return motifs, stdout, stderr

	def parse(self, fo):
		motifs = []
		nucs = {"A":0,"C":1,"G":2,"T":3}

		lines = fo.readlines()
		for i in range(0, len(lines), 5):
			align = []
			pwm = []
			pfm = []
			id = ""
			line = lines[i].strip()
			id = line[1:]
			number = id.split("_")[0][1:]
			if os.path.exists("%s.seq" % number):
				for l in open("%s.seq" % number).readlines():
					if not "x" in l and not "n" in l:
						l = l.strip().upper()
						align.append(l)
						if not pfm:
							pfm = [[0 for x in range(4)] for x in range(len(l))]
						for p in range(len(l)):
							pfm[p][nucs[l[p]]] += 1

			m = [l.strip().split(" ")[1].split("\t") for l in lines[i + 1: i + 5]]

			pwm = [[float(m[x][y]) for x in range(4)] for y in range(len(m[0]))]


			motifs.append(Motif(pwm))
			motifs[-1].id = id
			#motifs[-1].pwm = pwm
			if align:
				pass
				motifs[-1].pfm = pfm
				motifs[-1].align = align

		return motifs


class Meme(MotifProgram):
	def __init__(self):
		self.name = "meme"
		self.cmd = "meme.bin"

	def _run_program(self, bin, fastafile, savedir, params={}):
		from subprocess import Popen, PIPE, STDOUT, call
		import os, tempfile, shutil, StringIO
		
		#EVT = 1.0
		default_params = {"width":10, "single":False, "number":10}
		default_params.update(params)
		
		fastafile = os.path.abspath(fastafile)
		savedir = os.path.abspath(savedir)
		tmp = tempfile.NamedTemporaryFile()
		tmpname = tmp.name
	
		strand = "-revcomp"
		if default_params["single"]:
			strand = ""

		width = default_params["width"]
		number = default_params["number"]
		
		cmd = (bin, fastafile, "-text","-dna","-nostatus","-mod", "zoops","-nmotifs", "%s" % number, "-w","%s" % width, "-maxsize", "10000000", "%s" % strand)
		return [], " ".join(cmd), ""
		p = Popen(cmd, bufsize=1, stderr=PIPE, stdout=PIPE) 
		stdout,stderr = p.communicate()

		motifs = []
		motifs = self.parse(StringIO.StringIO(stdout))
		
		# Delete temporary files
		tmp.close()
	 	
		return motifs, stdout, stderr

	def parse(self, fo):
		motifs = []
		nucs = {"A":0,"C":1,"G":2,"T":3}

		p = re.compile('BL   MOTIF (\d+) width=(\d+) seqs=(\d+)')
		pa = re.compile('\)\s+(\w+)')
		line = fo.readline()
		while line:
			m = p.search(line)
			align = []
			pfm = []	
			if m:
				id = "Meme_%s_w%s" % (m.group(1), m.group(2))
				while not line.startswith("//"):
					ma = pa.search(line)
					if ma:
						l = ma.group(1)
						align.append(l)
						if not pfm:
							pfm = [[0 for x in range(4)] for x in range(len(l))]
						for pos in range(len(l)):
							pfm[pos][nucs[l[pos]]] += 1
					
					line = fo.readline()
				
				motifs.append(Motif(pfm[:]))
				motifs[-1].id = id
				motifs[-1].align = align[:]
			line = fo.readline()

		return motifs
