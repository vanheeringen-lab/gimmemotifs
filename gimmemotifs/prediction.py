# Copyright (c) 2009-2010 Simon van Heeringen <s.vanheeringen@ncmls.ru.nl>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

""" Parallel prediction of sequence motifs """

# Python imports
import sys
import logging
import os
import subprocess 
import thread
from time import time
# External imports
import pp
# GimmeMotifs imports
from gimmemotifs.tools import *
from gimmemotifs.comparison import *
from gimmemotifs.nmer_predict import *
from gimmemotifs.config import *

class PredictionResult:
	def __init__(self, logger=None):
		self.lock = thread.allocate_lock()
		self.motifs = []
		self.finished = []
		self.logger = logger

	def add_motifs(self, job, args):
		# Callback function for motif programs
		motifs, stdout, stderr = args
		
		if self.logger:
			self.logger.info("%s finished, found %s motifs" % (job, len(motifs))) 
		if len(motifs) > 0:
			self.lock.acquire()
			self.motifs += motifs
			self.lock.release()
		self.logger.debug("stdout %s: %s" % (job, stdout))
		self.logger.debug("stdout %s: %s" % (job, stderr))

		self.finished.append(job)

def pp_predict_motifs(fastafile, analysis="small", organism="hg18", single=False, background="", tools={}, job_server="", n_cpu=8, logger=None, max_time=None):
	
	config = MotifConfig()

	if not tools:
		tools = dict([(x,1) for x in config.get_default_params["tools"].split(",")])
	
	#logger = logging.getLogger('prediction.pp_predict_motifs')

	wmin = 5 
	step = 1
	if analysis in ["large","xl"]:
		step = 2
		wmin = 6
	
	analysis_max = {"xs":5,"small":8, "medium":10,"large":14, "xl":20}
	wmax = analysis_max[analysis]

	if analysis == "xs":
		sys.stderr.write("Setting analysis xs to small")
		analysis = "small"

	if not job_server:
		job_server = pp.Server(n_cpu, secret='pumpkinrisotto')
	
	jobs = {}
	
	# Some program specific code
	weeder_organism = ""
	weeder_organisms = {
		"hg18":"HS", 
		"hg19":"HS", 
		"mm9":"MM", 
		"rn4":"RN",
		"dm3":"DM",
		"fr2": "FR",
		"danRer7": "DR",
		"galGal3": "GG",
		"ce3": "CE",
		"anoGam1": "AG",
		"yeast":"SC",
		"sacCer2":"SC",
		"xenTro2":"XT"}
	if weeder_organisms.has_key(organism):
		weeder_organism = weeder_organisms[organism]
	else:
		logger.info("Weeder not supported for this organism; skipping Weeder")
		del tools["Weeder"]	
	
	ms_background = os.path.join(config.get_bg_dir(), "%s.%s.bg" % (organism, "MotifSampler"))
	
	result = PredictionResult(logger=logger)
	
	# Start with longer running jobs
	if tools.has_key("MoAn") and tools["MoAn"]:
		logger.info("WARNING: MoAn can take a very long time!")
		moan = MoAn()
		job_name = "MoAn"
		logger.debug("Starting MoAn job")
		jobs[job_name] = job_server.submit(moan.run, (fastafile, ".",{"analysis": analysis, "background":background}), (MotifProgram,),("gimmemotifs.config",), result.add_motifs, (job_name,))
	else:
		logger.debug("Skipping MoAn")

	
	if tools.has_key("GADEM") and tools["GADEM"]:
		logger.info("WARNING: GADEM  can take a very long time!")
		gadem = Gadem()
		logger.debug("Starting gadem job")
		job_name = "GADEM"
		jobs[job_name] = job_server.submit(gadem.run, (fastafile, ".",{}), (MotifProgram,),("gimmemotifs.config",), result.add_motifs, (job_name,))
	else:
		logger.debug("Skipping gadem")

	
	if tools.has_key("Weeder") and tools["Weeder"]:
		if analysis == "xl":	
			logger.info("WARNING: Weeder with analysis 'xl' can take a very long time!")
		weeder = Weeder()
		logger.debug("Starting Weeder job, analysis %s" % analysis)
		job_name = "Weeder"
		jobs[job_name] = job_server.submit(weeder.run, (fastafile, ".",{"analysis":analysis, "single":single}), (MotifProgram,),("gimmemotifs.config",), result.add_motifs, (job_name,))
	else:
		logger.debug("Skipping Weeder")

	if tools.has_key("MEME") and tools["MEME"]:
		meme = Meme()

		# This check is necessary because meme errors and pp don't play nice together
		#p = subprocess.Popen(meme.bin(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		#stdout,stderr = p.communicate()
		#if not stderr:
		for i in range(wmin, wmax + 1, step):
			job_name = "MEME_width_%s" % i
			logger.debug("Starting Meme job, width %s" % i)
			jobs[job_name] = job_server.submit(meme.run, (fastafile, ".",{"width":i, "single":single}), (MotifProgram,),("gimmemotifs.config",), result.add_motifs, (job_name,))
		#else:
		#	print "Error running meme: %s" % stderr
	else:
		logger.debug("Skipping MEME")
	
	
	if  tools.has_key("MotifSampler") and tools["MotifSampler"]:	
		motifsampler = MotifSampler()
		for i in range(wmin, wmax + 1, step):
	 		logger.debug("Starting MotifSampler job, width %s" % (i))
			job_name = "MotifSampler_width_%s" % i
			jobs[job_name] = job_server.submit(motifsampler.run, (fastafile, ".",{'width':i, 'background':ms_background, "single":single},), (MotifProgram,),("gimmemotifs.config",),  result.add_motifs, (job_name,))
	else:
		logger.debug("Skipping MotifSampler")

	if tools.has_key("trawler") and tools["trawler"]:
		trawler = Trawler()
		logger.debug("Starting trawler job, analysis %s" % analysis)
		job_name = "trawler"
		jobs[job_name] = job_server.submit(trawler.run, (fastafile, ".",{"analysis":analysis, "background":background, "single":single}), (MotifProgram,),("gimmemotifs.config",), result.add_motifs, (job_name,))
	else:
		logger.debug("Skipping trawler")

	if tools.has_key("Improbizer") and tools["Improbizer"]:
		improbizer = Improbizer()
		logger.debug("Starting improbizer job")
		job_name = "Improbizer"
		jobs[job_name] = job_server.submit(improbizer.run, (fastafile, ".",{"analysis":analysis, "background":background}), (MotifProgram,),("gimmemotifs.config",), result.add_motifs, (job_name,))
	else:
		logger.debug("Skipping Improbizer")

	if tools.has_key("MDmodule") and tools["MDmodule"]:	
		mdmodule = MDmodule()
		for i in range(wmin, wmax + 1, step):
			job_name = "MDmodule_width_%s" % i
			logger.debug("Starting MDmodule job, width %s" % (i))
			jobs[job_name] = job_server.submit(mdmodule.run, (fastafile, ".",{'width':i},), (MotifProgram,),("gimmemotifs.config",), result.add_motifs, (job_name,))
	else:
		logger.debug("Skipping MDmodule")
	
	if tools.has_key("BioProspector") and tools["BioProspector"]:	
		prospector = BioProspector()
		for i in range(wmin, wmax + 1, step):
			job_name = "BioProspector_width_%s" % i
			logger.debug("Starting BioProspector job, width %s" % (i))
			jobs[job_name] = job_server.submit(prospector.run, (fastafile, ".",{'width':i, "background":background, "single":single},), (MotifProgram,),("gimmemotifs.config",), result.add_motifs, (job_name,))
	else:
		logger.debug("Skipping BioProspector")
#
#	if tools.has_key("WannaMotif") and tools["WannaMotif"]:
#		logger.debug("Starting WannaMotif job")
#		jobs.append(job_server.submit(nmer_predict, (fastafile,),(),()))

	start_time = time()	
	try:
		# Run until all jobs are finished
		while len(result.finished) < len(jobs.keys()) and (not(max_time) or time() - start_time < max_time):
			pass
		if len(result.finished) < len(jobs.keys()):
			logger.info("Maximum allowed running time reached, destroying remaining jobs")
			job_server.destroy()
	except KeyboardInterrupt, e:
		# Destroy all running jobs
		logger.info("Caught interrupt, destroying all running jobs")
		job_server.destroy()

	
	return result.motifs
	
