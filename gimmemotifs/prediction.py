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
# External imports
import pp
# GimmeMotifs imports
from gimmemotifs.tools import *
from gimmemotifs.comparison import *
from gimmemotifs.nmer_predict import *
from gimmemotifs.config import *

def pp_predict_motifs(fastafile, analysis="small", organism="hg18", single=False, background="", tools={}, job_server="", n_cpu=8):
	
	config = MotifConfig()

	if not tools:
		tools = dict([(x,1) for x in config.get_default_params["tools"].split(",")])

	logger = logging.getLogger('prediction.pp_predict_motifs')

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
		job_server = pp.Server(n_cpu, secret='pumpkinrisotto")
	
	result = []
	jobs = []
	
	# Some program specific code
	weeder_organism = {
		"hg18":"HS", 
		"hg19":"HS", 
		"mm9":"MM", 
		"yeast":"SC",
		"sacCer2":"SC",
		"xenTro2":"XT"}[organism]
	
	ms_background = os.path.join(config.get_bg_dir(), "%s.%s.bg" % (organism, "MotifSampler"))
	
	# Start with longer running jobs
	if tools.has_key("MoAn") and tools["MoAn"]:
		moan = MoAn()
		logger.debug("Starting MoAn job")
		jobs.append(job_server.submit(moan.run, (fastafile, ".",{"analysis": analysis, "background":background}), (MotifProgram,),("gimmemotifs.config",)))
	else:
		logger.debug("Skipping MoAn")

	
	if tools.has_key("gadem") and tools["gadem"]:
		gadem = Gadem()
		logger.debug("Starting gadem job")
		jobs.append(job_server.submit(gadem.run, (fastafile, ".",{}), (MotifProgram,),("gimmemotifs.config",)))
	else:
		logger.debug("Skipping gadem")

	
	if tools.has_key("Weeder") and tools["Weeder"]:
		weeder = Weeder()
		logger.debug("Starting Weeder job, analysis %s" % analysis)
		jobs.append(job_server.submit(weeder.run, (fastafile, ".",{"analysis":analysis, "single":single}), (MotifProgram,),("gimmemotifs.config",)))
	else:
		logger.debug("Skipping Weeder")

	if tools.has_key("meme") and tools["meme"]:
		meme = Meme()

		# This check is necessary because meme errors and pp don't play nice together
		p = subprocess.Popen(meme.bin(),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		stdout,stderr = p.communicate()
		if not stderr:
			for i in range(wmin, wmax + 1, step):
				logger.debug("Starting Meme job, width %s" % i)
				jobs.append(job_server.submit(meme.run, (fastafile, ".",{"width":i, "single":single}), (MotifProgram,),("gimmemotifs.config",)))
		else:
			print "Error running meme: %s" % stderr
	else:
		logger.debug("Skipping meme")
	
	
	if  tools.has_key("MotifSampler") and tools["MotifSampler"]:	
		motifsampler = MotifSampler()
		for i in range(wmin, wmax + 1, step):
	 		logger.debug("Starting MotifSampler job, width %s" % (i))
			jobs.append(job_server.submit(motifsampler.run, (fastafile, ".",{'width':i, 'background':ms_background, "single":single},), (MotifProgram,),("gimmemotifs.config",)))
	else:
		logger.debug("Skipping MotifSampler")

	if tools.has_key("trawler") and tools["trawler"]:
		trawler = Trawler()
		logger.debug("Starting trawler job, analysis %s" % analysis)
		jobs.append(job_server.submit(trawler.run, (fastafile, ".",{"analysis":analysis, "background":background, "single":single}), (MotifProgram,),("gimmemotifs.config",)))
	else:
		logger.debug("Skipping trawler")
	
	if tools.has_key("Improbizer") and tools["Improbizer"]:
		improbizer = Improbizer()
		logger.debug("Starting improbizer job")
		jobs.append(job_server.submit(improbizer.run, (fastafile, ".",{"analysis":analysis, "background":background}), (MotifProgram,),("gimmemotifs.config",)))
	else:
		logger.debug("Skipping Improbizer")


	if tools.has_key("MDmodule") and tools["MDmodule"]:	
		mdmodule = MDmodule()
		for i in range(wmin, wmax + 1, step):
			logger.debug("Starting MDmodule job, width %s" % (i))
			jobs.append(job_server.submit(mdmodule.run, (fastafile, ".",{'width':i},), (MotifProgram,),("gimmemotifs.config",)))
	else:
		logger.debug("Skipping MDmodule")
	
	if tools.has_key("BioProspector") and tools["BioProspector"]:	
		prospector = BioProspector()
		for i in range(wmin, wmax + 1, step):
			logger.debug("Starting BioProspector job, width %s" % (i))
			jobs.append(job_server.submit(prospector.run, (fastafile, ".",{'width':i, "background":background, "single":single},), (MotifProgram,),("gimmemotifs.config",)))
	else:
		logger.debug("Skipping BioProspector")

	if tools.has_key("WannaMotif") and tools["WannaMotif"]:
		logger.debug("Starting WannaMotif job")
		jobs.append(job_server.submit(nmer_predict, (fastafile,),(),()))

	motifs = []
	for job in jobs:
		ret = job()
		if ret:
			#print "pp_predict", ret
			(m, stdout, stderr) = ret
			if m and len(m) > 1:
				#print "pp_predict motif ", m

				motifs += m
			logger.debug(stdout)
			logger.debug(stderr)
	
	return motifs

