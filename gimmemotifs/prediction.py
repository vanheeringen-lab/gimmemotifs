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
import inspect

# External imports
import pp

# GimmeMotifs imports
from gimmemotifs import tools as tool_classes
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
	
	result = PredictionResult(logger=logger)
	
	# Dynamically load all tools
	toolio = [x[1]() for x in inspect.getmembers(tool_classes, lambda x: inspect.isclass(x) and issubclass(x, tool_classes.MotifProgram)) if x[0] != 'MotifProgram']
	
	# TODO:
	# Add warnings for running time: Weeder, MoAn, GADEM
	for t in toolio:
		if tools.has_key(t.name) and tools[t.name]:
			if t.use_width:
				for i in range(wmin, wmax + 1, step):
	 				logger.debug("Starting %s job, width %s" % (t.name, i))
					job_name = "%s_width_%s" % (t.name, i)
					jobs[job_name] = job_server.submit(t.run, (fastafile, ".",{'width':i, 'background':background, "single":single, "organism":organism},), (tool_classes.MotifProgram,),("gimmemotifs.config",),  result.add_motifs, (job_name,))
			else:
				logger.debug("Starting %s job" % t.name)
				job_name = t.name
				jobs[job_name] = job_server.submit(t.run, (fastafile, ".",{"analysis": analysis, "background":background, "organism":organism}), (tool_classes.MotifProgram,),("gimmemotifs.config",), result.add_motifs, (job_name,))
		else:
			logger.debug("Skipping %s" % t.name)
	
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
