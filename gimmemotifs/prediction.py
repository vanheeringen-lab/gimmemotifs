# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
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

# GimmeMotifs imports
from gimmemotifs import tools as tool_classes
from gimmemotifs.comparison import *
from gimmemotifs.config import *
from gimmemotifs.fasta import *
from gimmemotifs import mytmpdir

try:
    from gimmemotifs.mp import pool
except:
    pass

def _calc_motif_stats(motif, fg_fa, bg_fa):
    try:
        stats = motif.stats(fg_fa, bg_fa)
    except Exception as e:
        sys.stderr.write("ERROR: {}\n".format(e))
        stats = {}
    return motif, stats

def _run_tool(job_name, t, fastafile, params):
    try:
        result = t.run(fastafile, ".", params, mytmpdir())
    except Exception as e:
        result = ([], "", "{} failed to run: {}".format(job_name, e))
    
    return job_name, result

class PredictionResult:
    def __init__(self, outfile, logger=None, fg_file=None, bg_file=None, job_server=None, do_counter=True):
        self.lock = thread.allocate_lock()
        self.motifs = []
        self.finished = []
        self.logger = logger
        self.stats = {}
        self.outfile = outfile
        self.job_server = job_server
        self.counter = 0
        self.do_counter = do_counter

        if fg_file and bg_file:
            self.fg_fa = Fasta(fg_file)
            self.bg_fa = Fasta(bg_file)
            self.do_stats = True
        else:
            self.do_stats = False

    def add_motifs(self, args):
        # Callback function for motif programs
        if args is None or len(args) != 2 or len(args[1]) != 3:
            try:
                job = args[0]
                self.logger.warn("job {} failed".format(job)) 
                self.finished.append(job)
            except:
                self.logger.warn("job failed") 
            return
        
        job, (motifs, stdout, stderr) = args
        
        if self.logger:
            self.logger.info("%s finished, found %s motifs" % (job, len(motifs))) 
        
        for motif in motifs:
            self.lock.acquire()
            if self.do_counter:
                self.counter += 1    
                motif.id = "gimme_{}_".format(self.counter) + motif.id
            f = open(self.outfile, "a")
            f.write("%s\n" % motif.to_pfm())
            f.close()
            self.motifs.append(motif)
            self.lock.release()
            
            if self.do_stats:
                job_id = "%s_%s" % (motif.id, motif.to_consensus())
                if self.logger:
                    self.logger.debug("Starting stats job of motif %s" % motif.id)
                job = self.job_server.apply_async(
                                    _calc_motif_stats, 
                                    (motif, self.fg_fa, self.bg_fa), 
                                    callback=self.add_stats
                                    )
        
        if self.logger:
            self.logger.debug("stdout %s: %s" % (job, stdout))
            self.logger.debug("stdout %s: %s" % (job, stderr))
        self.finished.append(job)

    def add_stats(self, args):
        motif, stats = args
        if self.logger:
            self.logger.debug("Stats: %s %s" % (motif, stats))
        self.stats["{}_{}".format(motif.id, motif.to_consensus())] = stats

    def get_remaining_stats(self):
        for motif in self.motifs:
            n = "%s_%s" % (motif.id, motif.to_consensus())
            if not self.stats.has_key(n):
                
                self.logger.info("Adding %s again!" % n)
                job_id = "%s_%s" % (motif.id, motif.to_consensus())
                job = self.job_server.apply_async(
                                    _calc_motif_stats, 
                                    (motif, self.fg_fa, self.bg_fa), 
                                    callback=self.add_stats)
                

def pp_predict_motifs(fastafile, outfile, analysis="small", organism="hg18", single=False, background="", tools=None, job_server="", ncpus=8, logger=None, max_time=None, fg_file=None, bg_file=None):
    if tools is None:
        tools = {}

    config = MotifConfig()

    if not tools:
        tools = dict([(x,1) for x in config.get_default_params["tools"].split(",")])
    
    if not logger:
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
        job_server = pool
    
    jobs = {}
    
    result = PredictionResult(outfile, logger=logger, fg_file=fg_file, bg_file=bg_file, job_server=job_server)
    
    # Dynamically load all tools
    toolio = [x[1]() for x in inspect.getmembers(
                                                tool_classes, 
                                                lambda x: 
                                                        inspect.isclass(x) and 
                                                        issubclass(x, tool_classes.MotifProgram)
                                                ) if x[0] != 'MotifProgram']
    
    # TODO:
    # Add warnings for running time: Weeder, GADEM
        
    ### Add all jobs to the job_server ###
    params = {
            'analysis': analysis, 
            'background':background, 
            "single":single, 
            "organism":organism
            }

    for t in toolio:
        if tools.has_key(t.name) and tools[t.name]:
            if t.use_width:
                for i in range(wmin, wmax + 1, step):
                    logger.debug("Starting %s job, width %s" % (t.name, i))
                    job_name = "%s_width_%s" % (t.name, i)
                    params['width'] = i
                    jobs[job_name] = job_server.apply_async(
                        _run_tool,
                        (job_name, t, fastafile, params), 
                        callback=result.add_motifs)
            else:
                logger.debug("Starting %s job" % t.name)
                job_name = t.name
                jobs[job_name] = job_server.apply_async(
                        _run_tool,
                        (job_name, t, fastafile, params), 
                        callback=result.add_motifs)
        else:
            logger.debug("Skipping %s" % t.name)
    
    logger.info("all jobs submitted")
    ### Wait until all jobs are finished or the time runs out ###
    start_time = time()    
    try:
        # Run until all jobs are finished
        while len(result.finished) < len(jobs.keys()) and (not(max_time) or time() - start_time < max_time):
            pass
        if len(result.finished) < len(jobs.keys()):
            logger.info("Maximum allowed running time reached, destroying remaining jobs")
            job_server.terminate()
            result.get_remaining_stats()
    ### Or the user gets impatient... ###
    except KeyboardInterrupt, e:
        # Destroy all running jobs
        logger.info("Caught interrupt, destroying all running jobs")
        job_server.terminate()
        result.get_remaining_stats()
        
    logger.info("waiting for motif statistics")
    n = 0
    last_len = 0 
    while len(set(result.stats.keys())) < len(set([str(m) for m in result.motifs])):
        if n >= 30:
            logger.debug("waited long enough")
            logger.debug("motifs: %s, stats: %s", len(result.motifs), len(result.stats.keys()))
            for i,motif in enumerate(result.motifs):
                if "{}_{}".format(motif.id, motif.to_consensus()) not in result.stats:
                    logger.debug("deleting %s", motif)
                    del result.motifs[i]
            break
        sleep(30)
        if len(result.stats.keys()) == last_len:
            n += 1
        else:
            last_len = len(result.stats.keys())
            n = 0
    
    return result
