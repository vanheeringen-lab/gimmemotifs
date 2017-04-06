# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.
""" Parallel prediction of sequence motifs """

# Python imports
import os
import sys
import logging
import thread
from time import time,sleep
import inspect

# GimmeMotifs imports
from gimmemotifs import tools as tool_classes
from gimmemotifs.config import MotifConfig, parse_denovo_params
from gimmemotifs.fasta import Fasta
from gimmemotifs import mytmpdir
from gimmemotifs.stats import calc_stats
from multiprocessing import Pool

def mp_calc_stats(motifs, fg_fa, bg_fa, bg_name=None):
    try:
        stats = calc_stats(motifs, fg_fa, bg_fa, ncpus=1)
    except Exception as e:
        sys.stderr.write("ERROR: {}\n".format(e))
        stats = {}

    if not bg_name:
        bg_name = "default"

    return bg_name, stats

def _run_tool(job_name, t, fastafile, params):
    try:
        result = t.run(fastafile, params, mytmpdir())
    except Exception as e:
        result = ([], "", "{} failed to run: {}".format(job_name, e))
    
    return job_name, result

class PredictionResult(object):
    def __init__(self, outfile, logger=None, fg_file=None, background=None, do_counter=True):
        self.lock = thread.allocate_lock()
        self.motifs = []
        self.finished = []
        self.logger = logger
        self.stats = {}
        self.stat_jobs = []
        self.outfile = outfile
        self.job_server = Pool(1)
        self.counter = 0
        self.do_counter = do_counter

        if fg_file and background:
            self.fg_fa = Fasta(fg_file)
            self.background = dict([(bg,Fasta(fname)) for bg,fname in background.items()])
            self.do_stats = True
        else:
            self.do_stats = False

    def add_motifs(self, args):
        self.lock.acquire()
        # Callback function for motif programs
        if args is None or len(args) != 2 or len(args[1]) != 3:
            try:
                job = args[0]
                self.logger.warn("job {} failed".format(job)) 
                self.finished.append(job)
            except Exception:
                self.logger.warn("job failed") 
            return
        
        job, (motifs, stdout, stderr) = args

        if self.logger:
            self.logger.info("%s finished, found %s motifs" % (job, len(motifs))) 
        
        for motif in motifs:
            if self.do_counter:
                self.counter += 1    
                motif.id = "gimme_{}_".format(self.counter) + motif.id
            f = open(self.outfile, "a")
            f.write("%s\n" % motif.to_pfm())
            f.close()
            self.motifs.append(motif)
            
        if self.do_stats:
            #job_id = "%s_%s" % (motif.id, motif.to_consensus())
            if self.logger:
                self.logger.debug("Starting stats job of motif %s" % motif.id)
            for bg_name, bg_fa in self.background.items():
                job = self.job_server.apply_async(
                                    mp_calc_stats, 
                                    (motifs, self.fg_fa, bg_fa, bg_name), 
                                    callback=self.add_stats
                                    )
                self.stat_jobs.append(job)
        if self.logger:
            self.logger.debug("stdout %s: %s" % (job, stdout))
            self.logger.debug("stdout %s: %s" % (job, stderr))
        self.finished.append(job)
        self.lock.release()

    def wait_for_stats(self):
        for job in self.stat_jobs:
            job.get()

    def add_stats(self, args):
        bg_name, stats = args
        if self.logger:
            self.logger.debug("Stats: %s %s" % (bg_name, stats))
        
        for motif_id in stats.keys():
            if not motif_id in self.stats:
                self.stats[motif_id] = {}
        
            self.stats[motif_id][bg_name] = stats[motif_id]

    def submit_remaining_stats(self):
        for motif in self.motifs:
            n = "%s_%s" % (motif.id, motif.to_consensus())
            if n in  self.stats:
                
                self.logger.info("Adding %s again!" % n)
                #job_id = "%s_%s" % (motif.id, motif.to_consensus())
                self.job_server.apply_async(
                                    _calc_motif_stats, 
                                    (motif, self.fg_fa, self.bg_fa), 
                                    callback=self.add_stats)
                

def pp_predict_motifs(fastafile, outfile, analysis="small", organism="hg18", single=False, background="", tools=None, job_server=None, ncpus=8, logger=None, max_time=None, stats_fg=None, stats_bg=None):
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
    
    result = PredictionResult(
                outfile, 
                logger=logger, 
                fg_file=stats_fg, 
                background=stats_bg)
    
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
        if t.name in tools and tools[t.name]:
            if t.use_width:
                for i in range(wmin, wmax + 1, step):
                    logger.debug("Starting %s job, width %s" % (t.name, i))
                    print("Starting %s job, width %s" % (t.name, i))
                    job_name = "%s_width_%s" % (t.name, i)
                    my_params = params.copy()
                    my_params['width'] = i
                    jobs[job_name] = job_server.apply_async(
                        _run_tool,
                        (job_name, t, fastafile, my_params), 
                        callback=result.add_motifs)
            else:
                logger.debug("Starting %s job" % t.name)
                print("Starting %s job" % t.name)
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
            result.submit_remaining_stats()
    ### Or the user gets impatient... ###
    except KeyboardInterrupt:
        # Destroy all running jobs
        logger.info("Caught interrupt, destroying all running jobs")
        job_server.terminate()
        result.submit_remaining_stats()
        
    
    if stats_fg and stats_bg:
        logger.info("waiting for motif statistics")
        n = 0
        last_len = 0 
       
        result.wait_for_stats()
    
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

def predict_motifs(infile, bgfile, outfile, params=None, stats_fg=None, stats_bg=None):
    """ Predict motifs, input is a FASTA-file"""

    # Parse parameters
    required_params = ["tools", "available_tools", "analysis", 
                                "genome", "use_strand", "max_time"]
    if params is None:
        params = parse_denovo_params()
    else:
        for p in required_params:
            if p not in params:
                params = parse_denovo_params()
                break
    
    # Define all tools
    tools = dict(
            [
                (x.strip(), x in [y.strip() for y in params["tools"].split(",")]) 
                    for x in params["available_tools"].split(",")
            ]
            )

    # Predict the motifs
    analysis = params["analysis"]
    #logger.info("starting motif prediction (%s)", analysis)
    #logger.info("tools: %s", 
    #        ", ".join([x for x in tools.keys() if tools[x]]))
    result = pp_predict_motifs(
                    infile, 
                    outfile, 
                    analysis, 
                    params["genome"], 
                    params["use_strand"], 
                    bgfile, 
                    tools, 
                    None, 
                    #logger=logger, 
                    max_time=params["max_time"], 
                    stats_fg=stats_fg, 
                    stats_bg=stats_bg
                )

    motifs = result.motifs
    print "predicted {} motifs".format(len(motifs))
    #logger.info("predicted %s motifs", len(motifs))
    #logger.debug("written to %s", pfm_file)

    if len(motifs) == 0:
        #logger.info("no motifs found")
        sys.exit()
    
    return result

try:
    from gimmemotifs.mp import pool
except:
    pass
