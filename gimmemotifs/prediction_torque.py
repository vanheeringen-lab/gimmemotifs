# Copyright (c) 2009-2016 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under 
# the terms of the MIT License, see the file COPYING included with this 
# distribution.

""" Prediction of sequence motifs using PBS/torque """

# Python imports
import sys
import logging
import os
import subprocess 
import thread
from time import time
import inspect
from tempfile import NamedTemporaryFile
import glob

# External imports
import pbs
from PBSQuery import PBSQuery
import yaml

# GimmeMotifs imports
from gimmemotifs import tools as tool_classes
from gimmemotifs.comparison import *
from gimmemotifs.config import *
from gimmemotifs.fasta import *
from gimmemotifs.motif import read_motifs
from gimmemotifs import mytmpdir

class PredictionResult:
    def __init__(self, outfile, logger=None, fg_file=None, bg_file=None):
        self.lock = thread.allocate_lock()
        self.motifs = []
        self.finished = []
        self.logger = logger
        self.stats = {}
        self.outfile = outfile

        if fg_file and bg_file:
            self.fg_fa = Fasta(fg_file)
            self.bg_fa = Fasta(bg_file)
            self.do_stats = True
        else:
            self.do_stats = False

    def add_motifs(self, job, args):
        # Callback function for motif programs
        motifs, stdout, stderr = args
        
        if self.logger:
            self.logger.info("%s finished, found %s motifs" % (job, len(motifs))) 
        
        for motif in motifs:
            self.lock.acquire()
            f = open(self.outfile, "a")
            f.write("%s\n" % motif.to_pfm())
            f.close()
            self.motifs.append(motif)
            if self.do_stats:
                self.logger.debug("Starting stats job of motif %s" % motif.id)
                self.job_server.submit(
                                    motif.stats, 
                                    (self.fg_fa, self.bg_fa), 
                                    (), 
                                    (),
                                    self.add_stats, 
                                    ("%s_%s" % (motif.id, motif.to_consensus()),), 
                                    group="stats"
                                    )
            self.lock.release()
        
        self.logger.debug("stdout %s: %s" % (job, stdout))
        self.logger.debug("stdout %s: %s" % (job, stderr))
        self.finished.append(job)

    def add_stats(self, motif, stats):
        self.logger.debug("Stats: %s %s" % (motif, stats))
        self.stats[motif] = stats

    def get_remaining_stats(self):
        for motif in self.motifs:
            n = "%s_%s" % (motif.id, motif.to_consensus())
            if n not in self.stats:
                self.logger.info("Adding %s again!" % n)
                self.job_server.submit(motif.stats, (self.fg_fa, self.bg_fa), (), (), self.add_stats, ("%s_%s" % (motif.id, motif.to_consensus()),), group="stats")


def write_shell_script(tool, fastafile, rundir=os.path.abspath("."), params=None):  
    
    if params is None:
        params = {}

    local_params = dict([(k,v) for k,v in params.items()])
    
    tmp = NamedTemporaryFile(
                            dir = rundir,
                            prefix = "{0}.{1}.".format(tool, os.getpid()),
                            suffix = ".sh",
                            delete = False,
                            )
    infa = "input.fa"
    bgfa = "bg.fa"
   
    paramfile = tmp.name.replace(".sh", ".params.yaml")
    
    base = os.path.splitext(os.path.basename(tmp.name))[0]
    tmpdir = "$TMPDIR/{0}".format(base)
    
    # Write script
    tmp.write("mkdir {0}\n".format(tmpdir))
    tmp.write("cp {0} {1}/{2}\n".format(os.path.abspath(fastafile), tmpdir, infa))
    if params["background"]:
        tmp.write("cp {0} {1}/{2}\n".format(
                                        os.path.abspath(params["background"]), 
                                        tmpdir,
                                        bgfa,
                                        ))
        local_params["background"] = bgfa

    
    tmp.write("cp {0} {1}/params.yaml\n".format(paramfile, tmpdir))
    tmp.write("cd {0}\n".format(tmpdir))
    
    # If we're running in virtualenv, use the same virtualenv
    activate = os.path.join(os.path.split(sys.executable)[0], "activate")
    if os.path.exists(activate):
        tmp.write("source {0}\n".format(activate))
    
    # The actual motif prediction command
    tmp.write("gimme prediction {0} {1} motifs.pwm -p params.yaml\n".format(
                                                                     tool,
                                                                     infa,
                                                                     ))
    
    tmp.write("mv motifs.pwm {0}/{1}.pwm\n".format(rundir, base))
    tmp.write("sleep 1\n")
    tmp.close()

    # Save parameters to file
    with open(paramfile, "w") as f:
        f.write(yaml.dump(local_params))
    
    return tmp.name

def pp_predict_motifs(fastafile, outfile, analysis="small", organism="hg18", single=False, background="", tools=None, job_server="", ncpus=8, logger=None, max_time=None, fg_file=None, bg_file=None):
    if tools is None:
        tools = {}

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

    jobs = {}
    
    result = PredictionResult(outfile, logger=logger, fg_file=fg_file, bg_file=bg_file)
    
    # Dynamically load all tools
    toolio = [x[1]() for x in inspect.getmembers(
                                                tool_classes, 
                                                lambda x: 
                                                        inspect.isclass(x) and 
                                                        issubclass(x, tool_classes.MotifProgram)
                                                ) if x[0] != 'MotifProgram']
    
    # TODO:
    # Add warnings for running time: Weeder GADEM
        
    # Prepare PBS submission
    server = pbs.pbs_default()
    c = pbs.pbs_connect(server)
    q = PBSQuery()
    attropl = pbs.new_attropl(6)
    # Name
    attropl[0].name  = pbs.ATTR_N
    # Restartable
    attropl[1].name  = pbs.ATTR_r
    attropl[1].value = 'y'
    # Walltime
    attropl[2].name  = pbs.ATTR_l
    attropl[2].resource = 'walltime'
    attropl[2].value = '600'
    # Node requirements
    attropl[3].name  = pbs.ATTR_l
    attropl[3].resource = 'nodes'
    attropl[3].value = '1:ppn=1'   # 
    attropl[4].name  = pbs.ATTR_o
    attropl[5].name  = pbs.ATTR_e
   
    rundir = os.path.join(os.path.split(os.path.abspath(fastafile))[0], "torque")
    if not os.path.exists(rundir):
        os.mkdir(rundir)

    params = {
              'analysis': analysis, 
              'background':background, 
              "single":single, 
              "organism":organism
              }
    
    jobs = {}
    for t in toolio:
        if t.name in tools and tools[t.name]:
            if t.use_width:
                for i in range(wmin, wmax + 1, step):
                    logger.info("Starting %s job, width %s" % (t.name, i))
                    params['width'] = i
                    sh = write_shell_script(t.name, fastafile, rundir=rundir, params=params)
                    job_name = os.path.basename(os.path.splitext(sh)[0]) 
                    # submit
                    attropl[0].value = job_name
                    attropl[4].value = "{0}/{1}.stdout".format(rundir, job_name)
                    attropl[5].value = "{0}/{1}.stderr".format(rundir, job_name)
                    job_id = pbs.pbs_submit(c, attropl, sh, "batchq", 'NULL')
                    e, e_txt = pbs.error()
                    if e:
                        logger.error("Failed: {0}".format(e_txt))
                    else:
                        jobs[job_id] = job_name
            else:
                logger.debug("Starting %s job" % t.name)
                sh = write_shell_script(t.name, fastafile, rundir=rundir, params=params)
                job_name = os.path.basename(os.path.splitext(sh)[0]) 
                # submit
                attropl[0].value = job_name
                attropl[4].value = "{0}/{1}.stdout".format(rundir, job_name)
                attropl[5].value = "{0}/{1}.stderr".format(rundir, job_name)
                job_id = pbs.pbs_submit(c, attropl, sh, "batchq", 'NULL')
                e, e_txt = pbs.error()
                if e:
                    logger.error("Failed submission: {0}".format(e_txt))
                else:
                    jobs[job_id] = job_name
        else:
            logger.debug("Skipping %s" % t.name)
    
    ### Wait until all jobs are finished or the time runs out ###
    start_time = time()  
    try:
        # Run until all jobs are finished
        while len(jobs) > 0 and not(max_time) or time() - start_time < max_time:
            for job_id,job_name in jobs.items():
                job = q.getjob(job_id)
                
                if not job: # or not job.is_running():
                    motifs = []
                    if job:
                        name = job['Job_Name']
                        # Some error checking here!
                    else:
                        pwmfile = os.path.join(rundir, "{0}.pwm".format(job_name))
                        if os.path.exists(pwmfile):
                            motifs = read_motifs(open(pwmfile), fmt="pwm")
                        else:
                            logger.error("Job {0} finished, but couldn find {1}!".format(job_name, pwmfile))
                    stdout = open(os.path.join(rundir, "{0}.stdout".format(job_name))).read()
                    stderr = open(os.path.join(rundir, "{0}.stderr".format(job_name))).read()
                    
                    result.add_motifs(job_id, (motifs, stdout, stderr))
                    #for fname in glob.glob("{0}*".format(job_name)):
                    #    logger.debug("Deleting {0}".format(fname))
                    #    #os.unlink(fname)
                    
                    del jobs[job_id]
            sleep(5)

    ### Or the user gets impatient... ###
    except KeyboardInterrupt as e:
        # Destroy all running jobs
        logger.info("Caught interrupt, destroying all running jobs")
        
    
    logger.info("Waiting for calculation of motif statistics to finish")
    while len(result.stats.keys()) < len(result.motifs):
        sleep(5)
    
    return result
