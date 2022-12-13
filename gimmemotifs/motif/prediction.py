# Copyright (c) 2009-2019 Simon van Heeringen <simon.vanheeringen@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.
"""Parallel prediction of sequence motifs """
import _thread
import inspect
import logging
import warnings
from multiprocessing import Pool
from time import sleep

from gimmemotifs import mytmpdir
from gimmemotifs import tools as tool_classes
from gimmemotifs.config import MotifConfig, parse_denovo_params
from gimmemotifs.fasta import Fasta
from gimmemotifs.stats import calc_stats

logger = logging.getLogger("gimme.prediction")


def mp_calc_stats(motifs, fg_fa, bg_fa, zscore, gc, genome, bg_name=None):
    """Parallel calculation of motif statistics."""
    stats = calc_stats(
        motifs=motifs,
        fg_file=fg_fa,
        bg_file=bg_fa,
        ncpus=1,
        zscore=zscore,
        gc=gc,
        genome=genome,
    )

    if not bg_name:
        bg_name = "default"

    return bg_name, stats


def _run_tool(job_name, t, fastafile, params):
    """Parallel motif prediction."""
    try:
        result = t.run(fastafile, params, mytmpdir())
    except Exception as e:
        result = ([], "", f"{job_name} failed to run: {e}")

    return job_name, result


class PredictionResult(object):
    """Store predicted motifs and calculate statistics."""

    def __init__(
        self,
        outfile,
        genome=None,
        fg_file=None,
        background=None,
        gc=False,
        do_counter=True,
        ncpus=2,
    ):
        self.lock = _thread.allocate_lock()
        self.motifs = []
        self.finished = []
        self.stats = {}
        self.stat_jobs = []
        self.outfile = outfile
        self.genome = genome
        self.pool = Pool(processes=ncpus, maxtasksperchild=1000)
        self.counter = 0
        self.do_counter = do_counter

        open(outfile, "w").close()

        if fg_file and background:
            self.fg_fa = Fasta(fg_file)
            self.background = dict(
                [(bg, Fasta(fname)) for bg, fname in background.items()]
            )
            self.do_stats = True
            self.gc = gc
            self.zscore = self.gc
            if self.gc:
                if genome is None:
                    raise ValueError(
                        "Need a genome when calculating GC% zscores for motif statistics"
                    )
                else:
                    self.genome = genome
        else:
            self.do_stats = False

    def __del__(self):
        # Close the multiprocessing.Pool to release memory
        try:
            self.pool.close()
            self.pool.join()
        except (ValueError, AttributeError):
            pass

    def add_motifs(self, args):
        """Add motifs to the result object."""
        self.lock.acquire()
        # Callback function for motif programs
        if args is None or len(args) != 2 or len(args[1]) != 3:
            try:
                job = args[0]
                logger.warning(f"job {job} failed")
                self.finished.append(job)
            except Exception:
                logger.warning("job failed")
            return

        job, (motifs, stdout, stderr) = args

        logger.info(f"{job} finished, found {len(motifs)} motifs")

        for motif in motifs:
            if self.do_counter:
                self.counter += 1
                motif.id = f"gimme_{self.counter}_{motif.id}"
            f = open(self.outfile, "a")
            f.write(f"{motif.to_pfm()}\n")
            f.close()
            self.motifs.append(motif)

        if self.do_stats and len(motifs) > 0:
            logger.debug(f"Starting stats job of {len(motifs)} motifs")
            for bg_name, bg_fa in self.background.items():
                job = self.pool.apply_async(
                    mp_calc_stats,
                    (
                        motifs,
                        self.fg_fa,
                        bg_fa,
                        self.zscore,
                        self.gc,
                        self.genome,
                        bg_name,
                    ),
                    callback=self.add_stats,
                )
                self.stat_jobs.append(job)

        logger.debug(f"stdout {job}: {stdout}")
        logger.debug(f"stderr {job}: {stderr}")
        self.finished.append(job)
        self.lock.release()

    def wait_for_stats(self):
        """Make sure all jobs are finished."""
        logger.debug("waiting for statistics to finish")
        for job in self.stat_jobs:
            job.get()
        sleep(2)

    def add_stats(self, args):
        """Callback to add motif statistics."""
        bg_name, stats = args
        logger.debug(f"Stats: {bg_name} {stats}")

        for motif_id in stats.keys():
            if motif_id not in self.stats:
                self.stats[motif_id] = {}

            self.stats[motif_id][bg_name] = stats[motif_id]


def pp_predict_motifs(
    fastafile,
    outfile,
    analysis="small",
    organism="hg19",
    single=False,
    background="",
    tools=None,
    ncpus=8,
    max_time=-1,
    stats_fg=None,
    stats_bg=None,
    gc=True,
):
    """Parallel prediction of motifs.

    Utility function for gimmemotifs.denovo.gimme_motifs. Probably better to
    use that, instead of this function directly.
    """
    config = MotifConfig()

    if not tools:
        tools = dict([(x, 1) for x in config.get_default_params()["tools"].split(",")])

    wmin = 5
    step = 1
    if analysis in ["large", "xl"]:
        step = 2
        wmin = 6

    analysis_max = {"xs": 5, "small": 8, "medium": 10, "large": 14, "xl": 20}
    wmax = analysis_max[analysis]

    if analysis == "xs":
        logger.info("Setting analysis xs to small")
        analysis = "small"

    if not ncpus:
        ncpus = config.get_default_params()["ncpus"]
    ncpus = int(ncpus)

    jobs = {}

    result = PredictionResult(
        outfile,
        organism,
        fg_file=stats_fg,
        background=stats_bg,
        gc=gc,
        ncpus=ncpus,
    )

    # Dynamically load all tools
    toolio = [
        x[1]()
        for x in inspect.getmembers(
            tool_classes,
            lambda x: inspect.isclass(x)
            and issubclass(x, tool_classes.motifprogram.MotifProgram),
        )
        if x[0] != "MotifProgram"
    ]

    # TODO: Add warnings for running time: Weeder, GADEM

    # Add all jobs to the pool
    params = {
        "analysis": analysis,
        "background": background,
        "single": single,
        "organism": organism,
    }

    # Tools that don't use a specified width usually take longer
    # ie. GADEM, XXmotif, MEME
    # Start these first.
    for t in [tool for tool in toolio if not tool.use_width]:
        if t.name in tools and tools[t.name]:
            logger.debug(f"Starting {t.name} job")
            job_name = t.name
            jobs[job_name] = result.pool.apply_async(
                _run_tool, (job_name, t, fastafile, params), callback=result.add_motifs
            )
        else:
            logger.debug(f"Skipping {t.name}")

    for t in [tool for tool in toolio if tool.use_width]:
        if t.name in tools and tools[t.name]:
            for i in range(wmin, wmax + 1, step):
                logger.debug(f"Starting {t.name} job, width {i}")
                job_name = f"{t.name}_width_{i}"
                my_params = params.copy()
                my_params["width"] = i
                jobs[job_name] = result.pool.apply_async(
                    _run_tool,
                    (job_name, t, fastafile, my_params),
                    callback=result.add_motifs,
                )
        else:
            logger.debug(f"Skipping {t.name}")

    logger.debug("all jobs submitted")
    for job in jobs.values():
        job.get()

    result.wait_for_stats()

    return result


def predict_motifs(infile, bgfile, outfile, params=None, stats_fg=None, stats_bg=None):
    """Predict motifs, input is a FASTA-file"""

    # Parse parameters
    required_params = [
        "tools",
        "available_tools",
        "analysis",
        "genome",
        "use_strand",
        "max_time",
    ]

    default_params = parse_denovo_params()
    if params is None:
        params = default_params
    else:
        for p in required_params:
            if p not in params and p in default_params:
                if p == "genome":
                    logger.info(f"Using default genome {default_params[p]}")
                params[p] = default_params[p]

    if "genome" not in params:
        logger.error("Need a genome for de novo motif prediction")

    # Define all tools
    tools = dict(
        [
            (x.strip(), x in [y.strip() for y in params["tools"].split(",")])
            for x in params["available_tools"].split(",")
        ]
    )

    # Predict the motifs
    analysis = params["analysis"]
    logger.info(f"starting motif prediction ({analysis})")
    logger.info(f"tools: {', '.join([x for x in tools.keys() if tools[x]])}")
    warnings.filterwarnings("ignore")
    result = pp_predict_motifs(
        infile,
        outfile,
        analysis,
        params.get("genome", None),
        params["use_strand"],
        bgfile,
        tools,
        max_time=params["max_time"],
        stats_fg=stats_fg,
        stats_bg=stats_bg,
    )
    warnings.resetwarnings()

    motifs = result.motifs
    logger.info(f"predicted {len(motifs)} motifs")
    logger.debug(f"written to {outfile}")

    if len(motifs) == 0:
        logger.info("no motifs found")
        result.motifs = []

    return result
