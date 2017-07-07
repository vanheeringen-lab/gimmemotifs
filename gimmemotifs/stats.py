"""Calculate motif enrichment statistics."""
import sys
from multiprocessing import Pool
import logging

import numpy as np
from scipy.stats import rankdata

from gimmemotifs import rocmetrics
from gimmemotifs.scanner import scan_to_best_match
from gimmemotifs.motif import read_motifs, Motif
from gimmemotifs.config import MotifConfig

logger = logging.getLogger("gimme.stats")

def calc_stats(motifs, fg_file, bg_file, genome=None, stats=None, ncpus=None):
    """Calculate motif enrichment metrics.

    Parameters
    ----------
    motifs : str, list or Motif instance
        A file with motifs in pwm format, a list of Motif instances or a 
        single Motif instance.

    fg_file : str
        Filename of a FASTA, BED or region file with positive sequences.

    bg_file : str
        Filename of a FASTA, BED or region file with negative sequences.

    genome : str, optional
        Genome or index directory in case of BED/regions.
    
    stats : dict, optional
        Names of metrics to calculate. See gimmemotifs.rocmetrics.__all__ 
        for available metrics.

    ncpus : int, optional
        Number of cores to use.

    Returns
    -------
    result : dict
        Dictionary with results where keys are motif ids and the values are
        dictionary with metric name and value pairs.
    """
    if not stats:
        stats = rocmetrics.__all__
    
    if isinstance(motifs, Motif):
        all_motifs = [motifs]
    else:
        try:
            with open(motifs) as f:
                all_motifs = read_motifs(f, fmt="pwm")
        except TypeError:
            all_motifs = motifs
    
    if ncpus is None:
        ncpus = int(MotifConfig().get_default_params()["ncpus"])
    chunksize = 240

    result = {}
    for i in range(0, len(all_motifs), chunksize):
        logger.debug("chunk %s of %s",
            (i / chunksize) + 1, len(all_motifs) / chunksize + 1)
        motifs = all_motifs[i:i + chunksize]
       
        fg_total = scan_to_best_match(fg_file, motifs, ncpus=ncpus, genome=genome)
        bg_total = scan_to_best_match(bg_file, motifs, ncpus=ncpus, genome=genome)
     
        logger.debug("calculating statistics")
        
        if ncpus == 1:
            it = _single_stats(motifs, stats, fg_total, bg_total) 
        else:
            it = _mp_stats(motifs, stats, fg_total, bg_total, ncpus) 
        
        for motif_id, s, ret in it:
            if motif_id not in result:
                result[motif_id] = {}
            result[motif_id][s] = ret
    return result

def _single_stats(motifs, stats, fg_total, bg_total):
    # Initialize multiprocessing pool
    
    for motif in motifs:
        motif_id = motif.id
        fg_vals = fg_total[motif_id]
        bg_vals = bg_total[motif_id]
        for s in stats:
            func = getattr(rocmetrics, s)
            if func.input_type == "score":
                fg = [x[0] for x in fg_vals]
                bg = [x[0] for x in bg_vals]
            elif func.input_type == "pos":
                fg = [x[1] for x in fg_vals]
                bg = [x[1] for x in bg_vals]
            else:
                raise ValueError("Unknown input_type for stats") 
            
            ret = func(fg, bg)
            yield str(motif), s, ret

def _mp_stats(motifs, stats, fg_total, bg_total, ncpus):
    # Initialize multiprocessing pool
    pool = Pool(ncpus, maxtasksperchild=1000)
    
    jobs = []
    for motif in motifs:
        motif_id = motif.id
        fg_vals = fg_total[motif_id]
        bg_vals = bg_total[motif_id]
        for s in stats:
            func = getattr(rocmetrics, s)
            if func.input_type == "score":
                fg = [x[0] for x in fg_vals]
                bg = [x[0] for x in bg_vals]
            elif func.input_type == "pos":
                fg = [x[1] for x in fg_vals]
                bg = [x[1] for x in bg_vals]
            else:
                raise ValueError("Unknown input_type for stats") 
            
            j = pool.apply_async(func, 
                        (fg, bg))
            jobs.append([str(motif), s, j])
    pool.close()
    pool.join()
    
    for motif_id, s, job in jobs:
        ret = job.get() 
        yield motif_id, s, ret

def star(stat, categories):
    stars = 0
    for c in sorted(categories):
        if stat >= c:
            stars += 1
        else:
            return stars
    return stars

def add_star(stats):
    all_stats = {
            "mncp": [2, 5, 8],
            "roc_auc": [0.6, 0.75, 0.9],
            "max_enrichment": [10, 20, 30],
            "enr_at_fpr": [4, 8, 12],
            "fraction_fpr": [0.4, 0.6, 0.8],
            "ks_significance": [4, 7, 10],
            "numcluster": [3, 6, 9],
    }
    
    for motif,s2 in stats.items():
        for bg, s in s2.items():
            stats[motif][bg]["stars"] = int(
                        np.mean([star(s[x], all_stats[x]) for x in all_stats.keys() if x in s]) + 0.5
                    )
    return stats

def rank_motifs(stats, metrics=("roc_auc", "recall_at_fdr")):
    """Determine mean rank of motifs based on metrics."""
    rank = {}
    combined_metrics = []
    motif_ids = stats.keys()
    background = list(stats.values())[0].keys()
    for metric in metrics:
        mean_metric_stats = [np.mean(
            [stats[m][bg][metric] for bg in background]) for m in motif_ids]
        ranked_metric_stats = rankdata(mean_metric_stats)
        combined_metrics.append(ranked_metric_stats)
    
    for motif, val in zip(motif_ids, np.mean(combined_metrics, 0)):
        rank[motif] = val

    return rank

def write_stats(stats, fname, header=None):
    """write motif statistics to text file."""
    # Write stats output to file

    for bg in list(stats.values())[0].keys():
        f = open(fname.format(bg), "w")
        if header:
            f.write(header)
        
        stat_keys = sorted(list(list(stats.values())[0].values())[0].keys())
        f.write("{}\t{}\n".format("Motif", "\t".join(stat_keys)))

        for motif in stats:
            m_stats = stats.get(str(motif), {}).get(bg)
            if m_stats:
                f.write("{}\t{}\n".format(
                    "_".join(motif.split("_")[:-1]),
                    "\t".join([str(m_stats[k]) for k in stat_keys])
                    ))
            else:
                logger.warn("No stats for motif {0}, skipping this motif!".format(motif.id))
            #motifs.remove(motif)
        f.close()

    return

