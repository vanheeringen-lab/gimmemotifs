from gimmemotifs import rocmetrics
from gimmemotifs.scanner import Scanner
from gimmemotifs.motif import read_motifs
from gimmemotifs.config import MotifConfig
from multiprocessing import Pool

def calc_stats(motifs, fg_file, bg_file, stats=None):
    """Calculate motif enrichment metrics.

    Parameters
    ----------
    motifs : str or list
        A file with motifs in pwm format or a list of Motif instances.

    fg_file : str
        Filename of a FASTA file with positive sequences.

    bg_file : str
        Filename of a FASTA file with negative sequences.

    stats : dict, optional
        Names of metrics to calculate. See gimmemotifs.rocmetrics.__all__ 
        for available metrics.

    Returns
    -------
    result : dict
        Dictionary with results where keys are motif ids and the values are
        dictionary with metric name and value pairs.
    """
    if not stats:
        stats = rocmetrics.__all__
    
    s = Scanner()
    s.set_motifs(motifs)
    
    try:
        motifs = read_motifs(open(motifs), fmt="pwm")
    except TypeError:
        pass
    
    ids = [m.id for m in motifs]
    fg_total = dict([(m.id, []) for m in motifs])
    for scores in s.best_score(fg_file):
        for motif,score in zip(motifs, scores):
            fg_total[motif.id].append(score)

    bg_total = dict([(m.id, []) for m in motifs])
    for scores in s.best_score(bg_file):
        for motif,score in zip(motifs, scores):
            bg_total[motif.id].append(score)
 
    jobs = []
    result = {}

    # Initialize multiprocessing pool
    ncpus = int(MotifConfig().get_default_params()["ncpus"])
    pool = Pool(ncpus)
    
    for motif_id in ids:
        result[motif_id] = {}
        fg_vals = fg_total[motif_id]
        bg_vals = bg_total[motif_id]
        for s in stats:
            func = getattr(rocmetrics, s)
            result[motif_id][s] = None
            j = pool.apply_async(func, 
                    (fg_vals, bg_vals))
            jobs.append([motif_id, s, j])
   
    pool.close()
    pool.join()
    
    for motif_id, s, job in jobs:
        ret = job.get() 
        result[motif_id][s] = ret
    
    return result
