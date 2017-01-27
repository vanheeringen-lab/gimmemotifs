"""Calculate motif enrichment statistics."""
import sys
from gimmemotifs import rocmetrics
from gimmemotifs.scanner import scan_fasta_to_best_score
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
    
    try:
        all_motifs = read_motifs(open(motifs), fmt="pwm")
    except TypeError:
        all_motifs = motifs
    
    ncpus = int(MotifConfig().get_default_params()["ncpus"])
    chunksize = 240

    result = {}
    for i in range(0, len(all_motifs), chunksize):
        sys.stderr.write("{} of {}\n".format(
            (i / chunksize) + 1, len(all_motifs) / chunksize + 1))
        motifs = all_motifs[i:i + chunksize]
        ids = [m.id for m in motifs]
       
        fg_total = scan_fasta_to_best_score(fg_file, motifs)
        bg_total = scan_fasta_to_best_score(bg_file, motifs)
     
        jobs = []
    
        sys.stderr.write("calculating statistics\n")
        # Initialize multiprocessing pool
        pool = Pool(ncpus, maxtasksperchild=1000)
        
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
