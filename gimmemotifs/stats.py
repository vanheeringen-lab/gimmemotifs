"""Calculate motif enrichment statistics."""
import sys
from gimmemotifs import rocmetrics
from gimmemotifs.scanner import scan_fasta_to_best_match
from gimmemotifs.motif import read_motifs, Motif
from gimmemotifs.config import MotifConfig
from multiprocessing import Pool

def calc_stats(motifs, fg_file, bg_file, stats=None, ncpus=None):
    """Calculate motif enrichment metrics.

    Parameters
    ----------
    motifs : str, list or Motif instance
        A file with motifs in pwm format, a list of Motif instances or a 
        single Motif instance.

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
    
    if isinstance(motifs, Motif):
        all_motifs = [motifs]
    else:
        try:
            all_motifs = read_motifs(open(motifs), fmt="pwm")
        except TypeError:
            all_motifs = motifs
    
    if ncpus is None:
        ncpus = int(MotifConfig().get_default_params()["ncpus"])
    chunksize = 240

    result = {}
    for i in range(0, len(all_motifs), chunksize):
        sys.stderr.write("{} of {}\n".format(
            (i / chunksize) + 1, len(all_motifs) / chunksize + 1))
        motifs = all_motifs[i:i + chunksize]
       
        fg_total = scan_fasta_to_best_match(fg_file, motifs, ncpus=ncpus)
        bg_total = scan_fasta_to_best_match(bg_file, motifs, ncpus=ncpus)
     
        sys.stderr.write("calculating statistics\n")
        
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
    
    ids = [m.id for m in motifs]
    for motif_id in ids:
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
            yield motif_id, s, ret

def _mp_stats(motifs, stats, fg_total, bg_total, ncpus):
    # Initialize multiprocessing pool
    pool = Pool(ncpus, maxtasksperchild=1000)
    
    jobs = []
    ids = [m.id for m in motifs]
    for motif_id in ids:
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
            jobs.append([motif_id, s, j])
    pool.close()
    pool.join()
    
    for motif_id, s, job in jobs:
        ret = job.get() 
        yield motif_id, s, ret

