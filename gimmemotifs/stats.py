from gimmemotifs import rocmetrics
from gimmemotifs.scanner import Scanner
from gimmemotifs.motif import read_motifs
from gimmemotifs.mp import pool

def calc_stats(motifs, fg_file, bg_file, stats=None):
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
   
    for motif_id, s, job in jobs:
        ret = job.get() 
        result[motif_id][s] = ret

    return result

if __name__ == "__main__":
    motifs = "test/data/pwms//motifs.pwm"
    fg = "test/data/pwmscan/promoters.fa"
    bg = "test/data/pwmscan/random_sequences.fa"
    
    #calc_stats(motifs, fg, bg, ["roc_auc"])
    print calc_stats(motifs, fg, bg)

