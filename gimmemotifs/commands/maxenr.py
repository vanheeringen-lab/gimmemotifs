#!/usr/bin/env python
from numpy import array, hstack
import sys
import os
from gimmemotifs.motif import pwmfile_to_motifs
from gimmemotifs.mp import pool

def get_scores(motif, file):
    from gimmemotifs.fasta import Fasta
    result = motif.pwm_scan_score(Fasta(file), cutoff=0.0, nreport=1)
    vals = [sorted(x)[-1] for x in result.values()]
    return vals

def maxenr(args):

    if not os.path.exists(args.sample):
        print "File %s does not exist!" % args.sample
        exit(1)
    
    if not os.path.exists(args.background):
        print "File %s does not exist!" % args.background
        exit(1)
    
    pwmfile = args.pwmfile
    fg_file = args.sample
    bg_file = args.background
    
    motifs = dict([(x.id, x) for x in pwmfile_to_motifs(pwmfile)])
    
    ids = []
    if args.ids:
      ids = args.ids.split(",")
    else:
        ids = motifs.keys()
    
    fg_jobs = {}
    bg_jobs = {}
    
    for id in ids:
        if motifs.has_key(id):
            bg_jobs[id] = pool.apply_async(get_scores, (motifs[id],bg_file,))
            fg_jobs[id] = pool.apply_async(get_scores, (motifs[id],fg_file,))
        else:
            print "Wrong id: %s" % id
            sys.exit()
    
    print "Motif\t# matches\tMax. enrichment\tScore\tCutoff"
    
    for id in ids:
        pos = array(fg_jobs[id].get())
        neg = array(bg_jobs[id].get())
        factor = len(neg) / float(len(pos))
    
        scores = array([s for s in hstack((pos, neg)) if sum(neg >= s) > 1])
        enr = array([(sum(pos >= x) / float(sum(neg >= x))) * factor for x in scores])
    
        #print len(scores), len(enr)
        #for x,y in zip(enr, scores):
        #    print "%s\t%s" % (x,y)
    
        max_score = scores[enr.argmax()]
        cutoff = (max_score - motifs[id].pwm_min_score()) / (motifs[id].pwm_max_score() - motifs[id].pwm_min_score())
    
        print "%s\t%s\t%0.2f\t%0.2f\t%0.3f" % (id, sum(pos >= scores[enr.argmax()]), max(enr), scores[enr.argmax()], cutoff)
    #print len(pos), len(neg)
