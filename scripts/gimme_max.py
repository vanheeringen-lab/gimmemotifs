#!/usr/bin/env python
from numpy import *
import sys
import os
from gimmemotifs.motif import pwmfile_to_motifs
from optparse import *
import pp

parser = OptionParser()
parser.add_option("-p", "--pwmfile", dest="pwmfile", help="File with pwms", metavar="FILE")
parser.add_option("-s", "--sample", dest="sample", help="Fasta formatted sample file", metavar="FILE")
parser.add_option("-b", "--background", dest="background", help="Fasta formatted background file", metavar="FILE")
parser.add_option("-i", "--ids", dest="ids", help="Comma-seperated list of motif ids to use (default is all ids)", metavar="IDS")

(options, args) = parser.parse_args()

if not options.pwmfile or not options.sample or not options.background:
	parser.print_help()
	exit()

if not os.path.exists(options.sample):
	print "File %s does not exist!" % options.sample
	exit()

if not os.path.exists(options.background):
	print "File %s does not exist!" % options.background
	exit()

def get_scores(motif, file):
	from gimmemotifs.fasta import Fasta
	result = motif.pwm_scan_score(Fasta(file), cutoff=0.0, nreport=1)
	vals = [sorted(x)[-1] for x in result.values()]
	return vals

job_server = pp.Server(secret="pumpkinrisotto")

pwmfile = options.pwmfile
fg_file = options.sample
bg_file = options.background

motifs = dict([(x.id, x) for x in pwmfile_to_motifs(pwmfile)])

ids = []
if options.ids:
  ids = options.ids.split(",")
else:
	ids = motifs.keys()

fg_jobs = {}
bg_jobs = {}

for id in ids:
	if motifs.has_key(id):
		bg_jobs[id] = job_server.submit(get_scores, (motifs[id],bg_file,))
		fg_jobs[id] = job_server.submit(get_scores, (motifs[id],fg_file,))
	else:
		print "Wrong id: %s" % id
		sys.exit()

print "Motif\t# matches\tMax. enrichment\tScore\tCutoff"

for id in ids:
	pos = array(fg_jobs[id]())
	neg = array(bg_jobs[id]())
	factor = len(neg) / float(len(pos))

	scores = array([s for s in hstack((pos, neg)) if sum(neg >= s) > 1])
	enr = array([(sum(pos >= x) / float(sum(neg >= x))) * factor for x in scores])

	#print len(scores), len(enr)
	#for x,y in zip(enr, scores):
	#	print "%s\t%s" % (x,y)

	max_score = scores[enr.argmax()]
	cutoff = (max_score - motifs[id].pwm_min_score()) / (motifs[id].pwm_max_score() - motifs[id].pwm_min_score())

	print "%s\t%s\t%s\t%s\t%s" % (id, sum(pos >= scores[enr.argmax()]), max(enr), scores[enr.argmax()], cutoff)
#print len(pos), len(neg)
