#! /usr/bin/env python

# Read glam2 output, and write it in MEME's PSFM format.  This can be
# used as input to TOMTOM.

import fileinput, string

def do_letter_counts(line):
    words = line.split()[2:-1]
    letters = [i.split('=')[0] for i in words]
    counts = [i.split('=')[1] for i in words]
    counts = map((lambda x: int(x) + 1e-100), counts)	# avoid divide by zero
    tot = sum(counts)
    probs = [i * 1.0 / tot for i in counts]
    probs = ['%.3f' % i for i in probs]
    print
    print 'ALPHABET= ' + string.upper(''.join(letters))
    print
    print 'strands: + -'
    print
    print 'Background letter frequencies (from dataset without add-one prior applied):'
    # avoid "generator expression", for stupid old versions of python:
    print ' '.join([string.upper(i) + ' ' + j for i, j in zip(letters, probs)])
    print

def write_matrix(count_matrix, motif_number, tot_seq_num):
    width = len(count_matrix)
    nsites = sum(count_matrix[0])
    alength = len(count_matrix[0]) - 1
    print 'MOTIF  ' + str(motif_number)
    print
    print 'BL   MOTIF ' + str(motif_number),
    print 'width=' + str(width), 'seqs=' + str(tot_seq_num)
    print 'letter-probability matrix:',
    print 'alength=', str(alength), 'w=', str(width), 'nsites=', str(nsites), 'E= 1'
    for counts in count_matrix:
        counts.pop()  # remove the deletion count
        counts = map((lambda x: x + 1e-100), counts)	# avoid divide by zero
        tot = sum(counts)
        probs = [i * 1.0 / tot for i in counts]
        probs = ['%f' % i for i in probs]
        print ' ' + '  '.join(probs)
    print

matrices = []
state = 0

for line in fileinput.input():
    if state == 0:
        if line.startswith('Version'):
	    print 'MEME version 4.0'
            #print 'GLAM2 version ' + line.split()[1] + '\n'
        elif line.startswith('Sequences'):
            tot_seq_num = line.split()[1]
        elif line.startswith('Residue counts'):
            do_letter_counts(line)
        elif line.rstrip().endswith('Del Ins Score'):
            matrices.append([])
            state = 1
    else:
        words = line.split()
        if len(words) > 2:
            words.pop()  # remove the score
            counts = map(int, words)
            matrices[-1].append(counts)
        elif not words:
            state = 0

for i, m in enumerate(matrices):
    write_matrix(m, i+1, tot_seq_num)
