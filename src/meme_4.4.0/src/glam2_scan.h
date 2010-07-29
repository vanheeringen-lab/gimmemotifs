/* Data structures for GLAM2 scanner */
#ifndef GLAM2_SCAN_H
#define GLAM2_SCAN_H

#include <stdio.h>
#include "glam2_alphabet.h"
#include "glam2_fasta.h"
#include "glam2_dirichlet.h"
#include "glam2_motif.h"
#include "glam2_scan_args.h"

typedef struct {
  char *name;
  char strand;
  int coord1;
  int coord2;
  int aln_size;
  int *seq1;
  int *seq2;
  double score;
} alignment;

typedef struct {
  beta d_prior;  /* prior for deletions */
  beta i_prior;  /* prior for insertions */
  dirichlet_mix e_prior;  /* prior for emissions */
} glam2_scorer;  /* merge this with the one in glam2.h */

typedef struct {
  args a;
  FILE *out;
  alphabet alph;
  motif m;
  glam2_scorer scorer;
  double **match_scores;
  double *delete_scores;
  double *insert_scores;
  int dp_seqlen;  /* sequence length that dp_mat is sized for */
  double **dp_mat;
  char **forbidden;  /* forbidden match positions in dp_mat */
  fasta f;  /* the sequence */
  int *hit_positions;
  int *hit_matches;  /* 0=deletion, 1=match */
  alignment *hits;
  int hit_num;  /* how many alignments are stored in hits so far */
} data;

#endif
