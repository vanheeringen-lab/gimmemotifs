/* Data structures for GLAM2 */
#ifndef GLAM2_GLAM2_H
#define GLAM2_GLAM2_H

#include <stdio.h>
#include "glam2_alphabet.h"
#include "glam2_fasta.h"
#include "glam2_dirichlet.h"
#include "glam2_args.h"
#include "glam2_glam2_aln.h"
#ifdef FFT
#include "glam2_convolve.h"
#endif

typedef struct {
  int dim;
  int *counts;
  double *probs;
  double *log_probs;
} prob_vec;

typedef struct {
  prob_vec bg;  /* background residue probabilities */
  beta d_prior;  /* prior for deletions */
  beta i_prior;  /* prior for insertions */
  dirichlet_mix e_prior;  /* prior for emissions */
} glam2_scorer;

typedef struct {
  double **match_scores;
  double *delete_scores;
  double *insert_scores;
  double **dp_mat;  /* dynamic programming matrix */
  double **rc_mat;  /* dynamic programming matrix for reverse strand */
  double dp_rescale;  /* scale factor to prevent dp_mat overflow */
  double rc_rescale;  /* scale factor to prevent rc_mat overflow */
  double **rc_match_scores;  /* match scores for reverse strand */
  double *rc_delete_scores;  /* delete scores for reverse strand */
  double underflow_flag;  /* lowest value hit during traceback */

  /* Data for nonlinear insertion probabilities: */
  double *convolver;
  double *convolved;
#ifdef FFT
  fft_convolver fft;
#endif
} score_matrices;  /* used for site sampling */

typedef struct {
  int seq_num;
  int width;
  int *offsets;
  int *fits;
  glam2_col col;
  double *scores;
/* probability of choosing a deletion move as a function of alignment width: */
  double *del_probs;  /* for width w, use del_probs[w-1] */
} column_sampler;

typedef struct {
  args a;
  FILE *out;
  char *txt_filename;
  char *html_filename;
  char *psfm_filename;
  alphabet alph;
  mfasta seqs;
  glam2_scorer scorer;
  glam2_aln aln;
  score_matrices sm;
  column_sampler col_sampler;
  int *seq_order;  /* sequence order for starting alignments */
} data;

double insertion_score(glam2_scorer *s, int insert_count, int no_insert_count);
double column_score(glam2_scorer *s, const glam2_col *col);
double marginal_score(glam2_scorer *s, glam2_aln *aln,
		      int seq, const fasta *f);

#endif
