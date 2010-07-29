/* Initialization done once only at program startup */
#include <assert.h>
#include <stdlib.h>  /* srand */
#include <string.h>  /* strcmp */
#include "glam2_util.h"
#include "glam2_recode3_20comp.h"
#include "glam2_dna_prior.h"
#include "glam2_init.h"
#include "io.h"

void init_alph(alphabet *alph, const char *name) {
  if (strcmp(name, "n") == 0) {  /* nucleotide alphabet */
    alphabet_n(alph);
  } else if (strcmp(name, "p") == 0) {  /* protein alphabet */
    alphabet_p(alph);
  } else {  /* alphabet file */
    FILE *fp = xfopen(name, "r");
    alphabet_read(alph, fp);
    xfclose(fp);
  }
}

void init_seqs(mfasta *seqs, const data *d) {
  int i;
  FILE *fp = xfopen(d->a.seq_file, "r");
  mfasta_read(seqs, fp);
  xfclose(fp);
  if (seqs->seqnum < 1)
    die("%s: no sequences found in %s\n", prog_name, d->a.seq_file);
  for (i = 0; i < seqs->seqnum; ++i)
    first_word(seqs->f[i].title);
  tr_mfasta(seqs, d->alph.encode);
  if (d->a.two_strands)
    rc_mfasta(seqs, d->alph.size);
}

void init_background(prob_vec *bg, const data *d) {
  const int alph_size = d->alph.size;
  const double bg_pseudo = d->a.bg_pseudo * alph_size;
  int *counts = xcalloc(alph_size+1, sizeof(int));  /* zero fill */
  double *probs = xmalloc(alph_size * sizeof(double));
  double *log_probs = xmalloc(alph_size * sizeof(double));
  int tot;
  int i;

  count_mfasta(&d->seqs, counts);

  if (d->a.two_strands)
    for (i = 0; i < (alph_size+1)/2; ++i)
      counts[i] = counts[alph_size-i-1] = counts[i] + counts[alph_size-i-1];

  tot = sum_int(counts, alph_size);

  for (i = 0; i < alph_size; ++i)
    probs[i] = (counts[i] + bg_pseudo * d->alph.prob[i]) / (tot + bg_pseudo);

  for (i = 0; i < alph_size; ++i)
    log_probs[i] = xlog(probs[i]);

  bg->dim = alph_size;
  bg->counts = counts;
  bg->probs = probs;
  bg->log_probs = log_probs;
}

void init_d_prior(beta *b, const data *d) {
  beta_init(b, d->a.delete_pseudo, d->a.no_delete_pseudo, d->seqs.seqnum);
}

void init_i_prior(beta *b, const data *d) {
  beta_init(b, d->a.insert_pseudo, d->a.no_insert_pseudo, d->seqs.seqnum);
}

void init_e_prior(dirichlet_mix *m, const data *d) {
  if (d->a.dirichlet_file) {
    FILE *fp = xfopen(d->a.dirichlet_file, "r");
    dmix_read(m, fp);
    xfclose(fp);
    assert(m->comp_num > 0);
    if (m->components->dim != d->alph.size)
      die("%s: wrong alphabet size in %s\n", prog_name, d->a.dirichlet_file);
  } else if (strcmp(d->a.alph_name, "p") == 0)
    dmix_recode3(m);
  else if (strcmp(d->a.alph_name, "n") == 0)
    dmix_dna(m);
  else
    dmix_uniform(m, d->alph.size);

  dmix_precalc(m, d->seqs.seqnum);
}

void init_convolve(score_matrices *sm, const data *d) {
  const int max_seqlen = d->seqs.maxlen;
  XMALLOC(sm->convolver, max_seqlen + 1);
  if (d->a.algorithm == 2) {
#ifdef FFT
    fft_init(&sm->fft, max_seqlen+1);
#endif
  } else
    XMALLOC(sm->convolved, max_seqlen + 1);    
}

void init_sm(score_matrices *sm, const data *d) {
  const int alph_size = d->alph.size;
  const int max_width = d->a.max_width;
  const int max_seqlen = d->seqs.maxlen;
  int i;

  sm->match_scores = xmalloc2(max_width, (alph_size+1) * sizeof(double));
  XMALLOC(sm->delete_scores, max_width);
  XMALLOC(sm->insert_scores, max_width);

  for (i = 0; i < max_width; ++i)  /* zero score for mask character */
    sm->match_scores[i][alph_size] = 0;

  /* dp_mat needs one extra row and column: */
  sm->dp_mat = xmalloc2(max_width+1, (max_seqlen+1) * sizeof(double));
  /* boundary condition for short-in-long alignment: */
  set_dbl(sm->dp_mat[0], max_seqlen+1, 1);

  if (d->a.two_strands) {
    sm->rc_match_scores = xmalloc2(max_width, (alph_size+1) * sizeof(double));
    XMALLOC(sm->rc_delete_scores, max_width);

    for (i = 0; i < max_width; ++i)  /* zero score for mask character */
      sm->rc_match_scores[i][alph_size] = 0;

    sm->rc_mat = xmalloc2(max_width+1, (max_seqlen+1) * sizeof(double));
    /* boundary condition for short-in-long alignment: */
    set_dbl(sm->rc_mat[0], max_seqlen+1, 1);
  }
}

void init_del_probs(double *del_probs, const int max) {
  int i;
  assert(max > 0);
  del_probs[1] = xlog(2);
  for (i = 1; i < max; ++i)
    del_probs[i+1] = (i+1) * (1 - del_probs[i]) / i;
}

void init_col_sampler(column_sampler *s, const data *d) {
  const int seq_num = d->seqs.seqnum;
  const int max_width = d->a.max_width;
  const int max_seqlen = d->seqs.maxlen;
  s->seq_num = seq_num;
  XMALLOC(s->offsets, seq_num);
  XMALLOC(s->fits, max_width);
  col_init(&s->col, seq_num, d->alph.size);
  XMALLOC(s->scores, max_seqlen + max_width * 2);  /* ? */
  XMALLOC(s->del_probs, max_width + 1);
  init_del_probs(s->del_probs, max_width);
}

void init_seq_order(data *d) {
  const int seq_num = d->seqs.seqnum;
  int i;
  XMALLOC(d->seq_order, seq_num);
  for (i = 0; i < seq_num; ++i)
    d->seq_order[i] = i;
}

void init(data *d) {
  srand(d->a.seed);
  init_alph(&d->alph, d->a.alph_name);
  init_seqs(&d->seqs, d);
#if 0
  if (d->a.min_seqs == 0) 		/* OOPS model requested */
    d->a.min_seqs = d->seqs.seqnum;
#endif
  init_background(&d->scorer.bg, d);
  init_d_prior(&d->scorer.d_prior, d);
  init_i_prior(&d->scorer.i_prior, d);
  init_e_prior(&d->scorer.e_prior, d);
  aln_init(&d->aln, d->seqs.seqnum, d->a.max_width, d->alph.size);
  init_sm(&d->sm, d);
  if (d->a.algorithm > 0)
    init_convolve(&d->sm, d);
  init_col_sampler(&d->col_sampler, d);
  init_seq_order(d);
  // Create output directory.       
  if (create_output_directory(d->a.out_dir, d->a.clobber, !d->a.quiet)) {
    die("Unable to create output directory %s.\n", d->a.out_dir);
  }
  // Open text output file
  d->txt_filename = make_path_to_file(d->a.out_dir, "glam2.txt");
  d->html_filename = make_path_to_file(d->a.out_dir, "glam2.html");
  d->psfm_filename = make_path_to_file(d->a.out_dir, "glam2.meme");
  d->out = xfopen(d->txt_filename, "w");
}
