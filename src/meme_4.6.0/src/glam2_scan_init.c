/* Initialization done once only at program startup */
#include <assert.h>
#include <float.h>  /* DBL_MAX */
#include <string.h>  /* strcmp */
#include "glam2_util.h"
#include "glam2_recode3_20comp.h"
#include "glam2_dna_prior.h"
#include "glam2_scan_init.h"

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

void init_motif(motif *m, const data *d) {
  FILE *fp = xfopen(d->a.motif_file, "r");
  read_motif(m, d->alph.size, d->alph.encode, fp);
  xfclose(fp);
}

void init_d_prior(beta *b, const data *d) {
  beta_init(b, d->a.delete_pseudo, d->a.no_delete_pseudo, d->m.seq_num);
}

void init_i_prior(beta *b, const data *d) {
  beta_init(b, d->a.insert_pseudo, d->a.no_insert_pseudo, d->m.seq_num);
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

  dmix_precalc(m, d->m.seq_num);
}

void init_scores(data *d) {
  const int alph_size = d->m.alph_size;
  const int width = d->m.width;
  const double *bg = d->alph.prob;
  glam2_scorer *s = &d->scorer;
  int i, j;

  d->match_scores = xmalloc2(width, (alph_size+1) * sizeof(double));
  XMALLOC(d->delete_scores, width);
  XMALLOC(d->insert_scores, width);

  for (i = 0; i < width; ++i) {
    const int *residue_counts = d->m.residue_counts[i];
    const int match_count = sum_int(residue_counts, alph_size);
    const int delete_count = d->m.delete_counts[i];
    const int insert_count = d->m.insert_counts[i];
    const int aligned_seq = match_count + delete_count;

    const double ds = beta_ratio_a(&s->d_prior, delete_count, match_count);
    const double ms = beta_ratio_b(&s->d_prior, delete_count, match_count);
    const double is = beta_ratio_a(&s->i_prior, insert_count, aligned_seq);
    const double js = beta_ratio_b(&s->i_prior, insert_count, aligned_seq);

    dmix_ratios(&s->e_prior, d->match_scores[i], residue_counts);

    if (i+1 < width) {
      for (j = 0; j < alph_size; ++j)
        d->match_scores[i][j] = xlog(js * ms * d->match_scores[i][j] / bg[j]);
      d->delete_scores[i] = xlog(ds * js);
      d->insert_scores[i] = xlog(is);
    } else {  /* no insertions after last column: special case */
      for (j = 0; j < alph_size; ++j)
        d->match_scores[i][j] = xlog(ms * d->match_scores[i][j] / bg[j]);
      d->delete_scores[i] = xlog(ds);
      d->insert_scores[i] = -DBL_MAX;  /* overflow risk? */
    }

    /* overflow risk? */
    d->match_scores[i][alph_size] = -DBL_MAX;  /* mask character forbidden */
  }
}

void init_dp(data *d) {
  const int width = d->m.width;
  int i;
  d->dp_seqlen = 0;
  d->dp_mat = xmalloc2(width+1, 1 * sizeof(double));
  d->forbidden = xmalloc2(width+1, 1);
  /* boundary conditions for short-in-long alignment: */
  d->dp_mat[0][0] = 0;
  for (i = 0; i < width; ++i)
    d->dp_mat[i+1][0] = d->delete_scores[i] + d->dp_mat[i][0];
}

void init_hit(data *d) {
  const int width = d->m.width;
  XMALLOC(d->hit_matches, width);
  XMALLOC(d->hit_positions, width);
  XMALLOC(d->hits, d->a.hit_num);
  d->hit_num = 0;
}

void init(data *d) {
  init_alph(&d->alph, d->a.alph_name);
  init_motif(&d->m, d);
  init_d_prior(&d->scorer.d_prior, d);
  init_i_prior(&d->scorer.i_prior, d);
  init_e_prior(&d->scorer.e_prior, d);
  init_scores(d);
  init_dp(d);
  init_hit(d);
  d->out = xfopen(d->a.out_file, "w");
}
