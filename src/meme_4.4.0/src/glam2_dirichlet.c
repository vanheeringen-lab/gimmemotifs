#include <assert.h>
#include <errno.h>
#include "glam2_util.h"
#include "glam2_dirichlet.h"

static void
grow_lgamma_alpha(double *table, int old_max, int new_max, double alpha) {
  double tot = table[old_max];
  while (old_max < new_max) {
    tot += xlog(old_max + alpha);
    ++old_max;
    table[old_max] = tot;
  }
}

static void start_lgamma_alpha(double *table, int max, double alpha) {
  assert(alpha > 0);  /* not checked in dirichlet initialization */
  table[0] = 0;
  grow_lgamma_alpha(table, 0, max, alpha);
}

static double get_lgamma_alpha(lgamma_alpha *look, int value) {
  assert(look != NULL);
  assert(value >= 0);
  if (value > look->max) {
    look->table = XREALLOC(look->table, value+1);
    grow_lgamma_alpha(look->table, look->max, value, look->alpha);
    look->max = value;
  }
  return look->table[value];
}

static void init_lgamma_alpha(lgamma_alpha *look, double alpha, int max) {
  assert(look != NULL);
  assert(alpha > 0);
  assert(max >= 0);
  look->alpha = alpha;
  look->max = max;
  XMALLOC(look->table, max + 1);
  start_lgamma_alpha(look->table, max, alpha);
}

static void free_lgamma_alpha(lgamma_alpha *look) {
  assert(look != NULL);
  free(look->table);
}

void beta_init(beta *b, double alpha, double beta, int max_lookup) {
  assert(b != NULL);
  assert(alpha > 0);
  assert(beta > 0);
  assert(max_lookup >= 0);
  b->alpha = alpha;
  b->beta = beta;
  b->sum = b->alpha + b->beta;
  init_lgamma_alpha(&b->alpha_lookup, b->alpha, max_lookup);
  init_lgamma_alpha(&b->beta_lookup, b->beta, max_lookup);
  init_lgamma_alpha(&b->sum_lookup, b->sum, max_lookup);
}

double beta_score(beta *b, const int count1, const int count2) {
  assert(b != NULL);
  assert(count1 >= 0);
  assert(count2 >= 0);
  return get_lgamma_alpha(&b->alpha_lookup, count1)
    + get_lgamma_alpha(&b->beta_lookup, count2)
    - get_lgamma_alpha(&b->sum_lookup, count1 + count2);
}

double beta_ratio_a(const beta *b, const int a_count, const int b_count) {
  assert(b != NULL);
  assert(a_count >= 0);
  assert(b_count >= 0);
  return (a_count + b->alpha) / (a_count + b_count + b->sum);
}

double beta_ratio_b(const beta *b, const int a_count, const int b_count) {
  assert(b != NULL);
  assert(a_count >= 0);
  assert(b_count >= 0);
  return (b_count + b->beta) / (a_count + b_count + b->sum);
}

void beta_free(beta *b) {
  assert(b != NULL);
  free_lgamma_alpha(&b->alpha_lookup);
  free_lgamma_alpha(&b->beta_lookup);
  free_lgamma_alpha(&b->sum_lookup);
}

static void dirichlet_extra_init(dirichlet *d) {
  d->sum = sum_dbl(d->alpha, d->dim);
  d->lookup_size = 0;
  d->alpha_lookup = NULL;
  d->sum_lookup = NULL;
}

void dirichlet_init(dirichlet *d, int dim, const double *alpha) {
  assert(d != NULL);
  assert(dim > 0);
  assert(alpha != NULL);
  d->dim = dim;
  d->alpha = XDUP(alpha, dim);
  dirichlet_extra_init(d);
}

void dirichlet_uniform(dirichlet *d, int dim) {
  assert(d != NULL);
  assert(dim > 0);
  d->dim = dim;
  XMALLOC(d->alpha, dim);
  set_dbl(d->alpha, dim, 1);
  dirichlet_extra_init(d);
}

void dirichlet_precalc(dirichlet *d, int max_lookup) {
  int i;
  assert(d != NULL);
  assert(max_lookup >= 0);

  d->lookup_size = max_lookup + 1;
  XMALLOC(d->alpha_lookup, d->dim * d->lookup_size);
  XMALLOC(d->sum_lookup, d->lookup_size);

  for (i = 0; i < d->dim; ++i)
    start_lgamma_alpha(d->alpha_lookup + i * d->lookup_size,
		       max_lookup, d->alpha[i]);

  start_lgamma_alpha(d->sum_lookup, max_lookup, d->sum);
}

void dirichlet_free(dirichlet *d) {
  assert(d != NULL);
  free(d->alpha);
  free(d->alpha_lookup);
  free(d->sum_lookup);
}

/* Sum elements of lookup table indicated by offsets */
static double dirichlet_lookup_sum(const dirichlet *d, const int *offsets) {
  const int *end = offsets + d->dim;
  const double *a = d->alpha_lookup;
  double sum = 0;

  while (offsets != end) {  /* tried to make this loop fast */
    a += *offsets;
    sum += *a;
    ++offsets;
  }

  return sum;
}

static void dmix_extra_init(dirichlet_mix *m) {
  int i;
  XMALLOC(m->log_weights, m->comp_num);
  XMALLOC(m->scratch, m->comp_num);
  XMALLOC(m->counts, m->components->dim);
  XMALLOC(m->offsets, m->components->dim);
  for (i = 0; i < m->comp_num; ++i)
    m->log_weights[i] = xlog(m->weights[i]);
}

void dmix_init(dirichlet_mix *m, int comp_num, int dim,
	       const double *weights, const double *all_alpha) {
  int i;
  assert(m != NULL);
  assert(comp_num > 0);
  assert(dim > 0);
  assert(weights != NULL);
  assert(all_alpha != NULL);
  m->comp_num = comp_num;
  m->weights = XDUP(weights, comp_num);
  XMALLOC(m->components, comp_num);
  for (i = 0; i < comp_num; ++i)
    dirichlet_init(&m->components[i], dim, all_alpha + i * dim);
  dmix_extra_init(m);
}

void dmix_uniform(dirichlet_mix *m, int dim) {
  assert(m != NULL);
  assert(dim > 0);
  m->comp_num = 1;
  XMALLOC(m->weights, 1);
  m->weights[0] = 1;
  XMALLOC(m->components, 1);
  dirichlet_uniform(m->components, dim);
  dmix_extra_init(m);
}

void dmix_precalc(dirichlet_mix *m, int max_lookup) {
  int i;
  assert(m != NULL);
  assert(max_lookup >= 0);
  for (i = 0; i < m->comp_num; ++i)
    dirichlet_precalc(&m->components[i], max_lookup);
}

void dmix_free(dirichlet_mix *m) {
  int i;
  assert(m != NULL);
  for (i = 0; i < m->comp_num; ++i)
    dirichlet_free(&m->components[i]);
  free(m->weights);
  free(m->components);
  free(m->log_weights);
  free(m->scratch);
  free(m->counts);
  free(m->offsets);
}

/* Get offsets into lookup table from counts */
static void dmix_offsets(dirichlet_mix *m, const int *counts) {
  int extra_offset = 0;
  int i;
  for (i = 0; i < m->components->dim; ++i) {
    m->offsets[i] = extra_offset + counts[i];
    extra_offset = m->components->lookup_size - counts[i];
  }
}

/* if all the counts are 0, the result should be 0, but might be a bit off */
/* tried to make this FAST */
double dmix_score(dirichlet_mix *m, const int *counts) {
  const int tot = sum_int(counts, m->components->dim);
  int i;

  assert(tot < m->components->lookup_size);

  dmix_offsets(m, counts);  /* get alpha_lookup offsets for speed */

  for (i = 0; i < m->comp_num; ++i)
    m->scratch[i] = m->log_weights[i]
      + dirichlet_lookup_sum(&m->components[i], m->offsets)
      - m->components[i].sum_lookup[tot];

  return log_sum(m->scratch, m->comp_num);
}

void dmix_ratios(dirichlet_mix *m, double *scores, const int *counts) {
  const int tot = sum_int(counts, m->components->dim);
/* this call leaves values on m->scratch that we'll use: */
  const double log_norm = dmix_score(m, counts);
  int i, j;

  set_dbl(scores, m->components->dim, 0);

  for (i = 0; i < m->components->dim; ++i)
    m->counts[i] = counts[i];  /* convert ints to doubles here for speed */

  for (i = 0; i < m->comp_num; ++i) {
    const double part_calc = xexp(m->scratch[i] - log_norm) /
      (tot + m->components[i].sum);

    for (j = 0; j < m->components->dim; ++j)
      scores[j] += part_calc * (m->counts[j] + m->components[i].alpha[j]);
  }
}

void dmix_read(dirichlet_mix *m, FILE *stream) {
  char word[21];  /* allow for terminating NUL */

  int mixture_num = 0;
  int alpha_num = 0;
  double *weights = NULL;
  dirichlet *components = NULL;
  double x;

  while (fscanf(stream, "%20s", word) == 1) {
    if (strcmp(word, "Mixture=") == 0) {
      const int r = fscanf(stream, "%lf", &x);
      if (r != 1)
	die("%s: error reading Dirichlet mixture file: Mixture= ???\n", prog_name);
      if (x <= 0)
	die("%s: error reading Dirichlet mixture file: 'Mixture=' values should be > 0\n", prog_name);
      ++mixture_num;
      weights = XREALLOC(weights, mixture_num);
      weights[mixture_num-1] = x;
    } else if (strcmp(word, "Alpha=") == 0) {
      int dim = 0;
      double *alpha = NULL;
      fscanf(stream, "%lf", &x);  /* 1st number is sum of alphas: ignore */
      while (fscanf(stream, "%lf", &x) == 1) {
	if (x <= 0)
	  die("%s: error reading Dirichlet mixture file: 'Alpha=' values should be > 0\n");
	++dim;
	alpha = XREALLOC(alpha, dim);
	alpha[dim-1] = x;
      }
      if (dim == 0)
	die("%s: error reading Dirichlet mixture file: zero pseudocounts\n", prog_name);
      ++alpha_num;
      components = XREALLOC(components, alpha_num);
      dirichlet_init(&components[alpha_num-1], dim, alpha);
      free(alpha);
      if (dim != components->dim)
	die("%s: error reading Dirichlet mixture file: unequal number of values on 'Alpha=' lines\n", prog_name);
    }
  }
  if (ferror(stream))
    die("%s: error reading Dirichlet mixture file: %s\n", prog_name, strerror(errno));

  if (mixture_num != alpha_num)
    die("%s: error reading Dirichlet mixture file: unequal number of 'Mixture=' and 'Alpha=' lines\n", prog_name);

  if (alpha_num == 0)
    die("%s: error reading Dirichlet mixture file: zero components\n", prog_name);

  /* check that weights sum to 1? */

  m->comp_num = alpha_num;
  m->weights = weights;
  m->components = components;
  dmix_extra_init(m);
}
