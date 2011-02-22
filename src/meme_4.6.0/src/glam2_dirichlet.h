/* structs and functions for Dirichlet mixture and Beta calculations */
#ifndef GLAM2_DIRICHLET_H
#define GLAM2_DIRICHLET_H

#include <stdio.h>

/* struct for storing values of log[ gamma(x + alpha) / gamma(alpha) ] */
typedef struct {
  double alpha;
  int max;  /* maximum value in lookup table */
  double *table;
} lgamma_alpha;

/* Beta distribution (special case of Dirichlet) with growable lookup tables */
typedef struct {
  /* basic data: */
  double alpha;  /* the 1st pseudocount */
  double beta;  /* the 2nd pseudocount */
  /* stored calculations for fast lookup: */
  double sum;  /* alpha + beta */
  lgamma_alpha alpha_lookup;
  lgamma_alpha beta_lookup;
  lgamma_alpha sum_lookup;
} beta;

/* Dirichlet distribution with fixed lookup tables for speed */
typedef struct {
  /* basic data: */
  int dim;  /* dimensionality: number of pseudocounts */
  double *alpha;  /* the pseudocounts */
  /* stored calculations for fast lookup: */
  double sum;  /* sum of pseudocounts */
  int lookup_size;  /* number of entries in each lookup table */
  double *alpha_lookup;  /* 2D array stuffed into 1D array for speed(?) */
  double *sum_lookup;
} dirichlet;

/* Dirichlet mixture */
typedef struct {
  /* basic data: */
  int comp_num;  /* number of components in the mixture */
  double *weights;  /* the weight of each component */
  dirichlet *components;
  /* stored calculations for fast lookup: */
  double *log_weights;  /* probably unnecessary */
  double *scratch;  /* scratch space for doing log_sum calculations */
  double *counts;  /* holds counts converted from int to double for speed */
  int *offsets;  /* holds offsets into alpha_lookup for speed */
} dirichlet_mix;

/* Initialize a Beta distribution */
void beta_init(beta *b, double alpha, double beta, int max_lookup);

/* Calculate log posterior probability of the counts */
/* Allow for lookup tables to be extended */
double beta_score(beta *b, const int count1, const int count2);

/* Calculate 1st ratio with pseudocounts */
double beta_ratio_a(const beta *b, const int a_count, const int b_count);

/* Calculate 2nd ratio with pseudocounts */
double beta_ratio_b(const beta *b, const int a_count, const int b_count);

/* Free the memory buffers in a Beta distribution */
void beta_free(beta *b);

/* Initialize a Dirichlet distribution, without lookup tables */
void dirichlet_init(dirichlet *d, int dim, const double *alpha);

/* Set up a uniform Dirichlet distribution (pseudocounts = 1) */
void dirichlet_uniform(dirichlet *d, int dim);

/* Make lookup tables for the Dirichlet distribution */
void dirichlet_precalc(dirichlet *d, int max_lookup);

/* No functions yet for calculations directly with a dirichlet */

/* Free the memory buffers in a Dirichlet distribution */
void dirichlet_free(dirichlet *d);

/* Initialize a Dirichlet mixture, without lookup tables */
/* all_alpha has the alpha parameters for all components, concatenated */
/* the alphas are passed like this because passing 2D arrays is hard in C */
void dmix_init(dirichlet_mix *m, int comp_num, int dim,
	       const double *weights, const double *all_alpha);

/* Read a Dirichlet mixture from a stream */
void dmix_read(dirichlet_mix *m, FILE *stream);

/* Set up a uniform Dirichlet mixture (1 component, pseudocounts = 1) */
void dmix_uniform(dirichlet_mix *m, int dim);

/* Make lookup tables for each component */
void dmix_precalc(dirichlet_mix *m, int max_lookup);

/* Calculate log posterior probability of the counts */
double dmix_score(dirichlet_mix *m, const int *counts);

/* Calculate "posterior mean estimators" */
void dmix_ratios(dirichlet_mix *m, double *scores, const int *counts);

/* Free the memory buffers in a Dirichlet mixture */
void dmix_free(dirichlet_mix *m);

#endif
