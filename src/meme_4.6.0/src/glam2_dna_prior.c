#include "glam2_dna_prior.h"

void dmix_dna(dirichlet_mix *m) {
  double a[] = { 0.4, 0.4, 0.4, 0.4 };
  double w[] = { 1 };
  dmix_init(m, 1, 4, w, a);
}
