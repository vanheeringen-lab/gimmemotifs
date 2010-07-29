#ifndef meta_pssm_distr_h
#define meta_pssm_distr_h

extern double *calc_pssm_cdf(
  int w,                // width of PSSM
  int alen,             // length of alphabet
  int range,            // largest value in PSSM
  double **pssm,        // scaled, integer PSSM: pssm[i][j] is score for
                        // j_th letter in i_th column of motif;
                        // entries in PSSM are in range [0..R]
  double *prob          // 0-order Markov background model
);

#endif

