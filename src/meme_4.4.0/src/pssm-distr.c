/*
 * $Id: pssm-distr.c 4452 2010-04-14 03:45:52Z james_johnson $
 *
 * $Log$
 * Revision 1.2  2005/10/25 21:23:44  nadya
 * add Id and Log lines for keeping track of changes
 * change C++-style comments to C-style
 *
 *
 */

#include "macros.h"
#include "pssm-distr.h"
#include "user.h"
#include <assert.h>

static double *calc_pssm_pdf(
  int w,                /* width of PSSM */
  int alen,             /* length of alphabet */
  int range,            /* largest value in PSSM */
  double **pssm,        /* scaled, integer PSSM: pssm[i][j] is score for
                           j_th letter in i_th column of motif;
                           entries in PSSM are in range [0..range] */
  double *prob,         /* 0-order Markov background model */
  double psfm[MAXSITE][MAXALPH] /* psfm corresponding to the pssm */
);

/**********************************************************************/
/*
        calc_pssm_cdf

        Calculate the 1 minus the theoretical cumulative distribution
        function for an (integer-valued) PSSM given a 0-order Markov
        background model.  PSSM entries must be in range [0..range].

        Returns an array of 1-cdf values:
          pvalue[x] = Pr(score >= x)
        for 0 <= x <= range*w.
*/
/**********************************************************************/
extern double *calc_pssm_cdf(
  int w,                /* width of PSSM */
  int alen,             /* length of alphabet */
  int range,            /* largest value in PSSM */
  double **pssm,        /* scaled, integer PSSM: pssm[i][j] is score for
                           j_th letter in i_th column of motif;
                           entries in PSSM are in range [0..range] */
  double *prob          /* 0-order Markov background model */
)
{
  int i;
  int size = w*range+1;                 /* size of cdf array */
  double *pdf = NULL;                   /* probability distribution */

  pdf = calc_pssm_pdf(w, alen, range, pssm, prob, NULL);
  if (!pdf) return NULL;

  /*for (i=0; i<size; i++) fprintf(stdout, "pdf %d %17.5e\n", i, pdf[i]);*/

  /* compute 1-cdf from the pdf from the right to preserve right accuracy */
  for (i=size-2; i>=0; i--) {
    pdf[i] += pdf[i+1];
    pdf[i] = MIN(1.0, pdf[i]);
    /*
      if (isnan(pdf[i])) {fprintf(stderr, "cdf: i %d pdf %f\n", i, pdf[i]);
        abort();
      }
    */
  }

  /* return the cdf */
  /*for (i=0; i<size; i++) fprintf(stdout, "cdf %d %17.5e\n", i, pdf[i]);*/
  return pdf;
} /* calc_pssm_cdf */

/**********************************************************************/
/*
        calc_pdf

        Calculate the theoretical distribution function for an
        (integer-valued) PSSM.
        Returns an array of pdf values:
          pvalue[x] = Pr(score == x)
        for 0 <= x <= range*w.

        ALTERATION TO INTERFACE, on 17-08-06:
        If "prob" is non-null, then calculate the pdf using that
        background model. In this case, the psfm must be null. Otherwise,
        calculate the pdf under the specified psfm instead (in which case
        the psfm must be non-null).
*/
/**********************************************************************/
static double *calc_pssm_pdf(
  int w,                /* width of PSSM */
  int alen,             /* length of alphabet */
  int range,            /* largest value in PSSM */
  double **pssm,        /* scaled, integer PSSM: pssm[i][j] is score for
                           j_th letter in i_th column of motif;
                           entries in PSSM are in range [0..range] */
  double *prob,         /* 0-order Markov background model */
  double psfm[MAXSITE][MAXALPH]
                        /* psfm corresponding to the pssm */
)
{
  // Check validity of prob/psfm combination:
  assert((prob == NULL) || (psfm == NULL));
  assert(!((prob == NULL) && (psfm == NULL)));

  int i, j, k;
  int size = w*range+1;                 /* size of pdf array */
  double *pdf_old=NULL, *pdf_new=NULL;

  /* set up the two arrays to hold probability density functions */
  Resize(pdf_old, size, double);
  Resize(pdf_new, size, double);

  /* set probabilities of each new score to zero except for score 0 */
  pdf_new[0] = 1;
  for (i=1; i<size; i++) pdf_new[i] = 0;

  /* recursively compute the pdf */
  for (i=0; i < w; i++) {               /* loop over columns in motif */
    int max_score = i * range;          /* maximum possible cumulative score */
    SWAP(double *, pdf_new, pdf_old); /* new column: swap old and new pdfs */
    /* zero out the new pdf; new maximum score is old max + range */
    for (k=0; k<=max_score+range; k++) pdf_new[k] = 0;
    for (j=0; j < alen; j++) {          /* loop over letters */
      int score = (int) pssm[i][j];     /* get integer PSSM entry */
      for (k=0; k<=max_score; k++) {
        /* Calculate the probability under either the background frequency or
           the motif psfm, as appropriate: */
        double curr_prob;
        if (prob != NULL) {
          curr_prob = prob[j];
        }
        else {
          curr_prob = psfm[i][j];
        }
        if (pdf_old[k] != 0) {
          pdf_new[k + score] += pdf_old[k] * curr_prob;
        }

        if (k+score >= size) {
          fprintf(stderr,
            "calc_pssm_pdf error: i=%d j=%d k=%d max_score=%d score=%d size=%d\n",
            i, j, k, max_score, score, size);
          return NULL;
        }
      }
    }
  }

  /* clean up */
  myfree(pdf_old);

  /* return the pdf */
  return pdf_new;

} /* calc_pssm_pdf */
