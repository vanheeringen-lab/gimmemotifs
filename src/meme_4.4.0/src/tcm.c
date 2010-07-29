/*
 * $Id: tcm.c 3885 2009-07-16 10:02:35Z tbailey $
 * 
 */

/***********************************************************************
*								       *
*	MEME++							       *
*	Copyright 1994, The Regents of the University of California    *
*	Author: Timothy L. Bailey				       *
*								       *
***********************************************************************/

/* 9-30-99 tlb; remove computation of cum. background probability */
/* 7-02-99 tlb; remove clobbering of theta */
/* 6-23-99 tlb; rewrite to support alternate DNA strands */
/*	
	EM algorithm.

	Two component mixture model. 
*/
	
#include "meme.h"

static double smooth(
  int w,				/* width to smooth over */
  MODEL *model,				/* the model */
  DATASET *dataset			/* the dataset */
);

/**********************************************************************/
/*
	tcm_e_step	

	Do the E step of EM.

	Estimate the expectation of model 1 for each position in the data.
	In other words, calculate E[z_ij] in z.

	Updates z.  

	Returns log pr(X | theta, lambda).

	Time: O(n_samples*lseq*w)
*/
/**********************************************************************/
double tcm_e_step(
  MODEL *model,			/* the model */
  DATASET *dataset  		/* the dataset */
)
{
  int i, j, k, ii;
  THETA logtheta1 = model->logtheta;	/* motif log(theta) */
  THETA logtheta1_rc = model->logtheta_rc;// motif log(theta) reverse complement
  int w = model->w;			/* motif width */
  int n_samples = dataset->n_samples;	/* number of sequences */
  BOOLEAN invcomp = model->invcomp;     /* use reverse complement strand, too */
  double log_sigma = invcomp ? log(0.5) : 0;	/* log \sigma */
  double lambda = model->lambda;	/* \lambda of tcm model */
  double log_lambda = LOG(lambda);	/* log \lambda */
  double log_1mlambda = LOG(1-lambda);	/* log (1 - \lambda) */
  double log_pX;			/* log likelihood; no erase or smooth */
  double log_Pij = 0;			// position-specific prior

  log_Pij = log_sigma;

  /* E step */

  convert_theta_to_log(model, dataset);

  /* calculate all the posterior offset probabilities */
  log_pX = 0;

  for (i=0; i < n_samples; i++) {	/* sequence */
    SAMPLE *s = dataset->samples[i];
    int lseq = s->length;
    double *zi = s->z;			/* Pr(z_ij=1 | X_i, \theta) */
    double *not_o = s->not_o;		/* Pr(V_ij = 1) */
    double *lcb = s->logcumback;	/* cumulative background probability */
    double log_pXij;			// log Pr(X_ij | \phi)
    double log_pXi = 0;			// log Pr(X_i | \phi)

    if (lseq < w) continue;		/* sequence too short for motif */

    int m = lseq - w + 1;		/* number of possible sites */
    for (k=0; k<m; k++) {		// position in sequence
      int j = k + 1;			// Z_ij = 1

      // log ( Pr(X_ij | Z_ij=1, \theta0) \sigma (1-\lambda) )
      double log_pXijtheta0 = log_sigma + log_1mlambda;

      // log ( Pr(X_ij | Z_ij=1, \theta1) \lambda P_ij )
      double log_pXijtheta1 = log_lambda + log_Pij;
      double log_pXijtheta1n = 0;
      // FIXME: this will change if P_i,j != P_i,-j
      if (invcomp) log_pXijtheta1n = log_pXijtheta1;

      /* calculate the probability of positions in the site under the
	background and foreground models
      */
      // background: both strands
      log_pXijtheta0 += Log_back(lcb, k, w);	// Pr(site | \theta_0

      //
      // Z_ij = log Pr(X_ij | Z_ij=1, \theta_1) \sigma \lambda
      //
      char *res = s->res+k;
      if (invcomp) {
	for (ii=0; ii<w; ii++) {
          // foreground: positive strand
          log_pXijtheta1 += logtheta1(ii, (int)res[ii]);
          // foreground: negative strand
          log_pXijtheta1n += logtheta1_rc(ii, (int)res[ii]);
        }
	Zi(j) = log_pXijtheta1;
        Zi(-j) = log_pXijtheta1n;
      } else {
        // foreground: positive strand
	for (ii=0; ii<w; ii++) log_pXijtheta1 += logtheta1(ii, (int)res[ii]);
	Zi(j) = log_pXijtheta1;
      }

      // log_pXij = log Pr(X_i | Z_ij=1, \phi) \sigma \lambda
      log_pXij = LOGL_SUM(log_pXijtheta0, log_pXijtheta1);
      if (invcomp) { 
        double log_pXijn = LOGL_SUM(log_pXijtheta0, log_pXijtheta1n);
        log_pXij = LOGL_SUM(log_pXij, log_pXijn); 
      }

      /* Z_ij : normalize, delog and account for erasing
	Pr(Z_ij=1 | X_i, \phi) \approx
	     P(Z_ij=1, | X_i, \phi) P(V_ij = 1)
      */
      Zi(j) = MIN(1.0, exp(Zi(j) - log_pXij) * not_o[k]);	/* roundoff */
      if (invcomp) Zi(-j) = MIN(1.0, exp(Zi(-j) - log_pXij) * not_o[k]);

      // log_pXi = log Pr(X_i | \phi) \sigma \lamba)
      log_pXi = (k==0) ? log_pXij : LOGL_SUM(log_pXi, log_pXij);

      // log_pX = log Pr(X | \phi) = sum_i,j log(Pr(X_ij)) */
      log_pX += log_pXi;
    } // Z_ij = 1

    // set tail of sequence Zi to 0
    for (j=m+1; j<=lseq; j++) {      // Z_ij = 1
      Zi(j) = 0;			// tail of sequence
      if (invcomp) Zi(-j) = 0;
    } // Z_ij = 1

  } /* sequence */

  /* smooth so no window of size w has z_i which sum to greater than 1.0 */
  (void) smooth(w, model, dataset);

  return (log_pX/log(2.0));
} /* tcm_e_step */

/***********************************************************************/
/*
  smooth

  Normalize so that no local region w wide has z_ij sum of > 1.0.
  Winner-take-all: the largest value of z_ij is never reduced.

  Returns the total expected number of sites of motif.
*/ 
/***********************************************************************/
static double smooth(
  int w,				/* width to smooth over */
  MODEL *model,				/* the model */
  DATASET *dataset			/* the dataset */
)
{
  int i, j, p;
  int n_samples = dataset->n_samples;
  SAMPLE **samples = dataset->samples;
  double p_sum = 0.0;
  BOOLEAN invcomp = model->invcomp;     /* use reverse complement strand, too */

  for (i=0; i<n_samples; i++) {		/* sequence */
    int ioff;
    SAMPLE *sample = samples[i];	/* sample i */
    int lseq = sample->length;
    double *zi= sample->z;		/* z */
    int max_o = lseq - w + 1;		/* largest possible offset */

    if (lseq < w) continue;		/* sequence too short for motif */

    /* normalize adjacent windows of length w, then shift and repeat */
    for (ioff = 0; ioff < MIN(w, max_o); ioff+=2) {	/* window start */
      for (j=ioff; j<max_o; j += w) {		/* adjacent windows */
	double local_z = 0.0;
        double max_z = 0;			/* find largest z_ij */
        int max_p = 0;
        int last_p = MIN(j+w, max_o);
	for (p=j; p<last_p; p++) {		/* position */
          int pp = p + 1;			// pp in [1,...,m]
          double z = (invcomp ? MIN(1.0, Zi(-pp)+Zi(pp)) : Zi(pp));
	  local_z += z;				/* compute local motif z sum */
	  if (z > max_z) {		
	    max_z = z;				/* largest z in window */
	    max_p = p;				/* position with largest z */
	  }
	}
	/* normalize if necessary; leave largest z in window unchanged */
	if (local_z > 1.0) {			/* normalize */
	  double scale = (1 - max_z) / (local_z - max_z);
	  for (p=j; p<last_p; p++) {		/* position */
	    if (p != max_p) {
              int pp = p + 1;			// pp in [1,...,m]
              Zi(pp) *= scale;			/* normalize z */
	      if (invcomp) Zi(-pp) *= scale;/* normalize z on neg. strand */
	    }
	  } /* position */
	} /* normalize */
      } /* adjacent windows */
    } /* window start */

    /* calculate p_sum */
    for (j=0; j < max_o; j++) {
      int pp = p + 1;				// pp in [1,...,m]
      p_sum += Zi(pp);
      if (invcomp) p_sum += Zi(-pp);
    }

  } /* n_samples loop */

  return p_sum;
} /* smooth */
