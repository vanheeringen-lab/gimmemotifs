/*
 * $Id: oops.c 4278 2009-12-23 09:58:37Z james_johnson $
 * 
*/

/***********************************************************************
*								       *
*	MEME++							       *
*	Copyright 1994-2008, The Regents of the University of California    *
*	Author: Timothy L. Bailey				       *
*								       *
***********************************************************************/

/* tlb 3-10-00; remove updating of background model */
/* tlb 7-23-99; add nsites and nsites_obs to model; keep lambda <= 1 */
/* tlb 7-12-99; multiply pXijtheta by not_o_ij in e_step before calc of zij */
/* tlb 7-02-99; remove clobbering of theta by using logtheta[] in model */
/* tlb 6-28-99; add prior to lambda estimation using wnsites */
/* tlb 6-23-99; use same method for background counts for all models */
/* tlb 6-23-99; use average probability of strand/rc strand for background */
/* tlb 6-23-99; remove regularization of background */
/* tlb 6-22-99; precalculate background probability prior to inner loop */
/* tlb 6-22-99; combine zoops.c into this file */
/* tlb 6-22-99; change theta to obs in many places */
/* tlb 6-21-99; changed strand directions to compute 
    szik[j] = Pr(z_ij=1, s_ijk=1 | X_i, \theta)
    zi[j] = Pr(z_ij=1 | X_i, \theta) = sum_k=0^NDIR-1 szik[j]
*/
/* tlb 9-13-94; added different strand directions */

/**********************************************************************/
/*
	EM algorithm
*/
/**********************************************************************/

#include "meme.h"
#include "psp.h"

/**********************************************************************/
/*
  m_step 

	Do the M step of EM.
	Calculate the expected residue frequencies using the
	posterior probabilities of the site starts for each
	sequence.

	Time: O(n_samples*lseq*lseq)
*/
/**********************************************************************/

void m_step(
  MODEL *model,			/* the model */
  DATASET *dataset,		/* the dataset */
  PRIORS *priors,		/* the priors */
  double wnsites 		/* weight on prior on nsites */
)
{
  int i, j, k, ii, jj;
  THETA theta = model->theta;		/* theta of motif */
  THETA obs = model->obs;		/* observed frequencies of motif */
  int w = model->w; 			/* width of motif */
  BOOLEAN invcomp = model->invcomp;	/* use reverse complement strand, too */
  int alength = dataset->alength;	/* length of alphabet */
  int n_samples = dataset->n_samples;	/* number of sequences in dataset */
  SAMPLE **samples = dataset->samples;	/* samples */
  double *back = dataset->back;		/* background frequencies */
  double q=0;				/* $Q$; sum of motif expectations */
  PTYPE ptype = priors->ptype;		/* type of prior */
  PriorLib *plib = priors->plib;	/* library for mixture priors */
  double total, total_obs;		/* scratch */
  double loglike0 = 0;			/* log-likelihood of sites under null */
  double loglike1;			/* log-likelihood of sites undr motif */
  double *p_count = priors->prior_count;/* pseudo counts */

  /* M step */

  /* initialize the observed frequency matrix to zero */
  for (i=0; i<w; i++) {
    for (j=0; j < alength; j++) obs(i, j) = 0; 
  }

  /* calculate the expected frequencies of residues in the different 
     site positions 
  */
  /* get the expected COUNTS, $c_k$ */
  for (i=0; i<n_samples; i++) { 	/* sequence */
    SAMPLE *s = samples[i];		/* array of samples */
    int lseq = s->length;		/* length this seq */
    double sw = s->sw;			/* sequence weight */
    double *lcb = s->logcumback;	/* cumulative backgnd. prob. */
    double qi = 0;			/* sum of z_ij */
    int m = lseq - w + 1;		// number of positions for site
    double *zi = s->z;			// zi[j], j in [-lseq...+lseq]
    if (lseq < w) continue;		/* skip sequence; too short */

    /* get the counts on the two strands (DNA) */
    for (k=0, j=1; k<m; k++,j++) {	// k=seq. pos., j: Z_ij = 1 (Z_i = j)
      char *r;				// for speed
      double z;				// weighted Pr(Z_ij=1 | X_i, \phi)

      //
      // calculate: E[theta | X, z] 
      //

      // positive strand contribution
      int off = k;			// + strand offset from left
      z = Zi(j) * sw;			// weighted Pr(Z_ij=1 | X_i, \phi)
      for (ii=0,r=s->res+off; ii<w; ii++,r++) {		/* position */
	if (*r<alength) {		/* normal letter */
	  obs(ii, (int) *r) += z;	/* motif counts */
	} else { 			/* 'X': spread z over all */
	  for (jj=0; jj<alength; jj++) obs(ii, jj) += z * back[jj];
	}
      } /* sequence position */
      loglike0 += z * Log_back(lcb, k, w);	// Pr(sites | \theta0)
      qi += z;					/* sum of z_ij */

      // negative strand contribution
      if (invcomp) {
        int off = lseq-w-k;		// - strand offset from right
        z = Zi(-j) * sw;		// weighted Pr(Z_ij=1 | X_i, \phi)
	for (ii=0,r=s->resic+off; ii<w; ii++,r++) {	/* position */
	  if (*r<alength) {		/* normal letter */
	    obs(ii, (int) *r) += z;	/* motif counts */
	  } else { 			/* 'X': spread z over all */
	    for (jj=0; jj<alength; jj++) obs(ii, jj) += z * back[jj];
	  }
	} /* sequence position */
        loglike0 += z * Log_back(lcb, k, w);	// Pr(sites | \theta0)
        qi += z;				/* sum of z_ij */
      } // invcomp

    } /* Z_i = j */

    q += qi;					/* $Q$; sum of q_i */
  }
  model->mll_0 = loglike0;			/* site log like. background */

  /* 
     M step for theta: convert COUNTS $c_k$ to FREQUENCIES $f_k$ 
     and use frequencies as estimate for theta--- $p_k^{(t+1)} = f_k$ 
  */
  /* convert counts to frequencies and regularize using pseudo-counts */
  for (i=0; i<w; i++) {				/* width */
    total = total_obs = 0.0;			/* total count for position i */
    /* get the total observed counts */
    for (j=0; j<alength; j++) total_obs += obs(i, j);
    /* get pseudo counts using prior */
    if (ptype == Dmix || ptype == Mega || ptype == MegaP)
      mixture_regularizer(obs[i], plib, p_count);
    /* adjust counts and total count by prior counts */
    for (j=0; j<alength; j++) {
      total += theta(i, j) = obs(i, j) + p_count[j];
    }
    /* normalize counts to probabilities */
    for (j=0; j<alength; j++) {
      obs(i, j) /= (total_obs ? total_obs : 1);	/* normalize to frequencies */
      theta(i, j) /= total;			/* normalize to probability */ 
    }
  } /* w */

  /* palindrome: enforce constraints */
  if (model->pal) {
    palindrome(theta, theta, w, alength);
    palindrome(obs, obs, w, alength);
  }

  /* compute log likelihood of sites under motif model */
  loglike1 = 0;
  for (i=0; i<w; i++) {
    for (j=0; j<alength; j++) {
      double f = obs(i, j);			/* letter frequency */
      if (f) loglike1 += f * LOG(f);		/* motif likelihood */
    } /* letter */
  } /* position in site */
  model->mll_1 = loglike1 *= q;			/* site log like. under motif */

  /* M step for observed lambda */
  model->nsites_obs = q;
  model->lambda_obs = MIN(model->nsites_obs / wps(dataset, w), 1.0);

  /* M step for estimated lambda */
  model->nsites = q*(1-wnsites) + model->psites*wnsites;
  model->lambda = MIN(model->nsites / wps(dataset, w), 1.0);

} /* m_step */

/**********************************************************************/
/*
	e_step

	Do the E step of EM. 

	OOPS and ZOOPS models.

	Updates z array.

	Returns log Pr(X | theta).

	Time: O(n_samples*lseq*w)

	See notes 9/13/94
*/
/**********************************************************************/
double e_step(
  MODEL *model,                 /* the model */
  DATASET *dataset		/* the dataset */
)
{
  int i, j, k, ii;
  MOTYPE mtype = model->mtype;		/* type of model */
  THETA logtheta1 = model->logtheta;	/* motif log(theta) */
  THETA logtheta1_rc = model->logtheta_rc;/* motif log(theta) rev. comp. */
  int w = model->w;			/* width of motif */
  int n_samples = dataset->n_samples;	/* number of sequences */
  BOOLEAN invcomp = model->invcomp;	/* use reverse complement strand, too */
  double gamma = 1.0;			// Pr(Z_i != 0)
  double log_1mgamma;			// log (1-\gamma)
  double logpX;				/* log Pr(X | \phi) */

  /* E step */

  // set up values of sequence model parameters
  if (mtype==Oops) {
    gamma = 1.0;
  } else if (mtype==Zoops) {
    // value of gamma takes into account sequence weighting
    double weighted_m = wps(dataset,w)/n_samples;
    gamma = MIN(1.0,(model->lambda*weighted_m));
  }
  log_1mgamma = LOG(1.0 - gamma);

  convert_theta_to_log(model, dataset);	/* convert theta to log(theta) */

  /* calculate all the posterior offset probabilities */
  logpX = 0.0;				/* prob X given \phi */
  for (i=0; i < n_samples; i++){	/* sequence */
    SAMPLE *s = dataset->samples[i]; 	/* sequence struct */
    int lseq = s->length;		/* length of the sequence */
    int m = lseq - w + 1;		/* number of possible sites */
    double *zi = s->z;			/* log Pr(Z_i=j | X_i, \phi) */ 
    double *not_o = s->not_o;		/* Pr(V_ij = 1) */
    double *lcb = s->logcumback;	/* log cumulative bkg. probability */
    double *log_Pij = s->log_psp;	// log Pr(Z_i=j)
    double log_gamma = LOG(gamma);	// log Pr(Z_i != 0)
    double log_pXijtheta1 = 0;		// log Pr(X_i | Z_i > 0, \theta1)
    double log_pXijtheta1n = 0;		// log Pr(X_i | Z_i < 0, \theta1)
    double log_pXiphi;			// log Pr(X_i | \phi)

    if (lseq < w) continue;		/* skip sequence; too short */

    // initialize log Pr(X_i | \phi) to =Pr(X_i|Z_i=0,phi)PR(Z_i=0|phi)
    // log_pXiphi = (mtype==Oops) ? LITTLE : lcb[lseq] + log_1mgamma;
    log_pXiphi = (mtype==Oops) ? LITTLE : lcb[lseq] + 
      MEME_LOG_SUM(log_1mgamma, log_gamma+log_Pij[0]);

    /* calculate P(X_i | Z_i=j, \phi) */
    for (k=0, j=1; k<m; k++,j++) {	// k=seq. pos., j: Z_ij = 1 (Z_i = j)
      
      // Initialize log Pr(X_i | Z_i=j, \theta1):
      // with erasing, priors and probability outside of site
      double init = lcb[lseq] - Log_back(lcb, k, w) + LOG(not_o[k]) + log_gamma;
      log_pXijtheta1 = init + log_Pij[j];
      // Initialize log Pr(X_i | Z_i=-j, \theta1)
      if (invcomp) log_pXijtheta1n = init + log_Pij[-j];

      //
      // calculate the probability of positions in the site on both strands
      // and use to compute z_ij and Pr(X_i | \phi)
      // z_ij = Pr(X_i | j<0, Z_ij=1 (Z_i=j), \theta1) \gamma P_ij
      //
      char *res = s->res+k;
      // positive strand
      for (ii=0; ii<w; ii++) log_pXijtheta1 += logtheta1(ii, (int) res[ii]);
      Zi(j) = log_pXijtheta1; 
      log_pXiphi = MEME_LOG_SUM(log_pXiphi, Zi(j));
      // negative strand
      if (invcomp) {
	for (ii=0; ii<w; ii++) log_pXijtheta1n += logtheta1_rc(ii,(int)res[ii]);
        Zi(-j) = log_pXijtheta1n;
        log_pXiphi = MEME_LOG_SUM(log_pXiphi, Zi(-j));
      }

    } // position in sequence

    /* calculate logpX = \sum_{i=1}^n Pr(X_i | \phi) */
    logpX += log_pXiphi;

    //  change zi to Pr(Z_i=j | X_i, \phi) (from Pr(X_i | Z_i=j, \phi) )
    //  by normalizing, de-logging and accounting for erasing
    for (k=0, j=1; k<m; k++,j++) {	// k=seq. pos., j: Z_ij = 1 (Z_i = j)
      Zi(j) = exp(Zi(j) - log_pXiphi) * not_o[k];
      if (invcomp) Zi(-j) = exp(Zi(-j) - log_pXiphi) * not_o[k];
    } /* position in Z_i */  

    // Set z_ij for positions motif won't fit to 0
    for (j=m+1; j<=lseq; j++) {		// j: Z_ij = 1 (Z_i = j)
      Zi(j) = 0;
      if (invcomp) Zi(-j) = 0;
    }

  } /* sequence */

  // convert to base-2 logarithms
  return logpX/log(2.0);
} /* e_step */

/**********************************************************************/
/*
	convert_theta_to_log

	Convert theta to log(theta) and include Pr(X).
	Create a reverse-complement version, too.
*/
/**********************************************************************/
extern void convert_theta_to_log(
  MODEL *model,				/* the model */
  DATASET *dataset			/* the dataset */
)
{
  int i, j;
  THETA theta = model->theta;			/* theta */
  THETA logtheta = model->logtheta;		/* log(theta) */
  THETA logtheta1_rc = model->logtheta_rc;	/* log(theta) reverse compl. */
  int w = model->w;				/* width */
  int alength = dataset->alength;		/* length of alphabet */
  double *back = dataset->back;			/* background freqs */

  for (i=0; i<w; i++) {				/* position */
    for (j=0; j<alength; j++) {			/* letter */
      logtheta(i, j) = LOGL(theta(i,j));
      if (model->invcomp) {
        int ic = dnahash(comp_dna(unhash(j)));
        logtheta1_rc(w-i-1, ic) = logtheta(i,j);
      }
    }
    logtheta(i, j) = LOGL(back[j]);		/* Pr(X) */
    if (model->invcomp) {
      int ic = dnahash(comp_dna(unhash(j)));
      logtheta1_rc(w-i-1, ic) = logtheta(i,j);
    }
  }
} /* convert_theta_to_log */
