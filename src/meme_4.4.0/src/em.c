/*
 * $Id: em.c 4278 2009-12-23 09:58:37Z james_johnson $
 * 
 * $Log$
 * Revision 1.2  2005/10/25 19:06:39  nadya
 * rm old macro for Header, all info is taken care of by Id and Log.
 *
 * Revision 1.1.1.1  2005/07/29 00:19:08  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/***********************************************************************
*                                                                      *
*       MEME++                                                         *
*       Copyright 1994, The Regents of the University of California    *
*       Author: Timothy L. Bailey                                      *
*                                                                      *
***********************************************************************/
/* 6-22-99 tlb; combine OOPS and ZOOPS e_step */
/* 6-28-99 tlb; add weight on prior on nsites, wnsites */
/* 5-17-99 tlb; fix bug--strcpy cons to model->cons */
/* em.c */

#include "meme.h"
#include "psp.h"
 
static BOOLEAN check_convergence(
  THETA	old_theta,			/* before EM iteration */
  THETA	new_theta,			/* after EM iteration */
  int w,				/* width of motif */
  double distance,			/* convergence radius */
  int alength, 				/* alphabet length */
  int iter,				/* current iteration number */
  int maxiter				/* maximum iterations */
);

/**********************************************************************/
/*
	em

	Uses a version of the EM algorithm (expectation maximization).

*/
/**********************************************************************/
extern void em(
  MODEL *model,			/* the model */
  DATASET *dataset 		/* the dataset */
)
{
  MOTYPE mtype = model->mtype;		/* type of model */
  int max_w = model->w;			/* width of motif */
  int alength = dataset->alength;	/* length of alphabet */
  PRIORS *priors = dataset->priors;	/* the priors */
  double wnsites = dataset->wnsites; 	/* weight on prior on nsites */
  int maxiter = dataset->maxiter; 	/* maximum number of iterations */
  double distance = dataset->distance;	/* stopping criterion */
  THETA theta_save;
  int iter;				/* iteration number */
  double (*E_STEP)(MODEL *, DATASET *); 	/* expectation step */
  double (*E_STEP0)(MODEL *, DATASET *); 	/* expectation step */
  /* maximization step function */
  void (*M_STEP)(MODEL *, DATASET *, PRIORS *, double);
  BOOLEAN converged = FALSE;		/* EM has converged */

  /* create a place to save old value of theta */
  create_2array(theta_save, double, max_w, alength);

  /* set up the correct type of EM to run */
  M_STEP = m_step;
  E_STEP = e_step;
  E_STEP0 = e_step;
  switch (mtype) {
    case Oops:
    case Zoops:
      E_STEP = e_step;
      break;
    case Tcm:
      E_STEP = tcm_e_step;
      break;
    default:
      fprintf(stderr, "Unknown model type in em()! \n");
      exit(1);
      break;
  }
  /* use like_e_step to set z matrix on iteration 0 if motifs were given */
  if (dataset->nkmotifs > 0) {E_STEP0 = E_STEP; E_STEP = like_e_step;}

  /* get the probability that a site starting at position x_ij would
     NOT overlap a previously found motif; used in E_STEP.
  */
  get_not_o(dataset, model->w);

  // renormalize the PSP to the current motif width.
  if (model->mtype != Tcm) {
    psp_renormalize(dataset, model->w, model->invcomp, model->mtype);
  }

  /* Perform EM for number of iterations or until no improvement */
  for (iter=0; iter < maxiter; iter++) {
    int w = model->w;				/* width of motif */
    THETA theta = model->theta;			/* final theta of motif */

    if (iter > 0 && dataset->nkmotifs > 0) E_STEP = E_STEP0;

    if (PRINTALL) { printf("\niter %d\n", iter); }
    if ((!NO_STATUS) && ((iter % 10) == 0)) {
      fprintf(stderr, 
        "\rem: w=%4d, psites=%4.0f, iter=%4d ", w, model->psites, iter);
    }

    /* save current contents of theta */
    copy_theta(theta, theta_save, w, alength);

    /* expectation step */
    model->ll = E_STEP(model, dataset);
    if (PRINT_Z) print_zij(dataset, model);

    /* maximization step */
    M_STEP(model, dataset, priors, wnsites);

    /* print status if requested */
    if (PRINT_LL) {
      double nsites_obs = model->lambda_obs * wps(dataset, w);
      double ll0 = model->mll_0;
      double ll1 = model->mll_1;
      double llr = ll1 - ll0;
      printf("iter=%d w=%d llr=%8.2f nsites_obs=%6.1f\n",
        iter, model->w, llr, nsites_obs);
      printf("w %d ll1 = %f ll0 = %f\n", model->w, ll1, ll0);
    }
    if (PRINTALL) {
      int n = model->nsites_dis;            /* number of sites */
      printf("lambda= %8.6f (theta and obs)\n", model->lambda);
      printf("obs: \n");
      print_theta(0, 1, n, model->obs, model->w, 0, "", dataset, stdout);  
      printf("freq: \n");
      print_theta(0, 1, n, model->theta, model->w, 0, "", dataset, stdout);
    }

    /* see if EM has converged */
    converged = check_convergence(theta_save, theta, w, distance, alength,
      iter, maxiter);

    if (converged) {iter++; break;}		/* done */
  }

  /* save the number of iterations (counting from zero)*/
  model->iter += iter;

  /* discretize, 1 m_step, get relative entropy  */
  (void) discretize(model, dataset);		
  /* use b=0 if using MegaP heuristic */
  if (priors->ptype == MegaP) SWAP(PriorLib *, priors->plib, priors->plib0);
  M_STEP(model, dataset, priors, wnsites);
  if (priors->ptype == MegaP) SWAP(PriorLib *, priors->plib, priors->plib0);
  
  /* get the consensus of the model */
  {
    THETA theta = model->theta;
    int w = model->w;
    char *cons;
    cons = get_consensus(theta, w, dataset, 1, MINCONS); 
    strcpy(model->cons, cons);
    myfree(cons);   
  }

  free_2array(theta_save, max_w);
} /* em */

/**********************************************************************/
/*
	check_convergence
*/
/**********************************************************************/
static BOOLEAN check_convergence(
  THETA	old_theta,			/* before EM iteration */
  THETA	new_theta,			/* after EM iteration */
  int w,				/* width of motif */
  double distance,			/* convergence radius */
  int alength, 				/* alphabet length */
  int iter,				/* current iteration number */
  int maxiter				/* maximum iterations */
)
{
  int i, j;
  double euclid;		/* distance between old_theta and new_theta */
  BOOLEAN converged;

  /* calculate the euclidean change in theta */
  euclid = 0;
  for(i=0; i<w; i++) {
    for(j=0; j<alength; j++) {
      double diff = theta_ref(old_theta, i, j) - theta_ref(new_theta, i, j);
      euclid += diff * diff;
    }
  }
  euclid = sqrt(euclid);
  if (PRINTALL || PRINT_LL) { printf(" d_theta = %f\n", euclid); }

  if (euclid < distance) {		/* converged? */
    if (TRACE) printf("Converged to motif (< %g change) after %d iterations\n",
      distance, iter+1);
    converged = TRUE;
  } else if (maxiter > 1 && iter == maxiter - 1) {
    /* Use fprintf to print from all nodes in parallel. */
    if (TRACE) {
      fprintf(stdout, "Failed to converge after %d iterations!\n", maxiter);
    }
    converged = FALSE;
  } else {
    converged = FALSE;
  }

  return converged;
}
