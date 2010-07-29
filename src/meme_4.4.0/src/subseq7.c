/*
 * $Id: subseq7.c 4278 2009-12-23 09:58:37Z james_johnson $
 * 
 */

/***********************************************************************
*								       *
*	MEME++							       *
*	Copyright 1994-2006, the Regents of the University of California*
*	Author: Timothy L. Bailey				       *
*								       *
***********************************************************************/
/* subseq7.c */
/* 5-23-00 tlb; change objective function to be log likelihood ratio for
  all models */
/* 7-12-99 tlb; move not_o out of pY calculation and move to local/global_max */
/* 7-12-99 tlb; multiply sw * not_o in m_step of align_top_subsequences so
  that erased starts will have low likelihood */
/* 7-01-99 tlb; multiply background counts by motif width for Tcm model */
/* 6-29-99 tlb; add reverse complement DNA strand support */

#include "calculate_p_y.h"
#include "heap.h"
#include "meme.h"
#include "psp.h"
#include "seed.h"
#include "sp_matrix.h"
#include "macros.h"
#include <assert.h>
#include <math.h> // For the "log" function; temporary.

#define trace(Y)     fprintf(stderr,Y)

/* minimum probability of NOT overlapping a previous site for starting points */
#define MIN_NOT_O .1

/* use logs and integers */
#define LOG_THETA_TYPE(V)		int *V[2][MAXSITE]
#define LOG_THETAG_TYPE(V)		int *V[MAXSITE]

/* print a subsequence from its integer encoded form */
#define Print_res(out, res, w) 						\
{									\
  int i;								\
  for (i=0; i<(w); i++) { fprintf(out, "%c", unhash((res)[i])); }	\
  fprintf(out, "\n");							\
}

/* local functions */
static void next_pY(
  DATASET *dataset,			/* the dataset */
  LOG_THETAG_TYPE(theta_1),		/* integer log theta_1 */
  int w,				/* width of motif */
  int *theta_0,				/* first column of previous theta_1 */
  int pYindex				/* which pY array to use */
);
static void init_theta_1(
  int w,			/* width of site */
  char *res,			/* (encoded) letters of subsequence */
  LOG_THETAG_TYPE(theta_1),	/* theta_1 */
  int lmap[MAXALPH][MAXALPH] 	/* matrix of frequency vectors */ 
);
static int global_max(
  DATASET *dataset,	/* the dataset */
  int w,		/* length of sites */ 
  P_PROB maxima, 	/* array of encoded site starts of local maxima */
  BOOLEAN ic 		/* use reverse complement, too */
);
static int local_max(
  DATASET *dataset,	/* the dataset */
  int w,		/* length of sites */ 
  P_PROB maxima,  	/* array of encoded site starts of local maxima */
  BOOLEAN ic 		/* use reverse complement, too */
);

/**********************************************************************/
/*
	subseq7	

	Try subsequences as starting points and choose the
	one which yields the highest score.
	Score is computed by:
		1) computing log p(Y | theta_1)
		2) finding the (sorted) postions of maximum log pY
		3) aligning the top NSITES0 scores for each value of NSITES0
		4) computing the expected likelihood for each value of NSITES0

	The computing of p(Y | theta_1) for successive
	subsequences (values of theta_1) is optimized by
	using p_i to calculate p_{i+1}.

	Returns number of starting points updated in s_points array.

	Updates s_points, array of starting points, one
	for each value of NSITES0 tried-- finds one \theta for each
	value of nsites0 specified in the input.

        NOTE: (22/11/06)
        The approach used to facilitate dynamic programming was left unchanged.
        HOWEVER, it could be "unified" with the dynamic programming approach
        used in branching_search, in which "SEED_DIFF" objects are used.
        This has not yet been done because "if it ain't broke, don't fix it".
*/
/**********************************************************************/
extern void subseq7(
  MODEL *model,			// the model
  DATASET *dataset,		/* the dataset */
  int w,			// w to use
  int n_nsites0,		/* number of nsites0 values to try */
  S_POINT s_points[],           /* array of starting points: 1 per nsites0 */
  HASH_TABLE evaluated_seed_ht 	/* A hash table used for remembering which seeds
                                   have been evaluated previously */
)
{
  MOTYPE mtype = model->mtype;		/* type of model */
  BOOLEAN ic = model->invcomp;		/* use reverse complement strand of DNA, too */
  THETA map = dataset->map;		/* freq x letter map */
  LOG_THETA_TYPE(ltheta);		/* integer encoded log theta */
  int iseq, ioff;
  int alength = dataset->alength;	/* length of alphabet */
  int n_samples = dataset->n_samples;	/* number of samples in dataset */
  SAMPLE **samples = dataset->samples;	/* samples in dataset */
  int n_starts = 0;			/* number of sampled start subseq */
  int n_maxima = ps(dataset, w);	/* upper bound on # maxima */
  /* the local maxima positions */
  P_PROB maxima = (P_PROB) mymalloc(n_maxima * sizeof(p_prob));
  int lmap[MAXALPH][MAXALPH];	/* consensus letter vs. log frequency matrix */
  double col_scores[MAXSITE];		/* not used */
#ifdef PARALLEL
  int start_seq, start_off=0, end_seq, end_off=0;
#endif
  char *str_seed;                       // A string representation of a seed.

  // PRECONDITIONS:

  // 1. If the sequence model is oops, then n_nsites0 is exactly 1:
  if (mtype == Oops) {
    assert(n_nsites0 == 1);
  }

  convert_to_lmap(map, lmap, alength);

  if (TRACE) { printf("w= %d\n", w); }

  /* get the probability that a site starting at position x_ij would
     NOT overlap a previously found motif.
  */
  get_not_o(dataset, w);

  // Set up log_not_o: log_not_o[site] is:
  // log ( Pr(site not overlapped) * scaled_to_one_Pr(site) )
  if (model->mtype != Tcm) {
    add_psp_to_log_not_o(dataset, w, model->invcomp, model->mtype);
  }

  /* score all the sampled positions saving the best position for
     each value of NSITES0 */
#ifdef PARALLEL
  /* Retrieve the previously-calculated starting and ending points. */
  get_start_n_end(&start_seq, &start_off, &end_seq, &end_off);
  /* Divide the various samples among processors. */
  for (iseq = start_seq; iseq <= end_seq; iseq++) { /* sequence */
#else /* not PARALLEL */
  for (iseq = 0; iseq < n_samples; iseq++) {	/* sequence */
#endif /* PARALLEL */

    SAMPLE *s = samples[iseq];
    int lseq = s->length;
    char *res = s->res;				/* left to right */
    char *name = s->sample_name;
    double *not_o = s->not_o;
    int max_off, init_off;

    if (lseq < w) continue;			/* shorter than motif */

#ifdef PARALLEL
    if (mpMyID() == 0)
#endif
    if ((!NO_STATUS) && ((iseq % 5) == 0)) {
      fprintf(stderr, "starts: w=%d, seq=%d, l=%d          \r", w, iseq, lseq); 
    }
    /* Set the appropriate starting and ending points. */
#ifdef PARALLEL
    if (iseq == start_seq)
      init_off = start_off;
    else
#endif
      init_off = 0;

#ifdef PARALLEL
    if (iseq == end_seq)
      max_off = MIN(end_off, lseq - w);
    else
#endif
      max_off = lseq - w;

    /*
      Loop over all subsequences in the current sequence testing them
      each as "starting points" (inital values) for theta
    */
    for (ioff = init_off; ioff <= max_off; ioff++) {/* subsequence */ 
      /* warning: always do the next step; don't ever
         "continue" or the value of pY will not be correct since
         it is computed based the previous value 
      */

      /* convert subsequence in dataset to starting point for EM */
      init_theta_1(w, res+ioff, &ltheta[1][0], lmap);

      if (ioff == init_off) { 			/* new sequence */

        /* Compute p(Y_ij | theta_1^0) */
        if (!ic) {
          get_pY(dataset, &ltheta[1][0], w, 0);
        } else {
          get_pY(dataset, &ltheta[1][0], w, 1);
          get_pY(dataset, &ltheta[1][0], w, 2);
        }

      } else {					/* same sequence */
        
        /* get theta[0][0]^{k-1} */
        init_theta_1(1, res+ioff-1, &ltheta[0][0], lmap);

        /* compute p(Y_ij | theta_1^k) */
        if (!ic) {
          next_pY(dataset, &ltheta[1][0], w, &ltheta[0][0][0], 0);
        } else {
          next_pY(dataset, &ltheta[1][0], w, &ltheta[0][0][0], 1);
          next_pY(dataset, &ltheta[1][0], w, &ltheta[0][0][0], 2);
        }
      } /* same sequence */

      /* skip if there is a high probability that this subsequence
         is part of a site which has already been found 
      */
      if (not_o[ioff] < MIN_NOT_O) continue;

      /*fprintf(stderr, "subseq: %d %d\r", iseq+1, ioff+1);*/

      // Put highest pY into first scratch array if using both DNA strands:
      if (ic) {
        combine_strands(samples, n_samples, w);
      }

      /* get a sorted list of the maxima of pY */
      n_maxima = get_max(mtype, dataset, w, maxima, ic, TRUE);

      /* "fake out" align_top_subsequences by setting each of the scores in
         the s_points objects to LITTLE, thereby forcing
         align_top_subsequences to record the attributes for the current seed
         in the s_points, rather than the seed with the highest respective
         scores: */
      int sp_idx;
      for (sp_idx = 0; sp_idx < n_nsites0; sp_idx++) {
        s_points[sp_idx].score = LITTLE;
      }

      /* align the top nsites0 subsequences for each value
         of nsites0 and save the alignments with the highest likelihood 
      */
      n_starts += align_top_subsequences(
        mtype,
        w,
        dataset,
        iseq,
        ioff, 
        res+ioff,
        name,
        n_nsites0,
        n_maxima,
        maxima,
        col_scores,
        s_points
      );

      /* A string version of the current seed is required for updating the
         S_POINT heaps: */
      str_seed = to_str_seed(res+ioff, w);

      /* For each of the S_POINT objects, add the current seed to that
         S_POINT'S heap.
         Also, branching search will require a hash_table of all seeds that
         have been evaluated prior to when branching search is called. Hence
         also record the current seed (string, nsites0) combination in the
         hash_table, for all nsites0, unless that seed was already in the
         hash_table:
      */
      hash_insert_str(str_seed, evaluated_seed_ht);
      update_s_point_heaps(s_points, str_seed, n_nsites0);

      myfree(str_seed);
    } /* subsequence */
  } /* sequence */

#ifdef PARALLEL
  reduce_across_heaps(s_points, n_nsites0);
#endif // PARALLEL 

  // Print the sites predicted using the seed after subsequence search, for
  // each of the starting points, if requested:
  if (dataset->print_pred) {
    int sp_idx;
    for (sp_idx = 0; sp_idx < n_nsites0; sp_idx++) {
      // Retrieve the best seed, from the heap:
      HEAP *heap = s_points[sp_idx].seed_heap;
      // Only print sites for the s_point if its heap was non-empty:
      if (get_num_nodes(heap) > 0) {
        SEED *best_seed = (SEED *)get_node(heap, get_best_node(heap));
        char *seed = get_str_seed(best_seed);

        /* Print the sites predicted using the motif corresponding to that seed,
           according to the sequence model being used:
        */
        int nsites0 = s_points[sp_idx].nsites0;
        fprintf(stdout,
                "PREDICTED SITES AFTER SUBSEQUENCE SEARCH WITH W = %i "
                "NSITES = %i MOTIF = %i\n", w, nsites0, dataset->imotif);
        int n_maxima = ps(dataset, w); // upper bound on number of maxima
        P_PROB psites = (P_PROB) mymalloc(n_maxima * sizeof(p_prob));
        n_maxima = get_pred_sites(psites, mtype, w, seed, ltheta[1], lmap,
                                  dataset, ic);
        print_site_array(psites, nsites0, stdout, w, dataset);
        myfree(psites);
      } // get_num_nodes > 0
    } //sp_idx
  } // print_pred

  if (TRACE){
    printf("Tested %d possible starts...\n", n_starts);
    }

  myfree(maxima);
} // subseq7


/**********************************************************************/
/*
	next_pY

	Compute the value of p(Y_ij | theta_1^{k+1})
	from p(Y_ij | theta_1^{k} and the probability
	of first letter of Y_ij given theta_1^k,
	p(Y_ij^0 | theta_1^k).
*/
/**********************************************************************/
static void next_pY(
  DATASET *dataset,			/* the dataset */
  LOG_THETAG_TYPE(theta_1),		/* integer log theta_1 */
  int w,				/* width of motif */
  int *theta_0,				/* first column of previous theta_1 */
  int pYindex				/* which pY array to use */
) {
  int i, k;
  int *theta_last = theta_1[w-1];	/* last column of theta_1 */
  int n_samples = dataset->n_samples;
  SAMPLE **samples = dataset->samples;
  
  for (i=0; i < n_samples; i++) { 	/* sequence */
    SAMPLE *s = samples[i];		/* sequence */
    int lseq = s->length;		/* length of sequence */
    char *res = pYindex<2 ? s->res : s->resic;	/* integer sequence */
    int *pY = s->pY[pYindex];		/* log p(Y_j | theta_1) */
    char *r = res+lseq-1;		/* last position in sequence */
    char *r0 = res+lseq-w-1;	        /* prior to start of last subsequence */
    int j, p;

    if (lseq < w) continue;		/* skip if sequence too short */

    /* calculate p(Y_ij | theta_1) */
    int *pY_shifted_1 = pY - 1;
    for (j=lseq-w; j>0; j--) {
      pY[j] = pY_shifted_1[j] + theta_last[(int)(*r--)] - theta_0[(int)(*r0--)];
    }

    /* calculate log p(Y_i0 | theta_1) */
    p = 0;
    r = res;
    for (k=0; k<w; k++) {     		/* position in site */
      p += theta_1[k][(int)(*r++)];
    }
    pY[0] = p;
  }
}

/**********************************************************************/
/*
	get_max

	Find the erased maxima of pY.  
	If add_psp_to_log_not_o() has been called, "erasing" includes
	the scaled_to_one PSP.

	Returns the length of the (sorted) list of maxima.

	For the oops and zoops models, one maximum is found per sequence.

	For the tcm model, all non-overlapping local maxima are found.
*/
/**********************************************************************/
extern int get_max(
  MOTYPE mtype,		/* the model type */
  DATASET *dataset,	/* the dataset */
  int w,		/* width of sites */ 
  P_PROB maxima, 	/* array of encoded site starts of local maxima */
  BOOLEAN ic, 		/* use reverse complement, too */
  BOOLEAN sort		/* sort the maxima */
)
{
  int n_maxima;

  /* find the maxima */
  if (mtype == Oops || mtype == Zoops) {
    n_maxima = global_max(dataset, w, maxima, ic);
  } else {
    n_maxima = local_max(dataset, w, maxima, ic);
  }

  /* sort the local maxima of pY[1] */
  if (sort) qsort((char *) maxima, n_maxima, sizeof(p_prob), pY_compare);

  return n_maxima;
} /* get_max */

/**********************************************************************/
/*
	global_max	

	Find the position in each sequence with the globally maximal
	value of log pY + log_not_o.

	Returns the number of maxima found and
 	the updated array of maxima positions.
*/
/**********************************************************************/
static int global_max(
  DATASET *dataset,	/* the dataset */
  int w,		/* length of sites */ 
  P_PROB maxima, 	/* array of encoded site starts of local maxima */
  BOOLEAN ic 		/* use reverse complement, too */
)
{
  int i, j;
  SAMPLE **samples = dataset->samples;		/* datset samples */
  int n_samples = dataset->n_samples;		/* number samples in dataset */
  int n_maxima = 0;				/* number of maxima found */

  /* find the position with maximum pY in each sequence */
  for (i=0; i<n_samples; i++) {			/* sequence */
    SAMPLE *s = samples[i];
    int lseq = s->length;
    int last_j = lseq-w;			/* start of last subseq */
    int *pY = s->pY[0];				/* log p(Y_j | theta_1) */
    char *pYic = s->pYic;			/* site on - strand */
    int *log_not_o = s->log_not_o;		// Pr(site) * Pr(no overlap)
    int max_j = 0;				/* best offset */
    int max = pY[0] + log_not_o[0]; 		/* initial maximum */

    if (lseq < w) continue;			/* skip if too short */

    for (j=0; j<=last_j; j++) {			/* subsequence */
      if (pY[j] + log_not_o[j] > max) {		/* new maximum found */
        max_j = j;
        // FIXME: We are assumming that priors are always symmetrical here.
	max = pY[j] + log_not_o[j]; 		// log (pY * Pr(site) * Pr(no overlap))
      } /* new maximum */
    } /* subsequence */

    /* record the maximum for this sequence */
    maxima[n_maxima].x = i;
    maxima[n_maxima].y = max_j;
    maxima[n_maxima].ic = ic && pYic[max_j];	/* on - strand */
    maxima[n_maxima].prob = max;
    n_maxima++;
  } /* sequence */

  return n_maxima;
} /* global_max */

/**********************************************************************/
/*
	local_max

	Find the local maxima of pY * log_not_o 
	subject to the constraint that they are separated by at 
	least w positions. 

	Returns the number of local maxima found and the
	updated array of maxima positions.
*/
/**********************************************************************/
static int local_max(
  DATASET *dataset,	/* the dataset */
  int w,		/* length of sites */ 
  P_PROB maxima,  	/* array of encoded site starts of local maxima */
  BOOLEAN ic		/* use reverse complement, too */
)
{
  int i, j, k, next_j, n_maxima;
  SAMPLE **samples = dataset->samples;		/* datset samples */
  int n_samples = dataset->n_samples;		/* number samples in dataset */

  /* Find the non-overlapping local maxima of p(Y_ij | theta_1) */
  n_maxima = 0;
  for (i=0; i<n_samples; i++) {			/* sequence */
    SAMPLE *s = samples[i];
    int lseq = s->length;			/* length of sequence */
    int *pY = s->pY[0];				/* log p(Y_j | theta_1) */
    int *log_not_o = s->log_not_o;		// Pr(site) * Pr(no overlap)
    int last_j = lseq-w;			/* last possible site */
    int max = pY[0]+log_not_o[0]; 		/* initial maximum */

    if (lseq < w) continue;			/* skip if too short */

    maxima[n_maxima].x = i;			/* candidate */
    maxima[n_maxima].y = 0;			/* candidate site */
    maxima[n_maxima].prob = max;
    next_j = MIN(w, last_j+1);			/* next possible maximum */

    for (j=0; j<=last_j; j++) {			/* subsequence */
      // FIXME: We are assumming that priors are always symmetrical here.
      int prob = pY[j]+log_not_o[j]; 		/* log (pY * Pr(site) * Pr(no overlap)) */
      if (j==next_j) n_maxima++;		/* candidate not exceeded */
      if (j==next_j || prob>max) {		/* create/overwrite */
        max = prob;				/* new max */
        maxima[n_maxima].x = i;			/* overwrite the candidate */
        maxima[n_maxima].y = j;			/* site */
        maxima[n_maxima].prob = max;		/* max */
        next_j = MIN(j+w, last_j+1);		/* next possible maximum */
      } /* create/overwrite candidate */
    }
    n_maxima++;					/* record last maxima */
  }

  /* set the strand and position */
  for (k=0; k<n_maxima; k++) {
    int i = maxima[k].x;			/* site position */
    int j = maxima[k].y;			/* site position */
    SAMPLE *s = samples[i];			/* sequence record */
    maxima[k].ic = ic && s->pYic[j];		/* on - strand */
  } /* n_maxima */

  return n_maxima;
} /* local_max */

/**********************************************************************/
/*
        pY_compare

        Compare the pY of two start sequences.  Return <0 0 >0
        if the first pY is <, =, > the first pY.
*/
/**********************************************************************/
extern int pY_compare(
  const void *v1, 
  const void *v2 
)
{
  double result;

  const struct p_prob * s1 = (const struct p_prob *) v1; 
  const struct p_prob * s2 = (const struct p_prob *) v2; 

  if ((result = s2->prob - s1->prob) != 0) {
    return (result<0) ? -1 : +1;
  } else if ((result = s2->x - s1->x) != 0) {
    return result;
  } else {
    return s2->y - s1->y;
  }
}

/**********************************************************************/
/*
	init_theta_1

	Convert a subsequence to a motif.

	Uses globals:
*/
/**********************************************************************/
static void init_theta_1(
  int w,			/* width of site */
  char *res,			/* (encoded) letters of subsequence */
  LOG_THETAG_TYPE(theta_1),	/* theta_1 */
  int lmap[MAXALPH][MAXALPH]  	/* matrix of frequency vectors */ 
)
{
  int m;
  for (m=0; m<w; m++) {
    theta_1[m] = lmap[(int)res[m]];
  }
} /* init_theta_1 */

/**********************************************************************/
/*
	score_llr_pop

     	Align the top nsites0 subsequences for each value
	of nsites0 and save the alignments with the highest 
        product of p-values of log likelihood ratio of the columns.
	Saves -log_pop as the score for the start.

	Returns number of values of nsites0 tried.
*/ 
/**********************************************************************/
static int score_llr_pop(
  MOTYPE mtype,				/* type of model */
  int w,				/* width of motif */
  DATASET *dataset,			/* the dataset */
  int iseq,				/* sequence number of starting point */
  int ioff,				/* sequence offset of starting point */
  char *eseq,				/* integer encoded subsequence */
  char *name,				/* name of sequence */
  int n_nsites0,			/* number of nsites0 values to try */
  int n_maxima,				/* number of local maxima */
  P_PROB maxima,			/* sorted local maxima indices */
  double *col_scores,			/* column scores for last start point */
  S_POINT s_points[]			/* array of starting points */
)
{
  int i, j, k, i_nsites0;
  int next_seq;				/* index of next subsequence to align */
  int n_starts = 0;			/* number of nsites0 tried */
  int nsites0;				/* starting nsites rounded down */
  int alength = dataset->alength;	/* lenght of alphabet */
  double *back = dataset->back;		/* background frequencies */
  SAMPLE **samples = dataset->samples;	/* the sequences */
  double counts[MAXSITE][MAXALPH];	/* array to hold observed counts */
  double wN;				/* weighted number of sites */
  double log_pop;			/* log product of p-values */
  double min_ic = dataset->min_ic;	/* min. per-column IC */

  /* initialize letter counts to 0 */
  wN = 0;				/* weighted number of sites */
  for (i=0; i<w; i++) for (j=0; j<alength; j++) { counts[i][j] = 0; }

  /* calculate the product of p-values of information content
     of the top nsite0 probability positions 
  */
  for (i_nsites0=0, next_seq=0; i_nsites0 < n_nsites0; i_nsites0++) {

    /* don't score this start if not enough maxima found */
    nsites0 = (int) s_points[i_nsites0].nsites0;	/* round down */
    if (n_maxima < nsites0) {
      continue;
    }
    n_starts++;					/* number of nsites0 tried */

    /* Align the next highest probability sites 
	1) count the number of occurrences of each letter in each column 
	   of the motif and, 
        2) compute the log likelihood of the sites under the background model
    */
    for (k=next_seq; k<nsites0; k++) {		/* site */
      int jj;
      BOOLEAN ic = maxima[k].ic;		/* on - strand */
      int y = maxima[k].y;			/* position of site */
      SAMPLE *s = samples[maxima[k].x];		/* sequence */
      int off = ic ? s->length-w-y : y;		/* - strand offset from rgt. */
      char *res = ic ? s->resic+off : s->res+off;	/* integer sequence */
      double sw = s->sw;			/* sequence weight */
      //
      // TLB: Note that log_not_o contains Pr(site) scaled to have max=1
      // when called from subseq7() but not when called from discretize().
      //
      // Why not revert to not_o[y] here?  TLB: Because the other one works
      // much better, although its kind of a hack.
      //double esw = sw * s->not_o[y];		// Pr(site not overlapped)
      //
      // FIXME: We are assumming that priors are always symmetrical here.
      double esw = sw * INT_DELOG(s->log_not_o[y]);	// Pr(site not overlapped) * Pr(site) 
      wN += esw;				/* total sequence wgt */

      /* residue counts */
      for (j=0; j<w; j++) {			/* position in sequence */
        int c = res[j];
        if (c < alength) {			/* normal letter */
          counts[j][c] += esw;
	} else {				/* 'X' : esw * back[letter] */
          for (jj=0; jj<alength; jj++) counts[j][jj] += esw * back[jj];
	}
        
      } /* position */

    } /* site */
    next_seq = k;				/* next site to align */
    
    /* 
      For DNA palindromes, combine the counts in symmetrically opposing columns
    */
    if (dataset->pal) palindrome(counts, counts, w, alength);

    // Updated on 13-12-06: Only calculate objective function score if the
    // current s_point is supposed to be evaluated:
    if (s_points[i_nsites0].evaluate) {
      /* 
	convert COUNTS to FREQUENCIES and calculate log likelihood ratio
      */
      log_pop = 0;				/* product of p-values */
      for (i=0; i<w; i++) {			/* position in site */
	double llr = 0;				/* log-like-ratio of column */
	double log_pv;				/* log of column p-value */
	double ic;

	/* compute log likelihood for position i */
	for (j=0; j<alength; j++) {		/* letter */
	  double f = wN ? counts[i][j] / wN : 1; 	/* observed letter frequency */
	  double p = back[j];			/* backgrnd letter frequency */
	  double log_f = LOGL(f);
	  double log_p = LOGL(p);
	  double llr_ij = (f&&p) ? f*(log_f - log_p) : 0;
	  llr += llr_ij;
	} /* letter */
	RND(llr/0.6934, RNDDIG, ic);		/* info content in bits */
	llr *= wN;				/* convert entropy to ll */ 
	RND(llr, RNDDIG, llr);			/* round to RNDDIG places */
	log_pv = get_llr_pv(llr, wN, 1, LLR_RANGE, 1.0, alength, back); 
	if (ic < min_ic) log_pv = 0; 		/* ignore low ic columns */

	if (dataset->use_llr) {
	  // Using llr instead of pop:
	  col_scores[i] = log_pv;
	  log_pop -= llr;
	} else {
	  log_pop += col_scores[i] = log_pv;
	}
      } /* position in site */
      RND(log_pop, RNDDIG, log_pop);

      /* print the start sequence and other stuff */
      if (TRACE) {
	if (eseq) {
	  char seq[MAXSITE+1];
	  r2seq(seq, eseq, w);
	  fprintf(stdout, 
	    "( %3d %3d ) ( %*.*s ) %.*s logpop %8.3f nsites0 %6d\n",
	    iseq+1, ioff+1, MSN, MSN, name, w, seq, -log_pop, nsites0);
	} else {
	  fprintf(stdout, 
	    "l_off %3d w %d logpop %8.3f nsites0 %6d\n",
	    iseq, w, -log_pop, nsites0);
	}
      }

      /* save the best start */
      if (-log_pop > s_points[i_nsites0].score) {
	/* Save the starting point and offset so we can re-calculate
	   eseq later. */
	s_points[i_nsites0].iseq = iseq;
	s_points[i_nsites0].ioff = ioff;
	s_points[i_nsites0].e_cons0 = eseq;
	s_points[i_nsites0].wgt_nsites = wN;
	s_points[i_nsites0].score = -log_pop;
      }
    } // Evaluating only if told to do so.
  } /* nsites0 */

  return n_starts;
} /* score_llr_pop */


/**********************************************************************/
/*
	align_top_subsequences

     	Align the top nsites0 subsequences for each value
	of nsites0 and save the alignments with the best values
	according to the objective function.

	Returns number of values of nsites0 tried.
*/ 
/**********************************************************************/
extern int align_top_subsequences(
  MOTYPE mtype,				/* type of model */
  int w,				/* width of motif */
  DATASET *dataset,			/* the dataset */
  int iseq,				/* sequence number of starting point */
  int ioff,				/* sequence offset of starting point */
  char *eseq,				/* integer encoded subsequence */
  char *name,				/* name of sequence */
  int n_nsites0,			/* number of nsites0 values to try */
  int n_maxima,				/* number of local maxima */
  P_PROB maxima,			/* sorted local maxima indices */
  double *col_scores,			/* column scores for last start point */
  S_POINT s_points[]			/* array of starting points */
)
{
  return(score_llr_pop(mtype, w, dataset, iseq, ioff, eseq, name, 
      n_nsites0, n_maxima, maxima, col_scores, s_points));
} /* align_top_subsequences */

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
