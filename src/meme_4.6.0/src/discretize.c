/*
 * $Id: discretize.c 4278 2009-12-23 09:58:37Z james_johnson $
 * 
 */

/*#define DEBUG*/
/***********************************************************************
*								       *
*	MEME++							       *
*	Copyright 2000-2006, The Regents of the University of California    *
*	Author: Timothy L. Bailey				       *
*								       *
***********************************************************************/
#include "meme.h"
#include "dpalign.h"

/* rounding stuff */
#define RNDEPS 1e-12

static int ma_adjust(
  P_PROB sites,                                 /* the sites */
  int nsites,                                   /* the number of sites */
  int w,                                        /* width of sites */
  int flank,                                    /* add flank cols on left+rgt */
  int min_w,		 			/* minimum width allowed */
  DATASET *dataset,  				/* the dataset */
  int *off 					/* best offset */
);

/***********************************************************************/
/*
	get_best_nsites

	Find the value of nsites that minimizes the significance of the
	log likelihood ratio.

	Returns nsites.
*/
/***********************************************************************/
static int get_best_nsites (
  MODEL *model,					/* the model */
  DATASET *dataset,				/* the dataset */
  int min_nsites,				/* minimum sites */
  int max_nsites,				/* minimum sites */
  int w, 	 				/* width of motif */
  int n_maxima,					/* number of maxima */
  P_PROB maxima,				/* the maxima positions */
  double *col_scores,				/* column scores */
  double *best_log_ev				/* best E-value */
)
{
  int i;
  MOTYPE mtype = model->mtype;			/* type of model */
  int n_nsites; 				/* number of different nsites */
  int nsites;					/* optimal number of sites */
  S_POINT *s_points = NULL;			/* for use with align_top...*/

  /* limit maximum nsites */
  max_nsites = MIN(n_maxima, max_nsites);	
  n_nsites = max_nsites-min_nsites+1;	

  /* create array for use by align_top_subsequences */
  Resize(s_points, n_nsites, S_POINT);

  /* initialize the s_points */
  for (i=0, nsites=min_nsites; i<n_nsites; i++, nsites++) {
    s_points[i].nsites0 = nsites;		/* number of sites */
    s_points[i].score = LITTLE;			/* no score yet */
    s_points[i].evaluate = TRUE;                /* Evaluate at every s_point */
  }

  // Renormalize the prior.
  get_not_o(dataset, model->w);
  // TLB: Works better without the following line on Yeast examples.
  //add_psp_to_log_not_o(dataset, model->w, model->invcomp, model->mtype);

  /* 
     align the top nsites sorted subsequences and compute the objective 
     function on each alignment 
  */
  (void) align_top_subsequences(mtype, w, dataset, 0, 0, 0, 0, n_nsites, 
    n_maxima, maxima, col_scores, s_points);

  /* 
     determine the significance of the score for each number of 
     sites and chose the number of sites to minimize it 
  */
  *best_log_ev = BIG;
  nsites = 0;
  for (i=0; i<n_nsites; i++) {			/* nsites */
    double score = s_points[i].score;		/* -log_pop */
    int N = s_points[i].nsites0;
    double wN = s_points[i].wgt_nsites; 
    double log_ev = get_log_sig(score, model->mtype, w, wN, N, model->invcomp,
      model->pal, dataset);
    RND(log_ev, RNDDIG, log_ev);

    if (TRACE) printf("w %d N %d wN %f log_ev %f \n", w, N, wN, log_ev);

    /* save if best so far */
    if (RNDEPS < *best_log_ev - log_ev) {	/* best E-value */
      nsites = N;				/* number of sites */
      *best_log_ev = log_ev;
    }
  } /* nsites */

  myfree(s_points);

  return nsites;
} /* get_best_nsites */

/***********************************************************************/
/*
	set_z

	Set the z to 1/0 using list of sites (maxima).
*/
/***********************************************************************/
extern void set_z (
  MODEL *model,					/* the model */
  DATASET *dataset 				/* the dataset */
)
{
  int i, j;
  int nsites = model->nsites_dis;		/* new nsites */
  P_PROB maxima = model->maxima;		/* the maxima positions */
  int n_samples = dataset->n_samples;		/* number of sequences */
  SAMPLE **samples = dataset->samples;		/* sequences */
  BOOLEAN invcomp = model->invcomp;

  /* set all z to 0 */
  for (i=0; i<n_samples; i++) {			/* sequence */
    SAMPLE *s = samples[i];			/* sample i */
    int lseq = s->length;			/* length of sequence */
    int min_j = invcomp ? -lseq : 0;		// minimum Z_i
    int max_j = lseq;				// maximum Z_i
    double *zi = s->z;				// zi[j], j in [-lseq...+lseq] 
    for (j=min_j; j<=max_j; j++) {		// Z_i = j
      Zi(j) = 0;
    }
  }

  /* set z 1 for selected sites */
  for (i=0; i<nsites; i++) {
    SAMPLE *s = samples[maxima[i].x];		/* sample */
    int y = maxima[i].y;			/* position of site */
    double *zi = s->z;				// zi[j], j in [-lseq...+lseq] 
    int j = maxima[i].ic ? -(y+1) : y+1;	// value of Z_i
    Zi(j) = 1.0;
  }

} /* set_z */

/***********************************************************************/
/*
	set_pY

	Initialize pY from z with given offset.
*/
/***********************************************************************/
static void set_pY(
  int w, 				/* motif width */
  BOOLEAN invcomp, 			/* use reverse complement strand, too */
  BOOLEAN pal,				/* force palindrome */
  DATASET *dataset,			/* the dataset */
  int off				/* offset to shift motif */
)
{
  int i, j;
  int n_samples = dataset->n_samples;		/* number of sequences */
  SAMPLE **samples = dataset->samples;		/* sequences */

  /* 
    put integerized, weighted log z into pY array
  */
  for (i=0; i<n_samples; i++) {			/* sequence */
    SAMPLE *s = samples[i];			/* sequence i */
    int lseq = s->length;			/* sequence length */
    int last_j = lseq-w;			/* last start */
    double *zi = s->z;				// zi[j], j in [-lseq...+lseq] 
    double sw = s->sw;				/* weight of sequence */
    int *pY = s->pY[0];				/* p(Y_j | theta_1) both */
    char *pYic = s->pYic;			/* site on - strand */

    if (lseq < w) continue;			/* sequence too short */

    /* initialize pY and pYic */
    for (j=0; j<=last_j; j++) {
      pY[j] = INT_LOG(0.0);			/* z == 0 */
      pYic[j] = '\0';				/* site on + strand */
    }
    for (j=0; j<lseq; j++) {			/* site start */
      int jj = j + off;				/* new site start */
      int k = jj+1;				// Z_i = k

      if (jj<0 || jj>last_j) continue;

      /* no z available? */
      if (j > last_j) {				/* no z available */
        pY[jj] = 0;				/* no site */
        pYic[jj] = '\0';			/* strand doesn't matter */
        continue;
      } 

      /* not using inverse strand, too? */
      if (!invcomp) {
        pY[jj] = INT_LOG(sw * Zi(k));
        pYic[jj] = '\0';			/* site on + strand */
        continue;
      } 

      /* using inverse complement strand, too */
      if (pal) {				// use sum of Zi(-k)+Zi(k)
        pY[jj] = INT_LOG(sw * MIN(1.0,(Zi(-k)+Zi(k))));	// FIXME??
      } else if (Zi(-k) > Zi(k)) {		// - strand
        pY[jj] = INT_LOG(sw * Zi(-k));
      } else {					// + strand
        pY[jj] = INT_LOG(sw * Zi(k));
      }
      pYic[jj] = (Zi(-k) > Zi(k)) ? '\1' : '\0';	/* choose strand */

    } /* site start */

  } /* sequence */

} /* set_pY */

/***********************************************************************/
/*
	discretize	

	Search over width and offset of motif to minimize E-value of
	log likelihood ratio.  

		1) get best nsites using E-value
		2) calculate p-value of each column of motif
		3) shorten using p-value
		4) get best nsites using E-value

	Sets z to 1.0 for sites, 0 for non-sites in dataset.
	Sets w, nsites_dis, maxima and sig in model.

	Returns the optimum number of sites.
*/
/***********************************************************************/
extern double discretize(
  MODEL *model,				/* the model */
  DATASET *dataset			/* the dataset */
)
{
  int i;
  int n_samples = dataset->n_samples;		/* number of sequences */
  BOOLEAN invcomp = model->invcomp;     	/* reverse complement strand */
  BOOLEAN pal = model->pal;     		/* force DNA palindromes */
  MOTYPE mtype = model->mtype;			/* type of model */
  int w, best_w, min_w, max_w, ini_w, ma_w;	/* motif width */
  int min_nsites = model->min_nsites;		/* minimum nsites */
  int max_nsites = model->max_nsites;		/* maximum nsites */
  int n_maxima;		 			/* number of possible sites */
  double log_pop;				/* log product of p-values */
  double log_pv, best_log_pv;			/* log p-value of motif */
  double log_ev, best_log_ev;			/* log E-value of score */
  int ma_off=0;					/* multiple alignment offset */
  int best_off;					/* motif offset */
  int nsites=0, best_nsites=0;			/* number of sites */
  double col_scores[MAXSITE];			/* column scores */

  /* 
    set minimum and maximum widths
  */
  ini_w = best_w = model->w;			/* initial motif width */
  min_w = model->min_w;				/* minimum width */
  max_w = ini_w;				/* maximum width */
  w = 0;

  /* create space for the maxima */
  n_maxima = (mtype==Tcm) ? ps(dataset, min_w) : n_samples;
  Resize(model->maxima, n_maxima, p_prob);

  /* 
    initialize pY from z 
  */
  set_pY(model->w, invcomp, pal, dataset, 0);

  // Don't include the PSP probabilities in log_not_o (used in get_max())
  // because z already takes them into account.
  get_not_o(dataset, model->w);

  /* 
    get the possible sites using maximum width
  */
  n_maxima = get_max(mtype, dataset, max_w, model->maxima, invcomp, FALSE);
  if (n_maxima < min_nsites) return nsites;

  /* 
    sort them by z 
  */
  for (i=0; i<n_maxima; i++) {
    int x = model->maxima[i].x;
    SAMPLE *s = dataset->samples[x];
    double *zi = s->z;				// Zi[j], j in [-lseq...+lseq]
    int y = model->maxima[i].y;
    int j = y + 1;				// Z_i = j
    double z = invcomp ? MIN(1.0,Zi(-j)+Zi(j)) : Zi(j);
    RND(z, 11, model->maxima[i].prob);
  }
  qsort((char *) model->maxima, n_maxima, sizeof(p_prob), pY_compare);

#ifdef DEBUG
#ifdef PARALLEL
#undef printf
#endif
  //printf("nmax = %d start %d %d %3.0f %3.0f %s\n", n_maxima, 
  //   model->iseq, model->ioff, model->pw, model->psites, model->cons0);
  //print_theta(1000*model->iseq+model->ioff, 1, model->nsites_dis,
  //  model->theta, model->w, 0, "", dataset, stdout);
  for (i=0; i<n_maxima; i++) {
    int x = model->maxima[i].x;
    int y = model->maxima[i].y;
    double prob = model->maxima[i].prob;
    SAMPLE *s = dataset->samples[x];
    printf("max %d %15.12f x %d y %d psp %.4g\n", i, prob, x, y);
  }
  fflush(stdout);
#ifdef PARALLEL
#define printf if (mpMyID() == 0) printf
#endif
#endif

  /* 
    get the best number of sites using the full motif width
  */
  best_nsites = (min_nsites==max_nsites) ? min_nsites :
    get_best_nsites(model, dataset, min_nsites, max_nsites,
     max_w, n_maxima, model->maxima, col_scores, &log_ev);
  if (best_nsites < min_nsites) return nsites;

  /* trim the alignment to include the minimum possible number of gaps */
  ma_w = ini_w;
  if (dataset->ma_adj) {
    int flank = w/2;				/* amt. of flank for mult. a. */
    ma_w = ma_adjust(model->maxima, best_nsites, max_w, flank, min_w, dataset, 
      &ma_off);
    max_w = ma_w;				/* new maximum width */
    /*
      update the maxima positions by shifting them
    */
    for (i=0; i<best_nsites; i++) {
      BOOLEAN ic = model->maxima[i].ic;			/* on - strand */
      model->maxima[i].y += (ic ? ini_w-ma_w-ma_off : ma_off);
    }
  } /* ma_adj */

  /* 
    get the p-values of the columns given the best number of sites
    and the gap-trimmed width
  */
  (void) get_best_nsites(model, dataset, best_nsites, best_nsites,
    max_w, n_maxima, model->maxima, col_scores, &log_ev);

  /* 
    shorten the motif: find subsequence of columns with best p-value
  */
  best_log_pv = best_log_ev = BIG;
  best_off = 0;
  for (w=min_w; w<=max_w; w++) {		/* width */
    int l_off;					/* left offset */

    for (l_off=0; l_off<=max_w-w; l_off++) {	/* left edge */

      /* get the product of column p-values */
      for (i=log_pop=0; i<w; i++) log_pop += col_scores[i+l_off];

      /* get the p-value of the pop */
      log_pv =
        get_log_sig(-log_pop, mtype, w, best_nsites, 0, invcomp, pal, dataset);

      if (TRACE) 
        printf(
        "ini_w %d ma_w %d w %d ma_off %d off %d log_pv %f init cons %*.*s\n", 
        ini_w, ma_w, w, ma_off, l_off, log_pv, w, w, model->cons0+l_off+ma_off);

      /* save if best so far: better log_pv */
      if (RNDEPS < best_log_pv - log_pv) {
        if (TRACE) printf("better: w %d best log_pv %g\n", w, log_pv);
	best_w = w;
	best_off = l_off;
	best_log_pv = log_pv;
      }
    } /* l_off */
  } /* w  */

  /*
    update the maxima positions by shifting them
  */
  for (i=0; i<best_nsites; i++) {
    BOOLEAN ic = model->maxima[i].ic;			/* on - strand */
    model->maxima[i].y += (ic ? ma_w-best_w-best_off : best_off);
  }
  
  /* 
    get the best number of sites for the shortened motif and the final E-value 
  */
  best_nsites =
    get_best_nsites(model, dataset, min_nsites, best_nsites,
     best_w, n_maxima, model->maxima, col_scores, &best_log_ev);

  /* 
    set the best motif info in the model and
    discretize the sites in z using the best width, offset and nsites 
  */
  model->w = best_w;				/* best width */
  model->nsites_dis = best_nsites;		/* after discretization */
  model->logpv = best_log_pv;			/* p-value */
  model->logev = best_log_ev;			/* E-value */
  set_z(model, dataset);

  if (TRACE) 
   printf( 
   "ini_w %d ma_w %d w %d ma_off %d off %d nsites %d pv %9.2e EV %9.2e cons %*.*s %s\n",
    ini_w, ma_w, best_w, ma_off, best_off, best_nsites, exp(best_log_pv), 
    exp(best_log_ev), best_w, best_w, (model->cons0)+best_off, model->cons0);

  return best_nsites; 
} /* discretize */

/**********************************************************************/
/*
	ma_adj

	Shorten a motif to the longest g-alignment of width at least
	min_w.  A g-alignment is an alignment with no more than g
	gapped sequences per column.  Values of g in [0..] are tried
	until an aligment of width min_w or greater is found.

	Returns best width and offset.
*/
/**********************************************************************/
static int ma_adjust(
  P_PROB sites,                                 /* the sites */
  int nsites,                                   /* the number of sites */ 
  int w,                                        /* width of sites */ 
  int flank,                                    /* add flank cols on left+rgt */
  int min_w,		 			/* minimum width allowed */
  DATASET *dataset,  				/* the dataset */
  int *off					/* best offset */
)
{
  char **ma;					/* multiple aligment */
  int left = MIN(flank, sites[0].y);		/* left edge after algnmnt. */
  int right = left+w-1;				/* right edge after algnmnt. */

  /* get the multiple alignment */
  ma = dp_multi_align(sites, nsites, w, flank, dataset);

  /* get longest g-alignment of width at least min_w */
  (void) g_align(ma, nsites, strlen(ma[0]), left, right, min_w, off, &w);

  /* free space */
  free_2array(ma, nsites);

  return w;					/* return the width */
} /* ma_adjust */

