/*
 * $Id: likelihood.c 4278 2009-12-23 09:58:37Z james_johnson $
 * 
 * $Log$
 * Revision 1.2  2005/10/25 19:06:39  nadya
 * rm old macro for Header, all info is taken care of by Id and Log.
 *
 * Revision 1.1.1.1  2005/07/29 00:21:43  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/*#define DEBUG*/
/**********************************************************************/
/*
	Likelihood and entropy functions.
*/
/**********************************************************************/
/* 7-14-00 tlb; use back for rentropy */
/* 6-23-00 tlb; add Ev objective function */

#ifdef lrt_standalone
#define DEFINE_GLOBALS
#endif
#include "meme.h"

/*
  constants
*/
  static int *len = NULL;

/*
  functions
*/
extern int int_compare(
  const void *v1,
  const void *v2
);
#ifdef OBSOLETE
static double log_non_overlap_comb(
  int m,				/* length of sequence */
  int n,				/* number of segments */
  int w 				/* width of segment */
); 
#endif

/**********************************************************************/
/* 
	get_log_sig

  	Calculate the statistical significance of the alignment given
	its score and the type of objective function in force.  

	If N>0, returns log E-value.
	If N==0, returns log p-value.
*/
/**********************************************************************/
extern double get_log_sig(
  double score,					/* score of alignment */
  MOTYPE mtype,					/* type of model */
  int w,					/* width of motif */
  double wN,					/* weighted number of sites */
  int N,					/* number of sites */
  BOOLEAN invcomp,				/* inv. compl. strand, too */
  BOOLEAN pal, 					/* motif is DNA palindrome */
  DATASET *dataset 				/* the dataset */
)
{
  double log_pv;			/* log of p-value */
  double log_sig = 0;			/* return value */

  log_pv = log_qfast(w, -score);	/* p-value of product of p-values */

  if (N) {				/* use E-value of alignment */
    log_sig = log_pv + get_log_nalign(mtype, w, N, invcomp && !pal, dataset);
  } else { 				/* use p-value of alignment */
    log_sig = log_pv;
  }

  return(log_sig);
} /* get_log_sig */

/**********************************************************************/
/*
	calc_entropy

	Calculate the entropy of each column and relative entropy per
	column of a motif.  Stores results in model.
		model->rentropy
		model->rel
		model->logev
		model->ic
*/
/**********************************************************************/
extern void calc_entropy (
  MODEL *model,			/* the model */
  DATASET *dataset  		/* the dataset */
)
{
  int i, j;
  double *rentropy = model->rentropy;		/* IC of each column */
  double ent = 0;				/* entropy per column */
  double rel = 0;				/* relative entropy per col. */
  int alength = dataset->alength;		/* length of alphabet */
  double *back = dataset->back;			/* background model freqs */
  THETA obs = model->obs;			/* observed frequencies */
  int w = model->w;				/* width of motif */
  int N = model->nsites_dis;			/* number of sites */
  double log_pop;				/* log product of col p-value */
  double max_ic = LOG(alength)/LOG(2);		// maximum IC per column
  //double e = (alength-1) / (2 * LOG(2) * N);	// small sample correction
  double ic = 0;				// "corrected" IC

  /* calculate the relative entropy of each column in motif */
  model->llr = log_pop = 0;
  for (i=0; i<w; i++) {			/* position */
    double llr;				/* log likelihood ratio of column */
    rentropy[i] = 0.0; 
    double H = 0;			// negative entropy in this column
    for (j=0; j<alength; j++) {		/* alphabet letter */
      double f = obs(i, j);		/* motif freq */
      double p = back[j];		/* background freq */
      double h = f ? f * LOG2(f) : 0;	// entropy of current letter
      rel += p ? f * LOG2(p) : 0;	/* total relative entropy */
      ent += h;				/* total entropy */
      H += -h;				// negative entropy in this column
      rentropy[i] += (f && p) ? f * LOG(f/p) : 0;
    } /* alphabet letter */
    llr = N * rentropy[i];		/* log likelihood ratio */
    RND(llr, RNDDIG, llr);		/* round to RNDDIG places */
    model->llr += llr;                  /* llr for model */
    log_pop += get_llr_pv(llr, N, 1, LLR_RANGE, 1.0, alength, back);
    rentropy[i] /= LOG(2);
    // ic += MAX(0, (max_ic - (H + e)));
    ic += max_ic - H;
  } /* position in motif */

  /* compute the log E-value of the motif */
  model->logev = get_log_sig(-log_pop, model->mtype, w, N, N, model->invcomp,
    model->pal, dataset);

  model->rel = (ent - rel)/w;		/* compute rel. entropy/col */

  // LOGO total information content
  RND(ic, RNDDIG, ic);			// round to RNDDIG places
  model->ic = ic;

} /* calc_entropy */

/**********************************************************************/
/*
        log_comb
 
        Compute logarithm of m choose n.

*/
/**********************************************************************/
 
extern double log_comb(
  int m,
  int n
) 
{
  int i;
  double x = 0;
  int big, little;
  if (m-n > n) { big = m-n; little = n; } else { big = n; little = m-n; }
  for (i=m; i>big; i--) x += log((double) i);
  for (i=2; i<=little; i++) x -= log((double) i);
  return x;
}

#ifdef OBSOLETE
/**********************************************************************/
/*
	log_non_overlap_comb

	Compute the log of the (approximate) number of combinations of n 
	non-overlapping segments of width w chosen from a sequence of
	length m.
*/
/**********************************************************************/
static double log_non_overlap_comb(
  int m,				/* length of sequence */
  int n,				/* number of segments */
  int w 				/* width of segment */
) 
{
  int i;
  double x = 0;

  if (m <= w*n) return BIG;		/* too many segments */

  for (i=1; i<=n; i++) {
    x += log(m-(w*i)+1.0) - log(i);	/* each segment w covered */
  }

  return x;
} /* log_non_overlap_comb */
#endif

/**********************************************************************/
/*
	get_log_nalign

	Get an upper bound on the number of independent alignments
	of segments of length w.
*/
/**********************************************************************/
extern double get_log_nalign(
  MOTYPE mtype,					/* type of model */
  int w,	 				/* width of motif */
  int N,					/* number of occurrences */
  BOOLEAN invcomp,				/* inv. compl. seq allowed */
  DATASET *dataset  				/* the dataset */
)
{
  int i, t;
  int nseqs = dataset->n_samples;		/* number of sequences */
  double log_nalign = 0;			/* log number alignments */
  int icfactor = invcomp ? 2 : 1;		/* double the possible sites */

  /* 
    sort the sequence lengths in decreasing order in array len[] first time thru
  */
  if (len == NULL) { 				/* first time through */
    Resize(len, nseqs, int);
    for (i=0; i< nseqs; i++) len[i] = dataset->samples[i]->length;
    qsort((char *) len, nseqs, sizeof(int), int_compare);
  }

  /*
    get upper bound on number of alignments 
  */
  switch (mtype) {
    case Oops:
    case Zoops:
      if (w > len[N-1]) {				/* impossible w */
	log_nalign = BIG;
      } else {
	for (i=0; i<N; i++) log_nalign += log(icfactor*(len[i]-w+1.0));
        if (N < nseqs) log_nalign += log_comb(nseqs, N);
      }
      break;
    case Tcm:
      for (i=t=0; i<nseqs; i++) t += len[i] - w + 1;	/* # starts */
      if (N > t) {					/* impossible w & N */
	log_nalign = BIG;
      } else {					/* remove 1 site per site */
        for (i=0; i<N; i++) log_nalign += log((t-i)*icfactor/(i+1));
      }
      break;
  } /* mtype */

  return log_nalign;
} /* double get_log_nalign */

/**********************************************************************/
/*
	log_qfast
	
	Calculate the log p-value of the log of the 
	product of uniform [0,1] random variables.

*/
/**********************************************************************/
extern double log_qfast(
  int n,			/* number of random variables in product */
  double logk			/* product of random variables */
)
{
  int i;
  double term, phi;
 
  if (n == 0) return 0;			/* worst possible log p-value */

  phi = term = 1;
  for (i=1; i<n; i++) {
    term *= -logk/i;
    phi += term;
  }

  return(logk + log(phi));
} /* qfast */

/**********************************************************************/
/*
        int_compare

        Compare two integers.  Return <0 0 >0
        if the second int is <, =, > the first int.
	For sorting with qsort in decreasing order.
*/
/**********************************************************************/
extern int int_compare(
  const void *v1,
  const void *v2
)
{
  const int * s1 = (const int *) v1;
  const int * s2 = (const int *) v2;

  return(*s2 - *s1);
} /* int_compare */


#ifdef obsolete
/**********************************************************************/
/*
	get_n_unique_multisets

	Get the (log) number of unique multisets of size n containing 
	items of (up to) m different types.
*/
/**********************************************************************/
static double get_n_unique_multisets (
  int n,					/* size of multisets */
  int m 					/* number of item types */
)
{
  int i, j, k;
  double **f;					/* dynamic programming array */

  create_2array(f, double, n+1, m+1)

  /* initialize */
  for (i=0; i<=n; i++) {
    for (j=0; j<=m; j++) {
      if (i==0 || j==1) {
        f[i][j] = 1;
      } else if (i==1) {
        f[i][j] = j;
      } else {
        f[i][j] = 0;
      }
    }
  }

  /* fill in array */
    for (j=2; j<=m; j++) 
      for (i=2; i<=n; i++) 
        for (k=0; k<=i; k++) 
          f[i][j] += f[k][j-1];

  fprintf(stderr, "unique multisets (%d %d) %g\n", n, m, f[n][m]);
  free_2array(f, n);

  return f[n][m];
} /* get_n_unique_multisets */
#endif
