/*
 * $Id: llr.c 2417 2008-01-30 00:06:38Z cegrant $
 * 
 * $Log$
 * Revision 1.2  2006/03/08 20:50:11  nadya
 * merge chamges from v3_5_2 branch
 *
 * Revision 1.1.1.1.4.2  2006/01/26 09:16:27  tbailey
 * Rename local function getline() to getline2() to avoid conflict with
 * function defined in stdio.h.
 *
 * Revision 1.1.1.1.4.1  2006/01/24 20:44:08  nadya
 * update copyright
 *
 * Revision 1.1.1.1  2005/07/29 00:22:18  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/************************************************************************
*	Copyright							*
*	(1999-2006) The Regents of the University of California.	*
*	All Rights Reserved.						*
*	Author: Timothy L. Bailey
*									*
*	Permission to use, copy, modify, and distribute any part of 	*
*	this software for educational, research and non-profit purposes,*
*	without fee, and without a written agreement is hereby granted, *
*	provided that the above copyright notice, this paragraph and 	*
*	the following three paragraphs appear in all copies.		*
*									*
*	Those desiring to incorporate this software into commercial 	*
*	products or use for commercial purposes should contact the 	*
*	Technology Transfer Office, University of California, San Diego,*
*	9500 Gilman Drive, La Jolla, California, 92093-0910, 		*
*	Ph: (619) 534 5815.						*
*									*
*	IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO 	*
*	ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR 	*
*	CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF 	*
*	THE USE OF THIS SOFTWARE, EVEN IF THE UNIVERSITY OF CALIFORNIA 	*
*	HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 		*
*									*
*	THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE *
*	UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE 		*
*	MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  *
*	THE UNIVERSITY OF CALIFORNIA MAKES NO REPRESENTATIONS AND 	*
*	EXTENDS NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED, *
*	INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 	*
*	MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT 	*
*	THE USE OF THE MATERIAL WILL NOT INFRINGE ANY PATENT, 		*
*	TRADEMARK OR OTHER RIGHTS.  					*
************************************************************************/

/***********************************************************************/
/*
	Information Content Statistics Routines
*/
/***********************************************************************/

#ifdef SO
#define DEFINE_GLOBALS
#endif
#include "llr.h"
#include "macros.h"
#include "general.h"
#include "io.h"
#include <stdio.h>
#include "logs.h"
#include "message.h"

/* probability distribution */ 
typedef struct distr {
  int w;			/* maximum width of alignment */
  double alpha;			/* scale factor for distribution */
  int *range;			/* range of distributions w = [1..w] */
  int *offset;			/* offset of distributions w = [1..w] */
  double **d;			/* array of (log) distributions w = [1..w] */
  double **cdf;			/* array of llr (log) 1 - CDFs w = [1..w] */
  double mean;			/* mean of llr for w = 1 */
} DISTR;
static int ndistrs = -1;	/* largest N for which distr known */
static DISTR *distrs = NULL;	/* array of distributions for N = [1..ndistrs]*/

/************************************************************************/
/*
	init_llr_pv_tables

	Initialize the single-column log likelihood ratio p-value tables.

        Note: Palindromes have effectively double the apparent sites but half 
        the width.

*/
/************************************************************************/
extern void init_llr_pv_tables(
  int min,				/* minimum number of sites */
  int max,				/* maximum number of sites */
  int alength,				/* alphabet length */
  double *back,				/* background frequencies */
  BOOLEAN pal				/* sites are palindromes */
)
{
  int nsites;				/* number of sites */

  /* set effective number of sites to double if pal */
  if (pal) { min *= 2; max *= 2; }

  if (!NO_STATUS)
    fprintf(stderr,
      "Initializing the motif probability tables for %d to %d sites...\n",
	min, max);

  /* make sure the distr table gets initialized on all nodes */
  (void) get_llr_pv(0, 1, 1, LLR_RANGE, 1.0, alength, back);

  for (nsites=min; nsites<=max; nsites += pal ? 2 : 1) {   /* nsites */

    /* allocate space for table */
    (void) get_llr_pv(0, nsites, 0, LLR_RANGE, 1.0, alength, back);

    if (!load_balance_llr(nsites, pal)) {
      continue;	/* for parallel */
    } 

    /* create table */
    if (!NO_STATUS) { fprintf(stderr, "nsites = %d\r", nsites); }
    (void) get_llr_pv(0, nsites, 1, LLR_RANGE, 1.0, alength, back);

  } /* nsites */
  broadcast_llr(min, max, pal);		/* for parallel; collect the tables */

#ifdef DEBUG
  /* print results */
  int n;
  for (n=min; n<=max; n++) {
    int I;
    int w = 1;
    printf("# N    I    llr         1-cdf\n");
    for (I=0; I<=distrs[n].range[w]; I++) {		/* LLR */
      double m2, e2;
      if (distrs[n].cdf[w][I] == LOGZERO) {
        m2 = e2 = 0;
      } else {
        exp10_logx(distrs[n].cdf[w][I]/log(10.0), m2, e2, 1);
      }
      printf("%3d %3d %5.1f %3.1fe%+05.0f\n",
        n, I, (distrs[n].offset[w]+I)/distrs[n].alpha, m2, e2);
    } /* LLR */
  }
#endif

  if (!NO_STATUS)fprintf(stderr, "\nDone initializing\n");

} /* init_llr_pv_tables */

/******************************************************************************/
/*
	llr_distr

	Compute the probability distribution of the log-likelihood
	ratio (llr) statistic for a discrete distribution.
	Because of roundoff, LLR may be as small as -1 (instead of 0).

	Returns an array p where
		p[i] = log(Pr(LLR=i))

*/
/******************************************************************************/
double *llr_distr(
  int A,				/* dimension of discrete distribution */
  double *dd,				/* discrete distribution */
  int N,				/* number of samples */
  int desired_range,			/* desired range for scaled LLR */
  double frac,				/* fraction of scores to use */
  double *alpha,			/* scale factor for scaled LLR */
  int *offset,				/* prob[0] = prob(offset) */
  int *range				/* range for scaled LLR */
)
{
  int i; 				/* index over alphabet */
  int n;				/* index over samples */
  int k;				/* other index */
  int I;				/* LLR */
  double dd_sum;			/* sum of dd */
  int **IP;				/* I'_i[n] */
  int *minI=NULL;			/* minimum intermediate value of I */
  int *maxI=NULL;			/* maximum intermediate value of I */
  int Irange;				/* maxI-minI+1 */
  double logNfact;			/* log N! */
  double **logP;			/* log P_i[n] */
  double **logSP;			/* log script_P[# samples][LLR] */
  double *prob=NULL;			/* final probability distribution */
  double min, max, min_dd;

  /* create space for IP, P, minI and maxI */
  create_2array(IP, int, A, N+1);
  create_2array(logP, double, A, N+1);
  Resize(minI, N+1, int);
  Resize(maxI, N+1, int);

  /* make sure distribution sums to 1.0 and has no 0's */
  for (i=dd_sum=0; i<A; i++) dd_sum += dd[i] + EPSILON;
  for (i=0; i<A; i++) dd[i] = (dd[i]+EPSILON)/dd_sum;

  /* compute N! */
  logNfact = 0;
  for (i=2; i<=N; i++) logNfact += log(i); 

  /* get estimates of minimum and miximum values of llr */
  for (i=0, min_dd=1; i<A; i++) min_dd = MIN(min_dd, dd[i]);
  max = NINT(-N * log(min_dd));
  for (i=min=0; i<A; i++) min += dd[i]*N*(log(dd[i]) - log(dd[i]));
  min = NINT(min);
  /*printf("min = %f max = %f\n", min, max);*/

  /* set alpha to achieve the desired range */
  *alpha = desired_range/((max-min));
  /* *alpha = NINT(((int)desired_range)/((max-min)));
  if (*alpha < 1) *alpha = 1;*/
  /*fprintf(stderr, "range %d max %f min %f alpha = %f\n",desired_range, max, min, *alpha);*/

  /* compute I', P, minI and maxI */ 
  for (n=0; n<=N; n++) minI[n] = maxI[n] = 0;
  for (i=0; i<A; i++) {				/* index over alphabet */
    double logdd = LOG(dd[i]);			/* log(dd[i]) */
    IP[i][0] = 0; 				/* I'_i(0) */
    logP[i][0] = 0;				/* log P_i(0) */
    for (n=1; n<=N; n++) {			/* index over samples */
      IP[i][n] = NINT(*alpha*n*log(n/(N*dd[i]))); 	/* I'_i(n) */
      logP[i][n] = logP[i][n-1] + logdd - log(n);	/* log P_i(n) */
      for (k=1; k<=n; k++) {			/* index over samples of new */
	minI[n] = MIN(minI[n], minI[n-k] + IP[i][k]);
	maxI[n] = MAX(maxI[n], maxI[n-k] + IP[i][k]);
      }
    }
  }

  /* get overall minI and maxI */
  for (n=1; n<=N; n++) {
    /*printf("minI[%d] %d maxI[%d] %d\n", n, minI[n], n, maxI[n]);*/
    minI[0] = MIN(minI[0], minI[n]);		/* min for intermediates */
    maxI[0] = MAX(maxI[0], maxI[n]);		/* max for intermediates */
    minI[n] = LOGZEROI;
    maxI[n] = 0;
  }
  Irange = maxI[0] - minI[0] + 2;
  *offset = minI[0] - 1;			/* I offset: I=-1 is array 0 */
  /*printf("minI %d maxI %d Irange %d\n", minI[0], maxI[0], Irange);*/
  minI[0] = LOGZEROI;
  maxI[0] = 0;

  /* create script_P arrays with enough space for intermediate calculations */
  create_2array(logSP, double, N+1, Irange+1);
  
  /* clear intermediate probability array */
  for (n=0; n<=N; n++) for(I=0; I<Irange; I++) logSP[n][I] = LOGZERO;

  /* init probability array for first letter in alphabet */
  for (n=0; n<=N; n++) {
    I = IP[0][n] - *offset;			/* offset I */
    logSP[n][I] = logNfact + logP[0][n];	/* init */
    minI[n] = maxI[n] = I;
  }

  /* compute probabilities recursively */
  for (i=1; i<A; i++) {			/* index over (rest of) alphabet */
    for (n=N; n>=0; n--) {		/* index over samples */
      for (k=1; k<=n; k++) {		/* index over samples of new letter */
        int min = minI[n-k];
        int max = MAX(min, maxI[n-k] - (1-frac)*(maxI[n-k]-minI[n-k]+1));
        /*printf("min %d maxI %d max %d\n", min, maxI[n-k], max);*/
        for (I=min; I<=max; I++) {	/* index over I */
          if (logSP[n-k][I] > LOGZERO) {
	    /*printf("i %d old: %d %d new: %d %d\n", i, n-k, I, n,I+IP[i][k]);*/
            logSP[n][I+IP[i][k]] = 
              LOGL_SUM(logSP[n][I+IP[i][k]], logP[i][k] + logSP[n-k][I]);
          }
	}
	/* get current minimum and maximum I in intermediate arrays */
	minI[n] = MIN(minI[n], minI[n-k]+IP[i][k]);
	maxI[n] = MAX(maxI[n], maxI[n-k]+IP[i][k]);
      }
      if (n==N && i==A-1) break;	/* all done */
    }
  }

  /* compute range */
  /*printf("minI[N] %d maxI[N] %d\n", minI[N], maxI[N]);*/
  *range = maxI[N] - minI[N]; 

  /* move to probability array with prob(offset) in position 0 */
  *offset += minI[N];			/* prob[0] = prob(offset) */
  Resize(prob, *range+2, double);
  for (I=minI[N]; I<=maxI[N]; I++) prob[I-minI[N]] = logSP[N][I];
  /*fprintf(stderr, "N= %d range= %d offset= %d alpha= %f\n", N, *range, 
    *offset, *alpha);*/
      
  /* free up space */
  free_2array(IP, A);
  free_2array(logP, A);
  free_2array(logSP, N+1);
  myfree(minI);
  myfree(maxI);

  return prob;
} /* llr_distr */

/******************************************************************************/
/*
	sum_distr

	Compute the distribution of the sum of two integer-valued
	random variables whose ranges are [0..r1] and [0..r2], 
	respectively, given the log of their distributions.
*/
/******************************************************************************/
double *sum_distr(
  double *d1,				/* (log) distribution of RV1 */
  int r1,				/* range of RV1 */
  double *d2,				/* (log) distribution of RV2 */
  int r2, 				/* range of RV2 */
  int *r_sum				/* range of sum of RV1 and RV2 */
)
{
  int i, j, k;
  int range = r1 + r2;			/* potential range of sum */
  double *d_sum = NULL;			/* distribution of sum */

  Resize(d_sum, range+1, double);	/* space for distribution */

  for (i=0; i<=range; i++) {		/* value of sum */
    d_sum[i] = LOGZERO;
  }

  for (i=0; i<=r1; i++) {		/* range of RV1 */
    if (d1[i]==LOGZERO) continue;
    for (j=0, k=i; j<=r2; j++, k++) {	/* range of RV2 */
      if (d2[j]==LOGZERO) continue;
      d_sum[k] = LOGL_SUM(d_sum[k], d1[i]+d2[j]);
    } /* RV2 */
  } /* RV1 */

  for (i=range; i>=0; i--) {		/* value of sum */
    if (d_sum[i] > LOGZERO) break;
  }
  *r_sum = i;				/* non-zero range */

  return d_sum;
} /* sum_distr */

/******************************************************************************/
/*
	cdf

	Compute (log) 1-CDF of an integer-valued (log) distribution.
	Smooths the CDF by linear interpolation so that adjacent positions
	in the table will have different values if possible.
*/
/******************************************************************************/
static double *cdf(
  double *d,				/* integer valued distribution */
  int r					/* range [0..r] */
)
{
  double *cdf=NULL, slope=0;
  int I, i, j, k;

  Resize(cdf, r+1, double);
  cdf[r] = d[r];
  for (I=r-1; I>=0; I--) {
    cdf[I] = LOGL_SUM(cdf[I+1], d[I]);
  }

  /* smooth cdf by linear interpolation in logs */
  for (i=r; i>0; i=j) {
    for (j=i-1; j>0 && d[j]==LOGZERO; j--) ;	/* find next non-zero p */ 
    if (i!=j) slope = (cdf[i]-cdf[j])/(i-j);	/* slope */
    for (k=j+1; k<i; k++) cdf[k] = cdf[j] + (k-j)*slope;
  }

  return cdf;
} /* cdf */

/******************************************************************************/
/*
	get_scaled_llr

	Given a set of observed counts (o_i) and expected counts (e_i)
	and a scale factor alpha, return the scaled information content: 
		Sum_i=1..A NINT(alpha * o_i log(o_i/e_i))
*/
/******************************************************************************/
int get_scaled_llr(
  double alpha,			/* scale factor */
  int A,			/* size of alphabet */
  double *o,			/* observed counts */
  double *e 			/* expected counts */
) {
  int i;
  int llr = 0;

  for (i=0; i<A; i++) {
    if (o[i]) {
      llr += NINT(alpha * o[i] * log(o[i]/e[i]));
    }
  }
  return llr >= 0 ? llr : 0;
} /* get_scaled_llr */

/******************************************************************************/
/*
	brute

	Calculate the distribution of scaled LLR using a brute force
	algorithm (for small values of n only).  Used for checking
	get_llr.
*/
/******************************************************************************/
void brute1(
  int A,				/* dimension of discrete distribution */
  double *dd,				/* discrete distribution */
  int N,				/* number of samples of distribution */
  double alpha,				/* scale factor for llr */
  double *o,
  double *e,
  double *p,
  double prob
){
  int i;
  int llr;

  if (N==0) {
    llr = get_scaled_llr(alpha, A, o, e);
    p[llr] += prob;
  } else {
    for (i=0; i<A; i++) {
      o[i]++;
      brute1(A, dd, N-1, alpha, o, e, p, prob*dd[i]);
      o[i]--;
    }
  } 
} /* brute1 */

double *brute(
  int A,				/* dimension of discrete distribution */
  double *dd,				/* discrete distribution */
  int N,				/* number of samples of distribution */
  int range,				/* desired range for scaled LLR */
  double alpha  			/* scale factor for llr */
){
  int i;
  double *o=NULL;
  double *e=NULL;
  double *p=NULL;
  Resize(o, A, double);
  Resize(e, A, double);
  Resize(p, range+2, double);

  for (i=0; i<A; i++) { 
    o[i] = 0;
    e[i] = dd[i] * N;
  }
  for (i=0; i<range+2; i++) { 
    p[i] = 0;
  }

  brute1(A, dd, N, alpha, o, e, p, 1.0);
 
  return p;
} /* brute */

/******************************************************************************/
/*
	get_llr_pv

	Get the log p-value of a weighted log-likelihood ratio given:
		llr	log-likelihood ratio
		n	(weighted) number of sequences in the alignment
		w	width of the alignment
		range 	desired number of values in range
		frac	speedup factor
		alength	number of letters in alphabet
		dd[]	alphabet frequency distribution

	Computes the distributions on the fly only as needed.
	Range should be the same each time get_llr_pv is called.
	Call with integral n and w=0 to allocate table in parallel mode.

	Returns the p-value.
*/
/******************************************************************************/
extern double get_llr_pv(
  double llr,				/* log likelihood ratio */
  double n,				/* wgtd number sequences in alignment */
  int w,				/* width of alignment */
  int range,				/* desired range for resolution */
  double frac,				/* speedup factor */
  int alength,				/* length of alphabet */
  double *dd 				/* alphabet frequency distribution */
)
{
  int i, N;
  double I;				/* weighted log likelihood ratio */
  int I0, I1;				/* position of llr in table */
  double logpv;				/* log pvalue */
  double n0, n1;			/* floor and ceil of n */
  double alpha;				/* scale factor used */

  if (n<=1) return 0.0;			/* only one site p-value = 1.0 */

  /* return geometric mean if N is not integral */
  if ( (n0=floor(n)) != (n1=ceil(n)) )
    return ( 
      (n1-n)*get_llr_pv(llr, n0, w, range, frac, alength, dd) +
      (n-n0)*get_llr_pv(llr, n1, w, range, frac, alength, dd) 
    );

  N = (int) n;				/* make n an integer */

  /* N larger than any previous N? */
  if (ndistrs < N) {			/* first call */
    Resize(distrs, N+1, DISTR);		/* create array of distributions */
    for (i=ndistrs+1; i<=N; i++) {
      distrs[i].w = 0;
      distrs[i].offset = NULL;
      distrs[i].range = NULL;
      distrs[i].d = NULL;
      distrs[i].cdf = NULL;
      distrs[i].mean = 0;
    }
    ndistrs = N;				/* set maximum N */
  }

  /* done if w == 0 */
  if (w == 0) return 0.0;

  /* w larger than any previous w for this N? */
  if (distrs[N].w < w) {			/* larger w */
    Resize(distrs[N].d, w+1, double *);
    Resize(distrs[N].cdf, w+1, double *);
    Resize(distrs[N].offset, w+1, int);
    Resize(distrs[N].range, w+1, int);

    /* first time? */ 
    if (distrs[N].w == 0) {			/* get the w=1 distribution */
      distrs[N].d[1] = llr_distr(alength, dd, N, range, frac,
        &distrs[N].alpha, &distrs[N].offset[1], &distrs[N].range[1]);
      /* get mean of LLR for w = 1 */
      for (i=0; i<=distrs[N].range[1]; i++) {
        double llr = (i + distrs[N].offset[1]) / distrs[N].alpha;
        distrs[N].mean += exp(distrs[N].d[1][i])*llr;
      }
      distrs[N].cdf[1] = cdf(distrs[N].d[1], distrs[N].range[1]);
      distrs[N].w = 1;
    } /* first time */

    /* get the distributions for widths oldw .. maxw */
    /*fprintf(stderr, "enter cdf N= %d w= %d oldw= %d\n", N, w, distrs[N].w);*/
    for (i=distrs[N].w+1; i<=w; i++) {		/* width */
      distrs[N].d[i] = sum_distr(
        distrs[N].d[i-1], 
        distrs[N].range[i-1], 
        distrs[N].d[1], 
        distrs[N].range[1], 
        &distrs[N].range[i]
      );
      distrs[N].offset[i] = distrs[N].offset[i-1] + distrs[N].offset[1];
      distrs[N].cdf[i] = cdf(distrs[N].d[i], distrs[N].range[i]);
    } /* width */
    /*fprintf(stderr, "leave cdf\n");*/
    distrs[N].w = w; 				/* set maximum w */
  } /* new w */
  
  /* get position in table */
  alpha = distrs[N].alpha;
  I = (alpha * llr) - distrs[N].offset[w];		/* position in table */
  I0 = (int) I;						/* floor of position */
  I1 = I0 + 1;						/* ceil. of position */
  if (I < 0) {						/* lower bound I */
    logpv = distrs[N].cdf[w][0];
  } else if (I0 >= distrs[N].range[w]) {		/* upper bound I */
    logpv = distrs[N].cdf[w][distrs[N].range[w]];
  } else {						/* lin. interpolate */
    logpv =
      distrs[N].cdf[w][I0] + (I-I0)*(distrs[N].cdf[w][I1]-distrs[N].cdf[w][I0]);
  }

  /* return log p-value */
  return logpv;
} /* get_llr_pv */

/******************************************************************************/
/*
	get_llr_mean

	Returns the mean of the LLR distribution for w=1.
	Assumes get_llr_pv has already been called.
*/
/******************************************************************************/
extern double get_llr_mean(
  double n 				/* wgt. number sequences in alignment */
)
{
  return(distrs[NINT(n)].mean);
} /* get_llr_mean */


#ifdef SO
/************************************************************************/
/*
	llr <alength> <freq_file> <N>

	Compute the probability distribution for the scaled information
	content (LLR) of N letters.  Only the most probable <frac> LLR values
	are used to speed the calculation.
*/
/************************************************************************/
extern int main(
  int argc,
  char** argv
)
{
  int i, j, I, n;
  int A=0;				/* length of alphabet */
  char *ffreq=NULL;			/* name of frequency file */
  FILE *ffile;				/* frequency file */
  double *dd = NULL;			/* discrete distribution */
  char *alphabet = NULL;		/* alphabet */
  int N=2, minN=0, maxN=0;		/* number of observations */
  int w=1;				/* maximum width of alignment */
  double frac = 1;			/* fraction of scores to use */
  int range = 100;			/* desired range per N */
  double *p2=NULL;			/* LLR statistic distribution */
  double *pv2=NULL;
  int brute_maxN = 0;			/* maximum N for checking llr_distr */
  char *line;				/* input buffer line */
  double alpha;				/* scale factor used */

  i = 1;
  argv[0] = "llr";
  DO_STANDARD_COMMAND_LINE(2,
    USAGE(<alength> <freq_file> [options]);
    USAGE(\n\t<alength>\tnumber of letters in alphabet); 
    USAGE(\t<freq_file>\tfile containing alphabet and letter frequencies);
    USAGE(\t\t\twhere each line is: <letter> <freq>);
    NON_SWITCH(1,\r,
      switch (i++) {
        case 1: A = atoi(_OPTION_); break;
        case 2: ffreq = _OPTION_; break;
        default: COMMAND_LINE_ERROR;
      }
    );
    DATA_OPTN(1, N, <N>, \tnumber of observations; default=2,
      N = atof(_OPTION_));
    DATA_OPTN(1, minN, <minN>, \tminimum number of observations; default=2,
      minN = atoi(_OPTION_));
    DATA_OPTN(1, maxN, <maxN>, \tmaximum number of observations; default=2,
      maxN = atoi(_OPTION_));
    DATA_OPTN(1, w, <w>, \tmaximum width of alignment; default=1,
      w = atoi(_OPTION_));
    DATA_OPTN(1, frac, <frac>, \tfraction of possible scores to use; default=1,
      frac = atoi(_OPTION_));
    DATA_OPTN(1, range, <range>, scale llr to have <range> values; default=100,
      range = atoi(_OPTION_));
    USAGE(\n\tCompute the probability distribution for the log-likelihood);
    USAGE(\tratio (LLR) of N letters.  If a range for N is given);
    USAGE(\tthe distributions for each value in the range are computed.);
    USAGE(\tIf <frac> is given only the most probable <frac> fraction of);
    USAGE(\tvalues are used to speed the calculation.);
    USAGE(\n\tCopyright);
    USAGE(\t(2000-2006) The Regents of the University of California);
    USAGE(\tAll Rights Reserved.);
    USAGE(\tAuthor: Timothy L. Bailey);
  );

  init_log();
  init_exp();

  /* set up range for N */
  if (minN==0) minN = N;
  if (maxN==0) maxN = N;
  if (minN > maxN) {
    fprintf(stderr, "You must specify -maxN larger than -minN.\n");
    exit(1);
  }

  /* get alphabet and frequencies */
  ffile = fopen(ffreq, "r"); 
  Resize(dd, A, double);
  Resize(alphabet, A, char);
  i = 0;
  while ((line = getline2(ffile))) {
    char c;
    double f;
    if (line[0] == '#') continue;		/* skip comments */
    if (sscanf(line, "%c %lf", &c, &f) != 2) break;
    alphabet[i] = c;
    dd[i++] = f;
    myfree(line);
  }
  if (i != A) {
    fprintf(stderr, "%d frequencies were found; %d were expected.\n", i,A);
    for (j=0; j<i; j++) printf("%c %f\n", alphabet[j], dd[j]);
    exit(1);
  }
  printf("# Alphabet frequency distribution:\n");
  for (j=0; j<A; j++) printf("# %c %f\n", alphabet[j], dd[j]);

  /* get distribution and print it */
  for (n=minN; n<=maxN; n++) { 			/* number of observations */

    printf("# desired range %d\n", range);

    /* get the distribution set */
    get_llr_pv(0, (double) n, w, range, frac, A, d);
    printf("# alpha %f\n", alpha);

    /* get the check distribution */
    if (n <= brute_maxN) { 
      p2 = brute(A, dd, n, range, alpha); 
      pv2 = cdf(p2, range);
    }

    /* print results */
    printf("# N    I    llr         p     1-cdf\n");
    for (I=0; I<=distrs[n].range[w]; I++) {		/* LLR */
      double m1, e1, m2, e2;
      if (distrs[n].d[w][I] == LOGZERO) {
        m1 = e1 = 0;
      } else {
        exp10_logx(distrs[n].d[w][I]/log(10.0), m1, e1, 1);
      }
      if (distrs[n].cdf[w][I] == LOGZERO) {
        m2 = e2 = 0;
      } else {
        exp10_logx(distrs[n].cdf[w][I]/log(10.0), m2, e2, 1);
      }
      printf("%3d %3d %5.1f %3.1fe%+05.0f %3.1fe%+05.0f\n",
        n, I, (distrs[n].offset[w]+I)/distrs[n].alpha, m1, e1, m2, e2);
    } /* LLR */
  } /* n */

  return 0;
}
#endif /* SO */
