/* ---------------------------------------------------------------------- */
/*  E_value                                                               */
/*                                                                        */
/*  Calculated the entropy ala MEME                                       */
/*                                                                        */
/*  This code was lifted directly from MEME v. 3.5.4, lightly modified    */
/*  and consolidated for understanding and modularity                     */
/*                                                                        */
/* ---------------------------------------------------------------------- */
/*                                                                        */
/*    MEME++                                                              */
/*    Author: Timothy L. Bailey                                           */
/*                                                                        */
/*    Copyright  (c)  1994-2006  The  Regents  of the  University of      */ 
/*    California.  All  Rights  Reserved.                                 */                           
/*                                                                        */                                                                  
/*    Permission  to use,  copy,  modify,  and  distribute  any part      */
/*    of this  software for  educational,  research  and  non-profit      */ 
/*    purposes,  without  fee,  and  without a written  agreement is      */ 
/*    hereby  granted,  provided  that the  above  copyright notice,      */ 
/*    this paragraph  and the following  three  paragraphs appear in      */ 
/*    all copies.                                                         */                                                   
/*                                                                        */                                                                  
/*    Those  desiring to  incorporate this  software into commercial      */
/*    products  or use for  commercial  purposes  should contact the      */     
/*    Technology  Transfer  Office,  University of California,   San      */ 
/*    Diego,  9500 Gilman Drive,  La Jolla,  California, 92093-0910,      */         
/*    Phone: (858) 534-5815.                                              */                                        
/*                                                                        */                                                                  
/*    IN  NO  EVENT  SHALL THE  UNIVERSITY  OF CALIFORNIA  BE LIABLE      */ 
/*    TO  ANY  PARTY FOR  DIRECT,  INDIRECT, SPECIAL, INCIDENTAL, OR      */     
/*    CONSEQUENTIAL  DAMAGES,  INCLUDING  LOST PROFITS, ARISING  OUT      */ 
/*    OF  THE  USE  OF  THIS  SOFTWARE,  EVEN  IF THE UNIVERSITY  OF      */ 
/*    CALIFORNIA  HAS  BEEN  ADVISED  OF  THE  POSSIBILITY  OF  SUCH      */ 
/*    DAMAGE.                                                             */                                                       
/*                                                                        */                                                                  
/*    THE SOFTWARE  PROVIDED HEREUNDER IS ON N  "AS IS" BASIS,  AND       */ 
/*    THE  UNIVERSITY OF CALIFORNIA  HAS  NO OBLIGATIONS  TO PROVIDE      */         
/*    MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.      */  
/*    THE UNIVERSITY  OF CALIFORNIA  MAKES  NO  REPRESENTATIONS  AND      */     
/*    EXTENDS  NO  WARRANTIES  OF  ANY  KIND,  EITHER  EXPRESSED  OR      */ 
/*    IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES      */
/*    OF  MERCHANTABILITY  OR  FITNESS FOR A  PARTICULAR PURPOSE, OR      */ 
/*    THAT  THE USE  OF THE MATERIAL  WILL NOT  INFRINGE ANY PATENT,      */         
/*    TRADEMARK OR OTHER RIGHTS.                                          */                                    
/*                                                                        */
/* ---------------------------------------------------------------------- */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include <time.h>
#include "evalue_meme.h"
#include "alloc.h"

static BOOLEAN init_log_done = FALSE;
static BOOLEAN init_exp_done = FALSE;

static double log_table[2 * (int) log_precision + 2];
static double exp_table[(int) BITS * (int) exp_precision + 2];

static double _logl_x;
static double _expl_x;


int int_compare(const void *v1, const void *v2);

double *llr_distr(int A, double *dd, int N, int desired_range, double frac,
                  double *alpha, int *offset, int *range);

double *sum_distr(double *d1, int r1, double *d2, int r2, int *r_sum);

double get_llr_pv(double llr, double n, int w, int range, double frac,
                  int alength, double *dd);

double log_comb(int m, int n);

double get_log_nalign(MTYPE mtype, int w, int N, BOOLEAN invcomp,
                       int Num_seq, int *seqLen);

double get_log_sig(double score, MTYPE mtype, int w, double wN, int N,
                   BOOLEAN invcomp, BOOLEAN pal, int Num_seq, int *seqLen);

double log_qfast(int n, double logk);


/* ---------------------------------------------------------------------- */
/* From MEME: llr.c                                                       */
/* ---------------------------------------------------------------------- */


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
double get_llr_pv(
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


/* ---------------------------------------------------------------------- */
/* From MEME: likelihood.c                                                */
/* ---------------------------------------------------------------------- */


/**********************************************************************/
/*
        log_comb
 
        Compute logarithm of m choose n.

*/
/**********************************************************************/
 
double log_comb(
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


/**********************************************************************/
/*
	get_log_nalign

	Get an upper bound on the number of independent alignments
	of segments of length w.
*/
/**********************************************************************/
double get_log_nalign(
  MTYPE mtype,					/* type of model */
  int w,	 				/* width of motif */
  int N,					/* number of occurrences */
  BOOLEAN invcomp,				/* inv. compl. seq allowed */
  int Num_seq,  				/* the dataset */
  int* seqLen
)
{
  int i, t;
  int nseqs = Num_seq;	                	/* number of sequences */
  double log_nalign = 0;			/* log number alignments */
  int icfactor = invcomp ? 2 : 1;		/* double the possible sites */
  int *len=NULL;

  /* 
    sort the sequence lengths in decreasing order in array len[] first time thru
  */
  if (len == NULL) { 				/* first time through */
    Resize(len, nseqs, int);
    for (i=0; i< nseqs; i++) len[i] = seqLen[i];
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
        get_log_sig

        Calculate the statistical significance of the alignment given
        its score and the type of objective function in force.

        If N>0, returns log E-value.
        If N==0, returns log p-value.
*/
/**********************************************************************/
double get_log_sig(
  double score,                                 /* score of alignment */
  MTYPE mtype,                                  /* type of model */
  int w,                                        /* width of motif */
  double wN,                                    /* weighted number of sites */
  int N,                                        /* number of sites */
  BOOLEAN invcomp,                              /* inv. compl. strand, too */
  BOOLEAN pal,                                  /* motif is DNA palindrome */
  int Num_seq,                                  /* the dataset */
  int *seqLen
)
{
  double log_pv;                        /* log of p-value */
  double log_sig = 0;                   /* return value */

  log_pv = log_qfast(w, -score);        /* p-value of product of p-values */

  if (N) {                              /* use E-value of alignment */
    log_sig = log_pv + get_log_nalign(mtype, w, N, invcomp && !pal, Num_seq, seqLen);
  } else {                              /* use p-value of alignment */
    log_sig = log_pv;
  }

  return(log_sig);
} /* get_log_sig */




/**********************************************************************/
/*
	log_qfast
	
	Calculate the log p-value of the log of the 
	product of uniform [0,1] random variables.

*/
/**********************************************************************/
double log_qfast(
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
int int_compare(const void *v1, const void *v2)
{
  const int * s1 = (const int *) v1;
  const int * s2 = (const int *) v2;

  return(*s2 - *s1);
} /* int_compare */


void init_log(void);
void init_exp(void);

double E_value(double **obs, int numsites, double *back, int motifLen, 
             int Numseq, int *seqLen) {

  int i, j;
  double *rentropy;                             /* IC of each column */
  double ent = 0;                               /* entropy per column */
  int alength = 4;                              /* length of alphabet */
  int  w = motifLen;                            /* width of motif */
  int N = numsites;                             /* number of sites */
  double log_pop, logev, m, e, totalllr;        /* log product of col p-value */
  MTYPE mtype = Tcm;                            /* oops, zoops, anr            */
  BOOLEAN invcomp = 1;                          /* invcomp=1 use inverse complement DNA strand, too */
  BOOLEAN pal = 0;                              /* pal=1, palindrome*/

  rentropy = alloc_double(w);

  init_log();
  init_exp();
  
 /* calculate the relative entropy of each column in motif */
  totalllr = log_pop = 0;
  for (i = 0; i < w; i++) {
    double llr;                         /* position */
    
    rentropy[i] = 0.0;
    
    for (j=0; j<alength; j++) {         /* alphabet letter */
      double f = obs[i][j];             /* motif freq */
      double p = back[j];               /* background freq */
      
      ent += f ? f * LOG2(f) : 0;       /* entropy */
      rentropy[i] += (f && p) ? f * LOG(f/p) : 0;
      
    }                                   /* alphabet letter */
    
    llr = N * rentropy[i];              /* log likelihood ratio */
    RND(llr, RNDDIG, llr);              /* round to RNDDIG places */
    totalllr += llr;            
    log_pop += get_llr_pv(llr, N, 1, LLR_RANGE, 1.0, alength, back);
    rentropy[i] /= LOG(2);
    
  }
  /* compute the log E-value of the motif */
  
  logev = get_log_sig(-log_pop, mtype, w, N, N, invcomp, pal, Numseq, seqLen);
  exp10_logx((logev)/log(10.0), m, e, 1);
  
  //printf("Num of binding sites=%d\tmotif width=%d\t",N, w);
  //printf("llr=%.0f\tE-value=%3.1fe%+04.0f\n",totalllr,m, e);
  //fprintf(fp, "Num of binding sites=%d\tmotif width=%d\tllr=%.0f\tE_value=%3.1fe%+04.0f\n", N, w, totalllr, m, e);

  if (rentropy) { free(rentropy); rentropy=NULL; }

  //*mantissa = m; *exponent = e;
  return (logev);
  
} /* E_value (MEME: calc_entropy) */
  
  

/* ---------------------------------------------------------------------- */
/* From MEME: logs.c                                                      */
/* ---------------------------------------------------------------------- */


/**********************************************************************/
/*
	init_log

	Setup lookup table for log(x), 0 < x <= 2
*/
/**********************************************************************/
void init_log(void)
{
  int i;
  double x;

  for (i=0; !init_log_done && i <= 2 * log_precision + 1; i++) {
    x = (double) i/log_precision;
    log_table[i] = LOG(x);
  }

} /* init_log */

/**********************************************************************/
/*
	init_exp

	Setup lookup table for exp(x), -BITS <= x < 0 
*/
/**********************************************************************/
void init_exp(void)
{
  int i;
  double x;

  for (i=0; !init_exp_done && i <= BITS * exp_precision + 1; i++) {
    x = -i/exp_precision;
    exp_table[i] = exp(x);
  }

} /* init_exp */




