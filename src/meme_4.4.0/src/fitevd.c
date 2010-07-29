/************************************************************************/
/*									*
* 	fitevd								*
* 	Author: Timothy L. Bailey					*
*									*
* 	Copyright							*
* 	(2001-2002) Timothy L. Bailey					*
*									*
* 	All rights reserved.						* 
*       Permission to use, copy, modify, and distribute any part of     *
*       this software for educational, research and non-profit purposes,*
*       without fee, and without a written agreement is hereby granted, *
*       provided that the above copyright notice, this paragraph and    *
*       the following three paragraphs appear in all copies.            *
*									*
*	Those desiring to incorporate this software into commercial	*
*	products or use for commercial purposes should contact the 	*
*	AUTHOR.
*                                                                       *
*       IN NO EVENT SHALL THE AUTHOR BE LIABLE TO     			*
*       ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR         *
*       CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF   *
*       THE USE OF THIS SOFTWARE, EVEN IF THE AUTHOR 			*
*       HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.             *
*                                                                       *
*       THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE *
*       AUTHOR HAS NO OBLIGATIONS TO PROVIDE    		      	*
*       MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  *
*       THE AUTHOR MAKES NO REPRESENTATIONS AND       			*
*       EXTENDS NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED, *
*       INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF        *
*       MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT    *
*       THE USE OF THE MATERIAL WILL NOT INFRINGE ANY PATENT,           *
*       TRADEMARK OR OTHER RIGHTS.                                      *
*                                                                       *
************************************************************************/

/******************************************************************************/
/*

	Routines for fitting an extreme value distribution or an exponential
	distribution to a set of score-length pairs.

	EVD:
	Routines for fitting the parameters of an extreme value distribution
	(EVD) to a set, X, of score-length pairs.  X = {x_i : x_i = (s_i, t_i)},
 	where s_i is a sequence similarity score, and t_i is the length of 
	the target sequence.  The length of the query sequence, is q, and is 
	assumed constant for all x_i.  The entry point is subroutine
	fit_score_distribution().  You should include fitevd.h wherever that 
	subroutine is called.  If maxscores is specified to 
   	fit_score_distribution(), the given score pairs are randomly 
	subsampled, and approximately maxscores of them are used in computing 
	the EVD.

	EXP: see fitexp

	If the macro MAIN is defined at compile time, a stand-alone 
	program is created that reads the (score length [<id>]) tuples from 
	standard input, estimates the EVD parameters, and prints them to 
	standard output.  It then prints the estimated p-values for all the
	scores, sorted in increasing order. 

	The subroutine main() thus serves as an example of how to call ml_evd().
*/
/******************************************************************************/

#include <string.h>
#include "fitevd.h"

/* turn this on for debugging */
static int TRACE = 0;				/* turn off debug trace */

#define PI 3.141569

/******************************************************************************/
/*
	local prototypes 
*/
/******************************************************************************/

static void init_evd(EVD *evd) {
  evd->lambda = 0.0; // parameters
  evd->K = 0.0;
  evd->H = 0.0;
  evd->n = 0;				 // number of scores used
  // Other stuff for extreme values.
  evd->g = 0;				 // length group number
  evd->min_t = 0.0;  // min_t < t (or gc) <= max_t
  evd->max_t = 0.0;
  evd->mid_t = 0.0;		 // midpoint of length range
  // Exponential parameters.
  evd->mu1 = 0.0;			// Exponential mean.
  evd->mu2 = 0.0;			// Gaussian mean.
  evd->sigma2 = 0.0;	// Gaussian sd.	
  evd->c = 0.0;				// Mixing parameter.
  // Log likelihood (applies to both types of distributions).
  evd->L = 0.0;
}

static void init_score_set(SCORE_SET *score_set) {
  score_set->n = 0;         // number of scores in set
  score_set->scores = NULL; // array containing X = {(s_i,t_i)}
  score_set->max_t = 0;     // maximum target length (t_i)
  score_set->q = 0;         // length of query
  score_set->negatives_only = FALSE; // All negatives?
  score_set->max_gap = 0;		// Maximum gap allowed.
  score_set->total_length = 0.0; // Total size of database.
  score_set->min_e = 0.0;	// Mininum E-value of scores.
  score_set->avgw = 0.0;	// Average motif width.
  score_set->phit = 0.0;  // Probability of a false hit.
  score_set->wimp = 0.0;  // Probability of false hit within max_gap.
  score_set->egap = 0.0;  // E[gap length]
  score_set->ehits = 0.0; // E[hits/match]
  score_set->espan = 0.0; // E[length of match]
  score_set->e_hit_score = 0.0; // E[score of single hit]
  score_set->avg_min_pvalue = 0.0; // Geometric mean of min p-values for all motifs.
}

static int p_compare(
  const void *v1,
  const void *v2
);

static EVD_SET fitexp(
  SCORE_SET score_set,				/* the set of scores and lens */
  int maxiter,					/* maximum iterations */
  double eps1,					// Stopping criterion for EM.
  int min_scores,				// Minimum number of scores/range.
  double max_ratio				// Maximum GC ratio between adj. ranges
);

static EVD_SET ml_evd(
  SCORE_SET score_set,				/* the set of scores and lens */
  double H,					/* initial estimate for H */
  int maxiter1,					/* maximum iterations for ML */
  int maxiter2,					/* maximum iterations for N-R */
  double eps1,					/* error tolerance for ML */
  double eps2, 					/* error tolerance for N-R */
  int maxscores,				/* randomly select scores */
  int size,					/* size of length ranges */
  double lspan,					/* minimum length ratio/group */
  double ethresh				/* outlier threshold */
);

EVD ml_evd_nr(
  EVD evd,					/* initial guess */
  SCORE_SET score_set,				/* the set of scores and lens */
  int maxiter1,					/* maximum iterations for ML */
  int maxiter2,					/* maximum iterations for N-R */
  double eps1,					/* error tolerance for ML */
  double eps2, 					/* error tolerance for N-R */
  double ethresh				/* outlier threshold */
);

EVD get_est (
  SCORE_SET score_set				/* set of scores */
);

EVD ml_lambda(
  EVD evd,					/* current parameters of EVD */ 
  SCORE_SET score_set,				/* set of scores */
  double eps,					/* stopping criterion */
  int maxiter					/* max. allowed iterations */
);

EVD ml_H(
  EVD evd,					/* current parameters of EVD */ 
  SCORE_SET score_set,				/* set of scores */
  double eps,					/* stopping criterion */
  int maxiter					/* max. allowed iterations */
);

DERIV dlambda(
  EVD evd,					/* current parameters of EVD */ 
  SCORE_SET score_set				/* set of scores */
);

DERIV dH(
  EVD evd,					/* current parameters of EVD */ 
  SCORE_SET score_set				/* set of scores */
);

#ifdef MAIN
SCORE_SET read_scores ();
#endif /* MAIN */

/******************************************************************************/
/*
	subroutines	
*/
/******************************************************************************/

/**********************************************************************/
/*
	set_up_score_set

	Initialize a score_set to contain scores and lengths
	for use in computing the score distribution.
*/
/**********************************************************************/
extern SCORE_SET *set_up_score_set(
  PROB_T p_threshold,		// P-value threshold for motif hits.
  PROB_T dp_threshold, 		// Score threshold for DP.
  int    max_gap, 		// Maximum gap length to allow in matches.
  BOOLEAN_T negatives_only,	// Negatives only in score set?
  MHMM_T* the_hmm 		// The HMM itself.
)
{
  int j;
  int i_state;
  int n_motifs = the_hmm->num_motifs;	// number of motifs in model
  SCORE_SET *score_set = NULL;		// Scores for calculating distribution.
  BOOLEAN_T need_max_gap = max_gap < 0;	// Need to estimate max_gap?

  //
  // Initialize everything.
  //
  mm_resize(score_set, 1, SCORE_SET);	// Create score set.
  score_set->scores = NULL;		// array of scores & lengths
  score_set->n = 0;			// number of scores saved
  score_set->max_t = 0;			// length of longest sequence
  score_set->q = 1;			// assume query has length 1
  score_set->negatives_only = negatives_only;
  score_set->total_length = 0;		// Total size of database.
  score_set->phit = 0;			// Probability of false hit.
  score_set->wimp = 0;			// Probability of false hit within max-gap.
  score_set->ehits = 0;			// E[hits/match]
  score_set->espan = 0;			// E[length of match]
  score_set->avg_min_pvalue = 0;	// Geometric average of minimum p-values.
  if (need_max_gap) max_gap = 0;	// Going to estimate max_gap.

  // 
  // Here we calculate the wimp factor (probability of
  // a hit within max-gap positions), expected number of
  // hits per match and the expected length (span) of a match.
  //
  if (p_threshold > 0 && dp_threshold > 0) {
    double avg_width = 0;

    // 
    // Get the number of motifs and average motif width.
    //
    for (i_state=0; i_state<the_hmm->num_states; i_state++) {       // state
      MHMM_STATE_T* state = the_hmm->states + i_state;  // Current state.
      if (state->type == START_MOTIF_STATE) {
        avg_width += state->w_motif;
        score_set->avg_min_pvalue += log(state->min_pvalue);

      } else if (need_max_gap && state->type == SPACER_STATE) {
	double in_cost = BIG;		// (MIN) In-transition cost.
	double out_cost = BIG;		// (MIN) In-transition cost.
	double loop_cost = 1;		// Self-loop cost.
	int max_loops;			// Maximum number of loops.

        // Find the cheapest in-transition cost and the self-loop cost.
        for (j=0; j<state->ntrans_in; j++) {	// in transition
          int j_state = state->itrans_in[j];
          if (j_state == 0) continue;		// ignore start state
          if (i_state==j_state) {		// self loop
            loop_cost = -get_array_item(j, state->trans_in); 
          } else {
            in_cost = MIN(in_cost, -get_array_item(j, state->trans_in));
          }
        }

        // Find the cheapest out-transition cost.
        for (j=0; j<state->ntrans_out; j++) {	// out transition
          int j_state = state->itrans_out[j];
          if (j_state == 0) continue;		// ignore start state
          if (i_state != j_state) {
            out_cost = MIN(out_cost, -get_array_item(j, state->trans_out));
          }
        }

 	// Use minimum values of transition costs to estimate the
        // maximum gap allowed by this spacer.
        max_loops = floor((dp_threshold - in_cost - out_cost)/loop_cost);
        max_gap += (max_loops > 0) ? 1 + max_loops : 0;
      }
    } // state
    
    // Use the average maximum gap as max_gap if we had to esitmate it.
    if (need_max_gap) {
      max_gap /= the_hmm->num_spacers;		// Estimate maximum gap.
    }
   score_set->max_gap = max_gap;		// Record max_gap.

    // Mean width of motifs.
    avg_width /= n_motifs; 
 
    // Geometric mean of minimum p-values for motifs.
    score_set->avg_min_pvalue = exp(score_set->avg_min_pvalue/n_motifs);

    //
    // Estimate phit, wimp, egap, ehits, espan.
    //
    if (the_hmm->type == STAR_HMM || the_hmm->type == COMPLETE_HMM) {
      int i;
      double cum_p, egap, ehits, espan;
      double q = EV(p_threshold, max_gap*n_motifs);	// Pr(hit within max_gap)
      double phit = EV(p_threshold, n_motifs);		// Pr(hit)
      double e_hit_score = 1.442;			// E[hit score] = 1/log(2).

      //
      // Compute statistics assuming no gap costs.
      //
      for (i=egap=0, cum_p=phit; i<=max_gap; i++) {
	egap += i * cum_p;
        //fprintf(stderr, "i %d cum_p = %g egap = %g\n", i, cum_p, egap);
	cum_p *= (1-phit);
      }
      egap /= q;
      ehits = (q < 1) ? q/(1-q) + 1 : 1e6;
      espan = (ehits * avg_width) + ((ehits-1) * egap);
      score_set->avgw = avg_width;
      score_set->phit = phit;
      score_set->wimp = q;
      score_set->egap = egap;
      score_set->ehits = ehits;
      score_set->espan = espan;
      score_set->e_hit_score = e_hit_score;

    } else {					// LINEAR_HMM
      // FIXME: need to figure out E[span] for LINEAR
      // models but I haven't thought about it.
    }
 
    //
    // Check that wimp is not too big.
    //
    if (score_set->wimp > MAX_WIMP) {
      die("The probability of a false hit is too high (%.2f > %.2f).\nReduce -p-thresh and/or max-gap.", score_set->wimp, MAX_WIMP);
    } 

    if (TRACE) {
      fprintf(stderr, "p %f g %d m %d w %f wimp %f ehits %f espan %f\n", 
	p_threshold, max_gap, n_motifs, avg_width,
	score_set->wimp, score_set->ehits, score_set->espan);
    }

  } // repeated match with p-value scoring

  return(score_set);
} // set_up_score_set

/**********************************************************************/
/*
        s_compare
 
        Compare two scores in descending order.  Return <0 >0
        if the first is >, < the second .  If they are equal,
        resolves ties by returning <0 if the second has smaller address.
*/
/**********************************************************************/
int s_compare(
  const void *v1,
  const void *v2
)
{
  const SCORE *s1 = (const SCORE *) v1;
  const SCORE *s2 = (const SCORE *) v2;
  double diff = s1->s- s2->s;
  if (diff == 0) diff = (double) (s2 - s1);
  return ((diff > 0) ? -1 : ( (diff < 0) ? 1 : 0) );
} /* s_compare */


/**********************************************************************/
/*
        p_compare
 
        Compare two p-values in ascending order.  Return <0 >0
        if the second is >, < the first.  If they are equal,
        resolves ties by returning <0 if the second has smaller address.
*/
/**********************************************************************/
static int p_compare(
  const void *v1,
  const void *v2
)
{
  const SCORE *s1 = (const SCORE *) v1;
  const SCORE *s2 = (const SCORE *) v2;
  double diff = s2->pv - s1->pv;
  if (diff == 0) diff = (double) (s1 - s2);
  return ((diff > 0) ? -1 : ( (diff < 0) ? 1 : 0) );
} /* p_compare */

/**********************************************************************/
/*
        gc_compare
 
        Compare two score gc's in ascending order.  Return <0 >0
        if the second is >, < the first.  If they are equal,
        resolves ties by returning <0 if the second has smaller address.
*/
/**********************************************************************/
int gc_compare(
  const void *v1,
  const void *v2
)
{
  const SCORE *s1 = (const SCORE *) v1;
  const SCORE *s2 = (const SCORE *) v2;
  double diff = s2->gc - s1->gc;
  if (diff == 0) diff = (double) (s1 - s2);
  return ((diff > 0) ? -1 : ( (diff < 0) ? 1 : 0) );
} /* gc_compare */

/******************************************************************************/
/*
	em_exp

	Fit an exponential distribution to scores using
	expectation maximization (EM).  The scores are fit
	to a mixture of an exponential and a Gaussian distribution.
	The exponential distribution fits the low scores, and
	the Gaussian fits the large outliers.
	The distribution function is:
		pv = exp(-x/mu1)

	Returns evd.
*/
/******************************************************************************/
static EVD em_exp(
  SCORE_SET score_set,				/* the set of scores and GC contents */
  int maxiter,					/* maximum iterations */
  double eps1					// Stopping criterion for EM.
)
{
  int i, iter;
  int n = score_set.n;
  SCORE *scores = score_set.scores;		// Scores.
  double L, old_L;				// Log-likelihood.
  double deltaL;				// Change in log-likelihood.
  double x_max;					// Maximum score.
  double x_sum;					// Sum of scores.
  double x_mean;				// Average score.
  double x_sd;					// Score sd.
  double e_max; 				// Expected largest score.
  int n_big;					// Number scores > e_max.
  double sx;					// Sum of (valid) scores.
  double c, mu1, mu2, sigma2;			// Parameters of distributions.
  double *z = NULL;				// Missing information array.
  double x_i, z_i,f1, f2;
  EVD evd;					// Return;

  // Initializatio
  init_evd(&evd);

  // Create missing information array.
  mm_resize(z, n, double);

  // Get sum of scores
  x_max = x_sum = 0;
  scores = score_set.scores;
  for (i=0; i<n; i++, scores++) {	// Scores.
    x_i = scores->s;			// Score.
    if (x_i > x_max) x_max = x_i;
    x_sum += x_i;
  } // Scores.
  evd.n = n;

  //
  // Get mean and standard deviation of scores.
  //
  x_mean = x_sum/n;			// Mean score.
  e_max = -x_mean*log(1.0/n);		// E[max score].
  x_sd = n_big = 0;
  scores = score_set.scores;
  for (i=0; i<n; i++, scores++) {	// Scores.
    x_i = scores->s;			// Score.
    x_sd += (x_i - x_mean) * (x_i - x_mean);
    if (x_i > e_max) n_big++;		// Scores larger than e_max.
  } // Scores.
  x_sd = sqrt(x_sd/(n-1));	// Score sd.

  //
  // Initialize starting point for EM using guesses.
  //
  c = score_set.negatives_only ? 
    1.0 : 				// All scores are negatives.
    (n - n_big - 1) / (double)n;	// Fraction less than e_max.
  mu1 = x_mean;				// Mean of (valid) data.
  mu2 = e_max;				// Expected largest score.
  sigma2 = fabs(x_max-e_max)/2;		// Put largest out 2sd.
  sx = x_mean * n;			// Invariant over iterations.
  if (TRACE) fprintf(stderr, "c = %f mu1 = %f mu2 = %f sigma2 = %f\n", c, mu1, mu2, sigma2);
  
  //
  // EM.
  //
  L = -BIG;
  for (iter=0; !score_set.negatives_only && iter<maxiter; iter++) {
    double sz, sxz, szx2, k1, k2;

    //
    // E-step.
    // 
    old_L = L;
    L = 0;
    sz = sxz = szx2 = 0;
    k1 = 1/(sqrt(2*PI) * sigma2);			// 1/(sqrt(2pi)sigma)
    k2 = 2 * sigma2 * sigma2;				// 2 sigma^2
    scores = score_set.scores;
    for (i=0; i<n; i++, scores++) {			// Scores.
      x_i = scores->s;					// Score.
      f1 = exp(-x_i/mu1)/mu1;				// Exponential.
      f2 = k1 * exp(- (x_i - mu2) * (x_i - mu2) / k2);	// Gaussian.
      z_i = z[i] = (c * f1) / (c * f1 + (1-c) * f2);	// Missing information variable.
      sz += z_i;					// Sum z_i
      sxz += z_i * x_i;					// Sum x_i z_i
      szx2 += (1 - z_i) * (x_i - mu2) * (x_i - mu2);	// Sum (1-z_i)(x_i-mu2)^2
      L += log( c * f1 + (1-c) * f2 );
    } // scores
    // Done if likelihood didn't change.
    deltaL = L - old_L;
    if (deltaL < eps1) break;

    //
    // M-step.
    //
    c = sz/n;					// Mixing parameter.
    mu1 = sxz/sz;				// Exponential mean.
    if (n == sz) {
      mu2 = mu1; 
      sigma2 = 1;
    } else {
      mu2 = (sx - sxz)/(n - sz);		// Gaussian mean.
      sigma2 = sqrt(szx2/(n - sz));		// Gaussian sd.
    }
    if (sigma2 < eps1) break;			// Ridiculously small.
    if (TRACE) fprintf(stderr, "i %3d c %8.3f mu1 %8.3f mu2 %8.3f sigma2 %8.3g L %8.3f deltaL %f\n",
      iter, c, mu1, mu2, sigma2, L, deltaL);
  } // EM

  //
  // Set up rest of evd.
  //
  evd.mu1 = mu1;
  evd.mu2 = mu2;
  evd.sigma2 = sigma2;
  evd.c = c;
  evd.L = L;

  // Free space.
  myfree(z);

  return(evd);
} // em_exp

/******************************************************************************/
/*
	fitexp

	Fit an exponential distribution to score-GC_content pairs using
	expectation maximization (EM).  

	A set of exponential distributions is returned corresponding to
	scores from sequence regions with GC-contents in distinct ranges.
	The GC-content ranges are determined by:
		min_scores	-- minimum number of scores per range
		max_ratio	-- maximum ratio of GC thresholds between ranges
				   (use 1 range if max_ratio > 100)
	Ranges are set up so that each range has at least min_scores in it.
	A new range i+1 is created when range i contains at least min_scores
	and:
		GC_{i+1}/GC_{i} >= max_ratio, if GC_{i} <= 0.5 
		(1-GC_{i})/(1-GC_{i+1}) >= max_ratio, if GC_{i} > 0.5 

	The scores in each GC-content range are fit to a mixture of an 
	exponential and a Gaussian distribution.
	The exponential distribution fits the low scores, and
	the Gaussian fits the large outliers.
	The distribution function is:
		pv = exp(-x/mu1)

	Returns evd_set.  Sets evd_set.n = 0 in case there is an error.
*/
/******************************************************************************/
static EVD_SET fitexp(
  SCORE_SET score_set,				/* the set of scores and lens */
  int maxiter,					/* maximum iterations */
  double eps1,					// Stopping criterion for EM.
  int min_scores,				// Minimum number of scores/range.
  double max_ratio				// Maximum GC ratio between adj. ranges
)
{
  int i;
  int n = score_set.n;				// Number of scores.
  int valid_n;					// # Scores of matches with ok spans.
  int g;					// Current GC range.
  SCORE *scores;				// Scores.
  EVD *evds = NULL;				// Array of EVDs.
  int nranges = max_ratio < 100 ? n/min_scores + 1 : 1;	
						// Maximum number of GC ranges.
  int i_low;					// Start of current GC range.
  double gc_mid;				// Mean of previous GC range.
  double gc_hi;					// Max of previous GC range.
  EVD_SET evd_set;				// Return.
  evd_set.evds = NULL; // Initialize evds

  // In case of error.
  evd_set.n = 0;

  // Using exponential distribution, not EVD.
  evd_set.exponential = TRUE;

  // Create arrays.
  mm_resize(evds, nranges, EVD);

  // Copy valid scores into new array saving the old array.
  scores = score_set.scores;
  n = score_set.n;
  score_set.scores = NULL;
  for (i=valid_n=0; i<n; i++, scores++) {
    // Span is too long if no room for extending match on each end.
    if ( (scores->t - scores->span) < 2*score_set.max_gap) 
      continue;				// Sequence too short; skip.
    if (valid_n % BLKSIZ == 0) {	/* make space in scores */
      mm_resize(score_set.scores, valid_n + BLKSIZ, SCORE);
    }
    score_set.scores[valid_n++] = *scores;
  }
  score_set.n = valid_n;

  // Bail if not enough scores to do estimation.
  if (valid_n < min_scores) {
    (void) sprintf(
       evd_set.msg,
       "Too few (%d < %d) scores from matches with short enough spans to estimate E-values.", 
        valid_n, min_scores
    );
    fprintf(stderr, "%s\n", evd_set.msg);
    return(evd_set);
  }

  // Sort scores by GC-content in ascending order.
  if (nranges > 1) {
    qsort((char *)score_set.scores, valid_n, (int)sizeof(SCORE), gc_compare);
  }

  // Setup distribution for different GC-content ranges.
  scores = score_set.scores;		// Don't need original list anymore.
  i_low = 0;				// Start of current range of GC-content.
  gc_mid = scores[0].gc/2;		// mid-point GC of previous range.
  gc_hi = 0;				// Previous range maximum GC-content.
  for (g=0; g<nranges; g++) {		// Range.
    int count;				// Number of scores in range.
    double gc_sum = 0;			// Sum of GC-contents.
    EVD evd;				// Current EVD.
    for (i=i_low, count=1; i<valid_n; i++, count++) {	// Scores.
      double gc_mid_new = gc_sum/count;
      double ratio = gc_mid <= 0.5 ? gc_mid_new/gc_mid : (1-gc_mid)/(1-gc_mid_new);
      gc_sum += scores[i].gc;
      // OK to end range here? 
      if (
        i == valid_n - 1 || (			// Always end if last score.
          count >= min_scores &&		// Need at least min_scores
          (ratio >= max_ratio) &&		// and ratio reached 
          i < n-min_scores) 			// and at least min_scores are left.
        ) {
        // Set up score_set to contain just this range.
        score_set.scores = scores+i_low;
        score_set.n = count;
        // Fit scores to distribution.
        evd = em_exp(score_set, maxiter, eps1);
        evd.g = g; 
        evd.n = count;
        evd.min_t = gc_hi;		// Low GC is previous range hi GC.
        evd.max_t = scores[i].gc;	// Hi GC.
        evd.mid_t = gc_sum/count;	// Mean GC.
        gc_mid = evd.mid_t;		// Save new mean GC.
        gc_hi = evd.max_t;		// Save new hi GC.
	// Flag problem if exponential mean greater than Gaussian: fitting wrong tail.
	if (evd.mu1 > evd.mu2) {
	  (void) sprintf(
	     evd_set.msg,
	     "Unable to estimate E-values because they don't follow an exponential distribution (mu1 %.2f, mu2 %.2f).",
	      evd.mu1, evd.mu2
	  );
	  fprintf(stderr, "%s\n", evd_set.msg);
	  evd_set.n = 0;			// Failure.
          return(evd_set);
	}
        // Save this distribution.
        evds[g] = evd;
        i_low = i+1;
//fprintf(stderr, "g %d min_t %.2f mid_t %.2f max_t %.2f n %d mu1 %.2f\n", g, evd.min_t, evd.mid_t, evd.max_t, evd.n, evd.mu1);
        break;
      } // Set up range.
    } // Scores.
    if (i == valid_n - 1) break;		// All done.
  } // Range.
  nranges = g+1;

  // Set max_t of last group to 1.
  evds[nranges-1].max_t = 1;

  //
  // Set up return stuff.
  //
  evd_set.n = nranges;			// Success.
  evd_set.evds = evds;
  return(evd_set);
} // fitexp 

/******************************************************************************/
/*
	ml_evd

	Get maximum likelihood estimates for an EVD distribution given
	a set of (score, length) pairs.

	A set of extreme value distributions is returned corresponding
	to target sequences in distinct length ranges.  The length ranges
	are determined by:
		"size" 	-- the (approx.) number of sequences per range
		"lspan" -- the minimum ratio of longest to shortest sequence
			   in a single length range

	An initial value of H, 0 <= H <= 1000 must be given.
	If H=0, the length-adjustment will not be made to the effective
	search space size.  H is ignored for sequences longer than
	300 residues since it has little effect for such long sequences.

	If maxscores>0, approximately that many randomly selected scores 
	are used in estimating the EVD.  If maxscores==0, all scores are used.  

	Returns evd_set.  If an error occurs, evd_set.n = 0 is returned.
*/
/******************************************************************************/
static EVD_SET ml_evd(
  SCORE_SET score_set,				/* the set of scores and lens */
  double H,					/* initial estimate for H */
  int maxiter1,					/* maximum iterations for ML */
  int maxiter2,					/* maximum iterations for N-R */
  double eps1,					/* error tolerance for ML */
  double eps2, 					/* error tolerance for N-R */
  int maxscores,				/* randomly select scores */
  int size,					/* size of length ranges */
  double lspan,					/* minimum length ratio/group */
  double ethresh				/* outlier threshold */
)
{
  int i, g; 
  double lower;
  int count, total, left;		/* counts of scores */
  int n = score_set.n;			/* number of scores */
  int old_n = n;			/* original n */
  SCORE *scores = score_set.scores;	/* array of scores */
  int max_t = score_set.max_t;		/* maximum target length */
  int *t_hist = NULL;			/* histogram of target lengths */
  EVD_SET evd_set;			/* return value */
  evd_set.evds = NULL;   /* Initialize evds member */
  EVD *evds = NULL;			/* array of EVDs */
  EVD evd_init;				/* initial EVD */
  int nranges = n/size + 1;		/* maximum number of length ranges */

  /* create arrays */
  mm_resize(t_hist, max_t+1, int);
  mm_resize(evds, nranges, EVD);

  // Using EVD, not exponential distribution.
  evd_set.exponential = FALSE;

  // Set up error signal.
  evd_set.n = 0; 

  /* check that there are enough points */
  if (n < MIN_SCORES) {
    (void) sprintf(
      evd_set.msg,
      "Too few (%d < %d) scores to estimate E-values.", 
      n, MIN_SCORES
    );
    fprintf(stderr, "%s\n", evd_set.msg);
    myfree(t_hist);
    myfree(evds);
    return(evd_set);
  }

  /* subsample the scores if the size of the set exceeds the maximum
     number of scores allowed
  */
  if (old_n > maxscores && maxscores > 0) {	/* too many scores */
    double frac = 1.0-(maxscores+0.0)/n;	/* fraction to sample */
    score_set.n = 0;				/* empty list of scores */
    score_set.scores = NULL;			/* initialize score list */
    srand48(0);					/* init rn generator */
    for (i=0; i<n; i++) {	 		/* subsample scores */
      if (drand48() >= frac) {			/* use this score */
	if (score_set.n % BLKSIZ == 0) { 
	  mm_resize(score_set.scores, score_set.n + BLKSIZ, SCORE);
	}
        score_set.scores[score_set.n++] = scores[i];
      }
    } /* subsample scores */
    scores = score_set.scores;			/* toss original scores */
    n = score_set.n;
  } /* too many scores */
  if (TRACE) fprintf(stderr, "n = %d\n", n);

  /* 
    get initial estimate of lambda and K and tag outliers setting pv=1
  */
  evd_init = get_est(score_set);
  evd_init.H = H <= MAXH ? H : MAXH;		/* limit size of H */
  if (evd_init.lambda == 0) {
    (void) sprintf( evd_set.msg, "Too few scores to estimate E-values.");
    fprintf(stderr, "%s\n", evd_set.msg);
    myfree(t_hist);
    myfree(evds);
    return(evd_set);
  }
      
  /* 
    get length histogram of non-outlier scores 
  */
  for (i=0; i<=max_t; i++) t_hist[i] = 0;	/* 0 counts of each length */
  for (i=total=0; i<n; i++) {
    if (scores[i].pv != 1) {			/* score */
      t_hist[scores[i].t]++;			/* # of non-outliers */
      total++;					/* total non-outliers */
    }
  } /* score */

  /* 
    get length ranges with about 'size' non-outliers each, but always
    covering at least a factor of lspan in lengths 
  */
  left = total;					/* number of scores left */
  for (i=g=count=0,lower=0; i<=max_t; i++) {	/* get cutoffs */
    count += t_hist[i];				/* number of scores <= i */
    left -= t_hist[i];				/* remaining after this group */
    if (left < 0.9*size) {			/* not enough for two groups */
      i = max_t;
      count = count+left;			/* combine into 1 group */
    }
    if ((count>=size && i/(lower+1)>=lspan) || i==max_t) {	/* new group */
      EVD evd = evd_init;			/* initialize lambda etc */
      evd.g = g;				/* group */
      evd.n = count;				/* size of group */
      count = 0;
      evd.min_t = lower;			/* lower length bound */
      evd.max_t = i;				/* maximum length */
      evd.mid_t = (i+lower)/2.0;		/* midpoint of lengths */
      lower = evd.max_t;			/* new lower bound */
      evds[g++] = evd;
    } /* new group */
  } /* get cutoffs */
  nranges = g;					/* set actual number ranges */

  /* 
    perform ML estimation for each length range 
  */
  for (g=0; g<nranges; g++) {			/* group */
    EVD evd = evds[g];				/* initialize evd */
    double lower = evd.min_t;
    double max_t = evd.max_t;
    double mid_t = 0;				/* get average t */

    if (lower > 300) evd.H = 0;			/* no H for long sequences */
    score_set.n = 0;				/* empty list of scores */
    score_set.scores = NULL;			/* initialize score list */
    for (i=0; i<n; i++) {	 		/* scores */
      int t = scores[i].t;			/* target length */
      if (t>lower && t<=max_t) {		/* use this score */
	if (score_set.n % BLKSIZ == 0) { 	/* make space in scores */
	  mm_resize(score_set.scores, score_set.n + BLKSIZ, SCORE);
	}
        score_set.scores[score_set.n++] = scores[i];
        mid_t += t; 
      } /* in current group */
    } /* scores */

    evd.n = score_set.n;			/* number of scores */
    evd.mid_t = mid_t/score_set.n;		/* mean of scores */
    if (TRACE) {
      printf("# g %d max_t %.0f gsize %d\n", evd.g+1, evd.max_t, evd.n);
      fprintf(stderr, "min_t %.0f mid_t %f max_t %.0f count %d\n", evd.min_t, 
	evd.mid_t, evd.max_t, evd.n);
    }

    /* update estimate of EVD */
    evd = ml_evd_nr(evd, score_set, maxiter1, maxiter2, eps1, eps2,ethresh);

    if (TRACE) 
      fprintf(stderr, "ml_evd: max_t %.0f = %d lambda = %f K = %f H = %f\n", 
        max_t, score_set.n, evd.lambda, evd.K, evd.H);

    /* free temporary score set */
    free(score_set.scores);

    /* exit if error occured */
    if (evd.lambda == 0) {
      (void) sprintf(evd_set.msg, "Too few scores to estimate E-values.");
      fprintf(stderr, "%s\n", evd_set.msg);
      return(evd_set);
    }

    /* set new evd */
    evds[g] = evd;
  } /* group */

  /* free local space if it was created */
  if (old_n > maxscores && maxscores > 0) {		/* too many scores */
    free(scores);
  }
    
  /* return evd_set */
  evd_set.n = nranges;
  evd_set.evds = evds;

  // Free local dynamic memory.
  myfree(t_hist);

  return(evd_set);
} /* ml_evd */

/******************************************************************************/
/*
	evd_set_pvalue

	Compute the p-value of a (score, length) pair from a set of EVDs
	applicable to different length ranges.

  	P-values of scores are the weighted geometric mean of
     	nearest two length (or GC-content) ranges, weighting by 0.5 + fraction of 
	distance from the boundary to midpoint of the length group containing
	the target sequence length (or GC-content).
*/
/******************************************************************************/
double evd_set_pvalue(
  double s,					/* score */
  double t,					/* target length (or GC-content) */
  int q,					/* query length */
  EVD_SET evd_set				/* set of EVDs */
)
{
  double N, el, pv, w;
  int g, g1;
  int nranges = evd_set.n;			/* number of length ranges */
  EVD *evds = evd_set.evds;			/* EVDs */
  EVD evd;

  /* find length/GC-content group */
  for (g=0; g<nranges; g++) {			/* find group */
    if (t>=evds[g].min_t && t<=evds[g].max_t) break;	/* found group */
  }
  if (g==nranges) g--;				/* put in last group */
  evd = evds[g];

  if (evd_set.exponential) {
    EVD evd;
    // Get two midpoints t is between.  g is the lower, g1 is the upper.
    if (t > evds[g].mid_t) {			// above midpoint 
      g1 = (g<nranges-1) ? g+1 : g;
    } else if (t < evds[g].mid_t)  {		// below midpoint
      g1 = g;
      g = (g>0) ? g-1 : 0;
    } else {
      g1 = g;
    }
    // Linearly interpolate the value of mu.
    evd.mu1 = (g==g1) ? evds[g].mu1 : evds[g].mu1 + 
      (t - evds[g].mid_t)*((evds[g1].mu1 - evds[g].mu1)/(evds[g1].mid_t - evds[g].mid_t));
    // get p-value given local GC-content.
    pv = Exp_pvalue(evd, s);
    return(pv);
  } else {
    /* get p-value in correct length group */
    Evd_pvalue(pv, N, el, evds[g], s, t, q);
  }

  /* get the weight for the correct length range and identity of next
     closest length range
  */
  w = 1;
  if (t >= evd.mid_t) {				/* above midpoint */
    g1 = g+1;
    if (g1 < nranges) {
      w = 0.5 + (0.5*(evd.max_t-t)) / (evd.max_t-evd.mid_t);
    }
  } else {					/* below midpoint */
    g1 = g-1;
    if (g1 >= 0) {
      w = 0.5 + (0.5*(t-evd.min_t)) / (evd.mid_t-evd.min_t);
    }
  }

  /* check for error in w */
  if (w < 0.5 || w > 1.0) {
    fprintf(stderr, "evd_set_pvalue error: t %.3f g %d g1 %d min %.3f mid %.3f max %.3f w = %f nranges %d\n", 
	t, g, g1, evd.min_t, evd.mid_t, evd.max_t, w, nranges);
    exit(0);
  }

  /* compute the geometric mean */
  if (w >= 0.5 && w < 1.0) {
    double pv1;
    if (evd_set.exponential) {
      // get p-value given local GC-content.
      pv1 = Exp_pvalue(evds[g1], s);
    } else {
      /* get p-value in correct length group */
      Evd_pvalue(pv1, N, el, evds[g1], s, t, q);
    }
    pv = pow(pv, w) * pow(pv1, 1.0-w);		/* geometric mean */
  }

  return(pv);
} /* evd_set_pvalue */

/******************************************************************************/
/*
	get_est

	Bin scores by length, remove outliers, measure variance to get initial 
	estimate of lambda; tag outliers with 1 pv.
*/
/******************************************************************************/
EVD get_est (
  SCORE_SET score_set			/* set of scores */
)
{
  int i, n1, bin, bin_cnt;
  double avgt;				/* temp variable */
  int n = score_set.n;			/* number of scores */
  int max_t = score_set.max_t;		/* largest length */
  int max_bin = Score_bin(max_t);	/* last bin */
  double *mu = NULL;			/* bin means */
  double *var = NULL;			/* bin means */
  int *count = NULL;			/* bin counts */
  double sd, newsd;			/* residual standard deviation */
  EVD evd;				/* return value */
  double factor = 3.1415/sqrt(6);	/* pi/sqrt(6) */

  /* initialize */
  init_evd(&evd);
  mm_resize(mu, max_bin+1, double);
  mm_resize(var, max_bin+1, double);
  mm_resize(count, max_bin+1, int);
  for (i=0; i<=max_bin; i++) { mu[i] = 0; var[i] = 0; count[i] = 0; }

  /* get bin average lengths */
  for (i=1,avgt=bin_cnt=0,bin=Score_bin(i); i<=max_t; avgt+=i, bin_cnt++, i++) {
    int new_bin = Score_bin(i);
    if (new_bin > bin) {			/* new bin found */
      avgt = bin_cnt = 0; bin = Score_bin(i);	/* initialize next bin */
    } 
  }

  /* get bin means */
  for (i=0; i<n; i++) {				/* scores */
    int bin = Score_bin(score_set.scores[i].t);
    mu[bin] += score_set.scores[i].s;
    count[bin]++; 
  } /* scores */
  for (i=0; i<=max_bin; i++) { if (count[i]) mu[i] /= count[i]; }

  /* get bin variances */
  for (i=sd=0; i<n; i++) {			/* scores */
    int bin = Score_bin(score_set.scores[i].t);
    double d = mu[bin] - score_set.scores[i].s;
    var[bin] += d*d;				/* bin variance */
    sd += d*d;					/* residual variance */
  } /* scores */
  for (i=0; i<=max_bin; i++) { if (count[i]>1) var[i] /= count[i]-1; }
  if (n>1) sd /= (n-1);
  sd = sqrt(sd);

  /*
    recompute mean and variance after removing outliers (> 5sd) 
  */

  /* recompute bin means removing outliers */
  for (i=0; i<n; i++) {				/* scores */
    int bin = Score_bin(score_set.scores[i].t);
    double d = mu[bin] - score_set.scores[i].s;
    if (fabs(d/sd) >= 5) { 			/* remove outlier */
      if (count[bin] > 1) {
        mu[bin] = (mu[bin]*count[bin] - score_set.scores[i].s)/--count[bin];
      } else {
        mu[bin] = count[bin] = 0;
      }
    } /* remove outlier */
  } /* scores */

  /* recompute variance removing outliers and tagging them */
  for (i=n1=newsd=0; i<n; i++) {			/* scores */
    int bin = Score_bin(score_set.scores[i].t);
    double d = mu[bin] - score_set.scores[i].s;
    if (fabs(d/sd) < 5) { 			/* not an outlier */
      newsd += d*d;
      n1++;
      score_set.scores[i].pv = 0;
    } else {
      score_set.scores[i].pv = 1;
    }
  } /* scores */
  for (i=0; i<=max_bin; i++) { if (count[i]>1) var[i] /= count[i]-1; }
  sd = newsd;
  if(n1>1) sd /= (n1-1);
  sd = sqrt(sd);

  if (sd == 0) {
    evd.lambda = 0;			// Signal error.
    evd.K = 0;
    return(evd);
  }

  /* estimate lambda = pi/(sqrt(6) sigma) */
  evd.lambda = sd ? factor/sd : 0;
  evd.K = 0.01;				/* seems to work */

  /* fprintf(stderr, "initial lambda = %8.3f\n", evd.lambda);*/

  /* free */
  free(mu); free(var); free(count);

  return(evd);
} /* get_est */

/******************************************************************************/
/*
	ml_evd_nr

	Find the maximum likelihood estimates of the parameters of an EVD
	using Newton-Raphson iteratively on lambda then H.
*/
/******************************************************************************/
EVD ml_evd_nr(
  EVD evd,					/* initial guess */
  SCORE_SET score_set,				/* the set of scores and lens */
  int maxiter1,					/* maximum iterations for ML */
  int maxiter2,					/* maximum iterations for N-R */
  double eps1,					/* error tolerance for ML */
  double eps2, 					/* error tolerance for N-R */
  double ethresh				/* outlier threshold */
)
{
  int i, iter;
  SCORE_SET cur_score_set;			/* temporary array of scores */
  int q = score_set.q;				/* length of query */
  int n = score_set.n;				/* total number of scores */
  double old_n;

  /* create a temporary score set */
  cur_score_set.scores = NULL;
  mm_resize(cur_score_set.scores, score_set.n, SCORE);
  cur_score_set.q = q;
  
  /* initialize log-likelihood */
  evd.L = -BIG;					/* very small */

  /* initialize number of scores being used */
  evd.n = 0;					/* number of non-outliers */

  for (iter=0; iter<maxiter1; iter++) {
    double H = evd.H;				/* paramter H */
    SCORE *scores = score_set.scores;		/* scores */
    SCORE *cur_scores = cur_score_set.scores;	/* current scores */
    double L;					/* log likelihood */

    /* move non-outliers to the current score list */
    for (i=0; i<n; i++, scores++) {		/* score */
      double s = scores->s;			/* score of target */
      int t = scores->t;			/* length of target */
      if (iter) {				/* move all first iteration */
        double pv, N, el;
        Evd_pvalue(pv, N, el, evd, s, t, q);
        if (pv*evd.n < ethresh)  continue;	/* skip if E-value < ethresh */
      }
      cur_scores->s = s; cur_scores->t = t; cur_scores++;
    } /* score */


    /* number of scores in current set */
    cur_score_set.n = cur_scores - cur_score_set.scores;
    old_n = evd.n;				/* old number of non-outliers */
    evd.n = cur_score_set.n;			/* number of non-outliers */
    if (cur_score_set.n < MIN_SCORES) {
      fprintf(stderr, "Too few scores to estimate E-values.\n");
      evd.lambda = 0;				// Signal error.
      return(evd);
    }

    /* maximize the likelihood */
    L = evd.L;					/* current value of log-like */
    evd = ml_lambda(evd, cur_score_set, eps2, maxiter2);
    if (TRACE) fprintf(stderr, 
      "evd l lambda %8.3f K %8.3f H %8.3f L %15.2f L/n %9.6f %12d\n", 
      evd.lambda, evd.K, evd.H, evd.L, evd.L/evd.n, evd.n);
    if (H) {
      evd = ml_H(evd, cur_score_set, eps2, maxiter2);
      if (TRACE) fprintf(stderr, 
	"H lambda %8.3f K %8.3f H %8.3f L %15.2f L/n %9.6f %12d del %g\n", 
	evd.lambda, evd.K, evd.H, evd.L, evd.L/evd.n, evd.n, fabs((L-evd.L)/L));
    }
    if (evd.n == old_n && fabs((L-evd.L)/L) < eps1) break;	/* converged */
  } /* iterations */

  /*fprintf(stderr, "# iter1 = %d\n", iter);*/
  if (TRACE) fprintf(stderr, "L %15.2f\n", evd.L);

  /* free temporary space */
  free(cur_score_set.scores);

  return(evd);
} /* ml_evd_nr */

/******************************************************************************/
/*
	ml_lambda

	Find the maximum likelihood estimates of lambda and K holding
	H constant.

	Returns the estimates.
*/
/******************************************************************************/
EVD ml_lambda(
  EVD evd,					/* current parameters of EVD */ 
  SCORE_SET score_set,				/* set of scores */
  double eps,					/* stopping criterion */
  int maxiter					/* max. allowed iterations */
)
{
  int iter;
  DERIV d;					/* 1st 2 derivs, log-like */

  /* Newton-Raphson loop */ 
  d.L = -BIG;					/* very small */
  for (iter=0; iter<maxiter; iter++) {
    double lambda = evd.lambda;			/* current value of lambda */
    double L = d.L;				/* current log likelihood */

    d = dlambda(evd, score_set); 		/* get the first two derivs */
    /*fprintf(stderr, "ml_lambda: L= %f K= %f d1= %e\n", d.L, d.K, d.d1);*/

    if (fabs(d.d1) < eps) break;		/* converged */
    if (d.L < L) { d.L = L; break; }		/* failed; reset log-like */

    evd.lambda -= d.d1/d.d2;			/* Newton-Raphson step */
    if (evd.lambda <= 0) { evd.lambda = lambda/2.0; }	/* binary search */
  } /* Newton-Raphson */

  /*fprintf(stderr, "# iter2 (lambda) : %d\n", iter);*/

  evd.L = d.L;					/* log likelihood */
  evd.K = d.K;					/* optimum K */

  return(evd);
} /* ml_lambda */

/******************************************************************************/
/*
	ml_H

	Find the maximum likelihood estimates of H holding
	lambda and K holding constant.

	Returns the estimates.
*/
/******************************************************************************/
EVD ml_H(
  EVD evd,					/* current parameters of EVD */ 
  SCORE_SET score_set,				/* set of scores */
  double eps,					/* stopping criterion */
  int maxiter					/* max. allowed iterations */
)
{
  int iter;
  DERIV d;					/* 1st 2 derivs and log like */

  /* Newton-Raphson loop */
  d.L = -BIG;					/* very small */
  for (iter=0; iter<maxiter; iter++) {
    double L = d.L;				/* current log likelihood */
    double H = evd.H;				/* current value of H */

    d = dH(evd, score_set); 			/* get the first two derivs */
    /*fprintf(stderr, "ml_H: H= %f L= %f d1 %g d2 %g\n", 
      evd.H, d.L, d.d1, d.d2);*/

    if (fabs(d.d1) < eps) break;		/* converged */
    if (d.L < L) { d.L = L; break;}		/* failed; reset log-like */

    if (d.d2 > 0) {				/* binary search w/1st deriv.*/
      if (d.d1 > 0) { 
        evd.H *= 2; 
      } else {
        evd.H /= 2; 
      }
      continue;
    }
    evd.H -= d.d1/d.d2;				/* Newton-Raphson step */
    if (evd.H <= 0) { evd.H = H/2.0; }		/* binary search */
    if (evd.H > MAXH) { break; }		/* don't allow H too large */
  } /* Newton-Raphson */

  /*fprintf(stderr, "# iter2 (H     ) : %d\n", iter);*/

  evd.L = d.L;					/* log likelihood */

  return(evd);
} /* ml_H */

/******************************************************************************/
/*
	dlambda

	Compute the first and second derivatives of the log likelihood
	with respect to lambda, holding H fixed.
*/
/******************************************************************************/
DERIV dlambda(
  EVD evd,					/* current parameters of EVD */ 
  SCORE_SET score_set				/* set of scores */
)
{
  int i;
  double lambda = evd.lambda, K = evd.K, H = evd.H; 
  int n = score_set.n;
  int q = score_set.q;
  double N, l;
  DERIV d;					/* the results */
  double t2, t3, t4, t5, sum0, sum1, sum2, sum3, sum4;

  sum0 = sum1 = sum3 = sum4 = 0;
  sum2 = 1e-200;				/* avoid divide by zero */
  for (i=0; i<n; i++) {				/* score */
    double s = score_set.scores[i].s;
    int t = score_set.scores[i].t;
    Get_N(N, l, q, t, K, H);
    t2 = N * exp(-lambda*s);
    t3 = s * t2;
    t4 = s * t3;
    sum0 += log(N);
    sum1 += s;
    sum2 += t2;
    sum3 += t3;
    sum4 += t4;
  } /* score */
  
  /* compute the derivatives and the log-likelihood */
  t5 =  sum3/sum2;
  d.d1 = 1/lambda - sum1/n + t5;
  d.d2 = -(1/(lambda*lambda)) - sum4/sum2 + t5*t5;
  d.K = n/sum2;
  d.L = n*log(lambda*d.K) + sum0 - lambda*sum1 - d.K*sum2;

  return(d);
} /* dlambda */

/******************************************************************************/
/*
	dH

	Compute the first and second derivatives of the log likelihood
	with respect to H, holding lambda and K fixed.
*/
/******************************************************************************/
DERIV dH(
  EVD evd,					/* current parameters of EVD */ 
  SCORE_SET score_set				/* set of scores */
)
{
  int i;
  double lambda = evd.lambda, K = evd.K, H = evd.H; 
  int n = score_set.n;
  int q = score_set.q;
  double N, l=0;
  DERIV d;					/* the results */
  double t1, t2, t3, t4, t5, t6;

  d.d1 = d.d2 = d.L = 0;
  for (i=0; i<n; i++) {				/* score */
    double s = score_set.scores[i].s;
    int t = score_set.scores[i].t;
    Get_N(N, l, q, t, K, H);
    t1 = 2*l - q - t;
    t2 = K * exp(-lambda*s);
    t3 = 1/N - t2;
    t4 = -l/H;
    t5 = t1 * t3 * t4;
    t6 = t1 * t4 / N;
    d.d1 += t5;
    d.d2 += 2*t3*t4*t4 - t6*t6 - 2*t5/H;
    d.L += log(lambda*K*N) - lambda*s - N*t2;
  } /* score */
  d.K = K;
  /*fprintf(stdout, ": L= %f lambda= %f K= %f H= %f\n", d.L, lambda, K, H);*/
  
  return(d);
} /* dH */

/******************************************************************************/
/*
	get_n

	Get the expected number of scores for use in converting
	p-values to E-values.

	Returns n, the number of non-outliers.
*/
/******************************************************************************/
extern int get_n(
  SCORE_SET score_set,			/* the set of scores and lens/gcs */
  EVD_SET evd_set 			// Distribution.
) 
{
  int i;
  SCORE *scores = score_set.scores;
  int n = score_set.n;			// Number of scores.
  int n_out;				// Current estimate of # outliers.
  double thresh;			// p-value threshold for outliers.

  // Compute the p-value of each score.
  for (i=0; i<n; i++, scores++) {
    scores->pv = Evd_set_evalue(
      1, 
      scores->s, 
      evd_set.exponential ? scores->gc : scores->t,
      1,
      evd_set
    );
  }

  /* sort the scores by p-value ascending order */
  qsort((char *)score_set.scores, n, (int)sizeof(SCORE), p_compare);

  // Iteratively count the number of outliers and subtract from total
  // number of scores until it does not increase.
  n_out = 0;
  while (n_out < n) {
    thresh = 1.0/(n-n_out);		// p-value outlier threshold.
    scores = score_set.scores;
    for (i=0; i<n && scores->pv < thresh; i++, scores++) 
      ;
    if (i > n_out) {			// Est. of outliers increased.
      n_out = i;			// New estimate.
    } else {
      break;				// Converged.
    }
  }
  return(n-n_out);
} // get_n

/******************************************************************************/
/*
	fit_score_distribution

 	Fit an exponential or extreme value distribution to a set of scores.

	See fitexp() and ml_evd() for meanings of parameters.
*/
/******************************************************************************/
extern EVD_SET fit_score_distribution(
  BOOLEAN_T exponential,			// Fit exponential distribution?
  SCORE_SET score_set,				/* the set of scores and lens */
  double H,					/* initial estimate for H */
  int maxiter1,					/* maximum iterations for ML */
  int maxiter2,					/* maximum iterations for N-R */
  double eps1,					/* error tolerance for ML */
  double eps2, 					/* error tolerance for N-R */
  int size,					/* size of length ranges */
  double lspan					/* minimum length ratio/group */
) 
{
  EVD_SET evd_set;
  evd_set.evds = NULL; // Initialize evds member

  // Fit the appropriate distribution to scores and lengths.
  if (exponential) {
    int maxiter = 1000;				// Limit EM to this.
    evd_set = fitexp(score_set, maxiter, eps1, MIN_SCORES, MAX_RATIO);
  } else {
    evd_set = ml_evd(score_set, H, maxiter1, maxiter2, eps1, eps2, 0,
                     size, lspan, 1);
  }

  // Get the total database length.
  evd_set.total_length = score_set.total_length;

  // Get the total number of matches.
  evd_set.nscores = score_set.n;

  // Get the number of non-outliers.
  if (evd_set.n > 0) {				// Succeeded.
    evd_set.non_outliers = get_n(score_set, evd_set);	// Number with E > 1.

    // See if enough scores remain after removing outliers for the
    // estimate to be reliable.
    if (evd_set.non_outliers < MIN_SCORES) {
      (void) sprintf(
        evd_set.msg,
        "After removing outliers, too few (%d < %d) scores remain to estimate E-values.", 
	evd_set.non_outliers, MIN_SCORES
      ); 
      fprintf(stderr, "%s\n", evd_set.msg);
      evd_set.n = 0;				// Signal error.
    }

  } // Succeeded

  if (evd_set.n > 0) strcpy(evd_set.msg, "Distribution estimated successfully.");

  return(evd_set);
} // fit_score_distribution

#ifdef MAIN
/******************************************************************************/
/*
	main
*/
/******************************************************************************/
int main(int argc, char **argv) 
{
  int i;
  int q;					/* length of query */
  double s;
  double s_min = 50;				/* minimum score to trace */
  double s_max = 50;				/* maximum score to trace */
  int ptrace = 0;				/* don't do ptrace */
  int size = 10000;				/* size of length groups */
  double lspan = 1.5;				/* min. length ratio/group */
  double H = 1.0;				/* initial H estimate */
  int maxscores = 0;				/* maximum number to use; all */
  int pvalues = 0;				/* don't print p-values */
  int maxiter1 = 10;				/* max. iterations in L-H loop*/
  int maxiter2 = 20;				/* max. iterations in N-R loop*/
  SCORE_SET score_set;				/* scores/lengths */
  EVD_SET evd_set;				/* extreme value distribution */

  /* check command line arguments */
  if (argc <= 1) {
    usage:
    fprintf(stderr, "USAGE: fitevd <q> [-s <s>] [-lr <lr>] [-r <n>] [-h <H> [-p] [-e <ethresh>] [-t]\n");
    fprintf(stderr, "\t\t[-ilh <ilh>] [-inr <inr>] [-d <smin> <smax>]\n\n");
    fprintf(stderr, "\t<q>\t\tquery length\n");
    fprintf(stderr, "\t[-s <s>]\tsize of length groups; default: %d\n", 
      size);
    fprintf(stderr, "\t[-lr <lr>]\tminimum ratio of shortest to longest in group;\n"); 
    fprintf(stderr, "\t\t\t\tdefault %.2f\n", lspan);
    fprintf(stderr, "\t[-r <n>]\tuse <n> randomly selected scores; default: all\n");
    fprintf(stderr, "\t[-p]\t\tprint p-values; default: don't print p-values\n");
    fprintf(stderr, "\t[-h <H>]	initial value of H; default: %3.1f\n", H);
    fprintf(stderr, "\t\t\tno edge-effect adjustment if <H> is 0.\n");
    fprintf(stderr, "\t[-ilh <ilh>]\tturn max. iterations in outer L-H loop; default: %d\n", maxiter1);
    fprintf(stderr, "\t[-inr <inr>]\tturn max. iterations in inner N-R loop; default: %d\n", maxiter2);
    fprintf(stderr, "\t[-t]\t\tturn on debug tracing; default: tracing off\n");
    fprintf(stderr, "\t[-d <smin> <smax>] print length vs p-value v for <smin> <= score <= <smax>\n");
    fprintf(stderr, "\n\tInput is from standard input:\n\t\t[score length]+\n");
    fprintf(stderr, "\n\tThe outer loop alternates between estimating lambda and K with H fixed\n");
    fprintf(stderr, "\tand estimating H with lambda and K fixed.  The inner loop is the Newton-Raphson\n");
    fprintf(stderr, "\talgorithm, which is used for solving the estimation problems.\n");
    
    return(1);
  }

  /* get command line arguments */
  i = 1;
  q = atoi(argv[i++]);				/* get length of query */
  while (i < argc) { 				/* get switches */
    if (!strcmp(argv[i], "-h")) {
      H = atof(argv[++i]);			/* initial H estimate */
    } else if (!strcmp(argv[i], "-lr")) {
      lspan = atof(argv[++i]);			/* min. length ratio/group */
    } else if (!strcmp(argv[i], "-s")) {
      size = atoi(argv[++i]);			/* size of length groups */
    } else if (!strcmp(argv[i], "-r")) {
      maxscores = atoi(argv[++i]);		/* maximum number to use */
    } else if (!strcmp(argv[i], "-p")) {
      pvalues = 1;				/* print p-values */
    } else if (!strcmp(argv[i], "-ilh")) {
      maxiter1 = atoi(argv[++i]);
    } else if (!strcmp(argv[i], "-inr")) {
      maxiter2 = atoi(argv[++i]);
    } else if (!strcmp(argv[i], "-t")) {
      TRACE = 1; 				/* turn on debug trace */
    } else if (!strcmp(argv[i], "-d")) {
      ptrace = 1; 				/* turn on ptraceing trace */
      s_min = atof(argv[++i]);			/* minimum score to print */
      s_max = atof(argv[++i]);			/* minimum score to print */
    } else {
      goto usage;
    }
    i++;
  }

  /* check input arguments */
  if (H > MAXH) {				/* check that H is in range */
    fprintf(stderr, "H must be <= %d\n", MAXH);
    return(1);
  }

  /* read in a set of scores and lengths */
  score_set = read_scores(stdin);		/* read scores */
  score_set.q = q;				/* set query length */

  /* estimate the parameters of the EVD */
  evd_set = ml_evd(score_set,
		   H,
		   maxiter1,
		   maxiter2,
		   EPS1, 
		   EPS2,
		   maxscores, 
		   size, 
		   lspan,
		   1.0); // Use ethresh=1.0 -- WSN 8/13/03

  /* 
    print EVD parameters
  */
  for (i=0; i<evd_set.n; i++) {			/* length range */
    EVD evd = evd_set.evds[i];
    printf("# t %g lambda %8.4f K %8.4f H %8.4f q %d L %10.2f L/n %10.6f n %6d\n", 
      evd.max_t, evd.lambda, evd.K, evd.H, q, evd.L, evd.L/evd.n, evd.n);
  } /* length range */

  /* 
    print the p-values for fixed scores as t is varied
  */
  if (ptrace) {
    int t;
    double pv;
    //int g, t;
    //double pv, N, el;
    //EVD *evds = evd_set.evds;
    for (s=s_min; s<=s_max; s+=10) {
      printf ("\n# s = %f\n", s);
      for (t=20; t<score_set.max_t; t++) {      /* len */
	pv = evd_set_pvalue(s, t, q, evd_set);
	printf("%d %g\n", t, pv);
      }
    } /* s */
  } /* (t, pv) */

  /* 
    print the p-values of the input scores 
  */
  if (pvalues && !ptrace) {				
    SCORE *scores;

    /* compute the p-values from the scores */
    for (i=0, scores=score_set.scores; i<score_set.n; i++, scores++) {
      scores->pv = evd_set_pvalue(scores->s, scores->t, q, evd_set);
    }

    /* sort the scores by p-value ascending order */
    qsort((char *)score_set.scores, score_set.n, (int)sizeof(SCORE), p_compare);

    printf("# p-value  t-length  id\n");
    for (i=0, scores=score_set.scores; i<score_set.n; i++, scores++) {
      if (scores->id) {
        printf("%10.2e %6d %s\n", scores->pv, scores->t, scores->id);
      } else {
        printf("%10.2e %6d\n", scores->pv, scores->t);
      }
    }
  } /* print p-values */

  return(0);
} /* main */

/******************************************************************************/
/*
	read_scores

	Read a file of scores with the format:
		[<score> <target_length> <id>]+
*/
/******************************************************************************/
SCORE_SET read_scores ()
{
  SCORE_SET score_set;			/* set of scores to return */
  int len, nwords;			/* length of lines read */
  char *s=NULL, **words=NULL;		/* holds lines/words read */

  /* initialize */
  init_score_set(&score_set);

  while (1) {					/* loop over stdin */
    Getline(stdin, s, len);			/* read a line */
    if (s == NULL) break; 			/* at EOF when Getline called */
    Split(s, words, nwords); 
    if (score_set.n % BLKSIZ == 0) { 
      mm_resize(score_set.scores, score_set.n + BLKSIZ, SCORE);
    }
    score_set.scores[score_set.n].s = atof(words[0]);
    score_set.scores[score_set.n].t = atoi(words[1]);
    if (score_set.scores[score_set.n].t > score_set.max_t) { 
      score_set.max_t = score_set.scores[score_set.n].t; 
    }
    if (nwords == 3) {				/* id given */
      mm_resize(score_set.scores[score_set.n].id, strlen(words[2])+1, char);
      strcpy(score_set.scores[score_set.n].id, words[2]);
    } else {
      score_set.scores[score_set.n].id = NULL;
    } /* id */
    score_set.n++;
  } /* loop over stdin */

  return(score_set);
} /* read_scores */  

#endif /* MAIN */
