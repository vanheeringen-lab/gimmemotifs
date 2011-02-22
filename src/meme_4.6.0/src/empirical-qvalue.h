/*************************************************************************
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 15 September 2009
 * PROJECT: MEME
 * COPYRIGHT: 2009, UW
 * DESCRIPTION: q value calculation from an empirical null distribution
 *************************************************************************/

/*************************************************************************
 * Convert a set of scores to a set of q-values based on an empirical
 * null score distribution.  For a given threshold, the q-value is
 * just the number of null scores above the threshold, divided by the
 * number of observed scores above the threshold.  The observed scores
 * are replaced by their corresponding q-values.
 *************************************************************************/
void compute_qvalues_empirical
  (BOOLEAN_T compute_fdr,      // If true, only compute FDR, not q-value.
   BOOLEAN_T good_score_is_low,
   ARRAY_T*  null_scores,
   ARRAY_T*  observed_scores);
