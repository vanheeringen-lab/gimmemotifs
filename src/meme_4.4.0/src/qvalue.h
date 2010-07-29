/*************************************************************************
 * FILE: qvalue.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 13 November 2007
 * PROJECT: MEME
 * COPYRIGHT: 2007, UW
 * DESCRIPTION: q value calculation, a la Benjamini Hochberg
 *************************************************************************/
#ifndef QVALUE_H
#define QVALUE_H

#include "utils.h"
#include "array.h"

#define NUM_BOOTSTRAPS 100          // Number of bootstraps to perform.
#define NUM_BOOTSTRAP_SAMPLES 1000  // Number of p-values in each bootstrap.
#define NUM_LAMBDA 100             // Number of lambda values to consider.
#define MAX_LAMBDA 0.99             // Maximum lambda value to consider.
/*****************************************************************************
 * Given a set of p-values, estimate the percentage of the empirical
 * distribution that is drawn according to the null hypothesis.
 *
 * This code implements the method described in John Storey, "A direct
 * approach to false discovery rates."  Journal of the Royal
 * Statistical Society B, 2002.
 *
 * NOTE: This function is only exported for use by the unit testing procedures.
 *****************************************************************************/
double estimate_pi_zero
  (int      num_bootstraps,
   int      num_bootstrap_samples, // Number of p-values in each bootstrap.
   int      num_lambda,      // Number of lambda values to consider.
   float    max_lambda,      // Maximum lambda value to consider.
   ARRAY_T* p_values);

/*************************************************************************
 * Convert a set of p-values to a set of q-values.  The p-values must
 * be in sorted (ascending) order.  This function replaces them with
 * the corresponding q-values.
 *
 * The total_values parameter may be larger than the number of values
 * stored in the array of p-values, but then pi_zero cannot be estimated.
 *************************************************************************/
void compute_qvalues
  (BOOLEAN_T compute_fdr,      // If true, only compute FDR, not q-value.
   BOOLEAN_T use_pi_zero,      // Estimate pi_zero; else just use 1.0.
   char*     pi_zero_filename, // File to store pi_zero estimate in.
   int       num_bootstraps,   // How many bootstraps to perform.
   int       num_bootstrap_samples, // Number of p-values in each bootstrap.
   int       num_lambda,       // Number of lambda values to consider.
   float     max_lambda,       // Maximum lambda value to consider.
   int       total_values,     // Total number of values.
   ARRAY_T*  pvalues);

/*****************************************************************************
 * Convert scores to p-values using an empirical null.
 *
 * The p-value associated with an observed score x is simply a/b,
 * where a is the number of null scores > x, and b is the total number
 * of null scores.
 *
 * NOTE: This function is only exported for use by the unit testing procedures.
 *****************************************************************************/
void convert_scores_to_pvalues
  (BOOLEAN_T good_score_is_low,
   ARRAY_T* observed_scores, 
   ARRAY_T* null_scores);
#endif
