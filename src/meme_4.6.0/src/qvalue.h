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
#include "string-list.h"

#define NUM_BOOTSTRAPS 100          // Number of bootstraps to perform.
#define NUM_BOOTSTRAP_SAMPLES 1000  // Number of p-values in each bootstrap.
#define NUM_LAMBDA 100              // Number of lambda values to consider.
#define MAX_LAMBDA 0.5              // Maximum lambda value to consider.

/*************************************************************************
 * Read a set of floats from the specified column of a tab-delimited
 * file, skipping a specified number of header lines.
 *************************************************************************/
ARRAY_T* read_scores_from_column
  (char*          score_filename,    // Input filename
   int            num_header_lines,  // Number of lines to skip at start.
   int            score_column,      // Column index (indexed from 0)
   STRING_LIST_T* header_lines,      // Store header lines, if not NULL.
   STRING_LIST_T* input_lines        // Store file contents, if not NULL.
   );


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
  (char*    pi_zero_filename,      // Name of file to store intermediate values.
   int      num_bootstraps,
   int      num_bootstrap_samples, // Number of p-values in each bootstrap.
   int      num_lambda,      // Number of lambda values to consider.
   float    max_lambda,      // Maximum lambda value to consider.
   ARRAY_T* p_values);

/*************************************************************************
 * Verify that the given array contains scores in sorted order, with
 * good scores first.
 *************************************************************************/
void verify_sort
  (BOOLEAN_T good_score_is_low,
   char*     array_name,
   ARRAY_T*  my_array);

/*************************************************************************
 * Convert a set of p-values to a set of q-values.  The p-values must
 * be in sorted (ascending) order.  This function replaces them with
 * the corresponding q-values.
 *
 * In order to estimate pi0 we need to have the distribution of pvalues.
 * We can get the distribution from the pvalues array only if it contains
 * all of the observed pvalues. If we don't have all the observed p-values
 * we can use a uniformly sampled subset of pvalues 
 * If neither the full set of pvalues, or uniformly sampled pvalues
 * are avilable, then we can't estimate pi0, and simply use pi0 = 1.0.
 *************************************************************************/
void compute_qvalues(
  BOOLEAN_T compute_fdr,      // If true, only compute FDR, not q-value.
  BOOLEAN_T use_pi_zero,      // Estimate pi_zero; else just use 1.0.
  char*     pi_zero_filename, // Filename to store pi-zero estimate in.
  int       num_bootstraps,   // How many bootstraps to perform.
  int       num_bootstrap_samples, // Number of p-values in each bootstrap.
  int       num_lambda,       // Number of lambda values to consider.
  float     max_lambda,       // Maximum lambda value to consider.
  long      total_values,     // Total number of values.
  ARRAY_T*  pvalues,          // retained pvalues
  ARRAY_T*  sampled_pvalues   // uniformly sampled pvalues
);

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
