/*************************************************************************
 * FILE: qvalue.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 13 November 2007
 * PROJECT: MEME
 * COPYRIGHT: 2007, UW
 * DESCRIPTION: q value calculation, a la Benjamini Hochberg
 *************************************************************************/
#include "utils.h"
#include "qvalue.h"
#include "array.h"
#include "string-list.h"
#include <time.h>
#include <assert.h>
#include <stdio.h>

/*************************************************************************
 * Read a set of floats from the specified column of a tab-delimited
 * file, skipping a specified number of header lines.
 *************************************************************************/
#define MAX_SCORES 100000000         // Maximum allowed number of input lines.
#define MAX_LINE 1000                // Longest allowed input line.
ARRAY_T* read_scores_from_column
  (char*          score_filename,    // Input filename
   int            num_header_lines,  // Number of lines to skip at start.
   int            score_column,      // Column index (indexed from 0)
   STRING_LIST_T* header_lines,      // Store header lines, if not NULL.
   STRING_LIST_T* input_lines        // Store file contents, if not NULL.
   )
{

  // Read the p-values into a very large array.
  FILE* score_file = NULL;
  if (open_file(score_filename, "r", TRUE, "scores", "scores",
		&score_file) == 0) {
    exit(1);
  }

  // Skip the header lines.
  int i_header = 0;
  char line[MAX_LINE];    // Buffer for reading.
  while (i_header < num_header_lines) {
    if (fgets(line, MAX_LINE, score_file) == NULL) {
      die("End of file found while reading line %d of the header.\n",
	  i_header);
    }

    // Replace the EOL.
    int length = strlen(line);
    if (length > 0) {
      line[length - 1] = '\0';
    }

    // Print the line, if requested.
    if (header_lines != NULL) {
      add_string(line, header_lines);
    }

    i_header++;
  }

  // Read the p-values.
  ARRAY_T* initial_scores = allocate_array(MAX_SCORES);
  int i_score = 0;
  while (fgets(line, MAX_LINE, score_file) != NULL) {

    // Replace the EOL.
    int length = strlen(line);
    if (length > 0) {
      line[length - 1] = '\0';
    }

    // Skip over the first columns.
    char* start_of_column = &(line[0]);
    int i_column = 0;
    while (i_column < score_column) {

      // Check for the start of a new column.
      if (*start_of_column == '\t') {
	i_column++;
      }

      // Move to the next letter in the string.
      start_of_column++;

      // Check for the end of the line.
      if (*start_of_column == '\0') {
	die("Couldn't find column %d in line %d of %s.\nline=%s\n",
	    score_column + 1, i_score, score_filename, line);
      }

    }

    // At the p-value column, read the value.
    ATYPE score;
    if (sscanf(start_of_column, "%lf", &score) != 1) {
      die("Can't parse line %d in %s.\n(column=%s)\n(line=%s)\n",
	  i_score, score_filename, start_of_column, line);
    }
    if (i_score >= MAX_SCORES) {
      die("Can't read more than %d p-values.\n", MAX_SCORES);
    }
    set_array_item(i_score, score, initial_scores);

    // If requested, store the line.
    if (input_lines != NULL) {
      add_string_with_score(line, input_lines, score);
    }

    i_score++;
  }
  fclose(score_file);
  int num_scores = i_score;
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Read %d values from %s.\n", num_scores,
	    score_filename);
  }

  // Make a smaller array to contain just the scores.
  ARRAY_T* scores = allocate_array(num_scores);
  for (i_score = 0; i_score < num_scores; i_score++) {
    set_array_item(i_score,
		   get_array_item(i_score, initial_scores),
		   scores);
  }
  free_array(initial_scores);

  return(scores);
}

/*****************************************************************************
 * Estimate pi_zero using a fixed lambda.
 *
 * This is step (b) in Algorithm 3.1 of Storey (2002).
 *
 * The estimated pi_zero is
 *
 *          W(\lambda)
 *        ---------------
 *        (1 - \lambda) m
 *
 * where W(\lambda) is the number of p-values > lambda, and m is the
 * total number of p-values.
 *****************************************************************************/
static double estimate_pi_zero_fixed_lambda
  (double   lambda,
   ARRAY_T* p_values)
{
  int num_pvalues = get_array_length(p_values);

  // Find the index of the first p-value that is greater than lambda.
  int start;
  for (start = 0; start < num_pvalues; start++) {
    if (get_array_item(start, p_values) >= lambda) {
      break;
    }
  }
  int num_greater = num_pvalues - start;

  // Compute pi_zero.
  double pi_zero = (double)num_greater /
    ((1.0 - lambda) * (double)num_pvalues);

  // If no p-values exceeded the given lambda, we have no estimate for pi_zero.
  // In this case, substitute a dummy huge value.
  if (pi_zero == 0.0) {
    pi_zero = HUGE_VAL;
  }

  return(pi_zero);
}


/*****************************************************************************
 * Given a set of p-values, estimate the percentage of the empirical
 * distribution that is drawn according to the null hypothesis.
 *
 * This code implements the method described in John Storey, "A direct
 * approach to false discovery rates."  Journal of the Royal
 * Statistical Society B, 2002.
 *
 * Following is an illuminating excerpt from an email from Tim Bailey
 * regarding this estimation procedure.
 *
 * "The paper shows that the estimate of pi_0 will always be biased high.
 *
 * "A careful reading also will convince one that E[pFDR] in the
 * algorithm is proportional to pi_0, for a fixed lambda.  (This
 * explains why your algorithm works only in terms of pi_0, whereas
 * the one in the paper works in terms of E[pFDR].)
 *
 * "So, the algorithm uses the minimum pi_0(lambda) observed for any
 * of the tested lambda as the "plug-in estimate" for pi_0.  It then
 * picks the lambda that gives the lowest MSE w/r to the plug-in
 * estimate under bootstrap sampling.  Since the plug-in estimate will
 * be biased high, so will the final, bootstrapped estimate of pi_0.
 *
 * "The last bit of explanation is a bit of hand-waving--I don't
 * really understand the rationale behind the plug-in estimate trick
 * in any depth."
 *****************************************************************************/
#define MIN_SAMPLES 100            // Minimum number of p-values required.

double estimate_pi_zero
  (char*    pi_zero_filename,      // Name of file to store intermediate values.
   int      num_bootstraps,        // How many bootstraps to perform.
   int      num_bootstrap_samples, // Number of p-values in each bootstrap.
   int      num_lambda,            // Number of lambda values to consider.
   float    max_lambda,            // Maximum lambda value to consider.
   ARRAY_T* p_values)
{
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Estimating pi_0.\n");
  }
  
  int num_pvalues = get_array_length(p_values);

  // Sort the p-values in ascending order.
  sort_array(FALSE, p_values);

  // Allocate various arrays.
  ARRAY_T* pi_zeroes = allocate_array(num_lambda);
  ARRAY_T* mean_squared_errors = allocate_array(num_lambda);

  // Step through the requested lambdas.
  int i_lambda;
  for (i_lambda = 0; i_lambda < num_lambda; i_lambda++) {
    double lambda = ((double)(i_lambda + 1) / (double)num_lambda) * max_lambda;

    // Get a simple estimate of pi_zero with a fixed lambda.
    double pi_zero = estimate_pi_zero_fixed_lambda(lambda, p_values);
    set_array_item(i_lambda, pi_zero, pi_zeroes);
  }

  // Find the minimal pi_zero.
  double min_pi_zero = get_array_minimum(pi_zeroes);
  if (verbosity > NORMAL_VERBOSE) {
    fprintf(stderr, "Minimal pi_zero = %g\n", min_pi_zero);
  }

  // Choose the lambda which minimizes the mean-squared error of our estimates.
  if (verbosity >= DUMP_VERBOSE) {
    fprintf(stderr, "Performing %d bootstraps.\n", num_bootstraps);
  }
  int i_bootstrap;
  for (i_bootstrap = 0; i_bootstrap < num_bootstraps; i_bootstrap++) {
    if (verbosity >= DUMP_VERBOSE) {
      fprintf(stderr, "Bootstrap %d\n", i_bootstrap);
    }

    // Create an array of bootstrapped p-values.
    ARRAY_T* bootstrapped_pvalues = bootstrap_array(p_values,
						    num_bootstrap_samples);
    sort_array(FALSE, bootstrapped_pvalues);

    for (i_lambda = 0; i_lambda < num_lambda; i_lambda++) {
      double lambda = ((double)(i_lambda + 1)
		       / (double)num_lambda) * max_lambda;

      // Get a simple estimate of pi_zero with a fixed lambda.
      double pi_zero = estimate_pi_zero_fixed_lambda(lambda,
						     bootstrapped_pvalues);

      // Estimated mean-squared error.
      // See Equation (28) in Storey (2002).
      double error = pi_zero - min_pi_zero;
      double squared_error = error * error;
      incr_array_item(i_lambda, squared_error, mean_squared_errors);
    }
    free_array(bootstrapped_pvalues);
  }

  // Find the lambda with the minimal error.
  int best_lambda_index = 0;
  double minimal_error = get_array_item(0, mean_squared_errors);
  for (i_lambda = 1; i_lambda < num_lambda; i_lambda++) {
    double this_error = get_array_item(i_lambda, mean_squared_errors);
    if (this_error < minimal_error) {
      best_lambda_index = i_lambda;
      minimal_error = this_error;
    }
  }
  if (verbosity >= DUMP_VERBOSE) {
    fprintf(stderr, "Minimal error = %g\n", minimal_error);
    fprintf(stderr, "Index of best lambda = %d\n", best_lambda_index);
  }

  // Use the corresponding pi_zero.
  double pi_zero = get_array_item(best_lambda_index, pi_zeroes);

  // Ensure the pi_zero estimate is in the right range.
  pi_zero = MAX(MIN(pi_zero, 1.0), 0.0);

  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Estimated pi_0=%g\n", pi_zero);
  }

  // If requested, store the estimated pi-zero as a function of lambda.
  if (pi_zero_filename != NULL) {
    FILE* pi_zero_file = NULL;
    if (open_file(pi_zero_filename, "w", TRUE, "pi-zero", "pi-zero",
		&pi_zero_file) == 0) {
      exit(1);
    }
    fprintf(pi_zero_file, "p-value threshold\tlocal pi-zero\tfinal pi-zero\n");
    for (i_lambda = 0; i_lambda < num_lambda; i_lambda++) {
      double lambda = ((double)(i_lambda + 1) 
		       / (double)num_lambda) * max_lambda;
      fprintf(pi_zero_file, "%g\t%g\t%g\n", 
	      lambda,
	      get_array_item(i_lambda, pi_zeroes),
	      pi_zero);
    }
    fclose(pi_zero_file);
  }

  // Free local dynamic memory.
  free_array(pi_zeroes);
  free_array(mean_squared_errors);
  return(pi_zero);
}

/*************************************************************************
 * Convert a set of p-values to a set of q-values.  The p-values must
 * be in sorted (ascending) order.  This function replaces them with
 * the corresponding q-values.
 *
 * In order to estimate pi0 we need to have the distribution of pvalues.
 * We can get the distribution from the pvalues array only if it contains
 * all the observed pvalues. If we don't have all the observed p-values
 * we can use the uniformly sampled subset of pvalues stored in the reservoir 
 * sampler. If neither the full set of pvalues, or uniformly sampled pvalues
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
  long      total_values,     // Total number of p-values observed.
  ARRAY_T*  pvalues,         // retained pvalues
  ARRAY_T*  sampled_pvalues  // uniformly sampled pvalues
) {

  // Make sure we got at least one p-value.
  int num_values = get_array_length(pvalues);
  if (num_values == 0) {
    return;
  }

  // Verify that the p-values are in the right range.
  int i_pvalue;
  for (i_pvalue = 0; i_pvalue < num_values; i_pvalue++) {
    ATYPE this_pvalue = get_array_item(i_pvalue, pvalues);
    myassert(1, this_pvalue >= 0.0, "Negative p-value (%g)\n", this_pvalue);
    myassert(1, this_pvalue <= 1.0, "p-value > 1 (%g)\n", this_pvalue);
  }

  // Verify that the p-values are in sorted order.
  myassert(1, is_sorted(TRUE, pvalues), "The p-values are not sorted.");

  double pi_zero = 1.0;
  
  // Estimate pi_zero, if requested.
  if (use_pi_zero) {

    BOOLEAN_T have_all_pvalues = (total_values == (long) num_values);

    if (have_all_pvalues) {
      if (num_values > (long) MIN_SAMPLES) {
        // Estimate pi0 from complete set of pvalues
        if (verbosity >= NORMAL_VERBOSE) {
          fprintf(
            stderr, 
            "Estimating pi_0 from all %d observed p-values.\n",
            num_values
          );
        }
        pi_zero = estimate_pi_zero(
           pi_zero_filename,
           num_bootstraps,
           num_bootstrap_samples,
           num_lambda,
           max_lambda,
           pvalues
        );
      } 
      else {
        if (verbosity >= NORMAL_VERBOSE) {
          fprintf(
            stderr, 
            "Warning: Cannot estimate pi_0 accurately from fewer than %d p-values."
            "         total p-values = %d. Using pi_zero = 1.0.\n",
            MIN_SAMPLES,
            num_values
          );
        }
      }
    } 
    else if (sampled_pvalues != NULL) {
      int num_samples = get_array_length(sampled_pvalues);
      if (num_samples > MIN_SAMPLES) {
        // Estimate pi0 from uniformly sampled pvalues
        if (verbosity >= NORMAL_VERBOSE) {
          fprintf(
            stderr,
            "Estimating pi_0 from a uniformly sampled set "
            "of %d p-values.\n",
            num_samples
          );
        }
        pi_zero = estimate_pi_zero(
           pi_zero_filename,
           num_bootstraps,
           num_bootstrap_samples,
           num_lambda,
           max_lambda,
           sampled_pvalues
        );
      } 
      else {
        fprintf(
          stderr, 
          "Warning: Cannot estimate pi_0 accurately from fewer "
          "than %d sampled p-values.\n"
          "         only %d p-values were sampled. Using pi_zero = 1.0.\n",
          MIN_SAMPLES,
          num_samples
        );
      }
    }
    else {
      fprintf(
        stderr, 
        "Warning: Cannot estimate pi_zero without the complete range of p-values,\n"
        "         or a uniformly sampled subset of pvalues. Using pi_zero = 1.0.\n"
      );
    }
  }

  // Traverse p-values from lowest to highest.
  for (i_pvalue = 0; i_pvalue < num_values; i_pvalue++) {
    ATYPE this_pvalue = get_array_item(i_pvalue, pvalues);

    ATYPE this_fdr = (pi_zero * this_pvalue * (float)total_values)
      / (float)(i_pvalue + 1);

    // Can't have an FDR > 1.0.
    set_array_item(i_pvalue, MIN(1.0, this_fdr), pvalues);
  }

  // No conversion if we only got one value or if only computing FDR.
  if ((num_values == 1) || (compute_fdr)) {
    return;
  }

  // Convert to q-values.
  ATYPE prev_pvalue = get_array_item(num_values - 1, pvalues);
  for (i_pvalue = num_values - 2; i_pvalue >= 0; i_pvalue--) {
    ATYPE this_pvalue = get_array_item(i_pvalue, pvalues);
    if (prev_pvalue < this_pvalue) {
      set_array_item(i_pvalue, prev_pvalue, pvalues);
    }
    prev_pvalue = get_array_item(i_pvalue, pvalues);
  }
}

/*****************************************************************************
 * Convert scores to p-values using an empirical null.
 *
 * The p-value associated with an observed score x is simply
 * (a+1)/(b+1), where a is the number of null scores better than x,
 * and b is the total number of null scores.
 *
 * The +1 correction is recommended in the article, "A Note on the
 * Calculation of Empirical P Values from Monte Carlo Procedures."
 * http://www.pubmedcentral.nih.gov/articlerender.fcgi?artid=379178
 *****************************************************************************/
void convert_scores_to_pvalues
  (BOOLEAN_T good_score_is_low,
   ARRAY_T*  observed_scores,
   ARRAY_T*  null_scores)
{
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Converting scores to p-values.\n");
  }

  // Sort both arrays in order, best to worst.
  sort_array(good_score_is_low != TRUE, observed_scores);
  sort_array(good_score_is_low != TRUE, null_scores);

  int num_observed = get_array_length(observed_scores);
  int num_null = get_array_length(null_scores);

  // Get the best null score.
  int i_null = 0;
  ATYPE null_score = get_array_item(i_null, null_scores);

  // Traverse the observed scores, best to worst.
  int i_observed;
  for (i_observed = 0; i_observed < num_observed; i_observed++) {
    ATYPE observed_score = get_array_item(i_observed, observed_scores);

    // Find the first null score that is worse than the observed score.
    while ((((good_score_is_low == TRUE) && (observed_score > null_score)) ||
	    ((good_score_is_low == FALSE) && (observed_score < null_score))) &&
	   (i_null < num_null)) {
      i_null++;
      if (i_null < num_null) {
	null_score = get_array_item(i_null, null_scores);
      }
    }

    ATYPE pvalue = (ATYPE)(i_null + 1) / (ATYPE)(num_null + 1);
    set_array_item(i_observed, pvalue, observed_scores);
  }
}

/*****************************************************************************
 * MAIN
 *****************************************************************************/
#ifdef MAIN

#include "simple-getopt.h"
VERBOSE_T verbosity = INVALID_VERBOSE;

int main
  (int    argc,
   char * argv[])
{

  // Default parameter settings.
  char*     null_filename = NULL; // By default, assume p-value input.
  BOOLEAN_T good_score_is_low = TRUE; // By default, good score is low.
  BOOLEAN_T compute_fdr = FALSE;  // By default, compute p-value.
  BOOLEAN_T use_pi_zero = FALSE; // By default, compute q-value, not just FDR.
  char*     pi_zero_filename = NULL; //  Store pi-zero in the given file.
  int       num_bootstraps = NUM_BOOTSTRAPS;
  int       num_header_lines = 0;
  int       pvalue_column = 0; // Indexed from 0.
  BOOLEAN_T append_output = FALSE; // Append q-values to end of input line?
  long      seed = time(0);
  verbosity = NORMAL_VERBOSE;

  const int num_options = 11;
  cmdoption const options[] = {
    { "null", REQUIRED_VALUE },
    { "good-score", REQUIRED_VALUE },
    { "pi-zero", NO_VALUE },
    { "pi-zero-file", REQUIRED_VALUE },
    { "fdr", NO_VALUE},
    { "bootstraps", REQUIRED_VALUE },
    { "header", REQUIRED_VALUE },
    { "column", REQUIRED_VALUE },
    { "append", NO_VALUE },
    { "seed", REQUIRED_VALUE },
    { "verbosity", REQUIRED_VALUE }
  };

  // Define the usage message.
  char      usage[400] = "";
  strcat(usage, "USAGE: qvalue [options] <pvalues>\n");
  strcat(usage, "\n");
  strcat(usage, "   Options:\n");
  strcat(usage, "     --null <file>\n");
  strcat(usage, "     --good-score high|low\n");
  strcat(usage, "     --pi-zero\n");
  strcat(usage, "     --pi-zero-file <file>\n");
  strcat(usage, "     --fdr\n");
  strcat(usage, "     --bootstraps <int> (default=1000)\n");
  strcat(usage, "     --header <int> (default=0)\n");
  strcat(usage, "     --column <int> (default=1)\n");
  strcat(usage, "     --append\n");
  strcat(usage, "     --seed <int> (default from clock)\n");
  strcat(usage, "     --verbosity 1|2|3|4 (default = 2)\n");
  strcat(usage, "\n");

  // Parse the command line.
  int option_index = 0;
  char* option_name = NULL;
  char* option_value = NULL;
  const char *  message = NULL;
  simple_setopt(argc, argv, num_options, options);
  while(1) {
    // Read the next option, and break if we're done.
    int c = simple_getopt(&option_name, &option_value, &option_index);
    if (c == 0) {
      break;
    } else if (c < 0) {
      simple_getopterror(&message);
      die("Error processing command line options (%s)\n", message);
    }

    if (strcmp(option_name, "null") == 0) {
      null_filename = option_value;
    } else if (strcmp(option_name, "good-score") == 0) {
      if (strcmp(option_value, "high") == 0) {
	good_score_is_low = FALSE;
      } else if (strcmp(option_value, "low") == 0) {
	good_score_is_low = TRUE;
      } else {
	fprintf(stderr, "Invalid option (--good-score %s)\n",
		option_value);
	exit(1);
      }
    } else if (strcmp(option_name, "pi-zero") == 0) {
      use_pi_zero = TRUE;
    } else if (strcmp(option_name, "pi-zero-file") == 0) {
      use_pi_zero = TRUE;
      pi_zero_filename = option_value;
    } else if (strcmp(option_name, "fdr") == 0) {
      compute_fdr = TRUE;
    } else if (strcmp(option_name, "bootstraps") == 0) {
      num_bootstraps = atoi(option_value);
    } else if (strcmp(option_name, "header") == 0) {
      num_header_lines = atoi(option_value);
    } else if (strcmp(option_name, "column") == 0) {
      pvalue_column = atoi(option_value) - 1; // User inputs indexed from 1.
    } else if (strcmp(option_name, "append") == 0) {
      append_output = TRUE;
    } else if (strcmp(option_name, "seed") == 0) {
      seed = atoi(option_value);
    } else if (strcmp(option_name, "verbosity") == 0) {
      verbosity = atoi(option_value);
    }
  }

  // Read the single required argument.
  if (option_index + 1 != argc) {
    fprintf(stderr, "%s", usage);
    exit(1);
  }
  char* pvalue_filename = argv[option_index];

  // Initialize the random number generator.
  my_srand(seed);

  // Make a string list to store the file contents.
  STRING_LIST_T* input_lines = NULL;
  STRING_LIST_T* header_lines = NULL;
  if (append_output) {
    input_lines = new_string_list();
    if (num_header_lines != 0) {
      header_lines = new_string_list();
    }
  }

  // Read the p-values from the specified column.
  ARRAY_T* pvalues = read_scores_from_column(pvalue_filename,
					     num_header_lines,
					     pvalue_column,
					     header_lines,
					     input_lines);
  int num_pvalues = get_array_length(pvalues);

  // If an empirical null was given, use it.
  ARRAY_T* raw_scores = NULL;
  if (null_filename != NULL) {

    // Store a copy of the raw scores for later reference, best score first.
    raw_scores = allocate_array(num_pvalues);
    copy_array(pvalues, raw_scores);
    sort_array(good_score_is_low != TRUE, raw_scores);

    // Read the null scores.
    ARRAY_T* null_scores = read_scores_from_column(null_filename,
						   num_header_lines,
						   pvalue_column,
						   NULL, // Do not store
						   NULL); // Do not store

    // Compute p-values.
    convert_scores_to_pvalues(good_score_is_low, pvalues, null_scores);
  }

  // Sort the p-values, if necessary.
  if (!is_sorted(TRUE, pvalues)) {
    sort_array(FALSE, pvalues);
    if (append_output == TRUE) {
      sort_string_list_by_score(input_lines, good_score_is_low != TRUE);
    }
  }

  // Make a copy of the p-values.
  ARRAY_T* qvalues = allocate_array(num_pvalues);
  copy_array(pvalues, qvalues);

  // Compute q-values.
  compute_qvalues(
      compute_fdr,
		  use_pi_zero,
		  pi_zero_filename,
		  num_bootstraps,
		  NUM_BOOTSTRAP_SAMPLES, //get_array_length(qvalues),
		  NUM_LAMBDA,
		  MAX_LAMBDA,
		  get_array_length(qvalues),
		  qvalues,
      NULL
    );

  // Print header.
  if (append_output) {
    int i_header = 0;
    for (i_header = 0; i_header < num_header_lines; i_header++) {
      if (raw_scores != NULL) {
	printf("%s\tp-value\tq-value\n",
	       get_nth_string(i_header, header_lines));
      } else {
	printf("%s\tq-value\n",
	       get_nth_string(i_header, header_lines));
      }
    }
  }

  // Print to stdout.
  int i_pvalue;
  for (i_pvalue = 0; i_pvalue < num_pvalues; i_pvalue++) {
    if (append_output) {
      if (raw_scores != NULL) {
	printf("%s\t%g\t%g\n",
	       get_nth_string(i_pvalue, input_lines),
	       get_array_item(i_pvalue, pvalues),
	       get_array_item(i_pvalue, qvalues));
      } else {
	printf("%s\t%g\n",
	       get_nth_string(i_pvalue, input_lines),
	       get_array_item(i_pvalue, qvalues));
      }
    } else if (raw_scores != NULL) {
      printf("%g\t%g\t%g\n",
	     get_array_item(i_pvalue, raw_scores),
	     get_array_item(i_pvalue, pvalues),
	     get_array_item(i_pvalue, qvalues));
    } else {
      printf("%g\t%g\n",
	     get_array_item(i_pvalue, pvalues),
	     get_array_item(i_pvalue, qvalues));
    }
  }
  free_array(pvalues);
  free_array(qvalues);
  free_string_list(input_lines);

  return(0);
}

#endif
