/*************************************************************************
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 15 September 2009
 * PROJECT: MEME
 * COPYRIGHT: 2009, UW
 * DESCRIPTION: q value calculation from an empirical null distribution
 *************************************************************************/
#include "utils.h"
#include "qvalue.h"
#include "array.h"
#include "string-list.h"
#include "empirical-qvalue.h"
#include <time.h>
#include <assert.h>
#include <stdio.h>

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
   ARRAY_T*  observed_scores)
{
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Computing FDR estimates.\n");
  }

  // Make sure we got at least one of each type of score.
  int num_observed = get_array_length(observed_scores);
  int num_null = get_array_length(null_scores);
  if ((num_observed == 0) || (num_null == 0)) {
    return;
  }

  /* If we don't have the same number of observed and null scores, include
     this multiplicative factor in each FDR estimate. */
  ATYPE null_observed_factor = (ATYPE)num_observed / (ATYPE)num_null;

  // Verify that the scores are in sorted order.
  myassert(1, is_sorted(good_score_is_low, observed_scores), 
	   "Unsorted observed scores.");
  myassert(1, is_sorted(good_score_is_low, null_scores),
	   "Unsorted null scores.");

  // Traverse the list of observed scores.
  int i_null = 0;
  ATYPE current_null = get_array_item(i_null, null_scores);
  int i_observed;
  for (i_observed = 0; i_observed < num_observed; i_observed++) {
    ATYPE threshold = get_array_item(i_observed, observed_scores);

    // Check for a block of tied scores.
    int end_of_tied_block = i_observed + 1;
    while ((end_of_tied_block < num_observed) &&
	   (get_array_item(end_of_tied_block, observed_scores) == threshold)) {
      //      fprintf(stderr, "Found a tie at %d.\n", end_of_tied_block);
      end_of_tied_block++;
    }

    // Find a null score that is greater than this one.
    while ((current_null >= threshold) && (i_null < num_null)) {
      i_null++;
      current_null = get_array_item(i_null, null_scores);
    }

    ATYPE fdr = null_observed_factor 
      * ((ATYPE)(i_null) / (ATYPE)(end_of_tied_block));
    //    fprintf(stderr, "fdr at %g = %d / %d = %g\n", threshold,
    //    	    i_null, end_of_tied_block, fdr);

    // Can't have an FDR > 1.
    if (fdr > 1.0) fdr = 1.0;

    // Assign FDR to everyone in this block.
    int i;
    for (i = i_observed; i < end_of_tied_block; i++) {
      set_array_item(i, fdr, observed_scores);

      if ((((i + 1) % 100000) == 0) && (verbosity >= NORMAL_VERBOSE)) {
	fprintf(stderr, "%5.2f%% complete.\n",
		100.0 * ((float)i / num_observed));
      }
    }

    // Skip to the end of the block.
    i_observed = end_of_tied_block - 1;
  }

  // No conversion if we only got one value or if only computing FDR.
  if ((num_observed == 1) || (compute_fdr)) {
    return;
  }

  // Convert to q-values.
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Converting FDRs to q-values.\n");
  }
  ATYPE prev_fdr = get_array_item(num_observed - 1, observed_scores);
  for (i_observed = num_observed - 2; i_observed >= 0; i_observed--) {
    ATYPE this_fdr = get_array_item(i_observed, observed_scores);
    if (prev_fdr < this_fdr) {
      set_array_item(i_observed, prev_fdr, observed_scores);
    }
    prev_fdr = get_array_item(i_observed, observed_scores);
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
  BOOLEAN_T good_score_is_low = FALSE; // By default, good score is high.
  BOOLEAN_T compute_fdr = FALSE;  // By default, compute q-value.
  int       num_header_lines = 0;
  int       score_column = 0; // Indexed from 0.
  BOOLEAN_T append_output = FALSE; // Append q-values to end of input line?
  long      seed = time(0);
  verbosity = NORMAL_VERBOSE;

  const int num_options = 7;
  cmdoption const options[] = {
    { "good-score", REQUIRED_VALUE },
    { "fdr", NO_VALUE},
    { "header", REQUIRED_VALUE },
    { "column", REQUIRED_VALUE },
    { "append", NO_VALUE },
    { "seed", REQUIRED_VALUE },
    { "verbosity", REQUIRED_VALUE }
  };

  // Define the usage message.
  char      usage[400] = "";
  strcat(usage, "USAGE: empirical-qvalue [options] <observed scores> <null scores>\n");
  strcat(usage, "\n");
  strcat(usage, "   Options:\n");
  strcat(usage, "     --good-score high|low (default=high)\n");
  strcat(usage, "     --fdr\n");
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

    if (strcmp(option_name, "good-score") == 0) {
      if (strcmp(option_value, "high") == 0) {
	good_score_is_low = FALSE;
      } else if (strcmp(option_value, "low") == 0) {
	good_score_is_low = TRUE;
      } else {
	fprintf(stderr, "Invalid option (--good-score %s)\n",
		option_value);
	exit(1);
      }
    } else if (strcmp(option_name, "fdr") == 0) {
      compute_fdr = TRUE;
    } else if (strcmp(option_name, "header") == 0) {
      num_header_lines = atoi(option_value);
    } else if (strcmp(option_name, "column") == 0) {
      score_column = atoi(option_value) - 1; // User inputs indexed from 1.
    } else if (strcmp(option_name, "append") == 0) {
      append_output = TRUE;
    } else if (strcmp(option_name, "seed") == 0) {
      seed = atoi(option_value);
    } else if (strcmp(option_name, "verbosity") == 0) {
      verbosity = atoi(option_value);
    }
  }

  // Read the two required arguments.
  if (option_index + 2 != argc) {
    fprintf(stderr, usage);
    exit(1);
  }
  char* observed_filename = argv[option_index];
  char* null_filename = argv[option_index + 1];

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

  // Read the two input files.
  ARRAY_T* observed_scores = read_scores_from_column(observed_filename,
						     num_header_lines,
						     score_column,
						     header_lines,
						     input_lines);
  int num_observed = get_array_length(observed_scores);

  ARRAY_T* null_scores = read_scores_from_column(null_filename,
						 num_header_lines,
						 score_column,
						 NULL, // Do not store
						 NULL); // Do not store

  // Sort the scores, if necessary.
  if (!is_sorted(good_score_is_low, observed_scores)) {
    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(stderr, "Sorting observed scores.\n");
    }
    sort_array(!good_score_is_low, observed_scores);
    if (append_output == TRUE) {
      sort_string_list_by_score(input_lines, !good_score_is_low);
    }
  }
  if (!is_sorted(good_score_is_low, null_scores)) {
    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(stderr, "Sorting null scores.\n");
    }
    sort_array(!good_score_is_low, null_scores);
  }

  // Make a copy of the observed scores.
  ARRAY_T* qvalues = allocate_array(num_observed);
  copy_array(observed_scores, qvalues);

  // Compute q-values.
  sort_array(!good_score_is_low, observed_scores);
  compute_qvalues_empirical(compute_fdr, good_score_is_low, null_scores,
			    qvalues);
  free_array(null_scores);

  // Print header.
  if (append_output) {
    int i_header = 0;
    for (i_header = 0; i_header < num_header_lines; i_header++) {
      printf("%s\tq-value\n",
	     get_nth_string(i_header, header_lines));
    }
  }

  // Print to stdout.
  int i_score;
  for (i_score = 0; i_score < num_observed; i_score++) {
    if (append_output) {
      printf("%s\t%g\n",
	     get_nth_string(i_score, input_lines),
	     get_array_item(i_score, qvalues));
    } else {
      printf("%g\t%g\n",
	     get_array_item(i_score, observed_scores),
	     get_array_item(i_score, qvalues));
    }
  }
  free_array(observed_scores);
  free_array(qvalues);
  free_string_list(input_lines);

  return(0);
}

#endif
