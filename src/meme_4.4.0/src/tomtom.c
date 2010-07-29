/**************************************************************************
 * FILE: tomtom.c
 * AUTHOR: Shobhit Gupta, Timothy Bailey, William Stafford Noble
 * CREATE DATE: 02-21-06
 * PROJECT: TOMTOM
 * DESCRIPTION: Motif-Motif comparison using the score statistic.
 **************************************************************************/
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include "utils.h"       // Generic utilities.
#include "motif.h"       // Needed for complement_dna_freqs.
#include "mhmm-state.h"  // HMM.
#include "matrix.h"      // Routines for floating point matrices.
#include "array.h"       // Routines for floating point arrays.
#include "alphabet.h"    // The alphabet.
#include "fasta-io.h"
#include "meme-io.h"
#include "simple-getopt.h"
#include "pssm.h"
#include "fitevd.h"
#include "float.h"
#include "config.h"
#include "io.h"
#include "ceqlogo.h"
#include "qvalue.h"
#include "string-list.h"
#include "object-list.h"
#include "tomtom.h"

#define DEFAULT_QUERY_PSEUDO   0.0
#define DEFAULT_TARGET_PSEUDO   0.0
#define MAX_TARGETS 5000
#define MAX_QUERYS 5000
#define BINS 100
#define ITERATIONS 1000
#define LOGOHEIGHT 10		// height of an alignment of logos in cm.
#define MIN_TARGET_DATABASE_SIZE 50

/* Define some string limits (these include the EOL).*/
// FIXME: these are brittle and will cause bugs if motifs are too long
#define MAX_MOTIF_LENGTH 100
#define MAX_LINE_LENGTH 100

/**************************************************************************
 * An object for storing information about one query-target alignment.
 **************************************************************************/
typedef struct tomtom_match {
  char  query_id[MAX_MOTIF_ID_LENGTH];
  char  target_id[MAX_MOTIF_ID_LENGTH];
  char  query_consensus[MAX_MOTIF_LENGTH];
  char  target_consensus[MAX_MOTIF_LENGTH];
  char  orientation;
  int   offset;
  int   overlap;
  float pvalue;
  float evalue;
  float qvalue;
  int   query_index;
  int   target_index;
  int   offset_by_length;
} TOMTOM_MATCH_T;

/**************************************************************************
 * Create a new query-target alignment object.
 * Initially, there is no q-value and no orientation.
 **************************************************************************/
TOMTOM_MATCH_T* new_tomtom_match
  (char*   query_id,
   char*   target_id,
   char*   query_consensus,
   char*   target_consensus,
   int     offset,
   int     overlap,
   float   pvalue,
   float   evalue,
   int     query_index,
   int     target_index,
   int     offset_by_length)
{
  TOMTOM_MATCH_T* return_value = mm_malloc(sizeof(TOMTOM_MATCH_T));

  strcpy(return_value->query_id, query_id);
  strcpy(return_value->target_id, target_id);
  strcpy(return_value->query_consensus, query_consensus);
  strcpy(return_value->target_consensus, target_consensus);
  return_value->offset = offset;
  return_value->overlap = overlap;
  return_value->pvalue = pvalue;
  return_value->evalue = evalue;
  return_value->query_index = query_index;
  return_value->target_index = target_index;
  return_value->offset_by_length = offset_by_length;
  return(return_value);
}

/**************************************************************************
 *
 * Creates a logo of the aligned target and query motif, writes it to an
 * EPS file and calls convert to convert the EPS file to a PNG image.
 *
 * Returns name of logo file (minus the extension).
 *
 **************************************************************************/
static char *create_logo(
  MOTIF_T* target,             // Target motif
  MOTIF_T* query,              // Query motif
  int offset,                  // Offset of target rel. to query
  BOOLEAN_T ssc, 	       // Small sample correction
  char *output_dirname 	       // Output directory
) 
{
  char *logo_filename, *logo_path;
  double logo_height = LOGOHEIGHT;
  double logo_width = offset <= 0 ? MAX(query->length, target->length)-offset : 
    MAX(query->length+offset, target->length);

  // create name of logo file
  logo_filename = (char*) malloc(
    sizeof(char) * (5+strlen(query->id)+1+strlen(target->id)+1)
  );
  sprintf(logo_filename, "logo_%s_%s", query->id, target->id);

  // create output file path
  logo_path = make_path_to_file(output_dirname, logo_filename);

  CL_create2(
    target, 		// first motif
    target->id+1,	// title (skip +/- character)
    query, 		// second motif
    query->id, 		// title
    TRUE,               // error bars
    ssc,                // small sample correction
    logo_height,        // logo height (cm)
    logo_width,         // logo width (cm)
    get_alphabet(FALSE),// alphabet
    -offset, 		// offset of second motif
    logo_path,		// output file path
    "Tomtom"		// program name
  );
  myfree(logo_path);

  return(logo_filename);
}


/**************************************************************************
*
* get_pv_lookup_new
*
* Create a lookup table for the pv.
*
* FIXME: document what this really does!
*
**************************************************************************/
void get_pv_lookup_new
  (
   MATRIX_T* pssm_matrix,		// The PSSM.
   int alen,				// Length of alphabet.
   int range,				// Range of scores.
   ARRAY_T* background,			// Background model (uniform).
   MATRIX_T* reference_matrix, 		// reference_matrix
   MATRIX_T* pv_lookup_matrix 		// pv_lookup_matrix
   ) 
{

  int i, j, k;
  int w = get_num_rows(pssm_matrix);	// Width of PSSM.
  int size = w*range+1;
  int start;				// Starting query motif column
  int address_index = 0;		// next free position in matrix
  ARRAY_T* pdf_old = allocate_array(size);
  ARRAY_T* pdf_new = allocate_array(size);

  // Compute pv tables for groups of motif columns [start...]
  for (start=w-1; start>=0; start--) {

    // Compute the pdf recursively.
    init_array(0, pdf_new);
    set_array_item(0, 1, pdf_new);    // Prob(0)
    for (i=0; start+i<w; i++) {
      int max = i * range;
      SWAP(ARRAY_T*, pdf_new, pdf_old)
      for (k=0; k<=max+range; k++) set_array_item(k, 0, pdf_new);
      for (j=0; j<alen; j++) {
	int s = (int) get_matrix_cell(start+i, j, pssm_matrix);
	for(k=0; k<=max; k++) {
	  double old = get_array_item(k, pdf_old);
	  if (old != 0) {
	    double new = get_array_item(k+s, pdf_new) +
	      (old * get_array_item(j, background));
	    set_array_item(k+s, new, pdf_new);
	  } // old
	} // k
      } // j

      // Compute 1-cdf for motif consisting of columns [start, start+i]
      // This is the p-value lookup table for those columns.
      ARRAY_T* pv = allocate_array(get_array_length(pdf_new));
      copy_array(pdf_new, pv);
      int ii;
      for (ii=size-2; ii>=0; ii--) {
	double p = get_array_item(ii, pv) + get_array_item(ii+1, pv);
	set_array_item(ii, MIN(1.0, p), pv);
      }

      // copy the pv lookup table into the lookup matrix
      set_matrix_row(address_index, pv, pv_lookup_matrix);
      free_array(pv);

      // Store the location of this pv lookup table in the reference array
      set_matrix_cell(start,
		      start+i,
		      (double) address_index++,
		      reference_matrix);
    } // i
  } // start position

  // Free space.
  free_array(pdf_new);
  free_array(pdf_old);
} // get_pv_lookup_new

/**************************************************************************
 * Parse a list of Tomtom output strings, converting p-values to q-values.
 **************************************************************************/
static void convert_tomtom_p_to_q
  (OBJECT_LIST_T* output_list)
{
  int i_line;
  STRING_LIST_T* completed_ids;

  // Extract all of the p-values.
  ARRAY_T* pvalues = get_object_list_scores(output_list);

  // Convert p-values to q-values.
  compute_qvalues(FALSE, // Don't stop with FDR.
		  TRUE,  // Estimate pi-zero.
		  NULL, // Don't store pi-zero in a file.
		  NUM_BOOTSTRAPS,
		  NUM_BOOTSTRAP_SAMPLES,
		  NUM_LAMBDA,
		  MAX_LAMBDA,
		  get_array_length(pvalues), pvalues);
  ARRAY_T* qvalues = pvalues;

  // Traverse the list of matches.
  TOMTOM_MATCH_T* my_match
    = (TOMTOM_MATCH_T*)retrieve_next_object(output_list);
  int index = 0;
  while (my_match != NULL) {

    // Store the q-value.
    my_match->qvalue = get_array_item(index, qvalues);

    // Get the next match.
    my_match
      = (TOMTOM_MATCH_T*)retrieve_next_object(output_list);
    index++;
  }

} // Function convert_tomtom_p_to_q


/**************************************************************************
 * Given a sorted list of output lines corresponding to one query
 * motif, filter the list such that it contains only the best-scoring
 * hit to each target.
 *
 * As a side effect, re-format the lines slightly.
 **************************************************************************/
static void select_strand
  (OBJECT_LIST_T* output_list,
   OBJECT_LIST_T* final_output)
{
  STRING_LIST_T* completed_ids;
  completed_ids = new_string_list();

  int i_output = 0;
  
  // Traverse the list of matches.
  TOMTOM_MATCH_T* my_match
    = (TOMTOM_MATCH_T*)retrieve_next_object(output_list);
  int index = 0;
  while (my_match != NULL) {

    // Get the orientation of the target motif.
    my_match->orientation = '-';
    if (my_match->target_id[0] == '+') {
      my_match->orientation = '+';
    }

    // Make a new ID that removes the orientation.
    char new_target_id[MAX_MOTIF_ID_LENGTH+1];
    sprintf(new_target_id, "%s", my_match->target_id+1);
    strcpy(my_match->target_id, new_target_id);

    // Have we already seen a better scoring match?
    if (!have_string(new_target_id, completed_ids)) {

      // Make a copy on the heap.
      TOMTOM_MATCH_T* match_copy = mm_malloc(sizeof(TOMTOM_MATCH_T));
      *match_copy = *my_match;

      // Add this match to the new list.
      store_object((void*)match_copy, NULL, my_match->pvalue, final_output);

      // Remember that we've seen this target already.
      add_string(new_target_id, completed_ids);
    }

    // Get the next match.
    my_match
      = (TOMTOM_MATCH_T*)retrieve_next_object(output_list);
    index++;
  }

  /* Free memory */
  free_string_list(completed_ids);

} // Function select_strand

/**************************************************************************
 * Gets the consensus DNA sequence from a pssm.
 **************************************************************************/
static void get_cons(MATRIX_T* freqs, char** cons){
  int index_rows;
  int index_cols;
  char curr_cons[1000] = "";
  STRING_LIST_T* alph = new_string_list();

  if (get_num_cols(freqs) == 4) {
    add_string("A", alph);
    add_string("C", alph);
    add_string("G", alph);
    add_string("T", alph);
  } else {
    add_string("A", alph);
    add_string("C", alph);
    add_string("D", alph);
    add_string("E", alph);
    add_string("F", alph);
    add_string("G", alph);
    add_string("H", alph);
    add_string("I", alph);
    add_string("K", alph);
    add_string("L", alph);
    add_string("M", alph);
    add_string("N", alph);
    add_string("P", alph);
    add_string("Q", alph);
    add_string("R", alph);
    add_string("S", alph);
    add_string("T", alph);
    add_string("V", alph);
    add_string("W", alph);
    add_string("Y", alph);
  }

  for (index_rows = 0; index_rows < get_num_rows(freqs);
       index_rows++) {
    int curr_index;
    curr_index = -1;
    float max;
   max = 0;

    for (index_cols = 0; index_cols < get_num_cols(freqs);
   index_cols++) {
      if (get_matrix_cell(index_rows, index_cols, freqs) > max){
  max = get_matrix_cell(index_rows, index_cols, freqs);
  curr_index = index_cols;
      }
    }
    char * temp;
    temp = strcat(curr_cons, get_nth_string(curr_index, alph));
  }
  *cons = curr_cons;
  free_string_list(alph);
} /* Function get_cons */


/**************************************************************************
 * Computes the score for a particular configuration
 **************************************************************************/
static double compute_overlap_score 
  (MATRIX_T* pairwise_column_scores,
   int current_configuration,
   int *index_start,
   int *index_stop 
  )
{

  /* Initialize the columns from where the comparison needs to be started */
  int current_target_column = 0;
  int current_query_column = 0;
  if (current_configuration <= pairwise_column_scores->num_rows){
    current_target_column = 0;
    current_query_column = pairwise_column_scores->num_rows
      - current_configuration;
  }
  else{
    current_query_column = 0;
    current_target_column = current_configuration
      - pairwise_column_scores->num_rows;
  }

  // Sum the pairwise column scores.
  double score = 0;
  *index_start = current_query_column;
  do {
    score += get_matrix_cell(
                         current_query_column++, 
                         current_target_column++,
                         pairwise_column_scores);
  } while ((current_query_column < pairwise_column_scores->num_rows)
     &&(current_target_column < pairwise_column_scores->num_cols));
  *index_stop = current_query_column - 1;/* -1 for the last increment */

  return(score);
} /* Function compute_overlap_score */

/**************************************************************************
 * Computes the optimal offset and score corresponding to the smallest
 * pvalue.
 **************************************************************************/
static void compare_motifs (
   MATRIX_T* pairwise_column_scores,
   int*      optimal_offset,
   MATRIX_T* reference_matrix,
   MATRIX_T* pv_lookup_matrix,
   double*   optimal_pvalue,
   int *     overlapping_bases,
   BOOLEAN_T internal, // Is shorter motif contained in longer motif?
   int       min_overlap
   //int       address_index 
   ){

  int first_pass = 1;

  /*  Slide one profile over another */
  int total_configurations = pairwise_column_scores->num_rows
    + pairwise_column_scores->num_cols - 1;
  int current_configuration;
  int overlap = -1;
  for(current_configuration = 1;
      current_configuration <= total_configurations;
      current_configuration++){

    /* If requested, check whether one motif is contained within the other. */
    BOOLEAN_T process = TRUE;
    if (internal) {
      int max, min;
      max = pairwise_column_scores->num_rows;
      min = pairwise_column_scores->num_cols;
      if (max < min) {
	int temp;
	temp = min;
	min = max;
	max = temp;
      }
      if (!((current_configuration <= max)
	    && (current_configuration >= min))) {
	process = FALSE;
      }
    }

    /* Override min-overlap if query/target is smaller than min-overlap */
    int mo;
    mo = min_overlap;
    if (pairwise_column_scores->num_rows < mo) {
      mo = pairwise_column_scores->num_rows;
    }
    if (pairwise_column_scores->num_cols < mo) {
      mo = pairwise_column_scores->num_cols;
    }

    /* Check if this configuration should be processed (min_overlap). */
    if (mo > 1) {
      int max, min;
      min = mo;
      max = pairwise_column_scores->num_rows
	+ pairwise_column_scores->num_cols - mo;
      if (!((current_configuration <= max)
	    && (current_configuration >= min))) {
	process = FALSE;
      }
    }

    /* Check whether the requested percentage of information content is 
       contained within the matching portion of each motif. */
    /* FIXME: Stub to be filled in.  WSN 4 Apr 2008
    info_content = MIN(compute_info_content(start1, stop1, motif1),
		       compute_info_content(start2, stop2, motif2));
    if (info_content < target_info_content) {
      process = FALSE;
    }
    */

    if (process) {
      double pvalue = 1;
      double TEMP_SCORE = 0;
      int current_overlap = -1;
      int index_start, index_stop;

      /* Compute the score for the current configuration */
      TEMP_SCORE = compute_overlap_score(
					 pairwise_column_scores,
					 current_configuration,
				         &index_start,
					 &index_stop
					 );

      /* Compute the pvalue of the score */
      int curr_address_index = (int) get_matrix_cell(
                   index_start,
		   index_stop,
		   reference_matrix);
// FIXME: what does this code do?
      //if (curr_address_index == -1) {
//	curr_address_index = address_index - 1;
 //     }
      pvalue = get_matrix_cell(
        curr_address_index, (int) TEMP_SCORE, pv_lookup_matrix
      );

      current_overlap = index_stop - index_start + 1;
      if (current_overlap > 0) {
	/* Keep the configuration with the lowest pvalue.
	   If pvalue is the same, keep the largest overlap */
	if (first_pass == 1){
	  *optimal_offset = current_configuration;
	  *optimal_pvalue = pvalue;
	  overlap = current_overlap;
	  first_pass = 0;
	}
	else{
	  if (pvalue == *optimal_pvalue){
	    if (overlap < current_overlap){
	      *optimal_offset = current_configuration;
	      *optimal_pvalue = pvalue;
	      overlap = current_overlap;
	    }
	  }
	  else{
	    if(pvalue < *optimal_pvalue){
	      *optimal_offset = current_configuration;
	      *optimal_pvalue = pvalue;
	      overlap = current_overlap;
	    }
	  }
	}/* For the pvalue minimization condition */
      } /* Overlap > 0 (happens sometimes if min-overlap is overridden by
	   target motif length). */
    } /* Process */
  } /* Loop over configurations. */
  *overlapping_bases = overlap;
} /* Function compare_motifs */

/**************************************************************************
 * Computes ALLR scores for all columns of the query matrix against target
 * motif. This can then be used to by compute_allr.
 **************************************************************************/
void allr_scores
  (MOTIF_T* query_motif,
   MOTIF_T* target_motif,
   ARRAY_T* background,
   MATRIX_T** score)
{
  /*   Recurse pairwise over the columns */
  int index_query;
  int index_target;
  for (index_query = 0; index_query < query_motif->length;
       index_query++) {
    for (index_target = 0; index_target < target_motif->length;
	 index_target++) {

      double SCORE_F = 0;
      double SCORE_NUMERATOR = 0;
      double SCORE_DENOMINATOR = 0;
      int index; /* Index corresponds to the alphabets. */

      for (index = 0; index < query_motif->alph_size; index++) {
	/* The numerator of the ALLR */
	double nr_part1, nr_part2;
	/*  Part 1 of the numerator  */
	nr_part1
	  = (((target_motif->num_sites)*
	      (get_matrix_cell(index_target, index,
			       target_motif->freqs))) /* nb  for target */
	     * log((get_matrix_cell(index_query, index,
				    query_motif->freqs))
		   / get_array_item(index,
				    background))); /* Likelihood for query */
	/*  Part 2 of the numerator */
	nr_part2 
	  = (((query_motif->num_sites)*
	      (get_matrix_cell(index_query, index,
			       query_motif->freqs))) /* nb  for query */
	     * log((get_matrix_cell(index_target, index,
				    target_motif->freqs))
		   / get_array_item(index,
				    background))); /* Likelihood for target */
	/* Sum of the two parts. */
	SCORE_NUMERATOR += nr_part1 + nr_part2;
	/*  The denominoator of the ALLR */
	SCORE_DENOMINATOR
	  += (
	      ((query_motif->num_sites)
	       *  (get_matrix_cell(index_query,
				   index,
				   query_motif->freqs)))/* nb for query */
	      + ((target_motif->num_sites)
		 * (get_matrix_cell(index_target,
				    index,
				    target_motif->freqs)))/* nb for target */
	      ); /* Sum of nb for query and nb for target */

      } /* Nr and Dr of the score summed over the entire alphabet. */
      SCORE_F = SCORE_NUMERATOR / SCORE_DENOMINATOR;
      set_matrix_cell(index_query, index_target, SCORE_F, *score);

    } /* Target motif */
  } /* Query motif */
} /* Function ALLR scores */

/**************************************************************************
 * Computes euclidian distance between all columns of the query matrix
 * This score is not normalized on length and cannot be used without a
 * normalization criterion as our DP based p-value approach.
 **************************************************************************/
void ed_scores(MOTIF_T* query_motif,
	       MOTIF_T* target_motif,
	       ARRAY_T* background,
	       MATRIX_T** score) {

  /*   Recurse pairwise over the columns */
  int index_query;
  int index_target;
  for (index_query = 0; index_query < query_motif->length;
       index_query++) {
    for (index_target = 0; index_target < target_motif->length;
	 index_target++) {

      double SCORE_T = 0;

      int index; /* Index corresponds to the alphabets. */
      for (index = 0; index < query_motif->alph_size; index++) {

	double query_freq = get_matrix_cell(index_query,
					    index,
					    query_motif->freqs);
	double target_freq = get_matrix_cell(index_target,
					     index,
					     target_motif->freqs);
	/* Euclidean distance */
	double sq;
	sq = pow((query_freq - target_freq), 2);
	SCORE_T = SCORE_T + sq;
      }
      double SCORE_F;
      SCORE_F = - sqrt(SCORE_T);
      set_matrix_cell(index_query, index_target, SCORE_F, *score);
    } /* Target motif */
  } /* Query motif */
} /* Function ed_scores */

/**************************************************************************
 * Computes sandelin distance between all columns of the query matrix
 * This score is not normalized on length and cannot be used without a
 * normalization criterion as our DP based p-value approach.
 **************************************************************************/
void sandelin_scores(MOTIF_T* query_motif,
		     MOTIF_T* target_motif,
		     ARRAY_T* background,
		     MATRIX_T** score){

  /*   Recurse pairwise over the columns */
  int index_query;
  int index_target;
  for (index_query = 0; index_query < query_motif->length;
       index_query++) {
    for (index_target = 0; index_target < target_motif->length;
	 index_target++) {

      double SCORE_T = 0;

      int index; /* Index corresponds to the alphabets. */
      for (index = 0; index < query_motif->alph_size; index++) {

	double query_freq = get_matrix_cell(index_query,
					    index,
					    query_motif->freqs);
	double target_freq = get_matrix_cell(index_target,
					     index,
					     target_motif->freqs);
	/* Euclidean distance */
	double sq;
	sq = pow((query_freq - target_freq), 2);
	SCORE_T = SCORE_T + sq;
      }
      double SCORE_F;
      SCORE_F = 2 - SCORE_T;
      set_matrix_cell(index_query, index_target, SCORE_F, *score);
    } /* Target motif */
  } /* Query motif */
} /* Function sandelin*/

/**************************************************************************
 * Computes Kullback-Leiber for all columns of the query matrix against target
 * motif. Note we use
 * a modified Kullback-leibler for computing single column scores. We do not
 * divide the score by the width. This division is replaced by our DP based
 * p-value computation.
 **************************************************************************/
void kullback_scores
  (MOTIF_T* query_motif,
   MOTIF_T* target_motif,
   ARRAY_T* background,
   MATRIX_T** score)
{
  /*   Recurse pairwise over the columns */
  int index_query;
  int index_target;
  for (index_query = 0; index_query < query_motif->length;
       index_query++) {
    for (index_target = 0; index_target < target_motif->length;
	 index_target++) {
      double SCORE_F = 0;
      int index; /* Index corresponds to the alphabets. */
      for (index = 0; index < query_motif->alph_size; index++) {
	double query_freq = get_matrix_cell(index_query,
					    index,
					    query_motif->freqs);
	double target_freq = get_matrix_cell(index_target,
					     index,
					     target_motif->freqs);
	/* Average Kullback */
	double avg_kull;
	avg_kull = ((query_freq * log10(query_freq/target_freq))
		    + (target_freq * log10(target_freq/query_freq))) / 2;
	SCORE_F = SCORE_F + avg_kull;
      }
      SCORE_F = SCORE_F * -1;
      set_matrix_cell(index_query, index_target, SCORE_F, *score);
    } /* Target motif */
  } /* Query motif */
} /* Function Kullback*/


/**************************************************************************
 * Computes Pearson scores for all columns of the query matrix against target
 * motif. 
 **************************************************************************/
void pearson_scores
  (MOTIF_T* query_motif,
   MOTIF_T* target_motif,
   ARRAY_T* background,
   MATRIX_T** score)
{

  /*   Define the average values */
  double x_bar = 1 / query_motif->alph_size;
  double y_bar = 1 / target_motif->alph_size;

  /*   Recurse pairwise over the columns */
  int index_query;
  int index_target;
  for (index_query = 0; index_query < query_motif->length;
       index_query++) {
    for (index_target = 0; index_target < target_motif->length;
	 index_target++) {

      double SCORE_F = 0;
      double SCORE_NUMERATOR = 0;
      double SCORE_SQ_DENOMINATOR_1 = 0;
      double SCORE_SQ_DENOMINATOR_2 = 0;

      int index; /* Index corresponds to the alphabets. */
      for (index = 0; index < query_motif->alph_size; index++) {
	/* The numerator of the Pearson */
	SCORE_NUMERATOR
	  += (get_matrix_cell(index_query,
			      index,
			      query_motif->freqs) - x_bar)
	  * (get_matrix_cell(index_target,
			     index,
			     target_motif->freqs) - y_bar);

	/*  The denominoator of the Pearson */
	SCORE_SQ_DENOMINATOR_1
	  += pow((get_matrix_cell(index_query,
				  index,
				  query_motif->freqs) - x_bar),
		 2);
	SCORE_SQ_DENOMINATOR_2
	  += pow((get_matrix_cell(index_target,
				  index,
				  target_motif->freqs) - y_bar),
		 2);
      } /* Nr and Dr components summed over the entire alphabet. */

      /*       Compute the Denominator */
      double SCORE_SQ_DENOMINATOR = SCORE_SQ_DENOMINATOR_1
	* SCORE_SQ_DENOMINATOR_2;
      double SCORE_DENOMINATOR = sqrt(SCORE_SQ_DENOMINATOR);

      SCORE_F = SCORE_NUMERATOR / SCORE_DENOMINATOR;
      set_matrix_cell(index_query, index_target, SCORE_F, *score);
    } /* Target motif */
  } /* Query motif */
} /* Function pearson */

/**************************************************************************
 * Estimate one Dirichlet component.
 *
 * {\hat p}_i
 *     = \frac{n_i + \alpha_i}{\sum_{j \in \{A,C,G,T\}} (n_j + \alpha_j)}
 **************************************************************************/
ARRAY_T* one_dirichlet_component
  (int      hyperparameters[],
   ARRAY_T* counts)
{
  int alph_size = get_array_length(counts);
  ARRAY_T* return_value = allocate_array(alph_size);

  // Calculate the denominator.
  double denominator = 0.0;
  int i_alph;
  for (i_alph = 0; i_alph < alph_size; i_alph++) {
    denominator += get_array_item(i_alph, counts)
      + (int)(hyperparameters[i_alph]);
  }

  // Estimate the source distribution.
  for (i_alph = 0; i_alph < alph_size; i_alph++) {
    double numerator = get_array_item(i_alph, counts) 
      + hyperparameters[i_alph];
    set_array_item(i_alph, numerator / denominator, return_value);
  }
  return(return_value);
}

/**************************************************************************
 * Helper function to get un-logged gamma.
 **************************************************************************/
static double true_gamma (double x) { return(exp(gamma(x))); }

/**************************************************************************
 * Likelihood of the data, given a Dirichlet component.
 *
 * This is the final equation on p. 13 of Habib et al.
 *
 * It should probably be implemented in log space, rather than using
 * the true_gamma function above.
 **************************************************************************/
static double compute_likelihood
  (int      hypers[],
   ARRAY_T* counts)
{

  // Sum the counts and the hyperparameters.
  double hyper_sum = 0.0;
  double count_sum = 0.0;
  int i_alph;
  int alphabet_size = get_array_length(counts);
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    hyper_sum += (float)(hypers[i_alph]);
    count_sum += get_array_item(i_alph, counts);
  }

  // Compute the product term.
  double product = 1.0;
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    double hyper = (float)(hypers[i_alph]);
    double count = get_array_item(i_alph, counts);
    product *= true_gamma(hyper + count) / true_gamma(hyper);
  }

  return((true_gamma(hyper_sum) / true_gamma(hyper_sum + count_sum))
	 * product);
  
}

/**************************************************************************
 * Estimate the source distribution from counts, using a 1-component
 * or 5-component Dirichlet prior.
 *
 * The 5-component prior only works for DNA.
 **************************************************************************/
static ARRAY_T* estimate_source
 (int      num_dirichlets,
  ARRAY_T* counts)
{
  int dna_hypers[6][4] = {{1, 1, 1, 1}, // Single component
			  {5, 1, 1, 1}, // Five components ...
			  {1, 5, 1, 1},
			  {1, 1, 5, 1},
			  {1, 1, 1, 5},
			  {2, 2, 2, 2}};
  int protein_hypers[20] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
			    1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

  // Single component or 5-component mixture?
  int alph_size = get_array_length(counts);
  if (num_dirichlets == 1) {
    if (alph_size == 4) {
      return(one_dirichlet_component(dna_hypers[0], counts));
    } else if (alph_size == 20) {
      return(one_dirichlet_component(protein_hypers, counts));
    } else {
      die("Invalid alphabet size (%d) for BLiC score.\n", alph_size);
    }
  }

  // Die if trying to do proteins.
  if (alph_size != 4) {
    die("Sorry, BLiC5 is only implemented for DNA motifs.\n");
  }

  // Compute the denominator.
  double denominator = 0.0;
  int i_hyper;
  int n_hypers = 5;
  for (i_hyper = 1; i_hyper <= n_hypers; i_hyper++) {
    denominator 
      += (1.0 / (float) n_hypers)
      * compute_likelihood(dna_hypers[i_hyper], counts);
  }

  ARRAY_T* return_value = allocate_array(alph_size);
  for (i_hyper = 1; i_hyper <= n_hypers; i_hyper++) {
  
    /* Compute the posterior.
     *   \Pr(k|n) = \frac{\Pr(k) \Pr(n|k)}{\sum_j \Pr(j) \Pr(n|j)}
     */
    double posterior 
      = ((1.0 / (float) n_hypers) 
	 * compute_likelihood(dna_hypers[i_hyper], counts))
      / denominator;

    // Compute the Dirichlet component.
    ARRAY_T* one_component
      = one_dirichlet_component(dna_hypers[i_hyper], counts);

    // Multiply them together and add to the return value.
    scalar_mult(posterior, one_component);
    sum_array(one_component, return_value);

    free_array(one_component);
  }

  return(return_value);
}

/**************************************************************************
 * Computes the BLiC (Bayesian Likelihood 2-Component) score described
 * by Equation (2) in Habib et al., PLoS CB 4(2):e1000010.
 *
 * Thus far, BLiC is only implemented for DNA.  It would be relatively
 * straightforward to extend to amino acids.
 **************************************************************************/
static double one_blic_score
  (int      num_dirichlets,
   ARRAY_T* query_counts, 
   ARRAY_T* target_counts,
   ARRAY_T* background)
{

  // Make a vector that sums the two sets of counts.
  int alphabet_size = get_array_length(query_counts);
  assert(alphabet_size == 4);  // BLiC only works for DNA.
  ARRAY_T* common_counts = allocate_array(alphabet_size);
  copy_array(query_counts, common_counts);
  sum_array(target_counts, common_counts);

  // Estimate the source distributions.
  ARRAY_T* query_source = estimate_source(num_dirichlets, query_counts);
  ARRAY_T* target_source = estimate_source(num_dirichlets, target_counts);
  ARRAY_T* common_source = estimate_source(num_dirichlets, common_counts);

  // First term.
  double first_term = 0.0;
  int i_alph;
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    first_term += 2 * get_array_item(i_alph, common_counts) 
      * log(get_array_item(i_alph, common_source));
  }

  // Second term.
  double second_term = 0.0;
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    second_term += get_array_item(i_alph, query_counts)
      * log(get_array_item(i_alph, query_source));
  }

  // Third term.
  double third_term = 0.0;
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    third_term += get_array_item(i_alph, target_counts)
      * log(get_array_item(i_alph, target_source));
  }
  
  // Fourth term.
  double fourth_term = 0.0;
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    fourth_term += get_array_item(i_alph, common_counts)
      * log(get_array_item(i_alph, background));
  }

  return(first_term - (second_term + third_term + fourth_term));
  
}

/**************************************************************************
 * Computes BLiC scores for all pairs of positions in two motifs.
 **************************************************************************/
static void blic_scores
  (int      num_dirichlets, // Number of Dirichlets [IN]
   MOTIF_T* query_motif,  // Single query motif [IN].
   MOTIF_T* target_motif, // Single target motif [IN].
   ARRAY_T* background,   // Background distribution [IN].
   MATRIX_T** scores)     // BLiC score matrix, indexed by start. [OUT]
{
  // Traverse the query and target motifs
  int index_query;
  for (index_query = 0; index_query < query_motif->length; index_query++) {
    ARRAY_T* query_counts = get_motif_counts(index_query, query_motif);

    int index_target;
    for (index_target = 0; index_target < target_motif->length;
	 index_target++) {
      ARRAY_T* target_counts = get_motif_counts(index_target, target_motif);

      // Compute the BLiC score according to Equation (2).
      double my_score = one_blic_score(num_dirichlets,
				       query_counts, 
				       target_counts,
				       background);
      free_array(target_counts);
      
      // Store the computed score.
      set_matrix_cell(index_query, index_target, my_score, *scores);
    }
    free_array(query_counts);
  }
}

/**************************************************************************
 * Computes the llr score 
 **************************************************************************/
static double one_llr_score
  (int      num_dirichlets,
   ARRAY_T* query_counts, 
   ARRAY_T* target_counts,
   ARRAY_T* background)
{

  // Make a vector that sums the two sets of counts.
  int alphabet_size = get_array_length(query_counts);
  assert(alphabet_size == 4);  // BLiC only works for DNA.
  ARRAY_T* common_counts = allocate_array(alphabet_size);
  copy_array(query_counts, common_counts);
  sum_array(target_counts, common_counts);

  // Estimate the source distributions.
  ARRAY_T* query_source = estimate_source(num_dirichlets, query_counts);
  ARRAY_T* target_source = estimate_source(num_dirichlets, target_counts);
  ARRAY_T* common_source = estimate_source(num_dirichlets, common_counts);

  // First term.
  double first_term = 0.0;
  int i_alph;
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    first_term += get_array_item(i_alph, common_counts) 
      * log(get_array_item(i_alph, common_source));
  }

  // Second term.
  double second_term = 0.0;
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    second_term += get_array_item(i_alph, query_counts)
      * log(get_array_item(i_alph, query_source));
  }

  // Third term.
  double third_term = 0.0;
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    third_term += get_array_item(i_alph, target_counts)
      * log(get_array_item(i_alph, target_source));
  }

  return(first_term - (second_term + third_term));
  
}

/**************************************************************************
 * Computes log likelihood ratio scores for all pairs of positions in two motifs.
 **************************************************************************/
static void llr_scores
  (int      num_dirichlets, // Number of Dirichlets [IN]
   MOTIF_T* query_motif,  // Single query motif [IN].
   MOTIF_T* target_motif, // Single target motif [IN].
   ARRAY_T* background,   // Background distribution [IN].
   MATRIX_T** scores)     // BLiC score matrix, indexed by start. [OUT]
{
  // Traverse the query and target motifs
  int index_query;
  for (index_query = 0; index_query < query_motif->length; index_query++) {
    ARRAY_T* query_counts = get_motif_counts(index_query, query_motif);

    int index_target;
    for (index_target = 0; index_target < target_motif->length;
	 index_target++) {
      ARRAY_T* target_counts = get_motif_counts(index_target, target_motif);

      
      double my_score = one_llr_score(num_dirichlets,
				       query_counts, 
				       target_counts,
				       background);
      free_array(target_counts);
      
      // Store the computed score.
      set_matrix_cell(index_query, index_target, my_score, *scores);
    }
    free_array(query_counts);
  }
}


/**************************************************************************
 **************************************************************************/
// FIXME: currently prints HTML
void print_tomtom_xml_header(
  FILE *output_path,
  char *query_filename,
  char *target_filename,
  char *distance_measure,
  char *sig_type,
  float sig_thresh
) 
{
  fprintf(output_path,
    "<html>\n"
    "<head>\n"
    "  <title> Tomtom output</title>\n"
    "</head> <body>\n" 
    "<table border=1>\n"
    "  <tr> <th align=center colspan=2><h2>TOMTOM OUTPUT</h2></th></tr>\n"
    "  <tr> <td align=left colspan=2><b>Query File:</b> %s</td></tr>\n"
    "  <tr> <td align=left colspan=2><b>Target File:</b> %s</td></tr>\n"
    "  <tr> <td align=left colspan=2><b>Distance Measure:</b> %s</td></tr>\n"
    "  <tr> <td align=left colspan=2><b>All Motif Matches with <i>%s</i>-value at most:</b> %.2g.",
    query_filename,
    target_filename,
    distance_measure,
    sig_type,
    sig_thresh
  );
  if (strcmp(sig_type, "q") == 0) {
    fprintf(output_path,
      " The <i>q</i>-value is the estimated false discovery rate if the occurrence"
      " is accepted as significant."
      " See Storey JD, Tibshirani R. \"Statistical significance for genome-wide studies\"."
      " <i>Proc. Natl Acad. Sci. USA (2003) 100:9440--9445</i>"
    );
  }
  fprintf(output_path, "</td></tr>\n");
} // print_tomtom_xml_header

// FIXME: currently prints HTML
void print_tomtom_xml_line(
  FILE *output_path,		// path to output XML file
  TOMTOM_MATCH_T* my_match,	// stored results line
  char *output_dirname,		// name of output directory
  BOOLEAN_T ssc,   	        // Small sample correction
  MOTIF_T *motifs_query,	// query motifs
  MOTIF_T *motifs_target,	// target motifs
  char *query_url_format,	// format string for query URL
  char *target_url_format 	// format string for target URL
  ) 
{
  int query, target, offset, overlap;
  char query_id[MAX_MOTIF_ID_LENGTH+1], target_id[MAX_MOTIF_ID_LENGTH+1];
  double evalue, uncorrected_pvalue;
  char orientation;

  // create alignment logo
  char *logo_filename = create_logo( 
    &(motifs_target[my_match->target_index]), 
    &(motifs_query[my_match->query_index]),
    my_match->offset, 
    ssc,
    output_dirname
  );

  // create href for motifs by using the given format
  char *query_href = (char*) malloc(
    sizeof(char) * (strlen(query_url_format) + 2*strlen(my_match->query_id) + 1)
  );
  sprintf(query_href, query_url_format, my_match->query_id, my_match->query_id);
  char *target_href = (char*) malloc(
    sizeof(char) * (strlen(target_url_format) + 2*strlen(my_match->target_id) + 1)
  );
  sprintf(target_href, target_url_format, my_match->target_id, my_match->target_id);

  //FIXME: change to XML
  fprintf(output_path,
    "  <tr>\n"
    "    <td>\n"
    "      <table>\n"
    "        <tr><th align=left nowrap>Target Motif:</th><td align=right>%s</td></tr>\n"	// target id
    "        <tr><th align=left nowrap>Target Description:</th><td align=right>%s</td></tr>\n"	// target description
    "        <tr><th align=left nowrap>Query Motif:</th><td align=right>%s</td></tr>\n"		// query id
    "        <tr><th align=left nowrap>Query Description:</th><td align=right>%s</td></tr>\n"	// query description
    "        <tr><th align=left nowrap><i>p</i>-value:</th><td align=right><b>%.2g<b></td></tr>\n"	// p-value
    "        <tr><th align=left nowrap><i>E</i>-value:</th><td align=right><b>%.2g<b></td></tr>\n"	// E-value
    "        <tr><th align=left nowrap><i>q</i>-value:</th><td align=right><b>%.2g<b></td></tr>\n"	// q-value
    "        <tr><th align=left nowrap>Overlap:</th><td align=right>%d</td></tr>\n"		// overlap
    "        <tr><th align=left nowrap>Query Offset:</th><td align=right>%d</td></tr>\n"	// offset
    "        <tr><th align=left nowrap>Orientation:</th><td align=right>%c</td></tr>\n"		// orientation
    "        <tr><th align=left nowrap>Figures:</th>"						// links
    "<td align=right><a href='%s.eps'>[EPS]</a><a href='%s.png'>[PNG]</a></td></tr>\n"
    "      </table>\n"
    "    </td>\n"
    "    <td><img src='%s.png'></td>\n"							// logo
    "  </tr>\n",
    target_href, 
    motifs_target[my_match->target_index].id2,	// target description
    query_href, 
    motifs_query[my_match->query_index].id2,	// query description
    my_match->pvalue,
    my_match->evalue,
    my_match->qvalue,
    my_match->overlap,
    my_match->offset,
    my_match->orientation,
    logo_filename, logo_filename,		// links
    logo_filename 				// image
  );

  myfree(logo_filename);
  myfree(query_href);
  myfree(target_href);
} // print_tomtom_xml_line

// FIXME: currently prints HTML
void print_tomtom_xml_tailer(
  FILE *output_path,
  int argc,
  char **argv,
  int total_matches,
  char *sig_type,
  float sig_thresh
) 
{
  //FIXME: change to XML
  int i;

  fprintf(output_path, 
    "  <tr>\n"
    "    <td> <b>Total matches with <i>%s</i>-value &le; %g:</b> %d\n"
    "    </td>\n"
    "  </tr>\n", 
    sig_type,
    sig_thresh,
    total_matches);
  fprintf(output_path, "<tr><th>Command line:</th><td>\n");
  for (i=0; i<argc; i++) {
    fprintf(output_path, "%s\n", argv[i]);
  }
  fprintf(output_path, "</td></tr>\n");
  fprintf(output_path,
    "</table>\n\n"
    "</body>\n"
    "</html>\n"
  );
} // print_tomtom_xml_tailer

/*****************************************************************************
 * MAIN
 *****************************************************************************/
VERBOSE_T verbosity = NORMAL_VERBOSE;
#ifdef MAIN

int main(int argc, char *argv[])
{
  char *default_output_dirname = "tomtom_out"; // where to write HTML
  char *text_filename = "tomtom.txt";
  FILE *text_output;			// where to write legacy output
  char *xml_filename = "tomtom.html";	// FIXME: change to .xml and output XML
  FILE *xml_output = NULL;
  //char *html_filename = "tomtom.html";// FIXME: convert XML to HTML
  //FILE *html_output = NULL;

  /**********************************************
   * Initialize command line parameters
   **********************************************/
  char* query_filename = NULL;
  char* target_filename = NULL;
  char *output_dirname = default_output_dirname;
  BOOLEAN_T clobber = FALSE;		// don't clobber existing files
  BOOLEAN_T text_only = FALSE; 		// default is HTML
  char* distance_measure = "ed";
  char* sig_type = "q";
  float sig_thresh = 0.5;
  char *target_url_format = "%s"; 	// no URL format
  char *query_url_format = "%s";	// no URL format
  float query_pseudocount = DEFAULT_QUERY_PSEUDO; 
  float target_pseudocount = DEFAULT_TARGET_PSEUDO; 
  char* columnwise_scores_filename = NULL;
  BOOLEAN_T ssc = TRUE;			// small sample correction in logos

  /**********************************************
   * Command line parsing.
   **********************************************/
  int option_count = 18;     // Must be updated if list below is
  cmdoption const motif_scan_options[] = {
    {"query", REQUIRED_VALUE},
    {"target", REQUIRED_VALUE},
    {"o", OPTIONAL_VALUE},
    {"oc", OPTIONAL_VALUE},
    {"text", OPTIONAL_VALUE},
    {"thresh", OPTIONAL_VALUE},
    {"q-thresh", OPTIONAL_VALUE},
    {"evalue", OPTIONAL_VALUE},
    {"dist", OPTIONAL_VALUE},
    {"internal", OPTIONAL_VALUE},
    {"min-overlap", OPTIONAL_VALUE},
    {"query-pseudo", OPTIONAL_VALUE},
    {"target-pseudo", OPTIONAL_VALUE},
    {"query-url", OPTIONAL_VALUE},
    {"target-url", OPTIONAL_VALUE},
    {"column-scores", OPTIONAL_VALUE},
    {"no-ssc", OPTIONAL_VALUE},
    {"verbosity", REQUIRED_VALUE},
    {"column-scores", OPTIONAL_VALUE}
  };

  int option_index = 0;

  /**********************************************
   * Define the usage message.
   **********************************************/
  char * usage = "\n   USAGE: tomtom [options] -query <query file> -target <target file>\n"
  "\n"
  "   Options:\n"
  "     -o <output dir>\t\tname of directory for output files;\n\t\t\t\twill not replace existing directory\n"
  "     -oc <output dir>\t\tname of directory for output files;\n\t\t\t\twill replace existing directory\n"
  "     -thresh <float>\t\tsignificance threshold; default: 0.5\n"
  "     -evalue\t\tuse E-value threshold; default: q-value\n"
  "     -dist allr|ed|kullback|pearson|sandelin|blic1|blic5\n\t\t\t\tdistance metric for scoring alignments;\n\t\t\t\tdefault: ed (Euclidean distance)\n"
  "     -internal\t\t\tonly allow internal alignments;\n\t\t\t\tdefault: allow overhangs\n"
  "     -min-overlap <int>\t\tminimum overlap between query and target;\n\t\t\t\tdefault: 1\n"
  "     -query-pseudo <float>\tdefault: %.1f\n"
  "     -target-pseudo <float>\tdefault: %.1f\n"
  "     -text\t\t\toutput in text format (default is HTML)\n"
  "     -query-url <string>\tformat string for links to query motifs\n\t\t\t\tdefault: none\n"
  "     -target-url <string>\tformat string for links to target motifs\n\t\t\t\tdefault: none\n"
  "     -no-ssc\t\t\tdon't apply small-sample correction to logos\n\t\t\t\tdefault: use small-sample correction\n"
  "     -verbosity [1|2|3|4]\tdefault: %d\n"
  // The "--column-scores <file>" option is undocumented -- developers only.
  "\n";

  /**********************************************
   * Parse the command line.
   **********************************************/
  if (simple_setopt(argc, argv, option_count, motif_scan_options) != NO_ERROR)
    {
      die("Error: option name too long.\n");
    }
  BOOLEAN_T internal = FALSE;
  int min_overlap = 1;
  while (TRUE) {
    int c = 0;
    char* option_name = NULL;
    char* option_value = NULL;
    const char * message = NULL;

    /*  Read the next option, and break if we're done. */
    c = simple_getopt(&option_name, &option_value, &option_index);
    if (c == 0) {
      break;
    } else if (c < 0) {
      (void) simple_getopterror(&message);
      die("Error processing command line options (%s)\n", message);
    }

    if (strcmp(option_name, "thresh") == 0) {
      if (! option_value) {
        fprintf(stderr, "error: missing value after -thresh.\n");
        exit(EXIT_FAILURE);
      }
      sig_thresh = atof(option_value);
    } else if (strcmp(option_name, "q-thresh") == 0) {	// kept for backward compatibility
      if (! option_value) {
        fprintf(stderr, "error: missing value after -q-thresh.\n");
        exit(EXIT_FAILURE);
      }
      sig_thresh = atof(option_value);
    } else if (strcmp(option_name, "evalue") == 0) {
      sig_type = "E";
    } else if (strcmp(option_name, "query") == 0) {
      query_filename = option_value;
    } else if (strcmp(option_name, "target") == 0) {
      target_filename = option_value;
    } else if (strcmp(option_name, "o") == 0) {
      output_dirname = option_value;
    } else if (strcmp(option_name, "oc") == 0) {
      output_dirname = option_value;
      clobber = TRUE;		// OK to clobber
    } else if (strcmp(option_name, "text") == 0) {
      text_only = TRUE;
    } else if (strcmp(option_name, "no-ssc") == 0) {
      ssc = FALSE;
    } else if (strcmp(option_name, "dist") == 0) {
      distance_measure = option_value;
    } else if (strcmp(option_name, "internal") == 0) {
      internal = TRUE;
    } else if (strcmp(option_name, "min-overlap") == 0) {
      double temp_var;
      if (! option_value) {
        fprintf(stderr, "error: missing value after -min-overlap.\n");
        exit(EXIT_FAILURE);
      }
      temp_var = atof(option_value);
      min_overlap = rint(temp_var);
    } else if (strcmp(option_name, "query-pseudo") == 0) {
      if (! option_value) {
        fprintf(stderr, "error: missing value after -query-pseudo.\n");
        exit(EXIT_FAILURE);
      }
      query_pseudocount = atof(option_value);
    } else if (strcmp(option_name, "target-pseudo") == 0) {
      if (! option_value) {
        fprintf(stderr, "error: missing value after -target-pseudo.\n");
        exit(EXIT_FAILURE);
      }
      target_pseudocount = atof(option_value);
    } else if (strcmp(option_name, "query-url") == 0) {
      query_url_format = option_value;
    } else if (strcmp(option_name, "target-url") == 0) {
      target_url_format = option_value;
    } else if (strcmp(option_name, "verbosity") == 0) {
      verbosity = atoi(option_value);
    } else if (strcmp(option_name, "column-scores") == 0) {
      columnwise_scores_filename = option_value;
    }
  }
  if ((query_filename == NULL)||(target_filename == NULL)) {
    fprintf(stderr, usage, DEFAULT_QUERY_PSEUDO, DEFAULT_QUERY_PSEUDO, NORMAL_VERBOSE);
    exit(EXIT_FAILURE);
  }
  if (text_only == TRUE) {
    // Legacy: plain text output to standard out.
    text_output = stdout;
  } else {
    // allow clobbering of the default output directory
    if (output_dirname == default_output_dirname) {
      clobber = TRUE;
    }
    if (create_output_directory(output_dirname, clobber, (verbosity >= NORMAL_VERBOSE))) {
      // Failed to create output directory.       
      die("Unable to create output directory %s.\n", output_dirname);
    }
    // Create the name of the output files (text, XML, HTML) and open them
    char *path;
    path = make_path_to_file(output_dirname, text_filename);
    text_output = fopen(path, "w"); //FIXME CEG check for errors
    myfree(path);
    path = make_path_to_file(output_dirname, xml_filename);
    xml_output = fopen(path, "w"); //FIXME CEG check for errors
    myfree(path);
    //path = make_path_to_file(output_dirname, html_filename);
    //html_output = fopen(path, "w"); //FIXME CEG check for errors
    //myfree(path);
  }

  /**********************************************
   * Read the motifs.
   **********************************************/

  // General variables for reading the meme file
  int num_motifs_query = 0;
  BOOLEAN_T has_reverse_strand_query = FALSE;
  MOTIF_T motifs_query[MAX_QUERYS];
  STRING_LIST_T* motif_occurrences = NULL;
  ARRAY_T* bg_freqs_query;

  read_meme_file(
		 query_filename,
		 "motif-file", // background freq. filename
		 query_pseudocount,
     REQUIRE_PSPM,
		 &num_motifs_query,
		 motifs_query,
		 &motif_occurrences,
		 &has_reverse_strand_query,
		 &bg_freqs_query
		 );

  int num_motifs_target = 0;
  MOTIF_T motifs_target[MAX_TARGETS];
  ARRAY_T* bg_freqs_target = NULL;
  BOOLEAN_T has_reverse_strand_target = FALSE;
  read_meme_file(
		 target_filename,
		 "motif-file", // background freq. filename
		 target_pseudocount,
     REQUIRE_PSPM,
		 &num_motifs_target,
		 motifs_target,
		 &motif_occurrences,
		 &has_reverse_strand_target,
		 &bg_freqs_target
		 );
  add_reverse_complements(&num_motifs_target, motifs_target);


  // Check for small target databases (multiply by 2 for reverse complements).
  if (num_motifs_target < (MIN_TARGET_DATABASE_SIZE * 2)) {
    fprintf(stderr,
	    "\nWarning: Target database size too small (%d). Provide at least %d motifs for accurate p-value computation\n\n",
	    num_motifs_target/2, MIN_TARGET_DATABASE_SIZE);
  }
  int motif_index_query;
  int motif_index_target;
  int total_target_length = 0;

  for(motif_index_target = 0; motif_index_target < num_motifs_target;
      motif_index_target++){
    MOTIF_T* motif_target = &(motifs_target[motif_index_target]);
    total_target_length += motif_target->length;
    /*   Remove the ambiguous letter columns */
    while(get_num_cols(motif_target->freqs) > motif_target->alph_size)
      {
	remove_matrix_col(motif_target->alph_size, motif_target->freqs);
      }
  }

  /*   Remove the ambiguous letter columns for query */
  for(motif_index_query = 0; motif_index_query < num_motifs_query;
      motif_index_query++){
    MOTIF_T* motif_query = &(motifs_query[motif_index_query]);
    while(get_num_cols(motif_query->freqs) > motif_query->alph_size)
      {
	remove_matrix_col(motif_query->alph_size, motif_query->freqs);
      }
  }

  /* Uniform background distribution array for the infinite alphabet*/
  ARRAY_T* bg_freqs_mod = allocate_array(total_target_length);
  init_array((ATYPE)(1.0 / total_target_length), bg_freqs_mod);

  /*  Recurse over the set of query motifs and target motifs */
  int total_matches = 0;

  // Print output
  fprintf(text_output, "#Query ID\tTarget ID\tOptimal offset");
  fprintf(text_output, "\tp-value\tE-value\tq-value\tOverlap\t");
  fprintf(text_output, "Query consensus\tTarget consensus\tOrientation\n");
  if (! text_only) { 
    print_tomtom_xml_header(
			    xml_output, 
			    query_filename, 
			    target_filename,
			    distance_measure,
			    sig_type,
			    sig_thresh
			    ); 
  }

  // Decide the number of runs.
  int query_iterations = num_motifs_query;
  int num_motifs_database = num_motifs_target;

  /* Recurse over the set of queries. */
  for(motif_index_query = 0; motif_index_query < query_iterations;
      motif_index_query++){
    OBJECT_LIST_T* output_list = new_object_list(NULL, NULL, NULL, free);

    fprintf(stderr, "Processing query %d out of %d \n",
	    motif_index_query + 1, query_iterations);
    int curr_query_index = motif_index_query;
    int target_database_length = total_target_length;
    int num_motifs_database = num_motifs_target;

    MOTIF_T* motif_query = &(motifs_query[curr_query_index]);

    /*  Initialize an array to store the columnwise scores. */
    MATRIX_T* all_columnwise_scores;
    all_columnwise_scores = allocate_matrix(
					    motif_query->length,
					    target_database_length
					    );

    int index_col_all_columnwise_scores = 0;
    int *target_lengths; /* Array to store lengths of target motifs */
    target_lengths = malloc(num_motifs_database * sizeof(int));
    int target_lengths_index = 0;
    STRING_LIST_T* motif_id_target = new_string_list();
    STRING_LIST_T* motif_cons_target = new_string_list();

    /* Recurse over the set of targets. */
    for(motif_index_target = 0; motif_index_target < num_motifs_target;
	motif_index_target++){
      MOTIF_T* motif_target = &(motifs_target[motif_index_target]);  ;

	add_string(get_motif_id(motif_target), motif_id_target);

	/* Add target Consensus Sequence to string list*/
	char* target_cons;
	get_cons(motif_target->freqs, &target_cons);
	add_string(target_cons, motif_cons_target);

	/* Compute columnwise scores */
	MATRIX_T* pairwise_column_scores;
	pairwise_column_scores = allocate_matrix(
						 motif_query->length,
						 motif_target->length
						 );

	/* Compute score matrix based on the distance measure */
	if (strcmp(distance_measure, "pearson") == 0) {
	  pearson_scores(
			 motif_query,
			 motif_target,
			 bg_freqs_target,
			 &pairwise_column_scores
			 );
	} else if (strcmp(distance_measure, "allr") == 0) {
	  allr_scores(
		      motif_query,
		      motif_target,
		      bg_freqs_target,
		      &pairwise_column_scores
		      );
	} else if (strcmp(distance_measure, "ed") == 0) {
	  ed_scores(
		    motif_query,
		    motif_target,
		    bg_freqs_target,
		    &pairwise_column_scores
		    );
	} else if (strcmp(distance_measure, "kullback") == 0) {
	  kullback_scores(
			  motif_query,
			  motif_target,
			  bg_freqs_target,
			  &pairwise_column_scores
			  );
	} else if (strcmp(distance_measure, "sandelin") == 0) {
	  sandelin_scores(
			  motif_query,
			  motif_target,
			  bg_freqs_target,
			  &pairwise_column_scores
			  );
	} else if (strcmp(distance_measure, "blic1") == 0) {
	  blic_scores(
		      1,
		      motif_query,
		      motif_target,
		      bg_freqs_target,
		      &pairwise_column_scores
		      );
	} else if (strcmp(distance_measure, "blic5") == 0) {
	  blic_scores(
		      5,
		      motif_query,
		      motif_target,
		      bg_freqs_target,
		      &pairwise_column_scores
		      );
	} else if (strcmp(distance_measure,"llr1") == 0) {
	  llr_scores(
		     1,
		     motif_query,
		     motif_target,
		     bg_freqs_target,
		     &pairwise_column_scores
		     );
	} else if (strcmp(distance_measure,"llr5") == 0) {
	  llr_scores(
		     5,
		     motif_query,
		     motif_target,
		     bg_freqs_target,
		     &pairwise_column_scores
		     );
	} else {
	  fprintf(stderr,
		  "Invalid Distance measure [allr|ed|kullback|pearson|sandelin|blic1|blic5|llr]\n"
		  );
	}

	/* Update the background frequency array */
        // tlb added documentation; 
        // the "frequency array" is [w_q, \sum w_t] and contains
	// scores of query motif columns (the rows) vs all target columns
        // with each new target's scores appended as a new set of columns
        // Each row of the frequency array contains the scores of a single
        // query column vs. all possible target columns.
	int index_rows;
	for(index_rows = 0;index_rows < pairwise_column_scores->num_rows;
	    index_rows++){
	  int index_columns;
	  for(index_columns = 0;
	      index_columns < pairwise_column_scores->num_cols;
	      index_columns++){
	    set_matrix_cell(
			    index_rows,
			    index_columns + index_col_all_columnwise_scores,
			    get_matrix_cell(index_rows,
					    index_columns,
					    pairwise_column_scores),
			    all_columnwise_scores
			    );
	  }
	}

	index_col_all_columnwise_scores +=
	  get_num_cols(pairwise_column_scores);
	target_lengths[target_lengths_index] = motif_target->length;
	target_lengths_index++;
	free_matrix(pairwise_column_scores);

    }/* Parsing all targets for the query done. */

    // If requested, store the columnwise scores in an external file.
    if (columnwise_scores_filename != NULL) {
      fprintf(stderr, 
	      "Storing %d by %d matrix of column scores in %s.\n", 
	      get_num_rows(all_columnwise_scores),
	      get_num_cols(all_columnwise_scores),
	      columnwise_scores_filename);
      FILE* columnwise_scores_file;
      if (!open_file(columnwise_scores_filename, "w", FALSE, "column scores",
		     "column scores", &columnwise_scores_file)) {
	exit(1);
      }
      print_matrix(all_columnwise_scores, 10, 7, FALSE,
		   columnwise_scores_file);
      fclose(columnwise_scores_file);
    }

    /* Scale the pssm */
    PSSM_T* pssm = build_matrix_pssm(
      all_columnwise_scores,	// matrix
      NULL,		// don't compute pv lookup table
      BINS 		// range
    );
    // Bad--exposing internals of PSSM object.
    // Get the scaled pssm matrix and throw away the pssm object.
    SWAP (MATRIX_T *, pssm->matrix, all_columnwise_scores);
    free_pssm(pssm);

    /* Compute the pvalue distributions */
    int num_pv_lookup_array;		// w * (w+1) / 2
    num_pv_lookup_array = (pow(all_columnwise_scores->num_rows, 2)
			   + (all_columnwise_scores->num_rows))/2;

    //int address_index = 0;
    MATRIX_T* reference_matrix;
    reference_matrix = allocate_matrix(
				       all_columnwise_scores->num_rows,
				       all_columnwise_scores->num_rows
				       );
    init_matrix(-1, reference_matrix);

    int num_array_lookups = (BINS * motif_query->length) + 1;
    MATRIX_T *pv_lookup_matrix = allocate_matrix(
				       num_pv_lookup_array,
				       num_array_lookups
				       );
    init_matrix(0, pv_lookup_matrix);

    get_pv_lookup_new(
	  all_columnwise_scores,	// scores of each query column
	  target_database_length,	// number of scores 
	  BINS,				// maximum score in a column
	  bg_freqs_mod,			// Uniform
	  reference_matrix,
	  pv_lookup_matrix
	  );

    /* Compute the motif scores */
    int colid = 0; /* Stores the lengths  of already read columns */
    for(motif_index_target = 0;
	motif_index_target < num_motifs_database; motif_index_target++) {

      MATRIX_T* pairwise_subset;
      pairwise_subset = allocate_matrix(
					motif_query->length,
					target_lengths[motif_index_target]
					);

      /* Get a query-target pair from the overall score matrix. */
      int col_index;
      for(col_index = 0; col_index < target_lengths[motif_index_target];
	  col_index++){
	ARRAY_T* temp_array = NULL;
	temp_array = get_matrix_column(col_index + colid,
				       all_columnwise_scores);
	set_matrix_column(temp_array, col_index, pairwise_subset);
	free_array(temp_array);
      }

      /* Compute the score and offset wrt the minimum pvalue CORRECTED. */
      int optimal_offset = -1;
      double optimal_pvalue = 1;
      int overlapping_bases = -1;

      compare_motifs(
		     pairwise_subset,
		     &optimal_offset,
		     reference_matrix,		// index into pv_lookup_matrix
		     pv_lookup_matrix,		// based on start, stop columns
		     &optimal_pvalue,		// of query motif
		     &overlapping_bases,
		     internal,
		     min_overlap
		     //address_index 
		     );

      /* Override min-overlap if query is smaller than min-overlap */
      int mo;
      mo = min_overlap;
      if (pairwise_subset->num_rows < mo) {
        mo = pairwise_subset->num_rows;
      }
      if (pairwise_subset->num_cols < mo) {
        mo = pairwise_subset->num_cols;
      }
      /* Correct the motif pvalue (minimization across several pvalues).
	 Configurations mulitiplied by 2, because later the minimum of
	 the two strands will be taken LATER. Both stands have the same
	 number of total_configurations, thus muliplication by 2.*/
      int total_configurations = pairwise_subset->num_rows
        + pairwise_subset->num_cols - (2 * mo) + 1;

      /* Correct the total number of configurations for internal alignments */
      if (internal) {
        total_configurations = abs(pairwise_subset->num_rows
				   - pairwise_subset->num_cols) + 1;
      }
      double motif_pvalue = EV(optimal_pvalue, total_configurations * 2);

      // Compute E-value.
      double motif_evalue = motif_pvalue * num_motifs_target / 2.0;

      /* Consensus Sequence */
      char* query_cons;
      get_cons(motif_query->freqs, &query_cons);

      /* Offset_by_length */
      int offset_by_length = optimal_offset - motif_query->length;

      // Store this match.
      TOMTOM_MATCH_T* my_match
	= new_tomtom_match(
			   motif_query->id, 
			   get_nth_string(motif_index_target,
					  motif_id_target),
			   query_cons,
			   get_nth_string(motif_index_target,
					  motif_cons_target),
			   offset_by_length,
			   overlapping_bases,
			   motif_pvalue,
			   motif_evalue,
			   curr_query_index,
			   motif_index_target,
			   offset_by_length
			   );

      // Store this match.  Don't bother with the key.
      store_object((void*)my_match, NULL, motif_pvalue, output_list);

      // Ensure that create_logo works
      assert(motif_query == &(motifs_query[curr_query_index]));

      colid += target_lengths[motif_index_target];
      free_matrix(pairwise_subset);
    } // Scores computed for all targets.

    // Free all the matrices/arrays related to the current query
    free_matrix(reference_matrix);
    free_matrix(pv_lookup_matrix);
    free_matrix(all_columnwise_scores);
    free_string_list(motif_id_target);
    free_string_list(motif_cons_target);
    free(target_lengths);

    // Sort the stringlist.
    sort_objects(output_list);

    // Compute q-values.
    convert_tomtom_p_to_q(output_list);

    // Do some re-formatting of the output strings.
    OBJECT_LIST_T* output_final = new_object_list(NULL, NULL, NULL, free);
    select_strand(output_list, output_final);

    // Traverse the list of matches.
    TOMTOM_MATCH_T* my_match
      = (TOMTOM_MATCH_T*)retrieve_next_object(output_final);
    int index = 0;
    while (my_match != NULL) {

      // Check significance threshold.
      if ( ((strcmp(sig_type, "q") == 0) && my_match->qvalue <= sig_thresh) ||
        ((strcmp(sig_type, "E") == 0) && my_match->evalue <= sig_thresh) ) {
         total_matches++;

	// Print text output.
	fprintf(text_output,
		"%s\t%s\t%d\t%g\t%g\t%g\t%d\t%s\t%s\t%c\n",
		my_match->query_id, 
		my_match->target_id,
		my_match->offset,
		my_match->pvalue,
		my_match->evalue,
		my_match->qvalue,
		my_match->overlap, 
		my_match->query_consensus,
		my_match->target_consensus, 
		my_match->orientation
		);

	// Print XML output.
	if (! text_only) {
	  print_tomtom_xml_line(
				xml_output, 
				my_match,
				output_dirname,
  				ssc,
				motifs_query,
				motifs_target,
				query_url_format,
				target_url_format
				);
	}
      }

      // Get the next match.
      my_match
	= (TOMTOM_MATCH_T*)retrieve_next_object(output_final);
      index++;
    }

    free_object_list(output_final);
    free_object_list(output_list);

  }/*  The loop across the queries. */

  // print XML tailer
  if (! text_only) {
    print_tomtom_xml_tailer(xml_output, argc, argv, total_matches, sig_type, sig_thresh);
  }

  free_motif(motifs_query);
  free_motif(motifs_target);
  free_array(bg_freqs_query);
  free_array(bg_freqs_target);
  free_array(bg_freqs_mod);
  return (0);
}/* Main tomtom*/
#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
