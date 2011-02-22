/********************************************************************
 * FILE: mhmm.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 06-23-96
 * PROJECT: MHMM
 * COPYRIGHT: 1996, WNG
 * DESCRIPTION: Given MEME output, write a motif-based HMM.
 ********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "utils.h"           // Generic utilities. 
#include "matrix.h"          // Matrix manipulation routines. 
#include "string-list.h"     // List of strings (for motif IDs).
#include "metameme.h"        // Global Meta-MEME functions. 
#include "motif.h"           // Data structure for holding motifs. 
#include "order.h"           // Linear motif order and spacing. 
#include "meme-io.h"         // Read motifs from MEME output. 
#include "mhmm-state.h"      // Data structure for HMMs. 
#include "build-hmm.h"       // Construct an HMM in memory. 
#include "write-mhmm.h"      // Functions for producing MHMM output. 

#define SLOP 1E-5

char* program = "mhmm";
extern char* END_TRANSITION; // Marker for end transition.

/***********************************************************************
 * Convert transition counts to transition probabilities, and compute
 * average spacer lengths.
 *
 * Each matrix is indexed 0 ... n+1, where n is the number of motifs.
 * The entry at [i,j] corresponds to the transition from motif i to
 * motif j.  Hence, after normalization, each row in the transition
 * matrix should sum to 1.
 ***********************************************************************/
static void normalize_spacer_counts(
  double    trans_pseudo,
  double    spacer_pseudo,    // Pseudocount for self-loop.
  BOOLEAN_T keep_unused,
  MATRIX_T* transp_freq,
  MATRIX_T* spacer_ave
) {
  int i_row;
  int i_col;
  int num_rows;
  double total_spacer;
  double num_transitions;
  double ave_spacer;
  
  /* Divide the spacer lengths by the number of occurrences. */
  num_rows = get_num_rows(transp_freq);
  for (i_row = 0; i_row < num_rows; i_row++) {
    for (i_col = 0; i_col < num_rows; i_col++) {
      total_spacer = get_matrix_cell(i_row, i_col, spacer_ave) + spacer_pseudo;
      num_transitions = get_matrix_cell(i_row, i_col, transp_freq);
      if (spacer_pseudo > 0) num_transitions++;
      if (num_transitions != 0.0) {
        ave_spacer = total_spacer / num_transitions;
        set_matrix_cell(i_row, i_col, ave_spacer, spacer_ave);
      }
    }
  }

  // Add pseudocounts.
  for (i_row = 0; i_row < num_rows; i_row++) {
    for (i_col = 0; i_col < num_rows; i_col++) {

      // Force some transitions to zero.
      if (// No transitions to the start state.
        (i_col == 0) || 
        // No transitions from the end state.
        (i_row == num_rows - 1) ||
        // No transition from start to end.
        ((i_row == 0) && (i_col == num_rows - 1))) {
        set_matrix_cell(i_row, i_col, 0.0, transp_freq);
      }
      else {
        // Only increment the used transitions.
        if ((keep_unused) || 
            (get_matrix_cell(i_row, i_col, transp_freq) > 0.0)) {
          incr_matrix_cell(i_row, i_col, trans_pseudo, transp_freq);
        }
      }
    }
  }

  // Normalize rows.
  for (i_row = 0; i_row < num_rows - 1; i_row++) {
    if (array_total(get_matrix_row(i_row, transp_freq)) > 0.0) {
      normalize(SLOP, get_matrix_row(i_row, transp_freq));
    }
  }
}

/***********************************************************************
 * Add a "+" or "-" to the beginning of a given motif ID, if it is not
 * already there.
 *
 * This function is only necessary to allow backward compatibility
 * with MEME files that do not include strand information in the motif
 * occurrence section.
 ***********************************************************************/
static void add_strand(char* motif_id) {

  static BOOLEAN_T warned = FALSE;
  static char new_motif_id[MAX_MOTIF_ID_LENGTH+1];

  // Do we already have a strand indicator?
  if ((motif_id[0] == '+') || (motif_id[0] == '-')) {
    return;
  }

  // Tell the user about old MEME file.
  if (!warned) {
    fprintf(stderr, "Warning: The MEME file you provided is old and ");
    fprintf(stderr, "therefore lacks strand\ninformation in the motif ");
    fprintf(stderr, "occurrence section.  Meta-MEME will treat all\n");
    fprintf(stderr, "motif occurrences as if they appear on the forward ");
    fprintf(stderr, "strand.\n");
    warned = TRUE;
  }

  // Read the motif id.
  if (sscanf(motif_id, "%s", new_motif_id) != 1) {
    die("Can't read motif ID (%s).", motif_id);
  }

  // Convert it back to a string.
  sprintf(motif_id, "+%s", new_motif_id);
}

/***********************************************************************
 * Find the location of a given motif in the transition or spacer
 * matrix.
 *
 * Each matrix is indexed 0 ... n+1, where n is the number of motifs.
 * The entry at [i,j] corresponds to the transition from motif i to
 * motif j.
 *
 * For single-stranded motifs, the motif_id string is just a number
 * between 1 and n, inclusive.  However, if the motif_id string
 * contains a sign (+ or -), then the positive and negative strand
 * versions of the motif are in the matrix.  These are indexed 0
 * ... 2n+1, where 1 ... n correpond to the forward strand (+)
 * versions, and n+1 .. 2n correspond to the reverse strand (-)
 * versions.
 *
 * The parameter 'num_motifs' is the total number of motifs: n for the
 * single-stranded case and 2n for the double-stranded case.
 ***********************************************************************/
static int find_matrix_location(
  MOTIF_T* motifs,
  char* motif_id,
  int   num_motifs
) {

  int motif_index = 0;

  // Make sure we're not looking for an end transition.
  assert(strcmp(motif_id, END_TRANSITION) != 0);

  if (strncmp(motif_id, "0", MAX_MOTIF_ID_LENGTH)==0) {
    // A motif_id of "0" could only come from the BEGIN state
    motif_index = 0;
  } else {
    // Search for this motif_id in the array of retained motifs
    int i;
    for (i = 0; i < num_motifs; i++) {
      if(strncmp(motif_id, motifs[i].id, MAX_MOTIF_ID_LENGTH)==0) {
        motif_index = i + 1;
        break;
      }
    }
    // We should always have found the motif in the array
    assert(i < num_motifs);
  }
  return(motif_index);
}

/***********************************************************************
 * Record one motif occurrence.
 ***********************************************************************/
static void record_occurrence(
  char*     sequence_id,
  char*     motif_id,
  double    p_threshold,
  double    motif_p,
  char*     prev_motif,
  int*      prev_position,
  int       motif_position,
  MATRIX_T* transp_freq,
  MATRIX_T* spacer_ave,
  ORDER_T*  new_order,
  int       num_motifs,
  MOTIF_T*  motifs
) {

  int prev_motif_location;
  int motif_location;

  /* Always include transitions to the end state. */
  if (strcmp(motif_id, END_TRANSITION) == 0) {
  }

  /* Check to see if this motif was not given. */
  else if (!have_motif(motif_id, num_motifs, motifs)) {

    if (verbosity > NORMAL_VERBOSE) {
      fprintf(stderr, "Skipping motif %s in sequence %s.\n", motif_id,
        sequence_id);
    }
    return;
  }

  /* Check to see if the threshold was exceeded. */
  else if ((p_threshold > 0.0) && (motif_p > p_threshold)) {

    if (verbosity > NORMAL_VERBOSE) {
      fprintf(stderr, "Skipping motif %s in sequence %s (%g > %g).\n",
        motif_id, sequence_id, motif_p, p_threshold);
    }
    return;
  }

  if (verbosity > NORMAL_VERBOSE) {
    if (strcmp(motif_id, END_TRANSITION) == 0) {
      fprintf(stderr, "Adding transition to end in sequence %s.\n",
        sequence_id);
    } else {
      fprintf(stderr, "Adding motif %s in sequence %s.\n", motif_id,
        sequence_id);
    }
  }

  // If we're at the end, store in the last column of the matrix.
  if (strcmp(motif_id, END_TRANSITION) == 0) {
    motif_location = num_motifs + 1;
  } else {
    // It's somewhat tricky to figure out where this motif goes in the matrix.
    motif_location = find_matrix_location(motifs, motif_id, num_motifs);
  }
  prev_motif_location = find_matrix_location(motifs, prev_motif, num_motifs);

  // Increment the transition count and spacer length matrices.
  incr_matrix_cell(prev_motif_location, motif_location, 1, transp_freq);
  incr_matrix_cell(prev_motif_location, motif_location, 
       motif_position - *prev_position, spacer_ave);

  // Add the occurrence to the order list.
  add_occurrence(motif_id, motif_position - *prev_position, new_order);

  strcpy(prev_motif, motif_id);
}

/***********************************************************************
 * Should the given motif be inserted into the model?
 * FIXME: These tests needn't be mutually exclusive.
 ***********************************************************************/
static BOOLEAN_T retain_motif(
  STRING_LIST_T* requested_motifs, // IDs of motifs to include.
  double         e_threshold,      // E-value to include motifs. 
  double         complexity_threshold, // Complexity threshold to include.
  ORDER_T*       order_spacing,    // Motif order and spacing (linear HMM). 
  MOTIF_T*       motif             // The motif. 
) {
  int num_requested;
  int i;
  char* motif_id;

  /* Method 1: Select motifs by index. */
  num_requested = get_num_strings(requested_motifs);
  if (num_requested > 0) {
    motif_id = get_motif_id(motif);
    for (i = 0; i < num_requested; i++) {
      if (strcmp(get_nth_string(i, requested_motifs), motif_id) == 0) {
        return(TRUE);
      }
    }
    return(FALSE);
  }

  /* Method 2: Select motifs below a certain E-value threshold. */
  else if (e_threshold != 0.0) {
    return (get_motif_evalue(motif) <= e_threshold);
  }

  /* Method 3: Select motifs that are included in the order string. */
  else if (order_spacing != NULL) {
    return order_contains(get_motif_id(motif), order_spacing);
  }

  // Method 4: Select motifs by their complexity score.
  else if (complexity_threshold != 0.0) {
    return(motif->complexity >= complexity_threshold);
  }

  /* Default is to include all motifs. */
  return(TRUE);
}

/***********************************************************************
 * Eliminate motifs that don't meet thresholds or are unused
 ***********************************************************************/
static void filter_motifs(
  STRING_LIST_T* requested_motifs, // Indices of motifs to include. 
  double         e_threshold,      // E-value to include motifs. 
  double         complexity_threshold, // For eliminate low-complexity motifs.
  ORDER_T**      order_spacing,    // Motif order and spacing (linear HMM). 
  int*           num_motifs,       // Number of motifs retrieved. 
  MOTIF_T*       motifs            // The retrieved motifs. 
) {
  int i_motif = 0;
  int num_motifs_kept = 0;

  for(i_motif = 0, num_motifs_kept = 0; i_motif < *num_motifs; i_motif++) {
    if (retain_motif(requested_motifs, e_threshold, complexity_threshold, 
         *order_spacing, &(motifs[i_motif]))) {
      if (i_motif > num_motifs_kept) {
        motifs[num_motifs_kept] = motifs[i_motif];
      }
      num_motifs_kept++;
    }
  }
  
  *num_motifs = num_motifs_kept;

  /* Tell the user which motifs were retained. */
  if (verbosity >= NORMAL_VERBOSE && *num_motifs > 0) {
    fprintf(stderr, "Retaining motif");
    for (i_motif = 0; i_motif < *num_motifs; i_motif++) {
      fprintf(stderr, " %s", get_motif_id(&(motifs[i_motif])));
    }
    fprintf(stderr, ".\n");
  }

  /* No point in proceeding if no motifs were retained. */
  if (*num_motifs == 0) {
    die("No motifs satisfied filters");
  }
    
}

/***********************************************************************
 * Discard motifs that are not connected.
 ***********************************************************************/
static void throw_out_unused_motifs
  (MATRIX_T* transp_freq,
   MATRIX_T* spacer_ave,
   int*      num_motifs,
   MOTIF_T*  motifs)
{
  int i_motif, j_motif;
  ARRAY_T* row_sums;
  ARRAY_T* col_sums;

  // Extract the margins of the transition matrix.
  row_sums = get_matrix_row_sums(transp_freq);
  col_sums = get_matrix_col_sums(transp_freq);

  for (i_motif = 0; i_motif < *num_motifs; i_motif++) {

    // Is this row or column empty?
    if ((get_array_item(i_motif + 1, row_sums) == 0.0) ||
      (get_array_item(i_motif + 1, col_sums) == 0.0)) {

      if (verbosity >= NORMAL_VERBOSE) {
        fprintf(stderr, "Removing unused motif %s. No occurrences of this motif were found.\n",
        get_motif_id(&(motifs[i_motif])));
      }

      // Remove the row and column from the transition matrix.
      remove_matrix_row(i_motif + 1, transp_freq);
      remove_matrix_col(i_motif + 1, transp_freq);
      assert(get_num_rows(transp_freq) == get_num_cols(transp_freq));

      remove_matrix_row(i_motif + 1, spacer_ave);
      remove_matrix_col(i_motif + 1, spacer_ave);
      assert(get_num_rows(spacer_ave) == get_num_cols(spacer_ave));

      // Remove the motif from the array.
      for (j_motif = i_motif + 1; j_motif < *num_motifs; j_motif++) {
        free_motif(&(motifs[j_motif - 1]));
        copy_motif(&(motifs[j_motif]), &(motifs[j_motif - 1]));
      }
      free_motif(&(motifs[j_motif - 1]));
      (*num_motifs)--;
      i_motif--;

      // Recalculate the row and column sums.
      free_array(row_sums);
      free_array(col_sums);
      row_sums = get_matrix_row_sums(transp_freq);
      col_sums = get_matrix_col_sums(transp_freq);

    }
  }

  free_array(row_sums);
  free_array(col_sums);
}

/***********************************************************************
 * Parse the motif occurrences.
 *
 * Each motif occurence string contains the following items
 *  - sequence id,
 *  - sequence p-value,
 *  - number n of motif occurrences, and
 *  - length of sequence.
 *
 * This is followed by n triples containing
 *  - motif id,
 *  - occurrence position, and
 *  - occurrence p-value.
 *
 ***********************************************************************/
static void parse_motif_occurrences(
  STRING_LIST_T*  motif_occurrences, // List of motif occurences OUT
  BOOLEAN_T  has_reverse_strand, // File included both strands? IN
  double     p_threshold,        // P-value to include motif occurences. OUT
  ORDER_T**  order_spacing,      // Motif order and spacing (linear HMM) 
                                 // IN OUT.
  MATRIX_T** transp_freq,        // Motif-to-motif transitions. OUT
  MATRIX_T** spacer_ave,         // Average inter-motif distances. OUT
  int        num_motifs,         // Number of motifs retrieved. IN
  MOTIF_T*   motifs              // The retrieved motifs. IN
) {

  ORDER_T*  new_order;      // New order and spacing. 
  BOOLEAN_T find_order;     // Should we look for the motif order? 

  // If we already have a motif order and spacing, don't find any more. 
  if (*order_spacing == NULL) {
    find_order = TRUE;
  } else {
    find_order = FALSE;
  }
  new_order = NULL;
  
  // Allocate the matrices. 
  *transp_freq = allocate_matrix(num_motifs + 2, num_motifs + 2);
  *spacer_ave = allocate_matrix(num_motifs + 2, num_motifs + 2);
  init_matrix(0.0, *transp_freq);
  init_matrix(0.0, *spacer_ave);

  int num_occurrence_strings = get_num_strings(motif_occurrences);
  int i;
  for (i = 0; i < num_occurrence_strings; i++) {
    char*  sequence_id;       // ID of the current sequence. 
    float  sequence_p;        // pvalue of the entire sequence. 
    int    num_occurs;        // Number of motif occurences in this sequence. 
    int    seq_length;        // Length of the current sequence. 
    int    i_occur;           // Index of the current occurrence. 
    char   prev_motif[MAX_MOTIF_ID_LENGTH + 1]; // Index of the previous motif. 
    int    prev_position;     // Location of the right edge of previous motif. 
    float  motif_p;           // P-value of the current occurrence. 
    char *c;          // Dummy to hold return of strtok.

    char* line = get_nth_string(i, motif_occurrences);

    /* Read the sequence identifier, p-value, number of occurrences
       and length. */
    // tlb; sscanf crashes if strtok returns NULL so pass it "" then
    sequence_id = strtok(line, " ");
    if (sequence_id == NULL) {
      die("Error reading motif occurrences.\n%s", line);
    }
    if (sscanf((c=strtok(NULL, " "))?c:"", "%f", &sequence_p) != 1) {
      die("Can't read p-value of sequence %s.", sequence_id);
    }
    if (sscanf((c=strtok(NULL, " "))?c:"", "%d", &num_occurs) != 1) {
      die("Can't read number of motif occurences in sequence %s.",
    sequence_id);
    }
    if (sscanf((c=strtok(NULL, " "))?c:"", "%d", &seq_length) != 1) {
      die("Can't read length of sequence %s.", sequence_id);
    }

    if (verbosity > NORMAL_VERBOSE) {
      fprintf(stderr, "Reading motif occurrences for sequence %s.\n", 
        sequence_id);
    }

    // If requested, try to create an order string. 
    if (find_order) {
      new_order = create_empty_order(num_occurs, sequence_p);
    }

    // Accumulate motif occurence data. 
    sprintf(prev_motif, "%d", 0);
    prev_position = 0;
    for (i_occur = 0; i_occur < num_occurs; i_occur++) {
      char  motif_id[MAX_MOTIF_ID_LENGTH + 1]; // ID of the current motif. 
      int   motif_position;    // Position of the current motif occurrence. 
      char *c;       // Dummy to hold return of strtok.
      
      // Read the three values. 
      if (sscanf((c=strtok(NULL, " "))?c:"", "%s", motif_id) != 1) {
        die("Can't read index of occurrence %d in sequence %s.",
        i_occur, sequence_id);
      }
      if (sscanf((c=strtok(NULL, " "))?c:"", "%d", &motif_position) != 1) {
        die("Can't read position of occurrence %d in sequence %s.",
        i_occur, sequence_id);
      }
      if (sscanf((c=strtok(NULL, " "))?c:"", "%f", &motif_p) != 1) {
        die("Can't read p-value of occurrence %d in sequence %s.",
        i_occur, sequence_id);
      }

      // Only include motifs that have been retained
      if (have_motif(motif_id, num_motifs, motifs)) {
        // Make sure we have strand information in the ID.
        if (has_reverse_strand) {
          add_strand(motif_id);
        }
  
        // Record this occurrence.
        record_occurrence(sequence_id,
          motif_id,
          p_threshold,
          motif_p,
          prev_motif,
          &prev_position,
          motif_position,
          *transp_freq,
          *spacer_ave,
          new_order,
          num_motifs, 
          motifs);
  
        /* Motifs are stored in order of their motif IDs, but they are
            indexed from zero rather than one. */
        prev_position = motif_position +
          (motifs[find_matrix_location(motifs, motif_id, num_motifs) - 1]).length;
      }
    }
  
    assert(seq_length >= prev_position);
  
    // Record the transition to the end state.
    record_occurrence(sequence_id,
          END_TRANSITION,
          p_threshold,
          motif_p,
          prev_motif,
          &prev_position,
          seq_length,
          *transp_freq,
          *spacer_ave,
          new_order,
          num_motifs,
          motifs);

    // Decide whether to keep the new order object. 
    if (find_order) {
      if ((get_num_distinct(new_order) > get_num_distinct(*order_spacing)) ||
         (((get_num_distinct(new_order) == get_num_distinct(*order_spacing))
         && (get_pvalue(new_order) < get_pvalue(*order_spacing))))) {
            if (verbosity > NORMAL_VERBOSE) {
              fprintf(stderr, "Storing order from sequence %s (%g < %g).\n",
                sequence_id, get_pvalue(new_order), 
                get_pvalue(*order_spacing));
              print_order_and_spacing(stderr, new_order);
            }
            free_order(*order_spacing);
            *order_spacing = new_order;
      } else {
        free_order(new_order);
      }
    }
  }
}

/***********************************************************************
 *  Select the motifs used to build the model, parse any motif
 *  occurences, build the motif order object, and the motif
 *  and spacer frequency matrices.
 ***********************************************************************/
void process_raw_motifs_for_model(
     int* num_motifs,                  // Number of motifs. IN, OUT
     MOTIF_T* motifs,                  // Array of motifs IN, OUT
     STRING_LIST_T* motif_occurrences, // List of motif occurrences. OUT
     STRING_LIST_T* requested_motifs,  // Explicitly requested motifs. IN
     BOOLEAN_T has_reverse_strand,     // Did file contain both strands? IN
     BOOLEAN_T keep_unused,            // Retain unsed motifs? IN
     double p_threshold,               // Motif p-value threshold IN
     double e_threshold,               // Motif e-value threshold IN
     double complexity_threshold,      // Motif complexity threshold IN
     ORDER_T** order_spacing,          // Motif/spacer order IN, OUT
     MATRIX_T** transp_freq,           // Motif transition freqs OUT
     MATRIX_T** spacer_ave,            // Spacer transition freqs OUT
     double trans_pseudo,              // Motif transition pseudo-counts IN
     double spacer_pseudo              // Spacer transition pseudo-counts IN
) {

  // If both strands, make reverse complements.
  if (has_reverse_strand) {
    add_reverse_complements(num_motifs, motifs);
  }

  /* Remove motifs not allowed by the command line parameters */
  filter_motifs(
    requested_motifs, 
    e_threshold, 
    complexity_threshold, 
    order_spacing,  
    num_motifs, 
    motifs
  );

  /* Turn the raw motifs and motif occurences into the */
  /* elements of the model */
  if (motif_occurrences != NULL && get_num_strings(motif_occurrences) > 0) {
    parse_motif_occurrences(
       motif_occurrences,
       has_reverse_strand,
       p_threshold,
       order_spacing,
       transp_freq,
       spacer_ave,
       *num_motifs,
       motifs
    );
  }
  else {
    // If no occurrences are found, initialize matrices uniformly.
    compute_naive_transitions_and_spacers(
      *num_motifs, 
      transp_freq, 
      spacer_ave
    );
  }

  // Convert spacer info to probabilities.
  normalize_spacer_counts(
        trans_pseudo, 
        spacer_pseudo,
        keep_unused,
        *transp_freq, 
        *spacer_ave);

  // Throw out unused motifs.
  throw_out_unused_motifs(*transp_freq, *spacer_ave, num_motifs, motifs);
}

#ifdef MAIN
#include "simple-getopt.h"

VERBOSE_T verbosity = NORMAL_VERBOSE;


/*************************************************************************
 * int main
 *************************************************************************/
int main(int argc, char *argv[])
{
  /* Data structures. */
  int       num_motifs;         /* The number of motifs in the model. */
  MOTIF_T   motifs[2 * MAX_MOTIFS]; /* The motifs. */
  STRING_LIST_T* motif_occurrences = NULL; /* Strings describing occurrences of
                                              motifs */
  BOOLEAN_T has_reverse_strand = FALSE;    /* MEME file contained both strands */
  ARRAY_T*  background;         /* Background probs for alphabet. */
  ORDER_T*  order_spacing;      /* Linear HMM order and spacing. */
  MATRIX_T* transp_freq = NULL; /* Matrix of inter-motif transitions freqs. */
  MATRIX_T* spacer_ave = NULL;  /* Matrix of average spacer lengths. */
  MHMM_T *  the_hmm = NULL;     /* The HMM being constructed. */

  /* Command line parameters. */
  char *    meme_filename;      /* Input file containg motifs. */
  char *    hmm_type_str;       /* HMM type. */
  HMM_T     hmm_type;
  STRING_LIST_T* requested_motifs; /* Indices of requested motifs. */
  int       request_n;          /* The user asked for the first n motifs. */
  double    e_threshold;        /* E-value threshold for motif inclusion. */
  double    complexity_threshold; // For eliminating low-complexity motifs.
  double    p_threshold;        /* p-value threshold for motif occurences. */
  char*     order_string;       /* Motif order and spacing. */
  int       spacer_states;      /* Number of states in each spacer. */
  BOOLEAN_T fim;                /* Represent spacers as free insertion
				   modules? */
  BOOLEAN_T keep_unused;        // Drop unused inter-motif transitions?
  double    trans_pseudo;       /* Transition pseudocount. */
  double    spacer_pseudo;      // Spacer (self-loop) pseudocount. */
  char*     description;        // Descriptive text to be stored in model.
  BOOLEAN_T print_header;       /* Print file header? */
  BOOLEAN_T print_params;       /* Print parameter summary? */
  BOOLEAN_T print_time;         /* Print timing data (dummy: always false). */

  /* Local variables. */
  int       i_motif;

  /**********************************************
   * COMMAND LINE PROCESSING
   **********************************************/
  // Define command line options.
  cmdoption const options[] = {
    {"type", OPTIONAL_VALUE},
    {"description", REQUIRED_VALUE},
    {"motif", REQUIRED_VALUE},
    {"nmotifs", REQUIRED_VALUE},
    {"ethresh", REQUIRED_VALUE},
    {"lowcomp", REQUIRED_VALUE},
    {"pthresh", REQUIRED_VALUE},
    {"order", REQUIRED_VALUE},
    {"nspacer", REQUIRED_VALUE},
    {"fim", NO_VALUE},
    {"keep-unused", NO_VALUE},
    {"transpseudo", REQUIRED_VALUE},
    {"spacerpseudo", REQUIRED_VALUE},
    {"verbosity", REQUIRED_VALUE},
    {"noheader", NO_VALUE},
    {"noparams", NO_VALUE},
    {"notime", NO_VALUE},
    {"quiet", NO_VALUE},
  };
  int option_count = 18;
  int option_index = 0;

  // Define the usage message.
  char      usage[1000] = "";
  strcat(usage, "USAGE: mhmm [options] <MEME file>\n");
  strcat(usage, "\n");
  strcat(usage, "   Options:\n");
  strcat(usage, "     --type [linear|complete|star] (default=linear)\n");
  strcat(usage, "     --description <string> (may be repeated)\n");
  strcat(usage, "     --motif <motif #> (may be repeated)\n");
  strcat(usage, "     --nmotifs <#>\n");
  strcat(usage, "     --ethresh <E-value>\n");
  strcat(usage, "     --lowcomp <value>\n");
  strcat(usage, "     --pthresh <p-value>\n");
  strcat(usage, "     --order <string>\n");
  strcat(usage, "     --nspacer <spacer length> (default=1)\n");
  strcat(usage, "     --fim\n");
  strcat(usage, "     --keep-unused\n");
  strcat(usage, "     --transpseudo <pseudocount>\n");
  strcat(usage, "     --spacerpseudo <pseudocount>\n");
  strcat(usage, "     --verbosity 1|2|3|4|5 (default=2)\n");
  strcat(usage, "     --noheader\n");
  strcat(usage, "     --noparams\n");
  strcat(usage, "     --notime\n");
  strcat(usage, "     --quiet\n");
  strcat(usage, "\n");

  /* Make sure various options are set to NULL or defaults. */
  meme_filename = NULL;
  hmm_type_str = NULL;
  hmm_type = INVALID_HMM;
  requested_motifs = new_string_list();
  request_n = 0;
  e_threshold = 0.0;
  complexity_threshold = 0.0;
  p_threshold = 0.0;
  order_string = NULL;
  spacer_states = DEFAULT_SPACER_STATES,
  fim = FALSE;
  keep_unused = FALSE;
  trans_pseudo = DEFAULT_TRANS_PSEUDO;
  spacer_pseudo = DEFAULT_SPACER_PSEUDO;
  description = NULL;
  print_header = TRUE;
  print_params = TRUE;
  print_time = FALSE;

	simple_setopt(argc, argv, option_count, options);

  // Parse the command line.
  while (1) { 
    int c = 0;
    char* option_name = NULL;
    char* option_value = NULL;
    const char * message = NULL;


    // Read the next option, and break if we're done.
    c = simple_getopt(&option_name, &option_value, &option_index);
    if (c == 0) {
      break;
    } else if (c < 0) {
    	simple_getopterror(&message);
      die("Error processing command line options (%s)\n", message);
    }

    if (strcmp(option_name, "type") == 0) {
			if (option_value != NULL) {
      	hmm_type_str = option_value;
			}
    } else if (strcmp(option_name, "description") == 0) {
      description = option_value;
    } else if (strcmp(option_name, "motif") == 0) {
      add_string(option_value, requested_motifs);
    } else if (strcmp(option_name, "nmotifs") == 0) {
      request_n = atoi(option_value);
    } else if (strcmp(option_name, "ethresh") == 0) {
      e_threshold = atof(option_value);
    } else if (strcmp(option_name, "lowcomp") == 0) {
      complexity_threshold = atof(option_value);
    } else if (strcmp(option_name, "pthresh") == 0) {
      p_threshold = atof(option_value);
    } else if (strcmp(option_name, "order") == 0) {
      order_string = option_value;
    } else if (strcmp(option_name, "nspacer") == 0) {
      spacer_states = atoi(option_value);
    } else if (strcmp(option_name, "fim") == 0) {
      fim = TRUE;
    } else if (strcmp(option_name, "keep-unused") == 0) {
      keep_unused = TRUE;
    } else if (strcmp(option_name, "transpseudo") == 0) {
      trans_pseudo = atof(option_value);
    } else if (strcmp(option_name, "spacerpseudo") == 0) {
      spacer_pseudo = atof(option_value);
    } else if (strcmp(option_name, "verbosity") == 0) {
      verbosity = (VERBOSE_T)atoi(option_value);
    } else if (strcmp(option_name, "noheader") == 0) {
      print_header = FALSE;
    } else if (strcmp(option_name, "noparams") == 0) {
      print_params = FALSE;
    } else if (strcmp(option_name, "notime") == 0) {
      print_time = FALSE;
    } else if (strcmp(option_name, "quiet") == 0) {
      print_header = print_params = print_time = FALSE;
      verbosity = QUIET_VERBOSE;
    }
  }

  // Read the single required argument.
  if (option_index + 1 != argc) {
    fprintf(stderr, "%s", usage);
    exit(1);
  }
  meme_filename = argv[option_index];

  // Set up motif requests. 
  if (request_n != 0) {
    if (get_num_strings(requested_motifs) != 0) {
      die("Can't combine the -motif and -nmotifs options.\n");
    } else {
      for (i_motif = 0; i_motif < request_n; i_motif++) {
        char motif_id[MAX_MOTIF_ID_LENGTH + 1];
        sprintf(motif_id, "%d", i_motif + 1);
        add_string(motif_id, requested_motifs);
      }
    }
  }

  /* Set the model type. */
  hmm_type = convert_enum_type_str(hmm_type_str, LINEAR_HMM, HMM_STRS, 
				   NUM_HMM_T);

  /* Gotta have positive spacer length. */
  if (spacer_states <= 0) {
    die("Negative spacer length (%d).\n", spacer_states);
  }

  /* Make sure motifs weren't selected redundantly. */
  // FIXME: Add tests for complexity threshold.
  if ((get_num_strings(requested_motifs) != 0) && (e_threshold != 0.0)) {
    die("Can't use -motif or -nmotifs with -ethresh.");
  }
  if ((get_num_strings(requested_motifs) != 0) && (order_string != NULL)) {
    die("Can't use -motif or -nmotifs with -order.");
  }
  if ((order_string != NULL) && (e_threshold != 0.0)) {
    die("Can't use -ethresh and -order.");
  }

  /* Prevent trying to build a complete or star model with ordering. */
  if (order_string != NULL) {
    if (hmm_type == COMPLETE_HMM) 
      die("Can't specify motif order with a completely connected model.");
    else if (hmm_type == STAR_HMM)
      die("Can't specify motif order with a star model.");
  } 

  // Parse the order string. 
  order_spacing = create_order(order_string);

  /**********************************************
   * READING THE MOTIFS
   **********************************************/

  BOOLEAN_T read_file = FALSE;
  double pseudocount = 0;

  read_meme_file(
		 meme_filename,
		 "motif-file", // Take bg freq. from motif file.
		 pseudocount,
     REQUIRE_PSPM,
		 &num_motifs,
		 motifs,
		 &motif_occurrences,
		 &has_reverse_strand,
		 &background
		 );

  process_raw_motifs_for_model(
       &num_motifs,
       motifs,
       motif_occurrences,
       requested_motifs,
       has_reverse_strand,
       keep_unused,
       p_threshold,
       e_threshold, 
       complexity_threshold, 
       &order_spacing,
       &transp_freq,
       &spacer_ave,
       trans_pseudo,
       spacer_pseudo
  );

  /**********************************************
   * BUILDING THE HMM
   **********************************************/

  /* Build the motif-based HMM. */
  if (hmm_type == LINEAR_HMM) {

    if (order_spacing != NULL) {
      reorder_motifs(order_spacing, &num_motifs, motifs);
    }
    else {
      die("No order specified for the motifs.\n"
          "For the linear model the motif file must contain motif occurence\n" 
          "data or the motif order must be specified using "
          "the --order option.");
    }

    build_linear_hmm(
      background,
		  order_spacing,
		  spacer_states,
		  motifs,
		  num_motifs, 
		  fim,
		  &the_hmm
    );

  } else if (hmm_type == COMPLETE_HMM) {

    build_complete_hmm(
      background,
		  spacer_states,
		  motifs,
		  num_motifs,
		  transp_freq,
		  spacer_ave,
		  fim,
		  &the_hmm
    );

  } else if (hmm_type == STAR_HMM) {

    build_star_hmm(
      background,
		  spacer_states,
		  motifs,
		  num_motifs,
		  fim,
		  &the_hmm
    );

  }

  // Add some global information.
  copy_string(&(the_hmm->motif_file), meme_filename);

  /**********************************************
   * WRITING THE HMM
   **********************************************/

  /* Print the header. */
  if (print_header)
    write_header(
     program, 
     "",
		 description,
		 meme_filename,
		 NULL,
		 NULL, 
		 stdout
    );

  /* Write the HMM. */
  write_mhmm(verbosity, the_hmm, stdout);

  /* Print the program parameters. */
  if (print_params) {
    printf("Program parameters for mhmm\n");
    printf("  MEME file: %s\n", meme_filename);
    printf("  Motifs:");
    write_string_list(" ", requested_motifs, stdout);
    printf("\n");
    printf("  Model topology: %s\n",
	   convert_enum_type(hmm_type, HMM_STRS, NUM_HMM_T));
    printf("  States per spacer: %d\n", spacer_states);
    printf("  Spacers are free-insertion modules: %s\n",
	   boolean_to_string(fim));
    printf("\n");
  }

  free_array(background);
  free_string_list(requested_motifs);
  free_order(order_spacing);
  free_matrix(transp_freq);
  free_matrix(spacer_ave);
  for (i_motif = 0; i_motif < num_motifs; i_motif++)
    free_motif(&(motifs[i_motif]));
  free_mhmm(the_hmm);
  return(0);
}
#endif /* MAIN */

