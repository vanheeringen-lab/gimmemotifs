/**************************************************************************
 * FILE: match.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 06/01/2002
 * PROJECT: MHMM
 * COPYRIGHT: 2002-2008, WSN
 * DESCRIPTION: Functions for extracting a sequence-to-model match
 *              from a traceback matrix.
 **************************************************************************/
#include <assert.h>
#include "mhmm-state.h"
#include "match.h"

#ifndef MATCH_DEBUG
#define MATCH_DEBUG 0
#endif

#ifndef SHIFT_DEBUG
#define SHIFT_DEBUG 0
#endif

// Initialize global variables associated with enumerated type.
char* SEARCH_STRS[] = {"invalid", "global single", "local single",
		       "global all", "local all", "repeated match"};
int NUM_SEARCH_T = 6;

// Instantiate the match type.
struct match_t {
  int      seq_length;       // The total length of the sequence.
  int      start_match;      // Start of this match, indexed from 0.
  int      end_match;        // End of this match, indexed from 0.
  int      nhits;            // Number of hits (motif-scoring only).
  PROB_T   score;            // Score associated with this match.
  ARRAY_T* trace;            /* An array of length >= seq_length.  Each 
			      * entry is the ID of the corresponding model 
			      * state (or -1 if there is no match here). */
};
    
/**************************************************************************
 * Allocate or de-allocate a match object.
 **************************************************************************/
MATCH_T* allocate_match()
{
  MATCH_T* return_value;
  return_value = (MATCH_T*)mm_malloc(sizeof(MATCH_T));
  return_value->seq_length = -1;
  return_value->start_match = -1;
  return_value->end_match = -1;
  return_value->score = 0.0;
  return_value->trace = allocate_array(1);
  return(return_value);
}

void free_match
  (MATCH_T* this_match)
{
  if (this_match == NULL) {
    return;
  }
  free_array(this_match->trace);
  myfree(this_match);
}

/**************************************************************************
 * Retrieve various fields from the match object.
 **************************************************************************/
PROB_T get_score
  (MATCH_T* this_match)
{
  assert(this_match != NULL);
  return(this_match->score);
}

int get_start_match
  (MATCH_T* this_match)
{
  assert(this_match != NULL);
  return(this_match->start_match);
}

int get_match_seq_length
  (MATCH_T* this_match)
{
  assert(this_match != NULL);
  return(this_match->seq_length);
}

int get_end_match
  (MATCH_T* this_match)
{
  assert(this_match != NULL);
  return(this_match->end_match);
}

int get_nhits
  (MATCH_T* this_match)
{
  assert(this_match != NULL);
  return(this_match->nhits);
}

/**************************************************************************
 * For a given observation, get the index of the state that generated
 * it.
 **************************************************************************/
int get_trace_item
  (int      item,
   MATCH_T* this_match)
{
  assert(this_match != NULL);
  return((int)(get_array_item(item, this_match->trace)));
}

/**************************************************************************
 * The total likelihood of the sequence, given the model, is the sum
 * of the values in the last column of the forward matrix.
 *
 * Note that the specified number of columns is generally smaller
 * than the total number of columns in the DP matrix, because the
 * matrix is allocated to accommodate the longest sequence in the set.
 **************************************************************************/
static PROB_T compute_forward_score
  (BOOLEAN_T local,      // Was local scoring performed?
   int       seq_length, // Number of columns in the matrix.
   MATRIX_T* forward_matrix)
{
  int    start_position; // First column to consider.
  int    i_position;     // Index of the current column.
  PROB_T max_score;      // Maximum log sum.

  // If it is local scoring, consider all positions in the matrix.
  if (local) {
    start_position = 0;
  } else {
    start_position = seq_length - 1;
  }

  // Compute sums for each column.
  max_score = LOG_ZERO;
  for (i_position = start_position; i_position < seq_length; i_position++) {

    // Extract this column from the matrix.
    ARRAY_T* this_column = get_matrix_column(i_position, forward_matrix);

    // Compute the sum.
    PROB_T this_score = log_array_total(this_column);

    // Is this the highest score we've seen so far?
    if (this_score > max_score) {
      max_score = this_score;
    }

    free_array(this_column);
  }

  return(max_score);
}


/**************************************************************************
 * Find the maximum value in a given matrix.  Return that value, along
 * with its coordinates.
 *
 * The number of columns in the matrix are passed in, allowing the method to
 * find the maximum within a submatrix, which is necessary because the
 * dp_matrix has as many columns as the longest sequence in the set.
 **************************************************************************/
static void get_max_cell
  (MATRIX_T* matrix,
   int       num_cols,
   int*      i_max_row,
   int*      i_max_col,
   MTYPE*    max_value)
{
  int num_rows = get_num_rows(matrix);
  int i_row;
  int i_col;
  MTYPE value;

  // Initialize the values we're looking for.
  *i_max_row = 0;
  *i_max_col = 0;
  *max_value = get_matrix_cell(0, 0, matrix);

  // Search for the maximum.
  for (i_row = 0; i_row < num_rows; i_row++) {
    for (i_col = 0; i_col < num_cols; i_col++) {
      value = get_matrix_cell(i_row, i_col, matrix);
      if (value > *max_value) {
	*i_max_row = i_row;
	*i_max_col = i_col;
	*max_value = value;
      }
    }
  }
}

/**************************************************************************
 * Expand a given match to a given sequence length.
 **************************************************************************/
static void set_seq_length
  (int       seq_length,    // Total length of the sequence.
   MATCH_T*  this_match)    // [In/Out] Sequence-to-model match.
{
  // Make sure the object exists.
  myassert(TRUE, this_match != NULL, "Tried to access null match object.");

  // Make sure the match object is big enough to hold this sequence.
  if (get_array_length(this_match->trace) < seq_length+1) {
    free_array(this_match->trace);
    this_match->trace = allocate_array(seq_length+1);
  }
  
  // Store the sequence length.
  this_match->seq_length = seq_length;

}

/**************************************************************************
 * Given a traceback matrix, return the series of states corresponding
 * to the most likely path.
 *
 * Results are returned in 'match', which must be pre-allocated.
 **************************************************************************/
static void get_match
  (int       end_row,       // Row position of end of match.
   int       end_col,       // Column position of end of match.
   MHMM_T*   the_log_hmm,   // Needed for motif widths.
   BOOLEAN_T motif_scoring, // Is motif-scoring in effect?
   MATRIX_T* trace_matrix,  // Matrix of prev state indices.
   MATCH_T*  this_match)    // [In/Out] Sequence-to-model match.
{
  int i_seq;          // Index in the sequence.
  int current_state;  // Index in the model.
  int nhits = 0;      // Number of motif hits.

  // Start at the end of the match.
  current_state = end_row;

  // Move backwards along the sequence.
  for (i_seq = end_col; i_seq >= 0; /* Decremented in loop. */) {
    int i_motif;  // Index in the current "motif".
    int w_motif;  // Width of current "motif".
    int prev_state;

    // Check for termination of local match.
    if (current_state == -1)  break;

    //
    // Record traceback for current "motif" if motif-scoring in effect.
    // Otherwise, just record current state's traceback.
    //
    w_motif = motif_scoring ? the_log_hmm->states[current_state].w_motif : 1;
    for (i_motif = 0; i_motif < w_motif; i_motif++) {
      // Assume that motif states are consecutively numbered.
      set_array_item(i_seq-i_motif, current_state-i_motif, this_match->trace);
    }

    // Keep track of number of motif hits if doing motif scoring.
    if (motif_scoring && 
      the_log_hmm->states[current_state].type == END_MOTIF_STATE) nhits++;
    
    // Get next state.
    prev_state = current_state;
    current_state = (int)(get_matrix_cell(current_state, i_seq, trace_matrix));

    // Update pointer in sequence.
    i_seq -= w_motif;

    // Kludge (currently turned on)
    // Match just ended during motif scoring mode.  Handle zero-length gap
    // case here: Get best state in column preceding motif start by
    // looking at the row zero in the motif start column.  If a motif
    // ended in the immediately preceding column it will be added to
    // the growing match.  No gap penalty is charged.
    //
    if (motif_scoring && current_state == 0 && prev_state != 0) {
      current_state = 
        (int) (get_matrix_cell(current_state, i_seq+1, trace_matrix));
    }

  }

  // Store the dimensions of the match.
  this_match->start_match = i_seq + 1;
  this_match->end_match = end_col;
  this_match->nhits = nhits;

  // Print traceback.
  if (MATCH_DEBUG) { 
    int i;
    fprintf(stderr, "traceback:\n");
    for (i=0; i<=end_col; i++) {
      fprintf(stderr, "%3d ", (int)get_array_item(i, this_match->trace));
      if ((i+1) % 10 == 0) fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
  }
  
}

/**************************************************************************
 * Find the index of the rightmost unmatched position in the sequence.
 *
 * Each entry in the trace is the index of the state that generated
 * the corresponding observation.
 **************************************************************************/
int find_last_unmatched_position
  (MATCH_T* complete_match)
{
  int trace_entry = 0;
  int i_seq;

  // (Subtract 1 due to zero-indexing, and another 1 due to flanking X.)
  for (i_seq = complete_match->seq_length - 2; i_seq > 0; i_seq--) {
    trace_entry = (int)(get_array_item(i_seq, complete_match->trace));
    if (trace_entry == 0) {
      break;
    }
  }

  if (SHIFT_DEBUG) {
    fprintf(stderr, "Unmatched position %d in segment of length %d.\n",
	    i_seq, complete_match->seq_length);
  }
  return(i_seq);
}

/**************************************************************************
 * Extract from a match object the first local match that occurs to
 * the right of the given position.
 *
 * Returns: Boolean indicating whether match is found.
 **************************************************************************/
BOOLEAN_T find_next_match
  (BOOLEAN_T is_complete,       // This match is versus a complete sequence.
   int       start_pos,         // Start looking for match here.
   int       align_width,       // Width of alignment lines.
   MATRIX_T* dp_matrix,
   MHMM_T*   the_log_hmm,   	// Needed for counting motif hits.
   MATCH_T*  complete_match,    // [In] Match against entire sequence.
   MATCH_T*  partial_match)     // [Out] Local match extracted from complete.
{
  BOOLEAN_T in_match = FALSE;
  BOOLEAN_T match_found = FALSE;
  int       start_of_match = -1;
  int       end_of_match = -1;
  int       i_row;
  int       i_col;
  int       num_cols;
  int	    nhits = 0;

  if (is_complete) i_row=0;			// avoid compiler warning.

  // Traverse the global match, looking for the first local match.
  num_cols = get_match_seq_length(complete_match);
  for (i_col = start_pos; i_col < num_cols; i_col++) {
    i_row = get_trace_item(i_col, complete_match);

    // Start of match?
    if (i_row != 0 && !in_match) {
      in_match = TRUE;
      start_of_match = i_col;
    }

    // End of match?
    else if (in_match && 
	     (i_row == 0 || (is_complete && i_col == num_cols-1))) {
      end_of_match = i_col-1;
      match_found = TRUE;
      break;
    }

    // Count number of hits.
    if (the_log_hmm->states[i_row].type == START_MOTIF_STATE) {
      nhits++;
    }
  }

  // If no complete match was found, quit.
  if (!match_found) {
    if (MATCH_DEBUG && (in_match == TRUE)) {
      fprintf(stderr, "Found partial match from %d.\n", start_of_match);
    }
    return(FALSE);
  }
  myassert(TRUE, start_of_match != -1, "Start of match not found.");
  myassert(TRUE, end_of_match != -1, "End of match not found.");
  myassert(TRUE, start_of_match <= end_of_match,
	   "Start of match (%d) is after end of match (%d).",
	   start_of_match, end_of_match);

  // Initialize the local match.
  set_seq_length(get_match_seq_length(complete_match), partial_match);

  // Set the boundaries of the match.
  partial_match->start_match = start_of_match;
  partial_match->end_match = end_of_match;

  // Set the number of motif hits in the match.
  partial_match->nhits = nhits;

  // The local match score is the DP entry at the end, minus the DP entry
  // at the start of the match.
  {
    int start_row = 0;	// No score saved at motif start but saved in row 0.
    int end_row = get_array_item(end_of_match, complete_match->trace);
    double start_score = get_matrix_cell(start_row, start_of_match, dp_matrix);
    double end_score = get_matrix_cell(end_row, end_of_match, dp_matrix);
    partial_match->score = end_score - start_score;

    if (MATCH_DEBUG) {
      fprintf(stderr, "Found a match from %d (%d) to %d (%d).\n",
	      start_of_match, start_row, end_of_match, end_row);
      fprintf(stderr, "start_score=%g end_score=%g total_score=%g\n", 
	      start_score, end_score, partial_match->score);
    }
  }

  // Initialize to -1's only the part of the trace that might be printed.
  for (i_col = MAX(0, start_of_match - align_width); 
       i_col < start_of_match; 
       i_col++) {
    set_array_item(i_col, -1, partial_match->trace);
  }
  for (i_col = end_of_match + 1; 
       i_col <= MIN(end_of_match + align_width + 1, 
		    complete_match->seq_length - 2); 
       i_col++) {
    set_array_item(i_col, -1, partial_match->trace);
  }

  // Transfer the local match.
  for (i_col = start_of_match; i_col <= end_of_match; i_col++) {
    set_array_item(i_col, get_array_item(i_col, complete_match->trace),
		   partial_match->trace);
  }
  
  return(TRUE);
}

/**************************************************************************
 * Find the (local or global) score and, if requested, the
 * corresponding sequence-to-model match.
 *
 * The sequence score and match are returned in this_match, which must
 * be pre-allocated.
 *
 * Returns TRUE iff the score is non-zero.
 **************************************************************************/
BOOLEAN_T traceback_and_scoring
  (SEARCH_T  search_type,    // What kind of search was performed?
   BOOLEAN_T compute_trace,  // Save path?
   int       num_states,     // Number of rows in matrix.
   int       seq_length,     // Number of columns in matrix.
   MHMM_T*   the_log_hmm,    // Needed for motif widths.
   BOOLEAN_T motif_scoring,  // Is motif-scoring in effect?
   MATRIX_T* dp_matrix,      // DP matrix.
   MATRIX_T* trace_matrix,   // Corresponding traceback matrix.
   MATCH_T*  this_match)     // [Out] Trace for this match.
{
  int end_row = 0; // End of match row.
  int end_col = 0; // End of match column.

  // With repeated match, take the upper right corner.
  switch (search_type) {
  case REPEATED_MATCH_SEARCH :
    end_row = 0;
    end_col = seq_length; // One past the end of the sequence.
    this_match->score = get_matrix_cell(end_row, end_col, dp_matrix);
    break;

  // With local Viterbi scoring, find max score in the whole matrix.
  case LOCAL_VITERBI_SEARCH :
    get_max_cell(dp_matrix, seq_length, &end_row, &end_col, 
		 &(this_match->score));
    break;

  // With global Viterbi scoring, take the lower right corner.
  case GLOBAL_VITERBI_SEARCH :
    end_row = num_states - 1;
    end_col = seq_length - 1;
    this_match->score = get_matrix_cell(end_row, end_col, dp_matrix);
    break;

  // Compute column sums for total probability scoring.
  case GLOBAL_TOTAL_SEARCH :
    this_match->score = compute_forward_score(FALSE, seq_length, dp_matrix);
    break;

  case LOCAL_TOTAL_SEARCH :
    this_match->score = compute_forward_score(TRUE, seq_length, dp_matrix);
    break;

  case INVALID_SEARCH :
    die("Invalid search type.\n");
  }

  // No path was found.
  if (this_match->score  <= LOG_SMALL) {
    this_match->score = LOG_ZERO;
    return(FALSE);
  }

  // Store the sequence length.
  set_seq_length(seq_length, this_match);

  // Fill the trace with -1s since this doesn't match any real state.
  init_array(-1, this_match->trace);

  // Get the trace.
  if (compute_trace) {
    get_match(end_row, end_col, the_log_hmm, motif_scoring, trace_matrix, 
	      this_match);
    // get_match may identify a run of identical spacer states as a match
    // We don't want to score matches like this. It is enough to test
    // whether the first state in the match is identical to the last
    // state in the match, since this would never happen in a genuine
    // match.
    int first_state = get_array_item(
      this_match->start_match, 
      this_match->trace
    );
    int last_state = get_array_item(
      this_match->end_match, 
      this_match->trace
    );
    if (first_state == last_state) {
      this_match->score = LOG_ZERO;
      return(FALSE);
    }
  }

  // Otherwise, assume it's a global match.
  else {
    this_match->start_match = 1;
    this_match->end_match = seq_length - 2; // Subtract spaces on either end.
  }

  return(TRUE);
}
