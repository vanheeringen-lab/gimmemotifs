/**************************************************************************
 * FILE: match.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 06/01/2002
 * PROJECT: MHMM
 * COPYRIGHT: 2002-2008 WSN 
 * DESCRIPTION: Functions for extracting a sequence-to-model match
 * from a traceback matrix.
 **************************************************************************/
#ifndef MATCH_H
#define MATCH_H
#include "utils.h"
#include "matrix.h"

// Enumerated type for different search algorithms.
typedef enum {INVALID_SEARCH, GLOBAL_VITERBI_SEARCH, LOCAL_VITERBI_SEARCH,
	      GLOBAL_TOTAL_SEARCH, LOCAL_TOTAL_SEARCH, REPEATED_MATCH_SEARCH}
SEARCH_T;
extern char* SEARCH_STRS[];
extern int NUM_STATE_T;

// One sequence-to-model match.
typedef struct match_t MATCH_T;

/**************************************************************************
 * Allocate or de-allocate a match object.
 **************************************************************************/
MATCH_T* allocate_match();
void free_match
  (MATCH_T* this_match);

/**************************************************************************
 * Retrieve various fields from the match object.
 **************************************************************************/
PROB_T get_score
  (MATCH_T* this_match);
int get_match_seq_length
  (MATCH_T* this_match);
int get_start_match
  (MATCH_T* this_match);
int get_end_match
  (MATCH_T* this_match);
int get_nhits
  (MATCH_T* this_match);

/**************************************************************************
 * For a given observation, get the index of the state that generated
 * it.
 **************************************************************************/
int get_trace_item
  (int      item,
   MATCH_T* this_match);

/**************************************************************************
 * Find the index of the last unmatched position in the sequence.
 **************************************************************************/
int find_last_unmatched_position
  (MATCH_T* complete_match);

/**************************************************************************
 * Extract from a match object the first local match that occurs to
 * the right of the given position.  If no match is found, returns
 * FALSE.
 **************************************************************************/
BOOLEAN_T find_next_match
  (BOOLEAN_T is_complete,       // This match is versus a complete sequence.
   int start_pos,	 	// Start looking for match here.
   int align_width,		// Width of alignment lines.
   MATRIX_T* dp_matrix,
   MHMM_T*   the_log_hmm,   	// Needed for counting motif hits.
   MATCH_T*  complete_match,
   MATCH_T*  partial_match);
    
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
   MATCH_T*  this_match);    // [Out] Trace for this match.

#endif
