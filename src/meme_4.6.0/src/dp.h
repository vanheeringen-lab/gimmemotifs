/**************************************************************************
 * FILE: dp.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 5/20/02
 * PROJECT: MHMM
 * DESCRIPTION: Search a database of sequences using a motif-based HMM.
 **************************************************************************/
#ifndef DP_H
#define DP_H
#include "seq.h"
#include "mhmm-state.h"
#include "match.h"

/**************************************************************************
 * Compute the total likelihood of a sequence, given a model.
 **************************************************************************/
void forward_algorithm
  (BOOLEAN_T local,         // Perform local scoring? (else global)
   int*      int_sequence,  // Sequence to be scored, int format.
   int       seq_length,    // The length of the sequence.
   MHMM_T *  the_log_hmm,   // An HMM in log format.
   MATRIX_T* dp_matrix,     // Dynamic programming matrix.
   MATCH_T*  this_match);   // Only stores total likelihood of sequence.

/**************************************************************************
 * Compute the backward probability matrix.
 **************************************************************************/
void backward_algorithm
  (const int *       int_sequence,  /* Sequence to be scored, int format. */
   const int         seq_length,    /* The length of the sequence. */
   const MHMM_T *    the_log_hmm,   /* An HMM in which the distributions have 
				       already been converted to log format. */
   MATRIX_T*         dp_matrix);    /* Dynamic programming matrix. */

/**************************************************************************
 * Compute the most likely path through a model, given a sequence.
 **************************************************************************/
BOOLEAN_T viterbi_algorithm
  (BOOLEAN_T local,         // Perform local scoring? (else global)
   int*      int_sequence,  // The sequence to be scored, int format.
   int       seq_length,    // The length of the sequence.
   MHMM_T*   the_log_hmm,   // An HMM in which the distributions have
                            // been converted to scaled log format.
   BOOLEAN_T compute_trace, // Save Viterbi path?
   MATRIX_T* motif_score_matrix, // Position-specific scores for each motif.
   MATRIX_T* dp_matrix,     // Dynamic programming matrix.
   MATRIX_T* trace_matrix,  // Traceback matrix.
   MATCH_T*  this_match);   // Most probable path through model.

/**************************************************************************
 * If the forward and backward matrices are correct, then the dot
 * product of corresponding columns in the two matrices should be the
 * same for each state.  This function verifies that equality.
 **************************************************************************/
BOOLEAN_T verify_matrices
  (int       num_states,
   int       max_length, 
   PROB_T    forward_score,
   MATRIX_T* forward_matrix,
   MATRIX_T* backward_matrix);

/**************************************************************************
 * Find repeated matches of a model to a sequence.
 **************************************************************************/
#define NO_REPEAT -999 // If this is dp_threshold, don't do repeat algorithm.
void repeated_match_algorithm
  (PROB_T    dp_threshold,  // Threshold for start of new match.
   int*      int_sequence,  // The sequence to be scored, int format.
   int       seq_length,    // The length of the sequence.
   MHMM_T*   the_log_hmm,   // An HMM in scaled log format.
   MATRIX_T* motif_score_matrix, // Position-specific scores for each motif.
   MATRIX_T* dp_matrix,     // Dynamic programming matrix.
   MATRIX_T* trace_matrix,  // Traceback matrix.
   MATCH_T*  this_match);   // Most probable path through model.

#endif
