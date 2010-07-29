/**************************************************************************
 * FILE: dp.c
 * AUTHOR: William Stafford Noble, Timothy L. Bailey
 * CREATE DATE: 5/20/02
 * PROJECT: MHMM
 * VERSION: $Revision: 1.1.1.1 $
 * COPYRIGHT: 1998-2002, WNG, 2001, TLB
 * DESCRIPTION: Core dynamic programming algorithms for Meta-MEME.
 **************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>
#include "utils.h"       // Generic utilities. 
#include "matrix.h"      // Routines for floating point matrices. 
#include "array.h"       // Routines for floating point arrays. 
#include "metameme.h"    // Global metameme functions. 
#include "alphabet.h"    // The amino acid / nucleotide alphabet. 
#include "seq.h"         // Biosequence data structure.
#include "mhmm-state.h"  // HMM data structure. 
#include "read-mhmm.h"   // HMM input/output. 
#include "log-hmm.h"     // HMM log/log-odds conversion. 
#include "fitevd.h"   // Extreme value distribution routines. 
#include "match.h"       // Find a match in the traceback matrix.
#include "pssm.h"   // Needed for computing motif scores.
#include "dp.h"

#define FIRST_COL 1    // First "real" (non-padding)  sequence column.

/**************************************************************************
* Print the dp and traceback matrices for debug purposes.
**************************************************************************/
void print_matrices(
  int       num_hot_states,
  int*       hot_states,
  int seq_length,
  MATRIX_T* dp_matrix,
  MATRIX_T* trace_matrix
) {

  int i,j,k;
  for (k=-1; k<num_hot_states; k++) {
    i = k<0 ? 0 : hot_states[k]; 
    fprintf(stderr, "\ndp/tb %3d:", i);
    for (j=0; j<=seq_length; j++) {
      if (j % 10 == 0) fprintf(stderr, "\n%4d    ", j);
      fprintf(stderr, "%6.2f %2d   ", MAX(0, get_matrix_cell(i, j, dp_matrix)),
        (int)get_matrix_cell(i, j, trace_matrix));
    }
    fprintf(stderr, "\n");
  }
} // print_matrices

/**************************************************************************
 * Compute the most likely path through a model, given a sequence.
 **************************************************************************/
#ifndef VITERBI_DEBUG
#define VITERBI_DEBUG 0
#endif
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
   MATCH_T*  this_match)    // Most probable path through model.
{
  int       i_hot;           // Index of hot state.
  int       i_row;           // Index of current row in dp matrix.
  int       i_col;           // Index of current column in dp matrix.
  const int num_states = the_log_hmm->num_states;
  BOOLEAN_T motif_scoring;
  int       num_hot_states = the_log_hmm->num_hot_states;
  int*       hot_states = the_log_hmm->hot_states;
  double    best_score;
  int       best_row;
  int       in_index;

  motif_scoring = (motif_score_matrix != NULL);
  // The matrices should start life with safe values
  init_matrix(LOG_ZERO, dp_matrix); // Minimal score
  init_matrix(-1, trace_matrix); // Invalid state

  // Initialize the first column.
  set_matrix_cell(0, 0, 0.0, dp_matrix);
  set_matrix_cell(0, 0, -1, trace_matrix);
  for (i_row = 1; i_row < num_states; i_row++) {
    set_matrix_cell(i_row, 0, LOG_ZERO, dp_matrix);
    set_matrix_cell(i_row, 0, -1, trace_matrix);
  }

  // Consider each position in the sequence.
  for (i_col = 1; i_col < seq_length; i_col++) {

    // Start state has no predecessors.
    set_matrix_cell(0, i_col, local ? 0 : LOG_ZERO, dp_matrix);
    set_matrix_cell(0, i_col, -1, trace_matrix);

    // Consider each (hot) state in the model.
    for (i_hot = 0; i_hot < num_hot_states; i_hot++) {
      MHMM_STATE_T*  this_state;   // The current state in the model.
      PROB_T         emit_logprob; // Log prob of emitting the current
                                   // character at the current position.
      int       start_col;         // Index of starting column of 
                                   // "motif" in dp matrix.
  
      // Set current row.
      i_row = hot_states[i_hot];
 
      // Get the current model state.
      this_state = &(the_log_hmm->states[i_row]);

      // Ignore start and mid-motif states if doing motif-scoring.
      if (motif_scoring && ((this_state->type == START_MOTIF_STATE) ||
          (this_state->type == MID_MOTIF_STATE))) {
         continue;
      }

      //
      // Compute the emission cost at the current position.
      // Adjust this_state to point to start of motif if doing motif-scoring.
      // Skip this state if the motif would hang off left edge of sequence.
      //
      if (motif_scoring && this_state->type == END_MOTIF_STATE) {
        int w = this_state->w_motif;
        start_col = i_col-w+1;                // Start pos in sequence.
        if (start_col < FIRST_COL) continue;  // Motif doesn't fit.
        this_state -= (w - 1);                // Start of motif.
        emit_logprob = get_matrix_cell(this_state->i_motif, start_col,
          motif_score_matrix);
      } else {
        start_col = i_col;      // Start pos in sequence.
        emit_logprob = get_array_item(int_sequence[i_col], 
              this_state->emit_odds);
      }

      // Find the preceding state with the highest score + transition cost.
      best_score = LOG_ZERO;
      best_row = -1;
      for (in_index = 0; in_index < this_state->ntrans_in; in_index++) {
        double trans;         // Transition cost to get here.
        double prv_score;     // Best score in previous state.
        double score;         // Current score.
        int prv_row;          // Index of the previous state. 
        int prv_col;          // Index of the previous seq. pos. 

        //
        // Get previous state from list of in-transitions to current state.
        //
        // Previous state index.
        prv_row = this_state->itrans_in[in_index];

        // Previous sequence position.
        prv_col = start_col - 1;

        // Transition cost to get here.
        trans = get_array_item(in_index, this_state->trans_in);  

        // Previous state score.
        prv_score = get_matrix_cell(prv_row, prv_col, dp_matrix); 

        // Total score to get here.
        score = prv_score + trans; 

        //
        // Save best predecessor.
        //
        if (score > best_score) {
          best_score = score;
          best_row = prv_row;
        }
      } // preceding state 

      // Compute the score for this cell.
      if (local && best_score < 0.0) {  // Local and no positive predecessor.
        best_score = emit_logprob > 0 ? emit_logprob : 0;
        best_row = -1;
      } else {                          // Positive predecessor or non-local.
        best_score += emit_logprob;
      }
      set_matrix_cell(i_row, i_col, best_score, dp_matrix);
      set_matrix_cell(i_row, i_col, best_row, trace_matrix);

    } // i_row
  } // i_col

  // Traceback and scoring.
  return traceback_and_scoring(
    local ? LOCAL_VITERBI_SEARCH : GLOBAL_VITERBI_SEARCH,
    compute_trace,
    num_states,
    seq_length,
    the_log_hmm,
    motif_scoring,
    dp_matrix,
    trace_matrix,
    this_match
  );
}



/**************************************************************************
 * Compute the total likelihood of a sequence, given a model.
 **************************************************************************/
#ifndef FORWARD_DEBUG
#define FORWARD_DEBUG 0
#endif
void forward_algorithm
  (BOOLEAN_T local,         // Perform local scoring? (else global)
   int*      int_sequence,  // Sequence to be scored, int format.
   int       seq_length,    // The length of the sequence.
   MHMM_T *  the_log_hmm,   // An HMM in log format.
   MATRIX_T* dp_matrix,     // Dynamic programming matrix.
   MATCH_T*  this_match)    // Only stores total likelihood of sequence.
{
  int       i_row;   // Index of current row in dp matrix.
  int       i_col;   // Index of current column in dp matrix.
  const int num_states = the_log_hmm->num_states; // Number of states in HMM.
  //FIXME: handle motif scoring here!
  BOOLEAN_T motif_scoring = FALSE;

  // Initialize the first column to zeroes.
  for (i_row = 1; i_row < num_states; i_row++) {
    set_matrix_cell(i_row, 0, LOG_ZERO, dp_matrix);
  }

  // Initialize the upper left corner to one.
  set_matrix_cell(0, 0, 0.0, dp_matrix); // = LOG2(1.0)

  // Consider each position in the sequence.
  for (i_col = 0; i_col < seq_length - 1; i_col++) {

    // Index of current character in alphabet.
    int this_char_index = int_sequence[i_col];

    // Initialize the next column to small numbers.
    for (i_row = 0; i_row < num_states; i_row++) {
      set_matrix_cell(i_row, i_col + 1, LOG_ZERO, dp_matrix);
    }

    // Consider each state in the model.
    for (i_row = 0; i_row < num_states; i_row++) {
      MHMM_STATE_T* this_state = &(the_log_hmm->states[i_row]);
      
      // Don't bother with unreachable states.
      if (get_matrix_cell(i_row, i_col, dp_matrix) <= LOG_SMALL) {
        continue;
      }

      // If we're in a motif, there is no transition cost.
      if ((this_state->type == START_MOTIF_STATE) ||
        (this_state->type == MID_MOTIF_STATE)) {

         set_matrix_cell(
           i_row + 1, 
           i_col + 1, 
           get_matrix_cell(i_row, i_col, dp_matrix),
           dp_matrix
         );

      } else {
        int i_trans; // Index of current transition out of this state.

        // Added during optimization.
        const int  ntrans_out = this_state->ntrans_out;
        const int* itrans_out = this_state->itrans_out;
        ARRAY_T*   trans_out = this_state->trans_out;

        // Add up the contributions from this state to the next column.
        for (i_trans = 0; i_trans < ntrans_out; i_trans++) {
        
          // Figure out where we're going to transition to.
          const int next_state = itrans_out[i_trans];
          
          // Figure out the log probability of reaching that state.
          const PROB_T total_logprob 
            = get_matrix_cell(i_row, i_col, dp_matrix) + 
            get_array_item(i_trans, trans_out);

          // Store the cost in the next column.
          set_matrix_cell(next_state, i_col + 1, 
              LOG_SUM(get_matrix_cell(next_state, i_col + 1, 
                    dp_matrix), 
                total_logprob), dp_matrix);

          if (FORWARD_DEBUG) {
            fprintf(stderr, 
              "Adding %g (= %g + %g) for transition from %d to %d.\n",
              total_logprob, get_matrix_cell(i_row, i_col, dp_matrix),
              get_array_item(i_trans, trans_out), i_row, next_state);
          }
        }
      }
    }

    // Add the emission probabilities to the next column.
    this_char_index = int_sequence[i_col + 1];
    if (i_col + 1 != seq_length - 1) { // Skip flanking Xs
      for (i_row = 0; i_row < num_states; i_row++) {
        MHMM_STATE_T* this_state = &(the_log_hmm->states[i_row]);

        PROB_T this_emit = get_array_item(this_char_index, 
            this_state->emit_odds);

        incr_matrix_cell(i_row, i_col + 1, this_emit, dp_matrix);

        if (FORWARD_DEBUG) {
          fprintf(stderr, "Adding %g for char %d in [%d][%d] -> %g.\n", 
            this_emit, this_char_index, i_row, i_col + 1,
            get_matrix_cell(i_row, i_col + 1, dp_matrix));
        }

        // Set negative values to zero if we're doing local scoring.
        if ((local) && (get_matrix_cell(i_row, i_col + 1, dp_matrix) < 0.0)) {
          set_matrix_cell(i_row, i_col + 1, 0.0, dp_matrix);

          if (FORWARD_DEBUG) {
            fprintf(stderr, "Setting negative value to zero at [%d][%d].\n",
              i_row, i_col + 1);
          }
        }
      }
    }
  }

  // Compute the sum of the last column.
  traceback_and_scoring(
    local ? LOCAL_TOTAL_SEARCH : GLOBAL_TOTAL_SEARCH,
    FALSE, // Don't compute traceback.
    num_states,
    seq_length,
    the_log_hmm,
    motif_scoring,
    dp_matrix,
    NULL,
    this_match
  );
}

/**************************************************************************
 * Compute the backward probability matrix.
 **************************************************************************/
void backward_algorithm
  (const int *       int_sequence,  // Sequence to be scored, int format. 
   const int         seq_length,    // The length of the sequence. 
   const MHMM_T *    the_log_hmm,   /* An HMM in which the distributions have 
               already been converted to log format. */
   MATRIX_T*         dp_matrix)     // Dynamic programming matrix. 
{
  int            i_row;           // Index of current row in dp matrix. 
  int            i_col;           // Index of current column in dp matrix. 
  int            i_trans;         /* Index of current transition out of this 
             state. */
  int            this_char_index; // Index of current character in alphabet. 
  MHMM_STATE_T * this_state;      // The current state in the model. 
  int            prev_state;      // Index of the previous state.
  PROB_T         emit_logprob;    /* Log prob of emitting the current
             character at the current position. */
  PROB_T         before_logprob;  // emit_logprob + current_logprob. 
  PROB_T         total_logprob;   /* before_logprob + log prob of transition
             from current state to the next state. */

  // Added during optimization. 
  const int      num_states = the_log_hmm->num_states;
  int            ntrans_in;
  int *          itrans_in;
  ARRAY_T*       trans_in;

  // Initialize the last column to 1.0 (0.0 in log space). 
  for (i_row = 0; i_row < num_states; i_row++) {
    set_matrix_cell(i_row, seq_length-1, 0.0, dp_matrix);
  }

  // Consider each position in the sequence. 
  for (i_col = seq_length-1; i_col > 0; i_col--) {
    this_char_index = int_sequence[i_col];

    // Initialize the preceding column to small numbers. 
    for (i_row = 0; i_row < num_states; i_row++) {
      set_matrix_cell(i_row, i_col-1, LOG_ZERO, dp_matrix);
    }

    // Consider each state in the model. 
    for (i_row = 0; i_row < num_states; i_row++) {
      
      // Don't bother with unreachable states. 
      if (get_matrix_cell(i_row, i_col, dp_matrix) <= LOG_SMALL) {
  continue;
      }

      // Compute the emission cost at the current position. 
      this_state = &(the_log_hmm->states[i_row]);
      if (i_col == seq_length-1) {
  emit_logprob = 0.0; // Skip flanking Xs 
      } else {
  emit_logprob
      = get_array_item(this_char_index, this_state->emit_odds);
      }

      /* Compute the total cost at this position after emitting but
   before transitioning. */
      before_logprob = get_matrix_cell(i_row, i_col, dp_matrix) + 
  emit_logprob;

   // If we're in a motif, there is no transition cost. 
      if ((this_state->type == END_MOTIF_STATE) ||
    (this_state->type == MID_MOTIF_STATE)) {
  set_matrix_cell(i_row - 1, i_col - 1, before_logprob, dp_matrix);

      } else {

  // Add up the contributions from this state to the next column. 
  ntrans_in = this_state->ntrans_in;
  itrans_in = this_state->itrans_in;
  trans_in = this_state->trans_in;
  for (i_trans = 0; i_trans < ntrans_in; i_trans++) {
  
    // Figure out where we're transitioning from. 
    prev_state = itrans_in[i_trans];
    
    // Figure out the log probability of reaching that state. 
    total_logprob = before_logprob + get_array_item(i_trans, trans_in);

    // Store the cost in the previous column. 
    if (get_matrix_cell(prev_state, i_col-1, dp_matrix) 
        <= LOG_SMALL) {
      set_matrix_cell(prev_state, i_col-1, total_logprob, dp_matrix);
    } else {
      set_matrix_cell(prev_state, i_col-1, 
          LOG_SUM(get_matrix_cell(prev_state, i_col-1,
                dp_matrix),
            total_logprob), dp_matrix);
    }
  }
      }
    }
  }
}

/**************************************************************************
 * If the forward and backward matrices are correct, then the dot
 * product of corresponding columns in the two matrices should be the
 * same for each state.  This function verifies that equality.
 **************************************************************************/
#define SLOP 1E-6
BOOLEAN_T verify_matrices
  (int       num_states,
   int       max_length, 
   PROB_T    forward_score,
   MATRIX_T* forward_matrix,
   MATRIX_T* backward_matrix)
{
  PROB_T this_value;       // (Forward * backward) for a particular cell. 
  PROB_T this_total;       // Sum of (forward * backward) for one column. 
  int    i_seq;            // Loop index. 
  int    i_state;          // Loop index. 

  for (i_seq = 0; i_seq < max_length; i_seq++) {

    // Add up (forward * backward) for each state. 
    this_total = LOG_ZERO;
    for (i_state = 0; i_state < num_states; i_state++) {
      this_value = get_matrix_cell(i_state, i_seq, forward_matrix) +
  get_matrix_cell(i_state, i_seq, backward_matrix);
      this_total = LOG_SUM(this_total, this_value);
    }

    // Check to see if the result is close to what we expect. 
    if (!almost_equal(forward_score, this_total, SLOP)) {
      fprintf(stderr, 
        "Matrix error for position %d: difference=%.4g.\n", i_seq,
        forward_score - this_total);
      abort();
      return(FALSE);
    }
  }
  return(TRUE);
}

/**************************************************************************
 * Find repeated matches of a model to a sequence.
 *
 * Let model have m states.
 * Let w(i) be the width of state i (= 1 for gap, else motif width).
 * let ps(i) be the predecessors of state i.
 * Let t(i,j) be the transition COST from i to j.
 * Let e(i,j) be the emission SCORE for aligning state i ENDING at position j.
 *
 * F(i,0) = 0, i= 1 to m
 *
 * The following are applied column-by-column, left to right, sequentially:
 *
 * F(0,j) = max F(0,j-1)
 *              F(i,j-1) - R, i= 1 to m
 *
 * F(i,j) = max F(0,j-w(i)+1) + e(i,j)
 *              F(k,j-w(i)) - t(k,i) + e(i,j), where k in ps(i)
 *
 **************************************************************************/
void repeated_match_algorithm
  (PROB_T    dp_threshold,  // Threshold for start of new match.
   int*      int_sequence,  // The sequence to be scored, int format.
   int       seq_length,    // The length of the sequence.
   MHMM_T*   the_log_hmm,   // An HMM in scaled log format.
   MATRIX_T* motif_score_matrix, // Position-specific scores for each motif.
   MATRIX_T* dp_matrix,     // Dynamic programming matrix.
   MATRIX_T* trace_matrix,  // Traceback matrix.
   MATCH_T*  this_match)    // Most probable path through model.
{
  int      i_hot;       // Index of current hot state.
  int       i_row;           // Index of current row in dp matrix.
  int       i_col;           // Index of current column in dp matrix.
  int       in_index;
  BOOLEAN_T motif_scoring = (motif_score_matrix != NULL);
  int       num_hot_states = the_log_hmm->num_hot_states;
  int*       hot_states = the_log_hmm->hot_states;
  const int num_states = the_log_hmm->num_states;

  // Initialize the upper left corner to one.
  set_matrix_cell(0, 0, 0.0, dp_matrix); // = LOG2(1.0)

  // Initialize the first column to 0.
  // Note: this column corresponds to the padding "X" at a start of
  // sequence.
  for (i_row = 1; i_row < num_states; i_row++) {
    set_matrix_cell(i_row, 0, LOG_ZERO, dp_matrix);
    set_matrix_cell(i_row, 0, 0, trace_matrix);
  }

  //
  // Consider each position in the sequence, starting with the second
  // position.
  //
  for (i_col = FIRST_COL; i_col <= seq_length; i_col++) {
    PROB_T top_row_value; // Value in the top row.
    int    top_row_index; // Index of row from which top-row value came.
    const int prv_col = i_col - 1; // Previous column.

    //
    // Extend an unmatched region or end the current match if
    // it exceeds the DP threshold.
    //
    // F(0,j) = max F(0,j-1)
    //              F(i,j-1) - R, i= 1 to m

    // Best score if previous column not in a match, too.
    top_row_value = get_matrix_cell(0, prv_col, dp_matrix);
    top_row_index = 0;
    // Get best score if a match ended in the previous column.
    for (i_hot=0; i_hot<num_hot_states; i_hot++) { // state
      PROB_T this_value;

      // Get index of current state.
      i_row = hot_states[i_hot];
      
      this_value = get_matrix_cell(i_row, prv_col, dp_matrix) - dp_threshold;

      if (top_row_value < this_value) {
        top_row_value = this_value;
        top_row_index = i_row;
      }
    } // state

    //
    // Store this value in end-state row.
    //
    set_matrix_cell(0, i_col, top_row_value, dp_matrix);
    set_matrix_cell(0, i_col, top_row_index, trace_matrix);

    // 
    // Done after setting start-state entry if past end of sequence.
    //
    if (i_col == seq_length) break;

    //
    // Find best sum of match scores assuming x_{i_col} is in 
    // a matched region.
    //
    // Consider "hot" state in the model.
    //
    for (i_hot=0; i_hot<num_hot_states; i_hot++) { 

      // The current state in the model.      
      MHMM_STATE_T* this_state;      

      // Log prob of emitting the current character at the current position.
      PROB_T        emit_logprob;    

      // Index of starting column of "motif" in dp matrix.
      int           start_col;

      // Row zero score for being in a non-match region.
      double       no_match_score;
   
      // Get index of current state.
      i_row = hot_states[i_hot];

      // Get the current model state.
      this_state = &(the_log_hmm->states[i_row]);

      //
      // Compute the emission score at the current position.
      // If motif-scoring:
      //  1) Adjust this_state to point to start of motif state.
      //  2) Adjust start_col to point to start of motif in sequence.
      //  3) Skip this state if the motif would hang off left edge of sequence.
      //
      if (motif_scoring && this_state->type == END_MOTIF_STATE) {
        int w = this_state->w_motif;
        start_col = i_col-w+1;      // Start pos in sequence.

        // Check to be sure motif fits.
        if (start_col < FIRST_COL) continue;

        this_state -= (w - 1);      // Start of motif.
        emit_logprob = get_matrix_cell(this_state->i_motif, start_col,
               motif_score_matrix);
      } else {
        start_col = i_col;
        emit_logprob = get_array_item(int_sequence[i_col], 
              this_state->emit_odds);
      }

      //
      // F(i,j) = max F(0,j-w(i)+1) + e(i,j)
      //              F(k,j-w(i)) - t(k,i) + e(i,j), where k in ps(i)
      //

      no_match_score = get_matrix_cell(0, start_col, dp_matrix);
      top_row_value = no_match_score;   
      top_row_index = 0;

      for (in_index = 0; in_index < this_state->ntrans_in; in_index++) {
        double trans;            // Transition cost to get here.
        double prv_score;        // Best score in previous state.
        double score;            // Current score.
        int prv_row;             // Index of the previous state.
        int prv_col;             // Index of the previous seq. pos.

        //
        // Get previous state from list of in-transitions to current state.
        //

        // Previous state index.
        prv_row = this_state->itrans_in[in_index];

        // Previous sequence position.
        prv_col = start_col - 1;              

        // Previous state score.
        prv_score = get_matrix_cell(prv_row, prv_col, dp_matrix);

        // Transition cost to get here.
        trans = get_array_item(in_index, this_state->trans_in); 

        // Total score to get here.
        score = prv_score + trans;

        //
        // Save best predecessor.
        //
        if (score > top_row_value) {
          top_row_value = score;
          top_row_index = prv_row;
        }
      } // preceding state
      top_row_value += emit_logprob;

      // Set this cell's values.
      set_matrix_cell(i_row, i_col, top_row_value, dp_matrix);
      set_matrix_cell(i_row, i_col, top_row_index, trace_matrix);

    } // i_hot
  } // i_col

  // Traceback and scoring.
  traceback_and_scoring(
    REPEATED_MATCH_SEARCH,
    TRUE,   // Compute the trace.
    num_states,
    seq_length,
    the_log_hmm,
    motif_scoring,
    dp_matrix,
    trace_matrix,
    this_match
  );
} // repeated_match
