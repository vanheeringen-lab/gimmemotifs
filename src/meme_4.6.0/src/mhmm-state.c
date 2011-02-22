/***********************************************************************
 * FILE: mhmm-state.c
 * AUTHOR: William Stafford Noble
 * PROJECT: MHMMM
 * COPYRIGHT: 1997-2008, WSN
 * DESCRIPTION: Data structure for representing one state in an HMM.
 ***********************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "utils.h"
#include "matrix.h"
#include "array.h"
#include "params.h"
#include "alphabet.h"
#include "pssm.h"
#include "mhmm-state.h"

/* Initialize global variables associated with enumerated types. */
char * STATE_STRS[] = {"invalid", "start", "start_motif", "mid_motif",
		       "end_motif", "spacer", "end"};
int NUM_STATE_T = 7;

char * HMM_STRS[] = {"invalid", "linear", "complete", "star"};
int NUM_HMM_T = 4;

/************************************************************************
 * See .h file for description.
 ************************************************************************/
void allocate_mhmm
  (int      num_states,
   MHMM_T** an_mhmm)
{
  /* First allocate the struct itself. */
  *an_mhmm = (MHMM_T*)mm_malloc(sizeof(MHMM_T));

  /* Then allocate the states. */
  (*an_mhmm)->states
    = (MHMM_STATE_T *)mm_calloc(num_states, sizeof(MHMM_STATE_T));

  /* ... and the transition matrix. */
  (*an_mhmm)->trans = allocate_matrix(num_states, num_states);

  /* Allocate the alphabet. */
  (*an_mhmm)->alphabet
    = (char *)mm_calloc(sizeof(char), get_alph_size(ALPH_SIZE) + 1);
  (*an_mhmm)->alph_size = get_alph_size(ALPH_SIZE);

  // Allocate the background distribution.
  (*an_mhmm)->background = allocate_array(get_alph_size(ALL_SIZE));

  // Set other stuff to NULL.
  (*an_mhmm)->description = NULL;
  (*an_mhmm)->motif_file = NULL;
  (*an_mhmm)->sequence_file = NULL;
  (*an_mhmm)->hot_states = NULL;
  (*an_mhmm)->num_hot_states = 0;
}

/************************************************************************
 * Get the motif ID, with or without strand info.
 *
 * Returns a pointer to local static storage.
 ************************************************************************/
#define MAX_MOTIF_ID 50
char* get_state_motif_id
  (BOOLEAN_T     stranded,
   MHMM_STATE_T* this_state)
{
  static char return_value[MAX_MOTIF_ID];

  // Didn't ask for strand, but strand is in the motif ID, then strip it off.
  if ((!stranded) && 
      ((this_state->motif_id[0] == '+') ||
       (this_state->motif_id[0] == '-'))) {
    strcpy(return_value, &((this_state->motif_id)[1]));
  } else {
    strcpy(return_value, this_state->motif_id);
  }
  return(return_value);
}


/************************************************************************
 * Determine the strandedness of a given motif.  Returns '+' or '-' if
 * the motif has strand and '.' otherwise.
 *************************************************************************/
char get_strand
  (MHMM_STATE_T* this_state)
{
  if (this_state->motif_id[0] == '+') {
    return('+');
  } else if (this_state->motif_id[0] == '-') {
    return('-');
  }
  return('.');
}


/************************************************************************
 * Copy one state of an MHMM.
 ************************************************************************/
static void copy_state
  (MHMM_STATE_T *  a_state,
   MHMM_STATE_T *  new_state)
{
  new_state->type = a_state->type;

  /* Allocate memory. */
  new_state->itrans_out = (int *)mm_malloc(a_state->ntrans_out * sizeof(int));
  new_state->trans_out = allocate_array(a_state->ntrans_out);
  new_state->itrans_in = (int *)mm_malloc(a_state->ntrans_in * sizeof(int));
  new_state->trans_in = allocate_array(a_state->ntrans_in);
  new_state->emit = allocate_array(get_alph_size(ALL_SIZE));
  new_state->emit_odds = allocate_array(get_alph_size(ALL_SIZE));

  /* Outgoing transitions. */
  new_state->ntrans_out = a_state->ntrans_out;
  copy_int_array(a_state->ntrans_out, 
		 a_state->itrans_out, 
		 new_state->itrans_out);
  copy_array(a_state->trans_out, new_state->trans_out);

  /* Incoming transitions. */
  new_state->ntrans_in = a_state->ntrans_in;
  copy_int_array(a_state->ntrans_in,
		 a_state->itrans_in,
		 new_state->itrans_in);
  copy_array(a_state->trans_in, new_state->trans_in);

  /* Emissions. */
  copy_array(a_state->emit, new_state->emit);
  copy_array(a_state->emit_odds, new_state->emit_odds);
  new_state->num_sites = a_state->num_sites;
  
  // Descriptive information.
  new_state->i_motif = a_state->i_motif;
  new_state->w_motif = a_state->w_motif;
  strcpy(new_state->motif_id, a_state->motif_id);
  new_state->i_position = a_state->i_position;
  new_state->id_char = a_state->id_char;
}

/************************************************************************
 * Mix two HMM states, according to a given mixing parameter.
 *
 * Assumes both states are already allocated and that both are in log
 * form.
 ************************************************************************/
static void mix_states
  (const float    mixing,    /* Percent of old to be retained. */
   MHMM_STATE_T * new_state,
   MHMM_STATE_T * old_state)
{
  new_state->type = old_state->type;

  /* Emissions. */
  mix_log_arrays(mixing, new_state->emit, old_state->emit);
  mix_log_arrays(mixing, new_state->emit_odds, old_state->emit_odds);
  new_state->num_sites = (mixing * old_state->num_sites)
    + ((1.0 - mixing) * new_state->num_sites);

  /* Transitions will be mixed in the transition array and updated later. */
}

/************************************************************************
 * See .h file for description.
 ************************************************************************/
void copy_mhmm
  (MHMM_T*        an_mhmm,
   MHMM_T*        new_mhmm)
{
  int i_state;

  /* Copy the top-level data. */
  new_mhmm->type = an_mhmm->type;
  new_mhmm->log_odds = an_mhmm->log_odds;
  new_mhmm->num_motifs = an_mhmm->num_motifs;
  new_mhmm->num_states = an_mhmm->num_states;
  new_mhmm->num_spacers = an_mhmm->num_spacers;
  new_mhmm->spacer_states = an_mhmm->spacer_states;
  new_mhmm->alph_size = an_mhmm->alph_size;
  new_mhmm->ambigs = an_mhmm->ambigs;
  copy_string(&(new_mhmm->alphabet), an_mhmm->alphabet);
  new_mhmm->background = allocate_array(get_alph_size(ALL_SIZE));
  copy_array(an_mhmm->background, new_mhmm->background);
  copy_string(&(new_mhmm->description), an_mhmm->description);
  copy_string(&(new_mhmm->motif_file), an_mhmm->motif_file);
  copy_string(&(new_mhmm->sequence_file), an_mhmm->sequence_file);
  // FIXME: Copy hot states array.
  new_mhmm->num_hot_states = an_mhmm->num_hot_states;

  /* Copy each state. */
  for (i_state = 0; i_state < an_mhmm->num_states; i_state++) {
    copy_state(&(an_mhmm->states[i_state]),
	       &(new_mhmm->states[i_state]));
  }

  /* Copy the transition matrix. */
  copy_matrix(an_mhmm->trans, new_mhmm->trans);
}
  

/************************************************************************
 * See .h file for description.
 ************************************************************************/
void mix_mhmm
  (const float    mixing,   /* Percent of old HMM that is retained. */
   MHMM_T* const new_mhmm, /* New HMM to be mixed into the old one. */
   MHMM_T*       old_mhmm) /* Target HMM to be modified via mixing. */
{
  int i_state;

  /* Mix the emission distributions at each state. */
  for (i_state = 0; i_state < old_mhmm->num_states; i_state++) {
    mix_states(mixing, &(new_mhmm->states[i_state]), 
	       &(old_mhmm->states[i_state]));
  }

  /* Mix the transition matrices. */
  mix_log_matrices(mixing, new_mhmm->trans, old_mhmm->trans);

  /* Propagate the new transition matrix into the states. */
  compute_ins_and_outs(old_mhmm, TRUE);
}
  

/************************************************************************
 * Free the memory used by an mhmm state. 
 ************************************************************************/
static void free_state
  (MHMM_STATE_T *a_state)
{
  
  /* Don't bother with empty states. */
  if (a_state == NULL)
    return;
  
  myfree(a_state->itrans_out);
  free_array(a_state->trans_out);
  myfree(a_state->itrans_in);
  free_array(a_state->trans_in);
  free_array(a_state->emit);
  free_array(a_state->emit_odds);
  free_pssm(a_state->pssm);
  free_pssm(a_state->npssm);
  //free_array(a_state->pvalues);
}

/************************************************************************
 * See .h file for description.
 ************************************************************************/
void free_mhmm
  (MHMM_T* an_mhmm)
{
  int i_state;    /* Index of the current state. */

  /* Don't deallocate an empty struct. */
  if (an_mhmm == NULL)
    return;

  myfree(an_mhmm->alphabet);
  free_array(an_mhmm->background);
  myfree(an_mhmm->description);
  myfree(an_mhmm->motif_file);
  myfree(an_mhmm->sequence_file);

  for (i_state = 0; i_state < an_mhmm->num_states; i_state++)
    free_state(&(an_mhmm->states[i_state]));
  myfree(an_mhmm->states);
  free_matrix(an_mhmm->trans);
  myfree(an_mhmm->hot_states);
  myfree(an_mhmm);
}

/************************************************************************
 * Count the number of transitions entering or exiting a given state.
 *
 * RETURN: Total number of transitions to or from the given state.
 ************************************************************************/
#define TRANS_IN 0
#define TRANS_OUT 1
#ifdef oldway
static int count_trans
  (MATRIX_T* trans,      /* The transition matrix. */
   BOOLEAN_T log_form,   /* Is the transition matrix in log form? */
   int       num_states, /* Number of states in the (square) matrix. */
   int       state_num,  /* Index of the state we're interested in. */
   int       in_or_out)  /* Incoming or outgoing transitions? */
{
  int i_row;
  int i_col;
  int ntrans = 0;  /* The return value. */

  for (i_row = 0; i_row < num_states; i_row++) {
    for (i_col = 0; i_col < num_states; i_col++) {
      if (!is_zero(get_matrix_cell(i_row, i_col, trans), log_form)) {
	if ((in_or_out == TRANS_IN) && (i_col == state_num))
	  ntrans++;
	else if ((in_or_out == TRANS_OUT) && (i_row == state_num))
	  ntrans++;
      }
    }
  }
  return(ntrans);
} // count_trans
#endif

/************************************************************************
 * Compute the indices and values of transitions to or from a state.
 ************************************************************************/
void compute_ins_and_outs
  (MHMM_T*   the_hmm,
   BOOLEAN_T log_form) /* Is the transition matrix in log form? */
{
  int i_row, i_col;
  int n = the_hmm->num_states;
  MATRIX_T *trans = the_hmm->trans;

  //
  // Visit the transition matrix cells just once each
  // to update ntrans, itrans and trans arrays.
  // This is quadratic in n. 
  //
  for (i_row = 0; i_row < n; i_row++) {
    for (i_col = 0; i_col < n; i_col++) {
      double p;				// The transition probability.
      int old_n, new_n;			// Number of transitions.
      if (!is_zero((p = get_matrix_cell(i_row, i_col, trans)), log_form)) {
        MHMM_STATE_T * out_state = &(the_hmm->states[i_row]);
        MHMM_STATE_T * in_state = &(the_hmm->states[i_col]);
        // out
        old_n = out_state->ntrans_out; 
        new_n = ++out_state->ntrans_out;
        mm_resize(out_state->itrans_out, new_n, int);
        out_state->trans_out = resize_array(out_state->trans_out, new_n);
        out_state->itrans_out[old_n] = i_col;
        set_array_item(old_n, p, out_state->trans_out);
        // in
        old_n = in_state->ntrans_in; 
        new_n = ++in_state->ntrans_in;
        mm_resize(in_state->itrans_in, new_n, int);
        in_state->trans_in = resize_array(in_state->trans_in, new_n);
        in_state->itrans_in[old_n] = i_row;
        set_array_item(old_n, p, in_state->trans_in);
      }
    } // col
  } // row

} // compute_ins_and_outs

/*************************************************************************
 * Given an HMM state, choose the most frequently occuring character.
 *
 * RETURN: The most common character in the given state's log-odds 
 *         emission distribution.
 *************************************************************************/
char choose_consensus
  (BOOLEAN_T      log_space, /* Choose based on log-odds probs. */
   MHMM_STATE_T * a_state)   /* State to be analyzed. */
{
  PROB_T maximum;
  int    i_max;
  int    i_freq; 

  /* Don't bother with start, end and spacer states. */
  if ((a_state->type == SPACER_STATE) || (a_state->type == START_STATE)
      || (a_state->type == END_STATE)) {
    return('.');
  }

  /* Find out which character has maximum logodds frequency. */
  if (log_space) {
    maximum = LOG_ZERO;
  } else {
    maximum = 0.0;
  }
  i_max = 0;
  for (i_freq = 0; i_freq < get_alph_size(ALPH_SIZE); i_freq++) {
    if (log_space) {
      if (get_array_item(i_freq, a_state->emit_odds) > maximum) {
	i_max = i_freq;
	maximum = get_array_item(i_max, a_state->emit_odds);
      }
    } else {
      if (get_array_item(i_freq, a_state->emit) > maximum) {
	i_max = i_freq;
	maximum = get_array_item(i_max, a_state->emit);
      }
    }
  }

  /* Convert the index to a character. */
  return(get_alph_char(i_max));
}

/*************************************************************************
 * See .h file for description. 
 *************************************************************************/
#define SLOP 1E-5
BOOLEAN_T has_fims
  (MHMM_T* the_hmm)
{
  int num_states;      /* Number of states in the model. */
  int i_state;         /* Index of the current state. */
  MTYPE this_row_sum;  /* Sum of the transitions out of this state. */

  num_states = get_num_rows(the_hmm->trans);
  for (i_state = 0; i_state < num_states - 1; i_state++) {
    this_row_sum = array_total(get_matrix_row(i_state, the_hmm->trans));

    if (almost_equal(this_row_sum, 2.0, SLOP)) {
      return(TRUE);
    }
  }
  return(FALSE);
}

/************************************************************************
 * Count the number of motifs in an HMM.
 *************************************************************************/
int get_num_motifs
  (MHMM_T* the_hmm)
{
  int return_value = 0;
  int i_state;
  int num_states = the_hmm->num_states;

  for (i_state = 0; i_state < num_states; i_state++) {
    if ((the_hmm->states[i_state]).type == START_MOTIF_STATE) {
      return_value++;
    }
  }
  return(return_value);
}


