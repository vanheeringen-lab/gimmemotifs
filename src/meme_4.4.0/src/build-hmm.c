/********************************************************************
 * FILE: build-hmm.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 7-14-97
 * PROJECT: MHMM
 * COPYRIGHT: 1997-2008, WSN
 * DESCRIPTION: Build an HMM in memory.
 ********************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "matrix.h"
#include "utils.h"
#include "motif.h"
#include "alphabet.h"
#include "mhmm-state.h"
#include "order.h"
#include "build-hmm.h"
#include "array.h"

#define NON_MOTIF_INDEX -1
#define NON_MOTIF_ID "---"
#define NON_MOTIF_POSITION -1 
#define START_INDEX 0			// Index of start state in star model.
#define SPACER_INDEX 1			// Index of spacer state in star model.

#define SPACER_NUMSITES 10          /* Number of sites associated with
				       spacer emission distributions. */

/*************************************************************************
 * Given the expected length of a state with a self-loop, calculate
 * the associated transition probability.
 * 
 * Consider an HMM node where the transition probability in is x, the
 * transition probability out is 1-x, and the probability of a
 * self-loop is x.  What is the expected value of the number of visits
 * to the node?
 * 
 * By definition this is
 *
 *        mu = sum_0^infinity n(1-x)x^n
 *  
 * At first there are two possibilities: visit the node or skip it.
 * Skipping the node gives a path of length 0, while visiting it gives
 * a path of expected length 1 + the expected remaining path length,
 * call this nu.  So we have
 *
 *        mu = (1-x)0 + x(1+nu)
 *
 * Now because of the Markov property, whatever the path length so
 * far, if we reach this node again, then the expected path length
 * from it is simply mu.  So we have
 *
 *        mu = x(1+mu)
 *
 * Solving gives mu = x/(1-x) or x = mu/(1+mu). The latter is used in
 * this function.
 *
 * IN: length - The desired expected length.
 * OUT: ---
 * RETURN: The required transition probability to generate a sequence
 *         of the desired expected length. 
 *************************************************************************/
#define MIN_TRANS 0.1
static PROB_T self_trans
  (float length)
{
  /* Don't allow zero self transitions. */
  if (length == 0.0) {
    length = MIN_TRANS;
  }

  return (PROB_T)(length / (length + 1.0));
}

/*************************************************************************
 * Verify that the rows of a transition matrix sum to 1.0. 
 *************************************************************************/
#define SLOP 5E-6
BOOLEAN_T verify_trans_matrix
  (BOOLEAN_T log_form,    /* Is the transition matrix in log form? */
   int       num_states,  /* Number of states in the (square) matrix. */
   MATRIX_T* trans)       /* The matrix. */
{
  int    i_state;
  PROB_T total;

  for (i_state = 0; i_state < num_states - 1; i_state++) {

    /* Cf. Rabiner, formula (43b), p. 265. */
    if (log_form) {
      total = log_array_total(get_matrix_row(i_state, trans));
      if ((!almost_equal(total, 0.0, SLOP)) &&
	  (!almost_equal(total, 1.0, SLOP)) && // Allow for FIMS.
	  (!almost_equal(EXP2(total), 0.0, SLOP))) { 
	fprintf(stderr,
		"Warning: Row %d of transition matrix differs from 0.0 by %g.\n",
		i_state, EXP2(total));
	return(FALSE);
      }
    } else {
      total = array_total(get_matrix_row(i_state, trans));
      if ((!almost_equal(total, 1.0, SLOP)) &&
	  (!almost_equal(total, 2.0, SLOP)) && // Allow FIMs.
	  (!almost_equal(total, 0.0, SLOP))) { // Allow inaccessible motifs.
	fprintf(stderr,
		"Warning: Row %d of transition matrix differs from 1.0 by %g.\n",
		i_state, 1.0 - total);
	return(FALSE);
      }
    }

    /* All transitions from the end state must be zero. */
    if ((log_form) &&
	(get_matrix_cell(num_states - 1, i_state, trans) > LOG_SMALL)) {
      fprintf(stderr,
	      "Warning: Transition %d from end state is non-zero (%g).\n", 
	      i_state, get_matrix_cell(num_states - 1, i_state, trans));
      return(FALSE);
    } else if (!(log_form) &&
	       (!almost_equal(get_matrix_cell(num_states - 1, i_state, trans), 
			      0.0, SLOP))) {
      fprintf(stderr,
	      "Warning: Transition %d from end state is non-zero (%g).\n", 
	      i_state, get_matrix_cell(num_states - 1, i_state, trans));
      return(FALSE);
    }
  }
  return(TRUE);
}  

/*************************************************************************
 * Using information stored in the states of an HMM, fill in the HMM's
 * transition matrix.
 *************************************************************************/
static void build_transition_matrix
  (MHMM_T *the_hmm)
{
  int            i_state;    /* Indices into the matrix. */
  int            j_state; 
  MHMM_STATE_T * this_state; /* Pointer to the current state. */
  int            num_out;    /* No. of trans out of the current state. */
  int            i_out;      /* Index of outgoing transition. */

  /* First make sure the matrix is zeroed. */
  for (i_state = 0; i_state < the_hmm->num_states; i_state++) {
    for (j_state = 0; j_state < the_hmm->num_states; j_state++) {
      set_matrix_cell(i_state, j_state, 0.0, the_hmm->trans);
    }
  }

  /* Look at each state in the model. */
  for (i_state = 0; i_state < the_hmm->num_states; i_state++) {
    this_state = &(the_hmm->states[i_state]);
    
    /* Find out how many transitions out of this state there are. */
    num_out = this_state->ntrans_out;

    for (i_out = 0; i_out < num_out; i_out++) {
      /* Get the index of the state being transitioned to. */
      j_state = this_state->itrans_out[i_out];
      assert(j_state != 0);

      /* Fill in the matrix with the appropriate value. */
      set_matrix_cell(i_state, j_state, 
		      get_array_item(i_out, this_state->trans_out),
		      the_hmm->trans);
    }
  }
  assert(verify_trans_matrix(FALSE, the_hmm->num_states, the_hmm->trans));
}


/*************************************************************************
 * Convert spacer states in a given HMM to free-insertion modules.
 * All transitions out of spacers states are set to 1.0.
 *
 * Note that this function only affects the transitions stored in each
 * individual state.  It is assumed that build_transition_matrix will
 * be called subsequently.
 *************************************************************************/
static void convert_to_fims
  (MHMM_T *the_hmm)
{
  int            i_state;   /* Index of the current state. */
  MHMM_STATE_T * the_state; /* The current state. */
  int            i_trans;   /* Index of outgoing transition. */

  for (i_state = 0; i_state < the_hmm->num_states; i_state++) {
    the_state = &(the_hmm->states[i_state]);
  
    if (the_state->type == SPACER_STATE) {
      for (i_trans = 0; i_trans < the_state->ntrans_out; i_trans++)
	set_array_item(i_trans, 1.0, the_state->trans_out);
    }
  }
}

/*************************************************************************
 * Set up one state in a linear HMM, given the appropriate data.
 *************************************************************************/
static void build_linear_state
  (STATE_T  state_type,      /* Type of state (START, SPACER,...) */
   int      i_state,         /* The state index. */
   int      expected_length, /* For spacers, expected length of
				output. */
   ARRAY_T* freqs,           /* Emission probability distrib. */
   double   num_sites,       // Number of sites for this emission.
   int      alph_size,       /* Size of emission distribution. */
   int      i_motif,         /* Index of motif this state is in. */
   char*    motif_id,        // MEME ID of this motif.
   int      i_position,      /* Position of this state within a
				motif or spacer. */
   int	    w_motif,	     // Width of motif.
   MOTIF_T*  motifs,         // Motifs.
   MHMM_STATE_T * a_state)   /* State to be filled in (pre-allocated). */
{
  if (verbosity >= NORMAL_VERBOSE) {
    switch (state_type) {
    case START_STATE :
      fprintf(stderr, "Building HMM: 0 ");
      break;
    case SPACER_STATE :
    case END_MOTIF_STATE :
      fprintf(stderr, "%d ", i_state);
      break;
    case START_MOTIF_STATE :
    case MID_MOTIF_STATE :
      fprintf(stderr, "%d-", i_state);
      break;
    case END_STATE :
      fprintf(stderr, "%d\n", i_state);
      break;
    case INVALID_STATE :
      die("Invalid state.\n");
    }
  }

  /* Record what type of state this is. */
  a_state->type = state_type;

  // Record the motif width if this is a motif.
  if (state_type == START_MOTIF_STATE ||
      state_type == MID_MOTIF_STATE ||
      state_type == END_MOTIF_STATE) {
    a_state->w_motif = w_motif;
  } else {
    a_state->w_motif = 1;
  }
  
  /* Set up the emission distribution and a few other tidbits. */
  a_state->emit = allocate_array(alph_size);
  a_state->emit_odds = allocate_array(alph_size);
  if (state_type == START_STATE || state_type == END_STATE) { 
    /* Start and end don't have emissions. */
    int i_alph;
    for (i_alph = 0; i_alph < get_alph_size(ALL_SIZE); i_alph++) {
      set_array_item(i_alph, 1.0, a_state->emit);
    }
  }
  else {
    copy_array(freqs, a_state->emit);
  }
  a_state->num_sites = num_sites;

  /* Record the motif index and ID. */
  a_state->i_motif = i_motif;
  strcpy(a_state->motif_id, motif_id);
  if ((state_type == START_STATE) ||
      (state_type == END_STATE) ||
      (state_type == SPACER_STATE)) {
    a_state->id_char = NON_MOTIF_ID_CHAR;
  } else {
    a_state->id_char = get_motif_id_char(i_position, &(motifs[i_motif]));
  }
  a_state->i_position = i_position;

  /* First set up the transitions into this state. */
  switch (state_type) {
  case START_STATE :
    a_state->ntrans_in = 0;
    a_state->itrans_in = NULL;
    a_state->trans_in = NULL;
    break;
  case START_MOTIF_STATE :
  case END_STATE :
    a_state->ntrans_in = 2;
    a_state->itrans_in = (int *)mm_malloc(sizeof(int) * 2);
    a_state->itrans_in[0] = i_state - 2;
    a_state->itrans_in[1] = i_state - 1;
    a_state->trans_in = allocate_array(2);
    set_array_item(0, 1.0 - self_trans(expected_length), a_state->trans_in);
    set_array_item(1, 1.0 - self_trans(expected_length), a_state->trans_in);
    break;
  case MID_MOTIF_STATE :
  case END_MOTIF_STATE :
    a_state->ntrans_in = 1;
    a_state->itrans_in = (int *)mm_malloc(sizeof(int));
    a_state->itrans_in[0] = i_state - 1;
    a_state->trans_in = allocate_array(1);
    set_array_item(0, 1.0, a_state->trans_in);
    break;
  case SPACER_STATE :
    a_state->ntrans_in = 2;
    a_state->itrans_in = (int *)mm_malloc(sizeof(int) * 2);
    a_state->itrans_in[0] = i_state - 1;
    a_state->itrans_in[1] = i_state;
    a_state->trans_in = allocate_array(2);
    set_array_item(0, 1.0 - self_trans(expected_length), a_state->trans_in);
    set_array_item(1, self_trans(expected_length), a_state->trans_in);
    break;
  default:
    die("Invalid state type.\n");
  }

  /* Then set up the transitions out of this state. */
  switch (state_type) {
  case START_STATE :
  case END_MOTIF_STATE :
    a_state->ntrans_out = 2;
    a_state->itrans_out = (int *)mm_malloc(sizeof(int) * 2);
    a_state->itrans_out[0] = i_state + 1;
    a_state->itrans_out[1] = i_state + 2;
    a_state->trans_out = allocate_array(2);
    set_array_item(0, self_trans(expected_length), a_state->trans_out);
    set_array_item(1, 1.0 - self_trans(expected_length), a_state->trans_out);
    break;
  case START_MOTIF_STATE :
  case MID_MOTIF_STATE :
    a_state->ntrans_out = 1;
    a_state->itrans_out = (int *)mm_malloc(sizeof(int));
    a_state->itrans_out[0] = i_state + 1;
    a_state->trans_out = allocate_array(1);
    set_array_item(0, 1.0, a_state->trans_out);
    break;
  case SPACER_STATE :
    a_state->ntrans_out = 2;
    a_state->itrans_out = (int *)mm_malloc(sizeof(int) * 2);
    a_state->itrans_out[0] = i_state;
    a_state->itrans_out[1] = i_state + 1;
    a_state->trans_out = allocate_array(2);
    set_array_item(0, self_trans(expected_length), a_state->trans_out);
    set_array_item(1, 1.0 - self_trans(expected_length), a_state->trans_out);
    break;
  case END_STATE :
    a_state->ntrans_out = 0;
    a_state->itrans_out = NULL;
    a_state->trans_out = NULL;
    break;
  default:
    die("Invalid state type.\n");
  }
}

/*************************************************************************
 * Build a linear HMM.
 *************************************************************************/
void build_linear_hmm
  (ARRAY_T*  background,
   ORDER_T*  order_spacing,
   int       spacer_states, 
   MOTIF_T*  motifs,
   int       nmotifs,
   BOOLEAN_T fim,
   MHMM_T**  the_hmm)
{
  int       model_length; /* Total number of states in the model. */
  int       i_state;      /* Index of the current state. */
  int       i_motif;      /* Index of the current motif. */
  int       i_position;   /* Index within the current motif. */

  /* Calculate the total length of the model. */
  i_motif = 0;
  for (model_length = 2, i_motif = 0; i_motif < nmotifs; i_motif++)
    model_length += (motifs[i_motif]).length;
  model_length += (nmotifs + 1) * spacer_states;

  /* Allocate the model. */
  allocate_mhmm(model_length, the_hmm);

  /* Record that this is a linear model. */
  (*the_hmm)->type = LINEAR_HMM;

  /* Record the number of motifs in the model. */
  (*the_hmm)->num_motifs = nmotifs;

  /* Record the number of states in the model. */
  (*the_hmm)->num_states = model_length;
  (*the_hmm)->num_spacers = nmotifs + 1;
  (*the_hmm)->spacer_states = spacer_states;

  /* Put the alphabet into the model. */
  (*the_hmm)->alph_size = get_alph_size(ALPH_SIZE);
  (*the_hmm)->ambigs = get_alph_size(AMBIG_SIZE);
  strcpy((*the_hmm)->alphabet, get_alphabet(FALSE));

  // Put the background distribution into the model.
  copy_array(background, (*the_hmm)->background);

  /* Begin the model with a non-emitting state. */
  i_state = 0;
  build_linear_state(START_STATE,
		     i_state,
		     get_spacer_length(order_spacing, 0),
		     NULL, // Emissions.
		     0, // Number of sites.
		     get_alph_size(ALL_SIZE),
		     NON_MOTIF_INDEX,
		     NON_MOTIF_ID,
		     NON_MOTIF_POSITION,
		     0,
		     motifs,
		     &((*the_hmm)->states[i_state]));
  i_state++;

  /* Build the first spacer. */
  for (i_position = 0; i_position < spacer_states; i_position++) {
    build_linear_state(SPACER_STATE,
		       i_state, 
		       get_spacer_length(order_spacing, 0),
		       background, 
		       SPACER_NUMSITES,
		       get_alph_size(ALL_SIZE),
		       NON_MOTIF_INDEX,
		       NON_MOTIF_ID,
		       i_position,
		       0,
		       motifs,
		       &((*the_hmm)->states[i_state]));
    i_state++;
  }

  /* Build each motif and subsequent spacer. */
  for (i_motif = 0; i_motif < nmotifs; i_motif++) {
    MOTIF_T *this_motif = &(motifs[i_motif]);

    /* Build the motif. */
    build_linear_state(START_MOTIF_STATE,
		       i_state,
		       get_spacer_length(order_spacing, i_motif),
		       get_matrix_row(0, this_motif->freqs),
		       this_motif->num_sites,
		       get_alph_size(ALL_SIZE), 
		       i_motif,
		       get_motif_id(this_motif),
                       0, // First position in the motif.
		       this_motif->length,
		       motifs,
		       &((*the_hmm)->states[i_state]));
    i_state++;
    for (i_position = 1; i_position < this_motif->length - 1; i_position++) {
      build_linear_state(MID_MOTIF_STATE, 
			 i_state, 
			 0, // Expected spacer length.
			 get_matrix_row(i_position, this_motif->freqs),
			 this_motif->num_sites,
			 get_alph_size(ALL_SIZE),
			 i_motif,
			 get_motif_id(this_motif),
			 i_position,
		         this_motif->length,
			 motifs,
			 &((*the_hmm)->states[i_state]));
      i_state++;
    }
    build_linear_state(END_MOTIF_STATE, 
		       i_state, 
		       get_spacer_length(order_spacing, i_motif+1),
		       get_matrix_row(i_position, this_motif->freqs),
		       this_motif->num_sites,
		       get_alph_size(ALL_SIZE), 
		       i_motif,
		       this_motif->id,
		       this_motif->length - 1,
		       this_motif->length,
		       motifs,
		       &((*the_hmm)->states[i_state]));
    i_state++;

    /* Build the following spacer. */
    for (i_position = 0; i_position < spacer_states; i_position++) {
      build_linear_state(SPACER_STATE, 
			 i_state, 
			 get_spacer_length(order_spacing, i_motif+1),
			 background,
			 SPACER_NUMSITES,
			 get_alph_size(ALL_SIZE), 
			 NON_MOTIF_INDEX, 
			 NON_MOTIF_ID,
			 i_position,
			 0,
			 motifs,
			 &((*the_hmm)->states[i_state]));
      i_state++;
    }
  }

  /* Finish up the model with a non-emitting end state. */
  build_linear_state(END_STATE, 
		     i_state, 
		     get_spacer_length(order_spacing, i_motif),
		     NULL, // Emissions.
		     0, // Number of sites.
			   get_alph_size(ALL_SIZE), // Size of emissions.
		     NON_MOTIF_INDEX,
		     NON_MOTIF_ID,
		     NON_MOTIF_POSITION,
		     0,
		     motifs,
		     &((*the_hmm)->states[i_state]));
  i_state++;

  /* Convert spacers to FIMs if requested. */
  if (fim) {
    convert_to_fims(*the_hmm);
  }

  /* Fill in the transition matrix. */
  build_transition_matrix(*the_hmm);
}

/*************************************************************************
 * Find the index, in a complete MHMM, of a spacer state, given the
 * indices of the motifs between which the spacer appears.
 *
 * In a complete MHMM, spacers appear at the beginning of the model,
 * starting at state 1 (state 0 is the start state).  The first n
 * spacers are those between the start state and each of the n motifs.
 * Next, there are n+1 spacers, corresponding to those between motif 1
 * motif and every motif plus the end state.
 *
 * The spacers can thus be thought of as indexed according to a
 * matrix, where the row is the motif from which the spacer comes, and
 * the column is the motif where it's going.
 *
 * Here's an example of such a matrix for the case of 3 motifs:
 *
 *    S  1  2  3  E
 *  S x  1  2  3  x
 *  1 x  4  5  6  7
 *  2 x  8  9 10 11
 *  3 x 12 13 14 15
 *  E x  x  x  x  x
 *
 * More generally, the spacers are indexed as follows:
 * 
 *  Index  From  To
 * ------------------
 *    1     0     1    Spacers from the start state to each motif.
 *    2     0     2
 *   ...   ...   ...
 *    n     0     n
 * ------------------
 *   n+1    1     1    Spacer from the first motif to every other motif,
 *   n+2    1     2    including the end state.
 *   ...   ...   ...
 *  2n+1    1    n+1
 * ------------------
 *  2n+2    2     1    Spacers from the second motif.
 *  2n+3    2     2
 *   ...   ...   ...
 *  3n+2    2    n+1
 * ------------------
 *    .     .     .
 *    .     .     .
 *    .     .     .     
 * ------------------
 *  n*n+n   n     1    Spacers from the last motif.
 * n*n+n+1  n     2
 *   ...   ...   ...
 * n*n+2n   n    n+1
 * ------------------
 *
 * This picture is further complicated by the fact that every spacer
 * can contain multiple states.  Thus, every 'n' above should be
 * replaced by a 'pn', where p is the number of states per spacer.
 *
 * RETURN: Index of the beginning or ending of the requested spacer.
 *************************************************************************/
static int spacer_index
  (const int from_motif,      /* Index of motif preceding this spacer
				 (0=start,n+1=end). */
   const int to_motif,        /* Index of motif following this spacer
			         (0=start,n+1=end). */
   const BOOLEAN_T end_state, /* Is this an end state? */
   const int nmotifs,         /* Total number of motifs in the model. */
   const int spacer_states)   /* Number of states per spacer in the model. */
{
  int return_value;

  /* Make sure we're not trying something illegal. */
  if (to_motif == 0) {
    die("Requested spacer leading to the start state.\n");
  } else if (from_motif == nmotifs + 1) {
    die("Requested spacer leading from the end state.\n");
  }

  /* The spacer indices start at 1. */
  return_value = 1;

  /* Add in the lengths of complete rows. */
  return_value += from_motif * spacer_states * (nmotifs+1);

  /* Add the final row. */
  return_value += (to_motif - 1) * spacer_states;

  /* The first row of the matrix is missing one value, so subtract it. */
  if (from_motif > 0)
    return_value -= spacer_states;

  /* If the end of the spacer was requested, add. */
  if (end_state)
    return_value += spacer_states - 1;

  return(return_value);
}

/*************************************************************************
 * Find the index of the starting or ending state of a given motif in
 * a given HMM.
 *************************************************************************/
static int motif_index
  (const int       motif_num,
   const BOOLEAN_T start_or_end,
   const int       num_spacers,
   const int       spacer_states,
   const MOTIF_T*  motifs)
{
  int i_motif;   
  int return_value;

  /* Skip the spacer states. */
  return_value = (num_spacers * spacer_states) + 1;

  /* Add the lengths of the preceding motifs. */
  for (i_motif = 0; i_motif < motif_num - 1; i_motif++)
    return_value += motifs[i_motif].length;
  
  /* If we're looking for the end of this motif, add its length as well. */
  if (start_or_end)
    return_value += motifs[i_motif].length - 1;

  /* fprintf(stderr, "Motif %d -> %d\n", motif_num, return_value); */

  return(return_value);
}	    



/*************************************************************************
 * Set up one state in a complete HMM, given the appropriate data.
 *************************************************************************/
static void build_complete_state
  (int       state_type,      /* Type of state (START, SPACER,..) */
   STATE_T   i_state,         /* State index. */
   int       expected_length, /* For spacers, the expected length
				 of output. */
   ARRAY_T*  freqs,           /* Emission probability distrib. */
   double    num_sites,       // Number of sites for this emission.
   int       alph_size,       /* Size of emission distribution. */
   int       i_motif,         /* Index of motif this state is in. */
   char*     motif_id,        // MEME-assigned ID of the motif.
   int       i_position,      /* Position of this state within motif */
   int	     w_motif,	      // Width of motif.
   int       nmotifs,         /* Total number of motifs. */
   int       prev_motif,      /* Index of previous motif. */
   int       next_motif,      /* Index of next motif. */
   MATRIX_T* transp_freq,     /* Transition freq matrix. */
   int       spacer_states,   /* Number of HMM states per spacer. */
   int       num_spacers,     /* Total number of spacers in HMM. */
   MOTIF_T*  motifs,          /* Motifs. */
   MHMM_STATE_T* a_state)     /* State to be filled in (pre-allocated). */
{
  int j_motif;    /* Index of the current motif. */

  /* Tell the user what's up. */
  if (verbosity >= NORMAL_VERBOSE) {
    switch (state_type) {
    case START_STATE :
      fprintf(stderr, "Building HMM: (0) ");
      break;
    case SPACER_STATE :
      fprintf(stderr, "%d ", i_state);
      break;
    case END_MOTIF_STATE :
      fprintf(stderr, "%d | ", i_state);
      break;
    case START_MOTIF_STATE :
    case MID_MOTIF_STATE :
      fprintf(stderr, "%d-", i_state);
      break;
    case END_STATE :
      fprintf(stderr, "(%d)\n", i_state);
      break;
    }
  }

  /* Record what type of state this is. */
  a_state->type = state_type;

  // Record the motif width if this is a motif.
  if (state_type == START_MOTIF_STATE ||
      state_type == MID_MOTIF_STATE ||
      state_type == END_MOTIF_STATE) {
    a_state->w_motif = w_motif;
  } else {
    a_state->w_motif = 0;
  }
  

  /* Set up the emission distribution and a few other tidbits. */
  if (freqs != NULL) { /* Start and end states have no emissions. */
    a_state->emit = allocate_array(alph_size);
    copy_array(freqs, a_state->emit);
  }
  a_state->num_sites = num_sites;
  a_state->i_motif = i_motif;
  strcpy(a_state->motif_id, motif_id);
  a_state->i_position = i_position;

  // Record the motif ID character at this position.
  if ((state_type == START_STATE) ||
      (state_type == END_STATE) ||
      (state_type == SPACER_STATE)) {
    a_state->id_char = NON_MOTIF_ID_CHAR;
  } else {
    a_state->id_char = get_motif_id_char(i_position, &(motifs[i_motif]));
  }
  assert(a_state->id_char != '\0');

  /* First set up the transitions into this state. */
  switch (state_type) {
  case START_STATE :
    a_state->ntrans_in = 0;
    a_state->itrans_in = NULL;
    a_state->trans_in = NULL;
    break;
  case START_MOTIF_STATE :
    /* Transitions come from any motif or from the start state. */
    a_state->ntrans_in = nmotifs + 1;
    a_state->itrans_in = (int *)mm_malloc(sizeof(int) * (nmotifs + 1));
    a_state->trans_in = allocate_array(nmotifs + 1);
    for (j_motif = 0; j_motif < nmotifs + 1; j_motif++) {
      a_state->itrans_in[j_motif]
	= spacer_index(j_motif, i_motif + 1, TRUE, nmotifs, spacer_states);
      set_array_item(j_motif, 
		     get_matrix_cell(j_motif, i_motif + 1, transp_freq), 
		     a_state->trans_in);
    }
    break;
  case END_STATE :
    /* Transitions come from any motif. */
    a_state->ntrans_in = nmotifs;
    a_state->itrans_in = (int *)mm_malloc(sizeof(int) * nmotifs);
    a_state->trans_in = allocate_array(nmotifs);
    for (j_motif = 0; j_motif < nmotifs; j_motif++) {
      a_state->itrans_in[j_motif] = spacer_index(j_motif + 1,
						 nmotifs + 1, TRUE,
						 nmotifs, spacer_states);
      set_array_item(j_motif, 
		     get_matrix_cell(j_motif + 1, nmotifs + 1, transp_freq), 
		     a_state->trans_in);
    }
    break;
  case MID_MOTIF_STATE :
  case END_MOTIF_STATE :
    a_state->ntrans_in = 1;
    a_state->itrans_in = (int *)mm_malloc(sizeof(int));
    a_state->itrans_in[0] = i_state - 1;
    a_state->trans_in = allocate_array(1);
    set_array_item(0, 1.0, a_state->trans_in);
    break;
  case SPACER_STATE :
    a_state->ntrans_in = 2;
    a_state->itrans_in = (int *)mm_malloc(sizeof(int) * 2);
    a_state->trans_in = allocate_array(2);
    /* For multi-state spacers, incoming transition from previous state. */
    if (i_position != 0)
      a_state->itrans_in[0] = i_state - 1;
    else 
      a_state->itrans_in[0] = motif_index(prev_motif, TRUE, num_spacers,
					  spacer_states, motifs);
    /* The other transition is a self-transition. */
    a_state->itrans_in[1] = i_state;
    set_array_item(0, 1.0 - self_trans(expected_length / spacer_states),
		   a_state->trans_in);
    set_array_item(1, self_trans(expected_length / spacer_states),
		   a_state->trans_in);
    break;
  }

  /* Then set up the transitions out of this state. */
  switch (state_type) {
  case START_STATE :
    /* Transitions go to each motif. */
    a_state->ntrans_out = nmotifs;
    a_state->itrans_out = (int *)mm_malloc(sizeof(int) * nmotifs);
    a_state->trans_out = allocate_array(nmotifs);
    for (j_motif = 0; j_motif < nmotifs; j_motif++) {
      a_state->itrans_out[j_motif] = spacer_index(0, j_motif + 1, FALSE,
						  nmotifs, spacer_states);
      set_array_item(j_motif,
		     get_matrix_cell(0, j_motif + 1, transp_freq),
		     a_state->trans_out);
    }
    break;
  case END_MOTIF_STATE :
    /* Can go to any other motif or to the end state. */
    a_state->ntrans_out = nmotifs + 1;
    a_state->itrans_out = (int *)mm_malloc(sizeof(int) * (nmotifs + 1));
    a_state->trans_out = allocate_array(nmotifs + 1);
    for (j_motif = 0; j_motif < nmotifs + 1; j_motif++) {
      a_state->itrans_out[j_motif] = spacer_index(i_motif + 1,
						  j_motif + 1, FALSE,
						  nmotifs, spacer_states);
      set_array_item(j_motif,
		     get_matrix_cell(i_motif + 1, j_motif + 1, transp_freq),
		     a_state->trans_out);
    }
    break;
  case START_MOTIF_STATE :
  case MID_MOTIF_STATE :
    a_state->ntrans_out = 1;
    a_state->itrans_out = (int *)mm_malloc(sizeof(int));
    a_state->itrans_out[0] = i_state + 1;
    a_state->trans_out = allocate_array(1);
    set_array_item(0, 1.0, a_state->trans_out);
    break;
  case SPACER_STATE :
    a_state->ntrans_out = 2;
    a_state->itrans_out = (int *)mm_malloc(sizeof(int) * 2);
    a_state->trans_out = allocate_array(2);
    /* The first transition is a self-transition. */
    a_state->itrans_out[0] = i_state;
    /* For multi-state spacers, outgoing transition to next state. */
    if (i_position < spacer_states - 1)
      a_state->itrans_out[1] = i_state + 1;
    else 
      a_state->itrans_out[1] = motif_index(next_motif, FALSE, num_spacers,
					   spacer_states, motifs);
    set_array_item(0, self_trans(expected_length), a_state->trans_out);
    set_array_item(1, 1.0 - self_trans(expected_length), a_state->trans_out);
    break;
  case END_STATE :
    a_state->ntrans_out = 0;
    a_state->itrans_out = NULL;
    a_state->trans_out = NULL;
    break;
  }
}


/*************************************************************************
 * Build a completely connected HMM.
 *************************************************************************/
void build_complete_hmm
  (ARRAY_T*  background,
   int       spacer_states, 
   MOTIF_T*  motifs,
   int       nmotifs,
   MATRIX_T* transp_freq,
   MATRIX_T* spacer_ave,
   BOOLEAN_T fim,
   MHMM_T**  the_hmm)
{
  int       motif_states; /* Total length of the motifs. */
  int       num_spacers;  /* Total number of spacer states. */
  int       num_states;   /* Total number of states in the model. */
  int       i_motif;      /* Index of the current "from" motif. */
  int       j_motif;      /* Index of the current "to" motif. */
  int       i_position;   /* Index within the current motif or spacer. */
  int       i_state = 0;  /* Index of the current state. */

  /* Count the width of the motifs. */
  for (motif_states = 0, i_motif = 0; i_motif < nmotifs; i_motif++)
    motif_states += (motifs[i_motif]).length;
  /* Count the spacer states adjacent to begin and end. */
  num_spacers = nmotifs * 2;
  /* Add the spacer states between motifs. */
  num_spacers += nmotifs * nmotifs;
  /* Total states = motifs + spacer_states + begin/end */
  num_states = motif_states + (num_spacers * spacer_states) + 2;
  /* fprintf(stderr, "motif_states=%d num_spacers=%d num_states=%d\n",
	  motif_states, num_spacers, num_states); */

  /* Allocate the model. */
  allocate_mhmm(num_states, the_hmm);

  /* Record that this is a completely connected model. */
  (*the_hmm)->type = COMPLETE_HMM;

  /* Record the number of motifs in the model. */
  (*the_hmm)->num_motifs = nmotifs;

  /* Record the number of states in the model. */
  (*the_hmm)->num_states = num_states;
  (*the_hmm)->num_spacers = ((nmotifs + 1) * (nmotifs + 1)) - 1;
  (*the_hmm)->spacer_states = spacer_states;

  /* Put the alphabet into the model. */
  (*the_hmm)->alph_size = get_alph_size(ALPH_SIZE);
  (*the_hmm)->ambigs = get_alph_size(AMBIG_SIZE);
  strcpy((*the_hmm)->alphabet, get_alphabet(FALSE));

  // Put the background distribution into the model.
  copy_array(background, (*the_hmm)->background);

  /* Build the begin state. */
  build_complete_state(START_STATE,
		       i_state,
		       0, // expected length
		       NULL, // Emissions.
		       0, // Number of sites.
		       0, // alphabet size
		       NON_MOTIF_INDEX,
		       NON_MOTIF_ID,
		       NON_MOTIF_POSITION,
		       0,
		       nmotifs,
		       0, // previous motif
		       0, // next motif
		       transp_freq,
		       spacer_states,
		       num_spacers,
		       motifs,
		       &((*the_hmm)->states[i_state]));
  i_state++;

  /* Build the spacer states. No transitions from the end state. */
  for (i_motif = 0; i_motif <= nmotifs; i_motif++) {
    /* No transitions to the start state. */
    for (j_motif = 1; j_motif <= nmotifs+1; j_motif++) {
      /* No transitions from start to end. */
      if ((i_motif == 0) && (j_motif == nmotifs+1))
	continue;
      /* Allow multi-state spacers. */
      for (i_position = 0; i_position < spacer_states; i_position++) {
	build_complete_state(SPACER_STATE, 
			     i_state, 
			     get_matrix_cell(i_motif, j_motif, spacer_ave),
			     background,
			     SPACER_NUMSITES,
			     get_alph_size(ALL_SIZE),
			     NON_MOTIF_INDEX,
			     NON_MOTIF_ID,
			     i_position,
			     0, 
			     nmotifs,
			     i_motif,
			     j_motif,
			     transp_freq,
			     spacer_states,
			     num_spacers,
			     motifs,
			     &((*the_hmm)->states[i_state]));
	i_state++;
      }
    }
  }

  /* Build the motif states. */
  for (i_motif = 0; i_motif < nmotifs; i_motif++) {
    MOTIF_T *this_motif = &(motifs[i_motif]);
    i_position = 0;
    build_complete_state(START_MOTIF_STATE, 
			 i_state,
			 0, // Expected spacer length.
			 get_matrix_row(i_position, this_motif->freqs),
			 this_motif->num_sites,
			 get_alph_size(ALL_SIZE),
			 i_motif,
			 get_motif_id(this_motif),
			 i_position,
                         this_motif->length,
			 nmotifs,
			 0, // Previous motif index.
			 0, // Next motif index.
			 transp_freq,
			 spacer_states,
			 num_spacers,
			 motifs,
			 &((*the_hmm)->states[i_state]));
    i_state++;
    for (i_position = 1; i_position < this_motif->length - 1; i_position++) {
      build_complete_state(MID_MOTIF_STATE, 
			   i_state,
			   0, // Expected spacer length. 
			   get_matrix_row(i_position, this_motif->freqs),
			   this_motif->num_sites,
			   get_alph_size(ALL_SIZE),
			   i_motif,
			   get_motif_id(this_motif),
			   i_position, 
		           this_motif->length,
			   nmotifs,
			   0, // Previous motif index.
			   0, // Next motif index.
			   transp_freq,
			   spacer_states,
			   num_spacers,
			   motifs,
			   &((*the_hmm)->states[i_state]));
      i_state++;
    }
    build_complete_state(END_MOTIF_STATE, 
			 i_state,
			 0, // Expected spacer length.
			 get_matrix_row(i_position, this_motif->freqs),
			 this_motif->num_sites,
			 get_alph_size(ALL_SIZE),
			 i_motif,
			 get_motif_id(this_motif),
			 i_position,
		         this_motif->length,
			 nmotifs,
			 0, // Previous motif index.
			 0, // Next motif index.
			 transp_freq,
			 spacer_states,
			 num_spacers,
			 motifs,
			 &((*the_hmm)->states[i_state]));
    i_state++;
  }

  /* Build the end state. */
  build_complete_state(END_STATE, 
		       i_state,
		       0, // Expected spacer length.
		       NULL, // Emissions
		       0, // Number of sites.
		       0, // Size of emissions.
		       NON_MOTIF_INDEX,
		       NON_MOTIF_ID,
		       NON_MOTIF_POSITION,
		       0,
		       nmotifs,
		       0, // Previous motif index.
		       0, // Next motif index.
		       transp_freq,
		       spacer_states,
		       num_spacers,
		       motifs,
		       &((*the_hmm)->states[i_state]));
  i_state++;

  /* Convert spacers to FIMs if requested. */
  if (fim) {
    convert_to_fims(*the_hmm);
  }

  /* Fill in the transition matrix. */
  build_transition_matrix(*the_hmm);
}

/*************************************************************************
 * Set up one state in a star HMM, given the appropriate data.
 *************************************************************************/
static void build_star_state
  (int       state_type,      // Type of state (START, SPACER,..)
   STATE_T   i_state,         // State index.
   int       expected_length, /* For spacers, the expected length
				 of output. */
   ARRAY_T*  freqs,           // Emission probability distrib.
   double    num_sites,       // Number of sites for this emission.
   int       alph_size,       // Size of emission distribution.
   int       i_motif,         // Index of motif this state is in.
   char*     motif_id,        // MEME-assigned ID of the motif.
   int       i_position,      // Position of this state within motif
   int	     w_motif,	      // Width of motif.
   int       nmotifs,         // Total number of motifs.
   int       spacer_states,   // Number of HMM states per spacer.
   MOTIF_T*  motifs,          // Motifs.
   MHMM_STATE_T* a_state)     // State to be filled in (pre-allocated).
{
  int j_motif;    		// Index of the current motif.
  int num_spacers = 1;     	// Total number of spacers in HMM.
  double in_p;			// Probability of transition into a state 

  // Tell the user what's up.
  if (verbosity >= NORMAL_VERBOSE) {
    switch (state_type) {
    case START_STATE :
      fprintf(stderr, "Building HMM: (0) ");
      break;
    case SPACER_STATE :
      fprintf(stderr, "%d ", i_state);
      break;
    case END_MOTIF_STATE :
      fprintf(stderr, "%d | ", i_state);
      break;
    case START_MOTIF_STATE :
    case MID_MOTIF_STATE :
      fprintf(stderr, "%d-", i_state);
      break;
    case END_STATE :
      fprintf(stderr, "(%d)\n", i_state);
      break;
    }
  }

  // Record what type of state this is.
  a_state->type = state_type;

  // Record the motif width if this is a motif.
  if (state_type == START_MOTIF_STATE ||
      state_type == MID_MOTIF_STATE ||
      state_type == END_MOTIF_STATE) {
    a_state->w_motif = w_motif;
  } else {
    a_state->w_motif = 0;
  }

  // Set up the emission distribution and a few other tidbits.
  if (freqs != NULL) { // Start and end states have no emissions.
    a_state->emit = allocate_array(alph_size);
    copy_array(freqs, a_state->emit);
  }
  a_state->num_sites = num_sites;
  a_state->i_motif = i_motif;
  strcpy(a_state->motif_id, motif_id);
  a_state->i_position = i_position;

  // Record the motif ID character at this position.
  if ((state_type == START_STATE) ||
      (state_type == END_STATE) ||
      (state_type == SPACER_STATE)) {
    a_state->id_char = NON_MOTIF_ID_CHAR;
  } else {
    a_state->id_char = get_motif_id_char(i_position, &(motifs[i_motif]));
  }
  assert(a_state->id_char != '\0');

  // First set up the transitions into this state.
  switch (state_type) {
  case START_STATE :
    a_state->ntrans_in = 0;
    a_state->itrans_in = NULL;
    a_state->trans_in = NULL;
    break;
  case END_STATE :
  case START_MOTIF_STATE :
    // Transitions come from spacer state.
    a_state->ntrans_in = 1;
    a_state->itrans_in = (int *)mm_malloc(sizeof(int));
    a_state->trans_in = allocate_array(1);
    a_state->itrans_in[0] = SPACER_INDEX;
    // Distribute non-self loop probability evenly among motifs and end state.
    in_p = (1 - self_trans(expected_length / spacer_states))/(nmotifs+1);
    set_array_item(0, in_p, a_state->trans_in);
    break;
  case MID_MOTIF_STATE :
  case END_MOTIF_STATE :
    // Transitions come from previous state.
    a_state->ntrans_in = 1;
    a_state->itrans_in = (int *)mm_malloc(sizeof(int));
    a_state->itrans_in[0] = i_state - 1;
    a_state->trans_in = allocate_array(1);
    set_array_item(0, 1.0, a_state->trans_in);
    break;
  case SPACER_STATE :
    // Transitions come from start and each motif except for internal
    // multi-states.
    a_state->ntrans_in = (i_position != 0) ? 2 : nmotifs + 2;
    a_state->itrans_in = (int *)mm_malloc(sizeof(int) * a_state->ntrans_in);
    a_state->trans_in = allocate_array(a_state->ntrans_in);

    // First transition is a self-transition. 
    a_state->itrans_in[0] = i_state;
    set_array_item(0, self_trans(expected_length / spacer_states), 
		   a_state->trans_in);

    // Next the transitions from all the motifs (or the previous spacer).
    if (i_position != 0) {
      a_state->itrans_in[1] = i_state - 1;
      set_array_item(1, 1.0 - self_trans(expected_length / spacer_states),
		     a_state->trans_in);
    } else {
      a_state->itrans_in[1] = START_INDEX;	// From start state.
      // From each motif.
      for (j_motif = 0; j_motif < nmotifs; j_motif++) {	
        a_state->itrans_in[j_motif+2] = 
          motif_index(j_motif+1, TRUE, num_spacers, spacer_states, motifs);
	set_array_item(j_motif+2, 1.0, a_state->trans_in);
      }
    }
    break;
  }

  // Then set up the transitions out of this state.
  switch (state_type) {
  case START_STATE :
  case END_MOTIF_STATE :
    // Transition goes to spacer. 
    a_state->ntrans_out = 1;
    a_state->itrans_out = (int *)mm_malloc(sizeof(int));
    a_state->trans_out = allocate_array(1);
    a_state->itrans_out[0] = SPACER_INDEX;
    set_array_item(0, 1.0, a_state->trans_out);
    break;
  case START_MOTIF_STATE :
  case MID_MOTIF_STATE :
    a_state->ntrans_out = 1;
    a_state->itrans_out = (int *)mm_malloc(sizeof(int));
    a_state->itrans_out[0] = i_state + 1;
    a_state->trans_out = allocate_array(1);
    set_array_item(0, 1.0, a_state->trans_out);
    break;
  case SPACER_STATE :
    // Transitions go to self, motifs and end (except for beginning
    // multi-state spacers)
    a_state->ntrans_out = (i_position < spacer_states -1 ) ? 2 : nmotifs + 2;
    a_state->itrans_out = (int *)mm_malloc(sizeof(int) * a_state->ntrans_out);
    a_state->trans_out = allocate_array(a_state->ntrans_out);

    // The first transition is a self-transition.
    a_state->itrans_out[0] = i_state;
    set_array_item(0, self_trans(expected_length), a_state->trans_out);

    // For multi-state spacers, outgoing transition to next state.
    if (i_position < spacer_states - 1) {
      a_state->itrans_out[1] = i_state + 1;
      set_array_item(1, 1-self_trans(expected_length), a_state->trans_out);
    } else {
      double out_p = (1 - self_trans(expected_length))/(nmotifs+1);
      // Out to each motif start.
      for (j_motif = 0; j_motif < nmotifs; j_motif++) {	
        a_state->itrans_out[j_motif+1] = 
          motif_index(j_motif+1, FALSE, num_spacers, spacer_states, motifs);
	set_array_item(j_motif+1, out_p, a_state->trans_out);
      }
      // Out to end state.
      a_state->itrans_out[j_motif+1] = 
        motif_index(nmotifs, TRUE, num_spacers, spacer_states, motifs) + 1;
      set_array_item(j_motif+1, out_p, a_state->trans_out);
    }
    break;
  case END_STATE :
    a_state->ntrans_out = 0;
    a_state->itrans_out = NULL;
    a_state->trans_out = NULL;
    break;
  }
} // build_star_state

/*************************************************************************
 * Build a star topology HMM.
 *************************************************************************/
void build_star_hmm
  (ARRAY_T*  background,
   int       spacer_states, 
   MOTIF_T*  motifs,
   int       nmotifs,
   BOOLEAN_T fim,
   MHMM_T**  the_hmm)
{
  int       motif_states; /* Total length of the motifs. */
  int       num_spacers;  /* Total number of spacer states. */
  int       num_states;   /* Total number of states in the model. */
  int       i_motif;      /* Index of the current "from" motif. */
  int       i_position;   /* Index within the current motif or spacer. */
  int       i_state = 0;  /* Index of the current state. */

  /* Count the width of the motifs. */
  for (motif_states = 0, i_motif = 0; i_motif < nmotifs; i_motif++)
    motif_states += (motifs[i_motif]).length;
  // Only 1 spacer.
  num_spacers = 1;
  /* Total states = motifs + spacer_states + begin/end */
  num_states = motif_states + (num_spacers * spacer_states) + 2;
  /* fprintf(stderr, "motif_states=%d num_spacers=%d num_states=%d\n",
	  motif_states, num_spacers, num_states); */

  /* Allocate the model. */
  allocate_mhmm(num_states, the_hmm);

  /* Record that this is a completely connected model. */
  (*the_hmm)->type = STAR_HMM;

  /* Record the number of motifs in the model. */
  (*the_hmm)->num_motifs = nmotifs;

  /* Record the number of states in the model. */
  (*the_hmm)->num_states = num_states;
  (*the_hmm)->num_spacers = 1;
  (*the_hmm)->spacer_states = spacer_states;

  /* Put the alphabet into the model. */
  (*the_hmm)->alph_size = get_alph_size(ALPH_SIZE);
  (*the_hmm)->ambigs = get_alph_size(AMBIG_SIZE);
  strcpy((*the_hmm)->alphabet, get_alphabet(FALSE));

  // Put the background distribution into the model.
  copy_array(background, (*the_hmm)->background);

  /* Build the begin state. */
  build_star_state(START_STATE,
		   i_state,
		   0, // expected length
		   NULL,
		   0, // Number of sites.
		   0, // alphabet size
		   NON_MOTIF_INDEX,
		   NON_MOTIF_ID,
		   NON_MOTIF_POSITION,
		   0,
		   nmotifs,
		   spacer_states,
		   motifs,
		   &((*the_hmm)->states[i_state]));
  i_state++;

  // Build the spacer state (state 0).  Allow multi-state spacers.
  for (i_position = 0; i_position < spacer_states; i_position++) {
    build_star_state(SPACER_STATE, 
		     i_state, 
		     DEFAULT_SPACER_LENGTH,
		     background,
		     SPACER_NUMSITES,
		     get_alph_size(ALL_SIZE),
		     NON_MOTIF_INDEX,
		     NON_MOTIF_ID,
		     i_position,
		     0, 
		     nmotifs,
		     spacer_states,
		     motifs,
		     &((*the_hmm)->states[i_state]));
    i_state++;
  }

  /* Build the motif states. */
  for (i_motif = 0; i_motif < nmotifs; i_motif++) {
    MOTIF_T *this_motif = &(motifs[i_motif]);
    assert(this_motif->length != 0);
    i_position = 0;
    build_star_state(START_MOTIF_STATE, 
		     i_state,
		     0, // Expected spacer length.
		     get_matrix_row(i_position, this_motif->freqs),
		     this_motif->num_sites,
		     get_alph_size(ALL_SIZE),
		     i_motif,
		     get_motif_id(this_motif),
		     i_position,
		     this_motif->length,
		     nmotifs,
		     spacer_states,
		     motifs,
		     &((*the_hmm)->states[i_state]));
    i_state++;
    for (i_position = 1; i_position < this_motif->length - 1; i_position++) {
      build_star_state(MID_MOTIF_STATE, 
		       i_state,
		       0, // Expected spacer length. 
		       get_matrix_row(i_position, this_motif->freqs),
		       this_motif->num_sites,
		       get_alph_size(ALL_SIZE),
		       i_motif,
		       get_motif_id(this_motif),
		       i_position, 
		       this_motif->length,
		       nmotifs,
		       spacer_states,
		       motifs,
		       &((*the_hmm)->states[i_state]));
      i_state++;
    }
    build_star_state(END_MOTIF_STATE, 
		     i_state,
		     0, // Expected spacer length.
		     get_matrix_row(i_position, this_motif->freqs),
		     this_motif->num_sites,
		     get_alph_size(ALL_SIZE),
		     i_motif,
		     get_motif_id(this_motif),
		     i_position,
		     this_motif->length,
		     nmotifs,
		     spacer_states,
		     motifs,
		     &((*the_hmm)->states[i_state]));
    i_state++;
  }

  /* Build the end state. */
  build_star_state(END_STATE, 
		   i_state,
		   0, // Expected spacer length.
		   NULL, // Emissions
		   0, // Number of sites.
		   0, // Size of emissions.
		   nmotifs,
		   NON_MOTIF_ID,
		   NON_MOTIF_POSITION,
		   0,
		   nmotifs,
		   spacer_states,
		   motifs,
		   &((*the_hmm)->states[i_state]));
  i_state++;

  /* Convert spacers to FIMs if requested. */
  if (fim) {
    convert_to_fims(*the_hmm);
  }

  /* Fill in the transition matrix. */
  build_transition_matrix(*the_hmm);
} // build_star_hmm

