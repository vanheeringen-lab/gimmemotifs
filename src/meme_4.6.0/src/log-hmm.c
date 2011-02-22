/**************************************************************************
 * FILE: log-hmm.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 1-28-98
 * PROJECT: MHMM
 * DESCRIPTION: Convert between a normal HMM and a log HMM.
 **************************************************************************/
#include <math.h>
#include <string.h>
#include "log-hmm.h"
#include "utils.h"
#include "mhmm-state.h"
#include "array.h"
#include "matrix.h"
#include "write-mhmm.h"
#include "subst-matrix.h"	/* substitution matrix routines */
#include "alphabet.h"

#ifndef LOG_DEBUG
#define LOG_DEBUG 0
#endif

/**************************************************************************
 * See .h file for description.
 **************************************************************************/
static void convert_to_from_log_state(
  BOOLEAN_T to_log, 
  BOOLEAN_T zero_spacer_emit_lo, // Set spacer emission log-odds = 0?
  ARRAY_T*  background,          // Background emission distribution.
  MATRIX_T* target_freq,         // Target frequency matrix.
  double beta,	                 // Weight on pseudocounts.
  MHMM_STATE_T* source_state,         // The state to be converted.
  MHMM_STATE_T* target_state     // The converted form of the given state.
) {
  double alpha = source_state->alpha;
  int i_emit;

  // If there is no background, just copy from emit to odds or vice versa.
  if (background == NULL) {
    if (to_log) {
      copy_array(source_state->emit, source_state->emit_odds);
    } else {
      copy_array(source_state->emit_odds, source_state->emit);
    }
  }

  // Compute odds for the emission distribution in the source state.
  else if (to_log) {
    ARRAY_T *g;      // Pseudocount frequencies.

    // Get the pseudocount frequencies for modifying observed frequencies.
    g = get_pseudocount_freqs(source_state->emit, background, target_freq);

    for (i_emit = 0; i_emit < get_alph_size(ALPH_SIZE); i_emit++) {

      // Start and end states have unity odds; spacers, too, if requested.
      if ((source_state->type == START_STATE) ||
          (source_state->type == END_STATE) ||
          (zero_spacer_emit_lo && source_state->type == SPACER_STATE)) { 
         set_array_item(i_emit, 1.0, source_state->emit_odds);
      }

      // Compute the adjusted odds.
      else {
        double fi = get_array_item(i_emit, source_state->emit);  // observed freq
        double gi = get_array_item(i_emit, g);      // pseudocount frequency
        double q = (alpha * fi + beta * gi)/(alpha + beta);  // adjusted frequency
        double b = get_array_item(i_emit, background);    // background frequency
        set_array_item(i_emit, q / b, source_state->emit_odds);
      }
    }
    free_array(g);      // Free pseudocounts.

  } // to_log

  // Allocate dynamic space in the new state.
  if (to_log) {
    target_state->emit = allocate_array(get_alph_size(ALL_SIZE));
    target_state->emit_odds = allocate_array(get_alph_size(ALL_SIZE));
    target_state->itrans_out = 
      (int *)mm_malloc(sizeof(int) * source_state->ntrans_out);
    target_state->trans_out = allocate_array(source_state->ntrans_out);
    target_state->itrans_in = 
      (int *)mm_malloc(sizeof(int) * source_state->ntrans_in);
    target_state->trans_in = allocate_array(source_state->ntrans_in);
  }

  // Copy the state type.
  target_state->type = source_state->type;

  // Copy outgoing transitions.
  target_state->ntrans_out = source_state->ntrans_out;
  copy_int_array(
    source_state->ntrans_out,
    source_state->itrans_out,
    target_state->itrans_out
  );
  convert_to_from_log_array(
    to_log,
    source_state->trans_out, 
    target_state->trans_out
  );

  // Copy incoming transitions.
  target_state->ntrans_in = source_state->ntrans_in;
  copy_int_array(
    source_state->ntrans_in,
    source_state->itrans_in,
    target_state->itrans_in
  );
  convert_to_from_log_array(
    to_log,
    source_state->trans_in, 
    target_state->trans_in
  );

  // Copy the emission distribution.
  convert_to_from_log_array(
    to_log,
    source_state->emit,
    target_state->emit
  );
  convert_to_from_log_array(
    to_log,
    source_state->emit_odds,
    target_state->emit_odds
  );

  // Set the odds for the ambiguous characters as the (unweighted) average
  // of the odds for the characters they match.
  //FIXME: should be weighted average for statistics to really be kosher.
  fill_in_ambiguous_chars(FALSE, target_state->emit_odds);

  // Copy the remaining stuff.
  target_state->i_motif = source_state->i_motif;
  target_state->w_motif = source_state->w_motif;
  strcpy(target_state->motif_id, source_state->motif_id);
  target_state->i_position = source_state->i_position;
  target_state->id_char = source_state->id_char;

  if (LOG_DEBUG) {
    for (i_emit = 0; i_emit < get_alph_size(ALL_SIZE); i_emit++) {
      printf("%5.3f ", get_array_item(i_emit, target_state->emit_odds));
    }
    printf("\n");
  }
}

/**************************************************************************

	Set the gap opening and extension costs in the adjacency
	lists and transition matrix.

              extend
              --->
      open/2 |   |  open/2 
	----->spacer----->

**************************************************************************/
static void set_gap_costs(
   double        gap_open,     // Cost to open a gap.
   double        gap_extend,   // Cost to extend a gap.
   MHMM_T*       hmm)          // The hmm.
{
  int i_state, i;
  double cost;

  for (i_state = 0; i_state < hmm->num_states; i_state++) {
    MHMM_STATE_T* state = hmm->states + i_state;
 
    // Set spacer incoming transitions to -gap_open/2 if it is >= 0
    // and -gap_extend if it is a self-loop.
    for (i=0; i<state->ntrans_in && (gap_open >= 0 || gap_extend >= 0); i++) {
      int j_state = state->itrans_in[i];			// from state
      if (i_state != j_state && gap_open < 0) continue;		// only do self-loop?
      if (state->type == SPACER_STATE) {			// to state is spacer
        cost = (i_state==j_state) ? -gap_extend : -gap_open/2;	// self-loop?
      } else if (hmm->states[j_state].type == SPACER_STATE) {	// from state is spacer
        cost = -gap_open/2;
      } else {							// neither is spacer
        continue;
      }
      set_array_item(i, cost, state->trans_in);			// set adj list
      set_matrix_cell(j_state, i_state, cost, hmm->trans);	// set xition matrix
    } // incoming state

    // Set self-loop transitions to -gap_extend if gap_extend >= 0
    // and set rest of outgoing spacer transitions to -gap_open/2.
    for (i=0; i<state->ntrans_out && (gap_open >= 0 || gap_extend >= 0); i++) {
      int j_state = state->itrans_out[i];			// to state
      if (state->type == SPACER_STATE) {			// from state is spacer
        cost = (i_state==j_state) ? -gap_extend : -gap_open/2;	// self-loop?
      } else if (hmm->states[j_state].type == SPACER_STATE) {	// to state is spacer
        cost = -gap_open/2;					// rest of open cost
      } else {							// neither is spacer
        continue;
      }
      set_array_item(i, cost, state->trans_out);		// set adj list
      set_matrix_cell(i_state, j_state, cost, hmm->trans);	// set xition matrix
    } // incoming state

  } // i_state

} // set_gap_costs

/**************************************************************************
 * Convert an HMM to or from log form.
 **************************************************************************/
void convert_to_from_log_hmm
  (BOOLEAN_T to_log,        // Convert to log form?
   BOOLEAN_T zero_spacer_emit_lo,     // Set spacer emission log-odds = 0?
   double    gap_open,      // Cost to open a gap; ignored if < 0
   double    gap_extend,    // Cost to extend a gap; ignored if < 0
   ARRAY_T*  background,    // The background distribution.
   char *    sc_filename,   // Name of score file (ignored if background NULL).
   int       pam_dist,      // PAM distance (ignored if sc_filename not NULL).
   double    beta,	    // Weight on pseudocounts.
   MHMM_T*   source_hmm,    // The HMM to convert.
   MHMM_T**  target_hmm)    // The same HMM in log/log-odds form 
{
  int i_state; /* Index of the current state. */
  MATRIX_T *target_freq = NULL; /* The target frequency matrix */

  /* Allocate space for the new HMM. */
  if (to_log) {
    allocate_mhmm(source_hmm->num_states, target_hmm);
  }

  // Copy the header information.
  (*target_hmm)->type = source_hmm->type;
  (*target_hmm)->log_odds = to_log;
  (*target_hmm)->num_motifs = source_hmm->num_motifs;
  (*target_hmm)->num_states = source_hmm->num_states;
  (*target_hmm)->num_spacers = source_hmm->num_spacers;
  (*target_hmm)->spacer_states = source_hmm->spacer_states;
  (*target_hmm)->alph_size = source_hmm->alph_size;
  (*target_hmm)->ambigs = source_hmm->ambigs;
  strcpy((*target_hmm)->alphabet, source_hmm->alphabet);

  // Create the substitution target frequency matrix if background given.
  if (background != NULL) {
    int asize = source_hmm->alph_size;
    BOOLEAN_T dna = (which_alphabet() == DNA_ALPH);
    char *alphabet = source_hmm->alphabet;
    char save_char = alphabet[asize];
    alphabet[asize] = '\0';			// remove ambigs 
    target_freq = get_subst_target_matrix(sc_filename, dna, pam_dist, 
      source_hmm->alphabet, background);
    alphabet[asize] = save_char;		// replace ambigs
  } // target frequency matrix 


  // Copy the states, converting to logs along the way.
  for (i_state = 0; i_state < source_hmm->num_states; i_state++) {
    (source_hmm)->states[i_state].alpha = 20;
    if (LOG_DEBUG) {
      printf("i_motif: %d i_position: %d\n", (source_hmm)->states[i_state].i_motif,
	(source_hmm)->states[i_state].i_position);
    }
    convert_to_from_log_state(to_log, zero_spacer_emit_lo,
			      background, target_freq, beta,
			      &(source_hmm->states[i_state]),
			      &((*target_hmm)->states[i_state]));
  }

  // Copy the background distribution, and convert to logs.
  convert_to_from_log_array(to_log,
			    source_hmm->background,
			    (*target_hmm)->background);

  // Copy the transition matrix, converting to logs along the way.
  convert_to_from_log_matrix(to_log,
			     source_hmm->trans,
			     (*target_hmm)->trans);

  // Set the transition costs to -gap_open and -gap_extend if
  // these are >= 0 and converting to logs.
  if (to_log && (gap_open >= 0 || gap_extend >= 0))
    set_gap_costs(gap_open, gap_extend, *target_hmm);

  // Copy the hot_states.
  (*target_hmm)->hot_states = NULL;
  if (source_hmm->num_hot_states > 0) {
    mm_resize((*target_hmm)->hot_states, source_hmm->num_hot_states, int);
  }
  (*target_hmm)->num_hot_states = source_hmm->num_hot_states;

  if (LOG_DEBUG) {
    write_mhmm(QUIET_VERBOSE, *target_hmm, stderr);
  }

  // Free local dynamic memory.
  free_matrix(target_freq);
}


#ifdef MAIN
#include "simple-getopt.h"
#include "read-mhmm.h"

VERBOSE_T verbosity = INVALID_VERBOSE;

/**************************************************************************
 * int main
 **************************************************************************/
int main (int argc, char *argv[])
{
  char *    hmm_filename= NULL;  // File containing the HMM.
  FILE *    hmm_file;
  MHMM_T *  the_hmm;        /* The HMM itself. */
  MHMM_T *  the_log_hmm;    /* The HMM, with probs converted to logs. */

  // Define command line options.
  cmdoption const options[] = {
      {"verbosity", REQUIRED_VALUE}
  };
	int option_count = 1;
  int option_index = 0;

  // Define the usage message.
  char      usage[1000] = "";
  strcat(usage, "USAGE: log-hmm [options] <HMM file>\n");
  strcat(usage, "\n");
  strcat(usage, "   Options:\n");
  strcat(usage, "     --verbosity 1|2|3|4|5 (default=2)\n");
  strcat(usage, "\n");

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
      die("Error process command line options (%s)\n", message);
    }

    if (strcmp(option_name, "verbosity") == 0) {
      verbosity = (VERBOSE_T)atoi(option_value);
    }
  }

  // Read the single required argument.
  if (option_index + 1 != argc) {

    fprintf(stderr, "%s", usage);
    exit(1);
  }
  hmm_filename = argv[option_index];

  /* Read the model. */
  read_mhmm(hmm_filename, &the_hmm);

  /* Convert back and forth. */
  convert_to_from_log_hmm(TRUE,  // To log form.
			  FALSE, // Set spacer emissions to zero?
			  -1,    // Ignore gap_open.
			  -1,    // Ignore gap_extend.
			  the_hmm->background,
			  NULL,  // No score file.
			  250.0, // PAM distance.
			  10,     // beta ?
			  the_hmm,
			  &the_log_hmm);
  convert_to_from_log_hmm(FALSE, // From log form.
			  FALSE, // Set spacer emissions to zero?
			  -1,    // Ignore gap_open.
			  -1,    // Ignore gap_extend.
			  the_hmm->background,
			  NULL,  // No score file.
			  250.0, // PAM distance.
			  10,     // beta
			  the_log_hmm,
			  &the_hmm);

  /* Print the HMM. */
  write_mhmm(verbosity, the_hmm, stdout);

  /* Free dynamic memory. */
  free_mhmm(the_hmm);
  free_mhmm(the_log_hmm);
  return(0);
}
#endif /* MAIN */
