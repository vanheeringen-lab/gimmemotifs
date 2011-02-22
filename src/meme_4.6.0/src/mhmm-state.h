/***********************************************************************
 * FILE: mhmm-state.h
 * AUTHOR: William Stafford Noble
 * PROJECT: MHMMM
 * COPYRIGHT: 1997-2008, WSN
 * VERSION: $Revision: 1.2 $
 * DESCRIPTION: Data structure for representing an HMM.
 ***********************************************************************/
#ifndef MHMM_STATE_H
#define MHMM_STATE_H

// referred to by pssm.h
typedef struct mhmm MHMM_T;

#include "matrix.h"
#include "array.h"
#include "utils.h" 	// Define PROB_T. 
#include "motif.h"	// Needed for MAX_MOTIF_ID_LENGTH
#include "pssm.h"

// Default pseudocounts for emissions and transitions.
#define DEFAULT_SPACER_LENGTH  100 // Length of spacer region.
#define DEFAULT_SPACER_STATES  1   // Number of states in spacer.
#define DEFAULT_TRANS_PSEUDO   0.1
#define DEFAULT_SPACER_PSEUDO  0.0
#define DEFAULT_MOTIF_PSEUDO  0.0
#define DEFAULT_EMIT_PSEUDO    0.1

// Enumerated type for different kinds of states. 
typedef enum {INVALID_STATE, START_STATE, START_MOTIF_STATE, MID_MOTIF_STATE,
	      END_MOTIF_STATE, SPACER_STATE, END_STATE} STATE_T;
extern char * STATE_STRS[];
extern int NUM_STATE_T;

typedef struct mhmm_state {
  STATE_T  type;
  int      ntrans_out;   // Number of transitions from this state. 
  int *    itrans_out;   // Indices of states reachable from this state. 
  ARRAY_T* trans_out;    // Transition probabilities from this state. 
  int      ntrans_in;    // Number of transitions to this state. 
  int *    itrans_in;    // Indices of states that can reach state. 
  ARRAY_T* trans_in;     // Transition probabilities to this state. 
  ARRAY_T* emit;         // Emission probability distribution. 
  double   num_sites;   // Number of sites associated with this emission.
  /* i_motif goes incrementally from 0; motif_id is the motif number
   * as determined by MEME.  It is arbitrary. */
  int      i_motif;      // Which motif is this state in? (-1=spacer) 
  int      w_motif;	 // Motif width.
  char     motif_id[MAX_MOTIF_ID_LENGTH + 1]; /* Motif ID.  If
				       signed, then sign determines strand. */
  int      i_position;   // What position in that motif is this state? 
  char     id_char;      // Character above this state when printing
                         // the motif ID centered w.r.t. the entire motif.
  // The following fields are only valid if the state is in log-odds format.
  ARRAY_T* emit_odds;    // Odds of emission probability distribution. 
  double   alpha;	 // Effective number of sequences in this motif.
  // In addition, the following fields are only valid at START_MOTIF_STATE.
  PSSM_T* pssm;		// PSSM for motif.
  PSSM_T* npssm;	// PSSM for motif on negative strand.
  // In addition, the following fields are only valid if using p-value scoring.
  //ARRAY_T* pvalues;      // Lookup table for p-values.
  double min_sig_score;	 // Minimum significant score (pv < p_threshold). 
  double min_pvalue;	 // P-value of largest possible score.
} MHMM_STATE_T;

// Enumerated type for different types of HMMs. 
typedef enum {INVALID_HMM, LINEAR_HMM, COMPLETE_HMM, STAR_HMM} HMM_T;
extern char * HMM_STRS[];
extern int NUM_HMM_T;

struct mhmm {
  HMM_T          type;          // LINEAR_HMM or COMPLETE_HMM 
  BOOLEAN_T      log_odds;      // Is this HMM in log-odds form? 
  int		 num_motifs;	// Number of motifs in this HMM.
  int            num_states;    // Total number of states. 
  int            num_spacers;   // Number of spacer states. 
  int            spacer_states; // Number of states per spacer. 
  int            alph_size;     // Number of symbols in the alphabet. 
  int            ambigs;        // Number of ambiguous symbols. 
  char*          alphabet;      // The alphabet without ambigs.
  ARRAY_T*       background;    // Background probability distribution.
  char*          description;   // Descriptive text.
  char*          motif_file;    // Name of motif file used in creation.
  char*          sequence_file; // Name of sequence file used in training.
  MHMM_STATE_T * states;        // Individual states. 
  MATRIX_T*      trans;         // Transition probability matrix. 
  // trans[i][j] = probability of transitioning from state i to state j. 
  // The hot_states list is all the states that need to be considered
  // during dynamic programming.  Some states can be skipped in motif-scoring
  // mode, so using the hot_states list saves loads of time.
  int*		 hot_states;	// States used in dp.
  int		 num_hot_states;// Number of hot states.
};

/************************************************************************
 * Dynamically allocate memory for one motif-based HMM.
 ************************************************************************/
void allocate_mhmm
  (int      num_states,
   MHMM_T** an_mhmm);

/************************************************************************
 * Get the motif ID, with or without strand info.
 *
 * Returns a pointer to local static storage.
 ************************************************************************/
char* get_state_motif_id
  (BOOLEAN_T     stranded,
   MHMM_STATE_T* this_state);

/************************************************************************
 * Determine the strandedness of a given motif.  Returns '+' or '-' if
 * the motif has strand and '.' otherwise.
 *************************************************************************/
char get_strand
  (MHMM_STATE_T* this_state);

/************************************************************************
 * Make a copy of one HMM into a previously allocated HMM.
 ************************************************************************/
void copy_mhmm
  (MHMM_T*        an_mhmm,
   MHMM_T*        new_mhmm);

/************************************************************************
 * Mix two HMMs according to the given mixing parameter.
 *
 * The transition and emission probability distributions in the old
 * HMM are mixed with a portion of the new HMM's distributions.
 *
 * Assumes that both HMMs are completely allocated and both are in log
 * form.
 ************************************************************************/
void mix_mhmm
  (const float    mixing,    /* Percent of new HMM that is the old HMM. */
   MHMM_T* const new_mhmm,  /* New HMM to be mixed into the old one. */
   MHMM_T*       old_mhmm); /* Target HMM to be modified via mixing. */

/************************************************************************
 * Free the memory used by one mhmm. 
 ************************************************************************/
void free_mhmm
  (MHMM_T*an_mhmm);

/************************************************************************
 * Given an HMM with the transition matrix filled in, compute the
 * state-by-state transition arrays.
 ************************************************************************/
void compute_ins_and_outs
  (MHMM_T*the_hmm,
   const BOOLEAN_T log_form); /* Is the transition matrix in log form? */

/*************************************************************************
 * Given an HMM state, choose the most frequently occuring character.
 *
 * RETURN: The most common character in the given state's log-odds 
 *         emission distribution.
 *************************************************************************/
char choose_consensus
  (BOOLEAN_T      log_space, /* Choose based on log-odds probs. */
   MHMM_STATE_T * a_state);  /* State to be analyzed. */

/*************************************************************************
 * Determine whether the transition matrix of a given HMM contains
 * free-insertion modules; i.e., spacer states in which the
 * self-transition and the exit transition are 1.0.
 *************************************************************************/
BOOLEAN_T has_fims
  (MHMM_T* the_hmm);

/************************************************************************
 * Count the number of motifs in an HMM.
 *************************************************************************/
int get_num_motifs
  (MHMM_T* the_hmm);

#endif 
