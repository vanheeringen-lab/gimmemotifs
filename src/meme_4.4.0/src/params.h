/***************************************************************************
 * FILE: params.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 8/23/01
 * PROJECT: MHMM
 * COPYRIGHT: 2001-2008, WSN
 * DESCRIPTION: Type for keeping track of different kinds of model parameters.
 ***************************************************************************/
#include "utils.h"
#include "mhmm-state.h"

/* An enumerated type for specifying particular parameter sets. */
typedef enum {
    INVALID_TRAIN,
    NONE_TRAIN,     /* Don't train anything. */
    TRANS_TRAIN,    /* Transitions from end of each motif. */
    LENGTH_TRAIN,   /* Self-transitions at the spacers. */
    SPACER_TRAIN,   /* Emissions within spacers. */
    MOTIF_TRAIN,    /* Emissions within motifs. */
    ALL_TRAIN,      /* Train everything. */
} TRAIN_T;
extern char*  TRAIN_STRS[];
extern int NUM_TRAIN_T;

/***************************************************************************
 * Add to the global list of parameters to train.
 ***************************************************************************/
void add_train_type
  (TRAIN_T new_train_type);

/***************************************************************************
 * Are any of the emission distributions being trained?
 ***************************************************************************/
BOOLEAN_T train_emissions();

/***************************************************************************
 * Are any of the transition distributions being trained?
 ***************************************************************************/
BOOLEAN_T train_transitions();

/***************************************************************************
 * Is a particular parameter subset being trained?
 ***************************************************************************/
BOOLEAN_T do_train
  (TRAIN_T train_what);

/***************************************************************************
 * Get a string representing all the parameters currently being trained.
 ***************************************************************************/
char*  train_string();

/***************************************************************************
 * Determine whether transitions from a given state should be trained.
 ***************************************************************************/
BOOLEAN_T train_state_transition
  (MHMM_STATE_T * this_state);

/***************************************************************************
 * Determine whether emissions from a given state should be trained.
 ***************************************************************************/
BOOLEAN_T train_state_emissions
  (MHMM_STATE_T * this_state);

