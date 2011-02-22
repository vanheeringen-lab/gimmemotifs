/********************************************************************
 * FILE: build-hmm.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 7-14-97
 * PROJECT: MHMM
 * COPYRIGHT: 1997-2008, WSN
 * DESCRIPTION: Build an HMM in memory.
 ********************************************************************/
#ifndef BUILD_HMM_H
#define BUILD_HMM_H
#include "matrix.h"
#include "utils.h"
#include "motif.h"
#include "order.h"
#include "mhmm-state.h"


/*************************************************************************
 * Verify that the rows of a transition matrix sum to 1.0. 
 *************************************************************************/
BOOLEAN_T verify_trans_matrix
  (BOOLEAN_T log_form,    /* Is the transition matrix in log form? */
   int       num_states,  /* Number of states in the (square) matrix. */
   MATRIX_T* trans);      /* The matrix. */

/*************************************************************************
 * Build a linear HMM.
 *************************************************************************/
void build_linear_hmm
  (ARRAY_T*  spacer_dist,
   ORDER_T*  order_spacing,
   int       spacer_length, 
   MOTIF_T*  motifs,
   int       nmotifs,
   BOOLEAN_T fim,
   MHMM_T**  the_hmm);

/*************************************************************************
 * Build a completely connected HMM.
 *************************************************************************/
void build_complete_hmm
  (ARRAY_T*  spacer_dist,
   int       spacer_length, 
   MOTIF_T*  motifs,
   int       nmotifs,
   MATRIX_T* transp_freq,
   MATRIX_T* spacer_ave,
   BOOLEAN_T fim,
   MHMM_T**  the_hmm);

/*************************************************************************
 * Build a star topology HMM.
 *************************************************************************/
void build_star_hmm
  (ARRAY_T*  background,
   int       spacer_length,
   MOTIF_T*  motifs,
   int       nmotifs,
   BOOLEAN_T fim,
   MHMM_T**  the_hmm);

#endif

