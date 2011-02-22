/****************************************************************************
 * FILE: motiph_scoring.h
 * AUTHOR: William Stafford Noble, Charles E. Grant
 * CREATE DATE: 12/03/2004
 * PROJECT: EVOMCAST
 * DESCRIPTION: Evaluate the log-likelihood of a phylogenetic tree
 * given an alignment and an evolutionary model.
 * COPYRIGHT: 2004, UW
 ****************************************************************************/
#ifndef MOTIPH_SCORING_H
#define MOTIPH_SCORING_H

#ifndef MOTIPH
#define MOTIPH 0
#endif

#ifndef SHADOW
#define SHADOW 0
#endif

#ifndef SCAN
#define SCAN 0
#endif

#include "alignment.h"
#include "cisml.h"
#include "evomodel.h"
#include "motif.h"
#include "seq.h"
#include "substmatrix-table.h"
#include "tree.h"
#include "scored-sites.h"
#include "pssm.h"

/****************************************************************************
 *  Calculate the likelihood of a phylogenetic tree at one column 
 *  in an alignment.
 ****************************************************************************/
double site_likelihood(
  char* alignment_col,
  STRING_LIST_T* seq_names,
  TREE_T* tree,
  MODEL_TYPE_T model_type,
  ARRAY_T* priors,
  struct substmatrix_array* substmatrix_array,
  GAP_SUPPORT_T gap_support
);

/*************************************************************************
 * Calculate the log odds score for each motif-sized window at each 
 * site in the specified sequence of the alignment using the given model.
 * Returns number of sites scored (not skipped).
 *************************************************************************/
int score_sequence_in_alignment(
  int             ref_seq_index,
  ALIGNMENT_T*    alignment,
  char*           motif_id,
  TREE_T*         tree,
  int             window_size,
  EVOMODEL_T**    models,
  ARRAY_T*        bg_freqs,
  PSSM_T*         pssm,
  int*            coord_conv_table,
  GAP_SUPPORT_T   gap_support,
  double          gap_cost,
  double          pthresh,
  SCANNED_SEQUENCE_T* scanned_seq
);

/*************************************************************************
 * Calculate the BLS score for each motif-sized window at each 
 * site in the specified sequence of the alignment.
 * Returns number of sites scored (not skipped).
 *************************************************************************/

int bls_score_sequence_in_alignment(
  int             ref_seq_index,
  ALIGNMENT_T*    alignment,
  char*           motif_id,
  TREE_T*         tree,
  int             window_size,
  ARRAY_T*        bg_freqs,
  EVOMODEL_T**    models, 
  PSSM_T*         pssm,
  char*           inverse_motif_id,
  EVOMODEL_T**    inverse_models, 
  PSSM_T*         inverse_pssm,
  int*            coord_conv_table,
  GAP_SUPPORT_T   gap_support,
  double          gap_cost,
  double          pthresh,
  int          	  distThreshold,
  SCANNED_SEQUENCE_T* scanned_seq
);

/*************************************************************************
 *  Build a position specific scoring matrix (PSSM) for an alignment
 *  and an array of evolutionary models. The first model in the array
 *  is assumed to be the background model. The PSSM will have one row
 *  for each model, and one column for each possible alignment column.
 *  The elements of the PSSM are the log-odds scores for the corresponding
 *  model and alignment column.
 *
 *  Caller is responsible for freeing the matrix.
 *************************************************************************/
MATRIX_T* build_alignment_pssm_matrix(
  STRING_LIST_T* seq_names,
  int num_models,
  EVOMODEL_T** models,
  TREE_T* tree,
  GAP_SUPPORT_T gap_support
);

/*************************************************************************
 * Get the (index of) the base in the sequence in the alignment column
 * corresponding to the given tree node.
 *************************************************************************/
int get_base_from_node(
  TREE_T* node,
  STRING_LIST_T* seq_names,
  char* alignment_col,
  char* alphabet,
  int alph_size
);

#endif
