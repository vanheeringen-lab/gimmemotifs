/*
 * @file calculate_p_y.h
 *
 * DESCRIPTION HERE...
 *
 * $Id: calculate_p_y.h 1499 2006-11-27 06:01:08Z tom $
 *
 * $Log$
 * Revision 1.1.2.2  2006/02/09 23:51:17  twhitington
 * Removed compile-time errors from branching_search module.
 *
 * Revision 1.1.2.1  2006/02/07 02:02:53  twhitington
 * Added calculate_p_y module to repository.
 *
 *
 */

#ifndef CALCULATE_P_Y_H
#define CALCULATE_P_Y_H
#include "meme.h"
#include "seed_diffs.h"
#include "sp_matrix.h"

/**
 * get_pY
 * 
 * Compute log Pr(Y_ij | theta_1) for all positions (i, j) in the dataset.
 *
 */
void get_pY(
  DATASET *dataset,  ///< the dataset
  int *theta_1[MAXSITE], ///< integer log theta_1
  int w,             ///< width of motif
  int pYindex        ///< which pY array to use
);


/**
 * combine_strands
 *
 * This function operates on an array of SAMPLE objects. For each sample, the
 * first pY scratch array is instantiated to contain the best pY values from
 * the positive and negative strands.
 *
 */
void combine_strands(
  SAMPLE **samples,  ///< An array of pointers to SAMPLE objects, whose first
                     ///< scratch pY arrays will be modified.
  int n_samples,     ///< The length of the array of SAMPLE pointers.
  int w              ///< The length of the motifs under which the pY
                     ///< arrays were initially instaniated.
);


extern void create_lmotif (
  char *seed_str,    ///< The seed as a string
  int lmap[MAXALPH][MAXALPH], ///< Matrix representing a seq-to-theta mapping fn
  int *lmotif[MAXSITE], ///< The resulting lmotif - OUT
  int *mot_width     ///< The width of the motif - OUT
);


extern void init_lmotif(
  int w,             ///< width of site
  char *res,         ///< (encoded) letters of subsequence
  int *theta_1[MAXSITE], ///< theta_1
  int lmap[MAXALPH][MAXALPH] ///< matrix of frequency vectors
);


extern void evaluate_seed_DP (
  char *new_seed,    ///< The string representing the new seed object.
  SEED_DIFFS *s_diffs, ///< An object representing the differences between the
                     ///< seed of interest and the "previous" seed.
  int lmap[MAXALPH][MAXALPH], ///< Matrix representing a seq-to-theta mapping fn
  MOTYPE mtype,       ///< The type of sequence model evaluation is performed
                     ///< under.
  BOOLEAN ic,        ///< Consider inverse complement in DNA.
  DATASET *dataset,  ///< The dataset, containing the sequences used during
                     ///< evaluation of the current seed.
  SP_MATRIX *sp_mat  ///< A matrix of starting points, into which the evaluated
                     ///< new seed will be placed.
);


extern void next_pY_branching(
  int *lmotif_new[MAXSITE], ///< The lmotif corresponding to the new seed.
  SEED_DIFFS *s_diffs, ///< The differences between the new and old seeds.
  DATASET *dataset,  ///< The dataset of sequences. Contains the pY arrays.
  int pYindex        ///< Variable indicating which pY array to use.
);


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */


#endif


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
