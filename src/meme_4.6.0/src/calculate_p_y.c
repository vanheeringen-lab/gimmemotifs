/**
 * @file calculate_p_y.c
 *
 * This module relates to generic pY calculations. "pY" is Pr(Y_ij | theta).
 * The pY calculation functions contained can be applied in multiple scenarios.
 * In contrast, other modules (eg branching_search.c and subseq7.c) might include
 * pY calculation functions specific to a process. Eg. branching search uses
 * dynamic programming, which dictates a specific approach to some pY
 * calculations.
 *
 * $Id: calculate_p_y.c 4278 2009-12-23 09:58:37Z james_johnson $
 *
 * $Log$
 * Revision 1.1.2.2  2006/02/09 23:51:15  twhitington
 * Removed compile-time errors from branching_search module.
 *
 * Revision 1.1.2.1  2006/02/07 02:02:52  twhitington
 * Added calculate_p_y module to repository.
 *
 *
 */


#include "calculate_p_y.h"
#include <assert.h>
#include <limits.h>
#include "meme.h"
#include "seed.h"

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
)
{
  int i, j, k;
  int n_samples = dataset->n_samples;
  SAMPLE **samples = dataset->samples;
  
  for (i=0; i<n_samples; i++) {         /* sequence */
    SAMPLE *s = samples[i];
    int lseq = s->length;
    char *res = pYindex<2 ? s->res : s->resic;  /* integer sequence */
    int *pY = s->pY[pYindex];           /* log p(Y_j | theta_1) */

    if (lseq < w) continue;             /* skip if sequence too short */

    for (j=0; j<=lseq-w; j++) {         /* site start */
      char *r = res + j;
      int p = 0;
      for (k=0; k<w; k++) {             /* position in site */
        p += theta_1[k][(int) (*r++)];
      }
      pY[j] = p;
    }
    for (j=lseq-w+1; j<lseq; j++) pY[j] = 0;    /* impossible positions */
  }
} // get_pY


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
)
{
  int i;
  for (i=0; i<n_samples; i++) { /* sequence */
    SAMPLE *s = samples[i];
    int lseq = s->length;
    int *pY = s->pY[0];                 	/* log p(Y_j | theta_1) both */
    int *pY1 = s->pY[1];                        /* log p(Y_j | theta_1) + strand */
    int *pY2 = s->pY[2];                        /* log p(Y_j | theta_1) - strand */
    char *pYic = s->pYic;                       /* site on - strand */
    int last_j = lseq-w;                        /* last possible site */
    
    if (lseq < w) continue;             /* skip if sequence too short */
    
    int j, k;
    for (j=0, k=last_j; j<=last_j; j++, k--) {  /* site start */
      if (pY2[k] > pY1[j]) {            /* site on - strand */
        pYic[j] = '\1'; pY[j] = pY2[k];
      } else {                          /* site on + strand */
        pYic[j] = '\0'; pY[j] = pY1[j]; 
      }
    } /* site start */
  } /* sequence */
} // Combine strands


/* NOTE: The concept of an "lmotif" is somewhat redundant with the
 * already-existing concept of an "ltheta". Merging these concepts and
 * adapting them to a more "object-oriented" philosophy could lead to clearer
 * code. However, this would require substantial revision of the existing code
 * where "ltheta" was already used (before I started working on MEME).
 */


/**
 * create_lmotif
 *
 * Create a motif (of iteger-converted log letter probabilities) corresponding
 * to the specified seed (in string format). This function "fills in" the array
 * of int pointers, lmotif, so that each pointer points to an array in lmap.
 *
 */
extern void create_lmotif (
  char *seed_str,    ///< The seed as a string
  int lmap[MAXALPH][MAXALPH], ///< Matrix representing a seq-to-theta mapping fn
  int *lmotif[MAXSITE], ///< The resulting lmotif - OUT
  int *mot_width     ///< The width of the motif - OUT
) {
  // Retrieve the integer-encoded representation of the seed:
  char *seed_enc = NULL;
  int seed_len;
  seed_enc = to_e_seed(seed_str, &seed_len);

  // Allocate space for the lmotif data structure:
/*   int **columns = NULL; */
/*   Resize(columns, seed_len, int *); */
/*   int col_idx; */
/*   for (col_idx = 0; col_idx < seed_len; col_idx++) { */
/*     int *curr_col; */
/*     Resize(curr_col, ALPH, int); */
/*     columns[col_idx] = curr_col; */
/*   } */

  // Convert the encoded seed to an lmotif, thus filling in the memory
  // allocated for the "lmotif" data structure:
  init_lmotif(seed_len, seed_enc, lmotif, lmap);

  *mot_width = seed_len;

  myfree(seed_enc);
} // create_lmotif


/**
 * init_lmotif
 *
 * Convert an integer-encoded subsequence to a motif.
 *
 * Memory must be pre-allocated for theta_1.
 *
 */
extern void init_lmotif (
  int w,             ///< width of site
  char *res,         ///< Integer encoded letters of subsequence
  int *theta_1[MAXSITE], ///< theta_1
  int lmap[MAXALPH][MAXALPH] ///< matrix of frequency vectors
) {
  int m;
  for (m=0; m<w; m++) {
    theta_1[m] = lmap[(int)res[m]];
  }
} // init_lmotif


/**
 * evaluate_seed_DP
 *
 * Evaluate a new seed using dynamic programming, given the differences between
 * it and the seed for which the "pY" array is currently set. This function
 * calculates the value of the LLR_POP objective function for the seed under
 * each of the nsites values in the sp_matrix. The function then creates
 * a new "SEED" object for each of the S_POINTS with the same length as the
 * new seed, and then adds each of the new SEED object to the heap in that
 * S_POINT.
 *
 * NOTE: This function was developed in order to support the dynamic programming
 * performed in branching search. HOWEVER, it can equally be applied to
 * substring search, although this approach has not yet been implemented.
 *
 */
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
) {
  int old_seed_len;
  int new_seed_len;
  get_seed_lengths(s_diffs, &old_seed_len, &new_seed_len);

  // The lmotif corresponding to this seed will be required by
  // next_pY_branching => Obtain it now:
  int *lmotif[MAXSITE];
  int mot_width;
  create_lmotif(new_seed, lmap, lmotif, &mot_width);
  assert(mot_width == new_seed_len);

  // Update the array of pY values given the new seed. This is where the
  // dynamic programming occurs:
  if (!ic) {
    //fprintf(stderr, "NPY....\n");
    next_pY_branching(lmotif, s_diffs, dataset, 0);
    //get_pY(dataset, lmotif, s_diffs->length_new, 0);
  } else {
    next_pY_branching(lmotif, s_diffs, dataset, 1);
    next_pY_branching(lmotif, s_diffs, dataset, 2);
  }

  /* Combine the log Pr(Y_ij | theta_1) values from each strand if using
   * the inverse complement:
   */
  if (ic) {
    combine_strands(dataset->samples, dataset->n_samples, new_seed_len);
  }
  
  // Get a list of "maxima" under the specified sequence model. These maxima
  // will indicate where the pY values are highest:
  int n_maxima = ps(dataset, new_seed_len); // upper bound on # maxima
  P_PROB maxima = (P_PROB) mymalloc(n_maxima * sizeof(p_prob));
  n_maxima = get_max(mtype, dataset, new_seed_len, maxima, ic, TRUE);
  
  // Obtain the S_POINT array of the new seed's length, from the s_point matrix:
  S_POINT *sp_arr = get_sp_arr(sp_mat, new_seed_len);
  int n_nsites = sp_get_num_cols(sp_mat);
  
  /* Evaluate the objective function, LLR_POP, given the alignment of "maxima".
     This will result in a score being entered in each S_POINT from the above
     array... */
  int iseq = -1;     // Dummy
  int ioff = -1;     // Dummy
  char *eseq = NULL; // Dummy
  char *name = NULL; // Dummy
  double col_scores[MAXSITE]; // Dummy
  // Force align_top_subsequences to set the score of every s_point:
  int sp_idx;
  for (sp_idx = 0; sp_idx < n_nsites; sp_idx++) {
    sp_arr[sp_idx].score = LITTLE;
  }

  align_top_subsequences(
    mtype,
    new_seed_len,
    dataset,
    iseq,
    ioff,
    eseq,
    name,
    n_nsites,
    n_maxima,
    maxima,
    col_scores,
    sp_arr
  );

  // Maxima not required anymore => delete them:
  myfree(maxima);
  
  // Update the heap in each of the relevant S_POINTs, adding a new "SEED"
  // object with the appropriate score the the S_POINT's heap in each case:
  update_s_point_heaps(sp_arr, new_seed, n_nsites);

} // evaluate_seed_DP


/**
 * next_pY_branching
 *
 * Next_pY_branching computes log Pr(Y_ij | theta_1) for all i and j in the
 * dataset.
 * It uses dynamic program to perform this computation, reusing the pY results
 * obtained for the previous seed evaluated.
 *
 * This function updates the "pY" arrays by considering the DIFFERENCES of
 * the seed of interest with respect to the seed for which pY is currently
 * evaluated.
 *
 * NOTE: Readability of this function is slightly compromised in order
 * to optimise the speed of the function.
 *
 * NOTE: This function was developed in order to support the dynamic programming
 * performed in branching search. HOWEVER, it can equally be applied to
 * substring search, although this approach has not yet been implemented.
 *
 */
extern void next_pY_branching(
  int *lmotif_new[MAXSITE], ///< The lmotif corresponding to the new seed.
  SEED_DIFFS *s_diffs, ///< The differences between the new and old seeds.
  DATASET *dataset,  ///< The dataset of sequences. Contains the pY arrays.
  int pYindex        ///< Variable indicating which pY array to use.
) {
  // Retrieve the information about the differences that are specifically
  // required by the "guts" of this function:
  int n_diffs = get_n_diffs(s_diffs);
  // Currently this function only caters for 1 difference or 2 differences:
  assert ((n_diffs == 1) || (n_diffs == 2));
  int *diff_idxs = get_diff_idxs(s_diffs);
  int **diff_cols = get_diff_cols(s_diffs);
  int diff_idx1 = 0;
  int diff_idx2 = 0;
  int *diff_col1 = NULL;
  int *diff_col2 = NULL;
  if (n_diffs == 1) {
    diff_idx1 = diff_idxs[0];
    diff_col1 = diff_cols[0];
  } else {
    diff_idx1 = diff_idxs[0];
    diff_col1 = diff_cols[0];
    diff_idx2 = diff_idxs[1];
    diff_col2 = diff_cols[1];
  }
  int new_old_shift = get_seed_shift(s_diffs);
  int old_seed_len, new_seed_len;
  get_seed_lengths(s_diffs, &old_seed_len, &new_seed_len);

  /* Currently, if there is a shift between the new seed and the old seed, then
     ALL the differences between the old and new seed must be accounted for
     by gaps.

     Otherwise this function is unable to compute pY for the new
     seed. I don't yet know how to improve on this, but I doubt it is
     necessary to do so anyway: */
  if (new_old_shift != 0) {
    int n_gaps = abs(new_seed_len + new_old_shift - old_seed_len) +
      abs(new_old_shift);
    assert(n_gaps == n_diffs);
  }

  SAMPLE **samples = dataset->samples;
  int n_samples = dataset->n_samples;

  /* If the new and old seeds are not of equal length and with zero offset, then
     there will be some pY values at the end(s) of each sequence, which must be
     calculated "from scrath". Figure out how many such positions there are: */
  int n_start_pos = MAX(0, new_old_shift);
  int n_end_pos = MAX(0, (old_seed_len - new_seed_len - new_old_shift));

  /* If the new seed is unshifted with respect to the previous seed, then
   * set each pY by simply iterating through dataset and using
   * Pr(Yij|old theta) in order to determine Pr(Yij|new theta) in each case: */
  if (new_old_shift == 0) {
    int seq_idx;
    for (seq_idx=0; seq_idx < n_samples; seq_idx++) { // Sequences in dataset
      SAMPLE *curr_seq = samples[seq_idx]; // The current sequence
      int lseq = curr_seq->length; // Length of current sequence
      // Integer encoding of the sequence:
      char *res = pYindex<2 ? curr_seq->res : curr_seq->resic;
      int *pY = curr_seq->pY[pYindex]; // log p(Y_j | theta_1)
      
      if (lseq < new_seed_len) continue; // Skip if sequence is too short
      
      // Shift the sequence (adjusting according to the indeces at which the
      // current and previous motifs differ) now, rather than inside inner
      // loop...

      // Residues - for the current sequence position - at which the new and old
      // seeds differ:
      char *diff_1 = res+diff_idx1;
      char *diff_2 = NULL;
      if (n_diffs == 2) {
        diff_2 = res+diff_idx2;
      }

      int lett_idx;
      for (lett_idx=0; lett_idx <= (lseq - new_seed_len - n_end_pos);
           lett_idx++) {
        if (n_diffs == 1) {
          /* Only one column differs compared to the previous motif for which pY
           * was evaluated. Evaluate pY (for all i,j in the dataset)
           * appropriately. It is found by adjusting the previous
           * pY value by the amount indicated in the difference column for the
           * letter at that column in the current subsequence...
           */
          pY[lett_idx] = pY[lett_idx] + diff_col1[(int)(*diff_1++)];
        } else {
          /* There are two columns differing between this motif and the previous
           * motif evaluated. Evaluate each pY taking this into consideration.
           * It is found by adjusting the previous pY value by the amounts
           * indicated in the TWO difference columns for the letters at those
           * columns in the current subsequence...
           */
          pY[lett_idx] =
            pY[lett_idx] +
            diff_col1[(int)(*diff_1++)] +
            diff_col2[(int)(*diff_2++)];
        } // examining col_b_idx
      } // offset lett_idx
    } // sequence seq_idx
  } else if (new_old_shift < 0) {
    /* The start of the new seed is shifted "to the left" of the start of the
       old seed => iterate forwards from the beggining of the sequence, and
       then deal with remaining positions (that don't have a valid score from
       the old sequence) explicitly: */
    int seq_idx;
    for (seq_idx=0; seq_idx < n_samples; seq_idx++) {
      SAMPLE *curr_seq = samples[seq_idx]; // The current sequence
      int lseq = curr_seq->length; // Length of current sequence
      // Integer encoding of the sequence:
      char *res = pYindex<2 ? curr_seq->res : curr_seq->resic;
      int *pY = curr_seq->pY[pYindex]; // p(Y_j | theta_1)
      
      if (lseq < new_seed_len) continue; // Skip if sequence is too short
      
      // Shift the sequence (adjusting according to the indeces at which the
      // current and previous motifs differ) now, rather than inside inner
      // loop...

      // Residues - for the current sequence position - at which the new and old
      // seeds differ:
      char *diff_1 = res+diff_idx1;
      char *diff_2 = NULL;
      if (n_diffs == 2) {
        diff_2 = res+diff_idx2;
      }
      int *pY_shifted = pY - new_old_shift;
      
      int lett_idx;
      for(lett_idx=0; lett_idx < (lseq - old_seed_len + new_old_shift);
          lett_idx++) {
        if (n_diffs == 1) {
          /* Only one column differs compared to the previous motif for which pY
           * was evaluated. Evaluate pY (for all i,j in the dataset)
           * appropriately. It is found by adjusting the previous
           * pY value by the amount indicated in the difference column for the
           * letter at that column in the current subsequence...
           */
          // (n_diffs == 1) && (shift < 0) => The difference is the single
          // shift. pY[ij] is set using pY[i(j+1)] under the old seed:
          pY[lett_idx] = pY_shifted[lett_idx] + diff_col1[(int)(*diff_1++)];
        } else {
          /* There are two columns differing between this motif and the previous
           * motif evaluated. Evaluate each pY taking this into consideration.
           * It is found by adjusting the previous pY value by the amounts
           * indicated in the TWO difference columns for the letters at those
           * columns in the current subsequence...
           */
          // Remember we have assumed that if there is a shift, all differences
          // (between seeds) are accounted for by the shift:
          pY[lett_idx] =
            pY_shifted[lett_idx] +
            diff_col1[(int)(*diff_1++)] +
            diff_col2[(int)(*diff_2++)];
        } // examining n_diffs
      } // offset lett_idx
    } // sequence seq_idx
  } else {
    /* The start of the new seed is shifted "to the right" of the start of
       the old seed => iterate backwards from the end of the sequence, and then
       deal with remaining positions (that don't have a valid score from the
       old seed) explicitly: */
    int seq_idx;

    for (seq_idx=0; seq_idx < n_samples; seq_idx++) {
      SAMPLE *curr_seq = samples[seq_idx]; // The current sequence
      int lseq = curr_seq->length; // Length of current sequence
      // Integer encoding of the sequence:
      char *res = pYindex<2 ? curr_seq->res : curr_seq->resic;
      int *pY = curr_seq->pY[pYindex]; // p(Y_j | theta_1)
      
      if (lseq < new_seed_len) continue; // Skip if sequence is too short

      // Shift the sequence (adjusting according to the indeces at which the
      // current and previous motifs differ) now, rather than inside inner
      // loop...

      int first_char_pos = (lseq-new_seed_len-n_end_pos);
      // Residues - for the current sequence position - at which the new and old
      // seeds differ:
      char *diff_1 = res+diff_idx1+first_char_pos;
      char *diff_2 = NULL;
      if (n_diffs == 2) {
        diff_2 = res+diff_idx2+first_char_pos;
      }
      int *pY_shifted = pY - new_old_shift;

      // Scan backwards (because we don't want to "overwrite" the previous
      // pY values prematurely):
      int lett_idx;
      for(lett_idx=first_char_pos; lett_idx >= (n_start_pos);
          lett_idx--) {
        if (n_diffs == 1) {
          // Only 1 column differs compared with the previous motif...
          pY[lett_idx] = pY_shifted[lett_idx] + diff_col1[(int)(*diff_1--)];
        } else {
          pY[lett_idx] =
            pY_shifted[lett_idx] +
            diff_col1[(int)(*diff_1--)] +
            diff_col2[(int)(*diff_2--)];
        } // examining n_diffs
      } // lett_idx
    } // seq_idx
  } // Considering new_old_shift
  
  /* pY has now been calculated for all positions that dynamic programming
     can be used for. However, some positions will remain un-calculated if
     the new and old seeds were not of equal length with zero offset.
     These pY values will be calculated now... */

  // NOTE: Optimisation by pointer arithmetic is NOT carried out below, as
  // the number of iterations in the following loops should be low.

  int seq_idx;
  for (seq_idx=0; seq_idx < n_samples; seq_idx++) {
    SAMPLE *curr_seq = samples[seq_idx]; // The current sequence
    int lseq = curr_seq->length; // Length of current sequence
    // Integer encoding of the sequence:
    char *res = pYindex<2 ? curr_seq->res : curr_seq->resic;
    int *pY = curr_seq->pY[pYindex]; // p(Y_j | theta_1)
    
    if (lseq < new_seed_len) continue; // Skip if sequence is too short

    // Calculate pY values at the start of the sequence:
    int lett_idx;
    for (lett_idx = 0; lett_idx < n_start_pos; lett_idx++) {

      int extra_pY_value = 0;
      int mot_idx; // Index within the lmotif for the new seed.
      // Retrieve the log prob from the letter for each position in the site:
      for (mot_idx = 0; mot_idx < new_seed_len; mot_idx++) {
/*       for (mot_idx = 0; mot_idx < (lett_idx + new_seed_len); mot_idx++) { */
        char curr_lett = res[lett_idx + mot_idx]; // Current letter in seq
        extra_pY_value += lmotif_new[mot_idx][(int)curr_lett];
      }
      pY[lett_idx] = extra_pY_value;
    } // lett_idx
    
    // Calculate pY values at the end of the sequence:
    for (lett_idx = (lseq - new_seed_len - n_end_pos + 1);
         lett_idx <= (lseq - new_seed_len);
         lett_idx++) {

      int extra_pY_value = 0;
      int mot_idx; // Index within the lmotif for the new seed.
      // Retrieve the log prob from the letter for each position in the site:
      for (mot_idx = 0; mot_idx < new_seed_len; mot_idx++) {
        char curr_lett = res[lett_idx + mot_idx]; // Current letter in seq
        extra_pY_value += lmotif_new[mot_idx][(int)curr_lett];
      }
      pY[lett_idx] = extra_pY_value;
    } // lett_idx
  }// seq_idx    

  //trace("Exited npyb...\n");
} // next_pY_branching


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
