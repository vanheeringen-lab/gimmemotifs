/**
 * @file branching_search.c
 *
 * This module implements the branching search algorithm.
 *
 * $Id: branching_search.c 2070 2007-09-12 01:49:09Z eredhead $
 *
 */

#include "calculate_p_y.h"
#include "hash_table.h"
#include "heap.h"
#include "macros.h"
#include "matrix.h"
#include "meme.h"
#include "seed.h"
#include "seed_diffs.h"
#include "sp_matrix.h"
#include <assert.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

#define NODE_NO 0
//#define NODE_NO 1


// Local functions...
static void do_branch (
  BRANCH_PARAMS *branch_params, ///< The parameters controlling branching search
  HEAP *branch_heap, ///< The heap of SEEDs. Branching search will be
                     ///< performed around each of these seeds.
  BOOLEAN ic,        ///< Allow motifs on reverse complement strand of DNA too.
  int lmap[MAXALPH][MAXALPH],
                     ///< Log freq x letter map (for converting seed to theta).
  DATASET *dataset,  ///< The dataset of sequences
  MTYPE mtype,       ///< The type of sequence model
  HASH_TABLE evaluated_seed_ht, ///< A hash_table recording all seeds evaluated
                     ///< thus far.
  SP_MATRIX *sp_mat  ///< A matrix referencing heaps that will be updated
                     ///< By this branching
);
static void eval_mutant_seeds (
  BRANCH_PARAMS *branch_params, ///< The parameters controlling branching search
  char *init_str,    ///> The initial seed, which will be "branched" from
  BOOLEAN ic,        ///< Allow motifs on reverse complement strand of DNA too.
  int lmap[MAXALPH][MAXALPH],
                     ///< Log freq x letter map (for converting seed to theta).
  DATASET *dataset,  ///< The dataset of sequences
  MTYPE mtype,       ///< The type of sequence model
  HASH_TABLE evaluated_seed_ht, ///< A hash_table recording all seeds evaluated
                     ///< thus far.
  SP_MATRIX *sp_mat  ///< A matrix referencing heaps that will be updated
                     ///< by this branching
);
static void evaluate_ACGTX_mutants (
  char *init_str,    ///< The initial seed, which will be "branched" from
  BOOLEAN ic,        ///< Allow motifs on reverse complement strand of DNA too.
  int lmap[MAXALPH][MAXALPH],
                     ///< Log freq x letter map (for converting seed to theta).
  DATASET *dataset,  ///< The dataset of sequences
  MTYPE mtype,       ///< The type of sequence model
  HASH_TABLE evaluated_seed_ht, ///< A hash_table recording all seeds evaluated
                     ///< thus far.
  SP_MATRIX *sp_mat  ///< A matrix referencing heaps that will be updated
                     ///< by this branching
);
static void evaluate_width_mutants (
  char *init_str,    ///> The initial seed, which will be "branched" from
  BOOLEAN ic,        ///< Allow motifs on reverse complement strand of DNA too.
  int lmap[MAXALPH][MAXALPH],
                     ///< Log freq x letter map (for converting seed to theta).
  DATASET *dataset,  ///< The dataset of sequences
  MTYPE mtype,       ///< The type of sequence model
  HASH_TABLE evaluated_seed_ht, ///< A hash_table recording all seeds evaluated
                     ///< thus far.
  SP_MATRIX *sp_mat  ///< A matrix referencing heaps that will be updated
                     ///< by this branching
);
static void consider_seed (
  char *candidate_seed,///< The seed being considered for evaluation.
  char *init_str,    ///< The initial seed. The modified seed belongs to
                     ///< a sequence generated from this initial seed.
  HASH_TABLE evaluated_seed_ht, ///< A hash_table recording all seeds evaluated
                     ///< thus far.
  char *pY_str,      ///< The seed to which the pY arrays currently
                     ///< correspond.
  int lmap[MAXALPH][MAXALPH],
                     ///< Log freq x letter map (for converting seed to theta).
  MTYPE mtype,       ///< The type of model.
  BOOLEAN ic,        ///< Whether or not to consider inverse complement
  DATASET *dataset,  ///< Contains the pY arrays
  SP_MATRIX *sp_mat  ///< This matrix will be updated with the seed, if novel.
);


/**
 * branching_search
 *
 * NOTE: New algorithm introduced on 06-11-06.
 *
 * Performs branching search - a component of MEME's global search - starting
 * from seed heaps that have been previously filled with candidate seeds.
 *
 * THE ALGORITHM:
 * For each iteration (from 1 to bfactor):
 *  For each heap of SEED objects in the sp_matrix:
 *    do_branch() from the current heap, updating heaps in the sp_matrix
 *    with the novel seeds
 *
 */
extern void branching_search (
  BRANCH_PARAMS *branch_params, ///< The parameters controlling branching search
  MODEL *model,      ///< The motif model (nascent)
  DATASET *dataset,  ///< The dataset of sequences
  SP_MATRIX *sp_mat, ///< A matrix of S_POINT objects. The heaps in these
                     ///< objects will be branched from and updated with
                     ///< the resulting novel seeds.
  HASH_TABLE evaluated_seed_ht ///< A hashtable recording which seeds have been
                     ///< evaluated at some point previously.
                     ///< Freq x letter map (for converting seed to theta)
) {
  MTYPE mtype = model->mtype;
  int ic = model->invcomp;
  THETA map = dataset->map;     // letter by frequency mapping matrix
  int alength = dataset->alength;

  // Get the consensus letter vs. log frequency matrix:
  int lmap[MAXALPH][MAXALPH];
  convert_to_lmap(map, lmap, alength);

  int branch_iter;
  // Branch from each heap in the matrix once for each iteration:
  for (branch_iter = 1; branch_iter <= branch_params->bfactor; branch_iter++) {

    // add all the seeds to one large heap to use in branching
    HEAP *mega_heap = create_heap_from_sp_matrix(sp_mat);

    // branch from the mega_heap
    do_branch(
                branch_params,
                mega_heap,
                ic,
                lmap,
                dataset,
                mtype,
                evaluated_seed_ht,
                sp_mat
              );
    destroy_heap(mega_heap);

#ifdef PARALLEL
    // do a reduction for each row of the sp_matrix
    int curr_w;

    // get the number of rows in the sp_matrix
    int num_rows = sp_get_num_rows(sp_mat);

    // get the minimum and maximum widths
    int min_width = get_min_width(sp_mat);
    int max_width = get_max_width(sp_mat);
  
    // get the number of s_points in a row of the sp_matrix
    int n_nsites = sp_get_num_cols(sp_mat);

    // do a reduction for each row in the sp_matrix
    for (curr_w = min_width ; curr_w <= max_width; curr_w++){
      // Get the row from the matrix that has the current width:
      S_POINT *curr_sp_row = get_sp_arr(sp_mat, curr_w);
      // do a reduction for the sp_matrix row
      reduce_across_heaps(curr_sp_row, n_nsites);
    }
#endif

  } // Branching iteration

/*
// print the sp_matrix heaps after branching
#ifdef PARALLEL
if (mpMyID() == NODE_NO) {
  int rows_idx, cols_idx;
  for (rows_idx = 0; rows_idx < get_num_rows(sp_mat); rows_idx++) {
    for (cols_idx = 0; cols_idx < get_num_cols(sp_mat); cols_idx++) {
      S_POINT *current_sp = get_spoint(sp_mat, rows_idx, cols_idx);
      HEAP *current_heap = current_sp->seed_heap;
      HEAP *branching_heap = copy_heap(current_heap);
      if (get_num_nodes(branching_heap) > 0){
        print_heap(stdout, branching_heap);
      }
      destroy_heap(branching_heap);
    }
  }
}
#endif
*/

} // branching_search


/**
 * do_branch
 *
 * This function performs a single round of branching search from the seeds
 * in a given heap ("branch_heap"). Any novel seeds generated by this
 * process are placed on the heaps of the spoints with the corresponding
 * motif width.
 *
 */
static void do_branch (
  BRANCH_PARAMS *branch_params, ///< The parameters controlling branching search
  HEAP *branch_heap, ///< The heap of SEEDs. Branching search will be
                     ///< performed around each of these seeds.
  BOOLEAN ic,        ///< Allow motifs on reverse complement strand of DNA too.
  int lmap[MAXALPH][MAXALPH],
                     ///< Log freq x letter map (for converting seed to theta).
  DATASET *dataset,  ///< The dataset of sequences
  MTYPE mtype,       ///< The type of sequence model
  HASH_TABLE evaluated_seed_ht, ///< A hash_table recording all seeds evaluated
                     ///< thus far.
  SP_MATRIX *sp_mat  ///< A matrix referencing heaps that will be updated
                     ///< By this branching
) {

  /* PRECONDITIONS:
   * 1. branch_heap is *not* one of the heaps referenced by sp_mat (although it
   * may be a copy of one).
   */

  // A temporary seed which will be produced by "mutating" init_seed:
  int nseeds = get_num_nodes(branch_heap); // Number of seeds to branch from
  int seed;
#ifdef PARALLEL
  // get the number of nodes
  int num_nodes = mpNodes();
  // get the current node number
  int my_id = mpMyID();
#endif

  // branch from each seed
  for (seed = 1; seed <= nseeds; seed++){
    // get the root seed
    char *init_str = get_str_seed((SEED *)pop_heap_root(branch_heap));
#ifdef PARALLEL    
    // distribute the work over the nodes
    if (seed % num_nodes == my_id){
      // branch from the seed
      eval_mutant_seeds(branch_params, init_str, ic, lmap,
                        dataset, mtype, evaluated_seed_ht, sp_mat);
    }
#else
    // branch from the seed
    eval_mutant_seeds(branch_params, init_str, ic, lmap,
                      dataset, mtype, evaluated_seed_ht, sp_mat);
#endif
  }

/*
// OLD CODE!!
  // Branch from all SEEDs in branch_heap. Root is at index 1:
  int nseeds = get_num_nodes(branch_heap);
  for (seed_idx=1; seed_idx<=nseeds; seed_idx++) {
    // Retrieve the seed from the branch_heap, which will be branched from:
    char *init_str;    // An ascii representation of an initial seed
    init_str = get_str_seed((SEED *)get_node(branch_heap, seed_idx));

    eval_mutant_seeds(branch_params, init_str, ic, lmap,
                      dataset, mtype, evaluated_seed_ht, sp_mat);
  } // Initial seeds: branched from

// END OLD CODE
*/
} // do_branch


/**
 * eval_mutant_seeds
 *
 * "Branch" from an initial seed, generating a series of "mutant" seeds in a
 * specific order. Uses "next_pY_branching" in order to do evaluate successive
 * novel mutants via dynamic programming.
 *
 * The order of seed evaluations is:
 * 1. Evaluate the initial seed using get_pY(). This simplifies the
 *    process.
 * 2. If doing [ACGTX] branching, call "evaluate_ACGTX_mutants()":
 *    - Evaluate all *novel* seeds in the sequence using next_pY_branching.
 *    The first novel seed is evaluated "relative to" the initial seed.
 *    - THEN, if a novel seed was evaluated, use next_pY_branching to revert the
 *      pY arrays to correspond to the initial seed
 * 3. If doing width-branching, then call "evaluate_width_mutants()":
 *    - Evaluate all *novel* seeds in the sequence using next_pY_branching.
 *    The first novel seed is evaluated "relative to" the initial seed.
 *    - Currently evaluate_width_mutants() does not have to revert the pY
 *    arrays back to their original state, as no branching DP will follow the
 *    funciton.
 */
static void eval_mutant_seeds (
  BRANCH_PARAMS *branch_params, ///< The parameters controlling branching search
  char *init_str,    ///> The initial seed, which will be "branched" from
  BOOLEAN ic,        ///< Allow motifs on reverse complement strand of DNA too.
  int lmap[MAXALPH][MAXALPH],
                     ///< Log freq x letter map (for converting seed to theta).
  DATASET *dataset,  ///< The dataset of sequences
  MTYPE mtype,       ///< The type of sequence model
  HASH_TABLE evaluated_seed_ht, ///< A hash_table recording all seeds evaluated
                     ///< thus far.
  SP_MATRIX *sp_mat  ///< A matrix referencing heaps that will be updated
                     ///< by this branching
) {
  // 1. Evaluate the initial seed using get_pY()...

  // Generate an integer-encoded representation of the initial seed, and
  // a "lmotif" representing the seed (by using the mapping matrix):
  // Declare the log motif corresponding to the current seed
  int *lmotif[MAXSITE];
  int mot_width;
  create_lmotif(init_str, lmap, lmotif, &mot_width);
  
  // Set the pY arrays under the initial, unmutated seed:
  if (!ic) {
    get_pY(dataset, lmotif, mot_width, 0);
  } else {
    get_pY(dataset, lmotif, mot_width, 1);
    get_pY(dataset, lmotif, mot_width, 2);
  }
  
  // 2. Generate mutants via [ACGTX]-branching, and evaluate any novel seeds
  // that occur:
  if (branch_params->point_branch != NO_POINT_B) {
    evaluate_ACGTX_mutants(init_str, ic, lmap, dataset, mtype,
                           evaluated_seed_ht, sp_mat);
  }

  // At this point, we can assume that the pY arrays correspond to the initial
  // seed...

  // 3. Generate mutants via width-branching, evaluating any novel seed that
  // occur:
  if (branch_params->w_branch) {
    evaluate_width_mutants(init_str, ic, lmap, dataset, mtype,
                           evaluated_seed_ht, sp_mat);
  }
} // eval_mutant_seeds


/**
 * consider_seed
 *
 * Consider the specified candidate seed. If it is novel, and its width is
 * encompassed by the sp_matrix, then evaluate the seed, updating the pY arrays
 * in the process.
 * 
 */
static void consider_seed (
  char *candidate_seed,///< The seed being considered for evaluation.
  char *init_str,    ///< The initial seed. The modified seed belongs to
                     ///< a sequence generated from this initial seed.
  HASH_TABLE evaluated_seed_ht, ///< A hash_table recording all seeds evaluated
                     ///< thus far.
  char *pY_str,      ///< The seed to which the pY arrays currently
                     ///< correspond.
  int lmap[MAXALPH][MAXALPH],
                     ///< Log freq x letter map (for converting seed to theta).
  MTYPE mtype,       ///< The type of model.
  BOOLEAN ic,        ///< Whether or not to consider inverse complement
  DATASET *dataset,  ///< Contains the pY arrays
  SP_MATRIX *sp_mat  ///< This matrix will be updated with the seed, if novel.
) {

  /* Determine whether the new candidate seed is novel. It is not novel if
     it is the same as init_seed, or if it is already registered in
     evaluated_seed_ht: */
  BOOLEAN seed_is_novel = TRUE;
  if (strcmp(init_str, candidate_seed) == 0) {
    seed_is_novel = FALSE;
  }
  if (hash_lookup_str(candidate_seed, evaluated_seed_ht) != NULL) {
    seed_is_novel = FALSE;
  }

  // Determine whether the length of the new candidate seed is encompassed by
  // sp_mat:
  BOOLEAN length_ok = TRUE;
  if ((strlen(candidate_seed) < get_min_width(sp_mat)) ||
      (strlen(candidate_seed) > get_max_width(sp_mat))) {
    length_ok = FALSE;
  }

  // Evaluate the seed IFF the seed is novel and its length is valid:
  BOOLEAN eval_seed = (seed_is_novel && length_ok);

  if (eval_seed) {
    /* Compare this candidate seed against the seed for which pY is set.
       Generate an object representing the differences between the candidate
       seed and the previous seed: */
    SEED_DIFFS *s_diffs = get_seed_diffs(pY_str, candidate_seed, lmap);
    
    /* Evaluate the candidate seed using DP, given that we know which columns
       differ with respect to the current "pY" seed. This function updates
       the heaps in the sp_matrix for which the width is the length
       of the current seed:
    */
    evaluate_seed_DP(
      candidate_seed,
      s_diffs,
      lmap,
      mtype,
      ic,
      dataset,
      sp_mat
    );
    
    // The candidate_seed is now the "pY" seed, since evaluation has
    // been performed using it:
    strcpy(pY_str, candidate_seed);
    
    // We no longer need to consider the differences between the new seed
    // and the previous pY seed => destroy the "SEED_DIFFS" object:
    free_seed_diffs(s_diffs);

    // Record the fact that the current seed has now been evaluated:
    hash_insert_str(candidate_seed, evaluated_seed_ht);
  } // Evaluating the candidate seed
} // consider_seed


/**
 * evaluate_ACGTX_mutants
 *
 * Evaluate all *novel* seeds generated by performing ACGTX branching from a
 * given initial seed. THEN, if a novel seed was evaluated, use
 * next_pY_branching to revert the pY arrays to correspond to the initial seed
 *
 * PRECONDITIONS:
 * 1. init_str represents the seed under which the pY arrays (in dataset)
 *    are currently set.
 *
 * POSTCONDITIONS:
 * 1. At the end of this function, the pY arrays correspond to the specified
 * initial seed (represented by "init_str").
 *
 */
static void evaluate_ACGTX_mutants (
  char *init_str,    ///> The initial seed, which will be "branched" from
  BOOLEAN ic,        ///< Allow motifs on reverse complement strand of DNA too.
  int lmap[MAXALPH][MAXALPH],
                     ///< Log freq x letter map (for converting seed to theta).
  DATASET *dataset,  ///< The dataset of sequences
  MTYPE mtype,       ///< The type of sequence model
  HASH_TABLE evaluated_seed_ht, ///< A hash_table recording all seeds evaluated
                     ///< thus far.
  SP_MATRIX *sp_mat  ///< A matrix referencing heaps that will be updated
                     ///< by this branching
) {
  // Keep a copy of the string that currently corresponds to the pY arrays:
  char pY_str[MAXSITE];
  strcpy(pY_str, init_str);

  // Set aside memory for a copy of the initial seed, which will be mutated:
  int seed_len = strlen(init_str);
  char modified_seed[MAXSITE];       // String representation
  
  /* Now, calculate the objective function score for each NOVEL seed in the
     series of mutants... The sequence of mutants will be generated by
     sequentially "mutating" each column to the alternative characters. */

  // First and last letters to try mutating to...
  int nrows = dataset->alength;
  int lett_start=nrows;
  if (dataset->branch_params->point_branch == X_ONLY) {
    // X only branching:
    lett_start=nrows;
  } else {
    lett_start=0;
  }
  // If branching search is considering "mutate to X" as a move, then
  // make this happen by setting the final letter differently:
  int lett_end;
  if ((dataset->branch_params->point_branch == ALL) ||
      (dataset->branch_params->point_branch == X_ONLY)) {
    lett_end=nrows;
  } else {
    lett_end=nrows - 1;
  }

  // Try mutating each column:
  int col_idx;
  for (col_idx=0; col_idx<seed_len; col_idx++) {
    // Try each of the alternative letters in place of the letter in the
    // current column:
    int lett;
    for (lett=lett_start; lett<=lett_end; lett++) {
      // Generate a new candidate seed with the letter "lett" at the column
      // "col_idx":

      // Copy the initial seed:
      strcpy(modified_seed, init_str);

      // Set the letter at column "col_idx" to "lett":
      modified_seed[col_idx] = unhash(lett);

      // Consider the new candidate current seed; evaluate IFF novel:
      consider_seed(modified_seed, init_str, evaluated_seed_ht, pY_str, lmap,
                    mtype, ic, dataset, sp_mat);
    } // Considering next letter in current column
  } // Considering next column

  // Use get_seed_diffs and next_pY_branching to set the pY arrays back to how
  // they were at the start of the function. HOWEVER I only need to do this if
  // pY_str is not the same as modified_seed:
  if (strcmp(init_str, pY_str) != 0) {
    SEED_DIFFS *s_diffs = get_seed_diffs(pY_str, init_str, lmap);

    // Update the pY arrays to correspond to the init_str (as they did at the
    // start of this function):
    int *lmotif[MAXSITE];
    int mot_width;
    create_lmotif(init_str, lmap, lmotif, &mot_width);
    assert(mot_width == strlen(init_str));
    if (!ic) {
      next_pY_branching(lmotif, s_diffs, dataset, 0);
    } else {
      next_pY_branching(lmotif, s_diffs, dataset, 1);
      next_pY_branching(lmotif, s_diffs, dataset, 2);
    }
    
    free_seed_diffs(s_diffs);
  }
} // evaluate_ACGTX_mutants


/**
 * evaluate_width_mutants
 *
 * Evaluate all *novel* seeds generated by performing width branching from a
 * given initial seed. The pY arrays are *not* reverted to correspond to
 * init_str after this function has completed.
 *
 * PRECONDITIONS:
 * 1. init_str represents the seed under which the pY arrays (in dataset)
 *    are currently set.
 *
 */
static void evaluate_width_mutants (
  char *init_str,    ///> The initial seed, which will be "branched" from
  BOOLEAN ic,        ///< Allow motifs on reverse complement strand of DNA too.
  int lmap[MAXALPH][MAXALPH],
                     ///< Log freq x letter map (for converting seed to theta).
  DATASET *dataset,  ///< The dataset of sequences
  MTYPE mtype,       ///< The type of sequence model
  HASH_TABLE evaluated_seed_ht, ///< A hash_table recording all seeds evaluated
                     ///< thus far.
  SP_MATRIX *sp_mat  ///< A matrix referencing heaps that will be updated
                     ///< by this branching
) {
  // Keep a copy of the string that currently corresponds to the pY arrays:
  char pY_str[MAXSITE];
  strcpy(pY_str, init_str);

  // Set aside memory for a copy of the initial seed, which will be mutated:
  char modified_seed[MAXSITE];       // String representation
  
  /* Now, calculate the objective function score for each NOVEL seed in the
     series of mutants... */

  // Determine letters to consider when appending at either end of the seed:
  int nrows = dataset->alength;
  int lett_start = 0;
  int lett_end = nrows - 1;

  // Consider mutants generated by appending a letter to the end of the
  // seed...

  int lett;
  strcpy(modified_seed, init_str); // modified_seed is finished below
  for (lett = lett_start; lett <= lett_end; lett++) {
    // Set the letter at the end of the seed (and terminate string with '\0'):
    int curr_end_idx = strlen(init_str);
    modified_seed[curr_end_idx] = unhash(lett);
    modified_seed[curr_end_idx + 1] = '\0';

    // Consider the new candidate current seed; evaluate IFF novel:
    consider_seed(modified_seed, init_str, evaluated_seed_ht, pY_str, lmap,
                  mtype, ic, dataset, sp_mat);
  } // Appending at end.

  // Consider the mutant generated by deleting a character from the end of
  // the initial seed...

  // Copy the initial seed:
  strcpy(modified_seed, init_str);
  
  // Terminate the seed one letter earlier:
  int curr_end_idx = strlen(init_str);
  modified_seed[curr_end_idx - 1] = '\0';

  // Consider the new candidate current seed; evaluate IFF novel:
  consider_seed(modified_seed, init_str, evaluated_seed_ht, pY_str, lmap, mtype,
                ic, dataset, sp_mat);

  // Consider the mutant generated by adding a character at the start of the
  // initial seed...


  // Copy the initial seed into the array section immediately to the right of
  // that letter. Extra letter will be added in loop below.
  strcpy((modified_seed+1), init_str);
  for (lett = lett_start; lett <= lett_end; lett++) {
    // Place the extra letter at the start of the new seed:
    modified_seed[0] = unhash(lett);

    // Consider the new candidate current seed; evaluate IFF novel:
    consider_seed(modified_seed, init_str, evaluated_seed_ht, pY_str, lmap,
                  mtype, ic, dataset, sp_mat);
  } // Appending at start

  // Consider the mutant generated by deleting a character from the start of
  // the initial seed...

  // Need to generate this mutant explicitly; can't make use of strcpy:
  int seed_idx;
  for (seed_idx = 0; seed_idx < (strlen(init_str) - 1); seed_idx++) {
    modified_seed[seed_idx] = init_str[seed_idx + 1];
  }
  modified_seed[strlen(init_str) - 1] = '\0';
    
  // Consider the new candidate current seed; evaluate IFF novel:
  consider_seed(modified_seed, init_str, evaluated_seed_ht, pY_str, lmap, mtype,
                ic, dataset, sp_mat);
} // evaluate_width_mutants        


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
