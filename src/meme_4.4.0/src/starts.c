/*
 * $Id: starts.c 3885 2009-07-16 10:02:35Z tbailey $
 * 
 */

/***********************************************************************
*                                                                      *
*       MEME++                                                         *
*       Copyright 1994, The Regents of the University of California    *
*       Author: Timothy L. Bailey                                      *
*                                                                      *
***********************************************************************/

/* 7-29-99 tlb; change MAX_NSITES to be sites/w */
/* 6-29-99 tlb; add ic to subseq7 */
/*
	Routines to set up an array of starting points for EM.
*/

#include "meme.h"
#include "seed.h"
#include "sp_matrix.h"
#include <assert.h>

/* static int get_nsites0 ( */
/*   int w,             ///< width of motif */
/*   double min_nsites, ///< minimum nsites0 */
/*   double max_nsites, ///< maximum nsites0 */
/*   int heap_size,     ///< size of seed heap */
/*   S_POINT *s_points  ///< array of starting points */
/* ); */


static void global_search (
  SP_MATRIX *sp_matrix,///< Matrix of spoints each with different n and w
                       ///< values
  DATASET *dataset,    ///< Dataset containing (eg) substrings
  MODEL *model         ///< The motif model (nascent)
);


static void substring_search (
  SP_MATRIX *sp_matrix, ///< Matrix of spoints each with different n and w
                        ///< values
  DATASET *dataset,    ///< Dataset containing sequences and parameters
  HASH_TABLE evaluated_seed_ht, ///< A hash table for recording novel seeds
  MODEL *model         ///< The motif model (nascent)
);


static int *make_width_prog (
  int min_width,     ///< Minimum width value in the progression
  int max_width,     ///< Maximum integer value in the progression
  double factor,     ///< Factor specifying rate of increase in progression
  int *n_vals,       ///< The final number of values in the progression - OUT
  BOOLEAN pal,        ///< Indicates whether palindromes are being enforced
  BOOLEAN all_widths  /// all widths between min and max
);


static void use_input_cons(
  SP_MATRIX *sp_matrix, ///< The s_point matrix that will be filled up
  char *e_cons          ///< This seed will be placed at every central s_point
);

/**
 * get_starts
 *
 * Create a list of starting points for EM. 
 *
 * If dataset->prob > 0, samples subsequences in the input dataset for good
 * starting points. Otherwise, e_cons is specified as the starting point.
 *
 * In order to generate starting points for EM, get_starts either performs
 * a "global search", or else uses a specified consensus string. Global search
 * always includes substring search, but several other algorithms ("branching"
 * searches) can optionally be run after substring search in order to improve
 * the starting points from which EM shall be run.
 *
 * Returns number of starting points in list.
 */
extern S_POINT *get_starts(
  DATASET *dataset,  ///< the dataset
  MODEL *model,      ///< the model
  char *e_cons,      ///< encoded consensus sequence
  int *n_starts      ///< number of starting points
)
{
  S_POINT *s_points = NULL;     // array of starting points to perform EM from
  BOOLEAN pal = dataset->pal;	/* force DNA palindromes */
  int min_w = model->min_w;     // minimum width
  int max_w = model->max_w;     // maximum width
  double min_nsites = model->min_nsites;        // min. allowed sites
  double max_nsites = model->max_nsites;        // max. allowed sites
  double sample_prob = dataset->prob; // sampling probability for subsequences

  // Generate geometric progressions over width and nsites. Make these
  // progressions identical to those used in original meme:
  int n_ns;
  int *central_ns = make_geometric_prog(min_nsites, max_nsites, 2, &n_ns);

  int n_ws;
  int *central_ws =
    make_width_prog(min_w, max_w, sqrt(2), &n_ws, pal, model->all_widths);

  // Create a matrix of S_POINTS over which global search will be performed.
  // Each S_POINT corresponds to a separate combination of n and w:
  SP_MATRIX *sp_matrix =
    create_spoint_matrix(central_ws, n_ws, central_ns, n_ns, dataset);

  // Fill in the starting points with the subsequence to use
  if (!e_cons && sample_prob != 0) { // sample subsequences (and do branching)
    // Perform "global search" in order to find starting points for EM.
    // global_search() will establish a string, from which EM will be
    // performed, for a number of S_POINT objects in the matrix (one string per
    // S_POINT). This process results in S_POINTS in the matrix being filled
    // in:
    global_search(sp_matrix, dataset, model);
  } else { // don't sample subsequences
    // Use the specified consensus sequence, "e_cons", in order to fill in
    // the sp_matrix:
    use_input_cons(sp_matrix, e_cons);
  }

  // Choose (from the matrix) which S_POINTS em will be executed from:
  s_points = choose_em_starts(sp_matrix, dataset, model, n_starts);

  // Free memory associated with the sp_matrix, as we no longer need it:
  free_sp_matrix(sp_matrix);

  return s_points;
} // get_starts


/**
 * make_width_prog
 *
 * Generate a list of central motif widths for global search, rounding values
 * to integer (using a cast to int) when they are not integers. The sequence
 * is identical that originally used by MEME: a geometric progression starting
 * from min_width, except that it is terminated by a specified max_width, and
 * the minimum increment between widths is 2.
 * 
 * If palindromes are being enforced, then width "w+1" will be included in
 * addition to every width "w".
 *
 * \return the array of integers in the sequence
  */
static int *make_width_prog (
  int min_width,     ///< Minimum width value in the progression
  int max_width,     ///< Maximum integer value in the progression
  double factor,     ///< Factor specifying rate of increase in progression
  int *n_vals,       ///< The final number of values in the progression - OUT
  BOOLEAN pal,        ///< Indicates whether palindromes are being enforced 
  BOOLEAN all_widths  /// all widths between min and max
) {
  assert (max_width >= min_width); // Precondition

  int curr_width;
  int *progression = NULL;
  *n_vals = 0;

  for (curr_width = min_width; curr_width < max_width;
       curr_width = all_widths? curr_width+1 :
	 (MAX(curr_width+2, (int)(curr_width*factor)))) {
    (*n_vals)++; // Adding an extra value to the progression.
    Resize(progression, *n_vals, int);
    progression[*n_vals - 1] = curr_width;

    // If palindromes are being enforced, then make sure adjacent widths are
    // "central", so that they will be considered for starting points (we
    // don't need to do this is all widths are being considered):
    if (pal && !all_widths) {
      (*n_vals)++;
      Resize(progression, *n_vals, int);
      progression[*n_vals - 1] = curr_width+1;
    }
  }
  // Add the terminal value:
  (*n_vals)++;
  Resize(progression, *n_vals, int);
  progression[*n_vals - 1] = max_width;

  // If palindromes are being enforced, then add an extra starting point at
  // the end of the array (this time we do need to do this even if all
  // widths are considered):
  if (pal) {
    (*n_vals)++;
    Resize(progression, *n_vals, int);
    progression[*n_vals - 1] = max_width+1;
  }

  assert (progression != NULL); // Postcondition

  return progression;
}


/**
 * global_search
 *
 * Perform substring search, followed by optional branching search. Modifies
 * sp_matrix as a result.
 *
 */
static void global_search (
  SP_MATRIX *sp_matrix,///< Matrix of spoints each with different n and w
                       ///< values
  DATASET *dataset,    ///< Dataset containing (eg) substrings
  MODEL *model         ///< The motif model (nascent)
) {
  // Allocate memory for a hash table recording which seeds have already
  // been evaluated...
  HASH_TABLE evaluated_seed_ht = NULL;

  /* Maximum number of entries in hash table (after branching search has
     also been performed) is:
     (Number of seeds evaluated in substring search) +
     (Number of seeds evaluated during branching search).

     Calculate rough upper bound on the number of novel seeds that could be
     found by substring search and branching search...
   */

  /* NOTE: The same branching factor is currently used for width-branching and
    x-branching alike. This was done because w-branching is not a default
    parameter, and so its performance has not been optimised. */
  int bfactor = dataset->branch_params->bfactor;
  int n_central_w = get_num_central_widths(sp_matrix);
  int n_nsites = sp_get_num_cols(sp_matrix);
  // Get the number of letters that can be branched to moves tried in point
  // branching:
  int point_br_alen;
  if (dataset->branch_params->point_branch == NORMAL) {
    point_br_alen = dataset->alength;
  } else if (dataset->branch_params->point_branch == ALL) {
    point_br_alen = dataset->alength + 1;
  } else if (dataset->branch_params->point_branch == NO_POINT_B) {
    point_br_alen = 0;
  } else {
    point_br_alen = 1; // X only branching
  }
  // Calculate max novel seeds generated for each seed branched from:
  double avg_w =
    (((double)get_max_width(sp_matrix) + get_min_width(sp_matrix))/2);
  int seeds_per_branch = (avg_w*point_br_alen);
  if (dataset->branch_params->w_branch) {
    seeds_per_branch = seeds_per_branch + 2*dataset->alength+2;
  }

  // Calculate approximate upper bound on number of novel seeds produced:
  int max_novel =
    (dataset->total_res*n_central_w) +
    3*dataset->main_hs*bfactor*n_central_w*n_nsites*seeds_per_branch;
  /* Placing an upper limit on the amount of memory used by the hashtable bins.
   * Note that a dynamically-resized hashtable could also be used to deal with
   * memory usage concerns. This did not appear necessary in our studies. */
  int n_bins = MIN(max_novel, MAX_BINS);
  evaluated_seed_ht = hash_create(n_bins);

  // Perform substring search:
  substring_search(sp_matrix, dataset, evaluated_seed_ht, model);

  // Use a new BRANCH_PARAMS object to control branching:
  BRANCH_PARAMS branch_params;

  if (dataset->branch_params->w_branch == TRUE) {
    // User requested width-branching => Set up branching parameters accordingly
    // and commence branching:
    branch_params.point_branch = NO_POINT_B;
    branch_params.w_branch = TRUE;
    branch_params.bfactor = dataset->branch_params->bfactor;

    branching_search(&branch_params, model, dataset,
                     sp_matrix, evaluated_seed_ht);
    transfer_final_scores(sp_matrix);
  }

  if (dataset->branch_params->point_branch != NO_POINT_B) {
    // User requested point-wise branching => Set up branching parameters
    // accordingly and commence branching:
    branch_params.point_branch = dataset->branch_params->point_branch;
    branch_params.w_branch = FALSE;
    branch_params.bfactor = dataset->branch_params->bfactor;

    branching_search(&branch_params, model, dataset,
                     sp_matrix, evaluated_seed_ht);
    transfer_final_scores(sp_matrix);
  }

  // Branching has finished => free the memory associated with the hash table.
  hash_destroy(evaluated_seed_ht);

  /* At this stage, the top of each of the heaps of the s_points contains the
     (hopefully) best starting point for the respective combination of w and
     nsites0. Set the attributes in each of the s_point objects accordingly:
  */
  transfer_final_scores(sp_matrix);
} // global_search


/**
 * use_intput_cons
 * 
 * Use the specified consensus string as the EM seed, for all nsites values.
 * 
 */
static void use_input_cons(
  SP_MATRIX *sp_matrix, ///< The s_point matrix that will be filled up
  char *e_cons          ///< This seed will be placed at every central s_point
)
{
  int w_idx, n_idx, lett_idx;

  // get the number of central widths and nsites
  int n_widths = get_num_central_widths(sp_matrix);
  int n_nsites = get_num_central_nsites(sp_matrix);

  // get the central widths
  int *central_ws = get_central_ws(sp_matrix);

  // get the central nsites 
  int *central_ns = get_central_ns(sp_matrix);

  // only fill the central spoints
  for (w_idx = 0; w_idx < n_widths; w_idx++){
    // get the current central width
    int curr_w = central_ws[w_idx];
    for (n_idx = 0; n_idx < n_nsites; n_idx++){
      // get the current central nsites 
      int curr_ns = central_ns[n_idx];
      // get the spoint with width=curr_w nsites=curr_ns
      S_POINT *curr_sp = get_s_point(sp_matrix, curr_w, curr_ns);
      // set the spoint encoded consensus
      curr_sp->e_cons0 = e_cons;
      // set the score to an abitrary value
      curr_sp->score = 1.0;
      // set up human readable consensus sequence cons0
      for (lett_idx = 0; lett_idx < curr_w; lett_idx++){
        curr_sp->cons0[lett_idx] = 
          (e_cons ? unhash(e_cons[lett_idx]) : 'x');
      } // lett_idx
      // add the final character to the string
      curr_sp->cons0[lett_idx] = '\0';
    } // n_idx
  } // w_idx
} // use_input_cons


/**
 * substring_search
 *
 * Perform substring search. For each central width in the sp_matrix, consider
 * every substring of that width from the dataset, as a possible EM seed, and
 * place it in the relevant s_point in the starting point matrix.
 *
  */
static void substring_search (
  SP_MATRIX *sp_matrix, ///< Matrix of spoints each with different n and w values
  DATASET *dataset,     ///< Dataset containing sequences and parameters
  HASH_TABLE evaluated_seed_ht, ///< A hash table for recording novel seeds
  MODEL *model          ///< The motif model (nascent)
) {

  /* For each central width in the sp_matrix, perform substring search on the
     appropriate S_POINT array. Finish on curr_w == max_w: */

  int w_idx; // The index of the current "central" width.
  int *central_ws = get_central_ws(sp_matrix);
  int n_nsites = sp_get_num_cols(sp_matrix);
  for (w_idx = 0; w_idx < get_num_central_widths(sp_matrix); w_idx++) {
    int curr_w = central_ws[w_idx];
    
    // Get the row from the matrix that has the current width:
    S_POINT *curr_sp_row = get_sp_arr(sp_matrix, curr_w);

    // Perform substring search under the current width:
    subseq7(model, dataset, curr_w, n_nsites, curr_sp_row, evaluated_seed_ht);

    // set up human readable consensus sequence for starting point
    int sp_idx;
    for (sp_idx=0; sp_idx < n_nsites; sp_idx++) { // consensus
      char *e_cons0 = curr_sp_row[sp_idx].e_cons0;
      int lett_idx;
      for (lett_idx=0; lett_idx<curr_w; lett_idx++) {
        curr_sp_row[sp_idx].cons0[lett_idx] =
          (e_cons0 ? unhash(e_cons0[lett_idx]) : 'x');
      }
      curr_sp_row[sp_idx].cons0[lett_idx] = 0;
      if (PRINT_STARTS) {
        printf("s=%d, score=%6.0f, w0=%3d, nsites0=%5.0f, cons=%s\n",
               sp_idx, curr_sp_row[sp_idx].score, curr_sp_row[sp_idx].w0,
               curr_sp_row[sp_idx].nsites0, curr_sp_row[sp_idx].cons0);
      }
    } // consensus
  } // w_idx
} // substring_search


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
