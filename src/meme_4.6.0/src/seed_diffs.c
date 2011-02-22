/**
 * @file sp_matrix.c
 *
 * sp_matrix class. Represents a matrix of s_points.
 *
 **/

#include "meme.h"
#include "macros.h"
#include "seed_diffs.h"
#include "seed.h"
#include <assert.h>


/* Note: The concept of an "offset" is used here to allow us to represent a
 * gapless alignment of two sequences (ie two seeds). */

static void create_diff_cols (
  char *seed1,     ///< A seed's string representation
  char *seed2,     ///< Another seed's string representation
  int offset,      ///< The number of characters that seed1 is shifted to the
                   ///< right of seed2
  int lmap[MAXALPH][MAXALPH], ///< Used for generating the difference columns
  int *diff_cols[MAXALPH], ///< The resulting lmotif difference columns - OUT
  int *diff_idxs,  ///< The indexes in seed1 where seed2 differs (given offset)
                   ///< - OUT
  int **d_c[MAXALPH],
  ///< The address of the pointer indicating the start of the array. This
  ///< function will set the pointer at this specified address.
  int **d_i,
  ///< The address of the pointer indicating the start of the array. This
  ///< function will set the pointer at this specified address.
  int *n_diffs     ///< The number of differences between seed1 and seed2
                   ///< - OUT
);


/**
 * get_seed_diffs
 *
 * Align a pair of seeds, and then analyse the alignment, returning an
 * object representing the essential differences between the two seeds
 * under that alignment. These differences can then be used in order to
 * update a "pY" array, so that it corresponds to the *new* seed rather than
 * the *old* seed.
 *
 * The following assumptions are made regarding the two seeds:
 * - The seeds differ by TWO letters at the most under the alignment with
 *   the least number of differences.
 * - There are only two valid types of difference:
 *   1. A unit of offset between the old and the new seed
 *   2. A single character differing between the old and the new seed
 * These assumptions are employed when aligning the two seeds.
 *
 * This function dynamically allocates memory for the returned seed_diffs
 * object. The caller is responsible for freeing that memory at an appropriate
 * time.
 *
 * \return A SEED_DIFFS object representing the differences between the two
 * seeds under the implicit alignment determined by this function.
 */
extern SEED_DIFFS *get_seed_diffs (
  char *old_seed,    ///< String representation of the old seed
  char *new_seed,    ///< String representation of the new seed
  int lmap[MAXALPH][MAXALPH] ///< Log freq x letter map (used for generating
                     ///< "lmotif" difference columns).
) {
  /* PRECONDTIONS:
     1. There exists an alignment (without internal gaps) of the two sequences
     such that the number of differences between the two seeds is <= 2 units.
     2. Both seeds are of length >= 1.
  */

  assert(strlen(old_seed) >= 1);
  assert(strlen(new_seed) >= 1);

  // Find out which seed is smaller:
  char *smaller = NULL;
  char *larger = NULL;

  if (strlen(old_seed) < strlen(new_seed)) {
    smaller = old_seed;
    larger = new_seed;
  } else {
    smaller = new_seed;
    larger = old_seed;
  }

  // "Scan" the smaller seed over the larger seed, and find the offset that
  // yields the smallest difference between the two seeds...

  // Determine the boundaries of the scanning...
  int scan_min, scan_max; // Offset of smaller sequence relative to larger seq

  if (strlen(smaller) == strlen(larger)) {
    // Only consider "overhang" at both ends of the alignment if both sequences
    // are the same length:
    scan_min = -1;
    scan_max = 1;
  } else {
    /* Otherwise only allow the longer sequence to "overhang". If the sequences
       were of different length and we allowed the smaller to "overhang" the
       larger, then there could not be <= 2 differences under that alignment:
    */
    scan_min = 0;
    scan_max = strlen(larger) - strlen(smaller);
  }

  // Consider each possible offset of the smaller sequence relative to the
  // larger sequence:
  int scan_idx;
  int best_n_diffs = (int)BIG;
  int best_offset = 0;
  for (scan_idx = scan_min; scan_idx <= scan_max; scan_idx++) {
    /* Calculate the number of differences between the smaller and the
       larger seed under the alignment specified by the current offset.
       Note that the offset indicates the number of characters by which
       the smaller seed is shifted to the right of the larger seed: */
    int curr_n_diffs = get_n_strdiffs(smaller, larger, scan_idx);
    if (curr_n_diffs < best_n_diffs) {
      best_n_diffs = curr_n_diffs;
      best_offset = scan_idx;
    }
  }

  // Assess the precondition that there exists an alignment producing <= 2
  // diffs:
  assert(best_n_diffs <= 2);

  // best_offset needs to be converted if the new seed is the *larger* seed:
  int new_seed_best_off;
  if (new_seed == larger) {
    new_seed_best_off = -best_offset;
  } else {
    new_seed_best_off = best_offset;
  }

  /* Obtain columns representing the differences in the "lmotifs" for the
     two seeds given the optimal offset. Note that this *could* have been
     achieved above by "get_n_diffs", but we do it separately instead here,
     because I want to avoid performing unnecessary dynamic memory allocation,
     which is computationally intensive. */
  int **diff_cols;
  int *diff_idxs;
  int n_diffs;
  create_diff_cols(new_seed, old_seed, new_seed_best_off, lmap,
                   diff_cols, diff_idxs, &diff_cols, &diff_idxs, &n_diffs);

  assert(n_diffs == best_n_diffs);

  // Create the seed_diffs object:
  SEED_DIFFS *seed_diffs_obj = (SEED_DIFFS *) mymalloc(sizeof(SEED_DIFFS));
  seed_diffs_obj->old_new_shift = new_seed_best_off;
  
  seed_diffs_obj->diff_cols = diff_cols;
  seed_diffs_obj->diff_idxs = diff_idxs;
  seed_diffs_obj->n_diffs = n_diffs;
  seed_diffs_obj->length_old = strlen(old_seed);
  seed_diffs_obj->length_new = strlen(new_seed);

  return seed_diffs_obj;
} // get_seed_diffs


/**
 * free_seed_diffs
 *
 * "Destructor" method: Free the memory associated with this seed diffs object.
 */
extern void free_seed_diffs (
  SEED_DIFFS *s_diffs ///< "This" seed_diffs object
) {
  // Destroy the diff_columns:
  int diff_idx;
  int **diff_cols;
  diff_cols = s_diffs->diff_cols;
  for (diff_idx = 0; diff_idx < s_diffs->n_diffs; diff_idx++) {
    int *curr_diff_col = diff_cols[diff_idx];
    // Destroy the current column:
    myfree(curr_diff_col);
  } // diff_idx
  // Destroy the array that references the individual diff columns:
  myfree(diff_cols);

  // Destroy the array containing the diff_idxs:
  int *diff_idxs;
  diff_idxs = s_diffs->diff_idxs;
  myfree(diff_idxs);

  // Destroy the object itself:
  myfree(s_diffs);
} // free_seed_diffs


/**
 * get_seed_shift
 *
 * Get the offset between the start of the new seed and the start of the old
 * seed under an implicit alignment of the two seeds. Positive if the new seed
 * is to the "right" of the old seed, and negative if it is to the "left".
 *
 * \return The offset
 */
extern int get_seed_shift(
  SEED_DIFFS *diffs  ///< The seed_diffs object
) {
  return diffs->old_new_shift;
}


/**
 * get_diff_cols
 *
 * Get an array of the lmotif difference columns for the two seeds.
 * Each column is: (lmotif column from the new seed) - (lmotif column from
 * the old seed), where an lmotif column is the column obtained by performing
 * a sequence-to-theta mapping from a letter to a column of integer-converted
 * logs of probability values.
 *
 * \return The array of lmotif difference columns
 */
extern int **get_diff_cols(
  SEED_DIFFS *diffs  ///< The seed_diffs object
) {
  return diffs->diff_cols;
}


/**
 * get_diff_idxs
 *
 * Get an array of the indexes indicating where differences occur between the
 * two seeds. These are all relative to the *first* letter of the *new*
 * seed. Eg. if the old seed has a character one unit to the left of the new
 * seed, then the index for this difference would be "-1".
 *
 * \return An array of indices.
 */
extern int *get_diff_idxs (
  SEED_DIFFS *diffs  ///< The seed_diffs object
) {
  return diffs->diff_idxs;
}


/**
 * get_n_diffs
 *
 * Get the number of differences between the two seeds.
 *
 * \return Number of differing characters in the implicit alignment.
 */
extern int get_n_diffs (
  SEED_DIFFS *diffs  ///< The seed_diffs object
) {
  return diffs->n_diffs;
}


/**
 * get_seed_lengths
 *
 * Get the lengths of each of the seeds
 *
 */
extern void get_seed_lengths (
  SEED_DIFFS *diffs, ///< The seed_diffs object
  int *length_old,   ///< The length of the old seed - OUT
  int *length_new    ///< The length of the new seed - OUT
) {
  *length_old = diffs->length_old;
  *length_new = diffs->length_new;
}


/**
 * create_diff_cols
 *
 * Given a pair of aligned seeds and a sequence-to-theta mapping function,
 * find each difference between the two seeds and generate an lmotif 
 * "difference" column for each such difference.
 *
 * For each difference index, the lmotif column corresponding to seed2's letter
 * will be subtracted from the lmotif column corresponding to new_seed's letter.
 *
 * This function dynamically allocates memory for diff_cols and diff_idxs.
 * The calling code is responsible for freeing that memory at an appropriate
 * time.
 */
static void create_diff_cols (
  char *new_seed,    ///< A seed's string representation
  char *old_seed,       ///< Another seed's string representation
  int offset,        ///< The number of characters that new_seed is shifted to the
                     ///< right of old_seed
  int lmap[MAXALPH][MAXALPH], ///< Used for generating the difference columns
  int *diff_cols[MAXALPH], ///< The resulting lmotif difference columns - OUT
  int *diff_idxs,    ///< The indexes in new_seed where old_seed differs (given
                     ///< offset) - OUT
  int **d_c[MAXALPH],
  ///< The address of the pointer indicating the start of the array. This
  ///< function will set the pointer at this specified address.
  int **d_i,
  ///< The address of the pointer indicating the start of the array. This
  ///< function will set the pointer at this specified address.
  int *n_diffs       ///< The number of differences between new_seed and old_seed
                     ///< - OUT
) {
  /* We want to work with positive offsets only => Work with "seed_off"
     and "seed_no_off". "pos_offset" will indicate the offset relative to the
     alignment. */
  char *seed_off, *seed_no_off;
  int pos_offset;
  if (offset >= 0) {
    pos_offset = offset;
    seed_off = new_seed;
    seed_no_off = old_seed;
  } else {
    pos_offset = -offset;
    seed_off = old_seed;
    seed_no_off = new_seed;
  }
  diff_cols = NULL;
  diff_idxs = NULL;
  *n_diffs = 0;

  // Get the integer-encoded representation of each of the seeds (offset and
  // non-offset):
  int non_off_len, off_len;
  char *off_enc = to_e_seed(seed_off, &off_len);
  char *non_off_enc = to_e_seed(seed_no_off, &non_off_len);

  // An array of zeros is required for generating "difference columns":
  static int zero_vector[MAXALPH];
  int vect_idx;
  for (vect_idx = 0; vect_idx < MAXALPH; vect_idx++) {
    zero_vector[vect_idx] = 0;
  }

  // Scan the alignment of the offset and non-offset seed, looking for
  // differences:
  int align_idx; // Index in the seed with no offset
  int align_len = MAX(non_off_len, (off_len + pos_offset));
  for (align_idx = 0; align_idx < align_len; align_idx++) {
    int off_idx = align_idx - pos_offset; // Index in the seed with the offset
    int non_off_idx = align_idx; // Index in the seed with no offset
  
    // NOTE: The following code is as ugly as sin. Could be improved.
    if ((off_idx < 0) || (off_idx >= off_len)) {
      // Align idx falls outside seed with offset...
      *n_diffs = (*n_diffs) + 1; // Register a new difference.

      // Allocate memory for the new difference:
      Resize(diff_idxs, *n_diffs, int);
      Resize(diff_cols, *n_diffs, int *);

      // Record the placement of the current difference, RELATIVE TO THE
      // FIRST CHARACTER OF NEW_SEED:
      if (seed_off == new_seed) {
        // calculate index of this difference relative to the new seed:
        diff_idxs[*n_diffs - 1] = align_idx - pos_offset;

        char non_off_echar = non_off_enc[non_off_idx];

        // Allocate memory for the new diff_column:
        int *diff_col = NULL;
        Resize(diff_col, MAXALPH, int);
        vector_subtract(zero_vector, lmap[(int)non_off_echar],
                        diff_col, MAXALPH);
        diff_cols[*n_diffs - 1] = diff_col; // Point to the new diff column
      } else {
        // calculate index of this difference relative to the new seed:
        diff_idxs[*n_diffs - 1] = align_idx;
        char non_off_echar = non_off_enc[non_off_idx];

        int *diff_col = NULL;
        Resize(diff_col, MAXALPH, int);
        vector_subtract(lmap[(int)non_off_echar], zero_vector,
                        diff_col, MAXALPH);
        diff_cols[*n_diffs - 1] = diff_col;
      }
    } else if (non_off_idx >= non_off_len) {
      // Align idx falls outside seed WITHOUT offset...
      *n_diffs = (*n_diffs) + 1; // Register a new difference.

      // Allocate memory for the new difference:
      Resize(diff_idxs, *n_diffs, int);
      Resize(diff_cols, *n_diffs, int *);

      // Record the placement of the current difference, RELATIVE TO THE
      // FIRST CHARACTER OF NEW_SEED:
      char off_echar = off_enc[off_idx]; // Get char from offset seed
      if (seed_off == new_seed) {
        // calculate index of this difference relative to the new seed:
        diff_idxs[*n_diffs - 1] = align_idx - pos_offset;

        int *diff_col = NULL;
        Resize(diff_col, MAXALPH, int);
        vector_subtract(lmap[(int)off_echar], zero_vector,
                        diff_col, MAXALPH);
        diff_cols[*n_diffs - 1] = diff_col;
      } else {
        // calculate index of this difference relative to the new seed:
        diff_idxs[*n_diffs - 1] = align_idx;
        // Allocate memory for the new diff_column:
        int *diff_col = NULL;
        Resize(diff_col, MAXALPH, int);
        vector_subtract(zero_vector, lmap[(int)off_echar],
                        diff_col, MAXALPH);
        diff_cols[*n_diffs - 1] = diff_col; // Point to the new diff column
      }
    } else {
      // Both the seeds (with or without an offset) have a character at
      // the current index...
      assert((off_idx >= 0) && (off_idx < off_len));
      assert((non_off_idx >= 0) && (non_off_idx < non_off_len));

      char non_off_echar = non_off_enc[non_off_idx];
      char off_echar = off_enc[off_idx];

      if (non_off_echar != off_echar) {
        *n_diffs = (*n_diffs) + 1; // Register a new difference.
        
        // Allocate memory for the new difference:
        Resize(diff_idxs, *n_diffs, int);
        Resize(diff_cols, *n_diffs, int *);
        
        if (seed_off == new_seed) {
          // new_seed was seed with no offset.
          diff_idxs[*n_diffs - 1] = align_idx;

          // Allocate memory for the new diff_column:
          int *diff_col = NULL;
          Resize(diff_col, MAXALPH, int);
          vector_subtract(lmap[(int)off_echar], lmap[(int)non_off_echar],
                          diff_col, MAXALPH);
          diff_cols[*n_diffs - 1] = diff_col; // Point to the new diff column
        } else {
          // new_seed was the one with an offset:
          diff_idxs[*n_diffs - 1] = align_idx - offset;

          // Allocate memory for the new diff_column:
          int *diff_col = NULL;
          Resize(diff_col, MAXALPH, int);
          vector_subtract(lmap[(int)non_off_echar], lmap[(int)off_echar],
                          diff_col, MAXALPH);
          diff_cols[*n_diffs - 1] = diff_col; // Point to the new diff column
        }
      } // Difference in the 2 characters
    } // Inspecting current characters in offset and non-offset seeds.
  } // Inspecting each index in the alignment

  // Encodings of the offset and non-offset seeds (their int representations)
  // are no longer needed => Delete them:
  myfree(non_off_enc);
  myfree(off_enc);

  *d_c = diff_cols;
  *d_i = diff_idxs;
} // create_diff_cols


/**
 *
 * Print a representation of this seed_diffs object to the specified file.
 */
extern void print_seed_diffs (
  SEED_DIFFS *diffs, ///< Seed diffs object
  FILE *out
) {
  fprintf(out, "----------SEED DIFFS OBJECT INSTANCE:---------\n");
  fprintf(out, "diffs->old_new_shift = %d\n", diffs->old_new_shift);
  fprintf(out, "diffs->n_diffs = %d\n", diffs->n_diffs);
  fprintf(out, "diffs->length_old = %d\n", diffs->length_old);
  fprintf(out, "diffs->length_new = %d\n", diffs->length_new);
  fprintf(out, "diffs->diff_idxs:");
  int idx;
  for (idx = 0; idx < diffs->n_diffs; idx++) {
    fprintf(out, " %d", (diffs->diff_idxs)[idx]);
  }
  fprintf(out, "\n");
  fprintf(out, "diffs->diff_cols (as rows):\n");
  for (idx = 0; idx < diffs->n_diffs; idx++) {
    fprintf(out, "col %d:", idx);
    int lett_idx;
    int *curr_col = diffs->diff_cols[idx];
    for (lett_idx = 0; lett_idx < MAXALPH; lett_idx++) {
      fprintf(out, " %d", curr_col[lett_idx]);
    }
    fprintf(out, "\n");
  }
  fprintf(out, "----------------------------------------------\n");
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
