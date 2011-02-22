/**
 * @file seed_diffs.h
 *
 * sp_diffs class. Represents the essential differences between an old
 * seed and a new seed given a particular alignment of those two seeds.
 * This concept lies at the heart of the dynamic programming strategy used
 * during seed evaluation in global search. An sp_diffs object is the key
 * input to the "next_pY" function.
 *
 **/

#ifndef SEED_DIFFS_H
#define SEED_DIFFS_H
#include "meme.h"
#ifdef DMALLOC
#include "dmalloc.h"
#endif

/**
 * A seed_diffs type. Stores the differences between an old seed and
 * a new seed given an alignment of the two seeds.
 */
typedef struct seed_diffs {
  int old_new_shift; /* The offset between the start of the old seed and
                        the start of the new seed, in the given alignment
                        of the two seeds. */
  int **diff_cols;   /* An array of "difference columns" representing the
                        differences in the lmotifs of the two seeds under
                        the alignment. */
  int *diff_idxs;    /* An array of indeces indicating where the differences
                        between the seeds occur. Each index is relative to the
                        new seed. */
  int n_diffs;       // Number of letters (in the implicit alignment of the
                     // two seeds) at which the two seeds differ.
  int length_old;    // Length of the "old" seed.
  int length_new;    // Length of the "new" seed.
} SEED_DIFFS;


extern SEED_DIFFS *get_seed_diffs (
  char *old_seed,    ///< String representation of the old seed
  char *new_seed,    ///< String representation of the new seed
  int lmap[MAXALPH][MAXALPH] ///< Log freq x letter map (used for generating
                     ///< "lmotif" difference columns).
);


extern void free_seed_diffs (
  SEED_DIFFS *s_diffs ///< "This" seed_diffs object
);


extern int get_seed_shift(
  SEED_DIFFS *diffs  ///< The seed_diffs object
);


extern int **get_diff_cols(
  SEED_DIFFS *diffs  ///< The seed_diffs object
);


extern int *get_diff_idxs(
  SEED_DIFFS *diffs  ///< The seed_diffs object
);


extern int get_n_diffs(
  SEED_DIFFS *diffs  ///< The seed_diffs object
);


extern void get_seed_lengths(
  SEED_DIFFS *diffs, ///< The seed_diffs object
  int *length_old,   ///< The length of the old seed - OUT
  int *length_new    ///< The length of the new seed - OUT
);


extern void print_seed_diffs (
  SEED_DIFFS *diffs, ///< Seed diffs object
  FILE *out
);

#endif


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

