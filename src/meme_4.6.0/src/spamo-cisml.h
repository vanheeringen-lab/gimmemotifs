
#ifndef SPAMO_CISML_H
#define SPAMO_CISML_H

#include "spamo-matches.h"
#include "red-black-tree.h"

/**************************************************************************
 * Reads a CISML file containing scores for the primary motif for each
 * of the sequences and records the best match above the score threshold
 * for each sequence.
 **************************************************************************/
void load_spamo_primary(
  const char *file,           // CISML file to be read
  int margin,                 // edge area to exclude primary
  double score_threshold,     // minimum score considered a hit
  RBTREE_T *sequences,        // sequence names and where matches are stored
  MOTIF_T *primary_motif      // primary motif
);

/**************************************************************************
 * Reads a CISML file containing the scores for secondary motif database
 * and tallys the spacings of the best matches of primary motifs to secondary
 * motifs.
 **************************************************************************/
void load_spamo_secondary(
  const char *file,           // CISML file to be read
  int margin,                 // area around the primary to find the secondary
  double score_threshold,     // minimum score considered a hit
  int bin,                    // bin size for histogram and pvalue calculations
  double sigthresh,           // pvalue significance threshold
  int test_max,               // distance from the primary for pvalue calculation
  int total_tests,            // test count for Bonferroni correction
  BOOLEAN_T output_sequences, // should sequence match data be dumped?
  char *output_directory,     // directory to create sequence match dumps
  RBTREE_T *sequences,        // sequence names and primary match details
  MOTIF_T *primary_motif,     // primary motif
  RBTREE_T *secondary_motifs, // all secondary motifs and where tally is stored
  int db_id                   // id of current motif database
);

#endif
