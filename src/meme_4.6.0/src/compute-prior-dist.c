/********************************************************************
 * FILE: compute-prior-dist.c
 * AUTHOR: William Stafford Noble, Charles E. Grant, Timothy Bailey
 * CREATE DATE: 11/03/2010
 * PROJECT: MEME suite
 * COPYRIGHT: 2010 UW
 ********************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prior-reader-from-psp.h"
#include "pssm.h"


/********************************************************************
 * This program reads a MEME PSP file and computes the binned
 * distribution of priors. The distribution is writen to stdout.
 ********************************************************************/
int main(int argc, char *argv[]) {

  char *usage = "compute-prior-dist <num-bins> <psp-file>";

  if (argc != 3) {
    fprintf(stderr, "Usage: %s\n", usage);
    return -1;
  }

  int num_bins = atoi(argv[1]);
  if (num_bins <= 0) {
    fprintf(stderr, "Usage: %s\n", usage);
    return -1;
  }

  const char *filename = argv[2];

  // Read each prior, find max and min of distribution.
  DATA_BLOCK_READER_T *psp_reader = NULL;
  psp_reader = new_prior_reader_from_psp(filename);
  DATA_BLOCK_T *psp_block = new_prior_block();
  double min_prior = 1.0;
  double max_prior = 0.0;

  while (go_to_next_sequence_in_data_block_reader(psp_reader) != FALSE) {
    while (get_next_block_from_data_block_reader(psp_reader, psp_block) != FALSE) {
      double *prior = get_prior_from_data_block(psp_block);
      if (*prior == 0.0L) {
        // Skip priors that are exactly 0.0L
        continue;
      }
      if (*prior > max_prior) {
        max_prior = *prior;
      }
      if (*prior < min_prior) {
        min_prior = *prior;
      }
    }
  }
  // Print min and max
  printf("%6.5f\n", min_prior);
  printf("%6.5f\n", max_prior);

  // Special case if priors are exactly uniform.
  if (min_prior == max_prior) {
    printf("%6.5f\n", 1.0);
    return 0;
  }

  // With min and max in hand we reset the prior reader
  // and re-read all the priors assigning them to bins.

  // Create the array of bins, intialized to 0.
  double *prior_dist = mm_calloc(num_bins, sizeof(double));

  reset_data_block_reader(psp_reader);
  // Since we've reset the reader we need a new data block.
  free_data_block(psp_block);
  psp_block = new_prior_block();
  double scale = num_bins / (max_prior - min_prior);
  double offset = min_prior;
  size_t num_priors = 0;
  size_t dist_index = 0;
  while (go_to_next_sequence_in_data_block_reader(psp_reader) != FALSE) {
    while (get_next_block_from_data_block_reader(psp_reader, psp_block) != FALSE) {
      double *prior = get_prior_from_data_block(psp_block);
      dist_index = raw_to_scaled(*prior, 1, scale, offset);
      ++prior_dist[dist_index];
      ++num_priors;
    }
  }

  for (dist_index = 0; dist_index < num_bins; ++dist_index) {
    // Print normalized bin counts
    prior_dist[dist_index] /= num_priors;
    printf("%6.5f\n", prior_dist[dist_index]);
  }

  return 0;
}
