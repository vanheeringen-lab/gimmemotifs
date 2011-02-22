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


/********************************************************************
 * This program reads a MEME PSP file and computes a uniform prior
 * equal to the mean of the input priors.
 * The a PSP file containing the uniform priors is writen to stdout.
 ********************************************************************/
int main(int argc, char *argv[]) {

  char *usage = "compute-uniform-prior [prior-value] <psp file>";

  const char *filename = NULL;
  double uniform_prior = -1.0L;

  if (argc == 2) {
    filename = argv[1];
  }
  else if (argc == 3) {
    uniform_prior = strtod(argv[1], NULL);
    filename = argv[2];
  }
  else {
    fprintf(stderr, "Usage: %s\n", usage);
    return -1;
  }

  char *seq_name = NULL;

  DATA_BLOCK_READER_T *prior_reader 
    = new_prior_reader_from_psp(filename);
  DATA_BLOCK_T *prior_block = new_prior_block();

  if (uniform_prior < 0.0L) {
    // Set uniform_prior to mean of all priors in the file
    int num_priors = 0;
    double sum_priors = 0.0L;
    while (go_to_next_sequence_in_data_block_reader(prior_reader) != FALSE) {
      while (get_next_block_from_data_block_reader(prior_reader, prior_block) != FALSE) {
        num_priors += get_num_read_into_data_block(prior_block);
        sum_priors += *get_prior_from_data_block(prior_block);
      }
    }
    uniform_prior = sum_priors / num_priors;
  }

  reset_data_block_reader(prior_reader);

  // Print the same sequences, but with uniform priors
  while (go_to_next_sequence_in_data_block_reader(prior_reader) != FALSE) {
    get_seq_name_from_data_block_reader(prior_reader, &seq_name);
    printf(">%s\n", seq_name);
    int prior_idx = 0;
    while (get_next_block_from_data_block_reader(prior_reader, prior_block) != FALSE) {
      ++prior_idx;
      // Print uniform prior
      printf("%5.4f", uniform_prior);
      // Print no more then 10 priors on a line
      if (prior_idx % 10 == 0) {
        printf("\n");
      }
      else {
        printf(" ");
      }
    }
    printf("\n");
  }

  return 0;
}
