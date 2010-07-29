/* data structure to hold command line arguments */
#ifndef GLAM2_ARGS_H
#define GLAM2_ARGS_H

#include <stdio.h>

typedef struct {
  int stop_after;
  int runs;
  char *out_dir;
  int clobber;		/* allow files to be clobbered */
  int two_strands;
  int min_width;
  int max_width;
  int min_seqs;  /* minimum number of sequences in the alignment */
  int init_width;  /* initial width */
  double delete_pseudo;
  double no_delete_pseudo;
  double insert_pseudo;
  double no_insert_pseudo;
  int profile;
  double column_sample_rate;
  int algorithm;  /* 0 = fast, 1 = slow, 2 = FFT */
  double temperature;
  double cool;
  double frozen;  /* lowest allowed temperature, to avoid numerical problems */
  const char *dirichlet_file;
  double bg_pseudo;  /* "weight" of background pseudocounts */
  unsigned seed;  /* seed for random number generator */
  const char *alph_name;  /* "n" = nucleotide, "p" = protein */
  const char *seq_file;
  int quiet;
  const char *embed;
  const char *address;
  char *description;
} args;

void printargs(FILE *fp, int argc, char **argv);

void getargs(args *a, int argc, char **argv);

#endif
