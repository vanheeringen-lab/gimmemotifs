/* data structure to hold command line arguments for GLAM2 scanner */
#ifndef GLAM2_SCAN_ARGS_H
#define GLAM2_SCAN_ARGS_H

#include <stdio.h>

typedef struct {
  const char *out_file;
  int hit_num;  /* number of alignments to report */
  int two_strands;
  double delete_pseudo;
  double no_delete_pseudo;
  double insert_pseudo;
  double no_insert_pseudo;
  const char *dirichlet_file;
  const char *alph_name;  /* "n" = nucleotide, "p" = protein */
  const char *motif_file;
  const char *seq_file;
} args;

void printargs(FILE *fp, int argc, char **argv);

void getargs(args *a, int argc, char **argv);

#endif
