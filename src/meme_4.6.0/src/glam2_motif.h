/* Data structure representing a motif found by GLAM2 */
/* Should probably replace this with glam2_aln */
#ifndef GLAM2_MOTIF_H
#define GLAM2_MOTIF_H

#include <stdio.h>

typedef struct {
  int width;  /* number of aligned columns */
  int alph_size;
  int seq_num;  /* number of sequences in the alignment */
  int **residue_counts;
  int *delete_counts;
  int *insert_counts;
} motif;

/* Read a motif from a file */
void read_motif(motif *m, int alph_size, const int *encode, FILE *fp);

#endif
