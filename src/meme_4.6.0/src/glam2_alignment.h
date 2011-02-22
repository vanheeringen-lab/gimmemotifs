/* Structs and functions for reading glam2 alignments */
#ifndef GLAM2_ALIGNMENT_H
#define GLAM2_ALIGNMENT_H

#include <stdio.h>

typedef struct aligned_seq {
  char *name;  /* sequence name */
  char *seq;  /* the aligned sequence with gaps */
  char *start;  /* start coordinate */
  /* add end coordinate, strand, score if needed */
} aligned_seq;

typedef struct alignment {
  /* could put alignment score, etc. here too */
  size_t seq_num;
  aligned_seq *seqs;
  char *key_positions;
} alignment;

/* Return 1 if all sequences and key_positions have same length, else 0 */
int aln_same_lengths(const alignment *a);

/* Read an alignment from a stream */
void aln_read(alignment *a, FILE *stream);

/* Free all the memory allocated in an alignment */
void aln_free(alignment *a);

#endif
