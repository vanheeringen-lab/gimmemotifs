#ifndef GLAM2_ALPHABET_H
#define GLAM2_ALPHABET_H

#include <limits.h>
#include <stdio.h>

typedef struct {
  int size;  /* number of symbols in the alphabet, apart from wildcards */
  int encode[UCHAR_MAX+1];
  int decode[UCHAR_MAX+1];  /* could use unsigned char? */
  double *prob;  /* probability of each symbol */
} alphabet;

/* Set up a nucleotide alphabet */
/* Uniform abundances */
void alphabet_n(alphabet *a);

/* Set up a protein alphabet */
/* Abundances of AB Robinson & LR Robinson 1991: same as BLAST */
void alphabet_p(alphabet *a);

/* Read an alphabet file in vmatch-like format */
/* Unspecified abundances are set to 1 */
void alphabet_read(alphabet *a, FILE *fp);

/* Free the memory allocated in an alphabet */
void alphabet_free(alphabet *a);

#endif
