#ifndef GENDB_H
#define GENDB_H

#include "seq.h"

SEQ_T *gendb(
  FILE *out,                    // Output stream.
  int type,                     // Type of alphabet.
                                // 0: protein w/ambigs
                                // 1: dna w/ambigs
                                // 2: codons
                                // 3: dna w/o ambigs
                                // 4: protein w/o ambigs
  char *bfile,                  // Name of Markov model file.
  int use_order,                // Order of Markov model to use.
  double *f,			// 0-order model; used if bfile==NULL.
  int nseqs,                    // Number of sequences.
  int min,                      // Shortest sequence.
  int max,                      // Longest sequence.
  int seed                      // Random seed.
);


#endif
