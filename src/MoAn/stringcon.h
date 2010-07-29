/*---------------------------------------------------------------*
 * PSSM Searcher                                                 *
 *                                                               *
 * By: Jes Frellsen, Ida Moltke and Martin Thiim                 *
 * University of Copenhagen, spring 2006                         *
 *---------------------------------------------------------------*/

#ifndef STRINGCON_H_
#define STRINGCON_H_

#include "hittable.h"

struct HashEntry {
	int firstSequence;
	int numSequences;
};

// Structure for a concatenated strings
struct structConcatenated {
  char   *pStr;
  int    sizeStr;
  int    *offsets;
  int    numOffsets;
  int    hashSegmentLength;
  struct HashEntry* pHashTable;
};

typedef struct structConcatenated *concatenated;

// Prototypes for the functions
concatenated concat(char **pSrc, int *lengths, int numOfStr);
inline int reverseLookup(concatenated data, int pos);
void freeConcatenated(concatenated data);

#endif
