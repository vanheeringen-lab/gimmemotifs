#ifndef READ_TAMO_H
#define READ_TAMO_H

#include "macros.h"
#include "motif.h"

#define TCHUNK 10

/**********************************************************************
  read_tamo

  Reads in a tamo motif file
**********************************************************************/
BOOLEAN read_tamo(
  char*		tamo_file,			  // tamo-file containing motifs IN
  double 	pseudocount,		// pseudocount (in frequency)
  int*      num_motifs,           // Number of motifs retrieved  OUT
  MOTIF_T** motif               // The retrieved motifs OUT
);

#endif
