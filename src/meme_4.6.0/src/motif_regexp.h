/*
 * motif_regexp.h
 *
 *  Created on: 22/09/2008
 *      Author: rob
 */

#ifndef MOTIF_REGEXP_H_
#define MOTIF_REGEXP_H_

#include "matrix.h"
#include "motif.h"
#include "utils.h"
#include "alphabet.h"

void read_regexp_file(
   char*      filename,          // Name of MEME file  IN
   int*       num_motifs,             // Number of motifs retrieved  OUT
   MOTIF_T*   motifs                 // The retrieved motifs - NOT ALLOCATED!
);

#endif /* MOTIF_REGEXP_H_ */
