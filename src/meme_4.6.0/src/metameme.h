/********************************************************************
 * FILE: metameme.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 05-19-98
 * PROJECT: MHMM
 * COPYRIGHT: 2001-2008, WSN
 * DESCRIPTION: Global functions for the Meta-MEME toolkit.
 ********************************************************************/
#ifndef METAMEME_H
#define METAMEME_H

#include <stdio.h>
#include "config.h"
#include "projrel.h"

/*************************************************************************
 * Print the program header, including help info, if necessary.
 *************************************************************************/
void write_header
  (char*    program,     // Which program called this function? */
   char*    contents,    // What is in this file?
   char*    description, // Text description.
   char*    meme_file,   // Motif filename.
   char*    model_file,  // HMM filename.
   char*    seq_file,    // Name of the input sequence file.
   FILE*    outfile);    // Stream to write to.

#endif
