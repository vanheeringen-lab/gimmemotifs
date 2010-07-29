/*************************************************************************
 * FILE: read-mhmm.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 8-12-97
 * PROJECT: MHMM
 * COPYRIGHT: 1997-2008, WSN
 * DESCRIPTION: Functions for reading HMMs in MHMM format.
 *************************************************************************/
#ifndef READ_MHMM_H
#define READ_MHMM_H
#include "mhmm-state.h"

/*************************************************************************
 * void read_mhmm
 *
 * Read an HMM in MHMM format.
 *
 * IN: model_file - 
 * OUT: the_hmm - 
 *************************************************************************/
void read_mhmm(
  char*    filename,   // The name of the file to be read from.
  MHMM_T** the_hmm     /* The HMM (not pre-allocated). */
);

#endif
