/*************************************************************************
 * FILE: write-mhmm.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 8-12-97
 * PROJECT: MHMM
 * COPYRIGHT: 1997-2008, WSN
 * DESCRIPTION: Functions for writing HMMs in MHMM format.
 *************************************************************************/
#ifndef WRITE_MHMM_H
#define WRITE_MHMM_H
#include "utils.h"
#include "mhmm-state.h"

/*************************************************************************
 * void write_mhmm
 *
 * Write out a hidden markov model in MHMM format.
 *
 * IN: model_file - 
 *     the_hmm - the HMM
 * OUT: ---
 * SIDE EFFECT: An MHMM-style HMM is written to model_file.
 *************************************************************************/
void write_mhmm
  (VERBOSE_T local_verbosity,
   MHMM_T*   the_hmm,     /* The HMM. */
   FILE*     model_file);  /* The file to be written to. */

/*************************************************************************
 * void write_mhmm_to_file
 *
 * Write a hidden markov model in MHMM format to file with the given name.
 *
 * IN: filename - the name of the file to be created.
 *     the_hmm - the HMM
 * OUT: ---
 * SIDE EFFECT: An MHMM-style HMM is written to the named file.
 *************************************************************************/
void write_mhmm_to_file(
  MHMM_T* the_hmm, /* The HMM. */
  char* filename   /* The name of the file to be written to. */
);

/*************************************************************************
 * write_mhmm_xml
 *
 * Write out XML describing HMM
 *************************************************************************/
void write_mhmm_xml(
  VERBOSE_T local_verbosity,
  MHMM_T*   the_hmm,    /* The HMM. */
  FILE*     model_file  /* The file to be written to. */
);

/*************************************************************************
 * write_mhmm_xml_to_file
 *
 * Write out XML describing HMM to the named file.
 *************************************************************************/
void write_mhmm_xml_to_file(
  MHMM_T* the_hmm, /* The HMM. */
  char* filename   /* The name of the file to be written to. */
);

#endif
