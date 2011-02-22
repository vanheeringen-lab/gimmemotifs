/******************************************************************************
 * FILE: prior-reader-from-psp.h
 * AUTHOR: Charles Grant, Bill Noble, Tim Bailey
 * CREATION DATE: 2010-09-17
 * COPYRIGHT: 2010 UW
 * 
 * This file contains the public declarations for a data block reader UDT for 
 * reading priors from plain text PSP files originally support by MEME.
 *****************************************************************************/

#ifndef PRIOR_READER_FROM_PSP_H
#define PRIOR_READER_FROM_PSP_H

#include "data-block-reader.h"

/******************************************************************************
 * This function creates an instance of a data block reader UDT for reading
 * priors from a MEME PSP file.
 *****************************************************************************/
DATA_BLOCK_READER_T *new_prior_reader_from_psp(const char *filenae);

#endif
