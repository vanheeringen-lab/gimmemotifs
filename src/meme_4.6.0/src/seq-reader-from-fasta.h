/******************************************************************************
 * FILE: seq-reader-from-fasta.h
 * AUTHOR: Charles Grant, Bill Noble, Tim Bailey
 * CREATION DATE: 2010-09-17
 * COPYRIGHT: 2010 UW
 * 
 * This file contains the public declarations for the concrete implementation of a 
 * data-block-reader UDT for reading sequence segments from a FASTA file.
 *****************************************************************************/

#ifndef SEQ_READER_FROM_FASTA_H
#define SEQ_READER_FROM_FASTA_H

#include "data-block-reader.h"

/******************************************************************************
 * This function creates an instance of a data block reader UDT for reading
 * sequence segments from a FASTA file.
 *****************************************************************************/
DATA_BLOCK_READER_T *new_seq_reader_from_fasta(const char *filenae);

#endif
