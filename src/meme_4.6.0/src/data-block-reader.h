/******************************************************************************
 * FILE: data-block-reader.h
 * AUTHOR: Charles Grant, Bill Noble, Tim Bailey
 * CREATION DATE: 2010-09-16
 * COPYRIGHT: 2010 UW
 * 
 * This file contains the declarations for the data structures and functions used 
 * for the data_block_reader. This is the interface wrapper for reading position
 * specific priors and sequences from files. Concrete implemenations exist for 
 * reading priors from MEME PSP files (prior-reader-from-psp.c), and for
 * reading sequence segments from FASTA files (seq-reader-from-fasta.c).
 *****************************************************************************/

#ifndef DATA_BLOCK_READER_H
#define DATA_BLOCK_READER_H

#include "data-block.h"
#include "utils.h"

typedef struct data_block_reader DATA_BLOCK_READER_T;

/******************************************************************************
 * This function creates an instanace of the data block reader UDT.
 * It is only intended to be called by the actual implemenations of the data 
 * block reader interface.
 *
 * Returns a pointer to an instance of a data block reader UDT.
 *****************************************************************************/
DATA_BLOCK_READER_T *new_data_block_reader(
  void *data,
  void (*free_data_block_reader)(DATA_BLOCK_READER_T *reader),
  BOOLEAN_T (*close_data_block_reader)(DATA_BLOCK_READER_T *reader),
  BOOLEAN_T (*reset_data_block_reader)(DATA_BLOCK_READER_T *reader),
  BOOLEAN_T (*get_next_block_from_data_block_reader)(
    DATA_BLOCK_READER_T *reader, 
    DATA_BLOCK_T *data_block 
  ),
  BOOLEAN_T (*go_to_next_sequence_data_block_reader)(
    DATA_BLOCK_READER_T *reader
  ),
  BOOLEAN_T (*get_seq_name_from_data_block_reader)(
    DATA_BLOCK_READER_T *reader, 
    char **sequence
  )
);

/******************************************************************************
 * This function frees an instance of the data block reader UDT.
 *****************************************************************************/
void free_data_block_reader(DATA_BLOCK_READER_T *reader);

/******************************************************************************
 * This function closes the data block reader.
 *
 * Returns TRUE if successful, FALSE otherwise.
 *****************************************************************************/
BOOLEAN_T close_data_block_reader(DATA_BLOCK_READER_T *reader);

/******************************************************************************
 * This function resets the reader to the start of the input by calling
 * the implemenation provided when this reader was created.
 * Any existing data blocks created using this reader should be freed,
 * and new data blocks created.
 *
 * Returns TRUE if successful, FALSE otherwize.
 *****************************************************************************/
BOOLEAN_T reset_data_block_reader(DATA_BLOCK_READER_T *reader);

/******************************************************************************
 * This function directs the data block reader to advance to the next
 * sequence in the input. As a side effect, it sets the value of the seq_name for 
 * the reader.
 *
 * Returns TRUE if successful, FALSE otherwise.
 *****************************************************************************/
BOOLEAN_T go_to_next_sequence_in_data_block_reader(
  DATA_BLOCK_READER_T *reader
);

/******************************************************************************
 * This function gets the name of the current sequence from the data block
 * reader. The name of the sequence is passed using the sequence parameter.
 * The caller is responsible for freeing the memory for the sequence name.
 *
 * Returns TRUE if successful, FALSE if there is no current sequence, as 
 * at the start of the file.
 *****************************************************************************/
BOOLEAN_T get_seq_name_from_data_block_reader(
  DATA_BLOCK_READER_T *reader, 
  char **name // OUT
);

/******************************************************************************
 * This function directs the reader to read the next block of data from the
 * input. The results are filled into the data block structure passed in the
 * argument list.
 *
 * Returns TRUE if successful, FALSE otherwise.
 *****************************************************************************/
BOOLEAN_T get_next_block_from_data_block_reader(
  DATA_BLOCK_READER_T *reader,
  DATA_BLOCK_T *data_block // OUT
);

/******************************************************************************
 * This function gets the data member for the data block reader.
 * The data member contains information specific to the concrete implementation
 * of the reader. The data should not be freed by the caller. It will be freed
 * by when close_data_block_reader() is called.
 *
 * Returns a void pointer to the data. The underlying implementation routines
 * will cast it to the appropriate type.
 *****************************************************************************/
void *get_data_block_reader_data(DATA_BLOCK_READER_T *reader);

#endif
