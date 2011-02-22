#ifndef PRIOR_BLOCK_H
#define PRIOR_BLOCK_H

/******************************************************************************
 * FILE: data-block.h
 * AUTHOR: Charles Grant, Bill Noble, and Tim Bailey
 * CREATION DATE: 2010-09-16
 * COPYRIGHT: 2010 UW
 * 
 * This file contains the public declarations for data structures and functions used 
 * to store and access a block of data. The meaning of a "block of data" is
 * determined by the constructor function, and the implementation of the 
 * functions providing the data_block_reader interface. Current implementations
 * include priors from MEME PSP files and sequence segments from FASTA files.
 *****************************************************************************/

typedef struct data_block DATA_BLOCK_T;

/******************************************************************************
 * This function creates and initializes a data block structure
 * which will contain one prior as a double.
 *****************************************************************************/
DATA_BLOCK_T *new_prior_block();

/******************************************************************************
 * This function creates and initializes a data block structure
 * which will contain a sequence segment of the given length.
 *****************************************************************************/
DATA_BLOCK_T *new_sequence_block(int block_size);

/******************************************************************************
 * This function frees a data block structure
 *****************************************************************************/
void free_data_block(DATA_BLOCK_T *data_block);

/******************************************************************************
 * This function returns the coordinate of the first data element in 
 * this block of data.
 *****************************************************************************/
int get_start_pos_for_data_block(DATA_BLOCK_T *data_block);

/******************************************************************************
 * This function sets the coordinate of the first data in 
 * this block of data.
 *****************************************************************************/
void set_start_pos_for_data_block(DATA_BLOCK_T *data_block, int start_pos);

/******************************************************************************
 * This function returns the number of data items allowed in this data block
 *****************************************************************************/
int get_block_size_from_data_block(DATA_BLOCK_T *data_block);

/******************************************************************************
 * This function returns the number of data items obtained in the last read
 *****************************************************************************/
int get_num_read_into_data_block(DATA_BLOCK_T *data_block);

/******************************************************************************
 * This function sets the number of data items read in the last read
 *****************************************************************************/
void set_num_read_into_data_block(DATA_BLOCK_T *data_block, int num_read);

/******************************************************************************
 * This function returns an array of characters. The size of the array
 * can be found by calling get_num_read_into_data_block()
 *
 * The user should not free this array. It will be freed when the data block
 * is freed.
 *****************************************************************************/
char *get_sequence_from_data_block(DATA_BLOCK_T *data_block);

/******************************************************************************
 * This function returns a pointer to a single prior.
 *
 * The user should not free this array. It will be freed when the data block
 * is freed.
 *****************************************************************************/
double *get_prior_from_data_block(DATA_BLOCK_T *data_block);

#endif
