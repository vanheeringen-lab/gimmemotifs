/****************************************************************************
 * FILE: substmatrix-table.h
 * AUTHOR: William Stafford Noble, Charles E. Grant
 * CREATE DATE: 10/25/2006
 * PROJECT: EVOMCAST
 * DESCRIPTION: Datastructure for storing substitution matrices for
 *              evolutioniary models.
 * COPYRIGHT: 2006, UW
 ****************************************************************************/
#ifndef SUBSTMATRIX_TABLE_H
#define SUBSTMATRIX_TABLE_H

#include "matrix.h"
#include "tree.h"
#include "evomodel.h"

typedef struct substmatrix_array SUBSTMATRIX_ARRAY_T;
typedef struct substmatrix_table SUBSTMATRIX_TABLE_T;

/******************************************************************************
  * This function free the memory associated with a substmatrix_array.
******************************************************************************/
void free_substmatrix_array(SUBSTMATRIX_ARRAY_T *a);

/******************************************************************************
  * This function returns the length of a transistion matrix array.
******************************************************************************/
int get_substmatrix_array_length(SUBSTMATRIX_ARRAY_T *a);

/******************************************************************************
  * This function returns a copy of the matrix for time t of a transition 
  * matrix array. Returns NUL is no matrix is found matching time t.
  * Caller is responsible for free the returned copy.
******************************************************************************/
MATRIX_T* get_substmatrix_for_time(
  SUBSTMATRIX_ARRAY_T *a, 
  double t
);

/******************************************************************************
  * This function returns the time value for the ith element of a transition
  * matrix array.
******************************************************************************/
double get_substmatrix_time(
  SUBSTMATRIX_ARRAY_T *a, 
  int i
);

/******************************************************************************
  * This function allocates an empty substmatrix_array with space allocated
  * for 'length' rows.
******************************************************************************/
SUBSTMATRIX_TABLE_T* allocate_substmatrix_table(int length);

/******************************************************************************
  * This function free the memory associated with a substmatrix_table.
******************************************************************************/
void free_substmatrix_table(SUBSTMATRIX_TABLE_T *table);

/*************************************************************************
  This function builds a table of transistion probability matrices.
  The indexes of the table are the (depth-first) edge number of the tree
  the index of the evolutionary model in the array of models.
 *************************************************************************/
SUBSTMATRIX_TABLE_T* make_substmatrix_table(
  TREE_T* tree, 
  int num_models, 
  EVOMODEL_T** models
);

/******************************************************************************
  * This function returns the length of a transistion matrix table.
******************************************************************************/
int get_substmatrix_table_length(SUBSTMATRIX_TABLE_T *table);

/*************************************************************************
  This function returns a pointer to a copy of the of substmatrix from the
  table correspoinding to the provided time and model number.
  Returns NULL if no such matrix.
  The caller is responsible for freeing the returned substmatrix_array.
 *************************************************************************/
MATRIX_T* get_substmatrix_from_table(
  SUBSTMATRIX_TABLE_T* table,
  int model_index,
  double t
);

/*************************************************************************
  This function returns a copy of the of substmatrix_array correspoinding
  to the provided model number from the provided table.
  The caller is responsible for freeing the returned substmatrix_array.
 *************************************************************************/
SUBSTMATRIX_ARRAY_T* get_substmatrix_array_from_table(
  SUBSTMATRIX_TABLE_T* table,
  int model_index
);

#endif
