/****************************************************************************
 * FILE: substmatrix-table.c
 * AUTHOR: William Stafford Noble, Charles E. Grant
 * CREATE DATE: 10/25/2006
 * PROJECT: EVOMCAST
 * DESCRIPTION: Datastructure for storing substitution matrices for
 *              evolutioniary models.
 * COPYRIGHT: 2006, UW
 ****************************************************************************/

#include <assert.h>
#include "substmatrix-table.h"

/******************************************************************************
  * Data structure for storing an array of transtion probability matricies
  * and corresponding times for a single model and tree.
******************************************************************************/
struct substmatrix_array {
  int length;
  double* t;
  MATRIX_T** substmatrix;
};

/******************************************************************************
  * Data structure for storing a table of transtion probability matricies
  * for a single tree and an array of models.
******************************************************************************/
struct substmatrix_table {
  int length;
  SUBSTMATRIX_ARRAY_T** substmatrix_array;
};

/******************************************************************************
  * This function allocates an empty substmatrix_array with space allocated
  * for 'length' pointers to transition matricese and corresponding times.
******************************************************************************/
SUBSTMATRIX_ARRAY_T* allocate_substmatrix_array(int length) {

  SUBSTMATRIX_ARRAY_T* a = 
    (SUBSTMATRIX_ARRAY_T*) mm_malloc(sizeof(SUBSTMATRIX_ARRAY_T));
  a->length = length;
  a->t = (double *) mm_malloc(length * sizeof(double));
  a->substmatrix = (MATRIX_T**) mm_malloc(length * sizeof(MATRIX_T*));
  int i;
  // Initialize all array elments.
  for (i = 0; i < length; i++) {
    a->t[i] = 0.0;
    a->substmatrix[i] = NULL;
  }

  return a;
}

/******************************************************************************
  * This function free the memory associated with a substmatrix_array.
******************************************************************************/
void free_substmatrix_array(SUBSTMATRIX_ARRAY_T *a) {
  if (a != NULL) {
    // Free the transition matrices.
    int i;
    for (i = 0; i < a->length; i++) {
      free_matrix(a->substmatrix[i]);
    }
    // Free the array for transtion matrices
    myfree(a->substmatrix);
    // Free the array for time values
    myfree(a->t);
    // Free the container
    myfree(a);
  }
}

/******************************************************************************
  * This function returns the length of a transistion matrix array.
******************************************************************************/
int get_substmatrix_array_length(SUBSTMATRIX_ARRAY_T *a) {
  assert(a != NULL);
  return a->length;
}

/******************************************************************************
  * This function returns the time value for the ith element of a transition
  * matrix array.
******************************************************************************/
double get_substmatrix_time(
  SUBSTMATRIX_ARRAY_T *a, 
  int i
) {
  assert(a != NULL);
  assert(i < a->length);
  return a->t[i];
}

/******************************************************************************
  * This function sets the time value for the ith element of a transition
  * matrix array.
******************************************************************************/
static void set_substmatrix_time(
  SUBSTMATRIX_ARRAY_T *a, 
  int i,
  double t
) {
  assert(a != NULL);
  assert(i < a->length);
  a->t[i] = t;
}

/******************************************************************************
  * This function returns a copy of the ith matrix of a transition matrix array.
  * Caller is responsible for freeing the returned copy.
******************************************************************************/
static MATRIX_T* get_substmatrix_matrix(
  SUBSTMATRIX_ARRAY_T *a, 
  int i
) {
  assert(a != NULL);
  assert(i < a->length);
  // Create a copy of our matrix
  MATRIX_T* src = a->substmatrix[i];
  MATRIX_T* dst = NULL;
  if (src != NULL) {
    int num_rows = get_num_rows(src);
    int num_cols = get_num_cols(src);
    dst = allocate_matrix(num_rows, num_cols);
    copy_matrix(src, dst);
  }
  return dst;
}

/******************************************************************************
  * This function sets the transistion prob. matrix for the ith element of 
  a transition matrix array to a copy of the provided matrix.
******************************************************************************/
static void set_substmatrix_matrix(
  SUBSTMATRIX_ARRAY_T *a, 
  int i,
  MATRIX_T* src
) {
  assert(a != NULL);
  assert(i < a->length);
  MATRIX_T* dst = NULL;
  if (src != NULL) {
    int num_rows = get_num_rows(src);
    int num_cols = get_num_rows(src);
    dst = allocate_matrix(num_rows, num_cols);
    copy_matrix(src, dst);
  }
  a->substmatrix[i] = dst;
}

/******************************************************************************
  * This function returns a copy of the matrix for time t of a transition 
  * matrix array. Returns NUL is no matrix is found matching time t.
  * Caller is responsible for free the returned copy.
******************************************************************************/
MATRIX_T* get_substmatrix_for_time(
  SUBSTMATRIX_ARRAY_T *a, 
  double t
) {

  assert(a != NULL);
  // Create a copy of our matrix
  int i;
  MATRIX_T* src = NULL;
  for (i = 0; i < a->length; i++) {
    if (a->t[i] == t) {
      src = a->substmatrix[i];
      break;
    }
  }
  MATRIX_T* dst = NULL;
  if (src != NULL) {
    int num_rows = get_num_rows(src);
    int num_cols = get_num_cols(src);
    dst = allocate_matrix(num_rows, num_cols);
    copy_matrix(src, dst);
  }

  return dst;
}

/*************************************************************************
  This function populates a trans_matix_array with pointers to matrices 
  and the corresponding time values indexed by the edge number of the 
  phylogenetic tree.  The edges are numbered in depth-first order.

  The three parameters are an evolutinary model, a phylogentic
  tree, and a pointer to substmatrix_array structure.

  The function returns an integer containing the number of matrices
  added to the substmatrix_array structure.
 *************************************************************************/
static int populate_substmatrix_array(
  EVOMODEL_T* model, // IN
  TREE_T* tree, // IN
  int current_position, // IN
  SUBSTMATRIX_ARRAY_T* array // OUT
) {
  // Recursively descend the tree, depth first
  int num_children = get_num_children(tree);
  if (is_leaf(tree) != TRUE) {
    int c = 0;
    for (c = 0; c < num_children; c++) {
      TREE_T* child = get_nth_child(c, tree);
      double t = get_length(child);
      set_substmatrix_time(array, current_position, t);
      MATRIX_T* prob_matrix = get_model_prob_matrix(t, model);
      set_substmatrix_matrix(array, current_position, prob_matrix);
      free_matrix(prob_matrix);
      current_position = populate_substmatrix_array(
        model, 
        child, 
        current_position + 1,
        array
      );
    }
  }
  return current_position;
}

/******************************************************************************
  * This function copies a substmatrix_array. The target substmatrix_array
  * must have been previously allocated.
******************************************************************************/
static void copy_substmatrix_array(
    SUBSTMATRIX_ARRAY_T *src, 
    SUBSTMATRIX_ARRAY_T *dst
){
  assert(src->length == dst->length);
  int i;
  for (i = 0; i < src->length; i++) {
    dst->t[i] = src->t[i];
    set_substmatrix_matrix(dst, i, src->substmatrix[i]);
  }
}

/******************************************************************************
  * This function allocates an empty substmatrix_table with space allocated
  * for 'length' pointers to substmatrix_arrays.
******************************************************************************/
SUBSTMATRIX_TABLE_T* allocate_substmatrix_table(int length) {

  SUBSTMATRIX_TABLE_T* table = 
    (SUBSTMATRIX_TABLE_T*) mm_malloc(sizeof(SUBSTMATRIX_TABLE_T));
  table->length = length;
  table->substmatrix_array = (SUBSTMATRIX_ARRAY_T**) 
    mm_malloc(length * sizeof(SUBSTMATRIX_ARRAY_T*));
  int i;
  // Initialize all array elments.
  for (i = 0; i < length; i++) {
    table->substmatrix_array[i] = NULL;
  }

  return table;
}

/*************************************************************************
  This function builds a table of transistion probability matrices: 
  The table is indexed by the edge number of the tree and the evo model.
 *************************************************************************/
SUBSTMATRIX_TABLE_T* make_substmatrix_table(
  TREE_T* tree, 
  int num_models, 
  EVOMODEL_T** models
) {
  SUBSTMATRIX_TABLE_T* table = allocate_substmatrix_table(num_models);
  SUBSTMATRIX_ARRAY_T* array = NULL;
  int num_edges = get_num_edges(tree);
  int i;
  for (i = 0; i < num_models; i++) {
    array = allocate_substmatrix_array(num_edges);
    populate_substmatrix_array(models[i], tree, 0, array);
    table->substmatrix_array[i] = array;
  }

  return table;
}

/******************************************************************************
  * This function returns the length of a transistion matrix table.
******************************************************************************/
int get_substmatrix_table_length(SUBSTMATRIX_TABLE_T *table) {
   return table->length;
 }

/******************************************************************************
  * This function free the memory associated with a substmatrix_table.
******************************************************************************/
void free_substmatrix_table(SUBSTMATRIX_TABLE_T *table) {
  if (table != NULL) {
    // Free the transition matrices.
    int i;
    for (i = 0; i < table->length; i++) {
      free_substmatrix_array(table->substmatrix_array[i]);
    }
    // Free the array for substmatrix_array
    myfree(table->substmatrix_array);
    myfree(table);
  }
}

/*************************************************************************
  This function returns a copy of the of substmatrix_array corresponding
  to the provided model number from the table.
  The caller is responsible for freeing the returned substmatrix_array.
 *************************************************************************/
SUBSTMATRIX_ARRAY_T* get_substmatrix_array_from_table(
  SUBSTMATRIX_TABLE_T* table,
  int model_index
) {
    assert(model_index < table->length);
    // Pick out the array
    SUBSTMATRIX_ARRAY_T* src = table->substmatrix_array[model_index];
    // Allocate the copy
    SUBSTMATRIX_ARRAY_T* dst = allocate_substmatrix_array(src->length);
    copy_substmatrix_array(src, dst);

    return dst;
}

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
) {
    assert(model_index < table->length);
    // Pick out the array for the model
    MATRIX_T* matrix = NULL;
    SUBSTMATRIX_ARRAY_T* array = table->substmatrix_array[model_index];
    // Pick out the matrix for the time
    if (array != NULL) {
      matrix = get_substmatrix_for_time(array, t);
    }

    return matrix;
}
