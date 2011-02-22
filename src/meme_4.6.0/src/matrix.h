/**************************************************************************
 * FILE: matrix.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 2-4-97
 * PROJECT: shared
 * COPYRIGHT: 1999-2008, WSN
 * VERSION: $Revision: 1.3 $
 * DESCRIPTION: Some simple matrix manipulation routines.
 **************************************************************************/
#ifndef MATRIX_H
#define MATRIX_H

#ifdef ARRAY_H
#warning "array.h included before matrix.h"
#warning "Matrix type may be ignored"
#endif

#ifdef IMATRIX
#define IARRAY
#else
#ifdef SMATRIX
#define SARRAY
#else
#ifdef LMATRIX
#define LARRAY
#else
#define NOT_INT
#endif
#endif
#endif

#include "array.h"
#include <stdio.h>

/**************************************************************************
 * Uses floats by default.  Defining IMATRIX, SMATRIX, LMATRIX or
 * DMATRIX changes the type.
 **************************************************************************/

#define MTYPE ATYPE
#define MSCAN ASCAN

/***************************************************************************
 * Define a matrix type.
 ***************************************************************************/
typedef struct matrix_t {
  int       num_rows;
  int       num_cols;
  ARRAY_T** rows;
}MATRIX_T;

/**************************************************************************
 * Allocate a matrix.
 **************************************************************************/
MATRIX_T* allocate_matrix
  (int num_rows,
   int num_columns);

/**************************************************************************
 * Duplicates a matrix.
 **************************************************************************/
MATRIX_T* duplicate_matrix
  (MATRIX_T* matrix);

/**************************************************************************
 * Converts a theta score matrix of MEME to a MATRIX_T matrix
 * Note that a new matrix is allocated that must be freed.
 * theta   : the input matrix
 * w       : width of motif  (= rows of input/output matrix)
 * alength : alphabet length (= cols of input/output matrix)
 **************************************************************************/
MATRIX_T* convert_matrix
  (double **theta, int w, int alength);

/**************************************************************************
 * Grow a matrix by adding one row to it.  Signal an error if the row
 * does not have the same number of columns as the given matrix.
 **************************************************************************/
void grow_matrix
  (ARRAY_T*   one_row,
   MATRIX_T*  matrix);

/**************************************************************************
 * Remove one row or column from a matrix.
 **************************************************************************/
void remove_matrix_row
  (int        row_index,
   MATRIX_T*  matrix);

void remove_matrix_col
  (int        row_index,
   MATRIX_T*  matrix);

/**************************************************************************
 * Basic access routines.
 **************************************************************************/
int get_num_rows
  (MATRIX_T* matrix);

int get_num_cols
  (MATRIX_T* matrix);

ARRAY_T* get_matrix_row
 (int       row,
  MATRIX_T* matrix);

void set_matrix_row
 (int       row,
  ARRAY_T* one_row,
  MATRIX_T* matrix);

#ifdef BOUNDS_CHECK
#define get_matrix_cell(row,col,matrix) \
   get_matrix_cell_defcheck(row,col,matrix)
#define set_matrix_cell(row,col,value,matrix) \
   set_matrix_cell_defcheck(row,col,value,matrix)
#define incr_matrix_cell(row,col,value,matrix) \
   incr_matrix_cell_defcheck(row,col,value,matrix)
#else
#define get_matrix_cell(row,col,matrix) \
   get_array_item(col, ((MATRIX_T*)matrix)->rows[row])
#define set_matrix_cell(row,col,value,matrix) \
   set_array_item(col, value, ((MATRIX_T *)matrix)->rows[row])
#define incr_matrix_cell(row,col,value,matrix) \
   set_array_item(col, get_array_item(col, ((MATRIX_T*)matrix)->rows[row]) + \
                                      value,((MATRIX_T*)matrix)->rows[row])
#endif

MTYPE get_matrix_cell_defcheck
  (int       row,
   int       col,
   MATRIX_T* matrix);

void set_matrix_cell_defcheck
  (int       row,
   int       col,
   MTYPE     value,
   MATRIX_T* matrix);

void incr_matrix_cell_defcheck
  (int       row,
   int       col,
   MTYPE     value,
   MATRIX_T* matrix);

/***********************************************************************
 * Get a column from a matrix.
 *
 * Returns a newly allocated copy of the requested column.
 ***********************************************************************/
ARRAY_T* get_matrix_column
  (int       i_col,
   MATRIX_T* matrix);

/***********************************************************************
 * Set a column in a matrix.
 ***********************************************************************/
void set_matrix_column
  (ARRAY_T*  column,
   int       i_col,
   MATRIX_T* matrix);

/***********************************************************************
 * Turn an array into a matrix.
 ***********************************************************************/
MATRIX_T* array_to_matrix
  (BOOLEAN_T one_row, /* Put the array in one row, or in many. */
   ARRAY_T*  array);

/**************************************************************************
 * Copy a matrix.
 **************************************************************************/
void copy_matrix
  (MATRIX_T* source_matrix,
   MATRIX_T* target_matrix);

/**************************************************************************
 * Initialize all cells of a given matrix to a given value.
 **************************************************************************/
void init_matrix
  (MTYPE      value,
   MATRIX_T*  matrix);

/**************************************************************************
 * Fill a matrix with a given raw matrix of values.
 **************************************************************************/
void fill_matrix
  (MTYPE*     raw_matrix,
   MATRIX_T*  matrix);

/**************************************************************************
 * Compute the sum of two matrices, assuming they have the same
 * dimension.  The sum is stored in the second matrix.
 **************************************************************************/
void sum_matrices
  (MATRIX_T* matrix1,
   MATRIX_T* matrix2);

/**************************************************************************
 * Extract the diagonal from a square matrix and return it in a
 * newly-allocated array.
 **************************************************************************/
ARRAY_T* extract_diagonal
  (MATRIX_T* matrix);

/**************************************************************************
 * Determine whether a given matrix is symmetric.
 **************************************************************************/
BOOLEAN_T is_symmetric
  (BOOLEAN_T verbose,
   MTYPE     slop,
   MATRIX_T* matrix);

/**************************************************************************
 * Create a matrix M such that M[x,y] = (M1[x,y] + M2[y,x]) / 2.
 **************************************************************************/
MATRIX_T* average_across_diagonal
  (MATRIX_T* matrix1,
   MATRIX_T* matrix2);

/**************************************************************************
 * Add a given value to each element in the diagonal of a square matrix.
 **************************************************************************/
void add_to_diagonal
  (MTYPE     value,
   MATRIX_T* matrix);

/***********************************************************************
 * Determine whether two matrices are equal, within a given bound.
 ***********************************************************************/
BOOLEAN_T equal_matrices
  (ATYPE    close_enough,
   MATRIX_T* matrix1,
   MATRIX_T* matrix2);

/**************************************************************************
 * Multiply all items in a matrix by a given scalar.
 **************************************************************************/
void scalar_mult_matrix
  (MTYPE      value,
   MATRIX_T*  matrix);

/**************************************************************************
 * Add a scalar to all items in a matrix.
 **************************************************************************/
void scalar_add_matrix
  (MTYPE      value,
   MATRIX_T*  matrix);

/**************************************************************************
 * Multiply together corresponding values in two matrices of equal
 * dimension.  Store the result in the second matrix.
 **************************************************************************/
void mult_matrix
  (MATRIX_T* matrix1,
   MATRIX_T* matrix2);

/**************************************************************************
 * Convert a given matrix to or from logs.
 **************************************************************************/
void convert_to_from_log_matrix
  (BOOLEAN_T  to_log,
   MATRIX_T*  the_matrix,
   MATRIX_T*  the_log_matrix);

/**************************************************************************
 * Mix two matrices in log space.
 **************************************************************************/
void mix_log_matrices
  (float     mixing,      /* Percent of matrix2 to be retained. */
   MATRIX_T* matrix1,
   MATRIX_T* matrix2);

/**************************************************************************
 * Read a matrix from a file.
 *
 * Each row of the matrix must appear on its own line.
 **************************************************************************/
MATRIX_T* read_matrix
  (FILE * infile);

/**************************************************************************
 * Read a matrix from a file, when the dimensions are known in advance.
 **************************************************************************/
MATRIX_T* read_known_matrix
  (int    num_rows,
   int    num_cols,
   FILE * infile);

/**************************************************************************
 * Print the matrix, optionally with row and column indices.
 **************************************************************************/
void print_matrix
  (MATRIX_T* matrix,        /* The matrix to be printed. */
   int       width,         /* Width of each cell. */
   int       precision,     /* Precision of each cell. */
   BOOLEAN_T print_titles,  /* Include row and column indices? */
   FILE*     outfile);      /* File to which to write. */

/**************************************************************************
 * void free_matrix
 **************************************************************************/
void free_matrix
  (MATRIX_T* matrix);


/***********************************************************************
 * Fill a matrix with random values between 0 and a given number.
 *
 * Assumes that the random number generator is initialized.
 ***********************************************************************/
void randomize_matrix
  (MTYPE     max_value,
   MATRIX_T* matrix);

/***********************************************************************
 * Compute the sum of the elements in a matrix.
 ***********************************************************************/
MTYPE sum_of_matrix
  (MATRIX_T* matrix);

/***********************************************************************
 * Compute the sum of the squares of a matrix.
 ***********************************************************************/
MTYPE sum_of_squares_matrix
  (MATRIX_T* matrix);

/***********************************************************************
 * Compute the sum of the squares of the differences between two
 * matrices.
 ***********************************************************************/
MTYPE sum_of_square_diff_matrices
  (MATRIX_T* matrix1,
   MATRIX_T* matrix2);

/***********************************************************************
 * Subtract the mean from each row or column of a matrix.
 ***********************************************************************/
void zero_mean_matrix_rows
  (MATRIX_T* matrix);

void zero_mean_matrix_cols
  (MATRIX_T* matrix);

/***********************************************************************
 * Divide each matrix row by its standard deviation.
 ***********************************************************************/
void variance_one_matrix_rows
  (MATRIX_T* matrix);

/***********************************************************************
 * Multiply two matrices to get a third.
 ***********************************************************************/
MATRIX_T* matrix_multiply
  (MATRIX_T* matrix1,
   MATRIX_T* matrix2);

/***********************************************************************
 * Normalize the rows of a matrix.
 ***********************************************************************/
void normalize_rows
  (float tolerance,
   MATRIX_T*   matrix);

/***********************************************************************
 * void normalize_matrix
 *
 * Rangarajan et al. (Neural Comp. (8) 1041-1060) mention a proof by
 * Sinkhorn (1964) that a square matrix can be converted to a double
 * stochastic matrix, in which each row and each column sums to 1.0,
 * via an iterative procedure in which rows and columns are normalized
 * alternately.  This program implements that procedure.
 ***********************************************************************/
void normalize_matrix
  (float tolerance,
   MATRIX_T*   matrix);

/*****************************************************************************
 * Extract one margin of a matrix and return it as an array.
 *****************************************************************************/
ARRAY_T* get_matrix_row_sums
  (MATRIX_T* matrix);
ARRAY_T* get_matrix_col_sums
  (MATRIX_T* matrix);

/*****************************************************************************
 * Sort a given matrix by row, according to a given set of sort keys.
 *****************************************************************************/
void sort_matrix_rows
  (BOOLEAN_T reverse_sort,
   ARRAY_T*  keys,
   MATRIX_T* matrix);

/*****************************************************************************
 * Randomly shuffle the rows or columns of a matrix. If the boolean cols is TRUE * shuffle the columns else shuffle the rows
 *****************************************************************************/
void shuffle_matrix_cols
(MATRIX_T* matrix, BOOLEAN_T cols);

/*****************************************************************************
 * Randomly shuffle the entries of a given matrix.
 *****************************************************************************/
void shuffle_matrix
  (MATRIX_T* matrix);

#endif

/***********************************************************************
 * Get a selection of the matrix rows from index_start to index_stop
 * (including both indices).
  ***********************************************************************/
/* Gets a selection of columns from a given matrix */
void get_matrix_rows(
     int,//Start column index
     int,//End column index
     MATRIX_T*,//Input matrix
     MATRIX_T**//Output matrix
     );


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
