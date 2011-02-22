/**************************************************************************
 * FILE: matrix.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 2-4-97
 * PROJECT: shared
 * COPYRIGHT: 1999-2008, WSN
 * VERSION: $Revision: 1.9 $
 * DESCRIPTION: Some simple matrix manipulation routines.
 **************************************************************************/
#include "matrix.h"
#include "array.h"
#include "utils.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/**************************************************************************
 * Allocate a matrix.
 **************************************************************************/
MATRIX_T* allocate_matrix
  (int num_rows,
   int num_cols)
{
  MATRIX_T* new_matrix;
  int i;

  if ((new_matrix = (MATRIX_T *)mm_malloc(sizeof(MATRIX_T))) == NULL) {
    die("Error allocating matrix.\n");
  }

  if ((new_matrix->rows = (ARRAY_T**)mm_malloc(sizeof(ARRAY_T*) * num_rows))
      == NULL) {
    die("Error allocating matrix rows.\n");
  }
  for (i = 0; i < num_rows; i++) {
    new_matrix->rows[i] = allocate_array(num_cols);
  }

  /* Store the number of rows and columns. */
  new_matrix->num_rows = num_rows;
  new_matrix->num_cols = num_cols;

  /* Return the new matrix. */
  return(new_matrix);
}

/**************************************************************************
 * Duplicates a matrix.
 **************************************************************************/
MATRIX_T* duplicate_matrix
  (MATRIX_T* matrix)
{
  MATRIX_T* new_matrix = allocate_matrix(get_num_rows(matrix),get_num_cols(matrix));

  if(new_matrix != NULL)
    copy_matrix(matrix, new_matrix);

  return(new_matrix);
}

/**************************************************************************
 * Converts a theta score matrix of MEME to a MATRIX_T matrix
 * Note that a new matrix is allocated that must be freed.
 * theta   : the input matrix
 * w       : width of motif  (= rows of input/output matrix)
 * alength : alphabet length (= cols of input/output matrix)
 **************************************************************************/
MATRIX_T* convert_matrix
  (double **theta, int w, int alength)
{
  MATRIX_T* new_matrix = allocate_matrix(w,alength);

  if(new_matrix != NULL) {
    int r, c;
    for(r=0; r<w; r++)
      for(c=0; c<alength; c++)
        set_matrix_cell(r,c, theta[r][c], new_matrix);
  }

  return(new_matrix);
}


/**************************************************************************
 * Grow a matrix by adding one row to it.  Signal an error if the row
 * does not have the same number of columns as the given matrix.
 **************************************************************************/
void grow_matrix
  (ARRAY_T*   one_row,
   MATRIX_T*  matrix)
{
  int num_rows;
  int num_cols;

  /* Get the matrix dimensions. */
  num_rows = get_num_rows(matrix);
  num_cols = get_num_cols(matrix);

  /* Check to be sure the dimensions match. */
  if (get_array_length(one_row) != num_cols) {
    die("Adding row of length %d to matrix of width %d.\n",
  get_array_length(one_row), num_cols);
  }

  /* Reallocate the rows. */
  if ((matrix->rows = (ARRAY_T**)mm_realloc(matrix->rows,
             sizeof(ARRAY_T*) * (num_rows + 1)))
      == NULL) {
    die("Error allocating matrix rows.\n");
  }

  /* Allocate the new row. */
  matrix->rows[num_rows] = allocate_array(num_cols);

  /* Copy the data. */
  copy_array(one_row, matrix->rows[num_rows]);

  /* Increase the number of rows. */
  matrix->num_rows = num_rows + 1;
}

/**************************************************************************
 * Remove one row from a matrix.
 **************************************************************************/
void remove_matrix_row
  (int        row_index,
   MATRIX_T*  matrix)
{
  int i_row;
  int num_rows;

  // Copy all the other rows.
  num_rows = get_num_rows(matrix);
  for (i_row = row_index + 1; i_row < num_rows; i_row++) {
    // FIXME: CEG Do we really need to copy the arrays?
    // Can't we just move the row pointers?
    copy_array(matrix->rows[i_row], matrix->rows[i_row-1]);
  }

  // Free the last row.
  free_array(matrix->rows[num_rows - 1]);

  // Reallocate the rows.
  if ((matrix->rows = (ARRAY_T**)mm_realloc(matrix->rows,
             sizeof(ARRAY_T*) * (num_rows - 1)))
      == NULL) {
    die("Error allocating matrix rows.\n");
  }

  // Decrease the number of rows. */
  (matrix->num_rows)--;
}


/**************************************************************************
 * Remove one column from a matrix.
 **************************************************************************/
void remove_matrix_col
  (int        row_index,
   MATRIX_T*  matrix)
{
  int i_row;
  int num_rows;

  num_rows = get_num_rows(matrix);
  for (i_row = 0; i_row < num_rows; i_row++) {
    remove_array_item(row_index, get_matrix_row(i_row, matrix));
  }

  (matrix->num_cols)--;
}



/**************************************************************************
 * Error checking routine called by all access functions to avoid core
 * dump when attempting to access a null pointer.
 **************************************************************************/
static void check_null_matrix
  (MATRIX_T* matrix)
{
#ifdef BOUNDS_CHECK
  if (matrix == NULL) {
    die("Attempted to access a null matrix.\n");
  }
#else
  /* Avoid compiler warning. */
  matrix = NULL;
#endif
}

/**************************************************************************
 * Implement bounds checking.
 **************************************************************************/
static void matrix_row_check
  (int       row,
   MATRIX_T* matrix)
{
  check_null_matrix(matrix);

#ifdef BOUNDS_CHECK
  if (row < 0) {
    die("Invalid matrix row (%d).\n", row);
  } else if (row > get_num_rows(matrix)) {
    die("Matrix row out of bounds (%d > %d).\n",
        row, get_num_rows(matrix));
  }
#else
  /* Avoid compiler warning. */
  matrix += row;
#endif
}


/**************************************************************************
 * Basic access routines.
 **************************************************************************/
int get_num_rows
  (MATRIX_T* matrix)
{
  check_null_matrix(matrix);
  return(matrix->num_rows);
}

int get_num_cols
  (MATRIX_T* matrix)
{
  check_null_matrix(matrix);
  return(matrix->num_cols);
}

ARRAY_T* get_matrix_row
 (int       row,
  MATRIX_T* matrix)
{
  matrix_row_check(row, matrix);
  return(matrix->rows[row]);
}

void set_matrix_row
 (int       row,
  ARRAY_T*  one_row,
  MATRIX_T* matrix)
{
  int num_cols;

  num_cols = get_num_cols(matrix);

  /* Check to be sure the dimensions match. */
  if (get_array_length(one_row) != num_cols) {
    die("Adding row of length %d to matrix of width %d.\n",
  get_array_length(one_row), num_cols);
  }

  matrix_row_check(row, matrix);
  copy_array(one_row, matrix->rows[row]);
  //fill_matrix(raw_array(one_row), matrix->rows[row]);
}

MTYPE get_matrix_cell_defcheck
  (int       row,
   int       col,
   MATRIX_T* matrix)
{
  matrix_row_check(row, matrix);
  return(get_array_item(col, matrix->rows[row]));
}

void set_matrix_cell_defcheck
  (int       row,
   int       col,
   MTYPE     value,
   MATRIX_T* matrix)
{
  matrix_row_check(row, matrix);
  set_array_item(col, value, matrix->rows[row]);
}

void incr_matrix_cell_defcheck
  (int       row,
   int       col,
   MTYPE     value,
   MATRIX_T* matrix)
{
  matrix_row_check(row, matrix);
  set_array_item(col, get_array_item(col, matrix->rows[row]) + value,
     matrix->rows[row]);
}

/***********************************************************************
 * Get a copy of a column from a matrix.
 *
 * Returns a newly allocated copy of the requested column.
 ***********************************************************************/
ARRAY_T* get_matrix_column
  (int       i_col,
   MATRIX_T* matrix)
{
  int i_row;
  int num_rows;
  ARRAY_T* column;

  /* Allocate space for the copy of the column. */
  num_rows = get_num_rows(matrix);
  column = allocate_array(num_rows);

  /* Copy elements of the column into an array. */
  for (i_row = 0; i_row < num_rows; i_row++) {
    set_array_item(i_row, get_matrix_cell(i_row, i_col, matrix), column);
  }

  /* Return the result. */
  return(column);
}

/***********************************************************************
 * Set a column in a matrix.
 ***********************************************************************/
void set_matrix_column
  (ARRAY_T*  column,
   int       i_col,
   MATRIX_T* matrix)
{
  int i_row;
  int num_rows;

  /* Copy elements of the array into the column. */
  num_rows = get_num_rows(matrix);
  for (i_row = 0; i_row < num_rows; i_row++) {
    set_matrix_cell(i_row, i_col, get_array_item(i_row, column), matrix);
  }
}



/***********************************************************************
 * Turn an array into a matrix.
 ***********************************************************************/
MATRIX_T* array_to_matrix
  (BOOLEAN_T one_row, /* Put the array in one row, or in many. */
   ARRAY_T*  array)
{
  MATRIX_T* new_matrix;
  int i_item;
  int num_items = get_array_length(array);

  if (one_row) {
    new_matrix = allocate_matrix(1, num_items);

    /* Do a fast memcpy. */
    memcpy(raw_array(new_matrix->rows[0]), raw_array(array),
     sizeof(ATYPE) * num_items);
  } else {
    new_matrix = allocate_matrix(num_items, 1);

    for (i_item = 0; i_item < num_items; i_item++) {
      set_matrix_cell(i_item, 0, get_array_item(i_item, array), new_matrix);
    }
  }

  return(new_matrix);
}

/**************************************************************************
 * Copy a matrix.  Assumes that the target matrix is allocated.
 **************************************************************************/
void copy_matrix
  (MATRIX_T* source_matrix,
   MATRIX_T* target_matrix)
{
  int num_rows;
  int i_row;

  // Verify that the matrices have the same dimensions.
  num_rows = get_num_rows(source_matrix);
  if (num_rows != get_num_rows(target_matrix)) {
    die("Attempted to copy matrices with different numbers of rows (%d != %d)\n",
      num_rows, get_num_rows(target_matrix));
  }
  myassert(TRUE, get_num_cols(source_matrix) == get_num_cols(target_matrix),
     "Copying matrix with %d columns into matrix with %d columns.\n",
     get_num_cols(source_matrix), get_num_cols(target_matrix));

  for (i_row = 0; i_row < num_rows; i_row++) {
    copy_array(get_matrix_row(i_row, source_matrix),
         get_matrix_row(i_row, target_matrix));
  }
}


/**************************************************************************
 * Initialize all cells of a given matrix to a given value.
 **************************************************************************/
void init_matrix
  (MTYPE      value,
   MATRIX_T*  matrix)
{
  int num_rows;
  int i_row;

  num_rows = get_num_rows(matrix);
  for (i_row = 0; i_row < num_rows; i_row++) {
    init_array(value, get_matrix_row(i_row, matrix));
  }
}

/**************************************************************************
 * Fill a matrix with a given raw matrix of values.  The input is an
 * array of values stored columnwise; i.e.,
 *
 *   array[ n*i+j ] = matrix[i][j].
 *
 **************************************************************************/
void fill_matrix
  (MTYPE*     raw_matrix,
   MATRIX_T*  matrix)
{
  int num_rows;
  int i_row;
  int num_cols;
  int i_col;
  MTYPE value;

  num_rows = get_num_rows(matrix);
  num_cols = get_num_cols(matrix);
  for (i_row = 0; i_row < num_rows; i_row++) {
    for (i_col = 0; i_col < num_cols; i_col++) {
      value = raw_matrix[(i_row * num_rows) + i_col];
      set_matrix_cell(i_col, i_row, value, matrix);
    }
  }
}


/**************************************************************************
 * Compute the sum of two matrices, assuming they have the same
 * dimension.  The sum is stored in the second matrix.
 **************************************************************************/
void sum_matrices
  (MATRIX_T* matrix1,
   MATRIX_T* matrix2)
{
  int num_rows;
  int i_row;

  /* Make sure the two matrices have the same number of rows. */
  num_rows = get_num_rows(matrix1);
  if (num_rows != get_num_rows(matrix2)) {
    die("Attempted to add matrices with different dimensions (%d != %d).",
  num_rows, get_num_rows(matrix2));
  }

  /* Add 'em up. */
  for (i_row = 0; i_row < num_rows; i_row++) {
    sum_array(get_matrix_row(i_row, matrix1), get_matrix_row(i_row, matrix2));
  }
}

/**************************************************************************
 * Extract the diagonal from a square matrix and return it in a
 * newly-allocated array.
 **************************************************************************/
ARRAY_T* extract_diagonal
  (MATRIX_T* matrix)
{
  int i_row;
  int num_rows;
  ARRAY_T* diagonal;

  /* Check to be sure the matrix is square. */
  num_rows = get_num_rows(matrix);
  if (num_rows != get_num_cols(matrix)) {
    die("The given matrix is not square.");
  }

  /* Allocate the diagonal array. */
  diagonal = allocate_array(num_rows);

  /* Extract the diagonal value. */
  for (i_row = 0; i_row < num_rows; i_row++) {
    set_array_item(i_row, get_matrix_cell(i_row, i_row, matrix), diagonal);
  }
  return(diagonal);
}

/**************************************************************************
 * Determine whether a given matrix is symmetric.
 **************************************************************************/
BOOLEAN_T is_symmetric
  (BOOLEAN_T verbose,
   MTYPE     slop,
   MATRIX_T* matrix)
{
  int num_rows;
  int num_cols;
  int i_row;
  int i_col;
  MTYPE upper;
  MTYPE lower;

  /* Get the matrix dimensions. */
  num_rows = get_num_rows(matrix);
  num_cols = get_num_cols(matrix);

  /* Check for symmetric across the diagonal. */
  for (i_row = 0; i_row < num_rows; i_row++) {
    for (i_col = 0; i_col < i_row; i_col++) {
      upper = get_matrix_cell(i_row, i_col, matrix);
      lower = get_matrix_cell(i_col, i_row, matrix);

      if (!almost_equal(upper, lower, slop)) {
  if (verbose) {
    fprintf(stderr, "matrix[%d][%d]=%g matrix[%d][%d]=%g diff=%g\n",
      i_row, i_col, upper, i_col, i_row, lower, upper - lower);
  }
  return(FALSE);
      }
    }
  }
  return(TRUE);
}

/**************************************************************************
 * Create a matrix M such that M[x,y] = (M1[x,y] + M2[y,x]) / 2.
 *
 * The output matrix has the same dimensions as the first input matrix.
 **************************************************************************/
MATRIX_T* average_across_diagonal
  (MATRIX_T* matrix1,
   MATRIX_T* matrix2)
{
  int num_rows;
  int i_row;
  int num_cols;
  int i_col;
  double value;
  MATRIX_T* new_matrix;

  /* Allocate the new matrix. */
  num_rows = get_num_rows(matrix1);
  num_cols = get_num_cols(matrix1);
  new_matrix = allocate_matrix(num_rows, num_cols);

  for (i_row = 0; i_row < num_rows; i_row++) {
    for (i_col = 0; i_col <= i_row; i_col++) {

      /* Compute the average of this and its reflection. */
      value = (get_matrix_cell(i_row, i_col, matrix1) +
         get_matrix_cell(i_col, i_row, matrix2)) / 2.0;

      /* Set the two cells to their average. */
      set_matrix_cell(i_row, i_col, value, new_matrix);
      set_matrix_cell(i_col, i_row, value, new_matrix);
    }
  }

  return(new_matrix);
}

/**************************************************************************
 * Add a given value to each element in the diagonal of a square matrix.
 **************************************************************************/
void add_to_diagonal
  (MTYPE     value,
   MATRIX_T* matrix)
{
  int num_rows;
  int i_row;

  num_rows = get_num_rows(matrix);

  if (num_rows != get_num_cols(matrix)) {
    fprintf(stderr, "Can't reflect a non-square matrix (%d != %d).\n",
      num_rows, get_num_cols(matrix));
  }

  for (i_row = 0; i_row < num_rows; i_row++) {
    incr_matrix_cell(i_row, i_row, value, matrix);
  }
}



/***********************************************************************
 * Determine whether two matrices are equal, within a given bound.
 ***********************************************************************/
BOOLEAN_T equal_matrices
  (ATYPE    close_enough,
   MATRIX_T* matrix1,
   MATRIX_T* matrix2)
{
  int num_rows;
  int i_row;

  /* Verify that the matrices have the same number of rows. */
  num_rows = get_num_rows(matrix1);
  if (num_rows != get_num_rows(matrix2)) {
    die("Attempted to compare matrices with different numbers of rows (%d != %d)\n",
      num_rows, get_num_rows(matrix2));
  }

  /* Compare each pair of corresponding rows. */
  for (i_row = 0; i_row < num_rows; i_row++) {
    if (!equal_arrays(close_enough, get_matrix_row(i_row, matrix1),
          get_matrix_row(i_row, matrix2))) {
      fprintf(stderr, "Matrices differ in row %d.\n", i_row);
      return(FALSE);
    }
  }
  return(TRUE);
}

/**************************************************************************
 * Multiply all items in a matrix by a given scalar.
 **************************************************************************/
void scalar_mult_matrix
  (MTYPE      value,
   MATRIX_T*  matrix)
{
  int i_row;
  int num_rows = get_num_rows(matrix);

  for (i_row = 0; i_row < num_rows; i_row++) {
    scalar_mult(value, get_matrix_row(i_row, matrix));
  }
}

/**************************************************************************
 * Add a scalar to all items in a matrix.
 **************************************************************************/
void scalar_add_matrix
  (MTYPE      value,
   MATRIX_T*  matrix)
{
  int i_row;
  int num_rows = get_num_rows(matrix);

  for (i_row = 0; i_row < num_rows; i_row++) {
    scalar_add(value, get_matrix_row(i_row, matrix));
  }
}

/**************************************************************************
 * Multiply together corresponding values in two matrices of equal
 * dimension.  Store the result in the second matrix.
 **************************************************************************/
void mult_matrix
  (MATRIX_T* matrix1,
   MATRIX_T* matrix2)
{
  int i_row;
  int i_col;
  int num_rows = get_num_rows(matrix1);
  int num_cols = get_num_cols(matrix1);

  for (i_row = 0; i_row < num_rows; i_row++) {
    for (i_col = 0; i_col < num_cols; i_col++) {
      set_matrix_cell(i_row, i_col, get_matrix_cell(i_row, i_col, matrix1) *
          get_matrix_cell(i_row, i_col, matrix2), matrix2);
    }
  }
}

/**************************************************************************
 * Convert a given matrix to or from logs.
 **************************************************************************/
void convert_to_from_log_matrix
  (BOOLEAN_T  to_log,
   MATRIX_T*  source_matrix,
   MATRIX_T*  target_matrix)
{
  int num_rows;
  int i_row;

  num_rows = get_num_rows(source_matrix);
  for (i_row = 0; i_row < num_rows; i_row++) {
    convert_to_from_log_array(to_log,
            get_matrix_row(i_row, source_matrix),
            get_matrix_row(i_row, target_matrix));
  }
}

/**************************************************************************
 * Mix two matrices in log space.
 **************************************************************************/
void mix_log_matrices
  (float     mixing,      /* Percent of matrix2 to be retained. */
   MATRIX_T* matrix1,
   MATRIX_T* matrix2)
{
  int i_row;
  int num_rows;

  num_rows = get_num_rows(matrix1);
  for (i_row = 0; i_row < num_rows; i_row++) {
    mix_log_arrays(mixing, get_matrix_row(i_row, matrix1),
       get_matrix_row(i_row, matrix2));
  }
}

/**************************************************************************
 * Read a matrix from a file.
 **************************************************************************/
#define MAX_ROW 100000
MATRIX_T* read_matrix
  (FILE * infile)
{
  MATRIX_T* matrix;             /* The matrix to be read. */
  char      first_row[MAX_ROW]; /* The first row of the matrix. */
  char      one_row[MAX_ROW];   /* One row of the matrix. */
  char *    fgets_result;       /* Error indicator for 'fgets' function. */
  int       i_row;              /* Index of the current row. */
  int       i_column;           /* Index of the current column. */
  int       num_columns;        /* Total number of columns. */
  int       num_scanned;
  char*     string_ptr = NULL;
  MTYPE     one_value;          /* One value read from the file. */

  if (infile == NULL) {
    die("Attempted to read matrix from null file.");
  }

  /* Allocate the matrix. */
  matrix = (MATRIX_T*)mm_malloc(sizeof(MATRIX_T));

  /* Read the first row. */
  fgets_result = fgets(one_row, MAX_ROW, infile);
  if (one_row[strlen(one_row) - 1] != '\n') {
    die("Matrix lines too long.  Increase MAX_ROW.");
  }

  /* Count the entries in the row. */
  strcpy(first_row, one_row);
  num_columns = 0;
  for (string_ptr = strtok(first_row, " \t"); string_ptr != NULL;
       string_ptr = strtok(NULL, " \t")) {
    if (strcmp(string_ptr, "\n") != 0) {
      num_columns++;
    }
  }
  matrix->num_cols = num_columns;

  /* Allocate the first dimension of the matrix. */
  matrix->rows = (ARRAY_T**)mm_malloc(MAX_ROW * sizeof(ARRAY_T*));

  /* Read the matrix. */
  for (i_row = 0; fgets_result != NULL; i_row++) {

    /* Allocate this row. */
    matrix->rows[i_row] = allocate_array(num_columns);

    /* Read the first value. */
    string_ptr = strtok(one_row, " \t");

    /* Read the row. */
    for (i_column = 0; i_column < num_columns; i_column++) {

      /* Make sure we read the number properly. */
      if (string_ptr == NULL) {
  die("Error reading matrix at position (%d,%d). ", i_row, i_column);
      }
      num_scanned = sscanf(string_ptr, MSCAN, &one_value);
      if ((num_scanned == 0) || (num_scanned == EOF)) {
  die("Error reading matrix at position (%d,%d).", i_row, i_column);
      }

      /* Record the current number of rows. */
      matrix->num_rows = i_row + 1;

      /* Store the value. */
      set_matrix_cell(i_row, i_column, one_value, matrix);

      /* Read the next value. */
      string_ptr = strtok(NULL, " \t");
    }

    /* Read the next line. */
    fgets_result = fgets(one_row, MAX_ROW, infile);
  }

  if (verbosity > NORMAL_VERBOSE) {
    fprintf(stderr, "Read %d x %d matrix.\n", matrix->num_rows,
      matrix->num_cols);
  }

  return(matrix);
}

/**************************************************************************
 * Read a matrix from a file, when the dimensions are known in advance.
 **************************************************************************/
MATRIX_T* read_known_matrix
  (int    num_rows,
   int    num_cols,
   FILE * infile)
{
  MATRIX_T* matrix;      /* The matrix to be read. */
  int       scan_result; /* Error indicator for 'fgets' function. */
  int       i_row;       /* Index of the current row. */
  int       i_col;       /* Index of the current column. */
  MTYPE     value;       /* One value read from the file. */

  if (infile == NULL) {
    die("Attempted to read matrix from null file.");
  }

  /* Allocate the matrix. */
  matrix = allocate_matrix(num_rows, num_cols);

  /* Read the matrix. */
  for (i_row = 0; i_row < num_rows; i_row++) {
    for (i_col = 0; i_col < num_cols; i_col++) {
      scan_result = fscanf(infile, MSCAN, &value);
      if (scan_result != 1) {
  die("Error reading matrix at row %d, column %d (%s).\n", i_row, i_col,
      value);
      }
      set_matrix_cell(i_row, i_col, value, matrix);
    }
  }

  if (verbosity > NORMAL_VERBOSE) {
    fprintf(stderr, "Read %d x %d matrix.\n", matrix->num_rows,
      matrix->num_cols);
  }

  return(matrix);
}

/**************************************************************************
 * Print the matrix, optionally with row and column indices.
 **************************************************************************/
void print_matrix
  (MATRIX_T* matrix,        /* The matrix to be printed. */
   int       width,         /* Width of each cell. */
   int       precision,     /* Precision of each cell. */
   BOOLEAN_T print_titles,  /* Include row and column indices? */
   FILE*     outfile)       /* File to which to write. */
{
  int num_rows = get_num_rows(matrix);
  int num_cols = get_num_cols(matrix);
  int i_row;
  int i_col;
  ARRAY_T* column;
  if (matrix == NULL) {
    return;
  }

  SWAP(int,  num_rows, num_cols);

  /* Print the header row. */
  if (print_titles) {
    fprintf(outfile, "   ");
    for (i_col = 0; i_col < num_cols; i_col++) {
      fprintf(outfile, "%*d ", width, i_col);
    }
    fprintf(outfile, "\n");
  }

  /* Print the matrix values. */
  for (i_row = 0; i_row < num_rows; i_row++) {
    if (print_titles) {
      fprintf(outfile, "%2d ", i_row);
    }
    column = get_matrix_column(i_row, matrix);
    print_array(column, width, precision, TRUE, outfile);
    free_array(column);
  }
}

/**************************************************************************
 * Free a given matrix.
 **************************************************************************/
void free_matrix
  (MATRIX_T* matrix)
{
  int i_row;
  int num_rows;

  if (matrix == NULL) {
    return;
  }

  num_rows = get_num_rows(matrix);
  for (i_row = 0; i_row < num_rows; i_row++) {
    free_array(get_matrix_row(i_row, matrix));
  }
  myfree(matrix->rows);
  myfree(matrix);
}

/***********************************************************************
 * Fill a matrix with random values between 0 and a given number.
 *
 * Assumes that the random number generator is initialized.
 ***********************************************************************/
void randomize_matrix
  (MTYPE     max_value,
   MATRIX_T* matrix)
{
  int num_rows;
  int i_row;

  num_rows = get_num_rows(matrix);
  for (i_row = 0; i_row < num_rows; i_row++) {
    randomize_array(max_value, get_matrix_row(i_row, matrix));
  }
}

/***********************************************************************
 * Compute the sum of the elements in a matrix.
 ***********************************************************************/
MTYPE sum_of_matrix
  (MATRIX_T* matrix)
{
  int num_rows = get_num_rows(matrix);
  int i_row;
  MTYPE return_value;

  /* Add up the sum of squared errors for each row. */
  return_value = (MTYPE)0;
  for (i_row = 0; i_row < num_rows; i_row++) {
    return_value += array_total(get_matrix_row(i_row, matrix));
  }
  return(return_value);
}


/***********************************************************************
 * Compute the sum of the squares of a matrix.
 ***********************************************************************/
MTYPE sum_of_squares_matrix
  (MATRIX_T* matrix)
{
  int num_rows = get_num_rows(matrix);
  int i_row;
  MTYPE return_value;

  /* Add up the sum of squared errors for each row. */
  return_value = (MTYPE)0;
  for (i_row = 0; i_row < num_rows; i_row++) {
    return_value += sum_of_squares(get_matrix_row(i_row, matrix));
  }
  return(return_value);
}

/***********************************************************************
 * Compute the sum of the squares of the differences between two arrays.
 ***********************************************************************/
MTYPE sum_of_square_diff_matrices
  (MATRIX_T* matrix1,
   MATRIX_T* matrix2)
{
  int num_rows = get_num_rows(matrix1);
  int i_row;
  MTYPE return_value;

  /* Check to be sure the matrices are the same size. */
  /* FIX: Make generic function for doing this. */

  /* Add up the sum of squared errors for each row. */
  return_value = (MTYPE)0;
  for (i_row = 0; i_row < num_rows; i_row++) {
    return_value += sum_of_square_diffs(get_matrix_row(i_row, matrix1),
          get_matrix_row(i_row, matrix2));
  }
  return(return_value);
}

/***********************************************************************
 * Normalize the rows of a matrix.
 ***********************************************************************/
void normalize_rows
  (float tolerance,
   MATRIX_T*   matrix)
{
  int       num_rows = get_num_rows(matrix);
  int       i_row;

  for (i_row = 0; i_row < num_rows; i_row++) {
    normalize(tolerance, get_matrix_row(i_row, matrix));
  }
}


/***********************************************************************
 * Subtract the mean from each row or column of a matrix.
 ***********************************************************************/
void zero_mean_matrix_rows
  (MATRIX_T* matrix)
{
  int num_rows = get_num_rows(matrix);
  int i_row;

  for (i_row = 0; i_row < num_rows; i_row++) {
    sum_to_zero(get_matrix_row(i_row, matrix));
  }
}

void zero_mean_matrix_cols
  (MATRIX_T* matrix)
{
  int num_cols = get_num_cols(matrix);
  int i_col;
  ARRAY_T* this_column;

  for (i_col = 0; i_col < num_cols; i_col++) {
    this_column = get_matrix_column(i_col, matrix);
    sum_to_zero(this_column);
    set_matrix_column(this_column, i_col, matrix);
    free_array(this_column);
  }
}

/***********************************************************************
 * Divide each matrix row by its standard deviation.
 ***********************************************************************/
void variance_one_matrix_rows
  (MATRIX_T* matrix)
{
  int num_rows = get_num_rows(matrix);
  int i_row;

  for (i_row = 0; i_row < num_rows; i_row++) {
    variance_one_array(get_matrix_row(i_row, matrix));
  }
}

/***********************************************************************
 * Iteratively normalize the rows and columns in a matrix.
 ***********************************************************************/
void normalize_matrix
  (float tolerance,
   MATRIX_T*   matrix)
{
  ARRAY_T*  one_column; /* Storage for one column of the matrix. */
  ARRAY_T*  this_row;   /* The current row of the matrix. */
  int       iter;
  int       num_rows;
  int       num_cols;
  int       i_row;
  int       i_col;
  BOOLEAN_T keep_going; /* Are we tolerably close to a normalized matrix? */
  MTYPE     total;

  /* Get the matrix dimensions. */
  num_rows = get_num_rows(matrix);
  num_cols = get_num_cols(matrix);

  /* Make sure the matrix is square. */
  if (num_rows != num_cols) {
    die("Cannot normalize non-square matrix (%d != %d).", num_rows, num_cols);
  }

  /* Allocate space for the column. */
  one_column = allocate_array(num_rows);

  /* Iterate until we're told to stop. */
  keep_going = TRUE;
  for (iter = 0; keep_going; iter++) {

    /* First normalize the rows. */
    for (i_row = 0; i_row < num_rows; i_row++) {
      normalize(0.0, get_matrix_row(i_row, matrix));
    }

    /* Then normalize the columns. */
    for (i_row = 0; i_row < num_rows; i_row++) {
      this_row = get_matrix_row(i_row, matrix);

      /* Copy column elements into an array. */
      for (i_col = 0; i_col < num_cols; i_col++) {
  set_array_item(i_col, get_matrix_cell(i_row, i_col, matrix),
           one_column);
      }

      /* Normalize the column. */
      normalize(0.0, one_column);

      /* Re-copy it back into the original position. */
      for (i_col = 0; i_col < num_cols; i_col++) {
  set_matrix_cell(i_row, i_col, get_array_item(i_col, one_column),
      matrix);
      }
    }

    /* Now find out whether all rows are close enough to 1.0. */
    keep_going = FALSE;
    for (i_row = 0; i_row < num_rows; i_row++) {
      total = array_total(get_matrix_row(i_row, matrix));
      if (!almost_equal(1.0, total, tolerance)) {

  /* Tell the user what's up. */
  if (verbosity > NORMAL_VERBOSE) {
    fprintf(stderr, "Iteration %d: row %d has difference of %g\n",
      iter, i_row, 1.0 - total);
  }
  keep_going = TRUE;
  break;
      }

    }
  }

  /* Tell the user what's up. */
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Performed %d iteration", iter);
    if (iter > 1) fprintf(stderr, "s");
    fprintf(stderr, ".\n");
  }

  myfree(one_column);
}

/*****************************************************************************
 * Extract one margin of a matrix and return it as an array.
 *****************************************************************************/
ARRAY_T* get_matrix_row_sums
  (MATRIX_T* matrix)
{
  int      num_rows;
  int      i_row;
  ARRAY_T* this_row;
  ARRAY_T* matrix_sums;

  /* Allocate the margin array. */
  num_rows = get_num_rows(matrix);
  matrix_sums = allocate_array(num_rows);

  /* Get the sums. */
  for (i_row = 0; i_row < num_rows; i_row++) {
    this_row = get_matrix_row(i_row, matrix);
    set_array_item(i_row, array_total(this_row), matrix_sums);
  }

  return(matrix_sums);
}

ARRAY_T* get_matrix_col_sums
  (MATRIX_T* matrix)
{
  int      num_cols;
  int      i_col;
  ARRAY_T* this_column;
  ARRAY_T* matrix_sums;

  /* Allocate the margin array. */
  num_cols = get_num_cols(matrix);
  matrix_sums = allocate_array(num_cols);

  /* Get the sums. */
  for (i_col = 0; i_col < num_cols; i_col++) {
    this_column = get_matrix_column(i_col, matrix);
    set_array_item(i_col, array_total(this_column), matrix_sums);
    free_array(this_column);
  }

  return(matrix_sums);
}

/*****************************************************************************
 * Sort a given matrix by row or column, according to a given set of
 * sort keys.
 *****************************************************************************/
int array_compare
  (const void* elem1,
   const void* elem2)
{
  ARRAY_T* array1 = *((ARRAY_T**)elem1);
  ARRAY_T* array2 = *((ARRAY_T**)elem2);
  ATYPE key1 = get_array_key(array1);
  ATYPE key2 = get_array_key(array2);

  if (key1 < key2) {
    return(-1);
  } else if (key1 > key2) {
    return(1);
  }
  return(0);
}

int reverse_array_compare
  (const void* elem1,
   const void* elem2)
{
  ARRAY_T* array1 = *((ARRAY_T**)elem1);
  ARRAY_T* array2 = *((ARRAY_T**)elem2);
  ATYPE key1 = get_array_key(array1);
  ATYPE key2 = get_array_key(array2);

  if (key1 > key2) {
    return(-1);
  } else if (key1 < key2) {
    return(1);
  }
  return(0);
}

void sort_matrix_rows
  (BOOLEAN_T reverse_sort,
   ARRAY_T*  keys,
   MATRIX_T* matrix)
{
  int      num_rows;
  int      i_row;

  /* Make sure the dimensions match up. */
  num_rows = get_num_rows(matrix);
  if (num_rows != get_array_length(keys)) {
    die("Tried to sort a matrix with %d rows using an array of %d keys.",
  num_rows, get_array_length(keys));
  }

  /* Put the keys into the rows. */
  for (i_row = 0; i_row < num_rows; i_row++) {
    set_array_key(get_array_item(i_row, keys), get_matrix_row(i_row, matrix));
  }

  /* Sort 'em! */
  if (reverse_sort) {
    qsort(matrix->rows, num_rows, sizeof(ARRAY_T*), reverse_array_compare);
  } else {
    qsort(matrix->rows, num_rows, sizeof(ARRAY_T*), array_compare);
  }

  /* Sort the keys, too. */
  sort_array(reverse_sort, keys);
}

/*****************************************************************************
 * Randomly shuffle the entries of a given matrix.
 *****************************************************************************/
void shuffle_matrix
  (MATRIX_T* matrix)
{
  int   i_row;
  int   i_col;
  int   entry;
  int   i_row_rand;
  int   i_col_rand;
  int   num_rows;
  int   num_cols;
  int   num_entries;
  MTYPE temp;

  /* Get bounds */
  num_rows = get_num_rows(matrix);
  num_cols = get_num_cols(matrix);
  num_entries = num_rows * num_cols;

  /* Randomly shuffle the matrix in one pass. */
  for (i_row = 0; i_row < num_rows; i_row++)
    {
      for (i_col = 0; i_col < num_cols; i_col++)
  {
    /* Call pseudo random number generator only once per entry. */
    entry = my_drand() * num_entries;
    i_row_rand = entry / num_cols;
    i_col_rand = entry % num_cols;

    /* Swap entries */
    temp = get_matrix_cell(i_row, i_col, matrix);
    set_matrix_cell(i_row, i_col,
        get_matrix_cell(i_row_rand, i_col_rand, matrix),
        matrix);
    set_matrix_cell(i_row_rand, i_col_rand, temp, matrix);
  }
    }
}

/*****************************************************************************
 * Randomly shuffle the rows or columns of a matrix. If the boolean cols is TRUE * shuffle the columns else shuffle the rows
 *****************************************************************************/
void shuffle_matrix_cols
  (MATRIX_T* matrix, BOOLEAN_T cols)
{
  int   i_row;
  int   i_col;
  int   i_col_rand;
  int   i_row_rand;
  int   num_rows;
  int   num_cols;

  /* Get bounds */
  num_rows = get_num_rows(matrix);
  num_cols = get_num_cols(matrix);

  if (cols){
    /* Call pseudo random number generator only once per entry. */
    for (i_col = 0; i_col < num_cols; i_col++){
      i_col_rand = my_drand() * (num_cols-1);
      /* Swap the columns */
      ARRAY_T* temp_array = NULL;
      temp_array = get_matrix_column(i_col, matrix);
      ARRAY_T* rand_col = NULL;
      rand_col = get_matrix_column(i_col_rand, matrix);
      set_matrix_column(
      rand_col,
      i_col,
      matrix
      );
      set_matrix_column(temp_array, i_col_rand, matrix);
      free_array(temp_array);
      free_array(rand_col);
    }
  }
  else{
    /* Call pseudo random number generator only once per entry. */
    for (i_row = 0; i_row < num_rows; i_row++){
      i_row_rand = my_drand() * (num_rows-1);
      /*       Swap the rows */
      int temp_index;
      for (temp_index = 0; temp_index < num_cols; temp_index++){
  double temp_var = get_matrix_cell(i_row, temp_index, matrix);
  set_matrix_cell(i_row, temp_index,
      get_matrix_cell(i_row_rand, temp_index, matrix),
      matrix
      );
  set_matrix_cell(i_row_rand, temp_index,
      temp_var,
      matrix
      );
      }//temp_index
    }//i_row
  }//else
}//shuffle_matrix_cols


/***********************************************************************
 * Multiply two matrices to get a third.
  ***********************************************************************/
MATRIX_T* matrix_multiply
  (MATRIX_T* matrix1,
   MATRIX_T* matrix2)
{
  int i, j, p;
  int r1 = get_num_rows(matrix1);
  int r2 = get_num_rows(matrix2);
  int c1 = get_num_cols(matrix1);
  int c2 = get_num_cols(matrix2);
  MATRIX_T *product;

  /* Check to be sure the dimensions match. */
  if (c1 != r2) {
    die("Tried to multiply %d x %d matrix times a %d x %d matrix.\n",
       r1, c1, r2, c2);
  }

  /* create the product matrix */
  product = allocate_matrix(r1, c2);

  for (i=0; i<r1; i++) {
    for (j=0; j<c2; j++) {
      double sum;
      /* dot product of the row of matrix 1 and column of matrix 2 */
      for (sum = p = 0; p < c1; p++) {
        sum += get_matrix_cell(i, p, matrix1) * get_matrix_cell(p, j, matrix2);
      }
      set_matrix_cell(i, j, sum, product);
    }
  }

  return(product);
} /* matrix_multiply */

/***********************************************************************
 * Get a selection of the matrix rows from index_start to index_stop
 * (including both indices).
  ***********************************************************************/
void get_matrix_rows (int index_start, int index_stop, MATRIX_T* back_freqs, MATRIX_T** subset){
  *subset = allocate_matrix(
          index_stop-index_start+1,
          get_num_cols(back_freqs)
          );
  int temp_index;
  int i_row;
  i_row=0;
  for(temp_index=index_start;temp_index<=index_stop;temp_index++){
    set_matrix_row(i_row,get_matrix_row(temp_index,back_freqs),*subset);
    i_row++;
  }
}



/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
