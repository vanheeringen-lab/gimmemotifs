/*****************************************************************************
 * FILE: rdb-matrix.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 4/19/99
 * PROJECT: shared
 * COPYRIGHT: 1999-2008, WSN
 * VERSION: $Revision: 1.1.1.1 $
 * DESCRIPTION: A matrix data structure with labels on rows and columns.
 *****************************************************************************/
#include "rdb-matrix.h"
#include "string-list.h"
#include "matrix.h"
#include "array.h"
#include "utils.h"
#include <string.h>
#include <stdio.h>
#include <assert.h>

/***********************************************************************
 * Allocate memory for an RDB matrix.
 *
 * Optional third argument is a raw matrix to be stored in the new RDB
 * matrix.
 ***********************************************************************/
RDB_MATRIX_T* allocate_rdb_matrix
  (int num_rows,
   int num_cols,
   MATRIX_T* matrix)
{
  RDB_MATRIX_T* return_value;

  /* Allocate the new RDB matrix. */
  return_value = (RDB_MATRIX_T*)mm_malloc(sizeof(RDB_MATRIX_T));

  return_value->row_names = new_string_list();
  return_value->col_names = new_string_list();
  if (matrix == NULL) {
    return_value->matrix = allocate_matrix(num_rows, num_cols);
  } else {
    return_value->matrix = matrix;
  }

  return(return_value);
}

/***********************************************************************
 * Convert a given matrix to an RDB matrix.
 ***********************************************************************/
RDB_MATRIX_T* rdbize_matrix
  (char*          corner_string,
   STRING_LIST_T* row_names,
   STRING_LIST_T* col_names,
   MATRIX_T*      matrix)
{
  RDB_MATRIX_T* new_matrix;
  assert(get_num_rows(matrix) == get_num_strings(row_names));
  assert(get_num_cols(matrix) == get_num_strings(col_names));

  new_matrix = allocate_rdb_matrix(get_num_rows(matrix),
				   get_num_cols(matrix),
				   matrix);

  set_corner_string(corner_string, new_matrix);
  set_row_names(row_names, new_matrix);
  set_col_names(col_names, new_matrix);
  
  return(new_matrix);
}

/***********************************************************************
 * Get or set the various pieces of the RDB matrix.
 ***********************************************************************/
STRING_LIST_T* get_row_names
  (RDB_MATRIX_T* rdb_matrix)
{
  if (rdb_matrix == NULL) {
    die("Attempted to access null matrix.");
  }
  return(rdb_matrix->row_names);
}

STRING_LIST_T* get_col_names
  (RDB_MATRIX_T* rdb_matrix)
{
  if (rdb_matrix == NULL) {
    die("Attempted to access null matrix.");
  }
  return(rdb_matrix->col_names);
}

void set_row_names
  (STRING_LIST_T* row_names,
   RDB_MATRIX_T* rdb_matrix)
{
  int num_rows;
  int i_row;

  if (rdb_matrix == NULL) {
    die("Attempted to access null matrix.");
  }

  num_rows = get_num_rows(get_raw_matrix(rdb_matrix));
  if (get_num_strings(row_names) != num_rows) {
    die("Adding %d row names to a matrix with %d rows.",
	get_num_strings(row_names), num_rows);
  }

  /* Add the strings. */
  for (i_row = 0; i_row < num_rows; i_row++) {
    add_string(get_nth_string(i_row, row_names), rdb_matrix->row_names);
  }
}

void set_col_names
  (STRING_LIST_T* col_names,
   RDB_MATRIX_T*  rdb_matrix)
{
  int i_col;
  int num_cols;

  if (rdb_matrix == NULL) {
    die("Attempted to access null matrix.");
  }

  num_cols = get_num_cols(get_raw_matrix(rdb_matrix));
  if (get_num_strings(col_names) != num_cols) {
    die("Adding %d column names to a matrix with %d columns.",
	get_num_strings(col_names), get_num_cols(get_raw_matrix(rdb_matrix)));
  }

  /* Add the strings. */
  for (i_col = 0; i_col < num_cols; i_col++) {
    add_string(get_nth_string(i_col, col_names), rdb_matrix->col_names);
  }
}

char* get_corner_string
  (RDB_MATRIX_T* rdb_matrix)

{
  if (rdb_matrix == NULL) {
    die("Attempted to access null matrix.");
  }
  return(rdb_matrix->corner_string);
}

void set_corner_string
  (char* corner_string,
   RDB_MATRIX_T* rdb_matrix)
{
  if (rdb_matrix == NULL) {
    die("Attempted to access null matrix.");
  }
  copy_string(&(rdb_matrix->corner_string), corner_string);
}


MATRIX_T* get_raw_matrix
  (RDB_MATRIX_T* rdb_matrix)
{
  if (rdb_matrix == NULL) {
    die("Attempted to access null matrix.");
  }
  return(rdb_matrix->matrix);
}

void set_raw_matrix
  (MATRIX_T*     matrix,
   RDB_MATRIX_T* rdb_matrix)
{
  if (rdb_matrix == NULL) {
    die("Attempted to access null matrix.");
  }
  free_matrix(rdb_matrix->matrix);
  rdb_matrix->matrix = matrix;
}

/***********************************************************************
 * Add a column to a matrix.
 ***********************************************************************/
void add_rdb_column
  (char*         col_name,
   ARRAY_T*      new_column,
   RDB_MATRIX_T* rdb_matrix)
{
  int i_col;
  int i_row;
  MATRIX_T* raw_matrix;
  int num_rows;
  MTYPE this_value;

  if (rdb_matrix == NULL) {
    die("Attempted to access null matrix.");
  }

  /* Add the new name. */
  add_string(col_name, rdb_matrix->col_names);

  /* Find out how many columns we have. */
  i_col = get_num_strings(rdb_matrix->col_names) - 1;

  /* Get the raw matrix. */
  raw_matrix = get_raw_matrix(rdb_matrix);

  /* Copy from the array to the matrix. */
  num_rows = get_array_length(new_column);
  for (i_row = 0; i_row < num_rows; i_row++) {
    this_value = get_array_item(i_row, new_column);
    set_matrix_cell(i_row, i_col, this_value, raw_matrix);
  }
}

/***********************************************************************
 * Add a row or column name to a matrix.
 ***********************************************************************/
void add_row_name
  (char*         row_name,
   RDB_MATRIX_T* rdb_matrix)
{
  if (rdb_matrix == NULL) {
    die("Attempted to access null matrix.");
  }
  add_string(row_name, rdb_matrix->row_names);
}

void add_col_name
  (char*         col_name,
   RDB_MATRIX_T* rdb_matrix)
{
  if (rdb_matrix == NULL) {
    die("Attempted to access null matrix.");
  }
  add_string(col_name, rdb_matrix->col_names);
}

/***********************************************************************
 * Is the given line a comment?
 ***********************************************************************/
static BOOLEAN_T is_comment
  (char* one_line)
{
  int  i_char;
  char this_char;

  /* Find the first non-space character. */
  i_char = 0;
  this_char = one_line[i_char];
  while (this_char == ' ') {
    this_char = one_line[i_char++];
  }

  /* Is it a hash mark? */
  return(this_char == '#');
}

/***********************************************************************
 * Read one row from a file, checking for long lines.
 ***********************************************************************/
static char *read_one_row
  (FILE* infile,
   int   length,
   char* one_row)
{
  char*     fgets_result;       /* Error indicator for 'fgets' function. */

  fgets_result = fgets(one_row, length, infile);

  if (fgets_result != NULL && one_row[strlen(one_row) - 1] != '\n') {
    die("Matrix lines too long.  Increase length.");
  }

  return(fgets_result);
}

/***********************************************************************
 * Read an RDB file into a matrix.
 ***********************************************************************/
#define MAX_ROW 325000
RDB_MATRIX_T* read_rdb_matrix
  (char* sep,			/* token separator */
   BOOLEAN_T has_corner_string,	/* is there a corner string? */
   int num_cols,		// no column hdrs if >0 and this is number
				// of data columns (not counting row names)
   BOOLEAN_T format_line,
   char* missing_marker,
   FILE* infile)
{
  MATRIX_T* matrix;             /* The matrix to be read. */
  char*     corner_string;      /* Upper left corner of matrix. */
  STRING_LIST_T* row_names;     /* Row names in the matrix. */
  STRING_LIST_T* col_names;     /* Column names in the matrix. */
  char      one_row[MAX_ROW];   /* One row of the matrix. */
  int       i_row;              /* Index of the current row. */
  int       i_column;           /* Index of the current column. */
  int       num_rows;           /* Total number of rows. */
  int       num_scanned;        /* Error indicator for fscanf function. */
  char*     string_ptr = NULL;  /* Storage for strtok function. */
  MTYPE     one_value;          /* One value read from the file. */
  ARRAY_T*  this_row = NULL;    /* One row of the matrix. */
  RDB_MATRIX_T* return_value;   /* The RDB matrix being created. */
  char      missing[MAX_ROW];   /* Storage space for missing values. */
  BOOLEAN_T failed;             /* Flag for improper read. */
  BOOLEAN_T row_read = FALSE;	// Flag for line already read.

  /* Make sure the file is valid. */
  if (infile == NULL) {
    die("Attempted to read matrix from null file.");
  }

  /* Create the row names and column names lists. */
  row_names = new_string_list();
  col_names = new_string_list();

  /* Read the first row. */
  if (read_one_row(infile, MAX_ROW, one_row) == NULL || one_row[0] == '\0')
    die("No data in background file\n");

  /* Keep reading till we get past the comments. */
  while (is_comment(one_row)) {
    if (read_one_row(infile, MAX_ROW, one_row) == NULL || one_row[0] == '\0')
      die("Nothing but commments in background file\n");
  }
    
  // Read and store the names of the columns and the corner string.
  if (!num_cols) {
    /* Store the name of the first column or corner string. */
    string_ptr = strtok(one_row, sep);
    if (has_corner_string) {
      copy_string(&corner_string, string_ptr);
    } else {
      add_string(string_ptr, col_names);
    }

    /* Store the names of the columns. */
    for (string_ptr = strtok(NULL, sep); string_ptr != NULL;
	 string_ptr = strtok(NULL, sep)) {

      /* Remove EOL. */
      if (string_ptr[strlen(string_ptr) - 1] == '\n') {
	string_ptr[strlen(string_ptr) - 1] = '\0';
      }

      /* Store string if it is non-empty. */
      if (strcmp(string_ptr, "") != 0) {
	add_string(string_ptr, col_names);
      }
    }
    num_cols = get_num_strings(col_names);
  } // read column names

  // Add dummy column names and mark row as read already.
  else {
    row_read = TRUE;
    for (i_column=0; i_column<num_cols; i_column++) 
      add_string("dummy", col_names);
  }

  /* Allocate the matrix. */
  matrix = allocate_matrix(0, num_cols);

  /* Allocate one row. */
  this_row = allocate_array(num_cols);

  /* Skip the format line, if necessary. */
  if (format_line) {
    (void) read_one_row(infile, MAX_ROW, one_row);
  }

  /* Read the matrix. */
  for (i_row = 0; ; i_row++) {

    /* Read the next line, stopping if it's empty. */
    if (!row_read && read_one_row(infile, MAX_ROW, one_row) == NULL) break;

    // Skip embedded comments.
    while (is_comment(one_row)) {
      if (read_one_row(infile, MAX_ROW, one_row) == NULL) break;
    }
    //need to test if we actually have something
    if (one_row[0] == '\0') break;

    /* Read the row name and store it. */
    string_ptr = strtok(one_row, sep);
    add_string(string_ptr, row_names);

    /* Read the row. */
    failed = FALSE;
    for (i_column = 0; i_column < num_cols; i_column++) {

      /* Read the first value. */
      string_ptr = strtok(NULL, sep);

      /* Make sure we got a valid string. */
      if (string_ptr == NULL) {
	failed = TRUE;
        break;
      }

      /* Extract a number */
      num_scanned = sscanf(string_ptr, MSCAN, &one_value);
      if (num_scanned == EOF) {
	failed = TRUE;
        break;
      }

      /* If we didn't get a number, deal with missing values. */
      if (num_scanned == 0) {
	num_scanned = sscanf(string_ptr, "%s", missing);
	if ((missing_marker == NULL) ||
	    (num_scanned == 0) ||
            (strcmp(missing, missing_marker) != 0)) {
	  failed = TRUE;
          break;
        }
	if (strcmp(ATYPENAME, "double") != 0) {
          fprintf(stderr,
            "Missing values only permitted in double type array.\n");
          failed = TRUE;
          break;
        }
	one_value = NaN();
      }

      /* Store the value. */
      set_array_item(i_column, one_value, this_row);
    }

    if (failed) {
      die("Error reading matrix at position (%d,%d). ", i_row, i_column);
    }

    /* Add this row to the matrix. */
    grow_matrix(this_row, matrix);
  
    row_read = FALSE;
  }
  num_rows = i_row - 1;

  /* Assemble it all into an RDB matrix. */
  return_value = allocate_rdb_matrix(num_rows, num_cols, matrix);
  if (!has_corner_string) copy_string(&corner_string, "X");
  set_corner_string(corner_string, return_value);
  set_row_names(row_names, return_value);
  set_col_names(col_names, return_value);

  /* Free local dynamic memory. */
  myfree(corner_string);
  free_array(this_row);
  free_string_list(row_names);
  free_string_list(col_names);

  return(return_value);
}

/***********************************************************************
 * Write a named array to a file.
 ***********************************************************************/
void print_rdb_array
  (STRING_LIST_T* names,         /* Structure storing names */ 
   char*          label1,        /* Label for column 1. */
   char*          label2,        /* Label for column 2. */
   ARRAY_T*       array,         /* The array to be printed. */
   int            width,         /* Width of each cell. */
   int            precision,     /* Precision of each cell. */
   FILE*          outfile)       /* File to which to write. */
{
  
  int   i_item;
  int   num_items;
  double item;

  fprintf(outfile, "%s\t%s\n", label1, label2);
  fprintf(outfile, "5S\t5N\n");

  num_items = get_array_length(array);
  for (i_item = 0; i_item < num_items; i_item++) {
    item = get_array_item(i_item, array);
    fprintf(outfile, "%s\t%*.*f\n", get_nth_string(i_item, names), width, 
            precision, item);
  }
}

/***********************************************************************
 * Write a labeled matrix in RDB format.
 ***********************************************************************/
void print_rdb_matrix
  (RDB_MATRIX_T*  rdb_matrix,    /* Matrix to be printed. */
   BOOLEAN_T      format_line,   /* Print format line? */
   int            width,         /* Width of each cell. */
   int            precision,     /* Precision of each cell. */
   FILE*          outfile)       /* File to which to write. */
{
  int   num_cols;
  int   num_rows;
  int   i_col;
  int   i_row;
  MTYPE item;

  /* Print the column names. */
  fprintf(outfile, "%s", rdb_matrix->corner_string);
  num_cols = get_num_cols(rdb_matrix->matrix);
  for (i_col = 0; i_col < num_cols; i_col++) {
    fprintf(outfile, "\t%s", get_nth_string(i_col, rdb_matrix->col_names));
  }
  fprintf(outfile, "\n");

  /* If requested, print the column widths. */
  if (format_line) {
    fprintf(outfile, "10S");
    for (i_col = 0; i_col < num_cols; i_col++) {
      fprintf(outfile, "\t%dN", width);
    }
    fprintf(outfile, "\n");
  }

  /* Print the matrix. */
  num_rows = get_num_rows(rdb_matrix->matrix);
  for (i_row = 0; i_row < num_rows; i_row++) {
    fprintf(outfile, "%s", get_nth_string(i_row, rdb_matrix->row_names));
    for (i_col = 0; i_col < num_cols; i_col++) {
      item = get_matrix_cell(i_row, i_col, rdb_matrix->matrix);
      fprintf(outfile, "\t%*.*g", width, precision, item);
    }
    fprintf(outfile, "\n");    
  }
}


/***********************************************************************
 * Free an RDB matrix.
 ***********************************************************************/
void free_rdb_matrix
  (RDB_MATRIX_T* rdb_matrix)
{
  if (rdb_matrix != NULL) {
    myfree(rdb_matrix->corner_string);
    free_string_list(rdb_matrix->row_names);
    free_string_list(rdb_matrix->col_names);
    free_matrix(rdb_matrix->matrix);
    myfree(rdb_matrix);
  }
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
