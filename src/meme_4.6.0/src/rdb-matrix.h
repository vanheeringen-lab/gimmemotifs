/*****************************************************************************
 * FILE: rdb-matrix.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 4/19/99
 * PROJECT: shared
 * COPYRIGHT: 1999-2008, WSN
 * VERSION: $Revision: 1.1.1.1 $
 * DESCRIPTION: A matrix data structure with labels on rows and columns.
 *****************************************************************************/
#ifndef RDB_MATRIX_H
#define RDB_MATRIX_H

#include "matrix.h"
#include "array.h"
#include "string-list.h"
#include "utils.h"

/***************************************************************************
 * Define the RDB matrix type.
 ***************************************************************************/
struct rdb_matrix_t {
  char*          corner_string; /* Upper left corner. */
  STRING_LIST_T* row_names;
  STRING_LIST_T* col_names;
  MATRIX_T*      matrix;
};
typedef struct rdb_matrix_t RDB_MATRIX_T;


/***********************************************************************
 * Allocate memory for an RDB matrix.
 ***********************************************************************/
RDB_MATRIX_T* allocate_rdb_matrix
  (int num_rows,
   int num_cols,
   MATRIX_T* matrix);

/***********************************************************************
 * Convert a given matrix to an RDB matrix.
 ***********************************************************************/
RDB_MATRIX_T* rdbize_matrix
  (char*          corner_string,
   STRING_LIST_T* row_names,
   STRING_LIST_T* col_names,
   MATRIX_T*      matrix);

/***********************************************************************
 * Get or set the various pieces of the RDB matrix.
 ***********************************************************************/
STRING_LIST_T* get_row_names
  (RDB_MATRIX_T* rdb_matrix);

STRING_LIST_T* get_col_names
  (RDB_MATRIX_T* rdb_matrix);

void set_row_names
  (STRING_LIST_T* row_names,
   RDB_MATRIX_T* rdb_matrix);

void set_col_names
  (STRING_LIST_T* col_names,
   RDB_MATRIX_T*  rdb_matrix);

char* get_corner_string
  (RDB_MATRIX_T* rdb_matrix);

void set_corner_string
  (char* corner_string,
   RDB_MATRIX_T* rdb_matrix);

MATRIX_T* get_raw_matrix
  (RDB_MATRIX_T* rdb_matrix);

void set_raw_matrix
  (MATRIX_T*     matrix,
   RDB_MATRIX_T* rdb_matrix);

/***********************************************************************
 * Add a column to a matrix.
 ***********************************************************************/
void add_rdb_column
  (char*         col_name,
   ARRAY_T*      new_column,
   RDB_MATRIX_T* rdb_matrix);

/***********************************************************************
 * Add a column name to a matrix.
 ***********************************************************************/
void add_col_name
  (char*         col_name,
   RDB_MATRIX_T* rdb_matrix);

/***********************************************************************
 * Read an RDB file into a matrix.
 ***********************************************************************/
RDB_MATRIX_T* read_rdb_matrix
  (char* sep,			 /* token separator */
   BOOLEAN_T has_corner_string,  /* is there a corner string? */
   int num_cols,                // no column hdrs if >0 and this is number
                                // of data columns (not counting row names)
   BOOLEAN_T format_line,
   char* missing_marker,         /* String to identify missing value */
   FILE* infile);

/***********************************************************************
 * Write a labeled matrix in RDB format.
 ***********************************************************************/
void print_rdb_matrix
  (RDB_MATRIX_T*  rdb_matrix,    /* Matrix to be printed. */
   BOOLEAN_T      format_line,   /* Print format line? */
   int            width,         /* Width of each cell. */
   int            precision,     /* Precision of each cell. */
   FILE*          outfile);      /* File to which to write. */

/***********************************************************************
 * Free an RDB matrix.
 ***********************************************************************/
void free_rdb_matrix
  (RDB_MATRIX_T* rdb_matrix);

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
