#ifndef READ_CSV_H
#define READ_CSV_H

#define SCHUNK 100
#define RCHUNK 5

#include "macros.h"
#include "string-list.h"
 /*******************************************************************
 * 
 * 	Read a csv file.
 * 	Returns a list of String-Lists (one for each line).
 * 
 * 	Supports CSV formats.
 ********************************************************************/

BOOLEAN read_csv(
  FILE *data_file,		/* file containing sequences */
  char comment,			/* the character that define a comment line */
  char separator,		/* the character used for splitting a line */
  STRING_LIST_T ***content,	/* retrieved content of the csv file as list of string lists */
  int *length			/* retrieve length of the file (number of elements in content) */
);

#endif
