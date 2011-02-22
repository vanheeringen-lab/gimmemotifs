
/********************************************************************
 * FILE: read_csv.c
 * AUTHOR: Fabian Buske
 * CREATE DATE: 18/06/2008
 * PROJECT: MEME suite
 * COPYRIGHT: 2008, UQ
 * 
 * 	Read a csv file.
 * 	Returns a list of String-Lists (one for each line).
 * 
 * 	Supports CSV formats.
 ********************************************************************/

#include "read_csv.h"
#include "string-list.h"
#include <stdlib.h>
#include <stdio.h>

BOOLEAN read_csv(
  FILE *data_file,		/* file containing sequences */
  char comment,			/* the character that define a comment line */
  char separator,		/* the character used for splitting a line */
  STRING_LIST_T ***content,	/* retrieved content of the csv file as list of string lists */
  int *length			/* retrieve length of the file (number of elements in content) */
)
{
  int i;
//  long j;
  char c;
  c = ' ';
  char* entry = NULL;
  
  if (data_file == NULL ){            	/* Could not open file */
	  die("Error with file! \n");
  }
  
  STRING_LIST_T* line = new_string_list();
  /* get the sample name */
  /* read to first blank/tab/ or end of line/file */
  for (i=0,*length=0; (c=fgetc(data_file))!=EOF; ) {
    if (c==comment) {
    	Skip_eol(c, data_file);	/* comment ends line */
    } else if (c== '\n' || c== '\r') {
    	Resize(entry, i+1, char);
    	entry[i] = '\0';
    	add_string(entry,line); /* return ends entry & line */
    	myfree(entry);
    	entry = NULL;
    	i=0;
    	/* resize array of pointers if necessary */
    	if ((*length % SCHUNK) == 0) {
    		Resize(*content, *length+SCHUNK, STRING_LIST_T*);
    	}
    	(*content)[(*length)++] = line;
		line = new_string_list();
    	
    } else if (c==separator) {
    	Resize(entry, i+1, char);
    	entry[i] = '\0';
    	add_string(entry,line); /* separator ends entry */
    	myfree(entry);
    	entry = NULL;
    	i=0;
    } else {
      if ((i % RCHUNK) == 0) {
        Resize(entry, i+RCHUNK, char);
      }
      entry[i++] = c;			/* non-blank: add to name */
    }
  }
  /* finish the last entry and line */
  Resize(entry, i+1, char);
  entry[i] = '\0';
  add_string(entry,line); 
  
  if (get_num_strings(line)==1  && strcmp(get_nth_string(0,line),entry)==0){
	  Resize(*content, *length+1, STRING_LIST_T*); /* last line contains no information so skip */
	  free_string_list(line);
  } else {
	  Resize(*content, *length+1, STRING_LIST_T*); /* add last line */
	  (*content)[(*length)++] = line;
  }
  myfree(entry);
  
  /* file contains no data */
  if (*length <= 0) {
    if (line != NULL)
    	free_string_list(line);
    myfree(*content);
    return FALSE;
  }

  return TRUE;
} /* read_csv */
