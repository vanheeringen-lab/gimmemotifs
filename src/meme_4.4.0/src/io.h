/*
 * $Id: io.h 1048 2006-07-06 20:07:44Z cegrant $
 * 
 * $Log$
 * Revision 1.2  2006/03/08 20:50:12  nadya
 * merge chamges from v3_5_2 branch
 *
 * Revision 1.1.1.1.4.1  2006/01/26 09:16:27  tbailey
 * Rename local function getline() to getline2() to avoid conflict with
 * function defined in stdio.h.
 *
 * Revision 1.1.1.1  2005/07/29 18:40:41  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

#ifndef IO_H
#define IO_H

#include "utils.h"

extern char *getline2(
  FILE *stream 						/* input stream */
);

int create_output_directory(
  char *output_dirname,	/* Name of the output directory to create */
  BOOLEAN_T clobber,	/* Whether or not to overwrite an existing dir */
  BOOLEAN_T warn	/* Print warning/informative messages to stderr? */
);

#endif
