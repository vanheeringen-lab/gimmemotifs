/*
 * $Id: io.c 2417 2008-01-30 00:06:38Z cegrant $
 * 
 * $Log$
 * Revision 1.2  2006/03/08 20:50:11  nadya
 * merge chamges from v3_5_2 branch
 *
 * Revision 1.1.1.1.4.1  2006/01/26 09:16:27  tbailey
 * Rename local function getline() to getline2() to avoid conflict with
 * function defined in stdio.h.
 *
 * Revision 1.1.1.1  2005/07/29 00:21:03  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

#include "io.h"
#include "macros.h"
#include <sys/stat.h>
#include <errno.h>


/******************************************************************************/
/*
	getline2

	Read a newline or EOF-terminated line from a file.

	Returns a pointer to the line or NULL if at EOF when called.
*/
/******************************************************************************/
extern char *getline2(
  FILE *stream 				/* input stream */
) {
  char *s = NULL;			/* string to return */
  int c;
  int i = 0;				/* current position in string */
  while ((c=getc(stream)) != EOF) {
    if (i % GLBUFSIZ == 0) Resize(s, i+GLBUFSIZ, char); 
    s[i++] = c;
    if (c == '\n') break;		/* end of line */
  }
  if (feof(stream) && i==0) {
    return NULL;
  } else {
    // Add terminating null
    if (i % GLBUFSIZ == 0) Resize(s, i+1, char); 
    s[i] = 0;
    return s;
  }
} /* getline2 */

/**********************************************************************/
/*
   create_output_directory

   Returns 0 if successful, -1 if an error occurs.
*/
/**********************************************************************/
int create_output_directory(
  char *output_dirname,	/* Name of the output directory to create */
  BOOLEAN_T clobber,	/* Whether or not to overwrite an existing dir */
  BOOLEAN_T warn	/* Print warning/informative messages to stderr? */
) 
{

  int result = -1;
  BOOLEAN_T path_is_directory = FALSE;
  BOOLEAN_T path_exists = FALSE;
  struct stat stat_buffer;

  // Does the output directory alredy exist?
  if (stat(output_dirname, &stat_buffer)) {
    if (errno == ENOENT) {
      // stat failed because the path doesn't exist.
      path_exists = FALSE;
      path_is_directory = FALSE;
    }
    else {
      // stat failed for some other reason
      fprintf(
        stderr,
        "Unable to check for status of output directory '%s': %s.\n",
        output_dirname,
        strerror(errno)
      );
      result = -1;
    }
  }
  else {
    path_exists = TRUE;
    path_is_directory = S_ISDIR(stat_buffer.st_mode);
  }

  if (path_exists) {
    if (!path_is_directory) {
      fprintf(
        stderr,
        "A non-directory file named '%s' already exists,\n"
        "so that name can't be used for an output directory.\n",
        output_dirname
      );
      result = -1;
    }
    else {
      if (!clobber) {
        fprintf(
          stderr,
          "The output directory '%s' already exists.\nIts contents will not"
          " be overwritten.\n",
          output_dirname
        );
        result = -1;
      }
      else {
        if (warn) fprintf(
          stderr,
          "The output directory '%s' already exists.\nIts contents will"
          " be overwritten.\n",
          output_dirname
        );
        result = 0;
      }
    }
  }
  else {
    // The directory doesn't exist, so we can create it.
    // Does this accomodate the case where one or more of the
    // parent directories doesn't exit?
    if (mkdir(output_dirname, 0777)) {
      // mkdir failed
      fprintf(
        stderr,
        "Unable to create output directory '%s': %s.\n",
        output_dirname,
        strerror(errno)
      );
      result = -1;
    }
    else {
      result = 0;
      if (warn) fprintf(
        stderr,
        "Writing results to output directory '%s'.\n",
        output_dirname
      );
    }
  }
  return result;
} // create_output_directory
