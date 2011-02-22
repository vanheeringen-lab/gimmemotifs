/*
 * $Id: banner.c 3690 2009-02-23 22:40:50Z cegrant $
 * 
 * $Log$
 * Revision 1.2  2005/10/25 19:06:39  nadya
 * rm old macro for Header, all info is taken care of by Id and Log.
 *
 * Revision 1.1.1.1  2005/07/29 17:15:49  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/***********************************************************************
*                                                                      *
*       MEME++                                                         *
*       Copyright 1994, The Regents of the University of California    *
*       Author: Timothy L. Bailey                                      *
*                                                                      *
***********************************************************************/
/*
        Print name, version, date, reference.
*/

#include "meme.h"
#include "projrel.h"

extern void banner(
  char *program, /* name of program */
  FILE *outfile  /* destination for output */
) 
{
  const char *archive_date = ARCHIVE_DATE;
  int i = strlen(archive_date);

  /* announce the program */
  PSTARS(outfile); 
  fprintf(outfile, "MEME - Motif discovery tool\n"); 
  PSTARS(outfile);
  fprintf(
    outfile,
    "MEME version %s (Release date: %.*s)\n\n"
    "For further information on how to interpret these results or to get\n"
    "a copy of the MEME software please access http://meme.nbcr.net.\n"
    "\n"
    "This file may be used as input to the MAST algorithm for searching\n"
    "sequence databases for matches to groups of motifs.  MAST is available\n"
    "for interactive use and downloading at http://meme.nbcr.net.\n",
    VERSION, i, archive_date
  );
  PSTARS(outfile);

  /* print reference citation */
  fprintf(outfile, "\n\n");
  PSTARS(outfile); 
  fprintf(outfile, "REFERENCE\n"); 
  PSTARS(outfile);
  fprintf(
    outfile,
    "If you use this program in your research, please cite:\n"
    "\n"
    "Timothy L. Bailey and Charles Elkan,\n"
    "\"Fitting a mixture model by expectation maximization to discover\n"
    "motifs in biopolymers\", Proceedings of the Second International\n"
    "Conference on Intelligent Systems for Molecular Biology, pp. 28-36,\n"
    "AAAI Press, Menlo Park, California, 1994.\n"
  );
  PSTARS(outfile);
}
