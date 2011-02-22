/********************************************************************
 * FILE: metameme.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 05-19-98
 * PROJECT: MHMM
 * COPYRIGHT: 2001-2008, WSN
 * DESCRIPTION: Global functions for the Meta-MEME toolkit.
 ********************************************************************/
#include <stdio.h>
#include <time.h>
#include "utils.h"
#include "metameme.h"

/*************************************************************************
 * Print the program header, including help info, if necessary.
 *************************************************************************/
void write_header
  (char*    program,     // Which program called this function? */
   char*    contents,    // What is in this file?
   char*    description, // Text description.
   char*    meme_file,   // Motif filename.
   char*    model_file,  // HMM filename.
   char*    seq_file,    // Name of the input sequence file.
   FILE*    outfile)     // Stream to write to.
{
  time_t current_time;
  const int head_skip = 18; // Skip SVN keyword
  const int tail_skip = 2; // Skip close of SVN keyword
  const char *release_date = ARCHIVE_DATE + head_skip;

  fprintf(outfile, "# The MEME Suite of motif-based sequence analysis tools.\n");
  fprintf(outfile, "# \n");
  int i = strlen(release_date) - tail_skip;
  fprintf(outfile, "# version %s (Release date: %.*s)\n", VERSION,
	  i, release_date);
  fprintf(outfile, "# \n");
  fprintf(outfile, "# Program: %s\n", program);
  fprintf(outfile, "# File contents: %s\n", contents);
  if (description != NULL) {
    fprintf(outfile, "# Model description: %s\n", description);
  }
  if (meme_file != NULL) {
    fprintf(outfile, "# Motif file: %s\n", meme_file);
  }
  if (model_file != NULL) {
    fprintf(outfile, "# HMM file: %s\n", model_file);
  }
  if (seq_file != NULL) {
    fprintf(outfile, "# Sequence file: %s\n", seq_file);
  }
  current_time = time(0);
  fprintf(outfile, "# Create date: %s", ctime(&current_time));
  fprintf(outfile, "# \n");
  fprintf(outfile, "# For further information on how to interpret these ");
  fprintf(outfile, "results or to get a copy\n");
  fprintf(outfile, "# of the MEME Suite software, ");
  fprintf(outfile, "please access http://meme.sdsc.edu.\n");
  fprintf(outfile, "# \n");
  fprintf(outfile, "# If you use this program in your research, please cite\n");
  fprintf(outfile, "# \n");
  fprintf(outfile, "#  William N. Grundy, Timothy L. Bailey, Charles P. ");
  fprintf(outfile, "Elkan and Michael E.\n");
  fprintf(outfile, "#   Baker.  \"Meta-MEME: Motif-based hidden ");
  fprintf(outfile, "Markov models of protein\n");
  fprintf(outfile, "#   families.\" Computer Applications ");
  fprintf(outfile, "in the Biosciences.  13(4):397-406, 1997.\n");
  fprintf(outfile, "# \n");

  fprintf(outfile, "# *****************************************************************************\n");
}


