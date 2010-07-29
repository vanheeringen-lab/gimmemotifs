/*************************************************************************
 * FILE: reduce-alignment.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 1 July 2008
 * PROJECT: MEME
 * COPYRIGHT: 2008, UW
 * DESCRIPTION: Extract specified range of columns from an alignment.
 *************************************************************************/
#include "alignment.h"
#include "clustalw-io.h"
#include <stdio.h>

/*****************************************************************************
 * MAIN
 *****************************************************************************/
#ifdef MAIN

VERBOSE_T verbosity = NORMAL_VERBOSE;

int main
  (int    argc,
   char * argv[])
{

  // Parse the command line.
  if (argc != 4) {
    fprintf(stderr, "USAGE: reduce-alignment <start> <width> <alignment>\n");
    exit(1);
  }
  int start_position = atoi(argv[1]);
  int width = atoi(argv[2]);
  char* alignment_filename = argv[3];

  // Read the alignment.
  ALIGNMENT_T* big_alignment = read_alignment_from_file(
    alignment_filename,
		FALSE, 
    FALSE, 
    NULL // pointer to ref_seq_index, not used.
  );
  fprintf(stderr, "Read alignment of %d sequences and %d columns.\n",
	  get_num_aligned_sequences(big_alignment),
	  get_alignment_length(big_alignment));

  if (start_position + width > get_alignment_length(big_alignment)) {
    fprintf(stderr, "Invalid coordinates: %d + %d > %d.\n",
	    start_position, width, get_alignment_length(big_alignment));
    exit(1);
  }

  // Extract the smaller alignment.
  ALIGNMENT_T* small_alignment = extract_subalignment(start_position,
						      width,
						      big_alignment);
  fprintf(stderr, "Created alignment of %d sequences and %d columns.\n",
	  get_num_aligned_sequences(small_alignment),
	  get_alignment_length(small_alignment));

  // Print the alignment.
  print_clustalw(stdout, FALSE, small_alignment);

  // Free locally allocated memory.
  free_alignment(big_alignment);
  free_alignment(small_alignment);

  return(0);
}

#endif
