/*************************************************************************
 * FILE: remove-alignment-gaps.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 1 July 2008
 * PROJECT: MEME
 * COPYRIGHT: 2008, UW
 * DESCRIPTION: Remove from an alignment all columns that contain a gap
 *              in the specified species.
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
  if (argc != 3) {
    fprintf(stderr, "USAGE: remove-alignment-gaps <species> <alignment>\n");
    exit(1);
  }
  char* species = argv[1];
  char* alignment_filename = argv[2];

  // Read the alignment.
  ALIGNMENT_T* big_alignment = read_alignment_from_file(
    alignment_filename,
	  FALSE, // sort by species name
    FALSE, // remove gaps
    NULL   // pointer to ref_seq_index not used
  );
  fprintf(stderr, "Read alignment of %d sequences and %d columns.\n",
	  get_num_aligned_sequences(big_alignment),
	  get_alignment_length(big_alignment));

  // Remove columns containing gaps in the specified species.
  ALIGNMENT_T* small_alignment = remove_alignment_gaps(species, big_alignment);
  fprintf(stderr, "Created alignment of %d sequences and %d columns.\n",
	  get_num_aligned_sequences(small_alignment),
	  get_alignment_length(small_alignment));

  // Print the reduced alignment.
  print_clustalw(stdout, FALSE, small_alignment);

  // Free locally allocated memory.
  free_alignment(big_alignment);
  free_alignment(small_alignment);

  return(0);
}

#endif
