/****************************************************************************
 * FILE: clustalw2phylip.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 21 July 2008 (at ISMB)
 * DESCRIPTION: Convert an alignment from clustalw to phyip format.
 * COPYRIGHT: 2008, UW
 ****************************************************************************/
#include "utils.h"
#include "alignment.c"

VERBOSE_T verbosity = NORMAL_VERBOSE;

/*************************************************************************
 * int main
 *************************************************************************/
int main(int argc, char *argv[]) {

  if (argc != 2) {
    fprintf(stderr, "USAGE: clustalw2phyilip <clustalw file>\n");
    exit(1);
  }
  char* clustalw_filename = argv[1];

  // Open the file for reading.
  FILE* clustalw_file;
  if (open_file(clustalw_filename, "r", TRUE, "CLUSTAL", "alignment", 
        &clustalw_file) == 0) {
    exit(1);
  }

  // Read the alignment
  ALIGNMENT_T* alignment = NULL;
  BOOLEAN_T result = read_clustalw(clustalw_file, &alignment);
  if (result) {
    print_phylip_alignment(alignment, stdout);
    free_alignment(alignment);
  } else {
    fprintf(stderr, "Unable to parse an alignment from the file %s.\n",
        clustalw_filename);
  }

  fclose(clustalw_file);
  return(!result);
}
