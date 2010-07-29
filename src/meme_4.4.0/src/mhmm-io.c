/**************************************************************************
 * FILE: mhmm-io.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 1-28-98
 * PROJECT: MHMM
 * DESCRIPTION: A simple program to test I/O of MHMMs.
 **************************************************************************/
#include "utils.h"
#include "mhmm-state.h"
#include "read-mhmm.h" 
#include "write-mhmm.h"
#include "simple-getopt.h"
#include <string.h>

VERBOSE_T verbosity = INVALID_VERBOSE;

/**************************************************************************
 * int main
 **************************************************************************/
int main (int argc, char *argv[])
{
  char *    hmm_filename = NULL;  // File containing the HMM.
  FILE *    hmm_file;
  MHMM_T *  the_hmm;        // The HMM itself.

	// Define command line options.
  cmdoption options[] = {
    {"verbosity", REQUIRED_VALUE}
  };

  int option_count = 1;
  int option_index = 0; // index to next argument

  // Define the usage message.
  char      usage[1000] = "";
  strcat(usage, "USAGE: mhmm-io [options] <HMM>\n");
  strcat(usage, "\n");
  strcat(usage, "   Options:\n");
  strcat(usage, "     -verbosity 1|2|3|4|5 (default=2)\n");
  strcat(usage, "\n");
 
  simple_setopt(argc, argv, option_count, options);

  // Parse the command line.
  while (1) { 
    int c = 0;
    char* option_name = NULL;
    char* option_value = NULL;
    const char* message = NULL;

     // Read the next option, and break if we're done.
    c = simple_getopt(&option_name, &option_value, &option_index);
    if (c == 0) {
      break;
    } else if (c < 0) {
    	simple_getopterror(&message);
      die("Error processing command line options (%s)\n", message);
    }

    if (strcmp(option_name, "verbosity") == 0) {
      verbosity = (VERBOSE_T)atoi(option_value);
    } 
  }

  // Read the single required argument.
  if (option_index + 1 != argc) {
    fprintf(stderr, "%s", usage);
    exit(1);
  }
  hmm_filename = argv[option_index];

  // Read the model.
  read_mhmm(hmm_filename, &the_hmm);

  // Print the HMM.
  write_mhmm(verbosity, the_hmm, stdout);
  free_mhmm(the_hmm);

  return(0);
}
