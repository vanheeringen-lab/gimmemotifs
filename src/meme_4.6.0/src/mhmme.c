/********************************************************************
 * FILE: mhmme.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 10/03/01
 * PROJECT: MHMM
 * COPYRIGHT: 2001-2008, WSN
 * DESCRIPTION: Given Meta-MEME model, randomly emit sequences.
 ********************************************************************/
#include "utils.h"
#include "mhmm-state.h"
#include "read-mhmm.h"
#include "alphabet.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include "simple-getopt.h"
#include <string.h>

/**************************************************************************
 * Randomly select an element according to a given discrete
 * probability distribution.
 **************************************************************************/
static double choose
  (ARRAY_T* distribution)
{
  int return_value;      // The selected bin.
  int length;            // The number of bins in the distribution.
  double random_number;  // The random number.
  double total;          // Total of all bins so far.

  // Choose a random number between 0 and 1.
  random_number = my_drand();

  // Compute cumulative total of distribution.
  total = 0.0;
  length = get_array_length(distribution);
  for (return_value = 0; return_value < length; return_value++) {
    total += get_array_item(return_value, distribution);

    // Stop when we reach the randomly selected number.
    if (total >= random_number) {
      return return_value;
    }
  }
  return(length - 1);
}

/**************************************************************************
 * Randomly emit one sequence from a given HMM.  The sequence is
 * printed in FASTA format to the given file.
 **************************************************************************/
#define LINE_LENGTH 70
static void emit_sequence
  (FILE*     outfile,
   BOOLEAN_T print_motifs, // Print motif IDs rather than sequence.
   char*     seqname,
   MHMM_T*   the_hmm)
{
  int current_state_index;
  MHMM_STATE_T* current_state;
  int position;

  // Print the header line.
  fprintf(outfile, ">%s\n", seqname);

  // Transition out of the start state.
  current_state_index = choose(get_matrix_row(0, the_hmm->trans));
  current_state = &(the_hmm->states[current_state_index]);

  // Continue until we reach the end state.
  position = 0;
  while (current_state->type != END_STATE) {

    /* Emit a character.
     *
     * N.B. We must do this even with print_motifs=TRUE, so that the
     * number of calls to my_drand() is the same.  Otherwise, two
     * different runs with the same random seed but with and without
     * print_motifs turned on would not yield the same results.
     */
    char this_char = get_alph_char(choose(current_state->emit));

    // Emit the motif index, rather than the characters.
    if (print_motifs) {
      this_char = current_state->id_char;
    }	
    fprintf(outfile, "%c", this_char);

    // Perhaps emit an end-of-line.
    if (position == LINE_LENGTH) {
      fprintf(outfile, "\n");
      position = 0;
    } else {
      position++;
    }

    // Change to the next state.
    current_state_index = choose(get_matrix_row(current_state_index,
						the_hmm->trans));
    current_state = &(the_hmm->states[current_state_index]);
  }
  fprintf(outfile, "\n");
}


VERBOSE_T verbosity = INVALID_VERBOSE;

/*************************************************************************
 * int main
 *************************************************************************/
int main(int argc, char *argv[])
{
  /**********************************************
   * COMMAND LINE PROCESSING
   **********************************************/

  // Set default options.
  int       num_seqs = 1;
  BOOLEAN_T print_motifs = FALSE;
  long      seed = time(0);

  // Define command line options.
  int option_count = 3;
  cmdoption const options[] = {
    {"numseq", REQUIRED_VALUE},
    {"occurs", NO_VALUE},
    {"seed", REQUIRED_VALUE}
  };

  // Define the usage message.
  char      usage[1000] = "";
  strcat(usage, "USAGE: mhmme [options] <HMM file>\n");
  strcat(usage, "\n");
  strcat(usage, "   Options:\n");
  strcat(usage, "     --numseq <n>\n");
  strcat(usage, "     --occurs\n");
  strcat(usage, "     --seed <n>\n");
  strcat(usage, "\n");

  // Parse the command line.
  simple_setopt(argc, argv, option_count, options);
  int option_index = 0;
  while (1) { 
    int c = 0;
    char* option_name = NULL;
    char* option_value = NULL;
    const char * message = NULL;

    // Read the next option, and break if we're done.
    c = simple_getopt(&option_name, &option_value, &option_index);
    if (c == 0) {
      break;
    } else if (c < 0) {
      simple_getopterror(&message);
      die("Error processing command line options (%s)\n", message);
    }

    // Get the option name (we only use long options).
    if (strcmp(option_name, "numseq") == 0) {
      num_seqs = atoi(option_value);
    } else if (strcmp(option_name, "occurs") == 0) {
      print_motifs = TRUE;
    } else if (strcmp(option_name, "seed") == 0) {
      seed = atoi(option_value);
    }
  }

  // Read the single required argument.
  if (option_index + 1 != argc) {
    fprintf(stderr, "%s", usage);
    exit(1);
  }
  char* hmm_filename = argv[option_index];

  /**********************************************
   * READ THE MODEL.
   **********************************************/

  // Read the model.
  MHMM_T* the_hmm = NULL;
  read_mhmm(hmm_filename, &the_hmm);

  /**********************************************
   * EMIT SEQUENCES.
   **********************************************/

  // Initialize the random number generator.
  my_srand(seed);

  int i_seq;
  for (i_seq = 0; i_seq < num_seqs; i_seq++) {
    char seqname[100];
    sprintf(seqname, "random_sequence%d", i_seq);
    emit_sequence(stdout, print_motifs, seqname, the_hmm);
  }

  free_mhmm(the_hmm);
  return(0);
}
