#include <stdio.h>
#include <string.h>
#include "clustalw-io.h"
#include "fasta-io.h"
#include "seq.h"
#include "simple-getopt.h"
#include "utils.h"

VERBOSE_T verbosity = NORMAL_VERBOSE;

void print_seqs(BOOLEAN_T show_gaps, int num_seq, SEQ_T** sequences) {

  int i = 0;
  int j = 1;
  char* c = NULL;

  if (show_gaps) {
    for (i = 0; i < num_seq; i++) {
      print_fasta(sequences[i], stdout);
    }
  } else {
    for (i = 0; i < num_seq; i++) {
      fprintf(stdout, ">%s\n", get_seq_name(sequences[i]));
      c = get_raw_sequence(sequences[i]);
      j = 0;
      while (*c != '\0') {
        // Only print non-gap symbols
        if (*c != '-') {
          fputc(*c, stdout);
          j++;
          if (j >= BLOCKSIZE) {
            fputc('\n', stdout);
            j = 0;
          }
        }
        c++;
      }
      fputc('\n', stdout);
    }
  }
}

int main(int argc, char* argv[]) {
  
  BOOLEAN_T consensus = FALSE;
  BOOLEAN_T result = FALSE;
  BOOLEAN_T show_gap = TRUE;
  char*     clustalw_filename = NULL;
  int       num_sequences = 0;
  double    threshold = 0.0;
  FILE*     clustalw_file = NULL;
  STRING_LIST_T* seq_order_list = NULL;
  SEQ_T*    sequence = NULL;
  SEQ_T**   sequences = NULL;
  ALIGNMENT_T* alignment = NULL;

  // Define command line options.
  cmdoption const options[] = {
    {"nogap", NO_VALUE},
    {"consensus", REQUIRED_VALUE},
    {"seqorder", REQUIRED_VALUE},
  };
	int option_count = 3;
  int option_index = 0;
  int c = 0;
  char* option_name = NULL;
  char* option_value = NULL;
	const char *  message = NULL;

  // Define the usage message.
  char      usage[400] = "";
  strcat(usage, "USAGE: clustalw2fasta [options] <alignment file>\n");
  strcat(usage, "\n");
  strcat(usage, "   Options:\n");
  strcat(usage, "     -nogap\n");
  strcat(usage, "     -consensus <threshold>\n");
  strcat(usage, "     -seqorder <seqorder filename>\n");
  strcat(usage, "\n");

	simple_setopt(argc, argv, option_count, options);

  // Parse the command line.
  while (1) { 

    // Read the next option, and break if we're done.
    c = simple_getopt(&option_name, &option_value, &option_index);
    if (c == 0) {
      break;
    } else if (c < 0) {
    	simple_getopterror(&message);
      die("Error processing command line options (%s)\n", message);
    }

    if (strcmp(option_name, "nogap") == 0) {
      show_gap = FALSE;
    } else if (strcmp(option_name, "consensus") == 0) {
      consensus = TRUE;
      threshold = atof(option_value)/100.0;
    } else if (strcmp(option_name, "seqorder") == 0) {
      seq_order_list = read_string_list_from_file(option_value);
    }
  }

  // Read the single required argument.
  if (option_index + 1 != argc) {
    fprintf(stderr, "%s", usage);
    exit(1);
  }
  clustalw_filename = argv[option_index];

  // Open the file for reading.
  if (open_file(clustalw_filename, "r", TRUE, "CLUSTAL", "alignment", 
        &clustalw_file) == 0) {
    exit(1);
  }

  // Read the alignment
  result = read_clustalw(clustalw_file, &alignment);
  if (result) {
    if (consensus) {
      sequence = get_consensus_sequence(threshold, alignment);
      print_fasta(sequence, stdout);
      free_seq(sequence);
    } else {
      num_sequences = get_num_aligned_sequences(alignment);
      sequences = get_alignment_sequences(alignment);
      // Sort the sequences in the given order if requested.
      if (seq_order_list) {
        int num_seq_order = get_num_strings(seq_order_list);
        int i, j, next;
        for (i=next=0; i<num_seq_order && next<num_sequences; i++) {
          for (j=next; j<num_sequences; j++) {
            if (strcmp(
                 get_nth_string(i, seq_order_list),
                 get_seq_name(sequences[j]) ) == 0) {
              SWAP(SEQ_T*, sequences[next], sequences[j]);
              next++;
              break;
            }
          }
        }
      }
      print_seqs(show_gap, num_sequences, sequences);
    }
    free_alignment(alignment);
  } else {
    fprintf(stderr, "Unable to parse an alignment from the file %s.\n",
        clustalw_filename);
  }
  if (clustalw_file != NULL) fclose(clustalw_file);

  return 0;
}
