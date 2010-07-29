/**************************************************************************
 * FILE: mhmmscan.c
 * AUTHOR: William Stafford Noble, Timothy L. Bailey
 * CREATE DATE: 5/21/02
 * PROJECT: MHMM
 * COPYRIGHT: 1998-2002, WNG, 2001, TLB
 * DESCRIPTION: Search a database of sequences using a motif-based
 * HMM.  Similar to mhmms, but allows arbitrarily long sequences and
 * multiple matches per sequence.
 **************************************************************************/

#define DEFINE_GLOBALS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "matrix.h"      // Routines for floating point matrices. 
#include "array.h"       // Routines for floating point arrays. 
#include "alphabet.h"    // The amino acid / nucleotide alphabet. 
#include "dp.h"          // Dynamic programming routines.
#include "fasta-io.h"    // Read FASTA format sequences. 
#include "fcodon.h"
#include "fitevd.h"      // Extreme value distribution routines. 
#include "gendb.h"       // Synthetic sequence generation.
#include "hash_alph.h"
#include "log-hmm.h"     // HMM log/log-odds conversion. 
#include "meme-io.h"     // Read background freq file.
#include "metameme.h"    // Global metameme functions. 
#include "mhmm-state.h"  // HMM data structure. 
#include "mhmms.h"       // Framework for database search.
#include "pssm.h"        // Position-specific scoring matrix.
#include "rdb-matrix.h"  // For reading background files 
#include "read-mhmm.h"   // HMM input/output. 
#include "simple-getopt.h"  // Cmd line processing
#include "utils.h"       // Generic utilities. 

VERBOSE_T verbosity = NORMAL_VERBOSE;

// Total size of DP matrix.
#define MAX_MATRIX_SIZE 10000000

// Number of columns of overlap between shifted versions of the matrix.
#define OVERLAP_SIZE 1000

// Number of matches required before we generate synthetic scores if
// -synth given.
// If fewer matches than this are found, we assume that EM (or fitevd)
// will be unable to accurately estimate the distribution.
// The current value of this was guessed at by assuming that:
//   1) estimation works well if 10:1 ratio of negative:positive matches
//   2) rarely will there be more than 1000 true positives
#define MIN_MATCHES 10000

#ifndef SHIFT_DEBUG
#define SHIFT_DEBUG 0
#endif

// Static variables
int       dp_rows = 0;           // Size of the DP matrix.
int       dp_cols = 0; 
MATRIX_T* dp_matrix = NULL;      // Dynamic programming matrix.
MATRIX_T* trace_matrix = NULL;   // Traceback for Viterbi.
MATCH_T*  complete_match = NULL; // All repeated matches.
MATCH_T*  partial_match = NULL;  // One partial match.

/**************************************************************************
 * Strip the filename off a path+filename string.
 * This function modifies the given string.
 **************************************************************************/
static char* strip_filename(char* filename) {
  int i;

  // Don't bother if it's empty.
  if (strlen(filename) == 0) {
    return(filename);
  }

  // Look for the last slash.
  for (i = strlen(filename) - 1; i >= 0; i--) {
    if (filename[i] == '/') {
      break;
    }
  }

  // Replace the slash with a '\0'.
  filename[i] = '\0';
  return(filename);
}

/**************************************************************************
 * Tell the user how far we are.
 **************************************************************************/
static void user_update(
   int changed,
   int num_seqs,
   int num_segments,
   int num_matches,
   int progress_every
) {
  if (verbosity >= NORMAL_VERBOSE) {
    int num_stored = get_num_stored();
    if ((changed != 0) && (changed % progress_every == 0)) {
      fprintf(stderr, "Sequences: %d Segments: %d Matches: %d Stored: %d\n",
              num_seqs, num_segments, num_matches, num_stored);
    }
  }
} // user_update

/**************************************************************************
 * read_and_score
 *
 * Read and score sequences from the input stream.
 * Set the number of sequences read.
 * Return the score_set.
 *
 **************************************************************************/
static SCORE_SET *read_and_score(
  SCORE_SET *score_set,     // Set of scores.
  FILE*     seq_file,       // Open stream to read from.
  SEQ_T     *sequence,      // Used if seq_file is NULL.
  BOOLEAN_T negatives_only, // Stream contains only negative examples.
  int       max_chars,      // Number of chars to read at once.
  int       max_gap,        // Maximum gap length to allow in matches.
  MHMM_T*   the_log_hmm,    // The HMM, with probs converted to logs.
  BOOLEAN_T motif_scoring,  // Perform motif-scoring.
  BOOLEAN_T use_pvalues,    // Use p-value scoring?
  PROB_T    p_threshold,    // P-value threshold for motif hits.
  PROB_T    dp_threshold,   // Score threshold for DP.
  PROB_T    e_threshold,    // E-value threshold for scores.
  BOOLEAN_T got_evd,        // Have distribution?
  EVD_SET   evd_set,        // EVD data structure.
  BOOLEAN_T print_fancy,    // Print alignments?
  BOOLEAN_T store_gff,      // Store GFF entries?
  int       output_width,   // Width of output, in chars.
  int       align_width,    // Width of one alignment line.
  int       progress_every, // Show progress after every n iterations.
  int       *num_seqs       // Number of sequences read.
) {
  int num_segments = 0; // Number of sequence segments processed.
  int num_matches = 0;  // Number of matches found.
  int start_pos = 0;    // Search from here.

  // Allocate memory.
  complete_match = allocate_match();
  partial_match = allocate_match();

  // Set up for computing score distribution.
  if (score_set == NULL) {
    score_set = set_up_score_set(
      p_threshold, 
      dp_threshold, 
      max_gap, 
      negatives_only,
      the_log_hmm);
  }

  *num_seqs = 0;

  while (
    (seq_file && read_one_fasta_segment(seq_file, max_chars, &sequence)) 
    || (sequence != NULL)
  ) {
    // Number of motifs x sequence length.
    MATRIX_T* motif_score_matrix = NULL;

    // Keep track of total database size for E-value calculation.
    score_set->total_length += get_seq_length(sequence);

    // Let the user know what's going on.
    if (verbosity > NORMAL_VERBOSE) {
      fprintf(stderr, "Scoring %s (length=%d) at position %d.\n",
        get_seq_name(sequence), get_seq_length(sequence),
        start_pos);
    }

    // Convert the sequence to alphabet-specific indices. 
    prepare_sequence(sequence);

    /* Allocate the dynamic programming matrix. Rows correspond to
       states in the model, columns to positions in the sequence. */
    if ((dp_rows < the_log_hmm->num_states) 
      || (dp_cols < get_seq_length(sequence))
    ) {
      free_matrix(dp_matrix);
      free_matrix(trace_matrix);
      if (dp_rows < the_log_hmm->num_states)
        dp_rows = the_log_hmm->num_states;
      if (dp_cols < get_seq_length(sequence))
        dp_cols = get_seq_length(sequence);
      // (Add one column for repeated match algorithm.)
      dp_matrix = allocate_matrix(dp_rows, dp_cols + 1);
      trace_matrix = allocate_matrix(dp_rows, dp_cols + 1);
    }

    // Compute the motif scoring matrix.
    if (motif_scoring) {
      motif_score_matrix = allocate_matrix(
        the_log_hmm->num_motifs,
        get_seq_length(sequence)
      );
      compute_motif_score_matrix(
        use_pvalues,
        p_threshold,
        get_int_sequence(sequence),
        get_seq_length(sequence),
        the_log_hmm,
        &motif_score_matrix
      );
    }

    // Fill in the DP matrix.
    repeated_match_algorithm(
      dp_threshold,
      get_int_sequence(sequence),
      get_seq_length(sequence),
      the_log_hmm,
      motif_score_matrix,
      dp_matrix,
      trace_matrix,
      complete_match
    );

    // Find all matches in the matrix.
    while (find_next_match(
      is_complete(sequence),
      start_pos,
      align_width,
      dp_matrix, 
      the_log_hmm,
      complete_match,
      partial_match)
    ) {

      // If this match starts in the overlap region, put if off until
      // the next segment.
      if (!is_complete(sequence) 
        && (get_start_match(partial_match) 
          > (get_seq_length(sequence) - OVERLAP_SIZE))
      ) {
        break;
      }

      // If this match starts at the beginning, then it's part of a
      // longer match from the previous matrix.  Skip it.  
      if (get_start_match(partial_match) == start_pos) {
        // fprintf(stderr, "Skipping match at position %d.\n", start_pos);

        // Change the starting position so that next call to find_next_match
        // will start to right of last match.
        start_pos = get_end_match(partial_match) + 1;
        continue;
      }
      start_pos = get_end_match(partial_match) + 1;

      // Store the score, ID, length and comment for later printing. 
      store_sequence(
        motif_scoring, 
        TRUE, // mhmmscan format 
        output_width,
        sequence,
        partial_match, 
        e_threshold,
        dp_threshold,
        p_threshold,
        score_set,
        got_evd,
        evd_set,
        print_fancy,
        the_log_hmm,
        motif_score_matrix,
        store_gff,
        NULL
      );

      /* Calculate the initial E-value distribution if the required
       * number sequences has been saved.  This will allow the
       * descriptions of low-scoring sequences not to be stored.  The
       * distribution will be recomputed using all scores when all
       * sequences have been read.  */
      if (score_set->n == EVD_NUM_SEQS && got_evd == FALSE) {
        evd_set = calc_distr(
          *score_set,   // Set of scores.  
          use_pvalues,  // Use exponential distribution?
          TRUE          // Use match E-values.
        );
        if (evd_set.n > 0) {
          got_evd = TRUE;
        }
      }

      num_matches++;
      user_update(
        num_matches, 
        *num_seqs, 
        num_segments, 
        num_matches, 
        progress_every
      );
    }

    // Does the input file contain more of this sequence?
    if (!is_complete(sequence)) { 
      int size_to_remove;

      // Compute the width of the buffer.
      size_to_remove = get_seq_length(sequence) - (OVERLAP_SIZE + align_width);

      // Adjust the sequence accordingly.
      remove_flanking_xs(sequence);
      shift_sequence(size_to_remove, sequence);
      if (SHIFT_DEBUG) {
        fprintf(stderr, "Retained %d bases.\n", get_seq_length(sequence));
      }

      // Remove the left-over portion from the total database size.
      score_set->total_length -= get_seq_length(sequence);

      // Next time, start looking for a match after the buffer.
      start_pos -= size_to_remove;
      if (start_pos < 0) {
        start_pos = 0;
      }

    } else {

      // Free the memory used by this sequence. 
      free_seq(sequence);
      sequence = NULL;
      start_pos = 0;

      *num_seqs += 1;
      user_update(*num_seqs, *num_seqs, num_segments, num_matches, 
      progress_every);
    }

    // Free the motif score matrix.
    free_matrix(motif_score_matrix);
    num_segments++;
    user_update(
      num_segments, 
      *num_seqs, 
      num_segments, 
      num_matches,
      progress_every
    );

    if (seq_file == NULL) {
      break;
    }
  }

  //
  // Free structures.
  //
  if (0) { // FIXME
    free_matrix(dp_matrix);
    free_matrix(trace_matrix);
    free_match(complete_match);
    free_match(partial_match);
  }

  // Return the set of scores.
  return(score_set);
} // read_and_score

/**************************************************************************
 * int main
 **************************************************************************/
int main (int argc, char * const argv[]) {

  // Command line parameters. 
  char*     hmm_filename;        // File containing the HMM. 
  FILE*     hmm_file;
  char*     seq_filename;        // File containing the sequences. 
  FILE*     seq_file;
  char*     gff_filename;        // GFF output file.
  FILE*     gff_file;
  PROB_T    p_threshold;         // P-value threshold for motif hits. 
  int       maxseqs;             // Maximum number of sequences to print
  double    max_gap;             // Maximum gap length to allow in matches.
  double    egcost;              // Fraction of E[score] that E[gap] should cost.
  PROB_T    dp_threshold;        // Score threshold for DP.
  BOOLEAN_T motif_scoring;       // Perform motif-scoring.
  BOOLEAN_T both_strands;        // Score both DNA strands.
  PROB_T    e_threshold;         // E-value threshold for scores. 
  char*     bg_filename;         // File containing background frequencies. 
  int       pam_distance;        // PAM distance 
  char*     sc_filename;         // File containing substitution scores. 
  double    beta;                // Weight on pseudocounts.
  BOOLEAN_T allow_weak_motifs;   // Allow motifs with min p-value > p-thresh?
  BOOLEAN_T zero_spacer_emit_lo; // Set spacer emission log-odds = 0?
  double    gap_open;            // Cost to open a gap; ignore if < 0.
  double    gap_extend;          // Cost to extend a gap; ignore if < 0.
  int       output_width;        // Width of output, in chars.
  int       max_chars;           // Number of chars to read at once.
  BOOLEAN_T print_fancy;         // Print Viterbi alignments? 
  BOOLEAN_T sort_output;         // Sort output scores? 
  int       progress_every;      // Show progress after every n iterations.
  BOOLEAN_T print_header;        // Print header information? 
  BOOLEAN_T print_params;        // Print program parameters? 
  BOOLEAN_T print_time;          // Print timing info? 
  BOOLEAN_T no_synth;            // Do not use synthetic scores for distribution.
  BOOLEAN_T text_output;         // Produce text output?

  // Data structures. 
  MHMM_T*   the_hmm;                 // The HMM itself. 
  MHMM_T*   the_log_hmm;             // The HMM, with probs converted to logs. 
  EVD_SET   evd_set;                 // EVD data structure 
  BOOLEAN_T got_evd = FALSE;         // no EVD found yet 
  BOOLEAN_T use_pvalues = FALSE;     // Don't use p-value scoring.
  SCORE_SET *score_set = NULL;       // Set of scores for computing distribution.
  SCORE_SET *synth_score_set = NULL; // Synthetic sequence scores.

  // Local variables. 
  double    start_time;  // Time at start of sequence processing. 
  double    end_time;    // Time at end of sequence processing. 
  int       num_seqs;    // Number of sequences processed. 
  int       align_width; // Width of one alignment line.
  int       dummy;       // Throw away variable.
  BOOLEAN_T use_synth;   // Distribution based on synthetic scores.
  FILE*     out_stream;  // Output stream (possibly running mhmm2html).

  // Record CPU time. 
  myclock();

  /***********************************************
   * Parse the command line.
   ***********************************************/

  // Set defaults. 
  hmm_filename = NULL;
  seq_filename = NULL;
  gff_filename = NULL;
  max_gap = -1;    // No maximum gap length.
  egcost = 0;    // No egcost given.
  p_threshold = -1;  // Don't do p-value scoring.
  dp_threshold = -BIG;
  e_threshold = DEFAULT_E_THRESHOLD;
  maxseqs = NO_MAX_SEQS;
  bg_filename = NULL;
  pam_distance = -1;  // Use default PAM distance.
  sc_filename = NULL;
  beta = -1;    // illegal beta, use defaults
  motif_scoring = FALSE;
  both_strands = FALSE; // Score given DNA strand only.
  allow_weak_motifs = FALSE;  // Don't allow weak motifs.
  zero_spacer_emit_lo = FALSE;
  gap_open = -1;   // No gap open penalty.
  gap_extend = -1;   // No gap extension penalty.
  output_width = DEFAULT_OUTPUT_WIDTH;
  max_chars = -1;
  print_fancy = FALSE;
  sort_output = TRUE;
  progress_every = DEFAULT_PROGRESS_EVERY;
  print_header = TRUE;
  print_params = TRUE;
  print_time = TRUE;
  text_output = FALSE;
  use_synth = FALSE;    // Assume there will be enough scores.
  no_synth = TRUE;    // Disallow generation of synthetic sequences.

  {
    // Define command line options.
    cmdoption const options[] = {
      {"gff", REQUIRED_VALUE},
      {"p-thresh", REQUIRED_VALUE},
      {"maxseqs", REQUIRED_VALUE},
      {"max-gap", REQUIRED_VALUE},
      {"e-thresh", REQUIRED_VALUE},
      {"fancy", NO_VALUE},
      {"width", REQUIRED_VALUE},
      {"text", NO_VALUE},
      {"nosort", NO_VALUE},
      {"bg-file", REQUIRED_VALUE},
      {"allow-weak-motifs", NO_VALUE},
      {"blocksize", REQUIRED_VALUE},
      {"synth", NO_VALUE},
      {"progress", REQUIRED_VALUE},
      {"verbosity", REQUIRED_VALUE},
      {"noheader", NO_VALUE},
      {"noparams", NO_VALUE},
      {"notime", NO_VALUE},
      {"quiet", NO_VALUE},
      {"zselo", NO_VALUE},
      {"gap-open", REQUIRED_VALUE},
      {"gap-extend", REQUIRED_VALUE},
      {"min-score", REQUIRED_VALUE},
      {"eg-cost", REQUIRED_VALUE},
      {"motif-scoring", NO_VALUE},
      {"pseudo-weight", REQUIRED_VALUE},
      {"pam", REQUIRED_VALUE},
      {"score-file", REQUIRED_VALUE}
    };
    int option_count = 27;
    int option_index = 0;

    // Define the usage message.
    char *usage = 
      "USAGE: mhmmscan [options] <HMM> <database>\n"
      "\n"
      "   Options:\n"
      "     --gff <file>\n"
      "     --maxseqs <int>\n"
      "     --p-thresh <p-value>\n"
      "     --max-gap <int>\n"
      "     --e-thresh <E-value> (default=10.0)\n"
      "     --fancy\n"
      "     --width <int> (default=79)\n"
      "     --text\n"
      "     --nosort\n"
      "     --bg-file <file>\n"
      "     --allow-weak-motifs\n"
      "     --blocksize <int>\n"
      "     --synth\n"
      "     --progress <int>\n"
      "     --verbosity 1|2|3|4|5 (default=2)\n"
      "     --noheader\n"
      "     --noparams\n"
      "     --notime\n"
      "     --quiet\n"
      "\n"
      "   Advanced options:\n"
      "     --zselo\n"
      "     --gap-open <cost>\n"
      "     --gap-extend <cost>\n"
      "     --min-score <score>\n"
      "     --eg-cost <fraction>\n"
      "     --motif-scoring\n"
      "     --pseudo-weight <weight> (default=10)\n"
      "     --pam <distance> (default=250 [protein] 1 [DNA])\n"
      "     --score-file <file>\n"
      "\n";

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
        die("Error process command line options (%s)\n", message);
      }

      if (strcmp(option_name, "gff") == 0) {
        gff_filename = option_value;
      } else if (strcmp(option_name, "p-thresh") == 0) {
        p_threshold = atof(option_value);
      } else if (strcmp(option_name, "max-gap") == 0) {
        max_gap = atof(option_value);
      } else if (strcmp(option_name, "e-thresh") == 0) {
        e_threshold = atof(option_value);
      } else if (strcmp(option_name, "maxseqs") == 0) {
        maxseqs = atoi(option_value);
      } else if (strcmp(option_name, "fancy") == 0) {
        print_fancy = TRUE;
      } else if (strcmp(option_name, "width") == 0) {
        output_width = atoi(option_value);
      } else if (strcmp(option_name, "text") == 0) {
        text_output = TRUE;
      } else if (strcmp(option_name, "nosort") == 0) {
        sort_output = FALSE;
      } else if (strcmp(option_name, "bg-file") == 0) {
        bg_filename = option_value;
      } else if (strcmp(option_name, "allow-weak-motifs") == 0) {
        allow_weak_motifs = TRUE;
      } else if (strcmp(option_name, "blocksize") == 0) {
        max_chars = atoi(option_value);
      } else if (strcmp(option_name, "synth") == 0) {
        no_synth = FALSE;
      } else if (strcmp(option_name, "progress") == 0) {
        progress_every = atoi(option_value);
      } else if (strcmp(option_name, "verbosity") == 0) {
        verbosity = (VERBOSE_T)atoi(option_value);
      } else if (strcmp(option_name, "noheader") == 0) {
        print_header = FALSE;
      } else if (strcmp(option_name, "noparams") == 0) {
        print_params = FALSE;
      } else if (strcmp(option_name, "notime") == 0) {
        print_time = FALSE;
      } else if (strcmp(option_name, "quiet") == 0) {
        print_header = print_params = print_time = FALSE;
        verbosity = QUIET_VERBOSE;
      } else if (strcmp(option_name, "zselo") == 0) {
        zero_spacer_emit_lo = TRUE;
      } else if (strcmp(option_name, "gap-open") == 0) {
        gap_open = atof(option_value);
      } else if (strcmp(option_name, "gap-extend") == 0) {
        gap_extend = atof(option_value);
      } else if (strcmp(option_name, "min-score") == 0) {
        dp_threshold = atof(option_value);
      } else if (strcmp(option_name, "eg-cost") == 0) {
        egcost = atof(option_value);
      } else if (strcmp(option_name, "motif-scoring") == 0) {
        motif_scoring = TRUE;
      } else if (strcmp(option_name, "pseudo-weight") == 0) {
        beta = atof(option_value);
      } else if (strcmp(option_name, "pam") == 0) {
        pam_distance = atoi(option_value);
      } else if (strcmp(option_name, "score-file") == 0) {
        sc_filename = option_value;
      } 
    }

    // Read the two required arguments.
    if (option_index + 2 != argc) {
      fprintf(stderr, "%s", usage);
      exit(1);
    }
    hmm_filename = argv[option_index];
    seq_filename = argv[option_index+1];
  }
     
  /***********************************************
   * Verify the command line.
   ***********************************************/

  // Make sure we got the required files. 
  if (hmm_filename == NULL) {
    die("No HMM file given.\n");
  }
  if (seq_filename == NULL) {
    die("No sequence file given.\n");
  }

  // Force p-value scoring if p-value threshold given.
  if (p_threshold != -1) {
    use_pvalues = TRUE;
  }

  // Check p-threshold is in range [0<p<=1].
  if (use_pvalues && (p_threshold <= 0 || p_threshold > 1)) {
    die("You may only specify p-thresh in the range [0<p<=1].\n");
  }

  // Check that no dp_threshold given if max-gap given.
  if (max_gap != -1 && dp_threshold != -BIG) {
    die("You may not specify both max-gap and min-score\n");
  }

  if (egcost > 0 && dp_threshold != -BIG) {
    die("You may not specify both egcost and min-score.\n");
  }

  // Set dp_threshold to default if none given by user.
  if (dp_threshold == -BIG) {
    dp_threshold = DEFAULT_DP_THRESHOLD;
  }

  // Check that min-score is positive.
  if (dp_threshold <= 0) die("You may only specify a positive min-score.\n");

  // Force motif_scoring if using p-value scoring.
  if (use_pvalues) {
    motif_scoring = TRUE;
  }

  // FIXME: both-strands not implemented
  if (both_strands) {
    die("Sorry, -both-strands not yet implemented.");
  }

  // Check egcost ok.
  if (egcost < 0) {
    die("You may only specify a positive egcost.\n");
  }
  if (egcost > 0 && max_gap == -1) {
    die("You may only specify egcost if you specify max-gap.\n");
  }

  // Check that -bg-file given if -synth given
  if (!no_synth && bg_filename == NULL) {
    die("You may only use -synth if you specify -bg-file.\n");
  }

  // Compute the alignment width.
  align_width = compute_align_width(output_width);

  /***********************************************
   * Open files.
   ***********************************************/

  if (open_file(seq_filename, "r", TRUE, "sequence", "sequences", &seq_file) 
    == 0) {
    exit(1);
  }

  if (gff_filename != NULL) {
    if (open_file(gff_filename, "w", FALSE, "gff", "gff", &gff_file) == 0) {
      exit(1);
    }
  } else {
    gff_file = NULL;
  }


  /***********************************************
   * Set up the model.
   ***********************************************/

  // Read the model. 
  read_mhmm(hmm_filename, &the_hmm);

  //
  // Check gap switches and set dp_threshold.
  //
  if (max_gap != -1) {
    // Don't allow gap-open, gap-extend or min-score if max-gap given.
    if ((gap_open != -1) || (gap_extend != -1) 
      || (dp_threshold != DEFAULT_DP_THRESHOLD)
    ) {
      die("You may not specify gap-open, gap-extend or min-score AND max-gap.\n");
    } 
    // Check that length is legal.
    if (max_gap < 0) {
      die("You may not specify a negative -max-gap length.\n");
    }

    // 
    // Set dp_threshold.
    // If using egcost, get model statistics in order to set it.
    //
    if (egcost > 0) {    // Gap cost a fraction of expected hit score.
      score_set = set_up_score_set(
        p_threshold, 
        dp_threshold, 
        max_gap, 
        FALSE,
        the_hmm
      );
      dp_threshold = 
        (egcost * max_gap * score_set->e_hit_score) / score_set->egap;
    } else {      // (Approximately) zero gap costs.
      if (dp_threshold == DEFAULT_DP_THRESHOLD) {
        dp_threshold = 1e-6;  // Very small number.
      }
    } // Check gap switches and set dp_threshold.

    //
    // Set gap costs.
    //
    gap_open = gap_extend = dp_threshold/max_gap;
    zero_spacer_emit_lo = TRUE;
  }

  //
  // Prepare the model for recognition.
  //
  if (pam_distance == -1) {
    pam_distance = (which_alphabet() == PROTEIN_ALPH) 
      ? DEFAULT_PROTEIN_PAM : DEFAULT_DNA_PAM;
  }
  if (beta < 0)  {
    beta = (which_alphabet() == PROTEIN_ALPH)
      ? DEFAULT_PROTEIN_BETA : DEFAULT_DNA_BETA;
  }
  free_array(the_hmm->background);
  the_hmm->background = get_background(bg_filename);
  convert_to_from_log_hmm(
    TRUE, // Convert to logs.
    zero_spacer_emit_lo,
    gap_open,
    gap_extend,
    the_hmm->background,
    sc_filename,
    pam_distance,
    beta,
    the_hmm, 
    &the_log_hmm
  );

  // Set up PSSM matrices if doing motif_scoring
  // and pre-compute motif p-values if using p-values.
  // Set up the hot_states list.
  set_up_pssms_and_pvalues(
    motif_scoring,
    p_threshold,
    both_strands,
    allow_weak_motifs,
    the_log_hmm
  );

  // Decide how many characters to read.
  if (max_chars == -1) {
    max_chars = (int)(MAX_MATRIX_SIZE / the_log_hmm->num_states);
  }
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(
      stderr, 
      "Reading sequences in blocks of %d characters.\n", 
      max_chars
    );
  }

  start_time = myclock();

  //
  // Read and score the real sequences.
  //
  score_set = read_and_score(
    NULL,      // No score set yet.
    seq_file,
    NULL,      // No sequence.
    FALSE,      // Positives and negatives.
    max_chars,
    max_gap,
    the_log_hmm,
    motif_scoring,
    use_pvalues,
    p_threshold,
    dp_threshold,
    e_threshold,
    got_evd,
    evd_set,
    print_fancy,
    gff_file != NULL,
    output_width,
    align_width,
    progress_every, 
    &num_seqs
  );

  //
  // Generate synthetic sequences and score them if fewer than
  // MIN_MATCHES matches was found.  Don't do this if no matches
  // were found (of course).
  //
  // Each sequence is 100000bp long.
  // Continues to generate sequences until enough matches are found
  // so that the standard error of mu1 will be 5% (std err = 100 * 1/sqrt(n))
  // or the maximum database size is reached.
  //
  if (!no_synth && score_set->n>0 && score_set->n<MIN_MATCHES) {  
    int i;
    int want = 1e3;  // Desired number of matches (for mu1 std err= 3%).
    double need = MIN_SCORES;  // Minimum number of matches (for mu1 std err= 5%).
    int found = 0;  // Matches found.
    double maxbp = 1e7;  // Maximum number of bp to generate.
    int db_size = 0;  // Number of bp generated. 
    int len = 1e5;  // Desired sequence length.
    int nseqs = 0;  // Number of sequences generated.
    SEQ_T *seq;    // Save time by not using tmpfile.
    double min_gc, max_gc;   // Min and Max observed GC-contents.
    double f[4];  // 0-order Markov model

    use_synth = TRUE;  // Distribution based on synthetic scores.

    // Get the minimum and maximum GC-contents of match regions.
    min_gc = 1;
    max_gc = 0;
    for (i=0; i<score_set->n; i++) {
      double gc = score_set->scores[i].gc;
      if (gc < min_gc) min_gc = gc;
      if (gc > max_gc) max_gc = gc;
    }

    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(stderr, "Generating and scoring synthetic sequences...\n");
    }

    // Loop until enough matches found or we get tired.
    while (found < want && db_size < maxbp) {  
      FILE* synth_seq_file = tmpfile();    // Use file not seq variable.
      BOOLEAN_T gc_stratify = (MAX_RATIO < 100);
      if (gc_stratify) {   // Select a GC content randomly and set f[].
        double gc = min_gc + (drand48() * (max_gc-min_gc));
        f[0] = f[3] = (1 - gc)/2;
        f[1] = f[2] = gc/2;
        fprintf(stderr, "gc = %f\r", gc);
      }
      seq = gendb(
        synth_seq_file, // output
        (which_alphabet() == DNA_ALPH) ? 3 : 4,  // type
        (gc_stratify) ? NULL : bg_filename,  // Markov model or f?
        -1, // use order in bg_file
        f, // 0-order model.
        1, // # of sequences
        len, // min length
        len, // max length
        0 // no seed
      );
      nseqs++;          // seed
      db_size += len;        // Size of db so far.
      if (synth_seq_file) {
        rewind(synth_seq_file);
      }
      if (verbosity > NORMAL_VERBOSE) {
        fprintf(stderr, "seqs: %d size: %d\r", nseqs, db_size);
      }

      //
      // Read and score synthetic sequences.
      //
      synth_score_set = read_and_score(
        synth_score_set,
        synth_seq_file,
        seq,
        TRUE, // Negatives only.
        max_chars,
        max_gap,
        the_log_hmm,
        motif_scoring,
        use_pvalues,
        p_threshold,
        dp_threshold,
        -1, // e_threshold: save nothing
        got_evd,
        evd_set,
        FALSE,
        FALSE,
        output_width,
        align_width,
        progress_every, 
        &dummy
      );
      found = synth_score_set->n;
      fprintf(stderr, "seqs: %d size: %d matches: %d\r", nseqs, db_size, found);

      // Close file to delete sequences.
      if (synth_seq_file) {
        fclose(synth_seq_file);
      }
      // free_seq(seq);      // free'd by read_and_score

      // Quit if we're not getting there.
      if (found/need < db_size/maxbp) {
        fprintf(
          stderr, 
          "\nGiving up generating synthetic sequences: match probability too low!\n"
        );
        break;
      }

    } // Generate synthetic scores.

    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(
        stderr, 
        "Generated %d synthetic sequences (%d characters) with %d matches.\n",
        nseqs, 
        db_size, 
        found
      );
    }

  } else {
    synth_score_set = score_set;    // Use actual scores for estimation.
  } // Generate and score synthetic sequences.

  end_time = myclock();

  /***********************************************
   * Calculate the E-values and store them as the keys.
   ***********************************************/
  // Recalculate the score distribution using synthetic scores. 
  // If successful, calculate the E-values and store them as keys.
  evd_set = calc_distr(
    *synth_score_set, // Set of scores.
    use_pvalues,  // Use exponential distribution?
    TRUE // Use match E-values.
  );
  // Distribution based on synthetic scores?
  evd_set.negatives_only = use_synth; 
  if (evd_set.n >= 1) {      // Found a valid score distribution.
    int q, t, N;
    q = 1; // Ignore query length.
    // Get p-value multiplier.
    if (use_synth) {
      N = evd_set.non_outliers = get_n(*score_set, evd_set);
    } else {
      N = evd_set.non_outliers; // Use number of non-outliers.
    }
    // p-value multiplier for E-value.
    evd_set.N = N;      
    // Record number of real scores in evd_set.
    evd_set.nscores = score_set->n;  

    // Get sequence length.
    if (use_pvalues) { // Exponential.
      t = 0; // Ignore sequence length.
    } else { // EVD
      t = score_set->total_length/score_set->n;  // average length
    }

    calc_evalues(&evd_set, N, q, t);
    got_evd = TRUE;
  }

  /***********************************************
   * Start the mhmm2html process, if requested.
   ***********************************************/
  if (text_output) {
    out_stream = stdout;
  } else {
    out_stream = open_command_pipe(
      "mhmm2html", // Program
      strip_filename(argv[0]), // Search directory
      "-test -", // Test arguments.
      "mhmm2html", // Expected reply.
      "-", // Real arguments.
      TRUE, // Return stdout on failure.
      "Warning: mhmm2html not found.  Producing text output."
    );
  }

  /***********************************************
   * Print header information.
   ***********************************************/
  if (print_header) {
    write_header(
      "mhmmscan",
      "Database search results",
      the_hmm->description,
      the_hmm->motif_file,
      hmm_filename,
      seq_filename,
      out_stream
    );
    if (gff_file != NULL) {
      write_header(
        "mhmmscan",
        "Database search results",
        the_hmm->description,
        the_hmm->motif_file,
        hmm_filename,
        seq_filename,
        gff_file
      );
    }
  }

  /***********************************************
   * Sort and print the results.
   ***********************************************/
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "\nSorting the scores.\n");
  }
  sort_and_print_scores(
    print_fancy, 
    print_header, 
    got_evd,
    motif_scoring,
    TRUE, // Print in mhmmscan format.
    maxseqs, // Maximum number of sequences to print
    output_width,
    e_threshold,
    sort_output,
    gff_file,
    out_stream
  );

  if (print_params) {
    print_parameters(
      argv, 
      argc,
      "mhmmscan",
      hmm_filename,
      seq_filename,
      TRUE, // Viterbi search.
      dp_threshold,
      motif_scoring,
      use_pvalues,
      p_threshold,
      e_threshold,
      both_strands,
      bg_filename,
      sc_filename,
      pam_distance,
      beta,
      zero_spacer_emit_lo,
      max_gap,
      egcost,
      gap_open,
      gap_extend,
      print_fancy,
      print_time,
      start_time,
      end_time,
      num_seqs,
      evd_set,
      *score_set,
      out_stream
    );
  }

  // Tie up loose ends.
  free_mhmm(the_hmm);
  free_mhmm(the_log_hmm);
  fclose(seq_file);
  if (gff_filename != NULL) {
    fclose(gff_file);
  }    
  pclose(out_stream);

  return(0);
}
