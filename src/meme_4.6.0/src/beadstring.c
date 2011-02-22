/********************************************************************
 * FILE: beadstring.c
 * AUTHOR: William Stafford Noble Charles E. Grant
 * CREATE DATE: 06-23-2007
 * PROJECT: MHMM
 * COPYRIGHT: 2007, WNG
 * DESCRIPTION: Given MEME output, create a motif-based linear HMM and
 *              use it to search a sequence database for motif clusters.
 ********************************************************************/
#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "matrix.h"         // Matrix functions
#include "alphabet.h"       // Amino acid and nucleotide alphabets
#include "beadstring.h"
#include "beadstring-xml.h" // Functions for printing Beadstring XML output
#include "build-hmm.h"      // Construct an HMM in memory.
#include "cisml.h"          // Structure for recording results
#include "dir.h"            // PATHS to needed directories
#include "dp.h"             // Dynamic programming algorithms
#include "fasta-io.h"       // Reading and writing FASTA files
#include "fitevd.h"         // Functions for fitting an extreme value dist
#include "io.h"             // File and directory utilites
#include "log-hmm.h"        // Convert an HMM to or from log form.
#include "match.h"          // Find sequence to model match from traceback
#include "metameme.h"       // Global Meta-MEME functions.
#include "meme-io.h"        // Read motifs from MEME output.
#include "mhmm.h"           // Scan a sequence db with HMM.
#include "mhmm-state.h"     // Data structure for HMMs.
#include "mhmms.h"          // Scan a sequence db with HMM.
#include "motif.h"          // Data structure for holding motifs.
#include "order.h"          // Linear motif order and spacing.
#include "projrel.h"
#include "pssm.h"           // Tools for position specific scoring matrix
#include "read-mhmm.h"      // Functions for reading HMM from file.
#include "simple-getopt.h"  // Command line processing
#include "string-list.h"    // List of strings (for motif IDs).
#include "utils.h"          // Generic utilities.
#include "write-mhmm.h"     // Functions for producing MHMM output.

#define SLOP 1E-5

// Marker for end transition.
extern char* END_TRANSITION;

VERBOSE_T verbosity = NORMAL_VERBOSE;

enum path_values { SINGLE, ALL };

/************************************************************************
 * Fix list of transitions in
 ************************************************************************/
void fix_trans_in(MHMM_T* the_hmm) {

  int i_row, i_col;
  int n = 0;
  MATRIX_T *trans = the_hmm->trans;
  n = trans->num_rows;

  //
  // Visit the transition matrix cells just once each
  // to update trans_in arrays.
  // This is quadratic in n.
  //
  for (i_col = 0; i_col < n; i_col++) {
    for (i_row = 0; i_row < n; i_row++) {
      double p = get_matrix_cell(i_row, i_col, trans);	// The transition probability.
      int m = 0;
      if (!is_zero(p, FALSE)) {
        MHMM_STATE_T * in_state = &(the_hmm->states[i_col]);
        assert(m < in_state->ntrans_in);
        set_array_item(m, p, in_state->trans_in);
        m++;
      }
    } // col
  } // row

} // compute_ins_and_outs

/*************************************************************************
 *  Build a linear HMM model for the motifs
 *************************************************************************/
MHMM_T* build_linear_model(struct options* options) {

  int num_motifs = 0;              // The number of motifs in the model.
  MOTIF_T  motifs[2 * MAX_MOTIFS]; // The motifs.
  BOOLEAN_T has_reverse_strand = FALSE; // MEME file contained both strands
  ARRAY_T* background = NULL;      // Background probs for alphabet.
  ORDER_T* order_spacing = NULL;   // Linear HMM order and spacing.
  MATRIX_T* transp_freq = NULL;    // Matrix of inter-motif transitions freqs.
  MATRIX_T* spacer_ave = NULL;     // Matrix of average spacer lengths.
  MHMM_T* the_hmm = NULL;          // The HMM being constructed.
  STRING_LIST_T* motif_occurrences = NULL; /* Strings describing occurrences
                                             of motifs */

  // Gotta have positive number of spacer states.
  if (options->spacer_states <= 0) {
    die("Negative spacer length (%d).\n", options->spacer_states);
  }

  // Make sure motifs weren't selected redundantly.
  if ((get_num_strings(options->requested_motifs) != 0)
      && (options->order_string != NULL)) {
    die("Can't combine --motif and --order options.");
  }

  // Parse the order string.
  order_spacing = create_order(options->order_string);

  /**********************************************
   * READING THE MOTIFS
   **********************************************/
  read_meme_file(
		 options->meme_filename,
		 options->bg_filename,
		 options->motif_pseudo,
     REQUIRE_PSPM, //use pspm
		 &num_motifs,
		 motifs,
		 &motif_occurrences,
		 &has_reverse_strand,
		 &background
		 );

  // If motifs use protein alphabet we will not scan both strands
  ALPH_T alphabet_type = which_alphabet();
  if (alphabet_type == PROTEIN_ALPH) {
    options->alphabet = PROTEIN_ALPH;
  }
  else {
    options->alphabet = DNA_ALPH;
  }

  process_raw_motifs_for_model(
       &num_motifs,
       motifs,
       motif_occurrences,
       options->requested_motifs,
       has_reverse_strand,
       FALSE, // keep unused motifs
       options->motif_p_threshold,
       options->motif_e_threshold,
       0.0, // complexity threshold
       &order_spacing,
       &transp_freq,
       &spacer_ave,
       options->trans_pseudo,
       options->spacer_pseudo
  );

  /**********************************************
   * BUILDING THE HMM
   **********************************************/

  if (order_spacing != NULL) {
    reorder_motifs(order_spacing, &num_motifs, motifs);
  }
  else {
    die("No order specified for the motifs.\n"
        "For the linear model the motif file must contain motif occurence\n"
        "data or the motif order must be specified using the --order option.");
  }

  build_linear_hmm(
    background,
    order_spacing,
    options->spacer_states,
    motifs,
    num_motifs,
    options->fim,
    &the_hmm
  );
  fix_trans_in(the_hmm);

  free_string_list(motif_occurrences);
  free_array(background);
  free_order(order_spacing);
  free_matrix(transp_freq);
  free_matrix(spacer_ave);
  int i_motif = 0;
  for (i_motif = 0; i_motif < num_motifs; i_motif++) {
    free_motif(&(motifs[i_motif]));
  }

  // Add some global information.
  copy_string(&(the_hmm->motif_file), options->meme_filename);

  return the_hmm;

}


void search_with_model(
 MHMM_T* the_mhmm,
 BOOLEAN_T* got_evd,
 OPTIONS_T *options,
 PATTERN_T *pattern
) {

  BOOLEAN_T viterbi = FALSE;
  BOOLEAN_T use_pvalues = FALSE;
  BOOLEAN_T both_strands = FALSE;
  EVD_SET evd_set;      // EVD data structure
  evd_set.evds = NULL;  // Initialize pointer to evds.

  /***********************************************
   * Set up the model.
   ***********************************************/

  // Figure out what kind of scoring to do.
  if (options->paths == SINGLE) {
    viterbi = TRUE;
  } else {
    viterbi = FALSE;
  }

  // Force p-value scoring if p-value threshold given.
  if (options->p_threshold != 0) {
    use_pvalues = TRUE;
  }

  // Check p-threshold is in range [0<p<=1].
  if (use_pvalues && (options->p_threshold <= 0 || options->p_threshold > 1)) {
    die("You must specify a p-value threshold in the range [0<p<=1]\n");
  }

  // Set the PAM distance.
  if (options->pam_distance == -1) {
    options->pam_distance = (which_alphabet() == PROTEIN_ALPH)
      ? DEFAULT_PROTEIN_PAM : DEFAULT_DNA_PAM;
  }
  int beta = (which_alphabet() == PROTEIN_ALPH) ? 10 : 1;

  // Convert the model to log space.
  MHMM_T*   the_log_hmm;      // The HMM, with probs converted to logs.
  convert_to_from_log_hmm(
    TRUE, // Convert to logs.
    options->zselo,
    options->gap_open,
    options->gap_extend,
    the_mhmm->background,
    options->score_filename,
    options->pam_distance,
    beta,
    the_mhmm,
    &the_log_hmm
  );

  // Set up PSSM matrices if doing motif_scoring
  // and pre-compute motif p-values if using p-values.
  // Set up the hot_states list.
  set_up_pssms_and_pvalues(
    options->motif_scoring,
    options->p_threshold,
    both_strands,
    options->allow_weak_motifs,
    the_log_hmm
  );

  //
  // Set up for computing score distribution.
  //
  SCORE_SET *score_set = set_up_score_set(
    options->p_threshold,
    -1,    // dp_threshold,
    -1,    // max_gap
    FALSE, // negatives_only
    the_log_hmm
  );

  /***********************************************
   * Search the database one sequence at a time.
  ***********************************************/

  // Open the file for reading.
  FILE* seq_file = NULL;
  if (open_file(
        options->seq_filename,
        "r",
        TRUE,
        "sequence",
        "sequences",
        &seq_file
      ) == 0) {
    exit(1);
  }

  double start_time = myclock();
  double old_loop_time = start_time;
  double seconds_since_last_report = 0.0;

  SEQ_T* sequence = NULL;               // Sequence to search against.
  MATRIX_T* motif_score_matrix = NULL;  // Number of motifs x sequence length.
  MATRIX_T* dp_matrix = NULL;           // Dynamic programming matrix.
  MATRIX_T* trace_matrix = NULL;        // Traceback for Viterbi.
  MATCH_T*  this_match = NULL;          // This sequence-to-model match.
  int dp_rows = 0;                      // Size of the DP matrix.
  int dp_cols = 0;
  int num_seqs = 0;

  this_match = allocate_match(); // This sequence-to-model match.
  while (read_one_fasta(seq_file, MAX_SEQ, &sequence)) {

    num_seqs++;

    // Create a scanned_sequence record and record it in pattern.
    char* seq_name = get_seq_name(sequence);
    int seq_length = get_seq_length(sequence);
    SCANNED_SEQUENCE_T *scanned_seq =
      allocate_scanned_sequence(seq_name, seq_name, pattern);
    set_scanned_sequence_length(scanned_seq, seq_length);
    // Keep track of total database size for E-value calculation.
    score_set->total_length += seq_length;

    // Let the user know what's going on.
    if (verbosity > NORMAL_VERBOSE) {
      fprintf(
        stderr,
        "Scoring %s (length=%d).\n",
        seq_name,
        seq_length
      );
    }

    // Convert the sequence to alphabet-specific indices.
    prepare_sequence(sequence);
    // Sequeence length may have changed
    seq_length = get_seq_length(sequence);
    assert(get_seq_char(seq_length, sequence) == '\0');

    /* Allocate the dynamic programming matrix. Rows correspond to
       states in the model, columns to positions in the sequence. */
    if ((dp_rows < the_log_hmm->num_states)
      || (dp_cols < seq_length)) {
      free_matrix(dp_matrix);
      free_matrix(trace_matrix);
      if (dp_rows < the_log_hmm->num_states) {
        dp_rows = the_log_hmm->num_states;
      }
      if (dp_cols < seq_length) {
        dp_cols = seq_length;
      }
      // (Add one column for repeated match algorithm.)
      dp_matrix = allocate_matrix(dp_rows, dp_cols + 1);
      trace_matrix = allocate_matrix(dp_rows, dp_cols + 1);
    }

    // Compute the motif scoring matrix.
    if (options->motif_scoring) {
      motif_score_matrix = allocate_matrix(
        the_log_hmm->num_motifs,
        seq_length
      );
      compute_motif_score_matrix(
        use_pvalues,
        options->p_threshold,
        get_int_sequence(sequence),
        seq_length,
        the_log_hmm,
        &motif_score_matrix
      );
    }

    // Perform the appropriate type of dynamic programming.
    if (viterbi) {
      viterbi_algorithm(
        options->local_scoring,
        get_int_sequence(sequence),
        seq_length,
        the_log_hmm,
        TRUE,  // Save Viterbi path?
        motif_score_matrix,
        dp_matrix,
        trace_matrix,
        this_match
      );
    } else {
      forward_algorithm(
        options->local_scoring,
        get_int_sequence(sequence),
        seq_length,
        the_log_hmm,
        //motif_score_matrix, // FIXME
        dp_matrix,
        this_match
      );
    }

    // Store the score, ID, length and comment for later printing.
    //assert(get_match_seq_length(this_match) == get_seq_length(sequence));
    //assert(get_match_seq_length(this_match) == get_seq_length(sequence));
    store_sequence(
      options->motif_scoring,
      FALSE, // Don't print in mhmmscan format.
      DEFAULT_OUTPUT_WIDTH,
      sequence,
      this_match,
      options->e_threshold,
      0,    // no dp_threshold
      options->p_threshold,
      score_set,
      *got_evd,
      evd_set,
      FALSE, // print fancy
      the_log_hmm,
      motif_score_matrix,
      FALSE, // Store GFF?,
      scanned_seq
    );

    /* Calculate the initial E-value distribution if the required
     * number of sequences has been saved.  This will allow the
     * descriptions of low-scoring sequences to not be stored.  The
     * distribution will be recomputed using all saved scores when all
     * sequences have been read.  */
    if (score_set->n == EVD_NUM_SEQS && *got_evd == FALSE) {
      evd_set = calc_distr(
        *score_set,   // Set of scores.
        use_pvalues, // Use exponential distribution?
        FALSE         // Use sequence E-values.
      );
      if (evd_set.n > 0) *got_evd = TRUE;
    }

    // Free the memory used by this sequence.
    free_matrix(motif_score_matrix);
    free_seq(sequence);

    // If requested, periodically report progress
    if (options->progress > 0.0) {
      double new_loop_time = myclock();
      seconds_since_last_report = (new_loop_time - old_loop_time) / 1e6;
      if (seconds_since_last_report > options->progress) {
        double elapsed_seconds = (new_loop_time - start_time) / 1e6;
        old_loop_time = new_loop_time;
        fprintf(
          stderr,
          "%d sequences processed in %d seconds.\n",
          num_seqs,
          (int) elapsed_seconds
        );
      }
    }
  }

  //double end_time = myclock();

  /***********************************************
   * Calculate the E-values and store them as the keys.
   ***********************************************/
  // Recalculate the EVD using all scores.
  // If successful, calculate the E-values and
  // store them as keys.
  evd_set = calc_distr(
    *score_set,   // Set of scores.
    use_pvalues, // Use exponential distribution?
    FALSE         // Use sequence E-values.
  );
  if (evd_set.n > 0) {
    int q, t, N;
    q = 1;        // Ignore query "length".
    t = 0;        // Use stored target lengths.
    N = evd_set.non_outliers;    // Use number of non-outliers.
    evd_set.N = N;
    calc_evalues(&evd_set, N, q, t);
    *got_evd = TRUE;
  }
  evd_set.negatives_only = FALSE;  // Used real sequences.

  /* Tie up loose ends. */
  myfree(evd_set.evds);
  free_mhmm(the_log_hmm);
  free_matrix(dp_matrix);
  free_matrix(trace_matrix);
  free_match(this_match);
  free(score_set->scores);
  myfree(score_set);
  fclose(seq_file);
}

/***********************************************************************
 * Print beadstring specific information to an XML file
 ***********************************************************************/
static  void print_settings_xml(
  FILE *out,
  OPTIONS_T *options
) {
  fputs("<settings>\n",  out);
  fprintf(out, "<setting name=\"%s\">%s</setting>\n", "output directory", options->output_dirname);
  fprintf(out, "<setting name=\"%s\">%s</setting>\n", "MEME file name", options->meme_filename);
  fprintf(out, "<setting name=\"%s\">%s</setting>\n", "background file name", options->bg_filename);
  fprintf(out, "<setting name=\"%s\">%s</setting>\n", "score file name", options->score_filename);
  fprintf(out, "<setting name=\"%s\">%s</setting>\n", "sequence file name", options->seq_filename);
  fprintf(out, "<setting name=\"%s\">%s</setting>\n", "model file name", options->model_filename);
  fprintf(out, "<setting name=\"%s\">%s</setting>\n", "order string", options->order_string);
  fprintf(
    out,
    "<setting name=\"%s\">%s</setting>\n",
    "allow weak motifs",
    boolean_to_string(options->allow_weak_motifs)
  );
  fprintf(
    out,
    "<setting name=\"%s\">%s</setting>\n",
    "allow clobber",
    boolean_to_string(options->allow_clobber)
  );
  fprintf(
    out,
    "<setting name=\"%s\">%s</setting>\n",
    "fim",
    boolean_to_string(options->fim)
  );
  fprintf(
    out,
    "<setting name=\"%s\">%s</setting>\n",
    "local scoring",
    boolean_to_string(options->local_scoring)
  );
  fprintf(
    out,
    "<setting name=\"%s\">%s</setting>\n",
    "motif scoring",
    boolean_to_string(options->motif_scoring)
  );
  fprintf(
    out,
    "<setting name=\"%s\">%s</setting>\n",
    "synth",
    boolean_to_string(options->synth)
  );
  fprintf(
    out,
    "<setting name=\"%s\">%s</setting>\n",
    "zselo",
    boolean_to_string(options->zselo)
  );
  fprintf(out, "<setting name=\"%s\">%d</setting>\n", "max_seq", options->max_seq);
  fprintf(out, "<setting name=\"%s\">%d</setting>\n", "pam_distance", options->pam_distance);
  fprintf(out, "<setting name=\"%s\">%d</setting>\n", "paths", options->paths);
  fprintf(out, "<setting name=\"%s\">%d</setting>\n", "spacer_states", options->spacer_states);
  fprintf(out, "<setting name=\"%s\">%3.1g</setting>\n", "gap extend", options->gap_extend);
  fprintf(out, "<setting name=\"%s\">%3.1g</setting>\n", "gap open", options->gap_open);
  fprintf(out, "<setting name=\"%s\">%3.2g</setting>\n", "e-threshold", options->e_threshold);
  fprintf(out, "<setting name=\"%s\">%3.2g</setting>\n", "p-threshold", options->p_threshold);
  fprintf(
    out,
    "<setting name=\"%s\">%s</setting>\n",
    "perform search",
    boolean_to_string(options->perform_search)
  );
  fprintf(out, "<setting name=\"%s\">%3.2g</setting>\n", "progress", options->progress);
  fprintf(out, "<setting name=\"%s\">%3.2g</setting>\n", "motif e-threshold", options->motif_e_threshold);
  fprintf(out, "<setting name=\"%s\">%3.2g</setting>\n", "motif p-threshold", options->motif_p_threshold);
  fprintf(out, "<setting name=\"%s\">%3.2g</setting>\n", "motif pseudocount", options->motif_pseudo);
  fprintf(out, "<setting name=\"%s\">%3.2g</setting>\n", "spacer pseudocount", options->spacer_pseudo);
  fprintf(out, "<setting name=\"%s\">%3.2g</setting>\n", "transition pseudocount", options->trans_pseudo);
  fprintf(out, "<setting name=\"%s\">%d</setting>\n", "verbosity", verbosity);
  int i = 0;
  int num_strings = get_num_strings(options->requested_motifs);
  for(i = 0; i < num_strings; i++) {
    fprintf(
      out,
      "<setting name=\"%s\">%s</setting>\n", "requested motif",
      get_nth_string(i, options->requested_motifs)
    );
  }

  fputs("</settings>\n",  out);

}

/***********************************************************************
 * Print beadstring specific information to an XML file
 ***********************************************************************/
static void print_beadstring_xml_file(
  FILE *out,
  OPTIONS_T *options,
  ARRAY_T *background,
  char  *stylesheet
) {

  fputs("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n", out);
  if (stylesheet != NULL) {
    fprintf(out, "<?xml-stylesheet type=\"text/xsl\" href=\"%s\"?>\n", stylesheet);
  }
  fputs("<!-- Begin document body -->\n", out);

  const char *archive_date = ARCHIVE_DATE;
  int i = strlen(archive_date);
  fprintf(out,
    "<beadstring version=\"%s\" release=\"%.*s\">\n",
    VERSION,
    i,
    archive_date
  );
  fputs("  xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"", out);
  fputs("\n", out);
  fputs("  xsi:schemaLocation=", out);
  fputs("\"http://noble.gs.washsington.edu/schema/beadstring beadstring.xsd\"\n", out);
  fputs("  xmlns=\"http://zlab.bu.edu/schema/cisml\"\n", out);
  fputs("  xmlns:bstr=\"http://noble.gs.washington.edu/schema/beadstring\"\n>\n", out);
  fprintf(out, "<command-line>%s</command-line>\n", options->command_line);
  print_settings_xml(out, options);
  fprintf(
    out,
    "<alphabet>%s</alphabet>\n",
    options->alphabet == DNA_ALPH ? "nucleotide" : "protein"
  );
  fprintf(
    out,
    "<background source=\"%s\">\n",
    options->bg_filename ? options->bg_filename : "non-redundant database"
  );
  int alph_size = get_alph_size(ALPH_SIZE);
  for (i = 0; i < alph_size; i++) {
    fprintf(
      out,
      "<value letter=\"%c\">%1.3f</value>\n",
      get_alph_char(i),
      get_array_item(i, background)
    );
  }
  fputs("</background>\n", out);
  fputs("<hmm-file>mhmm.xml</hmm-file>\n", out);
  fputs("<hmm-ps>mhmm.ps</hmm-ps>\n", out);
  fputs("<hmm-png>mhmm.png</hmm-png>\n", out);
  fputs("<cisml-file>cisml.xml</cisml-file>\n", out);
  fputs("</beadstring>\n", out);
}

/**********************************************************************
 * This function saves the sequences containing a hit to a motif
 * which is greater than the output p/q-value threshold. If the sequence
 * is less then 10k in length, a copy of the complete sequence will be
 * saved.
 *********************************************************************/
void print_sequences(FILE *out, OPTIONS_T *options, CISML_T *cisml) {

  SEQ_T* sequence = NULL;
  const int max_size_recorded = 10000;

  fputs("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n", out);
  fputs("<matched-sequences>\n", out);

  // Open the FASTA file for reading.
  FILE* fasta_file = NULL;
  if (open_file(options->seq_filename, "r", TRUE, "FASTA", "sequences", &fasta_file) == 0) {
    die("Couldn't open the file %s.\n", options->seq_filename);
  }

  // Read the FASTA file one sequence at a time.
  while (read_one_fasta(fasta_file, MAX_SEQ, &sequence)) {
    char *seq_name = get_seq_name(sequence);
    int seq_length = get_seq_length(sequence);
    // FIXME Only print those sequences that contain a matched element
    if (seq_length <= max_size_recorded) {
      // If sequence is less then limit record raw sequence
      char *raw_sequence = get_raw_sequence(sequence);
      fprintf(
        out,
        "<sequence name=\"%s\" length=\"%d\">\n%s\n",
        seq_name,
        seq_length,
        raw_sequence
      );
      fprintf(out, "</sequence>\n");
    }
    else {
      // Just record name and length
      fprintf(
        out,
        "<sequence name=\"%s\" length=\"%d\"\\>\n",
        seq_name,
        seq_length
      );
    }
  }

  fputs("</matched-sequences>\n", out);
  fclose(fasta_file);
}

/***********************************************************************
  Create diagrams of HMM in Postscript and PNG format using
  'dot'  and  'convert' programs.
 ***********************************************************************/
static void create_hmm_diagram_files(char *mhmm_xml_path, char *mhmm_ps_path,  char *mhmm_png_path) {

  int cmd_buff_size = strlen(mhmm_ps_path) + strlen(mhmm_png_path) + 80;
  char *cmd = mm_malloc(cmd_buff_size * sizeof(char));

  // Create diagram of HMM in Postscript format
  int num_char_printed = snprintf(
    cmd,
    cmd_buff_size,
    "draw-mhmm -consensus %s |dot -Tps > %s",
    mhmm_xml_path,
    mhmm_ps_path
  );
  if (num_char_printed > cmd_buff_size) {
    die("Command buffer too small: %s\n", cmd);
  }
  int result = system(cmd);
  if (result < 0) {
    die("Command %s failed with status %d.\n", cmd, result);
  }

  // Convert Postscript to PNG
  num_char_printed = snprintf(cmd, cmd_buff_size, "convert %s %s", mhmm_ps_path, mhmm_png_path);
  if (num_char_printed > cmd_buff_size) {
    die("Command buffer too small: %s\n", cmd);
  }
  result = system(cmd);
  if (result < 0) {
    die("Command %s failed with status %d.\n", cmd, result);
  }

  myfree(cmd);

}

/***********************************************************************
  Write results to files
 ***********************************************************************/
void print_beadstring_results(
  MHMM_T* the_hmm,
  BOOLEAN_T got_evd,
  OPTIONS_T* options,
  CISML_T *cisml
) {

  // Setup output directory
  if (create_output_directory(
       options->output_dirname,
       options->allow_clobber,
       FALSE /* Don't print warning messages */
      )
    ) {
    // Failed to create output directory.
    die("Couldn't create output directory %s.\n", options->output_dirname);
  }

  // Create model file if one isn't specified.
  char *mhmm_xml_filename = "mhmm.xml";
  char *mhmm_xml_path = NULL;
  if (!options->model_filename) {
    mhmm_xml_path = make_path_to_file(options->output_dirname, mhmm_xml_filename);
    write_mhmm_xml_to_file(the_hmm, mhmm_xml_path);
  }

  // Generate diagram of HMM
  char*  mhmm_ps_path = make_path_to_file(options->output_dirname, "mhmm.ps");
  char*  mhmm_png_path = make_path_to_file(options->output_dirname, "mhmm.png");
  create_hmm_diagram_files(mhmm_xml_path, mhmm_ps_path,  mhmm_png_path);

  // If search peformed output results of search
  if (options->perform_search) {

    // Create the paths to the output files
    const char *HTML_STYLESHEET = "beadstring.xsl";
    const char *CSS_STYLESHEET = "cisml.css";
    const char *GFF_STYLESHEET = "cisml-to-gff.xsl";
    const char *TEXT_STYLESHEET = "cisml-to-text.xsl";
    const char *SEQUENCE_FILENAME = "matched_sequences.xml";
    const char *CISML_FILENAME = "cisml.xml";
    const char *HTML_FILENAME = "beadstring.html";
    const char *TEXT_FILENAME = "beadstring.txt";
    const char *GFF_FILENAME = "beadstring.gff";
    char *sequence_path = make_path_to_file(options->output_dirname, SEQUENCE_FILENAME);
    char *cisml_path = make_path_to_file(options->output_dirname, CISML_FILENAME);
    char *html_stylesheet_path = make_path_to_file(ETC_DIR, HTML_STYLESHEET);
    char *css_stylesheet_path = make_path_to_file(ETC_DIR, CSS_STYLESHEET);
    char *text_stylesheet_path = make_path_to_file(ETC_DIR, TEXT_STYLESHEET);
    char *gff_stylesheet_path = make_path_to_file(ETC_DIR, GFF_STYLESHEET);
    char *html_path = make_path_to_file(options->output_dirname, HTML_FILENAME);
    char *text_path = make_path_to_file(options->output_dirname, TEXT_FILENAME);
    char *gff_path = make_path_to_file(options->output_dirname, GFF_FILENAME);
    char *html_stylesheet_copy_path = make_path_to_file(options->output_dirname, HTML_STYLESHEET);
    char *css_stylesheet_copy_path = make_path_to_file(options->output_dirname, CSS_STYLESHEET);

    // Copy XML to HTML and CSS stylesheets to output directory
    copy_file(html_stylesheet_path, html_stylesheet_copy_path);
    copy_file(css_stylesheet_path, css_stylesheet_copy_path);

    // Open the beadstring XML file.
    const char *BEADSTRING_FILENAME = "beadstring.xml";
    char *beadstring_path = make_path_to_file(options->output_dirname, BEADSTRING_FILENAME);
    FILE *beadstring_file = fopen(beadstring_path, "w");
    if (!beadstring_file) {
      die("Couldn't open file %s for output.\n", beadstring_path);
    }
    print_beadstring_xml_file(
      beadstring_file,
      options,
      the_hmm->background,
      "beadstring.xsl"
    );
    fclose(beadstring_file);

    // Output matched sequence XML
    FILE *sequence_file = fopen(sequence_path, "w");
    if (!sequence_file) {
      die("Couldn't open file %s for output.\n", cisml_path);
    }
    print_sequences(sequence_file, options, cisml);
    fclose(sequence_file);

    // Open the cisml XML file.
    FILE *cisml_file = fopen(cisml_path, "w");
    if (!cisml_file) {
      die("Couldn't open file %s for output.\n", cisml_path);
    }
    print_cisml(cisml_file, cisml, TRUE, NULL, TRUE);
    fclose(cisml_file);

    // Output HTML
    print_xml_filename_to_filename_using_stylesheet(
      beadstring_path, 
      html_stylesheet_copy_path, 
      html_path
    );

    // Output text
    print_xml_filename_to_filename_using_stylesheet(cisml_path, text_stylesheet_path, text_path);

    // Output GFF
    print_xml_filename_to_filename_using_stylesheet(cisml_path, gff_stylesheet_path, gff_path);

    /***********************************************
     * Print header information.
     ***********************************************/
    /*
    write_header(
      "beadstring",
      "Database search results",
      the_hmm->description,
      the_hmm->motif_file,
      NULL, // HMM
      options->seq_filename,
      text_file
    );
    */

    /***********************************************
     * Sort and print the results.
     ***********************************************/
    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(stderr, "\nSorting the E-values.\n");
    }

    /*
    sort_and_print_scores(
      FALSE, // print fancy
      TRUE,  // print_header,
      got_evd,
      options->motif_scoring,
      FALSE, // Don't print in mhmmscan format.
      options->max_seq, // Maximum number of sequences to print
      DEFAULT_OUTPUT_WIDTH,
      options->e_threshold,
      TRUE, // Sort output
      NULL, // FIXME: GFF output?
      text_file
    );
    */

    myfree(html_stylesheet_path);
    myfree(html_stylesheet_copy_path);
    myfree(css_stylesheet_path);
    myfree(css_stylesheet_copy_path);
    myfree(text_stylesheet_path);
    myfree(gff_stylesheet_path);
    myfree(cisml_path);
    myfree(html_path);
    myfree(text_path);
    myfree(gff_path);
    myfree(beadstring_path);
    myfree(mhmm_xml_path);
    myfree(mhmm_ps_path);
    myfree(mhmm_png_path);

  }

}

/***********************************************************************
  Process command line options
 ***********************************************************************/
static void process_command_line(
  int argc,
  char* argv[],
  OPTIONS_T *options
) {

  // Define command line options.
  cmdoption const option_list[] = {
    {"allow-weak-motifs", NO_VALUE},
    {"bgfile", REQUIRED_VALUE},
    {"e-thresh", REQUIRED_VALUE},
    {"fim", NO_VALUE},
    {"gap-extend", REQUIRED_VALUE},
    {"gap-open", REQUIRED_VALUE},
    {"global", NO_VALUE},
    {"maxseqs", REQUIRED_VALUE},
    {"model-file", REQUIRED_VALUE},
    {"motif", REQUIRED_VALUE},
    {"motif-e-thresh", REQUIRED_VALUE},
    {"motif-p-thresh", REQUIRED_VALUE},
    {"motif-pseudo", REQUIRED_VALUE},
    {"no-search", NO_VALUE},
    {"nspacer", REQUIRED_VALUE},
    {"o", REQUIRED_VALUE},
    {"oc", REQUIRED_VALUE},
    {"order", REQUIRED_VALUE},
    {"pam", REQUIRED_VALUE},
    {"paths", REQUIRED_VALUE},
    {"p-score", REQUIRED_VALUE},
    {"progress", REQUIRED_VALUE},
    {"score-file", REQUIRED_VALUE},
    {"spacer-pseudo", REQUIRED_VALUE},
    {"synth", NO_VALUE},
    {"trans-pseudo", REQUIRED_VALUE},
    {"verbosity", REQUIRED_VALUE},
    {"zselo", NO_VALUE},
  };
  int option_count = 28;
  int option_index = 0;

  // Define the usage message.
  char usage[] =
    "USAGE: beadstring [options] <motif file> <sequence file>\n"
     "\n"
     "   Options related to input and output:\n"
     "     --bgfile <background file>\n"
     "     --e-thresh <float>\n"
     "     --maxseqs <int>\n"
     "     --model-file <model file>\n"
     "     --no-search\n"
     "     --oc <output dir> (default=fimo_out)\n"
     "     --output-pthresh <float> (default 1e-4)\n"
     "     --progress <float>\n"
     "     --score-file <score file>\n"
     "     --verbosity 1|2|3|4|5 (default=2)\n"
     "\n"
     "   Options related to selecting motifs for the model:\n"
     "     --motif <string>\n"
     "     --motif-e-thresh <float>\n"
     "     --motif-p-thresh <float>\n"
     "     --order <string>\n"
     "\n"
     "   Options related to building the model:\n"
     "     --fim\n"
     "     --gap-extend <float>\n"
     "     --gap-open <float>\n"
     "     --motif-pseudo <float>\n"
     "     --nspacer <int>\n"
     "     --spacer-pseudo <float>\n"
     "     --trans-pseudo <float>\n"
     "     --zselo\n"
     "\n"
     "   Options related to scoring:\n"
     "     --allow-weak-motifs\n"
     "     --global\n"
     "     --pam <int>\n"
     "     --paths single|all\n"
     "     --p-score <float>\n"
     "\n";

  /* Make sure various options are set to NULL or defaults. */
  options->output_dirname = "beadstring_out";
  options->meme_filename = NULL;
  options->bg_filename = NULL;
  options->score_filename = NULL;
  options->seq_filename = NULL;
  options->model_filename = NULL;
  options->order_string = NULL;
  options->allow_weak_motifs = FALSE;
  options->allow_clobber = TRUE;
  options->fim = FALSE;
  options->local_scoring = TRUE;
  options->motif_scoring = TRUE; // Always do this (not really an option).
  options->synth = FALSE;
  options->zselo = FALSE;
  options->max_seq = -1;
  options->pam_distance = -1;
  options->paths = SINGLE;
  options->spacer_states = DEFAULT_SPACER_STATES,
  options->gap_extend = -1.0;
  options->gap_open = -1.0;
  options->e_threshold = 0.01;
  options->p_threshold = 0.0;
  options->perform_search = TRUE;
  options->progress = 0.0;
  options->requested_motifs = new_string_list();
  options->motif_e_threshold = 0.0;
  options->motif_p_threshold = 0.0;
  options->motif_pseudo = DEFAULT_MOTIF_PSEUDO;
  options->spacer_pseudo = DEFAULT_SPACER_PSEUDO;
  options->trans_pseudo = DEFAULT_TRANS_PSEUDO;
  verbosity = 2;

  simple_setopt(argc, argv, option_count, option_list);

  // Parse the command line.
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

    if (strcmp(option_name, "bgfile") == 0) {
       options->bg_filename = option_value;
    } else if (strcmp(option_name, "e-thresh") == 0) {
      options->e_threshold = atof(option_value);
    } else if (strcmp(option_name, "maxseqs") == 0) {
      options->max_seq = atoi(option_value);
    } else if (strcmp(option_name, "model-file") == 0) {
      options->model_filename = option_value;
    } else if (strcmp(option_name, "no-search") == 0) {
      options->perform_search = FALSE;
    } else if (strcmp(option_name, "o") == 0) {
      options->output_dirname = option_value;
      options->allow_clobber = FALSE;
    } else if (strcmp(option_name, "oc") == 0) {
      options->output_dirname = option_value;
      options->allow_clobber = TRUE;
    } else if (strcmp(option_name, "progress") == 0) {
      options->progress = atof(option_value);
    } else if (strcmp(option_name, "score-file") == 0) {
      options->score_filename = option_value;
    } else if (strcmp(option_name, "verbosity") == 0) {
      verbosity = (VERBOSE_T)atoi(option_value);

    } else if (strcmp(option_name, "motif") == 0) {
      add_string(option_value, options->requested_motifs);
    } else if (strcmp(option_name, "moitf-e-thresh") == 0) {
      options->motif_e_threshold = atof(option_value);
    } else if (strcmp(option_name, "motif-p-thresh") == 0) {
      options->motif_p_threshold = atof(option_value);
    } else if (strcmp(option_name, "order") == 0) {
      options->order_string = option_value;

    } else if (strcmp(option_name, "fim") == 0) {
      options->fim = TRUE;
    } else if (strcmp(option_name, "gap-extend") == 0) {
      options->gap_extend = atof(option_value);
    } else if (strcmp(option_name, "gap-open") == 0) {
      options->gap_open = atof(option_value);
    } else if (strcmp(option_name, "motif-pseudo") == 0) {
      options->motif_pseudo = atof(option_value);
    } else if (strcmp(option_name, "nspacer") == 0) {
      options->spacer_states = atoi(option_value);
    } else if (strcmp(option_name, "spacer-pseudo") == 0) {
      options->spacer_pseudo = atof(option_value);
    } else if (strcmp(option_name, "trans-pseudo") == 0) {
      options->trans_pseudo = atof(option_value);
    } else if (strcmp(option_name, "zselo") == 0) {
      options->zselo = TRUE;

    } else if (strcmp(option_name, "allow-weak-motifs") == 0) {
      options->allow_weak_motifs =  TRUE;
    } else if (strcmp(option_name, "global") == 0) {
      options->local_scoring = FALSE;
    } else if (strcmp(option_name, "motif-scoring") == 0) {
      options->motif_scoring = TRUE;
    } else if (strcmp(option_name, "pam") == 0) {
      options->pam_distance = atoi(option_value);
    } else if (strcmp(option_name, "paths") == 0) {
      if (strcmp(option_value, "single") == 0) {
        options->paths = SINGLE;
      }
      else if (strcmp(option_value, "all") == 0) {
        options->paths = ALL;
      }
      else {
        fprintf(stderr, "Invalid path option: %s.\n", option_value);
        fprintf(stderr, "%s", usage);
        exit(1);
      }
    } else if (strcmp(option_name, "p-score") == 0) {
      options->p_threshold = atof(option_value);
    }
  }

  // Read the two required arguments.
  if (option_index + 2 > argc) {
    fprintf(stderr, "Missing required arguments.\n");
    fprintf(stderr, "%s", usage);
    exit(1);
  }
  options->meme_filename = argv[option_index];
  options->seq_filename = argv[option_index+1];

  //  Record full command line
  int i;
  int arg_length = strlen(argv[0]);
  int total_arg_length = arg_length;
  int buffer_length = 2 * arg_length + 1;
  options->command_line = mm_malloc(buffer_length * sizeof(char));
  strcpy(options->command_line, argv[0]);
  for (i = 1; i < argc; i++) {
    arg_length = strlen(argv[i]) + 2; // +1 for leading space, +1 for trailing null
    total_arg_length += arg_length;
    if (total_arg_length > buffer_length) {
      buffer_length = 2 * total_arg_length + 1;
      options->command_line =
        mm_realloc(options->command_line, buffer_length * sizeof(char));
    }
    strcat(options->command_line, " ");
    strcat(options->command_line, argv[i]);
  }
}

/*************************************************************************
 *  Entry point for program.
 *************************************************************************/
int main(int argc, char *argv[]) {

  OPTIONS_T options;
  MHMM_T *the_hmm = NULL;

  process_command_line(argc, argv, &options);

  if (options.model_filename != NULL ) {

    if (!file_exists(options.model_filename)) {
      die("Couldn't find the model file %s.\n", options.model_filename);
    }

    // Read the HMM.
    read_mhmm(options.model_filename, &the_hmm);

    // Check that the HMM has linear topology
    if (the_hmm->type != LINEAR_HMM) {
      die(
        "beadstring requires a linear topology HMM.\n"
        "The HMM in %s has a %s topology.\n",
        options.model_filename,
        HMM_STRS[the_hmm->type]
      );
    }
  }
  else {
    // If a model file was not specifed, we'll
    // build the HMM from motif occurrences.
    the_hmm = build_linear_model(&options);
  }

  // Create cisml data structure for recording results
  CISML_T *cisml = allocate_cisml(
    "beadstring",
    options.meme_filename,
    options.seq_filename
  );
  PATTERN_T *pattern = allocate_pattern("Beadstring linear HMM","Beadstring linear HMM" );
  add_cisml_pattern(cisml, pattern);

  // Perform the actual search.
  BOOLEAN_T got_evd = FALSE;
  if (options.perform_search) {
    search_with_model(the_hmm, &got_evd, &options, pattern);
  }

  print_beadstring_results(the_hmm, got_evd, &options, cisml);

  // Print the program parameters.
  printf("Program parameters for beadstring\n");
  printf("  MEME file: %s\n", options.meme_filename);
  printf("  Motifs:");
  write_string_list(" ", options.requested_motifs, stdout);
  printf("\n");
  printf("  Model topology: linear\n");
  printf("  States per spacer: %d\n", options.spacer_states);
  printf(
    "  Spacers are free-insertion modules: %s\n",
    options.fim ? "true" : "false"
  );
  printf("\n");

  free_mhmm(the_hmm);
  free_string_list(options.requested_motifs);
  myfree(options.command_line);

  return(0);
}

