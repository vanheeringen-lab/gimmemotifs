/********************************************************************
 * FILE: fimo.c
 * AUTHOR: William Stafford Noble, Charles E. Grant, Timothy Bailey
 * CREATE DATE: 12/17/2004
 * PROJECT: MEME suite
 * COPYRIGHT: 2004-2007, UW
 ********************************************************************/

#define DEFINE_GLOBALS
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "matrix.h"
#include "alphabet.h"
#include "cisml.h"
#include "dir.h"
#include "fasta-io.h"
#include "hash_alph.h"
#include "io.h"
#include "meme-io.h"
#include "projrel.h"
#include "pssm.h"
#include "simple-getopt.h"
#include "string-list.h"
#include "utils.h"
#include "xml-util.h"

char* program_name = "fimo";
VERBOSE_T verbosity = NORMAL_VERBOSE;

// Structure for tracking fimo command line parameters.
typedef struct options {

  BOOLEAN_T allow_clobber;      // Allow overwritting of files in output directory.
  BOOLEAN_T compute_qvalues;    // Compute q-values
  BOOLEAN_T output_pthresh_set; // p-value threshold has been set from command line.
  BOOLEAN_T output_qthresh_set; // q-value threshold has been set from command line.
  BOOLEAN_T text_only;          // Generate only plain text output
  BOOLEAN_T scan_both_strands;  // Scan forward and reverse strands

  char* bg_filename;                     // Name of file file containg background freq.
  char* command_line;                    // Full command line
  char* meme_filename;                   // Name of file containg motifs.
  char* motif_name;                      // Use this motif name in the output.
  char* output_dirname;                  // Name of the output directory
  char* seq_filename;                    // Name of file containg sequences.
  char* seq_name;                        // Use this sequence name in the output.

  int max_seq_length; // Maximum allowed sequence length.
  int max_stored_scores; // Maximum number of matches to store per pattern.

  double pseudocount; // Pseudocount added to Motif PSFM.
  double pthresh;     // Maximum p-value to report.
  double qthresh;     // Maximum q-value to report.

  ALPH_T alphabet;    // Alphabet specified by MEME file.
  STRING_LIST_T* selected_motifs; // Indices of requested motifs.

  int num_seqs;     // Number of sequences in FASTA file
  int num_residues; // Number of residues in FASTA file

  char* pval_lookup_filename;   // Print p-value lookup table.
  char *html_stylesheet_path; // Path to master copy of HTML stlesheet for CisML
  char *html_stylesheet_local_path; // Path to working copy of HTML stylesheet for CisML
  char *css_stylesheet_path; // Path to master copy of CSS style sheet for CisML
  char *css_stylesheet_local_path; // Path to working copy of CSS stylesheet for CisML
  char *text_stylesheet_path; // Path to plain-text stylesheet for CisML
  char *gff_stylesheet_path; // Path to GFF stylesheeet for CisML
  char *html_path; // Path to FIMO HTML output file
  char *text_path; // Path to FIMO plain-text output file
  char *gff_path; // Pathe to FIMO GFF output file
  char *fimo_path; // Pathe to FIMO XML output file
  char *cisml_path; // Path to CisML XML output file

  const char* HTML_STYLESHEET;  // Name of HTML XSLT stylesheet.
  const char* CSS_STYLESHEET;   // Name of CSS stylesheet.
  const char* GFF_STYLESHEET;   // Name of GFF XSLT stylesheet.
  const char* TEXT_STYLESHEET;  // Name of plain-text XSLT stylesheet.
  const char* HTML_FILENAME;    // Name of HTML output file.
  const char* TEXT_FILENAME;    // Name of plain-text output file.
  const char* GFF_FILENAME;     // Name of GFF output file.
  const char* FIMO_FILENAME;    // Name of FIMO XML output file.
  const char* CISML_FILENAME;   // Name of CisML XML output file.

  const char* usage; // Usage statment

} FIMO_OPTIONS_T;

/***********************************************************************
  Free memory allocated in options processing
 ***********************************************************************/
static void cleanup_options(FIMO_OPTIONS_T *options) {
  myfree(options->command_line);
  myfree(options->text_stylesheet_path);
  myfree(options->gff_stylesheet_path);
  myfree(options->html_path);
  myfree(options->text_path);
  myfree(options->gff_path);
  myfree(options->html_stylesheet_path);
  myfree(options->html_stylesheet_local_path);
  myfree(options->css_stylesheet_path);
  myfree(options->css_stylesheet_local_path);
  myfree(options->fimo_path);
  myfree(options->cisml_path);
}

/***********************************************************************
  Process command line options
 ***********************************************************************/
static void process_command_line(
  int argc,
  char* argv[],
  FIMO_OPTIONS_T *options
) {

  // Define command line options.
  const int num_options = 14;
  cmdoption const fimo_options[] = {
    {"bgfile", REQUIRED_VALUE},
    {"max-seq-length", REQUIRED_VALUE},
    {"max-stored-scores", REQUIRED_VALUE},
    {"motif", REQUIRED_VALUE},
    {"motif-pseudo", REQUIRED_VALUE},
    {"norc", NO_VALUE},
    {"o", REQUIRED_VALUE},
    {"oc", REQUIRED_VALUE},
    {"output-pthresh", REQUIRED_VALUE},
    {"output-qthresh", REQUIRED_VALUE},
    {"no-qvalue", NO_VALUE},
    {"text", NO_VALUE},
    {"verbosity", REQUIRED_VALUE},
    {"pval-lookup", REQUIRED_VALUE} // This option is hidden from users.
  };

  // Define the usage message.
  options->usage =
    "USAGE: fimo [options] <motif file> <sequence file>\n"
    "\n"
    "   Options:\n"
    "     --bgfile <background> (default from NR sequence database)\n"
    "     --max-seq-length <int> (default=2.5e8)\n"
    "     --max-stored-scores <int> (default=100000)\n"
    "     --motif <id> (default=all)\n"
    "     --motif-pseudo <float> (default=0.1)\n"
    "     --norc\n"
    "     --o <output dir> (default=fimo_out)\n"
    "     --oc <output dir> (default=fimo_out)\n"
    "     --output-pthresh <float> (default 1e-4)\n"
    "     --output-qthresh <float> (default 1.0)\n"
    "     --no-qvalue\n"
    "     --text\n"
    "     --verbosity [1|2|3|4] (default 2)\n"
    "\n"
    "   Use \'-\' for <sequence file> to read the database from standard input.\n"
    "   Use \'--bgfile motif-file\' to read the background from the motif file.\n"
    "\n";

  int option_index = 0;

  /* Make sure various options are set to NULL or defaults. */
  options->allow_clobber = TRUE;
  options->compute_qvalues = TRUE;
  options->output_pthresh_set = FALSE;
  options->output_qthresh_set = FALSE;
  options->text_only = FALSE;
  options->scan_both_strands = TRUE;

  options->bg_filename = NULL;
  options->command_line = NULL;
  options->meme_filename = NULL;
  options->motif_name = "motif";
  options->seq_name = NULL;
  options->output_dirname = "fimo_out";
  options->seq_filename = NULL;

  options->max_seq_length = MAX_SEQ;
  options->max_stored_scores = 100000;

  options->pseudocount = 0.1;
  options->pthresh = 1e-4;
  options->qthresh = 1.0;

  options->num_seqs = 0;
  options->num_residues = 0;

  options->selected_motifs = new_string_list();
  options->pval_lookup_filename = NULL;
  verbosity = 2;

  simple_setopt(argc, argv, num_options, fimo_options);

  // Parse the command line.
  while (TRUE) {
    int c = 0;
    char* option_name = NULL;
    char* option_value = NULL;
    const char * message = NULL;

    // Read the next option, and break if we're done.
    c = simple_getopt(&option_name, &option_value, &option_index);
    if (c == 0) {
      break;
    }
    else if (c < 0) {
      (void) simple_getopterror(&message);
      die("Error processing command line options (%s)\n", message);
    }
    if (strcmp(option_name, "bgfile") == 0){
      options->bg_filename = option_value;
    }
    else if (strcmp(option_name, "max-seq-length") == 0) {
      // Use atof and cast to be able to read things like 1e8.
      options->max_seq_length = (int)atof(option_value);
    }
    else if (strcmp(option_name, "max-stored-scores") == 0) {
      // Use atof and cast to be able to read things like 1e8.
      options->max_stored_scores = (int)atof(option_value);
    }
    else if (strcmp(option_name, "motif") == 0){
      if (options->selected_motifs == NULL) {
        options->selected_motifs = new_string_list();
      }
      add_string(option_value, options->selected_motifs);
    }
    else if (strcmp(option_name, "motif-pseudo") == 0){
      options->pseudocount = atof(option_value);
    }
    else if (strcmp(option_name, "norc") == 0){
      options->scan_both_strands = FALSE;
    }
    else if (strcmp(option_name, "o") == 0){
      // Set output directory with no clobber
      options->output_dirname = option_value;
      options->allow_clobber = FALSE;
    }
    else if (strcmp(option_name, "oc") == 0){
      // Set output directory with clobber
      options->output_dirname = option_value;
      options->allow_clobber = TRUE;
    }
    else if (strcmp(option_name, "output-pthresh") == 0){
      options->pthresh = atof(option_value);
      options->qthresh = 1.0;
      options->output_pthresh_set = TRUE;
    }
    else if (strcmp(option_name, "output-qthresh") == 0){
      options->qthresh = atof(option_value);
      options->pthresh = 1.0;
      options->output_qthresh_set = TRUE;
    }
    else if (strcmp(option_name, "no-qvalue") == 0){
      options->compute_qvalues = FALSE;
    }
    else if (strcmp(option_name, "text") == 0){
      options->text_only = TRUE;
    }
    else if (strcmp(option_name, "verbosity") == 0){
      verbosity = atoi(option_value);
    }
    else if (strcmp(option_name, "pval-lookup") == 0) {
      options->pval_lookup_filename = option_value;
    }
  }

  // Check that qvalue options are consistent
  if (options->compute_qvalues == FALSE && options->output_qthresh_set == TRUE) {
    die("The --no-qvalue option cannot be used with the --output-qthresh options");
  }

  // Turn off q-values if text only.
  if (options->text_only == TRUE) {
    if (options->compute_qvalues) {
      fprintf(stderr, "Warning: text mode turns off computation of q-values\n");
    }
    options->compute_qvalues = FALSE;
  }

  // Must have sequence and motif file names
  if (argc != option_index + 2) {
    fprintf(stderr, "%s", options->usage);
    exit(EXIT_FAILURE);
  }
  // Record the command line
  options->command_line = get_command_line(argc, argv);

  // Record the input file names
  options->meme_filename = argv[option_index];
  option_index++;
  options->seq_filename = argv[option_index];
  option_index++;

  // Set up path values for needed stylesheets and output files.
  options->HTML_STYLESHEET = "fimo.xsl";
  options->CSS_STYLESHEET = "cisml.css";
  options->GFF_STYLESHEET = "cisml-to-gff.xsl";
  options->TEXT_STYLESHEET = "cisml-to-text.xsl";
  options->HTML_FILENAME = "fimo.html";
  options->TEXT_FILENAME = "fimo.txt";
  options->GFF_FILENAME = "fimo.gff";
  options->FIMO_FILENAME = "fimo.xml";
  options->CISML_FILENAME = "cisml.xml";
  options->text_stylesheet_path = make_path_to_file(ETC_DIR, options->TEXT_STYLESHEET);
  options->gff_stylesheet_path = make_path_to_file(ETC_DIR, options->GFF_STYLESHEET);
  options->html_path = make_path_to_file(options->output_dirname, options->HTML_FILENAME);
  options->text_path = make_path_to_file(options->output_dirname, options->TEXT_FILENAME);
  options->gff_path = make_path_to_file(options->output_dirname, options->GFF_FILENAME);
  options->html_stylesheet_path = make_path_to_file(ETC_DIR, options->HTML_STYLESHEET);
  options->html_stylesheet_local_path
    = make_path_to_file(options->output_dirname, options->HTML_STYLESHEET);
  options->css_stylesheet_path = make_path_to_file(ETC_DIR, options->CSS_STYLESHEET);
  options->css_stylesheet_local_path
    = make_path_to_file(options->output_dirname, options->CSS_STYLESHEET);
  options->fimo_path = make_path_to_file(options->output_dirname, options->FIMO_FILENAME);
  options->cisml_path = make_path_to_file(options->output_dirname, options->CISML_FILENAME);

}

/*************************************************************************
 * Calculate the log odds score for each motif-sized window at each
 * site in the sequence using the given nucleotide frequencies.
 *************************************************************************/
static void fimo_score_sequence(
  BOOLEAN_T       text_only,
  SEQ_T*          seq,
  MOTIF_T*        motif,
  ARRAY_T*        bg_freqs,
  char*           feature_type,
  char*           seq_name,
  double          pthresh,
  PSSM_T*         pssm,
  SCANNED_SEQUENCE_T* scanned_seq
)
{
  assert(seq != NULL);
  assert(motif != NULL);
  assert(bg_freqs != NULL);

  ARRAY_T* pv_lookup = pssm->pv;
  MATRIX_T* pssm_matrix = pssm->matrix;

  // Get data for sequence
  if (seq_name == NULL) {
    seq_name = get_seq_name(seq);
  }
  char* raw_seq = get_raw_sequence(seq);
  int seq_length = get_seq_length(seq);

  // Prepare storage for the string representing the portion
  // of the reference sequence within the window.
  char* window_seq = (char *) mm_malloc(sizeof(char) * (motif->length + 1));
  window_seq[motif->length] = '\0';

  // Get a character for the strand.
  char strand = 0;
  if (motif->id!= NULL) {
    strand = *motif->id; /* FIXME: This will not work when we generalize
                            the motif naming convention. */
  }
  else {
    strand = '.';
  }

  int max_index = seq_length - motif->length + 1;
  char* alphabet = get_alphabet(FALSE);
  int alph_size = get_alph_size(ALPH_SIZE);

  // For each site in the sequence
  int seq_index;
  for (seq_index = 0; seq_index < max_index; seq_index++) {

    BOOLEAN_T process_window = TRUE;
    double scaled_log_odds = 0.0;

    if ((verbosity >= HIGH_VERBOSE) && (((seq_index + 1) % 10000000) == 0)) {
      fprintf(stderr, "Scoring position %d.\n", seq_index + 1);
    }

    // For each site in the motif window
    int motif_position;
    for (motif_position = 0; motif_position < motif->length; motif_position++) {

      char c = raw_seq[seq_index + motif_position];
      int alph_index = alphabet_index(c, alphabet);
      window_seq[motif_position] = c;

      // Check for gaps and ambiguity codes at this site
      if(c == '-' || c == '.' || alph_index >= alph_size) {
          process_window = FALSE;
          break;
      }

      scaled_log_odds += get_matrix_cell(motif_position, alph_index, pssm_matrix);

    }

    if (process_window == TRUE) {
      if ((int)scaled_log_odds >= get_array_length(pv_lookup)) {
        scaled_log_odds = (float)(get_array_length(pv_lookup) - 1);
      }
      double pvalue = get_array_item((int) scaled_log_odds, pv_lookup);

      // Skip matches with pvalue less then the threshold.
      int start = 0;
      int stop = 0;
      if (strand == '-') {
        // Reverse the sense of start/stop for the reverse strand
        start = seq_index + motif->length;
        stop = seq_index + 1;
      }
      else {
        start = seq_index + 1;
        stop = seq_index + motif->length;
      }
      if (text_only) {
        if (pvalue <= pthresh) {
          fprintf(stdout, "%s", get_motif_id(motif));
          fprintf(stdout, "\t%s", get_seq_name(seq));
          fprintf(stdout, "\t%d", start);
          fprintf(stdout, "\t%d", stop);
          fprintf(
            stdout,
            "\t%g",
            get_unscaled_pssm_score(scaled_log_odds, pssm)
          );
          fprintf(stdout, "\t%g", pvalue);
          fprintf(stdout, "\t%s\n", window_seq);
        }
      } else {
        double score = get_unscaled_pssm_score(scaled_log_odds, pssm);
        // Create matched element and add too pattern.
        MATCHED_ELEMENT_T *element = allocate_matched_element_with_score(
          start,
          stop,
          score,
          pvalue,
          scanned_seq
        );
        add_scanned_sequence_scanned_element(scanned_seq);
        PATTERN_T* pattern = get_scanned_sequence_parent(scanned_seq);
        BOOLEAN_T added = add_pattern_matched_element(pattern, element);
        if (added == TRUE) {
          if (get_pattern_has_all_pvalues(pattern)) {
            // Raw sequence is left in scanned sequence
          }
          else {
            // Sub-set of raw sequence copied to matched element
            set_matched_element_sequence(element, window_seq);
          }
        }
        else {
          // Element was rejected by pattern
          free_matched_element(element);
        }
      }

    }
  }

  myfree(window_seq);
}

/***********************************************************************
 * Print FIMO settings information to an XML file
 ***********************************************************************/
static  void print_settings_xml(
  FILE *out,
  FIMO_OPTIONS_T *options
) {

  fputs("<settings>\n",  out);
  fprintf(out, "<setting name=\"%s\">%s</setting>\n", "output directory", options->output_dirname);
  fprintf(out, "<setting name=\"%s\">%s</setting>\n", "MEME file name", options->meme_filename);
  fprintf(out, "<setting name=\"%s\">%s</setting>\n", "sequence file name", options->seq_filename);
  fprintf(out, "<setting name=\"%s\">%s</setting>\n", "background file name", options->bg_filename);
  if (options->motif_name) {
    fprintf(out, "<setting name=\"%s\">%s</setting>\n", "motif name", options->motif_name);
  }
  if (options->seq_name) {
    fprintf(out, "<setting name=\"%s\">%s</setting>\n", "sequence name", options->seq_name);
  }
  fprintf(
    out,
    "<setting name=\"%s\">%s</setting>\n",
    "allow clobber",
    boolean_to_string(options->allow_clobber)
  );
  fprintf(
    out,
    "<setting name=\"%s\">%s</setting>\n",
    "compute q-values",
    boolean_to_string(options->compute_qvalues)
  );
  fprintf(
    out,
    "<setting name=\"%s\">%s</setting>\n",
    "output p-threshold set",
    boolean_to_string(options->output_pthresh_set)
  );
  fprintf(
    out,
    "<setting name=\"%s\">%s</setting>\n",
    "output q-threshold set",
    boolean_to_string(options->output_qthresh_set)
  );
  fprintf(
    out,
    "<setting name=\"%s\">%s</setting>\n",
    "text only",
    boolean_to_string(options->text_only)
  );
  fprintf(
    out,
    "<setting name=\"%s\">%s</setting>\n",
    "scan both strands",
    boolean_to_string(options->scan_both_strands)
  );
  fprintf(
    out,
    "<setting name=\"%s\">%d</setting>\n",
    "max sequence length",
    options->max_seq_length
  );
  fprintf(
    out,
    "<setting name=\"%s\">%3.2g</setting>\n",
    "output q-value threshold",
    options->qthresh);
  fprintf(
    out,
    "<setting name=\"%s\">%3.2g</setting>\n",
    "output p-value threshold",
    options->pthresh
  );
  fprintf(
    out,
    "<setting name=\"%s\">%3.2g</setting>\n",
    "pseudocount",
    options->pseudocount
  );
  fprintf(
    out,
    "<setting name=\"%s\">%d</setting>\n",
    "verbosity",
    verbosity
  );
  int i = 0;
  int num_strings = get_num_strings(options->selected_motifs);
  for(i = 0; i < num_strings; i++) {
    fprintf(
      out,
      "<setting name=\"%s\">%s</setting>\n", "selected motif",
      get_nth_string(i, options->selected_motifs)
    );
  }

  fputs("</settings>\n",  out);

}
/***********************************************************************
 * Print FIMO specific information to an XML file
 ***********************************************************************/
static void print_fimo_xml_file(
  FILE *out,
  FIMO_OPTIONS_T *options,
  MOTIF_T *motifs,
  int num_motifs,
  ARRAY_T *bgfreq,
  char  *stylesheet
) {

  fputs("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n", out);
  if (stylesheet != NULL) {
    fprintf(out, "<?xml-stylesheet type=\"text/xsl\" href=\"%s\"?>\n", stylesheet);
  }
  fputs("<!-- Begin document body -->\n", out);
  const int head_skip = 18; // Skip SVN keyword
  const int tail_skip = 2; // Skip close of SVN keyword
  const char *archive_date = ARCHIVE_DATE + head_skip;
  int i = strlen(archive_date) - tail_skip;
  fprintf(out,
    "<fimo version=\"%s\" release=\"%.*s\">\n",
    VERSION,
    i,
    archive_date
  );
  fputs("  xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"", out);
  fputs("\n", out);
  fputs("  xsi:schemaLocation=", out);
  fputs("  xmlns:fimo=\"http://noble.gs.washington.edu/schema/fimo\"\n>\n", out);
  fprintf(out, "<command-line>%s</command-line>\n", options->command_line);
  print_settings_xml(out, options);
  fprintf(
    out,
    "<sequence-data num-sequences=\"%d\" num-residues=\"%d\"/>\n",
    options->num_seqs,
    options->num_residues
  );
  fprintf(
    out,
    "<alphabet>%s</alphabet>\n",
    options->alphabet == DNA_ALPH ? "nucleotide" : "protein"
  );
  for (i = 0; i < num_motifs; i++) {
    char *bare_motif_id = motifs[i].id;
    // Drop the strand indicator from the motif id
    if (*bare_motif_id == '+' || *bare_motif_id == '-') {
      bare_motif_id++;
    }
    char *best_possible_match = get_best_possible_match(&(motifs[i]));
    fprintf(
      out,
      "<motif name=\"\%s\" width=\"%d\" best-possible-match=\"%s\"/>\n",
      bare_motif_id,
      motifs[i].length,
      best_possible_match
    );
    myfree(best_possible_match);
  }
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
      get_array_item(i, bgfreq)
    );
  }
  fputs("</background>\n", out);
  fputs("<cisml-file>cisml.xml</cisml-file>\n", out);
  fputs("</fimo>\n", out);
}

/**********************************************************************
 * This function saves the sequences containing a hit to a motif
 * which is greater then the output p/q-value threshold. If the sequence
 * is less then 10k in length, a copy of the complete sequence will be
 * saved.
 *********************************************************************/
void print_matched_sequences(FILE *out, FIMO_OPTIONS_T *options, CISML_T *cisml) {

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
  while (read_one_fasta(fasta_file, options->max_seq_length, &sequence)) {
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

/**********************************************************************
 * This function saves CisML results as a set of files in a
 * directory. The file names are provided by the input parameters:
 *   cisml is a pointer to the ciml structure to be printed.
 *   html_filename will be the name of the HTML output
 *   text_filename will be the name of the plain text output
 *   gff_filename will be the name of the GFF output
 *   allow_clobber will determine whether or not existing files will
 *                 be overwritten.
 *********************************************************************/
void create_output_files(
  CISML_T *cisml,
  FIMO_OPTIONS_T *options,
  MOTIF_T *motifs,
  int num_motifs,
  ARRAY_T *bgfreqs
) {

  // Copy XML to HTML and CSS stylesheets to output directory
  copy_file(options->html_stylesheet_path, options->html_stylesheet_local_path);
  copy_file(options->css_stylesheet_path, options->css_stylesheet_local_path);


  FILE *fimo_file = fopen(options->fimo_path, "w");
  if (!fimo_file) {
    die("Couldn't open file %s for output.\n", options->fimo_path);
  }
  print_fimo_xml_file(
    fimo_file,
    options,
    motifs,
    num_motifs,
    bgfreqs,
    "fimo.xsl"
  );
  fclose(fimo_file);

  // Output HTML
  print_xml_filename_to_filename_using_stylesheet(
    options->fimo_path,
    options->html_stylesheet_local_path,
    options->html_path
  );

  // Output text
  print_xml_filename_to_filename_using_stylesheet(
    options->cisml_path,
    options->text_stylesheet_path,
    options->text_path
  );

  // Output GFF
  print_xml_filename_to_filename_using_stylesheet(
    options->cisml_path,
    options->gff_stylesheet_path,
    options->gff_path
  );
}

/*************************************************************************
 * Entry point for fimo
 *************************************************************************/
int main(int argc, char *argv[]) {

  FIMO_OPTIONS_T options;

  /**********************************************
   * COMMAND LINE PROCESSING
   **********************************************/
  process_command_line(argc, argv, &options);

  // Create cisml data structure for recording results
  CISML_T *cisml = NULL;
  if (!options.text_only) {
    cisml = allocate_cisml("fimo", options.meme_filename, options.seq_filename);
    set_cisml_site_pvalue_cutoff(cisml, options.pthresh);
    set_cisml_site_qvalue_cutoff(cisml, options.qthresh);
  }

  // Open the output file for the p-value lookup table, if requested.
  FILE* pval_lookup_file = NULL;
  if (options.pval_lookup_filename != NULL) {
    if (open_file(options.pval_lookup_filename, "w", FALSE,
      "p-value lookup table",
      "p-value lookup table", &pval_lookup_file) == 0) {
      exit(EXIT_FAILURE);
    }
  }

  /**********************************************
   * Read the motifs.
   **********************************************/
  int num_motifs = 0;
  int num_motifs_after_rc = 0;
  MOTIF_T motifs[2 * MAX_MOTIFS];
  BOOLEAN_T has_reverse_strand;
  ARRAY_T* bg_freqs = NULL;
  FILE *cisml_file = NULL;

  read_meme_file(
    options.meme_filename,
    options.bg_filename,
    options.pseudocount,
    REQUIRE_PSPM,  //need PSPMs
    &num_motifs,
    motifs,
    NULL, //motif occurrences, not used
    &has_reverse_strand,
    &bg_freqs
  );

  // Open the FASTA file for reading.
  FILE* fasta_file = NULL;
  if (open_file(options.seq_filename, "r", TRUE, "FASTA", "sequences",
    &fasta_file) == 0) {
    die("Couldn't open the file %s.\n", options.seq_filename);
  }

  // We don't make use of the motif occurence data

  // If motifs use protein alphabet we will not scan both strands
  ALPH_T alphabet_type = which_alphabet();
  if (alphabet_type == PROTEIN_ALPH) {
    options.alphabet = PROTEIN_ALPH;
    options.scan_both_strands = FALSE;
  }
  else {
    options.alphabet = DNA_ALPH;
  }

  if (options.scan_both_strands == TRUE) {
    // Set up hash tables for computing reverse complement
    setup_hash_alph(DNAB);
    setalph(0);
    // Correct background by averaging on freq. for both strands.
    average_freq_with_complement(bg_freqs);
    int alph_size = get_alph_size(ALPH_SIZE);
    normalize_subarray(0, alph_size, 0.0, bg_freqs);
    fill_in_ambiguous_chars(FALSE, bg_freqs);
    // Make reverse complement motifs.
    add_reverse_complements(&num_motifs, motifs);
  }

  // Print the text-only header line.
  if (options.text_only) {
    fprintf(stdout, "Motif");
    fprintf(stdout, "\tSeq");
    fprintf(stdout, "\tStart");
    fprintf(stdout, "\tStop");
    fprintf(stdout, "\tLog-odds");
    fprintf(stdout, "\tp-value");
    fprintf(stdout, "\tSite\n");
  }
  else {
    // Create output directory
    if (create_output_directory(
         options.output_dirname,
         options.allow_clobber,
         FALSE /* Don't print warning messages */
        )
      ) {
      // Failed to create output directory.
      die("Couldn't create output directory %s.\n", options.output_dirname);
    }
    // Open the CisML XML file for output
    cisml_file = fopen(options.cisml_path, "w");
    if (!cisml_file) {
      die("Couldn't open file %s for output.\n", options.cisml_path);
    }

    // Print the opening section of the CisML XML
    print_cisml_start(
      cisml_file,
      cisml,
      TRUE, // print header
      options.HTML_STYLESHEET, // Name of HTML XSLT stylesheet
      TRUE // print namespace in XML header
    );

    // Print the parameters section of the CisML XML
    print_cisml_parameters(cisml_file, cisml);
  }

  /**************************************************************
   * Score each of the sites for each of the selected motifs.
   **************************************************************/
  int motif_index;
  for (motif_index = 0; motif_index < num_motifs; motif_index++) {

    MOTIF_T* motif = &(motifs[motif_index]);
    char* motif_id = get_motif_id(motif);
    char* bare_motif_id = motif_id;

    // Drop the strand indicator from the motif id
    if (*bare_motif_id == '+' || *bare_motif_id == '-') {
      bare_motif_id++;
    }

    if ((get_num_strings(options.selected_motifs) == 0)
      || (have_string(bare_motif_id, options.selected_motifs) == TRUE)) {

      if (verbosity >= NORMAL_VERBOSE) {
        fprintf(
          stderr,
          "Using motif %s of width %d.\n",
          motif_id,
          motif->length
        );
      }

      // Build PSSM for motif and tables for p-value calculation.
      // FIXME: the non-averaged freqs should be used for p-values
      PSSM_T* pos_pssm = build_motif_pssm(motif, bg_freqs, bg_freqs, PSSM_RANGE, 0, FALSE);

      if (options.pval_lookup_filename != NULL) {
        // Print the cumulative lookup table.
        print_array(pos_pssm->pv, 8, 6, TRUE, pval_lookup_file);
        // Also print the non-cumulative version.
        int i;
        for (i = 1; i < get_array_length(pos_pssm->pv); i++) {
          fprintf(pval_lookup_file, "%g ",
          get_array_item(i-1, pos_pssm->pv) - get_array_item(i, pos_pssm->pv));
        }
        fprintf(pval_lookup_file, "\n");
      }

      // If required, do the same for the reverse complement motif.
      MOTIF_T* rev_motif = NULL;
      PSSM_T* rev_pssm = NULL;
      if (options.scan_both_strands) {
        rev_motif = motif + 1;
        ++motif_index;
        motif_id = get_motif_id(rev_motif);
        // FIXME: the non-averaged freqs should be used for p-values
        rev_pssm = build_motif_pssm(rev_motif, bg_freqs, bg_freqs, PSSM_RANGE, 0, FALSE);
        if (verbosity >= NORMAL_VERBOSE) {
          fprintf(
            stderr,
            "Using motif %s of width %d.\n",
            motif_id,
            motif->length
          );
        }
      }

      // Create cisml pattern
      PATTERN_T* pattern = NULL;
      if (!options.text_only) {
        pattern = allocate_pattern(bare_motif_id, bare_motif_id);
        set_pattern_max_stored_matches(pattern, options.max_stored_scores);
      }

      options.num_seqs = 0;
      options.num_residues = 0;
      SEQ_T* sequence = NULL;

      // Read the FASTA file one sequence at a time.
      while (read_one_fasta(fasta_file, options.max_seq_length, &sequence)) {

        int seq_length = get_seq_length(sequence);

        ++options.num_seqs;
        options.num_residues += seq_length;
        char * seq_name = NULL;
        if (options.seq_name) {
         seq_name = options.seq_name;
        }
        else {
          seq_name = get_seq_name(sequence);
        }

        // Create a scanned_sequence record and record it in pattern.
        SCANNED_SEQUENCE_T* scanned_seq = NULL;
        if (!options.text_only) {
          scanned_seq =
            allocate_scanned_sequence(seq_name, seq_name, pattern);
          if (get_pattern_has_all_pvalues(pattern)) {
            // If we are holding all matches, store sequence in scanned_sequence
            set_scanned_sequence_raw_seq(scanned_seq, get_raw_sequence(sequence));
          }
          set_scanned_sequence_length(scanned_seq, seq_length);
        }

        if (verbosity >= DUMP_VERBOSE) {
          fprintf(stderr, "Scoring forward strand.\n");
        }

        // Score the forward strand.
        fimo_score_sequence(
          options.text_only,
          sequence,
          motif,
          bg_freqs,
          options.motif_name,
          seq_name,
          options.pthresh,
          pos_pssm,
          scanned_seq
        );

        // Score the reverse strand.
        if (options.scan_both_strands) {
          if (verbosity >= DUMP_VERBOSE) {
            fprintf(stderr, "Scoring reverse strand.\n");
          }
          fimo_score_sequence(
            options.text_only,
            sequence,
            rev_motif,
            bg_freqs,
            options.motif_name,
            seq_name,
            options.pthresh,
            rev_pssm,
            scanned_seq
          );
        }

        free_seq(sequence);

      }  // All sequences parsed

      // The pattern is complete.
      if (!options.text_only) {
        set_pattern_is_complete(pattern);
      }

      // Compute q-values, if requested.
      if (options.compute_qvalues) {
        pattern_calculate_qvalues(pattern);
      }
      if (!options.text_only) {
        print_cisml_start_pattern(cisml, cisml_file, pattern);
        int num_scanned_sequences = get_pattern_num_scanned_sequences(pattern);
        SCANNED_SEQUENCE_T **scanned_sequences
          = get_pattern_scanned_sequences(pattern);
        print_cisml_scanned_sequences(
          cisml,
          cisml_file,
          num_scanned_sequences,
          scanned_sequences
        );
        print_cisml_end_pattern(cisml_file);
        free_pattern(pattern);
      }


      rewind(fasta_file);

      // Free memory associated with this motif.
      free_pssm(pos_pssm);
      free_pssm(rev_pssm);

    }
    else {
      if (verbosity >= NORMAL_VERBOSE) {
        fprintf(stderr, "Skipping motif %s.\n", motif_id);
      }
    }


  } // All motifs parsed

  if (options.text_only) {
    // Close text file
  }
  else {
    // print cisml end
    print_cisml_end(cisml_file);
    fclose(fasta_file);
    fclose(cisml_file);
    // Output other file formats
    create_output_files(cisml, &options, motifs, num_motifs, bg_freqs);
    free_cisml(cisml);
  }

  // Clean up.
  // The motifs are automatically allocated, but we 
  // do need to free the frequency matrix for each one.
  for (motif_index = 0; motif_index < num_motifs; motif_index++) {
    free_matrix(motifs[motif_index].freqs);
  }
  free_array(bg_freqs);
  free_string_list(options.selected_motifs);
  cleanup_options(&options);

  return 0;

}

