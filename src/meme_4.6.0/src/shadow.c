/********************************************************************
 * FILE: shadow.c
 * AUTHOR: William Stafford Noble, Charles E. Grant
 * CREATE DATE: 12/17/2004
 * PROJECT: EVOMCAST
 * COPYRIGHT: 2004, UW
 ********************************************************************/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "alignment.h"
#include "alphabet.h"
#include "cisml.h"
#include "clustalw-io.h"
#include "evomodel.h"
#include "fasta-io.h"
#include "meme-io.h"
#include "mhmm-state.h"
#include "motif.h"
#include "motiph-scoring.h"
#include "pssm.h"
#include "seq.h"
#include "simple-getopt.h"
#include "tree.h"

// Nucleotide alphabest order as in motif.h
extern char alphabet[];

char* program_name = NULL;
VERBOSE_T verbosity = NORMAL_VERBOSE;

static ARRAY_T* read_bg_from_file(char* bg_filename);
static ARRAY_T* read_bg_from_alignment(
  ALIGNMENT_T* alignment,
  MODEL_TYPE_T model_type
);

/*************************************************************************
 * Entry point for shadow
 *************************************************************************/
int main(int argc, char *argv[]) {

  BOOLEAN_T compute_pvalues = TRUE;
  BOOLEAN_T compute_qvalues = TRUE;
  BOOLEAN_T use_file_list = FALSE;
  BOOLEAN_T allow_clobber = TRUE;
  BOOLEAN_T text_only = FALSE;
  BOOLEAN_T output_pthresh_set = FALSE;
  BOOLEAN_T output_qthresh_set = FALSE;
  int ref_seq_index = 0;  // Which sequence provides reference coords.
  //FIXME: Eventually, above should be accessible via a command-line switch.
  char *output_dirname = "shadow_out";
  char* bg_filename = NULL;
  int max_stored_scores = 100000;
  int selected_motif_index = 0;
  double output_pthresh = 1e-4;
  double output_qthresh = 1.0;
  double fg_rate = 1.0;
  double bg_rate = 1.0;
  double gap_cost = 0.0;
  double purine_pyrimidine = 1.0; // r
  double transition_transversion = 0.5; // R
  GAP_SUPPORT_T gap_support = SKIP_GAPS;
  MODEL_TYPE_T model_type = F81_MODEL;
  program_name = "shadow";

  /**********************************************
   * COMMAND LINE PROCESSING
   **********************************************/

  // Define command line options.
  const int num_options = 18;
  cmdoption const shadow_options[] = {
    {"bg", REQUIRED_VALUE},
    {"bgfile", REQUIRED_VALUE},
    {"fg", REQUIRED_VALUE},
    {"gap", REQUIRED_VALUE},
    {"gap-cost", REQUIRED_VALUE},
    {"list", NO_VALUE},
    {"max-stored-scores", REQUIRED_VALUE},
    {"model", REQUIRED_VALUE},
    {"no-pvalue", NO_VALUE},
    {"no-qvalue", NO_VALUE},
    {"o", REQUIRED_VALUE},
    {"oc", REQUIRED_VALUE},
    {"output-pthresh", REQUIRED_VALUE},
    {"output-qthresh", REQUIRED_VALUE},
    {"pur-pyr", REQUIRED_VALUE},
    {"text", NO_VALUE},
    {"transition-transversion", REQUIRED_VALUE},
    {"verbosity", REQUIRED_VALUE}
  };

  int option_index = 0;

  // Define the usage message.
  char usage[1000] = "";
  strcat(usage, "USAGE: shadow [options] <alignment> <tree>\n");
  strcat(usage, "\n");
  strcat(usage, "   Options:\n");

  // Evolutionary model parameters.
  strcat(usage, "     --bg <float> (default=1.0)\n");
  strcat(usage, "     --fg <float> (default=1.0)\n");
  strcat(usage, "     --gap skip | fixed | wildcard | minimum");
  strcat(usage, " (default=skip)\n");
  strcat(usage, "     --gap-cost <float> (default=0.0)\n");
  strcat(usage, "     --list\n");
  strcat(usage, "     --model [single|average|jc|k2|f81|f84|hky|tn]");
  strcat(usage, " (default=f81)\n");
  strcat(usage, "     --pur-pyr <float> (default=1.0)\n");
  strcat(usage, "     --transition-transversion <float> (default=0.5)\n");

  // Miscellaneous parameters
  strcat(usage, "     --bgfile <background> (default from alignment)\n");
  strcat(usage, "     --max-stored-scores (Default 100,000)\n");
  strcat(usage, "     --no-pvalue (default false)\n");
  strcat(usage, "     --no-qvalue (default false)\n");
  strcat(usage, "     --o <output dir> (default=shadow_out)\n");
  strcat(usage, "     --oc <output dir> (default=shadow_out)\n");
  strcat(usage, "     --output-pthresh <float> (default 1e-4)\n");
  strcat(usage, "     --text\n");
  strcat(usage, "     --verbosity [1|2|3|4] (default 2)]\n");
  strcat(usage, "\n");


  // Parse the command line.
  if (simple_setopt(argc, argv, num_options, shadow_options) != NO_ERROR) {
    die("Error processing command line options: option name too long.\n");
  }

  while (TRUE) {
    int c = 0;
    char* option_name = NULL;
    char* option_value = NULL;
    const char * message = NULL;

    // Read the next option, and break if we're done.
    c = simple_getopt(&option_name, &option_value, &option_index);
    if (c == 0) {
      break;
    } else if (c < 0) {
      (void) simple_getopterror(&message);
      die("Error processing command line options (%s)\n", message);
    }

    if (strcmp(option_name, "bg") == 0){
      bg_rate = atof(option_value);
    }
    else if (strcmp(option_name, "bgfile") == 0){
      bg_filename = option_value;
    }
    else if (strcmp(option_name, "fg") == 0){
      fg_rate = atof(option_value);
    }
    else if (strcmp(option_name, "gap") == 0){
      if (strcmp(option_value, "skip") == 0) {
        gap_support = SKIP_GAPS;
      }
      else if (strcmp(option_value, "fixed") == 0) {
        gap_support = FIXED_GAP_COST;
      }
      else if (strcmp(option_value, "wildcard") == 0) {
        gap_support = WILDCARD_GAP;
      }
      else if (strcmp(option_value, "minimum") == 0) {
        gap_support = MIN_GAPS;
      }
      else {
        die("Unknown gap handling method: %s\n", option_value);
      }
    }
    else if (strcmp(option_name, "gap-cost") == 0) {
      gap_cost = atof(option_value);
    }
    else if (strcmp(option_name, "list") == 0){
      use_file_list = TRUE;
    }
    else if (strcmp(option_name, "model") == 0) {
      if (strcmp(option_value, "jc") == 0) {
        model_type = JC_MODEL;
      } else if (strcmp(option_value, "k2") == 0) {
        model_type = K2_MODEL;
      } else if (strcmp(option_value, "f81") == 0) {
        model_type = F81_MODEL;
      } else if (strcmp(option_value, "f84") == 0) {
        model_type = F84_MODEL;
      } else if (strcmp(option_value, "hky") == 0) {
        model_type = HKY_MODEL;
      } else if (strcmp(option_value, "tn") == 0) {
        model_type = TAMURA_NEI_MODEL;
      } else if (strcmp(option_value, "single") == 0) {
        model_type = SINGLE_MODEL;
      } else if (strcmp(option_value, "average") == 0) {
        model_type = AVERAGE_MODEL;
      } else {
        die("Unknown model: %s\n", option_value);
      }
    }
    else if (strcmp(option_name, "no-pvalue") == 0){
      compute_pvalues = FALSE;
      compute_qvalues = FALSE;
    }
    else if (strcmp(option_name, "no-qvalue") == 0){
      compute_qvalues = FALSE;
    }
    else if (strcmp(option_name, "o") == 0){
      // Set output directory with no clobber
      output_dirname = option_value;
      allow_clobber = FALSE;
    }
    else if (strcmp(option_name, "oc") == 0){
      // Set output directory with clobber
      output_dirname = option_value;
      allow_clobber = TRUE;
    }
    else if (strcmp(option_name, "output-pthresh") == 0){
      output_pthresh = atof(option_value);
      output_qthresh = 1.0;
      output_pthresh_set = TRUE;
    }
    else if (strcmp(option_name, "output-qthresh") == 0){
      output_qthresh = atof(option_value);
      output_pthresh = 1.0;
      output_qthresh_set = TRUE;
    }
    else if (strcmp(option_name, "pur-pyr") == 0){
        purine_pyrimidine = atof(option_value);
    }
    else if (strcmp(option_name, "text") == 0){
      text_only = TRUE;
    }
    else if (strcmp(option_name, "transition-transervsion") == 0){
        transition_transversion = atof(option_value);
    }
    else if (strcmp(option_name, "verbosity") == 0){
        verbosity = atoi(option_value);
    }
  }

  // Check that pvalue and qvalue options are consistent
  if (compute_pvalues == FALSE
    && (output_pthresh_set == TRUE || output_qthresh_set == TRUE)) {
    die(
      "The --no-pvalue option cannot be used with the\n"
      " --output-pthresh or the --output-qthresh options"
    );
  }
  if (compute_qvalues == FALSE && output_qthresh_set == TRUE) {
    die("The --no-qvalue option cannot be used with the --output-qthresh options");
  }

  // Must have tree and alignment file names
  if (argc != option_index + 2) {
    fprintf(stderr, "%s", usage);
    exit(EXIT_FAILURE);
  }

  set_alphabet(verbosity, "ACGT");

  /****************************************************
   * Read the names of the alignment files.
   ****************************************************/
  STRING_LIST_T* filenames = NULL;
  int num_filenames = 0;
  if (use_file_list) {
    filenames = read_string_list_from_file(argv[option_index]);
    num_filenames = get_num_strings(filenames);
  }
  else {
    filenames = new_string_list();
    add_string(argv[option_index], filenames);
    num_filenames = 1;
  }
  option_index++;

  /**********************************************
   * Read the phylogenetic tree.
   **********************************************/
  char* tree_filename = NULL;
  TREE_T* tree = NULL;
  if (MOTIPH || SHADOW) {
    tree_filename = argv[option_index];
    option_index++;
    tree = read_tree_from_file(tree_filename);
  }

  // Verify that the tree doesn't have duplicate species.
  STRING_LIST_T* tree_ids = make_leaf_list(tree);
  if (has_duplicates("Duplicate IDs in tree:", tree_ids)) {
    exit(1);
  }
  free_string_list(tree_ids);

  /**********************************************
   * Determine the background frequencies if
   * background file provided.
   **********************************************/
  ARRAY_T* bg_freqs = NULL;
  if (bg_filename != NULL) {
    bg_freqs = read_bg_from_file(bg_filename);
    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(stderr, "Read background frequencies from %s.\n", bg_filename);
    }
  }

  // Create cisml data structure for recording results
  CISML_T *cisml = allocate_cisml("shadow", "none", "clustal-w alignment");
  set_cisml_site_pvalue_cutoff(cisml, output_pthresh);
  set_cisml_site_qvalue_cutoff(cisml, output_qthresh);

  // Create cisml pattern and add to cisml record
  PATTERN_T *pattern = allocate_pattern("conservation", "conservation");
  set_pattern_max_stored_matches(pattern, max_stored_scores);
  add_cisml_pattern(cisml, pattern);


  /****************************************************
   * Do the scoring.
   ****************************************************/
  int window_size = 1;
  const int num_models = 2;
  EVOMODEL_T* models[num_models];

  // Score the reference sequence in the alignment using these models.
  int file_index;
  for(file_index = 0; file_index < num_filenames; file_index++) {

    int current_ref_seq_index = ref_seq_index;

    // Get the next alignment.
    char* filename = get_nth_string(file_index, filenames);
    ALIGNMENT_T* alignment = read_alignment_from_file(
      filename,
      TRUE, // Sort by species name
      FALSE, // Remove gaps
      &current_ref_seq_index // update index to reference sequence
    );
    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(stderr, "Read alignment from %s.\n", filename);
    }

    // Make sure there are no duplicates.
    STRING_LIST_T* alignment_species = get_species_names(alignment);
    if (has_duplicates("Duplicate IDs in alignment:", alignment_species)) {
      exit(1);
    }

    // Trim the tree, eliminating species not in this alignment.
    TREE_T* trimmed_tree = trim_tree(TRUE, tree, alignment_species);

    // Check that at least one species was found.
    if (trimmed_tree == NULL) {
       die("Your tree doesn't contain any of the species in your alignment.");
    }

    if (verbosity >= HIGH_VERBOSE) {
      fprintf(stderr, "Trimmed tree: ");
      write_tree(trimmed_tree, stderr);
    }
    STRING_LIST_T* tree_species = make_leaf_list(trimmed_tree);
    sort_string_list(tree_species); 	// keep species alphabetical

    // Trim the alignment, eliminating species not in this tree.
    ALIGNMENT_T* trimmed_alignment = remove_alignment_seqs(tree_species, alignment);
    free_alignment(alignment);
    alignment = NULL;
    free_string_list(tree_species);
    alignment_species = get_species_names(trimmed_alignment);
    sort_string_list(alignment_species);	// sort species alphabetically

    // Take background freq. from alignment if no background file specified.
    if (bg_freqs == NULL) {
      bg_freqs = read_bg_from_alignment(trimmed_alignment, model_type);
      if (verbosity >= NORMAL_VERBOSE) {
        fprintf(stderr, "Read background frequencies from alignment.\n");
      }
    }

    // Build a table for converting the index into the alignment
    // to an index into ungapped reference sequence.
    int* coord_conv_table = make_alignment_to_seq_table(current_ref_seq_index, trimmed_alignment);

    // Build the background model using the equilibrium frequencies
    // from the alignment and the background mutation rate.
    models[0] = make_model(
			model_type,
			bg_rate,
			transition_transversion,
			purine_pyrimidine,
			bg_freqs,
			FALSE
    );

    // Build the foreground model using the equilibrium frequencies
    // from the alignment and the foreground mutation rate.
    models[1] = make_model(
			model_type,
			fg_rate,
			transition_transversion,
			purine_pyrimidine,
			bg_freqs,
			FALSE
    );
    if (verbosity >= HIGH_VERBOSE) {
      fprintf(stderr, "Created foreground and background model.\n");
    }

    // Create a scanned_sequence record and record it in pattern.
    SCANNED_SEQUENCE_T *scanned_seq =
      allocate_scanned_sequence(filename, filename, pattern);
    set_scanned_sequence_length(scanned_seq, get_alignment_length(trimmed_alignment));

    // If we need p-values, use a lot more memory.
    if (compute_pvalues) {
      // Build PSSM for motiph
      MATRIX_T* pssm_matrix = build_alignment_pssm_matrix(
        alignment_species,
        window_size + 1,
        models,
        trimmed_tree,
        gap_support
      );
      if (verbosity >= HIGH_VERBOSE) {
        fprintf(stderr, "Created PSSM.\n");
      }

      // The first row of the matrix contains the  probabilities for
      // the alignment columns given the background model
      ARRAY_T* alignment_col_probs = allocate_array(get_num_cols(pssm_matrix));
      copy_array(get_matrix_row(0, pssm_matrix), alignment_col_probs);
      remove_matrix_row(0, pssm_matrix);

      // Build tables to translate log-odds scores to p-values
      PSSM_T* pssm = build_matrix_pssm(
	pssm_matrix,
	alignment_col_probs,
	PSSM_RANGE 
      );
      if (verbosity >= HIGH_VERBOSE) {
        fprintf(stderr, "Created p-value lookup table.\n");
      }

      score_sequence_in_alignment(
        current_ref_seq_index,
        trimmed_alignment,
        "shadow", // motif_id
        NULL, // tree
        window_size,
        models,
        NULL, // background frequencies
        pssm,
        coord_conv_table,
        gap_support,
        gap_cost,
        output_pthresh,
        scanned_seq
      );
      free_array(alignment_col_probs);
      free_matrix(pssm_matrix);
    }

    // No pvalues; don't precompute and store scores.  Much slower.
    else {
      score_sequence_in_alignment(
				ref_seq_index,
				trimmed_alignment,
				"shadow", // motif ID (FIXME?)
				trimmed_tree,
				1, // Window size
				models, // 0 = background, 1 = foreground
				bg_freqs,
				NULL, // pssm
				coord_conv_table,
				gap_support,
				gap_cost,
				0, // output_pthresh
				scanned_seq
      );
    }

    free_string_list(alignment_species);
    free_tree(TRUE, trimmed_tree);
    free_alignment(trimmed_alignment);
    if (models != NULL) {
      int model_index;
      for (model_index = 0; model_index < (window_size + 1); model_index++) {
        free_model(models[model_index]);
      }
    }
  }

  // All the positions have been scanned
  set_pattern_is_complete(pattern);

  // Compute q-values, if requested.
  if (compute_pvalues && compute_qvalues) {
    pattern_calculate_qvalues(pattern, NULL);
  }

  // Write out results
  if (text_only) {
    print_cisml_as_text(cisml);
  }
  else {
    print_full_results(
      cisml,
      output_dirname,
      "shadow.xml",
      "shadow.html",
      "shadow.txt",
      "shadow.gff",
      allow_clobber,
      TRUE
    );
  }

  free_array(bg_freqs);
  free_tree(TRUE, tree);
  free_string_list(filenames);
  free_cisml(cisml);

  return(0);
}

/*******************************************************************
  Read the background nucleotide frequenices from a file
 ********************************************************************/
ARRAY_T* read_bg_from_file(char* bg_filename) {

  FILE* bg_file = NULL;
  if (open_file(bg_filename, "r", TRUE, "background",
                "background nucleotide frequencies", &bg_file) == 0) {
    die("Couldn't open the file %s.\n", bg_filename);
  }

  // Read four nucleotide frequencies from file
  const int bufsize = 80;
  char line[bufsize];
  char word[bufsize];
  int num_read = 0;
  float a_freq = 0.0;
  float c_freq = 0.0;
  float t_freq = 0.0;
  float g_freq = 0.0;
  while (TRUE) {
    if (fgets(line, 80, bg_file) == NULL) {
      break;
    }
    if (sscanf(line, "# %s", word) == 1) {
      continue;
    }
    if (sscanf(line, "A %f", &a_freq) == 1) {
      num_read++;
      continue;
    }
    if (sscanf(line, "C %f", &c_freq)) {
      num_read++;
      continue;
    }
    if (sscanf(line, "G %f", &g_freq) == 1) {
      num_read++;
      continue;
    }
    if (sscanf(line, "T %f", &t_freq) == 1) {
      num_read++;
      continue;
    }
    // Once we've read the four frequences we don't
    // need to read any more.
    if (num_read == 4) {
      break;
    }
  }

  fclose(bg_file);

  if (num_read != 4) {
    die("Can't read background frequences from file %s.\n", bg_filename);
  }

  int alph_size = get_alph_size(ALPH_SIZE);
  ARRAY_T* bg_freqs = allocate_array(alph_size);
  int index = 0;
  index = alphabet_index('A', alphabet);
  set_array_item(index, a_freq, bg_freqs);
  index = alphabet_index('C', alphabet);
  set_array_item(index, c_freq, bg_freqs);
  index = alphabet_index('G', alphabet);
  set_array_item(index, g_freq, bg_freqs);
  index = alphabet_index('T', alphabet);
  set_array_item(index, t_freq, bg_freqs);

  return bg_freqs;
}

/*******************************************************************
  Determine the background nucleotide frequenices from an alignment.
 ********************************************************************/
static ARRAY_T* read_bg_from_alignment(
  ALIGNMENT_T* alignment,
  MODEL_TYPE_T model_type
)
{
  // Set up the frequency array.
  ARRAY_T* bg_freqs;
  if (model_type == SINGLE_MODEL) {
    SEQ_T* seq = get_alignment_sequence(0, alignment);
    bg_freqs = get_sequence_freqs(seq);
    free_seq(seq);
    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(stderr, "Computing background frequencies from %s.\n",
              get_seq_name(seq));
    }
  } else {
    bg_freqs = get_alignment_freqs(alignment);
    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(stderr, "Computing background frequencies from alignment.\n");
    }
  }
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Background frequencies: ");
    print_array(bg_freqs, 6, 3, TRUE, stderr);
  }

  return bg_freqs;
}

