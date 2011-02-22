/********************************************************************
 * FILE: pmp_bf.c
 * AUTHOR: Timothy L. Bailey
 * CREATE DATE: 12/12/2007
 * PROJECT: PMP
 * COPYRIGHT: 2007, UQ
 ********************************************************************/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "matrix.h"		// avoid complaint when including array.h
#include "alphabet.h"
#include "evomodel.h"
#include "meme-io.h"
#include "mhmm-state.h"
#include "motif.h"
#include "motiph-scoring.h"
#include "pssm.h"
#include "seq.h"
#include "simple-getopt.h"
#include "tree.h"
#include "object-list.h"

char* program_name = NULL;
VERBOSE_T verbosity = NORMAL_VERBOSE;

// Nucleotide alphabet order as in motif.h
extern char alphabet[];

static void print_col_frequencies(ARRAY_T* alignment_col_freqs);

/*************************************************************************
 * Entry point for pmp_bf
 *************************************************************************/
int main(int argc, char *argv[]) {

  char* bg_filename = NULL;
  char* motif_name = "motif"; // Use this motif name in the output.
  STRING_LIST_T* selected_motifs = NULL;
  double fg_rate = 1.0;
  double bg_rate = 1.0;
  double purine_pyrimidine = 1.0; // r
  double transition_transversion = 0.5; // R
  double pseudocount = 0.1;
  GAP_SUPPORT_T gap_support = SKIP_GAPS;
  MODEL_TYPE_T model_type = F81_MODEL;
  BOOLEAN_T use_halpern_bruno = FALSE;
  char* ustar_label = NULL;	// TLB; create uniform star tree
  int i;

  program_name = "pmp_bf";

  /**********************************************
   * COMMAND LINE PROCESSING
   **********************************************/

  // Define command line options. (FIXME: Repeated code)
  // FIXME: Note that if you add or remove options you
  // must change n_options.
  int n_options = 12;
  cmdoption const pmp_options[] = {
    {"hb", NO_VALUE},
    {"ustar", REQUIRED_VALUE},
    {"model", REQUIRED_VALUE},
    {"pur-pyr", REQUIRED_VALUE},
    {"transition-transversion", REQUIRED_VALUE},
    {"bg", REQUIRED_VALUE},
    {"fg", REQUIRED_VALUE},
    {"motif", REQUIRED_VALUE},
    {"motif-name", REQUIRED_VALUE},
    {"bgfile", REQUIRED_VALUE},
    {"pseudocount", REQUIRED_VALUE},
    {"verbosity", REQUIRED_VALUE}
  };

  int option_index = 0;

  // Define the usage message.
  char      usage[1000] = "";
  strcat(usage, "USAGE: pmp [options] <tree file> <MEME file>\n");
  strcat(usage, "\n");
  strcat(usage, "   Options:\n");

  // Evolutionary model parameters.
  strcat(usage, "     --hb\n");
  strcat(usage, "     --model single|average|jc|k2|f81|f84|hky|tn");
  strcat(usage, " (default=f81)\n");
  strcat(usage, "     --pur-pyr <float> (default=1.0)\n");
  strcat(usage, "     --transition-transversion <float> (default=0.5)\n");
  strcat(usage, "     --bg <float> (default=1.0)\n");
  strcat(usage, "     --fg <float> (default=1.0)\n");

  // Motif parameters.
  strcat(usage, "     --motif <id> (default=all)\n");
  strcat(usage, "     --motif-name <string> (default from motif file)\n");

  // Miscellaneous parameters
  strcat(usage, "     --bgfile <background> (default from motif file)\n");
  strcat(usage, "     --pseudocount <float> (default=0.1)\n");
  strcat(usage, "     --ustar <label>\n");	// TLB; create uniform star tree
  strcat(usage, "     --verbosity [1|2|3|4] (default 2)\n");
  strcat(usage, "\n    Prints the FP and FN rate at each of 10000 score values.\n");
  strcat(usage, "\n    Output format: [<motif_id> score <score> FPR <fpr> TPR <tpr>]+\n");

  // Parse the command line.
  if (simple_setopt(argc, argv, n_options, pmp_options) != NO_ERROR) {
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
    
    if (strcmp(option_name, "model") == 0) {
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
    } else if (strcmp(option_name, "hb") == 0){
        use_halpern_bruno = TRUE;
    } else if (strcmp(option_name, "ustar") == 0){	// TLB; create uniform star tree
        ustar_label = option_value;
    } else if (strcmp(option_name, "pur-pyr") == 0){
        purine_pyrimidine = atof(option_value);
    } else if (strcmp(option_name, "transition-transversion") == 0){
        transition_transversion = atof(option_value);
    } else if (strcmp(option_name, "bg") == 0){
      bg_rate = atof(option_value);
    } else if (strcmp(option_name, "fg") == 0){
      fg_rate = atof(option_value);
    } else if (strcmp(option_name, "motif") == 0){
        if (selected_motifs == NULL) {
          selected_motifs = new_string_list();
        }
       add_string(option_value, selected_motifs);
    } else if (strcmp(option_name, "motif-name") == 0){
        motif_name = option_value;
    } else if (strcmp(option_name, "bgfile") == 0){
      bg_filename = option_value;
    } else if (strcmp(option_name, "pseudocount") == 0){
        pseudocount = atof(option_value);
    } else if (strcmp(option_name, "verbosity") == 0){
        verbosity = atoi(option_value);
    }
  }

  // Must have tree and motif file names
  if (argc != option_index + 2) {
    fprintf(stderr, "%s", usage);
    exit(EXIT_FAILURE);
  } 

  /**********************************************
   * Read the phylogenetic tree.
   **********************************************/
  char* tree_filename = NULL;
  TREE_T* tree = NULL;
  tree_filename = argv[option_index];
  option_index++;
  tree = read_tree_from_file(tree_filename);

  // get the species names
  STRING_LIST_T* alignment_species = make_leaf_list(tree);
  char *root_label = get_label(tree);	// in case target in center
  if (strlen(root_label)>0) add_string(root_label, alignment_species);
  //write_string_list(" ", alignment_species, stderr);

  // TLB; Convert the tree to a uniform star tree with
  // the target sequence at its center.
  if (ustar_label != NULL) {
    tree = convert_to_uniform_star_tree(tree, ustar_label);
    if (tree == NULL) 
      die("Tree or alignment missing target %s\n", ustar_label);
    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(stderr, 
	"Target %s placed at center of uniform (d=%.3f) star tree:\n", 
          ustar_label, get_total_length(tree) / get_num_children(tree) 
      );
      write_tree(tree, stderr);
    }
  }

  /**********************************************
   * Read the motifs.
   **********************************************/
  char* meme_filename = argv[option_index];
  option_index++;
  int num_motifs = 0; 
  MOTIF_T motifs[2 * MAX_MOTIFS];
  STRING_LIST_T* motif_occurrences = NULL;
  BOOLEAN_T has_reverse_strand = FALSE;
  ARRAY_T* bg_freqs = NULL;

  read_meme_file(
		 meme_filename, 
		 bg_filename,
		 pseudocount,
     REQUIRE_PSPM,
		 &num_motifs, 
		 motifs, 
		 &motif_occurrences,
		 &has_reverse_strand, 
		 &bg_freqs
		 );

  // We don't make use of the motif occurence data
  free_string_list(motif_occurrences);

  // TLB; need to resize bg_freqs array to ALPH_SIZE items
  // or copy array breaks in HB mode.  This throws away
  // the freqs for the ambiguous characters;
  int alph_size = get_alph_size(ALPH_SIZE);
  resize_array(bg_freqs, alph_size);

  /**************************************************************
  * Compute probability distributions for each of the selected motifs.
  **************************************************************/
  int motif_index;
  for (motif_index = 0; motif_index < num_motifs; motif_index++) {

    MOTIF_T* motif = &(motifs[motif_index]);
    char* motif_id = get_motif_id(motif);
    char* bare_motif_id = motif_id;

    // We may have specified on the command line that
    // only certain motifs were to be used.
    if (selected_motifs != NULL) {
      if (*bare_motif_id == '+' || *bare_motif_id == '-') {
        // The selected  motif id won't included a strand indicator.
        bare_motif_id++;
      }
      if (have_string(bare_motif_id, selected_motifs) == FALSE) {
        continue;
      }
    }

    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(
        stderr, 
        "Using motif %s of width %d.\n",
        motif_id, motif->length
      );
    }

    // Build an array of evolutionary models for each position in the motif.
    EVOMODEL_T** models = make_motif_models(
      motif, 
      bg_freqs,
      model_type,
      fg_rate, 
      bg_rate, 
      purine_pyrimidine, 
      transition_transversion, 
      use_halpern_bruno
    );

    // Get the frequencies under the background model (row 0) 
    // and position-dependent scores (rows 1..w)
    // for each possible alignment column.
    MATRIX_T* pssm_matrix = build_alignment_pssm_matrix(
      alignment_species,
      motif->length + 1, 
      models, 
      tree, 
      gap_support
    );
    ARRAY_T* alignment_col_freqs = allocate_array(get_num_cols(pssm_matrix)); 
    copy_array(get_matrix_row(0, pssm_matrix), alignment_col_freqs);
    remove_matrix_row(0, pssm_matrix);		// throw away first row
    //print_col_frequencies(alignment_col_freqs);

    //
    // Get the position-dependent null model alignment column frequencies
    //
    int w = motif->length;
    int ncols = get_num_cols(pssm_matrix); 
    MATRIX_T* pos_dep_bkg = allocate_matrix(w, ncols);
    for (i=0; i<w; i++) {
      // get the evo model corresponding to this column of the motif
      // and store it as the first evolutionary model.
      myfree(models[0]);
      // Use motif PSFM for equilibrium freqs. for model.
      ARRAY_T* site_specific_freqs = allocate_array(alph_size);
      int j = 0;
      for(j = 0; j < alph_size; j++) {
	double value = get_matrix_cell(i, j, motif->freqs);
	set_array_item(j, value, site_specific_freqs);
      }
      if (use_halpern_bruno == FALSE) {
	models[0] = make_model(
	  model_type,
	  fg_rate,
	  transition_transversion,
	  purine_pyrimidine,
	  site_specific_freqs,
          NULL
	);
      } else {
        models[0] = make_model(
	  model_type,
	  fg_rate,
	  transition_transversion,
	  purine_pyrimidine,
	  bg_freqs,
	  site_specific_freqs
	);
      }
      // get the alignment column frequencies using this model
      MATRIX_T* tmp_pssm_matrix = build_alignment_pssm_matrix(
	alignment_species,
	2,				// only interested in freqs under bkg
	models, 
	tree, 
	gap_support
      );
      // assemble the position-dependent background alignment column freqs.
      set_matrix_row(i, get_matrix_row(0, tmp_pssm_matrix), pos_dep_bkg);
      // chuck the pssm (not his real name)
      free_matrix(tmp_pssm_matrix);
    }

    //
    // Compute and print the score distribution under the background model
    // and under the (position-dependent) motif model.
    //
    int range = 10000;	// 10^4 gives same result as 10^5, but 10^3 differs

    // under background model
    PSSM_T* pssm = build_matrix_pssm(pssm_matrix, alignment_col_freqs, range);

    // under position-dependent background (motif) model
    PSSM_T* pssm_pos_dep = build_matrix_pssm(pssm_matrix, alignment_col_freqs, range);
    get_pv_lookup_pos_dep(
      pssm_pos_dep, 
      pos_dep_bkg, 
      NULL // no priors used
    );

    // print FP and FN distributions
    int num_items = get_pssm_pv_length(pssm_pos_dep);
    for (i=0; i<num_items; i++) {
      double pvf = get_pssm_pv(i, pssm);
      double pvt = get_pssm_pv(i, pssm_pos_dep);
      double fpr = pvf;
      double fnr = 1 - pvt;
      if (fpr >= 0.99999 || fnr == 0) continue;
      printf("%s score %d FPR %.3g FNR %.3g\n", motif_id, i, fpr, fnr);
    }

    // free stuff
    free_pssm(pssm);
    free_pssm(pssm_pos_dep);
    if (models != NULL) {
      int model_index;
      int num_models = motif->length + 1;
      for (model_index = 0; model_index < num_models; model_index++) {
        free_model(models[model_index]);
      }
      myfree(models);
    }

  } // motif

  for (motif_index = 0; motif_index < num_motifs; motif_index++) {
    free_matrix(motifs[motif_index].freqs);
  }

  /**********************************************
   * Clean up.
   **********************************************/
  // TLB may have encountered a memory corruption bug here
  // CEG has not been able to reproduce it. valgrind says all is well.
  free_array(bg_freqs);
  free_tree(TRUE, tree);
  free_string_list(selected_motifs);

  return(0);
} // main

/*******************************************************************
  Print the column frequency distribution.
 ********************************************************************/
static void print_col_frequencies(
  ARRAY_T* alignment_column_freqs
)
{
  int i;
  int num_freqs = get_array_length(alignment_column_freqs);
  int alph_size = get_alph_size(ALPH_SIZE);
  int num_leaves = NINT(log(num_freqs)/log(alph_size));
  char* alignment_col = mm_malloc((num_leaves + 1) * sizeof(char));
  for (i=0; i<num_freqs; i++) {
    unhash_alignment_col(
      i,                              //col_index
      alignment_col,
      num_leaves
    );
    printf("%s %d %g\n", alignment_col, i+1, 
      get_array_item(i, alignment_column_freqs));
  }
} // print_col_freqs
