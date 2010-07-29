/****************************************************************************
 * FILE: motiph_scoring.c
 * AUTHOR: William Stafford Noble, Charles E. Grant, Timothy L. Bailey
 * CREATE DATE: 12/03/2004
 * PROJECT: EVOMCAST
 * DESCRIPTION: Define scores using the ratio of log-likelihoods
 * between foreground an background models.
 * COPYRIGHT: 2004, UW
 ****************************************************************************/
#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "motiph-scoring.h"
#include "alphabet.h"
#include "motif.h"
#include "mhmm-state.h"
#include "pssm.h"
#include "utils.h"


extern char* program_name;

/*************************************************************************
 *  Build a Position Specific Scoring Matrix (PSSM) for alignment columns
 *  given a phylogenetic tree and an array of evolutionary models. 
 *  The first model in the array is assumed to be the background model. 
 *  The PSSM will have one row for each model, and one column for each 
 *  possible alignment column.
 *
 *  The elements of the PSSM are the log-odds scores for the corresponding
 *  model and alignment column.
 *
 *  The caller is responsible for freeing the returned PSSM.
 *************************************************************************/
MATRIX_T* build_alignment_pssm_matrix(
  STRING_LIST_T* seq_names, 
  int num_models,
  EVOMODEL_T** models, 
  TREE_T* tree, 
  GAP_SUPPORT_T gap_support
) 
{
  assert(seq_names != NULL);
  assert(num_models >= 2);
  assert(models != NULL);
  assert(tree != NULL);

  SUBSTMATRIX_TABLE_T* substmatrix_table = NULL;
  SUBSTMATRIX_ARRAY_T* bg_substmatrices = NULL;
  SUBSTMATRIX_ARRAY_T* fg_substmatrices = NULL;

  // The number of possible alignment columns is (alph_size^column_height)
  // The column height is given by the number of labels in the tree.
  int alph_size = get_alph_size(ALPH_SIZE);
  // Use the first column, not the background column (0) because
  // AVERAGE_MODEL has JC_MODEL in background column.
  MODEL_TYPE_T model_type = get_model_type(models[1]);
  //int num_leaves;
  int num_labels;		// TLB; Count labels, not leaves to get size
				// of multiple alignment.
  if (model_type == SINGLE_MODEL) {
    num_labels = 1;
  }
  else {
    num_labels = label_count(tree);
    // Create a table of the transition probability matrices
    // for each model and each edge of the tree.
    // One model is needed for the background to compute column frequencies
    // for the AVERAGE_MODEL.
    substmatrix_table = make_substmatrix_table(
      tree, 
      (model_type == AVERAGE_MODEL) ? 1 : num_models, 
      models
    );
    // Get the array of substitution matrices for the
    // background model
    bg_substmatrices = get_substmatrix_array_from_table(substmatrix_table, 0);
  }

  int num_alignment_cols = (int) pow((double) alph_size, (double) num_labels);
  char* alignment_col = mm_malloc((num_labels + 1) * sizeof(char));
  MATRIX_T* pssm_matrix = allocate_matrix(num_models, num_alignment_cols);
  EVOMODEL_T* bg_model = models[0];

  // Process in column major order to avoid re-calculating bg_likelihood
  int col_index;
  double total_bg_likelihood = 0.0;
  double total_fg_likelihood = 0.0;
  for (col_index = 0; col_index < num_alignment_cols; col_index++) {
    // Translate the column number into a string giving the
    // bases for that column.
    unhash_alignment_col(
      col_index, 
      alignment_col, 
      num_labels
    );

    // Calculate the likelihood for this column using the background model.
    double bg_likelihood = site_likelihood(
      alignment_col, 
      seq_names, 
      tree,
      (model_type == AVERAGE_MODEL) ? JC_MODEL : model_type,
      get_model_equil_freqs(bg_model),
      bg_substmatrices,
      gap_support
    );

    // save background likelihood of this column for use as column frequency
    int row_index = 0;
    set_matrix_cell(row_index, col_index, bg_likelihood, pssm_matrix);

    // recompute likelihood if using AVERAGE_MODEL for use in the
    // denominator of the log-odds
    if (model_type == AVERAGE_MODEL) {
      bg_likelihood = site_likelihood(
	alignment_col, 
	seq_names, 
	tree,
	model_type,
	get_model_equil_freqs(bg_model),
	bg_substmatrices,
	gap_support
      );
    }

    total_bg_likelihood += bg_likelihood;

    for (row_index = 1; row_index < num_models; row_index++) {
      if (model_type != SINGLE_MODEL && model_type != AVERAGE_MODEL) {
        // Get the array of substitution matrices for the
        // current foreground model.
        fg_substmatrices =
          get_substmatrix_array_from_table(substmatrix_table, row_index);
      }
      // Calculate the log-odds for this column using the
      // foreground and background models.
      double fg_likelihood = site_likelihood(
        alignment_col, 
        seq_names, 
        tree, 
        model_type,
        get_model_equil_freqs(models[row_index]),
        fg_substmatrices,
        gap_support
      );
      total_fg_likelihood += fg_likelihood;

      double log_odds = my_log2(fg_likelihood / bg_likelihood);
      set_matrix_cell(row_index, col_index, log_odds, pssm_matrix);
      free_substmatrix_array(fg_substmatrices);
    }
  }

  free_substmatrix_array(bg_substmatrices);
  free_substmatrix_table(substmatrix_table);
  myfree(alignment_col);

  return pssm_matrix;
}

/****************************************************************************
 *  This function implements the Felsentein pruning algorithm for evaluating 
 *  the log likelihood of a phylogentic tree for a given site in an alignment. 
 *  It sums over all possible base assignments of the internal nodes 
 *  using the law of total probability and recursively calls itself until it 
 *  reaches the tree leaves. At the leaves paths through the tree whose final
 *  base doesn't match the appropriate sequence site are "pruned" 
 *  by returning a value of 0.0.
 *  TLB: Compute
 *    Pr(t | r(t)=root_base), where "t" is the tree and "r(t)" is the base
 *  assigned to the root of t.
 ****************************************************************************/
static double pruning_algorithm(
  char* alignment_col,
  STRING_LIST_T* seq_names,
  int root_base,			// base at "root" of current subtree
  TREE_T* tree,
  SUBSTMATRIX_ARRAY_T* substmatrix_array,
  GAP_SUPPORT_T gap_support
) 
{
  assert(alignment_col != NULL);
  assert(seq_names != NULL);
  assert(tree != NULL);
  assert(substmatrix_array != NULL);

  char* alphabet = get_alphabet(FALSE);
  int alph_size = get_alph_size(ALPH_SIZE);
  assert(root_base < alph_size);

  int new_root_base = 0;
  double sum = 0.0; // Sum for normal case.
  double min = 0.0; // Minimum for min_gap strategy.
  double result = 0.0;
  if (is_leaf(tree) == TRUE) {
    // End of recursion, get character at this leaf
    new_root_base = get_base_from_node(
      tree, 
      seq_names, 
      alignment_col, 
      alphabet, 
      alph_size
    );

    if (root_base == new_root_base) {
      // TLB; Pr(T | r(t)=observed_base)
      result = 1.0;
    }
    else if (
          gap_support == WILDCARD_GAP 
          && (new_root_base == '-' || new_root_base == '.')
        ) {
        result = 1.0;
    } 
    else {
        result = 0.0;
    }
  } 
  else {
    result = 1.0;
    int num_children = get_num_children(tree);
    int c = 0;
    // TLB; Pr(t | r(t) = root_base) = \prod_{children(t)} Pr(child | r(t) = root_base)
    for (c = 0; c < num_children; c++) {
      min = HUGE_VAL;
      sum = 0.0;
      TREE_T* child = get_nth_child(c, tree);
      MATRIX_T* prob_matrix = 
        get_substmatrix_for_time(substmatrix_array, get_length(child));
      assert(prob_matrix != NULL);

      // TLB
      // Pr(child | r(t) = root_base) = 
      //   \sum_{x \in \alphabet} Pr(r(t)->x | time=d) Pr(child | r(child)=x)
      for (new_root_base = 0; new_root_base < alph_size; new_root_base++) {

        double prior = get_matrix_cell(root_base, new_root_base, prob_matrix);
        double prob = pruning_algorithm(
          alignment_col, 
          seq_names, 
          new_root_base, 
          child,
          substmatrix_array, 
          gap_support
        );

        // Implement minimum gap-handling strategy: choose
        // the base that minimizes the score.
        if (gap_support == MIN_GAPS) {

          if (is_leaf(child)) {
            char child_char;
            // Get the character at this leaf.
            int i = get_index_in_string_list(get_label(child), seq_names);
            child_char = alignment_col[i];
            // Is it a gap character?
            if ((child_char == '-') || (child_char == '.')) {
              if (prior < min) {
                // FIXME: Add a user-specified multiplicative factor.
                min = prior;
              }
            }
          }
        }

        // Normal case: just sum all the values.
        sum = sum + (prior * prob);

      }

      if ((gap_support == MIN_GAPS) && (min != HUGE_VAL)) {
        result = result * min;
      } 
      else {
        result = result * sum;
      }
      if (prob_matrix != NULL) {
        free_matrix(prob_matrix);
      }

    }

  }

  return result;
}

/****************************************************************************
 *  This function is our first implementation of the Felsentein pruning 
 *  algorithm for evaluating the log likelihood of a phylogentic tree for a 
 *  given site in an alignment. 
 *  It sums over all possible base assignments of the internal nodes 
 *  using the law of total probability and recursively calls itself until it 
 *  reaches the tree leaves. At the leaves paths through the tree whose final
 *  base doesn't match the appropriate sequence site are "pruned" 
 *  by returning a value of 0.0.
 ****************************************************************************/
static double old_pruning_algorithm(char* alignment_col,
                                STRING_LIST_T* seq_names,
                                char old_base,
                                TREE_T* tree,
                                EVOMODEL_T* model,
                                GAP_SUPPORT_T gap_support) {
  char new_base = 0;
  char* name = NULL;
  int a, b; // Indices for old base and new base 
  int c; // Index for child nodes
  int i;
  int alph_size = 0;
  int num_children = 0;
  double prior = 0.0;
  double prob = 0.0;
  double result = 0.0;
  double sum = 0.0; // Sum for normal case.
  double min = 0.0; // Minimum for min_gap strategy.
  double t = 0.0;
  TREE_T* child = NULL;
  MATRIX_T* prob_matrix = NULL;

  if (is_leaf(tree) == TRUE) {
    // End of recursion, get character at this leaf
    name = get_label(tree);
    i = get_index_in_string_list(name, seq_names);
    if (i == -1) {
      // No sequence matches the leaf label
      die("No sequence in the alignment matches the tree "
          "node labeled %s.\n", name);
    }
    new_base = alignment_col[i];
    if (old_base == new_base || new_base == 'N' || old_base == 'N') {
      result = 1.0;
    } else {
        result = 0.0;
    }

  } 
  else {
    result = 1.0;
    num_children = get_num_children(tree);
    alph_size = get_alph_size(ALPH_SIZE);
    for (c = 0; c < num_children; c++) {

      min = HUGE_VAL;
      sum = 0.0;
      child = get_nth_child(c, tree);
      t = get_length(child);
      prob_matrix = old_get_model_prob_matrix(model, t);
      a = alphabet_index(old_base, get_alphabet(FALSE));

      // Call recursively for each child
      for (b = 0; b < alph_size; b++) {

        new_base = get_alph_char(b);
        prior = get_matrix_cell(a, b, prob_matrix);
        prob = old_pruning_algorithm(
          alignment_col, 
          seq_names, 
          new_base, 
          child,
          model, 
          gap_support
        );

        // Normal case: just sum all the values.
        sum = sum + (prior * prob);
      }

      result = result * sum;

      if (prob_matrix != NULL) {
        free_matrix(prob_matrix);
      }
    }
  }

  return result;
}

/****************************************************************************
 *  Calculate the log likelihood of a phylogenetic tree
 *  at one site in the alignment.
 ****************************************************************************/
double site_likelihood(
  char* alignment_col,
  STRING_LIST_T* seq_names,
  TREE_T* tree,
  MODEL_TYPE_T model_type,
  ARRAY_T* priors,
  SUBSTMATRIX_ARRAY_T* substmatrix_array,
  GAP_SUPPORT_T gap_support) 
{
  assert(alignment_col != NULL);
  assert(seq_names != NULL);
  assert(tree != NULL);

  char* alphabet = get_alphabet(FALSE);
  int alpha_size = get_alph_size(ALPH_SIZE);
  int i = 0;
  double likelihood = 0.0;

  if (model_type == SINGLE_MODEL) {
    // Calculate likelihood directly from PSFM
    i = alphabet_index(*alignment_col, alphabet);
    likelihood = get_array_item(i, priors);

  } else if (model_type == AVERAGE_MODEL) {
    // Calculate likelihood as the *product* of the likelihoods
    // of the bases in the column.  This corresponds to *summing*
    // the log-likelihoods.
    likelihood = 1.0;
    while (*alignment_col != '\0') {
      i = alphabet_index(*alignment_col, alphabet);
      likelihood *= get_array_item(i, priors);
      alignment_col++;
    }

  } else {
    // Calcluate likelihood using Felsenstein pruning algorithm
    int root_base = 0;
    int start_base = 0;
    int end_base = alpha_size - 1;

    // TLB; Check if the root is labeled with a sequence.
    // If it is, this is the target genome.
    // Use the observed base only as the ancestor base.
    if (has_label(tree)) {			// target genome at root
      // TLB; Only compute conditional probability given target at root
      start_base = end_base = get_base_from_node(
        tree, 
        seq_names, 
        alignment_col, 
        alphabet, 
        alpha_size);
    } 

    for (root_base = start_base; root_base <= end_base; root_base++) {
      double prob = pruning_algorithm(
	alignment_col, 
	seq_names, 
	root_base, 
	tree, 
	substmatrix_array, 
	gap_support
      );
      double prior = get_array_item(root_base, priors);
      likelihood += (prob * prior);
    }

  }

  free_array(priors);

  return likelihood;
}

/****************************************************************************
 *  Calculate the log likelihood of a phylogenetic tree
 *  at one site in the alignment.
 ****************************************************************************/
static double old_site_likelihood(char* alignment_col,
  STRING_LIST_T* seq_names,
  TREE_T* tree,
  ARRAY_T* priors,
  EVOMODEL_T* model,
  GAP_SUPPORT_T gap_support) {

  char* a = NULL;
  int alpha_size = 0;
  int i, num_seq;
  double prob = 0.0;
  double prior = 0.0;
  double likelihood = 0.0;
  MODEL_TYPE_T model_type;

  assert(model != NULL);

  a = get_alphabet(FALSE);
  model_type = get_model_type(model);
  if (model_type == SINGLE_MODEL) {
    // Calculate likelihood directly from PSFM
    i = alphabet_index(*alignment_col, a);
    likelihood = get_array_item(i, priors);
  } else if (model_type == AVERAGE_MODEL) {
    // Calculate likelihood directly from PSFM averaged over 
    // bases in alignment column.
    num_seq = 0;
    while (*alignment_col != '\0') {
      i = alphabet_index(*alignment_col, a);
      likelihood += get_array_item(i, priors);
      alignment_col++;
      num_seq++;
    }
    likelihood = likelihood / num_seq;
  } else {
    // Calculate likelihood using Felsenstein pruning algorithm
    alpha_size = get_alph_size(ALPH_SIZE);
    for (i = 0; i < alpha_size; i++) {
      prob = old_pruning_algorithm(alignment_col, seq_names, get_alph_char(i), 
                               tree, model, gap_support);
      prior = get_array_item(i, priors);
      likelihood += (prob * prior);
    }
  }

  return likelihood;
}

/*************************************************************************
 * Calculate the log odds score score at a single position.
 *************************************************************************/
double slowly_compute_logodds(
  int            i,
  char*          alignment_col,
  ALIGNMENT_T*   alignment,
  STRING_LIST_T* seq_names,
  TREE_T*        tree,
  ARRAY_T*       bg_freqs,
  EVOMODEL_T*    fg_model,
  EVOMODEL_T*    bg_model,
  GAP_SUPPORT_T  gap_support
  )
{
  // Compute the foreground and background scores.
  get_alignment_col(i, alignment_col, alignment);
  double fg_result = old_site_likelihood(alignment_col, seq_names, tree,
					 bg_freqs, fg_model, gap_support);
  double bg_result = old_site_likelihood(alignment_col, seq_names, tree,
					 bg_freqs, bg_model, gap_support);

  // Compute the log-odds score.
  double log_odds = my_log2(fg_result / bg_result);
  return(log_odds);
}


/*************************************************************************
 * Calculate the log odds score and p-value for each motif-sized window 
 * at each site in the specified sequence of the alignment using 
 * the given evolutionary models.
 * Returns number of sites scored (not skipped).
 *************************************************************************/
int score_sequence_in_alignment(
  int             ref_seq_index,
  ALIGNMENT_T*    alignment,
  char*           motif_id,
  TREE_T*         tree,
  int             window_size,
  EVOMODEL_T**    models, 
  ARRAY_T*        bg_freqs,
  PSSM_T*         pssm,
  int*            coord_conv_table,
  GAP_SUPPORT_T   gap_support,
  double          gap_cost,
  double          pthresh,
  SCANNED_SEQUENCE_T* scanned_seq
) 
{
  assert(alignment != NULL);
  assert(models != NULL);

  // This function has two modes: fast and memory-intensive, or slow and not memory-intensive.
  // The function will operate in fast mode, unless the pssm argument is null.
  if (pssm == NULL) {
    assert(tree != NULL);
    assert(bg_freqs != NULL);
  } else {
    assert(tree == NULL);
    assert(bg_freqs == NULL);
  }

  // keep track of how many positions were scored
  int n_scored = 0;

  // maximum legal score in pv table
  int max_score = (pssm != NULL) ? get_pssm_pv_length(pssm) : 0;

  // JCH
  // This has been added so that the single mode scan
  // will scan the entire reference sequence and ignore
  // all gaps in the aligment

  ALIGNMENT_T*  new_alignment;

  if (get_model_type(models[0]) == SINGLE_MODEL) {

	// We need a new alignment file that has just the reference sequence
	char* seqname = get_seq_name(get_alignment_sequence(ref_seq_index, alignment));
	new_alignment = remove_alignment_gaps(seqname, alignment);
	// now remove all other seqs
	STRING_LIST_T* list = new_string_list();
	add_string(seqname, list);
	new_alignment = remove_alignment_seqs(list, new_alignment);
	ref_seq_index = 0;

	// and we need a dummy table for converting coordinates
	coord_conv_table = make_alignment_to_seq_table(ref_seq_index, new_alignment);
	//ALIGNMENT_T* temp = alignment;
	alignment = new_alignment;

  }

  // Print the bloody thing out and have a look
  //print_phylip_alignment(alignment, fopen("Test_alignment_output_file", "w"));

  // Prepare storage for the string representing 
  // the current alignment column.
  int alignment_length = get_alignment_length(alignment);
  int num_sequences = get_num_aligned_sequences(alignment);
  char* alignment_col = mm_malloc(sizeof(char) * (num_sequences + 1));

  // Prepare storage for the string representing the portion  
  // of the reference sequence within the window.
  char* window_seq = (char *) mm_malloc(sizeof(char) * (window_size + 1));
  window_seq[window_size] = '\0';
  SEQ_T* ref_seq = get_alignment_sequence(ref_seq_index, alignment);


  // Build a descriptive name for the sequence being scored.
  // If the user provided a sequence ID, use that.
  char* align_name = NULL;
  char* full_name = NULL;
  char* seq_name = NULL;

  seq_name = get_seq_name(ref_seq);
  align_name = get_alignment_name(alignment);
  full_name = (char*) mm_malloc(
    sizeof(char) * (strlen(align_name) + strlen(seq_name) + 2)
  );
  strcpy(full_name, align_name);
  strcat(full_name, ".");
  strcat(full_name, seq_name);


  // For each site in the alignment 
  int i;
  for (i = 0; i < alignment_length - window_size + 1; i++) {
    int j;

    // Fill in the sequence for this window.
    for (j = 0; j < window_size; j++) {
      window_seq[j] = get_seq_char(i + j, ref_seq);
    }

    
    // Decide whether to skip this window.
    BOOLEAN_T skip_window = FALSE;
    for (j = 0; j < window_size; j++) {
      if ((alignment_site_has_gaps(i + j, alignment) == TRUE) && 
        (gap_support == SKIP_GAPS)) {
        skip_window = TRUE;
        if (verbosity >= DUMP_VERBOSE) {
         fprintf(stderr, "Skipping sequence %s position %d (%s) due to gap at %d.\n",
         full_name, i, window_seq, i+j);
        }
        break;
      }
      if (alignment_site_ambiguous(i + j, alignment) == TRUE) {
        skip_window = TRUE;
        if (verbosity >= DUMP_VERBOSE) {
          fprintf(
            stderr, 
            "Skipping sequence %s position %d (%s) due to ambiguity at %d.\n",
            full_name, i, window_seq, i+j
          );
        }
        break;
      }
    }
    if (skip_window == TRUE) {
      if (verbosity >= DUMP_VERBOSE) {
        fprintf(
          stderr, 
          "Skipping sequence %s position %d (%s) due to ambiguity or gap.\n",
          full_name, i, window_seq
        );
      }
      continue;
    }

    //fprintf(stderr, "HERE 2\n");
    // Compute the total log-odds score across the window.
    double log_odds = 0.0;
    for (j = 0; j < window_size; j++) {

	//fprintf(stderr, "HERE 2.1\n");
      // Check for gaps at this site
      if ((alignment_site_has_gaps(i + j, alignment) == TRUE) && 
         (gap_support == FIXED_GAP_COST)) {
         log_odds -= gap_cost;
      }
	
      // Use the lookup table to compute log-odds.
      else if (pssm != NULL) {
	//fprintf(stderr, "HERE 2.2\n");
	//fprintf(stderr, "Model is : %s \n", get_model_type(models[0]));
        int alignment_col_index;
        if (get_model_type(models[0]) == SINGLE_MODEL) {
	  //fprintf(stderr, "HERE 2.2.1\n");
          alignment_col_index = alphabet_index(window_seq[j], get_alphabet(FALSE));
        } else {
	  
	  //fprintf(stderr, "HERE 2.2.2\n");
          get_alignment_col(i + j, alignment_col, alignment);
          alignment_col_index = hash_alignment_col(alignment_col, num_sequences);
        }
        //log_odds += get_matrix_cell(j, alignment_col_index, pssm);
        log_odds += get_pssm_score(j, alignment_col_index, pssm);
      }
	
      // If no lookup table is provided, compute log-odds the slow way.
      else {
	//fprintf(stderr, "HERE 2.3\n");
        log_odds += slowly_compute_logodds(
          i, alignment_col, alignment,
          get_species_names( alignment), 
          tree, bg_freqs,
          models[1], models[0], gap_support
        );
      }

      if (verbosity >= DUMP_VERBOSE) {
        fprintf(stderr, "%s %d+%d %c log_odds=%g\n",
        full_name, i, j, window_seq[j], log_odds);
      }
    }
    
    // Retrieve the p-value corresponding to this log-odds.
    double pvalue = 0.0;
    if (pssm != NULL) {
      if ((int)log_odds >= max_score) log_odds = (float)(max_score - 1);
      pvalue = get_pssm_pv((int) log_odds, pssm);
    }

    // Store information about this site.
    if ((pssm != NULL)) { 
      char strand = *motif_id;
      if ((strand != '-') && (strand != '+')) {
        strand = '.';
      }
      int start = 0;
      int stop = 0;
      if (strand == '-') {
        // Reverse the sense of start/stop for the reverse strand
        start = coord_conv_table[i + window_size - 1];
        stop = coord_conv_table[i];
      }
      else {
        start = coord_conv_table[i];
        stop = coord_conv_table[i + window_size - 1];
      }
      MATCHED_ELEMENT_T *element = 
        allocate_matched_element_with_score(start, stop, log_odds, pvalue, scanned_seq);
      add_scanned_sequence_scanned_element(scanned_seq);
      PATTERN_T* pattern = get_scanned_sequence_parent(scanned_seq);
      BOOLEAN_T added = add_pattern_matched_element(pattern, element);
      if (added == TRUE) {
        // Sub-set of raw sequence copied to matched element
        set_matched_element_sequence(element, window_seq);
      }
      else {
         // Element was rejected by pattern
        free_matched_element(element);
      }
    }

    n_scored++;
  }

  // Clean up memory.
  if (full_name != seq_name) {
    myfree(full_name);
  }
  myfree(alignment_col);
  myfree(window_seq);
  return(n_scored);
} // score_sequence_in_alignment




//~/Projects/MEME/trunk/src/motiph --bgfile intergenic_markov0 --output-pthresh 0.001 --model single --text --column-freqs empirical --pseudocount 0.01  --list alignments.txt yeast.tree ABF1.meme.txt > testout


/*************************************************************************
 * Calculate the Branch Length Score for each motif-sized window 
 * at each site in the specified sequence of the alignment using 
 * the given evolutionary models.
 * Returns number of sites scored (not skipped).
 *************************************************************************/
int bls_score_sequence_in_alignment(
  int             ref_seq_index,
  ALIGNMENT_T*    alignment,
  char*           motif_id,
  TREE_T*         tree,
  int             window_size,
  ARRAY_T*        bg_freqs,
  EVOMODEL_T**    models, 
  PSSM_T*         pssm,
  char*           inverse_motif_id,
  EVOMODEL_T**    inverse_models, 
  PSSM_T*         inverse_pssm,
  int*            coord_conv_table,
  GAP_SUPPORT_T   gap_support,
  double          gap_cost,
  double          pthresh,
  int          	  distThreshold,
  SCANNED_SEQUENCE_T* scanned_seq
) 
{

  // The first thing we do is get a single sequence scan score for all sequences
  // within the alignment. We use the existing functions and data structures to do this.
  int numSeqs = get_num_aligned_sequences(alignment);
  int i, j, t, e, temp_scored;
  int n_scored = 0;
  //int distThreshold = 20;

  //fprintf(stderr, "We have %d sequences in the alignment.\n", numSeqs);

  // We need to create some CISML crap to store each of the single scans results
  SCANNED_SEQUENCE_T **single_seq_scans = (SCANNED_SEQUENCE_T **) mm_malloc(numSeqs * sizeof(SCANNED_SEQUENCE_T *));
  // Using this appears to be causing problems so I need to create a 
  // a set of temporary CSMl crap instead
  // PATTERN_T* pat = get_scanned_sequence_parent(scanned_seq);
  PATTERN_T *pat = allocate_pattern(motif_id, motif_id);
  CISML_T *cisml = allocate_cisml("BLS", motif_id, "clustal-w alignment");
  set_cisml_site_pvalue_cutoff(cisml, pthresh);

  // Duplicate the above if we need to do inverse scans as well
  SCANNED_SEQUENCE_T **inverse_single_seq_scans = NULL;
  PATTERN_T *inverse_pat = NULL;
  CISML_T *inverse_cisml = NULL;
  if(inverse_models != NULL) { 
  	inverse_single_seq_scans = (SCANNED_SEQUENCE_T **) mm_malloc(numSeqs * sizeof(SCANNED_SEQUENCE_T *));
  	inverse_pat = allocate_pattern(inverse_motif_id, inverse_motif_id);
  	inverse_cisml = allocate_cisml("BLS", inverse_motif_id, "clustal-w alignment");
  	set_cisml_site_pvalue_cutoff(inverse_cisml, pthresh);
	//fprintf(stderr, "Inverse Data Structures Created: \n");
  }

  MATCHED_ELEMENT_T *** posElements = (MATCHED_ELEMENT_T ***) mm_malloc(numSeqs * sizeof(MATCHED_ELEMENT_T **));
  MATCHED_ELEMENT_T *** negElements = (MATCHED_ELEMENT_T ***) mm_malloc(numSeqs * sizeof(MATCHED_ELEMENT_T **));
  MATCHED_ELEMENT_T *** allElements = (MATCHED_ELEMENT_T ***) mm_malloc(numSeqs * sizeof(MATCHED_ELEMENT_T **));
  int *numHits = malloc( sizeof(int) * numSeqs); //int[numSeqs];

  for(i=0;i<numSeqs; i++) {
	single_seq_scans[i] = allocate_scanned_sequence("TEMP", "TEMP", pat);
	set_scanned_sequence_length(single_seq_scans[i], get_alignment_length(alignment));

	// Call the standard score_sequence_in_alignment function 
	// using a different reference sequence.
	// This returns the single sequence scans
	//double new_threshold = 0.001;
	double new_threshold = pthresh;
	if(i == ref_seq_index)
		new_threshold = pthresh;

	temp_scored = score_sequence_in_alignment( 
		i, 
		alignment, 
		motif_id, 
		NULL, 
		window_size, 
		models, 
		bg_freqs, 
		pssm, 
		coord_conv_table, 
		gap_support, 
		gap_cost, 
		new_threshold, 
		single_seq_scans[i] );

	posElements[i] = get_scanned_sequence_matched_elements(single_seq_scans[i]);
	numHits[i] = get_scanned_sequence_num_matched_elements(single_seq_scans[i]);

	//fprintf(stderr, "Positive motif scored ... %d \n", numHits[i]);

	// We don't run the inverse models on the reference sequence
	// The whole routine will be run a second time for this.
	// I know it is stupidlly inefficient, but the notation of whether something
	// is on the + or - strand is locked away higher in the CISML data structure
	// and I can't change it for individual matches.
	// If this ever changes then we can improve efficiency here

	if(inverse_models != NULL && i != ref_seq_index) {
		int numPos = numHits[i];
		// SCAN WITH INVERSE MODEL
		inverse_single_seq_scans[i] = allocate_scanned_sequence("TEMP", "TEMP", inverse_pat);
		set_scanned_sequence_length(inverse_single_seq_scans[i], get_alignment_length(alignment));
		
		//fprintf(stderr, "About to run Negative motif...\n");
		temp_scored = score_sequence_in_alignment( 
			i, 
			alignment, 
			inverse_motif_id, 
			NULL, 
			window_size, 
			inverse_models, 
			bg_freqs, 
			inverse_pssm, 
			coord_conv_table, 
			gap_support, 
			gap_cost, 
			new_threshold, 
			inverse_single_seq_scans[i] );


		negElements[i] = get_scanned_sequence_matched_elements(inverse_single_seq_scans[i]);
		numHits[i] += get_scanned_sequence_num_matched_elements(inverse_single_seq_scans[i]);

		//fprintf(stderr, "Negative motif scored ...%d \n", temp_scored);
		//fprintf(stderr, "Total motif scored ...%d \n", numHits[i]);
		// Now Copy pos and neg into all
		
		//fprintf(stderr, "Allocate space for %d  elements in total \n", numHits[i]);
		allElements[i] = (MATCHED_ELEMENT_T **) mm_malloc(numHits[i] * sizeof(MATCHED_ELEMENT_T *));

		for(e=0; e<numHits[i]; e++) {
			if(e<numPos) {
				allElements[i][e] = posElements[i][e];
			} else {
				allElements[i][e] = negElements[i][e-numPos];
			}
		}

	} else {
		allElements[i] = posElements[i];
	}


	if(numHits[i] > 0)
	   sort_matched_elements(TRUE, numHits[i], allElements[i]);

	// fprintf(stderr, "Scored %d positions in sequence %d with %d matches \n", temp_scored, i, numHits[i]);
  }
  

  // Ave a geezer
  /*
  for(i=0;i<numSeqs; i++) {
	fprintf(stderr, "In sequence %d \n",  i);

	MATCHED_ELEMENT_T ** elements = get_scanned_sequence_matched_elements(single_seq_scans[i]);
	for(j=0;j<numHits[i]; j++) {
		int start = get_matched_element_start(elements[j]);
		int stop = get_matched_element_stop(elements[j]);
		fprintf(stderr, "   %d - %d \n",  start, stop);
	}
  }
  */


  // Now we iterate over the results in the target sequence
  // Vars for reuse
  STRING_LIST_T * list_of_matching_seqs;
  int numOfMatches;
  float conservedBL, BLScore;

  STRING_LIST_T* speciesLabels = get_species_names(alignment);

  MATCHED_ELEMENT_T ** targetElements = allElements[ref_seq_index];
  // MATCHED_ELEMENT_T ** targetElements = get_scanned_sequence_matched_elements(single_seq_scans[ref_seq_index]);

  float totalLength = get_total_length(tree);

  int* cumulativeGapCounts = get_cumulative_gap_count(ref_seq_index, alignment);


  // We need a set of tables that map positions in each sequence back to the
  // alignment structure so that we can compare the positions.
  
  int** seq_to_aln_tables = (int **) mm_malloc((numSeqs) * sizeof(int*));
  for(i=0;i<numSeqs; i++) {
     seq_to_aln_tables[i] = make_seq_to_alignment_table(i, alignment);
  }


  for( t=0; t<numHits[ref_seq_index]; t++) {

	// Get the coordinates of this hit
	int targStartSeq = get_matched_element_start(targetElements[t]);
	int targStopSeq = get_matched_element_stop(targetElements[t]);

	// Translate them to positions in the alignment
	int targStart = seq_to_aln_tables[ref_seq_index][targStartSeq];
	int targStop = seq_to_aln_tables[ref_seq_index][targStopSeq];

	// We store matching seq labels here
	list_of_matching_seqs = new_string_list();
	numOfMatches = 0;

	// Iterate over all the other sequences and look for 
	// hits within the threshold

	for(i=0;i<numSeqs; i++) {
		if(i != ref_seq_index) {
			//MATCHED_ELEMENT_T ** elements = get_scanned_sequence_matched_elements(single_seq_scans[i]);
			MATCHED_ELEMENT_T ** elements = allElements[i];

			for(j=0;j<numHits[i]; j++) {

				// Use the start position as a flag
				int daStartSeq = get_matched_element_start(elements[j]);

				if(daStartSeq > 0) {
				   // Now we need to find the shortest distance
				   int daStopSeq = get_matched_element_stop(elements[j]);

				   // Translate these into positions in the alignment
				   int daStart = seq_to_aln_tables[i][daStartSeq];
				   int daStop = seq_to_aln_tables[i][daStopSeq];

				   int numGaps1 = abs( cumulativeGapCounts[targStart] - cumulativeGapCounts[daStart] );
				   int diff1 = abs(targStart - daStart) - numGaps1;
				   int smallest = diff1;
				   int numGaps2 = abs( cumulativeGapCounts[targStop] - cumulativeGapCounts[daStop] );
				   int diff2 = abs(targStop - daStop) - numGaps2;
				   if(diff2 < smallest)
					smallest = diff2;
				   int numGaps3 = abs( cumulativeGapCounts[targStart] - cumulativeGapCounts[daStop] );
				   int diff3 = abs(targStart - daStop) - numGaps3;
				   if(diff3 < smallest)
					smallest = diff3;
				   int numGaps4 = abs( cumulativeGapCounts[targStop] - cumulativeGapCounts[daStart] );
				   int diff4 = abs(targStop - daStart) - numGaps4;
				   if(diff4 < smallest)
					smallest = diff4;

				   if(smallest < distThreshold) {
					char* label = get_nth_string(i, speciesLabels);
					//fprintf(stderr, " Match in seq [%d] %s : %d - %d \n", i, label, start, targStart);
					add_string(label, list_of_matching_seqs);
					numOfMatches++;
					// Remove that match from further consideration
					// we set the start to -1.0 as a flag
					set_matched_element_start(elements[j], -1);
					
					break;
				   }
				}
			}
		} 
	}

	// If the number of matches in comparative species is greater than 0
	// Then we add in the target genome to the list of labels
	if(numOfMatches > 0) {
		char* label = get_nth_string(ref_seq_index, speciesLabels);
		add_string(label, list_of_matching_seqs);
	}

	// Now we have a list of leaf labels for seqs that have the motif
	// We retrieve the branch length of the sub-tree they define

	float subtreeLength = get_subtree_length(tree, list_of_matching_seqs);
	
	//fprintf(stderr, " Subtree length %f  and total tree length %f \n", subtreeLength, totalLength);
	float score = (float) subtreeLength / totalLength;
	//fprintf(stderr, " Final Score %f \n", score);
	

	// Store information about this site.
	
	MATCHED_ELEMENT_T *element = allocate_matched_element_without_inversion(
		targStartSeq, 
		targStopSeq, 
		get_matched_element_sequence(targetElements[t]), 
		scanned_seq
  );
  set_matched_element_score(element, score);
	set_matched_element_pvalue(element, 1-score);
  add_scanned_sequence_matched_element(scanned_seq, element);

	/*
	float p_val_sub = get_matched_element_pvalue(targetElements[t]);
	
	if(score > 0)
		p_val_sub = (1/score) * get_matched_element_pvalue(targetElements[t]);
	else
		p_val_sub = (2*numSeqs) * get_matched_element_pvalue(targetElements[t]);
	
	if(p_val_sub > 1)
		p_val_sub = 1; 
	
	set_matched_element_pvalue(element,  p_val_sub );
	*/
	
     	n_scored++;

	// Free up some memory 
	free_string_list(list_of_matching_seqs);
	
  }


  // Clean up memory.
  for(i=0;i<numSeqs; i++) {
     myfree(single_seq_scans[i]);
     myfree(seq_to_aln_tables[i]);
  }
  myfree(single_seq_scans);
  myfree(seq_to_aln_tables);
  free_cisml(cisml);
  myfree(cumulativeGapCounts);

  return(n_scored);
} // bls_score_sequence_in_alignment

/*************************************************************************
 * Get the (index of) the base in the sequence in the alignment column
 * corresponding to the given tree node.
 *************************************************************************/
int get_base_from_node(
  TREE_T* node,
  STRING_LIST_T* seq_names,
  char* alignment_col, 
  char* alphabet, 
  int alph_size
)
{
  char* name = get_label(node);
  int i = get_index_in_string_list(name, seq_names);
  if (i == -1) {
    // No sequence matches the leaf label
    die("No sequence in the alignment matches the tree "
          "node labeled %s.\n", name);
  }
  char new_char = alignment_col[i];
  int base = alphabet_index(new_char, alphabet);
  assert(base < alph_size);

  return(base);
} // get_sequence_index_from_node
