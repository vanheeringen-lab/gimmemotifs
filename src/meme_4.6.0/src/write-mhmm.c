/*************************************************************************
 * FILE: write-mhmm.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 8-12-97
 * PROJECT: MHMM
 * COPYRIGHT: 1997-2008, WSN
 * DESCRIPTION: Functions for writing HMMs in MHMM format.
 *************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "metameme.h"
#include "mhmm-state.h"
#include "write-mhmm.h"
#include "mhmm-dtd.h"

/*************************************************************************
 * void write_mhmm_state
 *
 * Write a single MHMM state.
 *************************************************************************/
static void write_mhmm_state(
  const int i_state,              /* Index of this state. */
  const MHMM_STATE_T * the_state, /* One HMM state. */
  const int alph_size,            /* Size of the alphabet. */
  FILE *model_file                /* Output file. */
) {
  int i_alph;
  int i_column;

  fprintf(model_file, "State: %d\n", i_state);
  fprintf(
    model_file, 
    "    type: %s\n", 
	  convert_enum_type(the_state->type, STATE_STRS, NUM_STATE_T)
  );
  fprintf(model_file, "    i_motif: %d\n", the_state->i_motif);
  fprintf(model_file, "    motif_id: %s\n", the_state->motif_id);
  fprintf(model_file, "    i_position: %d\n", the_state->i_position);
  fprintf(model_file, "    id_char: %c\n", the_state->id_char);
  fprintf(model_file, "    num_sites: %g\n", the_state->num_sites);

  /* Start and end states are non-emitting. */
  if ((the_state->type == START_STATE) || (the_state->type == END_STATE))
    return;

  /* Print the emission distribution */
  fprintf(model_file, "    emit: ");
  for (i_alph = 0, i_column = -1; i_alph < alph_size; i_alph++) {
    /* Don't make the output too wide. */
    if (i_column == 5) {
      fprintf(model_file, "\n          ");
      i_column = 0;
    } else {
      i_column++;
    }

    fprintf(model_file, "%9.7f ", get_array_item(i_alph, the_state->emit));
  }
  fprintf(model_file, "\n");
}

/*************************************************************************
 * void write_mhmm_transp
 *
 * Write a transition probability matrix in MHMM format.
 *************************************************************************/
static void write_mhmm_transp(
  const MHMM_T * the_hmm,
  FILE * model_file
) {
  int    i_state;  /* Current indices in the matrix. */
  int    j_state;  
  double value;    /* Current value from the matrix. */

  fprintf(
    model_file, 
    "Transition probability matrix (%d x %d):\n",
	  the_hmm->num_states, 
    the_hmm->num_states
  );
  for (i_state = 0; i_state < the_hmm->num_states; i_state++) {
    fprintf(model_file, "    ");
    for (j_state = 0; j_state < the_hmm->num_states; j_state++) {
      value = get_matrix_cell(i_state, j_state, the_hmm->trans);

      if (almost_equal(value, 0.0, 0.000001)) {
	      fprintf(model_file, "0 ");
      } else if (almost_equal(value, 1.0, 0.000001)) {
	      fprintf(model_file, "1 ");
      } else {
	      fprintf(model_file, "%8.6g ", value);
      }
    }
    fprintf(model_file, "\n");
  }
}


/*************************************************************************
 * write_mhmm
 *
 * See .h file for description.
 *************************************************************************/
void write_mhmm(
  VERBOSE_T local_verbosity,
  MHMM_T*   the_hmm,    /* The HMM. */
  FILE*     model_file  /* The file to be written to. */
) {
  int i_state;   /* The index of the current state. */
  int i_alph;    /* Index into the alphabet. */

  if (local_verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Writing MHMM: ");
  }

  /* Print the header. */
  fprintf(model_file, "MHMM v%s\n", VERSION);
  fprintf(model_file, "type: %s\n", 
	  convert_enum_type(the_hmm->type, HMM_STRS, NUM_HMM_T));
  fprintf(model_file, "num_states: %d\n", the_hmm->num_states);
  fprintf(model_file, "num_spacers: %d\n", the_hmm->num_spacers);
  fprintf(model_file, "spacer_length: %d\n", the_hmm->spacer_states);
  fprintf(model_file, "alph_size: %d\n", the_hmm->alph_size);
    

  /* Print the alphabet, minus ambiguous characters. */
  fprintf(model_file, "alphabet: ");
  for (i_alph = 0; i_alph < the_hmm->alph_size; i_alph++)
    fprintf(model_file, "%c", the_hmm->alphabet[i_alph]);
  fprintf(model_file, "\n");

  // Print the background distribution.
  fprintf(model_file, "background: ");
  print_sub_array(0, the_hmm->alph_size, the_hmm->background, 8, 7, TRUE,
		  model_file);

  // Print the optional fields.
  if (the_hmm->description != NULL) {
    fprintf(model_file, "description: %s\n", the_hmm->description);
  }
  if (the_hmm->motif_file != NULL) {
    fprintf(model_file, "motif_file: %s\n", the_hmm->motif_file);
  }
  if (the_hmm->sequence_file != NULL) {
    fprintf(model_file, "sequence_file: %s\n", the_hmm->sequence_file);
  }

  /* Write each state. */
  for (i_state = 0; i_state < the_hmm->num_states; i_state++) {
    if (local_verbosity >= NORMAL_VERBOSE) {
      fprintf(stderr, "%d ", i_state);
    }
    write_mhmm_state(
      i_state, 
      &(the_hmm->states[i_state]),
		  the_hmm->alph_size, 
      model_file
    );
  }
  if (local_verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "\n");
  }
  
  /* Write the transition probability matrix. */
  write_mhmm_transp(the_hmm, model_file);
  fprintf(model_file, "End of MHMM\n\n");
}

/*************************************************************************
 * void write_mhmm_to_file
 *
 * Write a hidden markov model in MHMM format to file with the given name.
 *
 * IN:  the_hmm - the HMM
 *      filename - the name of the file to be created.
 * OUT: ---
 * SIDE EFFECT: An MHMM-style HMM is written to a file named filename.
 *              If the filename is "-", HMM is written to stdout.
 *************************************************************************/
void write_mhmm_to_file(
  MHMM_T* the_hmm,/* The HMM. */
  char* filename  /* The name of the file to be written to. */
) {

    FILE* model_file = NULL;
    BOOLEAN_T close_model_file = TRUE;
    if (strcmp(filename, "-") == 0) {
      model_file = stdout;
      close_model_file = FALSE;
    }
    else {
      model_file = fopen(filename, "w");
      if (model_file == NULL) {
        die("Unable to open %s for writing model\n", filename);
      }
    }
    write_mhmm(verbosity, the_hmm, model_file);
    if (close_model_file) {
      fclose(model_file);
    }
}

/*************************************************************************
 * write_states_xml
 *
 * Write out XML describing HMM state
 *************************************************************************/
static void write_states_xml(
  VERBOSE_T local_verbosity,
  MHMM_T*   the_hmm,    /* The HMM. */
  FILE*     model_file  /* The file to be written to. */
) {

  int i_state;   /* Row index of the current state. */
  int j_state;   /* Column index to state. */
  int i_alph;    /* Alphabet index. */

  // Begin states element
  fprintf(model_file, "<states");
  fprintf(model_file, " num_states='%d'", the_hmm->num_states);
  fprintf(model_file, " num_spacers='%d'", the_hmm->num_spacers);
  fprintf(model_file, " spacer_states='%d'", the_hmm->spacer_states);
  fprintf(model_file, ">\n");

  // Print out each state element
  for (i_state = 0; i_state < the_hmm->num_states; i_state++) {

    MHMM_STATE_T the_state = the_hmm->states[i_state];

    if (local_verbosity >= NORMAL_VERBOSE) {
      fprintf(stderr, "%d ", i_state);
    }
    fprintf(model_file, "<state");
    fprintf(model_file, " id='state_%d'", i_state);
    fprintf(model_file, " index='%d'", i_state);
    fprintf(
      model_file, 
      " type='%s'", 
      convert_enum_type(the_state.type, STATE_STRS, NUM_STATE_T)
    );
    fprintf(model_file, " motif_index='%d'", the_state.i_motif);
    fprintf(model_file, " motif_id='%s'", the_state.motif_id);
    fprintf(model_file, " position='%d'", the_state.i_position);
    fprintf(model_file, " id_char='%c'", the_state.id_char);
    fprintf(model_file, " num_sites='%g'", the_state.num_sites);
    fprintf(model_file, ">\n");

    // Print emission probabilities element.
    // Start and end states are non-emitting.
    if ((the_state.type != START_STATE) && (the_state.type != END_STATE)) {
      fprintf(model_file, "<emission_probabilities>\n");
      fprintf(model_file, "<alphabet_array>\n");
      for (i_alph = 0; i_alph < the_hmm->alph_size; i_alph++) {
        fprintf(
          model_file, 
          "<value letter_id='letter_%c'>%9.7f</value>\n", 
          the_hmm->alphabet[i_alph],
          get_array_item(i_alph, the_state.emit)
        );
      }
      fprintf(model_file, "</alphabet_array>\n");
      fprintf(model_file, "</emission_probabilities>\n");
    }

    // Print transition elements.
    for (j_state = 0; j_state < the_hmm->num_states; j_state++) {
      ATYPE value = get_matrix_cell(i_state, j_state, the_hmm->trans);
      if (almost_equal(value, 0.0, 0.000001)) {
        // Don't print elements for zero transtions.
      }
      else if (almost_equal(value, 1.0, 0.000001)) {
	      fprintf(
          model_file, 
          "<transition to_state_index='%d'>1</transition>\n",
          j_state
        );
      } else {
	      fprintf(
          model_file, 
          "<transition to_state_index='%d'>%8.6g</transition>\n",
          j_state,
          value
        );
      }
    }

    // Close the state element
    fprintf(model_file, "</state>\n");

  }

  // Close the states element
  fprintf(model_file, "</states>\n");

}

/*************************************************************************
 * write_mhmm_xml
 *
 * Write out XML describing HMM
 *************************************************************************/
void write_mhmm_xml(
  VERBOSE_T local_verbosity,
  MHMM_T*   the_hmm,    /* The HMM. */
  FILE*     model_file  /* The file to be written to. */
) {

  int i_alph;    /* Index into the alphabet. */

  if (local_verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Writing MHMM: ");
  }

  // Print the DTD.
  fprintf(model_file, "%s", mhmm_dtd);

  // Start the document boduy.
  const int head_skip = 18; // Skip SVN keyword
  const int tail_skip = 2; // Skip close of SVN keyword
  const char *release_date = ARCHIVE_DATE + head_skip;
  const int release_date_width = strlen(release_date) - tail_skip;
  fprintf(model_file, "<!-- Document Body -->\n");
  fprintf(model_file, "<MEME_HMM");
  fprintf(model_file, " version='%s'", VERSION);
  fprintf(model_file, " release_date='%.*s'", release_date_width, release_date);
  fprintf(
    model_file, 
    " type='%s'", 
    convert_enum_type(the_hmm->type, HMM_STRS, NUM_HMM_T)
  );
  fprintf(model_file, " motif_file='%s'", the_hmm->motif_file);
  if (the_hmm->sequence_file != NULL) {
    fprintf(model_file, " sequencce_file='%s'", the_hmm->sequence_file);
  }
  if (the_hmm->description != NULL) {
    fprintf(model_file, " description='%s'", the_hmm->description);
  }
  fprintf(model_file, " >\n");
  fprintf(model_file, "<command_line>\n");
  fprintf(model_file, "\n");
  fprintf(model_file, "</command_line>\n");

  // Start the alphabet element
  fprintf(model_file, "<alphabet");
  fprintf(
    model_file, 
    " id='%s' length='%d'>\n", 
    the_hmm->alph_size == 4 ? "nucleotide" : "amino_acid",
    the_hmm->alph_size
  );
  for (i_alph = 0; i_alph < the_hmm->alph_size; i_alph++) {
    fprintf(
      model_file, 
      "<letter id='letter_%c' symbol='%c' />\n", 
      the_hmm->alphabet[i_alph],
      the_hmm->alphabet[i_alph]
    );
  }
  // Close the alphabet element
  fprintf(model_file, "</alphabet>\n");

  // Start the background frequencies element
  fprintf(model_file, "<background_frequencies>\n");
  fprintf(model_file, "<alphabet_array>\n");
  for (i_alph = 0; i_alph < the_hmm->alph_size; i_alph++) {
    fprintf(
      model_file, 
      "<value letter_id='letter_%c'>", 
      the_hmm->alphabet[i_alph]
    );
    ATYPE item = get_array_item(i_alph, the_hmm->background);
    fprintf(model_file, "%9.7f", item);
    fprintf(model_file, "</value>\n");
  }
  fprintf(model_file, "</alphabet_array>\n");
  // Close the background frequencies element
  fprintf(model_file, "</background_frequencies>\n");

  // Print the states element.
  write_states_xml(local_verbosity, the_hmm, model_file);

  // Close the document.
  fprintf(model_file, "</MEME_HMM>\n");
}

/*************************************************************************
 * write_mhmm_xml_to_file
 *
 * Write out XML describing HMM to the named file.
 *************************************************************************/
void write_mhmm_xml_to_file(
  MHMM_T* the_hmm, /* The HMM. */
  char* filename   /* The name of the file to be written to. */
) {
    FILE* model_file = NULL;
    BOOLEAN_T close_model_file = TRUE;
    if (strcmp(filename, "-") == 0) {
      model_file = stdout;
      close_model_file = FALSE;
    }
    else {
      model_file = fopen(filename, "w");
      if (model_file == NULL) {
        die("Unable to open %s for writing model\n", filename);
      }
    }
    write_mhmm_xml(verbosity, the_hmm, model_file);
    if (close_model_file) {
      fclose(model_file);
    }
}
