/*************************************************************************
 * FILE: read-mhmm.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 8-12-97
 * PROJECT: MHMM
 * COPYRIGHT: 1997-2008, WSN
 * DESCRIPTION: Functions for reading HMMs in MHMM format.
 *************************************************************************/
#include <assert.h>
#include <err.h>
#include <errno.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "build-hmm.h"
#include "mhmm-state.h"
#include "alphabet.h"
#include "array.h"
#include "read-mhmm.h"
#include "utils.h"
#include "xml-util.h"

/*************************************************************************
 * Read a single HMM state in MHMM format.
 *************************************************************************/
#define BUF_LENGTH 100
static void read_mhmm_state(
   FILE * const   model_file,
   const int      i_state, 
   MHMM_STATE_T * this_state
) {
  int state_num;    /* State index read from file. */
  int alph_size;    /* Number of characters in the alphabet. */
  int i_alph;       /* Index of the current letter. */
  int num_scan;     /* Number of items scanned by fscanf(). */
  char buffer[BUF_LENGTH]; /* Buffer for reading strings from the file. */
  ATYPE this_emit;

  /* Read the state index and make sure it's the right one. */
  num_scan = fscanf(model_file, " State: %d", &state_num);
  assert(num_scan == 1);
  assert(state_num == i_state);
#ifdef NDEBUG
  num_scan = i_state; /* Avoid compiler warning. */
#endif

  /* Read other state characteristics. */
  num_scan = fscanf(model_file, " type: %s", &(buffer[0]));
  assert(num_scan == 1);
  this_state->type = 
    convert_enum_type_str(buffer, INVALID_STATE, STATE_STRS, NUM_STATE_T);
  num_scan = fscanf(model_file, " i_motif: %d", &(this_state->i_motif));
  assert(num_scan == 1);
  num_scan = fscanf(model_file, " motif_id: %s", this_state->motif_id);
  assert(num_scan == 1);
  num_scan = fscanf(model_file, " i_position: %d", &(this_state->i_position));
  assert(num_scan == 1);
  num_scan = fscanf(model_file, " id_char: %c", &(this_state->id_char));
  assert(num_scan == 1);
  assert(this_state->id_char != '\0');
  num_scan = fscanf(model_file, " num_sites: %lf", &(this_state->num_sites));
  assert(num_scan == 1);

  /* Allocate space for the emission distributions. */
  this_state->emit = allocate_array(get_alph_size(ALL_SIZE));
  this_state->emit_odds = allocate_array(get_alph_size(ALL_SIZE));

  /* Read emission distribution from non-start-or-end states. */
  alph_size = get_alph_size(ALPH_SIZE);
  if ((this_state->type != START_STATE) && (this_state->type != END_STATE)) {
    fscanf(model_file, " emit: ");
    for (i_alph = 0; i_alph < alph_size; i_alph++) {
      num_scan = fscanf(model_file, PROB_SCAN, &this_emit);
      set_array_item(i_alph, this_emit, this_state->emit);
      assert(num_scan == 1);
    }

    /* Normalize the probabilities. */
    normalize(0.00001, this_state->emit);

    /* Calculate probabilities for ambiguous characters. */
    fill_in_ambiguous_chars(FALSE, this_state->emit);
  }

  /* Assign 1.0 probabilities to start and end emissions. */
  else {
    for (i_alph = 0; i_alph < get_alph_size(ALL_SIZE); i_alph++) {
      set_array_item(i_alph, 1.0, this_state->emit);
    }
  }
}


/*************************************************************************
 * Read an HMM transition matrix in MHMM format.
 *************************************************************************/
static void read_mhmm_trans(
  FILE *    model_file, 
  const int num_states,
  MATRIX_T* trans
) {
  int       num_rows;      /* Number of rows and columns (equal). */
  int       num_cols;
  int       i_row;         /* Row and column indices. */
  int       i_col;
  int       num_scan;      /* Number of items scanned by fscanf(). */
#ifndef NO_RENORMALIZE
  int       one_col;       /* Index of the column in which 1.0 appears. */
  MTYPE     row_sum;       /* Sum of given row (except for the 1.0 entry). */
  BOOLEAN_T found_one;     /* Does this column contain an entry of 1.0? */
  BOOLEAN_T found_nonzero; /* ... or an entry other than 0.0? */
#endif
  MTYPE     matrix_item;


  num_scan = fscanf(model_file, " Transition probability matrix (%d x %d):",
		    &num_rows, &num_cols);
  assert(num_scan == 2);
  assert((num_states == num_rows) && (num_states == num_cols));

  /* Read the matrix. */
  for (i_row = 0; i_row < num_states; i_row++) {
    for (i_col = 0; i_col < num_states; i_col++) {
      num_scan = fscanf(model_file, PROB_SCAN, &matrix_item);
      set_matrix_cell(i_row, i_col, matrix_item, trans);
      assert(num_scan == 1);
    }
  }

#ifndef NO_RENORMALIZE
  /* Values close to 1.0 get easily pushed to 1.0, so fix 'em. */
  for (i_row = 0; i_row < num_states; i_row++) {

    /* Find out if there's a 1.0 and something nonzero in this row. */
    found_one = FALSE;
    found_nonzero = FALSE;
    one_col = 0; /* To prevent compiler warning. */
    for (i_col = 0; i_col < num_states; i_col++) {
      if (get_matrix_cell(i_row, i_col, trans) == 1.0) {
        found_one = TRUE;
        one_col = i_col;
      } else if (get_matrix_cell(i_row, i_col, trans) != 0.0) {
	      found_nonzero = TRUE;
      }
      if (found_one && found_nonzero) {
	      break;
      }
    }

    /* If so, add up everything else and subtract it from the 1.0 entry. */
    if (found_one && found_nonzero) {
      fprintf(
        stderr, 
        "Renormalizing row %d of the transition matrix: ",
	      i_row
      );
      row_sum = 0.0;
      for (i_col = 0; i_col < num_states; i_col++) {
        if (i_col != one_col) {
          row_sum += get_matrix_cell(i_row, i_col, trans);
        }
      }
      incr_matrix_cell(i_row, one_col, -row_sum, trans);
      fprintf(
        stderr, 
        "trans[%d,%d]=%g\n", 
        i_row, one_col, 
	      get_matrix_cell(i_row, one_col, trans)
      );
    }
  }
#endif

  num_scan = fscanf(model_file, " End of MHMM");
  assert(num_scan == 0);
}


/***********************************************************************
 * Read the emission probabilities for a state from XML.
 * Return value: pointer to ARRAY_T structure containing 
 * emission probabilities.
 ***********************************************************************/
static void read_state_emission_prob_from_xml(
  xmlXPathContextPtr xpath_ctxt,
  ARRAY_T* emission_probs,
  int state_index // chosen state
) {

  const int MAX_XPATH_EXPRESSION = 200;
  char xpath_expression[MAX_XPATH_EXPRESSION];

  int alph_size = get_alph_size(ALPH_SIZE);
  int i_alph = 0;

  // Build the XPATH expression to get emission probabilities for chosen state.
  snprintf(
    xpath_expression,
    MAX_XPATH_EXPRESSION,
    "/MEME_HMM/states/state[@index='%d']/"
    "emission_probabilities/alphabet_array/value",
    state_index
  );

  // Count the number of emission probabilities for the state.
  xmlXPathObjectPtr xpath_obj = xpath_query(xpath_ctxt, xpath_expression);
  int num_values = (xpath_obj->nodesetval ? xpath_obj->nodesetval->nodeNr : 0);
  xmlXPathFreeObject(xpath_obj);

  // Start and end states won't have any emissions in XML.
  if (num_values > 0) {

    // If not zero, the number of emission probabilities 
    // should match the alphabet size.
    assert(num_values == alph_size);

    // XML doesn't enforce any order on the emission probability values,
    // so force reading emission probability values in alphabet order.
    xmlNodePtr currValueNode = NULL;
    for (i_alph = 0; i_alph < alph_size; i_alph++) {
      // Build the XPATH expression to get emission probability for 
      // a character in the chosen state.
      snprintf(
        xpath_expression,
        MAX_XPATH_EXPRESSION,
        "/MEME_HMM/states/state[@index='%d']/"
        "emission_probabilities/alphabet_array/value[@letter_id='letter_%c']",
        state_index, 
        get_alph_char(i_alph)
      );
      // Decode from node set to numeric value for emission prob.
      xpath_obj = xpath_query(xpath_ctxt, xpath_expression);
      currValueNode = xpath_obj->nodesetval->nodeTab[0];
      ATYPE value = xmlXPathCastNodeToNumber(currValueNode);
      set_array_item(i_alph, value, emission_probs);
      xmlXPathFreeObject(xpath_obj);
    }

    /* Normalize the probabilities. */
    normalize(0.00001, emission_probs);

    /* Calculate probabilities for ambiguous characters. */
    fill_in_ambiguous_chars(FALSE, emission_probs);
  }
  else {
    /* Assign 1.0 probabilities to start and end emissions. */
    for (i_alph = 0; i_alph < get_alph_size(ALL_SIZE); i_alph++) {
      set_array_item(i_alph, 1.0, emission_probs);
    }
  }
}

/***********************************************************************
 * Read the transition probabilities from XML.
 * transition probabilities.
 ***********************************************************************/
static void read_mhmm_trans_from_xml(
  xmlXPathContextPtr xpath_ctxt,
  MATRIX_T* trans,
  int from_state_index // chosen state
) {

  const int MAX_XPATH_EXPRESSION = 200;
  char xpath_expression[MAX_XPATH_EXPRESSION];

  // Build the XPATH expression to get emission probabilities for chosen state.
  snprintf(
    xpath_expression,
    MAX_XPATH_EXPRESSION,
    "/MEME_HMM/states/state[@index='%d']/transition",
    from_state_index
  );

  // Count the number of transition probabilities for the state.
  xmlXPathObjectPtr xpath_obj = xpath_query(xpath_ctxt, xpath_expression);
  int num_transitions = 
    (xpath_obj->nodesetval ? xpath_obj->nodesetval->nodeNr : 0);

  // Some states won't have any transtion probabilities.
  if (num_transitions > 0) {

    xmlNodePtr transition_node = NULL;
    int i_transition = 0;
    for (i_transition = 0; i_transition < num_transitions; i_transition++) {
      transition_node = xpath_obj->nodesetval->nodeTab[i_transition];
      xmlChar* property = 
        read_xml_node_property(transition_node, "to_state_index");
      int to_state_index = atoi((char *) property);
      xmlFree(property);
      ATYPE value = xmlXPathCastNodeToNumber(transition_node);
      set_matrix_cell(from_state_index, to_state_index, value, trans);
    }

  }

  xmlXPathFreeObject(xpath_obj);

}

/***********************************************************************
 * Read the states from XML.
 ***********************************************************************/
void read_mhmm_states_from_xml(
  xmlXPathContextPtr xpath_ctxt, 
  MHMM_T* the_mhmm
) {

  xmlXPathObjectPtr xpath_obj = NULL;
  xmlChar* property = NULL;

  // Select all states
  xpath_obj = xpath_query(xpath_ctxt, "/MEME_HMM/states/state");

  xmlNodePtr curr_state = NULL;
  int i_state = 0;
  for(i_state = 0; i_state < the_mhmm->num_states; i_state++) {

    curr_state = xpath_obj->nodesetval->nodeTab[i_state];
    if (curr_state == NULL) {
      die("Error: missing state %d\n", i_state);
    }
  
    // Get the state index attribute.
    // This will almost certainly be the same as i_state, 
    // but the order of the states in the XML is not guaranteed.
    property = read_xml_node_property(curr_state, "index");
    int state_index = atoi((char *) property);
    xmlFree(property);

    // Get the state type attribute
    property = read_xml_node_property(curr_state, "type");
    (the_mhmm->states[state_index]).type = convert_enum_type_str(
      (char *) property, 
      INVALID_STATE, STATE_STRS,
			NUM_STATE_T
    );
    xmlFree(property);

    // Get the num_sites attribute
    property = read_xml_node_property(curr_state, "num_sites");
    (the_mhmm->states[state_index]).num_sites = atof((char *) property);
    xmlFree(property);

    // Get the i_motif attribute
    property = read_xml_node_property(curr_state, "motif_index");
    (the_mhmm->states[state_index]).i_motif = atoi((char *) property);
    xmlFree(property);

    // Get the position attribute
    property = read_xml_node_property(curr_state, "position");
    (the_mhmm->states[state_index]).i_position = atoi((char *) property);
    xmlFree(property);

    // Get the motif_id attribute
    property = read_xml_node_property(curr_state, "motif_id");
    strncpy(
      (the_mhmm->states[state_index]).motif_id, 
      (char *) property, 
      MAX_MOTIF_ID_LENGTH + 1
    );
    xmlFree(property);

    // Get the id_char attribute
    property = read_xml_node_property(curr_state, "id_char");
    (the_mhmm->states[state_index]).id_char = property[0];
    xmlFree(property);

    /* Allocate space for the emission distributions. */
    the_mhmm->states[state_index].emit = 
      allocate_array(get_alph_size(ALL_SIZE));
    the_mhmm->states[state_index].emit_odds = 
      allocate_array(get_alph_size(ALL_SIZE));

    // Read the emission probabilities for the state.
    read_state_emission_prob_from_xml(
      xpath_ctxt, 
      the_mhmm->states[state_index].emit,
      state_index
    );

    // Read the transtion probabilities for the state.
    read_mhmm_trans_from_xml(
      xpath_ctxt, 
      the_mhmm->trans,
      state_index
    );

    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(stderr, " %d", i_state);
    }
  }
  
  xmlXPathFreeObject(xpath_obj);

}

/*************************************************************************
 * Create an HMM struct from a parsed XML struct.
 *
 * IN: mhmm_doc
 * OUT: the_hmm
 *************************************************************************/
static void read_mhmm_from_xml(
  xmlXPathContextPtr xpath_ctxt, 
  MHMM_T** the_hmm
) {

  xmlXPathObjectPtr xpath_obj = NULL;
  xmlChar* property = NULL;

  xpath_obj = xpath_query(xpath_ctxt, "/MEME_HMM/states");
  property = read_xml_node_property(
    xpath_obj->nodesetval->nodeTab[0], 
    "num_states"
  );
  int num_states = atoi((char *) property);
  xmlFree(property);

  // Allocate the HMM.
  allocate_mhmm(num_states, the_hmm);
  (*the_hmm)->num_states = num_states;
  strncpy(
    (*the_hmm)->alphabet, 
    get_alphabet(FALSE /* don't include ambigs */), 
    (*the_hmm)->alph_size
  );

  // Fill in the number and length of spacers since we already
  // have the states node in hand.
  property = read_xml_node_property(
    xpath_obj->nodesetval->nodeTab[0], 
    "num_spacers"
  );
  (*the_hmm)->num_spacers = atoi((char *) property);
  xmlFree(property);

  property = read_xml_node_property(
    (xmlNodePtr) xpath_obj->nodesetval->nodeTab[0], 
    "spacer_states"
  );
  (*the_hmm)->spacer_states = atoi((char *) property);
  xmlFree(property);
  xmlXPathFreeObject(xpath_obj);

  // Set the alphabet information
  (*the_hmm)->alph_size = get_alph_size(ALPH_SIZE /* Don't include ambigs */);
  strcpy((*the_hmm)->alphabet, get_alphabet(FALSE /* Don't include ambigs */));
  (*the_hmm)->ambigs = get_alph_size(AMBIG_SIZE);

  // Read and set the HMM type
  xpath_obj = xpath_query(xpath_ctxt, "/MEME_HMM");
  property = read_xml_node_property(xpath_obj->nodesetval->nodeTab[0], "type");
  (*the_hmm)->type = 
    convert_enum_type_str((char *) property, INVALID_HMM, HMM_STRS, NUM_HMM_T);
  xmlFree(property);

  // Read the optional fields
  property = xmlGetProp(
    xpath_obj->nodesetval->nodeTab[0], 
    BAD_CAST "description"
  );
  if (property != NULL) {
    (*the_hmm)->description = mm_malloc(sizeof(char) * strlen((char *) property) + 1);
    strcpy((*the_hmm)->description, (char *) property);
    xmlFree(property);
  }
  property = xmlGetProp(
    xpath_obj->nodesetval->nodeTab[0], 
    BAD_CAST "motif_file"
  );
  if (property != NULL) {
    (*the_hmm)->motif_file = mm_malloc(sizeof(char) * strlen((char *) property) + 1);
    strcpy((*the_hmm)->motif_file, (char *) property);
    xmlFree(property);
  }
  property = xmlGetProp(
    xpath_obj->nodesetval->nodeTab[0], 
    BAD_CAST "sequence_file"
  );
  if (property != NULL) {
    (*the_hmm)->motif_file = mm_malloc(sizeof(char) * strlen((char *) property) + 1);
    strcpy((*the_hmm)->sequence_file, (char *) property);
    xmlFree(property);
  }

  xmlXPathFreeObject(xpath_obj);

  // Initialize the transition matrix
  init_matrix(0.0, (*the_hmm)->trans);

  // Read all the states into the HMM
  read_mhmm_states_from_xml(xpath_ctxt, *the_hmm);

  // Record the number of motifs.
  (*the_hmm)->num_motifs = get_num_motifs(*the_hmm);

  /* N.B. We can't read directly into the model's background because
     it contains the alphabet plus ambiguous characters.  Hence
     the following shenanigans.  */
  {
    ARRAY_T* background = read_bg_freqs_from_xml(xpath_ctxt);
    int i = 0;
    for (i = 0; i < get_alph_size(ALPH_SIZE); i++) {
      set_array_item(i, get_array_item(i, background), (*the_hmm)->background);
    }
    free_array(background);
  }

  // Fill in ambiguous characters.
  fill_in_ambiguous_chars(FALSE, (*the_hmm)->background);


}

/*************************************************************************
 * Read HMM from an XML file.
 *
 * IN: model_file - 
 * OUT: the_hmm - 
 *************************************************************************/
BOOLEAN_T read_mhmm_xml_file(
  char*    model_filename, /* The name of the file to be read from. */
  MHMM_T** the_hmm         /* The HMM (not pre-allocated). */
) {

  xmlParserCtxtPtr ctxt;         // The parser
  xmlDocPtr mhmm_doc;            // The document tree
  xmlXPathContextPtr xpath_ctxt; // XPath context

  ctxt = xmlNewParserCtxt();
  if (ctxt == NULL) {
    die("Failed to create XML parser.\n");
  }

  // Parse and validate the file.
  mhmm_doc = xmlCtxtReadFile(
    ctxt, 
    model_filename, 
    NULL,  // Encoding
    XML_PARSE_DTDVALID | XML_PARSE_NOERROR | XML_PARSE_NOWARNING
  );

  // Did it parse?
  if (mhmm_doc == NULL) {
    // It may be an old plain-text file.
    xmlFreeParserCtxt(ctxt);
    xmlCleanupParser();
    if (verbosity > NORMAL_VERBOSE) {
      fprintf(stderr, "Failed to parse %s as an XML document.\n", model_filename);
    }
    return FALSE;
  } else {
    // Did it validate as a MEME_HMM document?
    if (ctxt->valid == 0 || strcmp((char *) ctxt->intSubName, "MEME_HMM") != 0) {
      xmlFreeDoc(mhmm_doc);
      xmlFreeParserCtxt(ctxt);
      xmlCleanupParser();
	    die("%s is not a valid MEME HMM XML document.\n", model_filename);
    }
  }
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "File %s is a valid MEME HMM XML file.\n", model_filename);
  }

  // Set up XPath context from parsed XML
  xpath_ctxt = xmlXPathNewContext(mhmm_doc);

  // Read and set alphabet
  read_alphabet_from_xml(xpath_ctxt);

  // Create an HMM struct from the XML struct.
  read_mhmm_from_xml(xpath_ctxt, the_hmm);

  /* Compute the state-by-state incoming and outgoing transitions. */
  compute_ins_and_outs(*the_hmm, FALSE);

  /* free up the resulting document */
  xmlXPathFreeContext(xpath_ctxt);
  xmlFreeDoc(mhmm_doc);
  xmlFreeParserCtxt(ctxt);
  xmlCleanupParser();

  return TRUE;

}
/*************************************************************************
 * void read_mhmm
 *
 * Read an HMM in MHMM format.
 *
 * IN: model_file - 
 * OUT: the_hmm - 
 *************************************************************************/
void read_mhmm(
  char*    hmm_filename,  /* The file to be read from. */
  MHMM_T** the_hmm        /* The HMM (not pre-allocated). */
) {

  char buffer[BUF_LENGTH];    /* Buffer for reading from file. */
  char version[BUF_LENGTH];   /* Buffer for reading from file. */
  int  type;          /* Type of model (MHMM or HMMER). */
  int  num_states;    /* Number of states in the model. */
  int  num_spacers;   /* Number of spacers in the model. */
  int  spacer_states; /* Number of states per spacer. */
  int  alph_size;     /* Number of characters in the alphabet. */
  int  i_state;       /* Index of the current state. */
  ARRAY_T* motif_widths;	// Array to hold motif widths 

  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Reading HMM: ");
  }

  // Begin by assuming that it is an MHMM XML file.
  BOOLEAN_T result = read_mhmm_xml_file(hmm_filename, the_hmm);
  if (result) {
    // If that worked, we're done!
    return;
  }

  // Open the file for reading.
  FILE* hmm_file;
  if (open_file(hmm_filename, "r", TRUE, "HMM", "HMM", &hmm_file) == 0) {
    exit(1);
  }

  /* Skip the header. */
  int num_scan = 0; // Number of items scanned by fscanf().
  while (num_scan != 1) {
    if (fgets(buffer, BUF_LENGTH, hmm_file) == NULL) {
      die("Can't find Meta-MEME version number.\n");
    }
    num_scan = sscanf(buffer, "MHMM v%s", version);
  }
  assert(strlen(version) < BUF_LENGTH);

  // Read the header information.
  num_scan = fscanf(hmm_file, " type: %s", &(buffer[0]));
  assert(num_scan == 1);
  type = convert_enum_type_str(buffer, INVALID_HMM, HMM_STRS, NUM_HMM_T);
  num_scan = fscanf(hmm_file, " num_states: %d", &num_states);
  assert(num_scan == 1);
  num_scan = fscanf(hmm_file, " num_spacers: %d", &num_spacers);
  assert(num_scan == 1);
  num_scan = fscanf(hmm_file, " spacer_length: %d", &spacer_states);
  assert(num_scan == 1);
  num_scan = fscanf(hmm_file, " alph_size: %d", &alph_size);
  assert(num_scan == 1);

  // Set up the alphabet.
  num_scan = fscanf(hmm_file, " alphabet: %s", buffer);
  assert(num_scan == 1);
  set_alphabet(verbosity, buffer);

  // Allocate the HMM.
  allocate_mhmm(num_states, the_hmm);

  // Store information that we have read so far.
  (*the_hmm)->type = type;
  (*the_hmm)->num_states = num_states;
  (*the_hmm)->num_spacers = num_spacers;
  (*the_hmm)->spacer_states = spacer_states;
  (*the_hmm)->alph_size = alph_size;
  strcpy((*the_hmm)->alphabet, get_alphabet(FALSE));
  (*the_hmm)->ambigs = get_alph_size(AMBIG_SIZE);

  // Read the background distribution.
  num_scan = fscanf(hmm_file, " background: ");
  assert(num_scan == 0);

  /* N.B. We can't read directly into the model's background because
     it contains the alphabet plus ambiguous characters.  Hence
     the following shenanigans.  */
  {
    int i;
    ARRAY_T* background;

    background = allocate_array(get_alph_size(ALPH_SIZE));
    read_array(hmm_file, background);
    for (i = 0; i < get_alph_size(ALPH_SIZE); i++) {
      set_array_item(i, get_array_item(i, background), (*the_hmm)->background);
    }
    free_array(background);
  }

  // Extend the distribution to account for ambiguous characters.
  fill_in_ambiguous_chars(FALSE, (*the_hmm)->background);

  // Read the optional fields.
  if (fscanf(hmm_file, " description: %[^\n]", buffer) == 1) {
    copy_string(&((*the_hmm)->description), buffer);
  }
  if (fscanf(hmm_file, " motif_file: %s", buffer) == 1) {
    copy_string(&((*the_hmm)->motif_file), buffer);
  }
  if (fscanf(hmm_file, " sequence_file: %s", buffer) == 1) {
    copy_string(&((*the_hmm)->sequence_file), buffer);
  }

  /* Read in each state. */
  motif_widths = allocate_array(num_states);  // tlb; need to set widths
  for (i_state = 0; i_state < num_states; i_state++) {
    if (verbosity >= NORMAL_VERBOSE)
      fprintf(stderr, "%d\r", i_state);		// tlb; replace " " with \r
    read_mhmm_state(hmm_file, i_state, &((*the_hmm)->states[i_state]));
    // Save motif width when you know it.
    if ((*the_hmm)->states[i_state].type == END_MOTIF_STATE) {
      set_array_item(
        (*the_hmm)->states[i_state].i_motif,  // Motif index.
        (*the_hmm)->states[i_state].i_position + 1,	// Width of motif. 
        motif_widths
      );
    }
  }
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "\n");
  }

  // Record the motif widths in the states.
  for (i_state = 0; i_state < num_states; i_state++) {
    if ((*the_hmm)->states[i_state].type == START_MOTIF_STATE ||
      (*the_hmm)->states[i_state].type == MID_MOTIF_STATE || 
      (*the_hmm)->states[i_state].type == END_MOTIF_STATE) {
      (*the_hmm)->states[i_state].w_motif = 
        get_array_item((*the_hmm)->states[i_state].i_motif, motif_widths);
    } else {
      // All other state types have width 1.
      (*the_hmm)->states[i_state].w_motif = 1;	// start, end, spacer
    }
  }
  free_array(motif_widths);

  // Record the number of motifs.
  (*the_hmm)->num_motifs = get_num_motifs(*the_hmm);

  /* Read the transition matrix. */
  read_mhmm_trans(hmm_file, num_states, (*the_hmm)->trans);
  assert(verify_trans_matrix(FALSE, num_states, (*the_hmm)->trans));

  fscanf(hmm_file, " End of MHMM");

  /* Compute the state-by-state incoming and outgoing transitions. */
  compute_ins_and_outs(*the_hmm, FALSE);

  fclose(hmm_file);

} // read_mhmm

