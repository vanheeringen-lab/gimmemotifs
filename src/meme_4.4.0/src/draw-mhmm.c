/*************************************************************************
 * FILE: draw-mhmm.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 3/20/99
 * PROJECT: MHMM
 * COPYRIGHT: 1997-2008, WSN
 * DESCRIPTION: Produce a graphviz version of a Meta-MEME HMM.
 * 
 * For more information about graphviz, see
 *
 *      http://www.research.att.com/sw/tools/graphviz
 *
 *************************************************************************/
#include "utils.h"
#include "mhmm-state.h"
#include "read-mhmm.h" 
#include "simple-getopt.h"
#include <stdio.h>
#include <assert.h>
#include <string.h>

#define FONT "Helvetica"             /* Font family used for all labels. */
#define FONT_SIZE 24                 /* Font size used for all labels. */
#define BEGIN_END_COLOR "black"      /* Color of begin and end nodes. */
#define SPACER_COLOR "blue"          /* Color of spacer nodes. */
#define MOTIF_BRIGHTNESS 0.8         /* Brightness of all motif nodes. */
#define LOW_MOTIF_HUE 0.1            /* Range of colors for motif nodes. */
#define HIGH_MOTIF_HUE 0.9
#define LOW_MOTIF_SATURATION 0.5     /* Range of saturation for motifs. */
#define HIGH_MOTIF_SATURATION 1.0
#define MIN_TRANS 0.0                /* Don't print edges with
                                        transitions less than this. */
#define EDGE_WEIGHT 10.0             /* Multiplier for transition
                                        probabilities to get edge
                                        weights. */
#define NODE_SIZE 0.1                /* Size of one node. */

// Shapes of the nodes.
#define BEGIN_END_SHAPE "ellipse"
#define SPACER_SHAPE "circle"
#define MOTIF_SHAPE "box"

#define MAX_MOTIF_LENGTH 100

// ID used in place of a motif ID for non-motif states.
#define NON_MOTIF_ID "---"

/*************************************************************************
 * A simple helper function to find the index of the start state of a
 * motif, given the index of the end state.
 *************************************************************************/
static int find_start_of_motif
  (MHMM_T* the_hmm,
   int     end_of_motif)
{
  int i_state;

  for (i_state = end_of_motif; 
       (the_hmm->states[i_state]).type != START_MOTIF_STATE;
       i_state--) {
    assert(i_state != 0);
  }
  //fprintf(stderr, "Motif from %d to %d.\n", i_state, end_of_motif);
  return(i_state);
}

/*************************************************************************
 * Find the consensus of a given motif.
 *************************************************************************/
static void get_motif_consensus
  (MHMM_T* the_hmm,
   int     start_of_motif,
   char*   consensus)
{
  int           i_motif;
  MHMM_STATE_T* this_state;
  
  for (i_motif = 0; ; i_motif++) {
    this_state = &(the_hmm->states[start_of_motif + i_motif]);
    consensus[i_motif] = choose_consensus(FALSE, this_state);

    if (this_state->type == END_MOTIF_STATE) {
      break;
    }
  }
  consensus[i_motif+1] = '\0';
}

/*************************************************************************
 * Choose the color of a motif, based upon its index.
 *************************************************************************/
static char* choose_motif_color
  (int  i_motif,
   int  num_motifs)
{
  static char return_value[100];
  float hue;
  float saturation;
  
  hue = LOW_MOTIF_HUE
    + ((HIGH_MOTIF_HUE - LOW_MOTIF_HUE)
       * ((float)i_motif / (float)(num_motifs - 1)));
  saturation = LOW_MOTIF_SATURATION
    + ((HIGH_MOTIF_SATURATION - LOW_MOTIF_SATURATION)
       * ((float)i_motif / (float)(num_motifs - 1)));

  /* Put the current settings into the return string. */
  sprintf(return_value, "%g %g %g", hue, saturation, MOTIF_BRIGHTNESS);
  return(return_value);
}

/*************************************************************************
 * Determine whether a given state is unreachable.
 * 
 * Currently only valid for spacer states and start motif states.
 *************************************************************************/
static BOOLEAN_T is_used
  (MHMM_T* the_hmm,
   int     i_state)
{
  int i_prev; // Index of the previous state.
  MHMM_STATE_T* this_state = &(the_hmm->states[i_state]);

  myassert(TRUE, i_state >= 0, "Negative state index (%d).", i_state);
  myassert(TRUE, i_state < the_hmm->num_states, 
	   "State index too large (%d > %d).", i_state, the_hmm->num_states);

  switch (this_state->type) {
  case START_STATE :
  case END_STATE : 
    return(TRUE);
    break;
  case MID_MOTIF_STATE :
  case END_MOTIF_STATE :
    // Recursively call on the single preceding state.
    return(is_used(the_hmm, this_state->itrans_in[0]));
    break;

  case START_MOTIF_STATE :
    // Recursively call on all preceding states.
    for (i_prev = 0; i_prev < this_state->ntrans_in; i_prev++) {
      if (is_used(the_hmm, this_state->itrans_in[i_prev])) {
	return(TRUE);
      }
    }
    return(FALSE);
    break;
  case SPACER_STATE :
    // Handle this case below
    break;
  case INVALID_STATE :
    die("Invalid state found at %d.\n", i_state);
  }

  // Move to the beginning of the spacer.
  while (1) {

    if (this_state->itrans_in[0] != i_state) {
      i_prev = this_state->itrans_in[0]; // Avoid self-loops.
    } else if (this_state->ntrans_in > 1) {
      i_prev = this_state->itrans_in[1];
    } else {
      // If there is only an incoming self-loop, this state is not used.
      return(FALSE);
    }

    this_state = &(the_hmm->states[i_prev]);
    if (this_state->type == SPACER_STATE) {
      i_state = i_prev;
    } else {
      break;
    }
  }
  this_state = &(the_hmm->states[i_state]);

  // Does this state have any incoming transitions (besides self-loop)?
  return(this_state->ntrans_in > 1);
}

/*************************************************************************
 * Print a single state (node) of the graph.
 *************************************************************************/
static void print_state
  (MHMM_T* the_hmm,
   int     state_num,
   char*   motif_id,
   char*   label,
   char*   shape,
   char*   color,
   FILE*   outfile)
{
  if (is_used(the_hmm, state_num)) {
    if (strcmp(motif_id, NON_MOTIF_ID) != 0) {
      // Print the motif ID and the consensus.
      if (strlen(label) > 0) {
	fprintf(outfile,
		"  state%d [label=\"%s : %s\",shape=%s,color=\"%s\",size=\"%f,%f\"];\n",
		state_num, motif_id, label, shape, color, NODE_SIZE, 
		NODE_SIZE);
      }
      // Print only the motif ID.
      else {
	fprintf(outfile,
		"  state%d [label=\"%s\",shape=%s,color=\"%s\",size=\"%f,%f\"];\n",
		state_num, motif_id, shape, color, NODE_SIZE, NODE_SIZE);
      }
    }
    else {
	fprintf(outfile,
		"  state%d [label=\"%s\",shape=%s,color=\"%s\",size=\"%f,%f\"];\n",
		state_num, label, shape, color, NODE_SIZE, NODE_SIZE);
    }
  }
}

/*************************************************************************
 * Print a single edge of the graph.
 *************************************************************************/
static void print_trans
  (MHMM_T* the_hmm,
   int     i_from,
   int     i_to, 
   PROB_T  this_trans,
   FILE*   outfile)
{

  if (is_used(the_hmm, i_from) && is_used(the_hmm, i_to)) {
    fprintf(outfile, "   state%d -> state%d [weight=%f,label=\"%.2g\"];\n",
	    i_from, i_to, this_trans * EDGE_WEIGHT, this_trans);
  }
}

/*************************************************************************
 * Produce a graphviz graph representing the given HMM.
 *************************************************************************/
static void draw_mhmm
  (MHMM_T*   the_hmm,
   float     min_trans,       // Minimum transition to display.
   BOOLEAN_T print_consensus, // Print consensus along with motif ID?
   FILE*     outfile)
{
  int           i_state;
  int           num_states;
  int           num_motifs;
  MHMM_STATE_T* this_state;
  char          consensus[MAX_MOTIF_LENGTH];

  // Count the number of motifs in this model.
  num_motifs = get_num_motifs(the_hmm);

  if (verbosity > NORMAL_VERBOSE) {
    fprintf(stderr, "State ");
  }

  /* Print the header. */
  fprintf(outfile, "digraph G {\n");
  fprintf(outfile, "  rankdir=LR;\n");
  fprintf(outfile, "  node[fontsize=%d,font=%s];\n", FONT_SIZE, FONT);
  fprintf(outfile, "  edge[fontsize=%d,font=%s];\n", FONT_SIZE, FONT);

  /* Print the nodes. */
  num_states = the_hmm->num_states;
  for (i_state = 0; i_state < num_states; i_state++) {
    this_state = &(the_hmm->states[i_state]);

    if (verbosity > NORMAL_VERBOSE) {
      fprintf(stderr, " %d", i_state);
    }

    /* Start state. */
    switch (this_state->type) {
    case START_STATE :
      print_state(the_hmm, 
		  i_state, 
		  NON_MOTIF_ID,
		  "begin",
		  BEGIN_END_SHAPE,
		  BEGIN_END_COLOR,
		  outfile);
      break;

    /* End state. */
    case END_STATE :
      print_state(the_hmm,
		  i_state,
		  NON_MOTIF_ID,
		  "end",
		  BEGIN_END_SHAPE,
		  BEGIN_END_COLOR,
		  outfile);
      break;

    /* Spacer state. */
    case SPACER_STATE :
      print_state(the_hmm,
		  i_state,
		  NON_MOTIF_ID,
		  "",
		  SPACER_SHAPE,
		  SPACER_COLOR,
		  outfile);
      break;

    /* Motif state. */
    case START_MOTIF_STATE :

      // Make sure the motif is connected to something.
      if (is_used(the_hmm, i_state)) {

	if (print_consensus) {
	  get_motif_consensus(the_hmm, i_state, consensus);
	} else {
	  consensus[0] = '\0';
	}
	print_state(the_hmm,
		    i_state,
		    this_state->motif_id,
		    consensus,
		    MOTIF_SHAPE,
		    choose_motif_color(this_state->i_motif, num_motifs), 
		    outfile);
      }
      break;

    case MID_MOTIF_STATE :
    case END_MOTIF_STATE :
      break;

    default :
      die("Invalid state found at %d.\n", i_state);
    }
  }

  if (verbosity > NORMAL_VERBOSE) {
    fprintf(stderr, "\nEdge ");
  }


  /* Print the edges. */
  {
    int    i_from;
    int    i_to;
    int    this_from;
    PROB_T this_trans;

    
    for (i_from = 0; i_from < num_states; i_from++) {
      this_from = i_from;
      this_state = &(the_hmm->states[i_from]);

      switch (this_state->type) {

      case END_MOTIF_STATE :
	this_from = find_start_of_motif(the_hmm, i_from);

      case START_STATE :
      case END_STATE :
      case SPACER_STATE :

	for (i_to = 0; i_to < num_states; i_to++) {
	  this_trans = get_matrix_cell(i_from, i_to, the_hmm->trans);
	  if (this_trans > min_trans) {
	    if (verbosity > NORMAL_VERBOSE) {
	      fprintf(stderr, " %d->%d", this_from, i_to);
	    }
	    print_trans(the_hmm, this_from, i_to, this_trans, outfile);
	  }
	}
	break;

      case START_MOTIF_STATE :
      case MID_MOTIF_STATE :
	break;

      default :
	die("Invalid state found at %d.\n", i_state);
      }
    }
  }
  if (verbosity > NORMAL_VERBOSE) {
    fprintf(stderr, "\n");
  }


  fprintf(outfile, "}\n");
}

VERBOSE_T verbosity = INVALID_VERBOSE;

/*************************************************************************
 * int main
 *************************************************************************/
int main (int argc, char *argv[])
{
  char *    hmm_filename = NULL;     // File containing the HMM.
  FILE *    hmm_file;
  MHMM_T *  the_hmm;                 // The HMM itself.
  float     min_trans = MIN_TRANS;   // Minimum transition prob to print.
  BOOLEAN_T print_consensus = FALSE; // Print the consensus in each motif?
  
  // Define command line options.
  static cmdoption options[] = {
    {"mintrans", REQUIRED_VALUE},
    {"consensus", NO_VALUE},
    {"verbosity", REQUIRED_VALUE}
  };
  int option_count = 3;
  int option_index = 0;

  // Define the usage message.
  char      usage[1000] = "";
  strcat(usage, "USAGE: draw-mhmm [options] <HMM file>\n");
  strcat(usage, "\n");
  strcat(usage, "   Options:\n");
  strcat(usage, "     --mintrans <prob> (default=0.0)\n");
  strcat(usage, "     --consensus\n");
  strcat(usage, "     --verbosity 1|2|3|4|5 (default=2)\n");
  strcat(usage, "\n");

  simple_setopt(argc, argv, option_count, options);

  // Parse the command line.
  while (1) { 
    int c = 0;
    char * option_value = NULL;
    char * option_name = NULL;
    const char * message = NULL;

    // Read the next option, and break if we're done.
    c = simple_getopt(&option_name, &option_value, &option_index);
    if (c == 0) {
      break;
    } else if (c < 0) {
    	simple_getopterror(&message);
      die("Error processing command line options (%s)\n", message);
    }

    if (strcmp(option_name, "mintrans") == 0) {
      min_trans = atof(option_value);
    } else if (strcmp(option_name, "consensus") == 0) {
      print_consensus = TRUE;
    } else if (strcmp(option_name, "verbosity") == 0) {
      verbosity = (VERBOSE_T)atoi(option_value);
    }
  }

  // Read the single required argument.
  if (option_index + 1 != argc) {
    fprintf(stderr, "%s", usage);
    exit(1);
  }
  hmm_filename = argv[option_index];

  /* Read the model. */
  read_mhmm(hmm_filename, &the_hmm);

  /* Print the HMM. */
  draw_mhmm(the_hmm, min_trans, print_consensus, stdout);
  free_mhmm(the_hmm);

  return(0);
}
