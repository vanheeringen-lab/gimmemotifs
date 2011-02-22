/**************************************************************************
 * FILE: mhmms.c
 * AUTHOR: William Stafford Noble, Timothy L. Bailey
 * CREATE DATE: 8-13-97
 * PROJECT: MHMM
 * VERSION: $Revision: 1.3 $
 * COPYRIGHT: 1998-2002, WNG, 2001-2002, TLB
 * DESCRIPTION: Search a database of sequences using a motif-based HMM.
 **************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "utils.h"       // Generic utilities. 
#include "matrix.h"      // Routines for floating point matrices. 
#include "array.h"       // Routines for floating point arrays. 
#include "metameme.h"    // Global metameme functions. 
#include "alphabet.h"    // The amino acid / nucleotide alphabet. 
#include "fasta-io.h"    // Read FASTA format sequences. 
#include "meme-io.h"     // Read background freq. file.
#include "mhmm-state.h"  // HMM data structure. 
#include "read-mhmm.h"   // HMM input/output. 
#include "log-hmm.h"     // HMM log/log-odds conversion. 
#include "fitevd.h"	     // Extreme value distribution routines. 
#include "rdb-matrix.h"  // For reading background files 
#include "match.h"       // Find sequence-to-model match from traceback matrix.
#include "dp.h"          // Dynamic programming routines.
#include "pssm.h"	       // Motif-scoring routines.
#include "seq.h"         // Sequence manipulation routines
#include "mhmms.h"

#ifndef DEBUG
#define DEBUG 0
#endif

/*************************************************************************
 * The following block of global variables and #define's concerns
 * storing sequence scores and descriptions for later sorting and
 * printing.
 *************************************************************************/
#define ALLOCATE_BLOCK 1000 // Number of sequences to allocate at once.
#define PV_WIDTH 11	    // Number of digits in one p-value.
#define SCORE_WIDTH 7       // Number of digits in one score.
#define SCORE_PRECISION 2   // Number of digits after zero.
#define LENGTH_WIDTH 10     // Number of digits to print sequence length.
#define MAX_SEQ_ID_LENGTH 40    // Maximum allowed ID length.
#define GC_WIN_SIZE 500	    // Size of GC window left and right of match.

/* scaling and unscaling of scores for non-local EVDs */
#define ALPHA 2.6	    // Kludge factor for non-local EVDs.
#define BETA 1		    // Kludge factor for non-local EVDs.
#define SCALE(s, m)		(pow((s) - (m) + BETA, ALPHA))
#define UNSCALE(s, m)		(pow((s), 1/ALPHA) + (m) - BETA)

static int     num_stored = 0; // Number of sequences stored so far.
static char**  stored_IDs;     // Stored sequence IDs.
static PROB_T* stored_keys;    // Stored sequence keys.
static PROB_T* stored_scores;  // Stored sequence scores.
static int*    stored_lengths; // Stored sequence lengths.
static double* stored_gcs;     // Stored match local GC-contents.
static char**  stored_seqs;    // Stored sequence scores and descriptions.
static char**  stored_traces;  // Stored Viterbi traces.
// Store the first part of the match line and the subsequent motif lines
// separately, so we can insert the correct E-value later.
static char**  stored_match_gffs;    // GFF entry for 'match' line.
static char**  stored_motif_gffs;    // GFF entry for 'motif' lines.
static int     longest_ID = 0;       // Length of the longest ID.
static int     num_allocated = 0;    // Number of sequences allocated.
BOOLEAN_T      skipped_last = FALSE; // Did we skip the previous sequence?

/*************************************************************************
*	See .h file.
 *************************************************************************/
int get_num_stored() { return(num_stored); }
 
/*************************************************************************
 * Print parameters associated with the extreme value distribution.
 *************************************************************************/
static void print_distr_parameters(
  EVD_SET evd_set,			// Distributions.
  SCORE_SET score_set,			// Scores.
  FILE* out_stream                      // Output stream.
)
{
  int i;
  fprintf(out_stream, "\nDerived parameters for E-value calculations\n");
  fprintf(out_stream, " type of distribution: %s\n", evd_set.exponential ?
    "Exponential" : "Extreme Value");
  fprintf(out_stream, " source of distribution: %s\n", evd_set.negatives_only ?
    "Synthetic sequence scores" : "Actual sequence scores");
  fprintf(out_stream, " status: %s\n", evd_set.msg);

  if (evd_set.n > 0) {
    for (i=0; i<evd_set.n; i++) {
      if (evd_set.n > 1) {		// Length ranges.
	fprintf(out_stream, " min_t: %.3f\n", evd_set.evds[i].min_t);
	fprintf(out_stream, " max_t: %.3f\n", evd_set.evds[i].max_t);
      }
      if (evd_set.exponential) {
	fprintf(out_stream, " mu1: %.3g\n", evd_set.evds[i].mu1);
        if (!evd_set.negatives_only) {	// No second component?
	  fprintf(out_stream, " mu2: %.3g\n", evd_set.evds[i].mu2);
	  fprintf(out_stream, " sigma2: %.3g\n", evd_set.evds[i].sigma2);
	  fprintf(out_stream, " c: %.3g\n", evd_set.evds[i].c);
        }
	fprintf(out_stream, " number of scores used for estimation: %d\n", evd_set.evds[i].n);
      } else {			// Extreme value parameters.
	fprintf(out_stream, " lambda: %.2g\n", evd_set.evds[i].lambda);
	fprintf(out_stream, " K: %.2g\n", evd_set.evds[i].K);
        fprintf(out_stream, " number of non-outliers: %d\n", evd_set.non_outliers);
      }
    }
    fprintf(out_stream, " N (E-value = N * p-value): %d\n", evd_set.N);
    fprintf(out_stream, " mininum E-value: %.3g\n", evd_set.min_e);
    fprintf(out_stream, " number of E-values < 1.0: %d\n", evd_set.outliers);
    fprintf(out_stream, " sum(-log E-values < 1.0): %.3g\n", -evd_set.sum_log_e);
  }

  fprintf(out_stream, " number of matches: %d\n", score_set.n);
  fprintf(out_stream, " db length: %.0f\n", score_set.total_length);
  if (evd_set.exponential) {
    fprintf(out_stream, " average minimum p-value: %.3g\n",
	    score_set.avg_min_pvalue);
    fprintf(out_stream, " average motif width: %.3g\n", score_set.avgw);
    /*
    fprintf(out_stream, " probability of a false hit: %.3g\n", score_set.phit);
    fprintf(out_stream, " probability of a false hit within max-gap: %.3g\n",
	    score_set.wimp);
    fprintf(out_stream, " E[distance between hits]: %.3g\n", score_set.egap);
    fprintf(out_stream, " E[hits/match]: %.3g\n", score_set.ehits);
    fprintf(out_stream, " E[span]: %.3g\n", score_set.espan);
    */
  }

  fprintf(out_stream, "\n");
} // print_distr_parameters

/*************************************************************************
 * Print a program's parameters and CPU time.
 *************************************************************************/
void print_parameters(
  char* const argv[],				// command line
  int argc,				// number of fields in command line
  char* name,				// program name
  char* hmm_filename,			// hmm input file name
  char* seq_filename,			// sequence file name
  BOOLEAN_T viterbi,			// Viterbi search?
  double  dp_threshold,			// Repeated match threshold
  BOOLEAN_T motif_scoring,		// No partial motifs?
  BOOLEAN_T use_pvalues,		// Convert matches to p-values?
  double p_threshold,			// P-value threshold for motif hits.
  double e_threshold,			// E-value threshold for printing.
  BOOLEAN_T both_strands,		// Score both DNA strands?
  char* bg_filename,			// Background file name.
  char* sc_filename,			// Score file name.
  int pam_distance,			// PAM distance,
  double beta,				// Pseudo-weight.
  BOOLEAN_T zero_spacer_emit_lo,	// Zero emission for spacers?
  int max_gap,				// Maximum gap.
  double egcost,			// Cost factor for gaps.
  double gap_open,			// Gap open penalty.
  double gap_extend,			// Gap extend penalty.
  BOOLEAN_T print_fancy,		// Print alignments?
  BOOLEAN_T print_time,			// Print CPU time?
  double start_time,			// Starting time
  double end_time,			// Ending time
  int num_seqs,				// Number of sequences processed.
  EVD_SET evd_set, 			// EVD data structure.
  SCORE_SET score_set,   		// Scores.
  FILE* out_stream                      // Output stream.
)
{

  // Print parameters used in calculating p-values.
  print_distr_parameters(evd_set, score_set, out_stream);

  fprintf(out_stream, "Program parameters for %s\n", name);
  fprintf(out_stream, "  HMM file: ");
  if (hmm_filename == NULL) {
    fprintf(out_stream, "stdin\n");
  } else {
    fprintf(out_stream, "%s\n", hmm_filename);
  }
  fprintf(out_stream, "  Sequence file: ");
  if (seq_filename == NULL) {
    fprintf(out_stream, "stdin\n");
  } else {
    fprintf(out_stream, "%s\n", seq_filename);
  }
  if (viterbi) {
    fprintf(out_stream, "  Paths: single\n");
  } else {
    fprintf(out_stream, "  Paths: all\n");
  }
  if (dp_threshold != NO_REPEAT) {
    fprintf(out_stream, "  Minimum match score: %g\n", dp_threshold);
  }
  fprintf(out_stream, "  Partial motif matches are: %s\n", 
	 motif_scoring ? "not allowed" : "allowed");
  if (motif_scoring) fprintf(out_stream, "  Motif scores are: %s\n", 
	 use_pvalues ? "-log(p_value(log_odds)/p-thresh))" : "log_odds");
  if (use_pvalues)
    fprintf(out_stream, "  P-value threshold for motif hits: %g\n",
	    p_threshold);
  if (which_alphabet() == DNA_ALPH) {
    // FIXME: Add this if we ever implement scoring both strands.
    //fprintf(out_stream, "  Matches can be on: %s\n", 
    //	   both_strands ? "either strand" : "given strand only");
    if (both_strands) both_strands = TRUE;	// prevent compiler warning
  }
  fprintf(out_stream, "  E-value threshold for scores: %g\n", e_threshold);
  fprintf(out_stream, "  Background: %s\n", bg_filename ? bg_filename :
	  "nrdb");
  fprintf(out_stream, "  Score matrix for pseudocount frequencies: %s",
	 sc_filename ? sc_filename : "PAM");
  if (sc_filename) {
    fprintf(out_stream, "\n"); 
  } else { 
    fprintf(out_stream, " %d\n", pam_distance);
  }
  fprintf(out_stream, "  Beta (weight on pseudocount frequencies): %6.2f\n",
	  beta);
  fprintf(out_stream, "  Zero spacer emission log-odds: %s\n",
	 zero_spacer_emit_lo ? "true" : "false");
  if (max_gap > 0) {
    fprintf(out_stream, "  Maximum allowed gap: %d\n", max_gap);
  }
  if (egcost > 0) {
    fprintf(out_stream, "  Cost factor for gaps: %8.2f\n", egcost);
  }
  if (gap_open >= 0) {
    fprintf(out_stream, "  Gap open cost: %8.2g\n", gap_open);
  }
  if (gap_extend >= 0) {
    fprintf(out_stream, "  Gap extension cost: %8.2g\n", gap_extend);
  }
  fprintf(out_stream, "  Fancy output format: %s\n",
	  print_fancy ? "true" : "false");
  fprintf(out_stream, "\n");

  // Print the actual command line.
  {
    char *command = NULL;      // To hold the actual command line.
    int i = 1, pos = 0, len = 0;
    len += strlen(name)+1;                 // +1 for space following 
    mm_resize(command, len+2, char);      	// +1 for null
    strcpy(command+pos, name);
    command[len-1] = ' ';
    command[len] = '\0';
    pos = len;
    for (i=1; i<argc; i++) {
      len += strlen(argv[i])+1;                 // +1 for space following 
      mm_resize(command, len+2, char);      	// +1 for null
      strcpy(command+pos, argv[i]);
      command[len-1] = ' ';
      command[len] = '\0';
      pos = len;
    }
    fprintf(out_stream, "Actual command line:\n  %s\n\n", command);
    myfree(command);
  }

  /* Print the total CPU time and the CPU name. */
  if (print_time) {
    double total_time = myclock();
    fprintf(out_stream, "CPU: %s\n", hostname());
    fprintf(out_stream, "Total time %.2f secs.\n", (float)total_time/1E6);
    fprintf(out_stream, "Overhead %.2f secs.\n",
	   (float)(total_time - (end_time - start_time))/1E6);
    fprintf(out_stream, "%d sequences\n", num_seqs);
    fprintf(out_stream, "%.1f millisec/seq\n", 
	   (float)(end_time - start_time)/(num_seqs * 1E3));
    fprintf(out_stream, "%.1f microsec/character\n", 
	   (float)(end_time - start_time)/evd_set.total_length);
  }
} // print_parameters

/*************************************************************************
 * Given the desired output width, calculate the width of one line of
 * an alignment.
 *************************************************************************/
int compute_align_width
  (int output_width)
{
  int align_width;

  // Calculate the width of the alignment, minus margins.
  align_width = output_width - 2 * (LENGTH_WIDTH + 1);

  // Remove the 1's digit.
  align_width -= align_width % 10;

  return(align_width);
}

/*************************************************************************
 * Given a sequence, model and traceback, generate a sequence-to-model
 * match.
 *
 * Side effect: replaces the flanking Xs in the sequence with spaces.
 *************************************************************************/
void generate_match_sequence
  (MHMM_T*  the_log_hmm,    // Log-odds version of the HMM.
   int      trace_start,    // Position in sequence to start trace at.
   int      trace_end,	    // Position in sequence to end trace at.
   MATCH_T* traceback,      // The path through the model.
   MATRIX_T* motif_score_matrix, // Motif score matrix.
   PROB_T   p_threshold,    // P-value threshold for motif hits.
   SEQ_T*   sequence,       // The sequence as characters.
   char*    score_sequence, // Motif scores printed above the motif IDs.
   char*    id_sequence,    // Motif IDs printed above the sequence.
   char*    model_sequence, // Consensus chars along path.
   char*    match_sequence) // Indications of matches between sequence & model.
{
  int           i_trace;     // Index into the path.
  assert(trace_end < get_seq_length(sequence));

  // Find the consensus letter and match character at each position in path.
  for (i_trace = trace_start; i_trace <= trace_end; i_trace++) {
    MHMM_STATE_T* the_state;   // The current state in the model.
    int trace_entry = get_trace_item(i_trace, traceback);

    // Find the current state, if there is one.
    if (trace_entry != -1) {
      the_state = &(the_log_hmm->states[trace_entry]);
    } else {
      the_state = NULL;
    }

    // Mark unmatched section with spaces.
    if ((trace_entry == -1) ||
	(the_state->type == START_STATE) ||
	(the_state->type == END_STATE)) {
      score_sequence[i_trace] = ' ';
      id_sequence[i_trace] = ' ';
      model_sequence[i_trace] = ' ';
      match_sequence[i_trace] = ' ';

    } else {

      // Mark spacer states with dots in consensus, spaces elsewhere.
      if (the_state->type == SPACER_STATE) {
	score_sequence[i_trace] = ' ';
	id_sequence[i_trace] = ' ';
	model_sequence[i_trace] = '.';
	match_sequence[i_trace] = ' ';
      }

      // Mark motif sections ...
      else {

	// Find the motif score at this position.
        if (the_state->type == START_MOTIF_STATE && motif_score_matrix) {
          char tmp[400];
          int i;
          double s = get_matrix_cell(the_state->i_motif, i_trace, motif_score_matrix);
          double score = p_threshold > 0 ? pow(2.0, -s) * p_threshold : s;
	  sprintf(tmp, "%-*.1e/", the_state->w_motif, score);
          for (i=0; i<the_state->w_motif; i++) {
            score_sequence[i_trace+i] = tmp[i];
          }
        }

	// Find the ID character at this position.
	id_sequence[i_trace] = the_state->id_char;
	assert(the_state->id_char != '\0');

	// Find the consensus letter at this position.
	model_sequence[i_trace] = choose_consensus(TRUE, the_state);

	// Label exact matches with letters.
	if (model_sequence[i_trace] == get_seq_char(i_trace, sequence)) {
	  match_sequence[i_trace] = get_seq_char(i_trace, sequence);
	}

	// Label close matches with plus signs.
	else if (get_array_item(get_seq_int(i_trace, sequence),
				the_state->emit_odds) > 0.0) {
	  match_sequence[i_trace] = '+';
	} else {
	  match_sequence[i_trace] = ' ';
	}
      }
    }
  }

  // Add null terminators.
  score_sequence[trace_end+1] = '\0';
  id_sequence[trace_end+1] = '\0';
  model_sequence[trace_end+1] = '\0';
  match_sequence[trace_end+1] = '\0';

  // Replace flanking Xs with spaces.
  set_seq_char(0, ' ', sequence);
  set_seq_char(get_seq_length(sequence) - 1, ' ', sequence);
}

/*************************************************************************
 * A trivial function to check whether a given string consists of all
 * spaces.
 *************************************************************************/
static BOOLEAN_T all_spaces
  (char* a_string)
{
  int i;
  int length;

  length = strlen(a_string);
  for (i = 0; i < length; i++) {
    if (a_string[i] != ' ') {
      return(FALSE);
    }
  }
  return(TRUE);
}

/*************************************************************************
 * Figure out which positions to start and end the alignment.
 *
 * We are given the length of this sequence segment, its offset
 * relative to the entire sequence, the width of the desired
 * alignment, and the start and end positions of the match within this
 * segment.  The function converts the start and end positions from
 * segment-relative to sequence-relative coordinates.  Using these new
 * positions, the function computes the start position of the row of
 * the alignment that would contain the start of this match, and the
 * end position of the row of the alignment that would contain the end
 * of this match.  These coordinates are then converted back into
 * segment-relative coordinates.
 *************************************************************************/
void compute_alignment_start_end(
  int  align_width,
  int  seq_length,
  int  seq_offset,
  int  start_match,
  int  end_match,
  int* align_start,
  int* align_end
) {
  // Convert the start and end to absolute sequence coordinates.
  int absolute_start_match = start_match + seq_offset;
  int absolute_end_match = end_match + seq_offset;

  // The trace starts at first multiple of align_width prior
  // to the start of the match.
  int absolute_align_start 
    = ((absolute_start_match / align_width) * align_width);

  // Trace ends in the first whole multiple of align_width
  // containing the end of the match.
  int absolute_align_end
    = ((absolute_end_match / align_width) * align_width)
    + align_width;

  // The end cannot occur beyond the end of the sequence.
  absolute_align_end = MIN(seq_length + seq_offset - 2, absolute_align_end);

  // Revert to the relative coordinate for this segment.
  *align_start = MAX(1, absolute_align_start - seq_offset);
  *align_end = absolute_align_end - seq_offset;

  // Make sure everything is hunky dory.
  assert(start_match <= end_match);
  assert(*align_start < *align_end);
}

/*************************************************************************
 * Given a sequence, a model, and a match, print the sequence-to-model
 * alignment in BLAST-style format.
 *
 * The formatted output is returned as a single string, which must
 * be freed by the caller.
 *************************************************************************/
void print_alignment
  (int      align_width,   // Number of characters in one alignment line.
   MHMM_T*  the_log_hmm,   // Log-odds version of the HMM.
   MATCH_T* this_match,	   // The match.
   SEQ_T*   sequence,      // The sequence as characters.
   MATRIX_T* motif_score_matrix, // Motif score matrix.
   PROB_T    p_threshold,  // P-value threshold for motif hits.
   char**   outstring)     // Output string.
{
  char* score_sequence;	  // Motif scores.
  char* id_sequence;      // Motif ID numbers.
  char* model_sequence;   // Consensus chars along path.
  char* match_sequence;   // Indications of matches between
                          // sequence & model.
  int   i_seq;            // Index into the sequence.
  int   seq_length;       // Length of the sequence.
  char* raw_sequence;     // The sequence itself.
  int   end_position;     // End position of the current line.
  int   outstring_length; // Number of chars in output string.
  int   align_start;	  // Position in sequence to start alignment at.
  int   align_end;	  // Position in sequence to end alignment at.

  // Buffers for storing stuff before printing it.
  static char* buffer0 = NULL;
  static char* buffer1;
  static char* buffer2;	
  static char* buffer3;
  static char* buffer4;

  // Create local dynamic storage.
  if (buffer0 == NULL) {
    buffer0 = (char*)mm_malloc(sizeof(char) * (align_width + 1));
    buffer1 = (char*)mm_malloc(sizeof(char) * (align_width + 1));
    buffer2 = (char*)mm_malloc(sizeof(char) * (align_width + 1));
    buffer3 = (char*)mm_malloc(sizeof(char) * (align_width + 1));
    buffer4 = (char*)mm_malloc(sizeof(char) * (align_width + 1));
  }

  // Return the empty string if there's no trace.
  if (this_match == NULL) {
    *outstring = (char *)mm_malloc(sizeof(char));
    (*outstring)[0] = '\0';
    return;
  }

  // Set up the start and end points for the trace.
  compute_alignment_start_end(align_width,
			      get_seq_length(sequence),
			      get_seq_offset(sequence),
			      get_start_match(this_match),
			      get_end_match(this_match),
			      &align_start,
			      &align_end);

  // Allocate space for sequences.
  seq_length = get_seq_length(sequence);
  score_sequence = (char *)mm_malloc(sizeof(char) * (seq_length + 1));
  id_sequence = (char *)mm_malloc(sizeof(char) * (seq_length + 1));
  model_sequence = (char *)mm_malloc(sizeof(char) * (seq_length + 1));
  match_sequence = (char *)mm_malloc(sizeof(char) * (seq_length + 1));

  // Find the correspondence between sequence and model.
  generate_match_sequence(the_log_hmm,
			  align_start,
			  align_end,
			  this_match,
			  motif_score_matrix,
			  p_threshold,
			  sequence,
			  score_sequence,
			  id_sequence,
			  model_sequence,
			  match_sequence);

  // Allocate plenty of space for the output string.
  outstring_length = ((align_end - align_start) + align_width) * 10;
  *outstring = (char *)mm_malloc(sizeof(char) * outstring_length);

  // Display the alignment.
  sprintf(*outstring, "\n");
  buffer0[align_width] = '\0';
  buffer1[align_width] = '\0';
  buffer2[align_width] = '\0';
  buffer3[align_width] = '\0';
  buffer4[align_width] = '\0';

  raw_sequence = get_raw_sequence(sequence);
  for (i_seq = align_start; i_seq <= align_end; i_seq += align_width) {

    // Figure out whether we're nearing the end of the sequence.
    if (i_seq + align_width <= align_end) {
      strncpy(buffer0, &(score_sequence[i_seq]), align_width);
      strncpy(buffer1, &(id_sequence[i_seq]), align_width);
      strncpy(buffer2, &(model_sequence[i_seq]), align_width);
      strncpy(buffer3, &(match_sequence[i_seq]), align_width);
      strncpy(buffer4, &(raw_sequence[i_seq]), align_width);
      end_position = i_seq + align_width - 1;
    } else {
      strcpy(buffer0, &(score_sequence[i_seq]));
      strcpy(buffer1, &(id_sequence[i_seq]));
      strcpy(buffer2, &(model_sequence[i_seq]));
      strcpy(buffer3, &(match_sequence[i_seq]));
      strncpy(buffer4, &(raw_sequence[i_seq]), align_end - i_seq + 1);
      buffer4[align_end - i_seq + 1] = '\0'; // Null terminate the sequence.
      end_position = align_end;
    }

    // If there is no alignment in this part, then skip it.
    if ((all_spaces(buffer1)) && (all_spaces(buffer2)) && 
	(all_spaces(buffer3))) {
      continue;
    }

    if (motif_score_matrix) {		// scores above hits?
      sprintf(&((*outstring)[strlen(*outstring)]), "%*s %s %*s\n", 
	      LENGTH_WIDTH, "", buffer0, LENGTH_WIDTH, "");
    }
    sprintf(&((*outstring)[strlen(*outstring)]), "%*s %s %*s\n", 
	    LENGTH_WIDTH, "", buffer1, LENGTH_WIDTH, "");
    sprintf(&((*outstring)[strlen(*outstring)]), "%*s %s %*s\n",
	    LENGTH_WIDTH, "", buffer2, LENGTH_WIDTH, "");
    sprintf(&((*outstring)[strlen(*outstring)]), "%*s %s %*s\n",
	    LENGTH_WIDTH, "", buffer3, LENGTH_WIDTH, "");
    sprintf(&((*outstring)[strlen(*outstring)]), "%*d %-*s %*d\n",
	    LENGTH_WIDTH, MAX(1, i_seq) + get_seq_offset(sequence),
	    align_width, buffer4, 
	    LENGTH_WIDTH,
	    MIN(end_position, get_seq_length(sequence) - 2) + 
	    get_seq_offset(sequence));
    sprintf(&((*outstring)[strlen(*outstring)]), "\n");
  }
  sprintf(&((*outstring)[strlen(*outstring)]), "\n");

  // Check for array bounds write.
  assert((int)strlen(*outstring) < outstring_length);

  // Free locally allocated memory.
  myfree(score_sequence);
  myfree(id_sequence);
  myfree(model_sequence);
  myfree(match_sequence);
}

/*************************************************************************
 * Combine a bunch of information into a single GFF line.  The result
 * is returned in a local static string.
 *************************************************************************/
#define MAX_GFF 10000
static char* make_motif_gff_line 
  (char*  seq_id,
   char*  motif_id,
   int    start_match,
   int    width,
   double pvalue,
   char   strand)
{
  static char return_value[MAX_GFF];

  sprintf(return_value, "%s\tMeta-MEME\tmotif%s\t%d\t%d\t%g\t%c\t.\n",
	  seq_id,
	  motif_id,
	  start_match,
	  start_match + width - 1,
	  pvalue,
	  strand);

  return(return_value);
}

/*************************************************************************
 * Store a GFF entry for the given match object in the given string.
 * The gff_entry string is allocated by this function and must be
 * freed by the caller.
 *************************************************************************/
static void print_gff_entry
  (SEQ_T*   sequence,
   MATCH_T* this_match,
   MHMM_T*  the_log_hmm,
   MATRIX_T* motif_score_matrix, // Motif score matrix.
   char**   match_string,
   char**   motif_string)
{
  static char temp_match[MAX_GFF];  // Temporary storage.
  static char temp_motif[MAX_GFF];  // Temporary storage.
  int i_match;

  // Create a single line for the entire match, minus the evalue.
  sprintf(temp_match, "%s\t%s\t%s\t%d\t%d\t",
	  get_seq_name(sequence), 	  
	  "Meta-MEME",
	  "match",
	  // GFF indexes from 1...
	  get_start_match(this_match) + get_seq_offset(sequence) + 1,
	  get_end_match(this_match) + get_seq_offset(sequence) + 1);

  // Create additional lines for the motif matches.
  temp_motif[0] = '\0';
  for (i_match = get_start_match(this_match);
       i_match < get_end_match(this_match); i_match++) {
    int i_state = get_trace_item(i_match, this_match);
    MHMM_STATE_T* this_state = &((the_log_hmm->states)[i_state]);

    // Are we at the start of a motif?
    if (this_state->type == START_MOTIF_STATE) {
      strcat(temp_motif, 
	     make_motif_gff_line(get_seq_name(sequence),
				 get_state_motif_id(FALSE, this_state),
				 get_seq_offset(sequence) + i_match + 1,
				 this_state->w_motif,
				 motif_score_matrix ? 
				 get_matrix_cell(this_state->i_motif, 
						 i_match, 
						 motif_score_matrix) : 
				 0.0,
				 get_strand(this_state)
				 )
	     );
    }

  }

  // Allocate space and copy the gff entries.
  copy_string(match_string, temp_match);
  copy_string(motif_string, temp_motif);
}

/*************************************************************************
 * Store a sequence's scores, ID and description for later printing.
 * Also stores the score and length in score_set.
 *************************************************************************/
void store_sequence
  (BOOLEAN_T motif_scoring,  // Motif-scoring being used.
   BOOLEAN_T mhmmscan,       // Print mhmmscan output (else mhmms output)?
   int       output_width,   // Width in chars of the output.
   SEQ_T*    sequence,       // Info about the sequence. 
   MATCH_T*  this_match,     // Info about this match.
   PROB_T    e_threshold,    // E-value reporting threshold. 
   double    dp_threshold,   // Subtract from viterbi score.
   PROB_T    p_threshold,    // P-value threshold for motif hits.
   SCORE_SET *score_set,     // Scores and lengths for computing distribution.
   BOOLEAN_T got_evd,	     // EVD is available 
   EVD_SET   evd_set,	     // EVD
   BOOLEAN_T store_trace,    // Store the Viterbi trace?
   MHMM_T*   the_log_hmm,    // The HMM.
   MATRIX_T* motif_score_matrix, // Motif score matrix.
   BOOLEAN_T store_gff,      // Store a GFF entry?
   SCANNED_SEQUENCE_T *scanned_seq // CISML scanned sequence structure
)
{
  int    length;         // Length of this sequence.
  int	 span;		 // Span of match (end - start + 1).
  int	 nhits;		 // Number of motif hits in match.
  double viterbi_score;  // Score of this sequence.
  int    out_len = 0; 	 // Length printed so far.
  double gc = 0;	 // GC-content around match.
  double t_or_gc;	 // Length or GC-content of sequence.

  viterbi_score = get_score(this_match) - dp_threshold;
  length = get_match_seq_length(this_match);
  span = get_end_match(this_match) - get_start_match(this_match) + 1;
  nhits = get_nhits(this_match);

  //
  // Save the score and sequence length.  Keep track of maximum length
  // and smallest score.
  //
  if (viterbi_score > LOG_SMALL) {		/* don't save tiny scores */
    // Save unscaled score and sequence length for calculating EVD. 
    if (!(score_set->n % ALLOCATE_BLOCK)){
      mm_resize(score_set->scores, score_set->n + ALLOCATE_BLOCK, SCORE);
      //fprintf(stderr, "Resized score_set to %d \n", score_set->n + ALLOCATE_BLOCK);
    }
    score_set->scores[score_set->n].s = viterbi_score;	/* unscaled score */
			 get_alph_size(ALL_SIZE), 
    score_set->scores[score_set->n].t = length;
    score_set->scores[score_set->n].nhits = nhits;
    score_set->scores[score_set->n].span = span;
    if (length > score_set->max_t) score_set->max_t = length;
// FIXME:  Bill: do I need to worry about the sequence offset if
// the match doesn't contain the whole sequence?
    if (which_alphabet() == DNA_ALPH) {		// save GC content
      // the GC count is (gc[match_end] - gc[match_start-1])
      int start = MAX(0, get_start_match(this_match)-GC_WIN_SIZE);
      int end = MIN(get_seq_length(sequence)-1, get_end_match(this_match)+GC_WIN_SIZE);
      double n = end - start;
      gc = score_set->scores[score_set->n].gc = 
       ( get_seq_gc(end, sequence) - get_seq_gc(start, sequence) ) / n;
    }
    score_set->n++;
    //fprintf(stderr, "score_set->n %d %f\n", score_set->n, viterbi_score);
  } else {				/* skip tiny scores */
    return;
  }

  /* Use E-value to determine if score should be skipped. 
   * The threshold is multiplied by 10 because E-values may be fairly
   * inaccurate at this point.  */
  t_or_gc = (got_evd && evd_set.exponential) ? gc : length;
  if ((e_threshold < 0) || 
      ((got_evd) && 
       (Evd_set_evalue(score_set->n, viterbi_score, t_or_gc, 1, evd_set) > 
	10.0 * e_threshold))) {
    return;
  }
    
  /* Allocate memory, if necessary. */
  if (num_stored >= num_allocated) {
    if (num_allocated > 0) {
      num_allocated = 2 * num_allocated;
    }
    else {
      num_allocated = ALLOCATE_BLOCK;
    }
    stored_IDs =
      (char**)mm_realloc(stored_IDs, sizeof(char*) * num_allocated);
    stored_keys = 
      (PROB_T *)mm_realloc(stored_keys, sizeof(PROB_T) * num_allocated);
    stored_scores = 
      (PROB_T *)mm_realloc(stored_scores, sizeof(PROB_T) * num_allocated);
    stored_lengths = 
      (int *)mm_realloc(stored_lengths, sizeof(int) * num_allocated);
    stored_gcs = 
      (double *)mm_realloc(stored_gcs, sizeof(double) * num_allocated);
    stored_seqs = 
      (char**)mm_realloc(stored_seqs, sizeof(char*) * num_allocated);
    stored_traces = 
      (char**)mm_realloc(stored_traces, sizeof(char*) * num_allocated);
    stored_match_gffs = 
      (char**)mm_realloc(stored_match_gffs, sizeof(char*) * num_allocated);
    stored_motif_gffs = 
      (char**)mm_realloc(stored_motif_gffs, sizeof(char*) * num_allocated);
  }

  /* Check how long this ID is. */
  if ((signed)strlen(get_seq_name(sequence)) > longest_ID) {
    longest_ID = strlen(get_seq_name(sequence));
  }

  /* Allocate memory for the new sequence. */
  stored_IDs[num_stored] = (char*)mm_calloc(longest_ID + 1, sizeof(char));
  stored_seqs[num_stored] = (char*)mm_calloc(output_width + 1, sizeof(char));

  /* Store the ID, the key and the length. */
  strcpy(stored_IDs[num_stored], get_seq_name(sequence));
  stored_keys[num_stored] = viterbi_score;
  stored_lengths[num_stored] = length;
  stored_gcs[num_stored] = gc;
  
  // If scaned sequence is provided, add the match to it
  if (scanned_seq) {
    int start = MAX(0, get_start_match(this_match)-GC_WIN_SIZE);
    int end = MIN(get_seq_length(sequence)-1, get_end_match(this_match)+GC_WIN_SIZE);
    char *matched_sequence = get_raw_subsequence(start, end, sequence);
    MATCHED_ELEMENT_T *element = 
      allocate_matched_element(start, end, scanned_seq);
    set_matched_element_score(element, viterbi_score);
    set_matched_element_sequence(element, matched_sequence);
  }

  //
  // Store the scores and the description.
  //

  // Print score or "NaN" if too small.
  if (viterbi_score <= LOG_SMALL) {
    sprintf(stored_seqs[num_stored], "%*s ", SCORE_WIDTH, "NaN");
  } else {
    sprintf(stored_seqs[num_stored], 
	    "%*.*f ", SCORE_WIDTH, SCORE_PRECISION, 
	    viterbi_score + dp_threshold);
  }
  out_len += SCORE_WIDTH + 1;

  // Print GC content.
  if (0) {		// Disabled.
    sprintf(stored_seqs[num_stored] + out_len, 
	    "%*.*f ", SCORE_WIDTH, SCORE_PRECISION, 
	    score_set->scores[(score_set->n)-1].gc);
    out_len += SCORE_WIDTH + 1;
  }
 
  // Print number of hits and span if doing motif scoring.
  if (motif_scoring) {
    sprintf(stored_seqs[num_stored] + out_len,
	    "%*d %*d ", 
	    LENGTH_WIDTH, 
	    nhits,
	    LENGTH_WIDTH, 
            span
	    );
    out_len += 2 * (LENGTH_WIDTH + 1);
  } 

  // Print start and end of match.
  sprintf(stored_seqs[num_stored] + out_len,
	  "%*d %*d ",
	  LENGTH_WIDTH,
	  get_start_match(this_match) + get_seq_offset(sequence),
	  LENGTH_WIDTH,
	  get_end_match(this_match) + get_seq_offset(sequence)
	  );
  out_len += 2 * (LENGTH_WIDTH+1);

  // Print sequence length.
  if (!mhmmscan) {
    sprintf(stored_seqs[num_stored] + out_len,
	    "%*d ",
	    LENGTH_WIDTH,
	    get_seq_length(sequence) - 2	       // Chop off X's.
	    );
    out_len += LENGTH_WIDTH+1;
  }

  // Print description of sequence.
  if (output_width > out_len) {
    sprintf(stored_seqs[num_stored] + out_len,
	    "%-*.*s",
	    output_width - out_len,
	    output_width - out_len,
	    get_seq_description(sequence));
  }

  // Store the Viterbi alignment.
  if (store_trace) {
    print_alignment(compute_align_width(output_width),
		    the_log_hmm,
		    this_match,
		    sequence,
		    motif_score_matrix,
                    p_threshold,
		    &(stored_traces[num_stored]));
  }

  // Store the GFF entry.
  if (store_gff) {
    print_gff_entry(sequence, 
		    this_match,
		    the_log_hmm,
		    motif_score_matrix,
		    &(stored_match_gffs[num_stored]),
		    &(stored_motif_gffs[num_stored]));
  }

  /* Increment the number of stored sequences. */
  num_stored++;
} // store_sequence


#define SEPARATOR "-----------------------------------------------------------------------------------------------------------------------------------"

/* I define this as a struct only so that the two items can be passed as
   a single pointer to the 'qsort' routine. */
typedef struct tosort {
  PROB_T key;
  char*  out_string;
  char*  gff_string;
} TOSORT_T;


/*************************************************************************
 * A comparison function used in sorting the output before printing.
 * Sorts by keys in ascending order.
 *************************************************************************/
static int tosort_compare
  (const void * elem1,
   const void * elem2)
{
  const PROB_T key1 = ((TOSORT_T *)elem1)->key;
  const PROB_T key2 = ((TOSORT_T *)elem2)->key;
  const char*  string1 = ((TOSORT_T *)elem1)->out_string;
  const char*  string2 = ((TOSORT_T *)elem2)->out_string;

  if (key1 < key2) {
    return(-1);
  } else if (key1 > key2) {
    return(1);
  } else {
    return(strcmp(string2, string1));
  }
}

/*************************************************************************
 * Tiny function to find an EOL in a string.  Dies if none is found.
 *************************************************************************/
static int find_eol
  (char* this_entry)
{
  int i;
  int length;

  length = strlen(this_entry);
  for (i = 0; i < length; i++) {
    if (this_entry[i] == '\n') {
      return(i);
    }
  }
  die("No EOL found in %s.", this_entry);
  return(0);
}

/*************************************************************************
 * Sort the stored list of sequences, scores and IDs, and print them
 * out.
 *************************************************************************/
void sort_and_print_scores
  (BOOLEAN_T print_fancy,  // Include Viterbi alignment? 
   BOOLEAN_T print_header, // Print header? 
   BOOLEAN_T got_evd,      // EVD found? 
   BOOLEAN_T motif_scoring,// Motif-scoring being used.
   BOOLEAN_T mhmmscan,     // Print in mhmmscan format (else mhmms)?
   int       maxseqs,     // Maximum number of sequences to print
   int       output_width, // Width of output, in chars.
   PROB_T    threshold,    // Print sequences scoring below this. 
   BOOLEAN_T sort_output,  // Sort scores? 
   FILE*     gff_file,     // Auxiliary GFF file.
   FILE*     outfile)      // Print to this file. 
{
  int        i_stored;    // Index into the arrays of stored sequences. 
  TOSORT_T*  tosort;      // Array of output strings to be sorted. 
  int        num_tosort;  // Number of output strings to be sorted. 
  int        i_tosort;    // Index into array of to-be-sorted sequences. 
  int        i_toprint;   // Index for printing sequences. 
  int        num_toprint; // Number of sequences to print. 
  int        num_under = 0;/* Number of sequences stored with E-values
                               below given threshold. */

  /* If we get ridiculously long sequence IDs, truncate them. */
  if (longest_ID > MAX_SEQ_ID_LENGTH) {
    fprintf(stderr, "Warning: truncating absurdly long sequence IDs.\n");
    longest_ID = MAX_SEQ_ID_LENGTH;
  }

  // Get the accurate number of sequences below the E-value threshold. 
  // If there were too few scores for the E-value computation, print
  // the first "threshold" matches, sorted by decreasing score.
  if (got_evd) {
    num_under = 0;
    for (i_stored = 0; i_stored < num_stored; i_stored++) {
      double evalue = stored_keys[i_stored];
      if (evalue <= threshold) num_under++;
    }
  } else {
    // The sorted keys contain the raw scores.  Replace them
    // with -score so that the sort by increasing order will
    // put the best scores first.
    fprintf(stderr, "Warning: Unable to compute E-values.  ");
    fprintf(stderr, "Sorting by score instead.\n");
    for (i_stored = 0; i_stored < num_stored; i_stored++) {
      stored_keys[i_stored] *= -1;	// -score so sort works.
    }
    num_under = num_stored;
    threshold = HUGE_VAL;		// All keys less than this.
  }

  tosort = (TOSORT_T *)mm_malloc(sizeof(TOSORT_T) * num_stored);

  // Consider each stored sequence in turn. 
  i_tosort = 0;
  for (i_stored = 0; i_stored < num_stored; i_stored++) {
    double score = stored_keys[i_stored];  // E-value or raw score.
    double evalue;                         // E-value or NaN.

    // Skip low-scoring sequences 
    if (score > threshold) continue;

    // Allocate memory for this item. 
    if (print_fancy) {
      tosort[i_tosort].out_string =
        (char *)mm_malloc(sizeof(char) * (output_width + 2 +
       strlen(stored_traces[i_stored])));
    } else {
      tosort[i_tosort].out_string =
         (char *)mm_malloc(sizeof(char) * (output_width + 2));
    }
    if (gff_file != NULL) {
      tosort[i_tosort].gff_string =
        (char *)mm_malloc(sizeof(char) * (strlen(stored_match_gffs[i_stored])
          + strlen(stored_motif_gffs[i_stored])
          + 50)); // Space for the evalue, etc.
    }

    // Format the ID, scores and description. 
    evalue = got_evd ? score : NaN();
    sprintf(
      tosort[i_tosort].out_string,
      "%-*.*s %*.2e %*.*s\n",
      longest_ID, 
      longest_ID, 
      stored_IDs[i_stored], 
      PV_WIDTH,
      evalue,
      output_width - (longest_ID + PV_WIDTH + 2),
      output_width - (longest_ID + PV_WIDTH + 2),
      stored_seqs[i_stored]
    );

    // Add the trace, if requested. 
    if (print_fancy) {
      strcat(tosort[i_tosort].out_string, stored_traces[i_stored]);
    }

    // Add the GFF entry, if requested.
    if (gff_file != NULL) {
      sprintf(
        tosort[i_tosort].gff_string, "%s%g\t.\t.\n%s",
        stored_match_gffs[i_stored], 
        fabs(score),
        stored_motif_gffs[i_stored]
      );
    }

    // Store the associated key. 
    tosort[i_tosort].key = score;
    i_tosort++;
  }
  num_tosort = i_tosort;

  // Sort the sequences. 
  if (sort_output) {
    qsort((void *)tosort, num_tosort, sizeof(TOSORT_T), tosort_compare);
  }

  // Decide how many sequences to print. 
  num_toprint = got_evd ? num_under : num_stored;
  // If the maxseqs limit has been set, apply it
  if (maxseqs > NO_MAX_SEQS) {
    num_toprint = MIN(num_toprint, maxseqs);
  }
    
  // Print the header. 
  if (print_header) {
    fprintf(outfile, "%-*s ", longest_ID, "ID");
    fprintf(outfile, "%*.*s ", PV_WIDTH, PV_WIDTH, "E-value");
    fprintf(outfile, "%*s ", SCORE_WIDTH, "Score");
    if (motif_scoring) {
      fprintf(outfile, "%*s ", LENGTH_WIDTH, "Hits");
      fprintf(outfile, "%*s ", LENGTH_WIDTH, "Span");
    }
    fprintf(outfile, "%*s ", LENGTH_WIDTH, "Start");
    fprintf(outfile, "%*s ", LENGTH_WIDTH, "End");
    if (!mhmmscan) {
      fprintf(outfile, "%*s ", LENGTH_WIDTH, "Length");
    }
    fprintf(outfile, "Description\n");
    fprintf(outfile, "%*.*s\n", output_width, output_width, SEPARATOR);
  }

  // Print the sequences (only the first line).
  for (i_toprint = 0; i_toprint < num_toprint; i_toprint++) {
    char* this_entry = tosort[i_toprint].out_string;
    int end_of_line = find_eol(this_entry);
    this_entry[end_of_line] = '\0';
    fprintf(outfile, "%s\n", this_entry);
    this_entry[end_of_line] = '\n';
  }

  // Print the traces (if requested).
  if (print_fancy) {
    fprintf(outfile, "%*.*s\n\n", output_width, output_width, SEPARATOR);
    for (i_toprint = 0; i_toprint < num_toprint; i_toprint++) {
      fprintf(outfile, "%s", tosort[i_toprint].out_string);
    }
  }
  if (print_header) {
    fprintf(outfile, "%*.*s\n\n", output_width, output_width, SEPARATOR);
  }

  // Print the GFFs.
  if (gff_file != NULL) {
    for (i_toprint = 0; i_toprint < num_toprint; i_toprint++) {
      fprintf(gff_file, "%s", tosort[i_toprint].gff_string);
    }
  }

  // Free up memory for the sorted sequences.
  for (i_tosort = 0; i_tosort < num_tosort; i_tosort++) {
    myfree(tosort[i_tosort].out_string);
    if (gff_file != NULL) {
      myfree(tosort[i_tosort].gff_string);
    }
  }
  myfree(tosort);  

  // Free up memory for the stored sequences.
  for (i_stored = 0; i_stored < num_stored; i_stored++) {
    myfree(stored_seqs[i_stored]);
    myfree(stored_IDs[i_stored]);
    if (print_fancy) {
      myfree(stored_traces[i_stored]);
    }
    if (gff_file != NULL) {
      myfree(stored_match_gffs[i_stored]);
      myfree(stored_motif_gffs[i_stored]);
    }      
  }
  myfree(stored_seqs);
  myfree(stored_IDs);
  myfree(stored_traces);
  myfree(stored_match_gffs);
  myfree(stored_motif_gffs);
  myfree(stored_keys);
  myfree(stored_scores);
  myfree(stored_lengths);
  myfree(stored_gcs);
}

/***********************************************
 * Calculate the score distribution.
 ***********************************************/
EVD_SET calc_distr(
  SCORE_SET score_set,			// Set of scores and lengths.
  BOOLEAN_T exponential,		// Use exponential distribution?
  BOOLEAN_T match_evalues		// Use match E-values?
)
{
  double H = 0;				       // don't use H 
  int maxiter1 = 10;                           // max. iterations in L-H loop
  int maxiter2 = 20;                           // max. iterations in N-R loop
  int size = 1000;    	                       // size of length groups 
  double lspan = 1.5;                          // min. length ratio/group 

  //
  // If using match E-values instead of score E-values, set all
  // the sequence lengths to total_length/# matches.
  //
  if (!exponential && match_evalues) {
    int i;
    int length = score_set.n ? score_set.total_length/score_set.n : score_set.total_length;
    for (i=0; i<score_set.n; i++) score_set.scores[i].t = length;
    score_set.max_t = length;
  }

  // Fit the distribution to the scores and lengths.
  return(fit_score_distribution(exponential, score_set, H, maxiter1, maxiter2, EPS1, EPS2, 
      size, lspan));

} // calc_distr 

/*************************************************************************
* Calculate the E-values of sequence scores and store them as the keys.
*************************************************************************/
void calc_evalues (
  EVD_SET *evd_set, 	// EVD parameters.
  int n, 		// Number of sequences or matches.
  int q,		// Length of query.
  int t		 	// Use as length of targets if non-zero.
)    
{
  int i;
  double t_or_gc;

  if (evd_set->n <= 0) return;			// no EVD available 

  // Get min(E-value) and sum of log(E).
  evd_set->outliers = 0;			// Number with E < 1.
  evd_set->min_e = BIG;				// Smallest E-value.
  evd_set->sum_log_e = 0;			// Sum of log(E < 1).

  for (i=0; i<num_stored; i++) {
    double evalue;
    stored_scores[i] = stored_keys[i];
    // E-value is new key.
    t_or_gc = evd_set->exponential ? stored_gcs[i] : (t ? t : stored_lengths[i]); 
    evalue = stored_keys[i] = Evd_set_evalue(n, 
				    stored_keys[i], 
      				    t_or_gc,
      				    q, 
				    *evd_set);
    if (evalue < evd_set->min_e) evd_set->min_e = evalue;
    if (evalue < 1) {
      evd_set->outliers++;
      evd_set->sum_log_e += log(evalue);
    }
  }
} // calc_evalues 

#ifdef MAIN
#include "simple-getopt.h"

VERBOSE_T verbosity = NORMAL_VERBOSE;

/**************************************************************************
 * int main
 **************************************************************************/
int main (int argc, char *argv[])
{
  // Command line parameters. 
  char*     hmm_filename;     // File containing the HMM. 
  FILE*     hmm_file;
  char*     seq_filename;     // File containing the sequences. 
  FILE*     seq_file;
  char*     viterbi_string;   // Command line buffer.
  BOOLEAN_T viterbi;          // Use Viterbi scoring.
  BOOLEAN_T local_scoring;    // Use local scoring.
  BOOLEAN_T motif_scoring;    // Perform motif-scoring.
  BOOLEAN_T use_pvalues;      // Convert motif scores to p-values.
  PROB_T    p_threshold;      // P-value threshold for motif hits.
  int       maxseqs;          // Maximum number of sequences to print
  BOOLEAN_T both_strands;     // Score both DNA strands.
  PROB_T    e_threshold;      // E-value threshold for scores.
  char*     bg_filename;      // File containing background frequencies. 
  int       pam_distance;     // PAM distance 
  char*     sc_filename;      // File containing substitution scores. 
  double    beta;	      // Weight on pseudocounts.
  BOOLEAN_T allow_weak_motifs;  // Allow motifs with min p-value > p-thresh?
  BOOLEAN_T zero_spacer_emit_lo; // Set spacer emission log-odds = 0?
  double    gap_open;	      // Cost to open a gap; ignore if < 0.
  double    gap_extend;	      // Cost to extend a gap; ignore if < 0.
  int       output_width;     // Width of output, in chars.
  BOOLEAN_T print_fancy;      // Print Viterbi alignments? 
  BOOLEAN_T sort_output;      // Sort output scores? 
  int       progress_every;   // Show progress after every n iterations.
  BOOLEAN_T print_header;     // Print header information? 
  BOOLEAN_T print_params;     // Print program parameters? 
  BOOLEAN_T print_time;       // Print timing info? 

  // Data structures. 
  MHMM_T*   the_hmm;          // The HMM itself. 
  MHMM_T*   the_log_hmm;      // The HMM, with probs converted to logs. 
  SEQ_T*    sequence;         // Sequence to search against. 
  MATRIX_T* motif_score_matrix; // Number of motifs x sequence length.
  MATRIX_T* dp_matrix = NULL; // Dynamic programming matrix. 
  MATRIX_T* trace_matrix = NULL; // Traceback for Viterbi.
  MATCH_T*  this_match = allocate_match(); // This sequence-to-model match.
  int       dp_rows = 0;      // Size of the DP matrix. 
  int       dp_cols = 0;
  EVD_SET   evd_set;          // EVD data structure 
  BOOLEAN_T got_evd = FALSE;  // no EVD found yet 
  SCORE_SET *score_set = NULL;// Set of scores/lengths for computing EVD.

  // Local variables. 
  double    start_time;    // Time at start of sequence processing. 
  double    end_time;      // Time at end of sequence processing. 
  int       num_seqs;      // Number of sequences processed. 

  // Record CPU time. 
  myclock();

  /***********************************************
   * Parse the command line.
   ***********************************************/

  // Set defaults. 
  hmm_filename = NULL;
  seq_filename = NULL;
  viterbi_string = NULL;
  viterbi = TRUE;
  local_scoring = TRUE;
  motif_scoring = FALSE;
  use_pvalues = FALSE;
  p_threshold = -1;	// Don't do p-value scoring.
  maxseqs = NO_MAX_SEQS;
  both_strands = FALSE; // Score given DNA strand only.
  e_threshold = DEFAULT_E_THRESHOLD;
  bg_filename = NULL;
  pam_distance = -1;    // Use default PAM distance.
  sc_filename = NULL;
  beta = -1;		// Illegal beta, use defaults.
  allow_weak_motifs = FALSE;    // Don't allow weak motifs.
  zero_spacer_emit_lo = FALSE;
  gap_open = -1;        // No gap open penalty.
  gap_extend = -1;      // No gap extension penalty.
  output_width = DEFAULT_OUTPUT_WIDTH;
  print_fancy = FALSE;
  sort_output = TRUE;
  progress_every = DEFAULT_PROGRESS_EVERY;
  print_header = TRUE;
  print_params = TRUE;
  print_time = TRUE;

  {

    // Define the usage message.
    char      usage[1000] = "";
    // Define command line options.
    cmdoption const options[] = {
			{"paths", REQUIRED_VALUE},
			{"global", NO_VALUE},
      {"maxseqs", REQUIRED_VALUE},
			{"pthresh", REQUIRED_VALUE},
			{"both-strands", NO_VALUE},
			{"ethresh", REQUIRED_VALUE},
			{"fancy", NO_VALUE},
			{"width", REQUIRED_VALUE},
			{"nosort", NO_VALUE},
			{"bg-file", REQUIRED_VALUE},
			{"allow-weak-motifs", NO_VALUE},
			{"progress", REQUIRED_VALUE},
			{"verbosity", REQUIRED_VALUE},
			{"noheader", NO_VALUE},
			{"noparams", NO_VALUE},
			{"notime", NO_VALUE},
			{"quiet", NO_VALUE},
			{"zselo", NO_VALUE},
			{"gap-open", REQUIRED_VALUE},
			{"gap-extend", REQUIRED_VALUE},
			{"motif-scoring", NO_VALUE},
			{"pseudo-weight", REQUIRED_VALUE},
			{"pam", REQUIRED_VALUE},
			{"score-file", REQUIRED_VALUE}
    };

    int option_count = 24;
    int option_index = 0;

    strcat(usage, "USAGE: mhmms [options] <HMM file> <FASTA file>\n");
    strcat(usage, "\n");
    strcat(usage, "   Options:\n");
    strcat(usage, "     --paths [single|all] (default=single)\n");
    strcat(usage, "     --global (default=local)\n");
    strcat(usage, "     --maxseqs <int>\n");
    strcat(usage, "     --pthresh <p-value>\n");
    strcat(usage, "     --both-strands\n");
    strcat(usage, "     --ethresh <E-value>\n");
    strcat(usage, "     --fancy\n");
    strcat(usage, "     --width <int> (default=79)\n");
    strcat(usage, "     --nosort\n");
    strcat(usage, "     --bg-file <file>\n");
    strcat(usage, "     --allow-weak-motifs\n");
    strcat(usage, "     --progress <int>\n");
    strcat(usage, "     --verbosity 1|2|3|4|5 (default=2)\n");
    strcat(usage, "     --noheader\n");
    strcat(usage, "     --noparams\n");
    strcat(usage, "     --notime\n");
    strcat(usage, "     --quiet\n");
    strcat(usage, "\n");
    strcat(usage, "   Advanced options:\n");
    strcat(usage, "     --zselo\n");
    strcat(usage, "     --gap-open <cost>\n");
    strcat(usage, "     --gap-extend <cost>\n");
    strcat(usage, "     --motif-scoring\n");
    strcat(usage, "     --pseudo-weight <weight> (default=10)\n");
    strcat(usage, "     --pam <distance> (default=250 [protein] 1 [DNA])\n");
    strcat(usage, "     --score-file <file>\n");
    strcat(usage, "\n");

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

      if (strcmp(option_name, "paths") == 0) {
        viterbi_string = option_value;
      } else if (strcmp(option_name, "global") == 0) {
        local_scoring = FALSE;
      } else if (strcmp(option_name, "maxseqs") == 0) {
        maxseqs = atoi(option_value);
      } else if (strcmp(option_name, "pthresh") == 0) {
        p_threshold = atof(option_value);
      } else if (strcmp(option_name, "both-strands") == 0) {
        both_strands = TRUE;
      } else if (strcmp(option_name, "ethresh") == 0) {
        e_threshold = atof(option_value);
      } else if (strcmp(option_name, "fancy") == 0) {
        print_fancy = TRUE;
      } else if (strcmp(option_name, "width") == 0) {
        output_width = atoi(option_value);
      } else if (strcmp(option_name, "nosort") == 0) {
        sort_output = FALSE;
      } else if (strcmp(option_name, "bg-file") == 0) {
        bg_filename = option_value;
      } else if (strcmp(option_name, "allow-weak-motifs") == 0) {
        allow_weak_motifs = TRUE;
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
      } else if (strcmp(option_name, "motif-scoring") == 0) {
        motif_scoring = TRUE;
      } else if (strcmp(option_name, "pseudo-weight") == 0) {
        beta = atof(option_value);
      } else if (strcmp(option_name, "pam") == 0) {
        pam_distance = atof(option_value);
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

  // Make sure we got the required files. 
  if (hmm_filename == NULL) 
    die("No HMM file given.\n");
  if (seq_filename == NULL)
    die("No sequence file given.\n");

  // Figure out what kind of scoring to do.
  if ((viterbi_string == NULL) || (strcmp(viterbi_string, "single") == 0)) {
    viterbi = TRUE;
  } else if (strcmp(viterbi_string, "all") == 0) {
    viterbi = FALSE;
  } else {
    die("Illegal option (\"-paths %s\").", viterbi_string);
  }

  // Force p-value scoring if p-value threshold given.
  if (p_threshold != -1) use_pvalues = TRUE;

  // Check p-threshold is in range [0<p<=1].
  if (use_pvalues && (p_threshold <= 0 || p_threshold > 1))
    die("You must specify p-thresh in the range [0<p<=1]\n");

  // Force motif_scoring if using p-value scoring.
  if (use_pvalues) motif_scoring = TRUE;

  // FIXME: both-strands not implemented
  if (both_strands) {
    die("Sorry, -both-strands not yet implemented.");
  }

  // Force motif_scoring if using both strands.
  if (both_strands) motif_scoring = TRUE;

  // We can't do fancy output with total probability scoring yet.
  if ((print_fancy) && (!viterbi)) {
    die("Sorry: mhmms cannot yet produce alignments from all paths scoring.");
  }

  /***********************************************
   * Set up the model.
   ***********************************************/

  // Read the model. 
  read_mhmm(hmm_filename, &the_hmm);

  // Set the PAM distance.
  if (pam_distance == -1) {
    pam_distance = (which_alphabet() == PROTEIN_ALPH) ? DEFAULT_PROTEIN_PAM
      : DEFAULT_DNA_PAM;
  }
  if (beta < 0) beta = (which_alphabet() == PROTEIN_ALPH) ? 10 : 1;

  // Read the background distribution.
  free_array(the_hmm->background);
  the_hmm->background = get_background(bg_filename);

  // Convert the model to log space.
  convert_to_from_log_hmm(TRUE, // Convert to logs.
			  zero_spacer_emit_lo,
			  gap_open,
			  gap_extend,
			  the_hmm->background,
			  sc_filename,
			  pam_distance,
			  beta,
			  the_hmm, &the_log_hmm);

  // Set up PSSM matrices if doing motif_scoring
  // and pre-compute motif p-values if using p-values.
  // Set up the hot_states list.
  set_up_pssms_and_pvalues(motif_scoring,
                           p_threshold,
                           both_strands,
			   allow_weak_motifs,
                           the_log_hmm);

  //
  // Set up for computing score distribution.
  //
  score_set = set_up_score_set(
    			p_threshold, 	
    			-1, // dp_threshold, 
    			-1, // max_gap
			FALSE, // negatives_only
    			the_log_hmm);

  /***********************************************
   * Search the database one sequence at a time.
  ***********************************************/

  // Open the file for reading. 
  if (open_file(seq_filename, "r", TRUE, "sequence", "sequences", &seq_file) 
      == 0)
    exit(1);

  start_time = myclock();

  num_seqs = 0;
  while (read_one_fasta(seq_file, MAX_SEQ, &sequence)) {

    num_seqs++;

    // Keep track of total database size for E-value calculation.
    score_set->total_length += get_seq_length(sequence);

    // Let the user know what's going on.
    if (verbosity > NORMAL_VERBOSE) {
      fprintf(stderr, "Scoring %s (length=%d).\n", get_seq_name(sequence),
	      get_seq_length(sequence));
    }

    // Convert the sequence to alphabet-specific indices. 
    prepare_sequence(sequence);
    assert(get_seq_char(get_seq_length(sequence), sequence) == '\0');

    /* Allocate the dynamic programming matrix. Rows correspond to
       states in the model, columns to positions in the sequence. */
    if ((dp_rows < the_log_hmm->num_states) 
	|| (dp_cols < get_seq_length(sequence))) {
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
      motif_score_matrix = allocate_matrix(the_log_hmm->num_motifs,
        				   get_seq_length(sequence));
      compute_motif_score_matrix(use_pvalues,
				 p_threshold,
				 get_int_sequence(sequence),
				 get_seq_length(sequence),
				 the_log_hmm,
				 &motif_score_matrix);
    } else {
      //FIXME: remove 
      motif_score_matrix = NULL;
    }

    // Perform the appropriate type of dynamic programming.
    if (viterbi) {
      viterbi_algorithm(
			local_scoring,
			get_int_sequence(sequence),
			get_seq_length(sequence),
			the_log_hmm,
			TRUE,  // Save Viterbi path?
			motif_score_matrix,
			dp_matrix,
			trace_matrix,
			this_match
			);
    } else {
      forward_algorithm(
			local_scoring,
			get_int_sequence(sequence),
			get_seq_length(sequence),
			the_log_hmm,
			//motif_score_matrix, // FIXME
			dp_matrix,
			this_match
			);
    }      

    // Store the score, ID, length and comment for later printing. 
    //assert(get_match_seq_length(this_match) == get_seq_length(sequence));
    //assert(get_match_seq_length(this_match) == get_seq_length(sequence));
    store_sequence(motif_scoring,
		   FALSE, // Don't print in mhmmscan format.
		   output_width, 
		   sequence, 
		   this_match, 
		   e_threshold,
		   0,		// no dp_threshold
		   p_threshold,
		   score_set,
		   got_evd, 
		   evd_set,
		   print_fancy,
		   the_log_hmm,
		   motif_score_matrix,
		   FALSE, // Store GFF?
       NULL
		   );

    /* Calculate the initial E-value distribution if the required
     * number of sequences has been saved.  This will allow the
     * descriptions of low-scoring sequences to not be stored.  The
     * distribution will be recomputed using all saved scores when all
     * sequences have been read.  */
    if (score_set->n == EVD_NUM_SEQS && got_evd == FALSE) {
      evd_set = calc_distr(
        *score_set,			// Set of scores.
        use_pvalues, 			// Use exponential distribution?
        FALSE				// Use sequence E-values.
      );
      if (evd_set.n > 0) got_evd = TRUE;
    }

    // Free the memory used by this sequence. 
    free_matrix(motif_score_matrix);
    free_seq(sequence);
    if ((verbosity >= NORMAL_VERBOSE) && (num_seqs % progress_every == 0)) {
      fprintf(stderr, "\rSequences: %d", num_seqs);
    }
  }
  end_time = myclock();

  /***********************************************
   * Calculate the E-values and store them as the keys.
   ***********************************************/
  // Recalculate the EVD using all scores. 
  // If successful, calculate the E-values and
  // store them as keys.
  evd_set = calc_distr(
    *score_set,     			// Set of scores.
    use_pvalues, 			// Use exponential distribution?
    FALSE				// Use sequence E-values.
  );
  if (evd_set.n > 0) {
    int q, t, N;
    q = 1;				// Ignore query "length".
    t = 0;				// Use stored target lengths.
    N = evd_set.non_outliers;		// Use number of non-outliers.
    evd_set.N = N;
    calc_evalues(&evd_set, N, q, t);
    got_evd = TRUE;
  }
  evd_set.negatives_only = FALSE;	// Used real sequences.

  /***********************************************
   * Print header information.
   ***********************************************/
  if (print_header) {
    write_header(
		 "MIAO",
		 "Database search results",
		 the_hmm->description,
		 the_hmm->motif_file,
		 hmm_filename,
		 seq_filename,
		 stdout
		 );
  }

  /***********************************************
   * Sort and print the results.
   ***********************************************/
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "\nSorting the E-values.\n");
  }
  sort_and_print_scores(
    print_fancy, 
    print_header, 
    got_evd,
    motif_scoring,
    FALSE, // Don't print in mhmmscan format.
    maxseqs, // Maximum number of sequences to print
    output_width,
    e_threshold,
    sort_output, 
    NULL, // FIXME: GFF output?
    stdout
  );

  /***********************************************
   * Print the program parameters.
   ***********************************************/
  if (print_params) {
    print_parameters(
     argv,
     argc,
     "mhmms",
     hmm_filename,
     seq_filename,
     viterbi,
     NO_REPEAT,   // Repeated match threshold.
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
     -1,    // No max_gap.
     0,      // No egcost.
     gap_open,
     gap_extend,
     print_fancy,
     print_time,
     start_time,
     end_time,
     num_seqs,
     evd_set,
     *score_set,
     stdout
    );
  }

  /* Tie up loose ends. */
  // myfree(evd_set.evds);
  free_mhmm(the_hmm);
  free_mhmm(the_log_hmm);
  free_matrix(dp_matrix);
  free_matrix(trace_matrix);
  free_match(this_match);
  fclose(seq_file);

  return(0);
}
#endif /* MAIN */

/*
 * Local Variables:
 * mode: c++
 * c-basic-offset: 2
 * End:
 */
