/**************************************************************************
 * FILE: mhmms.h
 * AUTHOR: William Stafford Noble and Timothy L. Bailey
 * CREATE DATE: 8-13-97
 * PROJECT: MHMM
 * VERSION: $Revision 1.10$
 * DESCRIPTION: Search a database of sequences using a motif-based HMM.
 **************************************************************************/
#ifndef MHMMS_H
#define MHMMS_H

#include "cisml.h"
#include "match.h"
#include "seq.h"

#define DEFAULT_OUTPUT_WIDTH 132  // Number of chars in output lines.
#define DEFAULT_E_THRESHOLD 10.0  // Default E-value threshold.
#define DEFAULT_DP_THRESHOLD 10.0 // Threshold for repeated matches.

#define EVD_NUM_SEQS 10000        // Number of sequences required for EVD.
#define NO_MAX_SEQS -1     // Value indicating maxseqs not set

// Print progress report every N iterations.
#define DEFAULT_PROGRESS_EVERY 1000

/*************************************************************************
 * Return the value of num_stored.
 *************************************************************************/
int get_num_stored();

/*************************************************************************
 * Print a program's parameters and CPU time.
 *************************************************************************/
void print_parameters(
  char* const argv[],                   // command line
  int argc,                             // number of fields in command line
  char* name,                           // program name
  char* hmm_filename,                   // hmm input file name
  char* seq_filename,                   // sequence file name
  BOOLEAN_T viterbi,                    // Viterbi search?
  double dp_threshold,                  // Repeated match threshold
  BOOLEAN_T motif_scoring,              // No partial motifs?
  BOOLEAN_T use_pvalues,                // Convert matches to p-values?
  double p_threshold,                   // P-value threshold for motif hits.
  double e_threshold,                   // E-value threshold for scores.
  BOOLEAN_T both_strands,               // Score both DNA strands?
  char* bg_filename,                    // Background file name.
  char* sc_filename,                    // Score file name.
  int pam_distance,                     // PAM distance,
  double beta,                          // Pseudo-weight.
  BOOLEAN_T zero_spacer_emit_lo,        // Zero emission for spacers?
  int max_gap,				// Maximum gap.
  double egcost,			// Cost factor for gaps.
  double gap_open,                      // Gap open penalty.
  double gap_extend,                    // Gap extend penalty.
  BOOLEAN_T print_fancy,                // Print alignments?
  BOOLEAN_T print_time,                 // Print CPU time?
  double start_time,                    // Starting time
  double end_time,                      // Ending time
  int num_seqs,				// Number of sequences processed.
  EVD_SET evd_set,			// EVD data structure.
  SCORE_SET score_set,			// Scores.
  FILE* out_stream                      // Output stream.
);

/*************************************************************************
 * Given the desired output width, calculate the width of one line of
 * an alignment.
 *************************************************************************/
int compute_align_width
  (int output_width);

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
);

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
   char*    match_sequence);// Indications of matches between sequence & model.

/*************************************************************************
 * Store a sequence's scores, ID and description for later printing.
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
   BOOLEAN_T store_gff,     // Store a GFF entry?
   SCANNED_SEQUENCE_T *scanned_seq // CISML scanned sequence structure
);

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
   char**   outstring);    // Output string.

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
   FILE*     outfile);     // Print to this file. 

/*************************************************************************
 * Calculate the score distribution.
 *************************************************************************/
EVD_SET calc_distr(
  SCORE_SET score_set,     		// Set of scores and lengths.
  BOOLEAN_T exponential,		// Use exponential distribution?
  BOOLEAN_T match_evalues		// Use match E-values?
);

/*************************************************************************
 * Calculate the E-values of sequence scores and store them as the keys.
 *************************************************************************/
void calc_evalues (
  EVD_SET *evd_set,	// EVD parameters.
  int n,         	// Number of sequences or matches.
  int q,                // Length of query.
  int t                 // Use as length of targets if non-zero.
);

#endif
