/**************************************************************************
 * FILE: mhmm.h
 * AUTHOR: William Stafford Noble, Timothy L. Bailey, and Charles E. Grant
 * CREATE DATE: 07-27-2007
 * PROJECT: MHMM
 * DESCRIPTION: Create a HMM from MEME motif and motif occurrence information.
 **************************************************************************/
#ifndef MHMM_H
#define MHMM_H

/***********************************************************************
 *  Select the motifs used to build the model, parse any motif
 *  occurences, build the motif order object, and the motif
 *  and spacer frequency matrices.
 ***********************************************************************/
void process_raw_motifs_for_model(
     int* num_motifs,                  // Number of motifs. IN, OUT
     MOTIF_T* motifs,                  // Array of motifs IN, OUT
     STRING_LIST_T* motif_occurrences, // List of motif occurrences. OUT
     STRING_LIST_T* requested_motifs,  // Explicitly requested motifs. IN
     BOOLEAN_T has_reverse_strand,     // Did file contain both strands? IN
     BOOLEAN_T keep_unused,            // Retain unsed motifs? IN
     double p_threshold,               // Motif p-value threshold IN
     double e_threshold,               // Motif e-value threshold IN
     double complexity_threshold,      // Motif complexity threshold IN
     ORDER_T** order_spacing,          // Motif/spacer order IN, OUT
     MATRIX_T** transp_freq,           // Motif transition freqs OUT
     MATRIX_T** spacer_ave,            // Spacer transition freqs OUT
     double trans_pseudo,              // Motif transition pseudo-counts IN
     double spacer_pseudo              // Spacer transition pseudo-counts IN
);

#endif
