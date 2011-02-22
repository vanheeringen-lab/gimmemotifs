/**************************************************************************
 * FILE: meme-io.h
 * CREATE DATE: 3/5/2001
 * AUTHOR: William Stafford Noble
 * PROJECT: MHMM
 * COPYRIGHT: 2001-2008, WSN
 * DESCRIPTION: Read a collection of motifs from a MEME output file.
 **************************************************************************/
#ifndef MEME_IO_H
#define MEME_IO_H

#include "array-list.h"
#include "matrix.h"
#include "array.h"
#include "motif.h"
#include "order.h"
#include "string-list.h"

/**************************************************************************
 * Read an emission distribution from a file in MEME -bfile format.
 **************************************************************************/
ARRAY_T *read_background_file(
  const char *bg_filename 		// Name of the file to read from. */
);

/***************************************************************************
 * Set a background distribution 
 *  - by reading values from a file if filename is given, or
 *  - equal to the NRDB frequencies if filename is NULL.
 *
 ***************************************************************************/
ARRAY_T *get_background(char *bg_filename);

#define WANT_PSSM 1
#define REQUIRE_PSSM 3
#define WANT_PSPM 4
#define REQUIRE_PSPM 12
#define REQUIRE_EVALUE 16
#define REQUIRE_NSITES 32
/***********************************************************************
 * Read motifs from a MEME file.
 * can read xml, html and txt.
 * Note: because of the arbitrary limit on number of motifs
 * this method is now depreciated in favor of read_meme_file2
 ***********************************************************************/
void read_meme_file (
   char*      meme_filename,          // Name of MEME file or "-" for stdin  IN
   char*      bg_filename,            // Name of background freq. file IN
   double     pseudocount,            // Value of pseudocount IN
   int        loading_mode,           // Load PSSM or PSPM (bit field) IN
   int*       num_motifs,             // Number of motifs retrieved  OUT
   MOTIF_T*   motifs,                 // The retrieved motifs 
   STRING_LIST_T** motif_occurrences, // Strings desc. motif occurrences  OUT
   BOOLEAN_T* has_reverse_strand,     // Does this file have both strands? OUT
   ARRAY_T**  background              // Background emission distribution  OUT
);

/***********************************************************************
 * Read motifs from a MEME file.
 * can read xml, html and txt.
 * Note that the motifs array list may contain motifs already.
 ***********************************************************************/
void read_meme_file2 (
   char*      meme_filename,          // Name of MEME file or "-" for stdin  IN
   char*      bg_filename,            // Name of background freq. file IN
   double     pseudocount,            // Value of pseudocount IN
   int        loading_mode,           // Load PSSM or PSPM (bit field) IN
   ARRAYLST_T* motifs,                // The retrieved motifs 
   STRING_LIST_T** motif_occurrences, // Strings desc. motif occurrences  OUT
   BOOLEAN_T* has_reverse_strand,     // Does this file have both strands? OUT
   ARRAY_T**  background              // Background emission distribution  OUT
);


/***********************************************************************
 * Free an array list and the contained motifs
 ***********************************************************************/
void free_motifs(
    ARRAYLST_T *motif_list
);

/***********************************************************************
 * Create two copies of each motif.  The new IDs are preceded by "+"
 * and "-", and the "-" version is the reverse complement of the "+"
 * version.
 * Note: because of the arbitrary limit on number of motifs
 * this methods is now depreciated in favor of add_reverse_complements2
 ***********************************************************************/
void add_reverse_complements
  (int* num_motifs,
   MOTIF_T* motifs);

/***********************************************************************
 * Create two copies of each motif.  The new IDs are preceded by "+"
 * and "-", and the "-" version is the reverse complement of the "+"
 * version.
 *
 * The reverse complement is always listed directly after the original.
 ***********************************************************************/
void add_reverse_complements2
  (ARRAYLST_T* motifs);

/***********************************************************************
 * Sets the trim boundaries for all the motifs by scanning inward until 
 * a position with an information content larger or equal to the 
 * threshold is encountered.
 ***********************************************************************/
void trim_motifs_by_bit_threshold
  (ARRAYLST_T* motifs, double threshold_bits);

/*************************************************************************
 * Setup motif-to-motif occurrence and spacer length frequency
 * transition matrices in a naive fashion, without taking into account
 * any motif occurrence information.
 *************************************************************************/
void compute_naive_transitions_and_spacers
  (const int  nmotifs,     // The number of motifs.
   MATRIX_T** transp_freq, // Motif-to-motif occurrence matrix.
   MATRIX_T** spacer_ave);  // Average spacer length matrix.

/**************************************************************************
 * Replace the elements an array of frequences with the average
 * over complementary bases.
 **************************************************************************/
void average_freq_with_complement(ARRAY_T *freqs);

#endif /* MEME_IO_H */

/*
 * Local Variables:
 * mode: c++
 * c-basic-offset: 2
 * End:
 */
