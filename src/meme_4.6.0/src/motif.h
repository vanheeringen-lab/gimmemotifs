/***********************************************************************
 * FILE: motif.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 7-13-97
 * PROJECT: MHMM
 * COPYRIGHT: 1997-2008, WSN
 * DESCRIPTION: Data structure for representing one motif.
 ***********************************************************************/
#ifndef MOTIF_H
#define MOTIF_H

#include "matrix.h"
#include "utils.h"

// For reading and writing, store null ID chars as this.
#define NON_MOTIF_ID_CHAR '.'

// Maximum number of motifs.
#define MAX_MOTIFS 400

// Maximum number of characters in a motif ID.
#define MAX_MOTIF_ID_LENGTH 100
#define MAX_MOTIF_URL_LENGTH 500

typedef struct motif_t {
  char id[MAX_MOTIF_ID_LENGTH+1]; // The motif ID.
  char id2[MAX_MOTIF_ID_LENGTH+1]; // The motif secondary ID (for Tomtom).
  int length;         // The length of the motif.
  int alph_size;      // The number of letters in the alphabet.
  int ambigs;         // The number of ambiguous characters.
  double evalue;      // E-value associated with this motif.
  double num_sites;   // Number of sites associated with this motif.
  double complexity;  // Complexity of this motif.
  MATRIX_T* freqs;    // The frequency matrix. [length by alph_size 
                      //are you sure? meme-io makes this 
                      //length by (alph_size + ambigs) ]
  MATRIX_T* scores;   // Position specific scoring matrix
                      // added so MAST could use meme-io
  char *url;          // Optional url, empty string if not set.
  int trim_left;      // count of columns trimmed on left, 0 = untrimmed
  int trim_right;     // count of columns trimmed on right, 0 = untrimmed
} MOTIF_T;

/***********************************************************************
 * Get or set the index number of a motif.
 ***********************************************************************/
char* get_motif_id
  (MOTIF_T* motif);

void set_motif_id
  (char* id,
   MOTIF_T* motif);

char* get_motif_id2
  (MOTIF_T* motif);

void set_motif_id2
  (char* id,
   MOTIF_T* motif);


/***********************************************************************
 * Get motif id without leading +/- indicating strand.
 ***********************************************************************/
char *get_bare_motif_id
 (MOTIF_T* motif);

/***********************************************************************
 * allocate memmory for a motif
 ***********************************************************************/
MOTIF_T* allocate_motif(
  char *id,
  MATRIX_T* freqs);

/***********************************************************************
 * Get or set the frequencies
 ***********************************************************************/

MATRIX_T* get_motif_freqs
  (MOTIF_T* motif);

void set_motif_freqs
  (MATRIX_T* freqs,
   MOTIF_T* motif);

/***********************************************************************
 * Get the motif length
 ***********************************************************************/
double get_motif_length
  (MOTIF_T* motif);

/***********************************************************************
 * Get the motif length after trimming
 ***********************************************************************/
int get_motif_trimmed_length
  (MOTIF_T* motif);

/***********************************************************************
 * Get the motif alphabet size (non-ambiguous letters)
 ***********************************************************************/
double get_motif_alph_size
  (MOTIF_T* motif);

/***********************************************************************
 * Get the motif ambiguous alphabet size (only ambiguous letters)
 ***********************************************************************/
double get_motif_ambiguous_size
  (MOTIF_T* motif);

/***********************************************************************
 * Get the E-value of a motif.
 ***********************************************************************/
double get_motif_evalue
  (MOTIF_T* motif);

/***********************************************************************
 * Get the number of sites of a motif.
 ***********************************************************************/
double get_motif_nsites
  (MOTIF_T* motif);

/***********************************************************************
 * Get the position specific score for a letter
 ***********************************************************************/
double get_motif_score
  (MOTIF_T* motif, int position, int i_alph);

/***********************************************************************
 * Return one column of a motif, as a newly allocated array of counts.
 ***********************************************************************/
ARRAY_T* get_motif_counts
  (int      position,
   MOTIF_T* motif);

/***********************************************************************
 * Get the url of a motif
 ***********************************************************************/
char* get_motif_url
  (MOTIF_T* motif);

/***********************************************************************
 * Set the url of a motif
 ***********************************************************************/
void set_motif_url
  (char    *url,
   MOTIF_T *motif);

/***********************************************************************
 * Get the number of positions to trim from the left of the motif
 ***********************************************************************/
int get_motif_trim_left
  (MOTIF_T *motif);

/***********************************************************************
 * Get the number of positions to trim from the right of the motif
 ***********************************************************************/
int get_motif_trim_right
  (MOTIF_T *motif);

/***********************************************************************
 * Clear the motif trim
 ***********************************************************************/
void clear_motif_trim
  (MOTIF_T *motif);

/***********************************************************************
 * Determine whether a given motif is in a given list of motifs.
 ***********************************************************************/
BOOLEAN_T have_motif
  (char*    motif_id,
   int      num_motifs,
   MOTIF_T* motifs);

/***********************************************************************
 * Copy a motif from one place to another.
 ***********************************************************************/
void copy_motif
  (MOTIF_T* source,
   MOTIF_T* dest);

/***********************************************************************
 * Say that the motif ID is printed centered above a given motif.
 * This function returns the character that appears in the nth
 * position of that motif ID string.
 ***********************************************************************/
char get_motif_id_char
  (int      position,
   MOTIF_T* a_motif);

/***********************************************************************
 * Compute the reverse complement of one DNA frequency distribution.
 *
 * Assumes DNA alphabet in order ACGT.
 ***********************************************************************/
void complement_dna_freqs
  (ARRAY_T* source,
   ARRAY_T* dest);


/***********************************************************************
 * Convert array by compute the average of complementary dna frequencies.
 *
 * Assumes DNA alphabet in order ACGT.
 ***********************************************************************/
void balance_complementary_dna_freqs
  (ARRAY_T* source);

/***********************************************************************
 * Turn a given motif into its own reverse complement.
 ***********************************************************************/
void reverse_complement_motif
  (MOTIF_T* a_motif);

/***********************************************************************
 * Set the trimming bounds on the motif.
 *
 * Reads from the left and right until it finds a motif position with
 * an information content larger or equal to the threshold in bits.
 * 
 ***********************************************************************/
void trim_motif_by_bit_threshold(
  MOTIF_T *a_motif, 
  double threshold_bits
);

/***********************************************************************
 * Compute the complexity of a motif as a number between 0 and 1.
 ***********************************************************************/
double compute_motif_complexity
  (MOTIF_T* a_motif);

/***********************************************************************
 * Compute the number of positions from the start or end of a motif
 * that contain a given percentage of the information content.
 *
 * Information content is the same as relative entropy, and is computed
 * as
 *
 *  \sum_i p_i log(p_i/f_i)
 *
 ***********************************************************************/
int get_info_content_position
  (BOOLEAN_T from_start, // Count from start?  Otherwise, count from end.
   float     threshold,  // Information content threshold (in 0-100).
   ARRAY_T*  background, // Background distribution.
   MOTIF_T*  a_motif);

/***********************************************************************
 * Returns the string that is the best possible match to the given motif.
 ***********************************************************************/
char *get_best_possible_match(MOTIF_T *motif);

/***********************************************************************
 * Free dynamic memory used by a given motif.
 ***********************************************************************/
void free_motif
  (MOTIF_T * a_motif);

#endif

