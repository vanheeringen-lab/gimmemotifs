/****************************************************************************
 * FILE: seq.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 06/25/2002
 * PROJECT: MHMM
 * DESCRIPTION: Biological sequences.
 * COPYRIGHT: 1998-2008, UCSD, UCSC, UW
 ****************************************************************************/
#ifndef SEQ_H
#define SEQ_H

#include "array.h"
#include "utils.h"

// A sequence object.
typedef struct seq SEQ_T;

/****************************************************************************
 * Allocate one sequence object.
 ****************************************************************************/
SEQ_T* allocate_seq
  (char* name,
   char* description,
   int   offset,
   char* sequence);

/****************************************************************************
 * Get and set various fields.
 ****************************************************************************/
char* get_seq_name
  (SEQ_T* a_sequence);
char* get_seq_description
  (SEQ_T* a_sequence);
int get_seq_length
  (SEQ_T* a_sequence);
int get_seq_offset
  (SEQ_T* a_sequence);
float get_seq_weight
  (SEQ_T* a_sequence);
void set_seq_weight
  (float  weight,
   SEQ_T* a_sequence);
unsigned char get_seq_char
  (int index,
   SEQ_T* a_sequence);
void set_seq_char
  (int    index,
   char   a_char,
   SEQ_T* a_sequence);
int get_seq_int
  (int index,
   SEQ_T* a_sequence);
int get_seq_gc
  (int index,
   SEQ_T* a_sequence);
char* get_raw_sequence
  (SEQ_T* a_sequence);
int* get_int_sequence
  (SEQ_T* a_sequence);
void allocate_int_sequence
  (SEQ_T* a_sequence);
int* get_gc_sequence
  (SEQ_T* a_sequence);
double get_total_gc_sequence
  (SEQ_T* a_sequence);
void set_total_gc_sequence
  (SEQ_T* a_sequence, double gc);

/**************************************************************************
 * Copy a sequence object.  Memory must be freed by caller.
 **************************************************************************/
SEQ_T* copy_sequence
  (SEQ_T* source_sequence);

/**************************************************************************
 * Sometimes a sequence object contains only a portion of the actual
 * sequence.  This function tells whether or not the end of this
 * sequence object corresponds to the end of the actual sequence.
 **************************************************************************/
void set_complete
  (BOOLEAN_T is_complete,
   SEQ_T* a_sequence);
BOOLEAN_T is_complete
  (SEQ_T* a_sequence);

/***************************************************************************
 * Remove Xs from either side of the sequence.
 ***************************************************************************/
void remove_flanking_xs
  (SEQ_T* sequence);

/**************************************************************************
 * Prepare a sequence for recognition by
 *  - making sure it is uppercase,
 *  - making sure it doesn't contain illegal characters,
 *  - adding flanking Xs to match START/END states, and
 *  - converting it to an integer format
 *
 * In the integer form, each character in the sequence is replaced by
 * the index of that character in the alphabet array.  Thus, if the
 * alphabet is 'ACGT', every occurence of the letter 'G' in the
 * sequence will be represented by the index 2.
 **************************************************************************/
void prepare_sequence
  (SEQ_T* sequence);

/***************************************************************************
 * Remove the first N bases of a given sequence.
 ***************************************************************************/
void shift_sequence
  (int    offset,
   SEQ_T* sequence);

/***************************************************************************
 * Get the maximum sequence length from a set of sequences.
 ***************************************************************************/
int get_max_seq_length
  (int     num_seqs,
   SEQ_T** sequences);

/***************************************************************************
 * Get the maximum sequence ID length from a set of sequences.
 ***************************************************************************/
int get_max_seq_name
  (int     num_seqs,
   SEQ_T** sequences);

/***************************************************************************
 * Set the sequence weights according to an external file.
 *
 * If the filename is "none," "internal," or NULL, then the weights are
 * set uniformly.
 ***************************************************************************/
void set_sequence_weights
  (char*    weight_filename, // Name of file containing sequence weights.
   int      num_seqs,        // Number of sequences.
   SEQ_T**  sequences);      // The sequences.

/****************************************************************************
 *  Return an array containing the frequencies in the sequence for each
 *  character of the alphabet. Characters not in the alphabet are not
 *  counted.
 ****************************************************************************/
ARRAY_T* get_sequence_freqs
  (SEQ_T* seq);

/**********************************************************************
  shuffle_sequence()

  shuffle a given sequences based on their content
**********************************************************************/
void shuffle_sequence(
  SEQ_T* seq,			/* original sequence IN */
  unsigned int seed,	/* seed IN */
  SEQ_T** target		/* target sequence OUT */
);

/****************************************************************************
 * Extract a subsequence string from a seq object.
 ****************************************************************************/
char* get_raw_subsequence(int start, int stop, SEQ_T* a_sequence);

/****************************************************************************
 * Free one sequence object.
 ****************************************************************************/
void free_seq
  (SEQ_T* a_sequence);

#endif
