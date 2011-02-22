/****************************************************************************
 * FILE: alignment.h
 * AUTHOR: William Stafford Noble, Charles E. Grant
 * CREATE DATE: 10/27/2004
 * PROJECT: EVOMCAST
 * DESCRIPTION: Multiple alignment of biological sequences.
 * COPYRIGHT: 1998-2008, UCSD, UCSC, UW
 ****************************************************************************/
#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include "matrix.h"
#include "array.h"
#include "string-list.h"
#include "seq.h"
#include "object-list.h"

// An alignment object.
typedef struct alignment ALIGNMENT_T;

/****************************************************************************
 * Allocate one alignment object.
 ****************************************************************************/
ALIGNMENT_T* allocate_alignment(
  char* name,
  char* description,
  int num_sequences,
  SEQ_T** sequences,
  char* consensus_string
);

/****************************************************************************
 * Get and set various fields.
 ****************************************************************************/
void set_alignment_name(char* name, ALIGNMENT_T* an_alignment);
char* get_alignment_name(ALIGNMENT_T* an_alignment);
char* get_alignment_description(ALIGNMENT_T* an_alignment);
char* get_consensus_string(ALIGNMENT_T* an_alignment);
int get_alignment_length(ALIGNMENT_T* an_alignment);
int get_num_aligned_sequences(ALIGNMENT_T* an_alignment);
int get_num_identical_sites(ALIGNMENT_T* an_alignment);
int get_num_conserved_sites(ALIGNMENT_T* an_alignment);
int get_num_semiconserved_sites(ALIGNMENT_T* an_alignment);
int get_num_nonconserved_sites(ALIGNMENT_T* an_alignment);
SEQ_T* get_consensus_sequence(double threshold, ALIGNMENT_T* an_alignment);
SEQ_T* get_alignment_sequence(int index, ALIGNMENT_T* an_alignment);
SEQ_T* get_alignment_sequence_by_name(char* name, ALIGNMENT_T* an_alignment);
SEQ_T** get_alignment_sequences(ALIGNMENT_T* an_alignment);

/****************************************************************************
 * Fill in a null terminated string with the bases in one column of the 
 * alignment. The user must allocate the memory for the string, which should
 * be large enough to store one characters from each sequence in the alignment
 * plus the trailing null character. This is done for reasons of efficiency,
 * since in most cases the user will be making for this call iterively over
 * the length of the alignment.
 ****************************************************************************/
void get_alignment_col(
  int col, 
  char* alignment_col, 
  ALIGNMENT_T* an_alignment
);

/*************************************************************************
 * Convert the string representing an alignment column into an integer
 * which will be the column index for that alignment column in the PSSM.
 * If the alphabet has m characters, and the alignment columns have n entries,
 * the array of all alignment columns is conveniently numbered by the set of
 * consecutive n-digit base m numerals:
 *   AAAA = 0000, AAAC = 0001, ..., TTTG = 3332, TTTT = 3333.
 *************************************************************************/
int hash_alignment_col(char* alignment_col, int alignment_col_size); 

/*************************************************************************
 * Convert an integer representing a column in a PSSM into the
 * corresponding alignment column string.
 * If the alphabet has m characters, and the alignment columns have n entries,
 * the array of all alignment columns is conveniently numbered by the set of
 * consecutive n-digit base m numerals:
 *   AAAA = 0000, AAAC = 0001, ..., TTTG = 3332, TTTT = 3333.
 * The caller must allocate the memory for the alignment column string.
 * The memory required is the number of sequences in the alignment, plus one
 * for the terminating null.
 *************************************************************************/
void unhash_alignment_col(
  int alignment_col_index,
  char *alignment_col,
  int alignment_col_size
);

/****************************************************************************
 *  Return an array containing the frequencies in the alignment for each 
 *  character of the alphabet. Gaps not counted. 
 ****************************************************************************/
ARRAY_T* get_alignment_freqs(ALIGNMENT_T* an_alignment);

/****************************************************************************
 *  Return an list containing the empirical column frequency distributions
 *  for all alignments in the input.
 ****************************************************************************/
OBJECT_LIST_T* get_alignment_column_freqs_list
  (STRING_LIST_T* filenames,
  BOOLEAN_T remove_allgap_seqs);


/****************************************************************************
 *  Get a cumulative count of gaps within one sequence of the alignment
 ****************************************************************************/
int* get_cumulative_gap_count(int seqIndex, ALIGNMENT_T* alignment);

/****************************************************************************
 *  Does a column of an alignment contain gaps?
 ****************************************************************************/
BOOLEAN_T alignment_site_has_gaps(int index, ALIGNMENT_T* an_alignment);

/****************************************************************************
 *  Does a column of an alignment contain any ambiguity codes?
 ****************************************************************************/
BOOLEAN_T alignment_site_ambiguous(int index, ALIGNMENT_T* alignment);

/****************************************************************************
 * Create a lookup table for converting an index into an alignment to an index
 * into a gapless version of one of the sequences in the alignment.
 ****************************************************************************/
int* make_alignment_to_seq_table(int seq_index, ALIGNMENT_T* an_alignment); 

/****************************************************************************
 * Create a lookup table for converting an index into a sequence to an index
 * into the alignment. Note that because there are many alignment positions
 * that correspond to a sequence position we take the first occurence.
 * JCH: I have added this function for the sake of the BLS scan mode
 * so that single mode matches in each sequence can be mapped back
 * to positions in the alignment.
 ****************************************************************************/
int* make_seq_to_alignment_table(int ref_seq_index, ALIGNMENT_T* an_alignment);

/****************************************************************************
 * Get a list of the names of the species in the alignment.
 ****************************************************************************/
STRING_LIST_T* get_species_names(ALIGNMENT_T* an_alignment);

/****************************************************************************
 * Count the number of non-gap characters in a sequence.
 ****************************************************************************/
int count_residues(char* seq);

/****************************************************************************
 * Extract a small alignment out of the middle of a larger alignment.
 ****************************************************************************/
ALIGNMENT_T* extract_subalignment
  (int start,
   int width,
   ALIGNMENT_T* alignment);

/****************************************************************************
 * Remove from the alignment all columns that contain gaps for the
 * specified species.
 ****************************************************************************/
ALIGNMENT_T* remove_alignment_gaps
  (char*        species,
   ALIGNMENT_T* alignment);

/****************************************************************************
 * Remove from an alignment any sequence whose ID is not in a given list.
 *
 * N.B. It is NOT an error for the given list to contain sequence IDs that 
 * are not in the alignment.
 ****************************************************************************/
ALIGNMENT_T* remove_alignment_seqs
  (STRING_LIST_T* seqs_to_keep,
   ALIGNMENT_T*   alignment);

/****************************************************************************
 * Read an alignment from a file.  Sort the sequences by sequence name.
 ****************************************************************************/
ALIGNMENT_T* read_alignment_from_file
  (char *filename, 
   BOOLEAN_T sort,
   BOOLEAN_T remove_allgap_seqs,
   int* ref_seq_index
  );

/*************************************************************************
 * Print an alignment in PHYLIP format.
 *************************************************************************/
void print_phylip_alignment
  (ALIGNMENT_T* the_alignment,
   FILE* outfile);

/****************************************************************************
 * Free one alignment object.
 ****************************************************************************/
void free_alignment(ALIGNMENT_T* an_alignment);


#endif
