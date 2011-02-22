/********************************************************************
 * FILE: alphabet.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 4-17-97
 * PROJECT: MHMM
 * COPYRIGHT: 1997-2008, WSN
 * DESCRIPTION: Define the amino acid and nucleotide alphabets.
 ********************************************************************/
#ifndef ALPHABET_H
#define ALPHABET_H

#include "array.h"
#include "utils.h"

/* Define a type and a global variable for the alphabet type. */
typedef enum {INVALID_ALPH, PROTEIN_ALPH, DNA_ALPH, DNA_ALPH_WITH_GAP } ALPH_T;

/********************************************************************
  The nucleic acid codes supported are:

        A --> adenosine           M --> A C (amino)
        C --> cytidine            S --> G C (strong)
        G --> guanine             W --> A T (weak)
        T --> thymidine           B --> G T C
        U --> T                   D --> G A T
        R --> G A (purine)        H --> A C T
        Y --> T C (pyrimidine)    V --> G C A
        K --> G T (keto)          N --> A G C T (any)

 ********************************************************************/
#define NUM_BASES 4
#define NUM_BASE_AMBIGS 12
#define BASE_AMBIGS "URYKMSWBDHVN"
#define ANY_BASE 'N'


/********************************************************************
  The accepted amino acid codes are:

    A  alanine                         P  proline
    B  aspartate or asparagine         Q  glutamine
    C  cystine                         R  arginine
    D  aspartate                       S  serine
    E  glutamate                       T  threonine
    F  phenylalanine                   U  any
    G  glycine                         V  valine
    H  histidine                       W  tryptophan
    I  isoleucine                      Y  tyrosine
    K  lysine                          Z  glutamate or glutamine
    L  leucine                         X  any
    M  methionine
    N  asparagine
 ********************************************************************/
#define NUM_AMINOS 20
#define NUM_AMINO_AMBIGS 4
#define AMINO_AMBIGS "BUXZ"
#define ANY_AMINO 'X'

/* Set the alphabet according to the size. */
void set_alphabet
  (VERBOSE_T verbose,
   char*     given_alphabet);

/* Find out which alphabet is in use. */
ALPH_T which_alphabet
  (void);

/* Get the size of the alphabet. */
typedef enum {INVALID_SIZE, ALPH_SIZE, AMBIG_SIZE, ALL_SIZE} ALPH_SIZE_T;
int get_alph_size
  (const ALPH_SIZE_T which_size);

/* Get the actual alphabet, including ambiguous characters. */
char * get_alphabet
  (BOOLEAN_T include_ambigs);

/* Get the nth character from the alphabet. */
char get_alph_char
  (const int char_index);

/* Get the character that represents complete ambiguity for this alphabet. */
char get_any_char
  (void);

/* Test whether a given character is one of the ambiguity codes. */
BOOLEAN_T is_ambiguous
  (const char my_char);

/********************************************************************
 * int alphabet_index
 *
 * Given a character, returns the index of that character in the given
 * alphabet.  Assumes that the character is an uppercase alphabetic
 * character.  If the character does not appear in the alphabet, die.
 * 
 * This function assumes that the alphabet does not change from one
 * call to the next.
 *
 * RETURN: The index of the given letter in the given alphabet. 
 ********************************************************************/
int alphabet_index
  (const char   letter,    /* The letter we want to find an index for. */
   const char * alphabet); /* Alphabet in which we're looking for an index. */ 

/********************************************************************
 * Fill a given array with frequencies from the non-redundant database
 * 9/22/96.  These were lifted from MEME 2.0 (init.c).
 *
 * 2/3/98: I noticed that they don't sum to 1.0, so I normalized them.
 ********************************************************************/
void get_nrdb_frequencies
  (ARRAY_T* freqs);

/********************************************************************/
/*	get_uniform_frequencies
*
*	Fill in array with uniform frequencies for the current 
*	alphabet.
*/
/********************************************************************/
void get_uniform_frequencies
  (ARRAY_T* freqs);

/********************************************************************
 * Given an emission distribution, append additional frequencies to
 * the end to account for ambiguous characters in the alphabet.
 ********************************************************************/
void fill_in_ambiguous_chars
  (BOOLEAN_T log_space,
   ARRAY_T*  freqs); /* The emission distribution to be extended. 
		       (Must be pre-mallocked large enough to accept
		       ambiguous characters). */

/********************************************************************
 * Given an emission distribution containing log counts, distribute the
 * counts for ambiguous character evenly among the corresponding
 * concrete characters, leaving zeroes in the ambiguous slots.
 ********************************************************************/
void distribute_ambiguous_counts
  (ARRAY_T* freqs);

/********************************************************************
 * Set all the ambiguous characters to zero.
 ********************************************************************/
void zero_ambigs
  (ARRAY_T* freqs);

#endif
