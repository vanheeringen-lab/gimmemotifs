/********************************************************************
 * FILE: alphabet.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 4-17-97
 * PROJECT: MHMM
 * COPYRIGHT: 1997-2008, WSN
 * DESCRIPTION: Define the amino acid and nucleotide alphabets.
 ********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "alphabet.h"
#include "utils.h"

#ifndef AMBIG_DEBUG
#define AMBIG_DEBUG 0
#endif

/* Define global variables to hold alphabet info. */
ALPH_T alph = INVALID_ALPH;
int alph_size = 0;
int ambigs = 0;
char alphabet[NUM_AMINOS + NUM_AMINO_AMBIGS + 1];
char any_char = '\0';

/********************************************************************
 * See .h file for description.
 ********************************************************************/
void set_alphabet
  (VERBOSE_T verbose,
   char*     given_alphabet)
{
  if (strlen(given_alphabet) == NUM_AMINOS) {
    alph = PROTEIN_ALPH;
    alph_size = NUM_AMINOS;
    ambigs = NUM_AMINO_AMBIGS;
    any_char = ANY_AMINO;
    strcpy(alphabet, given_alphabet);
    strcat(alphabet, AMINO_AMBIGS);
    if (verbose >= NORMAL_VERBOSE) {
      fprintf(stderr, "Using amino acid alphabet (%s).\n", alphabet);
    }
  } else if (strlen(given_alphabet) == NUM_BASES) {
    alph = DNA_ALPH;
    alph_size = NUM_BASES;
    ambigs = NUM_BASE_AMBIGS;
    any_char = ANY_BASE;
    strcpy(alphabet, given_alphabet);
    strcat(alphabet, BASE_AMBIGS);
    if (verbose >= NORMAL_VERBOSE) {
      fprintf(stderr, "Using nucleotide alphabet (%s).\n", alphabet);
    }
  } else if (strlen(given_alphabet) == (NUM_BASES + 1)) {
    alph = DNA_ALPH_WITH_GAP;
    alph_size = NUM_BASES + 1;
    ambigs = NUM_BASE_AMBIGS;
    any_char = ANY_BASE;
    strcpy(alphabet, given_alphabet);
    strcat(alphabet, BASE_AMBIGS);
    if (verbose >= NORMAL_VERBOSE) {
      fprintf(stderr, "Using nucleotide alphabet (%s).\n", alphabet);
    }
  } else {
    die("Unrecognized alphabet (%s).\n", given_alphabet);
  }
}

/********************************************************************
 * See .h file for description.
 ********************************************************************/
ALPH_T which_alphabet
  (void)
{
  return(alph);
}

/********************************************************************
 * See .h file for description.
 ********************************************************************/
int get_alph_size
  (const ALPH_SIZE_T which_size)
{
  if (alph == INVALID_ALPH) {
    die("Alphabet uninitialized.\n");
  }

  switch (which_size) {
  case ALPH_SIZE :
    return(alph_size);
    break;
  case ALL_SIZE :
    return(alph_size + ambigs);
    break;
  case AMBIG_SIZE :
    return(ambigs);
  default :
    die("Illegal alphabet size request.\n");
  }
  /* Unreachable. */
  return(alph_size);
}

/********************************************************************
 * See .h file for description.
 ********************************************************************/
char * get_alphabet
  (BOOLEAN_T include_ambigs)
{
  static char short_alphabet[NUM_AMINOS + 1];
  static BOOLEAN_T first_time = TRUE;

  if (alph == INVALID_ALPH) {
    die("Alphabet uninitialized.\n");
  }

  if (!include_ambigs) {
    if (first_time) {
      strncpy(short_alphabet, alphabet, alph_size);
      short_alphabet[alph_size] = '\0';
      first_time = FALSE;
    }
    return(short_alphabet);
  }
  return(alphabet);
}

/********************************************************************
 * See .h file for description
 ********************************************************************/
char get_alph_char
  (const int char_index)
{
  if (alph == INVALID_ALPH) {
    die("Alphabet uninitialized.\n");
  } 
  if ((char_index < 0) || (char_index > get_alph_size(ALL_SIZE))) {
    die("Requested character outside of alphabet (%d).\n", char_index);
  }
  return(alphabet[char_index]);
}

char get_any_char
  (void)
{
  switch(alph) {
  case INVALID_ALPH :
    die("Alphabet uninitialized.\n");
    break;
  case PROTEIN_ALPH :
    return(ANY_AMINO);
    break;
  case DNA_ALPH :
    return(ANY_BASE);
    break;
  case DNA_ALPH_WITH_GAP :
    return(ANY_BASE);
    break;
  }
  /* Unreachable. */
  return(ANY_BASE);
}

/* Test whether a given character is one of the ambiguity codes. */
BOOLEAN_T is_ambiguous
  (const char my_char)
{
  int i;
  for (i = 0; i < ambigs; i++) {
    if (my_char == alphabet[alph_size + i]) {
      return(TRUE);
    }
  }
  return(FALSE);
}


/********************************************************************
 * See .h file for description
 ********************************************************************/
int alphabet_index
  (const char   letter,    /* The letter we want to find an index for. */
   const char * alphabet)  /* Alphabet in which we're looking for an index. */ 
{
  static BOOLEAN_T first_time = TRUE; /* Is this the first call? */
  static int alph_indices[27];        /* Inverse alphabet array. */
  int        return_value;

  // The first time the function is called, compute an inverse array.
  if (first_time) {
    first_time = FALSE;

    // First fill the array with -1s.
    int i_alph;
    for (i_alph = 0; i_alph < 26; i_alph++)
      alph_indices[i_alph] = -1;

    // Compute the position of each letter in the alphabet.
    int alph_size = strlen(alphabet);
    for (i_alph = 0; i_alph < alph_size; i_alph++)
      alph_indices[(int)(alphabet[i_alph] - 'A')] = i_alph;
  }

  // Make sure we got something in A-Z.
  int index = (int)(letter - 'A');
  if ((index < 0) || (index >= 26)) {
    die("Non-alphabetic character (%c).\n", letter);
  }

  // Retrieve the index from the inverse array.
  return_value = alph_indices[index];

  //  Die if the character is not in our alphabet.
  if (return_value == -1) {
    die("Non-alphabetic character (%c).\n", letter);
  }

  return return_value;
}

/********************************************************************
 * See .h file for description.
 ********************************************************************/
void get_nrdb_frequencies
  (ARRAY_T* freqs) 
{
  if (which_alphabet() == PROTEIN_ALPH) {
    set_array_item( 0, 0.073164, freqs); /* A */
    set_array_item( 1, 0.018163, freqs); /* C */
    set_array_item( 2, 0.051739, freqs); /* D */
    set_array_item( 3, 0.062340, freqs); /* E */
    set_array_item( 4, 0.040283, freqs); /* F */
    set_array_item( 5, 0.069328, freqs); /* G */
    set_array_item( 6, 0.022428, freqs); /* H */
    set_array_item( 7, 0.056282, freqs); /* I */
    set_array_item( 8, 0.058493, freqs); /* K */
    set_array_item( 9, 0.091712, freqs); /* L */
    set_array_item(10, 0.023067, freqs); /* M */
    set_array_item(11, 0.046077, freqs); /* N */
    set_array_item(12, 0.050674, freqs); /* P */
    set_array_item(13, 0.040755, freqs); /* Q */
    set_array_item(14, 0.051897, freqs); /* R */
    set_array_item(15, 0.073802, freqs); /* S */
    set_array_item(16, 0.059411, freqs); /* T */
    set_array_item(17, 0.064362, freqs); /* V */
    set_array_item(18, 0.013341, freqs); /* W */
    set_array_item(19, 0.032682, freqs); /* Y */
  } else if (which_alphabet() == DNA_ALPH) {
    set_array_item( 0, 0.281774, freqs); /* A */
    set_array_item( 1, 0.222020, freqs); /* C */
    set_array_item( 2, 0.228876, freqs); /* G */
    set_array_item( 3, 0.267330, freqs); /* T */
  } else {
    die("Illegal alphabet size (%d).\n", alph_size);
  }
  fill_in_ambiguous_chars(FALSE, freqs);
}

void get_uniform_frequencies
  (ARRAY_T* freqs) 
{
  int i, n=0;

  if (which_alphabet() == PROTEIN_ALPH) {
    n = 20;
  } else if (which_alphabet() == DNA_ALPH) {
    n = 4;
  } else {
    die("Illegal alphabet size (%d).\n", alph_size);
  }

  for (i=0; i<n; i++) { 
    set_array_item( i, 1.0/n, freqs); 
  }
  fill_in_ambiguous_chars(FALSE, freqs);

} // get_uniform_frequencies

/********************************************************************
 * Make one position in an alphabetic frequency array the sum of
 * a set of other positions.
 ********************************************************************/
static void set_ambiguity
  (BOOLEAN_T log_space,
   char      target,
   char*     sources,
   ARRAY_T*  freqs)
{
  int    i_source;
  int    num_sources;
  char   source;
  PROB_T source_count;
  PROB_T this_count;

  /* Add up the frequencies of the source characters. */
  source_count = 0.0;
  num_sources = strlen(sources);
  for (i_source = 0; i_source < num_sources; i_source++) {
    source = sources[i_source];
    this_count = get_array_item(alphabet_index(source, alphabet), freqs);

    if (log_space) {
      source_count = LOG_SUM(source_count, this_count);
    } else {
      source_count += this_count;
    }
  }

  /* Store the result. */
  set_array_item(alphabet_index(target, alphabet), source_count, freqs);
}

/********************************************************************
 * See .h file for description.
 ********************************************************************/
void fill_in_ambiguous_chars
  (BOOLEAN_T log_space,
   ARRAY_T* freqs)  /* The emission distribution to be extended. 
		       (Must be pre-mallocked large enough to accept
		       ambiguous characters). */
{
  switch (which_alphabet()) {
  case PROTEIN_ALPH :
    set_ambiguity(log_space, 'B', "DN", freqs);
    set_ambiguity(log_space, 'Z', "EQ", freqs);
    set_ambiguity(log_space, 'U', "ACDEFGHIKLMNPQRSTVWY", freqs);
    set_ambiguity(log_space, 'X', "ACDEFGHIKLMNPQRSTVWY", freqs);
    break;

  case DNA_ALPH :
    set_ambiguity(log_space, 'U', "T", freqs);
    set_ambiguity(log_space, 'R', "GA", freqs);
    set_ambiguity(log_space, 'Y', "TC", freqs);
    set_ambiguity(log_space, 'K', "GT", freqs);
    set_ambiguity(log_space, 'M', "AC", freqs);
    set_ambiguity(log_space, 'S', "GC", freqs);
    set_ambiguity(log_space, 'W', "AT", freqs);
    set_ambiguity(log_space, 'B', "GTC", freqs);
    set_ambiguity(log_space, 'D', "GAT", freqs);
    set_ambiguity(log_space, 'H', "ACT", freqs);
    set_ambiguity(log_space, 'V', "GCA", freqs);
    set_ambiguity(log_space, 'N', "ACGT", freqs);
    break;

  default :
    die("Alphabet uninitialized.\n");
  }
  
  if (AMBIG_DEBUG) {
    fprintf(stderr, "Emission distribution: ");
    print_array(freqs, 5, 3, TRUE, stderr);
  }
}      

/********************************************************************
 * Take the counts from an ambiguous character and evenly distribute
 * them among the corresponding concreate characters.
 *
 * This function operates in log space.
 ********************************************************************/
void distribute_one_count
  (char     ambig_char,
   char*    concrete_chars,
   ARRAY_T* freqs)
{
  int    ambig_index;      /* Index of the ambiguous character. */
  PROB_T ambig_count;      /* Count of the ambiguous character. */
  int    num_concretes;    /* Number of associated definite characters. */
  int    i_concrete;       /* Which definite character are we considering? */
  char   concrete_char;    /* Current definite character. */
  int    concrete_index;   /* Index of the definite character. */
  PROB_T concrete_count;   /* Count of the definite character. */

  /* Get the count to be distributed. */
  ambig_index = alphabet_index(ambig_char, alphabet);
  ambig_count = get_array_item(ambig_index, freqs);

  /* Divide it by the number of corresponding concrete characters. */
  num_concretes = strlen(concrete_chars);
  ambig_count -= my_log2((PROB_T)num_concretes);

  /* Distribute it in equal portions to the given concrete characters. */
  for (i_concrete = 0; i_concrete < num_concretes; i_concrete++) {
    concrete_char = concrete_chars[i_concrete];
    concrete_index = alphabet_index(concrete_char, alphabet);
    concrete_count = get_array_item(concrete_index, freqs);

    /* Add the ambiguous counts. */
    concrete_count = LOG_SUM(concrete_count, ambig_count);

    set_array_item(concrete_index, concrete_count, freqs);
  }

  /* Set the ambiguous count to zero. */
  set_array_item(ambig_index, LOG_ZERO, freqs);
}
    


/********************************************************************
 * See .h file for description.
 ********************************************************************/
void distribute_ambiguous_counts
  (ARRAY_T* freqs)
{
  switch (which_alphabet()) {
  case PROTEIN_ALPH :
    distribute_one_count('B', "DN", freqs);
    distribute_one_count('Z', "EQ", freqs);
    distribute_one_count('U', "ACDEFGHIKLMNPQRSTVWY", freqs);
    distribute_one_count('X', "ACDEFGHIKLMNPQRSTVWY", freqs);
    break;

  case DNA_ALPH :
    distribute_one_count('U', "T", freqs);
    distribute_one_count('R', "GA", freqs);
    distribute_one_count('Y', "TC", freqs);
    distribute_one_count('K', "GT", freqs);
    distribute_one_count('M', "AC", freqs);
    distribute_one_count('S', "GC", freqs);
    distribute_one_count('W', "AT", freqs);
    distribute_one_count('B', "GTC", freqs);
    distribute_one_count('D', "GAT", freqs);
    distribute_one_count('H', "ACT", freqs);
    distribute_one_count('V', "GCA", freqs);
    distribute_one_count('N', "ACGT", freqs);
    break;

  default :
    die("Alphabet uninitialized.\n");
  }
}
    


/********************************************************************
 * Set all the ambiguous characters to zero.
 ********************************************************************/
void zero_ambigs
  (ARRAY_T* freqs)
{
  int i_ambig;
  int num_chars = 0;
  int num_ambigs = 0;

  switch (which_alphabet()) {
  case PROTEIN_ALPH :
    num_chars = NUM_AMINOS;
    num_ambigs = NUM_AMINO_AMBIGS;
    break;

  case DNA_ALPH :
    num_chars = NUM_BASES;
    num_ambigs = NUM_BASE_AMBIGS;
    break;

  default :
    die("Alphabet uninitialized.\n");
  }
  
  for (i_ambig = 0; i_ambig < num_ambigs; i_ambig++) {
    set_array_item(num_chars + i_ambig, 0.0, freqs);
  }
}      


  
