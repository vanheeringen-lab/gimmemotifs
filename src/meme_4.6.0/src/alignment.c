/****************************************************************************
 * FILE: alignment.c
 * AUTHOR: William Stafford Noble, Charles E. Grant
 * CREATE DATE: 10/28/2004
 * PROJECT: EVOMCAST
 * DESCRIPTION: Multiple alignment of biological sequences.
 * COPYRIGHT: 1998-2008, UCSD, UCSC, UW
 ****************************************************************************/
#include <assert.h>
#include <string.h>
#include "alignment.h"
#include "alphabet.h"
#include "utils.h"
#include "string-list.h"
#include "object-list.h"
#include "clustalw-io.h"

#define MAX_ALIGNMENT_NAME 40     // Longest alignment ID.
#define MAX_ALIGNMENT_COMMENT 128 // Longest comment.

extern char alphabet[];

// Instantiate the ALIGNMENT_T type.
struct alignment {
  char  name[MAX_ALIGNMENT_NAME + 1];     // Alignment program name.
  char  desc[MAX_ALIGNMENT_COMMENT + 1];  // Alignment description.
  int   length;                 // Length of the Alignment
  int   num_sequences;          // Number of sequences in the alignment
  SEQ_T** sequences;             // The aligned sequences.
  char* consensus_string;       // The string representing the alignment
  // N.B. At some point, we'll need to assign an alignment ID.
};

/****************************************************************************
 * Allocate one alignment object. Name and description may be NULL.
 * 
 * Returns a pointer to the newly created alignment.
 ****************************************************************************/
ALIGNMENT_T* allocate_alignment(
   char* name,
   char* description,
   int num_sequences,
   SEQ_T** sequences,
   char* consensus_string
)
{
  assert(num_sequences > 0);
  assert(sequences != NULL);
  assert(consensus_string != NULL);

  // Allocate the alignment object.
  ALIGNMENT_T* new_alignment = (ALIGNMENT_T*)mm_malloc(sizeof(ALIGNMENT_T));
  if (new_alignment == NULL) {
    die("Error allocating alignment\n");
  }

  // Store the name, truncating if necessary.
  if (name != NULL) {
    strncpy(new_alignment->name, name, MAX_ALIGNMENT_NAME);
    new_alignment->name[MAX_ALIGNMENT_NAME] = '\0';
    if (strlen(new_alignment->name) != strlen(name)) {
      fprintf(stderr, "Warning: truncating alignment program name %s to %s.\n",
	      name, new_alignment->name);
    }
  }

  // Store the description, truncating if necessary.
  if (description != NULL) {
    strncpy(new_alignment->desc, description, MAX_ALIGNMENT_COMMENT);
    new_alignment->desc[MAX_ALIGNMENT_COMMENT] = '\0';
  }

  // Store the sequences.
  new_alignment->sequences = (SEQ_T**) mm_malloc(num_sequences * sizeof(SEQ_T*));
  if (new_alignment->sequences == NULL) {
    die("Error allocating sequences\n");
  }
  new_alignment->num_sequences = num_sequences;
  int seq_length = strlen(get_raw_sequence(sequences[0]));
  int i;
  for (i = 0; i < num_sequences; i++) {
    myassert(TRUE,
	     strlen(get_raw_sequence(sequences[i])) == seq_length,
	     "Sequence #1 (%s) is length=%d, but sequence #%d (%s) is length=%d.\n<%s>\n",
	     get_seq_name(sequences[0]), seq_length, i, 
	     get_seq_name(sequences[i]), strlen(get_raw_sequence(sequences[i])),
	     get_raw_sequence(sequences[i]));
    new_alignment->sequences[i] = 
      allocate_seq(get_seq_name(sequences[i]),
        get_seq_description(sequences[i]),
        get_seq_offset(sequences[i]), 
        get_raw_sequence(sequences[i])
      );
  }

  // Fill in the remaining fields.
  new_alignment->length = seq_length;
  copy_string(&(new_alignment->consensus_string), consensus_string);

  return(new_alignment);
}

/****************************************************************************
 * Get and set various fields.
 ****************************************************************************/
char* get_alignment_name(ALIGNMENT_T* alignment) {
  assert(alignment != NULL);
  return(alignment->name);
}

void set_alignment_name(char* name, ALIGNMENT_T* alignment) {
  assert(name != NULL);
  assert(alignment != NULL);
  alignment->name[0] = 0;
  strncat(alignment->name, name, MAX_ALIGNMENT_NAME); 
}

char* get_alignment_description(ALIGNMENT_T* alignment) {
  assert(alignment != NULL);
  return(alignment->desc);
}

char* get_consensus_string(ALIGNMENT_T* alignment) {
  assert(alignment != NULL);
  return(alignment->consensus_string);
}

SEQ_T* get_consensus_sequence(double threshold, ALIGNMENT_T* alignment) {
  char* seq_string = NULL;
  unsigned char c = 0;
  unsigned char most_freq_char = 0;
  #define NUM_CHARS 127
  char char_counts[NUM_CHARS];
  int i = 0;
  int j = 0;
  double max_char_freq = 0.0;
  SEQ_T* consensus;
  assert(alignment != NULL);
  
  seq_string = mm_malloc(alignment->length * sizeof(char) + 1);
  if (seq_string == NULL) {
    die("Error allocating consensus sequence string\n");
  }
  // For each column in the alignment
  for (i = 0; i < alignment->length; i++) {
    most_freq_char = 0;
    memset(char_counts, 0, NUM_CHARS * sizeof(char));
    // Count character occurances
    for (j = 0; j < alignment->num_sequences; j++) {
      c = get_seq_char(i, alignment->sequences[j]);
      char_counts[c]++;
    }
    // Find the index of the character that occurs the most frequently
    for (c = 0; c < NUM_CHARS; c++) {
      most_freq_char = char_counts[most_freq_char] >= char_counts[c] ? 
        most_freq_char : c;
    }
    // If the most frequent character exceeds the threshold
    // it will be the consensus character.
    max_char_freq = (double) char_counts[most_freq_char] / 
      (double) alignment->num_sequences;
    if (max_char_freq >= threshold) {
      seq_string[i] = most_freq_char;
    } else {
      // Otherwise the consensus is the gap character
      seq_string[i] = '-';
    }
  }
  seq_string[i] = '\0';
  consensus = allocate_seq("Consensus", "", 0, seq_string);
  if (seq_string != NULL) myfree(seq_string);
  return(consensus);
}

int get_alignment_length(ALIGNMENT_T* alignment) {
  assert(alignment != NULL);
  return(alignment->length);
}

int get_num_aligned_sequences(ALIGNMENT_T* alignment) {
  assert(alignment != NULL);
  return(alignment->num_sequences);
}

SEQ_T** get_alignment_sequences(ALIGNMENT_T* alignment) {
  assert(alignment != NULL);
  assert(alignment->sequences != NULL);
  return(alignment->sequences);
}

SEQ_T* get_alignment_sequence(int index, ALIGNMENT_T* alignment) {
  assert(alignment != NULL);
  assert(alignment->sequences != NULL);
  assert(index >= 0);
  assert(index < alignment->num_sequences);
  return(alignment->sequences[index]);
}

SEQ_T* get_alignment_sequence_by_name(char* name, ALIGNMENT_T* alignment) {
  int i = 0;
  SEQ_T* sequence = NULL;

  assert(alignment != NULL);
  assert(alignment->sequences != NULL);
  assert(name != NULL);

  for (i = 0; i < alignment->num_sequences; i++) {
    if (strcmp(name, get_seq_name(alignment->sequences[i])) == 0) {
      sequence = alignment->sequences[i];
      break;
    }
  }

  return sequence;
}

int get_num_identical_sites(ALIGNMENT_T* alignment) {
  char *c = NULL;
  int length = 0;
  int num_identical_sites = 0;

  assert(alignment != NULL);
  assert(alignment->consensus_string != NULL);

  c = alignment->consensus_string;
  while (*c != '\0') {
    length++;
    if (*c == '*') {
      num_identical_sites++;
    }
    c++;
  }
  assert(length == alignment->length);
  return(num_identical_sites);
}

int get_num_conserved_sites(ALIGNMENT_T* alignment) {
  char *c = NULL;
  int length = 0;
  int num_conserved_sites = 0;
  
  assert(alignment != NULL);

  c = alignment->consensus_string;
  while (*c != '\0') {
    length++;
    if (*c == ':') {
      num_conserved_sites++;
    }
    c++;
  }
  assert(length == alignment->length);
  return(num_conserved_sites);
}

int get_num_semiconserved_sites(ALIGNMENT_T* alignment) {
  char *c = NULL;
  int length = 0;
  int num_semiconserved_sites = 0;

  assert(alignment != NULL);

  c = alignment->consensus_string;
  while (*c != '\0') {
    length++;
    if (*c == '.') {
      num_semiconserved_sites++;
    }
    c++;
  }
  assert(length == alignment->length);
  return(num_semiconserved_sites);
}

int get_num_nonconserved_sites(ALIGNMENT_T* alignment) {
  char *c = NULL;
  int length = 0;
  int num_nonconserved_sites = 0;
  
  assert(alignment != NULL);

  c = alignment->consensus_string;
  while (*c != '\0') {
    length++;
    if (*c == '.') {
      num_nonconserved_sites++;
    }
    c++;
  }
  assert(length == alignment->length);
  return(num_nonconserved_sites);
}

int count_residues(char* seq) {

  int count = 0;

  while(*seq != '\0') {
    if (*seq != '-' && *seq != '.') count++;
    seq++;
  }

  return count;
}

/****************************************************************************
 * Fill in a null terminated string with the bases in one column of the 
 * alignment. The user must allocate the memory for the string, which should
 * be large enough to store one characters from each sequence in the alignment
 * plus the trailing null character. This is done for reasons of efficiency,
 * since in most cases the user will be making for this call iteratively over
 * the length of the alignment.
 ****************************************************************************/
void get_alignment_col(int col, char* alignment_col, ALIGNMENT_T* alignment) {

  int i;

  for (i = 0; i < alignment->num_sequences; i++, alignment_col++) {
    *alignment_col = get_seq_char(col, alignment->sequences[i]);
  }
  *alignment_col = '\0';
}

/*************************************************************************
 * Convert the string representing an alignment column into an integer
 * which will be the column index for that alignment column in the PSSM.
 * If the alphabet has m characters, and the alignment columns have n entries,
 * the array of all alignment columns is conveniently numbered by the set of
 * consecutive n-digit base m numerals: 
 *   AAAA = 0000, AAAC = 0001, ..., TTTG = 3332, TTTT = 3333.
 *************************************************************************/
int hash_alignment_col(char* alignment_col, int alignment_col_size) {

  // The alignment column must have at least
  // one character aside from the terminating null.
  assert(alignment_col_size >= 1);
  assert(alignment_col != NULL);

  char* alphabet = get_alphabet(FALSE);
  int alph_size = get_alph_size(ALPH_SIZE);
  int hash = 0;
  int place = 1;
  int i;
  for (i = alignment_col_size - 1; i >= 0; i--) {
    // Find alignment column character in alphabet
    int j;
    for (j = 0; j < alph_size; j++) {
      if (alignment_col[i] == alphabet[j]) {
        break;
      }
    }
    if (j == alph_size)  {
      die("Alignment character %c not found in alphabet\n", alignment_col[i]);
    }
    hash += j * place;
    place *= alph_size;
  }

  return hash;
} // hash_alignment_col

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
) {
  char* alphabet = get_alphabet(FALSE);
  int alph_size = get_alph_size(ALPH_SIZE);

  assert(alignment_col_index >= 0);
  assert(
    alignment_col_index < pow(
      (double) alph_size, 
      (double) alignment_col_index
    )
  );
  assert(alignment_col != NULL);
  assert(alignment_col_size >= 1);
  
  alignment_col[alignment_col_size] = '\0';
  int i, j;
  for (i = alignment_col_size - 1; i >= 0; i--) {
    j = alignment_col_index % alph_size;
    alignment_col_index -= j;
    alignment_col[i] = alphabet[j];
    alignment_col_index /= alph_size;
  }
} // unhash_alignment_col

/*************************************************************************
 *  Build array containing the counts of columns in the alignment
 *  Caller is responsible for freeing the returned array.
 *  If input parameter "freqs" is NULL, allocates the array.
 *  Otherwise, the counts are added to the existing counts in the counts
 *  array.  Ignores all columns containing gaps or ambiguity characters:
 *    [.-nNxX]
 *************************************************************************/
static ARRAY_T* build_alignment_column_counts(
  ALIGNMENT_T* alignment,
  ARRAY_T* counts 
) 
{

  assert(alignment != NULL);

  int alph_size = get_alph_size(ALPH_SIZE);

  // Calculate number of possible alignment columns
  // and create storage for counting occurences.
  int num_seqs = get_num_aligned_sequences(alignment);
  int num_alignment_cols = (int) pow((double) alph_size, (double) num_seqs);
  if (counts == NULL) {
    counts = allocate_array(num_alignment_cols);
  }

  // Count how many examples of each column occur in the alignment.
  // Skip columns that contain gaps or ambiguity characters.
  int alignment_length = get_alignment_length(alignment);
  char* alignment_col = mm_malloc(sizeof(char) * (num_seqs + 1));
  alignment_col[num_seqs] = 0;
  int i, h;
  for(i = 0; i < alignment_length; i++) {
    get_alignment_col(i, alignment_col, alignment);
    if (strchr(alignment_col, '-') != NULL) { continue; }
    if (strchr(alignment_col, '.') != NULL) { continue; }
    if (strchr(alignment_col, 'N') != NULL) { continue; }
    if (strchr(alignment_col, 'n') != NULL) { continue; }
    if (strchr(alignment_col, 'X') != NULL) { continue; }
    if (strchr(alignment_col, 'x') != NULL) { continue; }
    h = hash_alignment_col(alignment_col, num_seqs);
    incr_array_item(h, 1, counts);
  }

  return counts;
} // build_alignment_column_counts

/****************************************************************************
 *  Return a list containing the empirical column frequency distributions
 *  for all alignments in the input.
 *
 *  Each file in the list of filenames is read and the species list is
 *  determined.  The counts of each occurring column are tallied.  
 *  All files with the same species lists get their counts combined.
 *
 *  The returned list contains one distribution per species list that 
 *  occurs in some alignment.
 ****************************************************************************/
OBJECT_LIST_T* get_alignment_column_freqs_list
  (STRING_LIST_T* filenames,
  BOOLEAN_T remove_allgap_seqs) 
{
  int file_index;
  int num_filenames = get_num_strings(filenames);
  ARRAY_T* alignment_column_freqs = NULL;
  OBJECT_LIST_T* alignment_column_freqs_list
    = new_object_list(equal_string_lists,
		      (void*)copy_string_list,
		      free_string_list,
		      free_array);

  // Consider each alignment in turn.
  for(file_index = 0; file_index < num_filenames; file_index++) { 
    char* filename = get_nth_string(file_index, filenames);

    if (verbosity >= NORMAL_VERBOSE && !(file_index % 1)) {
      fprintf(
	stderr, 
	"Computing column freqs: alignment file number %d of %d total files.\n",
	file_index+1, num_filenames
      );
    }

    // Read the alignment
    int ref_seq_index = 0;
    ALIGNMENT_T* alignment = 
      read_alignment_from_file(filename, TRUE, remove_allgap_seqs, &ref_seq_index);
    STRING_LIST_T* alignment_species = get_species_names(alignment);

    // Try to retrieve the counts so far for this list of species.
    alignment_column_freqs = 
      (ARRAY_T*)retrieve_object(
	alignment_species, 
	alignment_column_freqs_list
      );
    
    // Found counts for current species list?
    if (alignment_column_freqs) {
      // Add counts from current alignment.
      (void) build_alignment_column_counts(alignment, alignment_column_freqs);
      // Note: objects in lists are references, so no need to re-store
      // after modification.
    } 
    // Didn't find counts for this species list, so create new array of counts.
    else {
      alignment_column_freqs = build_alignment_column_counts(alignment, NULL);
      store_object(
	(void*)alignment_column_freqs,
	(void*)alignment_species,
	0.0, // Score
	alignment_column_freqs_list
      );
    }
    // free space used by alignment
    free_alignment(alignment);
  } // each filename
  fprintf(stderr, "\n");

  // Convert counts to frequencies by retrieving each array of counts
  // and dividing by the total counts for that list of species.
  while ( 
    (alignment_column_freqs = retrieve_next_object(alignment_column_freqs_list)
    ) != NULL )
  {
    int i;
    int num_freqs = get_array_length(alignment_column_freqs);
    double total_counts;

    // Get total counts.
    for (i=total_counts=0; i<num_freqs; i++) {
      total_counts += get_array_item(i, alignment_column_freqs);
    }

    // Get frequencies.
    for (i=0; i<num_freqs; i++) {
      double f = get_array_item(i, alignment_column_freqs);
      set_array_item(i, f/total_counts, alignment_column_freqs);

#ifdef DEBUG
      int alph_size = get_alph_size(ALPH_SIZE);
      int num_leaves = NINT(log(num_freqs)/log(alph_size));
      char* alignment_col = mm_malloc((num_leaves + 1) * sizeof(char));
      unhash_alignment_col(
        i, 				//col_index
	alignment_col,
	num_leaves
      );
      printf("%s %g %g\n", alignment_col, f, f/total_counts);
      myfree(alignment_col);
#endif
    } // get frequencies
  } // while more species lists

  return(alignment_column_freqs_list);
} // get_alignment_column_freqs_list
  
/****************************************************************************
*  Return an array containing the frequencies in the alignment for each 
*  character of the alphabet. Gaps and ambiguity characters other then
*  ANY_BASE are not counted. The freq. of ANY_BASE characters is stored
*  in the last element of the array.
****************************************************************************/
ARRAY_T* get_alignment_freqs(ALIGNMENT_T* alignment) {
  char c = 0;
  int alph_index = 0;
  int alph_size = 0;
  int i = 0;
  int s = 0;
  int total_bases = 0;
  int* num_bases = NULL;
  ARRAY_T* freqs = NULL;
  
  // Initialize counts for each character in the alphabet
  alph_size = get_alph_size(ALPH_SIZE);
  num_bases = mm_malloc(alph_size * sizeof(int));
  for (i = 0; i < alph_size; i++) {
    num_bases[i] = 0;
  }

  for (s = 0; s < alignment->num_sequences; s++) {
    for (i = 0; i < alignment->length; i++) {
      c = get_seq_char(i, alignment->sequences[s]);
      if (c != '-' && c != '.') {
        alph_index = alphabet_index(c, alphabet);
        // c might be an ambiguity code. We don't count ambiguity codes.
        if (alph_index < alph_size) {
          num_bases[alph_index]++;
          total_bases++;
        } 
      }
    }
  }

  freqs = allocate_array(alph_size);
  for (i = 0; i < alph_size; i++) {
    set_array_item(i, (double) num_bases[i] / (double) total_bases, freqs);
  }

  // Clean up the count of characters
  myfree(num_bases);

  return freqs;
}

/****************************************************************************
 *  Get a cumulative count of gaps within one sequence of the alignment
 ****************************************************************************/
int* get_cumulative_gap_count(int seqIndex, ALIGNMENT_T* alignment) {

  int* results = (int *)mm_malloc( sizeof(int) * alignment->length );

  int i = 0;
  char c = get_seq_char(i, alignment->sequences[seqIndex]);
  if (c == '-' || c == '.') {
      results[i] = 1;
  } else {
      results[i] = 0;
  }
  for (i = 1; i < alignment->length; i++) {
    c = get_seq_char(i, alignment->sequences[seqIndex]);
    if (c == '-' || c == '.') {
      results[i] = results[i-1] + 1;
    } else {
      results[i] = results[i-1];
    }
  }

  return results;
}


/****************************************************************************
 *  Does a column of an alignment contain gaps?
 ****************************************************************************/
BOOLEAN_T alignment_site_has_gaps(int index, ALIGNMENT_T* alignment) {

  int i = 0;
  for (i = 0; i < alignment->num_sequences; i++) {
    char c = get_seq_char(index, alignment->sequences[i]);
    if (c == '-' || c == '.') {
      return(TRUE);
    }
  }
  return(FALSE);
}

/****************************************************************************
 *  Does a column of an alignment contain any ambiguity codes?
 ****************************************************************************/
BOOLEAN_T alignment_site_ambiguous(int index, ALIGNMENT_T* alignment) {

  int i = 0;
  for (i = 0; i < alignment->num_sequences; i++) {
    char c = get_seq_char(index, alignment->sequences[i]);
    if (is_ambiguous(c)) {
      return(TRUE);
    }
  }
  return(FALSE);
}

/****************************************************************************
 * Create a lookup table for converting an index into an alignment to an index
 * into a gapless version of one of the sequences in the alignment.
 ****************************************************************************/
int* make_alignment_to_seq_table(int ref_seq_index, ALIGNMENT_T* an_alignment) {

  char* raw_seq = NULL;
  int align_length = 0;
  int align_index = 0;
  int seq_index = 0;
  int *table = NULL;
  SEQ_T* seq = NULL;

  align_length = get_alignment_length(an_alignment);
  seq = get_alignment_sequence(ref_seq_index, an_alignment);
  raw_seq = get_raw_sequence(seq);

  // Table is indexed by position in the alignment
  // Table values are the corresponding position in the
  // reference sequence not including gaps.
  table = (int *) mm_malloc((align_length) * sizeof(int));
  seq_index = 0;
  for (align_index = 0; align_index < align_length; align_index++) {
    if (raw_seq[align_index] != '-' && raw_seq[align_index] != '.') {
      seq_index++;
    }
    table[align_index] = seq_index;
  }

  return table;

}


/****************************************************************************
 * Create a lookup table for converting an index into a sequence to an index
 * into the alignment. Note that because there are many alignment positions
 * that correspond to a sequence position we take the first occurence.
 * JCH: I have added this function for the sake of the BLS scan mode
 * so that single mode matches in each sequence can be mapped back
 * to positions in the alignment.
 ****************************************************************************/
int* make_seq_to_alignment_table(int ref_seq_index, ALIGNMENT_T* an_alignment) {

  char* raw_seq = NULL;
  int align_length = 0;
  int align_index = 0;
  int seq_index = 0;
  int *table = NULL;
  SEQ_T* seq = NULL;

  align_length = get_alignment_length(an_alignment);
  seq = get_alignment_sequence(ref_seq_index, an_alignment);
  raw_seq = get_raw_sequence(seq);

  // Table is indexed by position in the sequence
  // Table values are the first corresponding
  // position in the alignment.
  table = (int *) mm_malloc((align_length) * sizeof(int));
  seq_index = 0;
  table[seq_index] = 0;

  for (align_index = 0; align_index < align_length; align_index++) {
    if (raw_seq[align_index] != '-' && raw_seq[align_index] != '.') {
      seq_index++;
      table[seq_index] = align_index;
    }
  }

  return table;

}



/****************************************************************************
 * Get a list of the names of the species in the alignment.
 ****************************************************************************/
STRING_LIST_T* get_species_names(ALIGNMENT_T* an_alignment) {
  STRING_LIST_T* return_value;
  int i_seq;
  int num_seqs;

  // Allocate a new string list.
  return_value = new_string_list();

  // Extract all the sequence names and add them to the list.
  num_seqs = get_num_aligned_sequences(an_alignment);
  for (i_seq = 0; i_seq < num_seqs; i_seq++) {
    add_string(get_seq_name(get_alignment_sequence(i_seq, an_alignment)),
	       return_value);
  }

  return(return_value);
}

/****************************************************************************
 * Extract a small alignment out of the middle of a larger alignment.
 ****************************************************************************/
ALIGNMENT_T* extract_subalignment
  (int start,
   int width,
   ALIGNMENT_T* alignment)
{
  int num_sequences = get_num_aligned_sequences(alignment);
  SEQ_T** sequences = get_alignment_sequences(alignment);
  SEQ_T** subsequences = (SEQ_T**)mm_malloc(num_sequences * sizeof(SEQ_T*));

  // Extract the specified columns into a new list of sequences.
  int i_seq = 0;
  char* subsequence = mm_malloc((width + 1) * sizeof(char));
  for (i_seq = 0; i_seq < num_sequences; i_seq++) {
    SEQ_T* this_seq = sequences[i_seq];
    char* raw_seq = get_raw_sequence(this_seq);
    strncpy(subsequence, raw_seq + start, width);
    subsequence[width] = '\0';
    subsequences[i_seq] = 
      allocate_seq(get_seq_name(this_seq),
		   get_seq_description(this_seq),
		   get_seq_offset(this_seq), 
		   subsequence);
  }

  // Extract the consensus string in the specified columns.
  char* consensus = get_consensus_string(alignment);
  char* subconsensus = mm_malloc(sizeof(char) * (width + 1));
  strncpy(subconsensus, consensus + start, width);
  subconsensus[width] = '\0';

  // Allocate and return the new alignment.
  ALIGNMENT_T* subalignment 
    = allocate_alignment(get_alignment_name(alignment),
			 get_alignment_description(alignment),
			 num_sequences,
			 subsequences,
			 subconsensus);

  // Free local dynamic memory.
  for (i_seq = 0; i_seq < num_sequences; i_seq++) {
    free_seq(subsequences[i_seq]);
  }
  myfree(subsequences);
  myfree(subsequence);
  return(subalignment);
}

/****************************************************************************
 * Remove from the alignment all columns that contain gaps for the
 * specified species.
 ****************************************************************************/
ALIGNMENT_T* remove_alignment_gaps
  (char*        species,
   ALIGNMENT_T* alignment)
{
  // Locate this species in the alignment.
  int species_index = get_index_in_string_list(species, 
					       get_species_names(alignment));
  if (species_index == -1) {
    die("Can't find %s in alignment.\n", species);
  }
  SEQ_T* this_seq = get_alignment_sequence(species_index, alignment);

  // Get the dimensions of the original matrix.
  int num_sequences = get_num_aligned_sequences(alignment);
  int alignment_length = get_alignment_length(alignment);

  // Allocate memory for raw sequences that will constitute the new alignment.
  char** raw_sequences = (char**)mm_malloc(sizeof(char*) * num_sequences);
  int i_seq = 0;
  for (i_seq = 0; i_seq < num_sequences; i_seq++) {
    raw_sequences[i_seq] 
      = (char*)mm_calloc(alignment_length + 1, sizeof(char*));
  }
  char* consensus = get_consensus_string(alignment);
  char* new_consensus 
    = (char*)mm_calloc(alignment_length + 1, sizeof(char*));

  // Iterate over all columns.
  int i_column;
  int i_raw = 0;
  for (i_column = 0; i_column < alignment_length; i_column++) {

    // Is there a gap?
    char this_char = get_seq_char(i_column, this_seq);
    if ((this_char != '-') && (this_char != '.')) {

      // If no gap, then copy over this column.
      for (i_seq = 0; i_seq < num_sequences; i_seq++) {
	SEQ_T* this_sequence = get_alignment_sequence(i_seq, alignment);
	char this_char = get_seq_char(i_column, this_sequence);
				      
	raw_sequences[i_seq][i_raw] = this_char;
      }
      new_consensus[i_raw] = consensus[i_column];
      i_raw++;
    }
  }

  // Create new sequence objects.
  SEQ_T** new_sequences = (SEQ_T**)mm_malloc(num_sequences * sizeof(SEQ_T*));
  for (i_seq = 0; i_seq < num_sequences; i_seq++) {
    SEQ_T* this_sequence = get_alignment_sequence(i_seq, alignment);
    new_sequences[i_seq] = allocate_seq(get_seq_name(this_sequence),
					get_seq_description(this_sequence),
					get_seq_offset(this_sequence),
					raw_sequences[i_seq]);
  }

  // Allocate and return the new alignment.
  ALIGNMENT_T* new_alignment
    = allocate_alignment(get_alignment_name(alignment),
			 get_alignment_description(alignment),
			 num_sequences,
			 new_sequences,
			 new_consensus);
  
  // Free local dynamic memory.
  for (i_seq = 0; i_seq < num_sequences; i_seq++) {
    myfree(raw_sequences[i_seq]);
    free_seq(new_sequences[i_seq]);
  }
  myfree(raw_sequences);
  myfree(new_sequences);
  myfree(new_consensus);

  return(new_alignment);
}


/****************************************************************************
 * Remove from an alignment any sequence whose ID is not in a given list.
 *
 * N.B. It is NOT an error for the given list to contain sequence IDs that 
 * are not in the alignment.
 ****************************************************************************/
ALIGNMENT_T* remove_alignment_seqs
  (STRING_LIST_T* seqs_to_keep,
   ALIGNMENT_T*   alignment)
{
  // Extract the names of the sequences in the alignment.
  STRING_LIST_T* alignment_species = get_species_names(alignment);
  int num_species = get_num_strings(alignment_species);

  // Count how many sequences will be in the new alignment.
  int i_species;
  int num_final = 0;
  for (i_species = 0; i_species < num_species; i_species++) {
    char* this_species = get_nth_string(i_species, alignment_species);

    if (have_string(this_species, seqs_to_keep)) {
      num_final++;
    } else {
      if (verbosity >= NORMAL_VERBOSE) {
        fprintf(stderr, "Removing %s from alignment.\n", this_species);
      }
    }
  }

  // Allocate space for the new sequences.
  SEQ_T** new_sequences = (SEQ_T**)mm_malloc(num_final * sizeof(SEQ_T*));

  // Copy the sequences.
  int final_index = 0;
  num_species = get_num_strings(seqs_to_keep);
  for (i_species = 0; i_species < num_species; i_species++) {
    char* this_species = get_nth_string(i_species, seqs_to_keep);

    // If the requested ID is in the alignment, then copy over the sequence.
    if (have_string(this_species, alignment_species)) {
      SEQ_T* this_seq 
	= get_alignment_sequence_by_name(this_species, alignment);
      new_sequences[final_index] = copy_sequence(this_seq);
      final_index++;
    }
  }

  // Allocate and return the new alignment.
  
  char *consensus = NULL;
  copy_string(&consensus, get_consensus_string(alignment));
  ALIGNMENT_T* new_alignment
    = allocate_alignment(get_alignment_name(alignment),
			 get_alignment_description(alignment),
			 num_final,
			 new_sequences,
			 consensus);

  return(new_alignment);
}


/****************************************************************************
 * Create a new alignment with any sequence that contains nothing but 
 * gap ('-') characters removed. Returns the new alignment.  Does not 
 * change the old alignment.
 * If there are no all-gap sequences, the returned alignment is the
 * same object as the original alignment.
 ****************************************************************************/
static ALIGNMENT_T* remove_allgap_sequences(ALIGNMENT_T* alignment)
{
  ALIGNMENT_T* new_alignment;
  int i_aln;
  int l_aln = get_num_aligned_sequences(alignment);
  STRING_LIST_T* keeper_seqs = new_string_list();

  // Identify the all-gap sequences.
  for (i_aln=0; i_aln<l_aln; i_aln++) {
    SEQ_T* sequence = get_alignment_sequence(i_aln, alignment);
    int i_seq;
    int l_seq = get_seq_length(sequence);
    // Add sequence to keepers if it contains a non-gap.
    for (i_seq=0; i_seq<l_seq; i_seq++) {
      if (get_seq_char(i_seq, sequence) != '-') {           // not gap?
	add_string(get_seq_name(sequence), keeper_seqs);    // non-gap: keeper
	break;
      }
    }
  }

  // Remove any sequences not in keeper list.
  if (get_num_strings(keeper_seqs) < l_aln) {
    new_alignment = remove_alignment_seqs(keeper_seqs, alignment);
    free_string_list(keeper_seqs);
  } else {
    new_alignment = alignment;
  }

  return(new_alignment);
} // remove_allgap_sequences

/****************************************************************************
 * Read an alignment from a file.  Sort the sequences by sequence name if
 * requested.  Remove all gap sequences if requested.
 ****************************************************************************/
ALIGNMENT_T* read_alignment_from_file
  (char *filename, 
   BOOLEAN_T sort,
   BOOLEAN_T remove_allgap_seqs,
   int* ref_seq_index
  )
{
  int i;

  // Read the sequences.
  ALIGNMENT_T* alignment = read_alignment_from_clustalw_file(filename);

  if (sort) {
    // Create a temporary array to hold sorted sequence pointers.
    int num_sequences = get_num_aligned_sequences(alignment);
    SEQ_T** sequences = (SEQ_T**) mm_malloc(num_sequences * sizeof(SEQ_T*));

    // Sort the sequences by name.
    STRING_LIST_T* alignment_species = get_species_names(alignment);
    // Store the name of the reference sequence.
    char *ref_name = get_nth_string(*ref_seq_index, alignment_species);
    sort_string_list(alignment_species); 	// keep species alphabetical
    for (i=0; i<num_sequences; i++) { 
      char *name = get_nth_string(i, alignment_species);
      sequences[i] = get_alignment_sequence_by_name(name, alignment);
    }
    myfree(alignment->sequences);
    alignment->sequences = sequences;

    // Find the new index of the reference sequence.
    *ref_seq_index = get_index_in_string_list(ref_name, alignment_species);
  }

  if (remove_allgap_seqs) {
    ALIGNMENT_T* new_alignment = remove_allgap_sequences(alignment);
    if (new_alignment != alignment) {
      free_alignment(alignment);
      alignment = new_alignment;
    }
  }

  return(alignment);
} // read_alignment_from_file

/*************************************************************************
 * Print an alignment in PHYLIP format.
 *************************************************************************/
#define OUTPUT_WIDTH 60
void print_phylip_alignment
  (ALIGNMENT_T* the_alignment,
   FILE* outfile)
{
  int i_seq;
  int i_position;
  char buffer[OUTPUT_WIDTH+1];
  char* this_sequence;

  fprintf(outfile, "%d %d\n", the_alignment->num_sequences, 
	  the_alignment->length);

  /* Print the IDs and initial sequences. */
  for (i_seq = 0; i_seq < the_alignment->num_sequences; i_seq++) {
    
    /* Print the ID. */
    strncpy(buffer, get_seq_name(the_alignment->sequences[i_seq]), 10);
    buffer[10] = '\0';
    fprintf(outfile, "%-10s ", buffer);

    /* Print the first block of sequence. */
    this_sequence = get_raw_sequence(the_alignment->sequences[i_seq]);
    strncpy(buffer, &(this_sequence[0]), OUTPUT_WIDTH);
    buffer[OUTPUT_WIDTH] = '\0';
    fprintf(outfile, "%s\n", buffer);
  }

  /* Blank line between sequences. */
  fprintf(outfile, "\n");

  /* Print successive blocks. */
  for (i_position = OUTPUT_WIDTH; i_position < the_alignment->length;
       i_position += OUTPUT_WIDTH) {

    for (i_seq = 0; i_seq < the_alignment->num_sequences; i_seq++) {
      this_sequence = get_raw_sequence(the_alignment->sequences[i_seq]);
      strncpy(buffer, &(this_sequence[i_position]), OUTPUT_WIDTH);
      buffer[OUTPUT_WIDTH] = '\0';
      fprintf(outfile, "           %s\n", buffer);
    }

    /* Blank line between sequences. */
    fprintf(outfile, "\n");
  }
}

/****************************************************************************
 * Free one alignment object.
 ****************************************************************************/
void free_alignment(ALIGNMENT_T* alignment) {

  if (alignment == NULL) {
    return;
  }
  else {
    if (alignment->consensus_string != NULL) {
      myfree(alignment->consensus_string);
    }
    if (alignment->sequences != NULL) {
      int i; 
      for(i = 0; i < alignment->num_sequences; i++) {
        assert(alignment->sequences[i] != NULL);
        free_seq(alignment->sequences[i]);
      }
      myfree(alignment->sequences);
    }
    myfree(alignment);
  }
}
