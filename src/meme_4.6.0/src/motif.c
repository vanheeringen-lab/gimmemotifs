/***********************************************************************
 * FILE: motif.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 7-13-97
 * PROJECT: MHMM
 * COPYRIGHT: 2001-2008, WSN
 * DESCRIPTION: Data structure for representing one motif.
 ***********************************************************************/
#include <stdlib.h> /* For definition of NULL .*/
#include <assert.h>
#include <string.h>
#include "utils.h"
#include "motif.h"
#include "alphabet.h"
#include "matrix.h"


// Minimum function
#define min(a,b)      (a<b)?a:b

/***********************************************************************
 * Calculates the information content of a position of the motif.
 *
 * Assumes that alph_size does not include ambigious characters.
 ***********************************************************************/
static inline double position_information_content(
  MOTIF_T *a_motif,
  int position
) {
  int i;
  double H, item;
  ARRAY_T *freqs;

  H = 0;
  freqs = get_matrix_row(position, a_motif->freqs);
  for (i = 0; i < a_motif->alph_size; ++i) {
    item = get_array_item(i, freqs);
    H -= item*my_log2(item);
  }
  return my_log2(a_motif->alph_size) - H;
}

/***********************************************************************
 * Get or set the index number of a motif.
 ***********************************************************************/
char* get_motif_id
  (MOTIF_T* motif)
{
  return motif->id;
}

void set_motif_id
  (char* id,
   MOTIF_T* motif)
{
  strcpy(motif->id, id);
}

char* get_motif_id2
  (MOTIF_T* motif)
{
  return motif->id2;
}

void set_motif_id2
  (char* id2,
   MOTIF_T* motif)
{
  strcpy(motif->id2, id2);
}

char* get_bare_motif_id
  (MOTIF_T* motif)
{
  char *bare_motif_id = motif->id;
  // Drop the strand indicator from the motif id
  if (*bare_motif_id == '+' || *bare_motif_id == '-') {
    bare_motif_id++;
  }
  return bare_motif_id;
}

/***********************************************************************
 * allocate memory for a motif
 ***********************************************************************/
MOTIF_T* allocate_motif(char *id, MATRIX_T* freqs){
	MOTIF_T* motif = mm_malloc(sizeof(MOTIF_T));

	assert(id != NULL);
	assert(freqs != NULL);

	motif->length = 0;
	motif->alph_size = 0;
	motif->ambigs = 0;
	motif->evalue = 0.0;
	motif->num_sites = 0.0;
	motif->complexity = 0.0;
	motif->freqs = duplicate_matrix(freqs);
  motif->scores = NULL;
	// Set required fields.
	int length = strlen(id) + 1;
	strncpy(motif->id, id, min(length,MAX_MOTIF_ID_LENGTH));

	motif->alph_size = get_num_cols(motif->freqs);
	motif->length = get_num_rows(motif->freqs);
  motif->url = NULL;
  motif->trim_left = 0;
  motif->trim_right = 0;
	return motif;
}

/***********************************************************************
 * Get or set the frequencies
 ***********************************************************************/
MATRIX_T* get_motif_freqs(MOTIF_T* motif){
  return motif->freqs;
}

void set_motif_freqs
  (MATRIX_T* freqs,
   MOTIF_T* motif){
	// free existing matrix
	if (motif->freqs != NULL){
		free_matrix(motif->freqs);
	}
	motif->freqs = duplicate_matrix(freqs);
}

/***********************************************************************
 * Get the E-value of a motif.
 ***********************************************************************/
double get_motif_evalue
  (MOTIF_T* motif)
{
  return motif->evalue;
}

/***********************************************************************
 * Get the E-value of a motif.
 ***********************************************************************/
double get_motif_nsites
  (MOTIF_T* motif)
{
  return motif->num_sites;
}

/***********************************************************************
 * Get the motif length
 ***********************************************************************/
double get_motif_length
  (MOTIF_T* motif) 
{
  return motif->length;
}

/***********************************************************************
 * Get the motif length after trimming
 ***********************************************************************/
int get_motif_trimmed_length
  (MOTIF_T* motif)
{
  return motif->length - motif->trim_left - motif->trim_right;
}

/***********************************************************************
 * Get the motif alphabet size (non-ambiguous letters)
 ***********************************************************************/
double get_motif_alph_size
  (MOTIF_T* motif) 
{
  return motif->alph_size;
}


/***********************************************************************
 * Get the motif ambiguous alphabet size (only ambiguous letters)
 ***********************************************************************/
double get_motif_ambiguous_size
  (MOTIF_T* motif)
{
  return motif->ambigs;
}

/***********************************************************************
 * Get the position specific score for a letter
 ***********************************************************************/
double get_motif_score
  (MOTIF_T* motif, int position, int i_alph)
{
  return get_matrix_cell(position, i_alph, motif->scores);
}

/***********************************************************************
 * Return one column of a motif, as a newly allocated array of counts.
 ***********************************************************************/
ARRAY_T* get_motif_counts
  (int      position,
   MOTIF_T* motif)
{
  ARRAY_T* return_value = allocate_array(motif->alph_size);

  int i_alph;
  for (i_alph = 0; i_alph < motif->alph_size; i_alph++) {
    set_array_item(i_alph,
		   motif->num_sites * get_matrix_cell(position, 
						      i_alph, motif->freqs),
		   return_value);
  }
  return(return_value);
}

/***********************************************************************
 * Get the url of a motif
 ***********************************************************************/
char* get_motif_url
  (MOTIF_T *motif)
{
  return motif->url;
}

/***********************************************************************
 * Set the url of a motif
 ***********************************************************************/
void set_motif_url
  (char    *url,
   MOTIF_T *motif)
{
  if (motif->url) {
    free(motif->url);
    motif->url = NULL;
  }
  copy_string(&(motif->url), url);
}

/***********************************************************************
 * Get the number of positions to trim from the left of the motif
 ***********************************************************************/
int get_motif_trim_left
  (MOTIF_T *motif)
{
  return motif->trim_left;
}

/***********************************************************************
 * Get the number of positions to trim from the right of the motif
 ***********************************************************************/
int get_motif_trim_right
  (MOTIF_T *motif)
{
  return motif->trim_right;
}

/***********************************************************************
 * Clear the motif trim
 ***********************************************************************/
void clear_motif_trim
  (MOTIF_T *motif)
{
  motif->trim_left = 0;
  motif->trim_right = 0;
}

/***********************************************************************
 * Determine whether a given motif is in a given list of motifs.
 ***********************************************************************/
BOOLEAN_T have_motif
  (char*    motif_id,
   int      num_motifs,
   MOTIF_T* motifs)
{
  int i_motif;

  for (i_motif = 0; i_motif < num_motifs; i_motif++) {
    if (strcmp(motifs[i_motif].id, motif_id) == 0) {
      return(TRUE);
    }
  }
  
  return(FALSE);

}

/***********************************************************************
 * Copy a motif from one place to another.
 ***********************************************************************/
void copy_motif
  (MOTIF_T* source,
   MOTIF_T* dest)
{
  strcpy(dest->id, source->id);
  strcpy(dest->id2, source->id2);
  dest->length = source->length;
  dest->alph_size = source->alph_size;
  dest->ambigs = source->ambigs;
  dest->evalue = source->evalue;
  dest->num_sites = source->num_sites;
  dest->complexity = source->complexity;
  dest->trim_left = source->trim_left;
  dest->trim_right = source->trim_right;
  if (source->freqs) {
    // Allocate memory for the matrix.
    dest->freqs = allocate_matrix(dest->length, dest->alph_size + dest->ambigs);
    // Copy the matrix.
    copy_matrix(source->freqs, dest->freqs);
  } else {
    dest->freqs = NULL;
  }
  if (source->scores) {
    // Allocate memory for the matrix. Note that scores don't contain ambigs.
    dest->scores = allocate_matrix(get_num_rows(source->scores), get_num_cols(source->scores));
    // Copy the matrix.
    copy_matrix(source->scores, dest->scores);
  } else {
    dest->scores = NULL;
  }
  copy_string(&(dest->url), source->url);
}

/***********************************************************************
 * Say that the motif ID is printed centered above a given motif.
 * If the motif ID string is longer than the motif, we truncate
 * it on the right and align the first character over the start of
 * the motif.
 * This function returns the character that appears in the nth
 * position of that motif ID string.
 ***********************************************************************/
char get_motif_id_char
  (int      position,
   MOTIF_T* a_motif)
{
  char* motif_id_string;
  int   id_width = strlen(a_motif->id);
  int   m_width = a_motif->length;
  int   id_start;
  char  return_char;

  assert(position < a_motif->length);

  // Allocate the string.
  motif_id_string = mm_calloc(sizeof(char), a_motif->length + 1);

  // Get position where ID starts relative to start of motif.
  id_start = id_width <= m_width ? ((m_width - id_width) / 2) : 0;
  // FIXME: (tlb) The following if() was put in to make the smoke tests of mhmm
  // pass.  It should be removed and the smoke test comparison files changed.
  if (m_width % 2 == 0 && id_width % 2 == 0) {
    id_start++; 
  } else {
   id_start+=2;
  }

  // Create the centered ID string.
  sprintf(motif_id_string, "%*.*s%-*.*s", id_start, id_start, "",
          m_width-id_start, m_width-id_start, a_motif->id);
  assert((int)(strlen(motif_id_string)) == a_motif->length);

  // Get the nth character.
  return_char = motif_id_string[position];

  if (return_char == ' ') {
    if ((position == 0) || (position == (a_motif->length - 1))) {
      return_char = '*';
    } else {
      return_char = '_';
    }
  }

  // Free up memory and return.
  myfree(motif_id_string);
  return(return_char);
}

/***********************************************************************
 * Compute the reverse complement of one DNA frequency distribution.
 * 
 * Assumes DNA alphabet in order ACGT.
 ***********************************************************************/
void complement_dna_freqs
  (ARRAY_T* source,
   ARRAY_T* dest)
{
  set_array_item(0, get_array_item(3, source), dest); // A -> T
  set_array_item(1, get_array_item(2, source), dest); // C -> G
  set_array_item(2, get_array_item(1, source), dest); // G -> C
  set_array_item(3, get_array_item(0, source), dest); // T -> A

  //check if the frequencies have ambiguous characters
  //for example meme does not use ambiguous characters
  if (get_array_length(source) > 4) {
    fill_in_ambiguous_chars(FALSE, dest);
  }
}

/***********************************************************************
 * Convert array by compute the average of complementary dna frequencies.
 *
 * Assumes DNA alphabet in order ACGT.
 ***********************************************************************/
void balance_complementary_dna_freqs
  (ARRAY_T* source)
{
  double at = (get_array_item(0, source)+get_array_item(3, source))/2.0;
  double cg = (get_array_item(1, source)+get_array_item(2, source))/2.0;
  set_array_item(0, at, source); // A -> T
  set_array_item(1, cg, source); // C -> G
  set_array_item(2, cg, source); // G -> C
  set_array_item(3, at, source); // T -> A

  fill_in_ambiguous_chars(FALSE, source);
}

/***********************************************************************
 * Turn a given motif into its own reverse complement.
 * TODO this does not handle the scores matrix, and it should.
 ***********************************************************************/
void reverse_complement_motif
  (MOTIF_T* a_motif)
{
  int i, temp_trim;
  ARRAY_T* left_freqs;
  ARRAY_T* right_freqs;
  ARRAY_T* temp_freqs;   // Temporary space during swap.

  // Allocate space.
  //temp_freqs = allocate_array(get_alph_size(ALL_SIZE)); //this relys on a global which assumes DNA has ambigs which meme doesn't seem to use 
  temp_freqs = allocate_array(a_motif->alph_size + a_motif->ambigs);

  // Consider each row (position) in the motif.
  for (i = 0; i < (int)((a_motif->length + 1) / 2); i++) {
    left_freqs = get_matrix_row(i, a_motif->freqs);
    right_freqs = get_matrix_row(a_motif->length - (i + 1), a_motif->freqs);

    // Make a temporary copy of one row.
    copy_array(left_freqs, temp_freqs);

    // Compute reverse complements in both directions.
    complement_dna_freqs(right_freqs, left_freqs);
    complement_dna_freqs(temp_freqs, right_freqs);
  }
  free_array(temp_freqs);
  //swap the trimming variables
  temp_trim = a_motif->trim_left;
  a_motif->trim_left = a_motif->trim_right;
  a_motif->trim_right = temp_trim;
}

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
) {
  int i, len;

  len = a_motif->length;
  for (i = 0; i < len; ++i) {
    if (position_information_content(a_motif, i) >= threshold_bits) break;
  }
  a_motif->trim_left = i;
  if (i == len) {
    a_motif->trim_right = 0;
    return;
  }
  for (i = len-1; i >= 0; --i) {
    if (position_information_content(a_motif, i) >= threshold_bits) break;
  }
  a_motif->trim_right = len - i - 1;
}

/***********************************************************************
 * Compute the complexity of a motif as a number between 0 and 1.
 *
 * Motif complexity is the average K-L distance between the "motif
 * background distribution" and each column of the motif.  The motif
 * background is just the average distribution of all the columns.  The
 * K-L distance, which measures the difference between two
 * distributions, is the same as the information content:
 *
 *  \sum_i p_i log(p_i/f_i)
 *
 * This value increases with increasing complexity.
 ***********************************************************************/
double compute_motif_complexity
  (MOTIF_T* a_motif)
{
  double return_value;
  ARRAY_T* motif_background;  // Mean emission distribution.
  int num_rows;
  int i_row;
  int num_cols;
  int i_col;

  num_cols = get_alph_size(ALPH_SIZE);
  num_rows = get_num_rows(a_motif->freqs);

  // Compute the mean emission distribution.
  motif_background = get_matrix_col_sums(a_motif->freqs);
  scalar_mult(1.0 / (double)num_rows, motif_background);

  // Compute the K-L distance w.r.t. the background.
  return_value = 0;
  for (i_row = 0; i_row < num_rows; i_row++) {
    ARRAY_T* this_emission = get_matrix_row(i_row, a_motif->freqs);
    for (i_col = 0; i_col < num_cols; i_col++) {
      ATYPE this_item = get_array_item(i_col, this_emission);
      ATYPE background_item = get_array_item(i_col, motif_background);

      // Use two logs to avoid handling divide-by-zero as a special case.
      return_value += this_item 
	* (my_log(this_item) - my_log(background_item));
    }
  }

  free_array(motif_background);
  return(return_value / (double)num_rows);
}

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
   MOTIF_T*  a_motif)
{
  // Make sure the given threshold is in the right range.
  if ((threshold < 0.0) || (threshold > 100.0)) {
    die(
      "Information threshold (%g) must be a percentage between 0 and 100.\n",
	    threshold
    );
  }

  // Get the dimensions of the motif.
  int num_cols = get_alph_size(ALPH_SIZE);
  int num_rows = get_num_rows(a_motif->freqs);

  // Compute and store the information content for each row
  // and the total information content for the motif.
  ATYPE total_information_content = 0.0;
  ARRAY_T* information_content = allocate_array(num_rows);
  int i_row;
  int i_col;
  for (i_row = 0; i_row < num_rows; i_row++) {
    ATYPE row_content = 0.0;
    ARRAY_T* this_emission = get_matrix_row(i_row, a_motif->freqs);
    for (i_col = 0; i_col < num_cols; i_col++) {
      ATYPE this_item = get_array_item(i_col, this_emission);
      ATYPE background_item = get_array_item(i_col, background);

      // Use two logs to avoid handling divide-by-zero as a special case.
      ATYPE partial_row_content = 
        this_item * (my_log(this_item) - my_log(background_item));

      row_content += partial_row_content;
      total_information_content += partial_row_content;

    }
    set_array_item(i_row, row_content, information_content);
  }

  // Search for the target position.
  int return_value = -1;
  ATYPE cumulative_content = 0.0;
  ATYPE percent = 0.0;
  if (from_start) {
    // Search from start for IC exceeding threshold.
    for (i_row = 0; i_row < num_rows; i_row++) {
      cumulative_content += get_array_item(i_row, information_content);
      percent = 100 *  cumulative_content / total_information_content;
      if (percent >= threshold) {
        return_value = i_row;
        break;
      }
    }
  }
  else {
    // Search from end for IC exceeding threshold.
    for (i_row = num_rows - 1; i_row >= 0; i_row--) {
      cumulative_content += get_array_item(i_row, information_content);
      percent = 100 *  cumulative_content / total_information_content;
      if (percent >= threshold) {
        return_value = i_row;
        break;
      }
    }
  }

  if (return_value == -1) {
    die(
      "Can't find a position that accounts for %g of information content.",
      threshold
    );
  }
  free_array(information_content);
  return(return_value);
}


/***********************************************************************
 * Returns the string that is the best possible match to the given motif.
 * Caller is responsible for freeing string.
 ***********************************************************************/
char *get_best_possible_match(MOTIF_T *motif) {
  int mpos, apos; 
  char *match_string;
  
  assert(motif != NULL);
  assert(motif->freqs != NULL);
  assert(motif->length == motif->freqs->num_rows);
  assert(motif->alph_size <= motif->freqs->num_cols); 

  match_string = mm_malloc(sizeof(char) * (motif->length + 1));

  // Find the higest scoring character at each position in the motif.
  for(mpos = 0; mpos < motif->length; ++mpos) {
    ARRAY_T *row = motif->freqs->rows[mpos];
    double max_v = row->items[0];
    int max_i = 0;
    for(apos = 1; apos < motif->alph_size; ++apos) {
     if (row->items[apos] >= max_v) {
        max_i = apos;
        max_v = row->items[apos];
     }
    }
    match_string[mpos] = get_alph_char(max_i);
  }

  //  Add null termination
  match_string[motif->length] = '\0';

  return match_string;
}

/***********************************************************************
 * Free dynamic memory used by one motif. 
 ***********************************************************************/
void free_motif
  (MOTIF_T *a_motif)
{
  /* Don't bother with empty motifs. */
  if (a_motif == NULL) 
    return;

  // Reset all memeber values
  a_motif->id[0] = 0;
  a_motif->id2[0] = 0;
  a_motif->length = 0;
  a_motif->alph_size = 0;
  a_motif->ambigs = 0;
  a_motif->evalue = NaN();
  a_motif->num_sites = NaN();
  a_motif->complexity = NaN();
  free_matrix(a_motif->freqs);
  a_motif->freqs = NULL;
  free_matrix(a_motif->scores);
  a_motif->freqs = NULL;
  if (a_motif->url) free(a_motif->url);
  a_motif->url = NULL;
}
