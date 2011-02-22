/********************************************************************
 * FILE: cisml.c
 * AUTHOR: William Stafford Noble, Charles E. Grant
 * CREATE DATE: 9/25/2007
 * PROJECT: MEME suite
 * COPYRIGHT: 2007, UW
 ********************************************************************/

#include <assert.h>
#include <float.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include "array.h"
#include "cisml.h"
#include "dir.h"
#include "hash_alph.h"
#include "fasta-io.h"
#include "heap.h"
#include "io.h"
#include "qvalue.h"
#include "xml-util.h"
#include "cisml-dtd.h"

const int PATTERN_INCREMENT = 5;
const int SEQUENCE_INCREMENT = 50;
const int ELEMENT_INCREMENT = 500;

// Define data structures related to CisML format.

struct multi_pattern {
  double *score; // May be NULL
  double *pvalue; // May be NULL
  int num_patterns; // Count of child patterns in this multi-pattern
  int num_allocated_patterns; // Size of pattern array
  PATTERN_T **patterns; // Array of child patterns in this multi-pattern
};

struct pattern {

  char *accession; // Required.
  char *name; // Required.

  double *pvalue; // May be NULL.
  double *score; // May be NULL.
  char *db; // May be NULL.
  char *lsid; // May be NULL.

  int num_allocated_sequences; // Size of child scanned-sequence array.
  int num_allocated_elements; // Size of matched-element array.
  int num_sequences; // Count of child scanned-sequences.
  long num_scanned_positions; // Count of all positions scanned for
                            // matched_element.
  int num_stored_matches; // Count of the matched-elements stored in the elements array.
  int max_stored_matches;   // maximum number of matches to store for this pattern.
  double max_pvalue_retained; // Largest pvalue of retained records.

  BOOLEAN_T has_all_pvalues;  // Has retained all matched elements, not just those
                              // with best p-values.
  BOOLEAN_T qvalues_computed; // Have q-values been calcuatled for matched-elements
  BOOLEAN_T is_complete;     /// All matched elements have been added to pattern
  SCANNED_SEQUENCE_T **sequences; // Array of child scanned-sequence pointers.
  HEAP *element_heap; // Heap of matched elements ordered by ascending score.
  MATCHED_ELEMENT_T **elements; // Array of matched element pointers, ordered by p-value
};

struct scanned_sequence {

  char* accession; // Required.
  char* name; // Required.

  double *pvalue; // May be NULL.
  double *score; // May be NULL.
  int *length; // May be NULL.
  char* db; // May be NULL.
  char* lsid; // May be NULL.

  long num_scanned_positions; // Count of all positions scanned for
                            // matched_elements.
  int num_matched_elements; // Count of all matched_elements
                            // <= num_scanned_elements because of filtering.
  int num_allocated_elements; // Number of elements in elements array.

  MATCHED_ELEMENT_T **elements; // Array of all matched elements for this sequence
  PATTERN_T *parent_pattern; // Pointer to containing pattern.

};

struct matched_element {

  int start; // Required.
  int stop; // Required.

  double score;
  BOOLEAN_T has_score;
  double pvalue;
  BOOLEAN_T has_pvalue;
  double qvalue;
  BOOLEAN_T has_qvalue;
  char* clusterid; // May be NULL.
  char* sequence; // May be NULL.
  char strand;

  SCANNED_SEQUENCE_T *parent_sequence; // Pointer to containing scanned-sequence.

};

struct cisml {

  char *program_name; // Required
  char *pattern_file; // Required
  char *sequence_file; // Required

  char *background_file; // May return NULL
  double *pattern_pvalue_cutoff; // May return NULL
  double *sequence_pvalue_cutoff; // May return NULL
  double *site_pvalue_cutoff; // May return NULL
  double *site_qvalue_cutoff; // May return NULL
  char* sequence_filter; // May return NULL

  int num_allocated_multi_patterns; // Size of multi-pattern array.
  int num_allocated_patterns; //Size of patern array.
  int num_multi_patterns; // Count of multi-pattern objects.
  int num_patterns; // Count of pattern objects.

  MULTI_PATTERN_T **multi_patterns; // Array of pointers to multi-pattern objects.
  PATTERN_T **patterns; // Array of pointers to pattern objects.

};

/**********************************************************************
  compare_matched_elements

  Compare two objects, treating them as matched elements
**********************************************************************/
int compare_matched_elements(void *p1, void *p2) {

  MATCHED_ELEMENT_T *e1 = (MATCHED_ELEMENT_T *) p1;
  MATCHED_ELEMENT_T *e2 = (MATCHED_ELEMENT_T *) p2;

  if (e1->pvalue < e2->pvalue) {
    return 1;
  }
  else if (e1->pvalue > e2->pvalue) {
    return -1;
  }
  else {
    // If p-values are equal, compare staring postions
    // to break the time.
    return e1->start < e2->start ? 1: -1;
  }

}

/**********************************************************************
  destroy_matched_element

  Destroy an object, treating it as a matched element
  For use with HEAP.
**********************************************************************/
void destroy_matched_element(void *p) {

  MATCHED_ELEMENT_T *e = (MATCHED_ELEMENT_T *) p;
  free_matched_element(e);

}

/**********************************************************************
  copy_matched_element

  Copy an object, treating it as a matched element
  For use with HEAP.
**********************************************************************/
void *copy_matched_element(void *p) {

  MATCHED_ELEMENT_T *e = (MATCHED_ELEMENT_T *) p;
  MATCHED_ELEMENT_T *new_e = allocate_matched_element_with_score(
      e->start,
      e->stop,
      e->score,
      e->pvalue,
      e->parent_sequence
    );
  new_e->qvalue = e->qvalue;
  new_e->clusterid = e->clusterid;
  new_e->sequence = e->sequence;
  new_e->strand = e->strand;
  return new_e;

}

/**********************************************************************
  reduce_pattern_matched_elements

  To conserve memory, remove the least significant matched_elements.
  We want to remove at least PERCENT_ELEMENT_DISCARD * num_matched_elements,
  but we may remove more if there remain matched elments with a
  p-value equal to the delected elements.

**********************************************************************/
static void reduce_pattern_matched_elements(PATTERN_T *pattern);

/*************************************************************************
 * sort_matched_elements
 *
 * Sort a an array of pointers to matched-elements sites by pvalue
 * or by sequence name and position.
 *************************************************************************/
void sort_matched_elements(
  BOOLEAN_T sort_by_pvalue,
  int num_elements,
  MATCHED_ELEMENT_T **elements
);


/**********************************************************************
  allocate_cisml

  Constructor for the cisml data structure.
  Sets required fields to point to copies of the provided arguments.
  Other fields set to 0 or NULL.
**********************************************************************/
CISML_T *allocate_cisml(
  char *program_name,
  char *pattern_file,
  char *sequence_file
) {

  assert(program_name != NULL);
  assert(pattern_file != NULL);
  assert(sequence_file != NULL);

  // Allocate memory and initialze fields
  CISML_T* cisml = mm_malloc(sizeof(CISML_T));
  cisml->program_name = NULL;
  cisml->pattern_file = NULL;
  cisml->sequence_file = NULL;
  cisml->background_file = NULL;
  cisml->pattern_pvalue_cutoff = NULL;
  cisml->sequence_pvalue_cutoff = NULL;
  cisml->site_pvalue_cutoff = NULL;
  cisml->site_qvalue_cutoff = NULL;
  cisml->sequence_filter = NULL;
  cisml->num_multi_patterns = 0;
  cisml->num_allocated_multi_patterns = 0;
  cisml->multi_patterns = NULL;
  cisml->num_patterns = 0;
  cisml->num_allocated_patterns = 0;
  cisml->patterns = NULL;

  // Set required fields.
  int length = strlen(program_name) + 1;
  cisml->program_name = mm_malloc(length * sizeof(char));
  strncpy(cisml->program_name, program_name, length);

  length = strlen(pattern_file) + 1;
  cisml->pattern_file = mm_malloc(length * sizeof(char));
  strncpy(cisml->pattern_file, pattern_file, length);

  length = strlen(sequence_file) + 1;
  cisml->sequence_file = mm_malloc(length * sizeof(char));
  strncpy(cisml->sequence_file, sequence_file, length);

  return cisml;
}

/**********************************************************************
  free_cisml

  Destructor for the cisml data structure.
**********************************************************************/
void free_cisml(CISML_T *cisml) {

  assert(cisml != NULL);

  while(cisml->num_multi_patterns > 0) {
    free_multi_pattern(cisml->multi_patterns[--cisml->num_multi_patterns]);
  }
  myfree(cisml->multi_patterns);
  cisml->num_allocated_multi_patterns = 0;

  while(cisml->num_patterns > 0) {
    free_pattern(cisml->patterns[--cisml->num_patterns]);
  }
  myfree(cisml->patterns);
  cisml->num_allocated_patterns = 0;

  myfree(cisml->sequence_filter);
  myfree(cisml->site_pvalue_cutoff);
  myfree(cisml->site_qvalue_cutoff);
  myfree(cisml->sequence_pvalue_cutoff);
  myfree(cisml->pattern_pvalue_cutoff);
  myfree(cisml->background_file);
  myfree(cisml->sequence_file);
  myfree(cisml->pattern_file);
  myfree(cisml->program_name);

  myfree(cisml);

}

/**********************************************************************
  get_cisml_program_name

  Gets the program_name member from a cisml object.
**********************************************************************/
char *get_cisml_program_name(CISML_T *cisml) {
  assert(cisml != NULL);
  return cisml->program_name;
}

/**********************************************************************
  get_cisml_pattern_file

  Gets the pattern_file member in a cisml object.
**********************************************************************/
char *get_cisml_pattern_file(CISML_T *cisml) {
  assert(cisml != NULL);
  return cisml->pattern_file;
}

/**********************************************************************
  get_cisml_sequence_file

  Gets the sequence_file member from a cisml object.
**********************************************************************/
char *get_cisml_sequence_file(CISML_T *cisml) {
  assert(cisml != NULL);
  return cisml->sequence_file;
}

/**********************************************************************
  set_cisml_background_file

  Sets the background_file member in a cisml object.
**********************************************************************/
void set_cisml_background_file(CISML_T *cisml, char *background_file) {
  assert(cisml != NULL);
  if (background_file == NULL) {
    if (cisml->background_file != NULL) {
      myfree(cisml->background_file);
    }
    cisml->background_file = NULL;
  }
  else {
    int new_length = strlen(background_file) + 1;
    int old_length = 0;
    if (cisml->background_file != NULL) {
      old_length = strlen(cisml->background_file) + 1;
    }
    if (old_length < new_length) {
      cisml->background_file = realloc(cisml->background_file, new_length);
    }
    strncpy(cisml->background_file, background_file, new_length);
  }
}

/**********************************************************************
  get_cisml_background_file

  Gets the background_file member from a cisml object.
  Return value may be NULL.
**********************************************************************/
char *get_cisml_background_file(CISML_T *cisml) {
  assert(cisml != NULL);
  return cisml->background_file;
}

/**********************************************************************
  set_cisml_pattern_pvalue_cutoff

  Sets the pattern_pvalue_cutoff member in a cisml object.
**********************************************************************/
void set_cisml_pattern_pvalue_cutoff(
  CISML_T *cisml,
  double pattern_pvalue_cutoff
) {
  assert(cisml != NULL);
  if (cisml->pattern_pvalue_cutoff == NULL) {
    cisml->pattern_pvalue_cutoff = mm_malloc(sizeof(double));
  }
  *(cisml->pattern_pvalue_cutoff) = pattern_pvalue_cutoff;
}

/**********************************************************************
  clear_cisml_pattern_pvalue_cutoff

  Sets the pattern_pvalue_cutoff member in a cisml object to null.
**********************************************************************/
void clear_cisml_pattern_pvalue_cutoff(CISML_T *cisml) {
  assert(cisml != NULL);
  if (cisml->pattern_pvalue_cutoff != NULL) {
    myfree(cisml->pattern_pvalue_cutoff);
  }
  cisml->pattern_pvalue_cutoff = NULL;
}

/**********************************************************************
  has_cisml_pattern_pvalue_cutoff

  Does a cisml object have a pattern_pvalue_cutoff?
**********************************************************************/
BOOLEAN_T has_cisml_pattern_pvalue_cutoff(CISML_T *cisml) {
  return cisml->pattern_pvalue_cutoff != NULL ? TRUE : FALSE;

}

/**********************************************************************
  get_cisml_pattern_pvalue_cutoff

  Gets the pattern_pvalue_cutoff member from a cisml object.
  Defaults to 1.0 if the pattern p-value cutoff has not been set.
**********************************************************************/
double get_cisml_pattern_pvalue_cutoff(CISML_T *cisml) {
  assert(cisml != NULL);
  if (cisml->pattern_pvalue_cutoff) {
    return *(cisml->pattern_pvalue_cutoff);
  }
  else {
    return 1.0;
  }
}

/**********************************************************************
  set_cisml_sequence_pvalue_cutoff

  Sets the sequence_pvalue_cutoff member in a cisml object.
**********************************************************************/
void set_cisml_sequence_pvalue_cutoff(
  CISML_T *cisml,
  double sequence_pvalue_cutoff
) {
  assert(cisml != NULL);
  if (cisml->sequence_pvalue_cutoff == NULL) {
    cisml->sequence_pvalue_cutoff = mm_malloc(sizeof(double));
  }
  *(cisml->sequence_pvalue_cutoff) = sequence_pvalue_cutoff;
}

/**********************************************************************
  clear_cisml_sequence_pvalue_cutoff

  Sets the sequence_pvalue_cutoff member in a cisml object to null.
**********************************************************************/
void clear_cisml_sequence_pvalue_cutoff(CISML_T *cisml) {
  assert(cisml != NULL);
  if (cisml->sequence_pvalue_cutoff != NULL) {
    myfree(cisml->sequence_pvalue_cutoff);
  }
  cisml->sequence_pvalue_cutoff = NULL;
}

/**********************************************************************
  has_cisml_sequence_pvalue_cutoff

  Does a cisml object have a sequence_pvalue_cutoff?
**********************************************************************/
BOOLEAN_T has_cisml_sequence_pvalue_cutoff(CISML_T *cisml) {
  return cisml->sequence_pvalue_cutoff != NULL ? TRUE : FALSE;
}

/**********************************************************************
  get_cisml_sequence_pvalue_cutoff

  Gets the sequence_pvalue_cutoff member from a cisml object.
  Defaults to 1.0 if the sequence p-value cutoff has not been set.
**********************************************************************/
double get_cisml_sequence_pvalue_cutoff(CISML_T *cisml) {
  assert(cisml != NULL);
  if (cisml->sequence_pvalue_cutoff) {
    return *(cisml->sequence_pvalue_cutoff);
  }
  else {
    return 1.0;
  }
}

/**********************************************************************
  set_cisml_site_pvalue_cutoff

  Sets the site_pvalue_cutoff member in a cisml object.
**********************************************************************/
void set_cisml_site_pvalue_cutoff(CISML_T *cisml, double site_pvalue_cutoff) {
  assert(cisml != NULL);
  if (cisml->site_pvalue_cutoff == NULL) {
    cisml->site_pvalue_cutoff = mm_malloc(sizeof(double));
  }
  *(cisml->site_pvalue_cutoff) = site_pvalue_cutoff;
}

/**********************************************************************
  clear_cisml_site_pvalue_cutoff

  Sets the site_pvalue_cutoff member in a cisml object to null.
**********************************************************************/
void clear_cisml_site_pvalue_cutoff(CISML_T *cisml) {
  assert(cisml != NULL);
  if (cisml->site_pvalue_cutoff != NULL) {
    myfree(cisml->site_pvalue_cutoff);
  }
  cisml->site_pvalue_cutoff = NULL;
}

/**********************************************************************
  has_cisml_site_pvalue_cutoff

  Does a cisml object have a site_pvalue_cutoff?
**********************************************************************/
BOOLEAN_T has_cisml_site_pvalue_cutoff(CISML_T *cisml) {
  return cisml->site_pvalue_cutoff != NULL ? TRUE : FALSE;
}

/**********************************************************************
  get_cisml_site_pvalue_cutoff

  Gets the site_pvalue_cutoff member from a cisml object.
  Defaults to 1.0 if the site p-value cutoff has not been set.
**********************************************************************/
double get_cisml_site_pvalue_cutoff(CISML_T *cisml) {
  assert(cisml != NULL);
  if (cisml->site_pvalue_cutoff) {
    return *(cisml->site_pvalue_cutoff);
  }
  else {
    return 1.0;
  }
}

/**********************************************************************
  set_cisml_site_qvalue_cutoff

  Sets the site_qvalue_cutoff member in a cisml object.
**********************************************************************/
void set_cisml_site_qvalue_cutoff(CISML_T *cisml, double site_qvalue_cutoff) {
  assert(cisml != NULL);
  if (cisml->site_qvalue_cutoff == NULL) {
    cisml->site_qvalue_cutoff = mm_malloc(sizeof(double));
  }
  *(cisml->site_qvalue_cutoff) = site_qvalue_cutoff;
}

/**********************************************************************
  clear_cisml_site_qvalue_cutoff

  Sets the site_qvalue_cutoff member in a cisml object to null.
**********************************************************************/
void clear_cisml_site_qvalue_cutoff(CISML_T *cisml) {
  assert(cisml != NULL);
  if (cisml->site_qvalue_cutoff != NULL) {
    myfree(cisml->site_qvalue_cutoff);
  }
  cisml->site_qvalue_cutoff = NULL;
}

/**********************************************************************
  has_cisml_site_qvalue_cutoff

  Does a cisml object have a site_qvalue_cutoff?
**********************************************************************/
BOOLEAN_T has_cisml_site_qvalue_cutoff(CISML_T *cisml) {
  return cisml->site_qvalue_cutoff != NULL ? TRUE : FALSE;
}

/**********************************************************************
  get_cisml_site_qvalue_cutoff

  Gets the site_qvalue_cutoff member from a cisml object.
  Defaults to 1.0 if the site q-value cutoff has not been set.
**********************************************************************/
double get_cisml_site_qvalue_cutoff(CISML_T *cisml) {
  assert(cisml != NULL);
  if (cisml->site_qvalue_cutoff) {
    return *(cisml->site_qvalue_cutoff);
  }
  else {
    return 1.0;
  }
}

/**********************************************************************
  get_cisml_sequence_filter

  Gets the sequence_filter member from a cisml object.
  May return NULL.
**********************************************************************/
char *get_cisml_sequence_filter(CISML_T *cisml) {
  assert(cisml != NULL);
  return cisml->sequence_filter;
}

/**********************************************************************
  set_cisml_sequence_filter

  Sets the sequence_filter member in a cisml object.
**********************************************************************/
void set_cisml_sequence_filter(CISML_T *cisml, char *sequence_filter) {
  assert(cisml != NULL);
  if (sequence_filter == NULL) {
    if (cisml->sequence_filter != NULL) {
      myfree(cisml->sequence_filter);
    }
    cisml->sequence_filter = NULL;
  }
  else {
    int new_length = strlen(sequence_filter) + 1;
    int old_length = 0;
    if (cisml->sequence_filter != NULL) {
      old_length = strlen(cisml->sequence_filter) + 1;
    }
    if (old_length < new_length) {
      cisml->sequence_filter = realloc(cisml->sequence_filter, new_length);
    }
    strncpy(cisml->sequence_filter, sequence_filter, new_length);
  }
}

/**********************************************************************
  get_cisml_num_multi_patterns

  Gets the number of multi-patterns from a cisml object.
**********************************************************************/
int get_cisml_num_multi_patterns(CISML_T *cisml) {
  assert(cisml != NULL);
  return cisml->num_multi_patterns;
}

/**********************************************************************
  get_cisml_multi_patterns

  Gets the array of pointers to multi_patterns from a cisml object.
  May return NULL.
**********************************************************************/
MULTI_PATTERN_T **get_cisml_multi_patterns(CISML_T *cisml) {
  assert(cisml != NULL);
  return cisml->multi_patterns;
}

/**********************************************************************
  get_cisml_num_patterns

  Gets the number of patterns from a cisml object.
**********************************************************************/
int get_cisml_num_patterns(CISML_T *cisml) {
  assert(cisml != NULL);
  return cisml->num_patterns;
}

/**********************************************************************
  get_cisml_patterns

  Gets the array of pointers to patterns from a cisml object.
  May return NULL.
**********************************************************************/
PATTERN_T **get_cisml_patterns(CISML_T *cisml) {
  assert(cisml != NULL);
  return cisml->patterns;
}

/**********************************************************************
  add_cisml_multi_pattern

  Adds a pattern to the array of pointers to patterns in a cisml object.
**********************************************************************/
void add_cisml_multi_pattern(CISML_T *cisml, MULTI_PATTERN_T* multi_pattern) {
  assert(cisml != NULL);
  if (multi_pattern == NULL) {
    return;
  }
  assert(cisml->num_multi_patterns <= cisml->num_allocated_multi_patterns);
  if (cisml->num_multi_patterns == cisml->num_allocated_multi_patterns) {
    cisml->num_allocated_multi_patterns += PATTERN_INCREMENT;
    cisml->multi_patterns = mm_realloc(
      cisml->multi_patterns,
      cisml->num_allocated_multi_patterns * sizeof(MULTI_PATTERN_T *)
    );
  }
  cisml->multi_patterns[cisml->num_multi_patterns] = multi_pattern;
  cisml->num_multi_patterns++;
}

/**********************************************************************
  add_cisml_pattern

  Adds a pattern to the array of pointers to patterns in a cisml object.
**********************************************************************/
void add_cisml_pattern(CISML_T *cisml, PATTERN_T* pattern) {
  assert(cisml != NULL);
  if (pattern == NULL) {
    return;
  }
  assert(cisml->num_patterns <= cisml->num_allocated_patterns);
  if (cisml->num_patterns == cisml->num_allocated_patterns) {
    cisml->num_allocated_patterns += PATTERN_INCREMENT;
    cisml->patterns = mm_realloc(
      cisml->patterns,
      cisml->num_allocated_patterns * sizeof(PATTERN_T *)
    );
  }
  cisml->patterns[cisml->num_patterns] = pattern;
  cisml->num_patterns++;
}

/**********************************************************************
  allocate_multi_pattern

  Constructor for the cisml multi_pattern data structure.
**********************************************************************/
MULTI_PATTERN_T *allocate_multi_pattern() {

  // Allocate memory and initialze fields
  MULTI_PATTERN_T *multi_pattern = mm_malloc(sizeof(MULTI_PATTERN_T));
  multi_pattern->score = NULL;
  multi_pattern->pvalue = NULL;
  multi_pattern->num_patterns = 0;
  multi_pattern->num_allocated_patterns = 0;
  multi_pattern->patterns = 0;

  return multi_pattern;
}

/**********************************************************************
  free_multi_pattern

  Destructor for the cisml multi_pattern data structure.
**********************************************************************/
void free_multi_pattern(MULTI_PATTERN_T *multi_pattern) {

  assert(multi_pattern != NULL);

  while(multi_pattern->num_patterns > 0) {
    free_pattern(multi_pattern->patterns[--multi_pattern->num_patterns]);
  }
  myfree(multi_pattern->patterns);

  myfree(multi_pattern->score);
  myfree(multi_pattern->pvalue);

  myfree(multi_pattern);

}

/**********************************************************************
  set_multi_pattern_score

  Sets the score member in a cisml pattern object.
**********************************************************************/
void set_multi_pattern_score(MULTI_PATTERN_T *multi_pattern, double score) {
  assert(multi_pattern != NULL);
  if (multi_pattern->score == NULL) {
    multi_pattern->score = mm_malloc(sizeof(double));
  }
  *(multi_pattern->score) = score;
}

/**********************************************************************
  has_multi_pattern_score

  Does a multi_pattern object have a score?
**********************************************************************/
BOOLEAN_T has_multi_pattern_score(MULTI_PATTERN_T *multi_pattern) {
  assert(multi_pattern != NULL);
  return multi_pattern->score != NULL ? TRUE : FALSE;
}

/**********************************************************************
  clear_multi_pattern_score

  Sets the score member in a cisml multi_pattern object to null.
**********************************************************************/
void clear_multi_pattern_score(MULTI_PATTERN_T *multi_pattern) {
  assert(multi_pattern != NULL);
  if (multi_pattern->score != NULL) {
    myfree(multi_pattern->score);
  }
  multi_pattern->score = NULL;
}

/**********************************************************************
  get_multi_pattern_score

  Gets the score member from a cisml multi_pattern object.
**********************************************************************/
double get_multi_pattern_score(MULTI_PATTERN_T* multi_pattern) {
  assert(multi_pattern != NULL);
  return *(multi_pattern->score);
}

/**********************************************************************
  set_multi_pattern_pvalue

  Sets the pvalue member in a cisml pattern object.
**********************************************************************/
void set_multi_pattern_pvalue(MULTI_PATTERN_T *multi_pattern, double pvalue) {
  assert(multi_pattern != NULL);
  if (multi_pattern->pvalue == NULL) {
    multi_pattern->pvalue = mm_malloc(sizeof(double));
  }
  *(multi_pattern->pvalue) = pvalue;
}

/**********************************************************************
  has_multi_pattern_pvalue

  Does a multi_pattern object have a pvalue?
**********************************************************************/
BOOLEAN_T has_multi_pattern_pvalue(MULTI_PATTERN_T *multi_pattern) {
  assert(multi_pattern != NULL);
  return multi_pattern->pvalue != NULL ? TRUE : FALSE;
}

/**********************************************************************
  clear_multi_pattern_pvalue

  Sets the pvalue member in a cisml multi_pattern object to null.
**********************************************************************/
void clear_multi_pattern_pvalue(MULTI_PATTERN_T *multi_pattern) {
  assert(multi_pattern != NULL);
  if (multi_pattern->pvalue != NULL) {
    myfree(multi_pattern->pvalue);
  }
  multi_pattern->pvalue = NULL;
}

/**********************************************************************
  get_multi_pattern_pvalue

  Gets the pvalue member from a cisml multi_pattern object.
**********************************************************************/
double get_multi_pattern_pvalue(MULTI_PATTERN_T* multi_pattern) {
  assert(multi_pattern != NULL);
  return *(multi_pattern->pvalue);
}

/**********************************************************************
  get_multi_pattern_num_patterns

  Gets the number of patterns from a cisml multi_pattern object.
**********************************************************************/
int get_multi_pattern_num_patterns(MULTI_PATTERN_T *multi_pattern) {
  assert(multi_pattern != NULL);
  return multi_pattern->num_patterns;
}

/**********************************************************************
  get_multi_pattern_patterns

  Gets the array of pointers to patterns from a cisml object.
**********************************************************************/
PATTERN_T **get_multi_pattern_patterns(MULTI_PATTERN_T *multi_pattern) {
  assert(multi_pattern != NULL);
  return multi_pattern->patterns;
}

/**********************************************************************
  add_multi_pattern_pattern

  Adds a pattern to the array of pointers to patterns in a multi_pattern
  object.
**********************************************************************/
void add_multi_pattern_pattern(
  MULTI_PATTERN_T *multi_pattern,
  PATTERN_T* pattern
) {
  assert(multi_pattern != NULL);
  if (pattern == NULL) {
    return;
  }
  assert(multi_pattern->num_patterns <= multi_pattern->num_allocated_patterns);
  if (multi_pattern->num_patterns == multi_pattern->num_allocated_patterns) {
    multi_pattern->num_allocated_patterns += PATTERN_INCREMENT;
    multi_pattern->patterns = mm_realloc(
      multi_pattern->patterns,
      multi_pattern->num_allocated_patterns * sizeof(PATTERN_T *)
    );
  }
  multi_pattern->patterns[multi_pattern->num_patterns] = pattern;
  multi_pattern->num_patterns++;
}

/**********************************************************************
  allocate_pattern

  Constructor for the cisml pattern data structure.
  Sets required fields to point to copies of the provided arguments.
  Other fields set to NULL.
**********************************************************************/
PATTERN_T *allocate_pattern(char *accession, char *name) {

  assert(accession != NULL);
  assert(name != NULL);

  // Allocate memory and initialze fields
  PATTERN_T *pattern = mm_malloc(sizeof(PATTERN_T));
  pattern->accession = NULL;
  pattern->name = NULL;
  pattern->pvalue = NULL;
  pattern->score = NULL;
  pattern->db = NULL;
  pattern->lsid = NULL;
  pattern->sequences = NULL;
  pattern->elements = NULL;

  pattern->num_allocated_sequences = 0;
  pattern->num_allocated_elements = 0;

  pattern->num_sequences = 0;
  pattern->num_scanned_positions = 0L;
  pattern->max_stored_matches = 100000;
  pattern->num_stored_matches = 0;
  pattern->max_pvalue_retained = 1.0;
  pattern->qvalues_computed = FALSE;
  pattern->has_all_pvalues = TRUE;
  pattern->is_complete = FALSE;

  // Set required fields
  int length = strlen(accession) + 1;
  pattern->accession = mm_malloc(length * sizeof(char));
  strncpy(pattern->accession, accession, length);
  length = strlen(name) + 1;
  pattern->name = mm_malloc(length * sizeof(char));
  strncpy(pattern->name, name, length);
  pattern->element_heap = create_heap(
    pattern->max_stored_matches,
    compare_matched_elements,
    copy_matched_element,
    destroy_matched_element,
    NULL, // Key function
    NULL  // Print function
  );
  pattern->elements = NULL;

  return pattern;
}

/**********************************************************************
  free_pattern

  Destructor for the cisml pattern data structure.
**********************************************************************/
void free_pattern(PATTERN_T *pattern) {

  assert(pattern != NULL);

  myfree(pattern->lsid);
  myfree(pattern->db);
  myfree(pattern->score);
  myfree(pattern->pvalue);
  myfree(pattern->name);
  myfree(pattern->accession);

  while (pattern->num_sequences > 0) {
    free_scanned_sequence(pattern->sequences[--pattern->num_sequences]);
  }
  pattern->num_allocated_sequences = 0;
  pattern->num_sequences = 0;
  myfree(pattern->sequences)

  pattern->num_scanned_positions = 0L;
  destroy_heap(pattern->element_heap);
  // Free any elements stored in the elements array
  if (pattern->num_stored_matches > 0) {
    int i = 0;
    for (i = 0; i < pattern->num_stored_matches; i++) {
      free_matched_element(pattern->elements[i]);
    }
  }
  myfree(pattern->elements)
  myfree(pattern);

}

/**********************************************************************
  get_pattern_has_all_pvalues

  Returns the flag indicating whether or not all matched elements have
  been retained, or only those with the smallest p-values.
**********************************************************************/
BOOLEAN_T get_pattern_has_all_pvalues(PATTERN_T *pattern) {

  assert(pattern != NULL);

  return pattern->has_all_pvalues;

}

/**********************************************************************
  set_pattern_has_all_pvalues

  Sets the flag indicating whether or not all matched elements have
  been retained, or only those with the smallest p-values.
**********************************************************************/
void set_pattern_has_all_pvalues(PATTERN_T *pattern, BOOLEAN_T has_all_pvalues) {
  assert(pattern != NULL);
  pattern->has_all_pvalues = has_all_pvalues;
}

/**********************************************************************
  get_pattern_is_complete

  Returns a flag indicating whether or not all matched elements have
  been added to the pattern. If flag is true, the element heap is no
  longer available and all matched elements are stored in an array
  of matched element pointers, sorted by p-value.
**********************************************************************/
BOOLEAN_T get_pattern_is_complete(PATTERN_T *pattern) {
  assert(pattern != NULL);
  return pattern->is_complete;
}

/**********************************************************************
  set_pattern_is_complete

  Sets the flag indicating that all matched elements have
  been added to the pattern to true.  Moves  all matched elments out of
  the element heap into an array of matched element pointers, sorted by
  p-value. No further elements can be added to the pattern once
  this function has been called.
**********************************************************************/
void set_pattern_is_complete(PATTERN_T *pattern) {

  assert(pattern != NULL);
  assert(pattern->is_complete == FALSE);

  pattern->is_complete = TRUE;

  // Now that the pattern is complete, move elements from heap to array
  // sorted by p-value
  int num_elements = pattern->num_stored_matches;
  int i_element = 0;
  MATCHED_ELEMENT_T *element = NULL;
  pattern->elements = mm_malloc(sizeof(MATCHED_ELEMENT_T *) * num_elements);
  // Elements come off the heap in descending p-value order
  // Need to have array in ascending p-value order.
  for (i_element = num_elements - 1; i_element >= 0; --i_element) {
    element = (MATCHED_ELEMENT_T *) pop_heap_root(pattern->element_heap);
    pattern->elements[i_element] = element;
  }


  // Now that the pattern is complete, update the scanned sequeences
  // with the matched elements
  add_pattern_elements_to_scanned_seq(pattern);

}

/**********************************************************************
  get_pattern_max_pvalue_retained

  Returns the maximum p-value of the matched elements retained by
  the pattern.
**********************************************************************/
double get_pattern_max_pvalue_retained(PATTERN_T *pattern) {
  assert(pattern != NULL);
  return pattern->max_pvalue_retained;
}

/**********************************************************************
  set_pattern_max_pvalue_retained

  Sets the maximum p-value of the matched elements retained by
  the pattern.
**********************************************************************/
void set_pattern_max_pvalue_retained(PATTERN_T *pattern, double max_pvalue) {
  assert(pattern != NULL);
  pattern->max_pvalue_retained = max_pvalue;
}

/**********************************************************************
  get_pattern_name

  Gets the program_name member from a cisml pattern object.
**********************************************************************/
char *get_pattern_name(PATTERN_T *pattern) {
  assert(pattern != NULL);
  return pattern->name;
}

/**********************************************************************
  get_pattern_accession

  Gets the accession member from a cisml pattern object.
**********************************************************************/
char *get_pattern_accession(PATTERN_T* pattern) {
  assert(pattern != NULL);
  return pattern->accession;
}

/**********************************************************************
  set_pattern_pvalue

  Sets the pvalue member in a cisml pattern object.
**********************************************************************/
void set_pattern_pvalue(PATTERN_T *pattern, double pvalue) {
  assert(pattern != NULL);
  if (pattern->pvalue == NULL) {
    pattern->pvalue = mm_malloc(sizeof(double));
  }
  *(pattern->pvalue) = pvalue;
}

/**********************************************************************
  clear_pattern_pvalue

  Sets the pvalue member in a cisml pattern object to null.
**********************************************************************/
void clear_pattern_pvalue(PATTERN_T *pattern) {
  assert(pattern != NULL);
  if (pattern->pvalue != NULL) {
    myfree(pattern->pvalue);
  }
  pattern->pvalue = NULL;
}

/**********************************************************************
  has_pattern_pvalue

  Does a pattern object have a pvalue?
**********************************************************************/
BOOLEAN_T has_pattern_pvalue(PATTERN_T *pattern) {
  return pattern->pvalue != NULL ? TRUE : FALSE;
}

/**********************************************************************
  get_pattern_pvalue

  Gets the pvalue member from a cisml pattern object.
**********************************************************************/
double get_pattern_pvalue(PATTERN_T* pattern) {
  assert(pattern != NULL);
  return *(pattern->pvalue);
}

/**********************************************************************
  has_pattern_qvalues

  Does the matched-elements for the pattern have qvalues available?
**********************************************************************/
BOOLEAN_T has_pattern_qvalues(PATTERN_T *pattern) {
  return pattern->qvalues_computed;
}

/**********************************************************************
  set_pattern_score

  Sets the score member in a cisml pattern object.
**********************************************************************/
void set_pattern_score(PATTERN_T *pattern, double score) {
  assert(pattern != NULL);
  if (pattern->score == NULL) {
    pattern->score = mm_malloc(sizeof(double));
  }
  *(pattern->score) = score;
}

/**********************************************************************
  clear_pattern_score

  Sets the score member in a cisml pattern object to null.
**********************************************************************/
void clear_pattern_score(PATTERN_T *pattern) {
  assert(pattern != NULL);
  if (pattern->score != NULL) {
    myfree(pattern->score);
  }
  pattern->score = NULL;
}

/**********************************************************************
  has_pattern_score

  Does a pattern object have a score?
**********************************************************************/
BOOLEAN_T has_pattern_score(PATTERN_T *pattern) {
  return pattern->score != NULL ? TRUE : FALSE;
}

/**********************************************************************
  get_pattern_score

  Gets the score member from a cisml pattern object.
**********************************************************************/
double get_pattern_score(PATTERN_T* pattern) {
  assert(pattern != NULL);
  return *(pattern->score);
}

/**********************************************************************
  set_pattern_db

  Sets the db member in a cisml pattern object.
**********************************************************************/
void set_pattern_db(PATTERN_T *pattern, char *db) {
  assert(pattern != NULL);
  if (db == NULL) {
    if (pattern->db != NULL) {
      myfree(pattern->db);
    }
    pattern->db = NULL;
  }
  else {
    int old_length = 0;
    int new_length = strlen(db) + 1;
    if (pattern->db != NULL) {
      old_length = strlen(pattern->db) + 1;
    }
    if (old_length < new_length) {
      pattern->db = mm_realloc(pattern->db, new_length);
    }
    strncpy(pattern->db, db, new_length);
  }
}

/**********************************************************************
  get_pattern_db

  Gets the db member from a cisml pattern object.
  May return NULL.
**********************************************************************/
char *get_pattern_db(PATTERN_T* pattern) {
  assert(pattern != NULL);
  return pattern->db;
}

/**********************************************************************
  set_pattern_lsid

  Sets the lsid member in a cisml pattern object.
**********************************************************************/
void set_pattern_lsid(PATTERN_T *pattern, char *lsid) {
  assert(pattern != NULL);
  if (lsid == NULL) {
    if (pattern->lsid != NULL) {
      myfree(pattern->lsid);
    }
    pattern->lsid = NULL;
  }
  else {
    int old_length = 0;
    int new_length = strlen(lsid) + 1;
    if (pattern->lsid != NULL) {
      old_length = strlen(pattern->lsid) + 1;
    }
    if (old_length < new_length) {
      pattern->lsid = mm_realloc(pattern->lsid, new_length);
    }
    strncpy(pattern->lsid, lsid, new_length);
  }
}

/**********************************************************************
  get_pattern_lsid

  Gets the lsid member from a cisml pattern object.
  May return NULL.
**********************************************************************/
char *get_pattern_lsid(PATTERN_T* pattern) {
  assert(pattern != NULL);
  return pattern->lsid;
}

/**********************************************************************
  get_pattern_num_scanned_sequences

  Gets the number of scanned_sequence objects in a cisml pattern object.
**********************************************************************/
int get_pattern_num_scanned_sequences(PATTERN_T *pattern) {
  assert(pattern != NULL);
  return pattern->num_sequences;
}

/**********************************************************************
  get_pattern_max_stored_matches

  Gets the maximum number of matched elements that will be stored in a cisml
  pattern object.
**********************************************************************/
int get_pattern_max_stored_matches(PATTERN_T *pattern) {
  assert(pattern != NULL);
  return pattern->max_stored_matches;
}

/**********************************************************************
  set_pattern_max_stored_matches

  Sets the maximum number of matched elements that will be stored in a cisml
  pattern object. This requires creating a new element heap using the new
  maximum number of elements, copying the elements from the existing heap,
  freeing the existing heap, and pointing the pattern to the new heap.

  It fails if the new max is smaller than the current number of nodes.

  Returns TRUE if successful, FALSE otherwise.
**********************************************************************/
BOOLEAN_T set_pattern_max_stored_matches(PATTERN_T *pattern, int max) {

  assert(pattern != NULL);

  HEAP *current_heap = pattern->element_heap;
  if (pattern->max_stored_matches == max) {
    // No point in doing anything, the max isn't changing.
    return TRUE;
  }
  else if (max > get_num_nodes(current_heap)) {
    HEAP *new_heap = create_heap(
      max,
      compare_matched_elements,
      copy_matched_element,
      destroy_matched_element,
      NULL, // Key function
      NULL  // Print function
    );
    MATCHED_ELEMENT_T *element = NULL;
    while ((element = pop_heap_root(current_heap))) {
      add_node_heap(new_heap, element);
    }
    destroy_heap(current_heap);
    pattern->element_heap = new_heap;
    pattern->max_stored_matches = max;
    return TRUE;
  }
  else {
    // Can't make max size of heap smaller than the current number of nodes.
    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(
        stderr,
        "Warning: The maximum size of the heap cannot be decreased.\n"
      );
    }
    return FALSE;
  }
}

/**********************************************************************
  get_pattern_num_scanned_positions

  Gets the number of sites scanned with a cisml pattern object.
**********************************************************************/
long get_pattern_num_scanned_positions(PATTERN_T *pattern) {
  assert(pattern != NULL);
  assert(pattern->num_scanned_positions >= 0L);
  return pattern->num_scanned_positions;
}

/**********************************************************************
  get_pattern_num_stored_matches

  Gets the total number of matched element objects stored in a cisml
  pattern object.
**********************************************************************/
int get_pattern_num_stored_matches(PATTERN_T *pattern) {
  assert(pattern != NULL);
  return pattern->num_stored_matches;
}

/**********************************************************************
  add_pattern_matched_element

  Adds a pointer to a matched element to the heap of pointers to
  matched_elements in a cisml pattern object.

  Before adding to the heap, check the pvalue for the element.
  if it's greater the pvalues we've already thrown away, we refuse to
  add it. Also, check max_stored_matches. If we've
  hit the limit, purge the heap of the least significian matched elements.
  Set to FALSE the Boolean that indicates that set of matched elements
  is complete.

  Returns TRUE if the element was added, FALSE otherwise
**********************************************************************/
BOOLEAN_T add_pattern_matched_element(
  PATTERN_T *pattern,
  MATCHED_ELEMENT_T *element
) {

  assert(element != NULL);
  assert(pattern != NULL);

  if (pattern->is_complete == TRUE) {
    // Don't add matched elements if pattern is marked as complete.
    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(
        stderr,
        "Warning: trying to add matched elements to pattern marked as complete.\n"
      );
    }
    return FALSE;
  }

  if (element->pvalue > pattern->max_pvalue_retained) {
    // Don't add the matched element if its pvalue is less
    // significant then pvalues we've already thrown away.
    return FALSE;
  }

  if (pattern->max_stored_matches > 0
      && (long) pattern->max_stored_matches <= pattern->num_stored_matches) {
    // Max element storage has been reached. We have to drop
    // the least significant elements.
    reduce_pattern_matched_elements(pattern);
  }

  // The least significant pvalue may have changed
  // during the reduction.
  if (element->pvalue <= pattern->max_pvalue_retained) {
    // Add the element to the pattern.
    add_node_heap(pattern->element_heap, (void *) element);
    ++pattern->num_stored_matches;
    return TRUE;
  }
  else {
    // Don't bother adding elements who's pvalue execeeds the
    // least sig. matched element retained.
    return FALSE;
  }
}

/**********************************************************************
  add_pattern_elements_to_scanned_seq

  Updates each scanned sequence belonging to pattern with the
  matched_elements assocaited with that sequence.
  Should  not be called until pattern is complete.

**********************************************************************/
void add_pattern_elements_to_scanned_seq(PATTERN_T *pattern) {

  assert(pattern != NULL);
  assert(pattern->is_complete == TRUE);
  assert(pattern->elements != NULL);

  int element_index = 0;
  for (element_index = 0; element_index < pattern->num_stored_matches; ++element_index) {
    MATCHED_ELEMENT_T *element = (pattern->elements)[element_index];
    SCANNED_SEQUENCE_T *seq = element->parent_sequence;
    add_scanned_sequence_matched_element(seq, element);
  }

}

/**********************************************************************
  reduce_pattern_matched_elements

  To conserve memory, remove the least significant matched_elements.
  We want to remove at least PERCENT_ELEMENT_DISCARD * num_matched_elements,
  but we may remove more if there remain matched elments with a
  p-value equal to the delected elements.

**********************************************************************/
static void reduce_pattern_matched_elements(PATTERN_T *pattern) {

  assert(pattern != NULL);
  assert(pattern->is_complete == FALSE);

  const float PERCENT_ELEMENT_DISCARD = 0.5;
  static BOOLEAN_T have_moved_sequences = FALSE;

  HEAP *heap = pattern->element_heap;

  // Delete PERCENT_ELEMENT_DISCARD with highest pvalue
  int num_elements_to_delete = PERCENT_ELEMENT_DISCARD * get_num_nodes(heap);
  if (verbosity > NORMAL_VERBOSE) {
	  fprintf(
      stderr,
      "Deleting at least %d matched elements from pattern %s.\n",
      num_elements_to_delete,
      pattern->name
    );
  }
  // Delete least significant matched elements.
  double min_pvalue_discarded = 1.0;
  MATCHED_ELEMENT_T *victim = NULL;
  int deletion_count = 0;
  for (deletion_count = 0; deletion_count  < num_elements_to_delete; ++deletion_count) {
    victim = (MATCHED_ELEMENT_T *) pop_heap_root(heap);
    min_pvalue_discarded = victim->pvalue;
    --pattern->num_stored_matches;
    free_matched_element(victim);
  }

  // Keep deleting matched elements until we find an element more 
  // significant then the elements we've already deleted.
  while (((MATCHED_ELEMENT_T *) get_node(heap, 1))->pvalue 
         >= min_pvalue_discarded) {
    victim = (MATCHED_ELEMENT_T *) pop_heap_root(heap);
    assert(victim != NULL);
    --pattern->num_stored_matches;
    free_matched_element(victim);
    if (get_num_nodes(heap) == 0) {
      // All the matched elements have been deleted!
      break;
    }
  }

  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(
      stderr, 
      "Warning: Reached max stored scores (%d).\n" 
      "Motif matches with p-value >= %3.2g have been "
      "deleted to reclaim memory.\n", 
      pattern->max_stored_matches,
      min_pvalue_discarded
    );
  }

  if (get_num_nodes(heap) > 0) {
    // Get the largest p-value retained from the top element of the heap.
    pattern->max_pvalue_retained = ((MATCHED_ELEMENT_T *) get_node(heap, 1))->pvalue;
  }
  else {
    // All items have been deleted!
    fprintf(
      stderr, 
      "Warning: there are no motif matches with p-value < %3.2g.\n"
      "Use --max-stored-scores to allocate more space for "
      "storing motif matches.\n", 
     min_pvalue_discarded 
    );
    // Set the largest p-value retained to something
    // slightly less the smallest p-value discarded.
    pattern->max_pvalue_retained 
      = get_next_smaller_double(min_pvalue_discarded);
  }

  set_pattern_has_all_pvalues(pattern, FALSE);
}

/**********************************************************************
  add_pattern_scanned_sequence

  Adds a pointer to a scanned_sequence to the array of pointers to
  scanned_sequences in a cisml pattern object.
**********************************************************************/
static void add_pattern_scanned_sequence(
  PATTERN_T *pattern,
  SCANNED_SEQUENCE_T *sequence
) {

  assert(pattern != NULL);
  assert(sequence != NULL);
  assert(pattern->num_sequences <= pattern->num_allocated_sequences);

  sequence->parent_pattern = pattern;

  if (pattern->num_sequences == pattern->num_allocated_sequences) {
    pattern->num_allocated_sequences += SEQUENCE_INCREMENT;
    pattern->sequences = mm_realloc(
      pattern->sequences,
      pattern->num_allocated_sequences * sizeof(SCANNED_SEQUENCE_T *)
    );
  }
  pattern->sequences[pattern->num_sequences] = sequence;
  pattern->num_sequences++;

}

/**********************************************************************
  get_pattern_scanned_sequences

  Gets the array of pointers to scanned_sequence objects in a cisml
  pattern object.  May return NULL.
**********************************************************************/
SCANNED_SEQUENCE_T **get_pattern_scanned_sequences(PATTERN_T *pattern) {
  assert(pattern != NULL);
  return pattern->sequences;
}

/**********************************************************************
  allocate_scanned_sequence

  Constructor for the cisml scanned_sequence data structure.
  Sets required fields to point to copies of the provided arguments.
  Other fields set to NULL.
**********************************************************************/
SCANNED_SEQUENCE_T *allocate_scanned_sequence(
  char *accession,
  char *name,
  PATTERN_T *parent_pattern
) {

  assert(accession != NULL);
  assert(name != NULL);

  // Allocate memory and initialze fields
  SCANNED_SEQUENCE_T *scanned_sequence = mm_malloc(sizeof(SCANNED_SEQUENCE_T));
  scanned_sequence->accession = NULL;
  scanned_sequence->name = NULL;
  scanned_sequence->pvalue = NULL;
  scanned_sequence->score = NULL;
  scanned_sequence->length = NULL;
  scanned_sequence->db = NULL;
  scanned_sequence->lsid = NULL;
  scanned_sequence->num_scanned_positions = 0L;
  scanned_sequence->num_matched_elements = 0;
  scanned_sequence->num_allocated_elements = 0;

  // Set required fields
  int length = strlen(accession) + 1;
  scanned_sequence->accession = mm_malloc(length * sizeof(char));
  strncpy(scanned_sequence->accession, accession, length);
  length = strlen(name) + 1;
  scanned_sequence->name = mm_malloc(length * sizeof(char));
  strncpy(scanned_sequence->name, name, length);
  add_pattern_scanned_sequence(parent_pattern, scanned_sequence);
  scanned_sequence->elements = NULL;

  return scanned_sequence;
}

/**********************************************************************
  free_scanned_sequence

  Destructor for the cisml scanned_sequence data structure.
**********************************************************************/
void free_scanned_sequence(SCANNED_SEQUENCE_T *scanned_sequence) {

  assert(scanned_sequence != NULL);

  myfree(scanned_sequence->lsid);
  myfree(scanned_sequence->db);
  myfree(scanned_sequence->length);
  myfree(scanned_sequence->score);
  myfree(scanned_sequence->pvalue);
  myfree(scanned_sequence->name);
  myfree(scanned_sequence->accession);
  myfree(scanned_sequence->elements); // Don't free individual elements
                                      // they belong to pattern.

  myfree(scanned_sequence);

}

/**********************************************************************
  get_scanned_sequence_accession

  Gets the accession member from a cisml scanned_sequence object.
**********************************************************************/
char *get_scanned_sequence_accession(SCANNED_SEQUENCE_T* scanned_sequence) {

  assert(scanned_sequence != NULL);

  return scanned_sequence->accession;

}

/**********************************************************************
  get_scanned_sequence_name

  Gets the program_name member from a cisml scanned_sequence object.
**********************************************************************/
char *get_scanned_sequence_name(SCANNED_SEQUENCE_T *scanned_sequence) {

  assert(scanned_sequence != NULL);

  return scanned_sequence->name;

}

/**********************************************************************
  get_scanned_sequence_parent

  Get a pointer to the pattern pointing to this scanned_sequence object.
**********************************************************************/
PATTERN_T *get_scanned_sequence_parent(
  SCANNED_SEQUENCE_T *scanned_sequence
) {

  assert(scanned_sequence != NULL);

  return scanned_sequence->parent_pattern;

}

/**********************************************************************
  set_scanned_sequence_pvalue

  Sets the pvalue member in a cisml scanned_sequence object.
**********************************************************************/
void set_scanned_sequence_pvalue(
  SCANNED_SEQUENCE_T *scanned_sequence,
  double pvalue
) {

  assert(scanned_sequence != NULL);

  if (scanned_sequence->pvalue == NULL) {
    scanned_sequence->pvalue = mm_malloc(sizeof(double));
  }
  *(scanned_sequence->pvalue) = pvalue;

}

/**********************************************************************
  clear_scanned_sequence_pvalue

  Sets the pvalue member in a cisml scanned_sequence object to null.
**********************************************************************/
void clear_scanned_sequence_pvalue(SCANNED_SEQUENCE_T *scanned_sequence) {

  assert(scanned_sequence != NULL);

  if (scanned_sequence->pvalue != NULL) {
    myfree(scanned_sequence->pvalue);
  }

  scanned_sequence->pvalue = NULL;

}

/**********************************************************************
  has_scanned_sequence_pvalue

  Does a scanned_sequence object have a pvalue?
**********************************************************************/
BOOLEAN_T has_scanned_sequence_pvalue(SCANNED_SEQUENCE_T *scanned_sequence) {
  assert(scanned_sequence != NULL);
  return scanned_sequence->pvalue != NULL ? TRUE : FALSE;

}

/**********************************************************************
  get_scanned_sequence_pvalue

  Gets the pvalue member from a cisml scanned_sequence object.
**********************************************************************/
double get_scanned_sequence_pvalue(SCANNED_SEQUENCE_T* scanned_sequence) {
  assert(scanned_sequence != NULL);
  return *(scanned_sequence->pvalue);

}

/**********************************************************************
  set_scanned_sequence_score

  Sets the score member in a cisml scanned_sequence object.
**********************************************************************/
void set_scanned_sequence_score(
  SCANNED_SEQUENCE_T *scanned_sequence,
  double score
) {
  assert(scanned_sequence != NULL);
  if (scanned_sequence->score == NULL) {
    scanned_sequence->score = mm_malloc(sizeof(double));
  }
  *(scanned_sequence->score) = score;
}

/**********************************************************************
  clear_scanned_sequence_score

  Sets the score member in a cisml scanned_sequence object to null.
**********************************************************************/
void clear_scanned_sequence_score(SCANNED_SEQUENCE_T *scanned_sequence) {
  assert(scanned_sequence != NULL);
  if (scanned_sequence->score != NULL) {
    myfree(scanned_sequence->score);
  }
  scanned_sequence->score = NULL;
}

/**********************************************************************
  has_scanned_sequence_score

  Does a scanned_sequence object have a score?
**********************************************************************/
BOOLEAN_T has_scanned_sequence_score(SCANNED_SEQUENCE_T *scanned_sequence) {
  assert(scanned_sequence != NULL);
  return scanned_sequence->score != NULL ? TRUE : FALSE;
}

/**********************************************************************
  get_scanned_sequence_score

  Gets the score member from a cisml scanned_sequence object.
**********************************************************************/
double get_scanned_sequence_score(SCANNED_SEQUENCE_T* scanned_sequence) {
  assert(scanned_sequence != NULL);
  return *(scanned_sequence->score);
}

/**********************************************************************
  set_scanned_sequence_length

  Sets the length member in a cisml scanned_sequence object.
**********************************************************************/
void set_scanned_sequence_length(
  SCANNED_SEQUENCE_T *scanned_sequence,
  int length
) {
  assert(scanned_sequence != NULL);
  if (scanned_sequence->length == NULL) {
    scanned_sequence->length = mm_malloc(sizeof(int));
  }
  *(scanned_sequence->length) = length;
}

/**********************************************************************
  clear_scanned_sequence_length

  Sets the length member in a cisml scanned_sequence object to null.
**********************************************************************/
void clear_scanned_sequence_length(SCANNED_SEQUENCE_T *scanned_sequence) {
  assert(scanned_sequence != NULL);
  if (scanned_sequence->length != NULL) {
    myfree(scanned_sequence->length);
  }
  scanned_sequence->length = NULL;
}

/**********************************************************************
  has_scanned_sequence_length

  Does a scanned_sequence object have a length?
**********************************************************************/
BOOLEAN_T has_scanned_sequence_length(SCANNED_SEQUENCE_T *scanned_sequence) {
  return scanned_sequence->length != NULL ? TRUE : FALSE;
}

/**********************************************************************
  get_scanned_sequence_length

  Gets the length member from a cisml scanned_sequence object.
**********************************************************************/
int get_scanned_sequence_length(SCANNED_SEQUENCE_T* scanned_sequence) {
  assert(scanned_sequence != NULL);
  return *(scanned_sequence->length);
}

/**********************************************************************
  set_scanned_sequence_db

  Sets the db member in a cisml scanned_sequence object.
**********************************************************************/
void set_scanned_sequence_db(SCANNED_SEQUENCE_T *scanned_sequence, char *db) {
  assert(scanned_sequence != NULL);
  if (db == NULL) {
    if (scanned_sequence->db != NULL) {
      myfree(scanned_sequence->db);
    }
    scanned_sequence->db = NULL;
  }
  else {
    int old_length = 0;
    int new_length = strlen(db) + 1;
    if (scanned_sequence->db != NULL) {
      old_length = strlen(scanned_sequence->db) + 1;
    }
    if (old_length < new_length) {
      scanned_sequence->db = mm_realloc(scanned_sequence->db, new_length);
    }
    strncpy(scanned_sequence->db, db, new_length);
  }
}

/**********************************************************************
  get_scanned_sequence_db

  Gets the db member from a cisml scanned_sequence object.
  May return NULL.
**********************************************************************/
char *get_scanned_sequence_db(SCANNED_SEQUENCE_T* scanned_sequence) {
  assert(scanned_sequence != NULL);
  return scanned_sequence->db;
}

/**********************************************************************
  set_scanned_sequence_lsid

  Sets the lsid member in a cisml scanned_sequence object.
**********************************************************************/
void set_scanned_sequence_lsid(
  SCANNED_SEQUENCE_T *scanned_sequence,
  char *lsid
) {
  assert(scanned_sequence != NULL);
  if (lsid == NULL) {
    if (scanned_sequence->lsid != NULL) {
      myfree(scanned_sequence->lsid);
    }
    scanned_sequence->lsid = NULL;
  }
  else {
    int old_length = 0;
    int new_length = strlen(lsid) + 1;
    if (scanned_sequence->lsid != NULL) {
      old_length = strlen(scanned_sequence->lsid) + 1;
    }
    if (old_length < new_length) {
      scanned_sequence->lsid = mm_realloc(scanned_sequence->lsid, new_length);
    }
    strncpy(scanned_sequence->lsid, lsid, new_length);
  }
}

/**********************************************************************
  get_scanned_sequence_lsid

  Gets the lsid member from a cisml scanned_sequence object.
  May return NULL.
**********************************************************************/
char *get_scanned_sequence_lsid(SCANNED_SEQUENCE_T* scanned_sequence) {
  assert(scanned_sequence != NULL);
  return scanned_sequence->lsid;
}

/**********************************************************************
  add_scanned_sequence_matched_element

  Adds a matched element to the array of matched elements in
  the scanned sequence.
**********************************************************************/
void add_scanned_sequence_matched_element(
  SCANNED_SEQUENCE_T *sequence,
  MATCHED_ELEMENT_T *element
) {

  assert(sequence != NULL);
  assert(element != NULL);
  assert(sequence->num_matched_elements <= sequence->num_allocated_elements);

  if (sequence->num_matched_elements == sequence->num_allocated_elements) {
    sequence->num_allocated_elements += ELEMENT_INCREMENT;
    sequence->elements = mm_realloc(
      sequence->elements,
      sequence->num_allocated_elements * sizeof(MATCHED_ELEMENT_T *)
    );
  }
  sequence->elements[sequence->num_matched_elements] = element;
  sequence->num_matched_elements++;
}

/**********************************************************************
  add_scanned_sequence_scanned_element

  Increments the count of scanned elements in a scanned_sequence.
**********************************************************************/
void add_scanned_sequence_scanned_element( SCANNED_SEQUENCE_T *sequence) {
  assert(sequence != NULL);
  sequence->num_scanned_positions++;
  sequence->parent_pattern->num_scanned_positions++;
}

/**********************************************************************
  get_scanned_sequence_num_matched_elements

  Gets the number of matched_element objects in a cisml
  scanned_sequence object.
**********************************************************************/
int get_scanned_sequence_num_matched_elements(SCANNED_SEQUENCE_T *sequence) {
  assert(sequence != NULL);
  return sequence->num_matched_elements;
}

/**********************************************************************
  get_scanned_sequence_num_scanned_positions

  Gets the number of positions in the scanned_sequence where we
  have scanned for a matched_element.
**********************************************************************/
long get_scanned_sequence_num_scanned_positions(SCANNED_SEQUENCE_T *sequence) {
  assert(sequence != NULL);
  return sequence->num_scanned_positions;
}

/**********************************************************************
  get_scanned_sequences_matched_elements

  Gets the array of pointers to matched_element objects in a cisml
  scanned_sequence object.  May return NULL.
**********************************************************************/
MATCHED_ELEMENT_T **get_scanned_sequence_matched_elements(
    SCANNED_SEQUENCE_T *sequence
) {
  assert(sequence != NULL);
  return sequence->elements;
}

/**********************************************************************
  allocate_matched_element

  Constructor for the cisml matched_element data structure.
  Sets required fields to the provided arguments.
  Other fields set to NULL.
**********************************************************************/
MATCHED_ELEMENT_T *allocate_matched_element(
  int start,
  int stop,
  SCANNED_SEQUENCE_T *parent_sequence
) {

  // Allocate memory and set required fields
  MATCHED_ELEMENT_T *element = mm_malloc(sizeof(MATCHED_ELEMENT_T));
  element->start = start;
  element->stop = stop;
  element->parent_sequence = parent_sequence;

  // Initialze optional fields
  element->score = 0.0;
  element->has_score = FALSE;
  element->pvalue = 0.0;
  element->has_pvalue = FALSE;
  element->qvalue = 0.0;
  element->has_qvalue = FALSE;
  element->clusterid = NULL;
  element->sequence = NULL;
  element->strand = '\0';

  return element;
}

/**********************************************************************
  allocate_matched_element_without_inversion

  Alternative Constructor for the cisml matched_element data structure.
  Sets required fields to the provided arguments.
  Other fields set to NULL.

  JH: I had to add this because the existing constructor inverts the
      DNA sequence if the start site is greater than the stop site. I
      could not simply copy an existing matched element's content and
      feed it to the constructor without the content being changed.

**********************************************************************/
MATCHED_ELEMENT_T *allocate_matched_element_without_inversion(
  int start,
  int stop,
  const char *seq,
  SCANNED_SEQUENCE_T *parent
) {

  // Allocate memory and set required fields
  MATCHED_ELEMENT_T *element = mm_malloc(sizeof(MATCHED_ELEMENT_T));
  element->start = start;
  element->stop = stop;
  int length = strlen(seq) + 1;
  element->sequence = mm_malloc(length * sizeof(char));
  strncpy(element->sequence, seq, length);

  element->parent_sequence = parent;

  // Initialze optional fields
  element->score = 0.0;
  element->has_score = FALSE;
  element->pvalue = 0.0;
  element->has_pvalue = FALSE;
  element->qvalue = 0.0;
  element->has_qvalue = FALSE;
  element->clusterid = NULL;
  element->strand = '\0';

  return element;
}


/**********************************************************************
  allocate_matched_element_with_score

  Constructor for the cisml matched_element data structure.
  Sets required fields to the provided arguments.
  Sets score and pvalue.
  Other fields set to NULL.

**********************************************************************/
MATCHED_ELEMENT_T *allocate_matched_element_with_score(
  int start,
  int stop,
  double score,
  double pvalue,
  SCANNED_SEQUENCE_T *parent_sequence
) {

  // Allocate memory and set required fields
  MATCHED_ELEMENT_T *element = allocate_matched_element(start, stop, parent_sequence);

  set_matched_element_score(element, score);
  set_matched_element_pvalue(element, pvalue);

  return element;
}

/**********************************************************************
  free_matched_element

  Destructor for the cisml matched_element data structure.
**********************************************************************/
void free_matched_element(MATCHED_ELEMENT_T *element) {

  assert(element != NULL);

  myfree(element->clusterid);
  myfree(element->sequence);

  myfree(element);

}

/**********************************************************************
  get_matched_element_start

  Gets the start member from a cisml matched_element object.
**********************************************************************/
int get_matched_element_start(MATCHED_ELEMENT_T* element) {

  assert(element != NULL);

  return element->start;

}
/**********************************************************************
  set_matched_element_start

  Sets the start member from a cisml matched_element object.
**********************************************************************/
void set_matched_element_start(MATCHED_ELEMENT_T* matched_element, int newstart) {
  assert(matched_element != NULL);
  matched_element->start = newstart;
}

/**********************************************************************
  get_matched_element_stop

  Gets the stop member from a cisml matched_element object.
**********************************************************************/
int get_matched_element_stop(MATCHED_ELEMENT_T* element) {
  assert(element != NULL);
  return element->stop;
}

/**********************************************************************
  set_matched_element_score

  Sets the score member in a cisml matched_element object.
**********************************************************************/
void set_matched_element_score(
  MATCHED_ELEMENT_T *element,
  double score
) {
  assert(element != NULL);
  element->score = score;
  element->has_score = TRUE;
}

/**********************************************************************
  has_matched_element_score

  Does a matched_element object have a score?
**********************************************************************/
BOOLEAN_T has_matched_element_score(MATCHED_ELEMENT_T *element) {
  return element->has_score;
}

/**********************************************************************
  get_matched_element_score

  Gets the score member from a cisml matched_element object.
**********************************************************************/
double get_matched_element_score(MATCHED_ELEMENT_T* element) {
  assert(element != NULL);
  return element->score;
}

/**********************************************************************
  set_matched_element_pvalue

  Sets the pvalue member in a cisml matched_element object.
**********************************************************************/
void set_matched_element_pvalue(
  MATCHED_ELEMENT_T *element,
  double pvalue
) {
  assert(element != NULL);
  element->pvalue = pvalue;
  element->has_pvalue = TRUE;
}

/**********************************************************************
  has_matched_element_pvalue

  Does a matched_element object have a pvalue?
**********************************************************************/
BOOLEAN_T has_matched_element_pvalue(MATCHED_ELEMENT_T *element) {
  return element->has_pvalue;
}

/**********************************************************************
  get_matched_element_pvalue

  Gets the pvalue member from a cisml matched_element object.
**********************************************************************/
double get_matched_element_pvalue(MATCHED_ELEMENT_T* element) {
  assert(element != NULL);
  return element->pvalue;
}

/**********************************************************************
  set_matched_element_qvalue

  Sets the qvalue member in a cisml matched_element object.
**********************************************************************/
void set_matched_element_qvalue(
  MATCHED_ELEMENT_T *element,
  double qvalue
) {
  assert(element != NULL);
  element->qvalue = qvalue;
  element->has_qvalue = TRUE;
}

/**********************************************************************
  has_matched_element_qvalue

  Does a matched_element object have a qvalue?
**********************************************************************/
BOOLEAN_T has_matched_element_qvalue(MATCHED_ELEMENT_T *element) {
  return element->has_qvalue;
}

/**********************************************************************
  get_matched_element_qvalue

  Gets the qvalue member from a cisml matched_element object.
**********************************************************************/
double get_matched_element_qvalue(MATCHED_ELEMENT_T* element) {
  assert(element != NULL);
  return element->qvalue;
}

/**********************************************************************
  set_matched_element_clusterid

  Sets the clusterid member in a cisml matched_element object.
**********************************************************************/
void set_matched_element_clusterid(
  MATCHED_ELEMENT_T *element,
  char *clusterid
) {

  assert(element != NULL);

  if (clusterid == NULL) {
    if (element->clusterid != NULL) {
      myfree(element->clusterid);
    }
    element->clusterid = NULL;
  }
  else {
    int old_length = 0;
    int new_length = strlen(clusterid) + 1;
    if (element->clusterid != NULL) {
      old_length = strlen(element->clusterid) + 1;
    }
    if (old_length < new_length) {
      element->clusterid =
        mm_realloc(element->clusterid, new_length);
    }
    strncpy(element->clusterid, clusterid, new_length);
  }

}

/**********************************************************************
  get_matched_element_clusterid

  Gets the clusterid member from a cisml matched_element object.
  May return NULL.
**********************************************************************/
char *get_matched_element_clusterid(MATCHED_ELEMENT_T* element) {

  assert(element != NULL);

  return element->clusterid;

}

/**********************************************************************
  get_matched_element_sequence

  Gets the sequence member from a cisml matched_element object.
  Caller is responsible for freeing the returned string.
**********************************************************************/
const char *get_matched_element_sequence(MATCHED_ELEMENT_T* element) {

  assert(element != NULL);

  return (const char *) element->sequence;

}

/**********************************************************************
  set_matched_element_sequence

  Sets the sequence member in a cisml matched_element object.
**********************************************************************/
void set_matched_element_sequence(MATCHED_ELEMENT_T* element, char *seq) {

  assert(element != NULL);

  if (element->sequence != NULL) {
    myfree(element->sequence);
  }
  element->sequence = strdup(seq);

}

/**********************************************************************
  set_matched_element_strand

  Sets the strand member in a cisml matched_element object.
**********************************************************************/
void set_matched_element_strand(MATCHED_ELEMENT_T* element, char strand) {

  assert(element != NULL);

  element->strand = strand;

}

/**********************************************************************
  get_matched_element_strand

  Gets the strand member in a cisml matched_element object.
**********************************************************************/
char get_matched_element_strand(MATCHED_ELEMENT_T* element, char strand) {

  assert(element != NULL);

  return element->strand;

}

/*********************************************************************
  print_cisml_matched_elements()

  Print the XML for the CisML matched elements under a scanned sequence.
**********************************************************************/
void print_cisml_matched_elements(
  CISML_T *cisml,
  FILE *out,
  int num_matched_elements,
  MATCHED_ELEMENT_T **elements
) {

  double qthresh = get_cisml_site_qvalue_cutoff(cisml);
  double pthresh = get_cisml_site_pvalue_cutoff(cisml);
  BOOLEAN_T have_output_sequence_start = FALSE;

  // We're going to output at least one matched element
  // for this sequence. Output the scanned sequence
  // starting XML tag.
  if (!have_output_sequence_start) {
    have_output_sequence_start = TRUE;
    // print_cisml_scanned_sequence_start(cisml, out, element->parent_sequence);
  }

  int i_element = 0;
  for(i_element = 0; i_element < num_matched_elements; i_element++) {

    MATCHED_ELEMENT_T *element = elements[i_element];

    if (element->pvalue > pthresh || element->qvalue > qthresh) {
      continue;
    }

    fprintf(
      out,
      "<matched-element start=\"%d\" stop=\"%d\"",
      get_matched_element_start(element),
      get_matched_element_stop(element)
    );
    if (has_matched_element_score(element)) {
      double score = get_matched_element_score(element);
      fprintf(out, " score=\"%g\"", score);
    }
    if (has_matched_element_pvalue(element)) {
      double pvalue = get_matched_element_pvalue(element);
      fprintf(out, " pvalue=\"%.3g\"", pvalue);
    }
    char *clusterid = get_matched_element_clusterid(element);
    if (clusterid != NULL) {
      fprintf(out, " clusterid=\"%s\"", clusterid);
    }
    fprintf(out, ">\n");
    const char *sequence = get_matched_element_sequence(element);
    if (sequence !=  NULL) {
      fprintf(out, "<sequence>%s</sequence>\n", sequence);
    }
    if (has_matched_element_qvalue(element)) {
      double qvalue = get_matched_element_qvalue(element);
      fprintf(out, "<mem:qvalue>%.3g</mem:qvalue>\n", qvalue);
    }
    fputs("</matched-element>\n", out);

  }

  if (have_output_sequence_start) {
    // We output the scanned sequence start XML tag.
    // Output the scanned sequence closing tag.
    // print_cisml_scanned_sequence_end(out);
  }

}

/**********************************************************************
  print_cisml_scanned_sequence_start()

  Print the starting XML tag for a CisML scanned sequence element.
**********************************************************************/
void print_cisml_scanned_sequence_start(
  CISML_T *cisml,
  FILE *out,
  SCANNED_SEQUENCE_T *seq
) {
    fprintf(
      out,
      "<scanned-sequence accession=\"%s\" name=\"%s\"",
      get_scanned_sequence_accession(seq),
      get_scanned_sequence_name(seq)
    );
    if (has_scanned_sequence_score(seq)) {
      double score = get_scanned_sequence_score(seq);
      fprintf(out, " score=\"%g\"", score);
    }
    if (has_scanned_sequence_pvalue(seq)) {
      double pvalue = get_scanned_sequence_pvalue(seq);
      fprintf(out, " pvalue=\"%g\"", pvalue);
    }
    if (has_scanned_sequence_length(seq)) {
      int length = get_scanned_sequence_length(seq);
      fprintf(out, " length=\"%d\"", length);
    }
    char *db = get_scanned_sequence_db(seq);
    if (db != NULL) {
      fprintf(out, " db=\"%s\"", db);
    }
    char *lsid = get_scanned_sequence_lsid(seq);
    if (lsid != NULL) {
      fprintf(out, " lsid=\"%s\"", lsid);
    }
    fprintf(out, ">\n");

}

/**********************************************************************
  print_cisml_scanned_sequence_end()

  Print the ending XML tag for a CisML scanned sequence element.
**********************************************************************/
void print_cisml_scanned_sequence_end(FILE *out) {

    fputs("</scanned-sequence>\n", out);

}

/**********************************************************************
  print_cisml_scanned_sequences()

  Print the XML for a CisML scanned sequence element.
**********************************************************************/
void print_cisml_scanned_sequences(
  CISML_T *cisml,
  FILE *out,
  int num_seqs,
  SCANNED_SEQUENCE_T **sequences
) {

  int i_seq = 0;
  for(i_seq = 0; i_seq < num_seqs; i_seq++) {

    SCANNED_SEQUENCE_T *seq = sequences[i_seq];

    print_cisml_scanned_sequence_start(cisml, out, seq);

    // Skip sequences with no matched elements
    if (seq->num_matched_elements == 0) {
      print_cisml_scanned_sequence_end(out);
      continue;
    }

    // We don't know if any matched elements pass the
    // p-value or q-value  threshold, so the first
    // matched elemenet printed will trigger the output
    // for the scanned sequence tags.
    print_cisml_matched_elements(
      cisml,
      out,
      seq->num_matched_elements,
      seq->elements
    );

    print_cisml_scanned_sequence_end(out);

  }

}

/**********************************************************************
  print_cisml_start_pattern

  Print the starting tag for a CisML pattern
**********************************************************************/
void print_cisml_start_pattern(
  CISML_T *cisml,
  FILE *out,
  PATTERN_T *pattern
) {

  fprintf(
    out,
    "<pattern accession=\"%s\" name=\"%s\"",
    get_pattern_accession(pattern),
    get_pattern_name(pattern)
  );
  if (has_pattern_score(pattern)) {
    double score = get_pattern_score(pattern);
    fprintf(out, " score=\"%g\"", score);
  }
  if (has_pattern_pvalue(pattern)) {
    double pvalue = get_pattern_pvalue(pattern);
    fprintf(out, " pvalue=\"%g\"", pvalue);
  }
  char *db = get_pattern_db(pattern);
  if (db != NULL) {
    fprintf(out, " db=\"%s\"", db);
  }
  char *lsid = get_pattern_lsid(pattern);
  if (lsid != NULL) {
    fprintf(out, " lsid=\"%s\"", lsid);
  }
  fputs(">\n", out);

}

/**********************************************************************
  print_cisml_end_pattern

  Print the ending tag for a CisML pattern
**********************************************************************/
void print_cisml_end_pattern(FILE *out) {
  fputs("</pattern>\n", out);
}

/**********************************************************************
  print_cisml_patterns

  Print pattern elements for CisML
**********************************************************************/
void print_cisml_patterns(
  CISML_T *cisml,
  FILE *out,
  int num_patterns,
  PATTERN_T **patterns
) {

  int i = 0;
  for(i = 0; i < num_patterns; i++) {
    int num_seq = 0;
    num_seq = get_pattern_num_scanned_sequences(patterns[i]);
    // only patterns with sequences can be printed to be conform with the DTD
    if (num_seq > 0) {
    	print_cisml_start_pattern(cisml, out, patterns[i]);
    	SCANNED_SEQUENCE_T **sequences
        = get_pattern_scanned_sequences(patterns[i]);
      print_cisml_scanned_sequences(cisml, out, num_seq, sequences);

      if (has_pattern_qvalues(patterns[i]) == TRUE) {
    	  fputs("<mem:has-qvalues>yes</mem:has-qvalues>\n", out);
      }

      print_cisml_end_pattern(out);
    }
  }
}

/**********************************************************************
  print_cisml_multi_patterns

  Print multi_patterns element for CisML
**********************************************************************/
void print_cisml_multi_patterns(
  CISML_T *cisml,
  FILE *out,
  int num_multi_patterns,
  MULTI_PATTERN_T **multi_patterns
) {

  int i = 0;
  for(i = 0; i < num_multi_patterns; i++) {
    fprintf(out, "<multi-pattern-scan");
    if (has_multi_pattern_score(multi_patterns[i])) {
      double score = get_multi_pattern_score(multi_patterns[i]);
      fprintf( out, " score=\"%g\"", score);
    }
    if (has_multi_pattern_pvalue(multi_patterns[i])) {
      double pvalue = get_multi_pattern_pvalue(multi_patterns[i]);
      fprintf( out, " pvalue=\"%g\"", pvalue);
    }
    fprintf(out, ">\n");
    int num_patterns = get_multi_pattern_num_patterns(multi_patterns[i]);
    if (num_patterns > 0) {
      PATTERN_T **patterns = get_multi_pattern_patterns(multi_patterns[i]);
      print_cisml_patterns(cisml, out, num_patterns, patterns);
    }
    fprintf(out, "</multi-pattern-scan>\n");
  }
}

/**********************************************************************
  print_cisml_parmeters

  Print parameters element for CisML
**********************************************************************/
void print_cisml_parameters(FILE *out, CISML_T *cisml) {

  char *sequence_filename = get_cisml_sequence_file(cisml);
  char *pattern_filename = get_cisml_pattern_file(cisml);

  fprintf(out, "<parameters>\n");
  fprintf(
    out,
    "<pattern-file>%s</pattern-file>\n",
    get_cisml_pattern_file(cisml)
  );
  fprintf(
    out,
    "<sequence-file>%s</sequence-file>\n",
    sequence_filename
  );
  char *background_file = get_cisml_background_file(cisml);
  if (background_file != NULL) {
    fprintf(
      out,
      "<background-seq-file>%s</background-seq-file>\n",
      background_file
    );
  }
  if (has_cisml_pattern_pvalue_cutoff(cisml)) {
    double pvalue = get_cisml_pattern_pvalue_cutoff(cisml);
    fprintf(
      out,
      "<pattern-pvalue-cutoff>%g</pattern-pvalue-cutoff>\n",
      pvalue
    );
  }
  if (has_cisml_sequence_pvalue_cutoff(cisml)) {
    double pvalue = get_cisml_sequence_pvalue_cutoff(cisml);
    fprintf(
      out,
      "<sequence-pvalue-cutoff>%g</sequence-pvalue-cutoff>\n",
      pvalue
    );
  }
  if (has_cisml_site_pvalue_cutoff(cisml)) {
    double pvalue = get_cisml_site_pvalue_cutoff(cisml);
    fprintf(
      out,
      "<site-pvalue-cutoff>%g</site-pvalue-cutoff>\n",
      pvalue
    );
  }
  char *filter = get_cisml_sequence_filter(cisml);
  if (filter != NULL) {
    fprintf(
      out,
      "<sequence-filtering on-off=\"on\" type=\"%s\" />\n",
      filter
    );
  }
  else {
    fprintf(
      out,
      "<sequence-filtering on-off=\"off\"/>\n"
    );
  }

  // Print MEME extensions for CisML


  fprintf(out, "</parameters>\n" );
}

/**********************************************************************
  print_cisml_xml_header

  Print the DTD and XML header for CisML
**********************************************************************/
static void print_cisml_xml_header(FILE* out, const char *stylesheet) {
  if (stylesheet == NULL){
	  fputs(cisml_dts,out);
  } else {
	  fputs("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n", out);
	  fputs("<?xml-stylesheet type=\"text/xsl\" ", out);
	  fprintf(out, "href=\"%s\"?>\n", stylesheet);
	  fputs("<!-- Begin document body -->\n", out);
  }
}

/**
 * register_namespaces:
 * @xpathCtx:		the pointer to an XPath context.
 * @nsList:		the list of known namespaces in
 *			"<prefix1>=<href1> <prefix2>=href2> ..." format.
 *
 * Registers namespaces from @nsList in @xpathCtx.
 *
 * Returns 0 on success and a negative value otherwise.
 */
int register_namespaces(xmlXPathContextPtr xpathCtx, const xmlChar* nsList) {
    xmlChar* nsListDup;
    xmlChar* prefix;
    xmlChar* href;
    xmlChar* next;

    assert(xpathCtx);
    assert(nsList);

    nsListDup = xmlStrdup(nsList);
    if(nsListDup == NULL) {
	fprintf(stderr, "Error: unable to strdup namespaces list\n");
	return(-1);
    }

    next = nsListDup;
    while(next != NULL) {
		/* skip spaces */
		while((*next) == ' ') next++;
		if((*next) == '\0') break;

		/* find prefix */
		prefix = next;
		next = (xmlChar*)xmlStrchr(next, '=');
		if(next == NULL) {
			fprintf(stderr,"Error: invalid namespaces list format\n");
			xmlFree(nsListDup);
			return(-1);
		}
		*(next++) = '\0';

		/* find href */
		href = next;
		next = (xmlChar*)xmlStrchr(next, ' ');
		if(next != NULL) {
			*(next++) = '\0';
		}

		fprintf(stderr,"%s %s\n",prefix,href);
		if(xmlXPathRegisterNs(xpathCtx, prefix, href) != 0) {
			fprintf(stderr,"Error: unable to register NS with prefix=\"%s\" and href=\"%s\"\n", prefix, href);
			xmlFree(nsListDup);
			return(-1);
		}
    }

    xmlFree(nsListDup);
    return(0);
}


/**********************************************************************
  read_cisml

  Reads in a CisML XML file and create the cisml_t data structure
  Attention: only the standard cisml file format is supported so far
  			 extension as <mem:qvalue> etc are not supported
  			 thus resulting in an error when a file containing such
  			 lines is given
**********************************************************************/
CISML_T* read_cisml(
	char* 		cisml_filename // name of the cisml file
){
  xmlParserCtxtPtr ctxt = NULL;         // The parser context
  xmlDocPtr cisml_doc = NULL;            // The resulting document tree
  xmlXPathContextPtr xpath_ctxt = NULL; // XPath context.

  ctxt = xmlNewParserCtxt();
  if (ctxt == NULL) {
	  die("Failed to create XML parser.\n");
  }

  // Parse and validate the file.
  cisml_doc = xmlCtxtReadFile(
		  ctxt,
		  cisml_filename,
		  NULL,  // Encoding
		  XML_PARSE_DTDVALID | XML_PARSE_NOWARNING
  );

  // Did it parse?
  if (cisml_doc == NULL) {
	  fprintf(stderr, "Failed to parse %s as cisml XML document.\n", cisml_filename);
	  xmlFreeParserCtxt(ctxt);
	  xmlCleanupParser();
	  return FALSE;
  } else {
	  // Did it validate?
	  if (ctxt->valid == 0) {
		  die("%s is not a valid cisml XML document.\n", cisml_filename);
	  }
  }
  if (verbosity >= HIGH_VERBOSE) {
	  fprintf(stderr, "File %s is a valid cisml XML file.\n", cisml_filename);
  }

  // Set up XPath context from parsed XML
  xpath_ctxt = xmlXPathNewContext(cisml_doc);

  xmlChar* prefix = BAD_CAST "mem";
  xmlChar* href = BAD_CAST "http://noble.gs.washington.edu/meme";
  xmlXPathRegisterNs(xpath_ctxt, prefix, href);

  // read program name
  xmlXPathObjectPtr xpath_obj = xpath_query(xpath_ctxt, "/cis-element-search/program-name");
  assert(xpath_obj->nodesetval->nodeNr == 1);
  xmlNodePtr currValueNode = xpath_obj->nodesetval->nodeTab[0];
  xmlXPathFreeObject(xpath_obj);
  char* program_name = (char *) xmlXPathCastNodeToString(currValueNode);

  // read meme file
  xpath_obj = xpath_query(xpath_ctxt, "/cis-element-search/parameters/pattern-file");
  assert(xpath_obj->nodesetval->nodeNr == 1);
  currValueNode = xpath_obj->nodesetval->nodeTab[0];
  xmlXPathFreeObject(xpath_obj);
  char* meme_filename = (char *) xmlXPathCastNodeToString(currValueNode);

  // read fasta file
  xpath_obj = xpath_query(xpath_ctxt, "/cis-element-search/parameters/sequence-file");
  assert(xpath_obj->nodesetval->nodeNr == 1);
  currValueNode = xpath_obj->nodesetval->nodeTab[0];
  xmlXPathFreeObject(xpath_obj);
  char* fasta_filename = (char *) xmlXPathCastNodeToString(currValueNode);

  // Create cisml data structure for recording results
  CISML_T* cisml = allocate_cisml(program_name, meme_filename, fasta_filename);

  // read pattern cutoff
  xpath_obj = xpath_query(xpath_ctxt, "/cis-element-search/parameters/pattern-pvalue-cutoff");
  if (xpath_obj->nodesetval->nodeNr == 1){
    currValueNode = xpath_obj->nodesetval->nodeTab[0];
    xmlXPathFreeObject(xpath_obj);
    double cutoff = xmlXPathCastNodeToNumber(currValueNode);
    if (!xmlXPathIsNaN(cutoff))
    	set_cisml_pattern_pvalue_cutoff(cisml, cutoff);
  }

  // read sequence cutoff
  xpath_obj = xpath_query(xpath_ctxt, "/cis-element-search/parameters/sequence-pvalue-cutoff");
  if (xpath_obj->nodesetval->nodeNr == 1){
    currValueNode = xpath_obj->nodesetval->nodeTab[0];
    xmlXPathFreeObject(xpath_obj);
    double cutoff = xmlXPathCastNodeToNumber(currValueNode);
    if (!xmlXPathIsNaN(cutoff))
    	set_cisml_sequence_pvalue_cutoff(cisml, cutoff);
  }

  // read site cutoff
  xpath_obj = xpath_query(xpath_ctxt, "/cis-element-search/parameters/site-pvalue-cutoff");
  if (xpath_obj->nodesetval->nodeNr == 1){
	  currValueNode = xpath_obj->nodesetval->nodeTab[0];
	  xmlXPathFreeObject(xpath_obj);
	  double cutoff = xmlXPathCastNodeToNumber(currValueNode);
	  if (!xmlXPathIsNaN(cutoff))
	  	set_cisml_site_pvalue_cutoff(cisml, cutoff);
  }

  int num_patterns = 0;
  read_patterns_from_cisml(xpath_ctxt, &num_patterns, cisml);

  // either pattern or multiple-pattern scan is possible
  if (num_patterns == 0){
	  read_multiple_patterns_from_cisml(xpath_ctxt, &num_patterns, cisml);
  }

  /*free the document */
  xmlFreeDoc(cisml_doc);

  /* Free the global variables that may have been allocated by the parser. */
  xmlXPathRegisteredNsCleanup(xpath_ctxt);
  xmlCleanupParser();

  return (cisml);
}

/**********************************************************************
  read_multiple_patterns_from_cisml

  reads a multi-pattern-scan from a CisML xml file
**********************************************************************/
void read_multiple_patterns_from_cisml(
		xmlXPathContextPtr 	xpath_ctxt,    	// cisml XPath context.
		int* 				num_patterns,	// cisml number of pattern in the file
		CISML_T* 			cisml			// the retrieved cisml struct
) {
  xmlXPathObjectPtr xpath_obj = NULL;
  xmlChar* property = NULL;
  // Use XPATH to get the set of patterns.
  xpath_obj = xpath_query(xpath_ctxt, "/cis-element-search/multi-pattern-scan");
  *num_patterns = (xpath_obj->nodesetval ? xpath_obj->nodesetval->nodeNr : 0);

  if (verbosity >= NORMAL_VERBOSE) {
	  fprintf(stderr, "Reading multi-pattern-scan.");
  }
  xmlNodePtr currMultiPatternNode = NULL;
  int i = 0;
  for (i = 0; i < *num_patterns; i++) {

	  currMultiPatternNode = xpath_obj->nodesetval->nodeTab[i];
  	  if (currMultiPatternNode == NULL) {
  		  die("Error: missing pattern %d\n", i);
  	  }

  	  MULTI_PATTERN_T* multipattern = allocate_multi_pattern();

  	  // Get the multi pattern pvalue
  	  if (check_xml_node_property(currMultiPatternNode, "pvalue")){
		  property = read_xml_node_property(currMultiPatternNode, "pvalue");
		  set_multi_pattern_pvalue(multipattern, atof((char *) property));
		  xmlFree(property);
  	  }

	  // Get the multi pattern score
  	  if (check_xml_node_property(currMultiPatternNode, "score")){
		  property = read_xml_node_property(currMultiPatternNode, "score");
		  set_multi_pattern_score(multipattern, atof((char *) property));
		  xmlFree(property);
  	  }

  	  xmlNodePtr currPatternNode = currMultiPatternNode->children;
  	  while (currPatternNode != NULL){
  		  if (currPatternNode->type == XML_ELEMENT_NODE){
	  		  // retrieve the pattern
	  		  PATTERN_T* pattern =  get_pattern_from_xmlnode(currPatternNode);
	  		  add_multi_pattern_pattern(multipattern, pattern);
  		  }
  		  currPatternNode = currPatternNode->next;
  	  }
  	  add_cisml_multi_pattern(cisml, multipattern);
  }

}
/**********************************************************************
  read_patterns_from_cisml

  reads the patterns from a CisML xml file
**********************************************************************/
void read_patterns_from_cisml(
		xmlXPathContextPtr 	xpath_ctxt,    	// cisml XPath context.
		int* 				num_patterns,	// cisml number of pattern in the file
		CISML_T* 			cisml			// the retrieved cisml struct
) {
  xmlXPathObjectPtr xpath_obj = NULL;
  // Use XPATH to get the set of patterns.
  xpath_obj = xpath_query(xpath_ctxt, "/cis-element-search/pattern");

  *num_patterns = (xpath_obj->nodesetval ? xpath_obj->nodesetval->nodeNr : 0);

  if (verbosity >= HIGH_VERBOSE) {
    fprintf(stderr, "Reading patterns ... ");
  }

  xmlNodePtr currPatternNode = NULL;
  int i_pattern = 0;
  for (i_pattern = 0; i_pattern < *num_patterns; i_pattern++) {
	  // get the pattern node
	  currPatternNode = xpath_obj->nodesetval->nodeTab[i_pattern];
	  // retrieve the pattern
	  PATTERN_T* pattern = get_pattern_from_xmlnode(currPatternNode);
	  // add pattern
	  add_cisml_pattern(cisml, pattern);
	  // get scanned sequences
	  read_scanned_sequences_from_cisml(cisml, pattern, currPatternNode->children);

	  // check if at least one sequence has been assigned
	  int num_seq = 0;
	  num_seq = get_pattern_num_scanned_sequences(pattern);
	  if (num_seq == 0){
		  die("CisML file not valid! At least one sequence needs to be assigned to each pattern.");
	  }
  }

  if (verbosity >= HIGH_VERBOSE) {
    fprintf(stderr, "finished\n");
  }
}

/**********************************************************************
  get_pattern_from_xmlnode

  returns a pattern from a specified xml Node
**********************************************************************/
PATTERN_T* get_pattern_from_xmlnode(xmlNodePtr currPatternNode){
  xmlChar* property = NULL;
  xmlChar* property2 = NULL;

  if (currPatternNode == NULL || strcmp((char *) currPatternNode->name,"pattern")!=0) {
	  die("Error: missing pattern or pattern not valid\n");
  }

  // Get the pattern accession attribute
  property = read_xml_node_property(currPatternNode, "accession");
  property2 = read_xml_node_property(currPatternNode, "name");
  PATTERN_T *pattern = allocate_pattern((char *) property, (char *) property2);

  xmlFree(property);
  xmlFree(property2);

  // Get the pattern pvalue
  if (check_xml_node_property(currPatternNode, "pvalue")){
	  property = read_xml_node_property(currPatternNode, "pvalue");
	  set_pattern_pvalue(pattern, atof((char *) property));
	  xmlFree(property);
  }

  // Get the pattern score
  if (check_xml_node_property(currPatternNode, "score")){
	  property = read_xml_node_property(currPatternNode, "score");
	  set_pattern_score(pattern, atof((char *) property));
	  xmlFree(property);
  }

  // Get the pattern db
  if (check_xml_node_property(currPatternNode, "db")){
	  property = read_xml_node_property(currPatternNode, "db");
	  set_pattern_db(pattern, (char *) property);
	  xmlFree(property);
  }

  // Get the pattern lsid
  if (check_xml_node_property(currPatternNode, "lsid")){
	  property = read_xml_node_property(currPatternNode, "lsid");
	  set_pattern_lsid(pattern, (char *) property);
	  xmlFree(property);
  }
  return pattern;
}

/**********************************************************************
  read_scanned_sequences_from_cisml

  Reads in the sequences for a pattern
**********************************************************************/
void read_scanned_sequences_from_cisml(
  CISML_T* 				cisml,			// the retrieved cisml struct
  PATTERN_T*			pattern, 		// the pattern of interest
  xmlNode*				first_sequence	// first child of the pattern node
) {

  xmlChar* property = NULL;
  xmlChar* property2 = NULL;

  // add all sequences to the pattern
  xmlNode *cur_node = NULL;
  for (cur_node = first_sequence; cur_node; cur_node = cur_node->next) {
	  if (cur_node->type != XML_ELEMENT_NODE) {
		  continue;
	  }

//	  printf("node type: Element, name: %s %s \n", cur_node->name,read_xml_node_property(cur_node, "accession"));
	  property = read_xml_node_property(cur_node, "accession");
	  property2 = read_xml_node_property(cur_node, "name");
	  SCANNED_SEQUENCE_T *scanned_seq =
		  allocate_scanned_sequence((char *) property, (char *) property2, pattern);

	  xmlFree(property);
	  xmlFree(property2);

	  // Get the pattern length
	  if (check_xml_node_property(cur_node, "length")){
		  property = read_xml_node_property(cur_node, "length");
		  set_scanned_sequence_length(scanned_seq, atoi((char *) property));
		  xmlFree(property);
	  }

	  // Get the pattern pvalue
	  if (check_xml_node_property(cur_node, "pvalue")){
	  	  property = read_xml_node_property(cur_node, "pvalue");
	  	  set_scanned_sequence_pvalue(scanned_seq, atof((char *) property));
	  	  xmlFree(property);
	  }

  	  // Get the pattern score
	  if (check_xml_node_property(cur_node, "score")){
	  	  property = read_xml_node_property(cur_node, "score");
	  	  set_scanned_sequence_score(scanned_seq, atof((char *) property));
	  	  xmlFree(property);
	  }

  	  // Get the pattern db
  	  if (check_xml_node_property(cur_node, "db")){
  		  property = read_xml_node_property(cur_node, "db");
  		  set_scanned_sequence_db(scanned_seq, (char *) property);
  		  xmlFree(property);
  	  }

  	  // Get the pattern lsid
  	  if (check_xml_node_property(cur_node, "lsid")){
  		  property = read_xml_node_property(cur_node, "lsid");
  		  set_scanned_sequence_lsid(scanned_seq, (char *) property);
  		  xmlFree(property);
  	  }

  	  // get matches
  	  read_matched_elements_from_cisml(cisml, scanned_seq, cur_node->children);
  }
}

/**********************************************************************
  read_matched_elements_from_cisml

  Reads in sequences matches from a cisml file
**********************************************************************/
void   read_matched_elements_from_cisml(
  CISML_T* 				cisml,			// the retrieved cisml struct
  SCANNED_SEQUENCE_T*	scanned_seq,	// the pattern of interest
  xmlNode*				firstmatch		// first child of the pattern node
) {

  // add all sequences to the pattern
  xmlNode *cur_node = NULL;
  for (cur_node = firstmatch; cur_node; cur_node = cur_node->next) {
	  if (cur_node->type != XML_ELEMENT_NODE) {
		  continue;
	  }

	  xmlChar* property = NULL;
	  xmlChar* window_seq = NULL;

	  int start = atoi((const char *) read_xml_node_property(cur_node, "start"));
	  int stop  = atoi((const char *) read_xml_node_property(cur_node, "stop"));

	  xmlNode *child_node = cur_node->children;
	  while (child_node != NULL){
		  if (strcmp((const char *) child_node->name, (const char *) "sequence")==0){
			  if (child_node->children->type == XML_TEXT_NODE){
				  window_seq = child_node->children->content;
				  child_node = NULL;
			  } else {
				  die("CisML document invalid!");
			  }
		  } else {
			  child_node = child_node->next;
		  }
	  }
	  MATCHED_ELEMENT_T *element =
		  allocate_matched_element(start, stop, scanned_seq);
    set_matched_element_sequence(element, (char *) window_seq);

	  // Get the match pvalue
	  if (check_xml_node_property(cur_node, "pvalue")){
	  	  property = read_xml_node_property(cur_node, "pvalue");
	  	  set_matched_element_pvalue(element, atof((char *) property));
	  	  xmlFree(property);
	  }

  	  // Get the match score
	  if (check_xml_node_property(cur_node, "score")){
	  	  property = read_xml_node_property(cur_node, "score");
	  	  set_matched_element_score(element, atof((char *) property));
	  	  xmlFree(property);
	  }

    // Get the match clusterid
    if (check_xml_node_property(cur_node, "clusterid")){
      property = read_xml_node_property(cur_node, "clusterid");
      set_matched_element_clusterid(element, (char *) property);
      xmlFree(property);
    }

  }
}

/**********************************************************************
  print_cisml_start

  Print the opening section of the CisML XML
**********************************************************************/
void print_cisml_start(FILE* out, CISML_T *cisml, BOOLEAN_T print_header,
		const char *stylesheet, BOOLEAN_T print_namespace) {

  if (print_header == TRUE) {
    print_cisml_xml_header(out, stylesheet);
  }
  fputs("<cis-element-search\n", out);
  fputs("  xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"", out);
  fputs("\n", out);
  fputs("  xsi:schemaLocation=", out);
  fputs("\"http://zlab.bu.edu/schema/cisml cisml.xsd\"\n", out);
  // FIXME the next line causes problems for the cisml reader
  if (print_namespace){
	  fputs("  xmlns=\"http://zlab.bu.edu/schema/cisml\"\n", out);
  }
  fputs("  xmlns:mem=\"http://noble.gs.washington.edu/meme\"\n>\n", out);

  fprintf(
    out,
    "<program-name>%s</program-name>\n",
    get_cisml_program_name(cisml)
  );
}

/**********************************************************************
  print_cisml_end

  Print the closing section of the CisML XML
**********************************************************************/
void print_cisml_end(FILE* out) {
  fprintf(out, "</cis-element-search>\n");
}

/**********************************************************************
  print_cisml

  Print the cisml data structure as CisML XML
**********************************************************************/
void print_cisml(FILE* out, CISML_T *cisml, BOOLEAN_T print_header,
		const char *stylesheet, BOOLEAN_T print_namespace) {

  print_cisml_start(out, cisml, print_header, stylesheet, print_namespace);
  print_cisml_parameters(out, cisml);

  int num_multi_patterns = get_cisml_num_multi_patterns(cisml);
  if (num_multi_patterns > 0) {
    MULTI_PATTERN_T **multi_patterns = get_cisml_multi_patterns(cisml);
    print_cisml_multi_patterns(cisml, out, num_multi_patterns, multi_patterns);
  }

  int num_patterns = get_cisml_num_patterns(cisml);
  if (num_patterns > 0) {
    PATTERN_T **patterns = get_cisml_patterns(cisml);
    print_cisml_patterns(cisml, out, num_patterns, patterns);
  }

  print_cisml_end(out);

}

/*************************************************************************
 * Compare two matched-elements by pvalue for 'qsort'.
 *************************************************************************/
static int matched_elements_compare_by_pvalue
  (const void* elem1,
   const void* elem2)
{
  const double key1 = ((MATCHED_ELEMENT_T *) elem1)->pvalue;
  const double key2 = ((MATCHED_ELEMENT_T *) elem2)->pvalue;
  const int start1 = ((MATCHED_ELEMENT_T *)elem1)->start;
  const int start2 = ((MATCHED_ELEMENT_T *)elem2)->start;

  if (key1 < key2) {
    return(-1);
  } else if (key1 > key2) {
    return(1);
  } else if (start1 < start2) {
    return(-1);
  } else if (start1 > start2) {
    return(1);
  }
  return(0);
}

/*************************************************************************
 * sort_matched_elements
 *
 * Sort a an array of pointers to matched-elements sites by pvalue
 * or by sequence name and position.
 *************************************************************************/
void sort_matched_elements(
  BOOLEAN_T sort_by_pvalue,
  int num_elements,
  MATCHED_ELEMENT_T **elements
) {

  assert(elements != NULL);

  // Tell the user what's up.
  if (verbosity > NORMAL_VERBOSE) {
    fprintf(stderr, "Sorting %d matched elements ", num_elements);
    if (sort_by_pvalue) {
      fprintf(stderr, "by p-value.\n");
    } else {
      fprintf(stderr, "by sequence name and start position.\n");
    }
  }

  if (sort_by_pvalue) {
    qsort(
	    (void *) elements,
	    num_elements,
	    sizeof(MATCHED_ELEMENT_T *),
	    matched_elements_compare_by_pvalue
	  );
  } else {
    /*
    qsort(
	    (void *) elements,
	    num_elements,
	    sizeof(MATCHED_ELEMENT_T *),
	    matched_elements_compare_by_position
	  );
    */
  }
}

/*************************************************************************
 * Calculate the q-values corresponding from the the p-values of the
 * matched elements.
 *************************************************************************/
void pattern_calculate_qvalues(PATTERN_T *pattern, ARRAY_T *sampled_pvalues) {

  assert(pattern != NULL);
  assert(pattern->is_complete == TRUE);

  int num_stored_matches = get_pattern_num_stored_matches(pattern);
  long num_scanned_positions = get_pattern_num_scanned_positions(pattern);
  if (verbosity >= HIGH_VERBOSE) {
    fprintf(stderr, "Num stored matches %d\n", num_stored_matches);
    fprintf(stderr, "Num scanned positions %ld\n", num_scanned_positions);
  }
  // Tell the user what's up.
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Computing q-values.\n");
  }

  if (num_stored_matches) {

    // Extract the p-values into an array.
    int i_element = 0;
    ARRAY_T* pvalues = allocate_array(num_stored_matches);
    MATCHED_ELEMENT_T *element = NULL;
    for (i_element = 0; i_element < num_stored_matches; i_element++) {
      element = pattern->elements[i_element];
      set_array_item(i_element, element->pvalue, pvalues);
    }

    // Convert them to q-values.
    compute_qvalues(
      FALSE, // Don't stop with FDR
      TRUE, // Try to esimate pi0
      NULL, // Don't store pi-zero in a file.
      NUM_BOOTSTRAPS,
      NUM_BOOTSTRAP_SAMPLES,
      NUM_LAMBDA,
      MAX_LAMBDA,
      num_scanned_positions,
      pvalues,
      sampled_pvalues
    );

    // Update the matched elements with the q-values.
    for (i_element = 0; i_element < num_stored_matches; i_element++) {
      set_matched_element_qvalue(
        pattern->elements[i_element],
        get_array_item(i_element, pvalues)
      );
    }

    free_array(pvalues);

    // Sort by sequence ID and position.
    // Since we are putting out XML there is no need to sot
    // matched-elements back into positon order.
  }

  pattern->qvalues_computed = TRUE;

}

/**********************************************************************
 * This function saves CisML results as a set of files in a
 * directory. The file names are provided by the input parameters:
 *   cisml is a pointer to the ciml structure to be printed.
 *   ouput_dirname will be the name of the output directory
 *   xml_filename will be the CisML output
 *   html_filename will be the name of the HTML output
 *   text_filename will be the name of the plain text output
 *   gff_filename will be the name of the GFF output
 *   allow_clobber will determine whether or not existing files will
 *                 be overwritten.
 *   print_namespace will determine whether or not the standard name
 *   				space will be written
 *********************************************************************/
void print_full_results(
  CISML_T *cisml,
  char *output_dirname,
  char *xml_filename,
  char *html_filename,
  char *text_filename,
  char *gff_filename,
  BOOLEAN_T allow_clobber,
  BOOLEAN_T print_namespace
) {

  // Create output directory
  if (create_output_directory(
       output_dirname,
       allow_clobber,
       FALSE /* Don't print warning messages */
      )
    ) {
    // Failed to create output directory.
    die("Couldn't create output directory %s.\n", output_dirname);
  }
  // Create the paths to the output files
  const char* HTML_STYLESHEET = "cisml-to-html.xsl";
  const char* CSS_STYLESHEET = "cisml.css";
  const char* GFF_STYLESHEET = "cisml-to-gff.xsl";
  const char* TEXT_STYLESHEET = "cisml-to-text.xsl";
  char *html_stylesheet_path = make_path_to_file(ETC_DIR, HTML_STYLESHEET);
  char *css_stylesheet_path = make_path_to_file(ETC_DIR, CSS_STYLESHEET);
  char *text_stylesheet_path = make_path_to_file(ETC_DIR, TEXT_STYLESHEET);
  char *gff_stylesheet_path = make_path_to_file(ETC_DIR, GFF_STYLESHEET);
  char *xml_path = make_path_to_file(output_dirname, xml_filename);
  char *html_path = make_path_to_file(output_dirname, html_filename);
  char *text_path = make_path_to_file(output_dirname, text_filename);
  char *gff_path = make_path_to_file(output_dirname, gff_filename);
  char *html_stylesheet_copy_path = make_path_to_file(output_dirname, HTML_STYLESHEET);
  char *css_stylesheet_copy_path = make_path_to_file(output_dirname, CSS_STYLESHEET);
  FILE *xml_file = fopen(xml_path, "w");
  if (!xml_file) {
    die("Couldn't open file %s for output.\n", xml_path);
  }

  // Output XML
  print_cisml(xml_file, cisml, TRUE, HTML_STYLESHEET,print_namespace);
  fclose(xml_file);

  // Output HTML
  print_xml_filename_to_filename_using_stylesheet(xml_path, html_stylesheet_path, html_path);

  // Output text
  print_xml_filename_to_filename_using_stylesheet(xml_path, text_stylesheet_path, text_path);

  // Output GFF
  print_xml_filename_to_filename_using_stylesheet(xml_path, gff_stylesheet_path, gff_path);

  // Copy XML to HTML and CSS stylesheets to output directory
  copy_file(html_stylesheet_path, html_stylesheet_copy_path);
  copy_file(css_stylesheet_path, css_stylesheet_copy_path);

  myfree(html_stylesheet_path);
  myfree(html_stylesheet_copy_path);
  myfree(css_stylesheet_path);
  myfree(css_stylesheet_copy_path);
  myfree(text_stylesheet_path);
  myfree(gff_stylesheet_path);
  myfree(xml_path);
  myfree(html_path);
  myfree(text_path);
  myfree(gff_path);
}

/**********************************************************************
 * This function saves the CisML results as a plain text to stdout.
 *********************************************************************/
void print_cisml_as_text(CISML_T *cisml) {

  const char* etc_dir = get_meme_etc_dir();

  // Create temp file for CisML output.
  char tmp_filename[] = "CISMLXXXXXX";
  int fd = mkstemp(tmp_filename);
  if (fd == -1) {
    die("Couldn't create temporary file for text results\n");
  }
  FILE *xml_file = fdopen(fd, "w");
  if (!xml_file) {
    die("Couldn't open file %s for output.\n", "fimo.xml");
  }

  // Output CisML
  print_cisml(xml_file, cisml, FALSE, NULL, TRUE);

  fclose(xml_file);


  // Output text using stylesheet that converts CisML temp file to text.
  const char* TEXT_STYLESHEET = "cisml-to-text.xsl";
  char *text_stylesheet_path = make_path_to_file(etc_dir, TEXT_STYLESHEET);
  print_xml_filename_to_file_using_stylesheet(tmp_filename, text_stylesheet_path, stdout);


  myfree(text_stylesheet_path);
  close(fd);
  int result = remove(tmp_filename);
  if (result == -1) {
    fprintf(stderr, "Couldn't remove temporary file %s.\n", tmp_filename);
  }

}

