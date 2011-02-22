/**********************************************************************
 * FILE: cisml.h
 * AUTHOR: Charles Grant
 * PROJECT: MEME
 * COPYRIGHT: 2007-2008, WSN
 * VERSION: $Revision: 1.1.1.1$
 * DESCRIPTION: Print routines for the cisml format.
 *
 * "CisML: an XML-based format for sequence motif detection software."
 * Haverty and Weng.  Bioinformatics.  20(11):1815-1817, 2004.
 *
 * http://bioinformatics.oxfordjournals.org/cgi/reprint/20/11/1815.pdf
 * http://zlab.bu.edu/~phaverty/SupplementalData/CisML/
 *
 * In CisML, the parent-child relationships go like this:
 *
 *  [ mulit-pattern -> ] pattern -> scanned-sequence -> matched-element
 *
 **********************************************************************/
#ifndef CISML_H
#define CISML_H

#include "utils.h"
#include "xml-util.h"

typedef struct cisml CISML_T;
typedef struct multi_pattern MULTI_PATTERN_T;
typedef struct pattern PATTERN_T;
typedef struct scanned_sequence SCANNED_SEQUENCE_T;
typedef struct matched_element MATCHED_ELEMENT_T;

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
);

/**********************************************************************
  free_cisml

  Destructor for the cisml data structure.
**********************************************************************/
void free_cisml(CISML_T *cisml);

/**********************************************************************
  get_cisml_program_name

  Gets the program_name member from a cisml object.
**********************************************************************/
char *get_cisml_program_name(CISML_T *cisml);

/**********************************************************************
  get_cisml_pattern_file

  Gets the pattern_file member from a cisml object.
**********************************************************************/
char *get_cisml_pattern_file(CISML_T *cisml);

/**********************************************************************
  get_cisml_sequence_file

  Gets the sequence_file member from a cisml object.
**********************************************************************/
char *get_cisml_sequence_file(CISML_T *cisml);

/**********************************************************************
  set_cisml_background_file

  Sets the background_file member in a cisml object.
**********************************************************************/
void set_cisml_background_file(CISML_T *cisml, char *background_file);

/**********************************************************************
  get_cisml_background_file

  Gets the background_file member from a cisml object.
  Return value may be NULL.
**********************************************************************/
char *get_cisml_background_file(CISML_T *cisml);

/**********************************************************************
  set_cisml_pattern_pvalue_cutoff

  Sets the pattern_pvalue_cutoff member in a cisml object.
  Defaults to 1.0 if not set.
**********************************************************************/
void set_cisml_pattern_pvalue_cutoff(
  CISML_T *cisml,
  double pattern_pvalue_cutoff
);

/**********************************************************************
  clear_cisml_pattern_pvalue_cutoff

  Sets the pattern_pvalue_cutoff member in a cisml object to null.
**********************************************************************/
void clear_cisml_pattern_pvalue_cutoff(CISML_T *cisml);

/**********************************************************************
  has_cisml_pattern_pvalue_cutoff

  Does a cisml object have a pattern_pvalue_cutoff?
**********************************************************************/
BOOLEAN_T has_cisml_pattern_pvalue_cutoff(CISML_T *cisml);

/**********************************************************************
  get_cisml_pattern_pvalue_cutoff

  Gets the pattern_pvalue_cutoff member from a cisml object.
  Defaults to 1.0 if the pattern p-value cutoff has not been set.
**********************************************************************/
double get_cisml_pattern_pvalue_cutoff(CISML_T *cisml);

/**********************************************************************
  set_cisml_sequence_pvalue_cutoff

  Sets the sequence_pvalue_cutoff member in a cisml object.
**********************************************************************/
void set_cisml_sequence_pvalue_cutoff(
  CISML_T *cisml,
  double sequence_pvalue_cutoff
);

/**********************************************************************
  clear_cisml_sequence_pvalue_cutoff

  Sets the sequence_pvalue_cutoff member in a cisml object to null.
**********************************************************************/
void clear_cisml_sequence_pvalue_cutoff(CISML_T *cisml);

/**********************************************************************
  has_cisml_sequence_pvalue_cutoff

  Does a cisml object have a sequence_pvalue_cutoff?
**********************************************************************/
BOOLEAN_T has_cisml_sequence_pvalue_cutoff(CISML_T *cisml);

/**********************************************************************
  get_cisml_sequence_pvalue_cutoff

  Gets the sequence_pvalue_cutoff member from a cisml object.
  Defaults to 1.0 if the sequence p-value cutoff has not been set.
**********************************************************************/
double get_cisml_sequence_pvalue_cutoff(CISML_T *cisml);

/**********************************************************************
  set_cisml_site_pvalue_cutoff

  Sets the site_pvalue_cutoff member in a cisml object.
**********************************************************************/
void set_cisml_site_pvalue_cutoff(CISML_T *cisml, double site_pvalue_cutoff);

/**********************************************************************
  clear_cisml_site_pvalue_cutoff

  Sets the site_pvalue_cutoff member in a cisml object to null.
**********************************************************************/
void clear_cisml_site_pvalue_cutoff(CISML_T *cisml);

/**********************************************************************
  has_cisml_site_pvalue_cutoff

  Does a cisml object have a pattern_pvalue_cutoff?
**********************************************************************/
BOOLEAN_T has_cisml_site_pvalue_cutoff(CISML_T *cisml);

/**********************************************************************
  get_cisml_site_pvalue_cutoff

  Gets the site_pvalue_cutoff member from a cisml object.
  Defaults to 1.0 if the site p-value cutoff has not been set.
**********************************************************************/
double get_cisml_site_pvalue_cutoff(CISML_T *cisml);

/**********************************************************************
  set_cisml_site_qvalue_cutoff

  Sets the site_qvalue_cutoff member in a cisml object.
**********************************************************************/
void set_cisml_site_qvalue_cutoff(CISML_T *cisml, double site_qvalue_cutoff);

/**********************************************************************
  clear_cisml_site_qvalue_cutoff

  Sets the site_qvalue_cutoff member in a cisml object to null.
**********************************************************************/
void clear_cisml_site_qvalue_cutoff(CISML_T *cisml);

/**********************************************************************
  has_cisml_site_qvalue_cutoff

  Does a cisml object have a site_qvalue_cutoff?
**********************************************************************/
BOOLEAN_T has_cisml_site_qvalue_cutoff(CISML_T *cisml);

/**********************************************************************
  get_cisml_site_qvalue_cutoff

  Gets the site_qvalue_cutoff member from a cisml object.
  Defaults to 1.0 if the site q-value cutoff has not been set.
**********************************************************************/
double get_cisml_site_qvalue_cutoff(CISML_T *cisml);

/**********************************************************************
  set_cisml_sequence_filter

  Sets the sequence_filter member in a cisml object.
**********************************************************************/
void set_cisml_sequence_filter(CISML_T *cisml, char *sequence_filter);

/**********************************************************************
  get_cisml_sequence_filter

  Gets the sequence_filter member from a cisml object.
  May return NULL.
**********************************************************************/
char *get_cisml_sequence_filter(CISML_T *cisml);

/**********************************************************************
  get_cisml_num_multi_patterns

  Gets the number of multi-patterns from a cisml object.
**********************************************************************/
int get_cisml_num_multi_patterns(CISML_T *cisml);

/**********************************************************************
  get_cisml_multi_patterns

  Gets the array of pointers to multi_patterns from a cisml object.
  May return NULL.
**********************************************************************/
MULTI_PATTERN_T **get_cisml_multi_patterns(CISML_T *cisml);

/**********************************************************************
  get_cisml_num_patterns

  Gets the number of patterns from a cisml object.
**********************************************************************/
int get_cisml_num_patterns(CISML_T *cisml);

/**********************************************************************
  get_cisml_patterns

  Gets the array of pointers to patterns from a cisml object.
**********************************************************************/
PATTERN_T **get_cisml_patterns(CISML_T *cisml);

/**********************************************************************
  add_cisml_multi_pattern

  Adds a pattern to the array of pointers to patterns in a cisml object.
**********************************************************************/
void add_cisml_multi_pattern(CISML_T *cisml, MULTI_PATTERN_T* multi_pattern);

/**********************************************************************
  add_cisml_pattern

  Adds a pattern to the array of pointers to patterns in a cisml object.
**********************************************************************/
void add_cisml_pattern(CISML_T *cisml, PATTERN_T* pattern);

/**********************************************************************
  allocate_multi_pattern

  Constructor for the cisml multi_pattern data structure.
**********************************************************************/
MULTI_PATTERN_T *allocate_multi_pattern();

/**********************************************************************
  free_multi_pattern

  Destructor for the cisml multi_pattern data structure.
**********************************************************************/
void free_multi_pattern(MULTI_PATTERN_T *multi_pattern);

/**********************************************************************
  set_multi_pattern_score

  Sets the score member in a cisml pattern object.
**********************************************************************/
void set_multi_pattern_score(MULTI_PATTERN_T *multi_pattern, double score);

/**********************************************************************
  has_multi_pattern_score

  Does a multi_pattern object have a score?
**********************************************************************/
BOOLEAN_T has_multi_pattern_score(MULTI_PATTERN_T *multi_pattern);

/**********************************************************************
  clear_multi_pattern_score

  Sets the score member in a cisml multi_pattern object to null.
**********************************************************************/
void clear_multi_pattern_score(MULTI_PATTERN_T *multi_pattern);

/**********************************************************************
  get_multi_pattern_score

  Gets the score member from a cisml multi_pattern object.
**********************************************************************/
double get_multi_pattern_score(MULTI_PATTERN_T *multi_pattern);

/**********************************************************************
  set_multi_pattern_pvalue

  Sets the pvalue member in a cisml pattern object.
**********************************************************************/
void set_multi_pattern_pvalue(MULTI_PATTERN_T *multi_pattern, double pvalue);

/**********************************************************************
  has_multi_pattern_pvalue

  Does a multi_pattern object have a pvalue?
**********************************************************************/
BOOLEAN_T has_multi_pattern_pvalue(MULTI_PATTERN_T *multi_pattern);

/**********************************************************************
  clear_multi_pattern_pvalue

  Sets the pvalue member in a cisml multi_pattern object to null.
**********************************************************************/
void clear_multi_pattern_pvalue(MULTI_PATTERN_T *multi_pattern);

/**********************************************************************
  get_multi_pattern_pvalue

  Gets the pvalue member from a cisml multi_pattern object.
**********************************************************************/
double get_multi_pattern_pvalue(MULTI_PATTERN_T *multi_pattern);

/**********************************************************************
  get_multi_pattern_num_patterns

  Gets the number of patterns from a cisml multi_pattern object.
**********************************************************************/
int get_multi_pattern_num_patterns(MULTI_PATTERN_T *multi_pattern);

/**********************************************************************
  get_multi_pattern_patterns

  Gets the array of pointers to patterns from a cisml object.
**********************************************************************/
PATTERN_T **get_multi_pattern_patterns(MULTI_PATTERN_T *multi_pattern);

/**********************************************************************
  add_multi_pattern_pattern

  Adds a pattern to the array of pointers to patterns in a multi_pattern
  object.
**********************************************************************/
void add_multi_pattern_pattern(
  MULTI_PATTERN_T *multi_pattern,
  PATTERN_T* pattern
);

/**********************************************************************
  allocate_pattern

  Constructor for the cisml pattern data structure.
  Sets required fields to point to copies of the provided arguments.
  Other fields set to NULL.
**********************************************************************/
PATTERN_T *allocate_pattern(char *accession, char *name);

/**********************************************************************
  free_pattern

  Destructor for the cisml pattern data structure.
**********************************************************************/
void free_pattern(PATTERN_T *pattern);

/**********************************************************************
  add_pattern_elements_to_scanned_seq

  Updates each scanned sequence belonging to pattern with the 
  matched_elements assocaited with that sequence.

**********************************************************************/
void add_pattern_elements_to_scanned_seq(PATTERN_T *pattern);

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

**********************************************************************/
BOOLEAN_T add_pattern_matched_element(
  PATTERN_T *pattern,
  MATCHED_ELEMENT_T *element
);

/**********************************************************************
  set_pattern_max_stored_matches

  Sets the maximum number of matched elements that will be stored in a cisml
  pattern object. This requires creating a new element heap using the new
  maximum number of elements, copying the elements from the existing heap,
  freeing the existing heap, and pointing the pattern to the new heap.

  It fails if the new max is smaller than the current max_stored_matches.

  Returns TRUE if successful, FALSE otherwise.
**********************************************************************/
BOOLEAN_T set_pattern_max_stored_matches(PATTERN_T *pattern, int max);

/**********************************************************************
  set_pattern_max_pvalue_retained

  Sets the maximum p-value of the matched elements retained by
  the pattern.
**********************************************************************/
void set_pattern_max_pvalue_retained(PATTERN_T *pattern, double max_pvalue);

/**********************************************************************
  get_pattern_max_stored_matches

  Gets the maximum number of elements that will be stored in a cisml
  pattern object.
**********************************************************************/
int get_pattern_max_stored_matches(PATTERN_T *pattern); 

/**********************************************************************
  get_pattern_num_stored_matches

  Gets the maximum number of elements that will be stored in a cisml
  pattern object.
**********************************************************************/
int get_pattern_num_stored_matches(PATTERN_T *pattern); 

/**********************************************************************
  get_pattern_has_all_pvalues

  Returns the flag indicating whether or not all matched elements have 
  been retained, or only those with the smallest p-values.
**********************************************************************/
BOOLEAN_T get_pattern_has_all_pvalues(PATTERN_T *pattern);

/**********************************************************************
  set_pattern_has_all_pvalues

  Sets the flag indicating whether or not all matched elements have 
  been retained, or only those with the smallest p-values.
**********************************************************************/
void set_pattern_has_all_pvalues(PATTERN_T *pattern, BOOLEAN_T has_all_pvalues);

/**********************************************************************
  get_pattern_is_complete

  Returns a flag indicating whether or not all matched elements have 
  been added to the pattern. If flag is true, the element heap is no
  longer available and all matched elements are stored in an array
  of matched element pointers, sorted by p-value.
**********************************************************************/
BOOLEAN_T get_pattern_is_complete(PATTERN_T *pattern);

/**********************************************************************
  set_pattern_is_complete

  Sets the flag indicating that all matched elements have 
  been added to the pattern to true.  Moves  all matched elments out of 
  the element heap into an array of matched element pointers, sorted by 
  p-value. No further elements can be added to the pattern once
  this function has been called.
**********************************************************************/
void set_pattern_is_complete(PATTERN_T *pattern);

/**********************************************************************
  get_pattern_max_pvalue_retained

  Returns the maximum p-value of the matched elements retained by
  the pattern.
**********************************************************************/
double get_pattern_max_pvalue_retained(PATTERN_T *pattern);

/**********************************************************************
  get_pattern_name

  Gets the name member from a cisml pattern object.
**********************************************************************/
char *get_pattern_name(PATTERN_T* pattern);

/**********************************************************************
  get_pattern_accession

  Gets the accession member from a cisml pattern object.
**********************************************************************/
char *get_pattern_accession(PATTERN_T* pattern);

/**********************************************************************
  set_pattern_pvalue

  Sets the pvalue member in a cisml pattern object.
**********************************************************************/
void set_pattern_pvalue(PATTERN_T *pattern, double pvalue);

/**********************************************************************
  clean_pattern_pvalue

  Sets the pvalue member in a cisml pattern object to null.
**********************************************************************/
void clear_pattern_pvalue(PATTERN_T *pattern);

/**********************************************************************
  has_pattern_pvalue

  Does a pattern object have a pvalue?
**********************************************************************/
BOOLEAN_T has_pattern_pvalue(PATTERN_T *pattern);

/**********************************************************************
  get_pattern_pvalue

  Gets the pvalue member from a cisml pattern object.
**********************************************************************/
double get_pattern_pvalue(PATTERN_T* pattern);

/**********************************************************************
  has_pattern_qvalues

  Does the matched-elements for the pattern have qvalues available?
**********************************************************************/
BOOLEAN_T has_pattern_qvalues(PATTERN_T *pattern);

/**********************************************************************
  set_pattern_score

  Sets the score member in a cisml pattern object.
**********************************************************************/
void set_pattern_score(PATTERN_T *pattern, double score);

/**********************************************************************
  clear_pattern_score

  Sets the score member in a cisml pattern object to null.
**********************************************************************/
void clear_pattern_pvalue(PATTERN_T *pattern);

/**********************************************************************
  has_pattern_score

  Does a pattern object have a score?
**********************************************************************/
BOOLEAN_T has_pattern_score(PATTERN_T *pattern);

/**********************************************************************
  get_pattern_score

  Gets the score member from a cisml pattern object.
**********************************************************************/
double get_pattern_score(PATTERN_T* pattern);

/**********************************************************************
  set_pattern_db

  Sets the db member in a cisml pattern object.
**********************************************************************/
void set_pattern_db(PATTERN_T *pattern, char *db);

/**********************************************************************
  get_pattern_db

  Gets the db member from a cisml pattern object.
  May return NULL.
**********************************************************************/
char *get_pattern_db(PATTERN_T* pattern);

/**********************************************************************
  set_pattern_lsid

  Sets the lsid member in a cisml pattern object.
**********************************************************************/
void set_pattern_lsid(PATTERN_T *pattern, char *lsid);

/**********************************************************************
  get_pattern_lsid

  Gets the lsid member from a cisml pattern object.
  May return NULL.
**********************************************************************/
char *get_pattern_lsid(PATTERN_T* pattern);

/**********************************************************************
  get_pattern_num_scanned_sequences

  Gets the number of scanned_sequence objects in a cisml pattern object.
**********************************************************************/
int get_pattern_num_scanned_sequences(PATTERN_T *pattern);

/**********************************************************************
  get_pattern_num_scanned_positions

  Gets the number of sites scanned with a cisml pattern object.
**********************************************************************/
long get_pattern_num_scanned_positions(PATTERN_T *pattern);

/**********************************************************************
  get_pattern_num_matched_elements

  Gets the total number of matched_element objects contained in a cisml
  pattern object.
**********************************************************************/
int get_pattern_num_matched_elements(PATTERN_T *pattern);

/**********************************************************************
  get_pattern_scanned_sequences

  Gets the array of pointers to scanned_sequence objects in a cisml
  pattern object.  May return NULL.
**********************************************************************/
SCANNED_SEQUENCE_T **get_pattern_scanned_sequences(PATTERN_T *pattern);

/**********************************************************************
  allocate_scanned_sequence

  Constructor for the cisml scanned_sequence data structure.
  Sets required fields to point to copies of the provided arguments.
  Other fields set to NULL.
**********************************************************************/
SCANNED_SEQUENCE_T *allocate_scanned_sequence(
  char *accession,
  char *name,
  PATTERN_T *parent
);

/**********************************************************************
  free_scanned_sequence

  Destructor for the cisml scanned_sequence data structure.
**********************************************************************/
void free_scanned_sequence(SCANNED_SEQUENCE_T *scanned_sequence);

/**********************************************************************
  free_scanned_sequence_from_matched_elements

  Removes the matches from a scanned sequences only
**********************************************************************/
void free_scanned_sequence_from_matched_elements(
		SCANNED_SEQUENCE_T *scanned_sequence
);

/**********************************************************************
  get_scanned_sequence_parent

  Get a pointer to the pattern pointing to this scanned_sequence object.
**********************************************************************/
PATTERN_T *get_scanned_sequence_parent(
  SCANNED_SEQUENCE_T *scanned_sequence
);

/**********************************************************************
  get_scanned_sequence_name

  Gets the program_name member from a cisml scanned_sequence object.
**********************************************************************/
char *get_scanned_sequence_name(SCANNED_SEQUENCE_T *scanned_sequence);

/**********************************************************************
  get_scanned_sequence_accession

  Gets the accession member from a cisml scanned_sequence object.
**********************************************************************/
char *get_scanned_sequence_accession(SCANNED_SEQUENCE_T* scanned_sequence);

/**********************************************************************
  set_scanned_sequence_pvalue

  Sets the pvalue member in a cisml scanned_sequence object.
**********************************************************************/
void set_scanned_sequence_pvalue(
  SCANNED_SEQUENCE_T *scanned_sequence,
  double pvalue
);

/**********************************************************************
  clear_scanned_sequence_pvalue

  Sets the pvalue member in a cisml scanned_sequence object to null.
**********************************************************************/
void clear_scanned_sequence_pvalue(SCANNED_SEQUENCE_T *scanned_sequence);

/**********************************************************************
  has_scanned_sequence_pvalue

  Does a scanned_sequence object have a pvalue?
**********************************************************************/
BOOLEAN_T has_scanned_sequence_pvalue(SCANNED_SEQUENCE_T *scanned_sequence);

/**********************************************************************
  get_scanned_sequence_pvalue

  Gets the pvalue member from a cisml scanned_sequence object.
**********************************************************************/
double get_scanned_sequence_pvalue(SCANNED_SEQUENCE_T* scanned_sequence);

/**********************************************************************
  set_scanned_sequence_score

  Sets the score member in a cisml scanned_sequence object.
**********************************************************************/
void set_scanned_sequence_score(
  SCANNED_SEQUENCE_T *scanned_sequence,
  double score
);

/**********************************************************************
  clear_scanned_sequence_score

  Sets the score member in a cisml scanned_sequence object to null.
**********************************************************************/
void clear_scanned_sequence_score(SCANNED_SEQUENCE_T *scanned_sequence);

/**********************************************************************
  has_scanned_sequence_score

  Does a scanned_sequence object have a score?
**********************************************************************/
BOOLEAN_T has_scanned_sequence_score(SCANNED_SEQUENCE_T *scanned_sequence);

/**********************************************************************
  get_scanned_sequence_score

  Gets the score member from a cisml scanned_sequence object.
**********************************************************************/
double get_scanned_sequence_score(SCANNED_SEQUENCE_T* scanned_sequence);

/**********************************************************************
  set_scanned_sequence_length

  Sets the length member in a cisml scanned_sequence object.
**********************************************************************/
void set_scanned_sequence_length(
  SCANNED_SEQUENCE_T *scanned_sequence,
  int length
);

/**********************************************************************
  clear_scanned_sequence_length

  Sets the length member in a cisml scanned_sequence object to null.
**********************************************************************/
void clear_scanned_sequence_length(SCANNED_SEQUENCE_T *scanned_sequence);

/**********************************************************************
  has_scanned_sequence_length

  Does a scanned_sequence object have a length?
**********************************************************************/
BOOLEAN_T has_scanned_sequence_length(SCANNED_SEQUENCE_T *scanned_sequence);

/**********************************************************************
  get_scanned_sequence_length

  Gets the length member from a cisml scanned_sequence object.
**********************************************************************/
int get_scanned_sequence_length(SCANNED_SEQUENCE_T* scanned_sequence);

/**********************************************************************
  set_scanned_sequence_db

  Sets the db member in a cisml scanned_sequence object.
**********************************************************************/
void set_scanned_sequence_db(SCANNED_SEQUENCE_T *scanned_sequence, char *db);

/**********************************************************************
  get_scanned_sequence_db

  Gets the db member from a cisml scanned_sequence object.
  May return NULL.
**********************************************************************/
char *get_scanned_sequence_db(SCANNED_SEQUENCE_T* scanned_sequence);

/**********************************************************************
  set_scanned_sequence_lsid

  Sets the lsid member in a cisml scanned_sequence object.
**********************************************************************/
void set_scanned_sequence_lsid(
  SCANNED_SEQUENCE_T *scanned_sequence, char *lsid
);

/**********************************************************************
  get_scanned_sequence_lsid

  Gets the lsid member from a cisml scanned_sequence object.
  May return NULL.
**********************************************************************/
char *get_scanned_sequence_lsid(SCANNED_SEQUENCE_T* scanned_sequence);

/**********************************************************************
  add_scanned_sequence_matched_element

  Adds a matched element to the array of matched elements in
  the scanned sequence.
**********************************************************************/
void add_scanned_sequence_matched_element(
  SCANNED_SEQUENCE_T *sequence,
  MATCHED_ELEMENT_T *element
);

/**********************************************************************
  add_scanned_sequence_scanned_element

  Increments the count of scanned elements in a scanned_sequence.
**********************************************************************/
void add_scanned_sequence_scanned_element( SCANNED_SEQUENCE_T *sequence);

/**********************************************************************
  get_scanned_sequence_num_matched_elements

  Gets the number of matched_element objects in a cisml 
  scanned_sequence object.
**********************************************************************/
int get_scanned_sequence_num_matched_elements(SCANNED_SEQUENCE_T *sequence);

/**********************************************************************
  get_scanned_sequence_num_scanned_positions

  Gets the number of positions in the scanned_sequence where we
  have scanned for a matched_element.
**********************************************************************/
long get_scanned_sequence_num_scanned_positions(SCANNED_SEQUENCE_T *sequence);

/**********************************************************************
  get_scanned_sequence_num_matched_elements

  Gets the number of matched_element objects in a cisml
  scanned_sequence object.
**********************************************************************/
int get_scanned_sequence_num_matched_elements(SCANNED_SEQUENCE_T *sequence);

/**********************************************************************
  get_scanned_sequence_num_scanned_positions

  Gets the number of positions in the scanned_sequence where we
  have scanned for a matched_element.
**********************************************************************/
long get_scanned_sequence_num_scanned_positions(SCANNED_SEQUENCE_T *sequence);

/**********************************************************************
  get_scanned_sequences_elements

  Gets the array of pointers to matched_element objects in a cisml
  scanned_sequence object.  May return NULL.
**********************************************************************/
MATCHED_ELEMENT_T **get_scanned_sequence_matched_elements(
    SCANNED_SEQUENCE_T *sequence
);

/**********************************************************************
  get_scanned_sequences_matched_elements

  Gets the array of pointers to matched_element objects in a cisml 
  scanned_sequence object.  May return NULL.
**********************************************************************/
MATCHED_ELEMENT_T **get_scanned_sequence_matched_elements(
    SCANNED_SEQUENCE_T *sequence
);

/**********************************************************************
  allocate_matched_element

  Constructor for the cisml matched_element data structure.
  Sets required fields to the provided arguments.
  Other fields set to NULL.
**********************************************************************/
MATCHED_ELEMENT_T *allocate_matched_element(
  int start,
  int stop,
  SCANNED_SEQUENCE_T *parent
);

/**********************************************************************
  allocate_matched_element_without_inversion

  Alternative Constructor for the cisml matched_element data structure. 
  Sets required fields to the provided arguments. 
  Other fields set to NULL.

  JH: I had to add this because the existing constructor inverts the
      DNA sequence if the start site is greater than the stop site. I
      could not simply copy an existing matched element's content and
      feed it to the constructor without the content being changed.
      For the record I think this is a poor way to design a constructor.
**********************************************************************/
MATCHED_ELEMENT_T *allocate_matched_element_without_inversion(
  int start, 
  int stop, 
  const char *seq,
  SCANNED_SEQUENCE_T *parent
);


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
  SCANNED_SEQUENCE_T *parent
);

/**********************************************************************
  free_matched_element

  Destructor for the cisml matched_element data structure.
**********************************************************************/
void free_matched_element(MATCHED_ELEMENT_T *element);

/**********************************************************************
  get_matched_element_start

  Gets the start member from a cisml matched_element object.
**********************************************************************/
int get_matched_element_start(MATCHED_ELEMENT_T* matched_element);

/**********************************************************************
  set_matched_element_start

  Sets the start member from a cisml matched_element object.
**********************************************************************/
void set_matched_element_start(MATCHED_ELEMENT_T* matched_element, int newstart);



/**********************************************************************
  get_matched_element_stop

  Gets the stop member from a cisml matched_element object.
**********************************************************************/
int get_matched_element_stop(MATCHED_ELEMENT_T* matched_element);

/**********************************************************************
  set_matched_element_score

  Sets the score member in a cisml matched_element object.
**********************************************************************/
void set_matched_element_score(
  MATCHED_ELEMENT_T *element,
  double score
);

/**********************************************************************
  has_matched_element_score

  Does a matched_element object have a score?
**********************************************************************/
BOOLEAN_T has_matched_element_score(MATCHED_ELEMENT_T *element);

/**********************************************************************
  get_matched_element_score

  Gets the score member from a cisml matched_element object.
**********************************************************************/
double get_matched_element_score(MATCHED_ELEMENT_T* element);

/**********************************************************************
  set_matched_element_pvalue

  Sets the pvalue member in a cisml matched_element object.
**********************************************************************/
void set_matched_element_pvalue(
  MATCHED_ELEMENT_T *element,
  double pvalue
);

/**********************************************************************
  has_matched_element_pvalue

  Does a matched_element object have a pvalue?
**********************************************************************/
BOOLEAN_T has_matched_element_pvalue(MATCHED_ELEMENT_T *element);

/**********************************************************************
  get_matched_element_pvalue

  Gets the pvalue member from a cisml matched_element object.
**********************************************************************/
double get_matched_element_pvalue(MATCHED_ELEMENT_T* element);

/**********************************************************************
  set_matched_element_clusterid

  Sets the clusterid member in a cisml matched_element object.
**********************************************************************/
void set_matched_element_clusterid(
  MATCHED_ELEMENT_T *element,
  char *clusterid
);

/**********************************************************************
  get_matched_element_clusterid

  Gets the clusterid member from a cisml matched_element object.
  May return NULL.
**********************************************************************/
char *get_matched_element_clusterid(MATCHED_ELEMENT_T* element);

/**********************************************************************
  get_matched_element_sequence

  Gets the sequence member from a cisml matched_element object.
  May return NULL.
**********************************************************************/
const char *get_matched_element_sequence(MATCHED_ELEMENT_T* element);

/**********************************************************************
  set_matched_element_strand

  Sets the strand member in a cisml matched_element object.
**********************************************************************/
void set_matched_element_strand(MATCHED_ELEMENT_T* element, char strand);

/**********************************************************************
  get_matched_element_strand

  Gets the strand member in a cisml matched_element object.
**********************************************************************/
char get_matched_element_strand(MATCHED_ELEMENT_T* element, char strand);

/*************************************************************************
 * sort_matched_elements
 *
 * Sort an array of pointers to matched-elements sites by pvalue 
 * or by sequence name and position.
 *
 * JH - I opened this up to make it accessible outside of the CISML code
 * because I needed the matched elements sorted for the BLS implementation
 *************************************************************************/

void sort_matched_elements(
  BOOLEAN_T sort_by_pvalue,
  int num_elements,
  MATCHED_ELEMENT_T **elements
);

/**********************************************************************
  set_matched_element_sequence

  Sets the sequence member in a cisml matched_element object.
**********************************************************************/
void set_matched_element_sequence(MATCHED_ELEMENT_T* element, char *seq);

/**********************************************************************
  print_cisml_start

  Print the opening section of the CisML XML
**********************************************************************/
void print_cisml_start(
  FILE* out, 
  CISML_T *cisml, 
  BOOLEAN_T print_header,
  const char *stylesheet, 
  BOOLEAN_T print_namespace
);

/**********************************************************************
  print_cisml_end

  Print the closing section of the CisML XML
**********************************************************************/
void print_cisml_end(FILE* out);

/**********************************************************************
  print_cisml_parmeters

  Print parameters element for CisML
**********************************************************************/
void print_cisml_parameters(FILE *out, CISML_T *cisml);

/**********************************************************************
  print_cisml

  Print the cisml data structure as CisML XML
**********************************************************************/
void print_cisml(
		FILE* out,                /* The target file IN */
		CISML_T *cisml,           /* The cisml data structure to save IN */
		BOOLEAN_T print_header,   /* FLAG indicates if the header should be printed IN */
		const char *stylesheet,   /* the stylesheet file IN */
		BOOLEAN_T print_namespace /* FLAG indicates if the namespace should be printed IN */
);

/*************************************************************************
 * Calculate the q-values for each matched-element in this pattern
 * from the p-value.
 *************************************************************************/
void pattern_calculate_qvalues(
  PATTERN_T *pattern,
  ARRAY_T *sampled_pvalues
);

/**********************************************************************
  print_cisml_start_pattern

  Print the starting tag for a CisML pattern 
**********************************************************************/
void print_cisml_start_pattern(
  CISML_T *cisml,
  FILE *out,
  PATTERN_T *pattern
);

/**********************************************************************
  print_cisml_end_pattern

  Print the ending tag for a CisML pattern 
**********************************************************************/
void print_cisml_end_pattern(FILE *out);

/**********************************************************************
  print_cisml_scanned_sequences()

  Print the XML for of a CisML scanned sequence element.
**********************************************************************/
void print_cisml_scanned_sequences(
  CISML_T *cisml,
  FILE *out,
  int num_seqs,
  SCANNED_SEQUENCE_T **sequences
);

/**********************************************************************
  print_cisml_scanned_sequence_start()

  Print the starting XML tag for a CisML scanned sequence element.
**********************************************************************/
void print_cisml_scanned_sequence_start(
  CISML_T *cisml,
  FILE *out,
  SCANNED_SEQUENCE_T *seq
);

/**********************************************************************
  print_cisml_scanned_sequence_end()

  Print the ending XML tag for a CisML scanned sequence element.
**********************************************************************/
void print_cisml_scanned_sequence_end(FILE *out);

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
);

/**********************************************************************
 * This function saves the CisML results as a set of files in a
 * directory. The files include:
 *   fimo.xml which uses the CisML schema
 *   fimo.txt which uses plain text with tab-delimited columns
 *   fimo.gff which uses the gff file format.
 *********************************************************************/
void print_cisml_as_text(CISML_T *cisml);

/**********************************************************************
  read_cisml

  Reads in a CisML XML file and create the cisml_t data structure
**********************************************************************/
CISML_T* read_cisml(
	char* 		cisml_filename // name of the cisml file
);

/**********************************************************************
  read_multiple_patterns_from_cisml

  reads a multi-pattern-scan from a CisML xml file
**********************************************************************/
void read_multiple_patterns_from_cisml(
		xmlXPathContextPtr 	xpath_ctxt,    	// cisml XPath context.
		int* 				num_patterns,	// cisml number of pattern in the file
		CISML_T* 			cisml			// the retrieved cisml struct
);

/**********************************************************************
  read_patterns_from_cisml

  reads the patterns from a CisML xml file
**********************************************************************/
void read_patterns_from_cisml(
		xmlXPathContextPtr 	xpath_ctxt,    	// cisml XPath context.
		int* 				num_patterns,	// cisml number of pattern in the file
		CISML_T* 			cisml			// the retrieved cisml struct
);

/**********************************************************************
  get_pattern_from_xmlnode

  returns a pattern from a specified xml Node
**********************************************************************/
PATTERN_T* get_pattern_from_xmlnode(xmlNodePtr currPatternNode);

/**********************************************************************
  read_scanned_sequences_from_cisml

  Reads in the sequences for a pattern
**********************************************************************/
void read_scanned_sequences_from_cisml(
  CISML_T* 				cisml,			// the retrieved cisml struct
  PATTERN_T*			pattern, 		// the pattern of interest
  xmlNode*				first_sequence	// first child of the pattern node
);

/**********************************************************************
  read_matched_elements_from_cisml

  Reads in sequences matches from a cisml file
**********************************************************************/
void  read_matched_elements_from_cisml(
  CISML_T* 				cisml,			// the retrieved cisml struct
  SCANNED_SEQUENCE_T*	scanned_seq,	// the pattern of interest
  xmlNode*				firstmatch		// first child of the pattern node
);

#endif
