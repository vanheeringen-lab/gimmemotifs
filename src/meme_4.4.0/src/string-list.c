/**************************************************************************
 * FILE: string-list.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 12-22-98
 * PROJECT: shared
 * COPYRIGHT: 1999-2008, WSN 
 * VERSION: $Revision: 1.5 $
 * DESCRIPTION: Data structure for manipulating a list of strings. 
 **************************************************************************/
#include "string-list.h"
#include "utils.h"
#include <string.h>
#include <assert.h>

/*************************************************************************
 * Primary data structure.
 *************************************************************************/
struct string_list_t {
  int     num_strings;    /* Number of strings in the array. */
  int     max_strings;    /* Total amount of memory allocated. */
  int     longest_string; /* Length of the longest allowed string. */
  char** strings;        /* The strings themselves. */
  double* scores;        /* Scores for optional sorting of strings */ 
};

/*************************************************************************
 * Allocate dynamic memory for a string list.
 *************************************************************************/
#define DEFAULT_MAX_STRINGS 10
#define DEFAULT_STRING_LENGTH  1000
STRING_LIST_T* new_string_list
  ()
{
  STRING_LIST_T* new_list;  /* The list being created. */
  int             i_string;  /* Index into the list. */

  new_list = (STRING_LIST_T*)mm_calloc(1, sizeof(STRING_LIST_T));
  new_list->num_strings = 0;
  new_list->max_strings = DEFAULT_MAX_STRINGS;
  new_list->longest_string = DEFAULT_STRING_LENGTH;
  new_list->strings = (char**)mm_calloc(DEFAULT_MAX_STRINGS, sizeof(char*));
  for (i_string = 0; i_string < DEFAULT_MAX_STRINGS; i_string++) {
    new_list->strings[i_string] = (char*)mm_calloc(DEFAULT_STRING_LENGTH + 1,
						   sizeof(char));
  }
  new_list->scores = (double*)mm_calloc(DEFAULT_MAX_STRINGS, sizeof(double));
  return(new_list);
}

/*************************************************************************
 * Is the given list a null list?
 *************************************************************************/
static void check_null_list
  (STRING_LIST_T*  a_list)
{
  if (a_list == NULL) {
    die("Attempted to access null string list.\n");
  }
}
/*************************************************************************
 * Get the length of the longest string in a list.
 *************************************************************************/
int max_string_length
  (STRING_LIST_T*  a_list)
{
  check_null_list(a_list);

  return(a_list->longest_string);
}

/*************************************************************************
 * Get the total number of strings in a list.
 *************************************************************************/
int get_num_strings
  (STRING_LIST_T*  a_list)
{
  check_null_list(a_list);

  return(a_list->num_strings);
}

/*************************************************************************
 * Make sure each string in the list is long enough.
 *************************************************************************/
static void resize_string_list
  (int            new_length,
   STRING_LIST_T* a_list)
{
  int i_string;

  if (new_length > a_list->longest_string) {
    a_list->longest_string = new_length + 1;

    for (i_string = 0; i_string < a_list->max_strings; i_string++) {
      a_list->strings[i_string] 
	= (char*)mm_realloc(a_list->strings[i_string], 
			    a_list->longest_string * sizeof(char));
    }
  }
}					    

/*************************************************************************
 * Get the nth string from a list.
 *************************************************************************/
char* get_nth_string
  ( int             n,
   STRING_LIST_T*  a_list)
{
  check_null_list(a_list);

  if (n > a_list->max_strings) {
    die("Attempted to access string beyond end of list.\n");
  } else if (n > a_list->num_strings) {
    die("Attempted to access uninitialized string.\n");
  }

  assert(a_list->strings[n] != NULL);
  return(a_list->strings[n]);
}

void set_nth_string
  (char*          new_string,            
   int            n,
   STRING_LIST_T* a_list)
{
  /* Check bounds. */
  get_nth_string(n, a_list);

  /* Make all the strings longer if this one is too long. */
  resize_string_list((int)strlen(new_string), a_list);

  /* Copy the string into the list. */
  strcpy(a_list->strings[n], new_string);
}

/*************************************************************************
 * Get nth score in a scored string list.
 *************************************************************************/
double get_nth_score
  (int            n,
   STRING_LIST_T* a_list)
{
  check_null_list(a_list);

  if (n > a_list->max_strings) {
    die("Attempted to access string beyond end of list.\n");
  } else if (n > a_list->num_strings) {
    die("Attempted to access uninitialized string.\n");
  }

  return(a_list->scores[n]);
}

void set_nth_score
  (double         new_score,
   int            n,
   STRING_LIST_T* a_list)
{
  // Check bounds.
  get_nth_string(n, a_list);

  // Copy the score into the list.
  a_list->scores[n] = new_score;
}

/*************************************************************************
 *  Extract an array of scores.  Memory must be freed by the caller.
 *************************************************************************/
ARRAY_T* get_string_list_scores
  (STRING_LIST_T* a_list)
{
  ARRAY_T* return_value = allocate_array(get_num_strings(a_list));

  int num_strings = get_num_strings(a_list);
  int i_string;
  for (i_string = 0; i_string < num_strings; i_string++) {
    set_array_item(i_string, a_list->scores[i_string], return_value);
  }
  return(return_value);
}

/*************************************************************************
 * Append to the nth string.
 *************************************************************************/
void append_to_nth_string
  (char*          new_string,            
   int            n,
   STRING_LIST_T* a_list)
{
  /* Check bounds. */
  get_nth_string(n, a_list);

  /* Make all the strings longer if this one is too long. */
  resize_string_list(a_list->longest_string + strlen(new_string), a_list);

  /* Copy the string into the list. */
  strcat(a_list->strings[n], new_string);
}


/*************************************************************************
 * Add a string to the end of a given list.
 *************************************************************************/
void add_string
  (char*     a_string,
   STRING_LIST_T* a_list)
{
  int i_string;

  check_null_list(a_list);

  /* Make sure we're not adding a null string. */
  if (a_string == NULL) {
    die("Adding null string to string list.");
  }

  /* Reallocate space if there isn't any. */
  if (a_list->num_strings >= a_list->max_strings) {
    a_list->strings = (char**)mm_realloc(a_list->strings, 
					 (a_list->max_strings 
					  + DEFAULT_MAX_STRINGS)
					 * sizeof(char*));
    for (i_string = 0; i_string < DEFAULT_MAX_STRINGS; i_string++) {
	a_list->strings[a_list->max_strings + i_string] 
	    = (char*)mm_calloc(a_list->longest_string + 1, sizeof(char));
    }
    a_list->max_strings += DEFAULT_MAX_STRINGS;



    a_list->scores
      = (double*)mm_realloc(a_list->scores,
			   (a_list->max_strings 
			    + DEFAULT_MAX_STRINGS) * sizeof(double));


  }

  /* Make all the strings longer if this one is too long. */
  resize_string_list((int)strlen(a_string), a_list);

  /* Put the string in the list. */
  strcpy(a_list->strings[a_list->num_strings], a_string);
  a_list->scores[a_list->num_strings] = 0;

  (a_list->num_strings)++;
}

/*************************************************************************
 * Add a string with score to the end of a given list.
 *************************************************************************/
void add_string_with_score
  (char*     a_string,
   STRING_LIST_T* a_list,
   double score)
{
  int index = a_list->num_strings;
  add_string(a_string, a_list);
  a_list->scores[index] = score;
}

/*************************************************************************
 * Add a list of strings to the end of a given list.
 *************************************************************************/
void add_strings
  (STRING_LIST_T*  source_list,
   STRING_LIST_T*       target_list)
{
  int i_string;

  check_null_list(source_list);
  check_null_list(target_list);

  for (i_string = 0; i_string < get_num_strings(source_list); i_string++) {
    add_string(get_nth_string(i_string, source_list), target_list);
  }
}

/*************************************************************************
 * Remove a string from a given list.
 *
 * It is an error if the string does not appear in the list.
 *************************************************************************/
void remove_string
  (char*     a_string,
   STRING_LIST_T* a_list)
{
  BOOLEAN_T found_it;    /* Have we found the desired string? */
  char*    this_string; /* The current string. */
  int       i_string;    /* Index of the current string.*/

  check_null_list(a_list);

  /* Make sure the string exists in the list. */
  if (!have_string(a_string, a_list)) {
    die("Attempted to remove string %s from list that doesn't contain it.\n",
	a_string);
  }

  found_it = FALSE;
  for (i_string = 0; i_string < get_num_strings(a_list) - 1; i_string++) {
    this_string = get_nth_string(i_string, a_list);

    /* Is this the string we are looking for? */
    if (strcmp(this_string, a_string) == 0) {
      found_it = TRUE;
    }

    /* Copy the next string into the current location. */
    if (found_it) {
      strcpy(this_string, get_nth_string(i_string + 1, a_list));
    }
  }
  (a_list->num_strings)--;
}

/*************************************************************************
 * Remove a set of strings from a given list.
 *
 * It is an error for the source list to contain any strings that
 * don't appear in the target list.
 *************************************************************************/
void remove_strings
  (STRING_LIST_T*  source_list,
   STRING_LIST_T*       target_list)
{
  int i_string;

  for (i_string = 0; i_string < get_num_strings(source_list); i_string++) {
    remove_string(get_nth_string(i_string, source_list), target_list);
  }
}

/*************************************************************************
 * Does a given string appear in a given list?
 *************************************************************************/
BOOLEAN_T have_string
  (char*           a_string,
   STRING_LIST_T*  a_list)
{
  char* this_string; /* The current string in the list. */
  int    i_string;    /* Index of the current string. */
  
  check_null_list(a_list);

  for (i_string = 0; i_string < get_num_strings(a_list); i_string++) {
    this_string = get_nth_string(i_string, a_list);
    if (strcmp(this_string, a_string) == 0) {
      return(TRUE);
    }
  }
  return(FALSE);
}

/*************************************************************************
 * Return a index of a string in the given list.
 * Returns -1 if string is not in list.
 *************************************************************************/
int get_index_in_string_list(char* a_string,
                             STRING_LIST_T* a_list) {
  char* this_string; /* The current string in the list. */
  int   i_string;    /* Index of the current string. */
  
  check_null_list(a_list);

  for (i_string = 0; i_string < get_num_strings(a_list); i_string++) {
    this_string = get_nth_string(i_string, a_list);
    if (strcmp(this_string, a_string) == 0) {
      return(i_string);
    }
  }
  // String was not found in list.
  return(-1);
}

/**************************************************************************
 * Make a copy of a string list.
 **************************************************************************/
STRING_LIST_T* copy_string_list
  (STRING_LIST_T* string_list)
{
  int i_string;
  char* this_string;
  STRING_LIST_T* new_list;

  /* Create the new list. */
  new_list = new_string_list();

  /* Add stuff to it. */
  for (i_string = 0; i_string < string_list->num_strings; i_string++) {
    this_string = get_nth_string(i_string, string_list);
    add_string(this_string, new_list);
  }

  /* Return it. */
  return(new_list);
}

/*************************************************************************
 * A version of 'strcmp' suitable for passing to 'qsort'.
 *************************************************************************/
static int string_compare
  (const void * elem1,
   const void * elem2)
{
  const char** string_ptr1 = (const char**)elem1;
  const char** string_ptr2 = (const char**)elem2;

  return(strcmp(*string_ptr1, *string_ptr2));
}

/*************************************************************************
 * Compare the scores for 'qsort'.
 *************************************************************************/

/* I define this as a struct only so that the two items can be passed as
   a single pointer to the 'qsort' routine. */
typedef struct forqsort {
  double score;
  char*  string;
} FORQSORT_T;

static int score_compare
  (const void * elem1,
   const void * elem2)
{
  const double key1 = ((FORQSORT_T *)elem1)->score;
  const double key2 = ((FORQSORT_T *)elem2)->score;
  const char* string1 = ((FORQSORT_T *)elem1)->string;
  const char* string2 = ((FORQSORT_T *)elem2)->string;

  if (key1 < key2) {
    return(-1);
  } else if (key1 > key2) {
    return(1);
  }
  return(strcmp(string1, string2));
}

static int score_compare_reverse
  (const void * elem1,
   const void * elem2)
{
  return score_compare(elem2, elem1);
}

/*************************************************************************
 * Sort a list of strings in place.
 *************************************************************************/
void sort_string_list 
  (STRING_LIST_T* a_list)
{
  check_null_list(a_list);

  qsort((void *)(a_list->strings), a_list->num_strings, 
	sizeof(char*), string_compare);
}

/*************************************************************************
 * Sort a list of strings by score.
 *************************************************************************/
void sort_string_list_by_score
  (STRING_LIST_T* a_list,
   BOOLEAN_T      reverse) //Descending order
{
  check_null_list(a_list);
  
  FORQSORT_T*  forqsort; // Array of strings to be sorted.
  int num_stored = a_list->num_strings; // Number of strings to be sorted.
  int i_stored = 0;

  forqsort = (FORQSORT_T *)mm_malloc(sizeof(FORQSORT_T) * 
				    num_stored);

  for (i_stored = 0; i_stored < num_stored; i_stored++) {

    forqsort[i_stored].string =
      (char *)mm_malloc(sizeof(char) * 
		       (strlen(get_nth_string(i_stored, 
					      a_list)) + 1));
    
    sprintf(forqsort[i_stored].string, "%s", get_nth_string(i_stored, 
							    a_list));

    forqsort[i_stored].score = a_list->scores[i_stored];
  }
  if (reverse) {
    qsort((void *)(forqsort), i_stored, 
	  sizeof(FORQSORT_T), score_compare_reverse);
  }
  else {
    qsort((void *)(forqsort), i_stored, 
	  sizeof(FORQSORT_T), score_compare);
  }
  // Transfer values in the original array
  for (i_stored = 0; i_stored < num_stored; i_stored++) {
    set_nth_string(forqsort[i_stored].string, i_stored, a_list);
    a_list->scores[i_stored] = forqsort[i_stored].score;
  }

  // Free up memory for the sorted sequences.
  for (i_stored = 0; i_stored < num_stored; i_stored++) {
    myfree(forqsort[i_stored].string);
  }
  myfree(forqsort);  

}
/***************************************************************************
 * Given two sets, A and B, find 
 *  - A intersect B,
 *  - A minus B, and
 *  - B minus A.
 *
 * Allocates memory for the three new string lists, which must be
 * freed by the caller.
 *
 * If any of the return parameters is supplied as NULL, then that
 * parameter will not be filled in
 ***************************************************************************/
void overlap_string_lists
  (STRING_LIST_T*  list_a,
   STRING_LIST_T*  list_b,
   STRING_LIST_T**      intersection,
   STRING_LIST_T**      a_minus_b,
   STRING_LIST_T**      b_minus_a)
{
  int             i_all;         /* Index of current string. */
  char*          this_string;   /* The current string. */

  check_null_list(list_a);
  check_null_list(list_b);

  /* Allocate memory for the new list. */
  if (intersection != NULL) {
    *intersection = new_string_list();
  }
  if (a_minus_b != NULL) {
    *a_minus_b = new_string_list();
  }
  if (b_minus_a != NULL) {
    *b_minus_a = new_string_list();
  }

  /* Enumerate the first list. */
  for (i_all = 0; i_all < get_num_strings(list_a); i_all++) {
    this_string = get_nth_string(i_all, list_a);

    if (have_string(this_string, list_b)) {
      if (intersection != NULL) {
	add_string(this_string, *intersection);
      }
    } else {
      if (a_minus_b != NULL) {
	add_string(this_string, *a_minus_b);
      }
    }
  }

  /* Enumerate the second list. */
  for (i_all = 0; i_all < get_num_strings(list_b); i_all++) {
    this_string = get_nth_string(i_all, list_b);

    if (!have_string(this_string, list_a)) {
      if (b_minus_a != NULL) {
	add_string(this_string, *b_minus_a);
      }
    }
  }
}

/*************************************************************************
 * Combine all the strings in a given list into a single string
 * separated by a given separator.
 *
 * Memory allocated for return value must be freed by caller.
 *************************************************************************/
char* combine_string_list
  (STRING_LIST_T*  a_list,
   char*          separator)
{
  int    string_length;   /* Length of the combined string. */
  char*  combined_string; /* String to be returned. */
  int    i_string;        /* Index of the current string. */

  check_null_list(a_list);

  /* Compute the length of the return string. */
  string_length = (get_num_strings(a_list) * 
		   (max_string_length(a_list) + strlen(separator))) + 2;

  /* Allocate memory for the string. */
  combined_string = (char*)mm_calloc(string_length, sizeof(char));

  /* Accumulate all the strings, with spaces between. */
  for (i_string = 0; i_string < get_num_strings(a_list); i_string++) {
    strcat(combined_string, get_nth_string(i_string, a_list));
    strcat(combined_string, separator);
  }

  /* Remove the final space. */
  combined_string[strlen(combined_string) - strlen(separator)] = '\0';

  /* Return the string. */
  return(combined_string);
}


/*****************************************************************************
 * Read a list of strings from a given file.
 *****************************************************************************/
STRING_LIST_T* read_string_list
  (FILE* infile)
{
  char           this_line[DEFAULT_STRING_LENGTH];
  char*          fgets_result;
  STRING_LIST_T* return_value;

  /* Allocate dynamic memory for a maximal list of names. */
  return_value = new_string_list();

  /* Read the first name. */
  fgets_result = fgets(this_line, DEFAULT_STRING_LENGTH, infile);
  this_line[strlen(this_line)-1] = 0; // chop

  /* Make sure we got at least one name. */
  if (fgets_result == NULL) {
      die("Couldn't read a single name from the given file.");
  }

  while (fgets_result != NULL) {

     /* Store the name in the name list. */
    add_string(this_line, return_value);

     /* Read the next name. */
     fgets_result = fgets(this_line, DEFAULT_STRING_LENGTH, infile);
     this_line[strlen(this_line)-1] = 0; // chop
  }

  /* Die if we didn't read anything. */
  if (get_num_strings(return_value) == 0) {
    die("Failed to read any names.");
  }

  /* Tell the user what's up. */
  if (verbosity > NORMAL_VERBOSE) {
    fprintf(stderr, "Read %d names.\n", get_num_strings(return_value));
  }
  
  return(return_value);
}

/*****************************************************************************
 * Read a list of strings from the named file.
 *****************************************************************************/
STRING_LIST_T* read_string_list_from_file
  (char* filename) 
{
    STRING_LIST_T* string_list;
    FILE* string_list_file = NULL;

    if (open_file(
          filename, 
          "r", 
          1, 
          "string list", 
          "string list", 
          &string_list_file
        ) == 0) {
      die("Couldn't open the file %s.\n", filename);
    }
    string_list = read_string_list(string_list_file);
    (void) fclose(string_list_file);
    if (verbosity > HIGH_VERBOSE) {
      fprintf(stderr, "Read string list: ");
      write_string_list(" ", string_list, stderr);
    }

    return string_list;
}

/*************************************************************************
 * Write out a list of strings.
 *************************************************************************/
void write_string_list
  (char*           separator, /* String to separate strings in list. */
   STRING_LIST_T*  a_list,    /* The list to be printed. */
   FILE*           outfile)
{
  int i_string;

  check_null_list(a_list);

  /* Don't do anything if we go no strings. */
  if (get_num_strings(a_list) == 0) {
      return;
  }

  for (i_string = 0; i_string < get_num_strings(a_list) - 1; i_string++) {
    fprintf(outfile, "%s%s", get_nth_string(i_string, a_list), separator);
  }
  fprintf(outfile, "%s\n", get_nth_string(i_string, a_list));
}


/*************************************************************************
 * Free dynamic memory used by a list of strings.
 *************************************************************************/
void free_string_list
  (STRING_LIST_T* a_list)
{
  int i_string;

  if (a_list == NULL) {
    return;
  }

  for (i_string = 0; i_string < a_list->max_strings; i_string++) {
    myfree(a_list->strings[i_string]);
  }
  myfree(a_list->strings);
  myfree(a_list->scores);
  myfree(a_list);
}

/*************************************************************************
 * Test two string lists for equality.  Order matters.
 *************************************************************************/
BOOLEAN_T equal_string_lists
 (STRING_LIST_T* a_list,
  STRING_LIST_T* b_list)
{
  int i_string;

  check_null_list(a_list);
  check_null_list(b_list);

  // If there are different numbers of strings, then they cannot be equal.
  if (get_num_strings(a_list) != get_num_strings(b_list)) {
    return(FALSE);
  }

  for (i_string = 0; i_string < get_num_strings(a_list); i_string++) {
    if (strcmp(get_nth_string(i_string, a_list),
	       get_nth_string(i_string, b_list))) {
      return(FALSE);
    }
  }
  return(TRUE);
}

/*************************************************************************
 * Test whether a given string list contains duplicates.
 *
 * Optionally, print the duplicated strings to stderr.
 *************************************************************************/
BOOLEAN_T has_duplicates
  (char*          message, // If message == "", then don't print the list to stderr.
   STRING_LIST_T* my_list)
{
  BOOLEAN_T return_value = FALSE;
  BOOLEAN_T printed_message = FALSE;

  check_null_list(my_list);

  // Consider each string in the list.
  int i_string;
  int num_strings = get_num_strings(my_list);
  for (i_string = 0; i_string < num_strings; i_string++) {
    char* first_string = get_nth_string(i_string, my_list);

    // If a string appears more than twice, just print it the first time.
    BOOLEAN_T already_printed = FALSE;

    // Look for another occurrence of this string.
    int j_string;
    for (j_string = 0; j_string < num_strings; j_string++) {
      char* second_string = get_nth_string(j_string, my_list);

      // Are the strings the same, but with different indices?
      if ((i_string != j_string) && 
	  (strcmp(first_string, second_string) == 0)) {
	return_value = TRUE;

	// If requested, print the string.
	if ((strcmp(message, "") != 0)
	    && (already_printed == FALSE)) {

	  // Print the initial message, if one is given.
	  if (printed_message == FALSE) {
	    fprintf(stderr, "%s", message);
	    printed_message = TRUE;
	  }

	  fprintf(stderr, " %s", first_string);
	  already_printed = TRUE;
	}
      }
    }
  }

  // Print the final EOL.  
  if (printed_message == TRUE) {
    fprintf(stderr, "\n");
  }

  return(return_value);
}


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */











