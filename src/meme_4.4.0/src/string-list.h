/**************************************************************************
 * FILE: string-list.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 12-22-98
 * PROJECT: shared
 * COPYRIGHT: 1999-2008, WSN
 * VERSION: $Revision: 1.4 $
 * DESCRIPTION: Data structure for manipulating a list of strings. 
 **************************************************************************/
#ifndef STRING_LIST_H
#define STRING_LIST_H

#include "utils.h"
#include "array.h"

/*************************************************************************
 * Primary data structure.
 *************************************************************************/
typedef struct string_list_t STRING_LIST_T;

/*************************************************************************
 * Allocate dynamic memory for a string list.
 *************************************************************************/
STRING_LIST_T* new_string_list
  ();

/*************************************************************************
 * Get the length of the longest string in a list.
 *************************************************************************/
int max_string_length
  (STRING_LIST_T* a_list);

/*************************************************************************
 * Get the total number of strings in a list.
 *************************************************************************/
int get_num_strings
  (STRING_LIST_T* a_list);

/*************************************************************************
 * Get or set the nth string in a list.
 *************************************************************************/
char* get_nth_string
  (int            n,
   STRING_LIST_T* a_list);

void set_nth_string
  (char*          new_string,            
   int            n,
   STRING_LIST_T* a_list);

/*************************************************************************
 * Get or set the nth score in a scored string list.
 *************************************************************************/
double get_nth_score
  (int            n,
   STRING_LIST_T* a_list);

void set_nth_score
  (double         new_score,
   int            n,
   STRING_LIST_T* a_list);

/*************************************************************************
 *  Extract an array of scores.  Memory must be freed by the caller.
 *************************************************************************/
ARRAY_T* get_string_list_scores
  (STRING_LIST_T* a_list);

/*************************************************************************
 *  Append a string onto the nth string in the list
 *************************************************************************/
void append_to_nth_string
  (char*          new_string,            
   int            n,
   STRING_LIST_T* a_list);

/*************************************************************************
 * Does the given string list have any space left in it?
 *************************************************************************/
BOOLEAN_T have_space
  (STRING_LIST_T* a_list);

/*************************************************************************
 * Add a string to the end of a given list.
 *************************************************************************/
void add_string
  (char*          a_string,
   STRING_LIST_T* a_list);

/*************************************************************************
 * Add a string with score to the end of a given list.
 *************************************************************************/
void add_string_with_score
  (char*          a_string,
   STRING_LIST_T* a_list,
   double score);

/*************************************************************************
 * Add a list of strings to the end of a given list.
 *************************************************************************/
void add_strings
  (STRING_LIST_T* source_list,
   STRING_LIST_T* target_list);

/*************************************************************************
 * Remove a string from a given list.
 *
 * It is an error if the string does not appear in the list.
 *************************************************************************/
void remove_string
  (char*    a_string,
   STRING_LIST_T* a_list);

/*************************************************************************
 * Remove a set of strings from a given list.
 *
 * It is an error for the source list to contain any strings that
 * don't appear in the target list.
 *************************************************************************/
void remove_strings
  (STRING_LIST_T* source_list,
   STRING_LIST_T* target_list);

/*************************************************************************
 * Does a given string appear in a given list?
 *************************************************************************/
BOOLEAN_T have_string
  (char*          a_string,
   STRING_LIST_T* a_list);

/*************************************************************************
 * Return index of a string in the given list.
 * Returns -1 if string is not in list.
 *************************************************************************/
int get_index_in_string_list(char* a_string,
                             STRING_LIST_T* a_list);

/*************************************************************************
 * Return a copy of the given list.
 *************************************************************************/
STRING_LIST_T* copy_string_list
  (STRING_LIST_T* string_list);

/*************************************************************************
 * Sort a list of strings in place.
 *************************************************************************/
void sort_string_list 
  (STRING_LIST_T* a_list);

/*************************************************************************
 * Sort a list of strings by score.
 *************************************************************************/
void sort_string_list_by_score
  (STRING_LIST_T* a_list,
   BOOLEAN_T      reverse); //Descending order

/***************************************************************************
 * Given two sets, A and B, find 
 *  - A union B,
 *  - A minus B, and
 *  - B minus A.
 ***************************************************************************/
void overlap_string_lists
  (STRING_LIST_T* list_a,
   STRING_LIST_T* list_b,
   STRING_LIST_T**      intersection,
   STRING_LIST_T**      a_minus_b,
   STRING_LIST_T**      b_minus_a);

/*************************************************************************
 * Combine all the strings in a given list into a single string
 * separated by a given separator.
 *************************************************************************/
char* combine_string_list
  (STRING_LIST_T* a_list,
   char*          separator);

/*****************************************************************************
 * Read a list of strings from a given file.
 *****************************************************************************/
STRING_LIST_T* read_string_list
  (FILE* infile);

/*****************************************************************************
 * Read a list of strings from the named file.
 *****************************************************************************/
STRING_LIST_T* read_string_list_from_file
  (char* filename);

/*************************************************************************
 * Write out a list of strings.
 *************************************************************************/
void write_string_list
  (char*          separator, /* String to separate strings in list. */
   STRING_LIST_T* a_list,    /* The list to be printed. */
   FILE*          outfile);

/*************************************************************************
 * Free dynamic memory used by a list of strings.
 *************************************************************************/
void free_string_list
  (STRING_LIST_T* a_list);

/*************************************************************************
 * Test two string lists for equality.  Order matters.
 *************************************************************************/
BOOLEAN_T equal_string_lists
  (STRING_LIST_T* a_list,
   STRING_LIST_T* b_list);

/*************************************************************************
 * Test whether a given string list contains duplicates.
 *
 * Optionally, print the duplicated strings to stderr.
 *************************************************************************/
BOOLEAN_T has_duplicates
  (char*          message, // If message == "", then don't print the list to stderr.
   STRING_LIST_T* my_list);

#endif



