/********************************************************************
 * FILE: object-list.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: May Day, 2006 (Alaska #8, SEA -> Newark)
 * PROJECT: META-MEME
 * COPYRIGHT: 2006, UW
 ********************************************************************/

#ifndef OBJECT_LIST_H
#define OBJECT_LIST_H

#include "array.h"


/*************************************************************************
 * Primary data structure.
 *************************************************************************/
typedef struct object_list_t OBJECT_LIST_T;

/*************************************************************************
 * Create an empty object list.
 *************************************************************************/
OBJECT_LIST_T* new_object_list
  (BOOLEAN_T (*compare_keys)(),    // Function to compare two keys.
   void*     (*copy_key)(),        // Function to copy a key.
   void      (*free_key)(),        // Function to free one key.
   void      (*free_object)());    // Function to free one object.

/*************************************************************************
 * Store an object in a list.
 *************************************************************************/
void store_object
 (void*          object,
  void*          key,       // Key for use in retrieval.
  double         score,     // Score for use in sorting.
  OBJECT_LIST_T* an_object_list);

/*************************************************************************
 * Sort the objects in the list.
 *************************************************************************/
void sort_objects
  (OBJECT_LIST_T* an_object_list);

/*************************************************************************
 * Store a object in a list.
 *************************************************************************/
void* retrieve_object
 (void*          key,
  OBJECT_LIST_T* an_object_list);

/*************************************************************************
 * Get the next object in a list; saves state.
 * First call retrieves first object; successive calls retrieve next ones.
 * Returns NULL if at end of list.  
 * The next call after that will retrieve the first object again.
 *************************************************************************/
void* retrieve_next_object
 (OBJECT_LIST_T* an_object_list);

/*************************************************************************
 *  Extract an array of scores.  Memory must be freed by the caller.
 *************************************************************************/
ARRAY_T* get_object_list_scores
  (OBJECT_LIST_T* a_list);

/*************************************************************************
 * Free dynamic memory used by a list of objects.
 *************************************************************************/
void free_object_list
  (OBJECT_LIST_T* a_list);

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
