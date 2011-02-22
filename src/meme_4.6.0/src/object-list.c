/********************************************************************
 * FILE: object-list.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: May 3, 2006
 * PROJECT: META-MEME
 * COPYRIGHT: 2006, UW
 ********************************************************************/

#include "object-list.h"

/*************************************************************************
 * Primary data structure.
 *************************************************************************/
// This second level struct allows easy sorting.
typedef struct object_t {
  void*     key;        // Used in retrieval.
  double    score;      // Used in sorting.
  void*     value;      // The object itself.
} OBJECT_T;

struct object_list_t {
  int        num_objects; // Number of objects in the array.
  int        max_objects; // Total amount of memory allocated.
  int        next_object; // Next object to retrieve.
  OBJECT_T*  objects;     // The list of objects.
  BOOLEAN_T (*compare_keys)();    // Function to compare two keys.
  void*     (*copy_key)();        // Function to copy a key.
  void      (*free_key)();        // Function to free one key.
  void      (*free_object)();     // Function to free one object.
};

/*************************************************************************
 * Create an empty object list.
 *
 * If you will not be retrieving by key, you can give NULL for
 * compare_keys, copy_key and free_key.
 *************************************************************************/
#define DEFAULT_MAX_OBJECTS 10
OBJECT_LIST_T* new_object_list
  (BOOLEAN_T (*compare_keys)(),    // Function to compare two keys.
   void*     (*copy_key)(),        // Function to copy a key.
   void      (*free_key)(),        // Function to free one key.
   void      (*free_object)())     // Function to free one object.
{
  OBJECT_LIST_T* new_object_list;

  // Allocate memory.
  new_object_list = (OBJECT_LIST_T*)mm_calloc(1, sizeof(OBJECT_LIST_T));
  new_object_list->num_objects = 0;
  new_object_list->max_objects = DEFAULT_MAX_OBJECTS;
  new_object_list->next_object = 0;
  new_object_list->objects = 
    (OBJECT_T*)mm_calloc(DEFAULT_MAX_OBJECTS, sizeof(OBJECT_T));

  // Store the functions for freeing and comparing.
  new_object_list->compare_keys = compare_keys;
  new_object_list->copy_key = copy_key;
  new_object_list->free_key = free_key;
  new_object_list->free_object = free_object;

  return(new_object_list);
}

/*************************************************************************
 * Store an object in a list.
 *
 * It is legal to give a NULL key; just don't try to retrieve by key later.
 *************************************************************************/
void store_object
 (void*          object,
  void*          key,       // Key for use in retrieval.
  double         score,     // Score for use in sorting.
  OBJECT_LIST_T* an_object_list)
{
  // Make sure none of the inputs is null.
  if (object == NULL) {
    die("Adding null OBJECT to list.");
  }
  if (an_object_list == NULL) {
    die("Adding OBJECT to null list.");
  }

  // If it's already there, don't do anything.
  if ((key != NULL) && (retrieve_object(key, an_object_list) != NULL)) {
    fprintf(stderr, "Adding OBJECT that was already in the list.\n");
    return;
  }

  // Reallocate space if there isn't any.
  if (an_object_list->num_objects >= an_object_list->max_objects) {
    an_object_list->max_objects += DEFAULT_MAX_OBJECTS;
    an_object_list->objects = (OBJECT_T*)mm_realloc(an_object_list->objects,
						    an_object_list->max_objects 
						    * sizeof(OBJECT_T));
  }

  // Add to the list.
  if (key != NULL) {
    an_object_list->objects[an_object_list->num_objects].key
      = an_object_list->copy_key(key);
  } else {
    an_object_list->objects[an_object_list->num_objects].key = NULL;
  }
  an_object_list->objects[an_object_list->num_objects].score = score;
  an_object_list->objects[an_object_list->num_objects].value = object;
  (an_object_list->num_objects)++;
}

/*************************************************************************
 * Sort the objects in the list.
 *************************************************************************/
static int object_compare
  (const void* elem1,
   const void* elem2)
{
  const double key1 = ((OBJECT_T *)elem1)->score;
  const double key2 = ((OBJECT_T *)elem2)->score;

  if (key1 <= key2) {
    return(-1);
  } else {
    return(1);
  }
  return(0);
}

void sort_objects
  (OBJECT_LIST_T* an_object_list)
{
  qsort((void *)(an_object_list->objects), 
	an_object_list->num_objects, 
	sizeof(OBJECT_T),
	object_compare);

  // Reset the counter to zero.
  an_object_list->next_object = 0;
}

/*************************************************************************
 * Retrieve an object from a list.
 *************************************************************************/
void* retrieve_object
 (void*          key,
  OBJECT_LIST_T* an_object_list)
{
  int i_object;

  if (key == NULL) {
    die("Tried to retrieve by NULL key.");
  }

  for (i_object = 0; i_object < an_object_list->num_objects; i_object++) {
    void* this_key = an_object_list->objects[i_object].key;
    if (this_key == NULL) {
      die("Tried to retrieve by key from a list containing a NULL key.");
    }
    if ((an_object_list->compare_keys)(key, this_key)) {
      return(an_object_list->objects[i_object].value);
    }
  }
  // If we didn't find it, return null pointer.
  return(NULL);
}

/*************************************************************************
 * Get the next object in a list; saves state.  
 * First call retrieves first object; successive calls retrieve next ones.
 * Returns NULL if at end * of list.  
 * The next call after that will retrieve the first object again.
 *************************************************************************/
void* retrieve_next_object
 (OBJECT_LIST_T* an_object_list)
{
  if (an_object_list->next_object == an_object_list->num_objects) {
    an_object_list->next_object = 0;		// rewind
    return(NULL);
  }
  else {
    return(an_object_list->objects[an_object_list->next_object++].value);
  }
}

/*************************************************************************
 *  Extract an array of scores.  Memory must be freed by the caller.
 *************************************************************************/
ARRAY_T* get_object_list_scores
  (OBJECT_LIST_T* a_list)
{
  ARRAY_T* return_value = allocate_array(a_list->num_objects);

  int i_object;
  for (i_object = 0; i_object < a_list->num_objects; i_object++) {
    set_array_item(i_object, a_list->objects[i_object].score, return_value);
  }
  return(return_value);
}

/*************************************************************************
 * Free dynamic memory used by a list of objects.
 *************************************************************************/
void free_object_list
  (OBJECT_LIST_T* a_list)
{
  int i_object;

  if (a_list == NULL) {
    return;
  }

  for (i_object = 0; i_object < a_list->num_objects; i_object++) {
    if (a_list->free_key != NULL) {
      (a_list->free_key)(a_list->objects[i_object].key);
    }
    (a_list->free_object)(a_list->objects[i_object].value);
  }
  myfree(a_list->objects);
  myfree(a_list->objects);
  myfree(a_list);
}

/*
 * Local Variables:
 * mode: c++
 * c-basic-offset: 2
 * End:
 */
