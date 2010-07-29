/**************************************************************************
 * FILE: merger.h
 * AUTHOR: James Johnson 
 * CREATE DATE: 19-August-2009 
 * PROJECT: shared
 * COPYRIGHT: TBA 
 * VERSION: $Revision: 1.0 $
 * DESCRIPTION: Data structure for iteration through multiple sorted lists 
 **************************************************************************/
#ifndef MERGER_H
#define MERGER_H
#include <sys/types.h>

typedef struct merger_t MERGER_T;

/*
 * merger_create
 * takes as a parameter the number of lists to merge,
 * the size of a list item and a comparison function
 */
MERGER_T *merger_create(int num_lists, size_t elem_size, int(*compar)(const void *, const void *));

/*
 * merger_destroy
 * destroys the passed merger
 */
void merger_destroy(MERGER_T *merger);

/*
 * merger_add
 * add a list to merge and the number of items in the list
 */
void merger_add(void *list, int list_elems, MERGER_T *merger);

/*
 * merger_peek
 * returns a pointer to the next item but doesn't step past it
 */
void *merger_peek(MERGER_T *merger);

/*
 * merger_next
 * returns a pointer to the next item or
 * null if all items have been returned
 */
void *merger_next(MERGER_T *merger);

#endif
