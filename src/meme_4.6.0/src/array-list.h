/**************************************************************************
 * FILE: array-list.h
 * AUTHOR: James Johnson 
 * CREATE DATE: 26-August-2009 
 * PROJECT: shared
 * COPYRIGHT: 2009, UQ 
 * VERSION: $Revision: 1.0 $
 * DESCRIPTION: Data structure for manipulating a list
 **************************************************************************/
#ifndef array_list_h
#define array_list_h

#include "utils.h"

typedef struct arraylst_t ARRAYLST_T;

/*
 * arraylst_create_sized
 * Creates an array list with an initial allocated
 * size as passed.
 */
ARRAYLST_T *arraylst_create_sized(int size);

/*
 * arraylst_create
 * Create an array list with an initial allocated 
 * size of ten.
 */
ARRAYLST_T *arraylst_create();

/*
 * arraylst_destroy
 * Destroys the array list. Will optionally destroy the contained items if
 * passed a destructor function. For simple things passing free will surfice.
 */
void arraylst_destroy(void (*optional_item_destructor)(void*), ARRAYLST_T *arraylst);

/*
 * arraylst_preallocate
 * Ensures that space up to the passed size is avaliable so
 * that add operations won't cause unnecessary expansions. 
 * If this has been called than remove opperations won't 
 * cause a shrink smaller than this size until this method 
 * is called again or arraylst_fit is called.
 */
void arraylst_preallocate(int size, ARRAYLST_T *arraylst);

/*
 * arraylst_fit
 * Ensures that the allocated size is exactly the currently
 * utilized size. If arraylst_preallocate has been called
 * then this cancels its requirement.
 */
void arraylst_fit(ARRAYLST_T *arraylst);

/*
 * arraylst_is_empty
 * Returns true if the array list is empty.
 */
BOOLEAN_T arraylst_is_empty(ARRAYLST_T *arraylst);

/*
 * arraylst_size
 * Gets the number of items added to the list. 
 */
int arraylst_size(ARRAYLST_T *arraylst); 

/*
 * arraylst_put_n
 * The passed item is put at the passed index the passed
 * number of times. The index must be from 0 to the size
 * of the list inclusive or the opperation will die.
 */
void arraylst_put_n(int times, int index, void *item, ARRAYLST_T *arraylst); 

/*
 * arraylst_add_n
 * The passed item is added at the end of the list the
 * passed number of times.
 */
void arraylst_add_n(int times, void *item, ARRAYLST_T *arraylst); 

/*
 * arraylst_put
 * The passed item is put at the passed index.
 * The passed index must be from 0 to the size of the list
 * inclusive or the opperation will die.
 */
void arraylst_put(int index, void *item, ARRAYLST_T *arraylst); 

/*
 * arraylst_add
 * The passed item is added to the end of the list.
 */
void arraylst_add(void *item, ARRAYLST_T *arraylst);

/*
 * arraylst_set
 * The passed item is set to be at the passed index.
 */
void arraylst_set(int index, void *item, ARRAYLST_T *arraylst);

/*
 * arraylst_swap
 * The items at the passed indexes are swapped
 */
void arraylst_swap(int index1, int index2, ARRAYLST_T *arraylst);

/*
 * arraylst_get
 * Gets the item at the passed index. Dies if the index is
 * out of bounds.
 */
void *arraylst_get(int index, ARRAYLST_T *arraylst);

/*
 * arraylst_remove_range
 * The items in the range are removed from the list
 * and the item in the smallest index is returned.
 */
void *arraylst_remove_range(int index1, int index2, void (*optional_item_destructor)(void*), ARRAYLST_T *arraylst);

/*
 * arraylst_clear
 * Clear the array list of any items.
 */
void arraylst_clear(void (*optional_item_destructor)(void*), ARRAYLST_T *arraylst);

/*
 * arraylst_remove
 * The item at the passed index is removed from the list
 * and returned. Dies if the index is out of bounds.
 */
void *arraylst_remove(int index, ARRAYLST_T *arraylst);

/*
 * arraylst_take
 * The last item is removed from the list and returned.
 * Dies if there are no more items.
 */
void *arraylst_take(ARRAYLST_T *arraylst);

/*
 * arraylst_map_range
 * Applies a function to each item in a range of the array list and updates 
 * the item to what was returned from the function.
 */
void arraylst_map_range(void* (*map_fun)(void*), int index, int count, ARRAYLST_T *arraylst);

/*
 * arraylst_map
 * Applies a function to each item of the array list and updates it to what 
 * was returned.
 */
void arraylst_map(void* (*map_fun)(void*), ARRAYLST_T *arraylst);

/*
 * arraylst_apply_range
 * Applies a function to each item in a range of the array list.
 */
void arraylst_apply_range(void (*fun)(void*), int index, int count, ARRAYLST_T *arraylst);

/*
 * arraylst_apply
 * Applies a function to each item of the array list.
 */
void arraylst_apply(void (*fun)(void*), ARRAYLST_T *arraylst);

/*
 * arraylst_accumulate_range
 * Applies a function to each item in a range of the array list combining it 
 * with an accumulated value. The accumulated value is returned.
 */
void* arraylst_accumulate_range(void (*accumulator_fun)(void*, void*), void *initval, int index, int count, ARRAYLST_T *arraylst);

/*
 * arraylst_accumulate
 * Applies a function to each item of the array list combining it with an 
 * accumulated value. The accumulated value is returned.
 */
void* arraylst_accumulate(void (*accumulator_fun)(void*, void*), void *initval, ARRAYLST_T *arraylst);

/*
 * arraylst_qsort
 * Sorts the array list using quick sort with the passed comparator.
 * Due to implementation (with qsort function) the comparator 
 * must take pointers to pointers to items.
 */
void arraylst_qsort(int(*compar)(const void *, const void *), ARRAYLST_T *arraylst);

/*
 * arraylst_bsearch
 * Does a binary search on the sorted array with the passed key and
 * comparison function. Returns the index of the an item that
 * equals the key according to the comparison function. Otherwise returns
 * the -(insert_position + 1) of the item.
 * The first item to compar will be a pointer to the address of the key 
 * the second item will be a pointer to the address of an item
 */
int arraylst_bsearch(const void *key, int (*compar)(const void *, const void *), ARRAYLST_T *arraylst);

/*
 * arraylst_compar_txt
 * Can be passed to arraylst_qsort or arraylst_bsearch as a comparator when
 * the items in the array and key are char*.
 */
int arraylst_compar_txt(const void *v1, const void *v2);

#endif
