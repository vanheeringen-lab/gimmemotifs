/**************************************************************************
 * FILE: array-list.c
 * AUTHOR: James Johnson 
 * CREATE DATE: 26-August-2009 
 * PROJECT: shared
 * COPYRIGHT: 2009, UQ 
 * VERSION: $Revision: 1.0 $
 * DESCRIPTION: Data structure for manipulating a list
 **************************************************************************/
#include <assert.h>
#include <stdlib.h>

#include "array-list.h"
#include "binary-search.h"

#define min(a,b)      (a<b)?a:b
#define max(a,b)      (a<b)?b:a

/*
 * arraylst_t
 * Stores an array of pointers to items.
 * Member 'array' will have an allocated size of
 * 'cur_alloc' void pointers. Member 'min_alloc'
 * contains the lower bound for 'cur_alloc' as the
 * array must not allocate less than 'min_alloc'
 * pointers even if the actual utilized pointer
 * count is much less.
 */
struct arraylst_t 
{
  void **array;
  int cur_alloc;
  int min_alloc;
  int size;
};

/*
 * arraylst_create_sized
 * Creates an array list with an initial allocated
 * size as passed.
 */
ARRAYLST_T *arraylst_create_sized(int size) 
{
  //test parameters
  if (size < 0) die("arraylst_create_sized: size must be zero or larger.\n");
  ARRAYLST_T *arraylst = (ARRAYLST_T*)mm_malloc(sizeof(ARRAYLST_T));
  if (size > 0) arraylst->array = (void**)mm_malloc(sizeof(void*)*size);
  arraylst->cur_alloc = size;
  arraylst->min_alloc = size;
  arraylst->size = 0;
  return arraylst;
}

/*
 * arraylst_create
 * Create an array list with an initial allocated 
 * size of ten.
 */
ARRAYLST_T *arraylst_create() {
  return arraylst_create_sized(10);
}


/*
 * arraylst_destroy
 * Destroys the array list. Will optionally destroy the contained items if
 * passed a destructor function. For simple things passing free will surfice.
 */
void arraylst_destroy(void (*optional_item_destructor)(void*), ARRAYLST_T *arraylst) {
  //test parameters
  if (arraylst == NULL) die("arraylst_destroy: arraylst is null!\n");

  if (optional_item_destructor != NULL) arraylst_apply(optional_item_destructor, arraylst);

  if (arraylst->cur_alloc > 0) free(arraylst->array);
  free(arraylst);
}

/*
 * arraylst_preallocate
 * Ensures that space up to the passed size is avaliable so
 * that add operations won't cause unnecessary expansions. 
 * If this has been called than remove opperations won't 
 * cause a shrink smaller than this size until this method 
 * is called again or arraylst_fit is called.
 */
void arraylst_preallocate(int size, ARRAYLST_T *arraylst) {
  //test parameters
  if (arraylst == NULL) die("arraylst_preallocate: arraylst is null!\n");
  if (size < 0)         die("arraylst_preallocate: size is smaller than zero.\n");

  arraylst->min_alloc = size;
  if (arraylst->min_alloc > arraylst->cur_alloc) {
    //expand
    if (arraylst->cur_alloc == 0) {
      arraylst->array = (void**)mm_malloc(sizeof(void*)*arraylst->min_alloc);
    } else {
      arraylst->array = (void**)mm_realloc(arraylst->array, sizeof(void*)*arraylst->min_alloc);
    }
    arraylst->cur_alloc = arraylst->min_alloc;
  } else if (arraylst->size < (arraylst->cur_alloc / 4) && arraylst->min_alloc < arraylst->cur_alloc) {
    //shrink
    int require = max(arraylst->size * 2, arraylst->min_alloc);
    if (require == 0) {
      free(arraylst->array);
    } else {
      arraylst->array = (void**)mm_realloc(arraylst->array, sizeof(void*)*require);
    }
    arraylst->cur_alloc = require;
  }
}

/*
 * arraylst_fit
 * Ensures that the allocated size is exactly the currently
 * utilized size. If arraylst_preallocate has been called
 * then this cancels its requirement.
 */
void arraylst_fit(ARRAYLST_T *arraylst) {
  //test parameters
  if (arraylst == NULL) die("arraylst_fit: arraylst is null!\n");

  arraylst->min_alloc = 0;
  if (arraylst->cur_alloc > arraylst->size) {
    if (arraylst->size == 0) {
      free(arraylst->array);
    } else {
      arraylst->array = (void**)mm_realloc(arraylst->array, sizeof(void*)*arraylst->size);
    }
    arraylst->cur_alloc = arraylst->size;
  }
}

/*
 * arraylst_is_empty
 * Returns true if the array list is empty.
 */
BOOLEAN_T arraylst_is_empty(ARRAYLST_T *arraylst) {
  //test parameters
  if (arraylst == NULL) die("arraylst_is_empty: arraylst is null!\n");

  return arraylst->size == 0;
}

/*
 * arraylst_size
 * Gets the number of items added to the list. 
 */
int arraylst_size(ARRAYLST_T *arraylst) {
  //test parameters
  if (arraylst == NULL) die("arraylst_size: arraylst is null!\n");

  return arraylst->size;
}

/*
 * arraylst_put_n
 * The passed item is put at the passed index the passed
 * number of times. The index must be from 0 to the size
 * of the list inclusive or the opperation will die.
 */
void arraylst_put_n(int times, int index, void *item, ARRAYLST_T *arraylst) {
  //test parameters
  if (arraylst == NULL)                     die("arraylst_put_n: arraylst is null!\n");
  if (index < 0 || index > arraylst->size)  die("arraylst_put_n: index must be between zero and size inclusive.\n");
  if (times <= 0)                           die("arraylst_put_n: times must be larger than zero.\n");
  //ensure we have enough space
  int new_size = arraylst->size + times;
  if (new_size > arraylst->cur_alloc) {
    int alloc = max(arraylst->cur_alloc, 1);
    do {
      alloc *= 2;
    } while (alloc < new_size);
    if (arraylst->cur_alloc == 0) {
      arraylst->array = mm_malloc(sizeof(void*)*alloc);
    } else {
      arraylst->array = mm_realloc(arraylst->array, sizeof(void*)*alloc);
    }
    arraylst->cur_alloc = alloc;
  }
  //move anything in the way
  if (index < arraylst->size) {
    memmove(arraylst->array+(index+times), arraylst->array+(index), sizeof(void*)*(arraylst->size - index));
  }
  //set the item
  int i;
  for (i = 0; i < times; ++i) {
    arraylst->array[index + i] = item;
  }
  arraylst->size = new_size;
}

/*
 * arraylst_add_n
 * The passed item is added at the end of the list the
 * passed number of times.
 */
void arraylst_add_n(int times, void *item, ARRAYLST_T *arraylst) {
  arraylst_put_n(times, arraylst->size, item, arraylst);
}

/*
 * arraylst_put
 * The passed item is put at the passed index.
 * The passed index must be from 0 to the size of the list
 * inclusive or the opperation will die.
 */
void arraylst_put(int index, void *item, ARRAYLST_T *arraylst) {
  arraylst_put_n(1, index, item, arraylst);
}

/*
 * arraylst_add
 * The passed item is added to the end of the list.
 */
void arraylst_add(void *item, ARRAYLST_T *arraylst) {
  arraylst_put_n(1, arraylst->size, item, arraylst);
}

/*
 * arraylst_set
 * The passed item is set to be at the passed index.
 */
void arraylst_set(int index, void *item, ARRAYLST_T *arraylst) {
  // test parameters
  if (arraylst == NULL)                     die("arraylst_set: arraylst is NULL!\n");
  if (index < 0 || index >= arraylst->size) die("arraylst_set: index must be within bounds.\n");
  arraylst->array[index] = item;
}

/*
 * arraylst_swap
 * The items at the passed indexes are swapped
 */
void arraylst_swap(int index1, int index2, ARRAYLST_T *arraylst) {
  void *temp;
  // test parameters
  if (arraylst == NULL)                       die("arraylst_set: arraylst is NULL!\n");
  if (index1 < 0 || index1 >= arraylst->size) die("arraylst_set: index1 must be within bounds.\n");
  if (index2 < 0 || index2 >= arraylst->size) die("arraylst_set: index2 must be within bounds.\n");
  temp = arraylst->array[index1];
  arraylst->array[index1] = arraylst->array[index2];
  arraylst->array[index2] = temp;
}

/*
 * arraylst_get
 * Gets the item at the passed index. Dies if the index is
 * out of bounds.
 */
void *arraylst_get(int index, ARRAYLST_T *arraylst) {
  // test parameters
  if (arraylst == NULL)                     die("arraylst_get: arraylst is NULL!\n");
  if (index < 0 || index >= arraylst->size) die("arraylst_get: index must be within bounds.\n");
  return arraylst->array[index];
}

/*
 * arraylst_remove_range
 * The items in the range are removed from the list and if a destructor is specified they are destroyed.
 * If a destructor is not specified the item that was at index is returned otherwise NULL is returned.
 */
void *arraylst_remove_range(int index, int count, void (*optional_item_destructor)(void*), ARRAYLST_T *arraylst) {
  void *item; 
  // test parameters
  if (arraylst == NULL)                       die("arraylst_remove_range: arraylst is NULL!\n"); 
  if (count < 1)                              die("arraylst_remove_range: count must be one or more elements.\n"); 
  if (index < 0 || index >= arraylst->size)   die("arraylst_remove_range: index must be within bounds.\n");
  if(index + count > arraylst->size)          die("arraylst_remove_range: index + count is larger than size!\n");
  
  // destroy the items if we can
  if (optional_item_destructor != NULL) {
    item = NULL;
    arraylst_apply(optional_item_destructor, arraylst);
  } else {
    // get the item to return
    item = arraylst->array[index];
  }
  // move the pointers after the removed range in to fill the deleted spot
  if (index + count < arraylst->size) {
    memmove(arraylst->array+(index), arraylst->array+(index + count), sizeof(void*)*(arraylst->size - index - count));
  }
  // update the size
  arraylst->size -= count;
  // shrink the array if the utilized portion drops to less than a quarter of the allocated size
  // unless we've already reached the minimum allowed allocation
  if (arraylst->size < arraylst->cur_alloc / 4 && arraylst->min_alloc < arraylst->cur_alloc) {
    int alloc = max(arraylst->size * 2, arraylst->min_alloc);
    if (alloc == 0) {
      free(arraylst->array);
    } else {
      arraylst->array = (void**)mm_realloc(arraylst->array, sizeof(void*)*alloc);
    }
    arraylst->cur_alloc = alloc;
  }
  return item;
}

/*
 * arraylst_clear
 * Clear the array list of any items.
 */
void arraylst_clear(void (*optional_item_destructor)(void*), ARRAYLST_T *arraylst) {
  // test parameters
  if (arraylst == NULL) die("arraylst_clear: arraylst is NULL!\n"); 
  if (arraylst->size > 0) arraylst_remove_range(0, arraylst->size, optional_item_destructor, arraylst);
}

/*
 * arraylst_remove
 * The item at the passed index is removed from the list
 * and returned. Dies if the index is out of bounds.
 */
void *arraylst_remove(int index, ARRAYLST_T *arraylst) {
  return arraylst_remove_range(index, 1, NULL, arraylst);
}

/*
 * arraylst_take
 * The last item is removed from the list and returned.
 * Dies if there are no more items.
 */
void *arraylst_take(ARRAYLST_T *arraylst) {
  if (arraylst->size == 0) die("arraylst_take: No more items to take from arraylst!\n");
  return arraylst_remove_range(arraylst->size-1, 1, NULL, arraylst);
}

/*
 * arraylst_map_range
 * Applies a function to each item in a range of the array list and updates 
 * the item to what was returned from the function.
 */
void arraylst_map_range(void* (*map_fun)(void*), int index, int count, ARRAYLST_T *arraylst) {
  int i;
  if (arraylst == NULL)                       die("arraylst_map_range: arraylst is NULL!\n"); 
  if (map_fun == NULL)                        die("arraylst_map_range: map_fun is NULL!\n");
  if (count < 0)                              die("arraylst_map_range: count must be zero or more elements.\n"); 
  if (index < 0 || index > arraylst->size)    die("arraylst_map_range: index must be within bounds.\n");
  if(index + count > arraylst->size)          die("arraylst_map_range: index + count is larger than size!\n");
  for (i = index; i < count; ++i) {
    arraylst->array[i] = map_fun(arraylst->array[i]);
  }
}

/*
 * arraylst_map
 * Applies a function to each item of the array list and updates it to what 
 * was returned.
 */
void arraylst_map(void* (*map_fun)(void*), ARRAYLST_T *arraylst) {
  int i;
  if (arraylst == NULL)                       die("arraylst_map: arraylst is NULL!\n"); 
  if (map_fun == NULL)                        die("arraylst_map: map_fun is NULL!\n");
  for (i = 0; i < arraylst->size; ++i) {
    arraylst->array[i] = map_fun(arraylst->array[i]);
  }
}

/*
 * arraylst_apply_range
 * Applies a function to each item in a range of the array list.
 */
void arraylst_apply_range(void (*fun)(void*), int index, int count, ARRAYLST_T *arraylst) {
  int i;
  if (arraylst == NULL)                       die("arraylst_map_range: arraylst is NULL!\n"); 
  if (fun == NULL)                            die("arraylst_map_range: fun is NULL!\n");
  if (count < 0)                              die("arraylst_map_range: count must be zero or more elements.\n"); 
  if (index < 0 || index > arraylst->size)    die("arraylst_map_range: index must be within bounds.\n");
  if(index + count > arraylst->size)          die("arraylst_map_range: index + count is larger than size!\n");
  for (i = index; i < count; ++i) {
    fun(arraylst->array[i]);
  }
}

/*
 * arraylst_apply
 * Applies a function to each item of the array list.
 */
void arraylst_apply(void (*fun)(void*), ARRAYLST_T *arraylst) {
  int i;
  if (arraylst == NULL)                       die("arraylst_apply: arraylst is NULL!\n"); 
  if (fun == NULL)                            die("arraylst_apply: fun is NULL!\n");
  for (i = 0; i < arraylst->size; ++i) {
    fun(arraylst->array[i]);
  }
}

/*
 * arraylst_accumulate_range
 * Applies a function to each item in a range of the array list combining it 
 * with an accumulated value. The accumulated value is returned.
 */
void* arraylst_accumulate_range(void (*accumulator_fun)(void*, void*), void *initval, int index, int count, ARRAYLST_T *arraylst) {
  int i;
  void *val;
  if (arraylst == NULL)                       die("arraylst_accumulate_range: arraylst is NULL!\n"); 
  if (accumulator_fun == NULL)                die("arraylst_accumulate_range: accumulator_fun is NULL!\n");
  if (count < 0)                              die("arraylst_accumulate_range: count must be zero or more elements.\n"); 
  if (index < 0 || index > arraylst->size)    die("arraylst_accumulate_range: index must be within bounds.\n");
  if(index + count > arraylst->size)          die("arraylst_accumulate_range: index + count is larger than size!\n");
  val = initval;
  for (i = index; i < count; ++i) {
    accumulator_fun(arraylst->array[i], val);
  }
  return val;
}

/*
 * arraylst_accumulate
 * Applies a function to each item of the array list combining it with an 
 * accumulated value. The accumulated value is returned.
 */
void* arraylst_accumulate(void (*accumulator_fun)(void*, void*), void *initval, ARRAYLST_T *arraylst) {
  int i;
  void *val;
  if (arraylst == NULL)                       die("arraylst_accumulate: arraylst is NULL!\n"); 
  if (accumulator_fun == NULL)                die("arraylst_accumulate: accumulator_fun is NULL!\n");
  val = initval;
  for (i = 0; i < arraylst->size; ++i) {
    accumulator_fun(arraylst->array[i], val);
  }
  return val;
}

/*
 * arraylst_qsort
 * Sorts the array list using quick sort with the passed comparator.
 * Due to implementation (with qsort function) the comparator must take pointers to pointers to items.
 */
void arraylst_qsort(int(*compar)(const void *, const void *), ARRAYLST_T *arraylst) {
  qsort(arraylst->array, arraylst->size, sizeof(void*), compar);
}

/*
 * arraylst_bsearch
 * Does a binary search on the sorted array with the passed key and
 * comparison function. Returns the index of the an item that
 * equals the key according to the comparison function. Otherwise returns
 * the -(insert_position + 1) of the item.
 * The first item to compar will be a pointer to the address of the key 
 * the second item will be a pointer to the address of an item
 */
int arraylst_bsearch(const void *key, int (*compar)(const void *, const void *), ARRAYLST_T *arraylst) {
  return binary_search(&key, arraylst->array, arraylst->size, sizeof(void*), compar); 
}

/*
 * arraylst_compar_txt
 * Can be passed to arraylst_qsort or arraylst_bsearch as a comparator when
 * the items in the array and key are char*.
 */
int arraylst_compar_txt(const void *v1, const void *v2) {
  char* s1 = *((char**)v1);
  char* s2 = *((char**)v2);
  return strcmp(s1, s2);
}


