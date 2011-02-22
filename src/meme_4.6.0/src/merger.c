/**************************************************************************
 * FILE: merger.c
 * AUTHOR: James Johnson 
 * CREATE DATE: 19-August-2009 
 * PROJECT: shared
 * COPYRIGHT: TBA 
 * VERSION: $Revision: 1.0 $
 * DESCRIPTION: Data structure for iteration through multiple sorted lists 
 **************************************************************************/
#include <assert.h>
#include "merger.h"
#include "utils.h"

/*
 * this points to a list with 'remaining' items
 * it is only used in association with the merger_t
 * which stores the information about the size of
 * elements.
 */
typedef struct merger_list_t MERGER_LIST_T;
struct merger_list_t {
  int remaining;
  void *list;
};

/*
 * merger_t maintains a heap of 'num_lists' sorted lists.
 * If 'num_lists' > 'list_to_insert' then not all lists
 * have been added and 'lists' is not yet heap ordered.
 * 'list_to_insert' maintains the next position a list
 * should be added in 'lists'.
 * The lists in 'lists' contain elements of 'elem_size' bytes.
 * A list is considered smaller than another list
 * if the first element of the smaller list is smaller
 * than the first element of the larger list.
 * Elements of the list are compared using the 'compar'
 * function which works exactly the same as the compar
 * to qsort. 
 */
struct merger_t {
  int num_lists;
  int list_to_insert;
  size_t elem_size;
  int(*compar)(const void *, const void *);
  MERGER_LIST_T *lists;
};
/*
 * Creates a merger
 */
MERGER_T *merger_create(int num_lists, size_t elem_size, int(*compar)(const void *, const void *)) {
  MERGER_T *merger = mm_malloc(sizeof(MERGER_T));
  merger->num_lists = num_lists;
  merger->list_to_insert = 0;
  merger->elem_size = elem_size;
  merger->compar = compar;
  merger->lists = (MERGER_LIST_T*)mm_malloc(sizeof(MERGER_LIST_T)*num_lists);
  return merger;
}

/*
 * Destroys a merger
 */
void merger_destroy(MERGER_T *merger) {
  assert(merger != NULL);
  free(merger->lists);
  free(merger);
}

/*
 * Takes a merger representing a heap and an index to check for
 * the heap property. It will fix the heap property ensuring that
 * the values of the heap elements are always larger than their
 * child elements.
 */
#define IN_BOUNDS(INDEX) ((INDEX) < merger->num_lists)
#define CMP(L_INDEX, R_INDEX) (merger->compar((merger->lists[L_INDEX].list),merger->lists[R_INDEX].list) < 0)
#define SWAP_LSTS(L_INDEX, R_INDEX) { \
    MERGER_LIST_T temp = merger->lists[L_INDEX]; \
    merger->lists[L_INDEX] = merger->lists[R_INDEX]; \
    merger->lists[R_INDEX] = temp; \
  }
void sift_down(int index, MERGER_T *merger) 
{
  assert(IN_BOUNDS(index));
  int leftI = index*2 + 1, rightI = index*2 + 2;
  if (IN_BOUNDS(rightI)) {
    //node has both children, find the smallest to consider a swap
    if (CMP(leftI, rightI)) { //left is smaller than right
      if (CMP(leftI, index)) { //swap with left
        SWAP_LSTS(leftI, index); 
        sift_down(leftI, merger); //fix left heap
      }
    } else { //either right is smaller than left or they are equal
      if (CMP(rightI, index)) { //swap with right
        SWAP_LSTS(rightI, index); 
        sift_down(rightI, merger); //fix right heap
      }
    }
  } else if (IN_BOUNDS(leftI)) { //only left child
      if (CMP(leftI, index)) { //swap with left
        SWAP_LSTS(leftI, index);
        sift_down(leftI, merger); //fix left heap
      }
  } 
}

/*
 * Adds a list to the merger, dies if too many lists are added
 */
void merger_add(void *list, int list_elems, MERGER_T *merger) {
  assert(merger != NULL);
  assert(merger->list_to_insert < merger->num_lists);
  if (list_elems > 0) {
    assert(list != NULL);
    MERGER_LIST_T *mlst = merger->lists+(merger->list_to_insert);
    mlst->list = list;
    mlst->remaining = list_elems;
    merger->list_to_insert += 1;
  } else {
    // they added an empty list
    merger->num_lists -= 1;
  }
  if (merger->list_to_insert == merger->num_lists) {
    int i; //order the lists so they have heap property
    for (i = ((merger->num_lists) >> 1) - 1; i >= 0; --i) {
      sift_down(i, merger);
    }
  }
}


/*
 * Gets the next item, but doesn't change the next item
 */
void *merger_peek(MERGER_T *merger) {
  assert(merger != NULL);
  assert(merger->list_to_insert >= merger->num_lists);
  if (merger->num_lists == 0) return NULL;
  return merger->lists[0].list;
}

/*
 * Gets the next item from the merged lists
 */
void *merger_next(MERGER_T *merger) {
  assert(merger != NULL);
  assert(merger->list_to_insert >= merger->num_lists);
  if (merger->num_lists == 0) return NULL;
  void *item = merger->lists[0].list;
  // update the list pointer to the next item
  merger->lists[0].list = (void *)(((char *)merger->lists[0].list)+(merger->elem_size));
  // deincrement the remaining items count
  merger->lists[0].remaining -= 1;
  // check if the list is now empty
  if (merger->lists[0].remaining == 0) {
    // deincrement the number of lists
    merger->num_lists -= 1;
    // move the last list to take its place
    merger->lists[0] = merger->lists[merger->num_lists];
    // restore the heap property
    if (merger->num_lists > 0) sift_down(0, merger);
  }
  return item;
}
