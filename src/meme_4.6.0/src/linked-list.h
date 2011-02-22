/**************************************************************************
 * FILE: linked-list.h
 * AUTHOR: James Johnson 
 * CREATE DATE: 17-August-2009 
 * PROJECT: shared
 * COPYRIGHT: TBA 
 * VERSION: $Revision: 1.0 $
 * DESCRIPTION: Data structure for manipulating a list
 **************************************************************************/

#ifndef LINKED_LIST_H
#define LINKED_LIST_H

/*
 * Types
 */
typedef struct linklst_t LINKLST_T;
typedef struct link_t LINK_T;


/*
 * linklst_create
 * creates and returns an empty linked list
 */
LINKLST_T *linklst_create();

/*
 * linklst_destroy_all
 * destroys a linked list using the passed function 
 * to deallocate all contained items
 */
void linklst_destroy_all(LINKLST_T *linklst, void(*free_item)(void *));

/*
 * linklst_destroy
 * destroys a linked list
 * caller must deallocate items themselves
 */
void linklst_destroy(LINKLST_T *linklst);

/*
 * linklst_size
 * returns the number of items in the list
 */
int linklst_size(LINKLST_T *linklst);

/*
 * linklst_first
 * gets the first link
 */
LINK_T *linklst_first(LINKLST_T *linklst);

/*
 * linklst_last
 * gets the last link
 */
LINK_T *linklst_last(LINKLST_T *linklst);

/*
 * linklst_add_after
 * adds an item after the specified link
 */
LINK_T *linklst_add_after(void *item, LINK_T *before, LINKLST_T *linklst);

/*
 * linklst_add_before
 * adds an item before the specified link
 */
LINK_T *linklst_add_before(void *item, LINK_T *after, LINKLST_T *linklst);

/*
 * linklst_add
 * adds an item at the end of the list
 */
LINK_T *linklst_add(void *item, LINKLST_T *linklst);

/*
 * linklst_take
 * removes an item from the end of the list
 * and returns a pointer to the item
 */
void *linklst_take(LINKLST_T *linklst);

/*
 * linklst_push
 * adds an item to the front of the list
 */
LINK_T *linklst_push(void *item, LINKLST_T *linklst);

/*
 * linklst_pop
 * removes an item from the front of the list
 * and returns a pointer to the item
 */
void *linklst_pop(LINKLST_T *linklst);

/*
 * linklst_remove
 * removes the passed link from the list
 * and returns a pointer to the item.
 */
void *linklst_remove(LINK_T *link, LINKLST_T *linklst);

/*
 * linklst_sort
 * sorts the list using the passed
 * comparator to compare items.
 */
void linklst_sort(int (*comparator)(void*,void*), LINKLST_T *linklst);

/*
 * linklst_next
 * gets the next link
 */
LINK_T *linklst_next(LINK_T *link);

/*
 * linklst_prev
 * gets the previous link
 */
LINK_T *linklst_prev(LINK_T *link);

/*
 * linklst_get
 * gets the item contained in the link
 */
void *linklst_get(LINK_T *link);

/*
 * linklst_set
 * sets the item contained in the link
 */
void linklst_set(void *item, LINK_T *link);

#endif

