/**************************************************************************
 * FILE: red-black-tree.h
 * AUTHOR: James Johnson 
 * CREATE DATE: 08-September-2009 
 * PROJECT: shared
 * COPYRIGHT: UQ, 2009
 * VERSION: $Revision: 1.0 $
 * DESCRIPTION: A red-black semi-balanced binary search tree. Usable as a set
 * or a map. Provides O(log(N)) lookup, insert and delete operations.
 **************************************************************************/

#ifndef RED_BLACK_TREE_H
#define RED_BLACK_TREE_H

#include "utils.h"

typedef struct rbtree_t RBTREE_T;
typedef struct rbnode_t RBNODE_T;

/*
 * rbtree_check
 * Debugging function. Checks that the tree is consistent with a red-black tree's structure.
 */
void rbtree_check(RBTREE_T *tree);

/*
 * rbtree_create
 * create a red-black tree.
 * takes a key comparator and optional key/value copy/free functions
 */
RBTREE_T *rbtree_create(
    int (*key_compare) (void*, void*),  // key comparison function, works just like the comparison function to qsort 
    void* (*key_copy)(void*),           // OPTIONAL key copying function, or NULL 
    void (*key_free) (void*),           // OPTIONAL key freeing function, or NULL
    void* (*value_copy) (void*),        // OPTIONAL value copying function or NULL
    void (*value_free) (void*)          // OPTIONAL value freeing function or NULL
);

/*
 * rbtree_destroy
 * destroys a red-black tree.
 * If any key or value free functions have been set then it will run 
 * them on the keys and values.
 */
void rbtree_destroy(RBTREE_T *tree);

/*
 * rbtree_alter_key_copy
 * Changes the function used to copy the key. Might be useful if you are
 * loading from multiple sources one of which you want to copy the key
 * and the other you just want to use it as is. If NULL is passed then
 * no attempt will be made to copy the keys.
 */
void rbtree_alter_key_copy(RBTREE_T *tree, void* (*key_copy)(void*));

/*
 * rbtree_alter_key_free
 * Changes the function used to free the key.
 * If NULL is passed then no attempt will be made to free the keys.
 */
void rbtree_alter_key_free(RBTREE_T *tree, void (*key_free)(void*));

/*
 * rbtree_alter_value_copy
 * Changes the function used to copy the value. Might be useful if you are
 * loading from multiple sources, one of which you want to copy the value
 * and the other you just want to use it as is. If NULL is passed then
 * no attempt will be made to copy the values.
 */
void rbtree_alter_value_copy(RBTREE_T *tree, void* (*value_copy)(void*));

/*
 * rbtree_alter_value_free
 * Changes the function used to free the value. If NULL is passed then
 * no attempt will be made to free the values.
 */
void rbtree_alter_value_free(RBTREE_T *tree, void (*value_free)(void*));

/*
 * rbtree_size
 * number of nodes in a red-black tree
 */
int rbtree_size(RBTREE_T *tree);

/*
 * rbtree_lookup
 * lookup a node in the tree. If the node doesn't exist and create is true then a new
 * node will be created using the key. If created is non-null then it will be set to true
 * when a new node is created.
 */
RBNODE_T* rbtree_lookup(RBTREE_T *tree, void *key, BOOLEAN_T create, BOOLEAN_T *created);

/*
 * rbtree_get
 * Gets the value for a key. If the key doesn't exist then NULL is returned.
 */
void *rbtree_get(RBTREE_T *tree, void *key);

/*
 * rbnode_get
 * Gets the value from a node.
 */
void *rbnode_get(RBNODE_T *node);

/*
 * rbtree_set
 * Updates the value. If the new value equals the previous value by pointer comparison then
 * nothing is done. Otherwise, if there was a previous value it will be freed using value_free 
 * (if it was passed to the constructor) and the new value will be copied using value_copy
 * (again only if it was passed to the constructor).
 */
void rbtree_set(RBTREE_T *tree, RBNODE_T *node, void *value);

/*
 * rbtree_put
 * Puts a key, value combination into the tree. If the key already exists the value is set to
 * the passed value as described in rbtree_set. Returns the node in the tree that represents
 * the key, value combination.
 */
RBNODE_T* rbtree_put(RBTREE_T *tree, void *key, void *value);

/*
 * rbtree_delete
 * removes the passed node from the passed tree and optionally returns the removed content. 
 * If removed_key is non-null then the free_key (from the constructor) will not be run on the
 * key, but it will instead be returned in the removed_key pointer. If removed_value is 
 * non-null then free_value will not be run on the value, but it will instad be returned in
 * the removed_value pointer. Destroys the node.
 */
void rbtree_delete(RBTREE_T *tree, RBNODE_T *node, void **removed_key, void **removed_value);

/*
 * rbtree_remove
 * Removes the key (and any value it has) from the tree. The key and value stored by the tree
 * are freed by free_key and free_value if these functions were set in the constructor.
 * Returns true if the key existed.
 */
BOOLEAN_T rbtree_remove(RBTREE_T *tree, void *key);

/*
 * rbtree_first
 * Returns the smallest node in the tree
 */
RBNODE_T *rbtree_first(RBTREE_T *tree);

/*
 * rbtree_last
 * Returns the largest node in the tree
 */
RBNODE_T *rbtree_last(RBTREE_T *tree);

/*
 * rbtree_next
 * Returns the next larger node in the tree or null if node is the largest
 */
RBNODE_T *rbtree_next(RBNODE_T *node);

/*
 * rbtree_prev
 * Returns the previous smaller node in the tree or null if node is the smallest
 */
RBNODE_T *rbtree_prev(RBNODE_T *node);

/*
 * rbtree_key
 * Returns the key. This key is used internally by the red black tree so do not
 * modify it as it will cause undefined results.
 */
void *rbtree_key(RBNODE_T *node);

/*
 * rbtree_value
 * Returns the value.
 */
void *rbtree_value(RBNODE_T *node);

/*
 * rbtree_strcmp
 * Utility function for using the red-black tree with normal strings.
 * Returns the result of comparing two strings using the strcmp function.
 */
int rbtree_strcmp(void *p1, void *p2);

/*
 * rbtree_strcpy
 * Utility function for using the red-black tree with normal strings.
 * Note that the free function can be used as the counterpart.
 * Returns a copy of the passed string.
 */
void* rbtree_strcpy(void *p);

/*
 * rbtree_intcmp
 * Utility function for using the red-black tree with pointers to integers.
 * Returns -1, 0 or 1 respectively if ints p1 < p2, p1 == p2 or p1 > p2.
 */
int rbtree_intcmp(void *p1, void *p2);

/*
 * rbtree_intcpy
 * Utility function for using the red-black tree with pointers to integers.
 * Note that the free function can be used as the counterpart.
 * Returns a malloc'ed copy of the passed int.
 */
void* rbtree_intcpy(void *p);


#endif

