/**************************************************************************
 * FILE: binary-search.h
 * AUTHOR: James Johnson 
 * CREATE DATE: 24-August-2009 
 * PROJECT: shared
 * COPYRIGHT: 2009, UQ 
 * VERSION: $Revision: 1.0 $
 * DESCRIPTION: Finds the index of an item, or place to insert in a sorted list.
 * Implements binary search in way that is slightly more useful than bsearch.
 **************************************************************************/


#ifndef BINARY_SEARCH_H
#define BINARY_SEARCH_H

#include <sys/types.h>

/*
 * binary_search
 * returns the index of the item equal to the key as compared using the compar function. If an equal item does not exist
 * then return -(indx + 1) where indx is the index that the item should be inserted. If the list is not sorted then the result
 * is undefined though it will return something...
 * The compar function works exactly like in bsearch so look that up if you need more details.
 */
int binary_search(const void *key, const void *base, size_t nmemb, size_t size, int (*compar)(const void *, const void *));


#endif
