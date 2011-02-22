/* A heap is a partially sorted array, where, in one-based coordinates:
   array[i] >= array[i * 2] and
   array[i] >= array[i * 2 + 1] */
#ifndef HEAP_H
#define HEAP_H

#include <stddef.h>

#define PUSH_HEAP(base, n, cmp) push_heap(base, n, sizeof *(base), cmp)
#define POP_HEAP(base, n, cmp) pop_heap(base, n, sizeof *(base), cmp)

/* cmp should return negative if the 1st argument is less than the
   2nd, zero if equal, and positive if greater */

/* Add base[n-1] into the heap, by re-ordering base[0] ... base[n-1] */
/* Assumes base[0] ... base[n-2] is initially a heap */
void push_heap(void *base, size_t n, size_t size,
               int (*cmp)(const void *, const void *));

/* Remove base[0] from the heap */
/* i.e. move it to base[n-1], leaving a heap in base[0] ... base[n-2] */
/* Assumes base[0] ... base[n-1] is initially a heap */
void pop_heap(void *base, size_t n, size_t size,
	      int (*cmp)(const void *, const void *));

#endif
