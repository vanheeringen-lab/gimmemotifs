/*---------------------------------------------------------------*
 * PSSM Searcher                                                 *
 *                                                               *
 * By: Jes Frellsen, Ida Moltke and Martin Thiim                 *
 * University of Copenhagen, spring 2006                         *
 *---------------------------------------------------------------*/
#include <stdlib.h>

#ifdef MEMCHECK
#ifndef MEMCHECK_h
#define MEMCHECK_h

#include <stdio.h>
#include <string.h>
#include <math.h>

static void *MEMmalloc(size_t size) {
  void *ptr = malloc(size);
  printf("\nMalloc:\t(%p)\n", ptr);
  return ptr;
}

static void MEMfree(void *ptr) {
  printf("\nFree:\t(%p)\n", ptr);
  free(ptr);
}

#define malloc(size) MEMmalloc(size)
#define free(ptr)    MEMfree(ptr)

#endif
#endif
