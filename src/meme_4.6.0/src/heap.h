/* heap.h */

#ifndef HEAP_H
#define HEAP_H

#include "hash_table.h"
#ifdef DMALLOC
#include "dmalloc.h"
#endif

#define HEAP_CHUNK 100

// Heap 
typedef struct heap HEAP;
struct heap {
  int max_size;                         // max size of heap; no max size if < 0
  int cur_size;                         // current size of malloced array
  int height;                           // the height of the tree 
  int next_node;                        // index in node list of the 
                                        // next available node slot
  void **node_list;                     // array of heap node pointers
  int (*compare)(void *p1, void *p2);   // function to compare two objects:
                                        // returns:     1  if obj1 > obj2
                                        //              -1 if obj1 < obj2
                                        //              0  if obj1 = obj2
  void *(*copy)(void *p);               // function to deep copy an object
  void (*destroy)(void *p);             // function to destroy an object 
  char *(*get_key)(void *p);            // function to get object key
  void *(*print)(FILE *f, void *p);     // function to print nodes
  HASH_TABLE ht;                        // hash table for object keys
};

HEAP *create_heap(
  int max_size,                         // max size of heap; no max size if < 0
  int (*compare)(void *p1, void *p2),   // function to compare two objects
  void *(*copy)(void *p),               // function to deep copy an object
  void (*destroy)(void *p),             // function to destroy an object
  char *(*get_key)(void *p),            // function to get an object string key 
  void *(*print)(FILE *f, void *p)      // function to print nodes
);

void *add_node_heap(
  HEAP *heap,           // heap
  void *node            // node to add to heap
);

void *pop_heap_root(
  HEAP *h		// heap
);

HEAP *copy_heap(
  HEAP *h                 // the heap to be copied
);

void destroy_heap(
  HEAP *h               // the heap to be destroy
);

int get_num_nodes(
  HEAP *h               // heap
);

void *get_node(
  HEAP *h,              // heap containing the node of interest
  int i                 // position of the node in the heap
);

int get_max_size(
  HEAP *h
);

int get_best_node(      // get the index of the node with the highest score
  HEAP *h
);

extern void print_heap(
  FILE *outfile,     ///< The file to print to
  HEAP *heap         ///< The heap being printed
);
#endif
