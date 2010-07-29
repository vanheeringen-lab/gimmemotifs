/**
 * @file heap.c
 *
 * This module the fixed-size heap datatype.
 *
 * $Id: heap.c 2070 2007-09-12 01:49:09Z eredhead $
 * $Log$
 * Revision 1.2.2.5  2006/03/10 06:58:00  twhitington
 * print_heaps and "eval_best_seed" command line switches added.
 *
 * Revision 1.2.2.4  2006/01/30 06:37:35  twhitington
 * Updated documentation of get_key parameter of create_heap, in heap.c.
 * Introduced myfree for str_seed in seed.c function free_seed().
 *
 * Revision 1.2.2.3  2006/01/16 08:31:33  tbailey
 * Add comments to help with branching search implementation.
 *
 * Revision 1.2.2.2  2006/01/15 23:59:24  twhitington
 * Added log line to heap.c. Also added branching factor parameter to init.c.
 *
*/

#include "macros.h"
#include "heap.h"
#include "hash_table.h"
#include <assert.h>

//
//	create_heap
//
HEAP *create_heap(
  int max_size,                         // max size of heap; can grow if max_size < 0
  int (*compare)(void *p1, void *p2), 	// function to compare two objects
  void *(*copy)(void *p),      		// function to copy an objects
  void (*destroy)(void *p),      	// function to destroy an object
  /* The function get_key must return an ascii character string. Note that this
   * module (heap.c) does not perform an memory managment on that string
   * pointer; at no point in this module is the memory associated with that
   * string freed.
   */
  char *(*get_key)(void *p),     	// function to get an object string key
  void *(*print)(FILE *f, void *p)       // function to print nodes
)
{
  int n;
  HEAP *heap = NULL;
  Resize(heap, 1, HEAP);
  heap->node_list = NULL;  	// a list containing nodes of the heap
  heap->max_size = max_size;  
  n = (max_size < 0) ? HEAP_CHUNK : max_size + 1; // heap size is predefined or can grow
  Resize(heap->node_list, n, void *);
  heap->cur_size = n; 
  heap->next_node = 1;   	// index of root node
  heap->compare = compare;  	// heap order function
  heap->copy = copy;    	// function to copy two ojects
  heap->destroy = destroy; 	// function to destroy ojects
  heap->get_key = get_key;
  heap->print = print;          // function to print nodes
#define HASH_SIZE 100
  if (heap->get_key != NULL) {
    heap->ht = hash_create(HASH_SIZE);
  }
  else {
    heap->ht = NULL;
  }
  return (heap);
} // create_heap

//
//	add_node_heap
//
/*
  Adds a new object to the heap.
  The heap is implemented using an array. Properties of the heap:
    The root node is at index position 1
    The root has the smallest score in the heap
    The score for a child node is >= the score of its parent
    The left child of a node at index p is at index position 2*p
    The right child of a node at index p is at index positon 2*p+1
  
  If the heap is not full the object is added to the next available 
  index in the node list. The new object is then compared with its 
  parent. If the new object is < the parent object/s the nodes are 
  swapped until the heap order is re-established. 

  If the heap is full the new object is compared to the root object. 
  If the new object is > the root object, the root is "bumped"
  from the heap and the new object replaces the root object.  
  The new object is then compared with its children nodes and is swapped
  with the smallest child node until the heap order is re-established.

  If the heap is not full when add_node_heap() is called, it returns NULL.
  Otherwise, if the new node is successfully added to the heap, it
  returns a pointer to the bumped node (the old root node); if the
  node cannot be added (it is smaller than the old root), it returns
  a pointer to the new node.

  If a get_key function was defined, no objects will be added
  to the heap that have duplicate keys.
*/
void *add_node_heap(
  HEAP *heap,           // heap
  void *node            // object to add to heap
)
{
  int i = heap->next_node; 	// index for the next node
  int max = heap->max_size; 	// node_list indexes from 0 to max_size
  void *parent;			// parent node objects
  int p_node; 	 		// index of parent
  void *left_node, *right_node;	// objects in the left and right child nodes
  int l_idx, r_idx;		// index of left and right child nodes
  int comp_n;			//
  void *bumped_node = NULL; 	//
  char *node_key = NULL;	// string used as key of node in hash table

  // get info for new node
  void *new = node;

  // Check if node with same key already is in heap if we would add this node.
  // Don't add node to heap if one does. 
  if (heap->ht && ((i<=max) || (heap->compare(node, heap->node_list[1]) > 0))) {
    node_key = heap->get_key(new);
    if (hash_lookup_str(node_key, heap->ht)) return new;
  }

  // add nodes until the heap is full
  if (i <= max) {				// heap NOT full yet, so add the node
    // add to heap
    heap->node_list[i] = new;
    heap->next_node++;
    if (heap->ht) hash_insert_str(node_key, heap->ht);	// log new node's key
    
    // move new node up while the parent node is greater than the new node
    while (
      p_node = i/2, 				// index of parent node
      parent = heap->node_list[p_node], 	// parent object
      i > 1 && heap->compare(parent, new) > 0
    ) {
      // parent > new node therefore swap
      heap->node_list[p_node] = new;
      heap->node_list[i] = parent;
      // update the index for the new node
      i = p_node;
    }
    return NULL; 				// heap not full yet

  } else if (heap->compare(node, heap->node_list[1]) > 0) { 
    // heap full and new node is larger than the root node

    // the new node becomes the root and the "percolates" down
    i = 1; 		// update the index of the new node
    bumped_node = heap->node_list[1];
    heap->node_list[1] = new;
    if (heap->ht) {
      hash_insert_str(node_key, heap->ht);			// log new node's key
      hash_remove_str(heap->get_key(bumped_node), heap->ht);	// remove bumped node's key
    }

    //
    // percolate nodes down to reestablish heap 
    //
    while (i <= max/2) {		// max_size/2 is the last parent node
      // get children nodes
      l_idx = 2*i;  			// index of left child
      r_idx = l_idx+1;  		// index of right child
      left_node = heap->node_list[l_idx];
      // it is possible that there is only a left node
      right_node = ((r_idx) > max) ? NULL: heap->node_list[r_idx]; 

      // get the smallest child node 
      comp_n = (right_node) ? heap->compare(left_node, right_node): -1;
      
      // swap with smallest child
      if (comp_n < 0) {  			// left node is the smallest
        if (heap->compare(new, left_node) > 0) {
          heap->node_list[l_idx] = new;  	// push new node down a level
          heap->node_list[i] = left_node;  	// becomes parent node
          i = l_idx; 				// update the index for new node
        } else {  				// leave the nodes where they are 
          break;
        }
      } else {  				// the right node is the smallest
        if (heap->compare(new, right_node) > 0) {
          heap->node_list[r_idx] = new;  	// push down new node
          heap->node_list[i] = right_node;  	// becomes parent node
          i = r_idx;
        } else { 				// leave nodes where they are
          break;
        }
      }

    } // end while

    return bumped_node; 			// RETURN the pointer to the removed node

  } else {				// heap full but new node smaller than root node
    // the new node is not added to the heap
    return new;  			// RETURN the pointer to the new node
  }
} // add_node_heap

//
//      pop_heap_root
//
//      This function removes and then returns the root node 
//      from a heap. Following the removal of the root node,
//      the heap-order properties are restored.
//
//      Returns the root node for non-empty heaps, returns
//      a NULL pointer otherwise.
//
void *pop_heap_root(HEAP *h)
{

  void *root_node, *tmp_root, *left_node, *right_node;
  int l_idx, r_idx, comp_n, i;

//  fprintf(stdout, "ENTERING pop_heap_root\n");

  // Return NULL if the heap is empty
  if (get_num_nodes(h) <= 0) {
   return NULL;
  } 

  // Remove the root node from the heap and hash table
  root_node = h->node_list[1];
  if (h->ht) {
    hash_remove_str(h->get_key(root_node), h->ht); 
  }

  // Restore the heap. Put the last node at the root, then 
  // restablish the heap order property.

  // Move the last node in the heap to the root
  h->node_list[1] = tmp_root = get_node(h, get_num_nodes(h));
  h->node_list[get_num_nodes(h)] = NULL;
  h->next_node--; // there is one less node in the heap

  // Return the root node if there is only one node in the heap
  if (get_num_nodes(h) == 0){
    return root_node;
  }

  // Restore the heap-order property by comparing the new root 
  // node with all it's child nodes.

  i = 1; // start at the new root node
  while (i <= get_num_nodes(h)/2) {   // get_num_nodes(h)/2 is the last parent node
    // get children nodes
    l_idx = 2*i;                      // index of left child
    r_idx = l_idx+1;                  // index of right child
    left_node = h->node_list[l_idx];
    // it is possible that there is only a left node
    right_node = ((r_idx) > get_num_nodes(h)) ? NULL: h->node_list[r_idx];

    // get the smallest child node 
    comp_n = (right_node) ? h->compare(left_node, right_node): -1;

    // swap with smallest child
    if (comp_n < 0) {                         // left node is the smallest
      if (h->compare(tmp_root, left_node) > 0) {
        h->node_list[l_idx] = tmp_root;       // push new node down a level
        h->node_list[i] = left_node;          // becomes parent node
        i = l_idx;                            // update the index for new node
      } else {                                // leave the nodes where they are 
        break;
      }
    } else {                                  // the right node is the smallest
      if (h->compare(tmp_root, right_node) > 0) {
        h->node_list[r_idx] = tmp_root;       // push down new node
        h->node_list[i] = right_node;         // becomes parent node
        i = r_idx;
      } else {                                // leave nodes where they are
        break;
      }
    }
  } // end while

//  fprintf(stdout, "EXITING pop_heap_root\n");

  // return the old root node
  return root_node;
} // pop_heap_root


//
// 	copy_heap
//
// 	Create a new heap with the same maxsize and compare functions as h.
//
// 	Fill the new heap with deep copies of the objects in h if a copy 
// 	function is given, otherwise return an empty heap.
//
HEAP *copy_heap(HEAP *h)
{
  int i;
  void *temp, *temp2;
  HEAP *new = create_heap(h->max_size, h->compare, h->copy, h->destroy, h->get_key, h->print); 

  if (h->copy != NULL) {  	// check if the heap has an object copy function
    // add deep copies of the objects to the new heap
    for (i=1; i<= get_num_nodes(h); i++){
      temp = get_node(h, i);
      temp2 = h->copy(temp);
      temp2 = add_node_heap(new, temp2);
    }
  } else {
    fprintf(stderr, "Heap cannot be copied\n");
    exit(1); 
  }
  return (new);
} // copy_heap


//
// 	destroy_heap
//
// 	Free the heap and all of the objects it contains.
//
void destroy_heap(
  HEAP *h
)
{
  int i;
  // Free each of the nodes in the heap:
  for (i = 1; i < h->next_node; i++){
    h->destroy(h->node_list[i]);
  } 
  myfree(h->node_list); // Free the node list itself.
  if (h->ht) {
    hash_destroy(h->ht); // Free the hash table.
  }
  free(h);
} // destroy_heap
  
//
// 	get_num_nodes
//
// 	Returns the number of nodes in the heap
//
int get_num_nodes(
  HEAP *h
)
{
  // nodes are found at index positions 1, ..., (next_node-1)
  return (h->next_node - 1);
} // get_num_nodes

//
// 	get_node
//
// 	Returns a pointer to the object at position i in the heap
// 	(where the root is at position 1 and nodes are labeled according 
// 	to a level order traversal).
//
void *get_node(
  HEAP *h,
  int i
)
{
  assert(get_num_nodes(h) > 0); // Heap must not be empty

  return (h->node_list[i]);
} // get_node

int get_max_size(
  HEAP *h
)
{
  return (h->max_size);
} // get_max_size

//
// get the node with the highest score in the heap
//
int get_best_node(
  HEAP *h
)
{
  int i;
  int best_node_index = 1;
  void *best_node = get_node(h, best_node_index);

  for (i=2; i<=get_num_nodes(h); i++) {
    if (h->compare(best_node, get_node(h, i)) < 0) {
      best_node_index = i;
      best_node = get_node(h, best_node_index);
    }
  }
  return (best_node_index);
} // get_best_node


/**
 * print_heap
 *
 * Print the contents of the heap to the specified file.
 *
 */
extern void print_heap(
  FILE *outfile,     ///< The file to print to
  HEAP *heap         ///< The heap being printed
)
{
  fprintf(outfile, "##################################################\n"
                  "                       HEAP                       \n\n");

  if (heap->print != NULL) { // Check if the heap has an object print function
    // Print every node in the heap:
    int i;
    for (i=1; i<= get_num_nodes(heap); i++){
      void *curr_node = get_node(heap, i);
      fprintf(outfile, "NODE %i:\n", i);
      heap->print(outfile, curr_node);
    }
  } else {
    fprintf(outfile, "Heap cannot be printed.\n");
    exit(1); 
  }

  fprintf(outfile, "\n                   END HEAP                     \n"
                  "##################################################\n");
}
/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

