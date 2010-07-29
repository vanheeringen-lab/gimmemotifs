/*************************************************************************
 * FILE: tree.h
 * AUTHOR: William Noble Grundy
 * CREATE DATE: 3-9-98
 * PROJECT: PHYLO
 * COPYRIGHT: 1998-1999, Regents of the University of California
 * DESCRIPTION: Data structure for manipulating phylogenetic trees.
 *************************************************************************/
#ifndef TREE_H
#define TREE_H

#include "rdb-matrix.h"
#include "string-list.h"
#include "utils.h"

#define MAX_LABEL 100   /* Length of the longest allowed label or key. */
#define MAX_DEGREE 100  /* Maximum number of children from one node. */

typedef struct tree_t TREE_T;

/*************************************************************************
 * Allocate memory for a new tree.
 *************************************************************************/
TREE_T* allocate_tree();

/*************************************************************************
 * Is the given tree a leaf node?
 *************************************************************************/
BOOLEAN_T is_leaf
  (TREE_T * const a_tree);

/*************************************************************************
 * What is the label at the root of the given tree?
 *************************************************************************/
char * get_label
  (TREE_T * const a_tree);

/*************************************************************************
 * Is there a label at the root of the given tree?
 *************************************************************************/
BOOLEAN_T has_label
  (TREE_T * const a_tree);

/*************************************************************************
 * What is the key at the root of the given tree?
 *************************************************************************/
char * get_key
  (TREE_T * const a_tree);

/*************************************************************************
 * Get or set the branch length at the root of the given tree?
 *************************************************************************/
float get_length
  (TREE_T * const a_tree);

void set_length
  (const float length,
   TREE_T *    a_tree);

/*************************************************************************
 * Get the number of children at the root of a tree.
 *************************************************************************/
int get_num_children
  (TREE_T * const a_tree);

/*************************************************************************
 * Get the total number of leaf descendants under a given root node.
 *************************************************************************/
int get_num_descendants
  (TREE_T * const a_tree, BOOLEAN_T force);

/*************************************************************************
 * Get the total number of edges under in a tree
 * by recursively summing all of the children.
 *************************************************************************/
int get_num_edges
  (TREE_T * const a_tree);

/*************************************************************************
 * Compute the maximum depth of a tree.
 *************************************************************************/
int compute_depth
  (TREE_T * const a_tree);

/*************************************************************************
 * Calculate pairwise distances between distinct leaves in a tree.
 * Returns a string labeled matrix. Caller is responsible for freeing matrix.
 *************************************************************************/
RDB_MATRIX_T* get_leaf_distances(TREE_T*  const the_tree);

/*************************************************************************
 * Add a given child to a tree.
 *************************************************************************/
void add_child
  (const VERBOSE_T verbosity,
   TREE_T * const  a_child,
   TREE_T *        a_tree);

/*************************************************************************
 * Retrieve the nth child of a given tree.
 *************************************************************************/
TREE_T * get_nth_child
  (const int      n,
   TREE_T * const a_tree);

/*************************************************************************
 * Remove the nth child of a given tree.
 *************************************************************************/
void remove_nth_child
  (const VERBOSE_T verbosity,
   const BOOLEAN free_child,    // free child node?
   const BOOLEAN_T free_children, /* Boolean: Free children as well? */
   const int       n,
   TREE_T *        a_tree);

/*************************************************************************
 * Read a tree from in New Hampshire parentheses format from an open
 * FILE structure.
 *************************************************************************/
void read_tree
  (FILE *    tree_file,
   TREE_T ** a_tree);

/*************************************************************************
 * Create a tree from the named file which should be in New Hampshire format.
 *************************************************************************/
TREE_T* read_tree_from_file
  (char* filename);

/*************************************************************************
 * Copy a tree. Assume memory is already allocated.
 *************************************************************************/
void copy_tree
  (TREE_T * source_tree,
   TREE_T * dest_tree);
  
/*************************************************************************
 * Write a tree to a file in New Hampshire parentheses format.
 *************************************************************************/
void write_tree
  (TREE_T * a_tree,
   FILE *   outfile);

/***************************************************************************
 * Find the length of the longest label in a tree.
 ***************************************************************************/
int longest_label
  (TREE_T * the_tree);

/***************************************************************************
 * Count the number of leaves in a tree.
 ***************************************************************************/
int leaf_count
  (TREE_T * the_tree);

/***************************************************************************
 * Count the number of labeled nodes in a tree.
 * TLB: This is needed when we are putting the target genome at the
 * root of the tree.  In that case, the number of labels gives us
 * the number of sequences in the alignment, and the leaf_count is one short.
 ***************************************************************************/
int label_count
  (TREE_T * the_tree);

/***************************************************************************
 * Accumulate into an array of strings the list of all leaf labels
 * that appear at or below a given tree.
 ***************************************************************************/
STRING_LIST_T * make_leaf_list
 (TREE_T *         the_tree);

/*************************************************************************
 * Remove leaves from a given tree.
 *************************************************************************/
TREE_T* trim_tree
(BOOLEAN_T      allocate,
 TREE_T*        the_tree,
 STRING_LIST_T* the_leaves);

/*************************************************************************
 * Convert a tree into a star tree with the given label
 * as the center.  Each leaf (except the given label)
 * becomes a leaf in the new tree, with distance equal
 * to total branch length of the old tree divided by N-1.
 *
 * This function modifies the tree in place and returns
 * the new root.  The old tree is lost.
 *************************************************************************/
TREE_T* convert_to_uniform_star_tree(
  TREE_T*       the_tree,
  char*         the_label
);

/*************************************************************************
 * Re-root a tree at the leaf node with the given label.
 * The old root becomes an internal node.  The new root
 * keeps its label (and associated sequence).
 *
 * This function modifies the tree in place and returns
 * the new root.
 *************************************************************************/
TREE_T* reroot_tree(
  TREE_T*       the_tree,
  char*         the_label,
  TREE_T*       the_parent              // Call with NULL
);

/*************************************************************************
  What is the total branch length of a tree?
*************************************************************************/
float get_total_length
  (TREE_T * const a_tree);


/*************************************************************************
  What is the branch length of the sub-tree specified by a list of labels?
*************************************************************************/


float get_subtree_length (
  TREE_T*  const the_tree,
  STRING_LIST_T* leaf_labels
);

/*************************************************************************
 * Free dynamic memory used by a tree.
 *************************************************************************/
void free_tree
  (const BOOLEAN_T free_children, /* Free children as well? */
   TREE_T * a_tree);

#endif
