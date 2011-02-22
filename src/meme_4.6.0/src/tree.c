/*************************************************************************
 * FILE: tree.c
 * AUTHOR: William Stafford Noble, Charles E. Grant
 * CREATE DATE: 3-9-98
 * PROJECT: EVOMCAST
 * COPYRIGHT: 1998-2004, W. S. N., UW
 *************************************************************************/
#include "tree.h"
#include "utils.h"
#include "string-list.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

struct tree_t {
  char      key[MAX_LABEL+1];   /* The key of a leaf node is just the
				   label, and for an interior node is the 
				   lexicographically lowest key of the 
				   children. */
  char      label[MAX_LABEL+1]; /* The label of this leaf. Interior nodes
				   are unlabelled. 
				   TLB; The exception is when the target 
			           genome is placed at the root of the tree.
			        */
  int       num_children;       /* How many children does this leaf have? */
  BOOLEAN_T has_length;         /* Does this node have a branch length? */
  float     length;             /* Length of branch leading to this node. */
  BOOLEAN_T has_descendants;    /* Has the number of descendents been 
				                           calculated? */
  int       num_descendants;    /* Number of nodes below this in the tree. */
  struct tree_t* children[MAX_DEGREE];    /* Children of the current node. */
};

static void convert_to_uniform_star_tree_helper(TREE_T* the_tree, TREE_T* new_tree, double length);

/*************************************************************************
 * Allocate memory for a new tree.
 *************************************************************************/
TREE_T* allocate_tree()
{
  TREE_T* a_tree = (TREE_T *) mm_malloc(sizeof(TREE_T));
  a_tree->label[0] = '\0';
  a_tree->key[0] = '\0';
  a_tree->num_children = 0;
  a_tree->has_length = FALSE;
  a_tree->length = 0.0;
  a_tree->has_descendants = FALSE;
  a_tree->num_descendants = 0;
  int i;
  for(i = 0; i < MAX_DEGREE; i++) {
    a_tree->children[i] = NULL;
  }

  return a_tree;
}

/*************************************************************************
 * Is the given tree a null tree?
 *************************************************************************/
static void check_null_tree
  (TREE_T * const a_tree)
{
  if (a_tree == NULL) {
    die("Attempted to access null tree.\n");
  }
}


/*************************************************************************
 * Is the given tree a leaf node?
 *************************************************************************/
BOOLEAN_T is_leaf
  (TREE_T * const a_tree)
{
  check_null_tree(a_tree);

  if (a_tree->num_children == 0) {
    return(TRUE);
  }
  return(FALSE);
}

/*************************************************************************
 * What is the label at the root of the given tree?
 *************************************************************************/
char * get_label
  (TREE_T * const a_tree)
{
  check_null_tree(a_tree);

  return(a_tree->label);
}

/*************************************************************************
 * Is there a label at the root of the given tree?
 *************************************************************************/
BOOLEAN_T has_label
  (TREE_T * const a_tree)
{
  check_null_tree(a_tree);

  return((a_tree->label[0] != '\0'));
}

/*************************************************************************
 * What is the key at the root of the given tree?
 *************************************************************************/
char * get_key
  (TREE_T * const a_tree)
{
  check_null_tree(a_tree);

  return(a_tree->key);
}

/*************************************************************************
 * What is the branch length at the root of the given tree?
 *************************************************************************/
float get_length
  (TREE_T * const a_tree)
{
  check_null_tree(a_tree);

  return(a_tree->length);
}

void set_length
  (const float length,
   TREE_T *    a_tree)
{
  check_null_tree(a_tree);

  if (length < 0.0) {
    die("Bad branch length (%g) for tree %s.\n",
	length, get_key(a_tree));
  }

  if (length > 0) {
    a_tree->has_length = TRUE;
  } else {
    a_tree->has_length = FALSE;
  }
  a_tree->length = length;
}

/*************************************************************************
 * Get the number of children at the root of a tree.
 *************************************************************************/
int get_num_children
  (TREE_T * const a_tree)
{
  check_null_tree(a_tree);
  return(a_tree->num_children);
}


/**************************************************************************
 * Compute and store the number of descendants of each node in a tree.
 **************************************************************************/
static void compute_descendants
  (TREE_T * a_tree, BOOLEAN_T force) 
{
  int i_child;
  int num_descendants;

  /* Don't do anything if we've already computed the descendants. */
  if (!force && a_tree->has_descendants) {
    return;
  }

  /* Base case: A leaf has only itself as a descendant. */
  if (is_leaf(a_tree)) {
    num_descendants = 1;
    a_tree->has_descendants = FALSE;
  }
  else {
    /* Recursive case: Descendants of each child. */
    num_descendants = 0;
    a_tree->has_descendants = TRUE;
    for (i_child = 0; i_child < a_tree->num_children; i_child++) {
      compute_descendants(get_nth_child(i_child, a_tree), force);
      num_descendants += get_num_descendants(get_nth_child(i_child, a_tree), force);
    }
  }
  a_tree->num_descendants = num_descendants;
}

/*************************************************************************
 * Get the total number of leaf descendants under a given root node.
 *************************************************************************/
int get_num_descendants
  (TREE_T * const a_tree, BOOLEAN_T force)
{
  check_null_tree(a_tree);
  compute_descendants(a_tree, force);
  return(a_tree->num_descendants);
}

/*************************************************************************
 * What is the total branch length of the given tree?
 *************************************************************************/
float get_total_length
  (TREE_T * const a_tree)
{
  check_null_tree(a_tree);
  int num_children = get_num_children(a_tree);
  double length = get_length(a_tree);
  int c = 0;

  if (!is_leaf(a_tree)) {
    for (c = 0; c < num_children; c++) {
      TREE_T* child = get_nth_child(c, a_tree);
      length += get_total_length(child);
    }
  }

  return length;
} // get_total_length

/*************************************************************************
  What is the branch length of the sub-tree specified by a list of labels?
*************************************************************************/
float get_subtree_length (
  TREE_T*  const the_tree,
  STRING_LIST_T* leaf_labels
) {

	float length = 0.0;
	float tempLength = 0.0;
	int i = 0;
  	int j = 0;

	// For each child, recursively calculate distance to leaves.
	// But only includes leaves with labels in the list.

	for (i = 0; i < the_tree->num_children; i++) {
	   tempLength = 0.0;
	   if (is_leaf(the_tree->children[i])) {
      		int leaf_index = get_index_in_string_list(the_tree->children[i]->label, leaf_labels);
      		if (leaf_index <  0) {
		   // This child is not in the list so we ignore it   
      		} else {
      		   length += the_tree->children[i]->length;
    		}
	   } else {
		// The child node is not a leaf so we recurse 
      		tempLength = get_subtree_length(the_tree->children[i], leaf_labels);
		// We check to see if anything came back from the recursion
		// before adding the length of the current child to total
		if (tempLength > 0.0) {
			length += tempLength;
			length += the_tree->children[i]->length;
		}
     	   }
	}

	return length;
}





/*************************************************************************
 * Get the total number of edges under in a tree
 * by recursively summing all of the children.
 *************************************************************************/
int get_num_edges
  (TREE_T * const a_tree)
{
  assert(a_tree != NULL);
  TREE_T* child = NULL;
  int num_children = get_num_children(a_tree);
  int num_edges = num_children;
  int c = 0;
  for (c = 0; c < num_children; c++) {
    child = get_nth_child(c, a_tree);
    num_edges += get_num_edges(child);
  }
  return num_edges;
}

/*************************************************************************
 * Compute the maximum depth of a tree.
 *************************************************************************/
int compute_depth
  (TREE_T * const a_tree)
{
  /* Base case: leaf. */
  if (is_leaf(a_tree)) {
    return(1);
  }

  /* Recursive case: internal node. */
  else {
    int max_depth = 0;
    int num_children = get_num_children(a_tree);
    int i_child;
    for (i_child = 0; i_child < num_children; i_child++) {
      int this_depth = compute_depth(get_nth_child(i_child, a_tree));
      if (this_depth > max_depth) {
	max_depth = this_depth;
      }
    }
    return(max_depth + 1);
  }
  /* Unreachable. */
  abort();
  return(0);
}

/*************************************************************************
 * Add a given child to a tree.
 *************************************************************************/
void add_child
  (const VERBOSE_T verbosity,
   TREE_T * const  a_child,
   TREE_T *        a_tree)
{
  check_null_tree(a_tree);

  /* Make sure we don't put too many children in the tree. */
  if (get_num_children(a_tree) >= MAX_DEGREE) {
    die("Attempted to add %s to tree (%s) with maximum degree (%d).\n",
	get_key(a_child), get_key(a_tree), get_num_children(a_tree));
  }

  if (verbosity >= HIGHER_VERBOSE) {
      fprintf(stderr, "Adding child %s to tree %s.\n", get_key(a_child),
	      get_key(a_tree));
  }

  /* Add the child to the tree. */
  a_tree->children[get_num_children(a_tree)] = a_child;
  (a_tree->num_children)++;

  /* Update the number of descendants, if necessary. */
  if (a_tree->has_descendants) {
    compute_descendants(a_child, FALSE);
    a_tree->num_descendants += a_child->num_descendants;
  }
  
}


/*************************************************************************
 * Retrieve the nth child of a given tree.
 *************************************************************************/
TREE_T * get_nth_child
  (const int      n,
   TREE_T * const a_tree)
{
  check_null_tree(a_tree);

  /* Make sure the child exists. */
  if (n >= get_num_children(a_tree)) {
    die("Attempted to retrieve child %d from a tree with %d children.\n", 
	n, get_num_children(a_tree));
  }

  /* Return the requested child. */
  return(a_tree->children[n]);
}

/*************************************************************************
 * Remove the nth child of a given tree.
 *
 * Does not free dynamic memory. 
 *************************************************************************/
void remove_nth_child
  (const VERBOSE_T verbosity,
   const BOOLEAN free_child,		// free child node?
   const BOOLEAN_T free_children, 	/* Boolean: Free children as well? */
   const int       n,
   TREE_T *        a_tree)
{
  int i_child;

  check_null_tree(a_tree);

  /* Make sure the child exists. */
  if (n >= get_num_children(a_tree)) {
    die("Attempted to remove child %d from a tree with %d children.\n", 
	n, get_num_children(a_tree));
  }

  if (verbosity >= HIGHER_VERBOSE) {
      fprintf(stderr, "Removing child %s from tree %s.\n", 
	      get_key(a_tree->children[n]), get_key(a_tree));
  }

  /* Free the child and its children. */
  if (free_child) free_tree(free_children, a_tree->children[n]);

  /* Move the remaining children left one slot. */
  for (i_child = n; i_child < get_num_children(a_tree) - 1; i_child++) {
    a_tree->children[i_child] = a_tree->children[i_child+1];
  }

  /* Decrement the number of children. */
  (a_tree->num_children)--;

} // remove_nth_child


/*************************************************************************
 * Read the branch length of a node.
 *************************************************************************/
static void process_length
  (char *   one_char,
   FILE *   tree_file,
   TREE_T * a_tree)
{
  long int  digit;
  long int  ordzero;
  double    value;
  double    divisor;
  BOOLEAN_T pointread; // Have we read a decimal point?
  BOOLEAN_T minusread; // Have we read a minus sign?

  ordzero = '0';
  pointread = FALSE;
  minusread = FALSE;
  value = 0.0;
  divisor = 1.0;
  get_non_blank(tree_file, one_char);
  digit = *one_char - ordzero;

  while (((unsigned long)digit <= 9) | (*one_char == '.') ||
         (*one_char == '-')) {
    if (*one_char == '.') {
      pointread = TRUE;
    } else if (*one_char == '-') {
      minusread = TRUE;
    } else {
      value = value * 10.0 + digit;
      if (pointread) {
	divisor *= 10.0;
      }
    }
    get_non_blank(tree_file, one_char);
    digit = *one_char - ordzero;
  }

  if (minusread) {
    a_tree->length = 0.0;
  } else {
    a_tree->length = value / divisor;
  }
}

/*************************************************************************
 * Recursively read a tree from a file.
 *************************************************************************/
static void read_node
(char *    one_char,
 FILE *    tree_file,
 TREE_T ** a_tree)
{
  BOOLEAN_T last_child;
  int i_child;
  int i_char;

  // Allocate memory for the new node.
  *a_tree = allocate_tree();

  // Are we starting a new subtree?
  if (*one_char == '(') {

    // Read in the children, one by one.
    last_child = FALSE;
    i_child = 0;
    while (!last_child) {
      get_non_blank(tree_file, one_char);

      // Read one child.
      read_node(one_char, tree_file, &(*a_tree)->children[i_child]);

      // Increment the number of children.
      i_child++;

      // Are we done with this set of children?
      if (*one_char == ')') {
        last_child = TRUE;

        // Allow labels on internal nodes.
        get_non_blank(tree_file, one_char);
	i_char = 0;
        while (*one_char != ':' && 
	       *one_char != ',' && 
	       *one_char != ')' && 
	       *one_char != '[' &&
	       *one_char != ';') {

	  // Replace special characters.
	  if ((*one_char & 255) == 255)
	    *one_char = '\'';
	  if ((*one_char & 255) > 175)
	    *one_char -= 48;
	  if ((*one_char & (~127)) != 0)
	    *one_char -= 64;

	  // Only store the first MAX_LABEL characters.
	  if (i_char < MAX_LABEL)
	    (*a_tree)->label[i_char] = *one_char;

	  // Skip line breaks.
	  if (is_eoln(tree_file)) {
	    fscanf(tree_file, "%*[^\n]");
	    getc(tree_file);
	  }
	  *one_char = getc(tree_file);
	  i_char++;
        } 
      }
    }
    (*a_tree)->num_children = i_child;
  }

  // Otherwise, this is a leaf node.
  else {
    (*a_tree)->num_children = 0;
    i_char = 0;
    do {

      // Replace special characters.
      if ((*one_char & 255) == 255)
        *one_char = '\'';
      if ((*one_char & 255) > 175)
        *one_char -= 48;
      if ((*one_char & (~127)) != 0)
        *one_char -= 64;

      // Only store the first MAX_LABEL characters.
      if (i_char < MAX_LABEL)
        (*a_tree)->label[i_char] = *one_char;

      // Skip line breaks.
      if (is_eoln(tree_file)) {
        fscanf(tree_file, "%*[^\n]");
        getc(tree_file);
      }
      *one_char = getc(tree_file);
      i_char++;
    } while (*one_char != ':' && *one_char != ',' && *one_char != ')');

    if (i_char > MAX_LABEL)
      i_char = MAX_LABEL;
    (*a_tree)->label[i_char] = '\0';
  }

  // Get the branch length.
  if (*one_char == ':') {
    process_length(one_char, tree_file, *a_tree);
    (*a_tree)->has_length = TRUE;
  }
}

/*************************************************************************
 * Read a tree from New Hampshire parentheses format from an open
 * FILE structure.
 *************************************************************************/
void read_tree
  (FILE *    tree_file,
   TREE_T ** a_tree)
{
  char one_char;

  get_non_blank(tree_file, &one_char);
  read_node(&one_char, tree_file, a_tree);
  fscanf(tree_file, "%*[^\n]");
  getc(tree_file);
}

TREE_T* read_tree_from_file
  (char* filename) 
{

    TREE_T* tree;
    FILE* tree_file = NULL;

    if (open_file(filename, "r", 1, "tree", "tree", &tree_file) == 0) {
      die("Couldn't open the file %s.\n", filename);
    }
    read_tree(tree_file, &tree);
    (void) fclose(tree_file);
    if (verbosity >= HIGH_VERBOSE) {
      fprintf(stderr, "Read tree: ");
      write_tree(tree, stderr);
    }

    return tree;
}

/*************************************************************************
 * Deep copy a tree. Assume memory is already allocated.
 *************************************************************************/
void copy_tree
  (TREE_T * source_tree,
   TREE_T * dest_tree)
{
  memcpy(dest_tree, source_tree, sizeof(TREE_T));
  int i;
  for (i = 0; i < dest_tree->num_children; i++) {
    dest_tree->children[i] = allocate_tree();
    copy_tree(source_tree->children[i], dest_tree->children[i]);
  }
}

/*************************************************************************
 * Write a tree to a file in New Hampshire parentheses format.
 *************************************************************************/

/* This auxiliary function is only necessary because the tree must end
   with a semi-colon. */
static void write_tree_aux
  (TREE_T * a_tree,
   FILE *   outfile)
{
  int i_child;

  /* If it's a leaf, print it. */
  if (is_leaf(a_tree)) {
    fprintf(outfile, "%s", a_tree->label);
    if (a_tree->has_length) {
      fprintf(outfile, ":%7.5f", a_tree->length);
    }
  }

  /* Otherwise, parenthesize and print the children with commas between. */
  else {
    fprintf(outfile, "(");
    for (i_child = 0; i_child < a_tree->num_children - 1; i_child++) {
      write_tree_aux(a_tree->children[i_child], outfile);
      fprintf(outfile, ",");
    }
    write_tree_aux(a_tree->children[i_child], outfile);
    fprintf(outfile, ")");
    if (a_tree->has_length) {
      fprintf(outfile, ":%7.5f", a_tree->length);
    }
  }
}

void write_tree
  (TREE_T * a_tree,
   FILE *   outfile)
{
  /* If the tree has only one node, print it in legal Newick format */
  if (is_leaf(a_tree)) {
    fprintf(outfile, "(%s);\n", a_tree->label);
  } else {
    write_tree_aux(a_tree, outfile);
    // Print label of root if there is one.
    fprintf(outfile, "%s;\n", a_tree->label); 
  }
}

/*************************************************************************
 * Compare two trees according to the lexicographical ordering of the
 * keys of their roots.
 *************************************************************************/
#ifndef TREE_DEBUG
#define TREE_DEBUG 0
#endif
static int key_compare
  (const void * t1,
   const void * t2)
{
  const TREE_T ** tree1 = (const TREE_T **) t1;
  const TREE_T ** tree2 = (const TREE_T **) t2;

  DEBUG_CODE(TREE_DEBUG,
	     printf("tree1->key=%s tree2->key=%s\n", (*tree1)->key,
		    (*tree2)->key);
	     );
  
  return(strcmp((*tree1)->key, (*tree2)->key));
}

/***************************************************************************
 * Find the length of the longest label in a tree.
 ***************************************************************************/
int longest_label
  (TREE_T * the_tree)
{
  int i_child;     /* Index of the current child of this node. */
  int this_length; /* Length of the label of this leaf. */
  int max_length;  /* Maximum leaf label length. */

  /* Base case: Just return the length of this label. */
  if (is_leaf(the_tree)) {
    max_length = strlen(the_tree->label);
  }

  /* Recursion: Find the longest label of the child subtrees. */
  else {
    max_length = 0;
    for (i_child = 0; i_child < the_tree->num_children; i_child++) {
      this_length = longest_label(the_tree->children[i_child]);
      if (this_length > max_length) {
	max_length = this_length;
      }
    }
  }
  return(max_length);
} // leaf_count

/***************************************************************************
 * Count the number of labeled nodes in a tree.
 * TLB: This is needed when we are putting the target genome at the
 * root of the tree.  In that case, the number of labels gives us
 * the number of sequences in the alignment, and the leaf_count is one short.
 ***************************************************************************/
int label_count
  (TREE_T * the_tree)
{
  int i_child;
  int num_labels;

  // See if current node has a label.
  num_labels = has_label(the_tree) ? 1 : 0;

  /* Base case: A single leaf. */
  if (is_leaf(the_tree)) {
    // NOOP
  }

  /* Recursion: Add up the labels of each child. */
  else {
    for (i_child = 0; i_child < the_tree->num_children; i_child++) {
      num_labels += label_count(the_tree->children[i_child]);
    }
  }
  return(num_labels);
} // label_count

/***************************************************************************
 * Count the number of leaves in a tree.
 ***************************************************************************/
int leaf_count
  (TREE_T * the_tree)
{
  int i_child;
  int num_leaves;

  /* Base case: A single leaf. */
  if (is_leaf(the_tree)) {
    num_leaves = 1;
  }

  /* Recursion: Add up the leaves of each child. */
  else {
    num_leaves = 0;
    for (i_child = 0; i_child < the_tree->num_children; i_child++) {
      num_leaves += leaf_count(the_tree->children[i_child]);
    }
  }
  return(num_leaves);
}

/***************************************************************************
 * Recursive version of 'make_leaf_list', this function does not
 * allocate memory.
 ***************************************************************************/
void get_leaf_names
  (TREE_T * const  the_tree,   /* The array of trees to be searched. */
   STRING_LIST_T * leaf_list)  /* The leaf labels. */
{
  int i_child;

  /* Base case: A single leaf. */
  if (is_leaf(the_tree)) {
    add_string(the_tree->label, leaf_list);
  }

  /* Recursion. */
  else {
    for (i_child = 0; i_child < the_tree->num_children; i_child++) {
      get_leaf_names(the_tree->children[i_child], leaf_list);
    }
  }
}

/***************************************************************************
 * Accumulate into an array of strings the list of all leaf labels
 * that appear at or below a given tree.
 ***************************************************************************/
STRING_LIST_T * make_leaf_list
 (TREE_T *         the_tree)
{
  STRING_LIST_T *  leaf_list;   /* The list of leaves. */

  /* Allocate the list of leaves. */
  leaf_list = new_string_list();

  /* Make a list of all the leaves in the tree. */
  get_leaf_names(the_tree, leaf_list);

  return(leaf_list);
}

/*************************************************************************
 * TLB;
 * Return a labeled tree node.
 *************************************************************************/
TREE_T* get_named_tree(
  TREE_T*       the_tree,
  char* 	the_label
)
{
  // Base case.
  if (is_leaf(the_tree)) {
    char* label = get_label(the_tree);
    if (!strcmp(label, the_label)) {
      return(the_tree);		// Found.
    } else {
      return(NULL);
    }
  } 

  // Recursive case.
  int num_children = get_num_children(the_tree);
  int i_child;
  for (i_child = 0; i_child < num_children; i_child++) {
    TREE_T* child = get_nth_child(i_child, the_tree);
    child = get_named_tree(child, the_label);
    if (child != NULL) return(child); // Found on this branch.
  }

  return(NULL);			// Not found on this branch.
} // get_named_tree

/*************************************************************************
 * TLB;
 * Convert a tree into a star tree with the given label
 * as the center.  Each leaf (except the given label) 
 * becomes a leaf in the new tree, with length equal
 * to total branch length of the old tree divided by N-1.
 *
 * This function modifies the tree in place and returns
 * the new root.  The old tree is lost.
 *************************************************************************/
TREE_T* convert_to_uniform_star_tree(
  TREE_T*       the_tree,
  char* 	the_label
)
{
  double d_tree = get_total_length(the_tree);
  int N = get_num_descendants(the_tree, FALSE);
  double d = d_tree/(N-1);		// new branch length
  // Find the node corresponding to the given label.
  TREE_T *new_tree = get_named_tree(the_tree, the_label);
  if (new_tree == NULL) die("Couldn't find label %s in tree.\n", the_label);
  // Set length to zero.
  set_length(0, new_tree);
  // Convert the tree.
  convert_to_uniform_star_tree_helper(the_tree, new_tree, d);
  // Make sure the tree is kosher.
  (void) get_num_descendants(the_tree, TRUE);

  return(new_tree);
} // convert_to_uniform_star_tree

/*************************************************************************
 * TLB;
 * Recursive helper routine for convert_to_uniform_star_tree
 *************************************************************************/
static void convert_to_uniform_star_tree_helper(
  TREE_T*       the_tree,
  TREE_T*       new_tree,
  double 	length
)
{
  if (the_tree == new_tree) return;
  
  if (is_leaf(the_tree)) {
    // add current leaf to new tree
    set_length(length, the_tree);
    add_child(verbosity, the_tree, new_tree);
  } else {
    int num_children = get_num_children(the_tree);
    int i_child;
    for (i_child = 0; i_child < num_children; i_child++) {
      TREE_T* child = get_nth_child(i_child, the_tree);
      convert_to_uniform_star_tree_helper(child, new_tree, length);
    }
    free_tree(FALSE, the_tree);		// free non-leaf
  }
  return;
} // convert_to_uniform_star_tree_helper

/*************************************************************************
 * TLB;
 * Re-root a tree at the leaf node with the given label.
 * The old root becomes an internal node.  The new root
 * keeps its label (and associated sequence).
 *
 * This function modifies the tree in place and returns
 * the new root.
 *
 * Note: this function appears to be pointless as it has 
 * no effect on the scoring function of Motiph.
 *************************************************************************/
TREE_T* reroot_tree(
  TREE_T*       the_tree,
  char* 	the_label,
  TREE_T*	the_parent		// Call with NULL
)
{
  // Base case.
  if (is_leaf(the_tree)) {
    char* label = get_label(the_tree);

    if (verbosity > NORMAL_VERBOSE) {
      fprintf(stderr, "Rerooting- at leaf (%s).\n", label);
    }

    if (!strcmp(label, the_label)) {
      // We have found the leaf we are looking for.
      // Make it the root.
      // Start by making its old parent its child.
      // Then return the non-null current node.
      add_child(verbosity, the_parent, the_tree);
      return(the_tree);
    } else {
      // This leaf is not the one we want.
      // Return NULL to signal this branch was a dead end.
      return(NULL);
    }
  } 

  // Recursive case.
  int num_children = get_num_children(the_tree);
  int i_child;
  for (i_child = 0; i_child < num_children; i_child++) {
    TREE_T* new_root = reroot_tree(
      get_nth_child(i_child, the_tree),
      the_label,
      the_tree
    );

    if (new_root != NULL) {
      // The desired leaf was found.  
      // Continue fixing the tree.
      // First, remove the current child.  
      if (verbosity > NORMAL_VERBOSE) {
	fprintf(stderr, "Rerooting- backing up from (%d).\n", i_child);
      }
      remove_nth_child(
        verbosity,
        FALSE,				// Don't free child node
        FALSE,				// Don't go recursive
        i_child,
        the_tree
      );
      // Then, make the current node's parent its child.
      // The new_root contains the branch-length for the current node.
      // Then, return new_root.
      // Old root has no parent!
      double old_length = get_length(the_tree); // save old length to pass back up
      double new_length = get_length(new_root);	// new length for node in new root
      set_length(new_length, the_tree);		// update current node's length
      if (the_parent != NULL) {
        // backing up towards old root
	add_child(verbosity, the_parent, the_tree);
        set_length(old_length, new_root);		// pass old length up to parent
      } else {
        // backed up to the old root
        new_root->has_length = FALSE;		// new root has no length
      }
      return(new_root);
    }

  } // loop over children
 
  return(NULL);

} // reroot_tree

/*************************************************************************
 * Helper function for trim_tree.
 *************************************************************************/
static TREE_T* trim_tree_recursive
(TREE_T*        the_tree,
 STRING_LIST_T* the_leaves)
{

  // Base case.
  if (is_leaf(the_tree)) {
    char* label = get_label(the_tree);

    // Should we prune this leaf?
    if (have_string(label, the_leaves)) {
      return(the_tree);
    } else {
      if (verbosity > NORMAL_VERBOSE) {
	fprintf(stderr, "Removing leaf (%s).\n", label);
      }
      free_tree(FALSE, the_tree);
      return(NULL);
    }
  } 

  // Recursive case.
  int num_children = get_num_children(the_tree);
  int i_child;
  int num_remaining_children = 0;
  for (i_child = 0; i_child < num_children; i_child++) {

    // Trim this child.
    TREE_T* this_child = trim_tree_recursive(get_nth_child(i_child, the_tree),
					     the_leaves);

    // If the child is non-NULL, keep it.
    if (this_child != NULL) {
      the_tree->children[num_remaining_children] = this_child;
      num_remaining_children++;
    } else {
      if (verbosity >= HIGHER_VERBOSE) {
	      fprintf(stderr, "Removed a child.\n");
      }
    }
  }

  // If we removed children, check for zero or one child case.
  if (num_remaining_children != num_children) {
    if (verbosity >= HIGHER_VERBOSE) {
      fprintf(stderr, "Changed from %d to %d children.\n", 
	      num_children, num_remaining_children);
    }

    // Update the number of children.
    the_tree->num_children = num_remaining_children;

    // Check for one remaining child.
    if (num_remaining_children == 1) {
      TREE_T* this_child = get_nth_child(0, the_tree);
      assert(this_child != NULL);

      if (verbosity >= HIGHER_VERBOSE) {
	      fprintf(stderr, "Collapsing a branch.\n");
      }

      // Compute the new branch length.
      float new_length = get_length(the_tree) + get_length(this_child);
      if (verbosity >= HIGHER_VERBOSE) {
	      fprintf(stderr, "New length: %g + %g = %g\n", 
		    get_length(the_tree), get_length(this_child), new_length);
      }

      // Copy the child into this location.
      copy_tree(this_child, the_tree);

      // Free the child.
      free_tree(TRUE, this_child);

      // Set the new branch length.
      set_length(new_length, the_tree);
    }

    // Check for no remaining children.
    else if (num_remaining_children == 0) {
      if (verbosity >= HIGHER_VERBOSE) {
	      fprintf(stderr, "No children.\n");
      }
      free_tree(TRUE, the_tree);
      return(NULL);
    }
  }

  return(the_tree);
}

/*************************************************************************
 * Remove leaves from a given tree.
 *************************************************************************/
TREE_T* trim_tree
(BOOLEAN_T      allocate,
 TREE_T*        the_tree,
 STRING_LIST_T* the_leaves)
{
  TREE_T* return_value;

  // If requested, allocate a copy of the tree.
  if (allocate) {
    return_value = allocate_tree();
    copy_tree(the_tree, return_value);
  } else {
    return_value = the_tree;
  }

  // Do the trimming.
  return_value = trim_tree_recursive(return_value, the_leaves);

  // Remove length from root (unless root is NULL).
  if (return_value) return_value->has_length = FALSE;

  return(return_value);
}

/*************************************************************************
 * Recursively calculate pairwise distances between distinct leaves.
 * Traverses tree in post-fix order.
 * Returns an array of doubles which are the distances from the leaves
 * to the root.
 *************************************************************************/
static double* compute_leaf_distances(
 TREE_T*  const the_tree,
 STRING_LIST_T* leaf_labels,
 MATRIX_T* matrix
) {

  double* distances[MAX_DEGREE]; // We will have an array of leaf distances
                                 // for each child.
  int i = 0;
  int j = 0;
  int num_leaves = get_num_strings(leaf_labels);

  // For each child, recursively calculate distance to leaves.
  for (i = 0; i < the_tree->num_children; i++) {
    if (is_leaf(the_tree->children[i])) {
      int leaf_index 
        = get_index_in_string_list(the_tree->children[i]->label, leaf_labels);
      if (leaf_index <  0) {
        die(
          "Leaf label %s not found in list of labels\n",
          the_tree->children[i]->label
        );
      }
      distances[i] = mm_calloc(num_leaves, sizeof(double));
      distances[i][leaf_index] = the_tree->children[i]->length;
    }
    else {
      distances[i] 
        = compute_leaf_distances(the_tree->children[i], leaf_labels, matrix);
      for (j = 0; j < num_leaves; j++) {
        if (distances[i][j] > 0.0) {
          distances[i][j] += the_tree->children[i]->length;
        }
      }
    }
  }
  
  // Calculate pairwise leaf distances.
  for (i = 0; i < the_tree->num_children; i++) {
    for (j = i + 1; j < the_tree->num_children; j++) {
      int p = 0;
      int q = 0;
      for (p = 0; p < num_leaves; p++) {
        for (q = 0; q < num_leaves; q++) {
          if (distances[i][p] > 0.0 && distances[j][q] > 0.0) {
            double pairwise_distance = distances[i][p] + distances[j][q];
            set_matrix_cell(p, q, pairwise_distance, matrix);
            set_matrix_cell(q, p, pairwise_distance, matrix);
          }
        }
      }
    }
  }

  // Consolidate distance arrays
  for (i = 1; i < the_tree->num_children; i++) {
    int p = 0;
    for (p = 0; p < num_leaves; p++) {
      distances[0][p] += distances[i][p];
    }
    myfree(distances[i]);
  }
  
  return distances[0];

}

/*************************************************************************
 * Calculate pairwise distances between distinct leaves in a tree.
 * Returns a string labeled matrix. Caller is responsible for freeing matrix.
 *************************************************************************/
RDB_MATRIX_T* get_leaf_distances(TREE_T*  const the_tree) {

  STRING_LIST_T *leaf_labels = make_leaf_list(the_tree);
  int num_leaves = get_num_strings(leaf_labels);
  MATRIX_T* raw_matrix = allocate_matrix(num_leaves, num_leaves);
  double *distances 
    = compute_leaf_distances(the_tree, leaf_labels, raw_matrix);
  RDB_MATRIX_T* pairwise_distances 
    = rdbize_matrix("", leaf_labels, leaf_labels, raw_matrix);
  myfree(distances); // We don't use the list of root-leaf distances.
  free_string_list(leaf_labels);

  return pairwise_distances;
}

/*************************************************************************
 * Free dynamic memory used by a tree.
 *************************************************************************/
void free_tree
  (const BOOLEAN_T free_children, /* Free children as well? */
   TREE_T * a_tree)
{
  int i_child;

  if (a_tree == NULL) {
    return;
  }

  /* Recursively free all the children in the tree. */
  if (free_children == TRUE) {
    for (i_child = 0; i_child < a_tree->num_children; i_child++) {
      free_tree(free_children, a_tree->children[i_child]);
    }
  }

  /* Free the tree struct. */
  myfree(a_tree);
}

/***************************************************************************
 * MAIN
 ***************************************************************************/
#ifdef MAIN

#include "simple-getopt.h"
VERBOSE_T verbosity = NORMAL_VERBOSE;

int main
  (int    argc,
   char * argv[])
{
  char *    tree_filename = NULL;
  FILE *    tree_file = NULL;
  int c = 0;
	int option_count = 1;
  int option_index = 0;
  char* option_name = NULL;
  char* option_value = NULL;
	const char *  message = NULL;
  TREE_T *  the_tree = NULL;
  VERBOSE_T verbosity = NORMAL_VERBOSE;
  cmdoption const options[] = {
    { "verbosity", REQUIRED_VALUE }
  };

  // Define the usage message.
  char      usage[400] = "";
  strcat(usage, "USAGE: tree [options] <tree file>\n");
  strcat(usage, "\n");
  strcat(usage, "   Options:\n");
  strcat(usage, "     -verbosity 1|2|3|4 (default = 2)\n");
  strcat(usage, "\n");

  /* Parse the command line. */
	simple_setopt(argc, argv, option_count, options);
  while(1) {
    // Read the next option, and break if we're done.
    c = simple_getopt(&option_name, &option_value, &option_index);
    if (c == 0) {
      break;
    } else if (c < 0) {
    	simple_getopterror(&message);
      die("Error processing command line options (%s)\n", message);
    }

    if (strcmp(option_name, "verbosity") == 0) {
      verbosity = FALSE;
    }
  }

  // Read the single required argument.
  if (option_index + 1 != argc) {
    fprintf(stderr, "%s", usage);
    exit(1);
  }
  tree_filename = argv[option_index];
  /* Read both trees. */
  if (open_file(tree_filename, "r", 1, "tree", "tree", &tree_file) == 0) {
    exit(1);
  }
  read_tree(tree_file, &the_tree);
  fclose(tree_file);

  /* Write the tree to stdout. */
  write_tree(the_tree, stdout);

  /* Tie up loose ends. */
  free_tree(TRUE, the_tree);
  return(0);
}

#endif
