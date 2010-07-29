/*************************************************************************
 * FILE: reconcile-tree-aligment.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 10 July 2008
 * PROJECT: MEME
 * COPYRIGHT: 2008, UW
 * DESCRIPTION: Remove from a tree all species not in the given alignment
 *              and vice versa.
 *************************************************************************/
#include "alignment.h"
#include "clustalw-io.h"
#include "tree.h"
#include <stdio.h>

/*****************************************************************************
 * MAIN
 *****************************************************************************/
#ifdef MAIN

VERBOSE_T verbosity = NORMAL_VERBOSE;

int main
  (int    argc,
   char * argv[])
{

  // Parse the command line.
  if (argc != 5) {
    fprintf(stderr, "USAGE: reconcile-tree-alignment <input tree> <input alignment> <output tree> <output alignment>\n");
    exit(1);
  }
  char* in_tree_filename = argv[1];
  char* in_alignment_filename = argv[2];
  char* out_tree_filename = argv[3];
  char* out_alignment_filename = argv[4];

  // Read the alignment.
  ALIGNMENT_T* in_alignment = read_alignment_from_file(
    in_alignment_filename, 
    FALSE, // sort by species name
    FALSE, // remove gaps
    NULL   // pointer to ref_seq_index not used
  );
  fprintf(stderr, "Read alignment of %d sequences and %d columns.\n",
	  get_num_aligned_sequences(in_alignment),
	  get_alignment_length(in_alignment));

  // Read the tree.
  FILE* in_tree_file;
  if (open_file(in_tree_filename, "r", FALSE, "tree", "tree", &in_tree_file) == 0) {
    exit(1);
  }
  TREE_T* in_tree;
  read_tree(in_tree_file, &in_tree);
  fprintf(stderr, "Read tree containing %d nodes.\n", get_num_descendants(in_tree, FALSE));

  // Trim the tree, eliminating species not in this alignment.
  STRING_LIST_T* alignment_species = get_species_names(in_alignment);
  TREE_T* out_tree = trim_tree(TRUE, in_tree, alignment_species);

  // Trim the alignment, eliminating species not in the tree.
  STRING_LIST_T* tree_species = make_leaf_list(out_tree);
  ALIGNMENT_T* out_alignment = remove_alignment_seqs(tree_species, in_alignment);

  // Print the alignment.
  FILE* out_alignment_file;
  if (open_file(out_alignment_filename, "w", FALSE, "alignment", "alignment", &out_alignment_file) == 0) {
    exit(1);
  }
  print_clustalw(out_alignment_file, FALSE, out_alignment);
  fprintf(stderr, "Printed alignment of %d sequences and %d columns.\n",
	  get_num_aligned_sequences(out_alignment),
	  get_alignment_length(out_alignment));

  // Print the tree.
  FILE* out_tree_file;
  if (open_file(out_tree_filename, "w", FALSE, "tree", "tree", &out_tree_file) == 0) {
    exit(1);
  }
  write_tree(out_tree, out_tree_file);
  fprintf(stderr, "Printed tree containing %d nodes.\n", get_num_descendants(out_tree, FALSE));

  // Free locally allocated memory.
  free_string_list(alignment_species);
  free_string_list(tree_species);
  free_tree(TRUE, in_tree);
  free_tree(TRUE, out_tree);
  free_alignment(in_alignment);
  free_alignment(out_alignment);

  return(0);
}

#endif
