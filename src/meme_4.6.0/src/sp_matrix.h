/**
 * @file sp_matrix.h
 *
 * sp_matrix class. Represents a matrix of s_points.
 *
 **/

#ifndef SP_MATRIX_H
#define SP_MATRIX_H
#include "meme.h"
#include "partition.h"
#ifdef DMALLOC
#include "dmalloc.h"
#endif

/**
 * An sp_matrix type. Stores the matrix, as an array of S_POINT pointers, as
 * well as information about the dimensions of the matrix and which rows
 * (ie width values) are "central", where substring search will occur.
 */
typedef struct sp_matrix {
  S_POINT **matrix;  // The matrix itself
  int *central_ws;   // Central width values. Smallest and largest indicate
                     // matrix bounds
  int n_ws;          // Number of central widths in the matrix.
  int *central_ns;   // Central nsites values. Smallest and largest indicate
                     // matrix bound.
  int n_ns;          // Number of central nsites values in the matrix
  PARTITION **partitions; // An array of partition objects for the matrix
  int n_parts;       // The number of partitions in the matrix
} SP_MATRIX;

extern SP_MATRIX *create_spoint_matrix (
  int *central_ws,   ///< A sorted list of central widths values for the matrix
  int n_ws,          ///< The number of central width values
  int *central_ns,   ///< A sorted list of central nsites values for the matrix
  int n_ns,          ///< The number of central nsites values
  DATASET *dataset   ///< Contains information on heap size and branching factor
);


extern S_POINT *create_spoint_row (
  int width,         ///< The width of all s_points in the current row
  int *central_ws,   ///< A sorted list of central widths values for the matrix
  int n_ws,          ///< The number of central width values
  int *central_ns,   ///< A sorted list of central nsites values for the matrix
  int n_ns,          ///< The number of central nsites values
  DATASET *dataset   ///< Contains information on heap size and branching factor
);


extern int get_hs(
  int curr_width,    ///< Current width under consideration
  int min_width,     ///< Minimum possible width value
  int max_width,     ///< Maximum possible width value
  int w_dist,        ///< Distance between width values
  int main_w_hs,     ///< Heap size at "main" width values
  DATASET *data      ///< Contains info used for determining heapsize
);


extern void transfer_final_scores (
  SP_MATRIX *sp_matrix ///< This object
);


extern void branching_search (
  BRANCH_PARAMS *branch_params, ///< The parameters controlling branching search
  MODEL *model,      ///< The motif model (nascent)
  DATASET *dataset,  ///< The dataset of sequences
  SP_MATRIX *sp_mat, ///< A matrix of S_POINT objects. The heaps in these
                     ///< objects will be branched from and updated with
                     ///< the resulting novel seeds.
  HASH_TABLE evaluated_seed_ht ///< A hashtable recording which seeds have been
                     ///< evaluated at some point previously.
                     ///< Freq x letter map (for converting seed to theta)
);


extern int sp_get_num_rows (
  SP_MATRIX *sp_mat
);


extern int sp_get_num_cols (
  SP_MATRIX *sp_mat
);


extern S_POINT **get_sp_matrix (
  SP_MATRIX *sp_mat
);


extern int get_num_central_widths (
  SP_MATRIX *sp_mat  ///< The sp_matrix object.
);


extern int get_num_central_nsites (
  SP_MATRIX *sp_mat  ///< The sp_matrix object.
);


extern int *get_central_ws (
  SP_MATRIX *sp_mat  ///< The sp_matrix object.
);

extern int *get_central_ns (
  SP_MATRIX *sp_mat  ///< The sp_matrix object.
);


extern int get_min_width (
  SP_MATRIX *sp_mat  ///< The sp_matrix object.
);


extern int get_max_width (
  SP_MATRIX *sp_mat  ///< The sp_matrix object.
);


extern int get_min_nsites (
  SP_MATRIX *sp_mat  ///< The sp_matrix object.
);


extern int get_max_nsites (
  SP_MATRIX *sp_mat  ///< The sp_matrix object.
);


extern S_POINT *get_sp_arr (
  SP_MATRIX *sp_mat, ///< The spoint matrix
  int row_idx        ///< The number of the desired row
);


extern S_POINT *get_spoint (
  SP_MATRIX *sp_mat, ///< The spoint matrix
  int row_idx,       ///< The row of the desired spoint
  int col_idx        ///< The column of the desired spoint
);

extern S_POINT *get_s_point (
  SP_MATRIX *sp_mat, ///< The spoint matrix
  int width,         ///< The width the desired spoint
  int nsites         ///< The nsites of the desired spoint
);

extern void update_s_point_heaps(
  S_POINT *s_points, ///< The s_points whose heaps are to be updated.
  char *str_seed,    ///< The ascii representation of the seed being added to
                     ///< each of the S_POINT heaps.
  int n_nsites0      ///< The length of array s_points.
);


extern S_POINT *choose_em_starts (
  SP_MATRIX *sp_matrix, ///< The matrix of S_POINT objects to choose from.
  DATASET *data,     ///< Sequence dataset.
  MODEL *model,      ///< The nascent motif model. Contains parameters needed
                     ///< in order to calculate significance of scores.
  int *n_starts      ///< The number of starting points in the array returned
                     ///< - OUT
);


extern void free_s_point (
  S_POINT *sp        ///< The item being destroyed.
);


extern void copy_s_point (
  S_POINT *s_point,  ///< The S_POINT being copied
  S_POINT *sp_copy   ///< The copy of the S_POINT - OUT
);


extern void print_s_point (
  S_POINT *s_point,  ///< The object being printed.
  FILE *out          ///< The output destination.
);


extern void free_sp_matrix (
  SP_MATRIX *sp_mat  ///< The object being destroyed.  
);


extern void print_sp_matrix (
  SP_MATRIX *sp_mat, ///< The matrix being printed
  FILE *out          ///< The output destination.
);

// MPI FUNTIONS
typedef struct seed_packet {
  double score;         // the score for the seed       
  int width;            // width of the seed
  int num_seed_packets; // the total number of seeds
  char seed[MAXSITE];   // the seed sequence
} SEED_PACKET;

extern void reduce_across_heaps(
  S_POINT *s_points,    // an array of S_POINTS
  int n_nsites0         // the number of S_POINTS in the s_points array
);

extern HEAP *create_heap_from_sp_matrix (
  SP_MATRIX *sp_mat	// the matrix of s_points 
);

#endif


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
