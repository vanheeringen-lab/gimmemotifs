/**
 * @file sp_matrix.c
 *
 * sp_matrix class. Represents a matrix of s_points.
 *
 **/

#include "partition.h"
#include "sp_matrix.h"
#include "seed.h"
#include "macros.h"
#include <assert.h>

#ifdef PARALLEL
#include <mpi.h>
#endif


#define NODE_NO 0
//#define NODE_NO 1


static void create_partitions (
  SP_MATRIX *sp_mat  ///< The sp_matrix object
  ///< Currently no parameters. May allow user manipulation later.
);


static S_POINT *search_matrix_partition (
  SP_MATRIX *sp_mat, ///< The matrix containing the partition of interest
  DATASET *data,     ///< Sequence dataset.
  MODEL *model,      ///< The nascent motif model. Contains parameters needed
                     ///< in order to calculate significance of scores
  PARTITION *part    ///< The coordinates of the partition within the matrix
);

#ifdef PARALLEL
void union_seed_packets(
  void *f_data,		  // sent seed packets
  void *f_result, 	  // received seed packets
  int *f_length,	  // number seed packets
  MPI_Datatype *datatype  // MPI data type
);
#endif


/**
 * create_spoint_matrix
 *
 * Creates a matrix of S_POINTS with P rows and Q columns in each row, where
 * P = (max_w - min_w + 1), and Q = (max_nsites - min_nsites + 1).
 *
 * The matrix is actually an array of array pointers. Each S_POINT in the matrix
 * is empty: it does not have a best string seed, and its heap is empty. The
 * heap at each S_POINT is created (and left empty). The size of a heap at a
 * given S_POINT is equal to H/F^^d, where d is the manhattan distance to
 * the closest "central" S_POINT (ie S_POINT with central w and nsites values).
 * H = main_hs, and F = (a constant specifying how quickly heap size
 * diminishes).
 *
 * NOTE: A width or nsites is "central" if its spoints have the maximum
 * possible heap size.
 *
 * \return An SP_MATRIX object containing the matrix of empty S_POINTS
 * and other information describing the matrix properties.
 */
extern SP_MATRIX *create_spoint_matrix (
  int *central_ws,   ///< A sorted list of central widths values for the matrix
  int n_ws,          ///< The number of central width values
  int *central_ns,   ///< A sorted list of central nsites values for the matrix
  int n_ns,          ///< The number of central nsites values
  DATASET *dataset   ///< Contains information on heap size and branching factor
) {
  // Matrix is implemented as an array of S_POINT array pointers. Allocate space
  // for this array now:
  S_POINT **matrix = NULL;
  // Calculate total number of rows in matrix:
  int min_w = central_ws[0];
  int max_w = central_ws[n_ws - 1];
  // number of widths in the range [min_w, max_w+1], where (max_w + 1) is only
  // considered if -pal is specified
  int n_widths = max_w - min_w + 2; 
  Resize(matrix, n_widths, S_POINT *);

  // Create a row for each value width from the minimum central width up to
  // the maximum central width:
  int curr_w;
  int w_idx = 0;
  for (curr_w = min_w; curr_w <= max_w; curr_w++) {
    // get list of nsites0 for this width and add it to the matrix:
    matrix[w_idx] =
      create_spoint_row(curr_w, central_ws, n_ws, central_ns, n_ns, dataset);
    w_idx++;
  }

  // Create an SP_MATRIX object containing all the relevant information:
  SP_MATRIX *sp_matrix = NULL;
  sp_matrix = (SP_MATRIX *)mymalloc(sizeof(SP_MATRIX));
  sp_matrix->matrix = matrix;
  sp_matrix->central_ws = central_ws;
  sp_matrix->n_ws = n_ws;
  sp_matrix->central_ns = central_ns;
  sp_matrix->n_ns = n_ns;

  // SP_MATRIX is almost finished; calculate and record partitions to finish
  // it off:
  sp_matrix->partitions = NULL;
  sp_matrix->n_parts = 0;
  create_partitions(sp_matrix);

  return sp_matrix;
} // create_spoint_matrix


/**      
 * create_spoint_row
 *
 * Creates an array of s_points with nsites0 values from the minimum up
 * to the maximum values specified, inclusive.
 *
 * \return The newly created row, which can then be added to a matrix.
 */
extern S_POINT *create_spoint_row (
  int width,         ///< The width of all s_points in the current row
  int *central_ws,   ///< A sorted list of central widths values for the matrix
  int n_ws,          ///< The number of central width values
  int *central_ns,   ///< A sorted list of central nsites values for the matrix
  int n_ns,          ///< The number of central nsites values
  DATASET *dataset   ///< Contains information on heap size and branching factor
) {
  // Allocate memory for the row:
  int min_nsites = central_ns[0];
  int max_nsites = central_ns[n_ns - 1];
  int n_nsites0 = max_nsites - min_nsites + 1; // Size of row
  S_POINT *s_point_row = NULL;
  Resize(s_point_row, n_nsites0, S_POINT);

  /* initialize the starting points, making sure everything is initalized
     so that MPI won't barf. */

  // Create each s_point:
  double curr_nsites;// Current number of sites
  int col_idx;       // Current index in this row
  for (col_idx=0, curr_nsites=min_nsites; curr_nsites <= max_nsites;
       col_idx++, curr_nsites+=1) {
    s_point_row[col_idx].score = LITTLE;
    s_point_row[col_idx].iseq = 0;
    s_point_row[col_idx].ioff = 0;
    s_point_row[col_idx].w0 = width;
    s_point_row[col_idx].nsites0 = curr_nsites<max_nsites ?
      curr_nsites : max_nsites;
    s_point_row[col_idx].wgt_nsites = 0;
    s_point_row[col_idx].e_cons0 = NULL;
    s_point_row[col_idx].cons0 = NULL;
    s_point_row[col_idx].sig = BIG;

    // Calculate the manhattan distance from the current (n,w) position to
    // the closest "central" (n,w) position...

    // Find the distances of the central w and n closest to the current w and
    // n:
    int w_idx;
    int curr_central_w;
    double lowest_w_dist = BIG;
    for (w_idx = 0; w_idx < n_ws; w_idx++) {
      curr_central_w = central_ws[w_idx];
      int curr_dist = (int)abs(curr_central_w - width);
      if (curr_dist < lowest_w_dist) {
        lowest_w_dist = curr_dist;
      }
    }

    int n_idx;
    int curr_central_n;
    double lowest_n_dist = BIG;
    for (n_idx = 0; n_idx < n_ns; n_idx++) {
      curr_central_n = central_ns[n_idx];
      int curr_dist = (int)abs(curr_central_n - curr_nsites); 
      if (curr_dist < lowest_n_dist) {
        lowest_n_dist = curr_dist;
      }
    }

    double manhattan_dist = lowest_w_dist + lowest_n_dist;

    // If branch_W search is to occur, then evaluation should occur at
    // every s_point. Otherwise it should only occur at central s_points:
    if (dataset->branch_params->w_branch && (lowest_n_dist == 0)) {
      s_point_row[col_idx].evaluate = TRUE;
    } else {
      // Initially only evaluate seeds at this s_point if the s_point is
      // central:
      s_point_row[col_idx].evaluate = (manhattan_dist == 0);
    }
    
    // Calculate the heap size for the current S_POINT:
    int max_hs = dataset->main_hs;
    double factor = dataset->hs_decrease;
    int heap_size = MAX((int)max_hs/(pow(factor, manhattan_dist)), 1);

    Resize(s_point_row[col_idx].cons0, width+1, char);
    s_point_row[col_idx].cons0[0] = '\0';
    s_point_row[col_idx].seed_heap = create_heap(
      heap_size,
      (int (*) (void *, void*))compare_seed,
      (void *)copy_seed,
      (void (*)(void*))free_seed,
      (char* (*)(void*))get_str_seed,
      (void *(*)(FILE *, void*))print_seed
    );
  } // col_idx

  return s_point_row;
} // create_spoint_row


/**
 * free_sp_matrix
 *
 */
extern void free_sp_matrix (
  SP_MATRIX *sp_mat  ///< The object being destroyed.  
) {
  // Free the internal representation of the matrix...
  S_POINT **matrix = get_sp_matrix(sp_mat);
  int row_idx, col_idx;
  for (row_idx = 0; row_idx < sp_get_num_rows(sp_mat); row_idx++) {
    S_POINT *curr_row = matrix[row_idx];
    // Free the current row...
    for (col_idx = 0; col_idx < sp_get_num_cols(sp_mat); col_idx++) {
      // Free the current s_point:
      free_s_point(&(curr_row[col_idx]));
    }
    myfree(curr_row);
  }
  myfree(matrix);

  // Free each of the matrix PARTITIONs, and the array referencing them:
  int part_idx;
  PARTITION **partitions = sp_mat->partitions;
  for (part_idx = 0; part_idx < sp_mat->n_parts; part_idx++) {
    myfree(partitions[part_idx]);
  }
  myfree(partitions);

  // Free the sp_matrix object:
  myfree(sp_mat);
} // free_sp_matrix


/**
 * free_s_point
 *
 * Destructor FUNCTION for s_points. Note that s_points were not
 * encapsulated as objects in the first place. Otherwise I would not
 * have to write this function and place it here. NOTE: This function
 * only frees the HEAP inside the spoint; it does not free the memory
 * allocated to the spoint struct itself. Nor does is free the
 * e_cons0 or the cons0 referenced by the S_POINT.
 *
 * It assumes that "destroy_heap" will free the s_point's e_cons0 and cons0.
 * This is a valid assumption if the seed representation (for the s_point)
 * is from one of the seeds in the heap (eg. the "best" seed in the heap).
 * However if this assumption is violated, then the e_cons0 and cons0 for
 * this spoint WILL NOT BE FREED. This could potentially result in a
 * memory leak.
 *
 */
extern void free_s_point (
  S_POINT *sp        ///< The item being destroyed.
) {
  destroy_heap(sp->seed_heap);
} // free_s_point


/**
 * copy_s_point
 *
 * Copies the fields of one spoint into the memory referenced by a second
 * specified S_POINT pointer. This function would go in the "S_POINT" class,
 * if it existed.
 *
 */
extern void copy_s_point (
  S_POINT *s_point,  ///< The S_POINT being copied
  S_POINT *sp_copy   ///< The copy of the S_POINT - OUT
) {
  assert(s_point != NULL);
  assert(sp_copy != NULL);

  sp_copy->score = s_point->score;
  sp_copy->iseq = s_point->iseq;
  sp_copy->ioff = s_point->ioff;
  sp_copy->w0 = s_point->w0;
  sp_copy->nsites0 = s_point->nsites0;
  sp_copy->wgt_nsites = s_point->wgt_nsites;
  sp_copy->evaluate = s_point->evaluate;
  sp_copy->sig = s_point->sig;
  char *e_cons_copy = NULL;
  Resize(e_cons_copy, s_point->w0, char);
  int char_idx;
  for(char_idx = 0; char_idx < s_point->w0; char_idx++) {
    e_cons_copy[char_idx] = (s_point->e_cons0)[char_idx];
  }
  sp_copy->e_cons0 = e_cons_copy;
  char *cons_copy = NULL;
  Resize(cons_copy, (s_point->w0 + 1), char);
  strcpy(cons_copy, s_point->cons0);
  sp_copy->cons0 = cons_copy;
  sp_copy->seed_heap = copy_heap(s_point->seed_heap);
} // copy_s_point


/**
 * print_s_point
 *
 * Print a representation of a specified s_point. This is yet another
 * function that should be a "method" of the s_point class.
 */
extern void print_s_point (
  S_POINT *s_point,  ///< The object being printed.
  FILE *out          ///< The output destination.
) {
  fprintf(out, "--------SPOINT OBJECT:--------\n");
  fprintf(out, "SCORE = %f\n", s_point->score);
  fprintf(out, "iseq = %d\n", s_point->iseq);
  fprintf(out, "ioff = %d\n", s_point->ioff);
  fprintf(out, "w0 = %d\n", s_point->w0);
  fprintf(out, "nsites0 = %f\n", s_point->nsites0);
  fprintf(out, "wgt_nsites = %f\n", s_point->wgt_nsites);
  fprintf(out, "cons0 = %s\n", s_point->cons0);
  fprintf(out, "evaluate = %i\n", s_point->evaluate);
  fprintf(out, "significance = %f\n", s_point->sig);
  fprintf(out, "seed_heap->max_size = %i\n", s_point->seed_heap->max_size);
  fprintf(out, "seed_heap->cur_size = %i\n", s_point->seed_heap->cur_size);
  fprintf(out, "------------------------------\n");
} // print_s_point


/**
 * transfer_final_scores
 *
 * Transfer the scores of the best seeds in the S_POINT heaps into the
 * S_POINTs themselves.
 */
extern void transfer_final_scores (
  SP_MATRIX *sp_matrix ///< This object
) {
  // Proceed through the entire matrix, transfering the details for each
  // S_POINT:
  int row_idx;
  int col_idx;
  for (row_idx = 0; row_idx < sp_get_num_rows(sp_matrix); row_idx++) {
    S_POINT *curr_row = sp_matrix->matrix[row_idx];
    for (col_idx = 0; col_idx < sp_get_num_cols(sp_matrix); col_idx++) {
      S_POINT *curr_sp = curr_row+col_idx;
      HEAP *sp_heap = curr_sp->seed_heap;
      
      if (get_num_nodes(sp_heap) >= 1) {
        SEED *best_seed = (SEED *)(get_node(sp_heap, get_best_node(sp_heap)));
        curr_sp->score = get_seed_score(best_seed);
        curr_sp->iseq = -1; // Seed does not correspond to a location in the dataset.
        curr_sp->ioff = -1; // Seed does not correspond to a location in the dataset.
        curr_sp->e_cons0 = get_e_seed(best_seed);
        curr_sp->cons0 = get_str_seed(best_seed);
      }

      /* If the seed heap of the current s_point is empty, then it could
         mean that no seeds added to the s_point had enough maxima to be
         evaluated by align_top_subsequences. Report this situation:
      */
      else if (TRACE) {
        fprintf(stderr,
                "Heap of spoint was empty, possibly because no seeds had"
                " enough local maxima. w = %i. nsites0 = %f.\n", curr_sp->w0,
                curr_sp->nsites0);
      }
    } // col_idx
  } // row_idx
} // transfer_final_scores


/**
 * get_min_width
 *
 * \return The minimum width encompassed by this sp_matrix.
 *
 */
extern int get_min_width (
  SP_MATRIX *sp_mat  ///< The sp_matrix object.
) {
  return sp_mat->central_ws[0];
}


/**
 * get_max_width
 *
 * \return The maximum width of any s_point in the matrix.
 *
 */
extern int get_max_width (
  SP_MATRIX *sp_mat  ///< The sp_matrix object.
) {
  return sp_mat->central_ws[sp_mat->n_ws - 1];
}


/**
 * get_central_ws
 *
 * \return An array of the central nsites-values in the matrix.
 *
 */
extern int *get_central_ws (
  SP_MATRIX *sp_mat  ///< The sp_matrix object.
) {
  return sp_mat->central_ws;
}


/**
 * get_central_ns
 *
 * \return An array of the central nsites values in the matrix.
 *
 */
extern int *get_central_ns (
  SP_MATRIX *sp_mat  ///< The sp_matrix object.
) {
  return sp_mat->central_ns;
}


/**
 * get_min_nsites
 *
 * \return The minimum nsites encompassed by this sp_matrix.
 *
 */
extern int get_min_nsites (
  SP_MATRIX *sp_mat  ///< The sp_matrix object.
) {
  return sp_mat->central_ns[0];
}


/**
 * get_max_nsites
 *
 * \return The maximum nsites of any s_point in the matrix.
 *
 */
extern int get_max_nsites (
  SP_MATRIX *sp_mat  ///< The sp_matrix object.
) {
  return sp_mat->central_ns[sp_mat->n_ns - 1];
}


/**
 *
 * get_num_central_widths
 *
 * \return The number of "central" width values in the matrix:
 *
 */
extern int get_num_central_widths (
  SP_MATRIX *sp_mat  ///< The sp_matrix object.
) {
  return sp_mat->n_ws;
}


/**
 *
 * get_num_central_nsites
 *
 * \return The number of "central" nsites values in the matrix:
 *
 */
extern int get_num_central_nsites (
  SP_MATRIX *sp_mat  ///< The sp_matrix object.
) {
  return sp_mat->n_ns;
}



/**
 * sp_get_num_rows
 *
 * \return The number of different widths encompassed by this matrix. This
 * is equal to the number of rows in the matrix.
 *
 */
extern int sp_get_num_rows (
  SP_MATRIX *sp_mat
) {
  return (get_max_width(sp_mat) - get_min_width(sp_mat) + 1);
}


/**
 * sp_get_num_cols
 *
 * \return The number of different nsites values encompassed by this matrix.
 * This is equal to the number of columns in the matrix.
 *
 */
extern int sp_get_num_cols (
  SP_MATRIX *sp_mat
) {
  return (get_max_nsites(sp_mat) - get_min_nsites(sp_mat) + 1);
}


/**
 * get_sp_matrix
 *
 * \return An array of pointers to the starts of S_POINT arrays; represents
 * the matrix itself
 *
 */
extern S_POINT **get_sp_matrix (
  SP_MATRIX *sp_mat
) {
  return sp_mat->matrix;
}


/**
 * get_sp_arr
 *
 * Retrieves a row from the spoint matrix. The row number is equal to the
 * width of the starting points/seeds considered in that row.
 *
 * \return An array of all S_POINTS with the specified seed width.
 *
 */
extern S_POINT *get_sp_arr (
  SP_MATRIX *sp_mat, ///< The spoint matrix
  int seed_width     ///< The number of the desired row
) {
  S_POINT **matrix = get_sp_matrix(sp_mat);

  // Index in the matrix is found by calculating index of the specified width
  // relative to the smallest width:
  int row_idx = seed_width - (get_min_width(sp_mat));

  return matrix[row_idx];
}


/**
 * get_spoint
 *
 * Retrieve the spoint located at the specified row and column of the
 * matrix
 *
 * \return A pointer to the apporpriate S_POINT object
 *
 */
extern S_POINT *get_spoint (
  SP_MATRIX *sp_mat, ///< The spoint matrix
  int row_idx,       ///< The row of the desired spoint
  int col_idx        ///< The column of the desired spoint
) {
  S_POINT **matrix = get_sp_matrix(sp_mat);
  return &(matrix[row_idx][col_idx]);
}

/**
 * get_s_point
 *
 * Retrieve the spoint with the given width and nsites.
 *
 * \return A pointer to the apporpriate S_POINT object
 *
 */
extern S_POINT *get_s_point (
  SP_MATRIX *sp_mat, ///< The s_point matrix
  int width,         ///< The width the desired s_point
  int nsites         ///< The nsites of the desired s_point
) {
  S_POINT **matrix = get_sp_matrix(sp_mat);
  // get the minimum width and nsites
  int min_w = get_min_width(sp_mat);
  int min_ns = get_min_nsites(sp_mat);
  // adjust width and nsites index the correct position in the sp_matrix
  return &(matrix[width-min_w][nsites-min_ns]);
}


/**
 * update_s_point_heaps updates the heaps in a specified array of s_points.
 * The specified seed is added to each S_POINT. The score for this seed is
 * possibly different for each S_POINT, and was stored in each of the S_POINTs
 * prior to the calling of this function.
 */
extern void update_s_point_heaps(
  S_POINT *s_points, ///< The s_points whose heaps are to be updated.
  char *str_seed,    ///< The ascii representation of the seed being added to
                     ///< each of the S_POINT heaps.
  int n_nsites0      ///< The length of array s_points.
)
{


  int sp_idx;
  for (sp_idx = 0; sp_idx < n_nsites0; sp_idx++) {
    /* The S_POINT object holds the score and the sequence for the current
       seed (rather than the best seed, as it did prior to Jan 2006)... */
    S_POINT curr_sp = s_points[sp_idx];

    /* If the score of the seed is recorded as LITTLE (which means, in the
     * current calling context of this function, that align_top_subsequences
     * did not evaluate the seed), then ignore this seed; don't add it to the
     * heap.
     */

    if (curr_sp.score > LITTLE) {
      // Construct a new seed from the score and the seed's string:
      SEED *seed_for_heap = new_seed(str_seed, curr_sp.score);

      // Add the node to the current heap:
      SEED *bumped_node = (SEED *)(add_node_heap(curr_sp.seed_heap,
                                                 seed_for_heap));

      // If the heap is full, then it will reject a node. If this happens,
      // then delete the "bumped" node, in order to avoid a memory leak:
      if (bumped_node != NULL) {
        free_seed(bumped_node);
      }
    }
  } // s_point
} // update_s_point_heaps


/**
 * create_partitions
 *
 * Generate a list of PARTITION objects, which store the coordinates of blocks
 * in the matrix. These blocks can then be considered in the future - eg
 * for the purpose of choosing a list of EM starts from the matrix. The
 * array of partitions objects is referenced by the sp_matrix.
 *
 */
static void create_partitions (
  SP_MATRIX *sp_mat  ///< The sp_matrix object
  ///< Currently no parameters. May allow user manipulation later.
) {
  /* For every central s_point (ie s_point with central w and n) in the
   * matrix, calculate the boundaries of the current partition and then
   * generate a partition that fits within those bounds... */
  int *central_widths = get_central_ws(sp_mat);
  int *central_nsites = get_central_ns(sp_mat);

  int w_idx;
  int n_partitions = 0;
  PARTITION **part_array = NULL;
  for (w_idx = 0; w_idx < get_num_central_widths(sp_mat); w_idx++) {
    int n_idx;
    for (n_idx = 0; n_idx < get_num_central_nsites(sp_mat); n_idx++) {
      int curr_w = central_widths[w_idx];
      int curr_n = central_nsites[n_idx];

      /* Calculate the boundaries within which the current partition must
       * be located. The boundaries are half way between the current
       * value and the adjacent value, unless there is no adjacent value
       * (in which case the boundary IS the current value):
       */
      double min_w_bounds, max_w_bounds, min_n_bounds, max_n_bounds;
      if (curr_w == get_min_width(sp_mat)) {
        min_w_bounds = curr_w;
      } else {
        int prev_w = central_widths[w_idx - 1];
        min_w_bounds = (curr_w + prev_w)/(double)2;
      }
      if (curr_w == get_max_width(sp_mat)) {
        max_w_bounds = curr_w;
      } else {
        int next_w = central_widths[w_idx + 1];
        max_w_bounds = (curr_w + next_w)/(double)2;
      }

      if (curr_n == get_min_nsites(sp_mat)) {
        min_n_bounds = curr_n;
      } else {
        int prev_n = central_nsites[n_idx - 1];
        min_n_bounds = (curr_n + prev_n)/(double)2;
      }
      if (curr_n == get_max_nsites(sp_mat)) {
        max_n_bounds = curr_n;
      } else {
        int next_n = central_nsites[n_idx + 1];
        max_n_bounds = (curr_n + next_n)/(double)2;
      }

      /* Calculate the minimum and maximum w and n values which define the
         current partition. The minimum value is the smallest integer that is
         >= the minimum boundary. The maximum value is the largest integer that
         is < the maximum boundary, unless the maximum boundary is the max
         value in the sp_matrix (in which case the maximum value is the
         maximum boundary itself). This ensures that every s_point in the
         matrix will be assigned to a partition.
       */
      int part_min_w, part_max_w, part_min_n, part_max_n;
      part_min_w = (int)ceil(min_w_bounds);
      part_min_n = (int)ceil(min_n_bounds);

      if (curr_w == get_max_width(sp_mat)) {
        part_max_w = curr_w;
      } else {
        // Largest integer LESS than max bounds:
        if (max_w_bounds == ceil(max_w_bounds)) {
          part_max_w = (int)max_w_bounds - 1;
        } else {
          part_max_w = (int)floor(max_w_bounds);
        }
      }

      if (curr_n == get_max_nsites(sp_mat)) {
        part_max_n = curr_n;
      } else {
        // Largest integer LESS than max bounds:
        if (max_n_bounds == ceil(max_n_bounds)) {
          part_max_n = (int)max_n_bounds - 1;
        } else {
          part_max_n = (int)floor(max_n_bounds);
        }
      }

      // Generate the current partition and add it to the growing array:
      PARTITION *curr_part = new_partition(part_min_w, part_max_w, curr_w,
                                           part_min_n, part_max_n, curr_n);
      (n_partitions)++;
      Resize(part_array, n_partitions, PARTITION *);
      part_array[(n_partitions) - 1] = curr_part;
    } // n_idx
  } // w_idx

  assert(sp_mat->partitions == NULL);
  sp_mat->partitions = part_array;
  sp_mat->n_parts = n_partitions;
} // create_partitions


/**
 * choose_em_starts
 *
 * Choose a set of S_POINT objects with the greatest significance (as
 * judged by E-value). Choose a single S_POINT from each partition in
 * the stored array of partitions.
 *
 * \return An array of the most significant s_points. Returns NULL if
 * no S_POINT had a score greater than LITTLE
 *
 */
extern S_POINT *choose_em_starts (
  SP_MATRIX *sp_matrix, ///< The matrix of S_POINT objects to choose from.
  DATASET *data,     ///< Sequence dataset.
  MODEL *model,      ///< The nascent motif model. Contains parameters needed
                     ///< in order to calculate significance of scores.
  int *n_starts      ///< The number of starting points in the array returned
                     ///< - OUT
) {
  assert(sp_matrix->partitions != NULL);
  assert(sp_matrix->n_parts > 0);

  // Array of best s_points:
  S_POINT *sp_arr = NULL;
  int n_spoints = 0;

  // For each partition in the SP_MATRIX, consider all S_POINTS in that
  // partition, and choose the S_POINT with the highest score:
  int part_idx;
  for (part_idx = 0; part_idx < sp_matrix->n_parts; part_idx++) {
    // Get the current partition:
    PARTITION *curr_part = (sp_matrix->partitions)[part_idx];

    S_POINT *best_sp =
      search_matrix_partition(sp_matrix, data, model, curr_part);

    // Check that the best S_POINT is non-null and has a score > LITTLE, and
    // if so, add to the list of chosen S_POINTS:
    if ((best_sp != NULL) && (best_sp->score > LITTLE)) {
      // sp_array needs to accomodate the new s_point pointer => Resize:
      n_spoints++;
      Resize(sp_arr, n_spoints, S_POINT);
      
      // Make a copy of the current S_POINT, as an element in the list of
      // chosen s_points:
      int curr_sp_idx = n_spoints - 1;
      copy_s_point(best_sp, &(sp_arr[curr_sp_idx]));
    }
  } // part_idx

  *n_starts = n_spoints;

  return sp_arr;
} // choose_em_starts


/**
 * search_matrix_partition
 *
 * Searches the specified partition of the matrix, for the S_POINT with the
 * smallest significance value.
 *
 * \return The S_POINT with the smallest significance value. Returns NULL if no
 * S_POINT (in the partition) had a significance less than BIG.
 */
static S_POINT *search_matrix_partition (
  SP_MATRIX *sp_mat, ///< The matrix containing the partition of interest
  DATASET *data,     ///< Sequence dataset.
  MODEL *model,      ///< The nascent motif model. Contains parameters needed
                     ///< in order to calculate significance of scores
  PARTITION *part    ///< The coordinates of the partition within the matrix
) {
  // Consider each S_POINT in the partition. Find the S_POINT with the
  // smallest significance value...
  double curr_best_sig = BIG;
  S_POINT *curr_best_sp = NULL;
  
  int curr_w, curr_n;
  for (curr_w = part->min_w; curr_w <= part->max_w; curr_w++) {
    assert (curr_w <= get_max_width(sp_mat));
    for (curr_n = part->min_n; curr_n <= part->max_n; curr_n++) {
      assert (curr_n <= get_max_nsites(sp_mat));

      int w_idx = curr_w - get_min_width(sp_mat);
      int n_idx = curr_n - get_min_nsites(sp_mat);
      // Get the current S_POINT from the SP_MATRIX:
      S_POINT *curr_sp = get_spoint(sp_mat, w_idx, n_idx);

      // Only consider the current s_point if it has been initialised. Note
      // that (spoint initialised <==> contains non-empty cons0 string):
      double curr_sig = BIG; // s_point is non-significant by default.
      if ((curr_sp->cons0 != NULL) && (strcmp(curr_sp->cons0, "") != 0)) {
        curr_sig = get_log_sig(
                     curr_sp->score,
                     model->mtype,
                     curr_sp->w0,
                     curr_sp->wgt_nsites,
                     curr_sp->nsites0,
                     model->invcomp,
                     model->pal,
                     data
                     );

        curr_sp->sig = curr_sig;

        if (curr_sig <= curr_best_sig) {
          curr_best_sig = curr_sig;
          curr_best_sp = curr_sp;
        }
      } // Only considering s_point if initialised
    } // curr_n
  } // curr_w

  return curr_best_sp;
} // search_matrix_partition


/**
 * print_sp_matrix
 *
 * Print a representation of this matrix.
 *
 */
extern void print_sp_matrix (
  SP_MATRIX *sp_mat, ///< The matrix being printed
  FILE *out          ///< The output destination.
) {
  fprintf(out, "-------------SP_MATRIX OBJECT-------------\n");
  fprintf(out, "Matrix:\n");
  S_POINT **matrix = get_sp_matrix(sp_mat);
  int row_idx, col_idx;
  for (row_idx = 0; row_idx < sp_get_num_rows(sp_mat); row_idx++) {
    for (col_idx = 0; col_idx < sp_get_num_cols(sp_mat); col_idx++) {
      fprintf(out, "row %i, col %i:\n", row_idx, col_idx);
      S_POINT *curr_sp = &(matrix[row_idx][col_idx]);
      print_s_point(curr_sp, out);
    } // col_idx
  } // row_idx

  fprintf(out, "Partitions:\n");
  int part_idx;
  PARTITION **parts = sp_mat->partitions;
  for (part_idx = 0; part_idx < sp_mat->n_parts; part_idx++) {
    print_partition(parts[part_idx], out);
  } // part_idx
  fprintf(out, "-----------END SP_MATRIX OBJECT-----------\n");
} // print_sp_matrix


#ifdef PARALLEL
/**
 * reduce_across_heaps
 *
 * Do a reduction across an array of S_POINT heaps. For each S_POINT in the
 * array, all the seeds from the heaps on each node are combinded (using a
 * union function). A heap containing the best seeds from every node is then
 * propogated to all nodes.
 *
 */
void reduce_across_heaps(
  S_POINT *s_points,     // an array of S_POINTS
  int n_nsites0          // the number of S_POINTS in the s_points array
) 
{
  static int init;
  static MPI_Datatype seed_packet_type;
  static MPI_Op union_seed_packets_op;
  int i_packet;

  // Initialise MPI stuff
  if (init==0){
    init = 1;
    SEED_PACKET seed_packet;
    int block_lengths[4];
    MPI_Aint displacements[4];
    MPI_Aint address[4]; 
    MPI_Datatype typelist[4];

    // Build the derived datatype
    // set the types
    typelist[0]=MPI_DOUBLE;
    typelist[1]=MPI_INT;
    typelist[2]=MPI_INT;
    typelist[3]=MPI_CHAR;

    // set number of elements of each type
    block_lengths[0] = block_lengths[1] = block_lengths[2] = 1;
    block_lengths[3] = MAXSITE;	// the maximum length of a seed

    // calculate the displacements
    MPI_Address(&seed_packet.score, &address[0]);
    MPI_Address(&seed_packet.width, &address[1]);
    MPI_Address(&seed_packet.num_seed_packets, &address[2]);
    MPI_Address(&seed_packet.seed, &address[3]);
    displacements[0]=0;
    displacements[1]=address[1]-address[0];
    displacements[2]=address[2]-address[0];
    displacements[3]=address[3]-address[0];

    // create the derived type
    MPI_Type_struct(4, block_lengths, displacements, typelist, &seed_packet_type);

    // commit the derived type
    MPI_Type_commit(&seed_packet_type);

    // set the MPI reduction operation
    MPI_Op_create(union_seed_packets, FALSE, &union_seed_packets_op);
  } // initialise MPI

  // do a reduction for each s_point in the s_point list
  int sp_idx;
  for (sp_idx = 0; sp_idx < n_nsites0; sp_idx++){
    // package the heap for the spoint at sp_idx in the s_points list
    HEAP *seed_heap = s_points[sp_idx].seed_heap;
    // get the maximum heap size and the number of seeds in the heap
    int max_heap_size = get_max_size(seed_heap);
    int num_seeds = get_num_nodes(seed_heap);

    // set the number of seed packets to the maximum heap size
    SEED_PACKET packets[max_heap_size], best_packets[max_heap_size];
    // set num_seed_packets to the number of filled nodes in the heap (in 
    // case the heap is empty)
    packets[0].num_seed_packets = num_seeds;
    // package each seed in the heap into a seed packet

    for (i_packet = 0; i_packet < num_seeds; i_packet++){
      // set the number of seed_packets that will be filled
      packets[i_packet].num_seed_packets = num_seeds;
      // get the seed at the root
      SEED *curr_seed = pop_heap_root(seed_heap);
      // set the seed packet score
      packets[i_packet].score = get_seed_score(curr_seed);
      // set the width of the string
      packets[i_packet].width = get_width(curr_seed);
      // set the seed
      char *seed_str = get_str_seed(curr_seed);
      strcpy(packets[i_packet].seed, seed_str); 
    }

/*
    // print the packets before the reduction
    if (mpMyID() == NODE_NO){
      fprintf(stdout, "BEFORE\n");
      for (i_packet = 0; i_packet < max_heap_size; i_packet++)
      fprintf(stdout, "node %d packet %d score= %g width= %i seed= %s\n",
                       mpMyID(), i_packet, packets[i_packet].score,
                       packets[i_packet].width, packets[i_packet].seed);
      fflush(stdout);
    }
*/

    // Do the reduction
    MPI_Allreduce((void *)&packets, (void *)&best_packets, max_heap_size,
                    seed_packet_type, union_seed_packets_op, MPI_COMM_WORLD);

/*
    // print the packets after the reduction
    if (mpMyID() == NODE_NO){
      fprintf(stdout, "AFTER\n");
      for (i_packet = 0; i_packet < max_heap_size; i_packet++)
      fprintf(stdout, "node %d packet %d score= %g width= %i seed= %s\n",
                       mpMyID(), i_packet, best_packets[i_packet].score,
                       best_packets[i_packet].width, best_packets[i_packet].seed);
      fflush(stdout);
    }
*/

    // Unpack the best seed packets into the heap

    // Get the number of filled packets
    int num_seed_packets = best_packets[0].num_seed_packets;

    // Add the best seeds to the heap
    for (i_packet = 0; i_packet < num_seed_packets; i_packet++){
      double score =  best_packets[i_packet].score;
      char *seed_str = best_packets[i_packet].seed;
      SEED *best_seed = new_seed(seed_str, score);
      //SEED *bumped_seed = (SEED *)(add_node_heap(seed_heap, best_seed));
      (void *)(add_node_heap(seed_heap, best_seed));
    }
  } // end n_nsites0

} // reduce_across_heaps


/**
 * union_seed_packets
 *
 * Find the union of two seed_packet arrays. Return an array containing 
 * the best seed_packets from this union.
 * 
 * This function is used in reduce_across_heaps.
 *
 */
void union_seed_packets(void *f_data, void *f_result, int *f_length,
                   MPI_Datatype *datatype)
{
  int i;
  int num_seed_packets; 	
  SEED *bumped_seed;
  
  // create a heap to do the heap union
  HEAP *heap = create_heap(
      *f_length, 
      (int (*) (void *, void*))compare_seed,
      (void *)copy_seed,
      (void (*)(void*))free_seed,
      (char* (*)(void*))get_str_seed,
      (void *(*)(FILE *, void*))print_seed
    );
  
  // get the number of seed_packets in f_data
  num_seed_packets = ((SEED_PACKET *)f_data + 0)->num_seed_packets; 

  // unpack the seeds from f_data and add them to the heap
  for (i = 0; i < num_seed_packets; i++){
    // get the data seed
    char *data_seed_str = ((SEED_PACKET *)f_data + i)->seed;
    double data_score = ((SEED_PACKET *)f_data + i)->score;
    SEED *data_seed = new_seed(data_seed_str, data_score);
    // add the seeds to the heap
    bumped_seed = (SEED *)(add_node_heap(heap, data_seed)); 
  } 

  // unpack the seeds from f_result and add them to the heap
  num_seed_packets = ((SEED_PACKET *)f_result + 0)->num_seed_packets;
  for (i = 0; i < num_seed_packets; i++){
    // get the result seed
    char *result_seed_str = ((SEED_PACKET *)f_result + i)->seed;
    double result_score = ((SEED_PACKET *)f_result + i)->score;
    SEED *result_seed = new_seed(result_seed_str, result_score);
    // add the seeds to the heap
    bumped_seed = (SEED *)(add_node_heap(heap, result_seed)); 
  }

  // pack the heap
  int num_seeds = get_num_nodes(heap);
  // set the number of filled packets (in case the heap is empty)
  ((SEED_PACKET *)f_result + 0)->num_seed_packets = num_seeds;
  for (i = 0; i < num_seeds; i++){
    // set the number of seed_packets
    ((SEED_PACKET *)f_result + i)->num_seed_packets = num_seeds;
    // get the index for the seed in the heap
    // (populated heap nodes are at index 1 to num_seeds)
    int heap_idx = i + 1;
    // get the node
    SEED *curr_seed = get_node(heap, heap_idx);
    //double score = get_seed_score(curr_seed);
    ((SEED_PACKET *)f_result + i)->score = get_seed_score(curr_seed);
    char *seed_str = get_str_seed(curr_seed);
    strcpy(((SEED_PACKET *)f_result + i)->seed, seed_str);
  }
} // union_seed_packets

#endif

/**      
 * create_heap_from_sp_matrix
 *
 * Creates a mega heap that contains all the seeds from every s_point in the 
 * sp_matrix.
 *
 * \return A heap containing the seeds from EVERY s_point heap in the matrix
 */

extern HEAP *create_heap_from_sp_matrix (
  SP_MATRIX *sp_mat    // the matrix of s_points 
) {

  int row_idx, col_idx, i;
  int num_seeds = 0;
  int num_rows = sp_get_num_rows(sp_mat);
  int num_cols = sp_get_num_cols(sp_mat);
  void *root, *temp;

  // iterate over the s_points in the sp_matrix to get the total number of seeds
  for (row_idx = 0; row_idx < num_rows; row_idx++) {
    for (col_idx = 0; col_idx < num_cols; col_idx++) {
      S_POINT *current_sp = get_spoint(sp_mat, row_idx, col_idx);
      HEAP *seed_heap = current_sp->seed_heap;
      num_seeds += get_num_nodes(seed_heap); 
    }
  }

  // create the heap
  HEAP *mega_heap = create_heap(
      num_seeds,
      (int (*) (void *, void*))compare_seed,
      (void *)copy_seed,
      (void (*)(void*))free_seed,
      (char* (*)(void*))get_str_seed,
      (void *(*)(FILE *, void*))print_seed
    );

  // add the seeds to the heap
  for (row_idx = 0; row_idx < num_rows; row_idx++) {
    for (col_idx = 0; col_idx < num_cols; col_idx++) {
      S_POINT *current_sp = get_spoint(sp_mat, row_idx, col_idx);
      HEAP *current_heap = current_sp->seed_heap;
      HEAP *seed_heap = copy_heap(current_heap);
      // add copies of the seeds to the mega_heap
      int num_nodes = get_num_nodes(seed_heap);
      for (i=1; i<= num_nodes; i++){
        root = pop_heap_root(seed_heap);
        temp = mega_heap->copy(root);
        temp = add_node_heap(mega_heap, temp);
      }
    }
  }

  // return the heap
  return mega_heap;

} // create_heap_from_sp_matrix



/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
