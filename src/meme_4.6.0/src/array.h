/********************************************************************
 * FILE: array.h
 * AUTHOR: William Stafford Noble
 * PROJECT: shared
 * COPYRIGHT: 1999-2008, WSN
 * VERSION: $Revision: 1.1.1.1 $
 * DESCRIPTION: Some simple array-handling routines.
 ********************************************************************/
#ifndef ARRAY_H
#define ARRAY_H

#include <stdio.h>
#include "utils.h"

/********************************************************************
 * By default use double arrays.  Defining IARRAY, SARRAY, or LARRAY
 * changes their type.
 *
 * Note that there would be no sense in defining a float array, since
 * all floats get converted to doubles when they are passed.
 ********************************************************************/
#ifdef IARRAY
#define ATYPE int
#define ATYPENAME "int"
#define ASCAN "%d"
#define APRINT "%*.*d "
#else
#ifdef SARRAY
#define ATYPE short
#define ATYPENAME "short"
#define ASCAN "%hd"
#define APRINT "%*.*d "
#else
#ifdef LARRAY
#define ATYPE long
#define ATYPENAME "long"
#define ASCAN "%ld"
#define APRINT "%*.*d "
#else
#define NOT_INT
#define ATYPE double
#define ATYPENAME "double"
#define ASCAN "%lf"
#define APRINT "%*.*f "
#endif
#endif
#endif

/***********************************************************************
 * Define the array type.  Include internal structure of the type here
 * in order to allow macro access.
 ***********************************************************************/
typedef struct array_t {
  int    num_items;
  ATYPE  key; /* Only used when sorting a matrix. */
  ATYPE* items;
}ARRAY_T;

/***********************************************************************
 * Allocate one array.
 ***********************************************************************/
ARRAY_T* allocate_array
  (const int length);

/***********************************************************************
 * Resize an array.  Allocates the array if array is NULL.
 ***********************************************************************/
ARRAY_T* resize_array
  (ARRAY_T* array,
   const int num_items);

/***********************************************************************
 * Basic access routines.
 ***********************************************************************/
int get_array_length
  (ARRAY_T* array);

#ifdef BOUNDS_CHECK
#define get_array_item(index,array) get_array_item_defcheck(index,array)
#define set_array_item(index,value,array) set_array_item_defcheck(index,value,array)
#else
#define get_array_item(index,array) ((ARRAY_T*)array)->items[index]
#define set_array_item(index,value,array) ((ARRAY_T *)array)->items[index] \
                                         = value
#endif 

/* The following two functions are only used when array bounds
   checking is turned on. Otherwise, they are replaced by macros
   above. */
ATYPE get_array_item_defcheck
  (int      index,
   ARRAY_T* array);

void set_array_item_defcheck
  (int      index,
   ATYPE    value,
   ARRAY_T* array);

void incr_array_item
  (int      index,
   ATYPE    value,
   ARRAY_T* array);

/***********************************************************************
 * Get the minimum and maximum values in an array.
 ***********************************************************************/
ATYPE get_array_minimum
  (ARRAY_T* array);

ATYPE get_array_maximum
  (ARRAY_T* array);

/***********************************************************************
 * Give out the innards of the array.  This function is included only
 * to allow optimizations.  It should be used sparingly.
 ***********************************************************************/
ATYPE* raw_array
  (ARRAY_T* array);

/***********************************************************************
 * Initialize a given array with a given value.
 ***********************************************************************/
void init_array
  (ATYPE    value,
   ARRAY_T* array);

/***********************************************************************
 * Remove one item from an array, shifting everything else left.
 ***********************************************************************/
void remove_array_item
  (int      item_index,
   ARRAY_T* array);

/***********************************************************************
 * Fill an array with a given raw array of values.
 ***********************************************************************/
void fill_array
  (ATYPE*   raw_array,
   ARRAY_T* array);

/***********************************************************************
 * Copy an array into another array of the same length.
 ***********************************************************************/
void copy_array
  (ARRAY_T* source_array,
   ARRAY_T* target_array);

/***********************************************************************
 * Determine whether two arrays are equal, within a given bound.
 ***********************************************************************/
BOOLEAN_T equal_arrays
  (ATYPE    close_enough,
   ARRAY_T* array1,
   ARRAY_T* array2);

/***********************************************************************
 * Add two arrays.
 ***********************************************************************/
void sum_array
  (ARRAY_T* array1,
   ARRAY_T* array2);

/***********************************************************************
 * Compute the element-by-element product of two arrays.
 *
 * Store the result in the second array.
 ***********************************************************************/
void element_product
  (ARRAY_T* array1,
   ARRAY_T* array2);

/***********************************************************************
 * Compute the dot product of two arrays.
 ***********************************************************************/
double dot_product
  (ARRAY_T* array1,
   ARRAY_T* array2);

/***********************************************************************
 * Compute the Euclidean distance between two arrays.
 ***********************************************************************/
double euclidean_distance
  (ARRAY_T* array1,
   ARRAY_T* array2);

/***********************************************************************
 * Add a scalar value to each element of an array.
 ***********************************************************************/
void scalar_add
  (ATYPE    value,
   ARRAY_T* array);

/***********************************************************************
 * Multiply each element of an array by a scalar value.
 ***********************************************************************/
void scalar_mult
  (ATYPE    value,
   ARRAY_T* array);

/***********************************************************************
 * Compute the sum of the elements in an array.
 ***********************************************************************/
ATYPE array_total
  (ARRAY_T* array);

/*************************************************************************
 * Test whether a given array is already sorted.
 *************************************************************************/
BOOLEAN_T is_sorted
  (BOOLEAN_T good_score_is_low,
   ARRAY_T*  my_array);

/***********************************************************************
 * Sort the elements in an array. 
 ***********************************************************************/
void sort_array
  (BOOLEAN_T reverse_sort,
   ARRAY_T*  array);

/***********************************************************************
 * Set and get the key used in sorting multiple arrays.
 ***********************************************************************/
void set_array_key
  (ATYPE key,
   ARRAY_T* array);

ATYPE get_array_key
  (ARRAY_T* array);

/***********************************************************************
 * Compute the sum of the squares of the elements in an array.
 ***********************************************************************/
ATYPE sum_of_squares
  (ARRAY_T* array);

/***********************************************************************
 * Compute the sum of the squares of the differences between two arrays.
 ***********************************************************************/
ATYPE sum_of_square_diffs
  (ARRAY_T* array1,
   ARRAY_T* array2);

/***********************************************************************
 * Compute the median value in an array.
 *
 * Sorts the array as a side effect.
 ***********************************************************************/
ATYPE compute_median
  (ARRAY_T* array);

/***********************************************************************
 * Read an array of known length from a file.
 ***********************************************************************/
void read_array
  (FILE *   infile,
   ARRAY_T* array);

/***********************************************************************
 * Read an array of unknown length from a file.
 * Caller is resposbile for freeing array.
 ***********************************************************************/
ARRAY_T *read_array_from_file(const char* filename);

/***********************************************************************
 * Write an array to a file.
 ***********************************************************************/
void print_array
  (ARRAY_T*  array,         /* The array to be printed. */
   int       width,         /* Width of each cell. */
   int       precision,     /* Precision of each cell. */
   BOOLEAN_T eol,           /* Include an EOL char at the end? */
   FILE*     outfile);      /* File to which to write. */

/***********************************************************************
 * Write an array to a file.
 ***********************************************************************/
void print_sub_array
  (int       start,         // First position to print.
   int       end,           // Last position to print.
   ARRAY_T*  array,         /* The array to be printed. */
   int       width,         /* Width of each cell. */
   int       precision,     /* Precision of each cell. */
   BOOLEAN_T eol,           /* Include an EOL char at the end? */
   FILE*     outfile);      /* File to which to write. */

/***********************************************************************
 * Free memory used by an array.
 ***********************************************************************/
void free_array
  (ARRAY_T* array);


#ifdef NOT_INT

/***********************************************************************
 * Divide corresponding elements in two arrays.
 ***********************************************************************/
void dot_divide
  (ARRAY_T* array1,
   ARRAY_T* array2);

/***********************************************************************
 * Create a bootstrapped copy of the given array.
 *
 * Both arrays must be pre-allocated to the same size.
 * Assumes that the random number generator is initialized. 
 ***********************************************************************/
ARRAY_T* bootstrap_array
  (ARRAY_T* source_array, 
   int      max_size);

/***********************************************************************
 * Compute the average of an array.
 ***********************************************************************/
ATYPE ave_array
  (ARRAY_T* array);

/***********************************************************************
 * Compute the variance of the elements in an array.
 ***********************************************************************/
ATYPE array_variance
  (ARRAY_T* array);

/***********************************************************************
 * Make an array sum to zero by subtracting the mean from each element.
 ***********************************************************************/
void sum_to_zero
  (ARRAY_T* array);

/***********************************************************************
 * Divide each element in an array by the standard deviation.
 ***********************************************************************/
void variance_one_array
  (ARRAY_T* array);

/***********************************************************************
 * Normalize the elements of an array to sum to 1.0.
 ***********************************************************************/
void normalize
  (ATYPE    close_enough, /* If the total is close to 1.0, don't bother. */
   ARRAY_T* array);

/***********************************************************************
 * Compute the Pearson correlation coefficient of two vectors.
 *
 * r(X,Y) = \frac{\sum X_i Y_i - \frac{\sum X_i \sum Y_i}{n}}
 *         {\sqrt{\left( \sum X_i^2 - \frac{(\sum X_i)^2}{n} \right)
 *                 \left( \sum Y_i^2 - \frac{(\sum Y_i)^2}{n}\right)}}
 ***********************************************************************/
ATYPE correlation
  (ARRAY_T* array1,
   ARRAY_T* array2);

/***********************************************************************
 * Convert an array to and from logs (base 2).
 ***********************************************************************/
void log_array
  (ARRAY_T* array);

void unlog_array
  (ARRAY_T* array);

/***********************************************************************
 * Compute the sum of an array in log space.
 ***********************************************************************/
ATYPE log_array_total
  (ARRAY_T* array);

/***********************************************************************
 * Normalize an array in log space.
 ***********************************************************************/
void log_normalize
  (ATYPE    close_enough,
   ARRAY_T* array);

/**************************************************************************
 * Convert a given array to or from logs.
 **************************************************************************/
void convert_to_from_log_array
  (BOOLEAN_T to_log,
   ARRAY_T*  source_array,
   ARRAY_T*  target_array);

/***********************************************************************
 * Mix two arrays in log space.
 ***********************************************************************/
void mix_log_arrays
  (float    mixing, /* Percent of array2 that will be retained. */
   ARRAY_T* array1,
   ARRAY_T* array2);

/***********************************************************************
 * Fill an array with random values between 0 and a given maximum.
 *
 * Assumes that the random number generator is initialized. 
 ***********************************************************************/
void randomize_array
  (ATYPE    magnitude,
   ARRAY_T* array);

/***********************************************************************
 * Add random noise to an array.
 ***********************************************************************/
void add_noise
  (ATYPE    magnitude,  /* Magnitude of the noise. */
   ARRAY_T* array);

/***********************************************************************
 * Make all the elements of an array positive by adding a constant to
 * each.
 ***********************************************************************/
void all_positive
  (ARRAY_T* array);

/***********************************************************************
 * A helper function to normalize a subarray.
 ***********************************************************************/
void normalize_subarray(
  int start_index,
  int length,
  double tolerance,
  ARRAY_T* array
);

#endif
#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
