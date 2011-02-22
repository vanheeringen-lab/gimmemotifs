/********************************************************************
 * FILE: array.c
 * AUTHOR: William Stafford Noble
 * PROJECT: shared
 * COPYRIGHT: 1999-2008, WSN
 * VERSION: $Revision: 1.1.1.1 $
 * DESCRIPTION: Some simple array-handling routines.
 ********************************************************************/
#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "array.h"


/***********************************************************************
 * Allocate one array.
 ***********************************************************************/
ARRAY_T* allocate_array
  (const int num_items)
{
  ARRAY_T* new_array = (ARRAY_T*)mm_malloc(sizeof(ARRAY_T));
  new_array->items = (ATYPE*)mm_calloc(num_items, sizeof(ATYPE));
  new_array->num_items = num_items;
  return(new_array);
}

/***********************************************************************
 * Resize an array.  Allocates the array if it is null.
 ***********************************************************************/
ARRAY_T* resize_array
  (ARRAY_T* array,
   const int num_items)
{
  ARRAY_T* new_array;

  if (array == NULL) {
    new_array = allocate_array(num_items);
  } else {
    new_array = array;
    mm_resize(new_array->items, num_items, ATYPE);
    new_array->num_items = num_items;
  }
  return(new_array);
}



/**************************************************************************
 * Error checking routine called by all access functions to avoid core
 * dump when attempting to access a null pointer.
 **************************************************************************/
static void check_null_array
  (ARRAY_T* array)
{
#ifdef BOUNDS_CHECK
  if (array == NULL) {
    die("Attempted to access a null array.\n");
  }
#else
  /* Avoid compiler warning. */
  array = NULL;
#endif
}

/**************************************************************************
 * Implement bounds checking.
 **************************************************************************/
static void array_bounds_check
  (int      index,
   ARRAY_T* array)
{

#ifdef BOUNDS_CHECK
  if (index < 0) {
    die("Invalid array index (%d).\n", index);
  } else if (index >= get_array_length(array)) {
    die("Array index out of bounds (%d >= %d).\n", 
	index, get_array_length(array));
  }
#else
  /* Avoid compiler warning. */
  array += index;
#endif
}


/***********************************************************************
 * Check to be sure that two arrays have the same length.
 ***********************************************************************/
static BOOLEAN_T check_array_dimensions
  (BOOLEAN_T die_on_mismatch,
   ARRAY_T*  array1,
   ARRAY_T*  array2)
{
  /* Check to see that the two arrays are of the same length. */
  if (get_array_length(array1) != get_array_length(array2)) {
    if (die_on_mismatch) {
      die("Arrays have differing lengths (%d != %d).\n", 
	  get_array_length(array1), get_array_length(array2));
    }
    return(FALSE);
  }
  return(TRUE);
}


/***********************************************************************
 * Basic access routines.
 ***********************************************************************/
int get_array_length
  (ARRAY_T* array)
{
  check_null_array(array);
  return(array->num_items);
}

/* The following two functions are only used when array bounds
   checking is turned on. Otherwise, they are replaced by macros in
   array.h. */
ATYPE get_array_item_defcheck
  (int      index,
   ARRAY_T* array)
{
  check_null_array(array);
  array_bounds_check(index, array);
  return(array->items[index]);
}

void set_array_item_defcheck
  (int      index,
   ATYPE    value,
   ARRAY_T* array)
{
  check_null_array(array);
  array_bounds_check(index, array);
  array->items[index] = value;
}

void incr_array_item
  (int      index,
   ATYPE    value,
   ARRAY_T* array)
{
  check_null_array(array);
  array_bounds_check(index, array);
  array->items[index] += value;
}

/***********************************************************************
 * Get the minimum and maximum values in an array.
 ***********************************************************************/
ATYPE get_array_minimum
  (ARRAY_T* array)
{
  check_null_array(array);

  int num_items = get_array_length(array);
  if (num_items == 0) {
    die("Attempted to retrieve minimum value from an empty array.\n");
  }

  ATYPE return_value = get_array_item(0, array);
  int i_item;
  for (i_item = 1; i_item < num_items; i_item++) {
    ATYPE item = get_array_item(i_item, array);
    if (item < return_value) {
      return_value = item;
    }
  }
  return(return_value);
}

ATYPE get_array_maximum
  (ARRAY_T* array)
{
  check_null_array(array);

  int num_items = get_array_length(array);
  if (num_items == 0) {
    die("Attempted to retrieve maximum value from an empty array.\n");
  }

  ATYPE return_value = get_array_item(0, array);
  int i_item;
  for (i_item = 1; i_item < num_items; i_item++) {
    ATYPE item = get_array_item(i_item, array);
    if (item > return_value) {
      return_value = item;
    }
  }
  return(return_value);
}

/***********************************************************************
 * Give out the innards of the array.  This function is included only
 * to allow optimizations.  It should be used sparingly.
 ***********************************************************************/
ATYPE* raw_array
  (ARRAY_T* array)
{
  return(array->items);
}

/***********************************************************************
 * Initialize a given array with a given value.
 ***********************************************************************/
void init_array
  (ATYPE    value,
   ARRAY_T* array)
{
  int i_item;
  int num_items;

  num_items = get_array_length(array);
  for (i_item = 0; i_item < num_items; i_item++) {
    set_array_item(i_item, value, array);
  }
}

/***********************************************************************
 * Remove one item from an array, shifting everything else left.
 ***********************************************************************/
void remove_array_item
  (int      item_index,
   ARRAY_T* array)
{
  int i_item;
  int num_items;

  // Shift everything left one position.
  num_items = get_array_length(array);
  for (i_item = item_index + 1; i_item < num_items; i_item++) {
    set_array_item(i_item - 1, get_array_item(i_item, array), array);
  }

  // Reallocate.
  if ((array->items = (ATYPE*)mm_realloc(array->items,
					sizeof(ATYPE) * (num_items - 1)))
      == NULL) {
    die("Error re-allocating array.\n");
  }
  (array->num_items)--;
}

/***********************************************************************
 * Fill an array with a given raw array of values.
 ***********************************************************************/
void fill_array
  (ATYPE*   raw_array,
   ARRAY_T* array)
{
  int i_item;
  int num_items;

  num_items = get_array_length(array);
  for (i_item = 0; i_item < num_items; i_item++) {
    set_array_item(i_item, raw_array[i_item], array);
  }
}


/***********************************************************************
 * Copy an array into another array of the same length.
 ***********************************************************************/
void copy_array
  (ARRAY_T* source,
   ARRAY_T* dest)
{
  int num_items;

  check_null_array(source);
  check_null_array(dest);
  check_array_dimensions(TRUE, source, dest);

  num_items = get_array_length(source);
  if (num_items != 0) {
    memcpy(dest->items, source->items, sizeof(ATYPE) * num_items);
  }
}

/***********************************************************************
 * Determine whether two arrays are equal, within a given bound.
 ***********************************************************************/
BOOLEAN_T equal_arrays
  (ATYPE    close_enough,
   ARRAY_T* array1,
   ARRAY_T* array2)
{
  int i_item;
  int num_items;

  check_null_array(array1);
  check_null_array(array2);

  /* Verify that the arrays are of the same length. */
  if (!check_array_dimensions(FALSE, array1, array2)) {
    return(FALSE);
  }

  /* Check to be sure that each value is the same. */
  num_items = get_array_length(array1);
  for (i_item = 0; i_item < num_items; i_item++) {
    if (!almost_equal(get_array_item(i_item, array1)
		     - get_array_item(i_item, array2), 0.0, close_enough)) {
      return(FALSE);
    }
  }

  return(TRUE);
}

/***********************************************************************
 * Add two arrays.
 ***********************************************************************/
void sum_array
  (ARRAY_T* array1,
   ARRAY_T* array2)
{
  int i_item;
  int num_items;
  
  check_null_array(array1);
  check_null_array(array2);
  check_array_dimensions(TRUE, array1, array2);

  num_items = get_array_length(array1);
  for (i_item = 0; i_item < num_items; i_item++) {
    incr_array_item(i_item, get_array_item(i_item, array1), array2);
  }
}

/***********************************************************************
 * Compute the element-by-element product of two arrays.
 ***********************************************************************/
void element_product
  (ARRAY_T* array1,
   ARRAY_T* array2)
{
  int i_item;
  int num_items;
  
  check_null_array(array1);
  check_null_array(array2);
  check_array_dimensions(TRUE, array1, array2);

  num_items = get_array_length(array1);
  for (i_item = 0; i_item < num_items; i_item++) {
    set_array_item(i_item, get_array_item(i_item, array1) *
		   get_array_item(i_item, array2), array2);
  }
}

/***********************************************************************
 * Compute the dot product of two arrays.
 ***********************************************************************/
double dot_product
  (ARRAY_T* array1,
   ARRAY_T* array2)
{
  int i_item;
  int num_items;
  double return_value;
  
  check_array_dimensions(TRUE, array1, array2);

  num_items = get_array_length(array1);
  return_value = 0.0;
  for (i_item = 0; i_item < num_items; i_item++) {
    return_value += get_array_item(i_item, array1) *
      get_array_item(i_item, array2);
  }
  return(return_value);
}

/***********************************************************************
 * Add a scalar value to each element of an array.
 ***********************************************************************/
void scalar_add
  (ATYPE    value,
   ARRAY_T* array)
{
  int i_item;
  int num_items;
  
  num_items = get_array_length(array);
  for (i_item = 0; i_item < num_items; i_item++) {
    incr_array_item(i_item, value, array);
  }
}

/***********************************************************************
 * Multiply each element of an array by a scalar value.
 ***********************************************************************/
void scalar_mult
  (ATYPE    value,
   ARRAY_T* array)
{
  int i_item;
  int num_items;
  
  check_null_array(array);
  num_items = get_array_length(array);

  for (i_item = 0; i_item < num_items; i_item++) {
    set_array_item(i_item, get_array_item(i_item, array) * value, array);
  }
}

/***********************************************************************
 * Compute the sum of the elements in an array.
 ***********************************************************************/
ATYPE array_total
  (ARRAY_T* array)
{
  int i_item;
  int num_items;
  ATYPE total = 0;
  
  check_null_array(array);
  num_items = get_array_length(array);

  for (i_item = 0; i_item < num_items; i_item++) {
    total += get_array_item(i_item, array);
  }
  return(total);
}

/*************************************************************************
 * Test whether a given array is already sorted.
 *************************************************************************/
BOOLEAN_T is_sorted
  (BOOLEAN_T good_score_is_low,
   ARRAY_T*  my_array)
{
  ATYPE prev_value = get_array_item(0, my_array);
  int index = 0;
  int length = get_array_length(my_array);
  for (index = 1; index < length; index++) {
    ATYPE this_value = get_array_item(index, my_array);
    if (good_score_is_low) {
      if (this_value < prev_value) {
	return(FALSE);
      }
    } else {
      if (this_value > prev_value) {
	return(FALSE);
      }
    }      
    prev_value = this_value;
  }
  return(TRUE);
}

/***********************************************************************
 * Sort the elements in an array. 
 ***********************************************************************/
static int sort_compare
  (const void* elem1,
   const void* elem2)
{
  ATYPE num1 = *((ATYPE*)elem1);
  ATYPE num2 = *((ATYPE*)elem2);

  if (num1 < num2) {		// less
    return(-1);
  } else if (num1 > num2) {	// greater
    return(1);
  } else if (elem1 < elem2) {	// break ties
    return(-1);
  } else if (elem1 > elem2) {	// break ties
    return(1);
  }
    
  return(0);			// identical
}

static int reverse_sort_compare
  (const void* elem1,
   const void* elem2)
{
  ATYPE num1 = *((ATYPE*)elem1);
  ATYPE num2 = *((ATYPE*)elem2);

  if (num1 < num2) {		// less
    return(1);
  } else if (num1 > num2) {	// greater
    return(-1);
  } else if (elem1 < elem2) {	// break ties
    return(1);
  } else if (elem1 > elem2) {	// break ties
    return(-1);
  }
  return(0);

  //if (num1 > num2) {
   // return(-1);
  //} else if (num1 < num2) {
   // return(1);
  //}
}

void sort_array
  (BOOLEAN_T reverse_sort,
   ARRAY_T* array)
{
  if (reverse_sort) {
    qsort(array->items, array->num_items, sizeof(ATYPE), reverse_sort_compare);
  } else {
    qsort(array->items, array->num_items, sizeof(ATYPE), sort_compare);
  }
}

/***********************************************************************
 * Set and get the key used in sorting multiple arrays.
 ***********************************************************************/
void set_array_key
  (ATYPE key,
   ARRAY_T* array)
{
  check_null_array(array);
  array->key = key;
}
ATYPE get_array_key
  (ARRAY_T* array)
{
  check_null_array(array);
  return(array->key);
}

/***********************************************************************
 * Compute the sum of the squares of the elements in an array.
 ***********************************************************************/
ATYPE sum_of_squares
  (ARRAY_T* array)
{
  int i_item;
  int num_items;
  ATYPE value;     /* One value from the array. */
  ATYPE total = 0; /* The sum of the squares. */
  
  check_null_array(array);
  num_items = get_array_length(array);

  for (i_item = 0; i_item < num_items; i_item++) {
    value = get_array_item(i_item, array);
    total += value * value;
  }
  return(total);
}

/***********************************************************************
 * Compute the sum of the squares of the differences between two arrays.
 ***********************************************************************/
ATYPE sum_of_square_diffs
  (ARRAY_T* array1,
   ARRAY_T* array2) 
{
  int i_item;
  int num_items;
  ATYPE diff;      /* The difference between two corresponding array values. */
  ATYPE total = 0; /* The sum of the squares. */
  
  check_null_array(array1);
  check_null_array(array2);
  check_array_dimensions(TRUE, array1, array2);

  num_items = get_array_length(array1);
  for (i_item = 0; i_item < num_items; i_item++) {
    diff = get_array_item(i_item, array1) - get_array_item(i_item, array2);
    total += diff * diff;
  }
  return(total);
}

/***********************************************************************
 * Compute the Euclidean distance between two arrays.
 ***********************************************************************/
double euclidean_distance
  (ARRAY_T* array1,
   ARRAY_T* array2)
{
  return(sqrt(sum_of_square_diffs(array1, array2)));
}

/***********************************************************************
 * Compute the median value in an array.
 *
 * Sorts the array as a side effect.
 ***********************************************************************/
ATYPE compute_median
  (ARRAY_T* array)
{
  int   num_items;
  ATYPE return_value;

  /* Sort the array. */
  sort_array(FALSE, array);

  /* If there are an odd number of elements, return the middle one. */
  num_items = get_array_length(array);
  if (num_items % 2 == 1) {
    return_value = get_array_item(num_items / 2, array);
  }

  /* Otherwise, return the average of the two middle ones. */
  else {
    return_value = get_array_item((num_items / 2) - 1, array);
    return_value += get_array_item(num_items / 2, array);
    return_value /= 2.0;
  }

  return(return_value);
}

/***********************************************************************
 * Read an array of known length from a file.
 ***********************************************************************/
void read_array
  (FILE *   infile,
   ARRAY_T* array)
{
  int   i_item;
  int   num_items;
  ATYPE value;
  
  check_null_array(array);

  num_items = get_array_length(array);
  for (i_item = 0; i_item < num_items; i_item++) {
    if (fscanf((FILE *)infile, ASCAN, &value) != 1) {
      die("Error reading array at position %d.\n", i_item);
    }
    set_array_item(i_item, value, array);
  }
}

/***********************************************************************
 * Read an array of unknown length from a file.
 * Caller is resposbile for freeing array.
 ***********************************************************************/
ARRAY_T *read_array_from_file(const char* filename) {

  int array_size = 100;
  int i_item = 0;
  int num_read = 0;
  ATYPE value;
  
  FILE *array_file = fopen(filename, "r");
  if (array_file == NULL) {
    die(
      "Unable to open file: %s.\nError message: %s.\n", 
      filename, 
      strerror(errno)
    );
  }

  ARRAY_T *array = allocate_array(array_size);
  while ((num_read = fscanf(array_file, ASCAN, &value)) == 1) {
    set_array_item(i_item, value, array);
    ++i_item;
    if (i_item >= array_size) {
      resize_array(array, 2 * array_size);
      array_size = 2 *array_size;
    }
  }

  if (num_read == 0) {
    die("Error reading array at position %d.\n", i_item);
  }
  fclose(array_file);

  resize_array(array, i_item);

  return array;
}

/***********************************************************************
 * Write an array to a file.
 ***********************************************************************/
void print_array
  (ARRAY_T*  array,         /* The array to be printed. */
   int       width,         /* Width of each cell. */
   int       precision,     /* Precision of each cell. */
   BOOLEAN_T eol,           /* Include an EOL char at the end? */
   FILE*     outfile)       /* File to which to write. */
{
  int   i_item;
  int   num_items;
  ATYPE item;
  
  check_null_array(array);

  num_items = get_array_length(array);
  for (i_item = 0; i_item < num_items; i_item++) {
    item = get_array_item(i_item, array);
    fprintf(outfile, APRINT, width, precision, item);
  }
  if (eol) {
    fprintf(outfile, "\n");
  }
}

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
   FILE*     outfile)       /* File to which to write. */
{
  int   i_item;
  int   num_items;
  ATYPE item;
  
  check_null_array(array);

  num_items = get_array_length(array);
  for (i_item = start; (i_item < end) && (i_item < num_items); i_item++) {
    item = get_array_item(i_item, array);
    fprintf(outfile, APRINT, width, precision, item);
  }
  if (eol) {
    fprintf(outfile, "\n");
  }
}


/***********************************************************************
 * Free memory used by an array.
 ***********************************************************************/
void free_array
  (ARRAY_T* array)
{
  if (array == NULL) {
    return;
  }

  myfree(array->items);
  myfree(array);
}

/***********************************************************************
 * Create a bootstrapped copy of the given array.
 *
 * Allocates the bootstrap array; must be freed by caller.
 * Assumes that the random number generator is initialized. 
 ***********************************************************************/
ARRAY_T* bootstrap_array
  (ARRAY_T* source_array, 
   int      num_items)
{
  ARRAY_T* return_array;

  check_null_array(source_array);

  // Allocate the bootstrap array.
  return_array = allocate_array(num_items);

  // Fill up the bootstrap array.
  int i_item;
  int num_inputs = get_array_length(source_array);
  for (i_item = 0; i_item < num_items; i_item++) {
    int random_index = (int)(my_drand() * (double)(num_inputs));
    ATYPE random_item = get_array_item(random_index, source_array);
    set_array_item(i_item, random_item, return_array);
  }
  return(return_array);
}

/* The following functions require non-integral array values. */
#ifdef NOT_INT

/***********************************************************************
 * Divide corresponding elements in two arrays.
 ***********************************************************************/
void dot_divide
  (ARRAY_T* array1,
   ARRAY_T* array2)
{
  int i_item;
  int num_items;
  
  check_null_array(array1);
  check_null_array(array2);
  check_array_dimensions(TRUE, array1, array2);

  num_items = get_array_length(array1);
  for (i_item = 0; i_item < num_items; i_item++) {
    set_array_item(i_item, get_array_item(i_item, array1) /
		   get_array_item(i_item, array2), array2);
  }
}


/***********************************************************************
 * Compute the average of an array.
 ***********************************************************************/
ATYPE ave_array
  (ARRAY_T* array)
{
  int num_items;

  check_null_array(array);  

  /* Check for zero length. */
  num_items = get_array_length(array);
  if (num_items == 0) {
    die("Attempting to average the elements of an empty array.\n");
  }

  /* Compute the average. */
  return(array_total(array) / (ATYPE)num_items);
}

/***********************************************************************
 * Compute the variance of the elements in an array.
 ***********************************************************************/
ATYPE array_variance
  (ARRAY_T* array)
{
  int num_items;
  int i_item;
  ATYPE average;
  ATYPE error;
  ATYPE total_error;

  /* Compute the average. */
  average = ave_array(array);

  total_error = 0.0;
  num_items = get_array_length(array);
  for (i_item = 0; i_item < num_items; i_item++) {
    error = get_array_item(i_item, array) - average;
    total_error += error * error;
  }
  return(total_error / (ATYPE)(num_items - 1));
}
  

/***********************************************************************
 * Make an array sum to zero by subtracting the mean from each element.
 ***********************************************************************/
void sum_to_zero
  (ARRAY_T* array)
{
  int num_items;
  int i_item;
  ATYPE ave;

  ave = ave_array(array);

  num_items = get_array_length(array);
  for (i_item = 0; i_item < num_items; i_item++) {
    incr_array_item(i_item, -ave, array);
  }
}

/***********************************************************************
 * Make an array have variance one by dividing by the standard
 * deviation.
 ***********************************************************************/
void variance_one_array
  (ARRAY_T* array)
{
  ATYPE variance;

  variance = array_variance(array);
  
  /* Avoid divide by zero. */
  if (variance == 0.0) {
    fprintf(stderr, "Warning: variance of zero.\n");
  } else {
    scalar_mult(1.0 / sqrt(array_variance(array)), array);
  }
}

/***********************************************************************
 * Normalize the elements of an array to sum to 1.0.
 ***********************************************************************/
void normalize
  (ATYPE    close_enough, /* If the total is close to 1.0, don't bother. */
   ARRAY_T* array)
{
  ATYPE total;

  /* Compute the sum of the elements in the array. */
  total = array_total(array);

  /* Only normalize if we're relatively far from 1.0. */
  if (almost_equal(total, 1.0, close_enough)) {
    return;
  }

  /* If there's nothing in the array, then return a uniform distribution. */
  if (total == 0.0) {
    init_array(1.0 / (ATYPE)get_array_length(array), array);
    return;
  }

  /* Divide by the total. */
  scalar_mult(1.0 / total, array);
}

/***********************************************************************
 * Compute the Pearson correlation coefficient of two vectors.
 *
 * r(X,Y) = \frac{\sum X_i Y_i - \frac{\sum X_i \sum Y_i}{n}}
 *         {\sqrt{\left( \sum X_i^2 - \frac{(\sum X_i)^2}{n} \right)
 *                 \left( \sum Y_i^2 - \frac{(\sum Y_i)^2}{n}\right)}}
 ***********************************************************************/
ATYPE correlation
  (ARRAY_T* array1,
   ARRAY_T* array2)
{
  ATYPE length;
  ATYPE sum1;
  ATYPE sum2;
  ATYPE dotproduct;
  ATYPE numerator;
  ATYPE variance1;
  ATYPE variance2;
  ATYPE denominator;
  ATYPE return_value;

  length = (ATYPE)get_array_length(array1);
  if (length != (ATYPE)get_array_length(array2)) {
    die("Computing correlation of vectors of different lengths.");
  }

  sum1 = array_total(array1);
  sum2 = array_total(array2);
  dotproduct = dot_product(array1, array2);

  numerator = dotproduct - ((sum1 * sum2) / length);
  
  variance1 = sum_of_squares(array1) - ((sum1 * sum1) / length);
  variance2 = sum_of_squares(array2) - ((sum2 * sum2) / length);

  denominator = sqrt(variance1 * variance2);

  return_value = numerator / denominator;

  /*fprintf(stderr, "return_value=%g\n", return_value);*/
  assert(return_value <= 1.0);
  assert(return_value >= -1.0);
  return(return_value);
}

/***********************************************************************
 * Convert an array to and from logs (base 2).
 ***********************************************************************/
void log_array
  (ARRAY_T* array)
{
  int i_item;
  int num_items;

  check_null_array(array);

  num_items = get_array_length(array);
  for (i_item = 0; i_item < num_items; i_item++) {
    set_array_item(i_item, my_log2(get_array_item(i_item, array)), array);
  }
}

void unlog_array
  (ARRAY_T* array)
{
  int i_item;
  int num_items;

  check_null_array(array);

  num_items = get_array_length(array);
  for (i_item = 0; i_item < num_items; i_item++) {
    set_array_item(i_item, EXP2(get_array_item(i_item, array)), array);
  }
}

/***********************************************************************
 * Compute the sum of an array in log space.
 ***********************************************************************/
ATYPE log_array_total
  (ARRAY_T* array)
{
  int   i_item;
  int   num_items;
  ATYPE total = LOG_ZERO;

  check_null_array(array);

  num_items = get_array_length(array);
  for (i_item = 0; i_item < num_items; i_item++) {
    total = LOG_SUM(total, get_array_item(i_item, array));
  }

  return(total);
}

/***********************************************************************
 * Normalize an array in log space.
 ***********************************************************************/
void log_normalize
  (ATYPE    close_enough,
   ARRAY_T* array)
{
  int i_item;
  int num_items;
  ATYPE total;
  ATYPE this_value;

  /* Get the sum of the elements. */
  total = log_array_total(array);

  /* If the array already sums to zero, don't bother. */
  if (almost_equal(total, 0.0, close_enough)) {
    return;
  }

  /* If there's nothing in the array, then return all zeroes. */
  if (total < LOG_SMALL) {
    init_array(LOG_ZERO, array);
    return;
  }

  num_items = get_array_length(array);
  for (i_item = 0; i_item < num_items; i_item++) {
    this_value = get_array_item(i_item, array) - total;

    /* If this value is small enough, just make it zero. */
    if (this_value < LOG_SMALL) {
      set_array_item(i_item, LOG_ZERO, array);
    } else {
      set_array_item(i_item, this_value, array);
    }
  }
}

/**************************************************************************
 * Convert a given array to or from logs.
 **************************************************************************/
void convert_to_from_log_array
  (BOOLEAN_T to_log,
   ARRAY_T*  source_array,
   ARRAY_T*  target_array)
{
  int num_items;
  int i_item;
  ATYPE new_value;

  // If the source is null, just return.
  if (source_array == NULL) 
    return;

  num_items = get_array_length(source_array);
  for (i_item = 0; i_item < num_items; i_item++) {
    if (to_log) {
      new_value = my_log2(get_array_item(i_item, source_array));
    } else {
      new_value = EXP2(get_array_item(i_item, source_array));
    }
    set_array_item(i_item, new_value, target_array);
  }
}

/***********************************************************************
 * Mix two arrays in log space.
 ***********************************************************************/
void mix_log_arrays
  (float    mixing, /* Percent of array2 that will be retained. */
   ARRAY_T* array1,
   ARRAY_T* array2)
{
  int   i_item;
  int   num_items;
  ATYPE mixed_value;

  check_null_array(array1);
  check_null_array(array2);

  /* Verify that the arrays are of the same length. */
  check_array_dimensions(TRUE, array1, array2);

  /* Verify that we've got a reasonable mixing parameter. */
  if ((mixing > 1.0) || (mixing < 0.0)) {
    die("Invalid mixing parameter (%g).\n", mixing);
  }

  num_items = get_array_length(array1);
  for (i_item = 0; i_item < num_items; i_item++) {
    mixed_value
      = LOG_SUM(my_log2(1.0 - mixing) + get_array_item(i_item, array1),
		my_log2(mixing) + get_array_item(i_item, array2));
    set_array_item(i_item, mixed_value, array2);
  }
}
  
/***********************************************************************
 * Fill an array with random values between 0 and a given maximum.
 *
 * Assumes that the random number generator is initialized. 
 ***********************************************************************/
void randomize_array
  (ATYPE    magnitude,
   ARRAY_T* array)
{
  int num_items;
  int i_item;

  check_null_array(array);

  num_items = get_array_length(array);
  for (i_item = 0; i_item < num_items; i_item++) {
    set_array_item(i_item, my_drand() * magnitude, array);
  }
}

/***********************************************************************
 * Add random noise to an array.
 ***********************************************************************/
void add_noise
  (ATYPE    magnitude,  /* Magnitude of the noise. */
   ARRAY_T* array)
{
  int   i_item;
  int   num_items;
  ATYPE noise;

  check_null_array(array);  

  num_items = get_array_length(array);
  for (i_item = 0; i_item < num_items; i_item++) {
    noise = magnitude * (2 * my_drand() - 1);
    incr_array_item(i_item, noise, array);
  }
}

/***********************************************************************
 * A helper function to normalize a subarray.
 ***********************************************************************/
void normalize_subarray(
  int start_index,
  int length,
  double tolerance,
  ARRAY_T* array
) {
  ATYPE total;
  int   i;

  /* Compute the total. */
  total = 0.0;
  for (i = start_index; i < start_index + length; i++) {
    total += get_array_item(i, array);
  }
  assert(total != 0.0);

  /* Don't bother if we're close enough. */
  if (almost_equal(1.0, total, tolerance)) {
    return;
  }

  /* Divide each element by the total. */
  for (i = start_index; i < start_index + length; i++) {
    set_array_item(i, get_array_item(i, array) / total, array);
  }
}

/***********************************************************************
 * Make all the elements of an array positive by adding a constant to
 * each.
 ***********************************************************************/
void all_positive
  (ARRAY_T* array)
{
  int   i_item;
  int   num_items;
  ATYPE min;

  check_null_array(array);  

  /* Find the minimum value. */
  min = get_array_item(0, array);
  num_items = get_array_length(array);
  for (i_item = 0; i_item < num_items; i_item++) {
    if (get_array_item(i_item, array) < min) {
      min = get_array_item(i_item, array);
    }
  }

  /* If there are negative elements ... */
  if (min < 0.0) {
    /* ... then subtract the minimum from all elements. */
    for (i_item = 0; i_item < num_items; i_item++) {
      incr_array_item(i_item, -min, array);
    }
  }
}
#endif /* NOT_INT */


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
