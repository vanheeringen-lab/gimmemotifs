/********************************************************************
 * FILE: prior-dist.c
 * AUTHOR: Charles Grant, Bill Noble, Tim Bailey
 * CREATION DATE: 2010-11-19
 * COPYRIGHT: 2010 UW
 *
 * This file contains the public interface for the data structures and
 * functions used to store and manipulate the distribution of priors.
 ********************************************************************/

#include "array.h"
#include "utils.h"

typedef struct psp_dist {
  double min;
  double max;
  double offset;
  double scale;
  double range;
  ARRAY_T *dist;
} PRIOR_DIST_T;

/***********************************************************************
 * Create a new PRIOR_DIST_T object by reading the distribution 
 * from a file.
 ***********************************************************************/
PRIOR_DIST_T *new_prior_dist(const char *filename) {
    
  PRIOR_DIST_T *prior_dist = mm_malloc(sizeof(PRIOR_DIST_T));
  prior_dist->dist = read_array_from_file(filename);
  // The first two items in the array contain
  // the minimum and maximium priors.
  prior_dist->min = get_array_item(0, prior_dist->dist);
  remove_array_item(0, prior_dist->dist);
  prior_dist->max = get_array_item(0, prior_dist->dist);
  remove_array_item(0, prior_dist->dist);
  prior_dist->range = get_array_length(prior_dist->dist);
  prior_dist->offset = prior_dist->min;
  prior_dist->scale = prior_dist->range / (prior_dist->max - prior_dist->min);

  return prior_dist;
}

/***********************************************************************
 * Free a PRIOR_DIST_T object.
 ***********************************************************************/
void free_prior_dist(PRIOR_DIST_T *prior_dist) {
  free_array(prior_dist->dist);
  myfree(prior_dist);
}

/***********************************************************************
 * Get minimum prior from PRIOR_DIST_T object.
 ***********************************************************************/
double get_prior_dist_minimum(PRIOR_DIST_T *prior_dist) {
  return prior_dist->min;
}

/***********************************************************************
 * Get maximum prior from PRIOR_DIST_T object.
 ***********************************************************************/
double get_prior_dist_maximum(PRIOR_DIST_T *prior_dist) {
  return prior_dist->max;
}

/***********************************************************************
 * Get the array containing distribution from PRIOR_DIST_T object.
 * The caller should not free this array. It will be freed when
 * the PRIOR_DIST_T is freed.
 ***********************************************************************/
ARRAY_T *get_prior_dist_array(PRIOR_DIST_T *prior_dist) {
  return prior_dist->dist;
}

/***********************************************************************
 * Get the length of the array containing the distribution from 
 * a PRIOR_DIST_T object.
 ***********************************************************************/
int get_prior_dist_length(PRIOR_DIST_T *prior_dist) {
  return get_array_length(prior_dist->dist);
}

/***********************************************************************
 * Get the offset for converting an index into the prior dist array
 * into a raw value and visa versa
 ***********************************************************************/
double get_prior_dist_offset(PRIOR_DIST_T *prior_dist) {
  return prior_dist->offset;
}
/***********************************************************************
 * Get the scale for converting an index into the prior dist array
 * into a raw value and visa verse
 * a PRIOR_DIST_T object.
 ***********************************************************************/
double get_prior_dist_scale(PRIOR_DIST_T *prior_dist) {
  return prior_dist->scale;
}
