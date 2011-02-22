/********************************************************************
 * FILE: prior-dist.h
 * AUTHOR: Charles Grant, Bill Noble, Tim Bailey
 * CREATION DATE: 2010-11-19
 * COPYRIGHT: 2010 UW
 *
 * This file contains the public interface for the data structures and
 * functions used to store and manipulate the distribution of priors.
 ********************************************************************/

#ifndef PSP_DIST_H
#define PSP_DIST_H

typedef struct prior_dist PRIOR_DIST_T;

/***********************************************************************
 * Create a new PRIOR_DIST_T object by reading the distribution 
 * from a file.
 ***********************************************************************/
PRIOR_DIST_T *new_prior_dist(const char *filename);

/***********************************************************************
 * Free a PRIOR_DIST_T object.
 ***********************************************************************/
void free_prior_dist(PRIOR_DIST_T *prior_dist);

/***********************************************************************
 * Get minimum prior from PRIOR_DIST_T object.
 ***********************************************************************/
double get_prior_dist_minimum(PRIOR_DIST_T *prior_dist);

/***********************************************************************
 * Get maximum prior from PRIOR_DIST_T object.
 ***********************************************************************/
double get_prior_dist_maximum(PRIOR_DIST_T *prior_dist);

/***********************************************************************
 * Get the array containing distribution from PRIOR_DIST_T object.
 * The caller should not free this array. It will be freed when
 * the PRIOR_DIST_T is freed.
 ***********************************************************************/
ARRAY_T *get_prior_dist_array(PRIOR_DIST_T *prior_dist);

/***********************************************************************
 * Get the length of the array containing the distribution from 
 * a PRIOR_DIST_T object.
 ***********************************************************************/
int get_prior_dist_length(PRIOR_DIST_T *prior_dist);

/***********************************************************************
 * Get the offset for converting an index into the prior dist array
 * into a raw value and visa versa
 ***********************************************************************/
double get_prior_dist_offset(PRIOR_DIST_T *prior_dist);

/***********************************************************************
 * Get the scale for converting an index into the prior dist array
 * into a raw value and visa verse
 * a PRIOR_DIST_T object.
 ***********************************************************************/
double get_prior_dist_scale(PRIOR_DIST_T *prior_dist);

#endif
