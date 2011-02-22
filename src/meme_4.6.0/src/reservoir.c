/******************************************************************************
 * FILE: reservoir.c
 * AUTHOR: Charles Grant and Bill Noble
 * CREATION DATE: 2010-10-22
 * COPYRIGHT: 2010 UW
 * 
 * This file contains the implementation for data structures and functions used 
 * for reservoir sampling.  Based on Knuth's "Semi-Numerical Algorithms", 
 * Algorithm R. pgs. 138-139
 *****************************************************************************/

#include "reservoir.h"
#include "utils.h"

/***********************************************************************
  Data structure for reservoir sampling.
 ***********************************************************************/
struct reservoir_sampler {
  size_t size;
  double *samples;
  size_t num_samples_seen;
  size_t num_samples_retained;
  size_t num_samples_swapped;
};

/***********************************************************************
  Creates a reservoir sampler.
 ***********************************************************************/
RESERVOIR_SAMPLER_T *new_reservoir_sampler(size_t size) {

  RESERVOIR_SAMPLER_T *reservoir = mm_malloc(sizeof(RESERVOIR_SAMPLER_T) * 1);
  reservoir->size = size;
  double *samples = mm_malloc(sizeof(double) * size);
  reservoir->samples = samples;
  reservoir->num_samples_seen = 0;
  reservoir->num_samples_retained = 0;
  reservoir->num_samples_swapped = 0;

  return reservoir;

}

/***********************************************************************
  Frees memory associated with a reservoir sampler
 ***********************************************************************/
void free_reservoir(RESERVOIR_SAMPLER_T *reservoir) {
  myfree(reservoir->samples);
  myfree(reservoir);
}

/***********************************************************************
  Resets the values in the reservoir sampler.
 ***********************************************************************/
void clear_reservoir(RESERVOIR_SAMPLER_T *reservoir) {
  reservoir->num_samples_seen = 0;
  reservoir->num_samples_retained = 0;
  reservoir->num_samples_swapped = 0;
}

/***********************************************************************
  Get the size of the reservoir sampler.
 ***********************************************************************/
size_t get_reservoir_size(RESERVOIR_SAMPLER_T *reservoir) {
  return reservoir->size;
}

/***********************************************************************
  Get the number of samples seen by the reservoir.
 ***********************************************************************/
size_t get_reservoir_num_samples_seen(RESERVOIR_SAMPLER_T *reservoir) {
  return reservoir->num_samples_seen;
}

/***********************************************************************
  Get the number of samples retained in the reservoir sampler.
 ***********************************************************************/
size_t get_reservoir_num_samples_retained(RESERVOIR_SAMPLER_T *reservoir) {
  return reservoir->num_samples_retained;
}

/***********************************************************************
  Retrieve the samples values form the reservoir sampler.
  The user should not free these values.
  They will be freed when the reservoir sampler is freed.
 ***********************************************************************/
double *get_reservoir_samples(RESERVOIR_SAMPLER_T *reservoir) {
  return reservoir->samples;
}

/***********************************************************************
  Submit one sample to the reservoir sampler.
 ***********************************************************************/
void reservoir_sample(RESERVOIR_SAMPLER_T *reservoir, double sample) {
  if (reservoir->num_samples_retained < reservoir->size) {
    // The first samples go directory into the reservoir
    // until it is filled.
    (reservoir->samples)[reservoir->num_samples_retained] = sample;
    ++reservoir->num_samples_seen;
    ++reservoir->num_samples_retained;
  }
  else {
    ++reservoir->num_samples_seen;
    int r = reservoir->num_samples_seen * my_drand();
    if (r < reservoir->size) {
      (reservoir->samples)[r] = sample;
      ++reservoir->num_samples_swapped;
    }
  }
}
