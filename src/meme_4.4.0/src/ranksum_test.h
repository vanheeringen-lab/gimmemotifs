/**
 * @file ranksum_test.h
 *
 * This module implements the statistical ranksum_test.c also known
 * as Mann-Whitney U-test or Mann-Whitney Wilcoxon test
 *
 * author: Fabian Buske
 * version: 1.0 (26.06.2008)
*/

#ifndef RANKSUM_TEST_H_
#define RANKSUM_TEST_H_

#include "utils.h"


typedef struct ranksum_result RSR_T;
typedef struct ranksum_datapoint RSDP_T;
typedef struct ranksum_data RSD_T;

/**********************************************************************
  ranksum_from_stats()

  Computes the ranksum test for a given statistic.
**********************************************************************/
RSR_T *ranksum_from_stats(
   int n,			/* number of samples */
   int na,  		/* number of positives (class a)*/
   double ta_obs	/* sum of positive ranks */
);

/**********************************************************************
  ranksum_from_sets()

  Calculates the ranksum test for two sets of samples and return a
  ranksum test result struct.
**********************************************************************/
RSR_T *ranksum_from_sets(
	double* sample_a,  	/* a set (a) of double values */
	int size_a, 		/* the size of sample set a */
	double* sample_b, 	/* a set (b) of double values */
	int size_b 			/* the size of sample set b */
);

/**********************************************************************
  ranksum_from_groups()

  Calculates the ranksum test for a set of samples which classes are
  specified in a corresponding array.
**********************************************************************/
RSR_T *ranksum_from_groups(
	double* samples,  	/* the sample set */
	BOOLEAN_T* group,	/* set of the group each samples belongs to */
	int size 			/* the size of sample set  */
);

/********************************************************************
  get_ranksum_rank()

  Gets the rank of the position
**********************************************************************/
double get_ranksum_rank(RSD_T *dataset, int position);

/**********************************************************************
  set_ranksum_groups()

  Assigns the given groups to the dataset structure. Both have to
  be of the same size.

  Most beneficial when used in conjuction with get_ranksum_dataset and
  run_ranksum_test_on_ordered_dataset i.e. when the ranking (ordering)
  does not change between the ranksum tests but the groups (classes) do.
**********************************************************************/
void set_ranksum_groups(
	RSD_T *dataset,		/* the dataset */
	BOOLEAN_T* groups	/* the groups to be assigned to the data*/
);

/**********************************************************************
  set_ranksum_group()

  Assigns the given position to the given group in the dataset structure. 
  For obvious reasons the position must exist.

  Most beneficial when used in conjuction with get_ranksum_dataset and
  run_ranksum_test_on_ordered_dataset i.e. when the ranking (ordering)
  does not change between the ranksum tests but the groups (classes) do.
**********************************************************************/
void set_ranksum_group(RSD_T *dataset, int position, BOOLEAN_T group);

/**********************************************************************
  get_ranksum_dataset()

  Creates a ranksum input data set that can be (re-)used with the
  function run_ranksum_test() or run_ranksum_test_on_ordered_dataset()

  Copies the list of samples and groups to the new data structure.
  Calculates the (adjusted) rank of each score.

  NOTE: The samples are expected to be sorted already. The arrays
        of samples and groups have to be of the specified size and each
        entry in group corresponds to the sample at the same index.
**********************************************************************/
RSD_T *get_ranksum_dataset(
	double* samples,  	/* the sample set */
	BOOLEAN_T* group,	/* set of the group each samples belongs to */
	int size 			/* the size of sample set  */
);

/*********************************************************************
  get_ranksum_dataset2()
 
  As above except uses accessor functions to allow copying of the
  samples and groups so you don't have to specially copy your data 
  into arrays of double and boolean.
 
 ********************************************************************/
RSD_T *get_ranksum_dataset2(
	void* samples,  	/* the sample set */
  double (*get_sample)(const void*, int), /* sample accessor */
	void* group,	/* set of the group each samples belongs to */
  BOOLEAN_T (*get_group)(const void*, int), /* group accessor */
	int size 			/* the size of sample set  */
);


/**********************************************************************
  run_ranksum_test()

  Runs the ranksum test on a given dataset. Use THIS function, if
  the dataset has not NOT ORDERED manually beforehand.
**********************************************************************/
RSR_T* run_ranksum_test(
	RSD_T *dataset 	// the dataset to work with
);

/**********************************************************************
  run_ranksum_test_on_ordered_dataset()

  Runs the ranksum test on a given already ordered dataset
**********************************************************************/
RSR_T* run_ranksum_test_on_ordered_dataset(
	RSD_T *dataset 	// the ordered dataset to work with
);

/**********************************************************************
  RSR_get_p_left()

  returns the p-value for the left side (lesser)
**********************************************************************/
double RSR_get_p_left(RSR_T* ranksum_result);

/**********************************************************************
  RSR_get_p_right()

  returns the p-value for the right side (greater)
**********************************************************************/
double RSR_get_p_right(RSR_T* ranksum_result);

/**********************************************************************
  RSR_get_p_onetailed()

  returns the p-value for the one tailed side (the smaller one of left
  and right)
**********************************************************************/
double RSR_get_p_onetailed(RSR_T* ranksum_result);

/**********************************************************************
  RSR_get_p_twotailed()

  returns the p-value for the twotailed test
**********************************************************************/
double RSR_get_p_twotailed(RSR_T* ranksum_result);

/**********************************************************************
  RSR_get_u()

  returns the u value for the class considered
**********************************************************************/
double RSR_get_u(RSR_T*  ranksum_result);

/**********************************************************************
  destroy_rsdt()

  Frees the memory, i.e. removes the pointers without deleting the
  final data points. Use destroy_rsd() if all memory should be freed
  instead.
**********************************************************************/
void destroy_rsd(RSD_T *datapoint);

/**********************************************************************
  destroy_rsr()

  Frees the memory of a RSR_T* data structure
**********************************************************************/
void destroy_rsr(RSR_T *dataset);

#endif /*RANKSUM_TEST_H_*/
