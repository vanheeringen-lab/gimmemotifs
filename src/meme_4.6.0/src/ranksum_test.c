/**
 * @file ranksum_test.c
 *
 * This module implements the statistical ranksum_test.c also known
 * as Mann-Whitney U-test or Mann-Whitney Wilcoxon test
 *
 * author: Fabian Buske
 * version: 1.0 (26.06.2008)
 *
*/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "ranksum_test.h"
#include "utils.h"

#define min(a,b)      (a<b)?a:b
#define DATA_INCREMENT  1000
#define sqrt2 sqrt(2.0)


/**********************************************************************
  ranksum_result

  this struct contains the results of a ranksum calculation
**********************************************************************/
struct ranksum_result{
	double p_left, p_right, p_twotailed; // p_values
	double u;			// u values
};

/**********************************************************************
  ranksum_datapoint()

  this struct contains one data point for the ranksum statistics
  consisting of a score an a group
**********************************************************************/
struct ranksum_datapoint{
	double score;   // measurement
	double rank;	// the rank
	BOOLEAN_T group;// group
};

/**********************************************************************
  ranksum_data()

  this struct contains a set of data points
**********************************************************************/
struct ranksum_data{
	int size;   	// number of data points
	int allocated_size;
	RSDP_T **data_points; // group
};

/**********************************************************************
  adjust_ranking()

  "Re-ranks" the list according to same value elements (later get the
  average rank assigned)
**********************************************************************/
void adjust_ranking(RSD_T *dataset);

/**********************************************************************
  allocate_ranksum_result()

  Constructor for the ranksum result data structure
**********************************************************************/
RSR_T *allocate_ranksum_result();

/**********************************************************************
  allocate_and_init_ranksum_data()

  Constructor for the ranksum data structure which allocates memory for
  up to size datapoints
**********************************************************************/
RSD_T *allocate_and_init_ranksum_data(unsigned int size);

/**********************************************************************
  allocate_ranksum_datapoint()

  Constructor for the ranksum data structure
**********************************************************************/
RSDP_T *allocate_ranksum_datapoint(double score, BOOLEAN_T group);

/**********************************************************************
  allocate_ranksum_datapoint()

  Constructor for the ranksum data structure
**********************************************************************/
RSD_T *allocate_ranksum_data();

/*************************************************************************
 * add_ranksum_datapoint()
 *
 * Adds a data point to the ranksum dataset
 *************************************************************************/
void add_ranksum_datapoint(RSD_T *data, RSDP_T* datapoint);

/**********************************************************************
  free_rsd()

  Frees the memory, i.e. removes the pointers without deleting the
  final data points. Use destroy_rsd() if all memory should be freed
  instead.
**********************************************************************/
void free_rsd(RSD_T *dataset);

/**********************************************************************
  allocate_ranksum_result()

  Constructor for the ranksum result data structure
**********************************************************************/
RSR_T *allocate_ranksum_result(){
	// Allocate memory and initialize fields
	RSR_T *data = mm_malloc(sizeof(RSR_T));
	data->u = 0.0;
	data->p_left = 0.0;
	data->p_right = 0.0;
	data->p_twotailed = 0.0;
	return data;
};

/**********************************************************************
  allocate_ranksum_datapoint()

  Constructor for the ranksum data structure
**********************************************************************/
RSDP_T *allocate_ranksum_datapoint(
		double score,
		BOOLEAN_T group
){
	// Allocate memory and initialize fields
	RSDP_T *datapoint = mm_malloc(sizeof(RSDP_T));
	datapoint->score = score;
	datapoint->group = group;
	datapoint->rank = -1.0;
	return datapoint;
}

/**********************************************************************
  allocate_ranksum_datapoint()

  Constructor for the ranksum data structure
**********************************************************************/
RSD_T *allocate_ranksum_data(){
	// Allocate memory and initialize fields
	RSD_T *data = mm_malloc(sizeof(RSD_T));

	data->size = 0;
	data->allocated_size = 0;
	data->data_points = NULL;
	return data;
}

/**********************************************************************
  allocate_and_init_ranksum_data()

  Constructor for the ranksum data structure which allocates memory for
  up to size datapoints
**********************************************************************/
RSD_T *allocate_and_init_ranksum_data(unsigned int size){
	// Allocate memory and initialize fields
	RSD_T *data = mm_malloc(sizeof(RSD_T));

	data->size = 0;
	data->allocated_size = size;
	data->data_points = (RSDP_T **) mm_malloc(data->allocated_size * sizeof(RSDP_T *));
	return data;
}

/**********************************************************************
  ranksum_from_stats()

  Computes the ranksum test for a given statistic.
**********************************************************************/
RSR_T *ranksum_from_stats(
   int n,			/* number of samples */
   int na,			/* number of positives */
   double ta_obs	/* sum of ranks of positives */
){
	RSR_T *mww = allocate_ranksum_result();

	int nb = n-na;

	double tab=n*(n+1)*0.5; 				// sum of n ranks in groups A and B combined
	double tb_obs = tab - ta_obs;

	double sd=sqrt((na*nb*(n+1.0))/12.0); 	// the standard deviation is the same in both sets
	double ta_null=na*(n+1.0)/2.0;			// the sum of the "null" case
	double tb_null=nb*(n+1.0)/2.0;			// the sum of the "null" case
	double ta_max=na*nb+(na*(na+1.0))/2.0;	// the max sum set A can take
	double tb_max=na*nb+(nb*(nb+1.0))/2.0;	// the max sum set B can take
	double ua=ta_max-ta_obs;				// the "U" value for A
	double ua_null=ta_max-ta_null;			// the U value for the null case
	double da=ta_obs>ta_null?-0.5:+0.5;		// a "continuity correction" for A
	double db=tb_obs>tb_null?-0.5:+0.5;		// a "continuity correction" for B
	double za=((ta_obs-ta_null)+da)/sd;		// the z value for A which is the mirror of ...
	double zb=((tb_obs-tb_null)+db)/sd;		// the z value for B (we only need one)
	double pa=0.5*(erfc(za/sqrt2));			// figure out the area of the normal distribution
	double pb=0.5*(erfc(zb/sqrt2));			// figure out the area of the normal distribution

	if (pa < pb){							// make sure we get the accuracy on
		mww->p_right = pa;					// the smaller of the p-values
		mww->p_left = 1.0-pa;				// due to numerical issues
		mww->p_twotailed = 2.0 * pa;
	} else {
		mww->p_left = pb;
		mww->p_right = 1.0-pb;
		mww->p_twotailed = 2.0 * pb;
	}
	mww->u = ua;
	return mww;
}

/*************************************************************************
 * ranksum_data_cmp()
 *
 * Compares two ranksum_data structures
 *************************************************************************/
int ranksum_data_cmp(const void *v1, const void *v2)
{
  assert(v1 != NULL);
  assert(v2 != NULL);
  RSDP_T* s1 = *(RSDP_T **)v1;
  RSDP_T* s2 = *(RSDP_T **)v2;
  double diff = s1->score - s2->score;
  if (diff < 0.0){
	  return 1;
  } else if (diff > 0.0){
  	  return -1;
  } else {
  	  return 0;
  }
}

/*************************************************************************
 * add_ranksum_datapoint()
 *
 * Adds a data point to the ranksum dataset
 *************************************************************************/
void add_ranksum_datapoint(
		RSD_T *data,
		RSDP_T* datapoint
){
  assert(data != NULL);
  if (datapoint == NULL) {
	  return;
  }
  assert(data->size <= data->allocated_size);
  if (data->size == data->allocated_size) {
	  data->allocated_size += DATA_INCREMENT;
	  data->data_points = mm_realloc(
			  data->data_points,
			  data->allocated_size * sizeof(RSDP_T *)
      );
  }
  data->data_points[data->size] = datapoint;
  data->size++;
}

/**********************************************************************
  get_ranksum_dataset2()

  Creates a ranksum input data set that can be (re-)used with the
  function run_ranksum_test() or run_ranksum_test_on_ordered_dataset()

  Copies the list of samples and groups to the new data structure.
  Calculates the (adjusted) rank of each score.

  NOTE: The samples are expected to be sorted already. The arrays
        of samples and groups have to be of the specified size and each
        entry in group corresponds to the sample at the same index.
**********************************************************************/
RSD_T *get_ranksum_dataset2(
	void* samples,  	/* the sample set, allows NULL */
  double (*get_sample)(const void*, int), /* sample accessor */
	void* group,	/* set of the group each samples belongs to, allows NULL */
  BOOLEAN_T (*get_group)(const void*, int), /* group accessor */
	int size 			/* the size of sample set  */
){
	assert(size > 0);
	int i;
	// Create data structure 
	RSD_T *dataset = allocate_and_init_ranksum_data(size);
	// Fill data structure
	for (i = 0; i < size; ++i) {
		RSDP_T* datapoint = allocate_ranksum_datapoint(get_sample(samples, i), get_group(group, i));
		add_ranksum_datapoint(dataset, datapoint);
	}
	// adjust rankings for same value elements
	adjust_ranking(dataset);
	return (dataset);
}

/*
 * Internal functions to get items from arrays of type double and boolean
 */
double double_array_get(const void *array, int index) {
  return ((double*)array)[index];
}
BOOLEAN_T boolean_array_get(const void *array, int index) {
  return ((BOOLEAN_T*)array)[index];
}

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
){
  assert(samples != NULL);
  assert(group != NULL);
	return get_ranksum_dataset2(samples, double_array_get, group, boolean_array_get, size);
}
/**********************************************************************
  ranksum_from_groups()

  Calculates the ranksum test for a set of samples which classes are
  specified in a corresponding array.
**********************************************************************/
RSR_T *ranksum_from_groups(
		double* samples,  	/* the sample set */
		BOOLEAN_T* group,	/* set of the group each samples belongs to */
		int size 			/* the size of sample set  */
){
	assert(samples != NULL);
	assert(group != NULL);
	assert(size > 0);

	int i;

	/* Create data structure */
	RSD_T *dataset = allocate_and_init_ranksum_data(size);

	/* Fill data structure */
	for (i=0;i<size;++i){
		RSDP_T* datapoint = allocate_ranksum_datapoint(samples[i], group[i]);
		add_ranksum_datapoint(dataset, datapoint);
	}

	/* calculate rank sum statistics*/
	RSR_T *mww = run_ranksum_test(dataset);

	destroy_rsd(dataset);
	return (mww);
}

/**********************************************************************
  ranksum_sample_sets()

  Calculates the ranksum test for two sets of samples and return a
  ranksum test result struct
**********************************************************************/
RSR_T *ranksum_from_sets(
		double* sample_a,  /* an array a of double values */
		int size_a, 		/* the size of sample set s */
		double* sample_b, 	/* an array b of double values */
		int size_b 			/* the size of sample set b */
){
	assert(sample_a != NULL);
	assert(sample_b != NULL);
	assert(size_a > 0);
	assert(size_b > 0);

	int i;

	/* Create data structure */
	RSD_T *dataset = allocate_ranksum_data();

	/* Fill data structure */
	for (i=0;i<size_a;++i){
		RSDP_T* datapoint = allocate_ranksum_datapoint(sample_a[i], TRUE);
		add_ranksum_datapoint(dataset, datapoint);
	}
	for (i=0;i<size_b;++i){
		RSDP_T* datapoint = allocate_ranksum_datapoint(sample_b[i], FALSE);
		add_ranksum_datapoint(dataset, datapoint);
	}

	/* calculate rank sum statistics*/
	RSR_T *mww = run_ranksum_test(dataset);
	destroy_rsd(dataset);
	return (mww);
}

/**********************************************************************
  run_ranksum_test()

  Runs the ranksum test on a given dataset. Use THIS function, if
  the dataset has not NOT ORDERED manually beforehand.
**********************************************************************/
RSR_T* run_ranksum_test(
		RSD_T *dataset 		// the dataset to work with
){
	// sort all entries according to the measurement
	qsort(dataset->data_points,dataset->size,sizeof(dataset->data_points[0]),ranksum_data_cmp);
	int i=0;
	adjust_ranking(dataset);
	RSR_T* mww = run_ranksum_test_on_ordered_dataset(dataset);
	return mww;
}

/**********************************************************************
  set_ranksum_groups()

  Assigns the given groups to the dataset structure. Both have to
  be of the same size.

  Most beneficial when used in conjuction with get_ranksum_dataset and
  run_ranksum_test_on_ordered_dataset i.e. when the ranking (ordering)
  does not change between the ranksum tests but the groups (classes) do.
**********************************************************************/
void set_ranksum_groups(RSD_T *dataset, BOOLEAN_T* groups){
	int i;
	for (i=0;i<dataset->size;++i){
		dataset->data_points[i]->group = groups[i];
	}
}

/**********************************************************************
  set_ranksum_group()

  Assigns the given position to the given group in the dataset structure. 
  For obvious reasons the position must exist.

  Most beneficial when used in conjuction with get_ranksum_dataset and
  run_ranksum_test_on_ordered_dataset i.e. when the ranking (ordering)
  does not change between the ranksum tests but the groups (classes) do.
**********************************************************************/
void set_ranksum_group(RSD_T *dataset, int position, BOOLEAN_T group){
  dataset->data_points[position]->group = group;
}

/********************************************************************
  get_ranksum_rank()

  Gets the rank of the position
**********************************************************************/
double get_ranksum_rank(RSD_T *dataset, int position) {
  return dataset->data_points[position]->rank;
}

/**********************************************************************
  adjust_ranking()

  "Re-ranks" the list according to same value elements by averaging
  their ranks
**********************************************************************/
void adjust_ranking(RSD_T *dataset){
	int first = 0,i=0, j=0, rank=0;
	double average;
	for (i=0;i<dataset->size;++i){
		if (dataset->data_points[first]->score != dataset->data_points[i]->score){ // data point differs from last one
			average = (rank - (first + 1))/2.0;
			for (j=first;j<i;++j){
				dataset->data_points[j]->rank = first + 1.0 + average;
			}
			first = i;
		}
		++rank;
	}
	// the last batch of entries is handled outside the loop
	average = (rank - (first + 1))/2.0;
	for (j=first;j<dataset->size;++j){
		dataset->data_points[j]->rank = first + 1.0 + average;
	}
}

/**********************************************************************
  run_ranksum_test()

  Runs the ranksum test on a given already ordered dataset
**********************************************************************/
RSR_T* run_ranksum_test_on_ordered_dataset(
		RSD_T *dataset 	// the ordered dataset to work with
){
	assert(dataset != NULL);

	int na=0; // number of samples in  class a
	double ta_obs = 0.0;	// sum of na ranks in group a
	int n = dataset->size;	// total number of measurements

	int i=0;
	for (i=0;i<dataset->size;++i){
		if (dataset->data_points[i]->group){
			++na;
			ta_obs += dataset->data_points[i]->rank; 	// sum the ranks
		}
	}

	RSR_T* mww = ranksum_from_stats(n,na,ta_obs);
	return mww;
}

/**********************************************************************
  RSR_get_p_left()

  returns the p-value for the left side (lesser)
**********************************************************************/
double RSR_get_p_left(RSR_T* ranksum_result){
	return ranksum_result->p_left;
}

/**********************************************************************
  RSR_get_p_right()

  returns the p-value for the right side (greater)
**********************************************************************/
double RSR_get_p_right(RSR_T* ranksum_result){
	return ranksum_result->p_right;
}

/**********************************************************************
  RSR_get_p_onetailed()

  returns the p-value for the one tailed side (the smaller one of left
  and right)
**********************************************************************/
double RSR_get_p_onetailed(RSR_T* ranksum_result){
	return min(ranksum_result->p_right,ranksum_result->p_left);
}

/**********************************************************************
  RSR_get_p_twotailed()

  returns the p-value for the twotailed side
**********************************************************************/
double RSR_get_p_twotailed(RSR_T* ranksum_result){
	return ranksum_result->p_twotailed;
}

/**********************************************************************
  RSR_get_u()

  returns the u value for the class considered
**********************************************************************/
double RSR_get_u(RSR_T*  ranksum_result){
	return ranksum_result->u;
}

/**********************************************************************
  free_rsd()

  Frees the memory, i.e. removes the pointers without deleting the
  underlying data points. Use destroy_rsd() if all memory should be freed
  instead.
**********************************************************************/
void free_rsd(RSD_T *dataset){
	myfree(dataset->data_points);
	dataset->size = 0;
	myfree(dataset);
}

/**********************************************************************
  destroy_rsd()

  Frees the memory, i.e. removes the pointers without deleting the
  final data points. Use destroy_rsd() if all memory should be freed
  instead.
**********************************************************************/
void destroy_rsd(RSD_T *dataset){
	while(dataset->size > 0) {
		RSDP_T *dp = dataset->data_points[--dataset->size] ;
		myfree(dp);
	}
	myfree(dataset->data_points);
	dataset->size = 0;

	myfree(dataset);
}

/**********************************************************************
  destroy_rsr()

  Frees the memory of a RSR_T* data structure
**********************************************************************/
void destroy_rsr(RSR_T *dataset){
	myfree(dataset);
}

/*****************************************************************************
 * MAIN
 *****************************************************************************/
#ifdef MAIN

#include "simple-getopt.h"
VERBOSE_T verbosity = INVALID_VERBOSE;
#define MAX_LINE 1000

int main
  (int    argc,
   char * argv[])
{

  const int P_LEFT = -1;
  const int P_ONETAILED = 0;
  const int P_RIGHT = 1;
  const int P_TWOTAILED = 2;
  const int P_ALL = 3;

  // Default parameter settings.
  verbosity = NORMAL_VERBOSE;
  int requested_p_value = P_ALL;

  const int num_options = 2;
  cmdoption const options[] = {
    { "verbosity", REQUIRED_VALUE },
    { "p-value", REQUIRED_VALUE }
  };

  // Define the usage message.
  char      usage[400] = "";
  strcat(usage, "USAGE: ranksum_test [options] <n> <p> <r>\n");
  strcat(usage, "\n");
  strcat(usage, "   <n> number of samples \n");
  strcat(usage, "   <p> number of positives\n");
  strcat(usage, "   <r> ranksum of positives (may be a real number)\n");
  strcat(usage, "\n");
  strcat(usage, "   Options:\n");
  strcat(usage, "     --p-value -1|0|1|2|3 (-1=left, 0=one-tailed,1=right,\n"
			    "                          	 2=two-tailed, 3=all (default))\n");
  strcat(usage, "     --verbosity 1|2|3|4 (default = 2)\n");
  strcat(usage, "\n");

  // Parse the command line.
  int option_index = 0;
  char* option_name = NULL;
  char* option_value = NULL;
  const char *  message = NULL;
  char line[MAX_LINE];    // Buffer for reading.
  simple_setopt(argc, argv, num_options, options);
  while(1) {
    // Read the next option, and break if we're done.
    int c = simple_getopt(&option_name, &option_value, &option_index);
    if (c == 0) {
      break;
    } else if (c < 0) {
      simple_getopterror(&message);
      die("Error processing command line options (%s)\n", message);
    }

    if (strcmp(option_name, "p-value") == 0) {
      requested_p_value = atoi(option_value);
      if (requested_p_value < -1 || requested_p_value > 3)
    	  die("Error requested p-value parameter unknown (%d)\n", requested_p_value);
	}
    else if (strcmp(option_name, "verbosity") == 0) {
      verbosity = atoi(option_value);
    }
  }

  // Read the single required argument.
  if (option_index + 3 != argc) {
    fprintf(stderr, "%s", usage);
    exit(1);
  }
  int n = atoi(argv[option_index]);option_index++;
  int p = atoi(argv[option_index]);option_index++;
  double r = atof(argv[option_index]);option_index++;

  RSR_T* mww = ranksum_from_stats(n,p,r);

  // Print to stdout.
  if (requested_p_value == P_RIGHT){
	  printf("%g\n",RSR_get_p_right(mww));
  } else if (requested_p_value == P_LEFT){
	  printf("%g\n",RSR_get_p_left(mww));
  } else if (requested_p_value == P_ONETAILED){
	  printf("%g\n",RSR_get_p_onetailed(mww));
  } else if (requested_p_value == P_TWOTAILED){
	  printf("%g\n",RSR_get_p_twotailed(mww));
  } else {
	  printf("p-left\t%g\tp-right\t%g\tone-tailed\t%g\ttwo-tailed\t%g\tU-value\t%g\n",
			  RSR_get_p_left(mww),
			  RSR_get_p_right(mww),
			  RSR_get_p_onetailed(mww),
			  RSR_get_p_twotailed(mww),
			  RSR_get_u(mww));
  }
  return(0);
}

#endif
