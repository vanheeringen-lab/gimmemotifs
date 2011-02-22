/*
 * $Id: spearman-rank-correlation.c $
 * 
 * $Log$
 *
 */

#include <math.h>

#include "macros.h"
#include "spearman-rank-correlation.h"

static int compare_spearman_rank_t_data(const void *a, const void *b);
static int compare_spearman_rank_t_orig_rank(const void *a, const void *b);

/*
 * 	spearman_rank_correlation:
 *		Returns the product moment correlation coefficient of the ranked input.
 */
double spearman_rank_correlation(
				int n,			/* number of points */
				double *x,			/* x values */
				double *y			/* y values */
				)
{
	int i;
	spearman_rank_t *sx, *sy; //Sorted-x and sorted-y in our struct type.
	double sumx, sumxsq, sumxy, sumy, sumysq;
	sumx = sumxsq = sumxy = sumy = sumysq = 0;

	sx = malloc(n*sizeof(spearman_rank_t));
	sy = malloc(n*sizeof(spearman_rank_t));

	//Copy into the array of structs and the pointer array
	for (i=0;i<n;i++) {
		sx[i].data = x[i];
		sy[i].data = y[i];
		sx[i].orig_rank = i;
		sy[i].orig_rank = i;
	}

	qsort(sx, n, sizeof(*sx), compare_spearman_rank_t_data);
	qsort(sy, n, sizeof(*sy), compare_spearman_rank_t_data);

	// Assign the ranks - we need to do this once per array due to the
	// fact that we're accounting for ties and may skip some.
	for (i=0;i<n;i++) {

		double rank = 0;
		int num_ties = 0;
		int j;

		// Look-ahead for ties.
		// This will include itself as a tie, deliberately
		// as a sanity check to make sure that it's working.
		for (j=i; j<n && sx[i].data == sx[j].data; j++) {
			rank += j+1; //ranks start at 1
			num_ties++;
		}
		//Set the ranks to the mean of the tied ranks.
		for (j=i; j<i+num_ties; j++) {
			sx[j].rank = rank / num_ties;
		}

		i += num_ties - 1; //skip those we've already assigned ranks to
	}

	for (i=0;i<n;i++) {

		double rank = 0;
		int num_ties = 0;
		int j;

		// Look-ahead for ties.
		// This will include itself as a tie, deliberately
		// as a sanity check to make sure that it's working.
		for (j=i; j<n && sy[i].data == sy[j].data; j++) {
			rank += j+1; //ranks start at 1
			num_ties++;
		}
		//Set the ranks to the mean of the tied ranks.
		for (j=i; j<i+num_ties; j++) {
			sy[j].rank = rank / num_ties;
		}

		i += num_ties - 1; //skip those we've already assigned ranks to
	}

	//Resort based on original ranks
	qsort(sx, n, sizeof(*sx), compare_spearman_rank_t_orig_rank);
	qsort(sy, n, sizeof(*sy), compare_spearman_rank_t_orig_rank);
	

	//Calculate the various components of the product moment correlation equ.
	for (i=0;i<n;i++) {
		sumxy += sx[i].rank*sy[i].rank;
		sumx += sx[i].rank;
		sumy += sy[i].rank;
		sumxsq += sx[i].rank*sx[i].rank;
		sumysq += sy[i].rank*sy[i].rank;
	}

	/*
	 * Complete the score:
	 */

	free(sx);
	free(sy);

	return (n*sumxy - sumx*sumy) / (sqrt(n*sumxsq - sumx*sumx) * sqrt(n*sumysq - sumy*sumy));

}

static int compare_spearman_rank_t_data(const void *a, const void *b)
{
	double temp = ((spearman_rank_t *) a)->data - ((spearman_rank_t *) b)->data;
	if (temp > 0)
		return 1;
	else if (temp < 0)
		return -1;
	else
		return 0;
}

static int compare_spearman_rank_t_orig_rank(const void *a, const void *b)
{
	double temp = ((spearman_rank_t *) a)->orig_rank - ((spearman_rank_t *) b)->orig_rank;
	if (temp > 0)
		return 1;
	else if (temp < 0)
		return -1;
	else
		return 0;
}
