/*
 * $Id: spearman-rank-correlation.h $
 * 
 * $Log$
 *
 */

/***********************************************************************
*                                                                      *
*       MEME++                                                         *
*       Copyright 1995, The Regents of the University of California    *
*       Author: Timothy L. Bailey                                      *
*                                                                      *
***********************************************************************/
 
#ifndef SPEARMANRC_H
#define SPEARMANRC_H

typedef struct {
	double data;
	double rank;
	int orig_rank;
} spearman_rank_t;

double spearman_rank_correlation(
  int n,			/* number of points */
  double *x,			/* x values */
  double *y			/* y values */
);

#endif
 
