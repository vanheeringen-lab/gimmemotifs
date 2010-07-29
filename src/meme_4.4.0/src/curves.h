/*
 * $Id: curves.h 1048 2006-07-06 20:07:44Z cegrant $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 18:36:29  nadya
 * Initial revision
 *
 */

/***********************************************************************
*                                                                      *
*       MEME++                                                         *
*       Copyright 1995, The Regents of the University of California    *
*       Author: Timothy L. Bailey                                      *
*                                                                      *
***********************************************************************/

#ifndef CURVES_H
#define CURVES_H

#define 	MAXPARMS	10

/* the type of the curve function */
typedef char (*FTYPE( 	float	X,	float	A[],
		float *YFit,	float	dYdA[],
		int mode));

extern FTYPE *get_curve_type(
  int mode,				/* type of curve */
  int *nparms				/* number of parameters */
);

extern char *negexp( 	float	X,	float	A[],
		float *YFit,	float	dYdA[],
		int mode);

extern char *logcurve(	float	X,	float	A[],
		float *YFit,	float	dYdA[],
		int mode);

extern char *logcurve2(	float	X,	float	A[],
		float *YFit,	float	dYdA[],
		int mode);

extern char *logexp( 	float	X,	float	A[],
		float *YFit,	float	dYdA[],
		int mode);

extern char *loglog(	float	X,	float	A[],
		float *YFit,	float	dYdA[],
		int mode);

extern char *negexp2( 	float	X,	float	A[],
		float *YFit,	float	dYdA[],
		int mode);

extern char *invsqrt(	float	X,	float	A[],
		float *YFit,	float	dYdA[],
		int mode);

extern char *logcurve3(	float	X,	float	A[],
		float *YFit,	float	dYdA[],
		int mode);

extern char *logcurve4(	float	X,	float	A[],
		float *YFit,	float	dYdA[],
		int mode);

extern char *gcurve(	float	X,	float	A[],
		float *YFit,	float	dYdA[],
		int mode);

extern char *rootlog(	float	X,	float	A[],
		float *YFit,	float	dYdA[],
		int mode);

extern char *gsdcurve(	float	X,	float	A[],
		float *YFit,	float	dYdA[],
		int mode);
#endif

