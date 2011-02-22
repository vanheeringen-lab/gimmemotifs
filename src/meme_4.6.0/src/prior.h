/*
 * $Id: prior.h 1048 2006-07-06 20:07:44Z cegrant $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 18:46:42  nadya
 * Initial revision
 *
 */

/* This are the data structures used in the adaptive prior method */

#include "user.h"
#define	Real		double	
#define RealPrec	double

/* This structure stores parameters of the prior distributions */
typedef struct{
	char	alphabet[MAXALPH+1];	/* The alphabet, in order */
	int	AlphaChar;	/* Number of alphabet characters */
	int	L;		/* Number of prior distributions */
	Real	*Mix;		/* Mixture coefficents for each prior */
	Real    *B;		/* strength of each prior */
	Real	**Distr;	/* Prior distributions. L Dirchlet's 
				   over AlphaChar positions:
				   Distribution[L][AlphaChar+1] */
	int	*FullUpdate;	/* !=0 re-estimate all alpha
				   ==0 re-estimate alpha 0 only */
	int	*QUpdate;	/* !=0 update mixture coefficents
				   ==0 do not update coefficents */
	char	**StructID;	/* Structure Tag */
	char	**Comment;
} PriorLib;

/* Allocate space for a prior library with L priors and Alpha characters */
PriorLib	*alloc_PriorLib( int L, int Alpha );
void free_PriorLib (PriorLib *lib);

/* This reads prior information from a file into a PriorLib */
PriorLib	*read_PriorLib( char *plib_name, double desired_beta );

/* This calculates the regularizer given the observed
   frequencies, the prior library, and a weight of the priors */
extern void mixture_regularizer(
  double *freq,			/* obs freq */
  PriorLib *Lib,		/* priors */ 
  double *reg			/* pseudo-counts */
);

Real logpajgy( Real *y, PriorLib *Lib, int j, int Calc);

RealPrec logpygaj( Real *y, Real *a, int AlphLength);

