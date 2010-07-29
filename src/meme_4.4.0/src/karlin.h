/*****************************************************************************
 * FILE: karlin.h
 * AUTHOR: Timothy Bailey (adapted from NCBI file karlin.h)
 * CREATE DATE: 4/23/2002
 * PROJECT: shared
 * COPYRIGHT: 2002, University of Queensland
 * VERSION: $Revision: 1.1.1.1 $
 * DESCRIPTION: Distribution parameters for scoring matrices
 *****************************************************************************/

#ifndef KARLIN_H
#define KARLIN_H

#include <math.h>
/*#include "stdinc.h"*/
#include "alphabet.h"
#include "array.h"

#define NEW(x,n,t)      (( (x = calloc(n,sizeof(t))) == NULL) ? \
                           (fprintf(stderr,"Out of Memory."), exit(1), NULL) : x)
#define valAlphaR(c,d,A)        ((A)->R[(c)][(d)])


#define Boolean char
/************************** alphabet datatype *****************************/
typedef struct {
        long     n;                     /* number of LETTERS */
        char    *alphabet;              /* ALPHABET */
        char    *code2let;              /* CODE2LETTER */
        char    *code2lower;            /* CODE2LETTER lower case */
        char    *let2code;              /* LETTER2CODE */
        char    **R;                    /* relatedness scoring matrix */
        char    *C;                     /* complementary bases */
        long    loR;                    /* lowest value in R */
        long    hiR;                    /* highest value in R */
        char    *prs;                   /* pairs string */
        char    **pairs;                /* residue pairs */
        Boolean *paired;                /* is residue r paired? */
        long    npairs;                 /* number of pairs */
} alphabet_type;
typedef alphabet_type *a_type;

/***************************************************************************
 * Define the KARLIN_INPUT type.
 ***************************************************************************/
typedef struct karlin_input_t {
  long low;			/* Lowest score (must be negative) */
  long high;			/* Highest score (must be positive) */
  double escore;		/* Expected score */
  ARRAY_T *prob;		/* Probabilities for all scores */
} KARLIN_INPUT_T;

#define KARLINMAXIT 50  /* Max. # iterations used in calculating K */
Boolean	karlin(long low,long high,double *pr,double *lambda,double *K,double *H);
double  ExpectedInformation(a_type A, double lambda, double *freq);
long 	karlin_gcd(long a,long b);

#endif

