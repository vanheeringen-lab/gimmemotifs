#if !defined(KARLIN)
#define KARLIN
#include <math.h>
#include "stdinc.h"
#include "alphabet.h"

#define KARLINMAXIT 50  /* Max. # iterations used in calculating K */
Boolean	karlin(long low,long high,double *pr,double *lambda,double *K,double *H);
double ExpectedInformation(a_type A, double lambda, double *freq);
long 	karlin_gcd(long a,long b);

#endif

