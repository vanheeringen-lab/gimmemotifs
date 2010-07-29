/*-----------------------------------------------------------------
        Random number generator by George Marsaglia algorithm

              From AMBER subroutines: amrset and amrand
   
   Testing:  Call amrset with iseed = 54185253.  This should result
   in is1 = 1802 and is2 = 9373.  Call amrand 20000 times, then six
   more times, printing the six random numbers * 2**24 (ie, 4096*4096)
   They should be: (6f12.1)
   6533892.0  14220222.0  7275067.0  6172232.0  8354498.0  10633180.0
 
   This random number generator originally appeared in "Toward a Universal
   Random Number Generator" by George Marsaglia and Arif Zaman.  Florida
   State University Report: FSU-SCRI-87-50 (1987)
 
   It was later modified by F. James and published in "A Review of Pseudo-
   random Number Generators"
 
   This is claimed to be the best known random number generator available.
   It passes ALL of the tests for random number generators and has a
   period of 2^144, is completely portable (gives bit identical results on
   all machines with at least 24-bit mantissas in the floating point
   representation).

   The algorithm is a combination of a Fibonacci sequence (with lags of 97
   and 33, and operation "subtraction plus one, modulo one") and an
   "arithmetic sequence" (using subtraction).


                         Date: March 26, 1999
            Last Modification: July 11, 2001
-------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>

#define  min(a,b) (((a)<(b))?(a):(b))
#define  max(a,b) (((a)>(b))?(a):(b))
#define  mod(a,b) ((a)-(b)*((a)/(b)))

#define IS1MAX 31328
#define IS2MAX 30081

static int i97,j97;
static double u[98],c,cd,cm;

void sgenrand_32bits(int seed) {

   register int num1,num2;
   register int i,j,k,l,m;
   double s, t;
   register int ii,jj;

   num1 = max((seed / IS2MAX)+1, 1);
   num1 = min(num1, IS1MAX);
   num2 = max(1, mod(seed, IS2MAX)+1);
   num2 = min(num2, IS2MAX);

   i = mod(num1/177, 177) + 2;
   j = mod(num1    , 177) + 2;
   k = mod(num2/169, 178) + 1;
   l = mod(num2    , 169);

   for( ii=1; ii< 98; ii++) {
       s = 0.0;
       t = 0.5;
       for ( jj = 1; jj< 25; jj++) {
          m = mod(mod(i*j, 179)*k, 179);
          i = j;
          j = k;
          k = m;
          l = mod(53*l+1, 169);
          if (mod(l*m, 64) >= 32)
             s += t;
          t = 0.5 * t;
       }
       *(u+ii) = s;
    }
    c  = 362436.0   / 16777216.0;
    cd = 7654321.0  / 16777216.0;
    cm = 16777213.0 / 16777216.0;

    i97 = 97;
    j97 = 33;
}

double genrand_32bits(void) {

   double uni;

   uni = *(u+i97) - *(u+j97);
   if ( uni < 0.0 ) uni += 1.0;
   u[i97] = uni;
   i97 -= 1;
   if (i97 == 0) i97 = 97;
   j97 -= 1;
   if (j97 == 0) j97 = 97;
   c -= cd;
   if ( c < 0.0 ) c += cm;
   uni -= c;
   if ( uni < 0.0 ) uni += 1.0;
   return (uni);
}

