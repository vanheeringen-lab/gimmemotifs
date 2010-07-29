/* 

Random number generation based on Mersenne twister
using implementations by Takuji Nishimura

*/

#ifndef __RANDOM_H__
#define __RANDOM_H__

#include "config.h"

#if SIZEOF_LONG == 8

   typedef unsigned long long seed_t;

   void sgenrand_64bits(unsigned long long seed);
   double genrand_64bits(void);
   
   #define sgenrand(x) sgenrand_64bits(x)
   #define genrand() genrand_64bits()
#else
   typedef unsigned long seed_t;

   void sgenrand_32bits(unsigned long seed);
   double genrand_32bits(void);
   
   #define sgenrand(x) sgenrand_32bits(x)
   #define genrand() genrand_32bits() 
#endif


#endif






