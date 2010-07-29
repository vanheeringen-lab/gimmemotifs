
#if !defined(RANDOM)
#define RANDOM
#include <stdio.h>
#include <math.h>

/* Random no. seeder and generator */
void sRandom(long x);
long Random(void);

#define DIM(A) (sizeof(A)/sizeof((A)[0]))

#endif 

