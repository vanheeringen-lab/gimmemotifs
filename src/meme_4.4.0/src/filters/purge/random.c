#include "random.h"

/*
	Additive random number generator

	Modelled after "Algorithm A" in
	Knuth, D. E. (1981). The art of computer programming, volume 2, page 27.

	7/26/90 WRG
*/

static long	state[33] = {
	(long)0xd53f1852,  (long)0xdfc78b83,  (long)0x4f256096,  (long)0xe643df7,
	(long)0x82c359bf,  (long)0xc7794dfa,  (long)0xd5e9ffaa,  (long)0x2c8cb64a,
	(long)0x2f07b334,  (long)0xad5a7eb5,  (long)0x96dc0cde,  (long)0x6fc24589,
	(long)0xa5853646,  (long)0xe71576e2,  (long)0xdae30df,  (long)0xb09ce711,
	(long)0x5e56ef87,  (long)0x4b4b0082,  (long)0x6f4f340e,  (long)0xc5bb17e8,
	(long)0xd788d765,  (long)0x67498087,  (long)0x9d7aba26,  (long)0x261351d4,
	(long)0x411ee7ea,  (long)0x393a263,  (long)0x2c5a5835,  (long)0xc115fcd8,
	(long)0x25e9132c,  (long)0xd0c6e906,  (long)0xc2bc5b2d,  (long)0x6c065c98,
	(long)0x6e37bd55 };

#define r_off	12

static long	*rJ = &state[r_off], *rK = &state[DIM(state)-1];

void sRandom(long x)
{
	register long	i;

	state[0] = x;
	/* linear congruential initializer */
	for (i=1; i<DIM(state); ++i) {
		state[i] = 1103515245*state[i-1] + 12345;
	}

	rJ = &state[r_off];
	rK = &state[DIM(state)-1];

	for (i=0; i<10*DIM(state); ++i) (void) Random();
}


/*
	Random --  return value in the range 0 <= x <= 2**31 - 1
*/
long Random(void)
{
	register long	r;

	r = *rK;
	r += *rJ--;
	*rK-- = r;

	if (rK < state)
		rK = &state[DIM(state)-1];
	else
		if (rJ < state)
			rJ = &state[DIM(state)-1];
	return (r>>1)&0x7fffffff; /* discard the least-random bit */
}


