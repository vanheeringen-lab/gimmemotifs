/* 

 64 Bit Mersenne Twister Random Number Generator 
 
 Reference:
 
 Nishimura, Takuji, "Tables of 64-bit Mersenne Twisters," ACM Transactions 
 on Modeling and Computer Simulation, Vol. 10, No. 4, October 2000, 
 ppg. 348-357.
 
*/


#define NN 312
#define M0 63
#define M1 151
#define M2 224

/* Constant vector a */
#define MATRIX_A 0xB3815B624FC82E2FULL

/* Most significant 33 bits */
#define UMASK    0xFFFFFFFF80000000ULL

/* Lease significant 31 bits */
#define LMASK    0x7FFFFFFFULL

/* Tempering parameters */
#define MASK_B 0x599CFCBFCA660000ULL
#define MASK_C 0xFFFAAFFE00000000ULL
#define UU 26
#define SS 17
#define TT 33
#define LL 39

/* The array for the state vector */
static unsigned long long mt[NN];

/* mti==NN+1 means mt[NN] is not initialized */
static int mti=NN+1;


void sgenrand_64bits(unsigned long long seed)
{
   unsigned long long ux, lx;
   
   for (mti = 0 ; mti < NN; mti++) {
      ux = seed & 0xFFFFFFFF00000000ULL;
      seed = 2862933555777941757ULL * seed + 1ULL;
      lx = seed >> 32;
      seed = 2862933555777941757ULL * seed + 1ULL;
      mt[mti] = ux | lx;
   }
}



double genrand_64bits(void) {

   int i;
   unsigned long long x;
   static unsigned long long mag01[2] = {0ULL, MATRIX_A};
   
   if (mti >= NN) {   /* Generate NN words at one time */
      
      /* If sgenrand() has not been called, a default initial seed is used */
      
      if (mti == NN+1) sgenrand_64bits(1ULL);
      
      for (i = 0 ; i < NN-M2; i++) {
         x = (mt[i] & UMASK) | (mt[i+1] & LMASK);
         mt[i]  = (x >> 1) ^ mag01[(int)(x & 1ULL)];
         mt[i] ^= mt[i+M0] ^ mt[i+M1] ^ mt[i+M2];
      }
      
      for (; i < NN-M1; i++) {
         x = (mt[i] & UMASK) | (mt[i+1] & LMASK);
         mt[i]  = (x >> 1) ^ mag01[(int)(x & 1ULL)];
         mt[i] ^= mt[i+M0] ^ mt[i+M1] ^ mt[i+M2-NN];
      }
      
      for (; i < NN-M0; i++) {
         x = (mt[i] & UMASK) | (mt[i+1] & LMASK);
         mt[i]  = (x >> 1) ^ mag01[(int)(x & 1ULL)];
         mt[i] ^= mt[i+M0] ^ mt[i+M1-NN] ^ mt[i+M2-NN];
      }
      
      for (; i < NN-1; i++) {
         x = (mt[i] & UMASK) | (mt[i+1] & LMASK);
         mt[i]  = (x >> 1) ^ mag01[(int)(x & 1ULL)];
         mt[i] ^= mt[i+M0-NN] ^ mt[i+M1-NN] ^ mt[i+M2-NN];
      }
      
      x = (mt[NN - 1] & UMASK) | (mt[0] & LMASK);
      mt[NN-1]  = (x >> 1) ^ mag01[(int)(x & 1ULL)];
      mt[NN-1] ^= mt[M0-1] ^ mt[M1-1] ^ mt[M2-1];
      
      mti = 0;
   }
   
   x = mt[mti++];
   x ^= (x >> UU);
   x ^= (x << SS) & MASK_B;
   x ^= (x << TT) & MASK_C;
   x ^= (x >> LL);
   
   return ((double) x / (double) 0xFFFFFFFFFFFFFFFFULL);
}


      
         
      
      
