
#include <stdlib.h>
#include "gadem.h"

/* sort by increasing order */
void sort_fitness(Fitness *fitness,int size) {

   int (*compar)(const void *,const void *);

   compar=Compare_fitness;
   qsort((void *)fitness,(size_t)size,sizeof(Fitness),compar);
}

/* sort by increasing order */
int Compare_fitness(const void *s1, const void *s2) {

   if (((Fitness *)s1)->value < ((Fitness *)s2)->value) { return -1; }
   if (((Fitness *)s1)->value > ((Fitness *)s2)->value) { return  1; }
      return 0;
}

/* descending */
void sort_double(double *data,int size) {

   int (*compar)(const void *,const void *);

   compar = Compare_double;
   qsort((void *)data,(size_t)size,sizeof(double),compar);
}

int Compare_double(const void *s1, const void *s2) {

   if (*((double *)s1)< (*(double *)s2)) { return 1; }
   if (*((double *)s1)> (*(double *)s2)) { return -1; }
      return 0;
}

void sort_llrDist(Pgfs *M1,int size) {

   int (*compar)(const void *,const void *);

   compar=Compare_M1;
   qsort((void *)M1,(size_t)size,sizeof(Pgfs),compar);
}

// sort in decending order
int Compare_M1(const void *s1, const void *s2) {

   if (((Pgfs *)s1)->score < ((Pgfs *)s2)->score) { return  1; }
   if (((Pgfs *)s1)->score > ((Pgfs *)s2)->score) { return -1; }
      return 0;
}

void sort_kmer_z(Ktuples *s2,int size) {

   int (*compar)(const void *,const void *);
   compar=Compare_z;
   qsort((void *)s2,(size_t)size,sizeof(Ktuples),compar);
}

int Compare_z(const void *s1,const void *s2) {

   if (((Ktuples *)s1)->z < ((Ktuples *)s2)->z) { return  1; }
   if (((Ktuples *)s1)->z > ((Ktuples *)s2)->z) { return -1; }
      return 0;
}

void sort_sites_score(Sites *site,int size) {

   int (*compar)(const void *,const void *);

   compar=Compare_score;
   qsort((void *)site,(size_t)size,sizeof(Sites),compar);
}

// sort in descending order
int Compare_score(const void *s1, const void *s2) {

   if (((Sites *)s1)->score < ((Sites *)s2)->score) { return  1; }
   if (((Sites *)s1)->score > ((Sites *)s2)->score) { return -1; }
      return 0;
}

void sort_sites_llr(Sites *site,int size) {

   int (*compar)(const void *,const void *);

   compar=Compare_llr;
   qsort((void *)site,(size_t)size,sizeof(Sites),compar);
}

// sort in descending order
int Compare_llr(const void *s1, const void *s2) {

   if (((Sites *)s1)->llr < ((Sites *)s2)->llr) { return  1; }
   if (((Sites *)s1)->llr > ((Sites *)s2)->llr) { return -1; }
      return 0;
}

/* sort by decreasing order */
void sort_vector(ColCons *data,int size) {

   int (*compar)(const void *,const void *);

   compar=Compare_vector;
   qsort((void *)data,(size_t)size,sizeof(ColCons),compar);
}

/* sort by decreasing order */
int Compare_vector(const void *s1, const void *s2) {

   if (((ColCons *)s1)->v < ((ColCons *)s2)->v) { return  1; }
   if (((ColCons *)s1)->v > ((ColCons *)s2)->v) { return -1; }
      return 0;
}

