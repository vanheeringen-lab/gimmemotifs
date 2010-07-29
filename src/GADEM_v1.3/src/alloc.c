#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include "gadem.h"

double **alloc_double_double(int size1,int size2) {

   double **tmp=NULL;
   register int i;

   tmp=(double **)calloc(size1,sizeof(double *));
   if (!tmp)   { printf("tmp calloc failed!\n"); exit(1); }

   tmp[0]=(double *)calloc(size1*size2,sizeof(double));
   if (tmp[0]==0) { printf("bit calloc failed!\n"); exit(1); }

   /* set up vector pointers */
   for (i=1; i<size1; i++) tmp[i]=tmp[0]+(size2 * i);

   return (tmp);
}

char **alloc_char_char(int size1,int size2) {

   char **tmp=NULL;
   register int i;

   tmp=(char **)calloc(size1,sizeof(char *));
   if (!tmp)   { printf("tmp calloc failed!\n"); exit(1); }

   tmp[0]=(char *)calloc(size1*size2,sizeof(char));
   if (tmp[0]==0) { printf("bit calloc failed!\n"); exit(1); }

   /* set up vector pointers */
   for (i=1; i<size1; i++)  tmp[i]=tmp[0]+(size2 * i);

   return (tmp);
}

int **alloc_int_int(int size1,int size2) {

   int **tmp=NULL;
   register int i;

   tmp=(int **)calloc(size1,sizeof(int *));
   if (!tmp)   { printf("tmp calloc failed!\n"); exit(1); }

   tmp[0]=(int *)calloc(size1*size2,sizeof(int));
   if (tmp[0]==0) { printf("bit calloc failed!\n"); exit(1); }

   /* set up vector pointers */
   for (i=1; i<size1; i++)  tmp[i]=tmp[0]+(size2 * i);

   return (tmp);
}

char *alloc_char(int size1) {

   char *tmp=NULL;
   tmp=(char *)calloc(size1,sizeof(char));
   if (!tmp)   { printf("tmp calloc failed!\n"); exit(1); }
   return (tmp);
}

int *alloc_int(int size1) {

   int *tmp=NULL;
   tmp=(int *)calloc(size1,sizeof(int));
   if (!tmp)   { printf("tmp calloc failed!\n"); exit(1); }
   return (tmp);
}

double *alloc_double(int size1) {

   double *tmp=NULL;
   tmp=(double *)calloc(size1,sizeof(double));
   if (!tmp)   { printf("tmp calloc failed!\n"); exit(1); }
   return (tmp);
}

double ***alloc_double_double_double(int num_niche,int pop_size,int size) {

   double ***tmp=NULL;
   register int i,j;

   tmp      =(double ***)calloc(num_niche,sizeof(double **));
   tmp[0]   =(double **) calloc(num_niche*pop_size,sizeof(double *));
   tmp[0][0]=(double *)  calloc(num_niche*pop_size*size,sizeof(double));

   for (i=1; i<num_niche; i++) tmp[i]=tmp[0]+i*pop_size;

   for (j=1; j<pop_size; j++) tmp[0][j]=tmp[0][0]+size*j;

   for (i=1; i<num_niche; i++) {
      tmp[i][0]=tmp[0][0]+size*pop_size*i;
      for (j=1; j<pop_size; j++)  tmp[i][j]=tmp[i][0]+size*j;
   }
   return (tmp);
}

int ***alloc_int_int_int(int num_niche,int pop_size,int size) {

   int ***tmp=NULL;
   register int i,j;

   tmp      =(int ***)calloc(num_niche,sizeof(int **));
   tmp[0]   =(int **) calloc(num_niche*pop_size,sizeof(int *));
   tmp[0][0]=(int *)  calloc(num_niche*pop_size*size,sizeof(int));

   for (i=1; i<num_niche; i++) tmp[i]=tmp[0]+i*pop_size;

   for (j=1; j<pop_size; j++) tmp[0][j]=tmp[0][0]+size*j;

   for (i=1; i<num_niche; i++) {
      tmp[i][0]=tmp[0][0]+size*pop_size*i;
      for (j=1; j<pop_size; j++)  tmp[i][j]=tmp[i][0]+size*j;
   }
   return (tmp);
}

Fitness *alloc_fitness(int size1) {

   Fitness *tmp=NULL;

   tmp=(Fitness *)calloc(size1,sizeof(Fitness));
   if (!tmp) { printf("calloc for fitness failed!\n"); exit(0); }
   return (tmp);
}

Wheel *alloc_wheel(int size1) {

   Wheel *tmp=NULL;

   tmp=(Wheel *)calloc(size1,sizeof(Wheel));
   if (!tmp) { printf("calloc for wheel failed!\n"); exit(0); }
   return (tmp);
}

Sites *alloc_site(int size1) {

   Sites *tmp=NULL;

   tmp=(Sites *)calloc(size1,sizeof(Sites));
   if (!tmp) { printf("calloc for wheel failed!\n"); exit(0); }
   return (tmp);
}

Sites **alloc_site_site(int size1,int size2) {

   Sites **tmp=NULL;
   register int i;

   tmp=(Sites **)calloc(size1,sizeof(Sites *));
   if (!tmp)   { printf("tmp calloc failed!\n"); exit(1); }

   tmp[0]=(Sites *)calloc(size1*size2,sizeof(Sites));
   if (tmp[0]==0) { printf("bit calloc failed!\n"); exit(1); }

   /* set up vector pointers */
   for (i=1; i<size1; i++)  tmp[i]=tmp[0]+(size2 * i);

   return (tmp);
}

Chrs **alloc_chrs(int size1,int size2) {

   Chrs **tmp=NULL;
   register int i;

   tmp=(Chrs **)calloc(size1,sizeof(Chrs *));
   if (!tmp)   { printf("tmp calloc failed!\n"); exit(1); }

   tmp[0]=(Chrs *)calloc(size1*size2,sizeof(Chrs));
   if (tmp[0]==0) { printf("bit calloc failed!\n"); exit(1); }

   /* set up vector pointers */
   for (i=1; i<size1; i++)  tmp[i]=tmp[0]+(size2 * i);

   return (tmp);
}

Words *alloc_word(int numCategory,int maxSize) {

   int i;
   Words *tmp=NULL;
   
   tmp=(Words *)calloc(numCategory,sizeof(Words));
   if (!tmp) { printf("calloc failed for Words.\n"); exit(0); }

   for (i=0; i<numCategory; i++) {
      tmp[i].s1=alloc_char_char(maxSize,10); 
      tmp[i].prob_sta=alloc_double(maxSize); 
      tmp[i].prob_end=alloc_double(maxSize); 
   }
   return (tmp);
}

void destroy_word(Words *word,int numCategory) {

   int i;
  
   for (i=0; i<numCategory; i++) {
      if (word[i].prob_sta) { free(word[i].prob_sta); word[i].prob_sta=NULL;  }
      if (word[i].prob_end) { free(word[i].prob_end); word[i].prob_end=NULL;  }
      if (word[i].s1[0])    { free(word[i].s1[0]);    word[i].s1[0]=NULL;     }
      if (word[i].s1)       { free(word[i].s1);       word[i].s1=NULL;        }
   }
   if (word) { free(word); word=NULL; }
}

Ktuples *alloc_ktuples(int numKmer,int kmerLen) {

   int i;
   Ktuples *tmp=NULL;

   tmp=(Ktuples *)calloc(numKmer,sizeof(Ktuples));

   for (i=0; i<numKmer; i++) tmp[i].seq=alloc_char(kmerLen+1);
   return (tmp);   
}

Pgfs *alloc_distr(int size) {

   Pgfs *llDist=NULL;
   llDist=(Pgfs *)calloc(size,sizeof(Pgfs));
   if (!llDist) { printf("calloc failed for llDist\n"); exit(0); }
   return (llDist);
}

void destroy_ktuples(Ktuples *ktuple,int numKmer) {

   int i;
   for (i=0; i<numKmer; i++) {
      if (ktuple[i].seq) { free(ktuple[i].seq); ktuple[i].seq=NULL; }
   }
   if (ktuple) { free(ktuple); ktuple=NULL; }
}
