#include "config.h"

#include <stdlib.h>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif 

#include <math.h>
#include <stdio.h>
#include "gadem.h"
#include "defines.h"
#include "random.h"

int Compare_double2(const void *s1, const void *s2);
void sort_double2(double *data,int size);

double vector_similarity(void) {

  double v1[4],v2[4],sum,dist;
  double *d;
  int numRandom;
  register int i,j;
 
  numRandom=100000;  
  d=alloc_double(numRandom);

  for (i=0; i<numRandom; i++) {
     for (j=0; j<4; j++) {
        v1[j]=genrand();
        v2[j]=genrand();
     }
     sum=0; for (j=0; j<4; j++) sum+=v1[j];
     if (sum!=0) { for (j=0; j<4; j++) v1[j]/=sum; }
     else        { for (j=0; j<4; j++) v1[j]=0.25; }

     sum=0; for (j=0; j<4; j++) sum+=v2[j];
     if (sum!=0) { for (j=0; j<4; j++) v2[j]/=sum; }
     else        { for (j=0; j<4; j++) v2[j]=0.25; }

     d[i]=0; for (j=0; j<4; j++) { d[i] +=fabs(v1[j]-v2[j]); } 
  }

  // sort in increasing order
  sort_double2(d,numRandom);
  dist=d[(int)(numRandom*SIMILARITY_ALPHA)];

  if (d) { free(d); d=NULL; }
  return (dist); 
}

/* increasing */
void sort_double2(double *data,int size) {

   int (*compar)(const void *,const void *);

   compar = Compare_double2;
   qsort((void *)data,(size_t)size,sizeof(double),compar);
}

int Compare_double2(const void *s1, const void *s2) {

   if (*((double *)s1)< (*(double *)s2)) { return -1; }
   if (*((double *)s1)> (*(double *)s2)) { return  1; }
      return 0;
}

int check_pwm_uniqueness_dist(double ***opwm,int *pwmLen,int populationSize,Fitness *fitness,
   double distCutoff,double EvalueCutoff,char *uniqMotif,int window) {

   register int ii,jj,i,j,l,m;
   int numUniq,similarPlus,similarMinus,u_id,p_id;
   int *uniqID;
   double dist;
   double ***rpwm; 

   uniqID=alloc_int(populationSize);
   rpwm=alloc_double_double_double(populationSize,MAX_PWM_LENGTH,4);

   // reverse complement
   for (i=0; i<populationSize; i++) {
      for (l=0; l<pwmLen[i]; l++) {
         for (m=0; m<4; m++) { rpwm[i][l][m]=opwm[i][pwmLen[i]-l-1][3-m]; }
      }
   }

   for (i=0; i<populationSize; i++) uniqMotif[i]='0';
   if (fabs(fitness[0].value-DUMMY_FITNESS)>0.1) { uniqMotif[0]='1'; numUniq=1; } // always designate the first motif unique
   else numUniq=0;

   uniqID[0]=fitness[0].index;

   for (i=1; i<populationSize; i++) {

      // fitness are sorted in increasing order
      if (fitness[i].value>EvalueCutoff) break;

      p_id=fitness[i].index;        // GA population member ID
      similarPlus=0; similarMinus=0;// plus and minus orientation similarity
      for (j=0; j<numUniq; j++) {
         u_id=uniqID[j];
         for (jj=0; jj<pwmLen[p_id]-window+1; jj++) {
            for (ii=0; ii<pwmLen[u_id]-window+1; ii++) {

               dist=0;
               for (l=0; l<window; l++) {
                  for (m=0; m<4; m++) { dist +=fabs(opwm[u_id][l+ii][m]-opwm[p_id][l+jj][m]); }
               }
               if (dist<=window*distCutoff) {
                  similarPlus=1; break;
               }
            }
            if (similarPlus==1) break;
         }
         if (similarPlus==1) break;
      }
      if (!similarPlus) {
         for (j=0; j<numUniq; j++) {
            u_id=uniqID[j];
            for (jj=0; jj<pwmLen[p_id]-window+1; jj++) {
               for (ii=0; ii<pwmLen[u_id]-window+1; ii++) {

                  dist=0;
                  for (l=0; l<window; l++) {
                     for (m=0; m<4; m++) { dist +=fabs(opwm[u_id][l+ii][m]-rpwm[p_id][l+jj][m]); }
                  }
                  if (dist<=window*distCutoff) { similarMinus=1; break; }
               }
               if (similarMinus==1) break;
            }
            if (similarMinus==1) break;
         }
      }
      if (!similarPlus && !similarMinus) {
         uniqMotif[i]='1';
         uniqID[numUniq]=fitness[i].index;
         numUniq++;
      }
   }
   printf("number of unique: %d\n",numUniq);

   if (uniqID)     { free(uniqID);     uniqID=NULL;     }
   if (rpwm[0][0]) { free(rpwm[0][0]); rpwm[0][0]=NULL; }
   if (rpwm[0])    { free(rpwm[0]);    rpwm[0]=NULL;    }
   if (rpwm)       { free(rpwm);       rpwm=NULL;       }
   
   return (numUniq);
}
