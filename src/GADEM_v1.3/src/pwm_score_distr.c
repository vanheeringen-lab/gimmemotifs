#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include "gadem.h"
#include "defines.h"

double find_pvalue(int score,Pgfs *M1,int dimension) {

    register int i;
    double p_result;

    p_result=0.0;
    for (i=0; i<dimension; i++) {
       if (M1[i].score>=score) p_result+=M1[i].prob;
       else break; 
    }
    return(p_result);
}

int ini_M(int pos,Pgfs *M0,int **pwm,double *bfreq) {

   register int i,j;
   int dim0;
   Pgfs *temp;

   temp=alloc_distr(4);

   for (i=0; i<4; i++) {
      temp[i].score=pwm[pos][i]; 
      temp[i].prob =bfreq[i];
   }
   sort_llrDist(temp,4);

   dim0=0;
   for (i=0; i<4; i++) {

      M0[dim0].score=temp[i].score; 
      M0[dim0].prob =temp[i].prob;

      for (j=i; j<4; j++) {
         if (j==i) continue;
         if (temp[i].score==temp[j].score) {
            M0[dim0].prob +=temp[j].prob;
            i=j;
         }
         else { break; }
      }
      dim0++;
   }
   if (temp) { free(temp); temp=NULL; }

   return (dim0);
}

int prod_M(Pgfs *llrDist,int dim0,Pgfs *M1,int dim1) {

   register int i,j;
   int total,dim2;
   Pgfs *temp;

   temp=alloc_distr(MAX_DIMENSION);

   // add up the scores/elements 
   // times the probabilities
   total=0;
   for (i=0; i<dim0; i++) {
       for (j=0; j<dim1; j++) {
           temp[total].score=llrDist[i].score+M1[j].score;
           temp[total].prob =llrDist[i].prob *M1[j].prob;
           total++;
           if (total==MAX_DIMENSION) { 
              printf("\nError: max p-table dimension reached\n");
              printf("  reset <MAX_DIMENSION> in defines.h\n");
              exit(0);
           }
       }
   }

   sort_llrDist(temp,total);
   // combine, if scores/elements are the same

   dim2=0; 
   for (i=0; i<total; i++) {
      llrDist[dim2].score=temp[i].score; 
      llrDist[dim2].prob =temp[i].prob;
      for (j=i; j<total; j++) {
         if (j==i) continue;
         if (temp[i].score==temp[j].score) {
            llrDist[dim2].prob +=temp[j].prob;
            i=j;
         }
         else { break; }
      }
      dim2++;
   }
   if (temp)  {free(temp); temp=NULL; }

   return (dim2);
}

int pwm_score_dist(int **pwm,int pwmLen,Pgfs *llrDist,double *bfreq) {

   register int i;
   int dim0,dim1,dim2,col;
   Pgfs *M1;

   M1=alloc_distr(4);

   col=0;
   dim0=ini_M(col,llrDist,pwm,bfreq);
   dim2=dim0;

   for (i=1; i<pwmLen; i++) {
      col=i;
      dim1=ini_M(col,M1,pwm,bfreq);
      dim2=prod_M(llrDist,dim0,M1,dim1);
      dim0=dim2;
   }

   if (M1)  { free(M1); M1=NULL; }
   return (dim2);
}

int determine_cutoff(Pgfs *scoreDist,int dimension,double p_cutoff_ipwm) {

   register int j;
   int scoreCutoff;
   double pvalue;

   scoreCutoff=scoreDist[0].score;
   pvalue=scoreDist[0].prob;
   for (j=1; j<dimension; j++) {
      pvalue +=scoreDist[j].prob;
      if (pvalue>p_cutoff_ipwm) { scoreCutoff=scoreDist[j-1].score; break; }
   }

   return (scoreCutoff);
}

