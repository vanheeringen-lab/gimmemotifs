#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "defines.h"
#include "gadem.h"

void assign_weight_uniform(int *seqLen,int numSeq,double **weight) {

   register int i,j;

   for (i=0; i<numSeq; i++) {
      for (j=0; j<seqLen[i]; j++)  weight[i][j]=1.0/(double)(2*seqLen[i]); // proximate 1/2*(L-w+0)
      // for (j=0; j<seqLen[i]; j++) printf("%8.7f\n",weight[i][j]); exit(0);
   }
}

void assign_weight_normal(int *seqLen,int numSeq,double **weight) {

   register int i,j;
   double sigma,mean,pi,sum;

   sigma=25.0; pi=3.1415926;
   for (i=0; i<numSeq; i++) {
      mean=seqLen[i]/2;
      // peakHeight=1.0/(sigma*sqrt(2.0*pi)); 
      for (j=0; j<seqLen[i]; j++) {
          // weight[i][j]=exp(-(j-mean)*(j-mean)/(2*sigma*sigma))/(sigma*sqrt(2.0*pi))/peakHeight; 
          weight[i][j]=exp(-(j-mean)*(j-mean)/(2*sigma*sigma))/(sigma*sqrt(2.0*pi)); 
      } 
      sum=0; for (j=0; j<seqLen[i]; j++) sum +=weight[i][j];
      for (j=0; j<seqLen[i]; j++) weight[i][j] /=(2.0*sum);
      // for (j=0; j<seqLen[i]; j++) printf("%8.7f\n",weight[i][j]); exit(0);
   }
}

void assign_weight_triangular(int *seqLen,int numSeq,double **weight) {

   register int i,j,k;
   double sum;

   for (i=0; i<numSeq; i++) {
      for (j=0; j<seqLen[i]/2; j++) {
          weight[i][j]=(double)(2*j)/(double)seqLen[i]; 
      } 
      for (k=seqLen[i]/2,j=seqLen[i]/2; j<seqLen[i]; j++,k--) {
          weight[i][j]=(double)(2*k)/(double)seqLen[i]; 
      } 
      sum=0; for (j=0; j<seqLen[i]; j++) sum +=weight[i][j];
      for (j=0; j<seqLen[i]; j++) weight[i][j]/=(2.0*sum);
      // for (j=0; j<seqLen[i]; j++) printf("%8.7f\n",weight[i][j]); exit(0);
   }
}

void assign_weight_triangular_uniform(int *seqLen,int numSeq,double **weight,int coreLen) {

   register int i,j,k;
   int halfWidth;
   double sum;

   halfWidth=coreLen/2;
 
   for (i=0; i<numSeq; i++) {
      if (2*halfWidth>=seqLen[i]) halfWidth=seqLen[i]/4;

      // exponential from say 0-175
      for (j=0; j<seqLen[i]/2-halfWidth; j++) {
          weight[i][j]=(double)(2*j)/(double)(seqLen[i]-2*halfWidth); 
      } 
      // uniform from say 175-225
      for (j=seqLen[i]/2-halfWidth; j<seqLen[i]/2+halfWidth; j++) {
         weight[i][j]=1.0;  
      }
      // exponential from say 225-400
      for (k=seqLen[i]/2-halfWidth,j=seqLen[i]/2+halfWidth; j<seqLen[i]; j++,k--) {
          weight[i][j]=(double)(2*k)/(double)(seqLen[i]-2*halfWidth); 
      } 
      sum=0; for (j=0; j<seqLen[i]; j++) sum +=weight[i][j];
      for (j=0; j<seqLen[i]; j++) weight[i][j]/=(2.0*sum);
      // for (j=0; j<seqLen[i]; j++) printf("%8.7f\n",weight[i][j]); exit(0);
   }
}

void assign_weight_rectangle(int *seqLen,int numSeq,double **weight,int coreLen) {

   register int i,j,k;
   int halfWidth;
   double lamda,sum;

   lamda=1.0;
   halfWidth=coreLen/2;
 
   for (i=0; i<numSeq; i++) {
      if (coreLen>=seqLen[i]) halfWidth=seqLen[i]/4;

      // zero weight from say 0-175
      for (k=seqLen[i]/2-halfWidth,j=0; j<seqLen[i]/2-halfWidth; j++,k--) {
          // weight[i][j]=lamda*exp(-lamda*(double)k); 
          weight[i][j]=0;
      }
      // uniform weight (1.0) from say 175-225
      for (j=seqLen[i]/2-halfWidth; j<seqLen[i]/2+halfWidth; j++) {
         weight[i][j]=lamda;  
      }
      // zero weight from say 225-400
      for (k=0,j=seqLen[i]/2+halfWidth; j<seqLen[i]; j++,k++) {
         // weight[i][j]=lamda*exp(-lamda*(double)k); 
         weight[i][j]=0;
      }
      sum=0; for (j=0; j<seqLen[i]; j++) sum+=weight[i][j];
      for (j=0; j<seqLen[i]; j++) weight[i][j]/=(2.0*sum);
      // for (j=0; j<seqLen[i]; j++) printf("%8.7f\n",weight[i][j]); exit(0);
   }
}
