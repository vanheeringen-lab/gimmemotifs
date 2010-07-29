#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defines.h"

double **alloc_double_double(int ,int );

void standardize_pwm(double **pwm,int pwmLen) {

   register int i,j;
   double sum;

   for (i=0; i<pwmLen; i++) {
      sum=0.0; for (j=0; j<4; j++)  sum +=pwm[i][j]; 
     
      if (sum>0.01) { 
         for (j=0; j<4; j++) pwm[i][j] /=sum;  
      }
      else { 
         for (j=0; j<4; j++) pwm[i][j]  =0.25; 
      } 
   }
}

void log_pwm(double **pwm,double **logpwm,int pwmLen){

   register int i,j;

   for (i=0; i<pwmLen; i++) {
      if (pwm[i][0]<PSEUDO_COUNT || pwm[i][1]<PSEUDO_COUNT || pwm[i][2]<PSEUDO_COUNT || pwm[i][3]<PSEUDO_COUNT) {
          for (j=0; j<4; j++) logpwm[i][j]=log((pwm[i][j]+PSEUDO_COUNT)/(1.0+PSEUDO_COUNT*4));
      }
      else {
         for (j=0; j<4; j++) logpwm[i][j]=log(pwm[i][j]);
      }
   }
}

void log_ratio_to_int(double **pwm,int **ipwm,int pwmLen,double *bfreq){

   register int i,j;
   double cellMin,cellMax;
   double **logratio;

   logratio=alloc_double_double(pwmLen,4);
 
   for (i=0; i<pwmLen; i++) {
      if (pwm[i][0]<PSEUDO_COUNT || pwm[i][1]<PSEUDO_COUNT || pwm[i][2]<PSEUDO_COUNT || pwm[i][3]<PSEUDO_COUNT) {
          for (j=0; j<4; j++) logratio[i][j]=log((pwm[i][j]+PSEUDO_COUNT)/(1.0+PSEUDO_COUNT*4)/bfreq[j]);
      }
      else {
          for (j=0; j<4; j++) logratio[i][j]=log(pwm[i][j]/bfreq[j]);
      }
   }
   //for (j=0; j<4; j++) {
   //   for (i=0; i<pwmLen;  i++) printf("%4.2f ",logratio[i][j]); printf("\n");
   //}

   cellMin=logratio[0][0]; cellMax=logratio[0][0]; 
   for (i=0; i<pwmLen; i++) {
      for (j=0; j<4; j++) {
        if (logratio[i][j]<cellMin) cellMin=logratio[i][j]; 
        if (logratio[i][j]>cellMax) cellMax=logratio[i][j]; 
      } 
   }
   //for (j=0; j<4; j++) {
   //   for (i=0; i<pwmLen;  i++) printf("%4.3f ",(logratio[i][j]-cellMin)/(cellMax-cellMin)); printf("\n");
   //}

   if (cellMax-cellMin<0.01) {
      for (i=0; i<pwmLen;  i++) {
         for (j=0; j<4; j++) ipwm[i][j]=DOUBLE_TO_INT_SCALE/4;
         // for (j=0; j<4; j++) ipwm[i][j]=0;
      } 
   }
   else {
      for (i=0; i<pwmLen;  i++) {
         for (j=0; j<4; j++) {
            ipwm[i][j]=(int)(DOUBLE_TO_INT_SCALE*((logratio[i][j]-cellMin)/(cellMax-cellMin))); 
         } 
      }
   }
   if (logratio[0]) { free(logratio[0]); logratio[0]=NULL; }
   if (logratio)    { free(logratio);    logratio=NULL;    }

   //for (j=0; j<4; j++) {
   //   for (i=0; i<pwmLen;  i++) printf("%4d",ipwm[i][j]); printf("\n");
   //}
}

