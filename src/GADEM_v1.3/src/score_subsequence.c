#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include "defines.h"

double *alloc_double(int );

void ll_score_motif_model(int numSeq,char **seq,char **rseq,int *seqLen,double **pwm,int pwmLen,
   double **score,double **rscore,char *Iseq,double *bfreq) {

   register int ii,i,m;
   double plusScore,minusScore;
   double *colAve;
   int posOnPlus;

   colAve=alloc_double(pwmLen);

   // weighted average for nucleotide 'n'
   for (m=0; m<pwmLen; m++) {
      colAve[m]=0; for (i=0; i<4; i++) { colAve[m] +=(pwm[m][i]*bfreq[i]); }
   }

   for (ii=0;  ii<numSeq; ii++) {
      if (Iseq[ii]=='0') continue;

      // plus strand
      for(i=0;  i<seqLen[ii]-pwmLen+1; i++) {
         plusScore=0;
         for (m=0; m<pwmLen;m++) {
            switch(seq[ii][i+m]) {
               case 'a': plusScore +=pwm[m][0]; break;
               case 'c': plusScore +=pwm[m][1]; break;
               case 'g': plusScore +=pwm[m][2]; break;
               case 't': plusScore +=pwm[m][3]; break;
               default:  plusScore +=colAve[m]; break;
            }
         }
         score[ii][i]=plusScore;
      }
 
      // minus strand
      for(i=0; i<seqLen[ii]-pwmLen+1; i++){
         posOnPlus=seqLen[ii]-i-pwmLen;
         minusScore=0;
         for (m=0; m<pwmLen;m++) {
            switch(rseq[ii][i+m]) {
               case 'a': minusScore +=pwm[m][0]; break;
               case 'c': minusScore +=pwm[m][1]; break;
               case 'g': minusScore +=pwm[m][2]; break;
               case 't': minusScore +=pwm[m][3]; break;
               default:  minusScore +=colAve[m]; break;
            }
         }
         rscore[ii][posOnPlus]=minusScore;
      }
   }

   if (colAve) { free(colAve); colAve=NULL; }
}

int llr_score(double *score,int numSeq,char **seq,int *seqLen,double **logepwm,
   int pwmLen,double *bfreq,double **bscore) {

   register int ii,i,m;
   int numKmer;
   double plusScore;
   double *colAve;

   colAve=alloc_double(pwmLen);
   for (m=0; m<pwmLen; m++) {
      colAve[m]=0; for (i=0; i<4; i++) { colAve[m] +=(logepwm[m][i]*bfreq[i]); } 
   }

   numKmer=0;
   for (ii=0; ii<numSeq; ii++) {
      for (i=0; i<seqLen[ii]-pwmLen+1; i++) {
         plusScore=0;
         for (m=0; m<pwmLen; m++) {
            switch(seq[ii][i+m]) {
               case 'a': plusScore +=logepwm[m][0]; break;
               case 'c': plusScore +=logepwm[m][1]; break;
               case 'g': plusScore +=logepwm[m][2]; break;
               case 't': plusScore +=logepwm[m][3]; break;
               default:  plusScore +=colAve[m];     break;
            }
         }
         score[numKmer]=plusScore-bscore[ii][i]; //llr
         numKmer++; 
      }
   }
   if (colAve) { free(colAve); colAve=NULL; }

   return (numKmer);
}
