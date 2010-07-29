#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include "defines.h"

int range(double score);

void normalize(double **score,double **rscore,int *seqLen,int pwmLen,int numSeq,char *Iseq,
   int maxp,double **weight,int weightType) {

   register int i,j,k,m;
   double totalOffset,localp,normalizeFactor,temp;
   int renormalize,totalSites;

   // convert log scale back to the original scale
   for(i=0; i<numSeq; i++){
      if (Iseq[i]=='0') continue;
      for(j=0; j<seqLen[i]-pwmLen+1; j++){
         score[i][j]=exp(score[i][j]); 
         if (score[i][j]==1)   score[i][j]=0.999999999;
         rscore[i][j]=exp(rscore[i][j]);
         if (rscore[i][j]==1) rscore[i][j]=0.999999999; 
      }
   }

   if (weightType!=0) {
      for (i=0; i<numSeq; i++) {
         for (j=0; j<seqLen[j]; j++) {
            score[i][j] *=weight[i][j]; rscore[i][j] *=weight[i][j];
         }
      }
   }

   totalOffset=0; 
   for (i=0; i<numSeq; i++) {
      if (Iseq[i]=='0') continue;
      for (j=0; j<seqLen[i]-pwmLen+1; j++) totalOffset += (score[i][j]+rscore[i][j]);
   }

   for (i=0; i<numSeq; i++) {
      if (Iseq[i]=='0') continue;
      for (j=0; j<seqLen[i]-pwmLen+1; j++) { 
         score[i][j] /=totalOffset; rscore[i][j]/=totalOffset; 
      }
   }

   // squash as in MEME
   totalSites=maxp;
   renormalize=1; totalOffset=1;
   while (renormalize) {
      renormalize=0;
      if (totalSites<1) break;

      // avoid dividing a small number by a small number
      if(totalOffset<10e-10 && totalSites>0){
         temp=range(totalOffset);
         normalizeFactor=(double)(totalOffset*pow(10,temp))/(double)(totalSites*pow(10,temp));
      }
      else {
         normalizeFactor=totalOffset/(double)totalSites;
      }

      totalOffset=0;
      for (i=0; i<numSeq; i++) {
         if (Iseq[i]=='0') continue;
         for (j=0; j<seqLen[i]-pwmLen+1; j++) {
            if (score[i][j]==0) continue;
            if (score[i][j]<1) {
               score[i][j]/=normalizeFactor; 
               if (score[i][j]>=1) {
                  score[i][j]=1.0; totalSites--;
                  renormalize=1; 
               }
            }
            if (rscore[i][j]==0) continue;
            if (rscore[i][j]<1) {
               rscore[i][j]/=normalizeFactor; 
               if (rscore[i][j]>=1) {
                  rscore[i][j]=1; totalSites--;
                  renormalize=1; 
               }
            }
            if (score[i][j]<1)  totalOffset += score[i][j]; 
            if (rscore[i][j]<1) totalOffset +=rscore[i][j]; 
         } 
      }

   };

  // smooth as in MEME
  for (i=0; i<numSeq; i++) {
      if (Iseq[i]=='0') continue;
      for(m=0; m<pwmLen; m++) {
         for(j=0; j<(seqLen[i]-m)/pwmLen; j++) {
            localp=0.0; 
            for(k=0; k<pwmLen; k++) localp+=(score[i][j*pwmLen+m+k]+rscore[i][j*pwmLen+m+k]);
             
            if(localp>1.0) {
               for(k=0; k<pwmLen; k++) {
                  score[i][j*pwmLen+m+k] = score[i][j*pwmLen+m+k]/localp;
                  rscore[i][j*pwmLen+m+k]=rscore[i][j*pwmLen+m+k]/localp;
               }
            }
         }
      }
   }
}

int range(double score) {

  register int i;
  double temp;

  for(i=1; i<21; i++){
     temp=score*pow(10,i);
     if(temp>=1.0) break;
  }
  if(i==20) printf("score is smaller than 10e-20\n");
  return(i);
}

