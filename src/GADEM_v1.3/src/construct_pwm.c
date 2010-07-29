
#include <stdlib.h>
#include <string.h>

void construct_pwm(double **pwm,double **score,double **rscore,char **seq,char **rseq,
   int *seqLen,int numseq,int motifLen,char *Iseq) {

   register int i,j,k;
   int posOnPlus;

   for (i=0; i<motifLen; i++) {
      for (j=0; j<4; j++) pwm[i][j]=0;
   }

   for(i=0; i<numseq; i++) {
      if (Iseq[i]=='0') continue;
      for(j=0; j<seqLen[i]-motifLen+1; j++){
         if (score[i][j]==0) continue;
         for(k=0; k<motifLen; k++) {
            switch (seq[i][j+k]) {
               case 'a': pwm[k][0] +=score[i][j]; break;
               case 'c': pwm[k][1] +=score[i][j]; break;
               case 'g': pwm[k][2] +=score[i][j]; break;
               case 't': pwm[k][3] +=score[i][j]; break;
               default: break;
            }
         }
      }

      for(j=0; j<seqLen[i]-motifLen+1; j++){
         if (rscore[i][j]==0) continue;
         posOnPlus=seqLen[i]-motifLen-j;
         for(k=0; k<motifLen; k++) {
            switch (rseq[i][posOnPlus+k]) {
               case 'a': pwm[k][0] +=rscore[i][j]; break;
               case 'c': pwm[k][1] +=rscore[i][j]; break;
               case 'g': pwm[k][2] +=rscore[i][j]; break;
               case 't': pwm[k][3] +=rscore[i][j]; break;
               default: break;
            }
         }
      }
   }
}

