#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "gadem.h"
#include "defines.h"

// check at most 6bp at each site of the alignment to see if the alignment
// can be extended/trimmed. This was done one base at a time

void extend_alignment(Sites *site,int numSeq,char **seq,char **rseq,int *seqLen,int nsites,
   int pwmLen,int *pwmnewLen) {

   register int i,j;
   int ii;
   int start,end,lshift,rshift,totalLen;
   double **pwm,*info;
   double colSum;

   lshift=BASE_EXT; rshift=BASE_EXT;
   totalLen=lshift+pwmLen+rshift;

   pwm=alloc_double_double(totalLen,4);
   info=alloc_double(totalLen);

   for (ii=0; ii<totalLen; ii++) {
      for (i=0; i<4; i++) pwm[ii][i]=PSEUDO_COUNT; 
   }

   for (ii=0; ii<lshift; ii++) {
      for (i=0; i<nsites; i++) {
         if (site[i].pos-lshift+ii>=0) {
            if (site[i].rev=='0') {
               switch (seq[site[i].seq][site[i].pos-lshift+ii]) {
                  case 'a': pwm[ii][0] +=1; break;
                  case 'c': pwm[ii][1] +=1; break;
                  case 'g': pwm[ii][2] +=1; break;
                  case 't': pwm[ii][3] +=1; break;
                  default: break;
               }
            }
            else {
               switch (rseq[site[i].seq][site[i].pos-lshift+ii]) {
                  case 'a': pwm[ii][0] +=1; break;
                  case 'c': pwm[ii][1] +=1; break;
                  case 'g': pwm[ii][2] +=1; break;
                  case 't': pwm[ii][3] +=1; break;
                  default: break;
               }
            }
         }
      }
   }

   for (i=0; i<nsites; i++) {
      if (site[i].rev=='0') {
         for (j=0; j<pwmLen; j++) {
            switch (seq[site[i].seq][site[i].pos+j]) {
               case 'a': pwm[lshift+j][0] +=1; break;
               case 'c': pwm[lshift+j][1] +=1; break;
               case 'g': pwm[lshift+j][2] +=1; break;
               case 't': pwm[lshift+j][3] +=1; break;
               default: break;
            }
         }
      }
      else {
         for (j=0; j<pwmLen; j++) {
            switch (rseq[site[i].seq][site[i].pos+j]) {
               case 'a': pwm[lshift+j][0] +=1; break;
               case 'c': pwm[lshift+j][1] +=1; break;
               case 'g': pwm[lshift+j][2] +=1; break;
               case 't': pwm[lshift+j][3] +=1; break;
               default: break;
            }
         }
      }
   }

   for (ii=0; ii<rshift; ii++) {
      for (i=0; i<nsites; i++) {
         if (site[i].pos+pwmLen+ii<seqLen[site[i].seq]-1) {
            if (site[i].rev=='0') {
               switch (seq[site[i].seq][site[i].pos+pwmLen+ii]) {
                  case 'a': pwm[lshift+pwmLen+ii][0] +=1; break;
                  case 'c': pwm[lshift+pwmLen+ii][1] +=1; break;
                  case 'g': pwm[lshift+pwmLen+ii][2] +=1; break;
                  case 't': pwm[lshift+pwmLen+ii][3] +=1; break;
                  default: break;
               }
            }
            else {
               switch (rseq[site[i].seq][site[i].pos+pwmLen+ii]) {
                  case 'a': pwm[lshift+pwmLen+ii][0] +=1; break;
                  case 'c': pwm[lshift+pwmLen+ii][1] +=1; break;
                  case 'g': pwm[lshift+pwmLen+ii][2] +=1; break;
                  case 't': pwm[lshift+pwmLen+ii][3] +=1; break;
                  default: break;
               }
            }
         }
      }
   }

   for (ii=0; ii<totalLen; ii++) {
      colSum=0; for (j=0; j<4; j++) colSum+=pwm[ii][j];
      if (colSum>PSEUDO_COUNT) { for (j=0; j<4; j++) pwm[ii][j]/=colSum; }
      else                     { for (j=0; j<4; j++) pwm[ii][j] =0.25;   }

      info[ii]=2.0; 
      for (j=0; j<4; j++) { 
         if (pwm[ii][j]>PSEUDO_COUNT) info[ii] +=(pwm[ii][j]*log(pwm[ii][j])/log(2.0)); 
      }
   }


   /*--------------------------------------------------------------------------------
   // debugging
   for (i=0; i<nsites; i++) {
      if (site[i].rev=='0') {
         for (j=0; j<pwmLen; j++) printf("%c",seq[site[i].seq][site[i].pos+j]); 
            printf("\t(+)\t%d\t%d\n",site[i].seq+1,site[i].pos); 
      }
      else {
         for (j=0; j<pwmLen; j++) printf("%c",rseq[site[i].seq][site[i].pos+j]);
            printf("\t(-)\t%d\t%d\n",site[i].seq+1,site[i].pos); 
      }
   }
   printf("\n");

   for (i=0; i<pwmLen+lshift+rshift; i++) printf("%4.2f ",info[i]);  printf("\n");
   --------------------------------------------------------------------------------*/
   // cutoffs are relative to a maximal information of 2 bits
   // for (i=0; i<pwmLen+lshift+rshift; i++) printf("%4.2f ",info[i]);  printf("\n");
   start=0;
   for (i=0; i<lshift+rshift+pwmLen-2; i++) {
      if ((info[i]>=  MIN_BITS1 && info[i+1]>=MIN_BITS1 && info[i+2]>=MIN_BITS1) ||
          (info[i]>=  MIN_BITS2 && info[i+1]>=MIN_BITS2) ||
          (info[i+1]>=MIN_BITS2 && info[i+2]>=MIN_BITS2) ||
          (info[i]>=  MIN_BITS3) ) { 
         for (j=0; j<nsites; j++) site[j].pos +=(i-lshift); 
         start=i; break; 
      }
   }
   // printf("start: %d\n",start);

   end=0;
   for (i=rshift+lshift+pwmLen-1; i>1; i--) {
      if ((info[i]>=  MIN_BITS1 && info[i-1]>=MIN_BITS1 && info[i-2]>=MIN_BITS1) ||
          (info[i-1]>=MIN_BITS2 && info[i-2]>=MIN_BITS2) ||
          (info[i]>=  MIN_BITS2 && info[i-1]>=MIN_BITS2) ||
          (info[i]>=  MIN_BITS3) ) { 
         // printf("to right %d\n",i);
         end=i; break; 
      }
   }
   // printf("end: %d\n",end);
   // new motif length
   if (end-start+1>0) *pwmnewLen=end-start+1;
   else *pwmnewLen=0;

   if (pwm[0])  { free(pwm[0]); pwm[0]=NULL; }
   if (pwm)     { free(pwm);    pwm=NULL;    }
   if (info)    { free(info);   info=NULL;   }
}

