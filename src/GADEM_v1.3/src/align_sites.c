
#include <stdlib.h>
#include <string.h>
#include "gadem.h"

void align_sites_count(Sites *site,char **seq,char **rseq,int nsites,int pwmLen,double **pwm) {

   register int i,j;

   for (j=0; j<pwmLen; j++) {
      for (i=0; i<4; i++) pwm[j][i]=0;
   }
   for (i=0; i<nsites; i++) {
      // plus strand
      if (site[i].rev=='0') {
         for (j=0; j<pwmLen; j++) {
            //printf("%c",seq[site[i].seq][site[i].pos+j]); 
            switch (seq[site[i].seq][site[i].pos+j]) {
               case 'a': pwm[j][0] +=1; break;
               case 'c': pwm[j][1] +=1; break;
               case 'g': pwm[j][2] +=1; break;
               case 't': pwm[j][3] +=1; break;
               default: break;
            }
         }
         //printf("\n");
      }
      // reverse strand
      else {
         for (j=0; j<pwmLen; j++) {
            //printf("%c",rseq[site[i].seq][site[i].pos+j]); 
            switch (rseq[site[i].seq][site[i].pos+j]) {
               case 'a': pwm[j][0] +=1; break;
               case 'c': pwm[j][1] +=1; break;
               case 'g': pwm[j][2] +=1; break;
               case 't': pwm[j][3] +=1; break;
               default: break;
            }
         }
         //printf("\n");
      }
   }
}

