
#include <stdlib.h>
#include "gadem.h"

//mask sites that have been found
void mask_sites(int predictedSiteCn,char **seq,char **rseq,int *seqLen,Sites *site,int pwmLen) {

   register int i,m;
   int pos,seqID;

   for (i=0; i<predictedSiteCn; i++) {
      pos=site[i].pos; seqID=site[i].seq;
      if (site[i].rev=='0') {
         for (m=0; m<pwmLen; m++) seq[seqID][pos+m]='n';
         for (m=0; m<pwmLen; m++) rseq[seqID][seqLen[seqID]-pwmLen+m-pos]='n';
      }
      else {
         for (m=0; m<pwmLen; m++) rseq[seqID][pos+m]='n';
         for (m=0; m<pwmLen; m++) seq[seqID][seqLen[seqID]-pwmLen+m-pos]='n';
      }
   }
}

