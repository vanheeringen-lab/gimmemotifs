
#include <stdlib.h>
#include <string.h>

void reverse_seq(char **seq,char **rSeq,int numSeq,int *seqLen) {

   register int i,j,k;

   for (i=0; i<numSeq; i++) {
      for (k=0,j=seqLen[i]-1; j>=0; j--,k++) {
         switch (seq[i][j]) {
            case 'a': rSeq[i][k]='t'; break;
            case 'c': rSeq[i][k]='g'; break;
            case 'g': rSeq[i][k]='c'; break;
            case 't': rSeq[i][k]='a'; break;
            case 'n': rSeq[i][k]='n'; break;
            default: break;
         }
      }
      rSeq[i][seqLen[i]]='\0';
   }
}

