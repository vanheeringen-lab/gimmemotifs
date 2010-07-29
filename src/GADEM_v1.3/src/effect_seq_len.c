#include <stdlib.h>

void effect_seq_length(char **seq,int numSeq,int *seqLen,char *Iseq,int *emSeqLen) {

   register int i,j;
   int emSeqCn,cn;

   emSeqCn=0;
   for (i=0; i<numSeq; i++) {
      if (Iseq[i]=='1') {
         cn=0;
         for (j=0; j<seqLen[i]; j++) {
            if (seq[i][j]!='n') cn++;
         }
         emSeqLen[emSeqCn]=cn; emSeqCn++;
      }
   }
}

void effect_seq_length_full(char **seq,int numSeq,int *seqLen,int *eSeqLen) {

   register int i,j;
   int cn;

   for (i=0; i<numSeq; i++) {
      cn=0;
      for (j=0; j<seqLen[i]; j++) {
         if (seq[i][j]!='n') cn++;
      }
      eSeqLen[i]=cn;
   }
}
