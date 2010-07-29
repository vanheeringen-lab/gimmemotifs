#include <stdlib.h>
#include <string.h>
#include "random.h"


void sample_without_replacement(char *which,int numReplacement,int numSeq) {

   register int i;
   int numFilled,dummy;
   double du;

   for (i=0; i<numSeq; i++) which[i]='0';

   numFilled=0;
   while (numFilled<numReplacement) {
      du=(double)numSeq*genrand();
      dummy=(int)du;
      if (dummy==numSeq) dummy--;

      if (which[dummy]=='0') { which[dummy]='1'; numFilled++; }
   }
}

void sample_without_replacement2(int *which,int numReplacement,int numSeq) {

   register int i;
   int used,numFilled,dummy;
   double du;

   for (i=0; i<numSeq; i++) which[i]=-1;

   numFilled=0;
   while (numFilled<numReplacement) {
      du=(double)numSeq*genrand();
      dummy=(int)du;
      if (dummy==numSeq) dummy--;

      used=0;
      for (i=0;i<numFilled; i++) {
         if (which[i]==dummy) { used=1; break; } 
      }
      if (!used) { 
         which[numFilled]=dummy; numFilled++; 
      }
   }
}

