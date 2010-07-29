#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include "gadem.h"
#include "defines.h"
#include "random.h"

void mutation(Chrs **dyad,int numWordGroup,Words *word,int minSpaceWidth,int maxSpaceWidth,
   Wheel *wheel,int populationSize,Fitness *fitness,char *uniqMotif,double *maxpFactor,double maxpMutationRate) {

   register int i,j;
   int popuCn;             // index of the number of dyads in population
   int whichToMutate;      // the spaced dyad for mutation
   int kmerGroup;          // pick a k-mer group (tetramer,pentamer,hexamer) in which a new word will be selected
   int whichPartDyad;      // pick a component of a spaced dyad for mutation
   int spaceWidth;
   double *tmpmaxpFactor;
   double rand,newFactor;
   Chrs **tmpDyad;

   tmpDyad=alloc_chrs(populationSize,4);
   tmpmaxpFactor=alloc_double(populationSize);

   popuCn=0; 
   for (i=0; i<populationSize; i++) {
      if (uniqMotif[i]=='1') {
         for (j=0; j<3; j++) {
            tmpDyad[popuCn][j].wordID   =dyad[fitness[i].index][j].wordID; 
            tmpDyad[popuCn][j].wordGroup=dyad[fitness[i].index][j].wordGroup; 
         }
         tmpmaxpFactor[popuCn]=maxpFactor[fitness[i].index];
         popuCn++;
      }
   }
   // printf("popuCn=%d\n",popuCn);

   while (popuCn<populationSize) {

      // sample with replacement with probability proportional to fitness score
      rand=(double)populationSize*genrand();
      whichToMutate=0;
      for (j=0; j<populationSize; j++) {
         if (rand>=wheel[j].start && rand<=wheel[j].end) { 
            whichToMutate=wheel[j].index; break; 
         }
      }
      tmpmaxpFactor[popuCn]=maxpFactor[whichToMutate];
      for (j=0; j<3; j++) { 
         tmpDyad[popuCn][j].wordID   =dyad[whichToMutate][j].wordID;
         tmpDyad[popuCn][j].wordGroup=dyad[whichToMutate][j].wordGroup; 
      }
    
      if (genrand()<maxpMutationRate) {
         // find a different maxp factor, keep dyad unchanged
         do {
            newFactor=MAXP_BASE+MAXP_FACTOR*((int)(MAXP_SCALE*genrand()));
         } while (fabs(newFactor-tmpmaxpFactor[popuCn])<0.001);
         tmpmaxpFactor[popuCn]=newFactor;
      }
      else {
         // change one component of the dyad, first find out how many components there are
         if (maxSpaceWidth==0) { 
            whichPartDyad=(int)(2*genrand());      // change w1
            if (whichPartDyad==1) whichPartDyad=2; // skip spacer and change w2
         }
         else {
            whichPartDyad=(int)(3*genrand());
            if (whichPartDyad==3) whichPartDyad--;
         }
         
         // if w1 of the dyad (w1-s-w2) is chosen
         if (whichPartDyad==0) {
            // decide which kmer kmerGroup (4-mer, 5-mer or 6-mer) from which w1 whill be replaced
            kmerGroup=(int)(numWordGroup*genrand());
            if (kmerGroup==numWordGroup) kmerGroup--;
   
            // select a k-mer from the k-group with prob proportional to z-score
            rand=genrand();
            tmpDyad[popuCn][0].wordID=0;
            for (j=0; j<word[kmerGroup].count; j++) {
               if (rand>=word[kmerGroup].prob_sta[j] && rand<=word[kmerGroup].prob_end[j]) { 
                  tmpDyad[popuCn][0].wordID=j; break; 
               }
            }
            tmpDyad[popuCn][0].wordGroup=kmerGroup;        // mark the kmer kmerGroup from which w1 is selected
         }
   
         // if the spacer is chosen
         else if (whichPartDyad==1) {
            // choose a gap that is different from the current one
            do {
               spaceWidth=minSpaceWidth+(int)((maxSpaceWidth-minSpaceWidth+1)*genrand());
            } while (spaceWidth==tmpDyad[popuCn][1].wordID); 
            tmpDyad[popuCn][1].wordID=spaceWidth;
            tmpDyad[popuCn][1].wordGroup=-1;                 // dummy - not used
         }
   
         // if w2 of the dyad (w1-space-w2) is chosen 
         else {
            kmerGroup=(int)(numWordGroup*genrand());
            if (kmerGroup==numWordGroup) kmerGroup--;
   
            rand=genrand();
            tmpDyad[popuCn][2].wordID=0;
            for (j=0; j<word[kmerGroup].count; j++) {
               if (rand>=word[kmerGroup].prob_sta[j] && rand<word[kmerGroup].prob_end[j]) { 
                  tmpDyad[popuCn][2].wordID=j; break; 
               }
            }
            tmpDyad[popuCn][2].wordGroup=kmerGroup;     // mark the kmer kmerGroup from which w2 is selected 
         }
      }
      popuCn++;
   }

   // update the population
   for (i=0; i<populationSize; i++) {
      for (j=0; j<3; j++) { 
         dyad[i][j].wordID   =tmpDyad[i][j].wordID; 
         dyad[i][j].wordGroup=tmpDyad[i][j].wordGroup; 
      }
      maxpFactor[i]=tmpmaxpFactor[i]; 
   }

   if (tmpDyad[0])    { free(tmpDyad[0]);    tmpDyad[0]=NULL;    }
   if (tmpDyad)       { free(tmpDyad);       tmpDyad=NULL;       }
   if (tmpmaxpFactor) { free(tmpmaxpFactor); tmpmaxpFactor=NULL; }
}
