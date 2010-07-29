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

void crossover(Chrs **dyad,int numWordGroup,Words *word,int minSpaceWidth,int maxSpaceWidth,
   Wheel *wheel,int populationSize,Fitness *fitness,char *motifUniq,double *maxpFactor,double maxpMutationRate) {

   register int i,j;
   int found,id1,id2,popuCn;
   double *tmpmaxpFactor;
   Chrs **tmpDyad;
   double rand,du,newFactor;

   tmpDyad=alloc_chrs(populationSize,4);
   tmpmaxpFactor=alloc_double(populationSize);

   // save a copy of all unique top-ranked spaced dyads
   popuCn=0;
   for (i=0; i<populationSize; i++) {
      if (motifUniq[i]=='1') {
         for (j=0; j<3; j++) {
            tmpDyad[popuCn][j].wordID   =dyad[fitness[i].index][j].wordID;
            tmpDyad[popuCn][j].wordGroup=dyad[fitness[i].index][j].wordGroup;
         }
         tmpmaxpFactor[popuCn]=maxpFactor[fitness[i].index];
         popuCn++;
      }
   }

   // crossover
   while (popuCn<populationSize) {
      do {
         // select the first dyad 
         du=(double)populationSize*genrand();
         found=0; id1=0; id2=0;
         for (j=0; j<populationSize; j++) {
            if (du>=wheel[j].start && du<=wheel[j].end) { id1=wheel[j].index; found=1; break; }
         }
         if (!found) {
            id1=(int)(populationSize*genrand());
            if (id1==populationSize) id1--;
         }
         // select the second dyad
         du=(double)populationSize*genrand();
         found=0;
         for (j=0; j<populationSize; j++) {
            if (du>=wheel[j].start && du<=wheel[j].end) { id2=wheel[j].index; found=1; break; }
         }
         if (!found) {
            id2=(int)(populationSize*genrand());
            if (id2==populationSize) id2--;
         }
         // printf("id1 id2: %d %d\n",id1,id2);
      } while (id1==id2); // do not pick the same "chromosome"

      if (genrand()<maxpMutationRate) {
         // if change maxp factor other than crossing over the dyads
         do {
            newFactor=MAXP_BASE+MAXP_FACTOR*((int)(MAXP_SCALE*genrand()));
         } while (fabs(newFactor-maxpFactor[id1])<0.001);
         tmpmaxpFactor[popuCn]=newFactor;

         tmpDyad[popuCn][0].wordID   =dyad[id1][0].wordID;
         tmpDyad[popuCn][1].wordID   =dyad[id1][1].wordID; 
         tmpDyad[popuCn][2].wordID   =dyad[id1][2].wordID;
         tmpDyad[popuCn][0].wordGroup=dyad[id1][0].wordGroup;
         tmpDyad[popuCn][1].wordGroup=dyad[id1][1].wordGroup; 
         tmpDyad[popuCn][2].wordGroup=dyad[id1][2].wordGroup;
         popuCn++;
         if (popuCn==populationSize) break;

         do {
            newFactor=MAXP_BASE+MAXP_FACTOR*((int)(MAXP_SCALE*genrand()));
         } while (fabs(newFactor-maxpFactor[id2])<0.001);
         tmpmaxpFactor[popuCn]=newFactor;

         tmpDyad[popuCn][0].wordID   =dyad[id2][0].wordID;
         tmpDyad[popuCn][1].wordID   =dyad[id2][1].wordID; 
         tmpDyad[popuCn][2].wordID   =dyad[id2][2].wordID;
         tmpDyad[popuCn][0].wordGroup=dyad[id2][0].wordGroup;
         tmpDyad[popuCn][1].wordGroup=dyad[id2][1].wordGroup; 
         tmpDyad[popuCn][2].wordGroup=dyad[id2][2].wordGroup;
         popuCn++;
         if (popuCn==populationSize) break;
      }

      else {
         rand=genrand();
         if (rand>=0 && rand<1/3) {
            // replace w1 of dyad1 by w1 of dyad2
            tmpDyad[popuCn][0].wordID   =dyad[id2][0].wordID;       // replacement
            tmpDyad[popuCn][1].wordID   =dyad[id1][1].wordID; 
            tmpDyad[popuCn][2].wordID   =dyad[id1][2].wordID;
            tmpDyad[popuCn][0].wordGroup=dyad[id2][0].wordGroup;    // replacement
            tmpDyad[popuCn][1].wordGroup=dyad[id1][1].wordGroup; 
            tmpDyad[popuCn][2].wordGroup=dyad[id1][2].wordGroup;
            tmpmaxpFactor[popuCn]=maxpFactor[id1];
            popuCn++;
            if (popuCn==populationSize) break;
            
            // replace w1 of dyad2 by w1 of dyad1
            tmpDyad[popuCn][0].wordID   =dyad[id1][0].wordID;       // replacement
            tmpDyad[popuCn][1].wordID   =dyad[id2][1].wordID; 
            tmpDyad[popuCn][2].wordID   =dyad[id2][2].wordID;
            tmpDyad[popuCn][0].wordGroup=dyad[id1][0].wordGroup;    // replacement
            tmpDyad[popuCn][1].wordGroup=dyad[id2][1].wordGroup; 
            tmpDyad[popuCn][2].wordGroup=dyad[id2][2].wordGroup;
            tmpmaxpFactor[popuCn]=maxpFactor[id2];
            popuCn++;
            if (popuCn==populationSize) break;
         }
         else if (rand>=1/3 && rand<2/3) {
            // replace spacer of dyad1 by spacer of dyad2
            tmpDyad[popuCn][0].wordID   =dyad[id1][0].wordID;
            tmpDyad[popuCn][1].wordID   =dyad[id2][1].wordID;       // replacement
            tmpDyad[popuCn][2].wordID   =dyad[id1][2].wordID;
            tmpDyad[popuCn][0].wordGroup=dyad[id1][0].wordGroup;
            tmpDyad[popuCn][1].wordGroup=dyad[id2][1].wordGroup;    // replacement 
            tmpDyad[popuCn][2].wordGroup=dyad[id1][2].wordGroup;
            tmpmaxpFactor[popuCn]=maxpFactor[id1];
            popuCn++;
            if (popuCn==populationSize) break;
   
            // replace spacer dyad2 by spacer of dyad1
            tmpDyad[popuCn][0].wordID   =dyad[id2][0].wordID;
            tmpDyad[popuCn][1].wordID   =dyad[id1][1].wordID;       // replacement
            tmpDyad[popuCn][2].wordID   =dyad[id2][2].wordID;      
            tmpDyad[popuCn][0].wordGroup=dyad[id2][0].wordGroup;
            tmpDyad[popuCn][1].wordGroup=dyad[id1][1].wordGroup;    // replacement 
            tmpDyad[popuCn][2].wordGroup=dyad[id2][2].wordGroup;
            tmpmaxpFactor[popuCn]=maxpFactor[id2];
            popuCn++;
            if (popuCn==populationSize) break;
         }
         else {
            // replace w2 of dyad1 by w2 of dyad2
            tmpDyad[popuCn][0].wordID   =dyad[id1][0].wordID;       
            tmpDyad[popuCn][1].wordID   =dyad[id1][1].wordID; 
            tmpDyad[popuCn][2].wordID   =dyad[id2][2].wordID;       // replacement
            tmpDyad[popuCn][0].wordGroup=dyad[id1][0].wordGroup;   
            tmpDyad[popuCn][1].wordGroup=dyad[id1][1].wordGroup; 
            tmpDyad[popuCn][2].wordGroup=dyad[id2][2].wordGroup;    // replacement
            tmpmaxpFactor[popuCn]=maxpFactor[id1];
            popuCn++;
            if (popuCn==populationSize) break;
   
            // replace w2 of dyad2 by w2 of dyad1
            tmpDyad[popuCn][0].wordID   =dyad[id2][0].wordID;       
            tmpDyad[popuCn][1].wordID   =dyad[id2][1].wordID; 
            tmpDyad[popuCn][2].wordID   =dyad[id1][2].wordID;       // replacement
            tmpDyad[popuCn][0].wordGroup=dyad[id2][0].wordGroup;
            tmpDyad[popuCn][1].wordGroup=dyad[id2][1].wordGroup; 
            tmpDyad[popuCn][2].wordGroup=dyad[id1][2].wordGroup;    // replacement
            tmpmaxpFactor[popuCn]=maxpFactor[id2];
            popuCn++;
            if (popuCn==populationSize) break;
         }
      }
      // printf("popuCn: %d\n",popuCn);
   }

   // update the population - both pwm lengths and members
   for (i=0; i<populationSize; i++) {
      for (j=0; j<3; j++) {
         dyad[i][j].wordID   =tmpDyad[i][j].wordID;
         dyad[i][j].wordGroup=tmpDyad[i][j].wordGroup;
      }
      maxpFactor[i]=tmpmaxpFactor[i];
   }

   if (tmpDyad[0])    { free(tmpDyad[0]);     tmpDyad[0]=NULL;    }
   if (tmpDyad)       { free(tmpDyad);        tmpDyad=NULL;       }
   if (tmpmaxpFactor) { free(tmpmaxpFactor);  tmpmaxpFactor=NULL; }
}

