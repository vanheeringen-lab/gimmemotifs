#include "config.h"

#include <stdio.h>
#include <stdlib.h>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include <math.h>
#include "gadem.h"
#include "defines.h"

// fitness scores have been arbitrarily set to either 2.0 (uniq motifs) or 1.0 (non-unique motifs)
// all unique motifs have the same weight (2.0) and all non-unique motifs have weights (1.0)

void roulett_wheel_weight(Fitness *fitness,int populationSize,Wheel *wheel) {

   register int i;
   double totalWeight;

   totalWeight=0; for (i=0; i<populationSize; i++) totalWeight += fitness[i].value;
 
   if (totalWeight==0) { printf("Error: \n"); exit(0); }

   wheel[0].start=0.0;
   wheel[0].end  =(double)populationSize*fitness[0].value/totalWeight;
   wheel[0].index=fitness[0].index;

   for (i=1; i<populationSize; i++) {
      wheel[i].start=wheel[i-1].end;
      wheel[i].end=(double)populationSize*fitness[i].value/totalWeight+wheel[i].start;
      wheel[i].index=fitness[i].index;
   }
   // for (i=0; i<populationSize; i++) { 
   //    printf("%4d\t%5.3f\t%3d\t%5.3f\t%5.3f\n",
   //       fitness[i].index,fitness[i].value,wheel[i].index,wheel[i].start,wheel[i].end);
   // }
}

// this subroutine assigns probabilities that are proportional to fitness scores
// Note that this minimizes rather than maximizes
// The minimum is assigned the largest probability
void roulett_wheel_fitness(Fitness *fitness,int populationSize,Wheel *wheel) {

   register int i;
   double totalScore,worstScore,range;
   double *scaledScore;

   worstScore=fitness[populationSize-1].value;
   range=worstScore-fitness[0].value;

   // for (i=0; i<populationSize; i++) printf("%d\t%5.3f\n",fitness[i].index,fitness[i].value);

   if ((range<0.0001) || 
       fitness[1].value==DUMMY_FITNESS || 
       fabs(fitness[1].value-fitness[populationSize-1].value)<0.0001) {
      printf("GA converged ...\n");
      for (i=0; i<populationSize; i++) {
         wheel[i].index=fitness[i].index;
         wheel[i].start=i;
         wheel[i].end=i+1;
      }
   }
   else {
      scaledScore=alloc_double(populationSize);

      totalScore=0;
      for (i=0; i<populationSize; i++) {
         // range scale 
         // make sure the minimum is assigned the largest area
         scaledScore[i]=1.0-(fitness[i].value-fitness[0].value)/range;
         totalScore += scaledScore[i];
      }
      for (i=0; i<populationSize; i++) scaledScore[i] /= totalScore;

      //for (i=0; i<populationSize; i++) printf("%6.5f\n",scaledScore[i]); exit(0);

      wheel[0].start=0;
      wheel[0].end  =(double)populationSize*scaledScore[0];
      wheel[0].index=fitness[0].index;
      //if (fitness[0].value<0) wheel[0].I_evalue='0';
      //else                    wheel[0].I_evalue='1';

      for (i=1; i<populationSize; i++) {
         wheel[i].start=wheel[i-1].end;
         wheel[i].end=(double)populationSize*scaledScore[i]+wheel[i].start;
         wheel[i].index=fitness[i].index;
         //if (fitness[i].value<0) wheel[i].I_evalue='0';
         //else                    wheel[i].I_evalue='1';
      }
      if (scaledScore) { free(scaledScore); scaledScore = NULL; }
   }
   //for (i=0; i<populationSize; i++) 
   //  printf("%4d\t%5.3f\t%3d\t%c\t%5.3f\t%5.3f\n",
   //  fitness[i].index,fitness[i].value,wheel[i].index,wheel[i].I_evalue,wheel[i].start,wheel[i].end);
}

void roulett_wheel_rank(Fitness *fitness,int populationSize,Wheel *wheel) {

   register int i;
   int sum;
   double *weight;
 
   weight=alloc_double(populationSize);

   sum=0; for (i=1; i<populationSize+1; i++) sum += i;
   for (i=0; i<populationSize; i++) weight[i]=(double)(populationSize-i)/(double)sum;
 
   wheel[0].start=0.0;
   wheel[0].end  =(double)populationSize*weight[0];
   wheel[0].index=fitness[0].index;

   for (i=1; i<populationSize; i++) {
      wheel[i].start=wheel[i-1].end;
      wheel[i].end=(double)populationSize*weight[i]+wheel[i].start;
      wheel[i].index=fitness[i].index;
   }
   //for (i=0; i<populationSize; i++) { 
   //   printf("%4d\t%8.3f\t%3d\t%8.6f\t%8.6f\n",
   //      fitness[i].index,fitness[i].value,wheel[i].index,wheel[i].start,wheel[i].end);
   //}
   if (weight) { free(weight); weight=NULL; }
}
