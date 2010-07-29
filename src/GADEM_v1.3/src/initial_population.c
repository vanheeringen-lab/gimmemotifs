#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include <math.h>
#include "defines.h"
#include "gadem.h"
#include "random.h"

void initialisation(Chrs **dyad,int populationSize,int numWordGroup,Words *word,
   int minSpaceWidth,int maxSpaceWidth,double *maxpFactor) {

   register int j;
   int group,popuCn;
   double rand;

   // initialize spaced dyads and the associated maxp factor - a population of GA "chromosomes"
   popuCn=0; 
   do { 
      // initialize w1 of dyad: w1-spacer-w2 
      group=(int)(numWordGroup*genrand());
      if (group==numWordGroup) group--;
      dyad[popuCn][0].wordGroup=group;

      rand=genrand();
      dyad[popuCn][0].wordID=0;  // update right below
      for (j=0; j<word[group].count; j++) {
         if (rand>=word[group].prob_sta[j] && rand<word[group].prob_end[j]) { 
            dyad[popuCn][0].wordID=j; break; 
         }
      }

      // choose a spacer width
      dyad[popuCn][1].wordGroup=-1;  // dummy
      dyad[popuCn][1].wordID=minSpaceWidth+(int)((maxSpaceWidth-minSpaceWidth+1)*genrand());

      // initialize w2 of dyad w1-spacer-w2 
      group=(int)(numWordGroup*genrand());
      if (group==numWordGroup) group--;
      dyad[popuCn][2].wordGroup=group;

      dyad[popuCn][2].wordID=0; 
      rand=genrand();
      for (j=0; j<word[group].count; j++) {
         if (rand>=word[group].prob_sta[j] && rand<word[group].prob_end[j]) { 
            dyad[popuCn][2].wordID=j; break; 
         }
      }

      // randomly assign a value to each maxp factor
      maxpFactor[popuCn]=MAXP_BASE+MAXP_FACTOR*((int)(MAXP_SCALE*genrand()));

      popuCn++;
   } while (popuCn<populationSize);
}

void dyad_to_pwm(Words *word,int populationSize,Chrs **dyad,double ***pwm,int *pwmLen) {

   register int i,j,k;
   int len0,len1,len2;
   double major,minor;

   minor=CELL_MINOR;
   major=CELL_MAJOR;
   for (i=0; i<populationSize; i++) {
      for (j=0; j<MAX_PWM_LENGTH; j++) {
         for (k=0; k<4; k++) pwm[i][j][k]=minor; 
      } 
   }

   for (i=0; i<populationSize; i++) {
      //word 1 of spaced dyad
      len0=strlen(word[dyad[i][0].wordGroup].s1[dyad[i][0].wordID]); 
      for (j=0; j<len0; j++) {
         switch (word[dyad[i][0].wordGroup].s1[dyad[i][0].wordID][j]) {
            case 'a': pwm[i][j][0]=major; break; 
            case 'c': pwm[i][j][1]=major; break; 
            case 'g': pwm[i][j][2]=major; break; 
            case 't': pwm[i][j][3]=major; break; 
            case 'n': { pwm[i][j][0]=1; pwm[i][j][1]=1; pwm[i][j][2]=1; pwm[i][j][3]=1; break; }
            default: break;
         } 
      }

      //spacer  of spaced dyad
      len1=dyad[i][1].wordID;
      for (j=0; j<len1; j++) {
         for (k=0; k<4; k++) pwm[i][len0+j][k]=1; 
      }

      //word 2 of spaced dyad
      len2=strlen(word[dyad[i][2].wordGroup].s1[dyad[i][2].wordID]);
      for (j=0; j<len2; j++) {
         switch (word[dyad[i][2].wordGroup].s1[dyad[i][2].wordID][j]) {
            case 'a': pwm[i][j+len0+len1][0]=major; break; 
            case 'c': pwm[i][j+len0+len1][1]=major; break; 
            case 'g': pwm[i][j+len0+len1][2]=major; break; 
            case 't': pwm[i][j+len0+len1][3]=major; break; 
            case 'n': { pwm[i][j+len0+len1][0]=1; pwm[i][j+len0+len1][1]=1; pwm[i][j+len0+len1][2]=1; pwm[i][j+len0+len1][3]=1; break; }
            default: break;
         } 
      }
      pwmLen[i]=len0+len1+len2;
   }

   /*---------------------------------------------------------------------------
   for (i=0; i<populationSize; i++) {
      for (k=0; k<4; k++) {
         for (j=0; j<pwmLen[i]; j++) printf("%4.2f ",pwm[i][j][k]); printf("\n");
      } printf("\n");
   }
   -----------------------------------------------------------------------------*/

   /*-----------------------for testing during development------------------------
   //ggtcannntgacc
   pwmLen[0]=13;
  
   for (i=0; i<pwmLen[0]; i++) { for (j=0; j<4; j++) pwm[0][i][j]=0; }

   pwm[0][0][2]=1; 
   pwm[0][1][2]=1; 
   pwm[0][2][3]=1;  
   pwm[0][3][1]=1; 
   pwm[0][4][0]=1; 
   pwm[0][5][0]=1; pwm[0][5][1]=1; pwm[0][5][2]=1; pwm[0][5][3]=1;
   pwm[0][6][0]=1; pwm[0][6][1]=1; pwm[0][6][2]=1; pwm[0][6][3]=1;
   pwm[0][7][0]=1; pwm[0][7][1]=1; pwm[0][7][2]=1; pwm[0][7][3]=1;
   pwm[0][8][3]=1; 
   pwm[0][9][2]=1;
   pwm[0][10][0]=1;
   pwm[0][11][1]=1;
   pwm[0][12][1]=1;

   for (j=0; j<4; j++) {
      for (i=0; i<pwmLen[0]; i++) printf("%1.0f ",pwm[0][i][j]); printf("\n");
   } 
   -------------------------------------------------------------------------*/


   /*-------------------------------for testing----------------------------------
   //atttgcatnwnaawr
   pwmLen[0]=15;
  
   for (i=0; i<pwmLen[0]; i++) { for (j=0; j<4; j++) pwm[0][i][j]=0; }

   pwm[0][0][0]=1;
   pwm[0][1][3]=1; 
   pwm[0][2][3]=1;  
   pwm[0][3][3]=1; 
   pwm[0][4][2]=1; 
   pwm[0][5][1]=1; 
   pwm[0][6][0]=1;
   pwm[0][7][3]=1;
   pwm[0][8][0]=1; pwm[0][8][1]=1; pwm[0][8][2]=1; pwm[0][8][3]=1;
   pwm[0][9][1]=1; pwm[0][9][3]=1;
   pwm[0][10][0]=1; pwm[0][10][1]=1; pwm[0][10][2]=1; pwm[0][10][3]=1;
   pwm[0][11][0]=1;
   pwm[0][12][0]=1;
   pwm[0][13][0]=1; pwm[0][13][3]=1;
   pwm[0][14][0]=1; pwm[0][14][2]=1;

   for (j=0; j<4; j++) {
      for (i=0; i<pwmLen[0]; i++) printf("%1.0f ",pwm[0][i][j]); printf("\n");
   } 
   -------------------------------------------------------------------------*/

   /*-------------------------------for testing----------------------------------
   //rrrcatgyyyrrrcatgyyy
   pwmLen[0]=20;
  
   for (i=0; i<pwmLen[0]; i++) { for (j=0; j<4; j++) pwm[0][i][j]=0; }

   pwm[0][0][0]=1;  pwm[0][0][2]=1;
   pwm[0][1][0]=1;  pwm[0][1][2]=1;
   pwm[0][2][0]=1;  pwm[0][2][2]=1;
   pwm[0][3][1]=1; 
   pwm[0][4][0]=1;
   pwm[0][5][3]=1; 
   pwm[0][6][2]=1;
   pwm[0][7][1]=1;  pwm[0][7][3]=1;
   pwm[0][8][1]=1;  pwm[0][8][3]=1; 
   pwm[0][9][1]=1;  pwm[0][9][3]=1;
   pwm[0][10][0]=1; pwm[0][10][2]=1;
   pwm[0][11][0]=1; pwm[0][11][2]=1;
   pwm[0][12][0]=1; pwm[0][12][2]=1;
   pwm[0][13][1]=1; 
   pwm[0][14][0]=1;
   pwm[0][15][3]=1; 
   pwm[0][16][2]=1;
   pwm[0][17][1]=1; pwm[0][17][3]=1;
   pwm[0][18][1]=1; pwm[0][18][3]=1; 
   pwm[0][19][1]=1; pwm[0][19][3]=1;

   for (j=0; j<4; j++) {
      for (i=0; i<pwmLen[0]; i++) printf("%1.0f ",pwm[0][i][j]); printf("\n");
   } 
   -------------------------------------------------------------------------*/

#ifdef DEBUG
   for (i=0; i<populationSize; i++) {
      for (k=0; k<4; k++) {
         for (j=0; j<pwmLen[i]; j++) printf("%2d",pwm[i][j][k]); printf("\n");
      } 
   }
#endif
}

