#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "defines.h"
#include "gadem.h"

double *alloc_double(int );

void pwm_profile(double **pwm,int pwmLen,char *dyad) {

   int i,j;
   for (i=0; i<pwmLen; i++) {
      if (pwm[i][0]==1 && pwm[i][1]==1 && pwm[i][2]==1 && pwm[i][3]==1) dyad[i]='n';
      else {
         for (j=0; j<4; j++) {
            if (pwm[i][j]==CELL_MAJOR) {
               switch (j) {
                  case 0:  dyad[i]='a'; break;
                  case 1:  dyad[i]='c'; break;
                  case 2:  dyad[i]='g'; break;
                  case 3:  dyad[i]='t'; break;
                  default: dyad[i]='n'; break;
               }
               break;
            }
         }
      }
   }
   dyad[pwmLen]='\0';
}

void consensus_pwm(double **pwm,int pwmLen,char *consensus) {

   register int i,j;
   ColCons *col;

   col=(ColCons *)calloc(4,sizeof(ColCons));

   for (i=0; i<pwmLen; i++) {
      for(j=0; j<4; j++) { col[j].v=pwm[i][j]; col[j].i=j; }
      sort_vector(col,4);
      if (col[0].v>0.6 && col[0].v>=2*col[1].v) {
         switch(col[0].i) {
            case 0: consensus[i]='A'; break;
            case 1: consensus[i]='C'; break;
            case 2: consensus[i]='G'; break;
            case 3: consensus[i]='T'; break;
            default: break;
         }
      }
      else if ((col[0].v+col[1].v>=0.70 && col[1].v<0.5) || (col[0].v+col[1].v>=0.5 && col[1].v>0.2 && col[2].v<0.1)) {
         if      ((col[0].i==0 && col[1].i==1)||(col[0].i==1 && col[1].i==0)) consensus[i]='m';
         else if ((col[0].i==0 && col[1].i==2)||(col[0].i==2 && col[1].i==0)) consensus[i]='r';
         else if ((col[0].i==0 && col[1].i==3)||(col[0].i==3 && col[1].i==0)) consensus[i]='w';
         else if ((col[0].i==1 && col[1].i==2)||(col[0].i==2 && col[1].i==1)) consensus[i]='s';
         else if ((col[0].i==1 && col[1].i==3)||(col[0].i==3 && col[1].i==1)) consensus[i]='y';
         else if ((col[0].i==2 && col[1].i==3)||(col[0].i==3 && col[1].i==2)) consensus[i]='k';
         else { }
      }
      else if (col[3].v<=0.075) {
         switch(col[3].i) {
            case 0: consensus[i]='b'; break; /* c, g, or t */
            case 1: consensus[i]='d'; break; /* a, g, or t */
            case 2: consensus[i]='h'; break; /* a, c, or t */
            case 3: consensus[i]='v'; break; /* a, c, or g */
            default: break;
         }
      }
      else consensus[i]='n';
   }
   consensus[pwmLen]='\0';
   if (col) { free(col); col=NULL; }
}

