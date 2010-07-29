#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gadem.h"
#include "defines.h"

int read_pwm0(char *fileName,double **pwm0) {

   register int m,n;
   int numCol,numRow;
   FILE *fp;

   fp=fopen(fileName,"r");
   if (!fp) { perror(fileName); exit(0); }

   fscanf(fp,"%d %d",&numRow,&numCol);
   if (numRow !=4 || numCol<2) { 
      printf("\n\nError:\n");
      printf("PMW format: numRow(4) numCol\n");
      printf("then the 4 x numCol matrix\n");
      printf("Note: gadem does not check if the number of actual columns matches the specified number.\n");
      printf("\nexample:\n"); 
      printf("      4  12\n");
      printf("      0  3       0       2       5       0       0       16      0       0       1       5\n");
      printf("      7  5       3       3       1       0       0       0       16      0       5       6\n");
      printf("      5  4       6       11      7       0       15      0       0       16      0       3\n");
      printf("      4  4       7       0       3       16      1       0       0       0       10      2\n");
      printf("\n                                  OR\n");
      printf("      4  12\n");
      printf("      0.000      0.188   0.000   0.125   0.312   0.000   0.000   1.000   0.000   0.000   0.062   0.312\n");
      printf("      0.438      0.312   0.188   0.188   0.062   0.000   0.000   0.000   1.000   0.000   0.312   0.375\n");
      printf("      0.312      0.250   0.375   0.688   0.438   0.000   0.938   0.000   0.000   1.000   0.000   0.188\n");
      printf("      0.250      0.250   0.438   0.000   0.188   1.000   0.062   0.000   0.000   0.000   0.625   0.125\n");
      printf("\n");

      exit(0); 
   }
 
   if (numCol>MAX_PWM_LENGTH) {
       printf("\nmaximal number of columns reached.\n"); exit(0);
   }
   for (m=0; m<4; m++) {
      for (n=0; n<numCol; n++) {
         fscanf(fp,"%lf",&(pwm0[n][m]));
         if (pwm0[n][m]<0) { 
            printf("\n\nError: cell must not be negative.\n"); exit(0); 
         } 
      }
   }
   fclose(fp);

   if (numCol<2) {
       printf("\npwm contains too few columns.\n"); exit(0);
   }
  
   return (numCol);
}

