
#include <stdlib.h>

void copy_pwm(double **pwm1,double **pwm2,int pwmLen) {

   register int i,j;

   for (i=0; i<pwmLen; i++) {
      for (j=0; j<4; j++) pwm2[i][j]=pwm1[i][j];
   }
}

