#include <stdlib.h>
#include <math.h>

// return the maximal difference between two pwms
double check_convergence(double **pwm1,double **pwm2,int pwmLen){

    register int i,j;
    double value,dif;

    dif=fabs(pwm1[0][0]-pwm2[0][0]);
    for (i=0; i<pwmLen; i++) {
       for (j=0; j<4; j++){
          value=fabs(pwm1[i][j]-pwm2[i][j]);
          if (value>dif) dif=value;
       }
    }
    return(dif);
}

