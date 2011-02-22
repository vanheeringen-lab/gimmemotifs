/********************************************************************
 * FILE: fisher_exact.c
 * AUTHOR: Robert McLeay
 * CREATE DATE: 9/10/2008
 * PROJECT: MEME suite
 * COPYRIGHT: 2008, Robert McLeay
 *
 * This is a log space version of the fisher exact test.
 * It does not use a gamma table.
 * Algorithmic order is O( max(b-a ; d-c)*n + n)
 ********************************************************************/

#include "fisher_exact.h"
#include <math.h>
#include <stdlib.h>

double* log_factorial;
double fet(int a, int b, int c, int d);

void fisher_exact_init(int len) {
	double* p;
	int i;

	p = malloc(sizeof(double)*(len+1));
	p[0] = 0;
	for (i=1;i<=len;i++) {
		p[i] = log(i) + p[i-1];
		//fprintf(stderr, "%i: %g\n", i, p[i]);
	}
	log_factorial = p;
}

double* fisher_exact_get_log_factorials() {
	return log_factorial;
}

double fet(int a, int b, int c, int d) {
	return  exp(log_factorial[a+b] + log_factorial[c+d] + log_factorial[a+c] + log_factorial[b+d]
	            - (log_factorial[a+b+c+d] + log_factorial[a] + log_factorial[b]
	               + log_factorial[c] + log_factorial[d]));
}

void fisher_exact(int a, //x[0,0]
				  int b, //x[0,1]
				  int c, //x[1,0]
				  int d, //x[1,1]
				  double* two, //two-tailed p-value (out)
				  double* left, //one-tailed left p-value (out)
				  double* right, //one-tailed right p-value (out)
				  double* p //exact p-value of this matrix (out)
				  ) {

		int tmpa,tmpb,tmpc,tmpd; //used for modifying the tables to be more 'extreme'.

		double tmpp;

		*two = *left = *right = *p = tmpp = 0;

		*p = fet(a,b,c,d);

		tmpa = a;
		tmpb = b;
		tmpc = c;
		tmpd = d;
		*left = *p; //for the current state
		while (tmpa != 0 && tmpd != 0) {
			tmpa--;
			tmpb++;
			tmpc++;
			tmpd--;
			tmpp = fet(tmpa,tmpb,tmpc,tmpd);
			if (tmpp <= *p)
				*left += tmpp;
		}

		//reset to go the other way (right)
		tmpa = a;
		tmpb = b;
		tmpc = c;
		tmpd = d;
		*two = *left;
		while (tmpb != 0 && tmpc != 0) {
			tmpa++;
			tmpb--;
			tmpc--;
			tmpd++;
			tmpp = fet(tmpa,tmpb,tmpc,tmpd);
			if (tmpp <= *p)
				*two += tmpp;
		}

		*right = 1 - *left + *p;
}

void fisher_exact_destruct() {
	free(log_factorial);
}
