/********************************************************************
 * FILE: fisher_exact.h
 * AUTHOR: Robert McLeay
 * CREATE DATE: 9/10/2008
 * PROJECT: MEME suite
 * COPYRIGHT: 2008, Robert McLeay
 *
 * This is a log space version of the fisher exact test.
 * It does not use a gamma table.
 * Algorithmic order is O( max(b-a ; d-c)*n + n)
 ********************************************************************/

#ifndef FISHER_EXACT_H_
#define FISHER_EXACT_H_

void fisher_exact_init(int len);
void fisher_exact_destruct();
void fisher_exact(int a, //x[0,0]
				  int b, //x[0,1]
				  int c, //x[1,0]
				  int d, //x[1,1]
				  double* two, //two-tailed p-value (out)
				  double* left, //one-tailed left p-value (out)
				  double* right, //one-tailed right p-value (out)
				  double* p //exact p-value of this matrix (out)
				  );

double* fisher_exact_get_log_factorials();

#endif /* FISHER_EXACT_H_ */
