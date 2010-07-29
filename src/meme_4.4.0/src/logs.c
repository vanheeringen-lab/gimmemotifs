/*
 * $Id: logs.c 355 2005-10-25 19:06:39Z nadya $
 * 
 * $Log$
 * Revision 1.2  2005/10/25 19:06:39  nadya
 * rm old macro for Header, all info is taken care of by Id and Log.
 *
 * Revision 1.1.1.1  2005/07/29 00:23:24  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/***********************************************************************
*                                                                      *
*       MEME++                                                         *
*       Copyright 1994, The Regents of the University of California    *
*       Author: Timothy L. Bailey                                      *
*                                                                      *
***********************************************************************/
/* logs.c */

#include "meme.h"
#include "logs.h"

static BOOLEAN init_log_done = FALSE;
static BOOLEAN init_exp_done = FALSE;

/**********************************************************************/
/*
	init_log

	Setup lookup table for log(x), 0 < x <= 2
*/
/**********************************************************************/
extern void init_log(void)
{
  int i;
  double x;

  for (i=0; !init_log_done && i<=2*log_precision+1; i++) {
    x = (double) i/log_precision;
    log_table[i] = LOG(x);
    /*printf("%d %f\n", i, log_table[i]);*/
  }
} /* init_log */

/**********************************************************************/
/*
	init_exp

	Setup lookup table for exp(x), -BITS <= x < 0 
*/
/**********************************************************************/
extern void init_exp(void)
{
  int i;
  double x;

  for (i=0; !init_exp_done && i<=BITS*exp_precision+1; i++) {
    x = -i/exp_precision;
    exp_table[i] = exp(x);
    /*printf("%d %f %g\n", i, x, exp_table[i]);*/
  }
} /* init_exp */

