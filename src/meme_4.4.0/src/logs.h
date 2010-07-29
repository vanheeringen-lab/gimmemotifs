/*
 * $Id: logs.h 2417 2008-01-30 00:06:38Z cegrant $
 * 
 * $Log$
 * Revision 1.2  2005/10/06 22:52:19  nadya
 * redefine cast because some SUN compilers cant reduce expression inside define
 *
 * Revision 1.1.1.1  2005/07/29 18:41:32  nadya
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

#ifndef logs_h
# define logs_h
#include "macros.h"
#ifdef DMALLOC
#include "dmalloc.h"
#endif

/* table of log(x) for 0 < x <= 2 */

#define log_precision 1.0e5
EXTERN double log_table[2*(int)log_precision+2];	/* leave room for n+1 */

/* log lookup function; use linear interpolation */
#define LOGL_Y(x) 	((x) * log_precision)
#define LOGL_I(x)  	((int) LOGL_Y(x))
#define LOGL_LOW(x) 	(log_table[LOGL_I(x)])
#define LOGL_HI(x) 	(log_table[LOGL_I(x) + 1])
#define LOGL(x)        (LOGL_LOW(x) + (LOGL_Y(x) - LOGL_I(x)) *        \
                               (LOGL_HI(x) - LOGL_LOW(x)))


#define exp_precision 1.0e3
EXTERN double exp_table[(int)BITS*(int)exp_precision+2];	/* leave room for n+1 */

/* exp lookup function; use linear interpolation */
#define EXPL_Y(x) 	(-(x) * exp_precision)
#define EXPL_I(x)  	((int) EXPL_Y(x))
#define EXPL_LOW(x) 	(exp_table[EXPL_I(x)])
#define EXPL_HI(x) 	(exp_table[EXPL_I(x) + 1])
#define EXPL(x)        (EXPL_LOW(x) + (EXPL_Y(x) - EXPL_I(x)) *        \
                               (EXPL_HI(x) - EXPL_LOW(x)))

/* used for summing logarithms:  log(x + y) where log(x) and log(y) are avail.*/
#define LOGL_SUM1(logx,logy)                                           \
  (((((logy)<=LOGZERO) || ((logx)-(logy))>BITS)) ?                     \
  (logx) : (logx) + LOGL( 1 + EXPL((logy) - (logx) ) ) )
#define LOGL_SUM(logx, logy) ( ( (logx) > (logy) ) ?                    \
  LOGL_SUM1( (logx), (logy) ) : LOGL_SUM1( (logy), (logx) ) )

/* function prototypes */
extern void init_log(void);
extern void init_exp(void);

#endif

