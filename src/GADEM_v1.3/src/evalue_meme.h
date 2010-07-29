/* ---------------------------------------------------------------------- */
/*  E_value                                                               */
/*                                                                        */
/*  Calculated the entropy ala MEME                                       */
/*                                                                        */
/*  This code was lifted directly from MEME v. 3.5.4 and consolidated     */
/*  for understanding and modularity                                      */
/*                                                                        */
/* ---------------------------------------------------------------------- */
/*                                                                        */
/*    MEME++                                                              */
/*    Author: Timothy L. Bailey                                           */
/*                                                                        */
/*    Copyright  (c)  1994-2006  The  Regents  of the  University of      */ 
/*    California.  All  Rights  Reserved.                                 */                           
/*                                                                        */                                                                  
/*    Permission  to use,  copy,  modify,  and  distribute  any part      */
/*    of this  software for  educational,  research  and  non-profit      */ 
/*    purposes,  without  fee,  and  without a written  agreement is      */ 
/*    hereby  granted,  provided  that the  above  copyright notice,      */ 
/*    this paragraph  and the following  three  paragraphs appear in      */ 
/*    all copies.                                                         */                                                   
/*                                                                        */                                                                  
/*    Those  desiring to  incorporate this  software into commercial      */
/*    products  or use for  commercial  purposes  should contact the      */     
/*    Technology  Transfer  Office,  University of California,   San      */ 
/*    Diego,  9500 Gilman Drive,  La Jolla,  California, 92093-0910,      */         
/*    Phone: (858) 534-5815.                                              */                                        
/*                                                                        */                                                                  
/*    IN  NO  EVENT  SHALL THE  UNIVERSITY  OF CALIFORNIA  BE LIABLE      */ 
/*    TO  ANY  PARTY FOR  DIRECT,  INDIRECT, SPECIAL, INCIDENTAL, OR      */     
/*    CONSEQUENTIAL  DAMAGES,  INCLUDING  LOST PROFITS, ARISING  OUT      */ 
/*    OF  THE  USE  OF  THIS  SOFTWARE,  EVEN  IF THE UNIVERSITY  OF      */ 
/*    CALIFORNIA  HAS  BEEN  ADVISED  OF  THE  POSSIBILITY  OF  SUCH      */ 
/*    DAMAGE.                                                             */                                                       
/*                                                                        */                                                                  
/*    THE SOFTWARE  PROVIDED HEREUNDER IS ON N  "AS IS" BASIS,  AND       */ 
/*    THE  UNIVERSITY OF CALIFORNIA  HAS  NO OBLIGATIONS  TO PROVIDE      */         
/*    MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.      */  
/*    THE UNIVERSITY  OF CALIFORNIA  MAKES  NO  REPRESENTATIONS  AND      */     
/*    EXTENDS  NO  WARRANTIES  OF  ANY  KIND,  EITHER  EXPRESSED  OR      */ 
/*    IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES      */
/*    OF  MERCHANTABILITY  OR  FITNESS FOR A  PARTICULAR PURPOSE, OR      */ 
/*    THAT  THE USE  OF THE MATERIAL  WILL NOT  INFRINGE ANY PATENT,      */         
/*    TRADEMARK OR OTHER RIGHTS.                                          */                                    
/*                                                                        */
/* ---------------------------------------------------------------------- */

#ifndef _EVALUE_MEME_H_
#define _EVALUE_MEME_H_

#include <math.h>

/* ---------------------------------------------------------------------- */
/* From MEME: user.h                                                      */
/* ---------------------------------------------------------------------- */

#define LLR_RANGE 200            /* range of scaled LLR statistic */

/* ---------------------------------------------------------------------- */
/* From MEME: mtype.h                                                     */
/* ---------------------------------------------------------------------- */

typedef enum {Tcm, Oops, Zoops} MTYPE;

/* ---------------------------------------------------------------------- */
/* From MEME: macros.h                                                    */
/* ---------------------------------------------------------------------- */


typedef void *malloc_t;


#define BIG HUGE_VAL

typedef int BOOLEAN;
#define FALSE 0
#define TRUE 1

#define EPSILON 1e-200         /* smallest LOG2(X) = -664.385 */

#define BITS 64.0
#define LOGZERO -1e100			/* log of number close to zero */
#define LOGZEROI -2147483647		/* log of number close  zero (integer)*/

#define exp10(X) pow(10.0, (X))
	
/* macro to round to the nearest int value, except halfway cases are rounded
   to the int value larger in magnitude.
*/
#define NINT(x) ((int)((x) >= 0 ? ((x) + 0.5) : ((x) - 0.5)))


/* max and min for all seasons */
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

/*
  round x to d significant digits and put in y
*/
#define RNDDIG 14
#define RND(x, d, y) {							               \
  if (x > 0) {                                           \
    double _z_ = exp10(ceil((d)-1-log10(x)));				\
    y = rint(_z_*(x))/_z_;						               \
  } else if (x < 0) {							               \
    double _z_ = exp10(ceil((d)-1-log10(-x)));				\
    y = -rint(_z_*(-x))/_z_;						\
  } else {								\
    y = 0;								\
  }									\
}


/* 
   get the exponent and mantissa of large numbers expressed as log10(x)
   for printing with prec digits after the decimal point in mantissa
*/
#define exp10_logx(logx, m, e, prec) { \
  (e) = floor(logx); \
  (m) = exp10((logx) - (e)); \
  if (m+(.5*pow(10,-prec)) >= 10) { (m) = 1; (e) += 1;} \
}


/* make a better free function */
#define myfree(P) {							\
  if (P) free((char *) (P));						\
}
#define dbmalloc(P, N, T) {						\
  (P) = (N) > 0 ? (T *) malloc((N)*sizeof(T)) : NULL;			\
}
#define myrealloc(P, OLDP, N, T) {					\
  (P) = (T *) realloc((malloc_t)(OLDP), (N)*sizeof(T));			\
}


/* dynamically create or grow an array; P = pointer, N = new size, T = type */
#define Resize(P,N,T) { 						\
  void *new_P = NULL; 							\
  int _n_ = (N);			/* avoid side effects */	\
  /*fprintf(stderr, "Resize(" #P ", " #N "(%ld), " #T ")\n", (long)_n_);*/\
  if (P) {								\
    myrealloc(new_P, (P), _n_, T);	 				\
  } else {								\
    dbmalloc(new_P, _n_, T);						\
  }									\
  if (!new_P || _n_<=0) {						\
    fprintf(stderr, "Resize(" #P ", " #N ", " #T ") failed!\n");	\
    fprintf(stderr, #N" = %ld\n", (long)_n_);				\
    /*Crash*/; 								\
    exit(1);								\
  } 									\
  (P) = (T *) new_P; 							\
}


/* macros to create and destroy an r x c array; 
   avoids problems with multiply subscripted C array elements 
   behaving differently when passed to subroutines or assigned locally 
*/
#define create_2array(v, typ, r, c) {					\
  int _i, _ir=(r), _ic=(c);		/* avoid side effects */	\
  /*fprintf(stderr, */							\
  /*  "create_2array(" #v ", " #typ ", " #r "(%d), " #c ")\n", r);*/	\
  dbmalloc((v), _ir,  typ *);						\
  if (!(v)) {fprintf(stderr, "malloc failed 1\n"); exit(1);}		\
  for (_i=0; _i<_ir; _i++) {						\
    typ *tmp;								\
    (v)[_i] = 0;							\
    dbmalloc(tmp, _ic, typ);						\
    (v)[_i] = tmp;							\
    if (!(v)[_i]) {printf("malloc failed 2\n");}			\
    if (!(v)[_i]) {fprintf(stderr, "malloc failed 2\n"); exit(1);}	\
  }									\
}

#define free_2array(v, r) {						\
  int _i, _ir=(r);			/* avoid side effects */	\
  /*fprintf(stderr, "free_2array(" #v ", " #r "(%d))\n", r);	*/	\
  for (_i=0; _i<_ir; _i++) myfree((v)[_i]);				\
  myfree(v);								\
}

/* ---------------------------------------------------------------------- */
/* From MEME: logs.h                                                      */
/* ---------------------------------------------------------------------- */

#define log_precision 1.0e5
#define exp_precision 1.0e3

#define LOG(X) log((double)((X) + EPSILON))  
#define LOG2(X) (LOG(X)/0.6931471805599452862)

/* log lookup function; use linear interpolation */
#define LOGL_Y(x) 	((x) * log_precision)
#define LOGL_I(x)  	((int) LOGL_Y(x))
#define LOGL_LOW(x) 	(log_table[LOGL_I(x)])
#define LOGL_HI(x) 	(log_table[LOGL_I(x) + 1])
#define LOGL(x)       (_logl_x = (x), LOGL_LOW(_logl_x) +     		\
                          (LOGL_Y(_logl_x) - LOGL_I(_logl_x)) *         \
                          (LOGL_HI(_logl_x) - LOGL_LOW(_logl_x)))


/* exp lookup function; use linear interpolation */
#define EXPL_Y(x) 	(-(x) * exp_precision)
#define EXPL_I(x)  	((int) EXPL_Y(x))
#define EXPL_LOW(x) 	(exp_table[EXPL_I(x)])
#define EXPL_HI(x) 	(exp_table[EXPL_I(x) + 1])
#define EXPL(x)       (_expl_x = (x), EXPL_LOW(_expl_x) +     		\
                          (EXPL_Y(_expl_x) - EXPL_I(_expl_x)) *         \
                          (EXPL_HI(_expl_x) - EXPL_LOW(_expl_x)))

/* used for summing logarithms:  log(x + y) where log(x) and log(y) are avail.*/
#define LOGL_SUM1(logx,logy)                                            \
  (((((logy)<=LOGZERO) || ((logx)-(logy))>BITS)) ?                      \
  (logx) : (logx) + LOGL( 1 + EXPL((logy) - (logx) ) ) )
#define LOGL_SUM(logx, logy) ( ( (logx) > (logy) ) ?                    \
  LOGL_SUM1( (logx), (logy) ) : LOGL_SUM1( (logy), (logx) ) )

#endif
