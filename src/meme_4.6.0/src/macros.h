/*
 * $Id: macros.h 1048 2006-07-06 20:07:44Z cegrant $
 *
 * $Log$
 * Revision 1.1  2005/07/29 18:43:23  nadya
 * Initial revision
 *
 */

/*
	Handy macros for general C programming
*/

/* 7-14-99 tlb; add LOGEV */
/* 7-12-99 tlb; add DELOG and INT_DELOG */
/* 7-09-99 tlb; fix documentation for exp10_logx; added LOGZERO check to LOG_SUM
*/
#include <unistd.h>

#ifndef macros_h
#define macros_h

/*#define MALLOC_DEBUG	*/	/* turn this on for debugging */

#ifdef sgi4d
#define _XOPEN_SOURCE		/* needed for isnan in math.h on sgi */
extern double   cbrt(double); /* _XOPEN_SOURCE turns of definition cbrt */
#endif
#include <math.h>
#ifdef sgi4d
  extern double rint();			/* missing from math.h ? */
#endif
#define exp10(X) pow(10.0, (X))
#ifndef ibmrs6000
  #define log10(X) (log(X)/2.30258509299405)
#endif

#ifdef UNIX

#ifndef SUNOS_CC
  /* this is defined in malloc.h which is in stdlib.h when compiler CC is
     used under SunOS; otherwise we must declare it
  */
    typedef void *malloc_t;
#endif

#else
  /* VAX stuff */
    typedef void *malloc_t;

#endif

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#ifdef sgi4d
  extern char *strdup(const char *);		/* missing for some reason */
#endif
#include <strings.h>
#include <stdlib.h>
#if defined(sgi4d) || defined(ibmrs6000) || defined(Linux)
  #include <time.h>
#else
  #include <sys/time.h>
#endif
#include <sys/resource.h>
#ifdef crayt3e
#include <fp.h>
#endif

/* random functions */
#ifdef __cplusplus
extern "C" {
#endif

/* flush functions */
extern int	filbuf(void);
extern int	flsbuf(void);

/* random number functions */
double drand48(void);
void srand48(long seedval);

#ifdef __cplusplus
}
#endif

/* EXTERN and DEXTERN can be used to allow global variables
   to be defined in a .h file.  The .h file can be included
   wherever the variables are used.  It should be preceeded
   by #define DEFINE_GLOBALS in exactly one file, usually the
   main program file.  This causes space to be allocated for
   the variable only once.
*/
#ifdef DEFINE_GLOBALS
#define EXTERN
#define DEXTERN(A,B,C) A B = C
#else
#define EXTERN extern
#define DEXTERN(A,B,C) extern A B
#endif

/*
  Cause a crash
*/
DEXTERN(char *, __crash_x__, NULL);
#define Crash (*__crash_x__ = 'x')

/*
  rounding
*/
/* macro to round to the nearest int value, except halfway cases are rounded
   to the int value larger in magnitude.
*/
#define NINT(x) ((int)((x) >= 0 ? ((x) + 0.5) : ((x) - 0.5)))
/* macro to round to the nearest int value, except halfway cases are rounded
   to the int value smaller in magnitude.
*/
#define _H_ 0.499999999999999
#define NINTL(x) ((int)((x) >= 0 ? ((x) + _H_) : ((x) - _H_)))
/*
  round x to d significant digits and put in y
*/
#define RNDDIG 14
#define RND(x, d, y) {							\
  if (x > 0) {								\
    double _z_ = exp10(ceil((d)-1-log10(x)));				\
    y = rint(_z_*(x))/_z_;						\
  } else if (x < 0) {							\
    double _z_ = exp10(ceil((d)-1-log10(-x)));				\
    y = -rint(_z_*(-x))/_z_;						\
  } else {								\
    y = 0;								\
  }									\
}

/* the largest and smallest double precision numbers */
#define BIG HUGE_VAL
#define LITTLE -BIG
#define MAXPOSLONG 0x7FFFFFFF
#define MAXNEGLONG 0x80000000

/* handy constants */
#define Log2 0.693147				/* log(2) */
#define MAXASCII 256 				/* number of ASCII codes */

/* max and min for all seasons */
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

/* swapping any two items x and y of type "type"*/
#define SWAP(type, x, y) {type tmp; tmp = (x); (x) = (y); (y) = (tmp);}

/* Used for preventing overflow with logarithms */
#define EPSILON 1e-200		/* smallest LOG2(X) = -664.385 */
#define LOG(X) log((double)((X) + EPSILON))
#define DELOG(X) (exp((double)(X)) - EPSILON)
#define LOG2(X) (LOG(X)/0.6931471805599452862)
/*#define LOG2(X) log2((double)((X) + EPSILON))*/

/* used for summing logarithms:  log(x + y) where log(x) and log(y) are avail.*/
#define BITS 64.0
#define LOGZERO -1e100			/* log of number close to zero */
#define LOGZEROI -2147483647		/* log of number close  zero (integer)*/
#define MEME_LOG_SUM1(logx, logy) (((((logy)<=LOGZERO) || ((logx)-(logy))>BITS)) ? \
    (logx) : (logx) + log( 1 + exp((logy) - (logx) ) ) )
#define MEME_LOG_SUM(logx, logy) ( ( (logx) > (logy) ) ? \
  MEME_LOG_SUM1( (logx), (logy) ) : MEME_LOG_SUM1( (logy), (logx) ) )

/* Largest X  L = e^(-1/SCALE_LOGS); larger X yeild INT_LOG(x) = 0 */
#define SCALE_LOGS 1e6			/* don't use larger! */
#define INT_LOG(X) ((int)(SCALE_LOGS * LOG(X)))
#define INT_DELOG(X) DELOG((X)/SCALE_LOGS)

/*
   get the exponent and mantissa of large numbers expressed as log10(x)
   for printing with prec digits after the decimal point in mantissa
*/
#define exp10_logx(logx, m, e, prec) { \
  (e) = floor(logx); \
  (m) = exp10((logx) - (e)); \
  if (m+(.5*pow(10,-prec)) >= 10) { (m) = 1; (e) += 1;} \
}

/*
  Compute the log p-value of the extreme value of n independent observations
  given the single-trial probability p of the observed extreme;
  	if (n*p < 1e-6) pv = n*p
        else if (n*p > A) pv = 1.0
        else pv = 1 - (1-A)**(n*p/A)
*/
#define _A_ 1e-6
#define LOGMINPROD -13.8155105579643	/* log(1e-6) */
#define LOGMAXPROD 4.60517018598809	/* log(100) */
#define LOGEV(logn, logp)						\
  ( ((logn)+(logp) < LOGMINPROD) ? (logn)+(logp) :			\
    ( ((logn)+(logp) > LOGMAXPROD) ? 0 :				\
     log(1 - pow((1-_A_), exp((logn)+(logp)-LOGMINPROD))) ) )

/*
	Malloc stuff for improved debugging
*/

/* make a better malloc function */
#define mymalloc(x) ((x) > 0 ? malloc(x) : NULL)

#ifdef MALLOC_DEBUG			/* debug mode */

/* make a better free function */
#define myfree(P) {							\
  fprintf(stdout, "myfree(" #P "(%x)(mhdr %x %x))\n",			\
    (unsigned)(P), ((unsigned *)(P))[-2], ((unsigned *)(P))[-1]);	\
  if (P) free((char *) (P)); (P) = NULL;				\
}

#define dbmalloc(P, N, T) {						\
  (P) = (N) > 0 ? (T *) malloc((N)*sizeof(T)) : NULL;			\
  fprintf(stdout, 							\
    "dbmalloc(" #P ", " #N "(%d), " #T ")=%x (mhdr %x %x)\n",		\
    (N), 								\
    (unsigned)(P), ((unsigned *)(P))[-2], ((unsigned *)(P))[-1]);	\
}

#define myrealloc(P, OLDP, N, T) {					\
  (P) = (T *) realloc((malloc_t)(OLDP), (N)*sizeof(T));			\
  fprintf(stdout, 							\
    "myrealloc(" #P ", " #OLDP ", " #N "(%d), " #T ")=%x (mhdr %x %x)\n",\
    (N), 								\
    (unsigned)(P), ((unsigned *)(P))[-2], ((unsigned *)(P))[-1]);	\
}

/* like free but just checks */
#define malloc_chk(P)							\
  fprintf(stdout, "myfree(" #P "(%X)(mhdr %X %X %X %X))\n",		\
    (unsigned)(P), ((unsigned *)(P))[-4], ((unsigned *)(P))[-3],	\
    ((unsigned *)(P))[-2], ((unsigned *)(P))[-1]);

#else				/* non-debug mode */

/* make a better free function */
#define myfree(P) {							\
  if (P) free((char *) (P)); (P) = NULL;						\
}
#define dbmalloc(P, N, T) {						\
  (P) = (N) > 0 ? (T *) malloc((N)*sizeof(T)) : NULL;			\
}
#define myrealloc(P, OLDP, N, T) {					\
  (P) = (T *) realloc((malloc_t)(OLDP), (N)*sizeof(T));			\
}

#endif	/* MALLOC_DEBUG */

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

/* extract a substring of length n from S2 starting at position i and place
   in new string S1
*/
#define Substr(S1, S2, i, n) { 						\
  S1 = NULL; Resize(S1, n+1, char);					\
  strncpy(S1, S2+i, n); 						\
  S1[n] = '\0';								\
}

/* Useful data types and constants */
typedef int BOOLEAN;
#define FALSE 0
#define TRUE 1

/* macros to create and destroy an r x c array;
   avoids problems with multiply subscripted C array elements
   behaving differently when passed to subroutines or assigned locally
*/
#define create_2array(v, typ, r, c) {					\
  int _i, _ir=(r), _ic=(c);		/* avoid side effects */	\
  /*fprintf(stderr, */							\
  /*  "create_2array(" #v ", " #typ ", " #r "(%d), " #c ")\n", r);*/	\
  dbmalloc((v), _ir,  typ *);						\
  if (!(v)) {fprintf(stderr, "malloc failed 1\n"); Crash;}		\
  for (_i=0; _i<_ir; _i++) {						\
    typ *tmp;								\
    (v)[_i] = 0;							\
    dbmalloc(tmp, _ic, typ);						\
    (v)[_i] = tmp;							\
    if (!(v)[_i]) {printf("malloc failed 2\n");}			\
    if (!(v)[_i]) {fprintf(stderr, "malloc failed 2\n"); Crash;}	\
  }									\
}

#define free_2array(v, r) {						\
  int _i, _ir=(r);			/* avoid side effects */	\
  /*fprintf(stderr, "free_2array(" #v ", " #r "(%d))\n", r);	*/	\
  for (_i=0; _i<_ir; _i++) myfree((v)[_i]);				\
  myfree(v);								\
}
#endif


/* I/O functions */
/* skip until end of line; start of next line to be read when done */
#define Skip_eol(c, file) while ((c) != '\n' && (c) != EOF) (c) = fgetc(file)
/* skip white space; start of next token to be read when done */
#define Skip_whi(c, file) { \
  while ((c) == ' ' || (c) == '\t') (c) = fgetc(file); \
  ungetc((c), (file)); \
}

/*
        cpyfile

        Copy from input file to output file until end of file is reached
*/
#define cpyfile(IN, OUT)						\
{									\
  int c;								\
  while ((c=getc(IN)) != EOF) putc(c, (OUT));				\
}

/*
	Perl2C stuff
*/

/*
  Getline(stream, string, len)
  Read a newline or EOF-terminated line from stream.
  string is a pointer to the line or NULL if at EOF when called.
  len is the length of the returned string including newline if present.
  The string is freed if it initially is non-NULL.
*/
#define GLBUFSIZ 1000
#define Getline(stream, string, len) \
{									\
  int c; 								\
  if (string) free(string); 						\
  string = NULL; 							\
  for (len=0; (c=getc(stream)) != EOF; len++) {				\
    if (len % GLBUFSIZ == 0) Resize(string, len+GLBUFSIZ+1, char); 	\
    string[len] = c;							\
    if (c == '\n') break;						\
  }									\
  if (string) {string[len] = '\0'; Resize(string,len+1,char);}		\
}									\

/*
  Split(string, array, len)
  Split a null-terminated string on whitespace.
  Place results into a array of strings.
  len is the length of the array.
  The array is freed if it initially is non-NULL.
  Don't change the values in array or it may cause crashes!
*/
#define SPBUFSIZ 100
#define Split(string, array, len) \
{									\
  int i, i1, i2;			/* position in string */	\
  for (i=0; array && array[i]; i++) {					\
    free(array[i]); 							\
  }									\
  if (array) free(array);						\
  array=NULL; 								\
  for (len=0, array=NULL, i1=0, i2=0; ; len++, i1=i2) {			\
    while (string[i1]==' ' || string[i1]=='\t' || string[i1]=='\n') i1++;\
    if (string[i1]=='\0') break;					\
    i2 = i1;								\
    while (string[i2]!=' ' && string[i2]!='\t' && string[i2]!='\n' && 	\
      string[i2]!='\0') i2++;						\
    if (len % SPBUFSIZ == 0) Resize(array, len+SPBUFSIZ, char *); 	\
    array[len] = NULL; Resize(array[len], i2-i1+1, char);		\
    for (i=0; i1<i2; i++, i1++) array[len][i] = string[i1];		\
    array[len][i] = '\0';						\
  }									\
  if (array) {								\
    Resize(array, len+1, char *);					\
    array[len] = NULL;							\
  }									\
}									\

/*
  Die(string)
  Print a string and exit 1.
*/
#define Die(string)							\
{									\
  fprintf(stderr, "%s", string);					\
  exit(1);								\
}

/* concatenate two strings into one new one */
#define Strcat(dst, s1, s2) {						\
  int w1=strlen(s1), w2=strlen(s2);					\
  dst = NULL; Resize(dst,w1+w2+1, char);				\
  strcpy(dst, (s1)); 							\
  strcpy(dst+w1, (s2));							\
}

/* turn a macro into a string: MAKE_STRING(x) turns macro "x" into a string */
#define MAKE_STRING(x) MAKE_STRING1(x)
#define MAKE_STRING1(x) #x


/* end of Perl2C stuff */

/*
  Compute the 1-CDF of the normal distribution:
	Pr(N(0,1) < z)
*/
#define SQRT2 1.4142135623731
#define NQ(z) (0.5*erfc((z)/SQRT2))

/* print a debug message once only */
#ifdef DEBUG
#define print_once(msg) {static BOOLEAN ft; if (!ft) printf(msg); ft=1;}
#else
#define print_once(msg) {}
#endif

#ifdef PARALLEL
#include "mp.h"
// #define printf printf("node %d: ", mpMyID()); fflush(stdout); printf
#define printf if (mpMyID() == 0) printf
//#define fprintf fprintf(stderr, "node %d: ", mpMyID()); fflush(stderr); fprintf //
#ifndef DEBUG_PARALLEL
#define fprintf if (mpMyID() == 0) fprintf
#else
#endif

#define exit(x) {mpFinalize(); exit(x);}
#endif
