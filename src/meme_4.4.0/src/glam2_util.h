/* Useful, general-purpose functions */
#ifndef GLAM2_UTIL_H
#define GLAM2_UTIL_H

#include <stdio.h>  /* FILE */
#include <stdlib.h>  /* qsort */
#include <string.h>  /* memcpy, memmove, memset */
#include "utils.h"              // needed for die() (meme version)

#define COPY(s, t, n) memcpy(s, t, (n) * sizeof *(s))
#define MOVE(s, t, n) memmove(s, t, (n) * sizeof *(s))
#define ZERO(s, n) memset(s, 0, (n) * sizeof *(s))
#define XDUP(s, n) xmemdup(s, (n) * sizeof *(s))
#define SORT(base, n, cmp) qsort(base, n, sizeof *(base), cmp)
#define ROL(s, n, k) memrol(s, (n) * sizeof *(s), (k) * sizeof *(s))
#define ROR(s, n, k) memror(s, (n) * sizeof *(s), (k) * sizeof *(s))
#define SHUFFLE(base, n) shuffle(base, n, sizeof *(base))
#define XREALLOC(p, n) xrealloc(p, (n) * sizeof *(p))
#define XMALLOC(p, n) ((p) = xmalloc((n) * sizeof *(p)))
#define COUNTOF(a) (sizeof (a) / sizeof *(a))

extern const char *prog_name;  /* set this to get nicer error messages */

/* Put some of these functions in the header so they may get inlined? */
/* Should some of these return their 1st argument rather than "void"? */

/* Return 1 if x and y can be added without overflow, 0 otherwise */
int can_add_int(const int x, const int y);

/* Return 1 if x and y can be multiplied without overflow, 0 otherwise */
int can_mul_int(const int x, const int y);

/* Swap memory blocks s and t of size n */
void memswap(void *s, void *t, size_t n);

/* Rotate memory k bytes to the left. Works for all values of k */
void memrol(void *s, size_t n, size_t k);

/* Rotate memory k bytes to the right. Works for all values of k */
void memror(void *s, size_t n, size_t k);

/* Set an array of ints all to "value" */
void set_int(int *ptr, const size_t size, const int value);

/* Set an array of doubles all to "value" */
void set_dbl(double *ptr, const size_t size, const double value);

/* Get the minimum element of a set of ints */
/* In case of ties, it will return the first one */
/* If size=0, it will return ptr */
int *min_int(const int *ptr, const size_t size);

/* Get the maximum element of a set of ints */
/* In case of ties, it will return the first one */
/* If size=0, it will return ptr */
int *max_int(const int *ptr, const size_t size);

/* Get the minimum element of a set of doubles */
/* In case of ties, it will return the first one */
/* If size=0, it will return ptr */
double *min_dbl(const double *ptr, const size_t size);

/* Get the maximum element of a set of doubles */
/* In case of ties, it will return the first one */
/* If size=0, it will return ptr */
double *max_dbl(const double *ptr, const size_t size);

/* Get the sum of a set of ints */
int sum_int(const int *ptr, const size_t size);

/* Get the sum of a set of doubles */
double sum_dbl(const double *ptr, const size_t size);

/* Add x to each member of a set of ints */
void add_int(int *ptr, const size_t size, const int x);

/* Add x to each member of a set of doubles */
void add_dbl(double *ptr, const size_t size, const double x);

/* Multiply each member of a set of doubles by x */
void mul_dbl(double *ptr, const size_t size, const double x);

/* Divide each member of a set of doubles by x */
void div_dbl(double *ptr, const size_t size, const double x);

/* Raise each member of a set of doubles to the power of x */
void pow_dbl(double *ptr, const size_t size, const double x);

/* Divide each member of a set of doubles by their sum, so they sum to 1 */
// Conflicts with function in array.c
//void normalize(double *ptr, const size_t size);

/* Get the sum of a set of doubles, working in log space. size must be >0 */
double log_sum(const double *ptr, const size_t size);

/* Random integer between 0 (inclusive) and n (exclusive) */
/* n must be > 0 and <= RAND_MAX+1 */
int rand_int(const unsigned n);

/* Random double between 0 (inclusive) and n (exclusive) */
/* n must be > 0 */
double rand_dbl(const double n);

/* Randomly pick a member of a set of doubles, weighted by their values */
double *pick_dbl(const double *ptr, const size_t size);

/* Randomly shuffle an array */
void shuffle(void *base, size_t n, size_t size);

/* Get the number of digits in an integer (must be >= 0) */
int digits(int n);

/* Print message to stderr and exit */
//void die(const char *fmt, ...);

/* malloc or die */
void *xmalloc(size_t size);

/* calloc or die */
void *xcalloc(size_t nmemb, size_t size);

/* realloc or die */
void *xrealloc(void *p, size_t size);

/* strdup or die */
char *xstrdup(const char *cs);

/* memdup or die */
void *xmemdup(const void *cs, size_t n);

/* If string s is longer than n, truncate it to length n */
char *strtrunc(char *s, size_t n);

/* Allocate memory for a two-dimensional array or die */
void *xmalloc2(size_t rows, size_t cols);

/* Allocate memory for a two-dimensional array with zero fill or die */
void *xcalloc2(size_t rows, size_t cols, size_t size);

/* Free memory allocated by xmalloc2 */
/* Does nothing if p is NULL */
void free2(void *p, size_t rows);

/* fopen or die */
/* if filename is "-", stdin or stdout is returned, and mode must be r|w|a */
FILE *xfopen(const char *filename, const char *mode);

/* fclose or die */
int xfclose(FILE *stream);

/* ungetc or die */
int xungetc(int c, FILE *stream);

/* getline or die: returns the number of characters read */
size_t xgetline(char **lineptr, size_t *n, FILE *stream);

/* Remove final newline, if any, from string s of length n */
/* Returns new length */
size_t chomp(char *s, size_t n);

/* Skip whitespace: return pointer to 1st non-white character */
char *skipws(const char *cs);

/* Skip non-white: return pointer to 1st whitespace character or '\0' */
/* (An alternative would be to return NULL if there's no whitespace) */
char *skipnw(const char *cs);

/* exp or die */
double xexp(double x);

/* log or die */
double xlog(double x);

/* pow or die */
double xpow(double x, double y);

/* The following functions convert strings to numbers, dying on failure */
/* Failure = non-numeric, extra chars after number, overflow, underflow */
/* Leading whitespace is allowed */
double xatof(const char *s);
int xatoi(const char *s);
/* Bizarrely, xatou allows -ve numbers, because so does strtoul */
unsigned xatou(const char *s);

/* Writes c to stream n times (n must be >= 0): returns c or EOF on error */
int put_pad(int c, int n, FILE *stream);

#endif
