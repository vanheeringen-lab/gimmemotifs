#include <assert.h>
#include <ctype.h>  /* isspace */
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdlib.h>
#include "glam2_util.h"

const char *prog_name = "";

int can_add_int(const int x, const int y) {
  if (x >= 0)
    return INT_MAX - x >= y;
  else
    return INT_MIN - x <= y;
}

/* Uses "div" because "/" is not well defined for negative operands */
int can_mul_int(const int x, const int y) {
  if (x > 0) {
    if (y >= 0)
      return y <= div(INT_MAX, x).quot;
    else
      return y >= div(INT_MIN, x).quot;
  } else if (x < 0) {
    if (y >= 0)
      return y <= div(INT_MIN, x).quot;
    else
      return y >= div(INT_MAX, x).quot;
  } else  /* x == 0 */
    return 1;
}

void memswap(void *s, void *t, size_t n) {
  char *a = s;
  char *b = t;

  while (n--) {  /* if n > 0, do-while might be faster? */
    const char tmp = *a;
    *a++ = *b;
    *b++ = tmp;
  }
}

/* Return the greatest common divisor of a and b: Euclidean algorithm */
static size_t gcd(size_t a, size_t b) {
  while (b) {
    size_t c = a % b;
    a = b;
    b = c;
  }
  return a;
}

void memrol(void *s, size_t n, size_t k) {
  size_t cycles;
  char *start;
  if (n == 0)  /* prevent % 0 */
    return;
  k %= n;  /* allow rotations of size > n */
  if (k == 0)
    return;
  assert(s != NULL);
  cycles = gcd(n, k);

  for (start = s; cycles != 0; ++start) {
    char *x = start + k;
    char tmp = *x;
    while (x != start) {
      char *next = x < (char *)s+n-k ? x+k : x-(n-k);
      *x = *next;
      x = next;
    }
    *start = tmp;
    --cycles;
  }
}

void memror(void *s, size_t n, size_t k) {
  size_t cycles;
  char *start;
  if (n == 0)  /* prevent % 0 */
    return;
  k %= n;  /* allow rotations of size > n */
  if (k == 0)
    return;
  assert(s != NULL);
  cycles = gcd(n, k);

  for (start = s; cycles != 0; ++start) {
    char *x = start + (n-k);
    char tmp = *x;
    while (x != start) {
      char *next = x < (char *)s+k ? x+(n-k) : x-k;
      *x = *next;
      x = next;
    }
    *start = tmp;
    --cycles;
  }
}

void set_int(int *ptr, const size_t size, const int value) {
  const int *end = ptr + size;
  assert(ptr != NULL || size == 0);
  for (; ptr != end; ++ptr)
    *ptr = value;
}

void set_dbl(double *ptr, const size_t size, const double value) {
  const double *end = ptr + size;
  assert(ptr != NULL || size == 0);
  for (; ptr != end; ++ptr)
    *ptr = value;
}

int *min_int(const int *ptr, const size_t size) {
  const int *end = ptr + size;
  const int *min = ptr;
  assert(ptr != NULL || size == 0);
  for (; ptr != end; ++ptr)
    if (*ptr < *min)
      min = ptr;
  return (int *)min;
}

int *max_int(const int *ptr, const size_t size) {
  const int *end = ptr + size;
  const int *max = ptr;
  assert(ptr != NULL || size == 0);
  for (; ptr != end; ++ptr)
    if (*ptr > *max)
      max = ptr;
  return (int *)max;
}

double *min_dbl(const double *ptr, const size_t size) {
  const double *end = ptr + size;
  const double *min = ptr;
  assert(ptr != NULL || size == 0);
  for (; ptr != end; ++ptr)
    if (*ptr < *min)
      min = ptr;
  return (double *)min;
}

double *max_dbl(const double *ptr, const size_t size) {
  const double *end = ptr + size;
  const double *max = ptr;
  assert(ptr != NULL || size == 0);
  for (; ptr != end; ++ptr)
    if (*ptr > *max)
      max = ptr;
  return (double *)max;
}

int sum_int(const int *ptr, const size_t size) {
  const int *end = ptr + size;
  int sum = 0;
  assert(ptr != NULL || size == 0);
  for (; ptr != end; ++ptr) {
    assert(can_add_int(sum, *ptr));
    sum += *ptr;
  }
  return sum;
}

double sum_dbl(const double *ptr, const size_t size) {
  const double *end = ptr + size;
  double sum = 0;
  assert(ptr != NULL || size == 0);
  for (; ptr != end; ++ptr)
    sum += *ptr;
  return sum;
}

void add_int(int *ptr, const size_t size, const int x) {
  const int *end = ptr + size;
  assert(ptr != NULL || size == 0);
  for (; ptr != end; ++ptr) {
    assert(can_add_int(x, *ptr));
    *ptr += x;
  }
}

void add_dbl(double *ptr, const size_t size, const double x) {
  const double *end = ptr + size;
  assert(ptr != NULL || size == 0);
  for (; ptr != end; ++ptr)
    *ptr += x;
}

void mul_dbl(double *ptr, const size_t size, const double x) {
  const double *end = ptr + size;
  assert(ptr != NULL || size == 0);
  for (; ptr != end; ++ptr)
    *ptr *= x;
}

void div_dbl(double *ptr, const size_t size, const double x) {
  const double *end = ptr + size;
  assert(ptr != NULL || size == 0);
  assert(x != 0 || size == 0);
  for (; ptr != end; ++ptr)
    *ptr /= x;
}


void pow_dbl(double *ptr, const size_t size, const double x) {
  const double *end = ptr + size;
  assert(ptr != NULL || size == 0);
  for (; ptr != end; ++ptr)
    *ptr = xpow(*ptr, x);
}

// conflicts with array.c function of same name
//void normalize(double *ptr, const size_t size) {
//  const double tot = sum_dbl(ptr, size);
//  div_dbl(ptr, size, tot);
//}

double log_sum(const double *ptr, const size_t size) {
  if (size == 1)  /* do special case faster */
    return *ptr;
  else {
    const double max = *max_dbl(ptr, size);
    const double *end = ptr + size;
    double sum = 0;
    assert(size > 0);
    for (; ptr != end; ++ptr)
      sum += xexp(*ptr - max);
    return xlog(sum) + max;
  }
}

int rand_int(const unsigned n) {  /* from C FAQ */
  unsigned d, threshold, r;
  assert(n > 0);
  assert(n <= RAND_MAX + 1u);
  d = (RAND_MAX + 1u) / n;
  threshold = d * n;
  do {
    r = rand();
  } while(r >= threshold);
  return r / d;
}

double rand_dbl(const double n) {
  double r;
  assert(n == n);  /* check for NaN */
  assert(n > 0);
  assert(n <= DBL_MAX);  /* check for inf */
  do {
    r = rand() / (RAND_MAX + 1.0) * n;  /* from C FAQ */
  } while (r >= n);  /* r can equal n when n is tiny (e.g. < DBL_MIN) */
  return r;
}

double *pick_dbl(const double *ptr, const size_t size) {
  const double tot = sum_dbl(ptr, size);
  const double r = rand_dbl(tot);
  const double *end = ptr + size;
  double s = 0;
  for (; ptr != end; ++ptr)
    if ((s += *ptr) > r)
      break;
  assert(ptr != end);
  return (double *)ptr;  /* strange that this cast is needed */
}

void shuffle(void *base, size_t n, size_t size) {
  for (; n > 1; --n) {
    const int r = rand_int(n);
    memswap((char *)base + (n-1) * size, (char *)base + r * size, size);
  }
}

int digits(int n) {
  int d = 1;
  assert(n >= 0);
  while (n /= 10)
    ++d;
  return d;
}

#ifdef VOID
void die(const char *fmt, ...) {
  va_list args;
  va_start(args, fmt);
  vfprintf(stderr, fmt, args);
  va_end(args);
  exit(EXIT_FAILURE);
}
#endif

void *xmalloc(size_t size) {
  void *p = malloc(size);
  if (p == NULL && size != 0)
    die("%s: error allocating %lu bytes: %s\n",
	prog_name, (unsigned long)size, strerror(errno));
  return p;
}

void *xcalloc(size_t nmemb, size_t size) {
  void *p = calloc(nmemb, size);
  if (p == NULL && nmemb != 0 && size != 0)
    die("%s: error calculated-allocating %lu * %lu bytes: %s\n",
        prog_name, (unsigned long)nmemb, (unsigned long)size, strerror(errno));
  return p;
}

void *xrealloc(void *p, size_t size) {
  void *q = realloc(p, size);
  if (q == NULL && size != 0)
    die("%s: error reallocating %lu bytes: %s\n",
	prog_name, (unsigned long)size, strerror(errno));
  return q;
}

char *xstrdup(const char *cs) {
  char *t;
  assert(cs != NULL);
  t = xmalloc(strlen(cs) + 1);
  return strcpy(t, cs);
}

void *xmemdup(const void *cs, size_t n) {
  char *t;
  assert(cs != NULL);
  t = xmalloc(n);
  return memcpy(t, cs, n);
}

char *strtrunc(char *s, size_t n) {
  size_t i;
  assert(s != NULL);
  for (i = 0; i != n; ++i)
    if (s[i] == 0)
      return s;
  s[n] = 0;
  return s;
}

void *xmalloc2(size_t rows, size_t cols) {
  char **ptrs;
  size_t i;

  if (rows > (size_t)-1 / sizeof(char *))
    die("%s: memory requirement too large: %lu * %lu bytes\n",
	prog_name, (unsigned long)rows, (unsigned long)sizeof(char *));

  XMALLOC(ptrs, rows);
  for (i = 0; i < rows; ++i)
    ptrs[i] = xmalloc(cols);
  return ptrs;
}

void *xcalloc2(size_t rows, size_t cols, size_t size) {
  char **ptrs;
  size_t i;

  if (rows > (size_t)-1 / sizeof(char *))
    die("%s: memory requirement too large: %lu * %lu bytes\n",
	prog_name, (unsigned long)rows, (unsigned long)sizeof(char *));

  XMALLOC(ptrs, rows);
  for (i = 0; i < rows; ++i)
    ptrs[i] = xcalloc(cols, size);
  return ptrs;
}

void free2(void *p, size_t rows) {
  char **ptrs = p;
  size_t i;
  if (p == NULL)  /* same behaviour as free (do we want this?) */
    return;

  for (i = 0; i < rows; ++i)
    free(ptrs[i]);
  free(ptrs);
}

FILE *xfopen(const char *filename, const char *mode) {
  FILE *fp;
  assert(filename != NULL);
  assert(mode != NULL);

  if (strcmp(filename, "-") == 0) {
    if (strcmp(mode, "r") == 0)
      return stdin;
    else if (strcmp(mode, "w") == 0 || strcmp(mode, "a") == 0)
      return stdout;
    else
      assert(0);  /* shouldn't get here (?) */
  }

  fp = fopen(filename, mode);
  if (fp == NULL)
    die("%s: error opening file %s: %s\n",
	prog_name, filename, strerror(errno));
  return fp;
}

int xfclose(FILE *stream) {
  int x = fclose(stream);
  if (x == EOF)
    die("%s: error closing file: %s\n", prog_name, strerror(errno));
  return x;
}

int xungetc(int c, FILE *stream) {
  int x = ungetc(c, stream);
  if (x == EOF)
    die("%s: error reading file: %s\n", prog_name, strerror(errno));
  return x;
}

size_t xgetline(char **lineptr, size_t *n, FILE *stream) {
  size_t i = 0;  /* number of characters read so far */
  int c;  /* this character */
  assert(lineptr != NULL);

  if (*lineptr == NULL || *n == 0) {
    *n = 1;  /* should probably use a bigger starting size */
    *lineptr = XREALLOC(*lineptr, *n);
  }

  while((c = getc(stream)) != EOF) {
    ++i;
    if (i == *n) {
      *n = *n * 2;
      if (*n <= i)  /* overflow */
        die("%s: line too long: >= %lu\n", prog_name, (unsigned long)i);
      *lineptr = XREALLOC(*lineptr, *n);
    }
    (*lineptr)[i-1] = c;  /* int -> char! */
    if (c == '\n')
      break;
  }

  if (ferror(stream))
    die("%s: error reading file: %s\n", prog_name, strerror(errno));
  assert(i < *n);
  (*lineptr)[i] = '\0';
  return i;
}

size_t chomp(char *s, size_t n) {
  assert(s != NULL);
  if (n != 0 && s[n-1] == '\n') {
    --n;
    s[n] = '\0';
  }
  return n;
}

char *skipws(const char *cs) {
  assert(cs != NULL);
  while (isspace((unsigned char)*cs))
    ++cs;
  return (char *)cs;
}

char *skipnw(const char *cs) {
  assert(cs != NULL);
  while (!isspace((unsigned char)*cs) && *cs != '\0')
    ++cs;
  return (char *)cs;
}

#if 0  /* is this function useful? */
char *next_word(char **stringp) {
  unsigned char *beg;  /* ??? isspace needs unsigned */
  unsigned char *end;  /* ??? isspace needs unsigned */

  assert(stringp != NULL);
  if (*stringp == NULL)
    return NULL;
  for (beg = *stringp; isspace(*beg); ++beg)
    ;
  for (end = beg, !isspace(*end) && *end != '\0'; ++end)
    ;
  if (*beg == '\0')
    beg = NULL;
  if (end == '\0')
    end = NULL;
  else {
    *end = '\0';
    ++end;  /* end might now point to '\0' */
  }
  *stringp = end;
  return beg;
}
#endif

double xexp(double x) {
  double z;
  errno = 0;
  z = exp(x);
  if (errno != 0 && !(errno == ERANGE && z == 0))  /* allow underflow */
    die("%s: error exponentiating %g: %s\n",
	prog_name, x, strerror(errno));
  return z;
}

double xlog(double x) {
  double z;
  errno = 0;
  z = log(x);
  assert(!(errno == ERANGE && z == 0));  /* underflow should be impossible */
  if (errno != 0)
    die("%s: error taking logarithm of %g: %s\n",
	prog_name, x, strerror(errno));
  return z;
}

double xpow(double x, double y) {
  double z;
  errno = 0;
  z = pow(x, y);
  if (errno != 0 && !(errno == ERANGE && z == 0))  /* allow underflow */
    die("%s: error raising %g to the power of %g: %s\n",
	prog_name, x, y, strerror(errno));
  return z;
}

double xatof(const char *s) {
  double out;
  char *endp;
  errno = 0;
  out = strtod(s, &endp);
  if (endp == s || *endp != '\0' || errno == ERANGE)  /* catches underflow */
    die("%s: error converting %s to double\n", prog_name, s);
  return out;
}

int xatoi(const char *s) {
  long out;
  char *endp;
  errno = 0;
  out = strtol(s, &endp, 10);
  if (endp == s || *endp != '\0' || errno == ERANGE ||
      out > INT_MAX || out < INT_MIN)
    die("%s: error converting %s to integer\n", prog_name, s);
  return out;
}

unsigned xatou(const char *s) {
  unsigned long out;
  char *endp;
  errno = 0;
  out = strtoul(s, &endp, 10);
  if (endp == s || *endp != '\0' || errno == ERANGE || out > UINT_MAX)
    die("%s: error converting %s to unsigned integer\n", prog_name, s);
  return (unsigned)out;
}

int put_pad(int c, int n, FILE *stream) {
  assert(n >= 0);
  while (n--)
    if (putc(c, stream) == EOF)
      return EOF;
  return (unsigned char)c;  /* like fputc(?) */
}
