#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include "glam2_util.h"
#include "glam2_alphabet.h"

/* almost general-purpose enough to go into util? */
static void counts2probs(double *probs, const int *counts, const size_t size) {
  const double tot = sum_int(counts, size);
  size_t i;
  for (i = 0; i < size; ++i)
    probs[i] = counts[i] / tot;
}

void alphabet_n(alphabet *a) {
  const char *string = "acgt";
  const int counts[] = { 1, 1, 1, 1 };
  int i;
  a->size = strlen(string);
  set_int(a->encode, UCHAR_MAX+1, a->size);

  for (i = 0; i < a->size; ++i) {
    const char x = string[i];
    a->encode[tolower(x)] = i;
    a->encode[toupper(x)] = i;
    a->decode[i] = tolower(x);
  }
  a->encode['u'] = a->encode['t'];  /* for RNA */
  a->encode['U'] = a->encode['T'];  /* for RNA */
  a->decode[a->size] = 'N';  /* wildcard symbol */

  XMALLOC(a->prob, a->size);
  counts2probs(a->prob, counts, a->size);
}

void alphabet_p(alphabet *a) {
  const char *string = "acdefghiklmnpqrstvwy";
  const int counts[] = { 35155, 8669, 24161, 28354, 17367, 33229, 9906, 23161, 25872, 40625, 10101, 20212, 23435, 19208, 23105, 32070, 26311, 29012, 5990, 14488 };
  int i;
  a->size = strlen(string);
  set_int(a->encode, UCHAR_MAX+1, a->size);

  for (i = 0; i < a->size; ++i) {
    const char x = string[i];
    a->encode[tolower(x)] = i;
    a->encode[toupper(x)] = i;
    a->decode[i] = toupper(x);
  }
  a->decode[a->size] = 'x';  /* wildcard symbol */

  XMALLOC(a->prob, a->size);
  counts2probs(a->prob, counts, a->size);
}

void alphabet_read(alphabet *a, FILE *fp) {
  int c;
  int state = 0;
  double tot;

  a->size = 0;
  set_int(a->encode, UCHAR_MAX+1, -1);
  a->prob = NULL;

  while((c = getc(fp)) != EOF) {
    if (c == '#') {  /* start of comment */
      fscanf(fp, "%*[^\n]");  /* skip line */
    } else if (c == '\n') {
      state = 0;
    } else if (isspace(c)) {
      if (state == 1) {
	double p = 0;
        const int r = fscanf(fp, "%lf", &p);
        if (r == 1 && p <= 0)
          die("%s: error reading alphabet file: abundance %g should be > 0\n",
              prog_name, p);
        a->prob[a->size - 1] = p;
        fscanf(fp, "%*[^\n]");  /* skip line */
      }
    } else {
      if (state == 0) {
	++a->size;
	a->decode[a->size - 1] = c;
        a->prob = XREALLOC(a->prob, a->size);
        a->prob[a->size - 1] = 1;  /* default abundance */
        state = 1;
      }
      if (a->encode[c] != -1 && a->encode[c] != a->size - 1)
	die("%s: error reading alphabet file: found %c twice\n",
	    prog_name, c);
      a->encode[c] = a->size - 1;
    }
  }
  if (ferror(fp))
    die("%s: error reading alphabet file: %s\n", prog_name, strerror(errno));

  if (a->size <= 1)  /* disallow size zero alphabet */
    die("%s: error reading alphabet file: no alphabet\n", prog_name);

  --a->size;  /* last line has the wildcard symbol */

  for (c = 0; c <= UCHAR_MAX; ++c)
    if (a->encode[c] == -1)
      a->encode[c] = a->size;

  tot = sum_dbl(a->prob, a->size);
  div_dbl(a->prob, a->size, tot);
}

void alphabet_free(alphabet *a) {
  assert(a != NULL);
  free(a->prob);
}
