#include <stdlib.h>  /* free */
#include <string.h>  /* strlen */
#include "glam2_util.h"
#include "glam2_alignment.h"

int aln_same_lengths(const alignment *a) {
  const size_t len = strlen(a->key_positions);
  size_t i;
  for (i = 0; i != a->seq_num; ++i)
    if (strlen(a->seqs[i].seq) != len)
      return 0;
  return 1;
}

static char *next_word(const char *cs) {
  return skipws(skipnw(cs));
}

void aln_read(alignment *a, FILE *stream) {
  char *line = NULL;
  size_t line_size = 0;
  int state = 0;
  a->seq_num = 0;
  a->seqs = NULL;
  a->key_positions = NULL;

  while (xgetline(&line, &line_size, stream)) {
    if (state == 0) {
      char *beg = skipws(line);
      char *end = skipnw(beg);
      if (*beg == '*' && *(end-1) == '*') {
        *end = 0;
        a->key_positions = xstrdup(beg);
        state = 1;
      }
    } else {
      char *name = line;
      char *start = next_word(name);
      char *seq = next_word(start);
      if (*seq == 0)
        break;
      *skipnw(name) = 0;
      *skipnw(start) = 0;
      *skipnw(seq) = 0;
      ++a->seq_num;
      /* slow but simple linear reallocation for now: */
      a->seqs = XREALLOC(a->seqs, a->seq_num);
      a->seqs[a->seq_num-1].name = xstrdup(name);
      a->seqs[a->seq_num-1].seq = xstrdup(seq);
      a->seqs[a->seq_num-1].start = xstrdup(start);
    }
  }

  free(line);
}

static void seq_free(aligned_seq *s) {
  free(s->name);
  free(s->seq);
  free(s->start);
}

void aln_free(alignment *a) {
  size_t i;
  for (i = 0; i != a->seq_num; ++i)
    seq_free(&a->seqs[i]);
  free(a->seqs);
  free(a->key_positions);
}
