/* Functions for handling FASTA format sequences */
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <stdlib.h>  /* free */
#include "glam2_util.h"
#include "glam2_fasta.h"

char *first_word(char *s) {
  const char *beg = skipws(s);
  const char *end = skipnw(beg);
  MOVE(s, beg, end - beg);  /* overflow danger? */
  s[end - beg] = '\0';  /* overflow danger? */
  return s;
}

void count_fasta(const fasta *f, int *counts) {
  int i;
  for (i = 0; i < f->seqlen; ++i)
    ++counts[f->seq[i]];
}

void count_mfasta(const mfasta *m, int *counts) {
  int i;
  for (i = 0; i < m->seqnum; ++i)
    count_fasta(&m->f[i], counts);
}

void rc_fasta(fasta *f, const int alph_size) {
  const int seqlen = f->seqlen;  /* might be zero */
  const int *seq = f->seq;  /* might be NULL */
  int *rcseq = xmalloc(seqlen * sizeof(int));  /* might be NULL */
  int i;
  for (i = 0; i < seqlen; ++i)
    rcseq[seqlen - i - 1] = seq[i] < alph_size ?
      alph_size - seq[i] - 1 : alph_size;
  f->rcseq = rcseq;
}

void rc_mfasta(mfasta *m, const int alph_size) {
  int i;
  for (i = 0; i < m->seqnum; ++i)
    rc_fasta(&m->f[i], alph_size);
}

void tr_fasta(fasta *f, const int *encode) {
  int i;
  for (i = 0; i < f->seqlen; ++i)
    f->seq[i] = encode[f->seq[i]];
}

void tr_mfasta(mfasta *m, const int *encode) {
  int i;
  for (i = 0; i < m->seqnum; ++i)
    tr_fasta(&m->f[i], encode);
}

static int get_maxlen(const mfasta *m) {
  int i;
  int max = 0;
  for (i = 0; i < m->seqnum; ++i)
    if (m->f[i].seqlen > max)
      max = m->f[i].seqlen;
  return max;
}

void put_fasta(const fasta *f, int line_size, FILE *fp) {
  int i;
  assert(line_size > 0);
  fprintf(fp, ">%s\n", f->title);

  for (i = 0; i < f->seqlen; ++i) {
    putc(f->seq[i], fp);
    if ((i+1) % line_size == 0 || i+1 == f->seqlen)
      putc('\n', fp);
  }
}

/* double memory allocation or die */
static void *moremem(void *buf, int *mem) {
  if (!can_mul_int(*mem, 2))
    die("%s: memory requirement too large: %d * 2 bytes\n", prog_name, *mem);
  *mem *= 2;
  return xrealloc(buf, *mem);
}

/* read FASTA title line */
static size_t get_title(char **title, FILE *fp) {
  char *buf = NULL;
  size_t mem = 0;
  size_t len = xgetline(&buf, &mem, fp);
  len = chomp(buf, len);
  buf = XREALLOC(buf, len+1);  /* shrink-wrap */
  *title = buf;
  return len;  /* might be zero */
}

/* read FASTA sequence */
static int get_seq(int **seq, FILE *fp) {
  int len = 0;  /* sequence length so far */
  int mem = sizeof(int);  /* bytes allocated to sequence so far */
  int *buf = xmalloc(mem);  /* sequence */
  register int c;

  while((c = getc(fp)) != EOF) {
    if (c == '>')
      break;
    else if (isspace(c))  /* we might want to read e.g. gap characters */
      continue;
    if (len * sizeof(int) == mem)
      buf = moremem(buf, &mem);
    buf[len] = c;
    ++len;
  }
  if (ferror(fp))
    die("%s: error reading fasta file: %s\n", prog_name, strerror(errno));

  if (c == '>')
    xungetc(c, fp);
  buf = XREALLOC(buf, len);  /* shrink-wrap */
  *seq = buf;  /* might be NULL */
  return len;  /* might be zero */
}

int fasta_read(fasta *f, FILE *fp) {
  char *title;
  int *seq;
  size_t titlen;
  int seqlen;
  register int c;

  /* skip any junk before '>': */
  while ((c = getc(fp)) != EOF)
    if (c == '>')
      break;
  if (ferror(fp))
    die("%s: error reading fasta file: %s\n", prog_name, strerror(errno));
  if (c == EOF)
    return EOF;

  titlen = get_title(&title, fp);
  seqlen = get_seq(&seq, fp);

  f->title = title;
  f->titlen = titlen;
  f->seq = seq;
  f->seqlen = seqlen;
  f->rcseq = NULL;  /* so it can be freed */
  return 0;
}

void mfasta_read(mfasta *m, FILE *fp) {
  int num = 0;  /* number of fastas read so far */
  int mem = sizeof(fasta);  /* bytes allocated for fastas so far */
  fasta *f = xmalloc(mem);  /* fastas */

  while (fasta_read(&f[num], fp) != EOF) {
    ++num;
    if (num * sizeof(fasta) == mem)
      f = moremem(f, &mem);
  }

  f = XREALLOC(f, num);  /* shrink-wrap */
  m->f = f;  /* might be NULL */
  m->seqnum = num;  /* might be zero */
  m->maxlen = get_maxlen(m);
}

void free_fasta(fasta *f) {
  free(f->title);
  free(f->seq);
  free(f->rcseq);
}

void free_mfasta(mfasta *m) {
  int i;
  for (i = 0; i < m->seqnum; ++i)
    free_fasta(&m->f[i]);
  free(m->f);
}
