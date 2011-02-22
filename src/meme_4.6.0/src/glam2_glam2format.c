/* Convert glam2 output to standard alignment formats */
#include <assert.h>
#include <ctype.h>  /* tolower */
#include <stdio.h>
#include <string.h>  /* strlen, strcmp */
#include <unistd.h>  /* non-ANSI */
#include "glam2_util.h"
#include "glam2_fasta.h"
#include "glam2_alignment.h"

/* Move to util? */
static char *lowercase(char *s) {  /* better to use strcasecmp? */
  unsigned char *t;
  assert(s != NULL);
  for (t = (unsigned char *)s; *t; ++t)
    *t = tolower(*t);
  return s;
}

/* Return count of characters in cs other than c */
static size_t strccnt(const char *cs, int c) {
  size_t count = 0;
  assert(cs != NULL);
  for (; *cs != '\0'; ++cs)
    count += *cs != (char)c;  /* same cast as strchr(?) */
  return count;
}

void seq_fasta_write(const aligned_seq *s, size_t line_size, FILE *stream) {
  const char *seq = s->seq;
  size_t i;
  assert(line_size > 0);
  fprintf(stream, ">%s\n", s->name);

  for (i = 0; seq[i] != 0; ++i) {
    putc(seq[i], stream);
    if ((i+1) % line_size == 0 || seq[i+1] == 0)
      putc('\n', stream);
  }
}

void aln_fasta_write(const alignment *a, FILE *stream) {
  static const size_t line_size = 60;
  size_t i;

  for (i = 0; i != a->seq_num; ++i)
    seq_fasta_write(&a->seqs[i], line_size, stream);
}

static size_t aln_max_len(const alignment *a) {
  size_t max = 0;
  size_t i;
  for (i = 0; i != a->seq_num; ++i) {
    const size_t len = strlen(a->seqs[i].seq);
    if (len > max)
      max = len;
  }
  return max;
}

void aln_msf_write(const alignment *a, FILE *stream) {
  static const size_t line_size = 50;
  static const size_t block_size = 10;
  const size_t max_len = aln_max_len(a);
  size_t i, j, k;

  fputs("PileUp\n\n", stream);  /* is this needed? */
  fprintf(stream, "   MSF: %lu Type: X Check: 0 ..\n\n",
	  (unsigned long)max_len);

  for (i = 0; i != a->seq_num; ++i) {
    const aligned_seq *s = &a->seqs[i];
    fprintf(stream, " Name: %-10.10s Len: %lu Check: 0 Weight: 1.00\n",
	    s->name, (unsigned long)strlen(s->seq));
  }

  fputs("\n//\n\n", stream);

  for (i = 0; i < max_len; i += line_size) {  /* theoretical overflow risk */
    for (j = 0; j != a->seq_num; ++j) {
      const aligned_seq *s = &a->seqs[j];
      const size_t len = strlen(s->seq);  /* slow */

      fprintf(stream, "%-10.10s", s->name);

      for (k = 0; k != line_size && i + k < len; ++k) {
	const int c = s->seq[i+k];
	if (k % block_size == 0)
	  putc(' ', stream);
	putc(c, stream);
      }

      putc('\n', stream);
    }
    putc('\n', stream);
  }
}

/* Calculate how many columns an alignment will have after extension */
static size_t aln_extended_size(const alignment *a) {
  size_t new_size = 0;
  size_t i, j;

  for (i = 0; a->key_positions[i] != 0; ++i)
    if (a->key_positions[i] == '*')
      ++new_size;
    else
      for (j = 0; j < a->seq_num; ++j)
	new_size += a->seqs[j].seq[i] != '.';

  return new_size;
}

/* Put extension of b in a, assuming sufficient memory has been allocated */
static void aln_do_extend(alignment *a, const alignment *b) {
  const char *kp = b->key_positions;
  const size_t seqnum = b->seq_num;
  size_t apos = 0;
  size_t s, t, bpos, bend, boff;

  for (bpos = 0; ; bpos = bend) {
    for (s = 0; s != seqnum; ++s)  /* fill in the key position */
      a->seqs[s].seq[apos] = b->seqs[s].seq[bpos];
    ++apos;
    ++bpos;

    if (kp[bpos] == 0)
      break;

    assert(strchr(kp + bpos, '*') != NULL);
    bend = strchr(kp + bpos, '*') - kp;  /* next key position */

    for (s = 0; s != seqnum; ++s) {  /* fill in the insertions */
      const char *bseq = b->seqs[s].seq;

      for (boff = bpos; boff != bend; ++boff)
	if (bseq[boff] != '.') {
	  for (t = 0; t != seqnum; ++t)  /* put gaps in all other sequences */
	    a->seqs[t].seq[apos] = t == s ? bseq[boff] : '.';
	  ++apos;
	}
    }
  }

  for (s = 0; s != seqnum; ++s)  /* terminate the strings */
    a->seqs[s].seq[apos] = 0;
}

/* Extend an alignment, by un-aligning inserted residues */
static void aln_extend(alignment *a) {
  const size_t new_size = aln_extended_size(a);
  alignment tmp;
  size_t s;

  XMALLOC(tmp.seqs, a->seq_num);

  for (s = 0; s != a->seq_num; ++s)
    tmp.seqs[s].seq = xmalloc(new_size + 1);  /* space for final NUL */

  aln_do_extend(&tmp, a);

  for (s = 0; s != a->seq_num; ++s) {
    free(a->seqs[s].seq);
    a->seqs[s].seq = tmp.seqs[s].seq;
  }

  free(tmp.seqs);
}

/* Copy an array of ints into an array of chars */
static void int2char(char *s, const int *v, const size_t size) {
  size_t i;
  for (i = 0; i != size; ++i)
    s[i] = v[i];
}

/* Get either sum or maximum of left flanking sequence lengths */
static int aln_tot_left(const alignment *a, int want_sum) {
  int tot = 0;
  size_t i;

  for (i = 0; i != a->seq_num; ++i) {
    const aligned_seq *s = &a->seqs[i];
    const int left = xatoi(s->start) - 1;  /* zero-based coordinate */
    if (want_sum)
      tot += left;
    else  /* want maximum */
      if (left > tot)
	tot = left;
  }

  return tot;
}

/* Get either sum or maximum of right flanking sequence lengths */
static int get_tot_right(const alignment *a, const mfasta *m, int want_sum) {
  int tot = 0;
  int i;
  size_t j = 0;

  for (i = 0; i < m->seqnum; ++i) {
    const fasta *f = &m->f[i];
    const aligned_seq *s = &a->seqs[j];

    if (j != a->seq_num && strcmp(f->title, s->name) == 0) {
      const int left = xatoi(s->start) - 1;  /* zero-based coordinate */
      const int seq_len = strccnt(s->seq, '.');
      const int seq_end = left + seq_len;
      const int right = f->seqlen - seq_end;
      if (want_sum)
	tot += right;
      else  /* want maximum */
	if (right > tot)
	  tot = right;
      ++j;
    }
  }

  return tot;
}

/* Stick flanking sequences onto the left and right of an alignment */
static void add_flanks(alignment *a, const mfasta *m, int want_extend) {
  const int left_tot = aln_tot_left(a, want_extend);
  const int right_tot = get_tot_right(a, m, want_extend);
  int left_cumu = 0;  /* only used for extended alignments */
  int right_cumu = 0;  /* only used for extended alignments */
  int i;
  size_t j = 0;

  for (i = 0; i < m->seqnum; ++i) {
    const fasta *f = &m->f[i];
    aligned_seq *s = &a->seqs[j];

    if (j != a->seq_num && strcmp(f->title, s->name) == 0) {
      const int left = xatoi(s->start) - 1;  /* zero-based coordinate */
      const int aln_len = strlen(s->seq);
      const int seq_len = strccnt(s->seq, '.');
      const int aln_end = left_tot + aln_len;
      const int seq_end = left + seq_len;
      const int right = f->seqlen - seq_end;
      const int left_extra = left_tot - left_cumu - left;
      const int right_extra = right_tot - right_cumu - right;
      const int left_beg = left_extra;
      const int right_beg = aln_end + right_cumu;
      const int left_end = left_beg + left;
      const int right_end = right_beg + right;
      const int tot_len = aln_end + right_tot;
      char *new_seq = xmalloc(tot_len+1);  /* space for the final 0 */

      memset(new_seq, '.', left_extra);
      int2char(new_seq + left_beg, f->seq, left);
      memset(new_seq + left_end, '.', left_cumu);
      memcpy(new_seq + left_tot, s->seq, aln_len);
      memset(new_seq + aln_end, '.', right_cumu);
      int2char(new_seq + right_beg, f->seq + seq_end, right);
      memset(new_seq + right_end, '.', right_extra);
      new_seq[tot_len] = 0;  /* terminate the string */

      free(s->seq);
      s->seq = new_seq;
      left_cumu += left * want_extend;
      right_cumu += right * want_extend;
      ++j;
    }
  }

  if (j != a->seq_num)
    die("%s: didn't find %s among the sequences\n",
	prog_name, a->seqs[j].name);
  assert(left_cumu == left_tot * want_extend);
  assert(right_cumu == right_tot * want_extend);
}

/* Mangle sequence names in the same way as glam2 */
static void mangle_names(mfasta *m) {
  static const int name_width = 12;  /* max name width for glam2 */
  int i;
  for (i = 0; i < m->seqnum; ++i) {
    first_word(m->f[i].title);
    strtrunc(m->f[i].title, name_width);
  }
}

/* Read sequences and stick flanking bits onto an alignment */
static void get_flanks(alignment *a, const char *seq_file, int want_extend) {
  mfasta m;
  FILE *fp = xfopen(seq_file, "r");
  mfasta_read(&m, fp);  /* read the sequences */
  xfclose(fp);
  mangle_names(&m);  /* make sequence names match alignment names */
  add_flanks(a, &m, want_extend);
  free_mfasta(&m);
}

static void usage(void) {
  die("\
Usage: glam2format [options] my_format my_motif.glam2\n\
Formats: fasta, msf\n\
Options (default settings):\n\
-o: output file (stdout)\n\
-c: make a compact alignment\n\
-f: sequence file for flanking sequences\n\
");
}

int main(int argc, char **argv) {
  const char *format;
  const char *in_file;
  const char *out_file = "-";
  const char *seq_file = NULL;
  int want_extend = 1;
  FILE *fp;
  alignment aln;
  int c;

  prog_name = "glam2format";  /* for error messages */

  while ((c = getopt(argc, argv, "o:cf:")) != -1) {  /* non-ANSI */
    switch (c) {
    case 'o':
      out_file = optarg;
      break;
    case 'c':
      want_extend = 0;
      break;
    case 'f':
      seq_file = optarg;
      break;
    case '?':
      usage();
    }
  }

  if (optind != argc-2)
    usage();
  format = lowercase(argv[optind++]);
  in_file = argv[optind++];

  fp = xfopen(in_file, "r");
  aln_read(&aln, fp);
  xfclose(fp);

  if (want_extend) {
    if (aln.key_positions == NULL)
      die("%s: no motif found in %s\n", prog_name, in_file);
    if (!aln_same_lengths(&aln))
      die("%s: unequal aligned lengths in %s\n", prog_name, in_file);
    aln_extend(&aln);
  }

  if (seq_file)
    get_flanks(&aln, seq_file, want_extend);

  fp = xfopen(out_file, "w");
  if (strcmp(format, "fasta") == 0)
    aln_fasta_write(&aln, fp); 
  else if (strcmp(format, "msf") == 0)
    aln_msf_write(&aln, fp);
  else
    usage();
  xfclose(fp);

  aln_free(&aln);
  return 0;
}
