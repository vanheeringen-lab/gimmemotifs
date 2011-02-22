/* GLAM2 scanner */
#include <assert.h>
#include <float.h>  /* DBL_MAX */
#include <stdlib.h>  /* free */
#include "glam2_util.h"
#include "glam2_heap.h"
#include "glam2_scan.h"
#include "glam2_scan_init.h"
#include "glam2_scan_output.h"

VERBOSE_T verbosity;            // needed by meme utilities

/* print vector of doubles: for debugging */
void print_vec(double *vec, int size) {
  int i;
  for (i = 0; i < size; ++i)
    printf("%.2g%c", vec[i], i < size-1 ? '\t' : '\n');
}

/* print matrix of doubles: for debugging */
void print_mat(double **mat, int rows, int cols) {
  int i, j;
  for (i = 0; i < rows; ++i)
    for (j = 0; j < cols; ++j)
      printf("%.2g%c", mat[i][j], j < cols-1 ? '\t' : '\n');
}

double max3(double x, double y, double z) {
  if (x > y)
    if (x > z)
      return x;
    else
      return z;
  else
    if (y > z)
      return y;
    else
      return z;
}

void init_alignment(alignment *aln, const data *d, int strand) {
  const int width = d->m.width;
  const int *seq = strand == '+' ? d->f.seq : d->f.rcseq;
  const int *hp = d->hit_positions;
  const int *hm = d->hit_matches;
  const int start = hp[0];
  const int end = hp[width-1] + hm[width-1];
  const int deletions = width - sum_int(hm, width);
  const int aln_size = end - start + deletions;
  char *name = xstrdup(d->f.title);
  int *seq1 = xmalloc(aln_size * sizeof(int));
  int *seq2 = xmalloc(aln_size * sizeof(int));
  int aln_pos = 0;
  int seq_pos = start;
  int mot_pos;

  for (mot_pos = 0; mot_pos < width; ++mot_pos) {
    while (seq_pos < hp[mot_pos]) {
      seq2[aln_pos] = seq[seq_pos++];
      seq1[aln_pos++] = 0;
    }
    if (hm[mot_pos] == 1) {  /* match */
      assert(seq[seq_pos] < d->alph.size);  /* ambiguous residue forbidden */
      seq2[aln_pos] = seq[seq_pos++];
    } else  /* deletion */
      seq2[aln_pos] = d->alph.size;
    seq1[aln_pos++] = 1;
  }

  aln->name = name;
  aln->strand = strand;
  aln->coord1 = strand == '+' ? start+1 : d->f.seqlen - start;
  aln->coord2 = strand == '+' ? end : d->f.seqlen - end + 1;
  aln->aln_size = aln_size;
  aln->seq1 = seq1;
  aln->seq2 = seq2;
  aln->score = d->dp_mat[width][end];
}

void free_alignment(alignment *aln) {
  free(aln->name);
  free(aln->seq1);
  free(aln->seq2);
}

int cmp_alignment(const void *a, const void *b) {
  const alignment *x = a;
  const alignment *y = b;
  if (x->score < y->score)
    return +1;
  if (x->score > y->score)
    return -1;
  /* break ties by sequence name, so results don't depend on order of input */
  return strcmp(x->name, y->name);
}

void reinit_dp(data *d, int seqlen) {
  const int width = d->m.width;
  int i;

  if (seqlen > d->dp_seqlen) {
    for (i = 0; i <= width; ++i) {
      d->dp_mat[i] = XREALLOC(d->dp_mat[i], seqlen+1);
      d->forbidden[i] = XREALLOC(d->forbidden[i], seqlen+1);
    }
    set_dbl(d->dp_mat[0], seqlen+1, 0);  /* boundary condition */
    d->dp_seqlen = seqlen;
  }

  for (i = 0; i <= width; ++i)
    ZERO(d->forbidden[i], seqlen+1);
}

void dyn_prog(data *d, int strand) {
  const int width = d->m.width;
  const int seqlen = d->f.seqlen;
  const int *seq = strand == '+' ? d->f.seq : d->f.rcseq;
  const int *end = seq + seqlen;
  int i;

  for (i = 0; i < width; ++i) {
    /* tried to optimize this inner loop as much as possible */
    /* it seems surprisingly easy to beat gcc -O3 */
    const int *s = seq;
    const double *ms = d->match_scores[i];
    const double ds = d->delete_scores[i];
    const double is = d->insert_scores[i];
    const double *mat1 = d->dp_mat[i];
    double *mat2 = d->dp_mat[i+1];
    double m1 = *mat1;
    double m2 = *mat2;

    while (s != end) {
      ++mat1;
      ++mat2;
      const double match = ms[*s] + m1;
      m1 = *mat1;
      const double delete = ds + m1;
      m2 += is;
      if (match > m2)
	m2 = match;
      if (delete > m2)
	m2 = delete;
      *mat2 = m2;
      ++s;
    }
  }
}

int traceback(data *d, int end, int strand) {
  const int width = d->m.width;
  const int *seq = strand == '+' ? d->f.seq : d->f.rcseq;
  int i = width - 1;
  int j = end - 1;

  while (i != -1 && j != -1) {
    const double match = d->forbidden[i][j] == 0 ?
      d->match_scores[i][seq[j]] + d->dp_mat[i][j] : -DBL_MAX;
    const double delete = d->delete_scores[i] + d->dp_mat[i][j+1];
    const double insert = d->insert_scores[i] + d->dp_mat[i+1][j];

    if (match >= delete && match >= insert) {  /* match */
      assert(seq[j] < d->alph.size);  /* ambiguous residue forbidden */
      d->hit_matches[i] = 1;
      d->hit_positions[i] = j;
      d->forbidden[i][j] = 1;
      --i;
      --j;
    } else if (delete >= insert) {  /* deletion */
      d->hit_matches[i] = 0;
      d->hit_positions[i] = j+1;
      --i;
    } else {  /* insertion */
      assert(i+1 != width);
      --j;
    }
  }

  for (; i != -1; --i) {
    d->hit_matches[i] = 0;
    d->hit_positions[i] = 0;
  }

  return j+1;
}

void recalculate(data *d, int start, int strand) {
  const int width = d->m.width;
  const int seqlen = d->f.seqlen;
  const int *seq = strand == '+' ? d->f.seq : d->f.rcseq;
  int end = start;
  int i, j;

  for (i = 0; i < width; ++i)
    for (j = start; j < seqlen; ++j) {
      const double match = d->forbidden[i][j] == 0 ?
	d->match_scores[i][seq[j]] + d->dp_mat[i][j] : -DBL_MAX;
      const double delete = d->delete_scores[i] + d->dp_mat[i][j+1];
      const double insert = d->insert_scores[i] + d->dp_mat[i+1][j];
      const double new = max3(match, delete, insert);

      if (new != d->dp_mat[i+1][j+1]) {
	d->dp_mat[i+1][j+1] = new;
	if (end < j)
	  end = j;
      } else
	if (j > end)
	  break;
    }
}

void scan_seq(data *d, int strand) {
  const int width = d->m.width;
  const int seqlen = d->f.seqlen;

  reinit_dp(d, d->f.seqlen);
  dyn_prog(d, strand);

  while (1) {
    /*
    print_mat(d->dp_mat, width+1, seqlen+1);
    */
    const double *last_row = d->dp_mat[width];
    const int end = max_dbl(last_row, seqlen+1) - last_row;
    const int start = traceback(d, end, strand);
    /* require at least 1 match, to prevent infinite loop: */
    if (sum_int(d->hit_matches, width) == 0)
      break;
    if (d->hit_num == d->a.hit_num) {
      if (d->hit_num == 0 || last_row[end] <= d->hits[0].score)
	break;
      POP_HEAP(d->hits, d->hit_num, cmp_alignment);
      --d->hit_num;
      free_alignment(d->hits + d->hit_num);
    }
    init_alignment(d->hits + d->hit_num, d, strand);
    ++d->hit_num;
    PUSH_HEAP(d->hits, d->hit_num, cmp_alignment);
    recalculate(d, start, strand);
  }
}

int main(int argc, char **argv) {
  data d;
  FILE *fp;

  prog_name = "glam2scan";  /* for error messages */
  getargs(&d.a, argc, argv);
  init(&d);

  fputs("GLAM2SCAN\nVersion "
#include "glam2_version.h"
	"\n\n", d.out);
  printargs(d.out, argc, argv);
  putc('\n', d.out);

  /*
  print_vec(d.delete_scores, d.m.width);
  print_vec(d.insert_scores, d.m.width);
  print_mat(d.match_scores, d.m.width, d.alph.size);
  */

  fp = xfopen(d.a.seq_file, "r");
  while (fasta_read(&d.f, fp) != EOF) {
    first_word(d.f.title);
    tr_fasta(&d.f, d.alph.encode);
    scan_seq(&d, '+');
    if (d.a.two_strands) {
      rc_fasta(&d.f, d.alph.size);
      scan_seq(&d, '-');
    }
    free_fasta(&d.f);
  }
  xfclose(fp);

  SORT(d.hits, d.hit_num, cmp_alignment);
  print_hits(d.out, d.hits, &d);
  return 0;
}
