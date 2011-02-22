#include <assert.h>
#include "glam2_util.h"
#include "glam2_glam2_aln.h"

void col_init(glam2_col *col, int seq_num, int alph_size) {
  col->seq_num = seq_num;
  XMALLOC(col->positions, seq_num);
  XMALLOC(col->matches, seq_num);
  col->alph_size = alph_size;
  XMALLOC(col->emission_counts, alph_size + 1);
}

void aln_init(glam2_aln *aln, int seq_num, int max_width, int alph_size) {
  int i;
  aln->seq_num = seq_num;
  XMALLOC(aln->cols, max_width);
  for (i = 0; i < max_width; ++i)
    col_init(&aln->cols[i], seq_num, alph_size);
  XMALLOC(aln->strands, seq_num);
  XMALLOC(aln->insert_counts, max_width);  /* back is unused */
}

static void col_zero(glam2_col *col) {
  ZERO(col->emission_counts, col->alph_size + 1);
  col->match_count = 0;
  col->delete_count = 0;
}

void aln_zero(glam2_aln *aln) {
  int i;
  ZERO(aln->strands, aln->seq_num);  /* set all sequences unaligned */
  aln->aligned_seq = 0;
  ZERO(aln->insert_counts, aln->width);  /* last element isn't really used */
  for (i = 0; i < aln->width; ++i)
    col_zero(&aln->cols[i]);
}

void col_copy(glam2_col *a, const glam2_col *b) {
  assert(a->seq_num == b->seq_num);
  assert(a->alph_size == b->alph_size);
  COPY(a->positions, b->positions, b->seq_num);
  COPY(a->matches, b->matches, b->seq_num);
  COPY(a->emission_counts, b->emission_counts, b->alph_size+1);
  a->match_count = b->match_count;
  a->delete_count = b->delete_count;
}

void aln_copy(glam2_aln *a, const glam2_aln *b) {
  int i;
  assert(a->seq_num == b->seq_num);
  a->width = b->width;
  for (i = 0; i < b->width; ++i)
    col_copy(&a->cols[i], &b->cols[i]);
  COPY(a->strands, b->strands, b->seq_num);
  a->aligned_seq = b->aligned_seq;
  COPY(a->insert_counts, b->insert_counts, b->width);
  a->score = b->score;
}
