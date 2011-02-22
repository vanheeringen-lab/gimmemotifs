/* This code seems a bit messy and complicated - can it be simplified? */
#include <assert.h>
#include <float.h>
#include "glam2_util.h"
#include "glam2_column_sample.h"

int sample(double *ptr, int size, const double temperature) {
  const double max = *max_dbl(ptr, size);
  const double *end = ptr + size;
  double *i;
  assert(size > 0);
  for (i = ptr; i != end; ++i)
    *i = xexp(*i - max);
  pow_dbl(ptr, size, 1 / temperature);
  return pick_dbl(ptr, size) - ptr;
}

/* Set positions of new column to the left of old column */
void col_put_left(column_sampler *s, const glam2_aln *aln, const int col) {
  int i;
  assert(col < aln->width);
  for (i = 0; i < aln->seq_num; ++i)
    if (ALIGNED(aln, i))
      s->col.positions[i] = LEFT_POS(&aln->cols[col], i) - s->offsets[i]
	- s->col.matches[i];
    else
      s->col.positions[i] = 0;  /* dummy value */
}

/* Set positions of new column to the right of old column */
void col_put_right(column_sampler *s, const glam2_aln *aln, const int col) {
  int i;
  assert(col < aln->width);
  for (i = 0; i < aln->seq_num; ++i)
    if (ALIGNED(aln, i))
      s->col.positions[i] = RIGHT_POS(&aln->cols[col], i) + s->offsets[i];
    else
      s->col.positions[i] = 0;  /* dummy value */
}

/* Count residues in one column of an alignment */
void count_residues(glam2_col *col, const glam2_aln *aln, const mfasta *m) {
  int i;
  ZERO(col->emission_counts, col->alph_size+1);

  for (i = 0; i < col->seq_num; ++i)
    if (ALIGNED(aln, i) && col->matches[i]) {
      const int *seq = aln->strands[i] == '+' ? m->f[i].seq : m->f[i].rcseq;
      ++col->emission_counts[seq[col->positions[i]]];
    }
}

/* Count insertions to the right of a column. Returns 0 for rightmost column */
int count_right_insert(const glam2_aln *aln, const int col) {
  int count = 0;
  int i;
  assert(col < aln->width);
  for (i = 0; i < aln->seq_num; ++i)
    if (ALIGNED(aln, i))
      count += RIGHT_INSERT(aln, col, i);
  return count;
}

int sum_offsets(const int *offsets, const glam2_aln *aln) {
  int sum = 0;
  int i;
  for (i = 0; i < aln->seq_num; ++i)
    if (ALIGNED(aln, i))
      sum += offsets[i];
  return sum;
}

/* Get offsets between columns col-1 and col */
void get_offsets(int *offsets, const glam2_aln *aln, const int col) {
  int i;
  int min = INT_MAX;
  assert(col > 0 && col < aln->width);
  for (i = 0; i < aln->seq_num; ++i)
    if (ALIGNED(aln, i)) {
      offsets[i] = LEFT_INSERT(aln, col, i);
      if (offsets[i] < min)
	min = offsets[i];
    } else
      offsets[i] = 0;  /* dummy value */
  add_int(offsets, aln->seq_num, -min);
}

/* Get the index of the first aligned sequence */
int first_aligned_seq(const glam2_aln *aln) {
  int i;
  for (i = 0; i < aln->seq_num; ++i)
    if (ALIGNED(aln, i))
      break;
  assert(i != aln->seq_num);
  return i;
}

/* Check whether matches/deletions and right offsets are identical */
int col_comp_left(const column_sampler *s, const glam2_aln *aln, int col) {
  int i = first_aligned_seq(aln);
  int off = RIGHT_INSERT(aln, col, i) - s->offsets[i];
  assert(col != aln->width-1);
  for (; i < aln->seq_num; ++i)
    if (ALIGNED(aln, i))
      if (s->col.matches[i] != aln->cols[col].matches[i] ||
	  RIGHT_INSERT(aln, col, i) - s->offsets[i] != off)
	return 0;
  return 1;
}

/* Check whether matches/deletions and left offsets are identical */
int col_comp_right(const column_sampler *s, const glam2_aln *aln, int col) {
  int i = first_aligned_seq(aln);
  int off = LEFT_INSERT(aln, col, i) - s->offsets[i];
  assert(col != 0);
  for (; i < aln->seq_num; ++i)
    if (ALIGNED(aln, i))
      if (s->col.matches[i] != aln->cols[col].matches[i] ||
	  LEFT_INSERT(aln, col, i) - s->offsets[i] != off)
	return 0;
  return 1;
}

/* Count matching columns */
int col_count_left(const column_sampler *s, const glam2_aln *aln) {
  int count = 0;
  int i;
  for (i = 0; i < aln->width-1; ++i)
    count += col_comp_left(s, aln, i);
  return count;
}

/* Count matching columns */
int col_count_right(const column_sampler *s, const glam2_aln *aln) {
  int count = 0;
  int i;
  for (i = 1; i < aln->width; ++i)
    count += col_comp_right(s, aln, i);
  return count;
}

/* How many columns can we fit between columns col-1 and col */
int col_fit(const column_sampler *s, const glam2_aln *aln, const int col) {
  int i;
  int min = INT_MAX;
  assert(col > 0 && col < aln->width);
  for (i = 0; i < aln->seq_num; ++i)
    if (ALIGNED(aln, i)) {
      const int fit = LEFT_INSERT(aln, col, i)
	- s->offsets[i] - s->col.matches[i] + 1;
      if (fit < min)
	min = fit;
    }
  return min > 0 ? min : 0;
}

/* How many columns can we fit before the first column? */
int col_fit_left(const column_sampler *s, const glam2_aln *aln) {
  int i;
  int min = INT_MAX;
  for (i = 0; i < aln->seq_num; ++i)
    if (ALIGNED(aln, i)) {
      const int fit = LEFT_POS(&aln->cols[0], i)
	- s->offsets[i] - s->col.matches[i] + 1;
      if (fit < min)
        min = fit;
    }
  return min > 0 ? min : 0;
}

/* How many columns can we fit after the last column? */
int col_fit_right(const column_sampler *s, const glam2_aln *aln,
		  const mfasta *m) {
  const int w = aln->width;
  int i;
  int min = INT_MAX;
  for (i = 0; i < aln->seq_num; ++i)
    if (ALIGNED(aln, i)) {
      const int fit = m->f[i].seqlen - RIGHT_POS(&aln->cols[w-1], i)
	- s->offsets[i] - s->col.matches[i] + 1;
      if (fit < min)
        min = fit;
    }
  return min > 0 ? min : 0;
}

void left_fits(column_sampler *s, const glam2_aln *aln) {
  int i;
  s->fits[0] = col_fit_left(s, aln);
  for (i = 1; i < aln->width; ++i)
    s->fits[i] = col_fit(s, aln, i);
}

void right_fits(column_sampler *s, const glam2_aln *aln, const mfasta *m) {
  int i;
  for (i = 1; i < aln->width; ++i)
    s->fits[i-1] = col_fit(s, aln, i);
  s->fits[aln->width-1] = col_fit_right(s, aln, m);
}

void col_scan_left(column_sampler *s, const glam2_aln *aln, data *d) {
  const int aligned_seq = aln->aligned_seq;
  const int off_tot = sum_offsets(s->offsets, aln);
  const int res_tot = off_tot + s->col.match_count;
  double *score_ptr = s->scores;
  int i, j;

  for (i = 0; i < aln->width; ++i) {
    const double old_score = i != 0 ?
      insertion_score(&d->scorer, aln->insert_counts[i-1], aligned_seq) : 0;
    int right_inserts = off_tot;
    int left_inserts = i != 0 ?
      aln->insert_counts[i-1] - res_tot : 0;
    col_put_left(s, aln, i);

    for (j = 0; j < s->fits[i]; ++j) {
      assert(score_ptr-s->scores < d->seqs.maxlen+d->a.max_width);
      count_residues(&s->col, aln, &d->seqs);
      *score_ptr = column_score(&d->scorer, &s->col)
	+ insertion_score(&d->scorer, right_inserts, aligned_seq);
      if (i != 0)
	*score_ptr += insertion_score(&d->scorer, left_inserts, aligned_seq)
	  - old_score;
      if (s->col.emission_counts[s->col.alph_size])
	*score_ptr = -DBL_MAX;  /* forbid ambiguous residues */
      add_int(s->col.positions, aln->seq_num, -1);  /* add -1 to positions */
      right_inserts += aligned_seq;
      left_inserts -= aligned_seq;
      ++score_ptr;
    }
  }
}

void col_scan_right(column_sampler *s, const glam2_aln *aln, data *d) {
  const int aligned_seq = aln->aligned_seq;
  const int off_tot = sum_offsets(s->offsets, aln);
  const int res_tot = off_tot + s->col.match_count;
  double *score_ptr = s->scores;
  int i, j;

  for (i = 0; i < aln->width; ++i) {
    const double old_score = i != aln->width-1 ?
      insertion_score(&d->scorer, aln->insert_counts[i], aligned_seq) : 0;
    int left_inserts = off_tot;
    int right_inserts = i != aln->width-1 ?
      aln->insert_counts[i] - res_tot : 0;
    col_put_right(s, aln, i);

    for (j = 0; j < s->fits[i]; ++j) {
      assert(score_ptr-s->scores < d->seqs.maxlen+d->a.max_width);
      count_residues(&s->col, aln, &d->seqs);
      *score_ptr = column_score(&d->scorer, &s->col)
	+ insertion_score(&d->scorer, left_inserts, aligned_seq);
      if (i != aln->width-1)
	*score_ptr += insertion_score(&d->scorer, right_inserts, aligned_seq)
	  - old_score;
      if (s->col.emission_counts[s->col.alph_size])
	*score_ptr = -DBL_MAX;  /* forbid ambiguous residues */
      add_int(s->col.positions, aln->seq_num, 1);  /* add 1 to positions */
      left_inserts += aligned_seq;
      right_inserts -= aligned_seq;
      ++score_ptr;
    }
  }
}

/* Rotate a column out of an alignment, and repair insert counts */
void rotate_out(glam2_aln *aln, const int col) {
  assert(col < aln->width);
  ROL(aln->cols + col, aln->width - col, 1);
  ROL(aln->insert_counts + col, aln->width - col, 1);
  --aln->width;
  if (col != 0)
    aln->insert_counts[col-1] = count_right_insert(aln, col-1);
}

/* Rotate a column back into an alignment, and repair insert counts */
void rotate_in(glam2_aln *aln, const int col) {
  ++aln->width;
  assert(col < aln->width);
  ROR(aln->cols + col, aln->width - col, 1);
  ROR(aln->insert_counts + col, aln->width - col, 1);
  if (col != 0)
    aln->insert_counts[col-1] = count_right_insert(aln, col-1);
  aln->insert_counts[col] = count_right_insert(aln, col);
}

int pick_new_col(column_sampler *s, int match_cols, double temperature) {
  const int tot_fit = sum_int(s->fits, s->width);
  int r, c;
  set_dbl(s->scores + tot_fit, match_cols, 0);
  /*
  printf("match_cols=%d tot_fit=%d scores=", match_cols, tot_fit);
  for (r = 0; r < tot_fit + match_cols; ++r) printf("%.3g ", s->scores[r]);
  printf("\n");
  */
  r = sample(s->scores, tot_fit + match_cols, temperature);
  for (c = 0; c < s->width; ++c) {
    if (r < s->fits[c])
      break;
    r -= s->fits[c];
  }
  if (c < s->width)
    add_int(s->offsets, s->seq_num, r);
  return c;
}

void sample_left(glam2_aln *aln, data *d, const double temperature) {
  column_sampler *s = &d->col_sampler;
  const int del_flag = rand_dbl(1) < s->del_probs[aln->width - 1];
  const int col_pick = rand_int(aln->width - 1);
  int match_cols, new_col;
  if (d->a.profile)
    fprintf(d->out, "column sample: side=left, delete=%d, col=%d\n",
	    del_flag, col_pick);
  get_offsets(s->offsets, aln, col_pick+1);
  col_copy(&s->col, &aln->cols[col_pick]);
  if (del_flag)
    rotate_out(aln, col_pick);  /* Remove column & decrement width */
  else if (aln->width == d->a.max_width)
    return;
  s->width = aln->width;
  left_fits(s, aln);
  col_scan_left(s, aln, d);
  match_cols = col_count_left(s, aln);
  if (aln->width < d->a.min_width)
    match_cols = 0;  /* force column re-insertion */
  new_col = pick_new_col(s, match_cols, temperature);
  if (new_col >= s->width)
    return;
  col_put_left(s, aln, new_col);
  count_residues(&s->col, aln, &d->seqs);
  col_copy(&aln->cols[aln->width], &s->col);
  rotate_in(aln, new_col);
}

void sample_right(glam2_aln *aln, data *d, const double temperature) {
  column_sampler *s = &d->col_sampler;
  const int del_flag = rand_dbl(1) < s->del_probs[aln->width - 1];
  const int col_pick = rand_int(aln->width - 1) + 1;
  int match_cols, new_col;
  if (d->a.profile)
    fprintf(d->out, "column sample: side=right, delete=%d, col=%d\n",
	    del_flag, col_pick);
  get_offsets(s->offsets, aln, col_pick);
  col_copy(&s->col, &aln->cols[col_pick]);
  if (del_flag)
    rotate_out(aln, col_pick);  /* Remove column & decrement width */
  else if (aln->width == d->a.max_width)
    return;
  s->width = aln->width;
  right_fits(s, aln, &d->seqs);
  col_scan_right(s, aln, d);
  match_cols = col_count_right(s, aln);
  if (aln->width < d->a.min_width)
    match_cols = 0;  /* force column re-insertion */
  new_col = pick_new_col(s, match_cols, temperature);
  if (new_col >= s->width)
    return;
  col_put_right(s, aln, new_col);
  count_residues(&s->col, aln, &d->seqs);
  col_copy(&aln->cols[aln->width], &s->col);
  rotate_in(aln, new_col+1);
}

void column_sample(glam2_aln *aln, data *d, double temperature) {
  assert(aln->width > 1);  /* need at least 2 columns */
  if (aln->aligned_seq == 0)  /* need at least 1 aligned sequence */
    return;
  if (rand_int(2))
    sample_left(aln, d, temperature);
  else
    sample_right(aln, d, temperature);
  assert(aln->width >= d->a.min_width && aln->width <= d->a.max_width);
}
