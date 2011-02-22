/* Site sampling functions for GLAM2 */
#include <assert.h>
#include <errno.h>
#include <float.h>
#include "glam2_util.h"
#include "glam2_site_sample.h"

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

/* Remove one sequence from an alignment's count matrices */
void unalign(glam2_aln *aln, const int seq_pick, const fasta *f) {
  const int *seq = aln->strands[seq_pick] == '+' ? f->seq : f->rcseq;
  int i;
  if (!ALIGNED(aln, seq_pick))  /* the sequence is unaligned */
    return;
  --aln->aligned_seq;

  for (i = 0; i < aln->width; ++i) {
    glam2_col *col = &aln->cols[i];
    if (col->matches[seq_pick]) {
      const int pos = col->positions[seq_pick];
      assert(seq[pos] < col->alph_size);
      --col->match_count;
      --col->emission_counts[seq[pos]];
    } else
      --col->delete_count;
    aln->insert_counts[i] -= RIGHT_INSERT(aln, i, seq_pick);
  }
}

/* Add one sequence to an alignment's count matrices */
void realign(glam2_aln *aln, const int seq_pick, const fasta *f) {
  const int *seq = aln->strands[seq_pick] == '+' ? f->seq : f->rcseq;
  int i;
  if (!ALIGNED(aln, seq_pick))  /* the sequence is unaligned */
    return;
  ++aln->aligned_seq;

  for (i = 0; i < aln->width; ++i) {
    glam2_col *col = &aln->cols[i];
    if (col->matches[seq_pick]) {
      const int pos = col->positions[seq_pick];
      assert(seq[pos] < col->alph_size);
      ++col->match_count;
      ++col->emission_counts[seq[pos]];
    } else
      ++col->delete_count;
    aln->insert_counts[i] += RIGHT_INSERT(aln, i, seq_pick);
  }
}

static void get_score_matrices(score_matrices *sm, const glam2_aln *aln,
			       glam2_scorer *s) {
  double **match_scores = sm->match_scores;
  double *delete_scores = sm->delete_scores;
  double *insert_scores = sm->insert_scores;
  int i, j;

  for (i = 0; i < aln->width; ++i) {  /* loop over aligned columns */
    const int alph_size = aln->cols[i].alph_size;  /* should be in sm */
    const int *residue_counts = aln->cols[i].emission_counts;
    const int match_count = aln->cols[i].match_count;
    const int delete_count = aln->cols[i].delete_count;
    const int insert_count = aln->insert_counts[i];
    const int aligned_seq = aln->aligned_seq;

    const double ds = beta_ratio_a(&s->d_prior, delete_count, match_count);
    const double ms = beta_ratio_b(&s->d_prior, delete_count, match_count);
    const double is = beta_ratio_a(&s->i_prior, insert_count, aligned_seq);
    const double js = beta_ratio_b(&s->i_prior, insert_count, aligned_seq);

    assert(residue_counts[alph_size] == 0);  /* no mask characters */

    dmix_ratios(&s->e_prior, match_scores[i], residue_counts);

    if (i+1 < aln->width) {
      for (j = 0; j < alph_size; ++j)
	match_scores[i][j] *= js * ms / s->bg.probs[j];
      delete_scores[i] = ds * js;
      insert_scores[i] = is;
    } else {  /* no insertions after last column: special case */
      for (j = 0; j < alph_size; ++j)
	match_scores[i][j] *= ms / s->bg.probs[j];
      delete_scores[i] = ds;
      insert_scores[i] = 0;
    }
  }
}

static void heat_scores(score_matrices *sm, int width, int alph_size,
			double temperature) {
  int i;
  for (i = 0; i < width; ++i)
    pow_dbl(sm->match_scores[i], alph_size, 1 / temperature);
  pow_dbl(sm->delete_scores, width, 1 / temperature);
  pow_dbl(sm->insert_scores, width, 1 / temperature);  /* back is 0 */
}

static void copy_scores(score_matrices *sm, int width, int alph_size) {
  int i;
  for (i = 0; i < width; ++i)
    COPY(sm->rc_match_scores[i], sm->match_scores[i], alph_size);
  COPY(sm->rc_delete_scores, sm->delete_scores, width);
}

/* Stochastically choose a strand (or "unaligned"), using the DP matrix */
int pick_strand(const glam2_aln *aln, const int seq_pick, const data *d) {
  const int width = aln->width;
  const int seq_len = d->seqs.f[seq_pick].seqlen;
  double **dp_mat = d->sm.dp_mat;
  double **rc_mat = d->sm.rc_mat;
  const double dp_rescale = d->sm.dp_rescale;
  const double rc_rescale = d->sm.rc_rescale;
  double tot[3];
  double log_scale[3];
  double max;
  int i, r;

  tot[0] = sum_dbl(dp_mat[width], seq_len+1);
  assert(tot[0] > 0);
  log_scale[0] = dp_rescale;

  if (d->a.two_strands) {
    tot[1] = sum_dbl(rc_mat[width], seq_len+1);
    assert(tot[1] > 0);
    log_scale[1] = rc_rescale;
  } else {
    tot[1] = 0;
    log_scale[1] = dp_rescale;
  }

  if (aln->aligned_seq >= d->a.min_seqs) {
    tot[2] = 1;
    log_scale[2] = 0;
  } else {
    tot[2] = 0;
    log_scale[2] = dp_rescale;
  }

  /*
  puts("tot:");
  print_vec(tot, 3);
  puts("log_scale:");
  print_vec(log_scale, 3);
  */

  max = *max_dbl(log_scale, 3);
  for (i = 0; i < 3; ++i)
    tot[i] *= xexp(log_scale[i] - max);
  r = pick_dbl(tot, 3) - tot;

  if (r == 0)
    return '+';
  else if (r == 1)
    return '-';
  else
    return 0;
}

/* Stochastically choose an alignment endpoint, using the DP matrix */
static int pick_endpoint(const score_matrices *sm, int width,
			 int seq_len, int strand) {
  double **dp_mat = strand == '+' ? sm->dp_mat : sm->rc_mat;

  if (strand == 0)  /* the sequence is unaligned */
    return -1;  /* dummy value, hopefully never used */

  return pick_dbl(dp_mat[width], seq_len+1) - dp_mat[width];
}

/* Stochastic traceback: */
void traceback(glam2_aln *aln, int seq_pick, int strand, int end, data *d) {
  const fasta *f = &d->seqs.f[seq_pick];
  const int *seq = strand == '+' ? f->seq : f->rcseq;
  glam2_col *cols = aln->cols;
  double **ms = strand == '+' ? d->sm.match_scores : d->sm.rc_match_scores;
  double *ds = strand == '+' ? d->sm.delete_scores : d->sm.rc_delete_scores;
  double **dp_mat = strand == '+' ? d->sm.dp_mat : d->sm.rc_mat;
  int i = aln->width - 1;
  int j = end - 1;

  aln->strands[seq_pick] = strand;
  if (strand == 0)  /* the sequence is unaligned */
    return;

  while (i != -1 && j != -1) {
    const double r = rand_dbl(dp_mat[i+1][j+1]);
    const double match = ms[i][seq[j]] * dp_mat[i][j];
    const double delete = ds[i] * dp_mat[i][j+1];

    if (dp_mat[i+1][j+1] < d->sm.underflow_flag)
      d->sm.underflow_flag = dp_mat[i+1][j+1];

    if (r < match) {  /* match */
      assert(seq[j] < d->alph.size);  /* ambiguous residue forbidden */
      cols[i].matches[seq_pick] = 1;
      cols[i].positions[seq_pick] = j;
      --i;
      --j;
    } else if (r < match + delete) {  /* deletion */
      cols[i].matches[seq_pick] = 0;
      cols[i].positions[seq_pick] = j+1;
      --i;
    } else {  /* insertion */
      assert(i+1 < aln->width);
      --j;
    }
  }

  for (; i != -1; --i) {
    if (dp_mat[i+1][j+1] < d->sm.underflow_flag)
      d->sm.underflow_flag = dp_mat[i+1][j+1];
    cols[i].matches[seq_pick] = 0;
    cols[i].positions[seq_pick] = 0;
  }
}

/* Stochastic traceback with nonlinear insertion scores */
void traceback_slow(glam2_aln *aln, int seq_pick, int strand, int end, data *d,
		    double temperature) {
  const fasta *f = &d->seqs.f[seq_pick];
  const int *seq = strand == '+' ? f->seq : f->rcseq;
  glam2_col *cols = aln->cols;
  double **ms = strand == '+' ? d->sm.match_scores : d->sm.rc_match_scores;
  double *ds = strand == '+' ? d->sm.delete_scores : d->sm.rc_delete_scores;
  double **dp_mat = strand == '+' ? d->sm.dp_mat : d->sm.rc_mat;
  int i = aln->width - 1;
  int j = end - 1;

  aln->strands[seq_pick] = strand;
  if (strand == 0)  /* the sequence is unaligned */
    return;

  while (i != -1 && j != -1) {
    double r = rand_dbl(dp_mat[i+1][j+1]);
    int a_count = aln->insert_counts[i];
    int b_count = aln->aligned_seq+1;  /* +1 */
    double insert_score = 1;

    if (dp_mat[i+1][j+1] < d->sm.underflow_flag)
      d->sm.underflow_flag = dp_mat[i+1][j+1];

    for (; j != -1; --j) {
      const double is = xpow(insert_score, 1 / temperature);
      const double match = ms[i][seq[j]] * dp_mat[i][j] * is;
      const double delete = ds[i] * dp_mat[i][j+1] * is;

      if (r < match) {  /* match */
	assert(seq[j] < d->alph.size);  /* ambiguous residue forbidden */
	cols[i].matches[seq_pick] = 1;
	cols[i].positions[seq_pick] = j;
	--i;
	--j;
	break;
      } else if (r < match + delete) {  /* deletion */
	cols[i].matches[seq_pick] = 0;
	cols[i].positions[seq_pick] = j+1;
	--i;
	break;
      }
      assert(i+1 < aln->width);
      r -= match + delete;
      insert_score *= beta_ratio_a(&d->scorer.i_prior, a_count, b_count);
      ++a_count;
    }
  }

  for (; i != -1; --i) {
    if (dp_mat[i+1][j+1] < d->sm.underflow_flag)
      d->sm.underflow_flag = dp_mat[i+1][j+1];
    cols[i].matches[seq_pick] = 0;
    cols[i].positions[seq_pick] = 0;
  }
}

/* Fill in the DP matrix, using the "forward algorithm" (sum of paths): */
void forward(const glam2_aln *aln, int seq_pick, int strand, data *d) {
  const fasta *f = &d->seqs.f[seq_pick];
  const int *seq = strand == '+' ? f->seq : f->rcseq;
  const int *end = seq + f->seqlen;
  double **ms = strand == '+' ? d->sm.match_scores : d->sm.rc_match_scores;
  double *ds = strand == '+' ? d->sm.delete_scores : d->sm.rc_delete_scores;
  double *is = d->sm.insert_scores;
  double **dp_mat = strand == '+' ? d->sm.dp_mat : d->sm.rc_mat;
  double *rescale = strand == '+' ? &d->sm.dp_rescale : &d->sm.rc_rescale;
  double r = 1;  /* rescaling factor for one row */
  int i;
  assert(is[aln->width-1] == 0);
  *rescale = 0;

  for (i = 0; i < aln->width; ++i) {
    if (r <= 0)
      die("%s: underflow in alignment i=%d\n", prog_name, i);
    div_dbl(ms[i], d->alph.size, r);
    ds[i] /= r;
    *rescale += xlog(r);

    /* tried to optimize this inner loop as much as possible */
    const int *s = seq;
    const double *msi = ms[i];
    const double dsi = ds[i];
    const double isi = is[i];
    const double *mat1 = dp_mat[i];
    double *mat2 = dp_mat[i+1];
    double m1 = *mat1;
    double m2 = dsi * m1;
    *mat2 = m2;
    r = m2;

    while (s != end) {
      ++mat1;
      ++mat2;
      m2 *= isi;
      m2 += msi[*s] * m1;
      m1 = *mat1;
      m2 += dsi * m1;
      *mat2 = m2;

      if (m2 > r)
	r = m2;
      if (m2 >= DBL_MAX)  /* this check doesn't seem to slow it down */
	die("%s: overflow %g in alignment i=%d isi=%.3g dsi=%.3g msi=%.3g m1=%.3g\n",
	    prog_name, m2, i, isi, dsi, msi[*s], m1);

      ++s;
    }
  }
}

/* Make vector of probabilities for different insertion sizes */
void beta_convolver_a(const beta *b, double *convolver, int size,
		      int a_count, int b_count) {
  double x = 1;
  int i;
  for (i = 0; i < size; ++i) {
    convolver[i] = x;
    x *= beta_ratio_a(b, a_count, b_count);
    ++a_count;
  }
}

/* z = convolution of x and y */
void convolve(double *z, const double *x, const double *y, const int size) {
  int i, j;
  for (i = 0; i < size; ++i) {
    z[i] = 0;
    for (j = 0; j <= i; ++j)
      z[i] += x[j] * y[i-j];
  }
}

/* Move this to util? */
void truncate_min_dbl(double *ptr, const size_t size, const double min) {
  const double *end = ptr + size;
  for (; ptr != end; ++ptr)
    if (*ptr < min)
      *ptr = min;
}

/* Fill in the DP matrix, using nonlinear insertion scores */
void forward_slow(const glam2_aln *aln, int seq_pick, int strand, data *d,
		  double temperature) {
  const fasta *f = &d->seqs.f[seq_pick];
  const int *seq = strand == '+' ? f->seq : f->rcseq;
  double **ms = strand == '+' ? d->sm.match_scores : d->sm.rc_match_scores;
  double *ds = strand == '+' ? d->sm.delete_scores : d->sm.rc_delete_scores;
  double **dp_mat = strand == '+' ? d->sm.dp_mat : d->sm.rc_mat;
  double *convolver = d->sm.convolver;
  double *convolved = d->sm.convolved;
  double *rescale = strand == '+' ? &d->sm.dp_rescale : &d->sm.rc_rescale;
  int i, j;
  *rescale = 0;

  for (i = 0; i < aln->width; ++i) {
    const double r = *max_dbl(dp_mat[i], f->seqlen+1);
    if (r <= 0)
      die("%s: underflow in alignment i=%d\n", prog_name, i);
    div_dbl(ms[i], d->alph.size, r);
    ds[i] /= r;
    *rescale += xlog(r);

    dp_mat[i+1][0] = ds[i] * dp_mat[i][0];

    for (j = 0; j < f->seqlen; ++j) {
      dp_mat[i+1][j+1] =
        ms[i][seq[j]] * dp_mat[i][j]
	+ ds[i] * dp_mat[i][j+1];

      if (dp_mat[i+1][j+1] >= DBL_MAX)
        die("%s: overflow %g in alignment i=%d j=%d M=%.3g*%.3g D=%.3g*%.3g\n",
            prog_name, dp_mat[i+1][j+1], i, j, ms[i][seq[j]],
            dp_mat[i][j], ds[i], dp_mat[i][j+1]);
    }

    if (i+1 < aln->width) {  /* no insertions after last column */
      beta_convolver_a(&d->scorer.i_prior, convolver, f->seqlen+1,
		       aln->insert_counts[i], aln->aligned_seq+1);  /* +1 */
      pow_dbl(convolver, f->seqlen+1, 1 / temperature);
      if (d->a.algorithm == 2) {
#ifdef FFT
	fft_convolve(dp_mat[i+1], dp_mat[i+1], convolver, f->seqlen+1, &d->sm.fft);
#endif
	/* set negative numbers produced by FFT to zero: big fat kludge */
	truncate_min_dbl(dp_mat[i+1], f->seqlen+1, 0);
      } else {
	convolve(convolved, dp_mat[i+1], convolver, f->seqlen+1);
	COPY(dp_mat[i+1], convolved, f->seqlen+1);
      }
    }
  }
}

#define FORWARD(aln, seq_pick, strand, d, temperature)\
        ((d)->a.algorithm == 0 ?\
         forward(aln, seq_pick, strand, d) :\
	 forward_slow(aln, seq_pick, strand, d, temperature))

#define TRACEBACK(aln, seq_pick, strand, end, d, temperature)\
        ((d)->a.algorithm == 0 ?\
         traceback(aln, seq_pick, strand, end, d) :\
         traceback_slow(aln, seq_pick, strand, end, d, temperature))

/* Stochastically re-align one sequence */
void site_sample(glam2_aln *aln, int seq_pick, data *d, double temperature) {
  int strand, end;
  const fasta *f = &d->seqs.f[seq_pick];
  unalign(aln, seq_pick, f);
  get_score_matrices(&d->sm, aln, &d->scorer);
  heat_scores(&d->sm, aln->width, d->alph.size, temperature);
  if (d->a.two_strands)
    copy_scores(&d->sm, aln->width, d->alph.size);
  FORWARD(aln, seq_pick, '+', d, temperature);
  if (d->a.two_strands)
    FORWARD(aln, seq_pick, '-', d, temperature);
  strand = pick_strand(aln, seq_pick, d);
  end = pick_endpoint(&d->sm, aln->width, f->seqlen, strand);
  /*
  puts("Delete scores:");
  print_vec(d->sm.delete_scores, aln->width);
  puts("Insert scores:");
  print_vec(d->sm.insert_scores, aln->width);
  puts("Match scores:");
  print_mat(d->sm.match_scores, aln->width, d->alph.size);
  puts("DP matrix:");
  print_mat(d->sm.dp_mat, aln->width+1, f->seqlen+1);
  printf("strand=%c end=%d\n", strand ? strand : '0', end);
  */
  TRACEBACK(aln, seq_pick, strand, end, d, temperature);
  realign(aln, seq_pick, f);
}
