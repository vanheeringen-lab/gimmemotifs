/* GLAM2 */
#include <assert.h>
#include <float.h>  /* DBL_MIN, DBL_EPSILON */
#include "glam2_util.h"
#include "glam2_glam2.h"
#include "glam2_init.h"
#include "glam2_output.h"
#include "glam2_site_sample.h"
#include "glam2_column_sample.h"
#include "dir.h"		// MEME specific distribution

VERBOSE_T verbosity;		// needed by meme utilities

/* log probability of counts given probabilities */
double prob_vec_score(const prob_vec *p, const int *counts) {
  int i;
  double score = 0;
  for (i = 0; i < p->dim; ++i)
    score += counts[i] * p->log_probs[i];
  return score;
}

double bg_score(const glam2_scorer *s, const int *counts) {
  return prob_vec_score(&s->bg, counts);
}

double emission_score(glam2_scorer *s, const int *counts) {
  return dmix_score(&s->e_prior, counts);
}

double deletion_score(glam2_scorer *s, int delete_count, int no_delete_count) {
  return beta_score(&s->d_prior, delete_count, no_delete_count);
}

double insertion_score(glam2_scorer *s, int insert_count, int no_insert_count) {
  return beta_score(&s->i_prior, insert_count, no_insert_count);
}

double column_score(glam2_scorer *s, const glam2_col *col) {
  return emission_score(s, col->emission_counts)
    + deletion_score(s, col->delete_count, col->match_count)
    - bg_score(s, col->emission_counts);
}

/* Calculate the score of an alignment */
double aln_score(glam2_scorer *s, const glam2_aln *aln) {
  double score = 0;
  int i;

  for (i = 0; i < aln->width; ++i)
    score += column_score(s, &aln->cols[i]);

  for (i = 1; i < aln->width; ++i)
    score += insertion_score(s, aln->insert_counts[i-1], aln->aligned_seq);

  return score;
}

/* Calculate one sequence's contribution to the alignment score */
double marginal_score(glam2_scorer *s, glam2_aln *aln,
		      int seq, const fasta *f) {
  double score = aln->score;
  unalign(aln, seq, f);
  score -= aln_score(s, aln);
  realign(aln, seq, f);
  return score;
}

/* get a random starting alignment */
void start_aln(glam2_aln *aln, data *d) {
  int i;
#if 0
  aln->width = d->a.min_width;  /* ?? initial number of columns */
  aln->width = sqrt(d->a.max_width * d->a.min_width);  /* geometric mean */
#endif
  aln->width = d->a.init_width;
  aln_zero(aln);
  SHUFFLE(d->seq_order, aln->seq_num);
  for (i = 0; i < aln->seq_num; ++i)
    site_sample(aln, d->seq_order[i], d, 1);
  aln->score = aln_score(&d->scorer, aln);
}

void update_aln(glam2_aln *aln, data *d, const double temperature) {
  assert(aln->seq_num > 0);
  if (rand_dbl(d->a.column_sample_rate + 1) < 1) {
    const int seq_pick = rand_int(aln->seq_num);
    if (d->a.profile)
      fprintf(d->out, "site sample: seq=%d\n", seq_pick);
    site_sample(aln, seq_pick, d, temperature);
  } else
    column_sample(aln, d, temperature);
  aln->score = aln_score(&d->scorer, aln);
}

void optimise_aln(glam2_aln *best, data *d) {
  glam2_aln *aln = &d->aln;
  int no_improvement = 0;
  int i;

  aln_copy(aln, best);
  if (d->a.profile)
    fputs("Temperature, Columns, Sequences, Score:\n", d->out);

  for (i = 0; no_improvement < d->a.stop_after; ++i) {
    double temperature = d->a.temperature /
      xpow(d->a.cool, (double)i / d->a.stop_after);
    if (temperature < d->a.frozen)
      temperature = d->a.frozen;
    if (d->a.profile)
      fprintf(d->out, "%g\t%d\t%d\t%g\n",
	      temperature, aln->width, aln->aligned_seq, aln->score / xlog(2));
    /*
    print_aln(d->out, aln, d);
    */
    update_aln(aln, d, temperature);
    if (aln->score > best->score) {
      aln_copy(best, aln);
      no_improvement = 0;
    } else
      ++no_improvement;
  }

  if (d->a.profile)
    putc('\n', d->out);
  if (!d->a.quiet) fprintf(stderr, "%d iterations\n", i);
}

void print_misc_info(FILE *fp, const data *d) {
  const int alph_size = d->alph.size;
  int i;
  fprintf(fp, "Sequences: %d\n", d->seqs.seqnum);
  fprintf(fp, "Greatest sequence length: %d\n", d->seqs.maxlen);
  fputs("Residue counts: ", fp);
  for (i = 0; i <= alph_size; ++i)
    fprintf(fp, "%c=%d%c", d->alph.decode[i], d->scorer.bg.counts[i],
	    i < alph_size ? ' ' : '\n');
}

/* Alignment comparison function for sorting */
int aln_cmp(const void *a, const void *b) {
  const double x = ((const glam2_aln *)a)->score;
  const double y = ((const glam2_aln *)b)->score;
  return x < y ? +1 : x > y ? -1 : 0;
}

int main(int argc, char **argv) {
  data d;
  glam2_aln *alns;
  int r;

  prog_name = "glam2";  /* for error messages */
  getargs(&d.a, argc, argv);
  init(&d);

  fputs("GLAM2: Gapped Local Alignment of Motifs\nVersion "
#include "glam2_version.h"
	"\n\n", d.out);
  printargs(d.out, argc, argv);
  print_misc_info(d.out, &d);
  putc('\n', d.out);
  XMALLOC(alns, d.a.runs);

  for (r = 0; r < d.a.runs; ++r) {
    glam2_aln *aln = &alns[r];
    if (!d.a.quiet) {
      fprintf(stderr, "Run %d... ", r+1);
      fflush(stderr);
    }
    aln_init(aln, d.seqs.seqnum, d.a.max_width, d.alph.size);
    d.sm.underflow_flag = 1;  /* do we care about underflow in start_aln? */
    start_aln(aln, &d);
    optimise_aln(aln, &d);
    if (d.sm.underflow_flag < (d.a.algorithm == 2 ? DBL_EPSILON : DBL_MIN))
      fprintf(stderr, "%s: accuracy loss due to numeric underflow (%g)\nIf the alignment looks suspect, try rerunning with higher -u, or maybe lower -b\n", prog_name, d.sm.underflow_flag);
    if (d.a.profile)
      print_aln_info(d.out, aln, &d);
  }

  if (!d.a.quiet) putc('\n', stderr);

  SORT(alns, d.a.runs, aln_cmp);
  if (!d.a.profile)
    print_alns(d.out, alns, &d);

  xfclose(d.out);			// close text output file

  // Create the HTML output and MEME format output
  // This depends on variable "MEME_DIR" containing the
  // path to the installation directory.
  int command_length = strlen(MEME_DIR) + strlen(d.txt_filename) + 
    strlen(d.html_filename) + 50;
  char *command = xmalloc(command_length);
  sprintf(command, "%s/bin/glam2html < %s > %s", 
    MEME_DIR, d.txt_filename, d.html_filename);
  system(command);
  sprintf(command, "%s/bin/glam2psfm < %s > %s", 
    MEME_DIR, d.txt_filename, d.psfm_filename);
  system(command);

  return 0;
}
