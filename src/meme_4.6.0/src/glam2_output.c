/* GLAM2 output functions */
#include <stdlib.h>  /* free */
#include "glam2_util.h"
#include "glam2_output.h"
#include "motif.h"
#include "user.h"		// for LOGOS
#include "ceqlogo.h"
#include "hash_alph.h"		// for DNA0

static const int name_width = 12;

/* Get maximum insertion between columns col and col+1 */
static int max_right_insert(const glam2_aln *aln, const int col) {
  int max = 0;
  int i;
  for (i = 0; i < aln->seq_num; ++i)
    if (ALIGNED(aln, i)) {
      int ins = RIGHT_INSERT(aln, col, i);
      if (ins > max)
	max = ins;
    }
  return max;
}

/* Print an alignment in a nice format */
void print_aln(FILE *fp, glam2_aln *aln, data *d) {
  const int width = aln->width;
  const glam2_col *cols = aln->cols;
  const fasta *seqs = d->seqs.f;
  const int coord_width = digits(d->seqs.maxlen);
  int *max_insert = xmalloc(width * sizeof(int));
  int i, j, k;

  /* get maximum insert size at each position: */
  for (i = 0; i < width; ++i)
    max_insert[i] = max_right_insert(aln, i);

  put_pad(' ', name_width + coord_width + 2, fp);
  for (j = 0; j < width; ++j) {
    putc('*', fp);
    put_pad('.', max_insert[j], fp);
  }
  putc('\n', fp);

  for (i = 0; i < aln->seq_num; ++i) {
    const int strand = aln->strands[i];
    const int *seq = strand == '+' ? seqs[i].seq : seqs[i].rcseq;
    const int start = LEFT_POS(&cols[0], i);  /* zero-based coordinate */
    const int end = RIGHT_POS(&cols[width-1], i);
    const int coord1 = strand == '+' ? start+1 : seqs[i].seqlen - start;
    const int coord2 = strand == '+' ? end : seqs[i].seqlen - end + 1;
    if (!ALIGNED(aln, i))
      continue;

    fprintf(fp, "%-*.*s", name_width, name_width, seqs[i].title);
    fprintf(fp, " %*d ", coord_width, coord1);

    for (j = 0; j < width; ++j) {
      int pos = LEFT_POS(&cols[j], i);
      if (cols[j].matches[i])
	putc(d->alph.decode[seq[pos++]], fp);
      else
	putc('.', fp);

      for (k = 0; k < RIGHT_INSERT(aln, j, i); ++k)
	putc(d->alph.decode[seq[pos++]], fp);
      put_pad('.', max_insert[j] - k, fp);
    }

    fprintf(fp, " %*d ", coord_width, coord2);
    fprintf(fp, "%c ", strand);
    fprintf(fp, "%#.3g\n",
	    marginal_score(&d->scorer, aln, i, &seqs[i]) / xlog(2));
  }

  free(max_insert);
}

/* Print a multilevel consensus sequence, like MEME */
void print_consensus(FILE *fp, const glam2_aln *aln, const data *d) {
  const int width = aln->width;
  const glam2_col *cols = aln->cols;
  const int alph_size = cols->alph_size;
  const int coord_width = digits(d->seqs.maxlen);
  int **counts = xmalloc2(width, alph_size * sizeof(int));
  int *max_insert = xmalloc(width * sizeof(int));
  int i, j;

  for (i = 0; i < width; ++i)
    COPY(counts[i], cols[i].emission_counts, alph_size);

  /* get maximum insert size at each position: */
  for (i = 0; i < width; ++i)
    max_insert[i] = max_right_insert(aln, i);

  for (i = 0; i < alph_size; ++i) {
    int end_flag = 1;
    for (j = 0; j < width; ++j) {
      const int *m = max_int(counts[j], alph_size);
      if (*m * 5 >= cols[j].match_count && *m > 0)
	  end_flag = 0;
    }
    if (end_flag)
      break;

    put_pad(' ', name_width + coord_width + 2, fp);

    for (j = 0; j < width; ++j) {
      int *m = max_int(counts[j], alph_size);
      if (*m * 5 >= cols[j].match_count && *m > 0)
	putc(d->alph.decode[m - counts[j]], fp);
      else
	putc(' ', fp);
      *m = 0;

      put_pad(' ', max_insert[j], fp);
    }

    putc('\n', fp);
  }

  free(max_insert);
  free2(counts, width);
}

/* Print residue counts, indel counts, and score for each column */
void print_col_scores(FILE *fp, glam2_aln *aln, data *d) {
  int i, j;

  for (j = 0; j < d->alph.size; ++j)
    fprintf(fp, "%2c ", d->alph.decode[j]);
  fputs("Del Ins Score\n", fp);  

  for (i = 0; i < aln->width; ++i) {
    if (i != 0) {
      put_pad(' ', 3 * d->alph.size, fp);
      fprintf(fp, "    %3d ", aln->insert_counts[i-1]);
      fprintf(fp, "% #.3g\n", insertion_score(&d->scorer, aln->insert_counts[i-1], aln->aligned_seq) / xlog(2));
    }
    for (j = 0; j < d->alph.size; ++j)
      fprintf(fp, "%2d ", aln->cols[i].emission_counts[j]);
    fprintf(fp, "%3d     ", aln->cols[i].delete_count);
    fprintf(fp, "% #.3g\n", column_score(&d->scorer, &aln->cols[i]) / xlog(2));
  }
}

/* Print LOGO for alignment */
void print_logo(glam2_aln *aln, data *d, int imotif) {
  int i, j;
  char* logodir = d->a.out_dir;
  int num_sites = aln->seq_num;
  int width = aln->width;
  int alph_size = aln->cols[0].alph_size;
  char *alphabet = alph_size == 4 ? DNA0 : PROTEIN0;
  double logo_height = LOGOHEIGHT;
  double logo_width =  width <= MAXLOGOWIDTH ? width : MAXLOGOWIDTH;

  /* convert alignment to motif struct */
  MOTIF_T motif;
  strcpy(motif.id, "0");
  motif.num_sites  = num_sites;
  motif.length     = width;
  motif.alph_size  = alph_size;
  motif.ambigs     = 0;
  motif.evalue     = 0.0;
  motif.complexity = 0.0;
  motif.trim_left = 0;
  motif.trim_right = 0;
  // convert counts to frequencies
  motif.freqs = allocate_matrix(width, alph_size);
  for (i = 0; i < width; ++i) {
    for (j = 0; j < alph_size; ++j) {
      set_matrix_cell(i, j, 
	// should give a different "n" for each row when making logo
        aln->cols[i].emission_counts[j]/((double) aln->cols[i].match_count),
        motif.freqs
      );
    }
  }

  // create the output path
  char *path = NULL;
  Resize(path, strlen(logodir)+29, char);     // room for "/logo_ssc<16digts>\0"

  // create logo without small sample correction
  sprintf(path, "%s/logo%d", logodir, imotif);
  CL_create2(
    &motif,                   // first motif
    "",                       // no title 
    NULL,                     // no second motif
    "",                       // no x-axis label
    FALSE,                    // no error bars
    FALSE,                    // ssc
    logo_height,              // logo height (cm)
    logo_width,               // logo width (cm)
    alphabet,		      // alphabet
    0,                        // no offset to second motif
    path,                     // output file path
    "GLAM2 (no SSC)"          // program name
  );

  // create logo with small sample correction and error bars
  sprintf(path, "%s/logo_ssc%d", logodir, imotif);
  CL_create2(
    &motif,                   // first motif
    "",                       // no title 
    NULL,                     // no second motif
    "",                       // no x-axis label
    TRUE,                     // error bars
    TRUE,                     // ssc
    logo_height,              // logo height (cm)
    logo_width,               // logo width (cm)
    alphabet,		      // alphabet
    0,                        // no offset to second motif
    path,                     // output file path
    "GLAM2 (with SSC)"        // program name
  );

  myfree(path);
} // print_logo

/* Print extended information about an alignment */
void print_aln_info(FILE *fp, glam2_aln *aln, data *d) {
  fprintf(fp, "Score: %#g  Columns: %d  Sequences: %d\n",
          aln->score / xlog(2), aln->width, aln->aligned_seq);
  putc('\n', fp);
  print_aln(fp, aln, d);
  putc('\n', fp);
  print_consensus(fp, aln, d);
  putc('\n', fp);
  print_col_scores(fp, aln, d);
  putc('\n', fp);
}

/* Print a list of alignments */
void print_alns(FILE *fp, glam2_aln *alns, data *d) {
  int i;
  for (i = 0; i < d->a.runs; ++i) {
    print_aln_info(fp, &alns[i], d);
    print_logo(&alns[i], d, i+1);
  }
}
