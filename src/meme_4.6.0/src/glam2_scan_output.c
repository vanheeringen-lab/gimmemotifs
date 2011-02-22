#include "glam2_util.h"
#include "glam2_scan_output.h"

void print_hit(FILE *fp, const alignment *aln, const alphabet *alph) {
  //const int name_width = 12;
  const int name_width = 24;		// TLB 
  const int coord1_width = digits(aln->coord1);
  int i;

  put_pad(' ', name_width + coord1_width + 2, fp);

  for (i = 0; i < aln->aln_size; ++i)
    if (aln->seq1[i] == 1)
      putc('*', fp);
    else
      putc('.', fp);

  putc('\n', fp);

  fprintf(fp, "%-*.*s", name_width, name_width, aln->name);
  fprintf(fp, " %d ", aln->coord1);

  for (i = 0; i < aln->aln_size; ++i)
    if (aln->seq2[i] < alph->size || aln->seq1[i] == 0)
      putc(alph->decode[aln->seq2[i]], fp);
    else
      putc('.', fp);

  fprintf(fp, " %d ", aln->coord2);
  fprintf(fp, "%c ", aln->strand);
  fprintf(fp, "%#.3g\n", aln->score / xlog(2));
}

void print_hits(FILE *fp, const alignment *alns, const data *d) {
  int i;
  for (i = 0; i < d->hit_num; ++i) {
    print_hit(fp, &alns[i], &d->alph);
    putc('\n', fp);
  }
}
