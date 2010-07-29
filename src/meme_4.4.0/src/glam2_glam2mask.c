/* Mask a glam2 motif in FASTA-format sequences */
#include <stdio.h>
#include <stdlib.h>  /* free */
#include <string.h>  /* strcmp */
#include <unistd.h>  /* non-ANSI */
#include "glam2_util.h"
#include "glam2_fasta.h"
#include "glam2_alignment.h"

void mask_seqs(const alignment *a, int mask_char, FILE *in, FILE *out) {
  static const int name_width = 12;  /* max name width for glam2 */
  static const int line_size = 60;  /* line size for FASTA-format output */
  fasta f;
  size_t i = 0;

  while (fasta_read(&f, in) != EOF) {
    const aligned_seq *s = &a->seqs[i];
    /* mangle the title in the same way as glam2: */
    char *title = xstrdup(f.title);
    first_word(title);
    strtrunc(title, name_width);

    if (i != a->seq_num && strcmp(title, s->name) == 0) {
      int pos = xatoi(s->start)-1;
      size_t j;

      for (j = 0; s->seq[j] != 0; ++j) {
	if (pos >= f.seqlen)
	  break;  /* die with error message? */
	if (s->seq[j] != '.') {
	  if (a->key_positions[j] == '*')
	    f.seq[pos] = mask_char;
	  ++pos;
	}
      }

      ++i;
    }

    free(title);
    put_fasta(&f, line_size, out);
    free_fasta(&f);
  }

  if (i != a->seq_num)
    die("%s: didn't find %s among the sequences\n",
	prog_name, a->seqs[i].name);
}

static void usage(void) {
  die("\
Usage: glam2mask [options] my_motif.glam2 my_seqs.fa\n\
Options (default settings):\n\
-o: output file (stdout)\n\
-x: mask character (x)\n\
");
}

int main(int argc, char **argv) {
  const char *mot_file;
  const char *seq_file;
  const char *out_file = "-";
  int mask_char = 'x';
  FILE *in;
  FILE *out;
  alignment aln;
  int c;

  prog_name = "glam2mask";  /* for error messages */

  while ((c = getopt(argc, argv, "o:x:")) != -1) {  /* non-ANSI */
    switch (c) {
    case 'o':
      out_file = optarg;
      break;
    case 'x':
      mask_char = *optarg;
      break;
    case '?':
      usage();
    }
  }

  if (optind != argc-2)
    usage();
  mot_file = argv[optind++];
  seq_file = argv[optind++];

  in = xfopen(mot_file, "r");
  aln_read(&aln, in);
  xfclose(in);

  if (aln.key_positions == NULL)
    die("%s: no motif found in %s\n", prog_name, mot_file);
  if (!aln_same_lengths(&aln))
    die("%s: unequal aligned lengths in %s\n", prog_name, mot_file);

  in = xfopen(seq_file, "r");
  out = xfopen(out_file, "w");
  mask_seqs(&aln, mask_char, in, out);
  xfclose(in);
  xfclose(out);

  aln_free(&aln);
  return 0;
}
