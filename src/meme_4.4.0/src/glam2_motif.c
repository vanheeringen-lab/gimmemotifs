#include <assert.h>
#include <limits.h>  /* INT_MAX */
#include "glam2_util.h"
#include "glam2_alignment.h"
#include "glam2_motif.h"

/* Move to util? */
static size_t strcnt(const char *cs, int c) {
  size_t count = 0;
  assert(cs != NULL);
  for (; *cs != '\0'; ++cs)
    count += *cs == (char)c;  /* same cast as strchr(?) */
  return count;
}

void aln2mot(motif *m, const alignment *a, int alph_size, const int *encode) {
  const int columns = strlen(a->key_positions);
  const int width = strcnt(a->key_positions, '*');
  int **residue_counts = xcalloc2(width, alph_size, sizeof(int));  /* zero */
  int *delete_counts = xcalloc(width, sizeof(int));  /* zero fill */
  int *insert_counts = xcalloc(width, sizeof(int));  /* zero fill */
  int i, j, k;

  for (i = 0; i != a->seq_num; ++i) {
    const char *seq = a->seqs[i].seq;
    k = 0;

    for (j = 0; j < columns; ++j) {
      const int c = (unsigned char)seq[j];  /* is this OK? */
      if (a->key_positions[j] == '*') {
	if (c != '.') {
	  if (encode[c] == alph_size)
	    die("%s: error reading motif file: ambiguous residue %c in aligned column\n", prog_name, c);
	  ++residue_counts[k][encode[c]];
	}
	delete_counts[k] += c == '.';
	++k;
      } else {
	assert(k > 0);
	assert(k < width);
	insert_counts[k-1] += c != '.';
      }
    }
  }  

  m->width = width;
  m->alph_size = alph_size;
  m->seq_num = a->seq_num;
  m->residue_counts = residue_counts;
  m->delete_counts = delete_counts;
  m->insert_counts = insert_counts;
}

void read_motif(motif *m, int alph_size, const int *encode, FILE *fp) {
  alignment aln;
  aln_read(&aln, fp);

  if (aln.seq_num == 0)
    die("%s: error reading motif file: no motif\n", prog_name);
  if (!aln_same_lengths(&aln))
    die("%s: error reading motif file: unequal aligned lengths\n", prog_name);
  assert(strlen(aln.key_positions) <= INT_MAX);  /* can fail */

  aln2mot(m, &aln, alph_size, encode);
  aln_free(&aln);
}
