/* structs for GLAM2 alignments */
#ifndef GLAM2_GLAM2_ALN_H
#define GLAM2_GLAM2_ALN_H

typedef struct {  /* one column in a glam2 alignment */
  /* fudamental data: */
  int seq_num;  /* number of sequences */
  int *positions;  /* zero-based coordinates */
  int *matches;  /* 0=deletion, 1=match */
  /* derived data: */
  int alph_size;
  int *emission_counts;
  int match_count;
  int delete_count;
} glam2_col;

typedef struct {  /* a glam2 alignment */
  int seq_num;  /* number of sequences */
  int width;  /* number of columns */
  glam2_col *cols;
  int *strands;  /* '+', '-', or 0 = unaligned */
  /* derived data: */
  int aligned_seq;  /* number of strands != 0 */
  int *insert_counts;
  double score;
} glam2_aln;

/* Is the s-th sequence of alignment "aln" aligned? */
#define ALIGNED(aln, s) ((aln)->strands[s])

#define LEFT_POS(col, seq) ((col)->positions[seq])

#define RIGHT_POS(col, seq) ((col)->positions[seq] + (col)->matches[seq])

#define INSERTS(col1, col2, seq) (LEFT_POS(col2, seq) - RIGHT_POS(col1, seq))

#define LEFT_INSERT(aln, col, seq)\
        ((col) > 0 ?\
         INSERTS(&(aln)->cols[(col)-1], &(aln)->cols[col], seq) : 0)

#define RIGHT_INSERT(aln, col, seq)\
        ((col) < (aln)->width-1 ?\
         INSERTS(&(aln)->cols[col], &(aln)->cols[(col)+1], seq) : 0)

/* Allocate memory for a glam2_col */
void col_init(glam2_col *col, int seq_num, int alph_size);

/* Allocate memory for a glam2_aln */
void aln_init(glam2_aln *aln, int seq_num, int max_width, int alph_size);

/* Set all sequences unaligned and all counts to zero */
/* Assumes width has been set */
void aln_zero(glam2_aln *aln);

/* Copy column b into column a */
/* Assumes a's memory has already been allocated */
void col_copy(glam2_col *a, const glam2_col *b);

/* Copy alignment b into alignment a */
/* Assumes a's memory has already been allocated */
void aln_copy(glam2_aln *a, const glam2_aln *b);

#endif
