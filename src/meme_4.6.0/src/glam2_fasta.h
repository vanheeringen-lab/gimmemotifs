/* Structs and functions for FASTA-format sequences */
#ifndef GLAM2_FASTA_H
#define GLAM2_FASTA_H

#include <stdio.h>

/* One FASTA sequence and its title */
typedef struct {
  char *title;  /* FASTA title line (without newline; NUL terminated) */
  size_t titlen;  /* length of title line */
  int *seq;  /* should probably be char or unsigned char */
  int seqlen;  /* maybe use unsigned or size_t? */
  int *rcseq;  /* reverse-complemented sequence */
} fasta;

/* Multiple FASTA sequences */
typedef struct {
  fasta *f;  /* the fastas */
  int seqnum;  /* the number of fastas */
  int maxlen;  /* maximum sequence length */
} mfasta;

/* Replace a string with its first word */
/* Often seems to be used with fasta titles */
char *first_word(char *s);

/* For now, these functions succeed or die, because error recovery is harder */

/* Read a FASTA-format sequence and title line from a file */
int fasta_read(fasta *f, FILE *fp);

/* Read multiple FASTA-format sequences from a file */
void mfasta_read(mfasta *m, FILE *fp);

/* Write a FASTA-format title & sequence to a file */
void put_fasta(const fasta *f, int line_size, FILE *fp);

/* Deallocate the memory buffers within the fasta object */
void free_fasta(fasta *f);

/* Deallocate all the memory within the mfasta */
void free_mfasta(mfasta *m);

/* Translate (encode) the sequence */
void tr_fasta(fasta *f, const int *encode);

/* Translate (encode) each sequence */
void tr_mfasta(mfasta *m, const int *encode);

/* Initialize rcseq */
void rc_fasta(fasta *f, const int alph_size);

/* Initialize rcseqs */
void rc_mfasta(mfasta *m, const int alph_size);

/* Count residue types */
/* Assumes counts has been suitably allocated and initialized */
void count_fasta(const fasta *f, int *counts);

/* Count residue types in all sequences */
/* Assumes counts has been suitably allocated and initialized */
void count_mfasta(const mfasta *m, int *counts);

#endif
