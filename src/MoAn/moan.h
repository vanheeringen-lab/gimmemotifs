/* MoAn is a motif discovery tool for DNA sequences */
/* Copyright (C) 2006 Eivind Valen */

/* This program is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU General Public License */
/* as published by the Free Software Foundation; either version 2 */
/* of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program; if not, write to the Free Software */
/* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. */

#ifndef _MOAN_H_
#define _MOAN_H_

#include <stdlib.h>
#include <stdio.h>

#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

#define VERSION "1.01"
#define TOTAL 100

/* Parameters*/
#define MAX_LENGTH 15
#define MIN_LENGTH 5
#define END_TEMP 0.1
#define NOOCCUR -10     /* log(0) */
#define MAX_SEQ 100000
#define ALPHLEN      4
#define STARTEMP_ITERATIONS 1000
/* Loglevels */
#define QUIET        0
#define NORMAL       1
#define VERBOSE      2
#define DBG          3
#define TRACE        4

/* The maximum and minimum cutoff for DYNAMIC */
/* #define CUT_MIN   0.8 */
#define CUT_MIN   0.8
#define CUT_MAX   2.0   
#define CUT_INT   (CUT_MAX - CUT_MIN)


/* Annealing schedule */
#define STATIC_TEMP
#define DYN_LIMIT 0.8

#ifdef  STATIC_TEMP
#define EXP_FACTOR search->step_ratio
#define LIN_FACTOR search->step_ratio
#endif

#ifdef  DYNAMIC_TEMP
#define EXP_FACTOR ((float)((float)search->accrej)/((float)(search->maxit - i)) > DYN_LIMIT) ? search->step_ratio : (2.0 - search->step_ratio)
#define LIN_FACTOR ((float)((float)search->accrej)/((float)(search->maxit - i)) > DYN_LIMIT) ? search->step_ratio : -search->step_ratio
#endif  

#define OUTPUT(lvl, text) if (search->loglevel >= lvl){ printf(text); }



/* RAND */
#define RAND_INT(max)   (int)   (((float) max) * rand() / (RAND_MAX + 1.0))
#define RAND_FLOAT(max) (float) (((float) max) * rand() / (RAND_MAX + 1.0))

#define MAX(x,y) ((x > y) ? x : y)
#define MIN(x,y) ((x > y) ? y : x)

#define LOG(x) (log(x)/log(2))
#define EXP(x) pow(2, x)
#define LOG2(x) (log(x)/log(2))
/* #define LOG(x) log(x) */
/* #define EXP(x) exp(x) */


typedef enum {FALSE, TRUE} bool;
typedef struct PSSM WM;
enum Alph {A,C,G,T};

typedef struct Sequences {
  char **sequence;                    /* The sequences */
  unsigned int pos_seq;               /* Positive sequences */
  unsigned int neg_seq;               /* Negative sequences */
  /* Length data */
  unsigned int *seq_len;              /* Length of each sequence */
  float *len_rec;                     /* Reciproc of the length of each sequence */
  unsigned int pos_len;               /* Combined length of positive sequences*/
  unsigned int neg_len;               /* Combined length of negative sequences*/
  /* Distribution */
  unsigned int freq[ALPHLEN];
  /* Files */
  FILE *pos_file;                     /* File containing positive set */
  FILE *neg_file;                     /* File containing negative set */
} Sequences;


inline extern int chr2ind(char i) {
  switch (i) {
  case 'A': return A;
  case 'C': return C;
  case 'G': return G;
  case 'T': return T;
  case 'a': return A;
  case 'c': return C;
  case 'g': return G;
  case 't': return T;
  }

  return -1;
}

inline extern char ind2chr(int i) {
  switch (i) {
  case 0: return 'A';
  case 1: return 'C';
  case 2: return 'G';
  case 3: return 'T';
  }

  return -1;
}



#endif

