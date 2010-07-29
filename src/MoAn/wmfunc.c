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

#include "pssm.h"
#include "wmfunc.h"

PSSM load_count_matrices(char *filename, int alphlen) {
  PSSM pssm;
  FILE *file;
  int i, j, k, order, length;

  file = fopen(filename, "r");
  fscanf(file, "# %d\n", &i);
  fscanf(file, "> %d %d\n", &order, &length);
/*   fprintf(stderr, "ORDER(%i) LEN(%i)\n", order, length); */
  
  pssm = initMatrix(order, length, alphlen);

  for (j = 0; j < alphlen; j++) {
    for (k = 0; k < length; k++) {
      fscanf(file, "%d", &pssm->counts[j * length + k]);
      pssm->scores[j * length + k] = (pssm->counts[j * length + k] > 0 ? log(((double) pssm->counts[j * length + k])) : LOG_ZERO);
/*       fprintf(stderr, "%5.2f\t", pssm->scores[j * length + k]); */
    }
/*     fprintf(stderr, "\n"); */
  }
  fclose(file);

  return pssm;
}

PSSM load_log_matrices(char *filename, int alphlen) {
  PSSM pssm;
  FILE *file;
  int i, j, k, order, length;

  file = fopen(filename, "r");
  fscanf(file, "# %d\n", &i);
  fscanf(file, "> %d %d\n", &order, &length);
/*   fprintf(stderr, "ORDER(%i) LEN(%i)\n", order, length); */
  
  pssm = initMatrix(order, length, alphlen);

  for (j = 0; j < alphlen; j++) {
    for (k = 0; k < length; k++) {
      fscanf(file, "%lf", &pssm->scores[k * alphlen + j]);
/*       fprintf(stderr, "%5.2f\t", pssm->scores[j * length + k]); */
    }
/*     fprintf(stderr, "\n"); */
  }
  fclose(file);

  return pssm;
}

void print_pssm(PSSM pssm) {
  int ltr, pos, ord;
  int alphlen = pssm->alphabetSize;

  for (ltr = 0; ltr < pssm->alphabetSize; ltr++) {
    ord = pssm->order;

    for (pos = 0; pos < pssm->length; pos++) {      
      printf("%7.2f", pssm->scores[pos * alphlen + ltr]);
    }
    printf("\n");
  }
}




/* 
 * Normalizes a column in a count matrix
 * pssm - The PSSM to normalize
 * pos  - The column number
 * sum  - The total sum for the column
 *
 * FIX: Only zeroth-order
 */
/* inline void normalize_col(PSSM pssm, unsigned int pos, unsigned int sum) { */
/*   unsigned int i, ind, total = 0; */
 
/*   i = pssm->alphabetSize; */
/*   while (i--) total += getScoreFast(pssm, i, pos) */

/*   if (total != 0) { */
/*     i = pssm->alphabetSize; */
/*     while (i--) { */
/*       ind = (int) (((float) (pssm->counts[ind] * sum)) / (float) total); */
/*       setCountFast(pssm, i, pos, ind); */
/*       setScoreFast(pssm, i, pos, logtable[ind]); */
/*     } */

/*   } else { */
/*     i = pssm->alphabetSize; */
/*     while (i--) { */
/*       ind = sum / pssm->alphabetSize; */
/*       setCountFast(pssm, i, pos, ind); */
/*       setScoreFast(pssm, i, pos, logtable[ind]); */
/*     } */

/*   } */
  
/* } */


/* 
 * Normalizes a count matrix
 * pssm - The PSSM to normalize
 * sum  - The total sum for each column
 *
 * FIX: Only zeroth-order
 */
/* inline void normalize_all(PSSM pssm, unsigned int sum) { */
/*   unsigned int pos;  */

/*   pos = pssm->max_length; */
/*   while (pos--) normalizeColumnCounts(pssm, pos, sum);   */
/* } */

/* 
 * Function for printing a PSSM
 * pssm - The PSSM to print
 *
 * FIX: Only zeroth-order
 */
