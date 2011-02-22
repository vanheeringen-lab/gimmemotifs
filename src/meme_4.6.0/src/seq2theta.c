/*
 * $Id: seq2theta.c 2220 2007-11-30 19:09:15Z cegrant $
 * 
 * $Log$
 * Revision 1.3  2006/03/08 20:50:11  nadya
 * merge chamges from v3_5_2 branch
 *
 * Revision 1.2.4.1  2006/01/24 20:44:08  nadya
 * update copyright
 *
 * Revision 1.2  2005/10/25 19:06:39  nadya
 * rm old macro for Header, all info is taken care of by Id and Log.
 *
 * Revision 1.1.1.1  2005/07/29 17:26:42  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/*#define DEBUG*/
/***********************************************************************
*                                                                      *
*       MEME++                                                         *
*       Copyright 1994-2006, The Regents of the University of California *
*       Author: Timothy L. Bailey                                      *
*                                                                      *
***********************************************************************/

/*
	Routines to map a sequence to a theta matrix.
*/

#include "meme.h"

/**********************************************************************/
/*
       	init_map

	Set up an |A| x |A| sequence_letter_to_frequency_vector
        or sequence_letter_to_logodds_vector
	matrix where each column is the 
	mapping from a given letter to a frequency/logodds vector.
	An extra row and column is added for the wildcard character 'X', and
	set to the background frequencies (or 0 if, logodds).

	Two types of mapping are possible:
		Uni -	add n prior
		Pam -	Dayhoff PAM for proteins, transversion/transition
			PAM for DNA

	Returns the map.

*/
/**********************************************************************/

extern THETA init_map(
  MAP_TYPE type,		/* type of mapping:
					Uni	- add n prior
					Pam	- pam matrix
				*/
  double scale,			/* degree of crispness, depends on type,
					Uni	- pseudo count n to add
					Pam	- pam distance
				*/
  int alength,			/* length of alphabet */
  double *back,			/* background frequencies for 'X' */
  BOOLEAN lo			/* setup a logodds mapping if true */
)
{
  int i, j, p;
  THETA map;			/* the map */

  /* dayhoff PAM 1 matrix; order of alphabet: ACDEFGHIKLMNPQRSTVWY */
  /* dayhoff_ij = Pr(amino acid j --> amino acid i | time=1) */
  double dayhoff[20][20] = {
    { 9867, 3, 10, 17, 2, 21, 2, 6, 2, 4, 6, 9, 22, 8, 2, 35, 32, 18, 0, 2},
    { 1, 9973, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 5, 1, 2, 0, 3},
    { 6, 0, 9859, 53, 0, 6, 4, 1, 3, 0, 0, 42, 1, 6, 0, 5, 3, 1, 0, 0},
    { 10, 0, 56, 9865, 0, 4, 2, 3, 4, 1, 1, 7, 3, 35, 0, 4, 2, 2, 0, 1},
    { 1, 0, 0, 0, 9946, 1, 2, 8, 0, 6, 4, 1, 0, 0, 1, 2, 1, 0, 3, 28},
    { 21, 1, 11, 7, 1, 9935, 1, 0, 2, 1, 1, 12, 3, 3, 1, 21, 3, 5, 0, 0},
    { 1, 1, 3, 1, 2, 0, 9912, 0, 1, 1, 0, 18, 3, 20, 8, 1, 1, 1, 1, 4},
    { 2, 2, 1, 2, 7, 0, 0, 9872, 2, 9, 12, 3, 0, 1, 2, 1, 7, 33, 0, 1},
    { 2, 0, 6, 7, 0, 2, 2, 4, 9926, 1, 20, 25, 3, 12, 37, 8, 11, 1, 0, 1},
    { 3, 0, 0, 1, 13, 1, 4, 22, 2, 9947, 45, 3, 3, 6, 1, 1, 3, 15, 4, 2},
    { 1, 0, 0, 0, 1, 0, 0, 5, 4, 8, 9874, 0, 0, 2, 1, 1, 2, 4, 0, 0},
    { 4, 0, 36, 6, 1, 6, 21, 3, 13, 1, 0, 9822, 2, 4, 1, 20, 9, 1, 1, 4},
    { 13, 1, 1, 3, 1, 2, 5, 1, 2, 2, 1, 2, 9926, 8, 5, 12, 4, 2, 0, 0},
    { 3, 0, 5, 27, 0, 1, 23, 1, 6, 3, 4, 4, 6, 9876, 9, 2, 2, 1, 0, 0},
    { 1, 1, 0, 0, 1, 0, 10, 3, 19, 1, 4, 1, 4, 10, 9913, 6, 1, 1, 8, 0},
    { 28, 11, 7, 6, 3, 16, 2, 2, 7, 1, 4, 34, 17, 4, 11, 9840, 38, 2, 5, 2},
    { 22, 1, 4, 2, 1, 2, 1, 11, 8, 2, 6, 13, 5, 3, 2, 32, 9871, 9, 0, 2},
    { 13, 3, 1, 2, 1, 3, 3, 57, 1, 11, 17, 1, 3, 2, 2, 2, 10, 9901, 0, 2},
    { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 9976, 1},
    { 1, 3, 0, 1, 21, 0, 4, 1, 0, 1, 0, 3, 0, 0, 0, 1, 1, 1, 2, 9945}
  };

  /* transition/transversion PAM 1 matrix; alphabet order "ACGT" */
  double trans[4][4] = {
    { 9900, 20, 60, 20}, 
    { 20, 9900, 20, 60},
    { 60, 20, 9900, 20},
    { 20, 60, 20, 9900}
  };

  /* special PAM mutation frequencies for proteins */
  double pam_dna_freq[] = { 
      0.25 /* A */,
      0.25 /* C */,
      0.25 /* G */,
      0.25 /* T */
  };

  /* special PAM mutation frequencies for proteins */
  double pam_prot_freq[] = { 
      0.096 /* A */,
      0.025 /* C */,
      0.053 /* D */,
      0.053 /* E */,
      0.045 /* F */,
      0.090 /* G */,
      0.034 /* H */,
      0.035 /* I */,
      0.085 /* K */,
      0.085 /* L */,
      0.012 /* M */,
      0.042 /* N */,
      0.041 /* P */,
      0.032 /* Q */,
      0.034 /* R */,
      0.057 /* S */,
      0.062 /* T */,
      0.078 /* V */,
      0.012 /* W */,
      0.030 /* Y */,
      0.0   /* X */
  };

  /* create the array for the sequence to theta map */
  create_2array(map, double, alength+1, alength+1);

  switch (type) {
    case Uni: {
      double main_letter;	/* probability of main letter */
      double other;		/* probability of each other letter */

      main_letter = (1.0 + scale)/(1.0 + alength * scale);
      other = scale/(1.0 + alength * scale);
      if (VERBOSE) {printf("main= %g\n\n", main_letter);}
      /* create a matrix of columns; each column gives mapping for a letter */
      for (i=0; i<alength; i++) 
	for (j=0; j<alength; j++)  
	  map[i][j] = (i == j) ? main_letter : other;
    } break;
    case Pam: {
      double mul[20][20]; 

      /* convert initial matrix to probabilities */ 
      for (i=0; i<alength; i++)
	for (j=0; j<alength; j++)
          map[i][j] = ((alength == 4) ? trans[i][j] : dayhoff[i][j]) / 10000;

      /* take pam matrix to desired power */ 

      /* copy: */
      for (i = 0; i < alength; i++)
	for (j = 0; j < alength; j++)
	  mul[i][j] = map[i][j];

      /* multiply: */
      while (--scale) {
        double result[20][20], sum;
	for (i = 0; i < alength; i++) {
	  for (j = 0; j < alength; j++) {
            for (sum = p = 0; p < alength; p++) {
              sum += mul[i][p] * map[p][j];
            }
            result[i][j] = sum;
          }
        }

        for (j = 0; j < alength; j++) {
          for (i = 0; i < alength; i++) {
            /*map[i][j] = result[i][j];*/
            RND(result[i][j], 8, map[i][j]);
          }
        }
      }
    }
  }

  /* add last row and column for mapping from the wildcard 'X' */
  for (i=0; i<alength; i++) map[alength][i] = map[i][alength] = back[i];

  /* convert to logodds matrix if requested */
  if (lo) {
    double *pfreq = alength==4 ? pam_dna_freq : pam_prot_freq;
    double avg, x_avg;

    for (i=0; i<alength; i++) {
      for (j=0; j<=i; j++) {
        map[i][j] = map[j][i] = NINT(2*LOG2(map[i][j]/pfreq[i]));
      }
    }

    /* do the last row and column for "X" matches: average of match to all chars */
    x_avg = 0;						/* average for X row */
    for (i=0; i<alength; i++) {
      avg = 0;						/* average for row */
      for (j=0; j<alength; j++) {
        avg += map[i][j];
      }
      x_avg += avg/alength;
      map[i][alength] = map[alength][i] = NINT(avg/alength);
    }
    map[alength][alength] = NINT(x_avg/alength);
#ifdef DEBUG
    for (i=0; i<=alength; i++) {
      for (j=0; j<=alength; j++) {
        printf("%3g ", map[i][j]);
      }
      printf("\n");
    }
#endif
  } /* lo */

  return map;
} /* init_map */


/**
 * convert_to_lmap
 *
 * Converts a matrix of sequence to theta probability mappings into a
 * matrix containing the int logs of those probabilities. Also sets up
 * a vector for mapping an "X" character to a vector of letter probabilities.
 * Those probabilities are uniform across all letters.
 */
extern void convert_to_lmap (
  THETA map,
  int lmap[MAXALPH][MAXALPH],
  int alength
)
{  
  /* 
    Set up the matrix of frequency vectors for each letter in the alphabet.
    Column and row for the "match-anything" character X are set to 1.0/alength 
    so that such matches are neither favored nor disfavored, and where 
    they match is irrelevant:
  */
  int i,j;
  for (i=0; i<alength+1; i++) {
    for (j=0; j<alength; j++) {
      lmap[i][j] = (i<alength) ? INT_LOG(map[j][i]) : INT_LOG(1.0/alength);
    }
    lmap[i][j] = INT_LOG(1.0/alength);			/* X */
  }
}


/**********************************************************************/
/*
	  init_theta

	  Set theta to represent a consensus sequence by copying
	  columns of the letter_to_theta_column_array map to theta.

	  For columns in sequence containing 'X' or if the sequence
	  is null, theta is set to the values in the extra (alength)
	  column in the map.  Each column is length alength+1, because
	  there is a row for 'X'.
*/
/**********************************************************************/
extern void init_theta(
  THETA theta,			/* theta */
  char *start,			/* integer encoded starting sequence */
  int w,			/* width of motif */
  THETA map,			/* frequency vectors for each letter */ 
  int alength			/* alphabet length */
)
{
  int i, j, c;

  /* initialize the insite frequencies */
  for (i=0; i < w; i++) {		/* column in theta */
    for (j=0; j <= alength; j++) {	/* row in theta */
      c = start ? start[i] : alength;	/* alength is index of 'X' */
      theta(i, j) = map[c][j];
    }
  }
} /* init_theta */


/**
 * convert_to_ltheta
 *
 * Convert the entries in the specified matrix motif model from doubles to INT
 * LOG values.
 */
extern void convert_to_ltheta (
  double matrix_ds[MAXSITE][MAXALPH], ///< The input matrix of doubles
  int matrix_il[MAXSITE][MAXALPH], ///< The output matrix of int logs
  int nrows,
  int ncols
)
{
  int row_idx;
  for (row_idx = 0; row_idx < nrows; row_idx++) {
    int col_idx;
    for (col_idx = 0; col_idx < ncols; col_idx++) {
      matrix_il[row_idx][col_idx] = INT_LOG(matrix_ds[row_idx][col_idx]);
      fprintf(stderr, "%i ", matrix_il[row_idx][col_idx]);
    }
    fprintf(stderr, "\n");
  }
}

