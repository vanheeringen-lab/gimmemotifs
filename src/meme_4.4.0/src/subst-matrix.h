/*****************************************************************************
 * FILE: subst_matrix.h
 * AUTHOR: Timothy Bailey
 * CREATE DATE: 4/23/2002
 * PROJECT: shared
 * COPYRIGHT: 2002, University of Queensland
 * VERSION: $Revision: 1.1.1.1 $
 * DESCRIPTION: Substitution matrices: scoring matrices and target frequencies
 *****************************************************************************/
#ifndef SUBST_MATRIX_H
#define SUBST_MATRIX_H

MATRIX_T *get_subst_target_matrix(
  char *score_filename,		/* name of score file */
  BOOLEAN_T dna,		/* DNA alphabet if true; default: protein */
  int dist,			/* PAM distance (ignored if score_filename != NULL) */
  char *alpha,			/* standard alphabet for final matrix */
  ARRAY_T *back			/* background frequencies of standard alphabet */
);

ARRAY_T *get_pseudocount_freqs(
   ARRAY_T *      f,            /* Foreground distribution. */
   ARRAY_T *      b,            /* Background distribution. */
   MATRIX_T *     target_freq   /* Target frequency matrix. */
);

#endif
