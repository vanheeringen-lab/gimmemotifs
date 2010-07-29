/*
 * $Id: llr.h 1048 2006-07-06 20:07:44Z cegrant $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 18:40:53  nadya
 * Initial revision
 *
 */
 
#ifndef LLR_H
#define LLR_H

#include "meme.h"

extern void init_llr_pv_tables(
  int min,				/* minimum number of sites */
  int max,				/* maximum number of sites */
  int alength,				/* alphabet length */
  double *back,				/* background frequencies */
  BOOLEAN pal				/* sites are palindromes */
);

extern double get_llr_pv(
  double llr,				/* log likelihood ratio */
  double n,				/* number of sequences in alignment */
  int w,				/* width of alignment */
  int range,				/* desired range of scaled LLR */
  double frac,				/* speedup factor */
  int alength,				/* length of alphabet */
  double *dd 				/* alphabet frequency distribution */
);

extern double get_llr_mean(
  double n 				/* number sequences in alignment */
);

#endif
