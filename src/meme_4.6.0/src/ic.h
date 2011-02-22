/*
 * $Id: ic.h 1048 2006-07-06 20:07:44Z cegrant $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 18:40:26  nadya
 * Initial revision
 *
 */
 
#ifndef IC_H
#define IC_H

double get_wllr_pv(
  double ll0,				/* likelihood under model 0 */
  double ll1,				/* likelihood under model 1 */
  double w0,				/* wllr = ll1 - w0*ll0 */
  int N,				/* number of sequences in alignment */
  int w,				/* width of alignment */
  double alpha,				/* scale factor for wllr */
  double frac,				/* speedup factor */
  int alength,				/* length of alphabet */
  double *dd 				/* alphabet frequency distribution */
);

#endif
