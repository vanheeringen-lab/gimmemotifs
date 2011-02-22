/*
 * $Id: regress.h 1048 2006-07-06 20:07:44Z cegrant $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 19:09:59  nadya
 * Initial revision
 *
 */

/***********************************************************************
*                                                                      *
*       MEME++                                                         *
*       Copyright 1995, The Regents of the University of California    *
*       Author: Timothy L. Bailey                                      *
*                                                                      *
***********************************************************************/
 
#ifndef REGRESS_H 
#define REGRESS_H

extern double regress(
  int n,                        /* number of points */
  double *x,                    /* x values */
  double *y,                    /* y values */
  double *m,                    /* slope */
  double *b                     /* y intercept */
);

#endif
 
