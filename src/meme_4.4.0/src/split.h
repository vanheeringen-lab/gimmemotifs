/*
 * $Id: split.h 1048 2006-07-06 20:07:44Z cegrant $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 19:11:18  nadya
 * Initial revision
 *
 */

#ifndef SPLIT_H
#  define SPLIT_H 

#include <macros.h>

/* longest line allowd in in file */
#define MAXLINE 10000

extern FILE *split(
  int nmode,		/* 0 - open; 1 - print */
  FILE *infile,
  int n
);

#endif
