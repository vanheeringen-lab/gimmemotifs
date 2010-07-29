/*
 * $Id: display_globals.h$
 *
 * $Log$
 *
 */


#ifndef DISPLAY_GLOBALS_H
#define DISPLAY_GLOBALS_H

// Shared by display.c and meme.c
EXTERN LO *los[MAXG];     /* logodds structure for each motif */
EXTERN double *pv[MAXG];  /* p-value tables for each motif */

#endif
