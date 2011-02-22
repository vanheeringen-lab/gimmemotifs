/*
 * $Id: message.h 1048 2006-07-06 20:07:44Z cegrant $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 18:44:39  nadya
 * Initial revision
 *
 */

/***********************************************************************
*								       *
*	MEME++							       *
*	Copyright 1996, The Regents of the University of California    *
*	Author: Bill Grundy					       *
*								       *
***********************************************************************/

#ifndef MESSAGE_H
#  define MESSAGE_H

#include "meme.h" /* We need the defn of S_POINT. */
/***************************************************************************/
/* Type declarations */
/***************************************************************************/
/* To do a reduction, we first pass small packets around containing
   the discriminating information. */
typedef struct reduction_packet {
  double data;				/* logev */
  double s_width;			/* starting width */
  double s_nsites;			/* starting nsites */
  double ID; /* Use a double so the MPI type handle is simple. */
} REDUCE_PACKET;

typedef struct s_point_packet {
  double score;
  double iseq; /* Use a double so the MPI type handle is simple. */
  double ioff; /* Use a double so the MPI type handle is simple. */
} S_POINT_PACKET;

/* This structure stores each processor's starting and ending points
   in the search for starting points. */
typedef struct seq_params {
  int start_seq;
  int start_off;
  int end_seq;
  int end_off;
} SEQ_PARAMS;

/***************************************************************************/
/* Global variables */
/***************************************************************************/

EXTERN SEQ_PARAMS *start_n_end;

/***************************************************************************/

/* Macros  */

/***************************************************************************/
/*
	llr_node(nsites, pal)

	True if this node does llr pvalue table for nsites
*/
/***************************************************************************/
#define llr_node(nsites, pal) ((nsites/((pal)?2:1)) % mpNodes())

/***************************************************************************/
/*
	balance_loop
*/
/***************************************************************************/
#ifdef PARALLEL
#define balance_loop(s, n) balance_loop1(s, n)
#else
#define balance_loop(s, n) 
#endif

/***************************************************************************/
/*
	load_balance_llr(nsites)

	Load balancing function for the llr pvalue tables.
	Returns true if the node computes table nsites, false otherwise.
*/
/***************************************************************************/

#ifdef PARALLEL
#define load_balance_llr(nsites, pal) (llr_node((nsites),(pal)) == mpMyID())
#else
#define load_balance_llr(nsites, pal) TRUE
#endif

/***************************************************************************/
/*
	Broadcast function for llr pvalue tables.
*/
/***************************************************************************/
#ifdef PARALLEL

#define broadcast_llr(min, max, pal) {					\
  int i;								\
  double hdr[5];				/* holds header */	\
  ndistrs = max;				/* set the global */	\
  Resize(distrs, ndistrs+1, DISTR);		/* set of tables */	\
  for (i=min; i<=max; i += (pal)?2:1) {		/* bcast from node */	\
    if (load_balance_llr(i, (pal))) {		/* pack header */	\
      hdr[0] = distrs[i].alpha; 					\
      hdr[1] = distrs[i].w ? distrs[i].range[1] : 0; 			\
      hdr[2] = distrs[i].offset[1]; 					\
      hdr[3] = distrs[i].mean;						\
      hdr[4] = distrs[i].w;						\
    } 									\
    /* broadcast header */						\
    mpBroadcast((void *)hdr, 5*sizeof(double), llr_node(i, (pal)));	\
    if (!load_balance_llr(i, (pal))) {		/* unpack header */	\
      distrs[i].range = NULL;						\
      distrs[i].offset = NULL;						\
      distrs[i].cdf = NULL;						\
      Resize(distrs[i].range, 2, int);					\
      Resize(distrs[i].offset, 2, int);					\
      Resize(distrs[i].cdf, 2, double *);				\
      distrs[i].alpha = hdr[0]; 					\
      distrs[i].range[1] = hdr[1]; 					\
      distrs[i].offset[1] = hdr[2]; 					\
      distrs[i].mean = hdr[3];						\
      distrs[i].w = hdr[4];						\
      distrs[i].cdf[1] = NULL;						\
      Resize(distrs[i].cdf[1], distrs[i].range[1]+1, double);		\
    } 									\
    /* broadcast the p-value table */					\
    if (distrs[i].w) 							\
      mpBroadcast((void *)distrs[i].cdf[1],				\
        sizeof(double)*(distrs[i].range[1]+1), llr_node(i, (pal))); 	\
  }									\
}

#else
#define broadcast_llr(min, max, pal)
#endif

/***************************************************************************/
/* Function prototypes */
/***************************************************************************/

extern void store_consensus(
  MODEL *model,
  CANDIDATE *candidates
);

extern void print_model(
   char *label,
   MODEL *model
);

extern void reduce_across_models(
   MODEL *model,
   int alength
);

extern void balance_loop1(
   SAMPLE **samples,
   int n_samples
);

extern void reduce_across_s_points(
  S_POINT *s_points,
  SAMPLE **samples,
  int n_nsites0,
  int n_starts
);			  

extern void get_start_n_end (
  int *start_seq,			/* starting sequence number */
  int *start_off,			/* offset in starting sequence */
  int *end_seq,				/* ending sequence number */
  int *end_off 				/* offset in ending sequence */
);

#endif
