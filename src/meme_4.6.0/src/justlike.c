/*
 * $Id: justlike.c 4046 2009-09-28 00:47:29Z james_johnson $
 * 
 * $Log$
 * Revision 1.2.2.1  2006/01/13 02:34:43  twhitington
 * Replace hash.(c,h) with new and improved hash_table.(c,h).
 *
 * Revision 1.2  2005/10/25 19:06:39  nadya
 * rm old macro for Header, all info is taken care of by Id and Log.
 *
 * Revision 1.1.1.1  2005/07/29 00:21:21  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/***********************************************************************
*								       *
*	MEME++							       *
*	Copyright 1994, The Regents of the University of California    *
*	Author: Timothy L. Bailey				       *
*								       *
***********************************************************************/
/* 7-13-99 tlb; remove dz */

/* justlike.c */
/*	
	Set the z-matrix to 1 for each motif start; to 0 otherwise.
	Motifs have been read in from a .motif file.
*/

#include "meme.h"

/**********************************************************************/
/*
	like_e_step

*/
/**********************************************************************/
extern double like_e_step(
  MODEL *model,					/* the model */
  DATASET *dataset 				/* the dataset */
)
{
  int i, j;
  int n_samples = dataset->n_samples;
  SAMPLE **samples = dataset->samples;
  int nsites = 0;				/* number of sites read in */
  int imotif = model->imotif;			/* motif number */
  MOTIF motif = dataset->motifs[imotif-1];	/* known motif */
  int w = motif.width;				/* width of known motif */

  /* set all z's to 0 except motif starts */
  for (i=0; i < n_samples; i++) {		/* sequence */
    SAMPLE *s = samples[i];
    double *zi = s->z;                          // zi[j], j in [-lseq...+lseq] 
    char *sample_name = s->sample_name;
    int lseq = s->length;
    for (j=0; j<=lseq-w; j++) {			/* offset */
      if (hash_lookup(sample_name, j+1, motif.ht) != NULL) {
        Zi(j) = 1;				/* motif start */
        nsites++;
      } else {
        Zi(j) = 0;				/* not motif start */
      }
    }
  }

  /* calculate lambda for motif */
  model->lambda = (double) nsites/wps(dataset, w); 

  return 0; 
} /* like_e_step */

