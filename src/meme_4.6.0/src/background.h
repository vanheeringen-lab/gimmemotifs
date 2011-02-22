/*
 * $Id: background.h 776 2006-05-10 17:27:13Z cegrant $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 18:34:37  nadya
 * Initial revision
 *
 */

/**********************************************************************/
/*
	MEME background models
*/
/**********************************************************************/

#ifndef BACKGROUND_H
#define BACKGROUND_H

/* maximum size of background model */
// 475254 is required for a third-order protein model
#define MAX_BACK_SIZE  475254

/* compute the log probability of a substring of a sequence given the log
   cumulative probability in lcb:
   log Pr(S_{i,...,i+w-1} | H_0) 
*/
#define Log_back(lcb, i, w) (lcb[(i)+w] - lcb[i])

extern double *read_markov_model(
  char *pfile,					/* name of probability file */
  double *freq,					/* letter frequencies */
  char *alpha,					/* alphabet expected */
  BOOLEAN add_x,				/* add x-tuples if TRUE */
  BOOLEAN rc,					/* average reverse complements*/
  int *order					/* order of model read */
);

extern char *i2s(
  int index
);

extern int s2i(
  char *string                                  /* the string to index */
);

extern int p2i(
  char *string,                                 /* the string to index */
  int w                                         // length of prefix
);

extern double *get_markov_from_sequence(
  char *seq,                    // the raw ASCII sequence
  char *alpha,                  // alphabet expected; should contain "X"
  BOOLEAN rc,                   // average reverse complements if TRUE
  int order,                    // order of Markov model to create
  double pseudocount            // the pseudocount
);

extern double log_cum_back(
  char *seq,					/* the sequence */
  double *a_cp,					/* the conditional probs */
  int order,					/* the order of the model */
  double *logcumback				/* the result */
);

#endif

