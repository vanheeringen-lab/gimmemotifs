/*
 * $Id: user.h 1048 2006-07-06 20:07:44Z cegrant $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 19:12:25  nadya
 * Initial revision
 *
 */

/*
	User settable parameters
*/

#ifndef user_h
#define user_h

#define MSN	 24 		/* maximum length of sample name */
				/* MSN + 40 < PAGEWIDTH (see meme.h) */
#define MAXALPH  28		/* maximum length of alphabet + 1 for 'X' */
#define MAXG 101		/* maximum number of motifs + 1 */
#define MAXSITE 300		/* maximum length of a site */
#define MINSITES 2		/* minimum number of sites in valid motif */
#define LLR_RANGE 200		/* range of scaled LLR statistic */

#define MINCONS 0.2		/* Display 'X' as consensus if no letter f > */ 
#define LOGOHEIGHT 7.5		// height of sequence logo in cm.
#define MAXLOGOWIDTH 30		// maximum width of sequence logo in cm.

/* minimum allowable motif width before shortening; 
   never make less than 2 or will crash! */
#define MIN_W 8
/* maximum allowable length before shortening */
#define MAX_W 50

#define MNAME 20		/* names of known motifs */
#define NMOTIFS MAXG		/* maximum number of known motifs */


/* default size of heap for branching search */
#define HSIZE 64
#define HS_DECREASE 2

/* default branching factor for branching search */
#define BFACTOR 3

/* Amount of error tolerated in probability column sums (they should sum to
   approximately 1): */
#define ERR_EPSILON 0.01


#endif
