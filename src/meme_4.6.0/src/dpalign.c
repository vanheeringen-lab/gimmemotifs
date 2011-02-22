/*
 * $Id: dpalign.c 61 2005-07-29 00:16:55Z nadya $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 00:18:53  nadya
 * Initial revision
 *
 */

/*#define DEBUG*/
/* dpalign.c */
/************************************************************************
*	Copyright							*
*	(1999) The Regents of the University of California.		*
*	All Rights Reserved.						*
*	Author: Timothy L. Bailey
*									*
*	Permission to use, copy, modify, and distribute any part of 	*
*	this software for educational, research and non-profit purposes,*
*	without fee, and without a written agreement is hereby granted, *
*	provided that the above copyright notice, this paragraph and 	*
*	the following three paragraphs appear in all copies.		*
*									*
*	Those desiring to incorporate this software into commercial 	*
*	products or use for commercial purposes should contact the 	*
*	Technology Transfer Office, University of California, San Diego,*
*	9500 Gilman Drive, La Jolla, California, 92093-09m, 		*
*	Ph: (619) 534 5815.						*
*									*
*	IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO 	*
*	ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR 	*
*	CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF 	*
*	THE USE OF THIS SOFTWARE, EVEN IF THE UNIVERSITY OF CALIFORNIA 	*
*	HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 		*
*									*
*	THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE *
*	UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE 		*
*	MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  *
*	THE UNIVERSITY OF CALIFORNIA MAKES NO REPRESENTATIONS AND 	*
*	EXTENDS NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED, *
*	INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 	*
*	MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT 	*
*	THE USE OF THE MATERIAL WILL NOT INFRINGE ANY PATENT, 		*
*	TRADEMARK OR OTHER RIGHTS.  					*
************************************************************************/

/* 10-3-99 tlb; replace nrdb with back */

#ifdef DPALIGN_SO
#define DEFINE_GLOBALS
#endif
#include "meme.h"
#include "general.h"

/*
  data structures
*/
typedef struct dp_table_entry {		/* entry in dp_table */
  double v;				/* score: max(e, f, g) */
  double e;				/* score when gap in first seq here */
  double f;				/* score when gap in second seq here */
  double g;				/* score when no gap at this position */
} DP_TABLE_ENTRY;
typedef DP_TABLE_ENTRY **DP_TABLE;	/* rows: seq1 colums: seq2 */

#define MATCH(i, j) logodds((i)-1, (int)eseq[(j)-1])


/******************************************************************************/
/*
	dp_align

	Align a sequence to a logodds matrix using dynamic programming.

	Returns a string containing the two aligned sequences each
	terminated by a null.

*/
/******************************************************************************/
extern char *dp_align(
  int alength,					/* length of alphabet */
  int n,					/* length of logodds matrix */
  LOGODDS logodds,				/* logodds matrix */
  char *seq1,					/* consensus of logodds */
  int m,					/* length of sequence */
  char *eseq,					/* integer encoded sequence */
  double wg,					/* gap cost (initialization) */ 
  double ws,					/* space cost (extension) */ 
  BOOLEAN endgaps 				/* penalize end gaps if TRUE */
)
{
  int i, j, k, ii;
  int gaptype;					/* used in traceback */
  DP_TABLE dpt;					/* dyn. prog. table */
  int best_i, best_j;				/* best endpoint */
  char *align1=NULL, *align2=NULL;

  /* create the dynamic programming table */
  create_2array(dpt, DP_TABLE_ENTRY, n+1, m+1);

  /* initialize the table */
  dpt[0][0].v = dpt[0][0].f = 0;
  for (i=0; i<=n; i++)dpt[i][0].v = dpt[i][0].f = dpt[i][0].g = dpt[i][0].e = 0;
  for (j=0; j<=m; j++)dpt[0][j].v = dpt[0][j].f = dpt[0][j].g = dpt[0][j].e = 0;
  if (endgaps) {				/* penalize end gaps */
    for (i=1; i<=n; i++) dpt[i][0].v = dpt[i][0].e = -wg - i*ws;
    for (j=1; j<=m; j++) dpt[0][j].v = dpt[0][j].f = -wg - j*ws;
  } else {					/* don't penalize end gaps */
    for (i=1; i<=n; i++) dpt[i][0].v = dpt[i][0].e = 0;
    for (j=1; j<=m; j++) dpt[0][j].v = dpt[0][j].f = 0;
  }

  /* fill up the table row by row; note: row in logodds matrix corresponds
     to a position in seq1 
  */
  best_i = 1; best_j = m;
  for (i=1; i<=n; i++) {			/* row: logodds matrix col */
    for (j=1; j<=m; j++) {			/* column: seq2 */
      dpt[i][j].e = MAX(dpt[i][j-1].e, dpt[i][j-1].v-wg)-ws;
      dpt[i][j].f = MAX(dpt[i-1][j].f, dpt[i-1][j].v-wg)-ws;
      dpt[i][j].g = dpt[i-1][j-1].v + MATCH(i,j);
      dpt[i][j].v = MAX(dpt[i][j].e, MAX(dpt[i][j].f, dpt[i][j].g));
      if (!endgaps && (j==m || i==n) && dpt[i][j].v>dpt[best_i][best_j].v) {
        best_i = i; best_j = j; 
      }
    } /* column */
  } /* row */

#ifdef DEBUG
  /* print the dp table */
  printf("dp table\n");
  for(i=-1; i<=m; i++) printf("%-3c %3.3s  ", i==0 ? '-' : i>0 ? unhash(eseq[i-1]) : ' ', " ");
  printf("\n");
  for (i=0; i<=n; i++) {
    printf("%3c ", i>0 ? seq1[i-1] : '-');
    for (j=0; j<=m; j++) {
      printf("%3.0f %3.0f  ", dpt[i][j].g, dpt[i][j].f);
    }
    printf("\n");
    printf("%3.3s ", " ");
    for (j=0; j<=m; j++) {
      printf("%3.0f %3.0f  ", dpt[i][j].e, dpt[i][j].v);
    }
    printf("\n\n");
  }
  printf("best_i = %d best_j = %d\n", best_i, best_j);
#endif

  /* 
    trace the optimal alignment and put into arrays align1 and align2 
  */
  Resize(align1, m+n+1, char);
  Resize(align2, m+n+1, char);
  k = 0;					/* pointer into alignment */
  gaptype = 0;
  if (!endgaps) {				/* print end gap */
    for (i=n; i>best_i;) {
      align1[k] = seq1[--i];
      align2[k++] = ' ';
      gaptype = (dpt[i][m].v-wg-ws == dpt[i+1][m].f) ? 0 : 2;
    }
    for (j=m; j>best_j;) {
      align1[k] = ' ';
      align2[k++] = unhash(eseq[--j]);
      gaptype = (dpt[n][j].v-wg-ws == dpt[n][j+1].e) ? 0 : 1;
    }
  } /* endgap */
  for (i=endgaps?n:best_i, j=endgaps?m:best_j; i>0 && j>0 && k<=m+n; ) {
#ifdef DEBUG
    printf("i %d j %d k %d e %f f %f g %f v %f gt %d %c %c\n", i,j,k,
      dpt[i][j].e, dpt[i][j].f, dpt[i][j].g, dpt[i][j].v, gaptype, seq1[i-1], 
      unhash(eseq[j-1]));
#endif
    switch (gaptype) {
      case 0:					/* not in gap */
	if (dpt[i][j].g > dpt[i][j].e) {
	  if (dpt[i][j].g > dpt[i][j].f) {	/* match */
	    align1[k] = seq1[--i];
	    align2[k++] = unhash(eseq[--j]);
	  } else {				/* gap in seq2 */
	    align1[k] = seq1[--i];
	    align2[k++] = ' ';
            gaptype = (dpt[i][j].v-wg-ws == dpt[i+1][j].f) ? 0 : 2;
	  }
	} else {
	  if (dpt[i][j].e > dpt[i][j].f) {	/* gap in seq1 */
	    align1[k] = ' ';
	    align2[k++] = unhash(eseq[--j]);
            gaptype = (dpt[i][j].v-wg-ws == dpt[i][j+1].e) ? 0 : 1;
	  } else {				/* gap in seq2 */
	    align1[k] = seq1[--i];
	    align2[k++] = ' ';
            gaptype = (dpt[i][j].v-wg-ws == dpt[i+1][j].f) ? 0 : 2;
	  }
	}
        break;
      case 1:					/* in seq1 gap */
	align1[k] = ' ';
	align2[k++] = unhash(eseq[--j]);
	gaptype = (dpt[i][j].v-wg-ws == dpt[i][j+1].e) ? 0 : 1;
        break;
      case 2:					/* in seq2 gap */
	align1[k] = seq1[--i];
	align2[k++] = ' ';
	gaptype = (dpt[i][j].v-wg-ws == dpt[i+1][j].f) ? 0 : 2;
        break;
    } /* gaptype */
  } /* i, j, k */

  /* fill in remainder of seq1 */
  for (ii=i; ii>0; ii--) {
    align1[k] = seq1[ii-1];
    align2[k++] = ' ';
  }
  for (ii=j; ii>0; ii--) {
    align2[k] = unhash(eseq[ii-1]);
    align1[k++] = ' ';
  }

  /* reverse the alignment */
  for (i=0; i<k/2; i++) {
    SWAP(char, align1[i], align1[k-1-i]);
    SWAP(char, align2[i], align2[k-1-i]);
  }
  align1[k] = align2[k] = '\0';

  /* combine the alignments into one string */
  Resize(align1, 2*k+3, char);
  for (i=0, j=k+1; i<=k; i++, j++) align1[j] = align2[i];
  myfree(align2);
  /*printf("%s\n%s\n", align1, align1+k+1);*/

  /* free the dynamic programming table */
  free_2array(dpt, n+1);

  return align1;

} /* dp_align */

/******************************************************************************/
/*
	dp_ma_motif_seqs

	Perform a multiple alignment of a motif versus a set of sequences.

	Returns and an array of sequences.    
*/
/******************************************************************************/
extern char **dp_ma_motif_seqs(
  int alength,				/* length of sequence alphabet */
  int w,				/* length of motif */
  LOGODDS logodds,			/* motif log-odds matrix */
  char *cons,				/* consensus sequence of motif */
  int n,				/* number of sequences */
  char **eseqs,				/* (encoded) sequences to align */
  int *lengths,				/* sequence lengths */
  char **names,				/* sequence names */
  double wg,				/* gap cost (initialization) */ 
  double ws,				/* space cost (extension) */ 
  BOOLEAN endgaps			/* penalize end gaps if TRUE */
)
{
  int i, j, k, l, nspaces;
  char c;
  char **motif=NULL, **seq2=NULL;	/* aligned motif/sequence pairs */
  int *spaces=NULL;			/* max. number of spaces in alignmement						   preceeding position i of motif */
  int *pos1=NULL;			/* position of char i of motif in 
					   alignment */
  char **alignment;			/* the alignment */
  int maxl;				/* maximum sequence length */

  /* get maximum sequence length */
  for (i=0, maxl=w; i<n; i++) maxl = MAX(maxl, lengths[i]);

  /* create arrays of aligned pairs */
  Resize(motif, n, char *);
  Resize(seq2, n, char *);

  /* 
    align motif to each sequence 
  */
  for (i=0; i<n; i++) {
    char *pa;				/* pairwise alignment */
    pa = dp_align(alength, w, logodds, cons, lengths[i], eseqs[i], wg, ws, 
      endgaps);
    motif[i] = pa; 			/* start of motif */
    while (*pa) pa++;			/* find end of motif */
    seq2[i] = ++pa;			/* start of seq2 */
  }

  /* 
    combine alignments into one multiple alignment 
  */

  /* compute maximum number of spaces preceeding each position in motif */
  Resize(spaces, maxl+1, int);
  for (i=0; i<=w; i++) spaces[i] = 0;
  for (i=0; i<n; i++) {				/* alignment */
    for (j=k=nspaces=0; ; j++) {		/* position in aligned motif */
      c = motif[i][j];				/* character in aligned motif */
      if (c == ' ') {				/* within a gap */
        nspaces++;				/* size of gap */
      } else {					/* new position in motif */
        spaces[k] = MAX(spaces[k], nspaces);	/* max spaces in any alignment*/
        k++;					/* pointer in aligned motif */
        nspaces = 0;
      }
      if (!c) break;				/* end of aligned motif */
    } /* position in aligned motif */
  } /* alignment */

  /* create the alignment */
  for (i=nspaces=0; i<=w; i++) nspaces += spaces[i];	/* total spaces */ 
  create_2array(alignment, char, n+1, w+nspaces+1);
  for (i=0; i<n+1; i++) {
    for (j=0; j<w+nspaces; j++) {
      alignment[i][j] = ' ';
    }
    alignment[i][j] = '\0';
  }

  /* put motif in first line of alignment */
  Resize(pos1, maxl+1, int);
  for (i=j=0; i<=w; i++, j++) {
    j += spaces[i];				/* insert spaces before */
    pos1[i] = j;				/* where motif i is in align */
    alignment[0][j] = cons[i];			/* place letter */
  }

  /* put sequences into alignment */
  for (i=0; i<n; i++) {				/* seq2 */
    for (j=k=l=0; seq2[i][j]; j++, l++) {	/* j: position in pairwise al. 
						   k: position in motif
						   l: position in final algmnt.
						*/
      if (motif[i][j] != ' ') l = pos1[k++];	/* map motif_pos->final_pos */
      alignment[i+1][l] = seq2[i][j];
    }
  } /* seq2 */

#ifdef DEBUG
  /* print alignment */
  printf("alignment:\n");
  printf("%10.10s .%s\n", "motif", alignment[0]);
  for (i=1; i<n+1; i++) printf("%10.10s .%s\n", names[i-1], alignment[i]);
#endif

  /* free up all space */
  free_2array(motif, n);
  myfree(seq2);					/* sequences freed by motif */
  myfree(pos1);
  myfree(spaces);

  return alignment;

} /* dp_ma_motif_seqs */

/******************************************************************************/
/*
	g_align

	Find the longest alignment of width at least minw with g or fewer gaps 
	per column.  Different values of g are tried [0..] until a g-alignment 
	of width minw or greater is found.  The legal alignments must
	start in position left or greater, and end in position right or less.

	Returns the starting point and width of the g-alignment, and sets
	off, the offset of the starting point of the g-alignment.
*/
/******************************************************************************/
extern int g_align(
  char **ma,				/* the multiple alignment */
  int nseq,				/* number of sequences in alignment 
					   (including the key sequence) */
  int ncol,				/* number of columns in alignment */
  int left,				/* leftmost allowed start */
  int right,				/* rightmost allowed end */
  int minw,				/* minimum width */
  int *off,				/* offset of g-alignment rel. to left */
  int *w				/* width of g-alignment */
)
{
  int i, j, pos;
  int *ngaps=NULL;			/* number of gaps in each column */
  int lalign;				/* length of g-aligment */
  int g;				/* allowed gaps per column */
  int best_w=0, best_end=left;		/* width and end of longest g-align */

  /* too few sequences? */
  if (nseq < 2) {
    *w = right-left+1;
    *off = 0;
    return(0);
  } 

  /* create arrays for gaps/col */
  Resize(ngaps, ncol, int);

  /* count the number of gapped sequences at each position; if there is
     a gap in the first sequence, subtract the number of gaps from nseq
  */ 
  for (i=0; i<ncol; i++) {			/* position in alignment */
    ngaps[i] = 0;
    for (j=1; j<nseq; j++) if (ma[j][i] == ' ') ngaps[i]++;
    if (ma[0][i] == ' ') ngaps[i] = nseq - ngaps[i] - 1;
  } /* position in alignment */

  /* 
    find the longest alignment with g or fewer gaps per column and width
    at least minw
  */
  for (g=0; g<nseq; g++) {			/* g */
    /* fill in the length of g-alignment ending at i */
    best_end = left;				/* assume start of alignm. */
    lalign = (!left && ma[0][0] != ' ' && ngaps[0] <= g) ? 1 : 0;
    lalign = 0;
    for (i=pos=best_w=0; i<ncol; i++) {		/* position in align & motif */
      if (pos>=left && ma[0][i] != ' ') lalign++;	/* alignment length */
      if (ngaps[i] > g) lalign = 0;		/* bad column--too many gaps */
      if (pos>=left && lalign>best_w) {		/* within motif and better */
	best_end = pos;
	best_w = lalign;
      }
      if (ma[0][i] != ' ') pos++;		/* next position in motif */
      if (pos>right) break;			/* at right edge */
    } /* position in alignment & motif */
    if (best_w >= minw) break;			/* found a legal one */
  } /* g */
  best_w = MAX(best_w, minw);			/* make sure at least min_w */

  /* set up return values */
  *w = best_w;
  *off = best_end-best_w+1-left;		/* relative to leftmost start */

/*
  printf("g_align: Old w %d nsites %d Best w = %d offset = %d with g %d\n", 
    right-left+1, nseq, *w, *off, g);
*/

  /* free space */
  myfree(ngaps);

  return(g);
} /* g_align */

/******************************************************************************/
/*
	dp_multi_align

	Perform a multiple alignment of a motif versus a set of sequences.
	Uses MEME data structures.  Logodds matrix uses back as
	background frequencies.

	Returns an array of null terminated, aligned sequences.
*/
/******************************************************************************/
extern char **dp_multi_align(
  P_PROB sites,					/* the sites */
  int nsites,					/* the number of sites */
  int w,					/* width of sites */
  int flank,					/* add flank cols on left+rgt */
  DATASET *dataset 				/* the dataset */
)
{
  int i;
  SAMPLE **samples = dataset->samples;		/* sequences */
  int alength = dataset->alength;		/* dataset */
  double wg = dataset->wg;			/* gap cost (initialization) */ 
  double ws = dataset->ws; 			/* space cost (extension) */ 
  BOOLEAN endgaps = dataset->endgaps; 		/* penalize end gaps if TRUE */
  THETA lomap = dataset->lomap;			/* seq to theta logodds map */
  char **eseqs=NULL;				/* encoded sequences */
  int *lengths=NULL;				/* lengths of sequences */
  char **names=NULL;				/* names of sequences */
  LOGODDS lo=NULL;				/* logodds matrix */
  char *seq = NULL;				/* ascii for first sequence  */
  char **ma;					/* multiple aligment */

  /* create the array of sequences */
  Resize(eseqs, nsites, char *);
  Resize(lengths, nsites, int); 
  Resize(names, nsites, char *);

  /* prepare input for dp_ma_motif_seqs */
  for (i=0; i<nsites; i++) {			/* site */
    int x = sites[i].x;				/* sequence of site */
    int y = sites[i].y;				/* position in sequence */
    BOOLEAN ic = sites[i].ic;			/* site on - strand? */
    SAMPLE *s= samples[x];			/* sample */
    int lseq = s->length;			/* length of sequence */
    char *eseq = ic ? s->resic : s->res;	/* encoded sequence */
    int left = MAX(0, y-flank);			/* left edge of subsequence */
    int right = MIN(lseq-1, y+w+flank-1);	/* right edge of subsequence */

    eseqs[i] = eseq+left;			/* subsequence */
    lengths[i] = right-left+1;			/* subsequence length */
    names[i] = s->sample_name;			/* sequence name */
  } /* site */

  /* convert the first site to a logodds matrix */
  create_2array(lo, double, lengths[0], alength+1);
  init_theta(lo, eseqs[0], lengths[0], lomap, alength);

  /* get the ascii version of the first site */
  Resize(seq, lengths[0]+1, char);
  for (i=0; i<lengths[0]; i++) seq[i] = unhash(eseqs[0][i]);
  seq[i] = '\0';

  /* get multiple alignment */
  ma = dp_ma_motif_seqs(alength, lengths[0], lo, seq, 
    nsites-1, eseqs+1, lengths+1, names+1, wg, ws, endgaps);

  /* free space */
  free_2array(lo, lengths[0]);
  myfree(eseqs);
  myfree(lengths);
  myfree(names);
  myfree(seq);

  /* return multiple alignment */
  return ma;

} /* dp_multi_align */


#ifdef DPALIGN_SO
/******************************************************************************/
/*
	dpalign
*/
/******************************************************************************/
int main(
  int argc,
  char** argv
)
{
  int i, j;
  double wg=11, ws=1, m=1;
  BOOLEAN endgaps = FALSE;
  int alen;
  int w;
  char *alphabet = PROTEINB;
  int length;
  char *eseq=NULL;
  LOGODDS logodds;
  char *seq1, *seq2;
  char *pa;

  i = 1;
  argv[0] = "dpalign";
  DO_STANDARD_COMMAND_LINE(2,
    USAGE(<seq1> <seq2> [options]);
    USAGE(\n\t<seq1>\t\tfirst sequence to align);
    USAGE(\t<seq2>\t\tsecond sequence to align);
    NON_SWITCH(1, \r,
      switch (i++) {
        case 1: seq1 = _OPTION_; break;
        case 2: seq2 = _OPTION_; break;
        default: COMMAND_LINE_ERROR;
      }
    );
    DATA_OPTN(1, wg, <wg>, \tgap cost, wg = atof(_OPTION_));
    DATA_OPTN(1, ws, <ws>, \tspace cost, ws = atof(_OPTION_));
    DATA_OPTN(1, m, <ws>, \tmatch score, m = atof(_OPTION_));
    FLAG_OPTN(1, eg, \tpenalize end gaps, endgaps = TRUE);
    USAGE(\tGlobally align two sequences.);
    USAGE(\n\tCopyright);
    USAGE(\t(1999) The Regents of the University of California);
    USAGE(\tAll Rights Reserved.);
    USAGE(\tAuthor: Timothy L. Bailey);
  );

  /*
    Set up the hashing functions for mapping between letters or codons and
    alphabet positions.
  */
  setup_hash_alph(PROTEINB);                    /* PROTEINB to position hash */
  alen = strlen(alphabet);

  /* convert seq1 to logodds matrix */
  w = strlen(seq1);
  create_2array(logodds, LOGODDSB, w, alen);
  for (i=0; i<w; i++) {
    for (j=0; j<alen; j++) {
      logodds(i,j) = hash(seq1[i]) == j ? m : -m; 
    }
  }

  /* convert seq2 to integer encoding */
  length = strlen(seq2);
  Resize(eseq, length, char);
  for (i=0; i<length; i++) eseq[i] = hash(seq2[i]);

  /* align the sequences */
  pa = dp_align(alen, w, logodds, seq1, length, eseq, wg, ws, endgaps);

} /* main */
#endif
