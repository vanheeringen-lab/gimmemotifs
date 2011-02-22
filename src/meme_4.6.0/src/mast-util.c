/*
 * $Id: mast-util.c 4837 2010-08-25 17:16:26Z cegrant $
 * 
 * $Log$
 * Revision 1.3  2006/03/08 20:50:11  nadya
 * merge chamges from v3_5_2 branch
 *
 * Revision 1.2.4.1  2006/01/24 20:44:08  nadya
 * update copyright
 *
 * Revision 1.2  2005/10/20 00:20:52  tbailey
 * *** empty log message ***
 *
 * Revision 1.1.1.1  2005/07/29 17:17:18  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/************************************************************************
*                                                                       *
*       MAST                                                            *
*       Author: Timothy L. Bailey                                       *
*                                                                       *
*       Copyright                                                       *
*       (1994 - 2006) The Regents of the University of California.      *
*       All Rights Reserved.                                            *
*                                                                       *
*       Permission to use, copy, modify, and distribute any part of     *
*       this software for educational, research and non-profit purposes,*
*       without fee, and without a written agreement is hereby granted, *
*       provided that the above copyright notice, this paragraph and    *
*       the following three paragraphs appear in all copies.            *
*                                                                       *
*       Those desiring to incorporate this software into commercial     *
*       products or use for commercial purposes should contact the      *
*       Technology Transfer Office, University of California, San Diego,*
*       9500 Gilman Drive, La Jolla, California, 92093-0910,            *
*       Ph: (858) 534 5815.                                             *
*                                                                       *
*       IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO     *
*       ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR         *
*       CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF   *
*       THE USE OF THIS SOFTWARE, EVEN IF THE UNIVERSITY OF CALIFORNIA  *
*       HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.             *
*                                                                       *
*       THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE *
*       UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE          *
*       MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  *
*       THE UNIVERSITY OF CALIFORNIA MAKES NO REPRESENTATIONS AND       *
*       EXTENDS NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED, *
*       INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF        *
*       MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT    *
*       THE USE OF THE MATERIAL WILL NOT INFRINGE ANY PATENT,           *
*       TRADEMARK OR OTHER RIGHTS.                                      *
************************************************************************/
#include "mast.h"

#define BCHUNK 100              /* maximum size of a block diagram entry */
#define BLIMIT 100              /* min. remainder before adding to block array */
#define SWATHSIZE 100000        // Size of swath for hit_list

/*
        Local Declarations
*/
SCORE **score_it(
  BOOLEAN xlate_dna,    /* database is DNA and motifs protein */
  BOOLEAN reverse_comp, /* same effect as reverse complementing the motifs */
  LO *los[],            /* array of pointers to log-odds matrices */
  int nmotifs,          /* number of motifs */
  char *sequence,       /* sequence */
  long length           /* length of the sequence */
);

static BOOLEAN best_hit(
  int index,            /* Position of the given motif. */
  int motif,            /* Motif number. */
  int *hits,            /* Motif indices. */
  double *pvalues,      /* Array of pvalues. */
  long length           /* Length of pvalue array. */
);

static TILING tile_sequence(
  STYPE stype,                          /* treatment of strands of DNA */
  double **pv,                          /* p-value tables for each motif */
  double thresh,                        /* mark hits with p-value < thresh */
  BOOLEAN use_seq_p,                    /* use sequence not position p-values */
  LO *los[],                            /* array of pointers to lo matrices */
  int nmotifs,                          /* number motifs read */ 
  long length,                          /* length of sample */
  SCORE **scores,                       /* scores for each motif vs seq. pos. */
  int *hits,                            // non-overlapping hits
  double *pvalues,                      // hit p-values
  int *svalues,                         // scaled hit scores
  int start,                            // begin tiling at this position in sequence
  BOOLEAN hit_list                      // just update hits and pvalues
);

/*
        Procedures
*/

/**********************************************************************/
/*
        score_sequence

        Compute the scores for each motif and a given sequence.

        Returns *scores[] : scores[i][j].score = score of motif i at position j
        Returns *scores[] : scores[i][j].ic = score on - strand
*/
/**********************************************************************/
extern SCORE **score_sequence(
  STYPE stype,          /* how to treat different strands of DNA */
  BOOLEAN neg_strand,   /* the strand to score */
  BOOLEAN xlate_dna,    /* database is DNA and motifs protein */
  char *sequence,       /* the sequence to score */
  long length,          /* length of the sequence */
  int nmotifs,          /* number of motifs */
  LO *los[]             /* array of pointers to log-odds matrices */
)
{
  long i, j, k;
  SCORE **scores;       /* scores[i][j].score = motif i at position j */

  /* 
    Score the sequence with each motif on positive strand
  */
  scores = score_it(xlate_dna, neg_strand, los, nmotifs, sequence, length);

  /* 
    Hash reverse complement sequence and score it; 
    combine scores saving the MAX score for both strands at each position
  */
  if (stype == Combine) {               /* combine scores from both strands */
    SCORE **icscores;                   /* scores on - strand */
    icscores = score_it(xlate_dna, !neg_strand, los, nmotifs, sequence, length); 

    /* combine scores */
    for (i = 0; i < nmotifs; i++) {                 /* motif */
      long last_j = length - los[i]->ws;                   /* last possible site */

      for (j = 0; j <= last_j; j++) { /* site */
        if (icscores[i][j].score > scores[i][j].score) {
          scores[i][j].score = icscores[i][j].score;
          scores[i][j].ic = icscores[i][j].ic;
        } 
      } /* site */
    } /* motif */
    free_2array(icscores, nmotifs);
  }

  return scores;
} /* score_sequence */

/**********************************************************************/
/*
        tile_sequence

        The sequence is tiled with motif occurrences with p-values > thresh
        such that:
          1) occurrences do not overlap
          2) smaller numbered motif occurrences place first
          3) new occurrence replaces old occurrences it overlaps if its p-value
             is less than the product of their p-values

        Returns a tiling structure containing
                hits =
                  +motif number,        hit on + strand
                  -motif number,        hit on reverse complement, 
                  0,                    no hit
                pvalues                 position p-value of hit
                pv                      p-value of the product of p-values
*/
/**********************************************************************/
static TILING tile_sequence(
  STYPE stype,                          /* treatment of strands of DNA */
  double **pv,                          /* p-value tables for each motif */
  double thresh,                        /* mark hits with p-value < thresh */
  BOOLEAN use_seq_p,                    /* use sequence not position p-values */
  LO *los[],                            /* array of pointers to lo matrices */
  int nmotifs,                          /* number motifs read */ 
  long length,                          /* length of sample */
  SCORE **scores,                       /* scores for each motif vs seq. pos. */
  int *hits,                            // non-overlapping hits
  double *pvalues,                      // hit p-values
  int *svalues,                         // scaled hit scores
  int start,                            // begin tiling at this position in sequence
  BOOLEAN hit_list                      // just update hits and pvalues
)
{
  long j, k;
  int imotif;
  int maxws = 0;                        /* maximum width in seq */
  TILING tiling;                        /* tiling record */
  double prod_of_p;                     /* product of sequence p-values */
  int smotifs;                          /* number of motifs scored */

  /*
    mark non-overlapping hits for all motifs
  */
  prod_of_p = 1;
  for (imotif=smotifs=0; imotif<nmotifs; imotif++) {    /* motif */
    LO *lo = los[imotif];               /* logodds matrix of motif */
    int ws = lo->ws;                    /* width of current motif in sequence */
    double pvalue;                      /* position p-value of score */
    double seq_pvalue;                  /* sequence p-value of score */
    double best_score = LITTLE;         /* best score for motif */
    long n = length - ws + 1;           /* possible positions for score */

    if (!scores[imotif]) continue;      /* skip motif if too long */
    smotifs++;                          /* motif was scored */

    maxws = MAX(ws, maxws);             /* maximum motif width in seq so far */

    for (j=start; j <= length - ws; j++) {      /* site start in sequence */
      int score = scores[imotif][j].score;      /* score of site */
      int ic = scores[imotif][j].ic;    /* site on reverse strand */

      /* update best score for motif */
      if (score > best_score) best_score = score;

      /* get position p-value of score */
      pvalue = pv[imotif][(int) score];
      // Don't correct even if min of two pvalues so will match FIMO p-value.
      // if (stype == Combine) EV(pvalue, 2, pvalue);

      /* get sequence p-value of score if using instead of position p-value */
      if (use_seq_p) EV(pvalue, n, pvalue);

      /* 
        create list of non-overlapping motifs for this motif
      */
      if (pvalue < thresh) {
        BOOLEAN ok_to_mark = TRUE;      /* make motif if true */
        double prod = 1;                /* product of overlapped p-values */
        long first = MAX(0,j-maxws+1);  /* last overlap from left */
        long last = MIN(j+ws, length);  /* past last overlap to right */

        /* 
          get product of p-values motifs this would overlap to right 
        */
        for (k=j; k<last && ok_to_mark; k++) {
          if (hits[k] != 0) {           /* motif already here? */
            prod *= pvalues[k];
            if (pvalue >= prod) ok_to_mark = FALSE;
          }
        }

        /* 
          get product of p-values motifs this would overlap from left 
        */
        for (k=first; k<j && ok_to_mark; k++) {
          int m = abs(hits[k])-1;
          if (m >= 0 && los[m]->ws > j-k) {/* motif already here? */
            prod *=  pvalues[k];
            if (pvalue >= prod) ok_to_mark = FALSE;
          }
        }

        /* 
          mark motif if ok 
        */
        if (ok_to_mark) {
          hits[j] = ic ? -(imotif+1) : imotif+1;
          pvalues[j] = pvalue;
          svalues[j] = scores[imotif][j].score;
          /* remove overlapped motifs on right */
          for (k=j+1; k<last; k++) {
            hits[k] = 0;
            pvalues[k] = 1;
            svalues[k] = 0;
          }

          /* 
            remove overlapped motifs on left 
          */
          for (k=first; k<j; k++) {
            int m = abs(hits[k])-1;
            if (m >= 0 && los[m]->ws > j-k) {
              hits[k] = 0;
              pvalues[k] = 1;
              svalues[k] = 0;
            }
          }
        } /* mark motif */

      } /* pvalue < thresh */
    } /* site start */

    /* get sequence p-value of best score and take product of p-values */
    pvalue = pv[imotif][(int) best_score];
    EV(pvalue, (stype == Combine ? 2*n : n), seq_pvalue);
    prod_of_p *= seq_pvalue;
  } /* imotif */

  /* 
    return the sequence tiling
  */
  tiling.hits = hits;
  tiling.pvalues = pvalues;
  tiling.svalues = svalues;
  tiling.pv = qfast(smotifs, prod_of_p);

  return(tiling);
} /* tile_sequence */

/**********************************************************************/
/*
        create_diagram

        Create a block diagram of the motifs in the sequence.

        Returns a block diagram showing the order and spacing of the hits
        where hits are either
                strong          p-value < thresh 
                weak            otherwise.
*/
/**********************************************************************/
extern char *create_diagram(
  BOOLEAN dna,                          /* database is DNA */
  STYPE stype,                          /* treatment of strands of DNA */
  BOOLEAN xlate_dna,                    /* database is DNA and motifs protein */
  BOOLEAN best_motifs,                  /* diagrams have only best motif */
  BOOLEAN print_p,                      /* print p-value in block */
  double thresh,                        /* strong hit threshold */
  int nmotifs,                          /* number of motifs */
  LO *los[],                            /* array of pointers to lo matrices */
  long length,                          /* length of sample */
  BOOLEAN hit_list,                     /* create hit list instead of diagram */
  char *name,                           /* name of sequence; only used if hit_list=TRUE */
  int skip,                             // skip first "skip" positions in sequence
  int offset,                           // offset of start of sequence to add to position
  TILING tiling                         /* tiling of sequence */
)
{
  long i;
  int c = 0;                                    /* current position */
  int spacer = 0;                               /* current spacer */
  char tmp[100];                                /* scratch space */
  char *blocks=NULL;                            /* motif diagram */
  int bsize = 0;                                /* current size of blocks array */
  int *hits = tiling.hits;                      /* non-overlapping hits */
  double *pvalues = tiling.pvalues;             /* hit p-values */
  int *svalues = tiling.svalues;                /* scaled hit scores */

  /* create the diagram */
  for (i=skip; i < length; i++) {
    int m = abs(hits[i])-1;                     /* motif of hit */
    BOOLEAN ic = hits[i] < 0;                   /* hit on reverse strand */
    char *strand = stype==Protein ? "" : (ic ? "-" : "+"); 
    int frame = xlate_dna ? i%3 + 1 : 0;        /* frame if translating DNA */

    if ((m >= 0) &&                             /* start of new motif */
      (!best_motifs || best_hit(i, m, hits, pvalues, length)) ) {
      if (spacer > 0) {                         /* spacer ending */
        if (bsize-c < BLIMIT) Resize(blocks, bsize+=BCHUNK, char);
        if (!hit_list) c += sprintf(blocks+c, "%d_", spacer);   /* spacer */
      }
      if (hit_list) {
        double scale = los[m]->scale;
        double off = los[m]->offset;
        int w = los[m]->w;
        if (bsize-c < BLIMIT) Resize(blocks, bsize+=BCHUNK, char);
        /* c += sprintf(blocks+c, " %s%d %ld %ld %8.2e,", */
        c += sprintf(blocks+c, "%s %s%d %ld %ld %8.2f %8.2e\n", 
          name, strand, los[m]->imotif, i+offset+1, i+offset+los[m]->ws, 
            scaled_to_bit(svalues[i], w, scale, off), pvalues[i]);
      } else {
        make_block(los[m]->imotif, strand, frame, thresh, pvalues[i], print_p, tmp);
        if (bsize-c < BLIMIT) Resize(blocks, bsize+=BCHUNK, char);
        c += sprintf(blocks+c, "%s_", tmp);     /* block */
      }
      spacer = -los[m]->ws + 1;                 /* account for width */
    } else {
      spacer++;                                 /* increase spacer */
    }
  }
  if (hit_list) {
    /* nothing to do */
  } else if (spacer > 0) {
    if (bsize-c < BLIMIT) Resize(blocks, bsize+=BCHUNK, char);
    sprintf(blocks+c, "%d", spacer);        /* final spacer */
  } else if (c>0) {
    blocks[c-1] = '\0';                     /* remove final dash or comma */
  } else {
    Resize(blocks, 1, char);
    blocks[0] = '\0';
  }

  return blocks;
} /* create_diagram */

/**********************************************************************/
/*
        print_diagram

        Print the motif block diagram.
*/
/**********************************************************************/
extern void print_diagram(
  char *dia,                            /* motif diagram string */
  char *hdr,                            /* prefix for each line of diagram */
  FILE *file                            /* destination file */
)
{
  int j;
  int dia_len = strlen(dia);                    /* length of diagram */
  int hlen = strlen(hdr);                       /* length of header */

  for (j=0; j < dia_len; ) {
    int remain = dia_len - j;                   /* left to print */
    int dlen;                                   /* room on line */
    char *h = ((j==0) ? hdr : " ");             /* current header */
    dlen = PAGEWIDTH - hlen - 6;
    if (remain <= PAGEWIDTH - hlen) dlen = remain;
    fprintf(file, "%-*.*s%.*s", hlen, hlen, h, dlen, dia+j);
    j += dlen;
    /* continue printing until a good breaking point */
    while (j < dia_len && dia[j-1] != '_' && dia[j-1] != ',') putc(dia[j++], file);
    putc('\n', file);
  }
} /* print_diagram */

/**********************************************************************/
/*
        score_it

        Score the sequence with each motif.

        Returns the array of scores:
                scores[i][j].score = score of motif i at position j
                scores[i][j].ic = score on - strand (set to reverse_comp)
*/
/**********************************************************************/
SCORE **score_it(
  BOOLEAN xlate_dna,    /* database is DNA and motifs protein */
  BOOLEAN reverse_comp, /* score opposing strand */
  LO *los[],            /* array of pointers to log-odds matrices */
  int nmotifs,          /* number of motifs */
  char *sequence,       /* sequence */
  long length           /* length of the sequence */
)
{
  int imotif;
  long i, j, k;
  SCORE **scores = NULL;        /* scores[i][j].score motif i at offset j */
  int *hash_seq = NULL;         /* hashed sequence */

  if (reverse_comp) {
    char *icseq = NULL;  
    /* reverse complement the sequence, note that the position of the motif
     * scores is adjusted to compensate later so the result looks like we've
     * RCd the motifs */                
    Resize(icseq, length+1, char);
    for (i = 0, j = length - 1; i < length; i++, j--) icseq[j] = comp_dna(sequence[i]);
    icseq[i] = '\0';
    //we don't need to refer to the original sequence
    sequence = icseq;
  }

  /* 
    Hash sequence to index in alphabet.
  */
  hash_seq = dhash_it(xlate_dna, los[0]->alen, sequence, length);

  /* 
    Create scores array. 
  */
  Resize(scores, nmotifs, SCORE *);

  /* 
    Score the sequence with each motif. 
  */
  for (imotif=0; imotif<nmotifs; imotif++) {    /* motif */
    LO *lo = los[imotif];               /* array of motifs */
    LOGODDS2 logodds2 = lo->logodds2;   /* motif matrix */
    int r = (lo->w+1)/2;                /* number of rows in the matrix */
    int ws = lo->ws;                    /* width the motif in the sequence */
    int inc1 = 1;                       /* increment to next site */
    int inc2 = 2*(xlate_dna ? 3 : 1);   /* increment to next hashed column */

    /* create array of scores unless sequence shorter than motif */
    scores[imotif] = NULL;
    if (ws > length) {
      continue;                         /* skip this motif */
    } else {
      Resize(scores[imotif], length, SCORE);
    }

    /* 
      Score each subsequence with the current motif.
    */
    for (j=0; j<=length-ws; j+=inc1) {  /* position */
      int *h;                           /* pointer in hash_seq */
      int score = 0;                    /* subsequence score */

      /* motif score (positive motif) */
      for (k=0, h=hash_seq+j; k<r; k++, h+=inc2) score += logodds2(k, *h);

      /* 3-class score */
      if (lo->pair) {
        int neg = 0;                    /* subsequence score, neg motif */
        double score3;                  /* 3-class bit score */
        for (k=r, h=hash_seq+j; k<2*r; k++, h+=inc2) neg += logodds2(k, *h);
        score3 = score3class(score, neg, lo);
        score = bit_to_scaled(score3, 1, lo->scale3, lo->offset3);
      }
      
      if (!reverse_comp) {
        scores[imotif][j].score = score;
        scores[imotif][j].ic = FALSE;
      } else {
        int real_pos = (length - ws) - j;//TODO double check this calculation
        scores[imotif][real_pos].score = score;
        scores[imotif][real_pos].ic = TRUE;
      }
    } /* position */

  } /* imotif */

  myfree(hash_seq);

  if (reverse_comp) myfree(sequence);
  return scores;
} /* score_it */

/**********************************************************************/
/*
        best_hit

        Determine whether a given motif occurrence has the lowest p-value
        for that motif.
*/
/**********************************************************************/
static BOOLEAN best_hit(
  int index,            /* Position of the given motif. */
  int motif,            /* Motif number. */
  int *hits,            /* Motif indices. */
  double *pvalues,      /* Array of pvalues. */
  long length           /* Length of pvalue array. */
)
{
  long i;
  double lowest_pvalue = 1.0;
  int best_index = 0;

  for (i = 0; i < length; i++) {
    int m = abs(hits[i])-1;                     /* motif of hit */
    if ((m == motif) && (pvalues[i] < lowest_pvalue)) {
      lowest_pvalue = pvalues[i];
      best_index = i;
    }
  }
  return (best_index == index);
} /* best_hit */
 
/**********************************************************************/
/*
        make_block

        Create a block string:

                [smf(p)] or <smf(p)>
                        s       strand (optional)
                        m       motif
                        f       frame (optional)
                        (p)     p-value (optional)      
*/
/**********************************************************************/
extern void make_block(
  int m,                        /* motif number */
  char *strand,                 /* strand */
  int f,                        /* frame number; f=0 not translating DNA */
  double thresh,                /* strong motif threshold */
  double p,                     /* p-value */
  BOOLEAN print_p,              /* print p-value in block if TRUE */
  char *block                   /* put block string here */
)       
{
  char left = p < thresh ? '[' : '<';
  char right = p < thresh ? ']' : '>';
  BOOLEAN xlate_dna = (f != 0);
  char *fnames = "abc";                                 /* frame 1=a, 2=b, 3=c*/

  if (print_p) {                                        /* print p-value */
    char *bfmt = f ? "%c%s%d%c(%8.2e)%c" : "%c%s%d(%8.2e)%c";
    if (xlate_dna) {                                    /* str., motif, frame */
      sprintf(block, bfmt, left, strand, m, fnames[f-1], p, right);
    } else {                                            /* strand, motif */
      sprintf(block, bfmt, left, strand, m, p, right);
    }
  } else {                                              /* don't print p-value */
    char *bfmt = f ? "%c%s%d%c%c" : "%c%s%d%c";
    if (xlate_dna) {                                    /* str., motif, frame */
      sprintf(block, bfmt, left, strand, m, fnames[f-1], right);
    } else {                                            /* strand, motif */
      sprintf(block, bfmt, left, strand, m, right);
    }
  }
} /* make_block */

/**********************************************************************/
/*
        score_tile_diagram

        Score a sequence, tile it with motif occurrences and create
        a block diagram string.

        Returns the tiling containing the hits, their position p-values,
        and the p-value of the product of p-values for the best hits.

        Caller is responsible for deallocating contents of TILING.
*/
/**********************************************************************/
extern TILING score_tile_diagram(
  FILE *mast_out,                       /* output */
  char *sequence,                       /* sequence to score and tile */
  long length,                          /* length of sequence */
  LO *los[],                            /* array of pointers to lo matrices */
  int nmotifs,                          /* number motifs read */
  BOOLEAN dna,                          /* database is DNA */
  STYPE stype,                          /* handling of different strands */
  BOOLEAN neg_strand,                   /* for dna sequences the negative strand may be scored and tiled */
  BOOLEAN xlate_dna,                    /* database is DNA and motifs protein */
  BOOLEAN best_motifs,                  /* only put the best motifs into the hit list */
  BOOLEAN print_p,                      /* print p-values in the block diagram */
  double **pv,                          /* p-value tables for each motif */
  double m_thresh,                      /* maximum motif p-value to print */
  double w_thresh,                      /* max. motif p-value for weak hits */
  BOOLEAN use_seq_p,                    /* use sequence not position p-values */
  BOOLEAN hit_list,                     /* create hit list instead of diagram */
  char *name                            /* name of sequence; only used if hit_list=TRUE */
)
{
  int i, j;
  SCORE **scores;                       /* scores for each motif vs seq. pos. */
  TILING tiling = {.hits=NULL, 
      .pvalues=NULL, .svalues=NULL, 
      .pv=0, .diagram=NULL};            /* tiling and diagram of sequence */
  int *hits = NULL;                     // non-overlapping hits
  double *pvalues = NULL;               // hit p-values
  int *svalues = NULL;                  // scaled hit scores
  int max_w;                            // maximum motif width
  int ssize;                            // size of swaths in hit_list mode

  // get maximum motif width
  for (i=max_w=0; i<nmotifs; i++) { if (los[i]->w > max_w) max_w = los[i]->w; }

  // define size of swath:
  //   !hit_list -- full sequence
  //   hit_list  -- SWATHSIZE + overlap but not longer than sequence length
  if (hit_list) {
    // check that SWATHSIZE is large enough
    if (SWATHSIZE <= 4*max_w) {
      die("Maximum motif width %d is too large; recompile MAST with larger SWATHSIZE\n", max_w);
    }
    ssize = MIN(SWATHSIZE + 2*max_w, length);
  } else {
    ssize = length;
  }

  // create storage for hits and pvalues
  Resize(hits, ssize, int);
  Resize(pvalues, ssize, double);
  Resize(svalues, ssize, int);

  // clear array of hits
  for (i=0; i<ssize; i++) hits[i] = 0;

  //
  // tile a swath
  //
  // Hits are printed immediately if hit_list == TRUE.
  //
  int offset = 0;                       // current position in sequence
  int prv_overlap = 0;                  // first position not yet printed
  while (length > 0) {                  // do a swath
    // figure out new size of overlap buffer
    int overlap = (ssize == length) ? 0 : 2*max_w; 
    
    // score the swath with each of the motifs
    scores = score_sequence(stype, neg_strand, xlate_dna, sequence+offset, ssize, nmotifs, los);

    // mark the non-overlapping motif positions
    tiling = tile_sequence(stype, pv, w_thresh, use_seq_p, los, nmotifs, ssize,
      scores, hits, pvalues, svalues, prv_overlap, hit_list);

    tiling.diagram = NULL;

    // add the block diagram to the tiling structure
    tiling.diagram = create_diagram(dna, stype, xlate_dna, best_motifs, print_p,
      m_thresh, nmotifs, los, ssize-overlap, hit_list, name, prv_overlap, offset, tiling);

    // hit_list: print hits in swath minus overlap region then shift all buffers
    if (hit_list) {                     // hit_list     
      
      // print hits
      if (tiling.diagram) {
        //printf(mast_out,"offset %d overlap %d prv_overlap %d ssize %d length %d\n",
         // offset, overlap, prv_overlap, ssize, length);
        if (mast_out) fprintf(mast_out, "%s", tiling.diagram);
        myfree(tiling.diagram);         // free diagram space
      }
      // The next position to be analyzed has index = ssize - overlap.
      // It requires seeing scores at positions with indices
      // as small as index = ssize - overlap - overlap.
      // So shift index ssize-2*overlap to position zero.
      int shift = ssize-2*overlap;
      for (i=0,j=shift; j<ssize; i++, j++) {    // shift hits, pvalues
        hits[i] = hits[j];
        pvalues[i] = pvalues[j];
      }
      offset += shift;                  // offset of start of swath
      length -= shift;                  // length of remaining sequence
      for (; i<MIN(ssize,length); i++) { // clear rest of hits
        hits[i] = 0;
        //pvalues[i] = 1;
      }
      ssize = MIN(SWATHSIZE + overlap, length);         // swath size
      prv_overlap = overlap;            // keep track to avoid reprinting hits
    } else {                            // ! hit_list
      length = 0;                       // all done
    }
    // free scores array
    free_2array(scores, nmotifs);
  } // swath

  return(tiling);
} /* score_tile_diagram */

/**********************************************************************/
/*
        qfast
        
        Calculate the p-value of the product of uniform [0,1] random
        variables.

*/
/**********************************************************************/
extern double qfast(
  int n,                        /* number of random variables in product */
  double k                      /* product of random variables */
)
{
  int i = 1;
  double mlnk, term, phi;
 
  if (n == 0) return 1.0;       /* worst possible p-value */
  if (k == 0) return 0.0;       /* p-value is 0 */

  mlnk = -log(k);
 
  phi = term = k;
  for (i=1; i<n; i++) {
    term *= mlnk/i;
    phi += term;
  }

  return phi;
} /* qfast */

/**********************************************************************/
/*
        init_mast_background

        Determine the BLAST alphabet corresponding to the input alphabet.
        Set up the hashing functions for the alphabet(s).
        Read in the background frequencies for letters (or use
        non-redundant database frequencies if no file name given).

        The background frequency file should contain frequencies
        for the DNA0 or PROTEIN0 alphabet in the format required
        by read_markov model.

        Sets
                blast_alphabet :        the corresponding BLAST alphabet
                p :                     permutation/substitution matrix
                                        from alphabet to blast_alphabet

        Returns the background frequencies for the BLAST alphabet.
*/
/**********************************************************************/
extern double *init_mast_background(
  char *bfile,                          /* name of background file */
  char *alphabet,                       /* motif alphabet */
  STYPE stype,                          /* handling of DNA strands */
  BOOLEAN translate_dna,                /* DNA sequences and protein motifs */
  char **blast_alphabet,                /* corresponding BLAST alphabet */
  int *p[MAXASCII]                      /* permutation/substitution matrix */
)
{
  int i;
  int alen, alen0;
  double *back=NULL, *tmp=NULL, *blast_freqs=NULL, *simple_freqs=NULL;
  char *bfile_alpha;                            /* alphabet in bfile */
  int order;                                    /* order of background */
  /* average reverse complement frequencies */
  // FIXME: tlb (17-jul-2008) Maybe we should use the reverse complement frequencies
  // for the background in Separate mode and, even, in Combine mode?
  BOOLEAN rc = (stype==Separate || stype==Combine);     

  /* 
    Set up the hashing functions for mapping between letters or codons and
    alphabet positions.
  */
  setup_hash_alph(DNAB);                        /* DNAB to position hashing */
  setup_hash_alph(PROTEINB);                    /* PROTEINB to position hash */

  /* 
    Determine the alphabet used by the motifs.
    Make a permutation/substitution mapping from that alphabet to one
    of the BLAST alphabets.
  */
  *blast_alphabet = get_blast_alphabet(alphabet, p);
  alen = strlen(*blast_alphabet);

  /*
    Set up to translate DNAB to PROTEINB if requested.
  */
  if (translate_dna) {
    if (strcmp(*blast_alphabet, PROTEINB)) {    /* motifs must be protein */
      die("\nThe -dna switch requires protein motifs.  "
          "Your motif(s) use a DNA alphabet.\n");
    }
    setup_hash_dnab2protb();                    /* DNAB to PROTEINB hashing */
  }

  /* 
    Determine the type of motif alphabet and set it as current
  */
  if (!strcmp(*blast_alphabet, DNAB)) { /* searching DNA with DNA motifs */
    blast_freqs =  ntfreq;
    bfile_alpha = DNA0;
    setalph(0);                         /* DNAB alphabet */
  } else {
    if (translate_dna) {                /* searching DNA with protein motifs */
      blast_freqs = frame0;
    } else {                            /* searching protein with protein mtfs*/
      blast_freqs = nrfreq;
    }
    rc = FALSE;                         /* don't average background freqs */
    bfile_alpha = PROTEIN0;
    setalph(1);                         /* PROTEINB alphabet */
  }

  // Convert from BLAST to simple alphabet if no bfile given.
  alen0 = strlen(bfile_alpha);          /* length of '0' alphabet */
  if (!bfile) {
    Resize(simple_freqs, alen, double);
    for (i=0; i<alen0; i++) simple_freqs[i] = blast_freqs[hash(bfile_alpha[i])];
  }

  /*
    Read in the background probabilities from a Markov model file or the
    tmp array using the basic alphabet.  Convert to the BLAST alphabet.
  */
  tmp = read_markov_model(bfile, simple_freqs, bfile_alpha, FALSE, rc, &order);
  Resize(back, alen, double);                   /* create results array */
  for (i=0; i<alen; i++) back[i] = 0;
  for (i=0; i<alen0; i++) back[hash(bfile_alpha[i])] = tmp[i];

  return(back);
} /* init_mast_background */

/**********************************************************************/
/*
        free_tiling

        Free a tiling structure.
*/
/**********************************************************************/
extern void free_tiling(
  TILING tiling
) 
{
    myfree(tiling.pvalues);
    myfree(tiling.svalues);
    myfree(tiling.diagram);

} /* free_tiling */

/**********************************************************************/
/*
        get_seq_comp

        Get the letter frequencies in a sequence.

        Doesn't take into account reverse complement.
*/
/**********************************************************************/
double *get_seq_comp(
  BOOLEAN xlate_dna,                    /* translate DNA */
  char *sequence,                       /* ASCII sequence */
  int alen                              /* length of alphabet */
)
{
  int i, n;
  double *freq = NULL;

  /* create the frequency array */ 
  Resize(freq, alen, double);
  for (i=0; i<alen; i++) freq[i] = 0;
  /*for (i=0; i<alen; i++) freq[i] = 1;*/               /* add one prior */

  /* count the number of letters of each type */
  for (n=0; sequence[n]; n++) {
    i = chash(xlate_dna, FALSE, sequence+n);    /* hash letter */
    freq[i]++;
  }

  /* convert counts to frequencies */
  //FIXME why doesn't the else return 1/alen or the background model value, wouldn't that be the prior with no other information?
  for (i=0; i<alen; i++) freq[i] = n ? freq[i]/n : 1.0;
  /*for (i=0; i<alen; i++) freq[i] = n ? freq[i]/(n+alength) : 1.0;*/   /* prior */

  /* return the frequency matrix */
  return(freq);
} /*  get_seq_comp */
