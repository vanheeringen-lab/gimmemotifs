/*
 * $Id: mast.c 5280 2011-01-05 11:10:34Z james_johnson $
 * 
 * $Log$
 * Revision 1.5  2006/03/08 20:50:11  nadya
 * merge chamges from v3_5_2 branch
 *
 * Revision 1.4.4.2  2006/01/31 20:51:38  nadya
 * rm initialization for 'name'
 *
 * Revision 1.4.4.1  2006/01/31 20:22:48  nadya
 * init name with zeros
 *
 * Revision 1.4  2005/10/25 19:02:26  nadya
 * change c++ style comment to proper c
 *
 * Revision 1.3  2005/10/20 00:20:52  tbailey
 * *** empty log message ***
 *
 * Revision 1.2  2005/10/01 23:58:04  nadya
 * update documentation comment
 *
 * Revision 1.1.1.1  2005/07/29 18:06:43  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/************************************************************************
*                                                                       *
*       MAST                                                            *
*       Author: Timothy L. Bailey                                       *
*                                                                       *
*       Copyright                                                       *
*       (1994 - 2001) The Regents of the University of California.      *
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
*       Ph: (619) 534 5815.                                             *
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

/**********************************************************************/
/*
  mast <memefile> <database> [optional arguments...]
        See <install-path>/bin/mast for documentation.
*/
/**********************************************************************/

#include <errno.h>
#include <stdarg.h>
#include <string.h>

#define DEFINE_GLOBALS
#include "mast.h"

#include "array-list.h"
#include "buffer.h"
#include "dir.h"
#include "io.h"
#include "meme-io.h"
#include "alphabet.h"
#include "diagram.h"
#include "projrel.h"
#include "utils.h"
#include "xml-out.h"
#include "xml-util.h"

/* error flag */
BOOLEAN errored = FALSE;
/* printing */
VERBOSE_T verbosity = QUIET_VERBOSE;

const char *program_name = "mast"; /* the program name */
static char *default_out_dir = "mast_out";  /* default name of output */
static const char *XML_FILENAME = "mast.xml";
static const char *HTML_STYLESHEET = "mast-to-html.xsl";
static const char *HTML_FILENAME = "mast.html";
static const char *TXT_FILENAME = "mast.txt";
static const char *MAST2TXT_FILENAME = "mast2txt";
    
/* MAST DTD */
/*{{{*/
char* MAST_DTD = 
"<!DOCTYPE mast[\n"
"<!ELEMENT mast (model, alphabet, motifs, sequences, runtime)>\n"
"<!ATTLIST mast version CDATA #REQUIRED release CDATA #REQUIRED>\n"
"<!ELEMENT model (command_line, max_correlation, remove_correlated, strand_handling, translate_dna, max_seq_evalue,\n"
"    adj_hit_pvalue, max_hit_pvalue, max_weak_pvalue, host, when)>\n"
"<!ELEMENT command_line (#PCDATA)>\n"
"<!ELEMENT max_correlation (#PCDATA)>\n"
"<!ELEMENT remove_correlated EMPTY>\n"
"<!ATTLIST remove_correlated value (y|n) #REQUIRED>\n"
"<!ELEMENT strand_handling EMPTY>\n"
"<!ATTLIST strand_handling value (combine|separate|norc|protein) #REQUIRED>\n"
"<!ELEMENT translate_dna EMPTY>\n"
"<!ATTLIST translate_dna value (y|n) #REQUIRED>\n"
"<!ELEMENT max_seq_evalue (#PCDATA)>\n"
"<!ELEMENT adj_hit_pvalue EMPTY>\n"
"<!ATTLIST adj_hit_pvalue value (y|n) #REQUIRED>\n"
"<!ELEMENT max_hit_pvalue (#PCDATA)>\n"
"<!ELEMENT max_weak_pvalue (#PCDATA)>\n"
"<!ELEMENT host (#PCDATA)>\n"
"<!ELEMENT when (#PCDATA)>\n"
"<!ELEMENT alphabet (letter+)>\n"
"<!ATTLIST alphabet type (amino-acid|nucleotide) #REQUIRED bg_source (preset|file|sequence_composition) #REQUIRED bg_file CDATA #IMPLIED>\n"
"<!ELEMENT letter EMPTY>\n"
"<!ATTLIST letter symbol CDATA #REQUIRED ambig (y|n) \"n\" bg_value CDATA #IMPLIED>\n"
"<!ELEMENT motifs (motif+,correlation*,nos*)>\n"
"<!ATTLIST motifs source CDATA #REQUIRED name CDATA #REQUIRED last_mod_date CDATA #REQUIRED>\n"
"<!ELEMENT motif EMPTY>\n"
"<!-- num is simply the loading order of the motif, it's superfluous but makes things easier for XSLT -->\n"
"<!ATTLIST motif id ID #REQUIRED num CDATA #REQUIRED name CDATA #REQUIRED width CDATA #REQUIRED\n"
"   best_f CDATA #REQUIRED best_r CDATA #IMPLIED bad (y|n) \"n\">\n"
"<!-- for n > 1 motifs there should be (n * (n - 1)) / 2 correlations, obviously there are none for only 1 motif -->\n"
"<!ELEMENT correlation EMPTY>\n"
"<!ATTLIST correlation motif_a IDREF #REQUIRED motif_b IDREF #REQUIRED value CDATA #REQUIRED>\n"
"<!-- nos: Nominal Order and Spacing diagram, a rarely used feature where mast can adjust pvalues for an expected motif spacing -->\n"
"<!ELEMENT nos (expect+)>\n"
"<!-- length is in the same unit as the motifs, which is not always the same unit as the sequence -->\n"
"<!ATTLIST nos length CDATA #REQUIRED>\n"
"<!-- the expect tags are expected to be ordered by pos ascending -->\n"
"<!ELEMENT expect EMPTY>\n"
"<!ATTLIST expect pos CDATA #REQUIRED gap CDATA #REQUIRED motif IDREF #REQUIRED>\n"
"<!ELEMENT sequences (database+, sequence*)>\n"
"<!-- the database tags are expected to be ordered in file specification order -->\n"
"<!ELEMENT database EMPTY>\n"
"<!ATTLIST database id ID #REQUIRED num CDATA #REQUIRED source CDATA #REQUIRED name CDATA #REQUIRED last_mod_date CDATA #REQUIRED \n"
"    seq_count CDATA #REQUIRED residue_count CDATA #REQUIRED type (amino-acid|nucleotide) #REQUIRED link CDATA #IMPLIED>\n"
"<!-- the sequence tags are expected to be ordered by best combined p-value (of contained score tags) ascending -->\n"
"<!ELEMENT sequence (score+,seg*)>\n"
"<!ATTLIST sequence id ID #REQUIRED db IDREF #REQUIRED num CDATA #REQUIRED name CDATA #REQUIRED comment CDATA \"\" length CDATA #REQUIRED>\n"
"<!ELEMENT score EMPTY>\n"
"<!-- frame is the starting offset for translation of dna sequences which gives the lowest pvalues for the provided protein motifs -->\n"
"<!ATTLIST score strand (both|forward|reverse) #REQUIRED frame (a|b|c) #IMPLIED combined_pvalue CDATA #REQUIRED evalue CDATA #REQUIRED>\n"
"<!-- within each sequence the seg tags are expected to be ordered by start ascending -->\n"
"<!ELEMENT seg (data,hit+)>\n"
"<!ATTLIST seg start CDATA #REQUIRED>\n"
"<!ELEMENT data (#PCDATA)>\n"
"<!-- within each seg the hit tags are expected to be ordered by pos ascending and then forward strand first -->\n"
"<!ELEMENT hit EMPTY>\n"
"<!-- gap, while superfluous, makes creating motif diagrams for the text version much easier when using XSLT -->\n"
"<!ATTLIST hit pos CDATA #REQUIRED gap CDATA #REQUIRED motif IDREF #REQUIRED pvalue CDATA #REQUIRED strand (forward|reverse) \"forward\" \n"
"    match CDATA #REQUIRED translation CDATA #IMPLIED>\n"
"<!ELEMENT runtime EMPTY>\n"
"<!ATTLIST runtime cycles CDATA #REQUIRED seconds CDATA #REQUIRED>\n"
"]>\n";
/*}}}*/

/* these two macros are used to convert a preprocessor define into a string
 * you would logicaly think that it would work with just the one macro but
 * for some odd reason it doesn't */
#define STRINGIZE2(s) #s
#define STRINGIZE(s) STRINGIZE2(s)

#ifndef EXP
#define EXP 0
#else
#define EXP 1
#endif

/* default e-thresh */
#define EXPECT 10

/* maximum pairwise motif correlation */
#define MAXCORR 0.6

/* number of different integral score values */
#define MAST_RANGE 100

/* the smallest size of a segment in the output */
#define SEG_CHUNK 75
#define SEG_FORMAT "%." STRINGIZE(SEG_CHUNK) "s\n"

/*
  Data Definitions (structs)
*/

typedef struct database_t DATABASE_T;
typedef struct strand_t STRAND_T;
typedef struct sseq_t SSEQ_T;

/*{{{*/
struct database_t {
  int index;            /* loading number of the database */
  char *source;         /* the source of the database "-" for stdin otherwise the file name */
  char *name;           /* the display name of the database */
  time_t last_mod;      /* the unix time the database was last modified */
  int sequence_count;   /* the count of sequences in this database */
  long residue_count;   /* the count of residues in this database */
  BOOLEAN is_dna;       /* is the database a nucleotide alphabet? */
  char *link;           /* the link to search for further information on sequences */
  FILE *file;           /* the file to read from, may be stdin */
  FILE *save;           /* an temporary file to save sequences read from stdin */
};

/* sortable sequence score record */
struct strand_t {
  SSEQ_T *sequence;     /* details of the sequence */
  double Pvalue;        /* p-value of product of p-values */
  int strand;           /* -1 neg. strand, 0 both/protein, +1 pos. strand */
  double *best_scores;  /* the best score for each motif */
  int *best_location;   /* the location of the best scores for each motif */
};

struct sseq_t {
  int db_index;         /* identifies the database */
  int index;            /* loading number of the sequence */
  long  fp;             /* file pointer to beginning of sequence rec. */
  long length;          /* length of sequence */
  double *comp;         /* actual sequence composition or NULL if  */
  BOOLEAN pv_allocated; /* was the pv variable allocated for this sequence */
  double **pv;          /* pvalue distribution for each motif */
  STRAND_T *pos_strand; /* forward strand */
  STRAND_T *neg_strand; /* reverse strand */
};
/*}}}*/

/************************************************************************/
/*
        error

        Prints an error message to stderr and sets a flag to exit on the next
        call to exit_on_error
*/
/************************************************************************/
static void error(char* format, ...) {
  va_list  argp;

  fprintf(stderr, "FATAL: ");
  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  fprintf(stderr, "\n");
  fflush(stderr);
  errored = TRUE;
}

/************************************************************************/
/*
        exit_on_error

        Checks if the function error has been called and quits if it has
*/
/************************************************************************/
static void exit_on_error() {
  if (errored) {
#ifdef DEBUG
    abort();
#else
    exit(1);
#endif
  }
}

/**********************************************************************/
/*
        calc_p_value

        Calculate the p-value of the product of the individual motif 
        p-values < observed

        Returns the combined p-value.
*/
/**********************************************************************/
static double calc_p_value(
  int sample_index,     /* index of sample */
  long length,          /* length of the sequence */
  STYPE stype,          /* handling of DNA strands */
  int nmotifs,          /* number of motifs */
  LO *los[],            /* array of pointers to log-odds matrices */
  double *best_score,   /* array of best score for each motif */
  int *best_loc,        /* position of best match for each motif */
  double range,         /* number of different score values */
  double **pv,          /* p-value tables for each motif */
  BOOLEAN sonly,        /* calculate p-value of observed spacings only */
  BOOLEAN lump,         /* combine spacings into one p-value */
  int norder           /* number of motifs in diag */
)
{
  int i;
  double pvalue;
  int nspaces;                  /* number of motif pairs spacings given for */
  double k_scores;              /* product of motif score p-values */
  double k_spacing;             /* product of spacing p-values */
  int smotifs;                  /* number of motifs that got scored */

  /* calculate the product of the motif score p-values */
  k_scores = 1.0;
  for (i=smotifs=0; i<nmotifs; i++) {
    int ws = los[i]->ws;        /* width of motif in sequence */
    int x = best_score[i];      /* observed EV (extreme value) */
    long n = length - ws + 1;   /* number of samples in EV */
    double p;                   /* p-value of EV */
    if (best_score[i] == LITTLE) continue;      /* motif wasn't scored */
    if (stype == Combine) n*=2; /* combining both DNA strands */
    EV(pv[i][x], n, p);         /* compute sequence p-value of x */
    k_scores *= p;              /* product of p-values */
    smotifs++;                  /* number of motifs scored */
  }

  /* multiply by p-values of motif spacings if provided (norder > 1) */
  k_spacing = 1.0;
  nspaces = 0;                  /* number of spaces p-valued */
  for (i=1; i<norder; i++) {    /* space to previous motif */
    if (space[i] >= 0) {        /* don't ignore this space */
      double p;                 /* to hold p-value of obs. spacing */
      int err;                  /* error in pos ith and i-1th motifs */
      int mi = order[i];        /* current motif */
      int mim1 = order[i-1];    /* previous motif */ 
      long npos;                /* number of positions for motif */
      int ws_avg;               /* average width of adjacent motifs */
      long mind;                /* minimum spacing that fits seq */
      long maxd;                /* maximum spacing that fits seq */

      /* get absolute error: diffence between observed and nominal spacing */
      err = abs(space[i] - (best_loc[mi] - best_loc[mim1]));

      /* get the average motif width and maximum number of placements */
      ws_avg = (int) (los[mim1]->ws + los[mi]->ws)/2.0; /* round down */
      npos = length - ws_avg + 1;

      mind = MAX(space[i] - err, ws_avg - length);
      maxd = MIN(space[i] + err, length - ws_avg);
      if (mind >= 0) {
        p = (maxd-mind+1) * (npos - (maxd+mind)/2.0);
      } else {
        p = -mind * (npos - (1-mind)/2.0) + (maxd+1) * (npos - maxd/2.0);
      }
      p /= npos * npos;
      if (p > 1.0 || p < 0.0) fprintf(stderr, 
        "\nerror in spacing p-value:seq_%d L=%8ld mind=%8ld E=%4d %-10g\n", 
         sample_index, length, mind, err, p);

      /* skip if no possible spacing */
      if (maxd < mind) {
        continue;
      }

      k_spacing *= p;           /* product of p-values */
      nspaces++;                /* number of spacings used */

    } /* don't ignore */
  } /* space to previous motif */

  /* finish the calculation */
  if (sonly) {                                  /* spacings only */
    pvalue = qfast(nspaces, k_spacing);
  } else if (lump && nspaces) {                 /* lump spacings */
    pvalue = qfast(nspaces, k_spacing);         /* spacing pvalue */
    pvalue = qfast(smotifs+1, k_scores*pvalue); /* spacing and motif p-value */
  } else {                                      /* spacings and motif scores */
    pvalue = qfast(smotifs+nspaces, k_scores*k_spacing);
  }

  return pvalue;
} /* calc_p_value */

/**********************************************************************/
/*
        best_frame

        Find the predominant frame of the matches to the motifs for
        a sequence.  The best frame is the frame with the minimum
        product of p-values.

        Returns an character a, b or c
*/
/**********************************************************************/
static char best_frame(
  int nmotifs,                          /* number of motifs */
  long length,                          /* length of sequence */
  TILING tiling                         /* tiling of sequence */
)
{
  int i, m;
  long j;
  double best_frame_p = 1;              /* p-value of best frame */
  int best = 0;                         /* best frame */

  /* find the product of p-values for each frame */
  for (i=0; i<3; i++) {                         /* frame */
    double best_p[MAXLO];                       /* best p-values for motifs */
    double prod_p = 1;                          /* product of p-values */

    /* clear array of best p-values for each motif */
    for (j=0; j<nmotifs; j++) best_p[j] = 1;    /* best p-value for motif j */

    /* find best p-value for each motif in this frame */
    for (j=i; j<length; j+=3) {                 /* position in sequence */
      if ((m = abs(tiling.hits[j])-1) >= 0) 
        best_p[m] = MIN(best_p[m], tiling.pvalues[j]);
    }

    /* form the product of p-values for this frame */
    for (j=0; j<nmotifs; j++) prod_p *= best_p[j];

    /* update best frame */
    if (prod_p < best_frame_p) {
      best_frame_p = prod_p;
      best = i;
    }
  } /* frame */
 
  switch (best) {
    case 0:
    return 'a';
    case 1:
    return 'b';
    case 2:
    return 'c';
    default:
    die("Illegal value for best.");
    return '\0';
  }
} /* best_frame */


/**********************************************************************/
/* 
 * best_scores
 * 
 * scores the sequence and pulls out the best scores for
 * each motif.
 * If there are existing scores in the best_score array it will consider
 * them in determining the best score.
 */
/**********************************************************************/
static void best_scores(
  BOOLEAN translate_dna,    /* database is DNA and motifs protein */
  LO *los[],            /* array of pointers to log-odds matrices */
  int nmotifs,          /* number of motifs */
  char *sequence,       /* sequence */
  long length,          /* length of the sequence */
  BOOLEAN rc_seq,       /* should the scores be found for the reverse complment */
  double **best_score,  /* OUT: the best scores for each motif, caller must deallocate */
  int **best_loc        /* OUT: the location for each best score, caller must deallocate */
  ) 
{
  SCORE **scores;
  double *score;            /* best score per motif for positive strand */
  int *loc, i, j, last_j;   /* best location per motif for positive strand */

  score = *best_score;
  loc = *best_loc;
  if (score == NULL) {
    Resize(score, nmotifs, double);
    for (i = 0; i < nmotifs; i++) score[i] = LITTLE;
  }
  if (loc == NULL) {
    Resize(loc, nmotifs, int);
    for (i = 0; i < nmotifs; i++) loc[i] = 0;
  }

  scores = score_it(translate_dna, rc_seq, los, nmotifs, sequence, length);
  

  for (i = 0; i < nmotifs; i++) { //motifs
    last_j = length - los[i]->ws;
    for (j = 0; j <= last_j; j++) {
      if (scores[i][j].score > score[i]) {
        score[i] = scores[i][j].score;
        loc[i] = j;
      }
    }
  }
  free_2array(scores, nmotifs);

  *best_score = score;
  *best_loc = loc;
}

/**********************************************************************/
/*
 * strand_destroy
 *
 * Destroys a strand score structure
 */
/**********************************************************************/
static void strand_destroy(STRAND_T *strand) {
    if (strand->best_scores) myfree(strand->best_scores);
    if (strand->best_location) myfree(strand->best_location);
    memset(strand, 0, sizeof(STRAND_T));
    myfree(strand);
}

/**********************************************************************/
/*
 * sseq_destroy
 *
 * Destroys a scored sequence
 * Note, also destroys connected
 * strand score structures.
 */
/**********************************************************************/
static void sseq_destroy(SSEQ_T *sseq, int nmotifs) {
  // deallocate connected strand objects
  if (sseq->neg_strand) strand_destroy(sseq->neg_strand);
  if (sseq->pos_strand) strand_destroy(sseq->pos_strand);
  // free up optional allocations
  if (sseq->comp) myfree(sseq->comp);
  //free expected allocations
  if (sseq->pv_allocated) {
    free_2array(sseq->pv, nmotifs);
  }
  //zero memory
  memset(sseq, 0, sizeof(SSEQ_T));
  //free struct
  myfree(sseq);
}

/**********************************************************************/
/*
 * add_strand
 *
 * Helper function to add a score for a strand (or both strands) to 
 * the list of scores
 */
/**********************************************************************/
static void add_strand(ARRAYLST_T *strand_scores, SSEQ_T *sseq, int strand, double pvalue, int nmotifs, double *best_score, int *best_loc) {
  STRAND_T *score = NULL;
  Resize(score, 1, STRAND_T);
  score->sequence = sseq;
  score->Pvalue = pvalue;
  score->strand = strand;
  if (sseq->comp) { // only needed if adjusting for sequence composition
    score->best_scores = mymalloc(nmotifs * sizeof(double));
    memcpy(score->best_scores, best_score, nmotifs * sizeof(double));
    score->best_location = mymalloc(nmotifs * sizeof(int));
    memcpy(score->best_location, best_loc, nmotifs * sizeof(int));
  } else {
    score->best_scores = NULL;
    score->best_location = NULL;
  }
  if (strand != -1) {
    sseq->pos_strand = score;
  } else {
    sseq->neg_strand = score;
  }
  arraylst_add(score, strand_scores);
}

/**********************************************************************/
/*
        get_scores

        Calculate the score and p-value for each sequence in the
        database. 

        Returns scores of sequences and sequence info in array lists
*/
/**********************************************************************/
static void get_scores(
  BOOLEAN dna,          /* database is DNA */
  STYPE stype,          /* handling of DNA strands */
  BOOLEAN translate_dna,    /* database is DNA and motifs protein */
  DATABASE_T *database, /* database */
  LO *los[],            /* array of pointers to log-odds matrices */
  int nmotifs,          /* number of log-odds matrices in los */
  int range,            /* set entries in matrices to [0..range] */
  double **pv,          /* p-value tables for each motif */
  BOOLEAN sonly,        /* calculate p-value of observed spacings only */
  BOOLEAN lump,         /* combine spacings into one p-value */
  BOOLEAN use_seq_comp, /* adjust E-value for actual sequence composition */
  double e_max,         /* maximum sequence E-value to print */
  int kmotifs,          /* number of known motifs*/
  MOTIF motif[],        /* data returned by read_motifs */
  BOOLEAN status,       /* show progress */
  int min_seqs,         /* lower bound on nseqs */
  int *nseqs,           /* total sequences in all databases */
  double *residues,     /* total number of residues in all seqs */
  ARRAYLST_T *sequences,/* IN/OUT: sequences saved because they have interesting scores */
  ARRAYLST_T *strand_scores /* IN/OUT: interesting strand scores */
)
{
  int i, j, imotif, sample_index;
  int est_seq; // estimated number of sequences
  long sample_len; // length of sample
  int file_pos;
  char *sample_name, *sample_seq, *sample_comment;
  
  double *best_score = NULL;            /* best score per motif for positive strand */
  int *best_loc = NULL;                 /* best location per motif for positive strand */

  STRAND_T *seqs = NULL;                /* sortable array of seq/score records*/
  double pvalue;                        /* expected number of sequences / total >= x */
  SCORE **scores;                       /* scores for each motif vs seq. pos. */
  int alen = los[0]->alen;              /* length of motif alphabet */
  BOOLEAN saved = FALSE;                /* sequence already saved */
  SSEQ_T *sseq;

  sample_index = 0;
  // loop through each sequence
  file_pos = (database->save != NULL) ? ftell(database->save) : ftell(database->file);
  while (read_sequence(database->file, &sample_name, &sample_comment, &sample_seq, &sample_len)) {
    sseq = (SSEQ_T*)mymalloc(sizeof(SSEQ_T));
    sseq->db_index = database->index;
    sseq->index = sample_index;
    sseq->length = sample_len;
    sseq->fp = file_pos;
    sseq->comp = use_seq_comp ? get_seq_comp(translate_dna, sample_seq, alen) : NULL;
    sseq->pv = pv;
    sseq->pv_allocated = FALSE;
    sseq->pos_strand = NULL;
    sseq->neg_strand = NULL;
    saved = FALSE;                                    /* new sequence */

    // update size of database 
    if (stype==Separate) { /* treat each DNA sequence as two */
      (*nseqs)+=2;                      /* number of sequences in database */
      *residues += 2*sample_len;        /* number of residues in database */
    } else if (stype!=Separate) {       /* treat each sequence as one */
      (*nseqs)++;                       /* number of sequences in database */
      *residues += sample_len;          /* number of residues in database */
    }
    database->sequence_count += 1;
    database->residue_count += sample_len;

    //clear the scores prior to use
    if (best_score) for (i = 0; i < nmotifs; ++i) best_score[i] = LITTLE;
    if (best_loc) for (i = 0; i < nmotifs; ++i) best_loc[i] = 0;

    est_seq = MAX(min_seqs, *nseqs);//get the current estimate sequence number

    //score the positive strand
    best_scores(translate_dna, los, nmotifs, sample_seq, sseq->length, FALSE, &best_score, &best_loc);
    if (stype == Combine) {//and the negative strand in combined mode
      best_scores(translate_dna, los, nmotifs, sample_seq, sseq->length, TRUE, &best_score, &best_loc);
    }

    // calculate the combined p-value for the positive (or maybe both) strand(s)
    pvalue = calc_p_value(sseq->index, sseq->length, stype, nmotifs, los, best_score, 
        best_loc, range, pv, sonly, lump, norder);
    
    if (pvalue * est_seq < e_max) {//save the pos strand
      saved = TRUE;
      add_strand(strand_scores, sseq, (stype == Combine || stype == Protein ? 0 : 1), pvalue, nmotifs, best_score, best_loc);
    }
    

    if (stype == Separate) {
      //clear the scores prior to use
      if (best_score) for (i = 0; i < nmotifs; ++i) best_score[i] = LITTLE;
      if (best_loc) for (i = 0; i < nmotifs; ++i) best_loc[i] = 0;
      //score negative strand
      best_scores(translate_dna, los, nmotifs, sample_seq, sseq->length, TRUE, &best_score, &best_loc);
      // calculate the combined p-value for the negative strand
      pvalue = calc_p_value(sseq->index, sseq->length, stype, nmotifs, los, best_score, 
          best_loc, range, pv, sonly, lump, norder);
      
      if (pvalue * est_seq < e_max) {//save the neg strand
        saved = TRUE;
        add_strand(strand_scores, sseq, -1, pvalue, nmotifs, best_score, best_loc);
      }
    }

    // save the sequence if a strand has an E-value that may be under threshold 
    if (saved) { 
      arraylst_add(sseq, sequences);
      // write the sequence to the save file if reading standard input
      if (database->file == stdin) {
        fprintf(database->save, ">%s %s\n%s\n", sample_name, sample_comment, sample_seq);
      }
    }

    // print progress report
    if (status && *nseqs % 100 == 0) 
      fprintf(stderr, "\rsequences: %6d ", *nseqs);

    // free up space
    if (!saved) sseq_destroy(sseq, nmotifs);
    myfree(sample_seq);
    myfree(sample_comment);
    myfree(sample_name);

    file_pos = (database->save != NULL) ? ftell(database->save) : ftell(database->file);
    sample_index++;
  } /* read_sequence */

  // free up space
  myfree(best_score);
  myfree(best_loc);

  /*
    recalculate the p-values based on actual sequence composition
  */
  for (i=0; use_seq_comp && i < arraylst_size(sequences); i++) {
    double **pv;
    sseq = arraylst_get(i, sequences);

    //skip if E-value too big, or no scores for this sequence
    if ((sseq->pos_strand == NULL || (sseq->pos_strand->Pvalue * (*nseqs) > e_max)) && 
        (sseq->neg_strand == NULL || (sseq->neg_strand->Pvalue * (*nseqs) > e_max))) continue;       

    // calculate p-values of all integer score values in range [0...w*MAST_RANGE] 
    pv = (double**)mymalloc(nmotifs * sizeof(double*));
    for (j=0; j<nmotifs; j++) {
      pv[j] = calc_pssm_cdf(los[j]->w, los[j]->alen, range, los[j]->logodds, sseq->comp);
    }

    // recalculate the combined p-value for this sequence using the actual
    // sequence composition as the background model
    if (sseq->pos_strand) {// for the positive strand
      STRAND_T *st = sseq->pos_strand;
      st->Pvalue = calc_p_value(sseq->index, sseq->length, stype, nmotifs, los, st->best_scores, 
          st->best_location, range, pv, sonly, lump, norder);
#ifdef DEBUG
      if (isnan(st->Pvalue)) {
        fprintf(stderr, "nan: i %d d-ind %d s-ind %d strand %d pvalue %f\n", i, sseq->db_index, 
            sseq->index, st->strand, st->Pvalue); 
        for (j=0; j<alen; j++) fprintf(stderr, "%c %f\n", unhash(j), sseq->comp[j]);
        abort();
      }
#endif
    }
    if (sseq->neg_strand) { // for the negative strand
      STRAND_T *st = sseq->neg_strand;
      st->Pvalue = calc_p_value(sseq->index, sseq->length, stype, nmotifs, los, st->best_scores, 
          st->best_location, range, pv, sonly, lump, norder);
#ifdef DEBUG
      if (isnan(st->Pvalue)) {
        fprintf(stderr, "nan: i %d d-ind %d s-ind %d strand %d pvalue %f\n", i, sseq->db_index, 
            sseq->index, st->strand, st->Pvalue); 
        for (j=0; j<alen; j++) fprintf(stderr, "%c %f\n", unhash(j), sseq->comp[j]);
        abort();
      }
#endif
    }
    /* 
      print progress report
    */
    if (status && i && i % 100 == 0) 
      fprintf(stderr, "\rrecalc p-value sequences: %6d ", i);

    /* 
      free pv space if E-value too big, otherwise
      save composition-based pv distribution
    */
    if ((sseq->pos_strand == NULL || (sseq->pos_strand->Pvalue * (*nseqs) > e_max)) &&
        (sseq->neg_strand == NULL || (sseq->neg_strand->Pvalue * (*nseqs) > e_max))) {
      free_2array(pv, nmotifs);
    } else {
      sseq->pv = pv;
      sseq->pv_allocated = TRUE;
    }
  } /* recalculate p-values */

} /* get_scores */


/**********************************************************************/
/*
        get_motifs

        Read in the log-odds matrices (motifs) from standard input, 
        remove unused motifs and interpret the ordering and spacing diagram.
        Exits if there is an error in the input.

        Returns number of motifs in the (compacted) los array.
        Updates los as well as globals order and space.
        Updates s2b.
*/
/**********************************************************************/
static int get_motifs(
  BOOLEAN translate_dna,            /* database is DNA and motifs protein */
  MOTIF_T *meme_io_motifs,      /* array of motifs loaded by meme-io */
  int meme_io_motif_num,        /* number of motifs loaded by meme-io */
  char *alphabet,               /* alphabet of logodds matrices */
  char *blast_alphabet,         /* corresponding BLAST alphabet */
  int *p[MAXASCII],             /* alphabet permutation/substitution matrix */
  BOOLEAN shuffle,              /* shuffle the motif columns */
  double min_ev,                /* minimum E-value of motifs to use */
  int umotifs,                  /* number motifs to be used */
  int mlist[MAXLO],             /* list of motifs given in -m option */
  BOOLEAN motifs[MAXLO],        /* list of motifs to output */
  LO *los[MAXLO],               /* array of logodds matrices */
  double *back,                 /* array of null letter probabilities */
  int range                     /* set logodds matrices to [0..range] */
)
{
  int i, imotif, cnt;
  int nmotifs;                  /* number of log-odds matrices in los */

  /* 
    read in log-odds matrices 
  */
  nmotifs = convert_2_log_odds(translate_dna, meme_io_motifs, meme_io_motif_num, alphabet, blast_alphabet,
    p, range, los, back, NULL);

  /*
    check that at least one motif was read
  */
  if (nmotifs == 0) {
    die("No scoring matrices found.\n");
  }

  /* 
    clear the list of flags showing motifs to use
  */
  for (i=0; i<MAXLO; i++) motifs[i] = FALSE;

  /* 
   set flags of motifs to use; check for valid motif numbers 
  */
  if (umotifs > 0) {                    /* using specified motifs */
    for (i=0; i<umotifs; i++) {
      int m = mlist[i];
      if (m < 1 || m > nmotifs) {
        die("Motif %d in -m option not in legal range: 1 to %d.\n",
          m, nmotifs); 
      }
      motifs[m-1] = TRUE;
    }
  } else {                              /* using motifs with E-value < min_ev */
    for (i=0; i<nmotifs; i++) {
      if (los[i]->ev < min_ev) { 
        motifs[i] = TRUE;
        umotifs++;
      }
    }
    if (!umotifs) {
      die("No input motifs pass the E-value (-mev) threshold %g.\n", min_ev);
    }
  }

  /*
    flag motifs with no information
  */
  for (i=0; i<nmotifs; i++) {
    if (motifs[i] && !los[i]->scale) {          /* motif has no info */
      motifs[i] = FALSE;
      umotifs--;
    }
  }

  /*
    check that valid motifs remain
  */
  if (!umotifs) {
    die("No scoring matrices contained any information.\n");
  }

  /* 
    Parse the ordering and spacing diagram.
    Variables norder, order and space are globals defined in diagram.h 
  */

  if (diagram == NULL) {
   norder = 0;                                  /* no diagram */
  } else {                                      /* have a diagram */
    BOOLEAN used[MAXLO];                        /* for checking use once */

    dptr = 0;                                   /* diag read pointer */
    norder = 0;                                 /* lnumber of motifs in diag */
    for (i=0; i<MAXLO; i++) space[i] = 0;       /* set all spacers to zero */ 
    if (yyparse()) exit(1);                     /* parse diagram; sets globals
                                                   norder, order, space */
    /* check that no unused motifs are in motif diagram */
    for (i=0; i<norder; i++) {
      int m = order[i];                         /* motif number (plus 1) */
      if (m < 1                                 /* number too small */
          || m > nmotifs                        /* number too large */
          || !motifs[m-1]                       /* motif not being used */
      ) {
        error("Unknown or unused motif %d given in motif diagram.\n", m);
      }
    }

    /* check that each motif only used once in diagram */
    for (i=0; i<MAXLO; i++) used[i] = FALSE; 
    for (i=0; i<norder; i++) {
      int m = order[i];
      if (used[m]) {
        error("Motif %d used more than once in motif diagram:\n  %s\n", m, diagram);
      } else {
        used[m] = TRUE;
      }
    }
    
    /* quit if there was a problem in the diagram */
    exit_on_error();

    /* change format for spacers to be relative to start of prior motif */
    space[0] = -1;                              /* ignore first spacer */
    for (i=1; i<norder; i++) {
      int m = order[i-1];                       /* previous motif in diagram  */
      space[i] += los[m-1]->w;                  /* internal format is -1 */
    }
  }
  
  /* 
    remove any motifs that are not to be used 
  */
  for (imotif=0, cnt=0; imotif<nmotifs; imotif++) {     /* motif */
    if (motifs[imotif]) {
      los[cnt++] = los[imotif];                 /* logodds matrix */
    }
  } /* motif */
  nmotifs = cnt;                                /* number of used motifs */

  /* 
    convert spacing diagram to internal motif names 
  */
  for (i=0; i<norder; i++) {
    for (imotif=0; imotif < nmotifs; imotif++) {
      if (order[i] == los[imotif]->imotif) break; /* found motif name */
    }
    order[i] = imotif;                          /* internal name for motif */
  }

  /* shuffle the columns of each motif matrix if requested */
  if (shuffle) shuffle_cols(los, nmotifs);

  /* compute the pairwise motif correlations (similarities) */
  motif_corr(nmotifs, los);

  /* return number of motifs remaining in (compacted) los array */
  return nmotifs;
} /* get_motifs */

/************************************************************************/
/*
        get_bad_motifs

        Get list of motifs to remove so that no pairs of motifs will
        have excessive correlation.

        Returns number of bad motifs and the list of bad motifs.
*/
/************************************************************************/
static int get_bad_motifs(
  double maxcorr,                       /* maximum allowed corellation */
  int nmotifs,                          /* number of motifs */ 
  LO *los[],                            /* array of logodds matrices */
  int bad_list[]                        /* list of highly correlated motifs */
)
{
  int i, j;
  int nbad = 0;                         /* number of highly correlated motifs */

  for (i=1; i<nmotifs; i++) {           /* from motif */
    for (j=0; j<i; j++) {               /* to motif */
      if (!los[j]->is_bad && los[i]->corr[j] > maxcorr) {
        bad_list[nbad++] = los[i]->imotif;
        los[i]->is_bad = TRUE;
        break; //no point in comparing this motif further
      }
    } /* to motif */
  } /* from motif */

  return nbad;
} /* get_bad_motifs */

/************************************************************************/
/*
        print_hit_list

        Create and print a list of non-overlapping hits.
*/
/************************************************************************/
void print_hit_list(
  FILE *mast_out,
  DATABASE_T *db,                       /* source database */
  LO *los[],                            /* array of pointers to lo matrices */
  int nmotifs,                          /* number motifs read */
  BOOLEAN dna,                          /* database is DNA */
  STYPE stype,                          /* handling of different strands */
  BOOLEAN translate_dna,                    /* database is DNA and motifs protein */
  BOOLEAN best_motifs,                  /* show only best motifs in diagrams */
  double **pv,                          /* p-value tables for each motif */
  double m_thresh,                      /* maximum motif p-value to print */
  double w_thresh,                      /* max. motif p-value for weak hits */
  BOOLEAN status                        // print sequence numbers to STDERR
) 
{
  long length;                  /* length of sample */
  char *sample_name;                    /* name of sample */
  char *sequence;                       /* sequence of sample */
  char *id;                             /* identifier text for sample */
  TILING tiling;                        /* tiling and diagram of sequence */
  int nseqs = 0;

  fprintf(mast_out,"# All non-overlapping hits in all sequences from \"%s\".\n", db->name);
  fprintf(mast_out,"# sequence_name motif hit_start hit_end score hit_p-value\n");

  while (read_sequence(db->file, &sample_name, &id, &sequence, &length)) {
    nseqs++;
    if (status && nseqs % 100 == 0) fprintf(stderr, "\rsequences: %6d ", nseqs);
    tiling = score_tile_diagram(mast_out, sequence, length, los, nmotifs, dna, stype, FALSE,
      translate_dna, best_motifs, FALSE, pv, m_thresh, w_thresh, FALSE, TRUE, sample_name);
    myfree(sample_name);
    myfree(sequence);
    myfree(id);
  } // read_sequence
  if (status) fprintf(stderr, "\n");
} // print_hit_list

/**********************************************************************/
/*
 * init_databases
 *
 * to cater for multiple databases, mast can read a tab delimited file with
 * one database record on each line.
 *
 * Returns array list of databases
 */
/**********************************************************************/
static ARRAYLST_T* init_databases(BOOLEAN is_list, BOOLEAN is_dna, 
    char *source, char *name, char *link) {
  int i;
  ARRAYLST_T *databases;
  DATABASE_T *db;
  time_t now;
  // check if multiple databases have been specified
  if (is_list) {
    //read the tab delimited file
    FILE *fp;
    BUF_T *buffer;
    int col, delim;
    char *src, *name, *link, *token;
    databases = arraylst_create();
    if (strcmp(source, "-") == 0) {
      fp = stdin;
    } else if (file_exists(source)) {
      fp = fopen(source, "r");
      if (!fp) die("Error reading file specified for database list \"%s\","
          " error was given as %s.\n", source, strerror(errno));
    } else {
      fp = NULL;//stop warnings about "uninitilized variables"...
      die("File specified for database list \"%s\" doesn't exist!\n", source);
    }
    buffer = buf_create(100);
    buf_flip(buffer);
    col = 0, src = NULL, name = NULL, link = NULL;
    //read the file token by token, only keep the first 3 tokens in a line
    while (!feof(fp) || buf_remaining(buffer)) {
      token = buf_fread_token(buffer, fp, buf_is_delim_letter, "\t\n\r", 
          FALSE, NULL, 0, NULL);
      delim = buf_getc(buffer);
      switch (col++) {
        case 0: src = token; break;
        case 1: name = token; break;
        case 2: link = token; break;
        default: myfree(token);
      }
      if (delim == -1 || delim == '\n' || delim == '\r') {
        buf_fread_consume(buffer, fp, buf_is_delim_letter, "\n\r", FALSE);
        if (!file_exists(src)) {
          error("File list contains a file \"%s\" "
              "that doesn't exist!\n", src);
        }
        if (name == NULL) copy_string(&name, src);
        db = mm_malloc(sizeof(DATABASE_T));
        memset(db, 0, sizeof(DATABASE_T));
        db->source = src;
        db->name = name;
        db->link = link;
        arraylst_add(db, databases);
        src = NULL, name = NULL, link = NULL;
        col = 0;
      }
    }
    if (src != NULL) myfree(src);
    if (name != NULL) myfree(name);
    if (link != NULL) myfree(link);
    arraylst_fit(databases);
    buf_destroy(buffer);
    fclose(fp);
    exit_on_error();
  } else {
    if (strcmp(source, "-") == 0) {
      if (name == NULL) name = "Sequences";
    } else {
      if (!file_exists(source)) die("File specified for database \"%s\""
          " doesn't exist!\n", source);
      if (name == NULL) name = source;
    }
    db = mm_malloc(sizeof(DATABASE_T));
    memset(db, 0, sizeof(DATABASE_T));
    copy_string(&(db->source), source);
    copy_string(&(db->name), name);
    copy_string(&(db->link), link);
    databases = arraylst_create_sized(1);
    arraylst_add(db, databases);
  }
  now = time(NULL);//for the default last mod date of databases
  for (i = 0; i < arraylst_size(databases); i++) {
    db = arraylst_get(i, databases);
    db->index = i;
    db->is_dna = is_dna;
    db->last_mod = now;
    db->sequence_count = 0;
    db->residue_count = 0;
    if (strcmp(db->source, "-") == 0) {
      db->file = stdin;
      db->save = tmpfile();
      if (db->save == NULL) die("Error creating temporary file for "
          "storing stdin, error given as: %s\n", strerror(errno));
    } else {
      db->file = fopen(db->source, "r");
      db->save = NULL;
      if (db->file == NULL) die("Error opening sequence file \"%s\", "
          "error given as: %s\n", db->source, strerror(errno));
#ifdef UNIX
      {//update the last mod date to the correct value
        struct stat stbuf;
        stat(db->source, &stbuf);
        db->last_mod = stbuf.st_mtime;
      }
#endif
    }
  }
  return databases;
}

/**********************************************************************/
/* 
 * score_asc
 *
 * Compares two strands to allow sorting in the array list.
 * Sorts primarily by combined p-value ascending and then breaks
 * ties by sorting on database, file position in database and then
 * on strand. This should mean that no two strands are considered equal.
 */
/**********************************************************************/
static int score_asc(const void *v1, const void *v2) {
  STRAND_T *hit1, *hit2;
  double diff;
  int db_diff;
  long fp_diff;
  hit1 = *((STRAND_T**)v1);
  hit2 = *((STRAND_T**)v2);
  //order by score
  diff = hit1->Pvalue - hit2->Pvalue;
  if (diff < 0) {
    return -1;
  } else if (diff > 0) {
    return 1;
  }
  //order by database
  db_diff = hit1->sequence->db_index - hit2->sequence->db_index;
  if (db_diff < 0) {
    return -1;
  } else if (db_diff > 0) {
    return 1;
  }
  //order by file position
  fp_diff = hit1->sequence->fp - hit2->sequence->fp;
  if (fp_diff < 0) {
    return -1;
  } else if (fp_diff > 0) {
    return 1;
  }
  //must be the same sequence
  //order by strand
  if (hit1->strand > hit2->strand) {
    return -1;
  } else if (hit1->strand < hit2->strand) {
    return 1;
  }
  //I don't think this is possible
  return 0;
}

/**********************************************************************/
/*
 * check_chunk
 *
 * checks a part of the hits array of size chunk starting at start
 * for motif hits. Checks if the hits found would overlap the end
 * of the chunk and makes a prediction on how many more chunks would
 * need to be checked.
 *
 * Returns 0 if there is nothing in this chunk, returns 1 if there is
 * something in this chunk but it does not overlap, returns n if there
 * is something in this chunk and at minimum it would overlap a further
 * n-1 chunks. Assumes that hits can not overlap the end of the sequence.
 */
/**********************************************************************/
static int check_chunk(int length, int *hits, LO **motifs, int start, int chunk) {
  int i, end, found;
  found = 0;
  end = MIN(length, start + chunk);
  for (i = start; i < end; ++i) {
    if (hits[i]) {
      found = 1;
      i += motifs[abs(hits[i]) - 1]->ws;
      if (i > end) {
        found = ((int)((i - end) / chunk)) + 2;
      }
      --i;//because of the loop inc
    }
  }
  return found;
}

/**********************************************************************/
/*
 * create_hit_match_patterns
 *
 * Creates the patterns showing which characters mast believes are good
 * matches and the translation from DNA to Protein to compare with the 
 * motif (only if translating dna).
 */
/**********************************************************************/
void create_hit_match_patterns(char *sequence, int pos, BOOLEAN neg_strand, LO *motif, BOOLEAN translate_dna, 
    char *matches, char *translation, int buflen) {
  char *s, *protalph;
  int i, inc, w, ws;
  double scale, offset;
  LOGODDS logodds;
  // double check that we can do what we want
  if (buflen <= motif->w) die("Buffer is too short to generate motif match patterns.\n");
  // protein alphabet (for translation)
  protalph = PROTEINB;
  // the width and width in sequence of the motif
  w = motif->w;
  ws = motif->ws; //should be w * inc
  // the motif log odds matrix
  logodds = motif->logodds;
  // the details required for converting log odds scores
  scale = motif->scale;
  offset = motif->offset;
  // increment in sequence for each entry in motif
  inc = (translate_dna ? 3 : 1);
  // for each entry in the motif, calculate if the sequence is a good match
  // and if the underlying sequence is a protein to the nucleotide motif then
  // translate the sequence into proteins
  for (i = 0, s = sequence+pos; i < w; i++, s += inc) {
    int cx = chash(translate_dna, neg_strand, s);
    double score = neg_strand ? logodds(w-i-1, cx) : logodds(i, cx); 
    matches[i] = (scaled_to_bit(score, 1, scale, offset) > 0) ? '+' : ' ';
    if (translate_dna) translation[i] = protalph[cx];
  }
  matches[w] = '\0';
  if (translate_dna) translation[w] = '\0';
}

/**********************************************************************/
/*
 * output_mast_xml_hit
 *
 * Outputs a hit for the current position if one exists.
 */
/**********************************************************************/
static inline void output_mast_xml_hit(
    FILE *out, char *sequence, int pos, TILING tiling, LO **motifs,
    char *matches_buffer, char *translation_buffer, 
    int buffer_len, BOOLEAN translate_dna, STYPE stype, int *last_gap_start) {
  int hit, gap;
  LO *motif;
  BOOLEAN neg_strand;
  double pvalue;
  hit = tiling.hits[pos];
  //check that a hit exists at that location
  if (!hit) return; //no motif, return
  gap = pos - *last_gap_start;
  motif = motifs[abs(hit) - 1];
  neg_strand = (hit < 0);
  pvalue = tiling.pvalues[pos];
  create_hit_match_patterns(sequence, pos, neg_strand, motif, translate_dna, matches_buffer, translation_buffer, buffer_len);
  fprintf(out, "\t\t\t\t<hit pos=\"%d\" gap=\"%d\" motif=\"motif_%d\" pvalue=\"%.1e\"", pos + 1, gap, motif->imotif, pvalue);
  if (stype != Protein && stype != Norc) {
    fprintf(out, " strand=\"%s\"", (!neg_strand ? "forward" : "reverse"));
  }
  fprintf(out, " match=\"%s\"", matches_buffer);
  if (translate_dna) fprintf(out, " translation=\"%s\"", translation_buffer);
  fprintf(out, "/>\n");
  *last_gap_start = pos + motif->ws;
}

/**********************************************************************/
/*
 * output_mast_xml
 *
 * Outputs mast xml results to out.
 */
/**********************************************************************/
void output_mast_xml(FILE *out, BOOLEAN translate_dna, 
    char *alphabet, BOOLEAN use_seq_comp, char *bg_file, double *back,
    BOOLEAN removed_correlated, double max_correlation, 
    int nmotifs, int nbad, LO **good_motifs, LO **all_motifs, 
    char *motifs_source, char *motifs_name, time_t motifs_last_mod,
    int nos_norder, int *nos_order, int *nos_space,
    int nseqs, ARRAYLST_T *databases, ARRAYLST_T *sequence_scores, 
    int argc, char **argv, STYPE strand_handling, double max_seq_evalue, 
    BOOLEAN adj_hit_pvalue, double max_hit_pvalue, double max_weak_pvalue) {
  int i, j, last_gap_start1, last_gap_start2, seg_start, seg_width, total_motifs, lpos, longest_motif_len;
  char *strand_handling_str, *bg_source_str, *letter, *matches, *translation;
  char xml_convert_buffer[500];
  time_t now;
  double cycles;
  fprintf(out, "<?xml version='1.0' encoding='UTF-8' standalone='yes'?>\n");
  fprintf(out, "%s", MAST_DTD);
  fprintf(out, "<mast version=\"" VERSION "\" release=\"" ARCHIVE_DATE "\">\n");
  
  //print the model first as otherwise anything trying to scan through this will have to make multiple passes
  //note that the runtime is excluded from the model because it has to be printed last to have any useful meaning
  fprintf(out, "\t<model>\n");
  // command line
  fprintf(out, "\t\t<command_line>mast");
  for (i = 1; i < argc; ++i) fprintf(out, " %s", replace_xml_chars(argv[i], xml_convert_buffer, 500, FALSE, NULL));
  fprintf(out, "</command_line>\n");
  // MOTIF RELATED MODEL PARAMETERS
  // max correlation
  fprintf(out, "\t\t<max_correlation>%.2f</max_correlation>\n", max_correlation);
  //remove correlated
  fprintf(out, "\t\t<remove_correlated value=\"%s\"/>\n", (removed_correlated ? "y" : "n"));
  // STRAND RELATED MODEL PARAMETERS
  // strand scoring
  switch (strand_handling) {
    case Combine:   strand_handling_str = "combine";   break;
    case Separate:  strand_handling_str = "separate";  break;
    case Norc:      strand_handling_str = "norc";      break;
    case Protein:   strand_handling_str = "protein";   break;
    default:        strand_handling_str = NULL;
      die("Impossible value for strand_handling.\n");
  }
  fprintf(out, "\t\t<strand_handling value=\"%s\"/>\n", strand_handling_str);
  fprintf(out, "\t\t<translate_dna value=\"%s\"/>\n", (translate_dna ? "y" : "n"));
  // maximum sequence evalue
  fprintf(out, "\t\t<max_seq_evalue>%.2g</max_seq_evalue>\n", max_seq_evalue);
  // HIT RELATED MODEL PARAMETERS
  // adjust hit pvalue to be relative to the sequence length (convert it into an evalue?)
  fprintf(out, "\t\t<adj_hit_pvalue value=\"%s\"/>\n", (adj_hit_pvalue ? "y" : "n"));
  // the threshold for considering a motif match to be a true hit.
  fprintf(out, "\t\t<max_hit_pvalue>%g</max_hit_pvalue>\n", max_hit_pvalue);
  // the threshold for considering a motif match to be a weak hit.
  fprintf(out, "\t\t<max_weak_pvalue>%g</max_weak_pvalue>\n", max_weak_pvalue);
  //output host and runtime information
  fprintf(out, "\t\t<host>%s</host>\n", hostname());
  now = time(NULL);
  fprintf(out, "\t\t<when>%s</when>\n", strtok(ctime(&now),"\n"));
  fprintf(out, "\t</model>\n");
  
  //get the alphabet background source
  if (use_seq_comp) {
    bg_source_str = "sequence";
  } else if (bg_file != NULL) {
    bg_source_str = "file";
  } else {
    bg_source_str = "preset";
  }
  //note that we should never get here if there are no good motifs...
  fprintf(out, "\t<alphabet type=\"%s\" bg_source=\"%s\"", (good_motifs[0]->dna ? "nucleotide" : "amino-acid"), bg_source_str);
  if (bg_file) fprintf(out, " bg_file=\"%s\"", replace_xml_chars(bg_file, xml_convert_buffer, 500, TRUE, NULL));
  fprintf(out, ">\n");
  for (letter = alphabet, i = 0; *letter != '\0'; letter++, i++) {
    fprintf(out, "\t\t<letter symbol=\"%c\"", *letter);
    if (!use_seq_comp) {
      lpos = translate_dna ? protbhash(alphabet[i]) : hash(alphabet[i]);
      fprintf(out," bg_value=\"%.3f\"", back[lpos]);
    }
    fprintf(out, "/>\n");
  }

  fprintf(out, "\t</alphabet>\n");
  fprintf(out, "\t<motifs source=\"%s\" ", replace_xml_chars(motifs_source, xml_convert_buffer, 500, TRUE, NULL));
  fprintf(out, "name=\"%s\" last_mod_date=\"%s\">\n", replace_xml_chars(motifs_name, xml_convert_buffer, 500, TRUE, NULL), 
      strtok(ctime(&motifs_last_mod), "\n"));
  total_motifs = nmotifs + (removed_correlated ? nbad : 0);
  for (i = 0; i < total_motifs; i++) {
    LO *motif = all_motifs[i];
    fprintf(out, "\t\t<motif id=\"motif_%d\" num=\"%d\" name=\"%s\" width=\"%d\""
        " best_f=\"%s\"", motif->imotif, motif->imotif, replace_xml_chars(motif->meme_name,
          xml_convert_buffer, 500, TRUE, NULL), motif->w, motif->best_seq);
    //output reverse complement of motif if nucleotide motif
    if (motif->dna) fprintf(out, " best_r=\"%s\"", motif->best_icseq);
    if (motif->is_bad) fprintf(out, " bad=\"y\"");
    fprintf(out, "/>\n");
  }
  for (i = 1; i < total_motifs; ++i) {
    LO *motifb = all_motifs[i];
    for (j = 0; j < i; ++j) {
      LO *motifa = all_motifs[j];
      fprintf(out, "\t\t<correlation motif_a=\"motif_%d\" motif_b=\"motif_%d\""
          " value=\"%.2f\"/>\n", motifa->imotif, motifb->imotif, motifb->corr[j]);
    }
  }
  if (nos_norder) {
    int nos_length, pos, gap;
    LO *prev_motif;
    //it appears from reading get_motifs that the spacer before the first 
    //motif is ignored and the spacer after the last motif isn't loaded
    //so the length will have to be from the first motif to the end motif. 
    //Additionally it appears that the length of the gaps between the motifs 
    //is assumed to be in the same alphabet as the motifs.
    for (i = 1, nos_length = 0; i < nos_norder; ++i) {
      nos_length += nos_space[i];//spacer plus width of previous motif
    }
    fprintf(out, "\t\t<nos length=\"%d\">\n", nos_length);
    for (i = 0, pos = 1, prev_motif = NULL, gap = 0; i < nos_norder; ++i) {
      if (i != 0) {
        pos += nos_space[i];
        gap = nos_space[i] - prev_motif->w;//w is not a bug, this diagram is not scaled to the units of the sequence
      }
      fprintf(out, "\t\t\t<expect motif=\"motif_%d\" pos=\"%d\" gap=\"%d\" />\n", nos_order[i], pos, gap);
      prev_motif = all_motifs[nos_order[i]-1];
    }
    fprintf(out, "\t\t</nos>\n");
  }
  fprintf(out, "\t</motifs>\n");
  fprintf(out, "\t<sequences>\n");
  //calculate longest motif used in scoring so I know what size buffers to allocate
  longest_motif_len = 0;
  for (i = 0; i < nmotifs; ++i) {
    LO *motif = good_motifs[i];
    if (motif->w > longest_motif_len) longest_motif_len = motif->w;
  }
  matches = mymalloc(sizeof(char) * (longest_motif_len + 1));
  translation = mymalloc(sizeof(char) * (longest_motif_len + 1));
  // output each database
  for (i = 0; i < arraylst_size(databases); ++i) {
    DATABASE_T *db;
    db = (DATABASE_T*)arraylst_get(i, databases);
    fprintf(out, "\t\t<database id=\"db_%d\" num=\"%d\" source=\"%s\"", i+1, i+1, replace_xml_chars(db->source, 
          xml_convert_buffer, 500, TRUE, NULL));
    //note that ctime puts a newline on the end of the time string, so I use strtok to null it
    fprintf(out, " name=\"%s\" last_mod_date=\"%s\" seq_count=\"%d\""
        " residue_count=\"%ld\" type=\"%s\"", replace_xml_chars(db->name, 
          xml_convert_buffer, 500, TRUE, NULL), strtok(ctime(&(db->last_mod)), "\n"), 
        db->sequence_count, db->residue_count, 
        (db->is_dna ? "nucleotide" : "amino-acid"));
    if (db->link) fprintf(out, " link=\"%s\"", replace_xml_chars(db->link, 
          xml_convert_buffer, 500, TRUE, NULL));
    fprintf(out, "/>\n");
  }
  //output each high scoring sequence
  for (i = 0; i < arraylst_size(sequence_scores); i++) {
    FILE *seq_file;
    char *sample_name, *sample_comment, *sample_sequence;
    long sample_length;
    TILING score_tiling, other_tiling;
    SSEQ_T *sequence;
    STRAND_T *score, *other;
    DATABASE_T *db;
    BOOLEAN neg_strand;

    //get the next score and its sequence.
    sequence = ((score = (STRAND_T*)arraylst_get(i, sequence_scores))->sequence);

    //check that this sequence hasn't already been printed
    if (strand_handling == Separate) {
      neg_strand = (score->strand == -1);
      // if there was the possibility of two strands then get both
      // and check that the sequence hasn't been printed already.
      if (!neg_strand) other = sequence->neg_strand;
      else other = sequence->pos_strand;
      if (other) {
        // check if this sequence has already been printed
        if (other->Pvalue < score->Pvalue || (other->Pvalue == score->Pvalue && score->strand == -1)) 
          continue;
        // check that the evalue is acceptable
        if (other->Pvalue * nseqs > max_seq_evalue)
          other = NULL;
      }
    } else {
      neg_strand = FALSE;
      other = NULL;
    }

    //get the database associated with the sequence
    db = (DATABASE_T*)arraylst_get(sequence->db_index, databases);
    //get the file that the sequence should be loaded from
    seq_file = (db->save ? db->save : db->file);

    // reload the sequence from file
    fseek(seq_file, sequence->fp, SEEK_SET);
    read_sequence(seq_file, &sample_name, &sample_comment, &sample_sequence, &sample_length);

    fprintf(out, "\t\t<sequence id=\"seq_%d_%d\" db=\"db_%d\" num=\"%d\" name=\"%s\" ", sequence->db_index+1, 
        sequence->index+1, sequence->db_index+1, sequence->index+1, replace_xml_chars(sample_name, xml_convert_buffer, 500, TRUE, NULL));

    fprintf(out, "comment=\"%s\" length=\"%ld\">\n", replace_xml_chars(sample_comment, xml_convert_buffer, 500, TRUE, NULL), 
      sample_length);

    // score, tile and diagram the sequence with each of the motifs 
    score_tiling = score_tile_diagram(NULL, sample_sequence, sample_length, 
        good_motifs, nmotifs, db->is_dna, strand_handling, neg_strand, translate_dna, 
        FALSE, FALSE, sequence->pv, max_hit_pvalue, max_weak_pvalue, adj_hit_pvalue, FALSE, sample_name);

    fprintf(out, "\t\t\t<score strand=\"%s\" combined_pvalue=\"%.2e\" evalue=\"%.2g\"", 
        (score->strand == 0 ? "both" : 
         (score->strand == 1 ? "forward" : "reverse")), 
        score->Pvalue, score->Pvalue * nseqs);
    if (translate_dna) fprintf(out, " frame=\"%c\"", best_frame(nmotifs, sample_length, score_tiling));
    fprintf(out, "/>\n");
    if (strand_handling == Separate) {
      // score, tile and diagram the sequence with each of the motifs 
      other_tiling = score_tile_diagram(NULL, sample_sequence, sample_length, 
          good_motifs, nmotifs, db->is_dna, strand_handling, !neg_strand, translate_dna, 
          FALSE, FALSE, sequence->pv, max_hit_pvalue, max_weak_pvalue, adj_hit_pvalue, FALSE, sample_name);
    } else { // avoid compilier complaints
      memset(&other_tiling, 0, sizeof(TILING));
    }

    if (other) {
      fprintf(out, "\t\t\t<score strand=\"%s\" combined_pvalue=\"%.2e\" evalue=\"%.2g\"", 
          (other->strand == 1 ? "forward" : "reverse"), other->Pvalue, other->Pvalue * nseqs);
      if (translate_dna) fprintf(out, " frame=\"%c\"", best_frame(nmotifs, sample_length, other_tiling));
      fprintf(out, "/>\n");
    } 
    //now output sequence segments in multiplies of SEG_CHUNK skipping segments that don't have hits
    seg_start = 0;
    last_gap_start1 = 0;
    last_gap_start2 = 0;
    while (seg_start < sample_length) {
      {
        int chunk, chunks, primary, secondary, start;
        for (chunk = 0, chunks = 0, start = seg_start; chunk == 0 || chunk < chunks; ++chunk, start += SEG_CHUNK) {
          primary = check_chunk(sample_length, score_tiling.hits, all_motifs, start, SEG_CHUNK);
          if (strand_handling == Separate) {
            secondary = check_chunk(sample_length, other_tiling.hits, all_motifs, start, SEG_CHUNK);
            chunks = MAX(chunks, chunk + MAX(primary, secondary));
          } else {
            chunks = MAX(chunks, chunk + primary);
          }
        }
        if (chunks) {
          seg_width = MIN(sample_length - seg_start, SEG_CHUNK * chunks);
        } else {
          seg_width = -SEG_CHUNK;
        }
      }
      if (seg_width > 0) {
        int offset, pos, seg_end;
        fprintf(out, "\t\t\t<seg start=\"%d\">\n", seg_start + 1);
        fprintf(out, "\t\t\t\t<data>\n");
        for (offset = 0; offset < seg_width; offset += SEG_CHUNK) {
          fprintf(out, SEG_FORMAT, sample_sequence+(seg_start+offset));
        }
        fprintf(out, "\t\t\t\t</data>\n");
        pos = seg_start;
        seg_end = seg_start + seg_width;
        for (pos = seg_start; pos < seg_end; ++pos) {
          output_mast_xml_hit(out, sample_sequence, pos, score_tiling, all_motifs, 
              matches, translation, longest_motif_len + 1, translate_dna, strand_handling, &last_gap_start1);
          if (strand_handling == Separate)
            output_mast_xml_hit(out, sample_sequence, pos, other_tiling, all_motifs, 
                matches, translation, longest_motif_len + 1, translate_dna, strand_handling, &last_gap_start2);
        }
        fprintf(out, "\t\t\t</seg>\n");
      }
      seg_start += abs(seg_width);
    }

    //clean up TILING objects, which have allocations in them.
    myfree(score_tiling.hits);
    myfree(score_tiling.pvalues);
    myfree(score_tiling.svalues);
    myfree(score_tiling.diagram);
    if (strand_handling == Separate) {
      myfree(other_tiling.hits);
      myfree(other_tiling.pvalues);
      myfree(other_tiling.svalues);
      myfree(other_tiling.diagram);
    }
    //clean up
    myfree(sample_name);
    myfree(sample_comment);
    myfree(sample_sequence);
    fprintf(out, "\t\t</sequence>\n");
  }
  myfree(matches);
  myfree(translation);
  fprintf(out, "\t</sequences>\n");
  //the run time has to be printed last, or it wouldn't be the run time
  //as it has to be in the xml then it can't include the time to create the html
  cycles = myclock();
  fprintf(out, "\t<runtime cycles=\"%.0f\" seconds=\"%.3f\"/>\n", cycles, cycles / CLOCKS_PER_SEC);
  fprintf(out, "</mast>\n");
}

/**********************************************************************/
/*
  main routine
*/
/**********************************************************************/
extern int main (int argc, char *argv[])
{
  /* declarations */
  /*{{{*/
  /* defaults */

  char *out_dir = default_out_dir;      /* the output directory */
  BOOLEAN clobber = TRUE;               /* should output directory be overwritten */
  char *xml_path = NULL;                /* the xml output path  */
  char *alphabet = "";                  /* alphabet as loaded from motifs file */
  char *mo_source = "";                 /* the file path or "-" indicating
                                           stdin. The file will be in a meme
                                           motif format as supported by 
                                           meme-io */
  char *db_source = "";                 /* the file path or "-" indicating 
                                           stdin. The file will be either a
                                           fasta file or if db_list is true
                                           then it will be a list of fasta
                                           files.
                                           */
  BOOLEAN db_list = FALSE;              /* when true, treat the db_source as 
                                           a file containing a list of db 
                                           files */
  ARRAYLST_T *databases;                /* the list of seq databases */
  char *mo_name = NULL;                 /* when not null, use to replace 
                                           motif file name in results */
  time_t mo_last_mod;                   /* last modification date of motifs */
  char *db_name = NULL;                 /* when not null, use to replace 
                                           sequence file name in results */
  char *db_link = NULL;                 /* when not null, use to create links 
                                           for searching sequence names in 
                                           results */
  BOOLEAN weak = FALSE;                 /* print weak hits */
  int first_motifs = 0;                 /* read only the first x motifs */
  double e_max = EXPECT;                /* maximum sequence p-value to output */
  double m_thresh = 1e-4;               /* maximum motif p-value to output */
  double w_thresh;                      /* maximum motif p-value--weak hits 
                                           default: m_thresh * 10 */
  double min_ev = BIG;                  /* minimum E-value of motifs to use */
  int min_seqs = 0;                     /* lower bound on total sequence 
                                           count, allows for earlier filtering
                                           of results with a large evalue */
  BOOLEAN adj_hit_pvalue = FALSE;            /* use sequence p-value for m_thresh */
  BOOLEAN shuffle = FALSE;              /* don't shuffle the motif columns */
  BOOLEAN sonly = FALSE;                /* print product of spacing p-values */
  BOOLEAN lump = FALSE;                 /* combine spacings into one p-value */
  BOOLEAN status = TRUE;                /* show progress */
  BOOLEAN html = TRUE;                  /* generate html output */
  BOOLEAN mast2txt = TRUE;              /* run mast2txt on the xml output */
  STYPE stype = Combine;                /* handling of DNA strands */
  BOOLEAN translate_dna = FALSE;            /* don't translate DNA sequences */
  BOOLEAN best_motifs = FALSE;          /* only print best motif in diagrams */
  BOOLEAN use_seq_comp = FALSE;         /* adjust E-values for actual compos. */
  BOOLEAN rem_corr = FALSE;             /* remove overly correlated motifs */ 
  BOOLEAN hit_list = FALSE;             /* print block diagrams */ 
  
  /* locals */
  int i, j, k;
  FILE *mast_out = NULL;                /* the output file */
  char *bfile = NULL;                   /* no background file */
  int nseqs = 0;                        /* total sequences in database */
  double residues = 0;                  /* total number of residues in seqs */
  int n_hits = 0;                       /* sequences less than e_max */
  BOOLEAN dna;                          /* true if sequences are DNAB */

  /* motifs and logodds matrices */
  LO **all_motif_los;                   /* all motifs, not just the good ones */
  int nmotifs;                          /* number motifs read */ 
  int umotifs=0;                        /* number motifs to be used */
  int kmotifs=0;                        /* number of known motifs */
  LO *good_motif_los[MAXLO];            /* array of logodds matrices for good motifs */
  BOOLEAN motifs[MAXLO];                /* list of motifs to output */
  int bad_list[MAXLO];                  /* list of highly correlated motifs */
  int nbad;                             /* number of bad (correlated) motifs */
  int mlist[MAXLO];                     /* list of motifs given in -m option */
  MOTIF motif[NMOTIFS];                 /* data returned by read_motifs */
  DATASET dummy;                        /* dummy arg for read_motifs */
  char *blast_alphabet;                 /* corresponding BLAST alphabet */
  int *perm[MAXASCII];                  /* permutation/substitution matrix */

  /* things for normalization */
  double **pv = NULL;           /* p-value tables for each motif */
  double *back;                 /* background letter probability distribution */
  /*}}}*/
#ifdef MALLOC_DEBUG
  int malloc_debug(int);
  malloc_debug(2);
#endif

  (void) myclock();                     /* record CPU time */

  argv[0] = "mast";

  /* get the command line arguments */
  /*{{{*/
  i = 1;
  DO_STANDARD_COMMAND_LINE(1,
    NON_SWITCH(1, <motif file> <sequence file> [options] \n,
      if (i == 1) mo_source = _OPTION_;
      else if (i == 2) db_source = _OPTION_;
      else COMMAND_LINE_ERROR;
      ++i;
    );
    // as stdin can be used it is indicated by a dash, which is picked up as an empty option
    FLAG_OPTN(2, , , 
      if (i == 1) mo_source = "-";
      else if (i == 2) db_source = "-";
      else COMMAND_LINE_ERROR;
      ++i;
      );
    //OPTIONAL INPUT
    DATA_OPTN(1, bfile, <bf>, \tread background frequencies from <bf>,
      bfile = _OPTION_);
    FLAG_OPTN(1, dblist, \tthe file specified as database contains a list of databases,
      db_list = TRUE);
    //OUTPUT OPTIONS
    DATA_OPTN(1, o, <dir>, \tdirectory to output mast results, out_dir = _OPTION_; clobber = FALSE);
    DATA_OPTN(1, oc, <dir>, \tdirectory to output mast results with overwriting allowed, out_dir = _OPTION_; clobber = TRUE);
    FLAG_OPTN(1, hit_list, \tprint only a list of non-overlapping hits to stdout, hit_list = TRUE);
    //WHICH MOTIFS TO USE
    FLAG_OPTN(1, remcorr, \tremove highly correlated motifs from query,
      rem_corr = TRUE);
    DATA_OPTN(1, m, <m>,+\tuse only motif(s) number <m> (overrides -mev), 
      mlist[umotifs++] = atoi(_OPTION_));
    DATA_OPTN(1, c, <count>,
      \tonly use the first <count> motifs or all motifs when <count> is zero (default: 0),
      first_motifs = atoi(_OPTION_));
    DATA_OPTN(1, mev, <mev>,+\tuse only motifs with E-values less than <mev>, 
      min_ev = atof(_OPTION_));
    DATA_OPTN(1, diag, <diag>, \tnominal order and spacing of motifs,
      diagram = _OPTION_);
    //DNA ONLY OPTIONS
    FLAG_OPTN(1, norc, \t\tdo not score reverse complement DNA strand, 
      stype = Norc);
    FLAG_OPTN(1, sep,
      \t\tscore reverse complement DNA strand as a separate \n\t\t\tsequence, 
      stype = Separate);
    FLAG_OPTN(1, dna, \t\ttranslate DNA sequences to protein, translate_dna = TRUE);
    FLAG_OPTN(1, comp, \t\tadjust p-values and E-values for sequence composition, 
      use_seq_comp = TRUE);
    //WHICH RESULTS TO PRINT
    DATA_OPTN(1, ev, <ev>,
      \tprint results for sequences with E-value < <ev>\n\t\t\t(default: 10),
      e_max = atof(_OPTION_)); 
    //APPEARANCE OF BLOCK DIAGRAMS
    DATA_OPTN(1, mt, <mt>, 
      \tshow motif matches with p-value < mt (default: 0.0001), 
      m_thresh = atof(_OPTION_)); 
    FLAG_OPTN(1, w, 
      \t\tshow weak matches (mt<p-value<mt*10) in angle brackets, weak = TRUE);
    FLAG_OPTN(1, best, \t\tinclude only the best motif in diagrams (hit_list mode only),
      best_motifs = TRUE);
    FLAG_OPTN(1, seqp, \t\tuse SEQUENCE p-values for motif thresholds\n\t\t\t(default: use POSITION p-values), adj_hit_pvalue = TRUE);
    //MISCELLANEOUS
    DATA_OPTN(1, mf, <mf>, \tin results use <mf> as motif file name, mo_name = _OPTION_);
    DATA_OPTN(1, df, <df>, \tin results use <df> as database name (ignored when -dblist specified), db_name = _OPTION_);
    DATA_OPTN(1, dl, <dl>, \tin results use <dl> as link to search sequence names (ignored when -dblist specified), db_link = _OPTION_);
    DATA_OPTN(1, minseqs, <ms>, \tlower bound on number of sequences in db,
      min_seqs = atoi(_OPTION_));
    FLAG_OPTN(1, nostatus, \tdo not print progress report, status = FALSE);
    FLAG_OPTN(1, notext, \tdo not generate text output, mast2txt = FALSE);
    FLAG_OPTN(1, nohtml, \tdo not generate html output, html = FALSE);
    // DEBUG AND EXPERIMENTAL OPTIONS
    FLAG_OPTN(EXP, shuffle, \tshuffle columns of motifs, shuffle = TRUE);
    FLAG_OPTN(EXP, sonly, \tuse only spacing p-values in product, sonly = TRUE);
    FLAG_OPTN(EXP, lump, \t\tcombine spacings into one p-value, lump = TRUE);
  );
  /*}}}*/
  // check that the two required parameters have been specified
  if (i != 3) {
    die("The meme file and sequences file are both required.\n");
  }
  
  if (strcmp(db_source, mo_source) == 0) {
    die("The sequence database and motif file have the same source.\n");
  }

  if (mo_name == NULL) {
    if (strcmp(mo_source, "-") == 0) {
      mo_name = "Motifs";
    } else {
      mo_name = mo_source;
    }
  }

  if (!hit_list && out_dir != NULL) {
    // configure output
    if (create_output_directory(out_dir, clobber, status)) {
      // Failed to create output directory.
      exit(1);
    }
    // create "<dir>/XML_FILENAME" and open it for writing
    xml_path = make_path_to_file(out_dir, XML_FILENAME);
    mast_out = fopen(xml_path, "w");
  } else {
    mast_out = stdout;
  }

  /* set e_max to BIG if printing hit-list */
  if (hit_list) e_max = BIG;

  /* set maximum p-value for weak hits */
  if (weak) {
    w_thresh = m_thresh * 10.0;                 /* 10 times less likely */
  } else {
    w_thresh = m_thresh;                        /* no weak threshold */
  }

  /* load motifs using meme-io, this is a bit of a hack since 
   * I don't want to break anything and so I'm making a minimum
   * of changes...*/
  {
    MOTIF_T meme_io_motifs[MAX_MOTIFS]; 
    int num_motifs;
    ARRAY_T *background;
    BOOLEAN_T has_reverse_strand;

    mo_last_mod = time(NULL);
#ifdef UNIX
    if (strcmp(mo_source, "-") != 0) {
      //update the last mod date to the correct value
      struct stat stbuf;
      stat(mo_source, &stbuf);
      mo_last_mod = stbuf.st_mtime;
    }
#endif

    read_meme_file(mo_source,
       bfile, // bg file name
       0, //pseudocounts
       REQUIRE_PSSM,
       &num_motifs,
       meme_io_motifs, 
       NULL,//motif occurrences
       &has_reverse_strand,
       &background);

    free_array(background);

    // note that set_alphabet is called by read_meme_file
    alphabet = get_alphabet(FALSE); 

    //limit the number of motifs
    if (first_motifs > 0 && first_motifs < num_motifs) {
      num_motifs = first_motifs;
    }
    // initialize the background frequencies and alphabets 
    back = init_mast_background(bfile, alphabet, stype, translate_dna, 
      &blast_alphabet, perm);

    // get the motifs and parse the ordering and spacing diagram
    nmotifs = get_motifs(translate_dna, meme_io_motifs, num_motifs, alphabet, 
        blast_alphabet, perm, shuffle, min_ev, umotifs, mlist, motifs, 
        good_motif_los, back, MAST_RANGE);
  }
  /* get list of motifs that should be removed */
  nbad = get_bad_motifs(MAXCORR, nmotifs, good_motif_los, bad_list);

  /* make a copy of all motifs in file order, used for printing */
  all_motif_los = (LO**)mymalloc(sizeof(LO*) * nmotifs);
  for (i = 0; i < nmotifs; i++) all_motif_los[i] = good_motif_los[i];

  /* remove highly correlated motifs from query */
  if (rem_corr && nbad) {
    for (i=j=k=0; i<nmotifs; i++) {
      if (good_motif_los[i]->imotif < bad_list[k] || k >= nbad) {
        good_motif_los[j++] = good_motif_los[i];
      } else {
        k++;
      }
    }
    nmotifs -= nbad;
  }

  /*
    Set current sequence alphabet to DNA if translating.
  */
  if (translate_dna) setalph(0);

  /* determine the type of database being searched */
  dna = (translate_dna || good_motif_los[0]->dna);

  /* get the type of handling of DNA strands */
  if (!dna) {                                   /* protein database */
    if (stype == Separate || stype == Norc) {
      die("You may not specify -sep or -norc with a protein database.\n");
    }
    stype = Protein;
  } /* database type */

  /*
    calculate p-values of all integer score values in range [0...w*MAST_RANGE]
  */
  Resize(pv, nmotifs, double *);
  for (i=0; i<nmotifs; i++) {
    pv[i] = calc_pssm_cdf(good_motif_los[i]->w, good_motif_los[i]->alen, MAST_RANGE, good_motif_los[i]->logodds, back);
    if (!pv[i]) {
      error("There is something wrong with motif %s\n", good_motif_los[i]->meme_name);
    }
  }
  exit_on_error();


  databases = init_databases(db_list, dna, db_source, db_name, db_link);

  /* Create and print a hit list and exit if requested. */
  if (hit_list) {
    for (i = 0; i < arraylst_size(databases); ++i) {
      DATABASE_T *db = (DATABASE_T*)arraylst_get(i, databases);
      print_hit_list(mast_out, db, good_motif_los, nmotifs, dna, stype, translate_dna, best_motifs, pv,
        m_thresh, w_thresh, status);
    }
    /* output arguments */
    fprintf(mast_out, "# mast");
    for (i = 1; i < argc; ++i) fprintf(mast_out, " %s", argv[i]);
    fprintf(mast_out, "\n");
    return 0;                           /* all done! */
  }

  // get the scores and p-values for all sequences in the database
  ARRAYLST_T *sequences = arraylst_create();
  ARRAYLST_T *hits = arraylst_create();
  for (i =  0; i < arraylst_size(databases); ++i) {
    DATABASE_T *db = (DATABASE_T*)arraylst_get(i, databases);
    get_scores(dna, stype, translate_dna, db, good_motif_los, nmotifs, 
      MAST_RANGE, pv, sonly, lump, use_seq_comp, e_max,
      kmotifs, motif, status, min_seqs, &nseqs, &residues, 
      sequences, hits);
  }

  // bail if no sequences were read succesfully
  if (nseqs == 0) {
    die("Quitting due to errors or empty database.\n"); 
  }

  //sort by score
  arraylst_qsort(score_asc, hits);
  //remove anything with an evalue worse than the threshold
  for (i = arraylst_size(hits)-1; i >= 0; --i) {
    STRAND_T *st = arraylst_get(i, hits);
    double evalue = st->Pvalue * nseqs;
    if (evalue <= e_max) break;
  }
  if (++i < arraylst_size(hits)) {//TODO FIXME check if the hits are freed somewhere
    arraylst_remove_range(i, arraylst_size(hits) - i, NULL, hits);
  }

  //output
  output_mast_xml(mast_out, translate_dna, 
      alphabet, use_seq_comp, bfile, back, 
      rem_corr, MAXCORR, 
      nmotifs, nbad, good_motif_los, all_motif_los, 
      mo_source, mo_name, mo_last_mod, 
      norder, order, space,
      nseqs, databases, hits, 
      argc, argv, stype, e_max, 
      adj_hit_pvalue, m_thresh, w_thresh);


  // finish xml output
  if (mast_out != stdout) fclose(mast_out);

  // finish status report
  if (status) fprintf(stderr, "\n");

  //cleanup
  for (i = 0; i < arraylst_size(sequences); ++i) {
    SSEQ_T *sseq = (SSEQ_T*)arraylst_get(i, sequences);
    sseq_destroy(sseq, nmotifs);
  }
  arraylst_destroy(NULL, sequences);
  for (i = 0; i < arraylst_size(databases); ++i) {
    DATABASE_T *db = (DATABASE_T*)arraylst_get(i, databases);
    myfree(db->source);
    myfree(db->name);
    myfree(db->link);
    if (db->save) fclose(db->save);
    if (db->file != stdin) fclose(db->file);
    myfree(db);
  }
  arraylst_destroy(NULL, databases);

  //output alternate formats
  if (xml_path && html) {
    char *stylesheet_path, *html_path;
    stylesheet_path = make_path_to_file(ETC_DIR, HTML_STYLESHEET);
    html_path = make_path_to_file(out_dir, HTML_FILENAME);
    if (file_exists(stylesheet_path)) {
      print_xml_filename_to_filename_using_stylesheet(
          xml_path,         /* path to XML input file IN */
          stylesheet_path,  /* path to MAST XSL stylesheet IN */
          html_path       /* path to HTML output file IN */
          );
    } else {
      if (verbosity >= NORMAL_VERBOSE) 
        fprintf(stderr, "Warning: could not find the stylesheet \"%s\" required for transformation of xml to html. Have you installed %s correctly?\n", 
          stylesheet_path, program_name);
    }
    myfree(stylesheet_path);
    myfree(html_path);
  }
  if (xml_path && mast2txt) {
    char *mast2txt_path, *txt_path;
    mast2txt_path = make_path_to_file(BIN_DIR, MAST2TXT_FILENAME);
    txt_path = make_path_to_file(out_dir, TXT_FILENAME);
    if (file_exists(mast2txt_path)) {
      char *cmd;
      int len;
      len = strlen(mast2txt_path) + strlen(xml_path) + strlen(txt_path) + 9;
      cmd = mm_malloc(sizeof(char) * len);
      snprintf(cmd, len, "\"%s\" \"%s\" \"%s\"", mast2txt_path, xml_path, txt_path);
      if (system(cmd)) {
        fprintf(stderr, "Warning: program mast2txt failed to run.\n");
      }
      myfree(cmd);
    } else {
      if (verbosity >= NORMAL_VERBOSE) 
        fprintf(stderr, "Warning: could not find the program \"%s\" required for transformation of xml to txt. Have you installed %s correctly?\n", 
          mast2txt_path, program_name);
    }
    myfree(mast2txt_path);
    myfree(txt_path);
  }

  myfree(xml_path);
  return 0;
} /* main */

