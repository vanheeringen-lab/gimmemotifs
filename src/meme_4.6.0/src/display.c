/*
 * $Id: display.c 5280 2011-01-05 11:10:34Z james_johnson $
 *
 * $Log$
 * Revision 1.4  2006/03/08 20:50:11  nadya
 * merge chamges from v3_5_2 branch
 *
 * Revision 1.3.4.1  2006/01/24 20:44:08  nadya
 * update copyright
 *
 * Revision 1.3.2.3  2006/03/02 03:59:25  twhitington
 * calc_pc2 changed name to calc_pc now that calculate_pc module is deleted.
 *
 * Revision 1.3.2.2  2006/02/22 07:10:46  twhitington
 * Added calculation and printing for performance coefficient.
 *
 * Revision 1.3.2.1  2006/02/16 00:02:06  twhitington
 * Introduced function calc_pc for calculating the performance coefficient.
 * Modified function print_results() and object MODEL for call to calc_pc.
 *
 * Revision 1.3  2005/10/25 19:00:37  nadya
 * change c++ style comment to proper c
 *
 * Revision 1.2  2005/10/20 00:20:52  tbailey
 * *** empty log message ***
 *
 * Revision 1.1.1.1  2005/07/29 23:35:53  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/***********************************************************************
*                      *
* MEME                     *
* Copyright 1994-2006, The Regents of the University of California *
* Author: Timothy L. Bailey              *
*                      *
***********************************************************************/
/**********************************************************************/
/*
  Display routines for the results of MEME
*/
/**********************************************************************/
/* 7-23-99 tlb; replace nsites() calls with model->nsites */
/* 7-16-99 tlb; move SITE to meme.h, get_sites to meme_util.c */
/* 7-14-99 tlb; change simplified prob. matrix to be observed frequencies,
  consensus to be based on observed frequencies */
/* 4-7-99 tlb; fix bug in get_sites: setting of best site in oops and zoops */

#include <assert.h>
#include "display.h"
#include "display_globals.h"		// isolate the global variables
#include "meme-dtd.h"
#include "mast.h"
#include "projrel.h"
#include "motif.h"
#include "config.h"
#include "ceqlogo.h"
#include "dir.h"
#include "alphabet.h"
#include <string.h>
#ifdef UNIX
#include <unistd.h>
#endif

/* FIXME ??? */
#ifndef PARALLEL
#define mpMyID() 0
#endif

static char *yesno[] = {"no", "yes"};
static char *stars = NULL;

/* number of different integral score values for computing p-values */
#define PVAL_RANGE 100

/* stepsize and size of smoothing window for get_q */
#define NSTEPS 100
/*#define WINDOW (NSTEPS/10)*/
#define WINDOW 0

/* round double to integer; round up if half way between */
/* only works for positive numbers */
#ifndef NINT
  #define NINT(x) ((int) ((x)+0.5))
#endif

/* maximum number of levels in consensus sequence */
#define MAXDEPTH ((int) (1.0/MINCONS))

/* encode and decode a (seq #, offset) pair as one integer */
/* this will bomb if seq. length > MAX_SEQ_LEN or
   if number of sequences * MAX_SEQ_LEN exceeds size of int
*/
#define MAX_SEQ_LEN 10000
#define ENCODE(I,J) (I) * MAX_SEQ_LEN + (J);
#define DECODE(N, SEQ, OFF) {\
  (OFF) = (N) % MAX_SEQ_LEN;\
  (SEQ) = ((N) - (OFF)) / MAX_SEQ_LEN;\
}

/* distance to indent start of RE histogram, consensus and simplified motif */
#define IND 13
#define IND2 6

/* record describing accuracy of motif */
typedef struct {
  double thresh;  /* optimal threshold */
  int err;    /* classification errors using thresh */
  double roc;   /* ROC */
} ACCURACY;

/* sortable sequence score record */
typedef struct {
  double score;   /* sequence score */
  int class;    /* class of sequence,  1=pos, 0=neg */
  char *id;
} SORTED_SCORE;

/* local functions */
static void print_sites(
  DATASET *dataset, /* the dataset */
  MODEL *model,     /* the model */
  int format,       /* 0=BLOCKS; 1=FASTA */
  char *com,        /* comment to append */
  FILE *outfile     /* where to print */
);
static void print_log_odds(
  int imotif,         /* motif number */
  DATASET *dataset,   /* the dataset */
  int w,              /* width of motif */
  BOOLEAN pair,       /* double matrix if true */
  LOGODDS logodds,    /* log-odds matrix */
  double bayes,       /* threshold */
  double log_ev,      /* log E-value of motif */
  FILE *outfile       /* where to print */
);
static void print_entropy(
  BOOLEAN logo,       /* true: prints the logo, false: otherwise */
  MODEL *model,       /* the model */
  DATASET *dataset,   /* the dataset */
  char *str_space,    /* space for printing strand direction */
  FILE *outfile       /* stream for output */
);
static void print_logo(
  MODEL *model,       /* the model */
  DATASET *dataset    /* the dataset */
);
static void print_candidates(
  CANDIDATE *candidates, /* list of candidate models IN */
  DATASET *dataset,      /* the dataset IN */
  int max_w,             /* maximum width for motifs IN */
  FILE *outfile          /* stream for output IN */
);
static ACCURACY *get_thresh(
  int w,        /* width of motif */
  LOGODDS logodds1,     /* log-odds matrix: LOG2(m_ij/b_j) */
  LOGODDS logodds2,     /* log-odds matrix: LOG2(m_ij/a_ij) */
  DATASET *pos,       /* positive examples */
  DATASET *neg,       /* negative examples */
  BOOLEAN print_scores      /* print sorted scores */
);
static double meme_score_sequence(
  char *eseq,   /* integer-coded sequence to score */
  int length,   /* length of the sequence */
  int w,    /* width of motif */
  LOGODDS logodds1,   /* log-odds matrix: LOG2(m_ij/b_j) */
  LOGODDS logodds2    /* log-odds matrix: LOG2(m_ij/n_ij) */
);
static int s_compare(
  const void *v1,
  const void *v2
);
static double get_q(
  int nsteps,         /* try nsteps from 0 to 1 */
  int window,         /* smoothing window radius */
  int w,          /* width of motif */
  THETA theta,          /* motif theta */
  THETA neg_theta,        /* anti-motif theta */
  double *back,         /* background motif */
  DATASET *dataset,       /* the dataset */
  DATASET *neg_dataset,       /* negative examples */
  char *str_space   /* space for printing strand direction */
);
static LO *create_lo(
  int imotif,       /* load number of motif */
  int w,        /* width of motif */
  int alen,       /* length of alphabet */
  LOGODDS logodds,      /* single-letter logodds matrix */
  BOOLEAN pair,
  double threshold      /* Bayes optimal threshold */
);
static void score_sites(
  DATASET *dataset,     /* the dataset */
  MODEL *model,       /* the model */
  LO *lo,       /* LO structure */
  double *pv        /* p-values for scores of this motif */
);
static void print_meme_header_xml(FILE *outfile /* output file */);
static void print_meme_training_set_xml(
  DATASET *dataset, /* the dataset */
  int nmotifs,      /* number of motifs */
  FILE* outfile     /* file for output */
);
static void print_meme_model_xml(
  MODEL *model,           /* the model */
  DATASET *dataset,       /* the dataset */
  char *stopping_reason,  /* LO structure */
  FILE* outfile           /* file for output */
);
static void print_meme_motifs_xml(
  MODEL *model,                   /* the model IN */
  DATASET *dataset,               /* the dataset */
  int nmotifs,                    /* number of motifs */
  MOTIF_SUMMARY *motif_summaries, /* List of final motif properties */
  FILE* outfile                   /* output file */
);
static void print_meme_pssm_xml(
  LOGODDS logodds, /* pointer to matrix of log-odds scores */
  int alength,     /* length of the alphabet */
  char* alphabet,  /* pointer to alphabet string */
  int width,       /* width of the motif */
  FILE* outfile    /* pointer to output file */
);
static void print_meme_psfm_xml(
  THETA theta,  /* pointer to matrix of frequencies */
  int alength,  /* length of the alphabet */
  char* alphabet,  /* pointer to alphabet string */
  int width,    /* width of the motif */
  FILE* outfile /* pointer to output file */
);
static void print_meme_regular_expression_xml(
  char* regexp,     /* regular expression  */
  FILE* outfile     /* pointer to output file */
);

static void print_meme_contributing_sites_xml(
  MODEL *model,				/* the model */
  MOTIF_SUMMARY *motif_summary, /* summary of the motif */
  DATASET* dataset, /* the dataset */
  FILE*  outfile    /* pointer to output file */
);

static void print_meme_scanned_sites_xml(
  MODEL *model,       /* the model */
  DATASET *dataset,     /* the dataset */
  LO *los[MAXG],      /* logodds structures for motifs */
  int nmotifs,        /* number of motifs */
  double *pv[MAXG],     /* p-value tables for each motif */
  FILE *outfile
);

static void print_site_diagrams(
  DATASET *dataset,     /* the dataset */
  MODEL *model,       /* the model */
  int nmotifs,        /* number of motifs in los */
  LO *los[MAXG],      /* logodds structure for motif */
  FILE *outfile       /* where to print */
);

void get_aligned_sequence_parts(
  int motif_start,
  int motif_width,
  BOOLEAN ic,
  int lseq,
  char *seq,
  char *pre,
  char *site,
  char *post
);
static void align_sites(
  DATASET *dataset,     /* the dataset */
  MODEL *model,       /* the model */
  LO *lo,       /* LO structure */
  double *pv,       /* pvalues for scores of this motif */
  FILE *outfile       /* stream for output */
);
static void print_block_diagrams(
  MODEL *model,       /* the model */
  DATASET *dataset,     /* the dataset */
  LO *los[MAXG],      /* logodds structures for motifs */
  int nmotifs,        /* number of motifs */
  double *pv[MAXG],     /* p-value tables for each motif */
  FILE *outfile
);
static void print_pllr (
  double pllr,       ///< LLR of the planted site alignment
  FILE *outfile      ///< The stream for output
);
static void calc_pllr(
  DATASET *data,     ///< The dataset of sequences containing the planted sites
  MODEL *model       ///< The model resulting from meme. The length of the
                     ///< planted sites is held here.
);
static char *get_regexp(
  MODEL *model,       /* the model */
  DATASET *dataset,   /* the dataset */
  int prosite         /* prosite ==1 [AK]-E-[EG]. ==0 [AK]E[EG] */
);


/**********************************************************************/
/*
  Record the results of EM
  Doesn't support negative datasets.
*/
/**********************************************************************/
extern void record_results(
  DATASET *dataset,              /* the dataset IN */
  MODEL *model,                  /* the model IN */
  MOTIF_SUMMARY *motif_summaries /* summaries of final motifs IN */
)
{
  int w = model->w;       /* width of last component */
  int nsites_dis = model->nsites_dis;   /* # of sites after discretiz.*/
  double m1, e1, m2, e2;      /* for printing significance */
  THETA theta = model->theta;
  THETA obs = model->obs;
  double lambda = model->lambda;
  LOGODDS logodds;
  int alength = dataset->alength;
  double *back = dataset->back;
  int imotif = model->imotif-1;     /* index of motif */
  double thresh;        /* Bayes optimal threshold */

  /* get p-value and E-value of motif */
  exp10_logx(model->logpv/log(10.0), m1, e1, 1);
  exp10_logx(model->logev/log(10.0), m2, e2, 1);

  /* Record the results for the model as a whole */
  motif_summaries[imotif].width = w;
  motif_summaries[imotif].num_sites = nsites_dis;
  motif_summaries[imotif].ic = model->ic;
  motif_summaries[imotif].re = w * model->rel;
  motif_summaries[imotif].llr = model->llr;
  motif_summaries[imotif].e_value_mant = m2;
  motif_summaries[imotif].e_value_exp = e2;
  // motif_summaries[imotif].sites = model->maxima; //FIXME CEG Copy this memory properly
  motif_summaries[imotif].sites = mymalloc(model->nsites_dis * sizeof(p_prob));
  memcpy(motif_summaries[imotif].sites, model->maxima, model->nsites_dis * sizeof(p_prob));

  /* Negative sample sets are not supported */
  /* make the log-odds matrices */
  logodds = make_log_odds(theta, NULL, back, 0, w, alength);
  /* calculate the optimal threshold (min classification error or Bayes' */
  thresh = LOG2((1-lambda)/lambda); /* Bayes' threshold */
  motif_summaries[imotif].bayes = thresh;
  create_2array(motif_summaries[imotif].pssm, double, w + 1, alength + 1);
  int i;
  int j;
  for (i = 0; i < w; i++) {
    for (j = 0; j < alength; j++) {
      motif_summaries[imotif].pssm[i][j] = logodds(i, j);
    }
  }
  create_2array(motif_summaries[imotif].psfm, double, w + 1, alength + 1);
  for (i = 0; i < w; i++) {
    for (j = 0; j < alength; j++) {
      motif_summaries[imotif].psfm[i][j] = obs(i, j);
    }
  }
  motif_summaries[imotif].regexp = get_regexp(model, dataset, 0);

  /* Record elapsed time */
  motif_summaries[imotif].elapsed_time = myclock()/1E6;

} /* Record results */

/**********************************************************************/
/*
  Print the results of EM
*/
/**********************************************************************/
extern void print_results(
  DATASET *dataset,      /* the dataset IN */
  DATASET *neg_dataset,  /* negative examples IN */
  MODEL *model,          /* the model */
  MODEL *neg_model,      /* negative model IN */
  CANDIDATE *candidates, /* candidate models found IN */
  FILE* outfile          /* file for text output IN */
)
{
  int i;
  int max_w = model->max_w;     /* maximum width tried */
  int nstrands = model->invcomp ? 2 : 1;  /* # of strands to use */
  int w = model->w;       /* width of last component */
  int nsites_dis = model->nsites_dis;   /* # of sites after discretiz.*/
  double m1, e1, m2, e2;      /* for printing significance */
  char *str_space = (nstrands == 1) ? "" : "       ";
  THETA theta = model->theta;
  THETA obs = model->obs;
  double lambda = model->lambda;
  LOGODDS logodds;
  char *cons = model->cons;
  int alength = dataset->alength;
  double *back = dataset->back;
  int imotif = model->imotif-1;     /* index of motif */
  double thresh;        /* Bayes optimal threshold */
  double *rounded_back = NULL;
  Resize(rounded_back, alength, double);

  /* get entropy and E-value of motif */
  calc_entropy(model, dataset);

  // Retrieve the array of sites predicted in the model:
  P_PROB pred_sites = model->maxima;
  int n_psites = model->nsites_dis;

  // If requested, print the final MEME site predictions (for the "best"
  // starting point after EM has been completed):
  if (dataset->print_pred) {
    fprintf(stdout, "PREDICTED SITES AFTER MEME COMPLETION MOTIF %i\n",
            model->imotif);
    print_site_array(pred_sites, n_psites, stdout, model->w, dataset);
    double sig = dataset->objfun==Pv ? model->logpv : model->logev;
    fprintf(stdout, "FINAL MODEL SIGNIFICANCE = %f\n", sig);
  }

  // Calculate the llr of the alignment of the planted sites:
  if (dataset->print_pllr) {
    calc_pllr(dataset, model);
  }

  /* get p-value and E-value of motif */
  exp10_logx(model->logpv/log(10.0), m1, e1, 1);
  exp10_logx(model->logev/log(10.0), m2, e2, 1);

  /* print the significant models */
  if (VERBOSE) {
    print_candidates(candidates, dataset, max_w, outfile);
  }

  /* print the results for the model as a whole */
  fprintf(
    outfile,
    "\n\n%s\nMOTIF%3d\twidth = %4d   sites = %3d   ",
    stars, model->imotif, w, nsites_dis
  );

  fprintf(
    outfile,
    "llr = %.0f   E-value = %3.1fe%+04.0f\n%s\n",
    model->llr, m2, e2, stars
  );

  if (VERBOSE) {
    fprintf(
      outfile,
      "p-value = %3.1fe%+04.0f   E-value = %3.1fe%+04.0f\n%s\n",
      m1, e1, m2, e2, stars
    );
  }

  /* print results for motif */
  if (VERBOSE) {
    char *cons0 = model->cons0;
    fprintf(outfile, "\n(best) %s --> %s\n", cons0, cons);
    fprintf(
      outfile,
      "(best) w %3d nsites %5.1f lambda %8.7f RE/col %6.3f\n",
      w, lambda*wps(dataset, w), lambda, model->rel
    );
    fprintf(outfile, "\n(best) RE %6.3f\n\n", w * model->rel);
  }

  /*
    print motif description
  */
  for (i=0; i<PAGEWIDTH; i++) {
    fputc('-', outfile);
  }
  fputc('\n', outfile);
  fprintf(outfile, "\tMotif %d Description\n", model->imotif);
  for (i=0; i<PAGEWIDTH; i++) {
    fputc('-', outfile);
  }
  fputc('\n', outfile);

  /* make the log-odds matrices */
  /* calculate the optimal threshold (min classification error or Bayes' */
  logodds = make_log_odds(theta, NULL, back, 0, w, alength);
  thresh = LOG2((1-lambda)/lambda); /* Bayes' threshold */
  print_theta(
    0, 2, model->nsites_dis, obs, w, model->logev,
    str_space, dataset, outfile
  );
  print_entropy(TRUE, model, dataset, str_space, outfile);
  for (i=0; i<PAGEWIDTH; i++) {
    fputc('-', outfile);
  }
  fprintf(outfile, "\n\n");

  /* create an LO structure and store it in the local array */
  los[imotif] = create_lo(model->imotif, w, alength, logodds, FALSE, thresh);

  // round background so p-values are same as for MAST
  for (i=0; i<alength; i++) RND(back[i], 3, rounded_back[i]);

  /* create a table of p-values and store it in the array */
  /*pv[imotif] = calc_cdf(los[imotif], PVAL_RANGE, dataset->back);*/
  pv[imotif] = calc_pssm_cdf(
    los[imotif]->w, los[imotif]->alen-1, PVAL_RANGE,
    los[imotif]->logodds, rounded_back
  );

  /* score the sites and sort by position p-value */
  score_sites(dataset, model, los[imotif], pv[imotif]);

  /* print alignment of the sites */
  align_sites(dataset, model, los[imotif], pv[imotif], outfile);

  /* print diagrams of the sites */
  print_site_diagrams(dataset, model, model->imotif, los, outfile);

  /* print the sites "making up" the model */
#ifdef PARALLEL
#ifdef DEBUG_PARALLEL
  fprintf(stderr, "%d: at print_print_sites\n", mpMyID()); fflush(stderr);
#endif
#endif
  print_sites(dataset, model, PRINT_FASTA, "", outfile);
#ifdef DEBUG_PARALLEL
  fprintf(stderr, "%d: past print_print_sites\n", mpMyID()); fflush(stderr);
#endif

  /* print the logodds matrix */
  print_log_odds(
    model->imotif, dataset, w, FALSE, logodds, thresh, model->logev, outfile
  );

  /* print the observed frequency matrix */
  print_theta(
    model->imotif, 1, model->nsites_dis, obs, w, model->logev,
    str_space, dataset, outfile
  );

  // Print the LLR of the alignment of the planted sites, if required:
  if (dataset->print_pllr) {
    print_pllr(dataset->pllr, stdout);
  }

  /* 
    print a regular expression corresponding to motif
  */
  for (i=0; i<PAGEWIDTH; i++) {
    fputc('-', outfile);
  }
  fputc('\n', outfile);
  fprintf(outfile, "\tMotif %d regular expression\n", model->imotif);
  for (i=0; i<PAGEWIDTH; i++) {
    fputc('-', outfile);
  }
  fputc('\n', outfile);
  fprintf(outfile, "%s\n", get_regexp(model, dataset, 0));
  for (i=0; i<PAGEWIDTH; i++) {
    fputc('-', outfile);
  }
  fprintf(outfile, "\n\n");

  /* display elapsed time */
  fprintf(outfile, "\n\n\nTime %5.2f secs.\n\n", myclock()/1E6);

  /* print line of stars */
  fprintf(outfile, "%s\n", stars);

  /* flush */
  fflush(outfile);

  myfree(rounded_back);

} /* print_results */

/**********************************************************************/
/*
  create_lo

  Create an an LO structure from the logodds matrix;
  include 'X' in the alphabet.
  PSSM is scaled by 100 rounded to the nearest integer,
  and then rescaled to PVAL_RANGE.
*/
/**********************************************************************/
static LO *create_lo(
  int imotif,			/* index of motif */
  int w,			/* width of motif */
  int alen,			/* length of alphabet */
  LOGODDS logodds,		/* single-letter logodds matrix */
  BOOLEAN pair,
  double threshold		/* Bayes optimal threshold */
)
{
  int i, j, len, tmp;
  LO *lo = NULL;      /* the LO structure */

  /* create a logodds structure */
  Resize(lo, 1, LO);

  /* initialize it */
  lo->alen = alen+1;      /* add 'X' column */
  lo->w = lo->ws = w;
  lo->pair = pair;
  lo->imotif = imotif;
  //convert the index into a string
  len = 1, tmp = imotif;
  while (tmp >= 1) { 
    tmp /= 10;
    ++len;
  }
  lo->meme_name = mymalloc(len*sizeof(char));
  snprintf(lo->meme_name, len, "%d", imotif);
  lo->thresh = threshold;

  /* make a copy of the logodds matrix and scale it to [0..range] */
  create_2array(lo->logodds, LOGODDSB, w, lo->alen);
  // make like when we print it out!
  for (i=0; i<w; i++) for(j=0; j<lo->alen; j++) lo->logodds(i,j) = NINT(100*logodds(i,j));
  //for (i=0; i<w; i++) for(j=0; j<lo->alen; j++) lo->logodds(i,j) = logodds(i,j);
  scale_lo(&lo, 1, PVAL_RANGE);        /* scale */
  make_double_lo(&lo, 1); /* make a double-letter logodds matrix */

  return(lo);
} /* create_lo */

/**********************************************************************/
/*
  print_block_diagrams

  Tile the dataset sequences with the motifs in los[] and print
  the block diagrams with the p-value of the product of p-values.
*/
/**********************************************************************/
static void print_block_diagrams(
  MODEL *model,       /* the model */
  DATASET *dataset,     /* the dataset */
  LO *los[MAXG],      /* logodds structures for motifs */
  int nmotifs,        /* number of motifs */
  double *pv[MAXG],     /* p-value tables for each motif */
  FILE *outfile
)
{
  int i;
  BOOLEAN dna = dataset->dna;   /* dataset is DNA */
  BOOLEAN invcomp = model->invcomp; /* use reverse complement strand */
  STYPE stype = dna ? (invcomp ? Combine : Norc) : Protein;
  BOOLEAN xlate_dna = FALSE;    /* don't translate DNA */
  BOOLEAN best_motifs = FALSE;    /* use all motifs */
  double m_thresh = 1e-4;   /* show motifs over p-value 0.0001 */
  double w_thresh = m_thresh;   /* show strong motifs only */
  BOOLEAN use_seq_p = FALSE;    /* use postion p-values */
  char *f = "%-*s%s %8s  %s\n";   /* format */

  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fputc('\n', outfile);
  fprintf(outfile, "\tCombined block diagrams:");
  fprintf(outfile, " non-overlapping sites with p-value < %6.4f\n", m_thresh);
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fputc('\n', outfile);
  fprintf(outfile, f, MSN, "SEQUENCE NAME", "", "COMBINED P-VALUE", "MOTIF DIAGRAM");
  fprintf(outfile, f, MSN, "-------------", "", "----------------", "-------------");

  for (i=0; i<dataset->n_samples; i++) {
    SAMPLE *s = dataset->samples[i];
    char *name = s->sample_name;
    int lseq = s->length;
    char *sequence = s->seq;
    TILING tiling = score_tile_diagram(NULL, sequence, lseq, los, nmotifs,
      dna, stype, FALSE, xlate_dna, best_motifs, TRUE, pv, m_thresh, w_thresh,
      use_seq_p, FALSE, NULL);
    fprintf(outfile, "%-*.*s %16.2e  %s\n", MSN, MSN, name, tiling.pv,
      tiling.diagram);
    free_tiling(tiling);
  } /* sequence */

  /* print a final line of hyphens */
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fprintf(outfile, "\n\n");

} /* print_block_diagrams */

/**********************************************************************/
/*
  print_meme_scanned_sites_xml

  Tile the dataset sequences with the motifs in los[] and print
  the XML for the motif occurence strings.
*/
/**********************************************************************/
static void print_meme_scanned_sites_xml(
  MODEL *model,     /* the model */
  DATASET *dataset, /* the dataset */
  LO *los[MAXG],    /* logodds structures for motifs */
  int nmotifs,      /* number of motifs */
  double *pv[MAXG], /* p-value tables for each motif */
  FILE *outfile     /* write XML to this file handle */
)
{
  int i_seq;
  BOOLEAN dna = dataset->dna;   /* dataset is DNA */
  BOOLEAN invcomp = model->invcomp; /* use reverse complement strand */
  STYPE stype = dna ? (invcomp ? Combine : Norc) : Protein;
  BOOLEAN xlate_dna = FALSE;    /* don't translate DNA */
  BOOLEAN best_motifs = FALSE;    /* use all motifs */
  double m_thresh = 1e-4;   /* show motifs over p-value 0.0001 */
  double w_thresh = m_thresh;   /* show strong motifs only */
  BOOLEAN use_seq_p = FALSE;    /* use postion p-values */

  fprintf(
    outfile,
    "<scanned_sites_summary p_thresh=\"%.4f\">\n",
    m_thresh
  );
  for (i_seq =0; i_seq <dataset->n_samples; i_seq ++) {
    SAMPLE *s = dataset->samples[i_seq ];
    char *name = s->sample_name;
    int lseq = s->length;
    char *sequence = s->seq;
    TILING tiling = score_tile_diagram(
      NULL,sequence, lseq, los, nmotifs,
      dna, stype, FALSE, xlate_dna, best_motifs, TRUE, pv, m_thresh, w_thresh,
      use_seq_p, FALSE, NULL
    );
    // Count the number of sites
    int i_pos = 0;
    int num_sites = 0;
    for (i_pos = 0; i_pos < lseq; i_pos++) {
      if (tiling.hits[i_pos] != 0) {
        num_sites++;
      }
    }
    // Print opening scanned sites tag
    // contains info about sequence.
    fprintf(
      outfile, 
      "<scanned_sites sequence_id=\"sequence_%d\""
      " pvalue=\"%.2e\" num_sites=\"%d\">",
      i_seq,
      tiling.pv,
      num_sites
    );
    for (i_pos = 0; i_pos < lseq; i_pos++) {
      // Print a scanned_site element for each site.
      // contains info about motif occurence.
      int motif_index = tiling.hits[i_pos]; 
      if (motif_index != 0) {
        char* strand = "none";
        if (dna == TRUE) {
          strand = motif_index > 0 ? "plus" : "minus";
        }
        fprintf(
          outfile, 
          "<scanned_site motif_id=\"motif_%d\" strand=\"%s\""
          " position=\"%d\" pvalue=\"%.2e\"/>\n",
          abs(tiling.hits[i_pos]),
          strand,
          i_pos,
          tiling.pvalues[i_pos]
        );
      }
    }
    // Print closing scanned sites tag
    fprintf( outfile, "</scanned_sites>\n");
    free_tiling(tiling);
  } /* sequence */
  fprintf(
    outfile,
    "</scanned_sites_summary>\n"
  );

} /* print_meme_scanned_sites_xml */

/**********************************************************************/
/*
  print_log_odds

  Print the log-odds matrix
*/
/**********************************************************************/
static void print_log_odds(
  int imotif,       /* motif number */
  DATASET *dataset, /* the dataset */
  int w,            /* width of motif */
  BOOLEAN pair,     /* double matrix if true */
  LOGODDS logodds,  /* log-odds matrix */
  double bayes,     /* threshold */
  double log_ev,    /* log E-value of motif */
  FILE *outfile     /* output file */
)
{
  int i, j;
  int alength = dataset->alength; /* length of alphabet */
  int n = wps(dataset, w);    /* weighted possible sites */
  char *type = pair ? "pair" : "";  /* type of matrix */
  double m1, e1;      /* for printing significance */

#ifdef PARALLEL
#ifdef DEBUG_PARALLEL
  fprintf(stderr, "%d: at print_print_logodds\n", mpMyID()); fflush(stderr);
#endif
#endif

  /* get E-value of motif */
  exp10_logx(log_ev/log(10.0), m1, e1, 1);

  /* double w if printing a matrix pair */
  if (pair) w *=2;

  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile);
  fprintf(outfile, "\n\tMotif %d position-specific scoring matrix\n", imotif);
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile);
  fprintf(outfile, "\n");
  fprintf(outfile,
    "log-odds matrix: alength= %d w= %d n= %d bayes= %g E= %3.1fe%+04.0f %s\n",
    alength, w, n, bayes, m1, e1, type);

  for (i=0; i < w; i++) {   /* site position */
    for (j=0; j < alength; j++) { /* letter */
      fprintf(outfile, "%6d ", NINT(100*logodds(i,j)));
    }
    fprintf(outfile, "\n");
  }
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile);
  fprintf(outfile, "\n\n");

#ifdef PARALLEL
#ifdef DEBUG_PARALLEL
  fprintf(stderr, "%d: leaving print_print_logodds\n", mpMyID()); fflush(stderr);
#endif
#endif
} /* print_log_odds */

/**********************************************************************/
/*
  print_entropy

  Displays the relative entropy of each column of the motif
  as a bar graph and as a logo.
*/
/**********************************************************************/
static void print_entropy(
  BOOLEAN logo,       /* true: prints the logo, false: otherwise */
  MODEL *model,       /* the model */
  DATASET *dataset,   /* the dataset */
  char *str_space,    /* space for printing strand direction */
  FILE *outfile       /* stream for output */
)
{
  int i, j;
  int w = model->w;       /* width of motif */
  THETA obs = model->obs;     /* observed frequencies */
  double *rentropy = model->rentropy;   /* RE of each column */
  double re = w *model->rel;     /* motif relative entropy */
  double *back = dataset->back;     /* background model */
  int alength = dataset->alength;
  char restring[15];        /* print string for re */
  char *consensus;        /* consensus strings */
  double min_freq;        /* minimum background freq */
  double maxre;         /* maximum relative entropy */
  int nsteps;         /* number of steps histogram */

  /* get minimum background frequency and maximum relative entropy */
  for (i=0, min_freq=1; i<alength; i++)
    if (back[i] < min_freq) min_freq = back[i];
  maxre = -LOG2(min_freq);    /* maximum relative entropy */

  /* create string containing RE */
  sprintf(restring, "(%.1f bits)", re);

  /* print the relative entropy of each column as a bar graph */
  nsteps = 10;
  for (i=0; i<nsteps; i++) {
    double level = maxre - (i * maxre/nsteps);
    fprintf(outfile, (i==0 ? "%*.*s %*.1f " : "%-*.*s %*.1f "), IND, IND,
      (i==0 ? "bits" : i==4 ? "Relative" : i==5 ? "Entropy" :
        i==6 ? restring : ""), IND2, level);
    for (j=0; j<w; j++) {
      if (NINT(nsteps * rentropy[j] / maxre) >= nsteps-i) {
        fputc('*', outfile);
      } else {
        fputc(' ', outfile);
      }
    }
    fputc('\n', outfile);
  }
  fprintf(outfile, "%-*.*s %*.1f ", IND, IND, "", IND2,0.0);
  for (i=0; i<w; i++) fputc('-', outfile);
  fprintf(outfile, "\n\n");
  /* get and print the consensus sequences */
  consensus = get_consensus(obs, w, dataset, MAXDEPTH, MINCONS);
  for (i=0; i<MAXDEPTH && i<alength; i++) {/* print next levels of consensus */
    fprintf(outfile, "%-*.*s %*.0s %*.*s\n", IND, IND,
      (i==0 ? "Multilevel" : i == 1 ? "consensus" : i == 2 ? "sequence" : ""),
      IND2, "", w, w, consensus+(i*w));
  }
  /* free up space */
  myfree(consensus);

  /* Prints a logo in EPS and PNG format to two files in the output directory */
  if(logo) print_logo(model, dataset);

} /* print_entropy */


/**********************************************************************/
/*
  print_logo

  Print the logo of a motif
*/
/**********************************************************************/
static void print_logo(
  MODEL *model,       /* the model */
  DATASET *dataset    /* the dataset */
)
{
  char* logodir = dataset->output_directory;
  double logo_height = LOGOHEIGHT;
  double logo_width =  model->w <= MAXLOGOWIDTH ? model->w : MAXLOGOWIDTH;

  if(logodir != NULL) {
    /* convert theta to motif */
    MOTIF_T motif;
    strcpy(motif.id, "0");
    motif.num_sites  = model->nsites_obs;
    motif.freqs      = convert_matrix(model->obs, model->w, dataset->alength); 
    motif.length     = model->w;
    motif.alph_size  = dataset->alength;
    motif.ambigs     = 0; //meme does not use ambigous characters it seems
    motif.evalue     = 0.0;
    motif.complexity = 0.0;
    motif.trim_left = 0;
    motif.trim_right = 0;

    // create the output path
    char *path = NULL;
    Resize(path, strlen(logodir)+29, char);	// room for "/logo_ssc<16digts>\0"

    // create logo without small sample correction
    sprintf(path, "%s/logo%d", logodir, model->imotif);
    CL_create2(
      &motif, 			// first motif
      "", 			// no title 
      NULL, 			// no second motif
      "", 			// no x-axis label
      FALSE, 			// no error bars
      FALSE,			// ssc
      logo_height,		// logo height (cm)
      logo_width,		// logo width (cm)
      dataset->alphabet, 	// alphabet
      0, 			// no offset to second motif
      path,			// output file path
      "MEME (no SSC)"		// program name
    );

    if (dataset->dna) {
      reverse_complement_motif(&motif);
      sprintf(path, "%s/logo_rc%d", logodir, model->imotif);
      CL_create2(
        &motif, 			// first motif
        "", 			// no title 
        NULL, 			// no second motif
        "", 			// no x-axis label
        FALSE, 			// no error bars
        FALSE,			// ssc
        logo_height,		// logo height (cm)
        logo_width,		// logo width (cm)
        dataset->alphabet, 	// alphabet
        0, 			// no offset to second motif
        path,			// output file path
        "MEME (no SSC)"		// program name
      );
    }
    /*
    // create logo with small sample correction and error bars
    sprintf(path, "%s/logo_ssc%d", logodir, model->imotif);
    CL_create2(
      &motif, 			// first motif
      "", 			// no title 
      NULL, 			// no second motif
      "", 			// no x-axis label
      TRUE, 			// error bars
      TRUE,			// ssc
      logo_height,		// logo height (cm)
      logo_width,		// logo width (cm)
      dataset->alphabet, 	// alphabet
      0, 			// no offset to second motif
      path,			// output file path
      "MEME (with SSC)"		// program name
    );
    */

    myfree(path);
    free_matrix(motif.freqs);
  }
} // print_logo

/**********************************************************************/
/*
  print_theta

    format=1    floating point; pos x letter
    format=2    1 digit; letter x pos

  Print the probability array.
*/
/**********************************************************************/
extern void print_theta(
  int imotif,   /* number of motif */
  int format,   /* 1 = floating point
                     2 = integer */
  int nsites,   /* number of sites (discrete) */
  THETA theta,    /* theta */
  int w,    /* width of motif */
  double log_ev,  /* log motif E-value */
  char *str_space,  /* space for printing strand direction */
  DATASET *dataset, /* the dataset */
  FILE *outfile   /* file to print to */
)
{
  int i, j;
  int alength = dataset->alength;

  if (format == 1) {
    double e1, m1;
    exp10_logx(log_ev/log(10.0), m1, e1, 3);
    for (i=0; i<PAGEWIDTH; i++) {
      fputc('-', outfile);
    }
    fprintf(
      outfile, "\n\tMotif %d position-specific probability matrix\n", imotif
    );
    for (i=0; i<PAGEWIDTH; i++) {
      fputc('-', outfile);
    }
    fprintf(
      outfile,
      "\nletter-probability matrix: alength= %d w= %d "
      "nsites= %d E= %3.1fe%+04.0f ",
      alength, w, nsites, m1, e1
    );
    fprintf(outfile, "\n");
    for (i=0; i < w; i++) {
      for (j=0; j < alength; j++) {
        fprintf(outfile, "%9.6f ", theta(i, j));
      }
      fprintf(outfile, "\n");
    }
    for (i=0; i<PAGEWIDTH; i++) {
      fputc('-', outfile);
    }
    fprintf(outfile, "\n");

  } else if (format == 2) {
    /* print theta: rows = letter; cols = position in motif; 1-digit integer */
    for (i=0; i < alength; i++) {
      /* print the letter */
      fprintf(
        outfile, "%-*.*s%*c  ", IND, IND,
        (i==0 ? "Simplified" : i==1 ? "pos.-specific" : i==2 ?
          "probability" : i==3 ? "matrix" : "" ), IND2, unhash(i)
      );
      for (j=0; j < w; j++) {
        int k = NINT(10.0 * theta(j,i));  /* round to 1 digit */
        if (k == 0) {
          fprintf(outfile, ":");    /* print 0 as colon */
        } else {
          fprintf(outfile, "%1x", k);   /* print 1 digit */
        }

      }
      fprintf(outfile, "\n");
    }
  }

  fprintf(outfile, "\n");
}

/**********************************************************************/
/*
  print_zij
*/
/**********************************************************************/
extern void print_zij(
  DATASET *dataset,     /* the dataset */
  MODEL *model        /* the model */
)
{
  int i, j;
  int n_samples = dataset->n_samples;
  SAMPLE **samples = dataset->samples;
  FILE *out=stdout;

  fprintf(out, "z_ij: lambda=%f ll=%f\n", model->lambda, model->ll);
  for (i=0; i<n_samples; i++) {     /* sequence */
    int lseq = samples[i]->length;
    double *zi = samples[i]->z;		// zi[j], j in [-lseq...+lseq]
    int w = model->w;
    fprintf(out, ">%s\nz : ", samples[i]->sample_name);
    for (j=0; j<lseq-w+1; j++) {    	/* position */
      int k = j+1;			// Z_i = k
      double z = model->invcomp ? MIN(1.0,Zi(-k)+Zi(k)) : Zi(k);
      int zij = NINT(10 * z);  /* round z */
      fprintf(out, "%1x", zij);
    } // position
    // print s0 and s1 for backwards compatibility
    if (model->invcomp) {
      fprintf(out, "\ns0: ");
      for (j=0; j<lseq-w+1; j++) {    	/* position */
	int k = j+1;			// Z_i = k
	double z = Zi(k);
	int zij = NINT(10 * z);  /* round z */
	fprintf(out, "%1x", zij);
      } // position
      fprintf(out, "\ns1: ");
      for (j=0; j<lseq-w+1; j++) {    	/* position */
	int k = j+1;			// Z_i = k
	double z = Zi(-k);
	int zij = NINT(10 * z);  /* round z */
	fprintf(out, "%1x", zij);
      } // position
    }
    fprintf(out, "\n");
  } /* sequence */
  printf("\n");
} /* print_zij */

/**********************************************************************/
/*
  print_wij
*/
/**********************************************************************/
extern void print_wij(
  DATASET *dataset      /* the dataset */
)
{
  int i,j;
  int n_samples = dataset->n_samples;
  SAMPLE **samples = dataset->samples;

  printf("w_ij:\n");
  for (i=0; i<n_samples; i++) {               /* sequence */
    int len = samples[i]->length;
    double *weights = samples[i]->weights;
    printf(">%s\n", samples[i]->sample_name);
    for (j=0; j<len; j++) {                   /* position */
      int w = NINT(10 * weights[j]);
      printf("%1x", w);
    }
    printf("\n");
  }
  printf("\n");
}

/**********************************************************************/
/*      get_consensus

        Get the consensus string from a motif.
  For each position, N consensus letters are found.
  If no letter has probability > min_prob,
        'x' is written for the first consensus
  letter and ' ' in the others.
        Otherwise, N letters are written in decreasing order of
  probability until one with min_prob is reached, and then ' ' is
  written for the remaining letters.
*/
/**********************************************************************/
extern char *get_consensus(
  THETA theta,      /* motif theta */
  int w,      /* width of motif */
  DATASET *dataset,   /* the dataset */
  int N,      /* number of letters for each position */
  double min_prob   /* minimum cumulative prob for N letters */
)
{
  int i, j, n;
  int alength = dataset->alength;
  char *alphabet = dataset->alphabet;
  char *string = NULL;

  Resize(string, w*N+2, char);

  for (i=0; i < w; i++) {   /* position in motif */
    int maxj[MAXDEPTH];     /* array of max indices in Theta */

    /* find N letters at position i with highest probability (in order) */
    for (n = 0; n < N; n++) {   /* current depth */
      double max = LITTLE;    /* current max probability */
      for (j=0; j < alength; j++) { /* letter */
  if (theta(i, j) > max) {
    max = theta(i, j);    /* maximum probability */
    maxj[n] = j;      /* current letter with n-th best prob */
  }
      }
      theta(i, maxj[n]) = -theta(i, maxj[n]); /* tag this position as used */
    }

    /* restore theta */
    for (n = 0; n < N; n++) {   /* current depth */
      theta(i, maxj[n]) = -theta(i, maxj[n]);   /* untag */
    }

    /* set up the consensus strings for position i */
    for (n = 0; n < N; n++) {       /* current depth */
      if (theta(i, maxj[n]) < min_prob) {
        string[(n*w)+i] = (n==0 ? 'x' : ' '); /* below cutoff */
      } else {
        string[(n*w)+i] = alphabet[maxj[n]];  /* set n'th consensus */
      }
    }
  }
  string[((N-1)*w)+i] = '\0';     /* terminate string */

  return string;
} /* get_consensus */


/**********************************************************************/
/*      get_regexp

  Get a regular expression using the same rules as for
  the multilevel consensus sequence.
*/
/**********************************************************************/
static char *get_regexp(
  MODEL *model,     	/* the model */
  DATASET *dataset,   	/* the dataset */
  int prosite           /* prosite ==1 [AK]-E-[EG]. ==0 [AK]E[EG] */
)
{
  THETA obs = model->obs; /* motif observed theta */
  int w = model->w;   	  /* width of motif */
  int i, j;
  char *pcons = get_consensus(obs, w, dataset, MAXDEPTH, MINCONS);
  char *re = NULL;    	  /* RE string to return */

  /* regular expression string */
  Resize(re, w*(MAXDEPTH+6), char);

  /* create the regular expression from the "packed" consensus */
  int pos = 0;
  for (i=0; i<w; i++) {
    /* see if there is more than one letter in consensus for this column */
    if (i>0 && prosite==1) re[pos++] = '-';
    if (MAXDEPTH > 1 && pcons[w + i] != ' ') re[pos++] = '[';
    for (j=0; j<MAXDEPTH; j++) {  /* copy consensus */
      char a = pcons[w*j + i];
      if (a == ' ') break;
      re[pos++] = a;
    }
    if (MAXDEPTH > 1 && pcons[w + i] != ' ') re[pos++] = ']';
  }
  if (prosite==1) re[pos++] = '.';
  re[pos++] = '\0';

  return re;
} /* get_regexp */

/**********************************************************************/
/*
  print_candidates

  Print the candidate models found.

*/
/**********************************************************************/
static void print_candidates(
  CANDIDATE *candidates, /* list of candidate models IN */
  DATASET *dataset,      /* the dataset IN */
  int max_w,             /* maximum width for motifs IN */
  FILE* outfile          /* file for text output IN */
)
{
  int i, w;
  int hdr = 35;
  int tail = PAGEWIDTH - hdr;
  double m, e;      /* for printing signficance */

  fprintf(outfile, "\nCandidate motifs:\n");
  fprintf(outfile, "width nsites  ll        signif     consensus\n");
  fprintf(outfile, "----- ------  --------- ---------- ---------\n");

  for (w=1; w<=max_w; w++) {

    if (candidates[w].sig > 1) continue;  /* skip unused candidates */

    exp10_logx(candidates[w].sig/log(10.0), m, e, 3);

    fprintf(
      outfile,
      "%5d %6.1f %c%9.0f %5.3fe%+04.0f ",
      w,
      candidates[w].lambda * wps(dataset, w),
      (candidates[w].pal ? 'P' : ' '),
      candidates[w].ll,
      m, e
    );
    fprintf(outfile, "%-*.*s\n", tail, tail, candidates[w].cons);
    for (i=tail; i < w; i+=tail) {
      fprintf(
        outfile, "%*.*s%-*.*s\n", hdr, hdr, "",
        tail, tail, candidates[w].cons+i
      );
    }
  }
} /* print_candidates */

/**********************************************************************/
/*
  print_dataset_summary
*/
/**********************************************************************/
extern void print_dataset_summary (
  DATASET *dataset, /* the dataset IN */
  FILE *outfile     /* where to print IN */
)
{
  int i, pcol;
  int w = MSN + 15;
  char *datafile = dataset->datafile; /* name of the training set file */
  char *alphabet = dataset->alphabet; /* alphabet of sequences */
  char *negfile = dataset->negfile; /* name of negative example file */

  /* set up printing spacers */
  Resize(stars, PAGEWIDTH+1, char);
  for (i=0; i<PAGEWIDTH; i++) {
    stars[i] = '*';
  }
  stars[PAGEWIDTH] = '\0';

  /* announce the training set */
  fprintf(outfile, "%s\nTRAINING SET\n%s\n", stars, stars);

  /* print name of file and alphabet */
  fprintf(
    outfile,
    "DATAFILE= %s\n"
    "ALPHABET= %s\n", datafile, alphabet
  );

  /* print name of negative dataset */
  if (negfile){
    fprintf(outfile, "NEGATIVES= %s\n", negfile);
  }

  /*
    print a table of sequence lengths
  */

  /*   print table header */
  for (pcol = w; pcol < 80; pcol += w) {
    fprintf(
      outfile,
      "%-*.*s %6s %6s%2s",
      MSN, MSN, "Sequence name", "Weight", "Length", " "
    );
  }
  fprintf(outfile, "\n");
  for (pcol = w; pcol < 80; pcol += w) {
    fprintf(
      outfile,
      "%-*.*s %6s %6s%2s",
      MSN, MSN, "-------------", "------", "------", " "
    );
  }
  fprintf(outfile, "\n");

  /*   print table columns */
  pcol = 0;
  for (i=0; i<dataset->n_samples; i++) {
    SAMPLE *sample = dataset->samples[i];
    char *sample_name = sample->sample_name;
    double wgt = sample->sw;
    int lseq = sample->length;
    /* print the sample name and its length */
    pcol += w;          /* new line for print out? */
    if (pcol >= 80) {
      fprintf(outfile, "\n");
      pcol = w;
    }
    fprintf(
      outfile,
      "%-*.*s %6.4f %6d%2s",
      MSN, MSN, sample_name, wgt, lseq, " "
    );
  }

  /* finish section */
  fprintf(outfile, "\n%s\n\n", stars);
} /* print_dataset_summary */

/**********************************************************************/
/*
  print_command_summary

  Print the command line summary
*/
/**********************************************************************/
extern void print_command_summary(
  MODEL *model,     /* the model IN */
  DATASET *dataset, /* the dataset IN */
  FILE *outfile     /* where to print IN */
)
{
  int i, pcol;
  char evt_string[12];

  if (dataset->evt == BIG) {
    strcpy(evt_string, "inf");
  } else {
    sprintf(evt_string, "%8g", dataset->evt);
  }

  fprintf(
    outfile,
    "%s\nCOMMAND LINE SUMMARY\n%s\n"
    "This information can also be useful in the event you wish to report a\n"
    "problem with the MEME software.\n\n"
    "command: %s\n\n"
    "model:  mod=      %8s    nmotifs=  %8d    evt=      %8s\n"
    "object function=  %s\n"
    "width:  minw=     %8d    maxw=     %8d    minic=    %8.2f\n",
    stars, stars,
    dataset->command,
    dataset->mod, dataset->nmotifs, evt_string,
    (dataset->objfun == Pv ?
       "P-value of product of p-values" :
       "E-value of product of p-values"),
    model->min_w, model->max_w, dataset->min_ic);
  if (dataset->ma_adj) {
    fprintf(
      outfile,
      "width:  wg=       %8g    ws=       %8g    endgaps=  %8s\n",
      dataset->wg, dataset->ws, yesno[dataset->endgaps]
    );
  }
  fprintf(outfile,
"nsites: minsites= %8g    maxsites= %8g    wnsites=  %8g\n"
"theta:  prob=     %8g    spmap=    %8s    spfuzz=   %8g\n", 
    model->min_nsites, model->max_nsites, dataset->wnsites,
    dataset->prob, dataset->mapname, dataset->map_scale);

  // print global search parameters
  fprintf(outfile,
"global: substring=%8s    branching=%8s    wbranch=  %8s\n",
    (dataset->p_point->c > 0) ? "no" : "yes", 
    (dataset->branch_params->point_branch != NO_POINT_B) ? "yes" : "no",
    yesno[dataset->branch_params->w_branch]);

  // print branching parameters
  if (dataset->branch_params->point_branch != NO_POINT_B){
    fprintf(outfile,
"        bfactor=  %8i    heapsize= %8i\n",
    dataset->branch_params->bfactor,
    dataset->main_hs);
  }

  // print local search parameters
  fprintf(outfile,
"em:     prior=   %9s    b=        %8g    maxiter=  %8d\n"
"        distance= %8g\n",
    dataset->priorname, dataset->beta, dataset->maxiter,
    dataset->distance);

  // print properties of the dataset
  fprintf(outfile,
"data:   n=        %8d    N=        %8d\n", 
    dataset->total_res, dataset->n_samples);

  if (!strcmp(dataset->alphabet, DNA0)) {
    fprintf(outfile, "strands: +");
    if (model->invcomp) {
      fprintf(outfile, " -");
    }
  }
  fprintf(outfile,
    "\n"
    "sample: seed=     %8d    seqfrac=  %8g\n",
    dataset->seed, dataset->seqfrac
  );
  if (dataset->plib_name) {
    fprintf(outfile, "Dirichlet mixture priors file: %s\n", dataset->plib_name);
  }
  /* print dataset frequencies of letters in alphabet */
  fprintf(outfile, "Letter frequencies in dataset:\n");
  for (i=0, pcol=0; i<dataset->alength; i++) {
    pcol += 8;          /* start of next printed thing */
    if (pcol >= PAGEWIDTH) {pcol=8; fprintf(outfile, "\n");}
    fprintf(outfile, "%c %5.3f ", dataset->alphabet[i], dataset->res_freq[i]);
  }
  /* print background frequencies of letters in alphabet */
  fprintf(
    outfile,
    "\nBackground letter frequencies (from %s):\n",
    dataset->bfile ? dataset->bfile : "dataset with add-one prior applied"
  );
  for (i=0, pcol=0; i<dataset->alength; i++) {
    pcol += 8;          /* start of next printed thing */
    if (pcol >= PAGEWIDTH) {
      pcol=8;
      fprintf(outfile, "\n");
    }
    fprintf(outfile, "%c %5.3f ", dataset->alphabet[i], dataset->back[i]);
  }
  fprintf(outfile, "\n%s\n", stars);

} /* print_command_summary */

/**********************************************************************/
/*
  meme_score_sequence

  Compute the sequence score for a motif.

  Returns the sequence score.
*/
/**********************************************************************/
static double meme_score_sequence(
  char *eseq,   /* integer-coded sequence to score */
  int length,   /* length of the sequence */
  int w,    /* width of motif */
  LOGODDS logodds1,   /* log-odds matrix: LOG2(m_ij/b_j) */
  LOGODDS logodds2    /* log-odds matrix: LOG2(m_ij/n_ij) */
)
{
  int i, j;
  double best = LITTLE;     /* sequence score */
  double score, sc1, sc2;
  double loge2 = log(2);

  /* score the sequence with motif */
  for (i=0; i <= length - w; i++) { /* site start */
    /* calculate score of subsequence */
    for (j=0, sc1=0, sc2=0; j<w; j++) { /* position in sequence */
      sc1 += logodds1(j, (int) eseq[i+j]);
      if (logodds2) sc2 += logodds2(j, (int) eseq[i+j]);
    } /* subsequence */
    score = logodds2 ? -LOGL_SUM(-sc1*loge2, -sc2*loge2)/loge2 : sc1;
    best = MAX(score, best);
  } /* site start */

  return best;

} /* meme_score_sequence */

/**********************************************************************/
/*
  get_thresh

  Get the optimal threshold for minimizing classification error
  by classifying positive and negative data using the motif,
  sorting and finding the minimum error.

  Returns optimal threshold, error rate and ROC.
*/
/**********************************************************************/
static ACCURACY *get_thresh(
  int w,        /* width of motif */
  LOGODDS logodds1,     /* log-odds matrix: LOG2(m_ij/b_j) */
  LOGODDS logodds2,     /* log-odds matrix: LOG2(m_ij/a_ij) */
  DATASET *pos,       /* positive examples */
  DATASET *neg,       /* negative examples */
  BOOLEAN print_scores      /* print sorted scores */
)
{
  int i, class, iseq;
  int err;        /* sum of false pos. and false neg. */
  int min_pos, max_pos;     /* best cutoff index in sorted list */
  int best_err;       /* best classification error */
  double thresh;      /* best threshold */
  SORTED_SCORE *scores=NULL;    /* array of class/score */
  int npos = pos->n_samples;
  int nneg = neg->n_samples;
  int nseqs = npos + nneg;    /* number of sequences */
  DATASET *dataset;
  /* for ROC */
  double roc;       /* receiver operating characteristic */
  double tpp, fpp;      /* true, false positive proportions */
  double tp, fp;      /* true, false positives so far */
  double newtpp, newfpp;
  ACCURACY *acc = NULL;
  double minposscore;     /* minimum score of a positive */
  double maxnegscore;     /* maximum score of a negative */

  /* allocate space for accuracy record */
  Resize(acc, 1, ACCURACY);

  /* allocate space for scores */
  Resize(scores, nseqs, SORTED_SCORE);

  /* score sequences */
  for (class=0, iseq = 0; class<2; class++) {
    if (class) dataset = pos; else dataset = neg;
    for (i=0; i<dataset->n_samples; i++) {
      SAMPLE *s = dataset->samples[i];  /* sample */
      char *eseq = s->res;      /* integer-coded sequence */
      int lseq = s->length;     /* length of sequence */
      if (lseq < w) continue;     /* sequence too short */
      scores[iseq].class = class;
      scores[iseq].score = meme_score_sequence(eseq, lseq, w, logodds1,
        logodds2);
      scores[iseq].id = s->sample_name;
      iseq++;
    }
  }

  /* sort sequences by score in descending order */
  qsort(scores, nseqs, sizeof(SORTED_SCORE), s_compare);

  /*
    get threshold with minimum classification error
    If a range of thresholds give the same error, choose the average.
  */
  roc = tpp = fpp = tp = fp = 0 ; /* for ROC */
  best_err = err = pos->n_samples;  /* all false negatives */
  min_pos = 0;        /* threshold must be below this */
  max_pos = 0;        /* threshold must be above this */
  minposscore = BIG;      /* smallest score of positives */
  maxnegscore = -BIG;     /* smallest score of negatives */
  for (i=0; i<nseqs; i++) {
    if (scores[i].class) {      /* positive */
      err--;          /* one fewer fn */
      tp++;
      minposscore = scores[i].score;
    } else {          /* negative */
      err++;          /* one more fp */
      fp++;
      maxnegscore = MAX(maxnegscore, scores[i].score);
    }
    if (err < best_err) {       /* start new range */
      best_err = err;
      min_pos = max_pos = i;
    } else if (err == best_err) {   /* extend current range */
      max_pos = i;
    }
    /* ROC trapezoidal rule : (y2 - y1)/2 dx */
    newtpp = tp / npos;
    newfpp = fp / nneg;
    roc += .5 * (newtpp + tpp) * (newfpp - fpp);
    tpp = newtpp;
    fpp = newfpp;
  }
  max_pos = MIN(max_pos+1, nseqs-1);
  thresh = (scores[min_pos].score + scores[max_pos].score)/2;

  /* normalize by fpp to get ROC */
  if (fpp == 0) {
    roc = 1.0;
  } else {
    roc /= fpp;
  }

  /* add difference between positives and negatives if ROC is 1.0 */
  if (roc == 1.0) roc += minposscore - maxnegscore;

  /* print the sorted list */
  if (print_scores) {
    printf("ROC= %f\n", roc);
    for (i=0; i<nseqs; i++) printf("%-*.*s %1d %g\n",
      MSN, MSN, scores[i].id, scores[i].class, scores[i].score);
  }

  acc->thresh = thresh;
  acc->err = best_err;
  acc->roc = roc;

  return acc;
} /* get_thresh */

/**********************************************************************/
/*
        s_compare

        Compare two scores in descending order.  Return <0 >0
        if the first is <, > the second.  If they are equal,
        resolves ties by returning <0 if the first has smaller class.
*/
/**********************************************************************/
static int s_compare(
  const void *v1,
  const void *v2
)
{
  const SORTED_SCORE * s1 = (const SORTED_SCORE *) v1;
  const SORTED_SCORE * s2 = (const SORTED_SCORE *) v2;
  double diff = s1->score - s2->score;
  if (diff == 0) diff = (double) (s1->class - s2->class);
  return ((diff > 0) ? -1 : ( (diff < 0) ? 1 : 0) );
} /* s_compare */


/**********************************************************************/
/*
  get_q

  Get the value of q which gives the optimal ROC on the training
  set.
*/
/**********************************************************************/
static double get_q(
  int nsteps,         /* try nsteps from 0 to 1 */
  int window,         /* smoothing window radius */
  int w,          /* width of motif */
  THETA theta,          /* motif theta */
  THETA neg_theta,        /* anti-motif theta */
  double *back,         /* background motif */
  DATASET *dataset,       /* the dataset */
  DATASET *neg_dataset,       /* negative examples */
  char *str_space     /* space for printing strand direction */
)
{
  int i, j;
  double *roc = NULL;       /* array to hold roc(q) */
  double smooth_roc;        /* smoothed value of roc */
  double best_roc;        /* maximum roc */
  double q=0;         /* mixing parameter */
  double incr = 1.0/nsteps;     /* increment for q */
  LOGODDS logodds;        /* motif log-odds matrix */
  ACCURACY *acc;        /* get_thresh struct */
  int alength = dataset->alength;   /* length of alphabet */

  /* create ROC array */
  Resize(roc, nsteps+1, double);

  /*
    get roc for different values of q
  */
  for (i=0; i<=nsteps; i++) {
    q = i * incr;
    logodds = make_log_odds(theta, neg_theta, back, q, w, alength);
    acc = get_thresh(w, logodds, NULL, dataset, neg_dataset, FALSE);
    roc[i] = acc->roc;        /* save roc for this q */
    myfree(acc);        /* free up space */
  } /* get roc */

  /*
    smooth ROC and find q that gives maximum
  */
  best_roc = 0;
  for (i=0; i<=nsteps; i++) {
    double avg = 0;
    int cnt = 0;
    for (j=MAX(0,i-window); j<=MIN(nsteps, i+window); j++) {
      avg += roc[j];
      cnt++;
    }
    smooth_roc = avg/cnt;
    if (smooth_roc > best_roc) {
      best_roc = smooth_roc;
      q = i * incr;
    }
    /*printf("q= %8.3f smoothed_roc= %8.5f\n", i*incr, smooth_roc);*/
  } /* smooth ROC and get max */

  myfree(roc);          /* release space */

  printf("Q= %8.3f ROC= %8.3f\n", q, best_roc);

  return q;
}

/**********************************************************************/
/*
  print_sites

  Print the sites making up the model.
*/
/**********************************************************************/
static void print_sites(
  DATASET *dataset,     /* the dataset */
  MODEL *model,		/* the model */
  int format,		/* 0=BLOCKS; 1=FASTA */
  char *com,		/* comment to append */
  FILE *outfile		/* where to print */
)
{
  int i, j;
  int w = model->w;     /* width of motif */
  P_PROB sites = model->maxima;   /* sites "defining" model */
  int n = model->nsites_dis;    /* number of sites */
  char *ftype = (format==0 ? "BLOCKS" : "FASTA");

  /* print header */
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fputc('\n', outfile);
  fprintf(outfile, "\tMotif %d in %s format%s\n", model->imotif, ftype, com);
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fputc('\n', outfile);

  /* print the sites */
  if (format == 0) fprintf(outfile, "BL   MOTIF %d width=%d seqs=%d\n",
    model->imotif, w, n);
  for (i=0; i<n; i++) {     /* print each site */
    int seqno = sites[i].x;   /* sequence number */
    SAMPLE *s = dataset->samples[seqno];/* sequence */
    BOOLEAN ic = sites[i].ic;   /* strand direction */
    int y = sites[i].y;     /* location of site */
    int off = ic ? s->length-w-y : y; /* - strand offset from rgt. */
    char *res = ic ? s->resic+off : s->res+off;       /* integer sequence */
    double weight = s->sw;    /* sequence weight */
    /*double weight = sites[i].prob;*/

    /* print sequence name and position of site */
    if (format == 0) {      /* BLOCKS format */
      fprintf(outfile, "%-*.*s ( %4d) ", MSN, MSN, s->sample_name, y+1);
    } else {        /* FASTA format */
      fprintf(outfile,">%-*.*s pos %4d\n", MSN, MSN, s->sample_name, y+1);
    }

    /* print site */
    for (j=0; j<w; j++) { fputc(unhash(res[j]), outfile); }
    if (format == 0) {      /* BLOCKS format */
      fprintf(outfile, "  %g ", weight);
    }
    fputc('\n', outfile);
  } /* print each site */
  if (format == 0) {
    fprintf(outfile, "//\n\n");
  }
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile);
  fprintf(outfile, "\n\n");

} /* print_sites */

/**********************************************************************/
/*
  print_summary

  Print the summary of all the motifs.
*/
/**********************************************************************/
extern void print_summary(
  MODEL *model,     /* the model IN */
  DATASET *dataset, /* the dataset IN */
  LO **los,         /* the LO structures IN */
  int nmotifs,      /* number of motifs IN */
  double **pv,      /* p-value of score distribution IN */
  FILE *outfile     /* where to print IN */
)
{
  /* print the motif block diagrams using all the motifs */
  fprintf(outfile, "\n\n%s\nSUMMARY OF MOTIFS\n%s\n\n", stars, stars);
  print_block_diagrams(model, dataset, los, nmotifs, pv, outfile);
  fprintf(outfile, "%s\n\n", stars);
} /* print_summary */

/**********************************************************************
  print_meme_file_xml

  Print MEME results in XML format. See DTD embeded in
  meme-dtd.h for description of MEME document.
 **********************************************************************/
extern void print_meme_file_xml(
  MODEL *model,                   /* the model IN */
  DATASET *dataset,               /* the dataset IN */
  LO *los[MAXG],                  /* logodds structures for motifs */
  int nmotifs,                    /* number of motifs IN */
  MOTIF_SUMMARY *motif_summaries, /* list of final motif properties IN */
  char *stopping_reason,          /* description of reason for stopping IN */
  char* xml_filename              /* full path to output file for xml IN */
)
{
  FILE* outfile = fopen(xml_filename, "w"); //FIXME CEG check for errors
  print_meme_header_xml(outfile);

  const char *archive_date = ARCHIVE_DATE;
  int i = strlen(archive_date);
  fprintf(outfile,
    "<MEME version=\"%s\" release=\"%.*s\">\n",
    VERSION,
    i,
    archive_date
  );
  print_meme_training_set_xml(dataset, nmotifs, outfile);
  print_meme_model_xml(model, dataset, stopping_reason, outfile);
  print_meme_motifs_xml(model, dataset, nmotifs, motif_summaries, outfile);
  if (nmotifs > 0) {
    print_meme_scanned_sites_xml(model, dataset, los, nmotifs, pv, outfile);
  }
  fprintf(outfile, "</MEME>\n");
  fclose(outfile);
} /* print_meme_file_xml */

/**********************************************************************
  print_meme_header_xml

  Print the DTD and XML header for MEME
**********************************************************************/
static void print_meme_header_xml(FILE *outfile) {
  fputs(meme_dts, outfile);
  fputs("<!-- Begin document body -->\n", outfile);
}


/**********************************************************************
  print_meme_training_set_xml

  Print XML elements for the training set. See DTD embeded in
  print_meme_header_xml for description of the training set XML syntax.
**********************************************************************/
static void print_meme_training_set_xml(
  DATASET *dataset, /* the dataset IN */
  int nmotifs,      /* number of motifs IN */
  FILE* outfile     /* file for output IN */
) {
  int i = 0;
  const char *alphabet_id = NULL;
  const char *alphabet = NULL;

  fprintf(
    outfile,
    "<training_set datafile=\"%s\" length=\"%d\">\n",
    dataset->datafile,
    dataset->n_samples
  );
  // Print alphabet
  if (dataset->dna) {
    alphabet_id = "nucleotide";
    alphabet = DNAB;
  }
  else {
    alphabet_id = "amino-acid";
    alphabet = PROTEINB;
  }
  fprintf(outfile,
    "<alphabet id=\"%s\" length=\"%d\">\n",
    alphabet_id,
    dataset->alength
  );
  // Print alphabet
  int c;
  for (i = 0; i < dataset->alength; i++) {
    c = dataset->alphabet[i];
    fprintf(outfile, "<letter id=\"letter_%c\" symbol=\"%c\"/>\n", c, c);
  }
  fprintf(outfile, "</alphabet>\n");
  // Print ambiguous characters
  fprintf(outfile, "<ambigs>\n");
  for (i = 0; (c=alphabet[i]) != 0; i++) {
    if (strchr(dataset->alphabet, c) == 0) {
      if (c == '*') {
        // '*' is not a valid character in an XML ID attribute
        // so we have to spell it out.
        fprintf(outfile, "<letter id=\"letter_star\" symbol=\"%c\"/>\n", c);
      }
      else {
        fprintf(outfile, "<letter id=\"letter_%c\" symbol=\"%c\"/>\n", c, c);
      }
    }
  }
  // FIXME HACK
  // The BLAST DNA alphabet does not include X but X can occur for DNA in the 
  // current implementation of MEME
  // This is the nicest of the possible hacks to get around that problem
  if (dataset->dna) {
    fprintf(outfile, "<letter id=\"letter_X\" symbol=\"X\"/>\n");
  }
  fprintf(outfile, "</ambigs>\n");
  // Print sequences in training set
  for(i = 0; i < dataset->n_samples; i++) {
    SAMPLE *s = dataset->samples[i];
    char *name = s->sample_name;
    fprintf(
      outfile,
      "<sequence "
      "id=\"sequence_%d\" "
      "name=\"%s\" "
      "length=\"%ld\" "
      "weight=\"%f\" "
      "/>\n",
      i,
      name,
      dataset->samples[i]->length,
      dataset->samples[i]->sw
    );
  }
  fprintf(outfile, "<letter_frequencies>\n");
  fprintf(outfile, "<alphabet_array>\n");
  for(i = 0; i < dataset->alength; i++) {
    fprintf(
      outfile,
      "<value letter_id=\"letter_%c\">%4.3f</value>\n",
      dataset->alphabet[i],
      dataset->res_freq[i]
    );
  }
  fprintf(outfile, "</alphabet_array>\n");
  fprintf(outfile, "</letter_frequencies>\n");
  fprintf(outfile, "</training_set>\n");
}

/**********************************************************************
  print_meme_model_xml

  Print XML elements for the model. See DTD embeded in
  print_meme_header_xml for the description of model XML syntax.
***********************************************************************/
static void print_meme_model_xml(
  MODEL *model,          /* the model IN */
  DATASET *dataset,      /* the dataset IN */
  char* stopping_reason, /* reason for stopping IN */
  FILE* outfile          /* output file IN */
) {
  char evt_string[12];

  if (dataset->evt == BIG) {
    strcpy(evt_string, "inf");
  } else {
    sprintf(evt_string, "%8g", dataset->evt);
  }

  #define MAX_HOST_NAME 100
  char hostname[MAX_HOST_NAME];
#ifdef UNIX
  int result = gethostname(hostname, MAX_HOST_NAME);
  if (result != 0) {
    // In most cases this just means that the hostname
    // has been truncated to fit in the buffer,
    // but may also be a complete failure.
    // We assume complete failure.
    snprintf(hostname, MAX_HOST_NAME, "CPU: unknown");
  }
#else
  // gethostname may be unavailable on non-unix systems
  snprintf(hostname, MAX_HOST_NAME, "CPU: unknown");
#endif /* UNIX */

  fprintf(
    outfile,
    "<model>\n"
    "<command_line>%s</command_line>\n"
    "<host>%s</host>\n"
    "<type>%s</type>\n"
    "<nmotifs>%d</nmotifs>\n"
    "<evalue_threshold>%s</evalue_threshold>\n"
    "<object_function>%s</object_function>\n"
    "<min_width>%d</min_width>\n"
    "<max_width>%d</max_width>\n"
    "<minic>%8.2f</minic>\n",
    dataset->command, hostname, dataset->mod, dataset->nmotifs, evt_string,
    (dataset->objfun == Pv ?
       "P-value of product of p-values" :
       "E-value of product of p-values"
    ),
    model->min_w, model->max_w, dataset->min_ic
  );
  if (dataset->ma_adj) {
    fprintf(
      outfile,
      "<wg>%g</wg>\n"
      "<ws>%g</ws>\n"
      "<endgaps>%s</endgaps>\n",
      dataset->wg, dataset->ws, yesno[dataset->endgaps]);
  }
  fprintf(
    outfile,
    "<minsites>%g</minsites>\n"
    "<maxsites>%g</maxsites>\n"
    "<wnsites>%g</wnsites>\n"
    "<prob>%g</prob>\n"
    "<spmap>%s</spmap>\n"
    "<spfuzz>%g</spfuzz>\n"
    "<prior>%s</prior>\n"
    "<beta>%g</beta>\n"
    "<maxiter>%d</maxiter>\n"
    "<distance>%g</distance>\n"
    "<num_sequences>%d</num_sequences>\n"
    "<num_positions>%d</num_positions>\n"
    "<seed>%d</seed>\n"
    "<seqfrac>%8g</seqfrac>\n",
    model->min_nsites, model->max_nsites, dataset->wnsites,
    dataset->prob, dataset->mapname, dataset->map_scale,
    dataset->priorname, dataset->beta, dataset->maxiter, dataset->distance,
    dataset->n_samples, dataset->total_res, dataset->seed, dataset->seqfrac
  );
  fprintf(outfile, "<strands>");
  if (dataset->dna == TRUE) {
    // Nucleotide sequence
    if(model->invcomp) {
      fprintf( outfile, "both");
    }
    else {
      fprintf( outfile, "forward");
    }
  }
  else {
    // Amino acid sequence
    fprintf( outfile, "none");
  }
  fprintf(outfile, "</strands>\n");
  fprintf(outfile, "<priors_file>");
  if (dataset->plib_name) {
    fprintf(outfile, "%s", dataset->plib_name);
  }
  fprintf(outfile, "</priors_file>\n");
  fprintf(
    outfile,
    "<reason_for_stopping>%s</reason_for_stopping>\n",
    stopping_reason
  );
  char *bfile = dataset->bfile ? dataset->bfile
                               : "dataset with add-one prior applied";
  fprintf(outfile, "<background_frequencies source=\"%s\">\n", bfile);
  fprintf(outfile, "<alphabet_array>\n");
  int i;
  for (i = 0; i < dataset->alength; i++) {
    fprintf(
      outfile,
      "<value letter_id=\"letter_%c\">%4.3f</value>\n",
      dataset->alphabet[i],
      dataset->back[i]
    );
  }
  fprintf(outfile, "</alphabet_array>\n");
  fprintf(outfile, "</background_frequencies>\n");
  fprintf(outfile, "</model>\n");
}

/**********************************************************************
  print_meme_motifs_xml

  Print XML elements for the motifs. See DTD embeded in
  print_meme_header_xml for the description of model XML syntax.

***********************************************************************/
static void print_meme_motifs_xml(
  MODEL *model,                   /* the model IN */
  DATASET *dataset,               /* the dataset IN */
  int nmotifs,                    /* number of motifs IN */
  MOTIF_SUMMARY *motif_summaries, /* List of final motif properties IN */
  FILE* outfile                   /* output file IN */
) {
  fprintf(outfile, "<motifs>\n");
  int i = 0;
  for (i = 0; i < nmotifs; i++) {
    fprintf(
      outfile,
      "<motif id=\"motif_%d\" name=\"%d\" width=\"%d\" sites=\"%d\""
      " ic=\"%.1f\" re=\"%.1f\""
      " llr=\"%.0f\" e_value=\"%3.1fe%+04.0f\" bayes_threshold=\"%g\""
      " elapsed_time=\"%f\">\n",
      i + 1,
      i + 1,
      motif_summaries[i].width,
      motif_summaries[i].num_sites,
      motif_summaries[i].ic,
      motif_summaries[i].re,
      motif_summaries[i].llr,
      motif_summaries[i].e_value_mant,
      motif_summaries[i].e_value_exp,
      motif_summaries[i].bayes,
      motif_summaries[i].elapsed_time
    );
    print_meme_pssm_xml(
      motif_summaries[i].pssm,
      dataset->alength,
      dataset->alphabet,
      motif_summaries[i].width,
      outfile
    );
    print_meme_psfm_xml(
      motif_summaries[i].psfm,
      dataset->alength,
      dataset->alphabet,
      motif_summaries[i].width,
      outfile
    );
    print_meme_regular_expression_xml(
      motif_summaries[i].regexp,
      outfile
    );
    print_meme_contributing_sites_xml(
      model,
      &(motif_summaries[i]),
      dataset,
      outfile
    );
    fprintf(outfile, "</motif>\n");
  }
  fprintf(outfile, "</motifs>\n");
}

/**********************************************************************
  print_meme_pssm_xml

  Print XML elements for the log-odds matrix. See DTD embeded in
  print_meme_header_xml for the description of model XML syntax.
**********************************************************************/
static void print_meme_pssm_xml(
  LOGODDS logodds, /* pointer to matrix of log-odds scores IN */
  int alength,     /* length of the alphabet IN */
  char* alphabet,  /* pointer to alphabet string IN */
  int width,       /* width of the motif IN */
  FILE* outfile    /* pointer to output file IN */
) {

  int i = 0;
  int j = 0;

  fprintf(outfile, "<scores>\n<alphabet_matrix>\n");
  for (i=0; i < width; i++) {   /* site position */
    fprintf(outfile, "<alphabet_array>\n");
    for (j=0; j < alength; j++) { /* letter */
      fprintf(outfile,
        "<value letter_id=\"letter_%c\">%d</value>\n",
        alphabet[j],
        NINT(100*logodds(i,j))
      );
    }
    fprintf(outfile, "</alphabet_array>\n");
  }
  fprintf(outfile, "</alphabet_matrix>\n</scores>\n");
}

/**********************************************************************
  print_meme_psfm_xml

  Print XML elements describing the frequency matrix. See DTD embeded in
  print_meme_header_xml for the description of model XML syntax.
**********************************************************************/
static void print_meme_psfm_xml(
  THETA theta,     /* pointer to matrix of frequencies */
  int alength,     /* length of the alphabet */
  char* alphabet,  /* pointer to alphabet string */
  int motif_width, /* width of the motif */
  FILE* outfile    /* pointer to output file */
) {
  int i = 0;
  int j = 0;

  fprintf(outfile, "<probabilities>\n<alphabet_matrix>\n");
  for (i=0; i < motif_width; i++) {   /* site position */
    fprintf(outfile, "<alphabet_array>\n");
    for (j=0; j < alength; j++) { /* letter */
      fprintf(
        outfile,
        "<value letter_id=\"letter_%c\">%f</value>\n",
        alphabet[j],
        theta(i, j)
      );
    }
    fprintf(outfile, "</alphabet_array>\n");
  }
  fprintf(outfile, "</alphabet_matrix>\n</probabilities>\n");
}

/**********************************************************************/
/*
  print_meme_regular_expression_xml

  Prints the XML element describing a motif's  regular expression.
  See DTD embeded in meme-dtd.h for the description of model
  XML syntax.
*/
/**********************************************************************/
static void print_meme_regular_expression_xml(
  char* regexp,     /* regular expression  */
  FILE *outfile
) {

  fprintf(
    outfile,
    "<regular_expression>\n%s\n</regular_expression>\n",
    regexp
  );
}

/**********************************************************************/
/*
  print_meme_contributing_sites_xml

  Print XML elements describing a motif's occurences. See DTD embeded in
  dtd.h for the description of model XML syntax.
*/
/**********************************************************************/
static void print_meme_contributing_sites_xml(
  MODEL *model,
  MOTIF_SUMMARY *motif_summary,
  DATASET *dataset,
  FILE* outfile
) {

    char site[MAXSITE+1], pre[10+1], post[10+1];
    fprintf(outfile, "<contributing_sites>\n");
    int i = 0;
    for(i = 0; i < motif_summary->num_sites; i++) {
      int seqno = motif_summary->sites[i].x;
      SAMPLE *s = dataset->samples[seqno];
      int lseq = s->length;   /* length of sequence */
      char *seq = s->seq;     /* the ascii sequence */
      BOOLEAN ic = motif_summary->sites[i].ic;
      int motif_start = motif_summary->sites[i].y;
      char *strand =  "none";
      if (dataset->dna == TRUE) {
        strand = ic ? "minus" : "plus";
      } 
      /* get the aligned sequence parts */
      get_aligned_sequence_parts(
        motif_start,
        motif_summary->width,
        ic,
        lseq,
        seq,
        pre,
        site,
        post
      );
      fprintf(
        outfile,
        "<contributing_site"
        " sequence_id=\"sequence_%d\""
        " position=\"%d\""
        " strand=\"%s\""
        " pvalue=\"%.2e\""
        " >\n",
        seqno,
        motif_start,
        strand,
        motif_summary->sites[i].prob
      );
      fprintf(outfile, "<left_flank>%s</left_flank>\n<site>\n", pre);
      int j = 0;
      for(j = 0; j < motif_summary->width; j++) {
        fprintf( outfile, "<letter_ref letter_id=\"letter_%c\"/>\n", site[j]);
      }
      fprintf(
        outfile,
        "</site>\n<right_flank>%s</right_flank>\n</contributing_site>\n",
        post
      );
    }
  fprintf(outfile, "</contributing_sites>\n");
}

/**********************************************************************/
/*
get_aligned_sequence_parts

Extract the sequence corresponding to the motif site and 10 bases
  on the left and right flanks.

*/
/**********************************************************************/
void get_aligned_sequence_parts(
  int motif_start,
  int motif_width,
  BOOLEAN ic,
  int lseq,
  char *seq,
  char *pre,
  char *site,
  char *post
) {
  int i = 0;
  int ii = 0;
  if (!ic) {        /* + strand */
    /* pre */
    for (i = motif_start - 10, ii = 0; i < motif_start; i++) {
      if (i<0) {
        continue;
      }
      pre[ii++] = seq[i];
    }
    pre[ii] = '\0';
    /* site */
    for (i = motif_start, ii = 0; ii < motif_width; i++) {
      site[ii++] = seq[i];
    }
    site[ii] = '\0';
    /* post */
    for (i = motif_start + motif_width, ii = 0; ii < 10 && i < lseq; i++) {
      post[ii++] = seq[i];
    }
    post[ii] = 0;
  }
  else {        /* - strand */
    /* pre */
    for (i = motif_start + motif_width + 9, ii = 0; i >= motif_start + motif_width; i--) {
      if (i>=lseq) {
        continue;
      }
      pre[ii++] = comp_dna(seq[i]);
    }
    pre[ii] = '\0';
    /* site */
    for (i = motif_start + motif_width - 1, ii = 0; ii < motif_width; i--) {
      site[ii++] = comp_dna(seq[i]);
    }
    site[ii] = '\0';
    /* post */
    for (i = motif_start - 1, ii = 0; ii < 10 && i >= 0; i--) {
      post[ii++] = comp_dna(seq[i]);
    }
    post[ii] = '\0';
  } /* strand */
}

/**********************************************************************/
/*
  score_sites

  Score and get the pvalues of the sites in the model.
  Sort in order of increasing p-value.
*/
/**********************************************************************/
static void score_sites(
  DATASET *dataset,     /* the dataset */
  MODEL *model,       /* the model */
  LO *lo,       /* LO structure */
  double *pv        /* p-values for scores of this motif */
)
{
  int isite;
  P_PROB sites = model->maxima;   /* sites "defining" model */
  int n = model->nsites_dis;    /* number of sites */
  BOOLEAN invcomp = model->invcomp; /* use reverse complement strand, too */
  STYPE stype = invcomp ? Combine : Protein;  /* Protein works for DNA too */
  SCORE **scores = NULL;    /* the site scores */
  int old_seqno = -1;

  for (isite=0; isite<n; isite++) { /* site */
    int seqno = sites[isite].x;   /* sequence number */
    int y = sites[isite].y;   /* location of site */
    SAMPLE *s = dataset->samples[seqno];/* sequence */
    int lseq = s->length;   /* length of sequence */
    char *seq = s->seq;     /* the ascii sequence */
    double pvalue;      /* score p-value */

    /* score the sequence if new */
    if (old_seqno != seqno) {
      BOOLEAN xlate_dna = FALSE;  /* not xlating */
      if (old_seqno >= 0) free_2array(scores, 1);
      scores = score_sequence(stype, FALSE, xlate_dna, seq, lseq, 1, &lo);
      old_seqno = seqno;
      s->minpv = 1.0;
    }

    pvalue = pv[(int) scores[0][y].score];  /* p-value */

    /* save MINUS the p-value in the .prob field of sites */
    sites[isite].prob = -pvalue;

    /* update minimum p-value of sites */
    if (pvalue < s->minpv) s->minpv = pvalue;

  } /* get p-values of sites */
  free_2array(scores, 1);               /* free space */

  /*
    sort the sites by p-value
  */
  qsort((char *) sites, n, sizeof(p_prob), pY_compare);

  /*
    change sign of p-values back
  */
  for (isite=0; isite<n; isite++) sites[isite].prob *= -1;

} /* score_sites */

/**********************************************************************/
/*
  print_site_diagrams

  Make block diagrams of the actual sites in the model
  and print them.
  Sequences are sorted by the minimum p-value of sites in them.
*/
/**********************************************************************/
static void print_site_diagrams(
  DATASET *dataset,     /* the dataset */
  MODEL *model,       /* the model */
  int nmotifs,        /* number of motifs in los */
  LO *los[MAXG],      /* logodds structure for motif */
  FILE *outfile       /* where to print */
)
{
  int i, j, isite;
  P_PROB sites = model->maxima;   /* sites "defining" model */
  int n = model->nsites_dis;    /* number of sites */
  int nseqs = dataset->n_samples; /* number of sequences in dataset */
  BOOLEAN dna = dataset->dna;   /* dataset is DNA if true */
  BOOLEAN invcomp = model->invcomp; /* use reverse complement strand, too */
  BOOLEAN xlate_dna = FALSE;    /* don't translate */
  BOOLEAN best_motifs = FALSE;    /* use all sites */
  double m_thresh = 1;      /* show all sites as strong */
  STYPE stype = dna ? (invcomp ? Combine : Norc) : Protein;
  int nseqs_with_sites;     /* number of sequences with sites */
  int *seqlist = NULL;      /* ordered list of sequences w/sites */
  int *hits = NULL;     /* store hits */
  double *pvalues = NULL;   /* store pvalues */
  char *f = "%-*s%s %8s  %s\n";   /* format */

  /*
    Create the list to contain sequences sorted by minimum p-value
  */
  Resize(seqlist, nseqs, int);

  /*
    Clear list of sites for each sequence
  */
  for (i=0; i<nseqs; i++) dataset->samples[i]->nsites = 0;

  /*
    Find which sequences have sites and create list of sites for each
  */
  for (isite=nseqs_with_sites=0; isite<n; isite++) {  /* site */
    int seqno = sites[isite].x;   /* sequence number */
    int y = sites[isite].y + 1;   /* location of site (plus 1) */
    int ic = sites[isite].ic;   /* site on reverse complement strand */
    SAMPLE *s = dataset->samples[seqno];/* sequence */

    /* record the sequence as containing a site */
    if (!s->nsites) seqlist[nseqs_with_sites++] = seqno;

    /* store the site in its list of sites */
    Resize(s->sites, s->nsites+1, int);
    s->sites[(s->nsites)++] = ic ? -y : y;  /* +/- site offset by 1 */
  } /* site */

  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fputc('\n', outfile);
  fprintf(outfile, "\tMotif %d block diagrams\n", nmotifs);
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fputc('\n', outfile);
  fprintf(outfile, f, MSN, "SEQUENCE NAME", "", "POSITION P-VALUE", "MOTIF DIAGRAM");
  fprintf(outfile, f, MSN, "-------------", "", "----------------", "-------------");

  /*
    create and print a block diagram for each sequence with sites
  */
  for (i=0; i<nseqs_with_sites; i++) {  /* sequence */
    int seqno = seqlist[i];   /* current sequence */
    SAMPLE *s = dataset->samples[seqno];/* sequence */
    int lseq = s->length;   /* length of sequence */
    char *name = s->sample_name;  /* name of sequence */
    double minpv = s->minpv;    /* minimum p-value of sites */
    char hdr[MSN+20];     /* header */
    TILING tiling;      /* tiling struct */
    tiling.diagram = NULL;	// silence compiler warnings

    /*
      create storage for hits and pvalues and clear them
    */
    Resize(hits, lseq, int);
    Resize(pvalues, lseq, double);
    for (j=0; j<lseq; j++) { hits[j] = 0; pvalues[j] = 0; }

    /* copy hits from s->nsites into hits array */
    for (j=0; j<s->nsites; j++) {
      int y = abs(s->sites[j]) - 1; /* position of site */
      int m = (s->sites[j] > 0) ? los[nmotifs-1]->imotif : -los[nmotifs-1]->imotif;
      hits[y] = m;      /* +/- motif */
    }

    /* put the hits in TILING struct */
    tiling.diagram = NULL;	// to prevent compiler warning
    tiling.hits = hits;
    tiling.pvalues = pvalues;

    /* create the block diagram */
    tiling.diagram = create_diagram(dna, stype, xlate_dna, best_motifs, FALSE,
      m_thresh, nmotifs, los, lseq, FALSE, NULL, 0, 0, tiling);

    /* print the diagram */
    sprintf(hdr, "%-*.*s %16.2g  ", MSN, MSN, name, minpv);
    print_diagram(tiling.diagram, hdr, outfile);

    myfree(tiling.diagram);   /* release space */
  } /* sequence */

  /* print a final line of hyphens */
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fprintf(outfile, "\n\n");

  myfree(seqlist);

} /* print_site_diagrams */

/**********************************************************************/
/*
  align_sites

  Align all sites that make up the model.
*/
/**********************************************************************/
static void align_sites(
  DATASET *dataset,     /* the dataset */
  MODEL *model,       /* the model */
  LO *lo,       /* LO structure */
  double *pv,       /* pvalues for scores of this motif */
  FILE *outfile       /* stream for output */
)
{
  int i, ii, isite;
  int w = model->w;     /* length of site */
  P_PROB sites = model->maxima;   /* sites "defining" model */
  int n = model->nsites_dis;    /* number of sites */
  BOOLEAN invcomp = model->invcomp; /* use reverse complement strand, too */
  int imotif = lo->imotif;    /* name of motif */
  char site[MAXSITE+1], pre[10+1], post[10+1];

  /* print header */
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fputc('\n', outfile);
  fprintf(outfile,
    "\tMotif %d sites sorted by position p-value\n", imotif);
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fputc('\n', outfile);
  fprintf(outfile, "%-*.*s%s ", MSN, MSN, "Sequence name",
    invcomp ? " Strand" : "");
  fprintf(outfile, "%6s %9s %10s %*sSite%*s\n",
    "Start", "P-value", "", w/2 - 2, "", w - w/2 - 4, "");
  fprintf(outfile, "%-*.*s%s ", MSN, MSN, "-------------",
    invcomp ? " ------" : "");
  fprintf(outfile, "%6s %8s %10s ", "-----", "---------", "");
  for (i=0; i<w; i++) fputc('-', outfile);
  fputc('\n', outfile);

  /*
    print sites that make up the model
  */
  for (isite=0; isite<n; isite++) { /* site */
    int seqno = sites[isite].x;   /* sequence number */
    int y = sites[isite].y;   /* location of site */
    BOOLEAN ic = sites[isite].ic; /* strand direction */
    double pvalue = sites[isite].prob;  /* position p-value */
    SAMPLE *s = dataset->samples[seqno];/* sequence */
    int lseq = s->length;   /* length of sequence */
    char *seq = s->seq;     /* the ascii sequence */
    char *sample_name = s->sample_name; /* name of sample */

    /* print name and strand */
    fprintf(outfile, "%-*.*s%s ", MSN, MSN, sample_name,
      invcomp ? (ic ? "     -" : "     +") : "");

    /* print position and position p-value */
    fprintf(outfile, "%6d %9.2e", y+1, pvalue);

    /* get the aligned sequence parts */
    if (!ic) {        /* + strand */
      /* pre */
      for (i=y-10, ii=0; i<y; i++) {
        if (i<0) continue;
        pre[ii++] = seq[i];
      }
      pre[ii] = '\0';
      /* site */
      for (i=y, ii=0; ii<w; i++)  site[ii++] = seq[i];
      site[ii] = '\0';
      /* post */
      for (i=y+w, ii=0; ii<10 && i<lseq; i++) post[ii++] = seq[i];
      post[ii] = 0;

    } else {        /* - strand */
      /* pre */
      for (i=y+w+9, ii=0; i>=y+w; i--) {
        if (i>=lseq) continue;
        pre[ii++] = comp_dna(seq[i]);
      }
      pre[ii] = '\0';
      /* site */
      for (i=y+w-1, ii=0; ii<w; i--) site[ii++] = comp_dna(seq[i]);
      site[ii] = '\0';
      /* post */
      for (i=y-1, ii=0; ii<10 && i>=0; i--) post[ii++] = comp_dna(seq[i]);
      post[ii] = '\0';
    } /* strand */

    /* print the alignment */
    if (pre[0] == '\0') {     /* print a dot in empty pre */
      fprintf(outfile, " %10s %-*s %-10s\n", ".", w, site, post);
    } else {
      fprintf(outfile, " %10s %-*s %-10s\n", pre, w, site, post);
    }

  } /* site */

  /* print line of hyphens */
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fprintf(outfile, "\n\n");

} /* align_sites */


/**
 * print_pllr
 *
 * Prints the specified LLR value for the alignment of the planted sites.
 */
static void print_pllr (
  double pllr,       ///< LLR of the planted site alignment
  FILE *outfile      ///< The stream for output
)
{
  fprintf(outfile, "\nPLANTED LLR = %0.4f\n", pllr);
}


/**
 * calc_pllr
 *
 * Compute the llr of the alignment of the planted sites. Stores the result
 * in dataset->pllr.
 */
static void calc_pllr (
  DATASET *dataset,  ///< The dataset of sequences containing the planted sites
  MODEL *model       ///< The model resulting from meme. The length of the
                     ///< planted sites is held here.
)
{
  // Generate the "alignment" of planted sites:
  int n_planted = dataset->n_samples;
  P_PROB planted_sites = (P_PROB) mymalloc(n_planted * sizeof(p_prob));
  int seq_idx;
  for (seq_idx = 0; seq_idx < n_planted; seq_idx++) {
    planted_sites[seq_idx].x = seq_idx;
    
    // Retrieve the expected position of the single site in the current seq:
    SAMPLE *curr_seq = dataset->samples[seq_idx];
    char *descript = curr_seq->descript;
    int kpos = get_first_siteloc(descript);

    planted_sites[seq_idx].y = kpos;
    planted_sites[seq_idx].ic = 0; // Fn currently only works without ic
  }

  // Use align_top_subsequences to calculate the LLR of that sequence alignment:
  int iseq = -1;     // Dummy
  int ioff = -1;     // Dummy
  char *eseq = NULL; // Dummy
  char *name = NULL; // Dummy
  int n_nsites0 = 1; // The s_points array will be of length 1.
  double col_scores[MAXSITE]; // Dummy
  // An array of length 1, used simply to retrieve the desired score.
  S_POINT s_points[1]; 
  // Force align_top_subsequences to set the score of the s_point:
  s_points[0].score = LITTLE;
  // The number of sequences in the alignment also needs to be specified:
  s_points[0].nsites0 = n_planted;
  s_points[0].evaluate = TRUE;
  /* Use align_top_subsequences (with dummy arguments) to calculate desired
   * score: */
  MTYPE mtype = Oops;// There is one planted site in each sequence.
  /* For the llr of the planted sites to be calculated, the user must have
   * specified the length w of the planted sites:
   */
  assert(model->min_w == model->max_w);
  int w = model->min_w;

  // Ensure that LLR is used within ats():
  BOOLEAN orig_use_llr = dataset->use_llr;
  dataset->use_llr = TRUE;

  align_top_subsequences(
    mtype,
    w,
    dataset,
    iseq,
    ioff,
    eseq,
    name,
    n_nsites0,
    n_planted,
    planted_sites,
    col_scores,
    s_points
  );
  dataset->use_llr = orig_use_llr;

  double pllr = s_points[0].score;
  dataset->pllr = pllr;
} // calc_pllr
