/* WidiPadi is a motif discovery tool for DNA sequences */
/* Copyright (C) 2006 Bioinformatics Centre */
/* Author: Eivind Valen and Ole Winther */

/* This program is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU General Public License */
/* as published by the Free Software Foundation; either version 2 */
/* of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program; if not, write to the Free Software */
/* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. */

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <gsl/gsl_sf_gamma.h>

#include "error.h"
#include "sesa.h"

#define USAGE "Usage: widipadi [OPTIONS] SEQUENCES"

enum bool {false, true};

/* FUNCTION PROTOTYPES */
PSSM hits_to_count_matrix(SESA sesa, struct HitTable *hits, const double threshold, int length, int alphabetSize);
int compare_hit_scores(const void *hit1, const void *hit2);
int parse_options(int argc, char **argv);
void print_help();

/* Acceptance */
int accept_all(SESA sesa, struct HitTable *hits);
int acceptIf_over_300_hits(SESA sesa, struct HitTable *hits);
int acceptIf_required_presence(SESA sesa, struct HitTable *hits);


/* DEFAULTS */
unsigned int length = 8;
unsigned int alphabetSize = 4;
unsigned int mismatches = 1;
int (*acceptanceFunc)(SESA, struct HitTable *) = acceptIf_required_presence;  
double required_presence = 0.7;
unsigned int exhaustive = false;

/**
 * Prints a help message to stdout.
 */
void print_help() {
  printf(USAGE);
  printf("\nOptions:\n");
/*   printf("  -a                 Size of alphabet\n"); */
  printf("  -h                 Prints this help and exits\n");
  printf("  -l LENGTH          Length of WM\n");
  printf("  -m MISMATCHES      Number of mismatches\n");
  printf("  -r RATIO           Sets the ratio of sequences where the PSSM must match\n");
  printf("  -x                 Do exhaustive search\n");
  printf("\n");
}


/**
 * Parses command line options
 */
int parse_options(int argc, char **argv) {
  int count = 1;
  int c;

  while ((c = getopt(argc, argv, "hl:m:r:x")) != -1) {
    count++;
    switch (c) {
    case 'h': 
      print_help(); 
      exit(EXIT_SUCCESS);
    case 'l':
      length = (unsigned int) strtol(optarg, (char **) NULL, 0);
      count++;
      break;
    case 'm':
      mismatches = (unsigned int) strtol(optarg, (char **) NULL, 0);
      count++;
      break;
    case 'r':
      acceptanceFunc = acceptIf_required_presence;  
      required_presence = strtod(optarg, (char **) NULL);
      count++;
      break;
    case 'x':
      exhaustive = true;
      break;
    default: 
      fprintf(stderr, USAGE);
      exit(EXIT_FAILURE);
    }
  }

  return count;
}


/*
 * Converts all hits above scoring threshold to a count matrix
 *
 * Only 0th order
 */
PSSM hits_to_count_matrix(SESA sesa, struct HitTable *hits, const double threshold, int length, int alphabetSize) {
  PSSM pssm = initMatrix(0, length, alphabetSize);
  ESA esa = sesa->esa;
  struct HitEntry *hit;
  int i, j;

  i = pssm->offsets[pssm->length];
  while (i--) pssm->counts[i] = 0;

  for (i = 0; i < hits->nScores; i++) {
    hit = &hits->pScores[i];

    if (hit->score < threshold) 
      break;

    for (j = 0; j < length; j++) {
      int ltr = esa->pStr[hit->position + j];
      /*       fprintf(stderr, "%c", ind2chr(ltr)); */
      pssm->counts[ pssm->offsets[j] + ltr ]++;
    }
/*     fprintf(stderr, " (%f)\n", hit->score); */
  }

  return pssm;
}


/**
 * Comparing function for sort. Sorts on scores from low to high.
 */
int compare_hit_scores(const void *hit1, const void *hit2) {
  return ((struct HitEntry *)hit2)->score - ((struct HitEntry *)hit1)->score;
}
 

/**
 * Accepts all PSSMs (for testing)
 */
int accept_all(SESA sesa, struct HitTable *hits) {
  return true;
}


/**
 * Accepts PSSMs with more than 300 hits (for testing)
 */
int acceptIf_over_300_hits(SESA sesa, struct HitTable *hits) {
  return (hits->nScores > 300 ? true : false);
}


/**
 * Accepts PSSMs with hits present in at least "required_presence" of
 * the sequences
 */
int acceptIf_required_presence(SESA sesa, struct HitTable *hits) {
  int present[sesa->nSeq];
  int i, j;

  i = sesa->nSeq;
  while(i--) present[i] = 0;
  
  i = hits->nScores;
  while (i--) present[ sesa->seq[ hits->pScores[i].position ] ]++;

  j = 0;
  i = sesa->nSeq;
  while (i--) if (present[i]) j++;


  return ( ((double) j / (double) sesa->nSeq)  > required_presence ? true : false);
}

/*                                                                                                
 * Log likelihood scoring
 */
double log_likelihood(ESA esa, PSSM pssm, int *count_bg_all, double *pcount, double *pcount_bg) {
  int pcount_tot = 0, pcount_bg_tot = 0, count_tot = 0, count_bg_tot = 0;
  int alphlen = pssm->alphabetSize;
  int count_bg[alphlen];
  int ccount, ltr, pos, ord, i;
  double ll = 0.0;

  i = alphlen;
  while (i--) count_bg[i] = count_bg_all[i];

  pos = pssm->length;
  while (pos--) {
    ltr = alphlen;
    while(ltr--) {
      ccount = pssm->counts[pos * alphlen + ltr];
      count_bg[ltr] -= ccount ;
      ll += gsl_sf_lngamma( ccount + pcount[ltr] );    
    }
  }

  ltr = alphlen;
  while (ltr--) {
    pcount_tot += pcount[ltr];
    pcount_bg_tot += pcount_bg[ltr];
    count_tot += pssm->counts[ltr];   
    count_bg_tot += count_bg[ltr];
    ll += gsl_sf_lngamma( count_bg[ltr] + pcount_bg[ltr] ); 
  }


  ll -= alphlen * gsl_sf_lngamma( count_tot + pcount_tot ) + gsl_sf_lngamma(count_bg_tot + pcount_bg_tot );
  
  return ll;
}



/*                                                                                                
 * Log likelihood scoring
 */
double higher_log_likelihood(ESA esa, PSSM pssm, double *pcount, double *pcount_bg) {
  int pos, lcpi, skip, curSuf, i, j, p;
  int numSufs = getSize(esa);
  int curEsaIndex = 0;
  int order = 2;
  int wcount[order+1];

  i = order+1;
  while (i--) wcount[i] = 0;

  /* Iterate all suffixes */
  while(curEsaIndex < numSufs){
    lcpi = getLcp(esa, curEsaIndex);
    curSuf = getSuf(esa, curEsaIndex);
    skip = getSkip(esa, curEsaIndex);


    /* We can skip if the lcp is longer than the order */
    if (lcpi > order) {

      /* FIX: HACK */
      if (skip == 0) {
	curEsaIndex++;
	continue;
      }

      /* Count all the skipped words */
      i = order+1;
      while (i--) wcount[i] += skip - curEsaIndex;

      curEsaIndex = skip;
      continue;

    } else {

      pos = pssm->length;
      while (pos--) { }

      /* Print the hits */
      p = 1;
      for (j = 0; j <= order; j++)
	if (ind2chr(esa->pStr[getSuf(esa, curEsaIndex-1)+j]) == '-')
	  p = 0;

      if (p) {
	for (i = 0; i <= order; i++) fprintf(stderr, "%c", ind2chr(esa->pStr[getSuf(esa, curEsaIndex-1)+i]));
	fprintf(stderr, ": %5i\n", wcount[2]);
      }

      wcount[order] = 1;

      /* For higher order models the first columns have fewer letters */
      i = order;
      while (i--) {
	if (i >= lcpi) {	  

	  p = 1;
	  for (j = 0; j <= i; j++)
	    if (ind2chr(esa->pStr[getSuf(esa, curEsaIndex-1)+j]) == '-')
		p = 0;

	  if (p) {
	    for (j = 0; j <= i; j++) fprintf(stderr, "%c", ind2chr(esa->pStr[getSuf(esa, curEsaIndex-1)+j]));
	    for (j = i; j < order; j++) fprintf(stderr, " ");
	    fprintf(stderr, ": %5i\n", wcount[i]);
	  }

	  wcount[i] = 1;
	}
      }


      curEsaIndex++;
    }

  }
  
  exit(0);
}

void print_esa (ESA esa) {
  int numSufs = getSize(esa);
  int i,j;

  for (i = 0; i < numSufs; i++) {
    fprintf(stderr, "%4i: SUF(%3i) LCP(%3i) SKIP(%4i) ", i, getSuf(esa, i), getLcp(esa, i), getSkip(esa, i));

    for (j = 0; j < 8; j++) {
      fprintf(stderr, "%c", ind2chr(esa->pStr[getSuf(esa, i)+j]));
    }
    fprintf(stderr, "\n");
  }
}

int main(int argc, char **argv) {
  PSSMScoreTable results;
  struct HitTable *hits = NULL;
  PSSM pssm, count;
  SESA sesa;
  unsigned int arg_count, j;
  int i;

  arg_count = parse_options(argc, argv);

  printf("Reading file: %s  (%i)\n", argv[arg_count], argc - arg_count);
  if (!argv[arg_count]) {
    fprintf(stderr, "Need an input file\n");
    exit(EXIT_FAILURE);
  }
  
  sesa = build_SESA(argc - arg_count, &argv[arg_count], "ACGT", "N", 0);

  if (!sesa) {
    fprintf(stderr, "%s\n", getError());
    exit(EXIT_FAILURE);
  }

  
  printf("Searching %s for motifs of length %i with %i mismatch%s\n", (exhaustive ? "exhaustively" : "exclusively"), length, mismatches, (mismatches > 1 ? "es" : ""));

  if (exhaustive) {
    /* Exhaustive search searches for all words og length n.  */
    results = SESA_exhaustive_search(sesa, length, mismatches, alphabetSize, acceptanceFunc);
  } else {
    /* Exclusive search only searches for PSSMs representing words present in the sequence.  */
    results = SESA_exclusive_search(sesa, length, mismatches, alphabetSize, acceptanceFunc);
  }

  /** 
   * Both the above methods return a list of integers that uniquely
   * identifies pssms.  We do not return the hits of the PSSMs even
   * though we have searched with them since this might be very(!)
   * memory demanding, depending on criteria used for acceptance.
   *
   * If we know the criteria is strict we might considering altering
   * this to avoid re-searching and thus speed up the process.
   */

  pssm = initMatrix(0, length, alphabetSize);

  for (i = 0; i < results->nPssm; i++) {
    /* Sets the pssm to the values corresponding to the ID*/
    set_mismatch_scores(pssm, results->pssmNr[i]);    
    calcAndSetThresholds(pssm, -((double)mismatches)-0.5);    

    /* Search with the PSSM again */
    init_hittable();                 /* Creates a hittable */
    hits = SESA_search(sesa, pssm);  /* Searches again */

    /* Sort the hits */
    qsort((void*)hits->pScores, hits->nScores, sizeof(struct HitEntry), compare_hit_scores);

    printf("\n\n");
    print_pssm(pssm);
    for (j = 0; j <= mismatches; j++) {
      printf("Count matrix for %i mismatches\n", j);
      count = hits_to_count_matrix(sesa, hits, -((double)j)-0.5, pssm->length, pssm->alphabetSize);

      double ost1[] = {1.0, 1.0, 1.0, 1.0};
      double ost2[] = {1.0, 1.0, 1.0, 1.0};

/*       double foo = log_likelihood(count, sesa->bg, ost1, ost2); */

      double foo = higher_log_likelihood(sesa->esa, count, ost1, ost2); 

/*       fprintf(stderr, "LL: %g\n", foo); */
/*       print_counts(count);  */
      printf("\n\n");
    }
    
    release_hits(hits); /* Frees the allocated hittable*/
  }


  exit(EXIT_SUCCESS);
}
