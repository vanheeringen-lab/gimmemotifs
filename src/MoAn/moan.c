/* MoAn is a motif discovery tool for DNA sequences */
/* Copyright (C) 2006 Eivind Valen */

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
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <math.h>
#include "error.h"
#include "pssm.h"
#include "search.h"
#include "score.h"
#include "anneal.h"
#include "moan.h"

#define USAGE "Usage: moan [options] POSITIVE NEGATIVE\n"
/* #define SINGLE score_single_likelihood */
#define SINGLE score_single_likelihood
#define DOUBLE score_single_likelihood
#define COOC score_cooc_likelihood

/* Default priors */
#define POSCOOC_C 0.6
#define POSONE_C  0.2
#define POSNONE_C 0.0
#define NEGCOOC_C 0.0
#define NEGONE_C  0.1
#define NEGNONE_C 0.8

#define POSCOOC_S 0.0
#define POSONE_S  0.95
#define POSNONE_S 0.05
#define NEGCOOC_S 0.0
#define NEGONE_S  0.2
#define NEGNONE_S 0.8

#define MAXIT 30000000

/* Globals */
char *wm_filename;
int wmcount; 
float cut;

FILE *open_file(char *filename);

FILE *open_file(char *filename) {
  FILE *file = (FILE *)NULL;

  if ( !filename || !(file = fopen(filename,"r"))){
    fprintf(stderr, "Error: Couldn't open file %s\n", filename);
    exit(EXIT_FAILURE);
  }

  return file;
}


void print_header(Annealing *search) {

  if (search->loglevel >= NORMAL) {
    printf("======================================================================\n");
    printf("-------------------------------- MoAn --------------------------------\n");
    printf("======================================================================\n");
    printf("Version: %s\n", VERSION);
    printf("Strands: %s\n", ( search->ds ? "double" : "single"));
    printf("WMs: %i\n", search->pssm_count);
    printf("Co-occurrence: %s\n", ( search->cooc ? "yes" : "no"));
    printf("Motif length: %i - %i\n", search->min_len, search->max_len);
    printf("Positive set: %i sequences \n", search->sesa->nSeqSet[1]); 
    printf("Negative set: %i sequences \n", search->sesa->nSeqSet[0]); 
    printf("Number of iterations: %i \n", search->maxit); 
  }

  if (search->loglevel >= VERBOSE) {
    printf("---------------------------------\n");
    printf("Priors     NONE   SINGLE    COOC \n");
    printf("---------------------------------\n");
    printf("Positive: %6.4f  %6.4f   %6.4f\n", search->posnone, search->posone, search->poscooc);
    printf("Negative: %6.4f  %6.4f   %6.4f\n", search->negnone, search->negone, search->negcooc);
    printf("---------------------------------\n");
    printf("End temperature: %f\n", END_TEMP);
    printf("Move probabilities:\n");
    printf("MOD : %5.2f\n", 1.0 - PROB_CUT);
    printf("INC : %5.2f\n", PROB_INC_R);
    printf("DEC : %5.2f\n", PROB_DEC_R - PROB_INC_R);
    printf("MOV : %5.2f\n", PROB_MOV_R - PROB_DEC_R);
    printf("CUT : %5.2f\n", PROB_CUT - PROB_MOV_R);

    if (search->ctype == ANNEALED) {
      fprintf(stderr, "Cutoff : %4.2f - %4.2f\n", CUT_MIN, CUT_MAX);
    }
  }
}

void report(Annealing *search) {
  unsigned int i;

  printf("\n======================================================================\n");
  printf("------------------------------- REPORT -------------------------------\n");
  printf("======================================================================\n");
  printf("Score: %f\n", search->pos_score + search->neg_score);
  printf("Positive single hits:   %i\n", search->pos_matches);
  printf("Positive coocurrences:  %i\n", search->pos_cooc_matches);
  printf("Negative single hits:   %i\n", search->neg_matches);
  printf("Negative coocurrences:  %i\n", search->neg_cooc_matches);
  printf("======================================================================\n\n");

  for (i = 0; i < search->pssm_count; i++) {
    printf("Position Specific Score Matrix no. %i   Cutoff: %f\n", i+1, (search->cut[i] * (float)WM(i)->length));
    printf("----------------------------------------------------------------------\n");
    print_pssm(WM(i));
    printf("----------------------------------------------------------------------\n");
    print_counts(WM(i));
    printf("----------------------------------------------------------------------\n\n");
    printf("Hits for PSSM %i : %i\n", i+1, HITS(i)->nScores);
    printf("----------------------------------------------------------------------\n");
    print_pos_hits(search->sesa, HITS(i), WM(i)->length);
    printf("----------------------------------------------------------------------\n\n");
  }

}

void help() {
  printf(USAGE);
  printf("\nOptions:\n");
  printf("  -c                 Search for co-occurrence of two motifs\n");
  printf("  -D                 Search both strands\n");
  printf("  -h                 Prints this help and exits\n");
  printf("  -i  ITERATIONS     Number of iterations\n");
  printf("  -l  LOGLEVEL       Change amount of output 0-4 (Default: 1)\n");
  printf("  -R  RANGE          Range of the size of the WM, i.e: 5,15 (Default: %i,%i)\n", MIN_LENGTH, MAX_LENGTH);
  printf("\n\n");
}

int parse_options(int argc, char **argv, Annealing *search) {
  int posnone, posone, negnone, negone, tmp;
  int pos_set = 0, neg_set = 0;
  int c, w = 0, count = 1;
  long int s = time(NULL) * getpid();  /* PID because when running multiple processes they are sometimes executed the same second */ 
  unsigned int seed = (unsigned int) s;

  /* Defaults */
  search->ctype      =  ANNEALED;
  search->end_temp   =  END_TEMP; 
  search->start_temp =  -1;
  search->step_ratio =  -1;
  search->loglevel   =  NORMAL;
  search->maxit      =  MAXIT;
  search->schedule   =  EXPONENTIAL;
  search->scorefunc  =  SINGLE;
  search->stepfunc   =  wm_steal;
  search->max_len    =  MAX_LENGTH;
  search->min_len    =  MIN_LENGTH;
  search->ds         =  FALSE;
  search->cooc       =  FALSE;
  cut                =  CUT_MIN + CUT_INT / 2;
  wmcount            =  1;


  /* Dafault priors */
  search->posnone = POSNONE_S;
  search->posone  = POSONE_S;
  search->poscooc = POSCOOC_S;
  search->negnone = NEGNONE_S;
  search->negone  = NEGONE_S;
  search->negcooc = NEGCOOC_S;


  while ((c = getopt(argc, argv, "a:AcC:De:E:f:hi:l:N:p:P:s:St:r:R:w:W:")) != -1) {
    switch (c) {
    case 'a': 
      switch (optarg[0]) {
      case 'E': search->schedule = EXPONENTIAL; break;
      case 'L': search->schedule = LINEAR; break;
      default: 	
	fprintf(stderr, "Invalid argument to a: must be E or L\n");
	exit(EXIT_FAILURE);
      }
      count++;
      break;
    case 'A': search->ctype = ANNEALED; cut = CUT_MIN + CUT_INT / 2; break;
/*     case 'c': search->ctype = RATIO; cut = (float) strtod(optarg, (char **) NULL); count++;  */
/*       if (cut <= 0 || cut > 2) { */
/* 	fprintf(stderr, "Invalid argument to c: must be in range <0,2]\n"); */
/* 	exit(EXIT_FAILURE); */
/*       } */
/*       break; */
    case 'c': 
      search->cooc = TRUE;
      wmcount = 2;
      search->scorefunc  =  COOC; 
      break;
    case 'C': search->ctype = DYNAMIC; cut = (float) strtod(optarg, (char **) NULL); count++; 
      if (cut <= 0 || cut  > 1) {
	fprintf(stderr, "Invalid argument to C: must be in range <0,1]\n");
	exit(EXIT_FAILURE);
      }
      break;
    case 'D': 
      search->ds = TRUE;
      break;
    case 'e': search->end_temp = (float) strtod(optarg, (char **)NULL); count++; break;
    case 'E': seed =  (unsigned int) strtol(optarg, (char **) NULL, 0); count++; break;
    case 'f': search->ctype = FIXED; cut = (float) strtod(optarg, (char **) NULL); count++; break;
    case 'h': help(); exit(EXIT_SUCCESS);
    case 'i': search->maxit = (unsigned int) strtol(optarg, (char **) NULL, 0); count++; break;
    case 'l': search->loglevel =  (unsigned int) strtol(optarg, (char **) NULL, 0); count++; break;
    case 'p': 
      switch (optarg[0]) {
      case 'S': search->stepfunc = wm_steal; break;
      default: 
	fprintf(stderr, "Unrecognized step function: %s\n", optarg);
	exit(EXIT_FAILURE);
      }
      count++;
      break; 
    case 'P':
      tmp = sscanf(optarg, "%i,%i", &posnone, &posone);
      if (tmp == 2) {
	search->posnone = (float) posnone / 100.0;
	search->posone =  (float) posone / 100.0;
	search->poscooc = 1.0 - search->posnone - search->posone;
      } else {
	fprintf(stderr, "Unrecognized prior format: %s\n", optarg);
	exit(EXIT_FAILURE);
      }
      pos_set = 1;
      count++;
      break;
    case 'N':
      tmp = sscanf(optarg, "%i,%i",  &negnone, &negone);
      if (tmp == 2) {
	search->negnone = (float) negnone / 100.0;
	search->negone =  (float) negone / 100.0;
	search->negcooc = 1.0 - search->negnone - search->negone;
      } else {
	fprintf(stderr, "Unrecognized prior format: %s\n", optarg);
	exit(EXIT_FAILURE);
      }
      neg_set = 1;
      count++;
      break;
    case 'r': search->step_ratio = (float) strtod(optarg, (char **) NULL); count++; break;
    case 'R': 
      sscanf(optarg, "%i,%i", &search->min_len, &search->max_len);
      count++;
      break;
    case 'S':
      search->shared_cutoff = TRUE;
      break;
    case 't': search->start_temp = (float) strtod(optarg, (char **) NULL); count++; break;
    case 'w': 
      if(w++ > 0 ) {
	fprintf(stderr, "Can't specify both -w and -W\n"); 
	exit(EXIT_FAILURE); 
      } 

      wmcount = (unsigned int) strtol(optarg, (char **) NULL, 0); 

/*       fprintf(stderr, "WMCOUNT (%i) \n", wmcount); */

      switch(wmcount) {
      case 1:
      fprintf(stderr, "SINGLE \n");
	search->scorefunc  =  SINGLE; 

	if (!pos_set) {
	  search->posnone = POSNONE_S;
	  search->posone  = POSONE_S;
	  search->poscooc = POSCOOC_S;
	} else {
	  
	}

	if (!neg_set) {
	  search->negnone = NEGNONE_S;
	  search->negone  = NEGONE_S;
	  search->negcooc = NEGCOOC_S;
	}

	break;
      default:
	search->scorefunc  =  DOUBLE; 
	break;
/* 	fprintf(stderr, "Number of WMs must be either 1 or 2\n");  */
/* 	exit(EXIT_FAILURE); */
      }

      count++; 
      break;
    case 'W': 
      if(w++ > 0 ) {
	fprintf(stderr, "Can't specify both -w and -W\n"); 
	exit(EXIT_FAILURE); 
      } 
      wm_filename = (char *) malloc(sizeof(char) * (strlen(optarg)+1));
      strcpy(wm_filename, optarg);
      count++;
      break;
    case '?': fprintf(stderr, "Unrecognized option: -%c\n", optopt); exit(EXIT_FAILURE);
    case ':': fprintf(stderr," Option -%c needs an argument\n", optopt); exit(EXIT_FAILURE);
    default: exit(EXIT_FAILURE);
    }
    count++;
  }


  /* Set co-occurrence default priors */
  if (search->cooc) {
    if (!pos_set) {
      search->posnone = POSNONE_C;
      search->posone  = POSONE_C;
      search->poscooc = POSCOOC_C;
    }
    
    if (!neg_set) {
      search->negnone = NEGNONE_C;
      search->negone  = NEGONE_C;
      search->negcooc = NEGCOOC_C;
    }
  }

  /* Make sure input files are present */
  if (argc - count !=  2) {
    fprintf(stderr, USAGE);
    exit(EXIT_FAILURE);
  }

  /* Init rand function */
  srand(seed);

  /* Switch: Positive set should be index 1, negative 0 */
  search->setfiles = (char **) malloc(sizeof(char *) * 2);
  search->setfiles[0] = (char *) malloc(sizeof(char *) * (strlen(argv[count+1]) + 1));
  search->setfiles[1] = (char *) malloc(sizeof(char *) * (strlen(argv[count]) + 1));
  strcpy(search->setfiles[0], argv[count+1]);
  strcpy(search->setfiles[1], argv[count]);


  return count;
}


/************* MAIN **************/
int main(int argc, char **argv) {
  Annealing *search = (Annealing *) malloc(sizeof(Annealing));
  PSSMSet pssmset;
  int ind;
  
  /******************** INIT ********************/
  ind = parse_options(argc, argv, search);
  create_logtable(TOTAL);


  /******************** SEQUENCES ********************/
  if (search->loglevel >= VERBOSE) fprintf(stderr, "Building suffix array...\n");
  search->sesa = build_SESA(2, search->setfiles, "ACGT", "N", 0);

  if (!search->sesa) {
    fprintf(stderr, "Could not build suffix array: %s\n", getError());
    exit(EXIT_FAILURE);
  }

  /******************** PSSMs ********************/
  if (wm_filename) {
    pssmset = load_matrices(wm_filename, MIN_LENGTH, MAX_LENGTH, ALPHLEN, TOTAL);

    if (pssmset->pssm_count == 1) {
      search->scorefunc  =  SINGLE;
    } else if (pssmset->pssm_count == 2){
      search->scorefunc  =  DOUBLE;
    } else {
      fprintf(stderr, "File contained %i Pssms. Number of Pssms must be either 1 or 2\n", wmcount);
      exit(EXIT_FAILURE);
    }

  } else {
    pssmset = random_pssmset(wmcount, 0, search->min_len, search->max_len, ALPHLEN, TOTAL); 
  }

  /******************** INITIALIZE ********************/
  init_annealing(search, pssmset->pssm_count, pssmset->pssms, TOTAL);

  /* Set the initial cutoff */
  ind = pssmset->pssm_count;
  while (ind--) {
    search->cut[ind] = cut;

    if (search->ctype == RATIO || search->ctype == ANNEALED) {
      search->cutoff[ind] = search->cut[ind] * (float) pssmset->pssms[2 * ind]->length;
    } else if (search->ctype == FIXED) {
      search->cutoff[ind] = search->cut[ind];
    }

    calcAndSetThresholds(pssmset->pssms[2 * ind], search->cutoff[ind]);
    calcAndSetThresholds(pssmset->pssms[2 * ind + 1], search->cutoff[ind]);

    if (search->loglevel >= VERBOSE) {
      printf("Setting threshold for %i to %f \n", ind, search->cutoff[ind]);
    }
  }

  /******************** LOG ********************/
  print_header(search);

  /******************** MAIN ********************/
  anneal(search);
  report(search); 


  /******************** CLEAN UP ********************/
  free_annealing(search); 
  releasePSSMSet(pssmset);

  return EXIT_SUCCESS;
}

