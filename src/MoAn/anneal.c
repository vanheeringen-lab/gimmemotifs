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
#include <string.h>
#include <unistd.h>
#include <float.h>
#include <math.h>
#include "pssm.h"
#include "search.h"
#include "moan.h"
#include "anneal.h"


/* DEFINITIONS */
#define LOG_STEP(pre, loglvl) if (search->loglevel >= loglvl) printf("%s %8i  Score: %8.2f   Temp: %6.3f\n", pre, search->maxit - i, score, temperature); 


/* Possible steps */
#define MOD           0
#define INC_L         1
#define INC_R         2
#define DEC_L         3
#define DEC_R         4
#define MOV_L         5
#define MOV_R         6
#define CUT           7

/* Not temperature dependent (LN) */
#define INCLEN  1 + RAND_INT((float)(pssm->max_length - pssm->length - 1))
#define DECLEN  1 + RAND_INT((float)(pssm->length - pssm->min_length - 1))
#define MOVLEN  1 + RAND_INT(4)
#define STEAL_AMOUNT  (1 + (int) RAND_INT(nwm->counts[WM_IND(nwm, pos,victim)] - 1))

/* Div */
#define WGT        0.001
#define MWGT       0.999
#define HIGH_TEMP FLT_MAX

/* Table of precalculated log values */
float *logtable;
float acc;


/*************** PROTOTYPES ***************/

/* Annealing */
float anneal_step_wm(Annealing *search, float temperature, float score, int wm_index);
float anneal_search(Annealing *search, int wm_index);
void backroll_wm(Annealing *search, int wm_index);
void save_wm(Annealing *search, int wm_index);

/* Anneal stepping */
void decrease_left(Annealing *search, const unsigned int wm_index, const unsigned int length); 
void decrease_right(Annealing *search, const unsigned int wm_index, const unsigned int length);
void increase_left(Annealing *search, const unsigned int wm_index, const unsigned int total, const unsigned int length);
void increase_right(Annealing *search, const unsigned int wm_index, const unsigned int total, const unsigned int length);
void move_left(Annealing *search, const unsigned int wm_index, const unsigned int total, const unsigned int length); 
void move_right(Annealing *search, const unsigned int wm_index, const unsigned int total, const unsigned int length); 
void step_cutoff(Annealing *search, int wmindex, float temperature);

/* Div */
void new_columns(Annealing *search, const int wm_index, const unsigned int first_col, unsigned int exp_length, const unsigned int total);
void init_temp(Annealing *search, int iter);
void init_starttemp(Annealing *search, int iter);
inline void normalize_pssm(PSSM pssm, unsigned int sum);
void create_logtable(int seq);

/* Random matrices */
PSSMSet random_pssms(unsigned int count, const unsigned char order,  const unsigned char min_length, const unsigned char max_length, const unsigned char alphlen, const int seq);
PSSM random_pssm(unsigned char order,  const unsigned char min_length, const unsigned char max_length, unsigned char alphlen, int seq);



/*************** FUNCTIONS ***************/

char *get_match(SESA sesa, int pos, int len) {
  char *match = (char *) malloc(sizeof(char) * (len+1));

  strncpy(match, &sesa->sequences[pos], len);
  match[len] = '\0';

  return match;
}


int compare_hits(const void *pEntry1, const void *pEntry2) {
  return ((struct HitEntry *)pEntry1)->position - ((struct HitEntry *)pEntry2)->position;
}


void print_pos_hits(SESA sesa, struct HitTable *hits, int wm_length) {
  struct HitEntry *end = hits->pScores + hits->nScores;
  struct HitEntry *hit = hits->pScores;
  int seq, last_seq, neg_seq = sesa->nSeqSet[0];

  qsort((void*)hits->pScores, hits->nScores, sizeof(struct HitEntry), compare_hits);

  last_seq = -1;

  for (hit = hits->pScores; hit < end; hit++) {
    seq = sesa->seq[hit->position];

    if (sesa->set[seq]) {
      printf("%3i:%5i  %s   %7.2f   %-10s\n", seq-neg_seq, sesa->pos[hit->position],  get_match(sesa, hit->position, wm_length), hit->score, sesa->seqnames[seq]);
    }
  } 
  
}

void print_hits(SESA sesa, struct HitTable *hits, int wm_length) {
  struct HitEntry *end = hits->pScores + hits->nScores;
  struct HitEntry *hit = hits->pScores;
  int seq, last_seq;

  qsort((void*)hits->pScores, hits->nScores, sizeof(struct HitEntry), compare_hits);

  last_seq = -1;

  for (hit = hits->pScores; hit < end; hit++) {
    seq = sesa->seq[hit->position];
    printf("%3i: %10s   %5i   %s   %7.2f\n", seq, sesa->seqnames[seq], sesa->pos[hit->position], get_match(sesa, hit->position, wm_length), hit->score);
  } 
  
}


/* 
 * Normalizes a column in a count matrix
 * pssm - The PSSM to normalize
 * pos  - The column number
 * sum  - The total sum for the column
 *
 * FIX: Only zeroth-order
 */
inline void normalize_col(PSSM pssm, unsigned int pos, unsigned int sum) {
  unsigned int i, ind, total = 0;
  unsigned int offset = pssm->offsets[pos];

  i = pssm->alphabetSize;
  while (i--) {
    total += pssm->counts[offset + i];
  }

  if (total != 0) {
    i = pssm->alphabetSize;
    while (i--) {
      ind = offset + i;
      pssm->counts[ind] = (int) (((float) (pssm->counts[ind] * sum)) / (float) total);
      pssm->scores[ind] = logtable[pssm->counts[ind]];
    }
  } else {
    i = pssm->alphabetSize;
    while (i--) {
      ind = offset + i;
      pssm->counts[ind] = sum / pssm->alphabetSize;
      pssm->scores[ind] = logtable[pssm->counts[ind]];
    }
  }

}


/* 
 * Normalizes a count matrix
 * pssm - The PSSM to normalize
 * sum  - The total sum for each column
 *
 * FIX: Only zeroth-order
 */
inline void normalize_pssm(PSSM pssm, unsigned int sum) {
  unsigned int pos;

  pos = pssm->length;
  while (pos--) normalize_col(pssm, pos, sum);
}

/********** ANNEALING **********/

/* FIX: DEPEND UPON COMPOSITION? */
void create_logtable(int seq) {
  logtable = (float *) malloc((seq + 1) * sizeof(float));
  int i;

  logtable[0] = NOOCCUR;

  for (i = 1; i <= seq; i++) 
    logtable[i] = (float) LOG((((float) i / seq) * 4.0));

}

float anneal(Annealing *search) {
  unsigned int wmcount = search->pssm_count;
  unsigned int i, k, wm;
  float temperature, max_score;;
  float score = 0.0;


  init_temp(search, STARTEMP_ITERATIONS);
  temperature = search->start_temp;

  acc = 0.5;

  PSSM *best = malloc(wmcount * sizeof(PSSM));

  i = wmcount;
  while (i--) {
    best[i] = clone_pssm(WM(i));   
    search->best_cut[i] = search->cut[i];
    score = anneal_search(search, i);
/*     struct HitTable *hits = HITS(i);   */
  }

  
  max_score = score; 
  search->best_step = 0;
  search->best_pos = search->pos_score;
  search->best_neg = search->neg_score;

  i = wmcount;
  while (i--) 
    save_wm(search, i);

  /* Annealing */
  if (search->loglevel >= NORMAL) {
    printf("\nAnnealing\n=============================================\n");
  }

  wm = 0;
  i = search->maxit;
  LOG_STEP("---", NORMAL);

  while (i--) {
    if (wmcount == 2) wm = !wm;
    score = anneal_step_wm(search, temperature, score, wm);

    /* New best score => save results */
    if (score > max_score) {
      max_score = score;
      search->best_step = search->maxit - i;
      search->best_pos = search->pos_score;
      search->best_neg = search->neg_score;
      
      k = wmcount;
      while (k--) {
	copy_pssm(WM(k), best[k]);
	search->best_cut[k] = search->cut[k];
      }

      LOG_STEP("***", NORMAL);
    }

    if (!(i % 10000)){
      LOG_STEP("---", NORMAL);
    } else {
      LOG_STEP("---", TRACE);
    }
    
    /* Dynamic ratio */
    if (search->schedule == EXPONENTIAL) {
      temperature *= EXP_FACTOR;
    } else {
      temperature -= LIN_FACTOR;
    }
  }

  k = wmcount;
  while (k--) {
    copy_pssm(best[k], WM(k));
    search->cut[k] = search->best_cut[k];
  }

  search->best_score = max_score;

  /* Clean up */
  i = wmcount;
  while (i--) {
    anneal_search(search, i);
    free(best[i]);
  }
  free(best);

  score = (*search->scorefunc)(search);

  return score;
}


float anneal_step_wm(Annealing *search, float temperature, float score, int wm_index) {
  float ra = RAND_FLOAT(1.0); 
  float newscore;
  PSSM pssm = WM(wm_index);

  search->last_step[wm_index] = -1;
    /* Modify the WM */

  do {
    ra = RAND_FLOAT(1.0);

    if (ra < PROB_INC_L) {
      if (pssm->length < pssm->max_length) {
	increase_left(search, wm_index, TOTAL, INCLEN);
	search->last_step[wm_index] = INC_L;
      }
    } else if (ra < PROB_INC_R) {
      if (pssm->length < pssm->max_length) {
	increase_right(search, wm_index, TOTAL, INCLEN);
	search->last_step[wm_index] = INC_R;
      }
    } else if (ra < PROB_DEC_L) {
      if (pssm->length > pssm->min_length) {
	decrease_left(search, wm_index, DECLEN);
	search->last_step[wm_index] = DEC_L;
      }
    } else if (ra < PROB_DEC_R) {
      if (pssm->length > pssm->min_length) {
	decrease_right(search, wm_index, DECLEN);
	search->last_step[wm_index] = DEC_R;
      }
    } else if (ra < PROB_MOV_L) {
      if (search->pos_matches > 0) {
	search->last_step[wm_index] = MOV_L;
	move_left(search, wm_index, TOTAL, MOVLEN);
      }
    } else if (ra < PROB_MOV_R) {
      if (search->pos_matches > 0) {
	search->last_step[wm_index] = MOV_R;
	move_right(search, wm_index, TOTAL, MOVLEN);
      }
    } else if (ra < PROB_CUT && search->ctype == ANNEALED) {
      step_cutoff(search, wm_index, temperature);
      search->last_step[wm_index] = CUT;
    } else {
      (*search->stepfunc)(search, wm_index, temperature);
      search->last_step[wm_index] = MOD;
    }
  } while (search->last_step[wm_index] == -1);

  pssm =  WM(wm_index);
  calcAndSetThresholds(pssm, search->cut[wm_index] * (float)pssm->length);  
  newscore = anneal_search(search, wm_index);

  double prob = exp(-(score - newscore)/temperature);

  if (RAND_FLOAT(1.0) < prob) {      /* Accept */   
    save_wm(search, wm_index);
    acc = WGT + MWGT * acc;
    return newscore;
  } else {              /* Reject */
    backroll_wm(search, wm_index);
    acc *= MWGT;
    return score;
  }

}



float anneal_search(Annealing *search, int wm_index) {
  double score;

  SWITCH_HITS(wm_index);

  /* Fetch the correct hittable */
  g_pHitTable = HITS(wm_index);
  g_pHitTable->nScores = 0;
  g_EndOfEntries = END(wm_index);
  g_curSize = SIZE(wm_index);
  g_NextHitEntry = g_pHitTable->pScores;

  SESA_search(search->sesa, WM(wm_index));
  search->fwd[wm_index] = g_pHitTable->nScores; /* Save FWD hit count */

  if (search->ds) {
    revcomp_copy_pssm(WM(wm_index), search->ds_tmp); 
    SESA_search(search->sesa, search->ds_tmp); 
  }

  /* May have been reallocated during search */
  HITS(wm_index) = g_pHitTable;
  END(wm_index) =   g_EndOfEntries;
  SIZE(wm_index) = g_curSize;

  score = (*search->scorefunc)(search);

  return score;
}


/** Stepping **/
void step_cutoff(Annealing *search, int wmindex, float temperature) {
  float ratio = temperature / search->start_temp;
  float interval = ratio * CUT_INT;
  float ra = RAND_FLOAT(1.0); 
  float res = interval * ra - (interval/2.0);
  search->cut[wmindex] += res;
  
  if (search->cut[wmindex] > CUT_MAX) {
    search->cut[wmindex] = CUT_MAX;
  } else if (search->cut[wmindex] < CUT_MIN) {
    search->cut[wmindex] = CUT_MIN;
  }
}



/* Initialize the start temperature and the step ratio */
void init_temp(Annealing *search, int iter){

  if (search->start_temp < 0) {
    if (search->loglevel >= VERBOSE) printf("Calculating start temperature \n");

    init_starttemp(search, iter);

    if (search->loglevel >= NORMAL) printf("Start temperature: %f \n", search->start_temp);
  }

  if (search->step_ratio < 0) {
    if (search->schedule == EXPONENTIAL) search->step_ratio = exp( log(search->end_temp/search->start_temp)/search->maxit);
    else                                 search->step_ratio =  (search->start_temp - search->end_temp) / (float)search->maxit;    
  }

  
}


/* Initializes the starting temperature */
void init_starttemp(Annealing *search, int iter) {
  float mean = 0.0, stddev = 0.0;
  float score[iter];
  unsigned int wmcount = search->pssm_count;
  int i, wm;
  PSSM wms[wmcount];

  /* Backup WMs*/
  i = wmcount;
  while (i--) {
    wms[i] = clone_pssm(WM(i));
  }

  search->start_temp = HIGH_TEMP;

  /* FIX search all */
  i = wmcount * 2;
  while (i--) {
    g_pHitTable = search->hits[i];
    g_EndOfEntries = search->hittable_last[i];
    g_curSize = search->hittable_size[i];
    reset_hits();

    SESA_search(search->sesa,  search->wm[i]);

    /* If they are reallocated */
    search->hits[i] = g_pHitTable;
    search->hittable_size[i] = g_curSize;
    search->hittable_last[i] = g_EndOfEntries;
  }

  if (search->loglevel >= VERBOSE) {
    printf("Starting temp search \n");
  }

  i = iter - 1;
  score[i] = (*search->scorefunc)(search);

  wm = 0;
  while (i--) {
    if (wmcount == 2) wm = !wm;
    score[i] = anneal_step_wm(search, HIGH_TEMP, score[i+1], wm);
    mean += fabs(score[i+1] - score[i]);
  }

  if (search->loglevel >= VERBOSE) {
    printf("Stopping temp search\n");
  }

  mean /= (float) iter;

  i = iter;
  while (i--) stddev += (score[i] - mean) * (score[i] - mean);


  stddev /= iter;

  if (search->loglevel >= DBG) {
    printf("MEAN  : %f\n", mean);
    printf("VAR   : %f\n", stddev);
  }

  stddev = sqrt(stddev);

  if (search->loglevel >= DBG) 
    printf("STDDEV: %f\n", stddev);

  search->start_temp = (mean + (stddev * 0.85)) * -log(0.8);

  /* Restore WMs */
  i = wmcount;
  while (i--) {
    copy_pssm(wms[i], WM(i));
    free(wms[i]);
  }

}


void save_wm(Annealing *search, int wm_index) {
  search->oldcut[wm_index] = search->cut[wm_index];
}

void backroll_wm(Annealing *search, int wm_index) {
  switch (search->last_step[wm_index]) {
  case MOD: SWITCH_WM(wm_index); break; 
  case CUT: search->cut[wm_index] = search->oldcut[wm_index]; break;
  case INC_R: WM(wm_index)->length -= WM(wm_index)->last_col; break;
  case INC_L: SWITCH_WM(wm_index); break;
  case DEC_R: WM(wm_index)->length += WM(wm_index)->last_col; break;
  case DEC_L: SWITCH_WM(wm_index); break;
  case MOV_R: SWITCH_WM(wm_index); break;
  case MOV_L: SWITCH_WM(wm_index); WM(wm_index)->length += WM(wm_index)->last_col; break; 
  }

  SWITCH_HITS(wm_index);
}


void init_annealing(Annealing *search, int wmcount, WM **wm, int total) {
  int i;

  search->wm = wm;
  search->pssm_count = wmcount;
  search->cut = (float *) calloc(wmcount+1, sizeof(float));
  search->oldcut = (float *) calloc(wmcount+1, sizeof(float));
  search->best_cut = (float *) calloc(wmcount+1, sizeof(float));
  search->cutoff = (float *) calloc(wmcount+1, sizeof(float));
  search->match_count = (unsigned int *) calloc(wmcount, sizeof(unsigned int));
  search->cur = (unsigned int *) calloc(wmcount, sizeof(unsigned int));
  search->cur_hits = (unsigned int *) calloc(wmcount, sizeof(unsigned int));
  search->last_step = (int *) malloc(wmcount * sizeof(int)); 
  search->shared_cutoff = FALSE;
  search->pos_cooc_matches = 0;
  search->neg_cooc_matches = 0;

  search->hittable_size = (int *) malloc(sizeof(int) * 2 * wmcount);
  search->hits = (struct HitTable **) malloc(sizeof(struct HitTable *) * 2 * wmcount);
  search->hittable_last = (struct HitEntry **) malloc(sizeof(struct HitEntry *) * 2 * wmcount);
  search->fwd = (unsigned int *) malloc(sizeof(unsigned int) * wmcount);
  
  if (search->ds) {
    search->ds_tmp = random_pssm(0, search->min_len, search->max_len, ALPHLEN, TOTAL);
  }

  i = 2 * wmcount;
  while (i--) {
    g_inited = 0;
    init_hittable();
    search->hits[i] = g_pHitTable;
    search->hittable_last[i] = g_EndOfEntries;
    search->hittable_size[i] = g_curSize;
  }

}

void free_annealing(Annealing *search) {
  int i;

  i = 2 * search->pssm_count;
  while (i--) {
    release_hits(search->hits[i]);
  }

  free_SESA(search->sesa);

  free(search->hits);
  free(search->cur);
  free(search->best_cut);
  free(search->cutoff);
  free(search->oldcut);
  free(search->cut);
  free(search->match_count);
  free(search->last_step);
  free(search->setfiles[1]);
  free(search->setfiles[0]);
  free(search->setfiles);

  if (search->ds) {
    releaseMatrix(search->ds_tmp);
  }

  free(search);
}



/************************** STEPPING **********************************/

/* 
 * Creates the contents of newly added columns
 * 
 * wm_index - The index of the pssm where columns are added
 * col - 
 * length -
 * total - The total that the count columns should sum to
 *
 */

void new_columns(Annealing *search, const int wm_index, const unsigned int first_col, unsigned int exp_length, const unsigned int total) {
  SESA sesa = search->sesa;
  PSSM pssm = WM(wm_index);
  unsigned int last = pssm->offsets[first_col + exp_length];
  unsigned int first = pssm->offsets[first_col];
  struct HitTable *hits = HITS(wm_index);
  unsigned int seq_first, seq_last, i, j;
	     

  i = last;
  while (i-- > first) pssm->counts[i] = 0;

  if (hits->nScores) {    

    i = hits->nScores;
    while (i--) {
      int pos = hits->pScores[i].position;
      int seq = sesa->seq[pos];

      if (sesa->set[seq]) {

	if (search->ds && i >= search->fwd[wm_index]) {
	  /* REVERSE STRAND */
	  seq_first = pos + WM(wm_index)->length - first_col - 1;
	  seq_last = seq_first - exp_length;
	
	  for (j = 0; j < exp_length; j++) {
	    int ltr = (int) search->sesa->esa->pStr[seq_first - j]; /* Get complement */
	    if (ltr == 254) break; /* End of sequence, jump to next hit */
	    ltr = 3 - ltr; /* Complement */

	    pssm->counts[ pssm->offsets[ first_col + j ] + ltr]++;
	  }
	} else {
	  /* FORWARD STRAND */
	  seq_first = pos + first_col;
	  seq_last = seq_first + exp_length;
	
	  for (j = 0; j < exp_length; j++) {
	    int ltr = (int) search->sesa->esa->pStr[seq_first + j];
	    if (ltr == 254) break; /* End of sequence, jump to next hit */

	    pssm->counts[ pssm->offsets[ first_col + j ] + ltr]++;
	  }
	}
      }
    }


    i = first_col + exp_length; 
    while (i--  > first_col) normalize_col(pssm, i, total);
    
  }  else {
    i = last;
    while (i-- > first) {
      /* FIX: Not higher order friendly */
      pssm->counts[i] = TOTAL/pssm->alphabetSize;
      pssm->scores[i] = logtable[TOTAL/pssm->alphabetSize];
    }
  }
}


/* Matrix length alteration and displacement */

void decrease_right(Annealing *search, const unsigned int wm_index, const unsigned int length) {
  PSSM wm = WM(wm_index);

  wm->last_col = length;
  wm->length -= wm->last_col;

}

void decrease_left(Annealing *search, const unsigned int wm_index, const unsigned int length) {
  SWITCH_WM(wm_index);
  PSSM wm = WM(wm_index);
  PSSM nwm = NWM(wm_index);
  int pos, old_pos;
  int offset = length * ALPHLEN;    

  wm->last_col = length; 
  wm->length = nwm->length - wm->last_col; 

  pos = wm->length * ALPHLEN;
  while (pos--) {
    old_pos = pos + offset;
    wm->counts[pos] = nwm->counts[old_pos];
    wm->scores[pos] = nwm->scores[old_pos];
  }

}


void increase_right(Annealing *search, const unsigned int wm_index, const unsigned int total, const unsigned int length) {
  PSSM wm = WM(wm_index);

  new_columns(search, wm_index, wm->length, length, total);
  wm->last_col = length;
  wm->length += wm->last_col; 

}

void increase_left(Annealing *search, const unsigned int wm_index, const unsigned int total, const unsigned int length) {

  SWITCH_WM(wm_index);
  PSSM wm = WM(wm_index);
  PSSM nwm = NWM(wm_index);
  int cell, old_cell;

  wm->last_col = length;
  wm->length = nwm->length + length;

  old_cell = nwm->offsets[nwm->length];
  cell = wm->offsets[nwm->length + length];

  while (old_cell--) {
    cell--;
    wm->counts[cell] = nwm->counts[old_cell];
    wm->scores[cell] = nwm->scores[old_cell];
  }

  new_columns(search, wm_index, 0, wm->last_col, total);

}

void move_left(Annealing *search, const unsigned int wm_index, const unsigned int total, const unsigned int length) {
  decrease_right(search, wm_index, length); 
  increase_left(search, wm_index, total, length); 


}

void move_right(Annealing *search, const unsigned int wm_index, const unsigned int total, const unsigned int length) {
  decrease_left(search, wm_index, length);
  increase_right(search, wm_index, total, length);

}


/*
 * Modifies a WM by stealing from one letter and given to another at each position
 */
void wm_steal(Annealing *search, int wm_index, float temperature) {
  unsigned int pos, thief, victim;
  unsigned int amount, i, nind, wind;
  PSSM wm, nwm;

  wm = WM(wm_index);
  nwm = NWM(wm_index);
  SWITCH_WM(wm_index);
  wm = WM(wm_index);
  nwm = NWM(wm_index);
  wm->length = nwm->length;

  
  pos = wm->length;
  while (pos--) {    
    
    /* Steal from one letter, and add to another */
    victim = RAND_INT(nwm->alphabetSize);
    while ( nwm->counts[ WM_IND(nwm, pos, victim)  ] == 0) victim = RAND_INT(4.0);

    thief = RAND_INT(nwm->alphabetSize - 1);
    if (thief >= victim) thief++; 
    amount = STEAL_AMOUNT;

    i = ALPHLEN;
    while (i--) {      
      nind = WM_IND(nwm, pos, i);
      wind = WM_IND(wm, pos, i);
      wm->counts[wind] = nwm->counts[nind];

      if (i == thief) 
	wm->counts[wind]  += amount;
      else if(i == victim)
	wm->counts[wind] -= amount;
      wm->scores[wind] = logtable[wm->counts[wind]];
    }
  }
}

/**
 * Loads a count matrix from file and converts it to log-odds scores
 */
PSSMSet load_matrices(char *filename, unsigned int min_length, unsigned int max_length, unsigned int alphlen, unsigned int total) {
  PSSMSet pssmset = (PSSMSet) malloc(sizeof(struct PSSMSet));
  unsigned int i, j, k, order, length;
  FILE *file;
  PSSM pssm;

  file = fopen(filename, "r");
  fscanf(file, "# %d\n", &i);
  pssmset->pssms = (PSSM *) malloc(2 * i * sizeof(PSSM));
  pssmset->pssm_count = i;

  fprintf(stderr, "PSSMs: %i\n", i);

  /* We need two copies of every matrix to speed up annealing */
  for (i = 0; i < 2 * pssmset->pssm_count; i+=2) {
    fscanf(file, "> %d %d\n", &order, &length);        
    pssmset->pssms[i] = initMatrix(order, max_length, alphlen);
    pssmset->pssms[i]->length = length;
    
    pssm = pssmset->pssms[i];
    for (j = 0; j < alphlen; j++) 
      for (k = 0; k < length; k++) 
	fscanf(file, "%d", &pssm->counts[alphlen * k + j]);


    fscanf(file, "\n");
    
    print_counts(pssm);
    normalize_pssm(pssm, total);
    print_counts(pssm);
    
  }
  fclose(file);
    
  for (i = 0; i < 2 * pssmset->pssm_count; i += 2) {
    calcAndSetThresholds(pssmset->pssms[i], 0.0); /* To avoid unitialized values */
    pssmset->pssms[i]->max_length = max_length;
    pssmset->pssms[i]->min_length = min_length;
    pssmset->pssms[i+1] = clone_pssm(pssmset->pssms[i]);  
    pssmset->pssms[i+1]->max_length = max_length;
    pssmset->pssms[i+1]->min_length = min_length;
  }

  return pssmset;
}


PSSMSet random_pssmset(unsigned int count, const unsigned char order,  const unsigned char min_length, const unsigned char max_length, const unsigned char alphlen, const int seq) {
  PSSMSet pssmset = (PSSMSet) malloc(sizeof(struct PSSMSet));
  pssmset->pssms = (PSSM *) malloc(2 * count * sizeof(PSSM));
  pssmset->pssm_count = count;

  while (count--) {
    pssmset->pssms[2 * count] = random_pssm(order, min_length, max_length, alphlen, seq);
    pssmset->pssms[2 * count + 1] = clone_pssm(pssmset->pssms[2 * count]);

    pssmset->pssms[2 * count]->min_length = min_length;
    pssmset->pssms[2 * count]->max_length = max_length;
    pssmset->pssms[2 * count + 1]->min_length = min_length;
    pssmset->pssms[2 * count + 1]->max_length = max_length;
  }

  return pssmset;
}

PSSM random_pssm(unsigned char order,  const unsigned char min_length, const unsigned char max_length, unsigned char alphlen, int seq) {
  PSSM pssm = initMatrix(order, max_length, alphlen);
  int i;

  pssm->max_length = max_length;

  i = pssm->max_length * pssm->alphabetSize;
  while (i--) pssm->counts[i] = 0 + RAND_INT(seq);

  normalize_pssm(pssm, seq);

  pssm->min_length = min_length;
  pssm->length = min_length + RAND_INT(max_length - min_length);
  
  calcAndSetThresholds(pssm, RAND_FLOAT((2.0 * pssm->length)));
  
  return pssm;
}



