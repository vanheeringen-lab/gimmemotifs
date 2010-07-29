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

#include <math.h>
#include <limits.h>


#include "score.h"
#include "pssm.h"
#include "search.h"
#include "moan.h"

#define MAX_HIT   100   /* Maximum number of hits in one sequence */

/* Functions */
inline bool overlap(unsigned int s1, unsigned int l1, unsigned int s2, unsigned int l2);

inline bool overlap(unsigned int s1, unsigned int l1, unsigned int s2, unsigned int l2) {
  if ( (s1 >= s2 && s1 <= s2+l2) || (s2 >= s1 && s2 <= s1+l1) ) {
    return TRUE;
  } else {
    return FALSE;
  }
}


float score_cooc_likelihood(Annealing *search) {
  SESA sesa = search->sesa;
  unsigned int count[search->pssm_count][sesa->nSeq];
  float scores[search->pssm_count][sesa->nSeq][MAX_HIT];
  int positions[search->pssm_count][sesa->nSeq][MAX_HIT];
  int pos_seqs = sesa->nSeqSet[1];
  int neg_seqs = sesa->nSeqSet[0];
  float ps[sesa->nSeq], ns[sesa->nSeq];
  float pc[sesa->nSeq], nc[sesa->nSeq];
  int wm, nwm, seq, pos, i, c;
  int wm_len, nwm_len;
  float score;

  search->pos_score = search->posnone;
  search->neg_score = search->negnone;
  search->pos_matches = 0;
  search->neg_matches = 0;
  search->pos_cooc_matches = 0;
  search->neg_cooc_matches = 0;

  /* Init count matrix */
  wm = search->pssm_count;
  while (wm--) {
    seq = sesa->nSeq;
    while (seq--) count[wm][seq] = 0;
  }
 
  seq = sesa->nSeq;
  while (seq--) {
    ps[seq] = 0.0;
    ns[seq] = 0.0;
    pc[seq] = 0.0;
    nc[seq] = 0.0;
  }

  wm = search->pssm_count;
  while (wm--) {
    struct HitTable *hits = HITS(wm);
    
    nwm = !wm;
    wm_len = WM(wm)->length;
    nwm_len = WM(nwm)->length;

    i = hits->nScores;
    while (i--) {
      struct HitEntry *hit = &hits->pScores[i];
      seq = sesa->seq[hit->position];
      pos = sesa->pos[hit->position];
      score = hit->score;

      /* 
       * More hits than MAX_HITS are ignored. The optimal solution
       * doesn't have hits all over the place.
       */
      c = count[wm][seq];
      if (count[wm][seq] < MAX_HIT) { 
	scores[wm][seq][c] = score;
	positions[wm][seq][c] = pos;
	count[wm][seq]++;
	
	if (sesa->set[seq]) {  /* POSITIVE SEQUENCES */
	  ns[seq] += EXP(score);
	  ps[seq] += EXP(score);
	  search->pos_matches++;

	  score *= search->poscooc;
	  c = count[nwm][seq];
	  while (c--) {
	    /* Check for overlap */
	    if (!overlap(pos, wm_len, positions[nwm][seq][c], nwm_len)) {
	      nc[seq] += EXP(scores[nwm][seq][c]) * EXP(score);
	      pc[seq] += EXP(scores[nwm][seq][c]) * EXP(score);
	      search->pos_cooc_matches++;
	    }
	  }	
	} else { /* NEGATIVE SEQUENCES */
	  ns[seq] += EXP(score);
	  ps[seq] += EXP(score);
	  search->neg_matches++;

	  score *= search->negcooc;
	  c = count[nwm][seq];
	  while (c--) {
	    /* Check for overlap */
	    if (!overlap(pos, wm_len, positions[nwm][seq][c], nwm_len)) {
	      nc[seq] += EXP(scores[nwm][seq][c]) * EXP(score);
	      pc[seq] += EXP(scores[nwm][seq][c]) * EXP(score);
	      search->neg_cooc_matches++;
	    }
	  }		  
	}
      }
    }
  }


  i = sesa->nSeq;
  while (i-- > neg_seqs) {
    if (!((count[0][i]) || count[1][i])) {
      search->pos_score += LOG((float)pos_seqs / (float)(pos_seqs + neg_seqs));
    } else {
      ps[i] *= pos_seqs * search->posone;
      ns[i] *= neg_seqs * search->negone;
      nc[i] *= neg_seqs * search->negcooc;
      pc[i] *= pos_seqs * search->poscooc;
      float p = ps[i] + pc[i];
      float n = ns[i] + nc[i];
      search->pos_score += LOG(p / (n + p));
    }
  }

  while (i--) {
    if (!(count[0][i] || count[1][i])) {
      search->neg_score += LOG((float)neg_seqs / (float)(pos_seqs + neg_seqs));
    } else {
      ns[i] *= neg_seqs * search->negone;
      ps[i] *= pos_seqs * search->posone;
      nc[i] *= neg_seqs * search->negcooc;
      pc[i] *= pos_seqs * search->poscooc;
      float p = ps[i] + pc[i];
      float n = ns[i] + nc[i];
      search->neg_score += LOG(n / (n + p));
    }
  }

  return search->pos_score + search->neg_score;
}


/* The log-likelihood cost for a single WM */
float score_single_likelihood(Annealing *search) {
  int pos_seq = search->sesa->nSeqSet[1];
  int neg_seq = search->sesa->nSeqSet[0];
  float scores[search->pssm_count][search->sesa->nSeq];
  unsigned int *seq = search->sesa->seq;
  int i, wm, s;

  float prior_nn = search->negnone;
  float prior_ns = search->negone;
  float prior_pn = search->posnone;
  float prior_ps = search->posone;

  search->pos_score = 0.0; 
  search->neg_score = 0.0; 
  search->pos_matches = 0;
  search->neg_matches = 0;


  wm = search->pssm_count;
  while (wm--) {
    s = search->sesa->nSeq;
    while(s--) scores[wm][s] = 0.0;
  }

  wm = search->pssm_count;
  while (wm--) {
    struct HitTable *hits = HITS(wm);

    i = hits->nScores;
    while (i--){         /* Iterate the leaves */
      struct HitEntry *hit = &hits->pScores[i];
      scores[wm][ seq[hit->position] ] += EXP(hit->score);


      if (search->sesa->set[seq[hit->position]]) 
	search->pos_matches++;
      else
	search->neg_matches++;
    } 
    
    /* Sum up the scores */  
    s = search->sesa->nSeq;
    while (s-- > neg_seq) {
      float p;
       if (scores[wm][s] > 1.0) { 
	 p = ((float)pos_seq) * scores[wm][s];
       } else { 
	 p = (float) pos_seq; 
       } 

       p = ((p * prior_ps) + (pos_seq * prior_pn));

       search->pos_score += LOG( p / (((float)neg_seq) * prior_nn + ((float)neg_seq) * prior_ns * scores[wm][s] + p));
    }


    while (s--) {
      float n = (((float)neg_seq) * prior_nn + ((float)neg_seq) * prior_ns * scores[wm][s]);

      if (scores[wm][s] > 1.0) { 
	search->neg_score += LOG(n / (n + ((float)pos_seq) * scores[wm][s] * prior_ps +  ((float)pos_seq) * prior_pn));
      } else { 
 	search->neg_score += LOG(((float)neg_seq) / (((float)neg_seq) + ((float)pos_seq))); 
      } 
    }
    
  } 
  
  return search->pos_score + search->neg_score;
}


float score_likelihood_lennorm(Annealing *search) {
  SESA sesa = search->sesa;
  unsigned int count[search->pssm_count][sesa->nSeq];
  float scores[search->pssm_count][sesa->nSeq][MAX_HIT];
  int positions[search->pssm_count][sesa->nSeq][MAX_HIT];
  int pos_seqs = sesa->nSeqSet[1];
  int neg_seqs = sesa->nSeqSet[0];
  float ps[sesa->nSeq], ns[sesa->nSeq];
  float pc[sesa->nSeq], nc[sesa->nSeq];
  int wm, nwm, seq, pos, i, c;
  int wm_len, nwm_len;
  float score, norm;

  search->pos_score = search->posnone;
  search->neg_score = search->negnone;
  search->pos_matches = 0;
  search->neg_matches = 0;
  search->pos_cooc_matches = 0;
  search->neg_cooc_matches = 0;

  /* Init count matrix */
  wm = search->pssm_count;
  while (wm--) {
    seq = sesa->nSeq;
    while (seq--) count[wm][seq] = 0;
  }
 
  seq = sesa->nSeq;
  while (seq--) {
    ps[seq] = 0.0;
    ns[seq] = 0.0;
    pc[seq] = 0.0;
    nc[seq] = 0.0;
  }

  wm = search->pssm_count;
  while (wm--) {
    struct HitTable *hits = HITS(wm);
    
    nwm = !wm;
    wm_len = WM(wm)->length;
    nwm_len = WM(nwm)->length;
    norm = 1.0/wm_len;

    i = hits->nScores;
    while (i--) {
      struct HitEntry *hit = &hits->pScores[i];
      seq = sesa->seq[hit->position];
      pos = sesa->pos[hit->position];
      score = hit->score;

      /* 
       * More hits than MAX_HITS are ignored. The optimal solution
       * doesn't have hits all over the place.
       */
      c = count[wm][seq];
      if (count[wm][seq] < MAX_HIT) { 
	scores[wm][seq][c] = score *  norm;
	positions[wm][seq][c] = pos;
	count[wm][seq]++;
	
	if (sesa->set[seq]) {  /* POSITIVE SEQUENCES */
	  ns[seq] += EXP(score) * norm;
	  ps[seq] += EXP(score) * norm;
	  search->pos_matches++;

	  score *= search->poscooc;
	  c = count[nwm][seq];
	  while (c--) {
	    /* Check for overlap */
	    if (!overlap(pos, wm_len, positions[nwm][seq][c], nwm_len)) {
	      nc[seq] += EXP(scores[nwm][seq][c]) * EXP(score) *  norm;
	      pc[seq] += EXP(scores[nwm][seq][c]) * EXP(score) *  norm;
	      search->pos_cooc_matches++;
	    }
	  }	
	} else { /* NEGATIVE SEQUENCES */
	  ns[seq] += EXP(score) * norm;
	  ps[seq] += EXP(score) * norm;
	  search->neg_matches++;

	  score *= search->negcooc;
	  c = count[nwm][seq];
	  while (c--) {
	    /* Check for overlap */
	    if (!overlap(pos, wm_len, positions[nwm][seq][c], nwm_len)) {
	      nc[seq] += EXP(scores[nwm][seq][c]) * EXP(score) * norm;
	      pc[seq] += EXP(scores[nwm][seq][c]) * EXP(score) * norm;
	      search->neg_cooc_matches++;
	    }
	  }		  
	}
      }
    }
  }


  i = sesa->nSeq;
  while (i-- > neg_seqs) {
    if (!((count[0][i]) || count[1][i])) {
      search->pos_score += LOG((float)pos_seqs / (float)(pos_seqs + neg_seqs));
    } else {
      ps[i] *= pos_seqs * search->posone;
      ns[i] *= neg_seqs * search->negone;
      nc[i] *= neg_seqs * search->negcooc;
      pc[i] *= pos_seqs * search->poscooc;
      float p = ps[i] + pc[i];
      float n = ns[i] + nc[i];
      search->pos_score += LOG(p / (n + p));
    }
  }

  while (i--) {
    if (!(count[0][i] || count[1][i])) {
      search->neg_score += LOG((float)neg_seqs / (float)(pos_seqs + neg_seqs));
    } else {
      ns[i] *= neg_seqs * search->negone;
      ps[i] *= pos_seqs * search->posone;
      nc[i] *= neg_seqs * search->negcooc;
      pc[i] *= pos_seqs * search->poscooc;
      float p = ps[i] + pc[i];
      float n = ns[i] + nc[i];
      search->neg_score += LOG(n / (n + p));
    }
  }  

  return search->pos_score + search->neg_score;
}


float score_single_aic(Annealing *search) {
  float lhood = score_single_likelihood(search);
  int params = 1 + 3 * WM(0)->length; /* 3 free parameters <=> DNA alphabet - 1 */

  return 2.0 * (lhood - ((float)params));
}

float score_aic(Annealing *search) {
  float lhood = score_cooc_likelihood(search);
  int params = 2 + 3 * (WM(0)->length + WM(1)->length); /* 3 free parameters <=> DNA alphabet - 1 */

  return 2.0 * (((float)params) - lhood);
}


/*
 * Log of the posterior distribution P(W|x,y) for a single matrix
 */
float score_single_bayesian(Annealing *search) {
  return score_single_likelihood(search) + LOG(search->wmpriors[WM(0)->length]);
}


/*
 * Log of the posterior distribution P(W|x,y) 
 */
float score_bayesian(Annealing *search) {
  return 1.0;
}
