/* SESA is a wrapper for the ESA implementation */
/* Copyright (C) 2006 Bioinformatics Centre */
/* Author: Eivind Valen */

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
#include <ctype.h>
#include <math.h>

#include "sesa.h"
#include "hittable.h"
#include "moan.h"
#include "build.h"
#include "shared.h"
#include "error.h"


#define TABLE_START    100
#define TABLE_INCREASE 2
#define SEP_CHAR '|'
#define notN(c) (toupper(c) == 'N' ? FALSE : TRUE)
#define isACGT(c) ((toupper(c) == 'A' || toupper(c) == 'C' || toupper(c) == 'G' || toupper(c) == 'T') ? TRUE : FALSE)

/* Globals */
unsigned int pssm_list_size;

/* Prototypes */
int reinit_hittable();
int compare_integers(const int *int1, const int *int2);
inline void set_mismatch_scores(PSSM pssm, unsigned int pssmNr);
int *get_pssmnr_list(unsigned char *seq, int seq_length, int motif_length);
void SESA_recur_exhaust_search(SESA sesa, PSSM pssm, int pos, unsigned int pssmNr, double threshold, PSSMScoreTable scoreTable, int (*acceptanceFunc)(SESA, struct HitTable *));
inline void SESA_pfxsc(double C[], int *d, int l, PSSM pssm, unsigned char *base, int offset);
void SESA_add_sequence(SESA sesa, char *seq, unsigned int set);
int create_reverse_strand(SESA sesa, char *sequences, int size, int start_size, int start_seq);
PSSMScoreTable initScoreTable(int size); 
void reverse(char *sequences, int start, int end);

inline char complement(char i) {
  switch (i) {
  case 'A': case 'a': return 'T';
  case 'C': case 'c': return 'G';
  case 'G': case 'g': return 'C';
  case 'T': case 't': return 'A';
  case 'N': case 'n': return 'N';
  }

  return SEP_CHAR;
}

inline int chr2ind(char i) {
  switch (i) {
  case 'A': case 'a': return A;
  case 'C': case 'c': return C;
  case 'G': case 'g': return G;
  case 'T': case 't': return T;
  }

  return -1;
}

inline char ind2chr(int i) {
  switch (i) {
  case 0: return 'A';
  case 1: return 'C';
  case 2: return 'G';
  case 3: return 'T';
  }

  return '-';
}

inline void set_mismatch_scores(PSSM pssm, unsigned int pssmNr) {
  unsigned int pos, ltr, cor, off;

  pos = (unsigned int) pssm->length;
  while (pos--) {
    cor = (pssmNr >> (unsigned int)(2 * pos)) & 3;
    ltr = (unsigned int) pssm->alphabetSize;
    off = (unsigned int) pssm->offsets[pos];
    while (ltr--) {
      pssm->scores[off + ltr] = (ltr == cor ? 0.0 : -1.0); 
    }
  }
}

SESA init_SESA() {
  int i;

  SESA sesa = (SESA) malloc(sizeof(struct SESAStruct));
  sesa->maxSeq = 1000; /* DEFAULT_MAX_SEQ; */
  sesa->seqnames = (char **) malloc(sizeof(char *) * sesa->maxSeq);
  sesa->set = (unsigned int *) malloc(sizeof(unsigned int) * sesa->maxSeq);
  sesa->nSeq = 0;

  i = 4;
  while (i--) sesa->bg[i] = 0;

  return sesa;
}



/* SEARCHING */

struct HitTable *SESA_search(SESA sesa, PSSM pssm) { /* , int nHits) { *//* , struct HitTable *pHitTable) { */
  // Declaration of variables
  ESA esa = sesa->esa;
  int i, j, total, dCur, dPrev, numSufs, special;
  int curEsaIndex, pssmLength, curSuf;
  int lcpi;
  struct HitTable *pHitTable;
  double Ci[MAXPSSMSIZE];
  unsigned char *pStr;

  g_inited = 1;
  total = 0;

  // Getting the string and stringsize from the ESA
  pStr = getStr(esa);
  numSufs = getSize(esa);

  // Initialisation of d, C and the scoring system
  // dCur and dPrev keeps track of the number of letters scored
  // before the threshold is/was met in the current and in
  // the previous suffix (= matrixlength, suffix is a hit)

  /*   init_hittable();  */
  pssmLength = getLengthFast(pssm);

  // First d and C are calculated for the 1st suffix in the array:
  curEsaIndex = 0;
    
  curSuf = getSuf(esa,curEsaIndex);
  SESA_pfxsc(Ci, &dCur, pssmLength-1, pssm, &pStr[curSuf], curEsaIndex);
    
  // if the suffix is a hit it is registered in the hittable
  if(dCur == pssmLength-1){
    register_hit(curSuf, Ci[pssmLength-1]);
    total++;
/*     fprintf(stderr, "--- %i\n", curSuf); */
  }
        

  // The d and C are calc for the rest of the suffixes: 
  curEsaIndex++;    
  dPrev = dCur;

  while(curEsaIndex < numSufs){
    // For each suffix lcp between suf i-1 and suf i in the 
    // array is checked
    lcpi = getLcp(esa,curEsaIndex);
    
    // NOTE: It a special case, if d[curEsaIndex-1]+1 == lcpi ==
    // pssmLength and then we can also skip.
    special = dPrev+1 == lcpi && lcpi == pssmLength;
    
    // If the lcp =< di-1+1 (the number of letters 
    // succesfully scored for suf i-1) then the local scores
    // are calculated from the point of lcp (as the local 
    // scores must be the same up until the lcp)
    if( dPrev+1 >= lcpi && !special ){
       curSuf = getSuf(esa,curEsaIndex);	
       SESA_pfxsc(Ci, &dCur, pssmLength-1, pssm, &pStr[curSuf], lcpi);
      
       if(dCur == pssmLength-1){
	 register_hit(curSuf, Ci[pssmLength-1]);
	 total++;
       } 
       curEsaIndex++;
       dPrev = dCur;
    }
    // If they the lcp is larger than (di-1)+1 then di, C and score
    // must be the samme as those for the previous suffix and
    // this must hold for all following suffixes that have lcp > 
    // with the previous suf 
    else {
      j = getSkip(esa,curEsaIndex);
      
      // NOTE: In the special case we can skip any letter at position d[cEI-1]+1 !
      while( j < numSufs && 
	     (getLcp(esa,j)>dPrev+1 || (special && getLcp(esa,j)>dPrev))) {
	j = getSkip(esa,j);
      }
      
      // If all those we skip pass was hits: register them
      if(dPrev == pssmLength-1) {
	for(i = curEsaIndex; i < j; i++) {
	  register_hit(getSuf(esa,i), Ci[pssmLength-1]);
	  total++;
	}
      }
      curEsaIndex = j;
    }
  }  

  pHitTable = finish();

  return pHitTable;
}	


PSSMScoreTable initScoreTable(int size) {
  PSSMScoreTable table = (PSSMScoreTable) malloc(sizeof(struct PSSMScoreTable));

  table->pssmNr = (int *) malloc(sizeof(int) * size);
  table->size = size;
  table->nPssm = 0;

  return table;
}


inline void add_pssm(PSSMScoreTable table, int pssmNr) {
  table->pssmNr[table->nPssm++] = pssmNr;

  if (table->size == table->nPssm) {
    table->size *= TABLE_INCREASE;
    table->pssmNr = (int *) realloc(table->pssmNr, sizeof(int) * table->size);
  }
}



PSSMScoreTable SESA_exhaustive_search(SESA sesa, int length, int mismatches, int alphabetSize, int (*acceptanceFunc)(SESA, struct HitTable *)) {
  PSSMScoreTable scoreTable = initScoreTable(TABLE_START);
  PSSM pssm = initMatrix(0, length, alphabetSize);
  int i;

  init_hittable();
  i = pssm->offsets[pssm->length];
  while (i--) pssm->scores[i] = -1;
  SESA_recur_exhaust_search(sesa, pssm, 0, 0, -((double)mismatches)-0.5, scoreTable, acceptanceFunc);
  releaseMatrix(pssm);
  release_hits(finish());

  return scoreTable;
}


/* 
 * Recursively descend to last position generating all possible
 * matrices along the way.
 */
void SESA_recur_exhaust_search(SESA sesa, PSSM pssm, int pos, unsigned int pssmNr, double threshold, PSSMScoreTable scoreTable, int (*acceptanceFunc)(SESA, struct HitTable *)) {
  struct HitTable *hits;
  unsigned int i, new;

  if (pos == pssm->length) {    
    set_mismatch_scores(pssm, pssmNr);
    calcAndSetThresholds(pssm, threshold);
    hits = SESA_search(sesa, pssm);
    if (acceptanceFunc(sesa, hits)) {
      add_pssm(scoreTable, pssmNr);
    }
    reset_hits();
  } else {
    for (i = 0; i < pssm->alphabetSize; i++) {
      new = i << (pos * 2);
      SESA_recur_exhaust_search(sesa, pssm, pos + 1, pssmNr | new, threshold, scoreTable, acceptanceFunc);
    }
  }

}


void free_SESA(SESA sesa) {
  free_ESA(sesa->esa);
  free(sesa->seq);
  free(sesa->pos);
  free(sesa->set);
  free(sesa->nSeqSet);

  while (sesa->nSeq--) {
    free(sesa->seqnames[sesa->nSeq]);
  }
  free(sesa->seqnames);
}

PSSMScoreTable SESA_exclusive_search(SESA sesa, int length, int mismatches, int alphabetSize, int (*acceptanceFunc)(SESA, struct HitTable *)) {
  PSSMScoreTable scoreTable = initScoreTable(TABLE_START);
  PSSM pssm = initMatrix(0, length, alphabetSize);
  int last = 0, size = sesa->esa->size;
  struct HitTable *hits = NULL;
  int *list;
  int i, s = 0;


  init_hittable();
  list = get_pssmnr_list(sesa->esa->pStr, size, length);

  /* Sort, so we can skip non-unique */
  qsort(list, pssm_list_size, sizeof(int), (void *) compare_integers);

  i = pssm_list_size;
  while (i--) {
    if (list[i] == last) continue;

    set_mismatch_scores(pssm, list[i]);
    calcAndSetThresholds(pssm, -((double)mismatches)-0.5);    

    hits = SESA_search(sesa, pssm);
    if (acceptanceFunc(sesa, hits)) {
      add_pssm(scoreTable, list[i]);
    }

    s++;
    reset_hits();
    last = list[i];
  }

  release_hits(hits);
  releaseMatrix(pssm);
  free(list);

  return scoreTable; 
}


int *get_pssmnr_list(unsigned char *seq, int seq_length, int motif_length) {
  int *results = (int *) malloc(sizeof(int) * seq_length);
  int last = seq_length - motif_length;
  unsigned int pssmNr = 0;
  unsigned int valid_fields = pow(2, 2*motif_length) - 1;
  int start, pos, i = 0;

  for (pos = 0; pos < motif_length; pos++)
    pssmNr = (pssmNr << 2) | seq[pos];


  do {
    results[i++] = pssmNr;
    
    /* Skip ignore symbols */
    if (seq[pos] > 250) {   
      while (seq[++pos] > 250); 

      start = pos + motif_length - 1;
      pssmNr = seq[pos];

      while (pos++ < start) pssmNr = (pssmNr << 2) | seq[pos];
    }


    pssmNr = ((pssmNr << 2) | seq[pos]) & valid_fields;
  } while(pos++ < last);

  pssm_list_size = i;

  return results;
}




/* HELPER FUNCTIONS */

// Function for calculation Ci[d] (the prefix score) and di.
//
// The function assumes that the C[offset-1] is correct and that pssm
// is a single PSSM.
// C_i is updated from offset and up to l, as long a the score is
// above the threshold and the string separator isn't met.
inline void SESA_pfxsc(double C[], int *d, int l, PSSM pssm, unsigned char *base, int offset) {
  int j;

  // If the offset is greater that l, we set d_i to l
  *d = l;

  // Iterate over the letters in the string, starting at base+offset
  for(j = offset; j <= l; j++) {
    // Stop if a string seperator is met the suffix is too short.
    if(base[j] == BREAKSYM) {
      *d = j - 1;
      break;
    }

    if(j == 0)
      C[j] =  getScoreFast(pssm, &base[j], j);
    else 
      C[j] = C[j-1] + getScoreFast(pssm, &base[j], j);
    

    // Stop if the local score is below the local threshold
    if(C[j] < getThresholdFast(pssm, j) ) {
      *d = j - 1;
      break;
    }
    *d = j;
  }
  return;
}



int compare_integers(const int *int1, const int *int2) {
  return (int) (*int1 - *int2);
}

		
int reinit_hittable() {
  g_pHitTable = malloc(sizeof(struct HitTable));	
  if(!g_pHitTable) {
    setError("Couldn't allocate memory for hit table.");
    return 0;
  }
  g_curSize = START_SCORES;
  g_pHitTable->pScores = malloc ( sizeof(struct HitEntry) * g_curSize);	
  if(!g_pHitTable->pScores) {
    setError("Couldn't allocate memory for hit table.");
    free(g_pHitTable);
    return 0;		
  }

  g_NextHitEntry = g_pHitTable->pScores;
  g_EndOfEntries = g_pHitTable->pScores + g_curSize;
  g_overload = 0;
  g_inited = 1;	
  return 1;
}


/**
 * Helper function for build_SESA and build_dsSESA.
 * Adds a sequence to the SESA.
 *
 * sesa - Pointer to the SESA.
 * seq - Pointer to the name of the sequence.
 * size - Length of the string pointed to by seq.
 * set - The set this sequence belongs to.
 */
void SESA_add_sequence(SESA sesa, char *seq, unsigned int set) {
  unsigned int nSeq = sesa->nSeq;

  if (nSeq == sesa->maxSeq-1) {
    sesa->maxSeq *= 2;
    sesa->seqnames = (char **) realloc((char **) sesa->seqnames, sizeof(char *) * sesa->maxSeq);
    sesa->set = (unsigned int *) realloc((unsigned int *) sesa->set, sizeof(unsigned int) * sesa->maxSeq);

    if (sesa->set == NULL) {
      fprintf(stderr, "Memory reallocation failed, out of memory? Sequences: %d\n", sesa->nSeq);
      exit(EXIT_FAILURE);
    }
  }

/*   if (nSeq > 0)  */
/*   fprintf(stderr, "-NSEQ(%i) of (%i)  === (%i)  (%X)\n", nSeq-1, sesa->maxSeq, sesa->set[nSeq-1], sesa->set);  */

  sesa->set[nSeq] = set;
/*   fprintf(stderr, "NSEQ(%i) of (%i)  === (%i)   (%X)\n", nSeq, sesa->maxSeq, sesa->set[nSeq], sesa->set); */
  sesa->nSeqSet[set]++;
  sesa->seqnames[nSeq] = (char *) malloc(sizeof(char) * (strlen(seq) + 1));
  strcpy(sesa->seqnames[nSeq], seq);
  sesa->nSeq++;
}


/**
 * Builds a wrapper for an ESA. 
 *
 * files - Number of files.
 * filenames - Array of filnames and paths.
 * pAlphabet - The alphabet the files consist of [ACGT].
 * pIgnore - The alphabet that should be ignored [N].
 * ds - Should the tree be double stranded.
 */
SESA build_SESA(int files, char **filenames, char *pAlphabet, char *pIgnore, int ds) {
  SESA sesa = init_SESA();
  char *sequences;
  FILE *file[files];
  int i, seq_size, size = 1;
  char namebuffer[MAX_NAME_SIZE+1];
  int c, n;
  int start_size = 0, start_seq = 0;
  char *newIgnore;

  sesa->nSeqSet = (unsigned int *) malloc(sizeof(unsigned int) * files);

  for (i = 0; i < files; i++) {
    sesa->nSeqSet[i] = 0;

    if(!filenames[i]) {
      setError("Need to supply a filename for each file");
      return NULL;
    }

    /* Open sequencefile */
    if(!(file[i] = fopen(filenames[i], "r"))) {
      char st[512];
      sprintf(st, "Could not open '%s' for reading.", filenames[i]);
      setError(st);
      return NULL;
    }

    /* Check file size. Will overestimate by size of labels, but who cares. */
    fseek(file[i], 0, SEEK_END);
    size += ftell(file[i]);
    rewind(file[i]);
  }

  if (ds) size *= 2;

/*   fprintf(stderr, "SIZE is %i\n",size); */

  sequences = (char *) malloc(sizeof(char) * size);
  sesa->seq = (unsigned int *) malloc(sizeof(unsigned int) * size);
  sesa->pos = (unsigned int *) malloc(sizeof(unsigned int) * size);


/*   unsigned int tmp = sesa->set; */
  /* Files to string */
  size = 0;
  seq_size = 0;

  for (i = 0; i < files; i++) {
    
    while ((c = fgetc(file[i])) != EOF) {     

      if (c == '\n') continue;

      /* Fasta header */
      if (c == '>') {
	if (sesa->nSeq) sequences[size++] = SEP_CHAR;
	
	/* Lasts till end of line (or til MAX_NAME_SIZE) */
	n = 0;
	while ((c = fgetc(file[i])) != '\n') { 
	  if (n < MAX_NAME_SIZE) namebuffer[n++] = c;
	}
	namebuffer[n++] = '\0';

	SESA_add_sequence(sesa, namebuffer, i);	
/* 	tmp = (unsigned int) sesa->set; */

	seq_size = 0;
	continue;
      }

      sesa->seq[size] = sesa->nSeq-1;

      sesa->pos[size] = seq_size;
      sequences[size] = toupper(c);      


      if (isACGT(c)) {
	sesa->bg[chr2ind(c)]++;
      } else if (notN(c)) {
	  fprintf(stderr, "Character not in alphabet: %c\n", c);
      }

      seq_size++;
      size++;
    } 
    
    fclose(file[i]);

    /* Create reverse complement */
/*     if (ds) size = create_reverse_strand(sesa, sequences, size, start_size, start_seq); */
    start_size = size;
    start_seq = sesa->nSeq-1;
  }


  sequences[size] = '\0';
  newIgnore = malloc(sizeof(char) * strlen(pIgnore) + 2);
  strcpy(newIgnore, pIgnore);
  strcat(newIgnore, "|");

  sesa->sequences = sequences;

  /* Build the ESA */
  sesa->esa = build_ESA(sequences, size, pAlphabet, newIgnore, 0); 
  free(newIgnore);
  
  return (sesa->esa ? sesa : NULL);
}


/* int create_reverse_strand(SESA sesa, char *sequences, int size, int start_size, int seq) { */
/*   int forward_size = size; */
/*   int i, start; */
/*   char *revcomp = "REVCOMP: "; */
/*   char *name; */

/*   sequences[size++] = SEP_CHAR;  */
/*   start = size; */

/*   for (i = start_size; i < forward_size; i++) { */
/*     if (sequences[i] == SEP_CHAR) { */
/*       name = (char *) malloc(strlen(sesa->seqnames[seq]) + strlen(revcomp) + 1); */
/*       strcpy(name, revcomp); */
/*       strcat(name, sesa->seqnames[seq]); */
/*       SESA_add_sequence(sesa, name, sesa->set[seq]);  */
/*       free(name);  */
/*       seq++; */

/*       /\* Reverse the string *\/ */
/*       reverse(sequences, start, size-1); */
/*       start = size + 1; */
/*     } */

/*     sesa->seq[size] = sesa->nSeq-1; */
/*     sesa->pos[size] = size-start; */
/*     sequences[size] = complement(sequences[i]); */
/*     size++; */
/*   } */

/*   name = (char *) malloc(strlen(sesa->seqnames[seq]) + strlen(revcomp) + 1); */
/*   strcpy(name, revcomp); */
/*   strcat(name, sesa->seqnames[seq]); */
/*   SESA_add_sequence(sesa, name, sesa->set[seq]); */
/*   reverse(sequences, start, size-1); */

/*   return size; */
/* } */


/* void reverse(char *sequences, int start, int end) { */
/*   char tmp; */

/*   while( start < end ) { */
/*     tmp = sequences[start]; */
/*     sequences[start] = sequences[end]; */
/*     sequences[end] = tmp; */
	  
/*     ++start; */
/*     --end; */
/*   } */
/* } */
