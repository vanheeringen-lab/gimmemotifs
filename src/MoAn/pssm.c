/*---------------------------------------------------------------*
 * PSSM Searcher                                                 *
 *                                                               *
 * By: Jes Frellsen, Ida Moltke and Martin Thiim                 *
 * University of Copenhagen, spring 2006                         *
 *---------------------------------------------------------------*/

#include "memcheck.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "pssm.h"
#include "error.h"

#define min(x, y) (x < y ? x : y)

// Function prototypes
inline int map(const unsigned char *letters, int length, int alphabetSize);
inline int powI(int x, int y);


// ******************************************************************
// Functions for creating and mutating a pssm 
// ******************************************************************

// Function that creates and allocates a matrix (except for the scores) 
PSSM initMatrix(int order, int length, int alphabetSize){

  // Check if longer than max-length of PSSMs
  if(length >= MAXPSSMSIZE) {
    setError("Matrix is too long.");
    return NULL;
  }

  // All attributes are initialised
  PSSM pssm = malloc(sizeof(struct PSSM));
  if(!pssm){setError("Couldn't allocate memory for the PSSM.");return NULL;};

  pssm->order = order;
  pssm->length = length;
  pssm->alphabetSize = alphabetSize;

  calcAndSetOffsets(pssm, order, length, alphabetSize);

  // Allocate space for the scores
  pssm->scores = malloc(pssm->offsets[length] * sizeof(double));
  if(!pssm->scores) {free(pssm); setError("Couldn't allocate memory for the PSSM.");return NULL;}

  /* Allocate space for the counts */
  pssm->counts = malloc( pssm->offsets[length] * sizeof(int));
  if(!pssm->counts) {free(pssm); setError("Couldn't allocate memory for the PSSM."); return NULL;}

  return pssm;
}


// Creates a matrix given an array of scores.
//
// The function also checks if the length is too high and if numbers of
// scores in '*scores' array matches order, length and alphabetSize
// (returns NULL is check fails)
PSSM initMatrixScore(int order, int length, int alphabetSize, double *scores, int nScores, double threshold){

  // Check if longer than max-length of PSSMs
  if(length >= MAXPSSMSIZE) {
    setError("Matrix is to long.");
    return NULL;
  }

  // All attributes are initialised
  PSSM pssm = malloc(sizeof(struct PSSM));
  if(!pssm){setError("Couldn't allocate memory for the PSSM."); return NULL;};
  pssm->order = order;
  pssm->length = length;
  pssm->alphabetSize = alphabetSize;
  pssm->scores = scores;

  calcAndSetOffsets(pssm,order,length,alphabetSize);

  // Check if length of scores is correct
  if(nScores != pssm->offsets[length]) {
    free(pssm);
    char errormsg[160];
    sprintf(errormsg, "Mismatch between list size (%i) and size calculated from order, length and alphabet size (%i).",
	    nScores, pssm->offsets[length]);
    setError(errormsg);
    return NULL;
  }

  calcAndSetThresholds(pssm, threshold);

  return pssm;
}


// Function that frees the memory of a matrix
void releaseMatrix(PSSM pssm){
  if(pssm){
    if(pssm->scores)
      free(pssm->scores);
    if(pssm->counts)
      free(pssm->counts);
    free(pssm);
  }
}


// Function that fills in a given pssms offset array
void calcAndSetOffsets(PSSM pssm, int order, int length, int alphabetSize){
  int i;
  int offset, realorder;

  // The offset for position 0 is set to 0
  offset = 0;
  pssm->offsets[offset] = 0;

  // The offset for position 'i' is set to the offset for posi-
  // tion 'i-1' plus the number of scores for position 'i-1' 
  for(i = 1; i <= length; i++){
    realorder = min(i,order+1);
    offset += powI(alphabetSize,realorder);
    pssm->offsets[i] = offset;
  }
}


// Function that fills in a given pssms threshold array
void calcAndSetThresholds(PSSM pssm, double threshold){
  int pos, key, offset, numScores;
  double last, curMax, curScore;

  // For each position 'pos' (starting with the last) 
  last = threshold;
  for(pos = getLength(pssm) - 1; pos >= 0; pos--){
    curMax = -HUGE_VAL;
    numScores = pssm->offsets[pos+1] - pssm->offsets[pos];
    offset = pssm->offsets[pos];

    // The highest score is found 
    for(key = 0; key < numScores; key++) {
      curScore = pssm->scores[offset+key];
      if(curScore > curMax) {
	curMax = curScore;
      }
    }

    // The threshold for 'pos' is set to 'threshold' minus
    // minus the sum of the max scores on all positions after
    pssm->thresholds[pos] = last;
    last -= curMax;     
  }
}


// ******************************************************************
// Accessors for pssms
// ******************************************************************

// Function that returns the score of a given letter in a suffix
double getScore(PSSM pssm, const unsigned char *baseLetter, int pos){
  int index = 0;
  int realorder = min(pos, pssm->order);
  unsigned char *letters = (unsigned char *)(baseLetter-realorder);

  // The index of the letters are calculated
  index =  pssm->offsets[pos];
  index += map(letters, realorder+1, pssm->alphabetSize);

  return pssm->scores[index];
}


// Function that returns the threshold of a given pos in a given pssm
double getThreshold(PSSM pssm, int pos){
  return pssm->thresholds[pos];
}


// Function that returns the global threshold for a given pssm
double getGlobalThreshold(PSSM pssm){
  return pssm->thresholds[pssm->length - 1];
}


// Function that returns the length of a given pssm
int getLength(PSSM pssm){
  return pssm->length;
}

/******************************************************************
 * Setters for PSSMs
 ******************************************************************/

// Function that sets a specific score in a given pssm
void setScore(PSSM pssm, const unsigned char *letters, int pos, double score){
  int index = 0;
  int realorder = min(pos, pssm->order);

  index =  pssm->offsets[pos];
  index += map(letters, realorder+1, pssm->alphabetSize);

  pssm->scores[index] = score;
}

void setCount(PSSM pssm, const unsigned char *baseLetter, int pos, unsigned int newValue){
  int index = 0;
  int realorder = min(pos, pssm->order);
  unsigned char *letters = (unsigned char *)(baseLetter-realorder);

  /* The index of the letters are calculated */
  index =  pssm->offsets[pos];
  index += map(letters, realorder+1, pssm->alphabetSize);

  pssm->counts[index] = newValue;
}


// ******************************************************************
// Additional functions
// ******************************************************************

// Function for calculation small powers of intergers
int powI(int x, int y) {
  int i;
  int res = 1;
  
  for(i = 1; i<=y; i++)
    res = res*x;

  return res;
}


// Fuctions for mapping from char array to integer
inline int map(const unsigned char *letters, int length, int alphabetSize) {
  int k;
  int value = 0;
  
  for(k = length-1; k >= 0; k--){
    value += (*(letters + k)) * powI(alphabetSize, (length - 1) - k);
  }

  return value;
}

/******************************************************* NEW STUFF *********************************************/


/**
 * Loads a count matrix from file and converts it to log-odds scores
 */
PSSMSet load_count_matrices(char *filename, int alphlen, double log_zero) {
  PSSMSet pssmset = (PSSMSet) malloc(sizeof(struct PSSMSet));
  PSSM pssm;
  FILE *file;
  int i, j, k, order, length;


  file = fopen(filename, "r");
  fscanf(file, "# %d\n", &i);
  fscanf(file, "> %d %d\n", &order, &length);

/*   fprintf(stderr, "%d %d %d\t",i, order, length); */

  pssmset->pssms = (PSSM *) malloc(2 * i * sizeof(PSSM));
  pssmset->pssm_count = i;

  pssmset->pssms[0] = initMatrix(order, length, alphlen);

  pssm = pssmset->pssms[0];
  for (j = 0; j < alphlen; j++) {
    for (k = 0; k < length; k++) {
      fscanf(file, "%d", &pssm->counts[j * length + k]);
/*       fprintf(stderr, "%d\t", pssm->counts[j * length + k]); */
      pssm->scores[j * length + k] = (pssm->counts[j * length + k] > 0 ? log(((double) pssm->counts[j * length + k])) : log_zero);
    }
/*       fprintf(stderr, "\n"); */
  }
  fclose(file);

  pssmset->pssms[1] = clone_pssm(pssmset->pssms[0]);

  return pssmset;
}

/**
 * Loads a log matrix from file
 */
PSSM load_log_matrices(char *filename, int alphlen) {
  PSSM pssm;
  FILE *file;
  int i, j, k, order, length;

  file = fopen(filename, "r");
  fscanf(file, "# %d\n", &i);
  fscanf(file, "> %d %d\n", &order, &length);
  
  pssm = initMatrix(order, length, alphlen);
  for (j = 0; j < alphlen; j++) {
    for (k = 0; k < length; k++) {
      fscanf(file, "%lf", &pssm->scores[k * alphlen + j]);
    }
  }
  fclose(file);

  return pssm;
}


/* 
 * Function for printing a PSSM
 * pssm - The PSSM to print
 *
 * FIX: Only zeroth-order
 */
void print_counts(PSSM pssm) {
  int ltr, pos, ord;
  int alphlen = pssm->alphabetSize;

/*   fprintf(stderr, "Order: %i  Length: %i  Size of alphabet: %i\n", pssm->order, pssm->length, pssm->alphabetSize); */

  for (pos = 0; pos < pssm->length; pos++) {      
    printf("%7i", pos+1);
  }
  printf("\n");

  for (ltr = 0; ltr < pssm->alphabetSize; ltr++) {
    ord = pssm->order;

    for (pos = 0; pos < pssm->length; pos++) {      
      printf("%7i", pssm->counts[pos * alphlen + ltr]);
    }
    printf("\n");
  }
}



/* 
 * Function for printing a PSSM
 * pssm - The PSSM to print
 *
 * FIX: Only zeroth-order
 */
void print_pssm(PSSM pssm) {
  int ltr, pos, ord;
  int alphlen = pssm->alphabetSize;

/*   fprintf(stderr, "Order: %i  Length: %i  Size of alphabet: %i\n", pssm->order, pssm->length, pssm->alphabetSize); */

  for (pos = 0; pos < pssm->length; pos++) {      
    printf("%7i", pos+1);
  }
  printf("\n");

  for (ltr = 0; ltr < pssm->alphabetSize; ltr++) {
    ord = pssm->order;

    for (pos = 0; pos < pssm->length; pos++) {      
      printf("%7.2f", pssm->scores[pos * alphlen + ltr]);
    }
    printf("\n");
  }

  printf("Thresholds: \n");
  for (pos = 0; pos < pssm->length; pos++) {      
    printf("%7.2f", pssm->thresholds[pos]);
  }
  printf("\n");
}


/**
 * Returns a pointer to a new PSSM with the same content as the
 * parameter PSSM.
 */
PSSM clone_pssm(PSSM pssm) {
  PSSM copy = initMatrix(pssm->order, pssm->max_length, pssm->alphabetSize);
  copy_pssm(pssm, copy);

  return copy;
}


/**
 * Copies the content of one PSSM to another.
 */
void copy_pssm(PSSM old, PSSM new) {
  int i;

  i = old->length; 
  while (i--) {
    new->thresholds[i] = old->thresholds[i]; 
/*     new->offsets[i] = old->offsets[i];  */
  }
  
  i = old->offsets[old->length];
  while (i--) {
    new->scores[i] = old->scores[i];
    new->counts[i] = old->counts[i];
  }

  new->alphabetSize = old->alphabetSize;
  new->order = old->order;
  new->length = old->length;
}

/**
 * Copies the reverse complement of one PSSM to another.
 */
void revcomp_copy_pssm(PSSM old, PSSM new) {
  int i, max;
  
  max = i = old->offsets[old->length];
  while (i--) {
    new->scores[max - i - 1] = old->scores[i];
    new->counts[max - i - 1] = old->counts[i];
  }

  new->alphabetSize = old->alphabetSize;
  new->order = old->order;
  new->length = old->length;

  calcAndSetThresholds(new, old->thresholds[old->length-1]);
}


void releasePSSMSet(PSSMSet pssmset) {
  while (pssmset->pssm_count--) 
    releaseMatrix(pssmset->pssms[pssmset->pssm_count]);
  free(pssmset->pssms);
  free(pssmset);
}

