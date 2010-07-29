/*---------------------------------------------------------------*
 * PSSM Searcher                                                 *
 *                                                               *
 * By: Jes Frellsen, Ida Moltke and Martin Thiim                 *
 * Modified by: Eivind Valen                                     *
 * University of Copenhagen, spring 2006                         *
 *---------------------------------------------------------------*/

#ifndef PSSM_H_
#define PSSM_H_

#include "shared.h"


/* -------------------------------------------------------------------
 *  Definition of the pssm used for scoring
 * -------------------------------------------------------------------
 *
 * The attributes are: 
 *
 * length       : the length of the pssm 
 * order        : the order of the pssm
 * alphabetSize : the size of the alphabet the pssm supports
 *                (the alphabet is assumed to consists of chars
 *                with ASCII code 0,1,2...)
 * scores       : the scores of the pssm. They are kept in a one
 *                dim array with all the scores for the first 
 *                position first then all ofg those for the second
 *                and so forth. 
 *                E.g if it is a second order pssm for a four letter 
 *                alphabet the scores will be in the order correspon-
 *                ding to
 *
 *                pos 
 *                0: 0,1,2,3,
 *                1: 00,01,02,03,10,11,12,13,20,21,22,23,30,31,32,33,
 *                2: 000,001,002,003,010,011,012,013,020,021,022,...
 *                3: 000,001,002,003,010,011,012,013,020,021,022,...
 *                4: 000,001,002,003,010,011,012,013,020,021,022,...
 *                ...
 *
 * offsets       : an array the for each position in the pssm contains
 *                 the index of the first spot in the score array that
 *                 contains scores for that position
 * thresholds    : an array that contains the threshold for each 
 *                 position. The one for position length-1 is of
 *                 course also the global threshold
 *
 * ------------------------------------------------------------------- */


typedef struct PSSM {
  char           *name;         /* Name if any of the WM */
  unsigned char  order;         /* Order of the WM, currently always 0 */
  unsigned char  max_length;    /* Maximum length for WM (allocated space) */
  unsigned char  min_length;    /* Minimum allowed length for WM */
  unsigned char  length;        /* Length of the WM */
  unsigned char  alphabetSize;  /* Length of the alphabet, always 4 */
  double         *scores;       /* The score matrix */
  unsigned int   *counts;       /* The count matrix */
  double         thresholds[MAXPSSMSIZE + 1];   /* Data for lookahead search */
  int            offsets[MAXPSSMSIZE];
  int            last_step;     /* Last operation performed on the WM */
  int            last_col;      /* Number of cols changed by last op if INC/DEC*/
} *PSSM;

/* typedef struct PSSM* PSSM;  */

typedef struct PSSMSet {
  unsigned int pssm_count;
  PSSM *pssms;
} *PSSMSet;

/* typedef struct PSSMSet* PSSMSet;  */

#define LOG_ZERO  -10


// Mutators
PSSM initMatrix(int order, int length, int alphabetSize);
PSSM initMatrixScore(int order, int length, int alphabetSize, double *scores, int nScores, double threshold);
void releaseMatrix(PSSM pssm);
void calcAndSetOffsets(PSSM pssm, int order, int length, int alphabetSize);
void calcAndSetThresholds(PSSM pssm, double threshold);
void setScore(PSSM pssm, const unsigned char *letters, int pos, double score);


// Acsessors
double getScore(PSSM pssm, const unsigned char *baseLetter, int pos);
double getThreshold(PSSM pssm, int pos);
double getGlobalThreshold(PSSM pssm);
int getLength(PSSM pssm);


// Fast accessors
#define getThresholdFast(pssm, pos)  (pssm)->thresholds[pos]
#define getGlobalThresholdFast(pssm) (pssm)->thresholds[(pssm)->length - 1]
#define getLengthFast(pssm)          (pssm)->length
#define getNumScores(pssm)           (pssm)->offsets[(pssm)->length]
#define getScoreFast(pssm, baseLetter, pos) ( (pssm)->order == 0 ? (pssm)->scores[(pssm)->offsets[pos] + (int)*(baseLetter)] : getScore(pssm, baseLetter, pos) )


/***** NEW *****/

void releasePSSMSet(PSSMSet pssmset);

/* Read from file */
PSSM load_log_matrices(char *filename, int alphlen);
PSSMSet load_count_matrices(char *filename, int alphlen, double log_zero);

/* Random matrices */
PSSMSet random_pssms(unsigned int count, const unsigned char order,  const unsigned char min_length, const unsigned char max_length, const unsigned char alphlen, const int seq);
PSSM random_pssm(unsigned char order,  const unsigned char min_length, const unsigned char max_length, unsigned char alphlen, int seq);

/* DIV */
void copy_pssm(PSSM old, PSSM new);
PSSM clone_pssm(PSSM pssm);
void revcomp_copy_pssm(PSSM old, PSSM new);
void print_pssm(PSSM pssm);
void print_counts(PSSM pssm);

/* Setters */
void setScore(PSSM pssm, const unsigned char *baseLetter, int pos, double newValue);
void setCount(PSSM pssm, const unsigned char *baseLetter, int pos, unsigned int newValue);

/* Fast setters */
#define setScoreFast(pssm, baseLetter, pos, new) ( (pssm)->order == 0 ? (pssm)->scores[(pssm)->offsets[pos] + (int)*(baseLetter)] = new : setScore(pssm, baseLetter, pos, new) )
#define setCountFast(pssm, baseLetter, pos, new) ( (pssm)->order == 0 ? (pssm)->counts[(pssm)->offsets[pos] + (int)*(baseLetter)] = new : setCount(pssm, baseLetter, pos, new) )

#endif
