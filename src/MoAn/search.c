/*---------------------------------------------------------------*
 * PSSM Searcher                                                 *
 *                                                               *
 * By: Jes Frellsen, Ida Moltke and Martin Thiim                 *
 * University of Copenhagen, spring 2006                         *
 *---------------------------------------------------------------*/

#include "memcheck.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "search.h"
#include "hittable.h"
#include "error.h"
#include "timer.h"

inline void pfxsc(double C[], int *d, int l, PSSM pssm, unsigned char *base, int offset);


struct HitTable* search(ESA esa, PSSM pssm) {
/*   initTimer(); */

  // Declaration of variables
  int i,j;
  int curEsaIndex, pssmLength, curSuf;
  int lcpi;
  struct HitTable *pHitTable;

  // Getting the string and stringsize from the ESA
  unsigned char *pStr = getStr(esa);
  int numSufs = getSize(esa);

  // Initialisation of d, C and the scoring system
  // dCur and dPrev keeps track of the number of letters scored
  // before the threshold is/was met in the current and in
  // the previous suffix (= matrixlength, suffix is a hit)
  int dCur;
  int dPrev; 

  // - C is an array; for suffix i the entry j contains the 
  //   score for the first j-1 letters -- the prefix scores
  //   (it is only filled for entries less than or equal d)
  
  double *Ci = (double*) malloc(sizeof(double) * MAXPSSMSIZE);
  if(!Ci) {setError("Couldn't allocate memory for scores in the search.");return NULL;}
  init_hittable(); 
  pssmLength = getLengthFast(pssm);


  // First d and C are calculated for the 1st suffix in the array:
  curEsaIndex = 0;
    
  curSuf = getSuf(esa,curEsaIndex);
  pfxsc(Ci, &dCur, pssmLength-1, pssm, &pStr[curSuf], curEsaIndex);
    
  // if the suffix is a hit it is registered in the hittable
  if(dCur == pssmLength-1){
    register_hit(curSuf, Ci[pssmLength-1]);
  }
        

  // The d and C are calc for the rest of the suffixes: 
  curEsaIndex++;    
  dPrev = dCur;
  while(curEsaIndex<numSufs){
    // For each suffix lcp between suf i-1 and suf i in the 
    // array is checked
    lcpi = getLcp(esa,curEsaIndex);
    
    // NOTE: It a special case, if d[curEsaIndex-1]+1 == lcpi ==
    // pssmLength and then we can also skip.
    int special = dPrev+1 == lcpi && lcpi == pssmLength;
    
    // If the lcp =< di-1+1 (the number of letters 
    // succesfully scored for suf i-1) then the local scores
    // are calculated from the point of lcp (as the local 
    // scores must be the same up until the lcp)
    if( dPrev+1 >= lcpi && !special ){
       curSuf = getSuf(esa,curEsaIndex);	
       pfxsc(Ci, &dCur, pssmLength-1, pssm, &pStr[curSuf], lcpi);
      
       if(dCur == pssmLength-1){
	      register_hit(curSuf, Ci[pssmLength-1]);
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
	     }
      }
      curEsaIndex = j;
    }
  }
  
  pHitTable = finish();
  
  // Deallocation 
  free(Ci);
  
  // Print the timer
/*   addTimer("ESA Search"); */
/*   printTimer(); */
/*   freeTimer(); */

  return pHitTable;
}	

// Function for calculation Ci[d] (the prefix score) and di.
//
// The function assumes that the C[offset-1] is correct and that pssm
// is a single PSSM.
// C_i is updated from offset and up to l, as long a the score is
// above the threshold and the string separator isn't met.
inline void pfxsc(double C[], int *d, int l, PSSM pssm, unsigned char *base, int offset) {
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




struct HitTable* searchNaively(unsigned char* pStr, int numSufs, PSSM pssm){
  initTimer();

  // Declaration of variables
  int pssmLength, curSuf, curLetterPos;
  double curThreshold, score;
  unsigned char *firstLetter, *curLetter;
  bool hit;
  struct HitTable *pHitTable;

  // Initialisation of scores to return
  init_hittable(); 

  // Score all suffixes are scored one at a time:
  pssmLength = getLengthFast(pssm);
  
  for(curSuf = 0; curSuf < numSufs; curSuf++){
    score = 0;
    firstLetter = pStr+curSuf;
    hit = true;

    for(curLetterPos = 0; curLetterPos < pssmLength; curLetterPos++){
      curLetter = firstLetter+curLetterPos;	
	      
      // If a string seperator is met the suffix is too short.
      // The calculations are therefore stopped and the suffix 
      // is given the score -1
      if(*curLetter == BREAKSYM){
	hit = false;
	break;
      }
      else{
	score += getScoreFast(pssm,curLetter,curLetterPos);
      }				
    }
    
    // If the score is not below the threshold it is registered    
    curThreshold = getGlobalThresholdFast(pssm);
    if(score>=curThreshold && hit){
      register_hit(curSuf, score);
    }    
  }

  pHitTable = finish();


  // Print timer
  addTimer("Search NA");
  printTimer();
  freeTimer();

  return pHitTable;
}




struct HitTable* searchNaivelyWithLA(unsigned char* pStr, int numSufs, PSSM pssm){

  initTimer();

  // Declaration of variables
  int pssmLength, curSuf, curLetterPos;
  double curThreshold, score;
  unsigned char *firstLetter, *curLetter;
  bool hit;
  struct HitTable *pHitTable;

  // Allocation
  init_hittable(); 
  
  // All suffixes are scored one at a time:
  pssmLength = getLengthFast(pssm);
        
  for(curSuf = 0; curSuf < numSufs; curSuf++){
    score = 0;
    firstLetter = pStr+curSuf;
    hit = true;

    for(curLetterPos = 0; curLetterPos < pssmLength; curLetterPos++){
      curLetter = firstLetter+curLetterPos;
      
      // If a string seperator is met the suffix is too short.
      // The calculations are therefore stopped and the suffix 
      // is given the score -1
      if(*curLetter == BREAKSYM){
	hit = false;
	break;
      }
      else{
	score += getScoreFast(pssm,curLetter,curLetterPos);
	curThreshold = getThresholdFast(pssm,curLetterPos);
	
	// If a local score below the local threshold is met 
	// the suffix is also given a score of -1 
	if(curThreshold > score){
	  hit = false;
	  break;
	}
      }
    }
    if (hit){ 
      register_hit(curSuf, score);
    }
  }

  pHitTable = finish();   

  // Print timer
  addTimer("Search NA with LA");
  printTimer();
  freeTimer();

  return pHitTable;
}
