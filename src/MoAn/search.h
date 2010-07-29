/*---------------------------------------------------------------*
 * PSSM Searcher                                                 *
 *                                                               *
 * By: Jes Frellsen, Ida Moltke and Martin Thiim                 *
 * University of Copenhagen, spring 2006                         *
 *---------------------------------------------------------------*/

#ifndef SEARCH_H_
#define SEARCH_H_

#include "shared.h"
#include "pssm.h"

/** This function searches a string with a given PSSM 
 * using a preconstructed enhanced suffix array
 */ 
struct HitTable* search(ESA esa, PSSM pssm);


/** This function searches a string with a given PSSM 
 * using the naïve method
 */ 
struct HitTable* searchNaively(unsigned char* pStr, int numSufs, PSSM pssm);


/** This function searches a string with a given PSSM 
 * using the naïve method with lookahead 
 */ 
struct HitTable* searchNaivelyWithLA(unsigned char* pStr, int numSufs, PSSM pssm);

#endif
