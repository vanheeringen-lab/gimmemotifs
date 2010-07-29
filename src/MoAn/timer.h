/*---------------------------------------------------------------*
 * PSSM Searcher                                                 *
 *                                                               *
 * By: Jes Frellsen, Ida Moltke and Martin Thiim                 *
 * University of Copenhagen, spring 2006                         *
 *---------------------------------------------------------------*/

#ifndef TIMER_H_
#define TIMER_H_

#ifndef TIMER
#define initTimer() 
#define addTimer(text) 
#define printTimer() 
#define printTimerWithoutTotal() 
#define freeTimer() 

#else

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void initTimer(void);
void addTimer(const char *);
void printTimer(void);
void printTimerWithoutTotal(void);
void freeTimer(void);

#endif
#endif
