/*---------------------------------------------------------------*
 * PSSM Searcher                                                 *
 *                                                               *
 * By: Jes Frellsen, Ida Moltke and Martin Thiim                 *
 * University of Copenhagen, spring 2006                         *
 *---------------------------------------------------------------*/

#define TIMER

#include "memcheck.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <sys/times.h>
#include <sys/resource.h>
#include "timer.h"

#define MAXTIMERS 100
#define CPUCLOCK  1000000.0


struct timer {
  double time;
  const char *name;
};

double getTime ( void )
{
   double usertime,systime;
   struct rusage usage;

   getrusage ( RUSAGE_SELF, &usage );

   usertime = (double)usage.ru_utime.tv_sec +
     (double)usage.ru_utime.tv_usec / CPUCLOCK;

   systime = (double)usage.ru_stime.tv_sec +
     (double)usage.ru_stime.tv_usec / CPUCLOCK;

   return(usertime+systime);
}



struct timer *globTimer = NULL;
int globTimerCount = 0;

void initTimer() {
  if(globTimer)
    freeTimer();
  
  globTimer = malloc(sizeof(struct timer)*MAXTIMERS);
  
  globTimer[globTimerCount].name = "Start";
  globTimer[globTimerCount].time = getTime();
}

void addTimer(const char *name) {
  if(globTimerCount < MAXTIMERS-1) {
    globTimerCount++;

    globTimer[globTimerCount].name = name;
    globTimer[globTimerCount].time = getTime();
  }
}


void printTimer() {
  int i;
  double diff, total;
  
  printf("--------------------------------\n");
  for(i = 1; i <= globTimerCount; i++) {
    diff = globTimer[i].time - globTimer[i-1].time;
    printf("%.2f s\t%s\n", diff, globTimer[i].name);
  }

  total = globTimer[globTimerCount].time - globTimer[0].time;
  printf("%.2f s\t%s\n",total,"TOTAL");
  printf("--------------------------------\n");
}


void printTimerWithoutTotal() {
  int i;
  double diff;

  for(i = 1; i <= globTimerCount; i++) {
    diff = globTimer[i].time - globTimer[i-1].time;
    printf("%s\t%.2f seconds\n", globTimer[i].name,diff);
 
  }
}

void freeTimer() {
  if(globTimer)
    free(globTimer);

  globTimerCount = 0;
  globTimer = NULL;
}

