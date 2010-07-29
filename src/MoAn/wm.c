#define TIMER

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "error.h"
#include "eesa.h"
#include "sesa.h"
#include "search.h"
#include "timer.h"

#define RAND_INT(max)   (int)   (((float) max) * rand() / (RAND_MAX + 1.0))

int main(int argc, char **argv) {
/*   Hits hits; */
/*   Hit *hit; */
  EESA eesa;
  PSSM pssm;
  unsigned int i;
/*   unsigned int p; */
  struct HitTable *s_hits;


  srand((unsigned)time(NULL));  

  pssm = initMatrix(0, 5, 4);
  pssm = load_log_matrices("../test/small-wm", 4);
  print_pssm(pssm); 


  /* BUILD ESA */
  eesa = load_fasta_DNA("../test/PR1.real.100.0", 0);
  calcAndSetThresholds(pssm, -0.5);

  if (!eesa->esa) {
    fprintf(stderr, "ERROR: %s\n", getError());
  }

  initTimer();
  for (i = 0; i < 100000; i++) {
    setMismatchScores(pssm,RAND_INT(255));
    s_hits = search(eesa->esa, pssm);
    release_hits(s_hits);
  }
  addTimer("ESA");
  printTimer();
  freeTimer();


  /* SEARCH */
/*   initTimer(); */
/*   for (i = 0; i < 1000000; i++) { */
/*     hits = single_search(eesa, pssm); */
/*   } */
/*   addTimer("EESA"); */
/*   printTimer(); */
/*   freeTimer(); */

/*   for (i = 0; i < eesa->nSeq; i++) { */
/*     printf("BORDER: %i is %u\n", i, eesa->seqborders[i]); */
/*   } */

/*   printf("HIT COUNT: %i\n", hits->nScores); */
/*   for (i = 0; i < hits->nScores; i++) { */
/*     hit = &hits->pScores[i]; */
/*     printf("HIT: %3i SEQ: %3i POS: %5i ==> %lf  \n", i, hit->seq, hit->pos, hit->score); */
/*   } */


  SESA sesa = EESA2SESA(eesa);
/*   p = (unsigned int) strtol(argv[1], (char **) NULL, 0); */
/*   shits = SESA_exhaustive_search(sesa, p, 0, 4, 50);  */

  initTimer();
  init_hittable();
  for (i = 0; i < 100000; i++) {
    setMismatchScores(pssm,RAND_INT(255));
    s_hits =  SESA_search(sesa, pssm);
    reset_hits();
  }
  release_hits(s_hits);
  addTimer("SESA");
  printTimer();
  freeTimer();

/*   printf("HIT COUNT: %i\n", shits->nScores); */
/*   for (i = 0; i < shits->nScores; i++) { */
/*     hit = &shits->pScores[i]; */
/*     printf("HIT: %3i SEQ: %3i POS: %5i ==> %lf  \n", i, hit->seq, hit->pos, hit->score); */
/*   } */

  exit(0);
}







/* int compare_pEntries(const void *pEntry1, const void *pEntry2){ */
/*   return ((struct HitEntry *)pEntry1)->position - ((struct HitEntry *)pEntry2)->position; */
/* } */

/* void writeResultsToFile(char *file, struct HitTable *pScoreData) { */
/*   int i, j, k; */
/*   struct HitEntry *pEntry; */
/*   FILE *f; */

/*   i = 0; */
/*   pEntry = pScoreData->pScores; */
/*   k = pScoreData->nScores; */
/*   printf("\nResult from using ESA search:\n"); */
/*   printf("Matrix %d has %d entries\n", i, k); */
 
/*   if(k<100){ */
/*     qsort((void*)pEntry,k,sizeof(struct HitEntry),compare_pEntries); */
/*     for(j = 0; j < k; j++){ */
/*       printf("Entry %d has pos %d ", j, pEntry[j].position); */
/*       printf("and score %f\n", (float)pEntry[j].score); */
/*     } */
/*   } */
/*   else{ */
/*     f = fopen(file, "wt"); */
    
/*     if(!f) { */
/*       printf("Could not write file %s\n", file); */
/*     } */
/*     else { */
/*       qsort((void*)pEntry,k,sizeof(struct HitEntry), compare_pEntries); */
/*       for(j = 0; j < k; j++){ */
/* 	fprintf(f,"pos %d ", pEntry[j].position); */
/* 	fprintf(f,"score %f\n", (float)pEntry[j].score); */
/*       } */
/*       fclose(f); */
/*       printf("Result can be seen in %s\n", file); */
/*     } */
/*   } */
/* } */


/*   exit(0); */

/*   fprintf(stderr, "Searching...\n"); */
/*   results = search(esa, pssm); */
/*   writeResultsToFile("esa.res-2", results); */
  

/*   fprintf(stderr, "Number of hits: %i\n", results->nScores); */
/*   i = results->nScores; */
/*   while (i--) { */
/*     fprintf(stderr, "%i ::: %f\n", results->pScores[i].position, results->pScores[i].score); */
/*   } */

/*   normalizeAllCounts(pssm, 1000);  */
/*   calc_sufmax(cur, cur->length - 1); */


/*   pssm = initMatrix(0, 4, 4); */

/*   unsigned char LETTER_A = 0; */
/*   unsigned char LETTER_C = 1; */
/*   unsigned char LETTER_G = 2; */
/*   unsigned char LETTER_T = 3; */
  
/*   setScore(pssm, &LETTER_A, 0 , -1.0 ); */
/*   setScore(pssm, &LETTER_A, 1 , -1.0 ); */
/*   setScore(pssm, &LETTER_A, 2 , -1.0 ); */
/*   setScore(pssm, &LETTER_A, 3 , -1.0 ); */
/*   setScore(pssm, &LETTER_C, 0 , -1.0 ); */
/*   setScore(pssm, &LETTER_C, 1 , -1.0 ); */
/*   setScore(pssm, &LETTER_C, 2 , -1.0 ); */
/*   setScore(pssm, &LETTER_C, 3 , -1.0 ); */
/*   setScore(pssm, &LETTER_G, 0 , -1.0 ); */
/*   setScore(pssm, &LETTER_G, 1 , -1.0 ); */
/*   setScore(pssm, &LETTER_G, 2 , -1.0 ); */
/*   setScore(pssm, &LETTER_G, 3 , -1.0 ); */
/*   setScore(pssm, &LETTER_T, 0 , 0.0 ); */
/*   setScore(pssm, &LETTER_T, 1 , 0.0 ); */
/*   setScore(pssm, &LETTER_T, 2 , 0.0 ); */
/*   setScore(pssm, &LETTER_T, 3 , 0.0 ); */

/*   calcAndSetThresholds(pssm, -0.5); */

/*   printPSSM(pssm);  */

/*   pssm->scores[ */


/*   results = search(esa, pssm); */
/*   writeResultsToFile("esa.res-1", results);  */

/*   exit(0); */
  

  /* LOAD MATRIX */
/*   fprintf(stderr, "Loading matrix...\n"); */
/*   file = fopen("../test/HNF1", "r"); */
/*   fscanf(file, "# %d\n", &i); */
/*   fscanf(file, "> %d %d\n", &order, &length); */
/*   fprintf(stderr, "ORDER(%i) LEN(%i)\n", order, length);   */
  
/*   pssm = initMatrix(0, length, 4); */

/*   for (j = 0; j < alphlen; j++) { */
/*     for (k = 0; k < length; k++) { */
/*       fscanf(file, "%d", &pssm->counts[j * length + k]);   */
/*       pssm->scores[j * length + k] = (pssm->counts[j * length + k] > 0 ? log(((double) pssm->counts[j * length + k])) : -10); */
/*       fprintf(stderr, "%5.2f\t", pssm->scores[j * length + k]);  */
/*     } */
/*     fprintf(stderr, "\n"); */
/*   } */
/*   fclose(file); */

/*   calcAndSetThresholds(pssm, 15.0); */

