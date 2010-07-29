/*---------------------------------------------------------------*
 * PSSM Searcher                                                 *
 *                                                               *
 * By: Jes Frellsen, Ida Moltke and Martin Thiim                 *
 * University of Copenhagen, spring 2006                         *
 *---------------------------------------------------------------*/

#include "memcheck.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "shared.h"
#include "build.h"
#include "hittable.h"
#include "search.h"
#include "pssm.h"
#include "error.h"

int compare_pEntries(const void *pEntry1, const void *pEntry2){
  return ((struct HitEntry *)pEntry1)->position - ((struct HitEntry *)pEntry2)->position;
}

void writeResultsToFile(char *file, struct HitTable *pScoreData);

int main() {
  int k,n;
  char *text;
  FILE *f;
  ESA esa;
  PSSM pssm;  
  
  // The pssm array is initialised
  pssm = initMatrix(0, 12, 4);
  #include "data/matrixinit.h"
  calcAndSetThresholds(pssm, 12.0);

  // Open file
  f = fopen("../test/testseq", "r");
  if(f == 0) {
    printf("File not found.\n");
    return 0;
  }
  // Obtain file size.
  fseek (f, 0 , SEEK_END);
  n = ftell (f);    
  rewind (f);
  printf("File size: %d\n", n);

  text = malloc(sizeof(char)*(n));
  printf("TEXT\n");
  k = fread(text, 1, n, f);
  printf("Read size: %d\n", k);
  fclose(f);

  // Build the ESA
  esa = build_ESA(text, n, "ACTG", "|", 0);
  
  /*esa = read_ESA_from_file("/tmp/esa.dat", &hej, &kk);
  
  printf("Extra Amount %d", kk);
  printf("Extra data %s", hej);
  
  if(esa == NULL)
  {
  	printf("Build error. Aborting.\n");
  	return 0;
  }*/
  
  printf("Finished build\n");
  
  
  // Testing the different search algorithms
  struct HitTable *pScoreData;	
  
  
  // The suffixes are scored ESA search
  pScoreData = search(esa, pssm);	
  writeResultsToFile("/tmp/esa.res", pScoreData);
  release_hits(pScoreData);
 
  // The suffixes are scored using the naive algorithm
  pScoreData = searchNaively(getStr(esa), getSize(esa), pssm);
  writeResultsToFile("/tmp/naiv.res", pScoreData);
  release_hits(pScoreData);

  // The suffixes are scored using the naive algorithm with LA
  pScoreData = searchNaivelyWithLA(getStr(esa), getSize(esa), pssm);
  writeResultsToFile("/tmp/naivLA.res", pScoreData);
  release_hits(pScoreData);

  writeESAToFile("/tmp/esa.dat", esa, "HEJSA", 6);
  
  free_ESA(esa);
  releaseMatrix(pssm);

  free(text);

  return 0;
}


void writeResultsToFile(char *file, struct HitTable *pScoreData) {
  int i, j, k;
  struct HitEntry *pEntry;
  FILE *f;

  i = 0;
  pEntry = pScoreData->pScores;
  k = pScoreData->nScores;
  printf("\nResult from using ESA search:\n");
  printf("Matrix %d has %d entries\n", i, k);
 
  if(k<100){
    qsort((void*)pEntry, k, sizeof(struct HitEntry), compare_pEntries);
    for(j = 0; j < k; j++){
      printf("Entry %d has pos %d ", j, pEntry[j].position);
      printf("and score %f\n", (float)pEntry[j].score);
    }
  }
  else{
    f = fopen(file, "wt");
    
    if(!f) {
      printf("Could not write file %s\n", file);
    }
    else {
      qsort((void*)pEntry,k,sizeof(struct HitEntry),compare_pEntries);
      for(j = 0; j < k; j++){
	fprintf(f,"pos %d ", pEntry[j].position);
	fprintf(f,"score %f\n", (float)pEntry[j].score);
      }
      fclose(f);
      printf("Result can be seen in %s\n", file);
    }
  }
}

