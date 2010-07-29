#include <stdlib.h>
#include <stdio.h>
#include "eesa.h"


int main() {
  EESA eesa;
  PSSM pssm;
  Hits hits;
  unsigned int i;

  /* 
   * NOTE: It can only load one matrix at the monment! 
   * Should be changed to return a PSSMSet 
   */
  pssm = load_log_matrices("small-wm", 4);
  print_pssm(pssm); 
  
  eesa = load_fasta_DNA("sequences", 0);
  calcAndSetThresholds(pssm, -0.5);


  /* 
   * This is stupid when you are doing many searches. Use SESA
   * instead.
   */
  hits = single_search(eesa, pssm);
  
  for (i = 0; i < hits->nScores; i++) {
    printf("HIT: %3i is in seq %3i at pos %5i with score %lf\n", i, hits->pScores[i].seq, hits->pScores[i].pos, hits->pScores[i].score);
  }
  
  exit(0);
}


/**
 *
 * MATRIX Format:
 *
 *    # NUMBER_OF MATRICES
 *    > ORDER LENGTH
 *    ..VALUES...
 *
 *
 * EXAMPLE:
 *
 *    # 1
 *    > 0 3
 *    1.0    2.1   1.0   
 *    -2.3   4.0   0.3
 *    -0.5   0.0   -0.3
 *    0.0    -10   -0.1
 *
 *
 *
 * Count matrices uses the same format. Use: load_count_matrices
 */
