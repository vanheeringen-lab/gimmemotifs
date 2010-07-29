#include <stdlib.h>
#include <stdio.h>
#include "sesa.h"


int main() {
  SESA sesa;
  PSSM pssm;
  unsigned int i;

  /* 
   * NOTE: It can only load one matrix at the monment! 
   * Should be changed to return a PSSMSet 
   */
  pssm = load_log_matrices("small-wm", 4);
  print_pssm(pssm); 
  
  sesa = EESA2SESA(load_fasta_DNA("sequences", 0));
  calcAndSetThresholds(pssm, -0.5);


  /* 
   * This is stupid when you are doing many searches. Use SESA
   * instead.
   */
  struct HitTable hits = SESA_search(sesa, pssm);
  
  for (i = 0; i < hits->nScores; i++) {
    int p = hits->pScores[i].position;
    printf("HIT: %3i is in seq %3i at pos %5i with score %lf\n", i, sesa->seq[p], sesa->pos[p], hits->pScores[i].score);
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
