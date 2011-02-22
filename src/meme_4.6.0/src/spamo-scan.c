
#include "matrix.h"
#include "alphabet.h"
#include "assert.h"
#include "pssm.h"
#include "spamo-matches.h"
#include "spamo-scan.h"

#include <stdlib.h>

/**************************************************************************
 * Scans specified regions of a sequence with a PSSM and its reverse 
 * complement and returns a single best match if it found a position that
 * scores above the specified score threshold. If there were multiple 
 * equally good matches then one is picked at random.
 **************************************************************************/
static inline int best_in_region(
  PSSM_T *pssm,             // the motif score matrix
  PSSM_T *pssm_rc,          // the revese complement motif score matrix
  char *seq,                // the sequence to score
  int *regions,             // the regions to scan
  double score_threshold,   // the minimum score considered a hit
  int *hits,                // the pre-allocated hits cache should be large enough to store a hit for every position
  int hits_size             // the size of the pre-allocatted hits cache.
) {
  int region_index, position_index, motif_offset, motif_length, region_end, alphabet_size, alph_index, success, hit_count;
  char c, *alphabet;
  double lo, lo_rc, best_lo;
  MATRIX_T *pssm_matrix, *pssm_rc_matrix;
  //basic checks on parameters
  assert(pssm != NULL);
  assert(pssm_rc != NULL);
  assert(get_pssm_w(pssm) == get_pssm_w(pssm_rc));
  assert(seq != NULL);
  assert(regions != NULL);
  //get the alphabet
  //oh how detest having to use globals... see alphabet.h if this comment doesn't make any sense
  alphabet = get_alphabet(FALSE);
  alphabet_size = get_alph_size(ALPH_SIZE);
  //initilize the best position: 0 means no match; +ve is match on 5` strand; -ve is match on 3` strand
  // to make this work an origin of 1 is used.
  best_lo = 0;
  hit_count = 0;
  // get the motif length
  motif_length = get_pssm_w(pssm);
  // get the pssm matrixes
  pssm_matrix = pssm->matrix;
  pssm_rc_matrix = pssm_rc->matrix;
  //scan regions, check for overlaps and those not in sorted order
  for (region_index = 0; regions[region_index] != -1; region_index += 2) {
    //error checks, so nothing is too suprising
    assert(regions[region_index] >= 0);
    assert(region_index == 0 || regions[region_index-1] < regions[region_index]);
    //scan a region
    region_end = regions[region_index+1] - motif_length + 1;
    for (position_index = regions[region_index]; position_index <= region_end; ++position_index) {
      //scan a position
      success = TRUE; // assume success
      lo = 0; //reset sums to zero
      lo_rc = 0;
      for (motif_offset = 0; motif_offset < motif_length; ++motif_offset) {
        c = seq[position_index+motif_offset];
        alph_index = alphabet_index(c, alphabet);

        // Check for gaps and ambiguity codes at this site
        if(c == '-' || c == '.' || alph_index >= alphabet_size) {
          success = FALSE;
          break;
        }
        //note these scores are scaled
        lo += get_matrix_cell(motif_offset, alph_index, pssm_matrix);
        lo_rc += get_matrix_cell(motif_offset, alph_index, pssm_rc_matrix); 
      }
      if (success) {
        //now revert the scaled lo scores to unscaled
        lo = get_unscaled_pssm_score(lo, pssm);
        lo_rc = get_unscaled_pssm_score(lo_rc, pssm);

        if (lo > lo_rc) {
          if (lo > best_lo) {
            best_lo = lo;
            hit_count = 1;
            hits[0] = position_index + 1;
          } else if (lo == best_lo) {
            hits[hit_count++] = position_index + 1;
          }
        } else if (lo == lo_rc) {
          if (lo > best_lo) {
            best_lo = lo;
            hit_count = 0;
            hits[hit_count++] = position_index + 1;
            hits[hit_count++] = -(position_index + 1);
          } else if (lo == best_lo) {
            hits[hit_count++] = position_index + 1;
            hits[hit_count++] = -(position_index + 1);
          }
        } else {
          if (lo_rc > best_lo) {
            best_lo = lo_rc;
            hit_count = 1;
            hits[0] = -(position_index + 1);
          } else if (lo_rc == best_lo) {
            hits[hit_count++] = -(position_index + 1);
          }
        }
      }
    }
  }
  if (best_lo >= score_threshold) {
    return hits[(int)(hit_count * ((double)rand()/RAND_MAX))];
  }
  return 0;
}

/**************************************************************************
 * Scan all the sequences with the primary motif and store the best match
 * for each in the sequence data structure.
 **************************************************************************/
void scan_spamo_primary(
  int margin, 
  double score_threshold, 
  ARRAY_T *background, 
  MOTIF_T *motif, 
  RBTREE_T *sequences
) {
  MOTIF_T *motif_rc;
  PSSM_T *pssm, *pssm_rc;
  RBNODE_T *node;
  SEQUENCE_T *sequence;
  int regions[3];
  int *hits, hits_size;
  if ((node = rbtree_first(sequences)) == NULL) return;
  sequence = (SEQUENCE_T*)rbtree_value(node);
  hits_size = (sequence->length - 2 * margin) * 2;
  hits = mm_malloc(sizeof(int) * hits_size);

  regions[0] = margin; // first valid index to scan
  regions[2] = -1; //terminate list with negative 1
  //prepare a reverse complement of the motif
  motif_rc = mm_malloc(sizeof(MOTIF_T));
  copy_motif(motif, motif_rc);
  reverse_complement_motif(motif_rc);
  //convert the motif and the reverse complement into PSSMs
  pssm =    build_motif_pssm(motif,    background, background, NULL, 0, PSSM_RANGE, 0, FALSE);
  pssm_rc = build_motif_pssm(motif_rc, background, background, NULL, 0, PSSM_RANGE, 0, FALSE);
  //scan each sequence
  for (node = rbtree_first(sequences); node != NULL; node = rbtree_next(node)) {
    sequence = (SEQUENCE_T*)rbtree_value(node);
    if (((sequence->length - 2 * margin) * 2) > hits_size) {
      hits_size = (sequence->length - 2 * margin) * 2;
      hits = mm_realloc(hits, sizeof(int) * hits_size);
    }
    regions[1] = sequence->length - margin - 1; // last valid index to include in scan
    sequence->primary_match = best_in_region(pssm, pssm_rc, sequence->data, regions, score_threshold, hits, hits_size);
  }
  //clean up the PSSMs
  free_pssm(pssm);
  free_pssm(pssm_rc);
  //clean up the motif
  free_motif(motif_rc);
  free(motif_rc);
  //clean up the hits
  free(hits);
}

/**************************************************************************
 * Scan all the sequences with the secondary motif and store the best match
 * for each in matches list.
 **************************************************************************/
void scan_spamo_secondary(
  int margin, 
  double score_threshold,
  ARRAY_T *background, 
  MOTIF_T *motif,
  RBTREE_T *sequences, 
  int *matches,
  int *hits,
  int hits_size
) {
  MOTIF_T *motif_rc;
  PSSM_T *pssm, *pssm_rc;
  RBNODE_T *node;
  SEQUENCE_T *sequence;
  int regions[5];
  if ((node = rbtree_first(sequences)) == NULL) return;
  sequence = (SEQUENCE_T*)rbtree_value(node);
  // set up the two regions flanking the primary motif
  regions[0] = 0; // first valid index to scan
  regions[1] = margin - 1;
  regions[2] = sequence->length - margin;
  regions[3] = sequence->length - 1;
  regions[4] = -1; // terminate list with negative 1
  // prepare a reverse complement of the motif
  motif_rc = mm_malloc(sizeof(MOTIF_T));
  copy_motif(motif, motif_rc);
  reverse_complement_motif(motif_rc);
  // convert the motif and the reverse complement into PSSMs
  pssm =    build_motif_pssm(motif,    background, background, NULL, 0, PSSM_RANGE, 0, FALSE);
  pssm_rc = build_motif_pssm(motif_rc, background, background, NULL, 0, PSSM_RANGE, 0, FALSE);
  // scan all the sequences
  for (node = rbtree_first(sequences); node != NULL; node = rbtree_next(node)) {
    sequence = (SEQUENCE_T*)rbtree_value(node);
    matches[sequence->index] = best_in_region(pssm, pssm_rc, sequence->data, regions, score_threshold, hits, hits_size);
  }
  //clean up the PSSMs
  free_pssm(pssm);
  free_pssm(pssm_rc);
  //clean up the motif
  free_motif(motif_rc);
  free(motif_rc);
}
