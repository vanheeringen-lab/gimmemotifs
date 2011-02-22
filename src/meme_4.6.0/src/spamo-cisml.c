
#include <string.h>

#include "cisml-sax.h"
#include "motif.h"
#include "red-black-tree.h"
#include "spamo-cisml.h"
#include "spamo-output.h"
#include "utils.h"

// datastructure to keep track of the loading of the primary cisml file
typedef struct PRIMARY_LOADER PRIMARY_LOADER_T;
struct PRIMARY_LOADER {
  // the primary motif
  MOTIF_T *motif;
  // the margin or minimum distance from the edge of a sequence that a primary can appear
  int margin;
  // the sequences
  RBTREE_T *sequences;
  // the current sequence for which scores are being read
  SEQUENCE_T *current_sequence;
  // the current best score for the sequence (0 means no hit)
  double current_score;
  // are we currently reading data for the primary motif?
  BOOLEAN_T in_motif;
  // a list of positions that have the best score
  int *hits;
  // the allocated size of the hits array
  int hits_size;
  // the number of hits currently set
  int hit_count;
  // the minimum score allowed for a hit
  double score_threshold;
};

// datastructure to keep track of the loading of the secondary cisml file
typedef struct SECONDARY_LOADER SECONDARY_LOADER_T;
struct SECONDARY_LOADER {
  // should sequence matches be dumped?
  BOOLEAN_T output_sequences;
  // the output directory
  char *output_directory;
  // the area on either side of the primary hit that a secondary can appear
  int margin;
  // the bin size used for graphing and calculating p-values
  int bin;
  // the largest pvalue considered significant
  double significant_pvalue;
  // the maximum distance from the primary that pvalues are calculated
  int test_max;
  // the total number of tests, used for Bonferroni correction
  int total_tests;
  // the current motif database being loaded
  int db_id;
  // the trimmed, centered, filtered sequences
  RBTREE_T *sequences;
  // the current sequence for which scores are being read
  SEQUENCE_T *current_sequence;
  // the primary motif
  MOTIF_T *primary_motif;
  // cached value for the leftmost position of the primary motif
  // in the original sequence
  int primary_lpos;
  // cached value for the rightmost position of the primary motif
  // in the original sequence
  int primary_rpos;
  // the secondary motifs
  RBTREE_T *secondary_motifs;
  // the current secondary motif for which scores are being read
  SECONDARY_MOTIF_T *secondary_motif;
  // a best match for each sequence
  int *secondary_matches;
  // the best score for the current sequence
  double secondary_score;
  // a list of positions that have the best score
  int *hits;
  // the allocated size of the list of hits
  int hits_size;
  // the current count of positions that have the best score
  int hit_count;
  // the minimum score allowed for a hit
  double score_threshold;
};

//parsing functions for reading the CISML file

/**************************************************************************
 * Callback invoked when matching an opening pattern tag in the CISML file 
 * for the primary motif. Checks to see if the pattern refers to the primary
 * motif and sets the flag in_motif to indicate this.
 **************************************************************************/
void motif_primary(void *ctx, char *accession, char *name, char *db, char *lsId, double *pvalue, double *score) {
  PRIMARY_LOADER_T *loader = (PRIMARY_LOADER_T*)ctx;
  loader->in_motif = (strcmp(get_motif_id(loader->motif), accession) == 0);
}

/**************************************************************************
 * Callback invoked when matching an opening scanned_sequence tag in the 
 * CISML file for the primary motif. Checks if the sequence is one we are 
 * scoring and if so records it as the current sequence as well as clearing
 * the hits list.
 **************************************************************************/
void sequence_primary(void *ctx, char *accession, char *name, char *db, char *lsId, double *score, double *pvalue, long *length) {
  PRIMARY_LOADER_T *loader = (PRIMARY_LOADER_T*)ctx;
  if (!(loader->in_motif)) {
    loader->current_sequence = NULL;
  } else {
    RBNODE_T *node = rbtree_lookup(loader->sequences, name, FALSE, NULL);
    if (node) {
      loader->current_sequence = rbtree_value(node);
      if (loader->current_sequence->primary_match) die("Already seen this sequence! We can't process this information "
          "because the scoring information from the previous sighting has already been discarded.\n");
      loader->current_score = 0; // reset the current score
      loader->hit_count = 0; //reset the hit count
    } else {
      loader->current_sequence = NULL;
    }
  }
}

/**************************************************************************
 * Callback invoked when matching a matched_element tag in the CISML file
 * for the primary motif. If we are recording scores for this motif and
 * sequence then it:
 * 1) Checks that a score was supplied
 * 2) Checks that the start and stop are correctly spaced for the expected 
 *    motif.
 * 3) Checks that the hit does not overlap the margin on each end of the
 *    sequence.
 * 4) If we don't have a best score, or this score is better: 
 *      - clear the list of best hits and add this one.
 * 5) Alternately if this score is equal to the existing one:
 *      - add the hit to the list of best hits.
 **************************************************************************/
void match_primary(void *ctx, long start, long stop, double *score, double *pvalue, char *clusterId) {
  PRIMARY_LOADER_T *loader = (PRIMARY_LOADER_T*)ctx;
  int lpos, rpos, rc;
  //check we're actually loading data
  if (loader->current_sequence == NULL) return;
  //check that this match is worth investigating further
  if (score == NULL) return;
  if (start <= 0 || stop <= 0) {
    die("Expected start and stop fields in cisml to be 1 or larger.\n");
  }
  if (start < stop) {
    lpos = start;
    rpos = stop;
    rc = FALSE;
  } else {
    lpos = stop;
    rpos = start;
    rc = TRUE;
  }
  //check that gap makes sense
  if ((rpos - lpos + 1) != get_motif_length(loader->motif)) {
    die("Motif %s has length %d but a match in a CISML file had a start of %ld and stop of %ld which evaluates to a length of %d\n", 
        get_motif_id(loader->motif), get_motif_length(loader->motif), start, stop, (rpos - lpos + 1) );
  }
  //check left margin
  // For example if we had a margin of 1 then the primary motif must start at 
  // 2 or larger which would allow a secondary motif of length 1 to fit at 
  // position 1
  if (lpos <= loader->margin) return;
  //check right margin
  // For example if we had a sequence of length 5 and a margin of 1 then the 
  // primary motif must finish at 4 or smaller which would allow a secondary 
  // motif of length 1 to fit at position 5
  if (rpos > (loader->current_sequence->length - loader->margin)) return;
  //now see if our existing best match is worse than this one
  if (loader->hit_count == 0 || *score > loader->current_score) {
    loader->current_score = *score;
    loader->hit_count = 1;
    loader->hits[0] = (rc ? -lpos : lpos);
  } else if (*score == loader->current_score) {
    if (loader->hit_count >= loader->hits_size) {
      loader->hits_size = loader->hit_count + 10;
      loader->hits = mm_realloc(loader->hits, sizeof(int) * loader->hits_size);
    }
    loader->hits[loader->hit_count++] = (rc ? -lpos : lpos);
  }
}

/**************************************************************************
 * Callback invoked when matching a scanned_sequence closing tag in the
 * CISML file for the primary motif. If there are no hits for the sequence 
 * then set the primary match to zero. If there is one hit then record it as 
 * the primary match. If there are many hits then select one at random to
 * record as the primary match.
 **************************************************************************/
void sequence_end_primary(void *ctx) {
  PRIMARY_LOADER_T *loader = (PRIMARY_LOADER_T*)ctx;
  int match;
  //check we're actually loading data
  if (loader->current_sequence == NULL) return;
  if (loader->hit_count == 0 || loader->current_score < loader->score_threshold) {
    match = 0;
  } else if (loader->hit_count == 1) {
    match = loader->hits[0];
  } else {
    //we had multiple matches at the best score so pick one at random
    match = loader->hits[(int)(loader->hit_count * ((double)rand()/RAND_MAX))];
  }
  loader->current_sequence->primary_match = match;
}

/**************************************************************************
 * Reads a CISML file containing scores for the primary motif for each
 * of the sequences and records the best match above the score threshold
 * for each sequence.
 **************************************************************************/
void load_spamo_primary(const char *file, int margin, double score_threshold, RBTREE_T *sequences, MOTIF_T *motif) {
  CISML_CALLBACKS_T callbacks;
  PRIMARY_LOADER_T data;
  memset(&callbacks, 0, sizeof(CISML_CALLBACKS_T));
  callbacks.start_pattern = motif_primary;
  callbacks.start_scanned_sequence = sequence_primary;
  callbacks.start_matched_element = match_primary;
  callbacks.end_scanned_sequence = sequence_end_primary;
  data.motif = motif;
  data.margin = margin;
  data.sequences = sequences;
  data.current_sequence = NULL;
  data.in_motif = FALSE;
  data.hits_size = 10;
  data.hits = mm_malloc(sizeof(int) * data.hits_size);
  data.hit_count = 0;
  data.score_threshold = score_threshold;
  parse_cisml(&callbacks, &data, file);
  free(data.hits);
}



/**************************************************************************
 * Callback invoked when matching an opening pattern tag for a CISML file 
 * of a secondary motif database. It checks that the motif should be scored,
 * clears out the list of sequence matches and stores the current motif.
 **************************************************************************/
void motif_secondary(void *ctx, char *accession, char *name, char *db, char *lsId, double *pvalue, double *score) {
  SECONDARY_LOADER_T *loader = (SECONDARY_LOADER_T*)ctx;
  SECONDARY_KEY_T key;
  RBNODE_T *node;
  int i, seq_count;
  key.db_id = loader->db_id;
  key.motif_id = accession;
  node = rbtree_lookup(loader->secondary_motifs, &key, FALSE, NULL);
  if (node != NULL) {
    loader->secondary_motif = (SECONDARY_MOTIF_T*)rbtree_value(node);
    if (!(loader->secondary_motif->loaded)) {
      seq_count = rbtree_size(loader->sequences);
      for (i = 0; i < seq_count; ++i) loader->secondary_matches[i] = 0;
    } else {
      die("Already seen CISML data for this motif!");
    }
  } else {
    loader->secondary_motif = NULL;
  }
}

/**************************************************************************
 * Callback invoked when matching an opening scanned_sequence tag for a
 * CISML file of a secondary motif database. It calcualtes and caches the
 * left and right bounds of the primary motif and stores the current 
 * sequence.
 **************************************************************************/
void sequence_secondary(void *ctx, char *accession, char *name, char *db, char *lsId, double *score, double *pvalue, long *length) {
  SECONDARY_LOADER_T *loader = (SECONDARY_LOADER_T*)ctx;
  RBNODE_T *node;
  int pmatch;
  if (loader->secondary_motif == NULL) return;
  node = rbtree_lookup(loader->sequences, accession, FALSE, NULL);
  if (node != NULL) {
    loader->current_sequence = (SEQUENCE_T*)rbtree_value(node);
    pmatch = loader->current_sequence->primary_match;
    loader->primary_lpos = (pmatch < 0 ? -pmatch : pmatch);
    loader->primary_rpos = loader->primary_lpos + get_motif_length(loader->primary_motif) - 1;
    if (loader->secondary_matches[loader->current_sequence->index] != 0) {
      die("Already seen this sequence!");
    }
    loader->secondary_score = 0;
    loader->hit_count = 0;
  } else {
    loader->current_sequence = NULL;
  }
}

/**************************************************************************
 * Callback invoked when matching an opening matched_element tag for a 
 * CISML file of a secondary motif database. A hit must pass the checks:
 * 1) The current match is for a sequence/motif that we're interested in.
 * 2) A score is supplied.
 * 3) The score supplied is better or equal to the existing best score.
 * 4) Consistant with CISML format so start and stop are larger than 0.
 * 5) The distance between the start and stop matches the motif.
 * 6) It fits within the margin region around the primary motif
 * 7) It does not overlap the primary motif. 
 * Provided all those checks pass then the hit is calculated relative to
 * the start of the matched region. If the score is equal to the current
 * best then the relative hit position is added to the list of best hits,
 * otherwise the list is cleared, the best score is updated and the hit
 * is added to the previously empty list.
 **************************************************************************/
void match_secondary(void *ctx, long start, long stop, double *score, double *pvalue, char *clusterId) {
  SECONDARY_LOADER_T *loader = (SECONDARY_LOADER_T*)ctx;
  int lpos, rpos, rc, relative_position, match;
  //check if we're loading this match
  if (loader->current_sequence == NULL) return;
  //check if this match has enough information to be considered
  if (score == NULL) return;
  //check to see if the existing match is better
  if (loader->hit_count > 0 && loader->secondary_score > *score) return;
  //convert the coordinates of the match into easier to use ones
  if (start <= 0 || stop <= 0) {
    die("Expected start and stop fields in cisml to be 1 or larger.\n");
  }
  if (start < stop) {
    lpos = start;
    rpos = stop;
    rc = FALSE;
  } else {
    lpos = stop;
    rpos = start;
    rc = TRUE;
  }
  //check that gap makes sense
  if ((rpos - lpos + 1) != get_motif_length(loader->secondary_motif->motif)) {
    die("Motif %s has length %d but a match in a CISML file had a start of %ld and stop of %ld which evaluates to a length of %d\n", 
        get_motif_id(loader->secondary_motif->motif), get_motif_length(loader->secondary_motif->motif), start, stop, (rpos - lpos + 1) );
  }
  //check for overlap with the primary match
  //and that the secondary motif fits within the margin
  if (rpos < loader->primary_lpos) { // left side (upstream)
    if ((loader->primary_lpos - lpos) > loader->margin) return;//outside margin    
  } else if (lpos > loader->primary_rpos) { // right side (downstream)
    if ((rpos - loader->primary_rpos) > loader->margin) return;//outside margin
  } else {
    return;//overlap
  }
  //match seems valid and better than anything we've seen previous so update
  //note that stored position is relative to the start of the margin, as if
  //this was scored on a trimmed sequence indexing from 1 
  //this has the advantage that we only need the width of the primary
  //motif and the size of the margin to calculate the offset
  relative_position = lpos - (loader->primary_lpos - loader->margin) + 1;
  //now make the scale pos/neg dependent on if the match is with a
  //reverse complement
  match = (rc ? -relative_position : relative_position);
  if (loader->hit_count == 0 || loader->secondary_score > *score) {
    loader->secondary_score = *score;
    loader->hit_count = 1;
    loader->hits[0] = match;
  } else if (loader->secondary_score == *score) {
    if (loader->hit_count >= loader->hits_size) {
      loader->hits_size = loader->hit_count + 10;
      loader->hits = mm_realloc(loader->hits, sizeof(int) * loader->hits_size);
    }
    loader->hits[loader->hit_count++] = match;
  }
}

/**************************************************************************
 * Callback invoked when matching a closing scanned_sequence tag for a CISML
 * file of a secondary motif database. If no hits were found then the 
 * secondary match is set to zero. If one hit was found then the secondary
 * match is set to the hit. If multiple best hits were found then one is 
 * picked at random to become the secondary match.
 **************************************************************************/
void sequence_end_secondary(void *ctx) {
  SECONDARY_LOADER_T *loader = (SECONDARY_LOADER_T*)ctx;
  int match;
  //check we're actually loading data
  if (loader->current_sequence == NULL) return;
  if (loader->hit_count == 0 || loader->secondary_score < loader->score_threshold) {
    match = 0;
  } else if (loader->hit_count == 1) {
    match = loader->hits[0];
  } else {
    //we had multiple matches at the best score so pick one at random
    match = loader->hits[(int)(loader->hit_count * ((double)rand()/RAND_MAX))];
  }
  loader->secondary_matches[loader->current_sequence->index] = match;
}

/**************************************************************************
 * Callback invoked when matching a closing pattern tag for a CISML file of
 * a secondary motif database. Processes the matches found (see 
 * spamo-matches.c for more details) and if the motif had a significant 
 * spacing it optionally dumps the sequence matches to file.
 **************************************************************************/
void motif_end_secondary(void *ctx) {
  SECONDARY_LOADER_T *loader = (SECONDARY_LOADER_T*)ctx;
  if (loader->secondary_motif == NULL) return;
  // process matches
  process_matches(loader->margin, loader->bin, loader->significant_pvalue, loader->test_max, loader->total_tests,
      loader->primary_motif, loader->sequences, loader->secondary_motif, loader->secondary_matches);
  // dump sequences (if requested)
  if (loader->output_sequences && loader->secondary_motif->sig_count > 0) {
    output_sequence_matches(loader->output_directory, loader->margin, loader->sequences, 
        loader->primary_motif, loader->secondary_motif, loader->secondary_matches);
  }
}

/**************************************************************************
 * Reads a CISML file containing the scores for secondary motif database
 * and tallys the spacings of the best matches of primary motifs to secondary
 * motifs.
 **************************************************************************/
void load_spamo_secondary(
    const char *file, int margin, double score_threshold, int bin, double sigthresh, 
    int test_max, int total_tests, BOOLEAN_T output_sequences, char *output_directory, 
    RBTREE_T *sequences, MOTIF_T *primary_motif, RBTREE_T *secondary_motifs, int db_id
  ) {
  CISML_CALLBACKS_T callbacks;
  SECONDARY_LOADER_T data;
  memset(&callbacks, 0, sizeof(CISML_CALLBACKS_T));
  callbacks.start_pattern = motif_secondary;
  callbacks.start_scanned_sequence = sequence_secondary;
  callbacks.start_matched_element = match_secondary;
  callbacks.end_scanned_sequence = sequence_end_secondary;
  callbacks.end_pattern = motif_end_secondary;
  data.output_sequences = output_sequences;
  data.output_directory = output_directory;
  data.margin = margin;
  data.bin = bin;
  data.significant_pvalue = sigthresh;
  data.test_max = test_max;
  data.total_tests = total_tests;
  data.sequences = sequences;
  data.primary_motif = primary_motif;
  data.primary_lpos = 0;
  data.primary_rpos = 0;
  data.secondary_motifs = secondary_motifs;
  data.secondary_motif = NULL;
  data.secondary_matches = mm_malloc(sizeof(int) * rbtree_size(sequences));
  data.hits_size = 10;
  data.hits = mm_malloc(sizeof(int) * data.hits_size);
  data.hit_count = 0;
  data.score_threshold = score_threshold;
  data.db_id = db_id;
  parse_cisml(&callbacks, &data, file);
  free(data.secondary_matches);
  free(data.hits);
}
