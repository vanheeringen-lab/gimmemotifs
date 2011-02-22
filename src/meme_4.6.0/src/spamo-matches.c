
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "spamo-matches.h"

/**************************************************************************
 * Compares two keys of a secondary motif for equality.
 * First compares based on the motif database and second compares
 * based on the name of the motif.
 **************************************************************************/
int secondary_key_compare(void *p1, void *p2) {
  SECONDARY_KEY_T *k1, *k2;
  k1 = (SECONDARY_KEY_T*)p1;
  k2 = (SECONDARY_KEY_T*)p2;
  if (k1->db_id == k2->db_id) {
    return strcmp(k1->motif_id, k2->motif_id);
  } else if (k1->db_id < k2->db_id) {
    return -1;
  } else {
    return 1;
  }
}

/**************************************************************************
 * Copies a key of a secondary motif
 **************************************************************************/
void* secondary_key_copy(void *p) {
  SECONDARY_KEY_T *source, *dest;
  source = (SECONDARY_KEY_T*)p;
  dest = mm_malloc(sizeof(SECONDARY_KEY_T));
  dest->db_id = source->db_id;
  dest->motif_id = source->motif_id;
  return dest;
}

/**************************************************************************
 * Creates a structure to hold the sequence database information.
 * All strings are copied.
 **************************************************************************/
SEQUENCE_DB_T* create_sequence_db(char *file) {
  SEQUENCE_DB_T* db;
  struct stat stbuf;
  char *c;
  int len, i, stat_result;
  //allocate memory
  db = (SEQUENCE_DB_T*)mm_malloc(sizeof(MOTIF_DB_T));
  //copy the source
  copy_string(&(db->source), file);
  //copy the name, assumes unix style file separator
  copy_string(&(db->name), ((c = strrchr(file, '/')) ? c+1 : file));
  //clean up the name
  len = strlen(db->name);
  //if it ends with .fasta or .fa then truncate the string
  if (len > 6 && strcmp(db->name+(len - 6), ".fasta") == 0) db->name[len-6] = '\0';
  else if (len > 3 && strcmp(db->name+(len - 3), ".fa") == 0) db->name[len-3] = '\0';
  //replace underscores
  for (c = db->name; *c != '\0'; c++) if (*c == '_') *c = ' ';
  //stat the file and get the last modified date
  if (stat(file, &stbuf) < 0) {
    //an error occured
    die("Failed to stat file \"%s\", error given as: %s\n", file, strerror(errno));
  }
  db->last_mod = stbuf.st_mtime;
  //initilize defaults
  db->loaded = 0;
  db->excluded_nomatch = 0;
  db->excluded_similar = 0;
  return db;
}

/**************************************************************************
 * Destroy the sequence db
 **************************************************************************/
void destroy_sequence_db(SEQUENCE_DB_T *db) {
  free(db->source);
  free(db->name);
  free(db);
}

/**************************************************************************
 * Creates a structure to hold the motif database information.
 * All strings are copied.
 **************************************************************************/
MOTIF_DB_T* create_motif_db(int id, char *file, char *cisml) {
  MOTIF_DB_T* db;
  struct stat stbuf;
  char *c;
  int len, i, stat_result;
  //allocate memory
  db = (MOTIF_DB_T*)mm_malloc(sizeof(MOTIF_DB_T));
  //copy the id
  db->id = id;
  //copy the source
  copy_string(&(db->source), file);
  //copy the cisml
  copy_string(&(db->cisml), cisml);
  //copy the name, assumes unix style file separator
  copy_string(&(db->name), ((c = strrchr(file, '/')) ? c+1 : file));
  //clean up the name
  len = strlen(db->name);
  //if it ends with .meme then truncate the string
  if (len > 5 && strcmp(db->name+(len - 5), ".meme") == 0) db->name[len-5] = '\0';
  //replace underscores
  for (c = db->name; *c != '\0'; c++) if (*c == '_') *c = ' ';
  //stat the file and get the last modified date
  if (stat(file, &stbuf) < 0) {
    //an error occured
    die("Failed to stat file \"%s\", error given as: %s\n", file, strerror(errno));
  }
  db->last_mod = stbuf.st_mtime;
  //initilize defaults
  db->loaded = 0;
  db->excluded = 0;
  return db;
}

/**************************************************************************
 * Destroy the motif db.
 **************************************************************************/
void destroy_motif_db(void *secondary_db) {
  MOTIF_DB_T *db;

  db = (MOTIF_DB_T*)secondary_db;
  free(db->source);
  free(db->name);
  free(db->cisml);
  free(db);
}

/**************************************************************************
 * Destroy the sequence (including the name attribute)
 * Note that the data attribute is destroyed seperately
 **************************************************************************/
void destroy_sequence(void *sequence) {
  SEQUENCE_T *seq;

  seq = (SEQUENCE_T*)sequence;
  free(seq->name);
  free(seq);
}

/**************************************************************************
 * Used by create_secondary_motif to simplify the initilization of the
 * spacings array for each of the 4 spacing types.
 **************************************************************************/
static void init_spacings(SPACING_T *frequencies, int bin_count) {
  int i;
  frequencies->bins = (int*)mm_malloc(sizeof(int) * bin_count);
  frequencies->pvalues = (double*)mm_malloc(sizeof(double) * bin_count);
  for (i = 0; i < bin_count; ++i) {
    frequencies->bins[i] = 0;
    frequencies->pvalues[i] = 1;
  }
}

/**************************************************************************
 * Used by destroy_secondary_motif to deallocate memory for the spacings
 * array of the 4 spacing types.
 **************************************************************************/
static void destroy_spacings(SPACING_T *frequencies) {
  free(frequencies->bins);
  free(frequencies->pvalues);
}

/**************************************************************************
 * Create a secondary motif. As the number of sequences is unknown at this
 * point the sequence_matches array is left unallocated. All pvalues are
 * initilized to 1.
 **************************************************************************/
SECONDARY_MOTIF_T* create_secondary_motif(int margin, int bin, 
    MOTIF_DB_T *db, MOTIF_T *motif) {
  int bin_count, i;
  SECONDARY_MOTIF_T *smotif;
  smotif = mm_malloc(sizeof(SECONDARY_MOTIF_T));
  smotif->db = db;
  smotif->motif = motif;
  //set loaded to false
  smotif->loaded = FALSE;
  //calculate the number of bins needed for this motif
  bin_count = (int)((margin - get_motif_trimmed_length(motif) + 1) / bin) + 1;
  //allocate spacings
  for (i = 0; i < 4; ++i) init_spacings((smotif->spacings)+i, bin_count);
  smotif->total_spacings = 0;
  smotif->max_in_one_bin = 0;
  //these will be allocated after we've filled the spacings tables
  //and calculated the most significant spacings
  smotif->sigs = NULL;
  smotif->sig_count = 0;
  smotif->seqs = NULL;
  smotif->seq_count = 0;
  return smotif;
}

/**************************************************************************
 * Destroys a secondary motif. 
 * It takes a void * pointer so it can be used in the collection objects.
 **************************************************************************/
void destroy_secondary_motif(void *p) {
  int i;
  SECONDARY_MOTIF_T *smotif = (SECONDARY_MOTIF_T*)p;
  free_motif(smotif->motif);
  free(smotif->motif);
  for (i = 0; i < 4; ++i) destroy_spacings((smotif->spacings)+i);
  if (smotif->sigs) free(smotif->sigs);
  if (smotif->seqs) free(smotif->seqs);
  free(smotif);
}

/**************************************************************************
 * Create a group with a good secondary motif and a group
 * of redundant secondary motifs.
 **************************************************************************/
GROUPED_MOTIF_T* create_grouped_motif(SECONDARY_MOTIF_T* best) {
  GROUPED_MOTIF_T *group;
  group = (GROUPED_MOTIF_T*)mm_malloc(sizeof(GROUPED_MOTIF_T));
  group->best = best;
  group->others = linklst_create();
  return group;
}

/**************************************************************************
 * Destroys a grouped motif
 * 
 * Relies on the secondary motifs begin destroyed elsewhere.
 **************************************************************************/
void destroy_grouped_motif(void *p) {
  GROUPED_MOTIF_T *group = (GROUPED_MOTIF_T*)p;
  linklst_destroy(group->others);
  free(group);
}

/**************************************************************************
 * Puts counts into the spacing bins.
 **************************************************************************/
void bin_matches(int margin, int bin_size, RBTREE_T *sequences, MOTIF_T *primary_motif, SECONDARY_MOTIF_T *secondary_motif, int *matches) {
  int primary_len, secondary_len, secondary, secondary_pos, primary_rc, secondary_rc, quad, distance, max_distance;
  RBNODE_T *node;
  SECONDARY_MOTIF_T *smotif;
  SEQUENCE_T *sequence;
  SPACING_T *spacing;

  primary_len = get_motif_trimmed_length(primary_motif);

  smotif = secondary_motif;
  secondary_len = get_motif_trimmed_length(smotif->motif);

  // Note that distance counts from zero
  max_distance = margin - secondary_len;

  // for each sequence
  for (node = rbtree_first(sequences); node != NULL; node = rbtree_next(node)) {
    sequence = (SEQUENCE_T*)rbtree_value(node);
    secondary = matches[sequence->index];
    // check for a match
    if (!secondary) continue;
    // convert the encoded form into easier to use form
    primary_rc = sequence->primary_match < 0;
    secondary_rc = secondary < 0;
    secondary_pos = (secondary_rc ? -secondary : secondary);
    // calculate the distance (counts from zero) and side
    if (secondary_pos <= margin) {
      distance = margin - secondary_pos - secondary_len + 1;
      if (primary_rc) {//rotate reference direction
        quad = RIGHT;
      } else {
        quad = LEFT;
      }
    } else {
      distance = secondary_pos - margin - primary_len - 1;
      if (primary_rc) {//rotate reference direction
        quad = LEFT;
      } else {
        quad = RIGHT;
      }
    }
    // check that we're within the acceptable range
    if (distance < 0 || distance > max_distance) {
      die("Secondary motif match not within margin as it should be due to prior checks!");
    }
    // calculate the strand
    if (secondary_rc == primary_rc) {
      quad |= SAME;
    } else {
      quad |= OPPO;
    }
    // add a count to the frequencies
    spacing = smotif->spacings+(quad);
    spacing->bins[(int)(distance / bin_size)] += 1;
    smotif->total_spacings += 1;
  }
}

/**************************************************************************
 * compute a pvalue by summing the tail of a binomial distribution
 **************************************************************************/
double binomial_exact(int succ, int trials, double prob) {
  double out;
  //note: 
  //gamma(x + 1) = x!
  //lgamma(x) = log(gamma(x)) = log((x-1)!)
  out = lgamma(trials + 1) - lgamma(trials - succ + 1) - lgamma(succ + 1) + log(prob) * succ + log(1.0-prob) * (trials - succ) ;

  return exp(out);
}

/**************************************************************************
 * compute a pvalue by summing the tail of a binomial distribution
 **************************************************************************/
double one_minus_binomial_cdf(double probability, int successful_trials, int total_trials) {
  double sum;
  int x;
  for (x = successful_trials, sum = 0; x <= total_trials; ++x) {
    sum += binomial_exact(x, total_trials, probability);
  }
  return sum;
}

/**************************************************************************
 * Compute a pvalue for a particular bin and quadrant. 
 * Store the location of the bin if the pvalue is significant.
 **************************************************************************/
static void compute_spacing_pvalue(int tests, double threshold, int quad, int bin, int test_bin_count, double prob, SECONDARY_MOTIF_T *smotif) {
  double pvalue;
  int counts;
  
  counts = smotif->spacings[quad].bins[bin];

  //keep track of the maximum count for the histogram.
  if (counts > smotif->max_in_one_bin) {
    smotif->max_in_one_bin = counts;
  }

  if (bin < test_bin_count) {
    pvalue = one_minus_binomial_cdf(prob, counts, smotif->total_spacings);

    //bonferoni correction
    pvalue = exp(LOGEV(log(tests), log(pvalue)));

    smotif->spacings[quad].pvalues[bin] = pvalue;

    if (pvalue <= threshold) {
      SIGSPACE_T *sig;
      smotif->sig_count += 1;
      smotif->sigs = mm_realloc(smotif->sigs, sizeof(SIGSPACE_T) * smotif->sig_count);
      sig = smotif->sigs+(smotif->sig_count - 1);
      sig->pvalue = pvalue;
      sig->bin = bin;
      sig->quad = quad;
    }
  } else {
    //don't test this bin, set pvalue to special value to indicate it is untested
    pvalue = 2;
    smotif->spacings[quad].pvalues[bin] = pvalue;
  }
}

/**************************************************************************
 * compare the pvalues for two significant spacings
 **************************************************************************/
int compare_sigs(const void *p1, const void *p2) {
  SIGSPACE_T *s1 = (SIGSPACE_T*)p1;
  SIGSPACE_T *s2 = (SIGSPACE_T*)p2;
  if (s1->pvalue < s2->pvalue) {
    return -1;
  } else if (s1->pvalue == s2->pvalue) {
    return 0;
  } else {
    return 1;
  }
}

/**************************************************************************
 * compute the pvalues for the frequencies of each spacing
 **************************************************************************/
void compute_spacing_pvalues(int margin, int bin_size, int independent_tests, int test_max, double threshold, SECONDARY_MOTIF_T *smotif) {
  int quad_opt_count, quad_bin_count, quad_leftover, total_opt_count, i, j;
  double general_prob, leftover_prob;
  //the number of possible values for spacings in one quadrant
  quad_opt_count = margin - get_motif_trimmed_length(smotif->motif) + 1;
  //the number of bins in one quadrant (excluding a possible leftover bin)
  quad_bin_count = (int)(quad_opt_count / bin_size);
  //the number of spacings that don't fit in the full bins (the number that would go into the leftover bin)
  quad_leftover = quad_opt_count % bin_size;
  //the total number of possible values for spacings
  total_opt_count = quad_opt_count * 4;
  //prior probability of a bin that has bin_size possible spacings that could go into it
  general_prob = (double)bin_size / total_opt_count;
  //prior probability of the final bin that has less than bin_size possible spacings that could go into it
  leftover_prob = (double)quad_leftover / total_opt_count;
  //calculate the significance of each bin
  for (i = 0; i < 4; ++i) {
    for (j = 0; j < quad_bin_count; ++j) {
      compute_spacing_pvalue(independent_tests, threshold, i, j, test_max, general_prob, smotif);
    }
    if (quad_leftover) { //bin only exists if quad_leftover is non-zero
      compute_spacing_pvalue(independent_tests, threshold, i, j, test_max, leftover_prob, smotif);
    }
  }
  //sort the significant finds
  qsort(smotif->sigs, smotif->sig_count, sizeof(SIGSPACE_T), compare_sigs);
}

/**************************************************************************
 * compute the list of ids for the most significant spacing
 **************************************************************************/
void compute_idset(int margin, int bin_size, RBTREE_T *sequences, MOTIF_T *primary_motif, SECONDARY_MOTIF_T *secondary_motif, int *matches) {
  int primary_len, secondary_len, secondary, secondary_pos, primary_rc, secondary_rc, quad, distance;
  RBNODE_T *node;
  SEQUENCE_T *sequence;

  if (secondary_motif->sig_count == 0) return;

  primary_len = get_motif_trimmed_length(primary_motif);
  secondary_len = get_motif_trimmed_length(secondary_motif->motif);

  // for each sequence
  for (node = rbtree_first(sequences); node != NULL; node = rbtree_next(node)) {
    sequence = (SEQUENCE_T*)rbtree_value(node);
    secondary = matches[sequence->index];
    // check for a match
    if (!secondary) continue;
    // convert the encoded form into easier to use form
    primary_rc = sequence->primary_match < 0;
    secondary_rc = secondary < 0;
    secondary_pos = (secondary_rc ? -secondary : secondary);
    // calculate the distance and side
    // note that distance can be zero meaning the primary is next to the secondary
    if (secondary_pos <= margin) {
      distance = margin - secondary_pos - secondary_len + 1;
      quad = LEFT;
    } else {
      distance = secondary_pos - margin - primary_len;
      quad = RIGHT;
    }
    // calculate the strand
    if (secondary_rc == primary_rc) {
      quad |= SAME;
    } else {
      quad |= OPPO;
    }
    // add the sequence id to the set if the bin matches    
    if (quad == secondary_motif->sigs->quad && (distance / bin_size) == secondary_motif->sigs->bin) {
      secondary_motif->seq_count += 1;
      secondary_motif->seqs = (int*)mm_realloc(secondary_motif->seqs, sizeof(int) * secondary_motif->seq_count);
      secondary_motif->seqs[secondary_motif->seq_count-1] = sequence->index;
    }
  }
}

/**************************************************************************
 * Calculate the total number of pvalue calculations that will be done
 * by the program. This number is used to correct the pvalues for multiple
 * tests using a bonferoni correction.
 **************************************************************************/
int calculate_test_count(int margin, int bin, int test_max, RBTREE_T *secondary_motifs) {
  int total_tests, test_bin_count, quad_opt_count, quad_bin_count;
  SECONDARY_MOTIF_T *smotif;
  RBNODE_T *node;

  total_tests = 0;
  for (node = rbtree_first(secondary_motifs); node != NULL; node = rbtree_next(node)) {
    smotif = (SECONDARY_MOTIF_T*)rbtree_value(node);
    //the number of possible values for spacings in one quadrant
    quad_opt_count = margin - get_motif_trimmed_length(smotif->motif) + 1;
    //the number of bins in one quadrant (excluding a possible leftover bin)
    quad_bin_count = (int)(quad_opt_count / bin) + (quad_opt_count % bin ? 1 : 0);
    //add the number of tested bins
    total_tests += (test_max < quad_bin_count ? test_max : quad_bin_count) * 4;
  }
  return total_tests;
}

/**************************************************************************
 * This does most of the calculation steps after the matching
 * positions of the motifs have been derived from either a scan or a cisml
 * file. The steps undertaken are:
 * 1) bin the matches
 * 2) compute the pvalues of the tested region around the primary
 * 3) compute the sequence id set of the most significant peak
 * 4) set the state of the motif to loaded.
 *
 **************************************************************************/
void process_matches(
  int margin, 
  int bin, 
  double significance_threshold, 
  int test_max,
  int total_tests,
  MOTIF_T *primary_motif, 
  RBTREE_T *sequences,
  SECONDARY_MOTIF_T *secondary_motif,
  int *secondary_matches
) {
  bin_matches(margin, bin, sequences, primary_motif, secondary_motif, secondary_matches);
  compute_spacing_pvalues(margin, bin, total_tests, test_max, significance_threshold, secondary_motif);
  compute_idset(margin, bin, sequences, primary_motif, secondary_motif, secondary_matches);
  secondary_motif->loaded = TRUE; 
}
