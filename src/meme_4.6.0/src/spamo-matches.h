#ifndef SPAMO_MATCHES_H
#define SPAMO_MATCHES_H

#include "red-black-tree.h"
#include "linked-list.h"
#include "motif.h"

typedef struct SEQUENCE_DB SEQUENCE_DB_T;
typedef struct SEQUENCE SEQUENCE_T;
typedef struct MOTIF_DB MOTIF_DB_T;
typedef struct SPACING SPACING_T;
typedef struct SIGSPACE SIGSPACE_T;
typedef struct SECONDARY_MOTIF SECONDARY_MOTIF_T;
typedef struct SECONDARY_KEY SECONDARY_KEY_T;
typedef struct GROUPED_MOTIF GROUPED_MOTIF_T;

#define SAME 0
#define OPPO 1
#define LEFT 0
#define RIGHT 2

#define STRAND(quad) ((quad) & OPPO)
#define SIDE(quad) ((quad) & RIGHT)

/*
 * To conserve space the different values have the following meanings:
 * 0            : there is no match
 * value > 0    : there is a match on the positive strand at value
 * value < 0    : there is a match on the negative strand at abs(value)
 */

struct SEQUENCE_DB {
  char *source;
  char *name;
  time_t last_mod;
  int loaded;
  int excluded_tooshort;
  int excluded_nomatch;
  int excluded_similar;
};

struct SEQUENCE {
  int index;
  int length;
  char *data;
  char *name;
  int primary_match;
};

struct MOTIF_DB {
  int id; //an id assigned to the database based on the order it was listed
  char *source;
  char *name;
  char *cisml;
  time_t last_mod;
  int loaded;
  int excluded;
};

struct SPACING {
  int *bins;
  double *pvalues;
};

struct SIGSPACE {
  double pvalue;
  int bin;
  int quad;
};

struct SECONDARY_MOTIF {
  MOTIF_DB_T *db;
  MOTIF_T *motif;
  BOOLEAN_T loaded;
  SPACING_T spacings[4];
  int total_spacings;
  int max_in_one_bin;
  SIGSPACE_T *sigs; //significant spacings
  int sig_count;
  int *seqs;       // the ids of the best sequences
  int seq_count;
};

struct SECONDARY_KEY {
  int db_id;
  char *motif_id;
};

struct GROUPED_MOTIF {
  SECONDARY_MOTIF_T *best;
  LINKLST_T *others;
};

/**************************************************************************
 * Compares two keys of a secondary motif for equality.
 * First compares based on the motif database and second compares
 * based on the name of the motif.
 **************************************************************************/
int secondary_key_compare(void *p1, void *p2);

/**************************************************************************
 * Copies a key of a secondary motif
 **************************************************************************/
void* secondary_key_copy(void *p);

/**************************************************************************
 * Creates a structure to hold the motif database information.
 * All strings are copied.
 **************************************************************************/
SEQUENCE_DB_T* create_sequence_db(char *file);

/**************************************************************************
 * Destroy the sequence db
 **************************************************************************/
void destroy_sequence_db(SEQUENCE_DB_T *db);

/**************************************************************************
 * Creates a structure to hold the secondary database information.
 **************************************************************************/
MOTIF_DB_T* create_motif_db(int id, char *file, char *cisml);

/**************************************************************************
 * Destroy the secondary db
 **************************************************************************/
void destroy_motif_db(void *secondary_db);

/**************************************************************************
 * Destroy the sequence (including the name attribute)
 **************************************************************************/
void destroy_sequence(void *sequence);

/**************************************************************************
 * Create a secondary motif. As the number of sequences is unknown at this
 * point the sequence_matches array is left unallocated. All pvalues are
 * initilized to 1.
 **************************************************************************/
SECONDARY_MOTIF_T* create_secondary_motif(int margin, int bin, 
    MOTIF_DB_T *db, MOTIF_T *motif);

/**************************************************************************
 * Destroys a secondary motif. 
 * It takes a void * pointer so it can be used in the collection objects.
 **************************************************************************/
void destroy_secondary_motif(void *);

/**************************************************************************
 * Create a group with a good secondary motif and a group
 * of redundant secondary motifs.
 **************************************************************************/
GROUPED_MOTIF_T* create_grouped_motif(SECONDARY_MOTIF_T* best);

/**************************************************************************
 * Destroys a grouped motif
 * 
 * Relies on the secondary motifs begin destroyed elsewhere.
 **************************************************************************/
void destroy_grouped_motif(void *p);

/**************************************************************************
 * Calculate the total number of pvalue calculations that will be done
 * by the program. This number is used to correct the pvalues for multiple
 * tests using a bonferoni correction.
 **************************************************************************/
int calculate_test_count(
  int margin, 
  int bin, 
  int test_max, 
  RBTREE_T *secondary_motifs
);

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
);

#endif
