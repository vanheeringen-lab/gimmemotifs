/**************************************************************************
 * FILE: tomtom.c
 * AUTHOR: Shobhit Gupta, Timothy Bailey, William Stafford Noble
 * CREATE DATE: 02-21-06
 * PROJECT: TOMTOM
 * DESCRIPTION: Motif-Motif comparison using the score statistic.
 **************************************************************************/
#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#include "array-list.h"
#include "dir.h"
#include "motif.h"       // Needed for complement_dna_freqs.
#include "matrix.h"      // Routines for floating point matrices.
#include "array.h"       // Routines for floating point arrays.
#include "alphabet.h"    // The alphabet.
#include "meme-io.h"
#include "simple-getopt.h"
#include "projrel.h"
#include "pssm.h"
#include "fitevd.h"
#include "config.h"
#include "ceqlogo.h"
#include "io.h"
#include "qvalue.h"
#include "string-list.h"
#include "object-list.h"
#include "red-black-tree.h"
#include "utils.h"
#include "xml-out.h"
#include "xml-util.h"

#define DEFAULT_QUERY_PSEUDO   0.0
#define DEFAULT_TARGET_PSEUDO   0.0
#define BINS 100
#define LOGOHEIGHT 10                // height of an alignment of logos in cm.
#define MIN_TARGET_DATABASE_SIZE 50

/* Define some string limits (these include the EOL).*/
// FIXME: these are brittle and will cause bugs if motifs are too long
#define MAX_MOTIF_LENGTH 100
#define MAX_LINE_LENGTH 100

static const char *XML_FILENAME = "tomtom.xml";
static const char *HTML_STYLESHEET = "tomtom-to-html.xsl";
static const char *HTML_FILENAME = "tomtom.html";
static const char *OLD_HTML_FILENAME = "old_tomtom.html";
static const char *TXT_FILENAME = "tomtom.txt";
const char* tomtom_dtd = 
"<?xml version='1.0' encoding='UTF-8' standalone='yes'?>\n"
"<!DOCTYPE tomtom[\n"
"<!ELEMENT tomtom (model, targets, queries, runtime)>\n"
"<!ATTLIST tomtom version CDATA #REQUIRED release CDATA #REQUIRED>\n"
"<!ELEMENT model (command_line, distance_measure, threshold, background, host, when)>\n"
"<!ELEMENT command_line (#PCDATA)>\n"
"<!ELEMENT distance_measure EMPTY>\n"
"<!ATTLIST distance_measure value (allr|blic1|blic5|ed|kullback|llr1|llr5|pearson|sandelin) #REQUIRED>\n"
"<!ELEMENT threshold (#PCDATA)>\n"
"<!ATTLIST threshold type (evalue|qvalue) #REQUIRED>\n"
"<!ELEMENT background EMPTY>\n"
"<!ATTLIST background from (first_target|file) #REQUIRED A CDATA #REQUIRED C CDATA #REQUIRED G CDATA #REQUIRED T CDATA #REQUIRED file CDATA #IMPLIED>\n"
"<!ELEMENT host (#PCDATA)>\n"
"<!ELEMENT when (#PCDATA)>\n"
"<!-- each target is listed in order that the target was specified to the command line\n"
"     the motifs are not listed in any particular order -->\n"
"<!ELEMENT targets (target_file+)>\n"
"<!ELEMENT target_file (motif*)>\n"
"<!ATTLIST target_file index CDATA #REQUIRED source CDATA #REQUIRED name CDATA #REQUIRED \n"
"  loaded CDATA #REQUIRED excluded CDATA #REQUIRED last_mod_date CDATA #REQUIRED>\n"
"<!-- currently there can only be one query file (but users should not assume this will always be true)\n"
"     the query motifs are specified in the order that they appear in the file\n"
"     the matches are ordered from best to worst -->\n"
"<!ELEMENT queries (query_file+)>\n"
"<!ELEMENT query_file (query+)>\n"
"<!ATTLIST query_file source CDATA #REQUIRED name CDATA #REQUIRED last_mod_date CDATA #REQUIRED>\n"
"<!ELEMENT query (motif,match*)>\n"
"<!ELEMENT match EMPTY>\n"
"<!ATTLIST match target IDREF #REQUIRED orientation (forward|reverse) \"forward\" \n"
"  offset CDATA #REQUIRED pvalue CDATA #REQUIRED evalue CDATA #REQUIRED qvalue CDATA #REQUIRED>\n"
"<!-- motif contains the probability of each of the nucleotide bases at each position;\n"
"     i starts at 1; A, C, G and T are probabilities that sum to 1 -->\n"
"<!ELEMENT motif (pos+)>\n"
"<!ATTLIST motif id ID #REQUIRED name CDATA #REQUIRED alt CDATA #IMPLIED length CDATA #REQUIRED \n"
"  nsites CDATA #IMPLIED evalue CDATA #IMPLIED url CDATA #IMPLIED>\n"
"<!ELEMENT pos EMPTY>\n"
"<!ATTLIST pos i CDATA #REQUIRED A CDATA #REQUIRED C CDATA #REQUIRED G CDATA #REQUIRED T CDATA #REQUIRED>\n"
"<!ELEMENT runtime EMPTY>\n"
"<!ATTLIST runtime cycles CDATA #REQUIRED seconds CDATA #REQUIRED>\n"
"]>\n";



typedef struct tomtom_match TOMTOM_MATCH_T;
typedef struct motif_db MOTIF_DB_T;

/**************************************************************************
 * An object for storing information about one query-target alignment.
 **************************************************************************/
struct tomtom_match {
  MOTIF_T *query;   // the query motif
  MOTIF_T *target;  // the target motif, may be the reverse complement
  MOTIF_T *original;// the non reverse complemented motif
  MOTIF_DB_T *db;
  char  orientation;
  int   offset;
  int   overlap;
  float pvalue;
  float evalue;
  float qvalue;
  int   offset_by_length;
};


/**************************************************************************
 * An object for storing information about a file containing a set of motifs
 **************************************************************************/
struct motif_db {
  int id; //an id assigned to the database based on the order it was listed
  char *source;
  char *name;
  time_t last_mod;
  int loaded;
  int excluded;
  int list_index;//the index that in that list of motifs that this db starts at
  int list_entries;//the number of entries in the list of motifs from this db
  RBTREE_T *matched_motifs;
};

/**************************************************************************
 * Create a new query-target alignment object.
 * Initially, there is no q-value and no orientation.
 **************************************************************************/
TOMTOM_MATCH_T* new_tomtom_match (
  MOTIF_T *query,
  MOTIF_T *target,
  MOTIF_T *original,
  MOTIF_DB_T *target_db,
  int     offset,
  int     overlap,
  float   pvalue,
  float   evalue,
  int     offset_by_length
) {
  TOMTOM_MATCH_T* return_value = mm_malloc(sizeof(TOMTOM_MATCH_T));

  return_value->query = query;
  return_value->target = target;
  return_value->original = original;
  return_value->db = target_db;
  return_value->offset = offset;
  return_value->overlap = overlap;
  return_value->pvalue = pvalue;
  return_value->evalue = evalue;
  return_value->offset_by_length = offset_by_length;
  return(return_value);
}

/**************************************************************************
 * Compare two strings.
 * Used by the red-black tree
 **************************************************************************/
static int compare_strings(void *v1, void *v2) {
  return strcmp((char*)v1, (char*)v2);
}

/**************************************************************************
 * Compare two ints.
 * Used by the red-black tree
 **************************************************************************/
static int compare_ints(void *v1, void *v2) {
  int i1 = *((int*)v1);
  int i2 = *((int*)v2);
  if (i1 < i2) return -1;
  else if (i1 == i2) return 0;
  else return 1;
}

void* copy_int(void *v) {
  int *i_data = (int*)mm_malloc(sizeof(int));
  *i_data = *((int*)v);
  return i_data;
}

/**************************************************************************
 * Create a new motif database entry.
 * Makes a copy of the strings so there is no confusion on what can be
 * freed later.
 **************************************************************************/
static MOTIF_DB_T* create_motif_db(int id, char *source, char *name, time_t now) {
  MOTIF_DB_T *db;

  assert(source != NULL);
  assert(name != NULL);

  db = mm_malloc(sizeof(MOTIF_DB_T));
  memset(db, 0, sizeof(MOTIF_DB_T));
  db->id = id;
  copy_string(&(db->source), source);
  copy_string(&(db->name), name);
  //if the input is from stdin, or we can't do any better than use now as the timestamp
  db->last_mod = now;
  if (strcmp(source, "-") != 0) {
    //update the last mod date to the correct value
    struct stat stbuf;
    int result;
    result = stat(source, &stbuf);
    if (result < 0) {
      //an error occured
      die("Failed to stat file \"%s\", error given as: %s\n", source, strerror(errno));
    }
    db->last_mod = stbuf.st_mtime;
  }
  // the strings that I will be using are fixed size buffers of the MOTIF_T object
  // and they won't be deallocated in the life of the database entry so I don't need
  // a key copy or key free function.
  db->matched_motifs = rbtree_create(compare_strings, NULL, NULL, NULL, NULL);

  return db;
}

/**************************************************************************
 * Destroys the motif database entry.
 **************************************************************************/
void destroy_motif_db(void *ptr) {
  MOTIF_DB_T *db = (MOTIF_DB_T*)ptr;
  rbtree_destroy(db->matched_motifs);
  free(db->source);
  free(db->name);
  memset(db, 0, sizeof(MOTIF_DB_T));
  free(db);
}

/**************************************************************************
 *
 * Creates a logo of the aligned target and query motif, writes it to an
 * EPS file and calls convert to convert the EPS file to a PNG image.
 *
 * Returns name of logo file (minus the extension).
 *
 * Caller is responsible for freeing memory.
 *
 **************************************************************************/
static char *create_logo(
  MOTIF_DB_T* target_db,       // Target Database
  MOTIF_T* target,             // Target motif
  MOTIF_T* query,              // Query motif
  int offset,                  // Offset of target rel. to query
  BOOLEAN_T ssc,               // Small sample correction
  char *output_dirname         // Output directory
) {
  char *logo_filename, *logo_path;
  int num_len, num;
  double logo_height = LOGOHEIGHT;
  double logo_width
    = offset <= 0
      ? MAX(query->length, (target->length - offset))
      : MAX((query->length + offset), target->length);

  // calculate the length of the id even though technically it should be 1 in a large number of cases
  num_len = 1;
  num = target_db->id;
  if (num < 0) {
    num_len++;
    num = -num;
  }
  while (num >= 10) {
    num_len++;
    num /= 10;
  }

  // create name of logo file
  logo_filename
    = (char*) malloc(
        sizeof(char) * (6 + strlen(query->id) + 1 + num_len + 1 + strlen(target->id) + 1)
      );
  sprintf(logo_filename, "align_%s_%d_%s", query->id, target_db->id, target->id);

  // create output file path
  logo_path = make_path_to_file(output_dirname, logo_filename);

  CL_create2(
    target,             // first motif
    target->id+1,       // title (skip +/- character)
    query,              // second motif
    query->id,          // title
    TRUE,               // error bars
    ssc,                // small sample correction
    logo_height,        // logo height (cm)
    logo_width,         // logo width (cm)
    get_alphabet(FALSE),// alphabet
    -offset,            // offset of second motif
    logo_path,          // output file path
    "Tomtom"            // program name
  );
  myfree(logo_path);

  return(logo_filename);
}


/**************************************************************************
*
* get_pv_lookup_new
*
* Create a lookup table for the pv.
*
* FIXME: document what this really does!
*
**************************************************************************/
void get_pv_lookup_new(
  MATRIX_T* pssm_matrix,      // The PSSM.
  int alen,                   // Length of alphabet.
  int range,                  // Range of scores.
  ARRAY_T* background,        // Background model (uniform).
  MATRIX_T* reference_matrix, // reference_matrix
  MATRIX_T* pv_lookup_matrix  // pv_lookup_matrix
) {

  int i, j, k;
  int w = get_num_rows(pssm_matrix);        // Width of PSSM.
  int size = w*range+1;
  int start;                                // Starting query motif column
  int address_index = 0;                // next free position in matrix
  ARRAY_T* pdf_old = allocate_array(size);
  ARRAY_T* pdf_new = allocate_array(size);

  // Compute pv tables for groups of motif columns [start...]
  for (start=w-1; start>=0; start--) {
    // Compute the pdf recursively.
    init_array(0, pdf_new);
    set_array_item(0, 1, pdf_new);    // Prob(0)
    for (i=0; start+i<w; i++) {
      int max = i * range;
      SWAP(ARRAY_T*, pdf_new, pdf_old)
      for (k=0; k<=max+range; k++) {
        set_array_item(k, 0, pdf_new);
      }
      for (j=0; j<alen; j++) {
        int s = (int) get_matrix_cell(start+i, j, pssm_matrix);
        for(k=0; k<=max; k++) {
          double old = get_array_item(k, pdf_old);
          if (old != 0) {
            double new = get_array_item(k+s, pdf_new) +
              (old * get_array_item(j, background));
            set_array_item(k+s, new, pdf_new);
          } // old
        } // k
      } // j

      // Compute 1-cdf for motif consisting of columns [start, start+i]
      // This is the p-value lookup table for those columns.
      ARRAY_T* pv = allocate_array(get_array_length(pdf_new));
      copy_array(pdf_new, pv);
      int ii;
      for (ii=size-2; ii>=0; ii--) {
        double p = get_array_item(ii, pv) + get_array_item(ii+1, pv);
        set_array_item(ii, MIN(1.0, p), pv);
      }

      // copy the pv lookup table into the lookup matrix
      set_matrix_row(address_index, pv, pv_lookup_matrix);
      free_array(pv);

      // Store the location of this pv lookup table in the reference array
      set_matrix_cell(
        start,
        start+i,
        (double) address_index++,
        reference_matrix
      );
    } // i
  } // start position

  // Free space.
  free_array(pdf_new);
  free_array(pdf_old);
} // get_pv_lookup_new

/**************************************************************************
 * Parse a list of Tomtom output strings, converting p-values to q-values.
 **************************************************************************/
static void convert_tomtom_p_to_q(OBJECT_LIST_T* output_list) {
  int i_line;
  STRING_LIST_T* completed_ids;

  // Extract all of the p-values.
  ARRAY_T* pvalues = get_object_list_scores(output_list);

  // Convert p-values to q-values.
  compute_qvalues(
    FALSE,  // Don't stop with FDR.
    TRUE,   // Estimate pi-zero.
    NULL,   // Don't store pi-zero in a file.
    NUM_BOOTSTRAPS,
    NUM_BOOTSTRAP_SAMPLES,
    NUM_LAMBDA,
    MAX_LAMBDA,
    get_array_length(pvalues), 
    pvalues,
    NULL    // No sampled p-values provided
  );
  ARRAY_T* qvalues = pvalues;

  // Traverse the list of matches.
  TOMTOM_MATCH_T* my_match = (TOMTOM_MATCH_T*)retrieve_next_object(output_list);
  int index = 0;
  while (my_match != NULL) {

    // Store the q-value.
    my_match->qvalue = get_array_item(index, qvalues);

    // Get the next match.
    my_match = (TOMTOM_MATCH_T*)retrieve_next_object(output_list);
    index++;
  }
  free_array(pvalues); // note: same array as qvalues
} // Function convert_tomtom_p_to_q


/**************************************************************************
 * Given a sorted list of output lines corresponding to one query
 * motif, filter the list such that it contains only the best-scoring
 * hit to each target.
 *
 * As a side effect, re-format the lines slightly.
 **************************************************************************/
static void select_strand(
  OBJECT_LIST_T* output_list,
  OBJECT_LIST_T* final_output
) {
  TOMTOM_MATCH_T* my_match;
  RBTREE_T *completed_ids;
  BOOLEAN_T is_new;

  completed_ids = rbtree_create(compare_strings, NULL, NULL, NULL, NULL);

  // Traverse the list of matches.
  my_match = (TOMTOM_MATCH_T*)retrieve_next_object(output_list);
  while (my_match != NULL) {

    // Get the orientation of the target motif.
    my_match->orientation = '-';
    if (get_motif_id(my_match->target)[0] == '+') {
      my_match->orientation = '+';
    }

    // Have we already seen a better scoring match? (note: skipping over +/- at the start)
    rbtree_lookup(completed_ids, get_motif_id(my_match->target)+1, TRUE, &is_new);

    if (is_new) {
      // Make a copy on the heap.
      TOMTOM_MATCH_T* match_copy = mm_malloc(sizeof(TOMTOM_MATCH_T));
      *match_copy = *my_match;

      // Add this match to the new list.
      store_object((void*)match_copy, NULL, my_match->pvalue, final_output);
    }

    // Get the next match.
    my_match = (TOMTOM_MATCH_T*)retrieve_next_object(output_list);
  }

  /* Free memory */
  rbtree_destroy(completed_ids);

} // Function select_strand

/**************************************************************************
 * Gets the consensus DNA sequence from a pssm.
 * caller is responsible for deallocating memory
 **************************************************************************/
static char* get_cons(MATRIX_T* freqs){
  char *nucleotide_syms = "ACGT";
  char *aminoacid_syms = "ACDEFGHIKLMNPQRSTVWY";
  char *syms, *cons;

  int i_row, i_col, size_alph, size_len;

  size_alph = get_num_cols(freqs);
  size_len = get_num_rows(freqs);

  cons = (char*)mm_malloc(sizeof(char) * (size_len + 1));

  if (size_alph == 4) {
    syms = nucleotide_syms;
  } else if (size_alph == 20) {
    syms = aminoacid_syms;
  } else {
    die("get_cons was passed a matrix with a size that didn't match one of the known alphabets.\n");
    syms = NULL; //silly compiler warning
  }

  for (i_row = 0; i_row < size_len; i_row++) {
    int best_i;
    float max;
    best_i = -1;
    max = 0;
    //find the column with the best score
    for (i_col = 0; i_col < size_alph; i_col++) {
      if (get_matrix_cell(i_row, i_col, freqs) > max) {
        max = get_matrix_cell(i_row, i_col, freqs);
        best_i = i_col;
      }
    }
    //append the corresponding letter to the consensus
    cons[i_row] = syms[best_i];
  }
  cons[size_len] = '\0';

  return cons;
} /* Function get_cons */


/**************************************************************************
 * Computes the score for a particular configuration
 **************************************************************************/
static double compute_overlap_score(
  MATRIX_T* pairwise_column_scores,
  int current_configuration,
  int *index_start,
  int *index_stop 
) {

  /* Initialize the columns from where the comparison needs to be started */
  int current_target_column = 0;
  int current_query_column = 0;
  if (current_configuration <= pairwise_column_scores->num_rows) {
    current_target_column = 0;
    current_query_column = pairwise_column_scores->num_rows
      - current_configuration;
  } else {
    current_query_column = 0;
    current_target_column = current_configuration
      - pairwise_column_scores->num_rows;
  }

  // Sum the pairwise column scores.
  double score = 0;
  *index_start = current_query_column;
  do {
    score += get_matrix_cell(
               current_query_column++,
               current_target_column++,
               pairwise_column_scores
             );
  } while (
    (current_query_column < pairwise_column_scores->num_rows)
      && (current_target_column < pairwise_column_scores->num_cols)
  );
  *index_stop = current_query_column - 1; // -1 for the last increment

  return(score);
} /* Function compute_overlap_score */

/**************************************************************************
 * Computes the optimal offset and score corresponding to the smallest
 * pvalue.
 **************************************************************************/
static void compare_motifs(
  MATRIX_T* pairwise_column_scores,
  int*      optimal_offset,
  MATRIX_T* reference_matrix,
  MATRIX_T* pv_lookup_matrix,
  double*   optimal_pvalue,
  int*      overlapping_bases,
  BOOLEAN_T internal, // Is shorter motif contained in longer motif?
  int       min_overlap
) {

  int first_pass = 1;

  /*  Slide one profile over another */
  int total_configurations = pairwise_column_scores->num_rows
      + pairwise_column_scores->num_cols - 1;
  int current_configuration;
  int overlap = -1;
  for (current_configuration = 1;
      current_configuration <= total_configurations;
      current_configuration++) {

    // If requested, check whether one motif is contained within the other.
    BOOLEAN_T process = TRUE;
    if (internal) {
      int max, min;
      max = pairwise_column_scores->num_rows;
      min = pairwise_column_scores->num_cols;
      if (max < min) {
        int temp;
        temp = min;
        min = max;
        max = temp;
      }
      if (!((current_configuration <= max)
            && (current_configuration >= min))) {
        process = FALSE;
      }
    }

    // Override min-overlap if query/target is smaller than min-overlap
    int mo;
    mo = min_overlap;
    if (pairwise_column_scores->num_rows < mo) {
      mo = pairwise_column_scores->num_rows;
    }
    if (pairwise_column_scores->num_cols < mo) {
      mo = pairwise_column_scores->num_cols;
    }

    // Check if this configuration should be processed (min_overlap).
    if (mo > 1) {
      int max, min;
      min = mo;
      max = pairwise_column_scores->num_rows
        + pairwise_column_scores->num_cols - mo;
      if (!((current_configuration <= max)
            && (current_configuration >= min))) {
        process = FALSE;
      }
    }

    /* Check whether the requested percentage of information content is 
       contained within the matching portion of each motif. */
    /* FIXME: Stub to be filled in.  WSN 4 Apr 2008
    info_content = MIN(compute_info_content(start1, stop1, motif1),
                       compute_info_content(start2, stop2, motif2));
    if (info_content < target_info_content) {
      process = FALSE;
    }
    */

    if (process) {
      double pvalue = 1;
      double TEMP_SCORE = 0;
      int current_overlap = -1;
      int index_start, index_stop;

      /* Compute the score for the current configuration */
      TEMP_SCORE 
        = compute_overlap_score(
            pairwise_column_scores,
            current_configuration,
            &index_start,
            &index_stop
          );

      /* Compute the pvalue of the score */
      int curr_address_index = (int) get_matrix_cell(
                   index_start,
                   index_stop,
                   reference_matrix);
      // FIXME: what does this code do?
      //if (curr_address_index == -1) {
      //  curr_address_index = address_index - 1;
      //}
      pvalue = get_matrix_cell(
        curr_address_index,
        (int) TEMP_SCORE,
        pv_lookup_matrix
      );

      current_overlap = index_stop - index_start + 1;
      if (current_overlap > 0) {
        /* Keep the configuration with the lowest pvalue.
           If pvalue is the same, keep the largest overlap */
        if (first_pass == 1) {
          *optimal_offset = current_configuration;
          *optimal_pvalue = pvalue;
          overlap = current_overlap;
          first_pass = 0;
        } else {
          if (pvalue == *optimal_pvalue) {
            if (overlap < current_overlap) {
              *optimal_offset = current_configuration;
              *optimal_pvalue = pvalue;
              overlap = current_overlap;
            }
          } else {
            if (pvalue < *optimal_pvalue) {
              *optimal_offset = current_configuration;
              *optimal_pvalue = pvalue;
              overlap = current_overlap;
            }
          }
        }/* For the pvalue minimization condition */
      } /* Overlap > 0 (happens sometimes if min-overlap is overridden by
           target motif length). */
    } /* Process */
  } /* Loop over configurations. */
  *overlapping_bases = overlap;
} /* Function compare_motifs */

/**************************************************************************
 * Computes ALLR scores for all columns of the query matrix against target
 * motif. This can then be used to by compute_allr.
 **************************************************************************/
void allr_scores(
  MOTIF_T* query_motif,
  MOTIF_T* target_motif,
  ARRAY_T* background,
  MATRIX_T** score
) {
  /*   Recurse pairwise over the columns */
  int index_query;
  int index_target;
  for (index_query = 0; index_query < query_motif->length; index_query++) {
    for (index_target = 0; index_target < target_motif->length; index_target++) {

      double SCORE_F = 0;
      double SCORE_NUMERATOR = 0;
      double SCORE_DENOMINATOR = 0;
      int index; /* Index corresponds to the alphabets. */

      for (index = 0; index < query_motif->alph_size; index++) {
        /* The numerator of the ALLR */
        double nr_part1, nr_part2;
        /*  Part 1 of the numerator  */
        nr_part1
          = (
              target_motif->num_sites
                * (
                    get_matrix_cell( index_target, index, target_motif->freqs)
                    /* nb for target */
                  )
                * log(
                    get_matrix_cell(index_query, index, query_motif->freqs)
                    / get_array_item(index, background)
                  )
          ); /* Likelihood for query */
        /*  Part 2 of the numerator */
        nr_part2 
          = (
              query_motif->num_sites
                * (
                    get_matrix_cell(index_query, index, query_motif->freqs)
                    /* nb for query */
                  )
                * log(
                    get_matrix_cell(index_target, index, target_motif->freqs)
                    / get_array_item(index, background)
                  )
          ); /* Likelihood for target */

        /* Sum of the two parts. */
        SCORE_NUMERATOR += nr_part1 + nr_part2;

        /*  The denominoator of the ALLR */
        SCORE_DENOMINATOR
          += (
               (
                 query_motif->num_sites
                   *  get_matrix_cell(
                        index_query, 
                        index, 
                        query_motif->freqs
                       )
                  // nb for query
                )
            + (
                target_motif->num_sites
                * get_matrix_cell(index_target, index, target_motif->freqs)
              ) // nb for target
           ); /* Sum of nb for query and nb for target */

      } /* Nr and Dr of the score summed over the entire alphabet. */
      SCORE_F = SCORE_NUMERATOR / SCORE_DENOMINATOR;
      set_matrix_cell(index_query, index_target, SCORE_F, *score);

    } /* Target motif */
  } /* Query motif */
} /* Function ALLR scores */

/**************************************************************************
 * Computes euclidian distance between all columns of the query matrix
 * This score is not normalized on length and cannot be used without a
 * normalization criterion as our DP based p-value approach.
 **************************************************************************/
void ed_scores(
  MOTIF_T* query_motif,
  MOTIF_T* target_motif,
  ARRAY_T* background,
  MATRIX_T** score
) {

  /*   Recurse pairwise over the columns */
  int index_query;
  int index_target;
  for (index_query = 0; index_query < query_motif->length; index_query++) {
    for (index_target = 0; index_target < target_motif->length; index_target++) {

      double SCORE_T = 0;

      int index; /* Index corresponds to the alphabets. */
      for (index = 0; index < query_motif->alph_size; index++) {

        double query_freq = get_matrix_cell(index_query, index, query_motif->freqs);
        double target_freq = get_matrix_cell(index_target, index, target_motif->freqs);
        /* Euclidean distance */
        double sq;
        sq = pow((query_freq - target_freq), 2);
        SCORE_T = SCORE_T + sq;
      }
      double SCORE_F;
      SCORE_F = - sqrt(SCORE_T);
      set_matrix_cell(index_query, index_target, SCORE_F, *score);
    } /* Target motif */
  } /* Query motif */
} /* Function ed_scores */

/**************************************************************************
 * Computes sandelin distance between all columns of the query matrix
 * This score is not normalized on length and cannot be used without a
 * normalization criterion as our DP based p-value approach.
 **************************************************************************/
void sandelin_scores(
  MOTIF_T* query_motif,
  MOTIF_T* target_motif,
  ARRAY_T* background,
  MATRIX_T** score
) {

  /*   Recurse pairwise over the columns */
  int index_query;
  int index_target;
  for (index_query = 0; index_query < query_motif->length; index_query++) {
    for (index_target = 0; index_target < target_motif->length; index_target++) {

      double SCORE_T = 0;

      int index; /* Index corresponds to the alphabets. */
      for (index = 0; index < query_motif->alph_size; index++) {

        double query_freq
          = get_matrix_cell(index_query, index, query_motif->freqs);
        double target_freq 
          = get_matrix_cell(index_target, index, target_motif->freqs);
        /* Euclidean distance */
        double sq;
        sq = pow((query_freq - target_freq), 2);
        SCORE_T = SCORE_T + sq;
      }
      double SCORE_F;
      SCORE_F = 2 - SCORE_T;
      set_matrix_cell(index_query, index_target, SCORE_F, *score);
    } /* Target motif */
  } /* Query motif */
} /* Function sandelin*/

/**************************************************************************
 * Computes Kullback-Leiber for all columns of the query matrix against target
 * motif. Note we use
 * a modified Kullback-leibler for computing single column scores. We do not
 * divide the score by the width. This division is replaced by our DP based
 * p-value computation.
 **************************************************************************/
void kullback_scores(
  MOTIF_T* query_motif,
  MOTIF_T* target_motif,
  ARRAY_T* background,
  MATRIX_T** score
) {
  /*   Recurse pairwise over the columns */
  int index_query;
  int index_target;
  for (index_query = 0; index_query < query_motif->length; index_query++) {
    for (index_target = 0; index_target < target_motif->length; index_target++) {
      double SCORE_F = 0;
      int index; /* Index corresponds to the alphabets. */
      for (index = 0; index < query_motif->alph_size; index++) {
        double query_freq
          = get_matrix_cell(index_query, index, query_motif->freqs);
        double target_freq
          = get_matrix_cell(index_target, index, target_motif->freqs);
        /* Average Kullback */
        double avg_kull;
        avg_kull = ((query_freq * log10(query_freq/target_freq))
          + (target_freq * log10(target_freq/query_freq))) / 2;
        SCORE_F = SCORE_F + avg_kull;
      }
      SCORE_F = SCORE_F * -1;
      set_matrix_cell(index_query, index_target, SCORE_F, *score);
    } /* Target motif */
  } /* Query motif */
} /* Function Kullback*/


/**************************************************************************
 * Computes Pearson scores for all columns of the query matrix against target
 * motif. 
 **************************************************************************/
void pearson_scores(
  MOTIF_T* query_motif,
  MOTIF_T* target_motif,
  ARRAY_T* background,
  MATRIX_T** score
) {

  /*   Define the average values */
  double x_bar = 1 / query_motif->alph_size;
  double y_bar = 1 / target_motif->alph_size;

  /*   Recurse pairwise over the columns */
  int index_query;
  int index_target;
  for (index_query = 0; index_query < query_motif->length; index_query++) {
    for (index_target = 0; index_target < target_motif->length; index_target++) {

      double SCORE_F = 0;
      double SCORE_NUMERATOR = 0;
      double SCORE_SQ_DENOMINATOR_1 = 0;
      double SCORE_SQ_DENOMINATOR_2 = 0;

      int index; /* Index corresponds to the alphabets. */
      for (index = 0; index < query_motif->alph_size; index++) {
        /* The numerator of the Pearson */
        SCORE_NUMERATOR
          += (get_matrix_cell(index_query, index, query_motif->freqs) - x_bar)
          * (get_matrix_cell(index_target, index, target_motif->freqs) - y_bar);

        /*  The denominoator of the Pearson */
        SCORE_SQ_DENOMINATOR_1
          += pow(
                (get_matrix_cell(index_query, index, query_motif->freqs)
                  - x_bar),
                2
             );
        SCORE_SQ_DENOMINATOR_2
          += pow(
                (get_matrix_cell(index_target, index, target_motif->freqs)
                  - y_bar),
                 2
             );
      } /* Nr and Dr components summed over the entire alphabet. */

      /*       Compute the Denominator */
      double SCORE_SQ_DENOMINATOR 
        = SCORE_SQ_DENOMINATOR_1 * SCORE_SQ_DENOMINATOR_2;
      double SCORE_DENOMINATOR = sqrt(SCORE_SQ_DENOMINATOR);

      SCORE_F = SCORE_NUMERATOR / SCORE_DENOMINATOR;
      set_matrix_cell(index_query, index_target, SCORE_F, *score);
    } /* Target motif */
  } /* Query motif */
} /* Function pearson */

/**************************************************************************
 * Estimate one Dirichlet component.
 *
 * {\hat p}_i
 *     = \frac{n_i + \alpha_i}{\sum_{j \in \{A,C,G,T\}} (n_j + \alpha_j)}
 **************************************************************************/
static ARRAY_T* one_dirichlet_component(
  int      hyperparameters[],
  ARRAY_T* counts
) {
  int alph_size = get_array_length(counts);
  ARRAY_T* return_value = allocate_array(alph_size);

  // Calculate the denominator.
  double denominator = 0.0;
  int i_alph;
  for (i_alph = 0; i_alph < alph_size; i_alph++) {
    denominator += get_array_item(i_alph, counts)
      + (int)(hyperparameters[i_alph]);
  }

  // Estimate the source distribution.
  for (i_alph = 0; i_alph < alph_size; i_alph++) {
    double numerator = get_array_item(i_alph, counts) 
      + hyperparameters[i_alph];
    set_array_item(i_alph, numerator / denominator, return_value);
  }
  return(return_value);
}

/**************************************************************************
 * Helper function to get un-logged gamma.
 **************************************************************************/
static double true_gamma (double x) { return(exp(gamma(x))); }

/**************************************************************************
 * Likelihood of the data, given a Dirichlet component.
 *
 * This is the final equation on p. 13 of Habib et al.
 *
 * It should probably be implemented in log space, rather than using
 * the true_gamma function above.
 **************************************************************************/
static double compute_likelihood(
  int      hypers[],
  ARRAY_T* counts
) {

  // Sum the counts and the hyperparameters.
  double hyper_sum = 0.0;
  double count_sum = 0.0;
  int i_alph;
  int alphabet_size = get_array_length(counts);
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    hyper_sum += (float)(hypers[i_alph]);
    count_sum += get_array_item(i_alph, counts);
  }

  // Compute the product term.
  double product = 1.0;
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    double hyper = (float)(hypers[i_alph]);
    double count = get_array_item(i_alph, counts);
    product *= true_gamma(hyper + count) / true_gamma(hyper);
  }

  return((true_gamma(hyper_sum) / true_gamma(hyper_sum + count_sum))
         * product);
  
}

/**************************************************************************
 * Estimate the source distribution from counts, using a 1-component
 * or 5-component Dirichlet prior.
 *
 * The 5-component prior only works for DNA.
 **************************************************************************/
static ARRAY_T* estimate_source
 (int      num_dirichlets,
  ARRAY_T* counts)
{
  int dna_hypers[6][4] 
    = {
        {1, 1, 1, 1}, // Single component
        {5, 1, 1, 1}, // Five components ...
        {1, 5, 1, 1},
        {1, 1, 5, 1},
        {1, 1, 1, 5},
        {2, 2, 2, 2}
      };
  int protein_hypers[20]
    = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

  // Single component or 5-component mixture?
  int alph_size = get_array_length(counts);
  if (num_dirichlets == 1) {
    if (alph_size == 4) {
      return(one_dirichlet_component(dna_hypers[0], counts));
    } else if (alph_size == 20) {
      return(one_dirichlet_component(protein_hypers, counts));
    } else {
      die("Invalid alphabet size (%d) for BLiC score.\n", alph_size);
    }
  }

  // Die if trying to do proteins.
  if (alph_size != 4) {
    die("Sorry, BLiC5 is only implemented for DNA motifs.\n");
  }

  // Compute the denominator.
  double denominator = 0.0;
  int i_hyper;
  int n_hypers = 5;
  for (i_hyper = 1; i_hyper <= n_hypers; i_hyper++) {
    denominator 
      += (1.0 / (float) n_hypers)
      * compute_likelihood(dna_hypers[i_hyper], counts);
  }

  ARRAY_T* return_value = allocate_array(alph_size);
  for (i_hyper = 1; i_hyper <= n_hypers; i_hyper++) {

    /* Compute the posterior.
     *   \Pr(k|n) = \frac{\Pr(k) \Pr(n|k)}{\sum_j \Pr(j) \Pr(n|j)}
     */
    double posterior 
      = ((1.0 / (float) n_hypers) 
         * compute_likelihood(dna_hypers[i_hyper], counts))
      / denominator;

    // Compute the Dirichlet component.
    ARRAY_T* one_component
      = one_dirichlet_component(dna_hypers[i_hyper], counts);

    // Multiply them together and add to the return value.
    scalar_mult(posterior, one_component);
    sum_array(one_component, return_value);

    free_array(one_component);
  }

  return(return_value);
}

/**************************************************************************
 * Computes the BLiC (Bayesian Likelihood 2-Component) score described
 * by Equation (2) in Habib et al., PLoS CB 4(2):e1000010.
 *
 * Thus far, BLiC is only implemented for DNA.  It would be relatively
 * straightforward to extend to amino acids.
 **************************************************************************/
static double one_blic_score(
  int      num_dirichlets,
  ARRAY_T* query_counts, 
  ARRAY_T* target_counts,
  ARRAY_T* background
) {

  // Make a vector that sums the two sets of counts.
  int alphabet_size = get_array_length(query_counts);
  assert(alphabet_size == 4);  // BLiC only works for DNA.
  ARRAY_T* common_counts = allocate_array(alphabet_size);
  copy_array(query_counts, common_counts);
  sum_array(target_counts, common_counts);

  // Estimate the source distributions.
  ARRAY_T* query_source = estimate_source(num_dirichlets, query_counts);
  ARRAY_T* target_source = estimate_source(num_dirichlets, target_counts);
  ARRAY_T* common_source = estimate_source(num_dirichlets, common_counts);

  // First term.
  double first_term = 0.0;
  int i_alph;
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    first_term += 2 * get_array_item(i_alph, common_counts) 
      * log(get_array_item(i_alph, common_source));
  }

  // Second term.
  double second_term = 0.0;
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    second_term += get_array_item(i_alph, query_counts)
      * log(get_array_item(i_alph, query_source));
  }

  // Third term.
  double third_term = 0.0;
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    third_term += get_array_item(i_alph, target_counts)
      * log(get_array_item(i_alph, target_source));
  }

  // Fourth term.
  double fourth_term = 0.0;
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    fourth_term += get_array_item(i_alph, common_counts)
      * log(get_array_item(i_alph, background));
  }

  return(first_term - (second_term + third_term + fourth_term));

}

/**************************************************************************
 * Computes BLiC scores for all pairs of positions in two motifs.
 **************************************************************************/
static void blic_scores(
  int num_dirichlets,     // Number of Dirichlets [IN]
  MOTIF_T* query_motif,   // Single query motif [IN].
  MOTIF_T* target_motif,  // Single target motif [IN].
  ARRAY_T* background,    // Background distribution [IN].
  MATRIX_T** scores)      // BLiC score matrix, indexed by start. [OUT]
{
  // Traverse the query and target motifs
  int query_i;
  for (query_i = 0; query_i < query_motif->length; query_i++) {
    ARRAY_T* query_counts = get_motif_counts(query_i, query_motif);

    int target_i;
    for (target_i = 0; target_i < target_motif->length; target_i++) {
      ARRAY_T* target_counts = get_motif_counts(target_i, target_motif);

      // Compute the BLiC score according to Equation (2).
      double my_score 
        = one_blic_score(
            num_dirichlets,
            query_counts,
            target_counts,
            background
          );
      free_array(target_counts);

      // Store the computed score.
      set_matrix_cell(query_i, target_i, my_score, *scores);
    }
    free_array(query_counts);
  }
}

/**************************************************************************
 * Computes BLiC scores for all pairs of positions in two motifs.
 * The number of dirichlets is set to 1.
 **************************************************************************/
static void blic1_scores(
  MOTIF_T* query_motif,   // Single query motif [IN].
  MOTIF_T* target_motif,  // Single target motif [IN].
  ARRAY_T* background,    // Background distribution [IN].
  MATRIX_T** scores       // BLiC score matrix, indexed by start. [OUT]
) {
  blic_scores(1, query_motif, target_motif, background, scores);
}

/**************************************************************************
 * Computes BLiC scores for all pairs of positions in two motifs.
 * The number of dirichlets is set to 5.
 **************************************************************************/
static void blic5_scores(
  MOTIF_T* query_motif,  // Single query motif [IN].
  MOTIF_T* target_motif, // Single target motif [IN].
  ARRAY_T* background,   // Background distribution [IN].
  MATRIX_T** scores      // BLiC score matrix, indexed by start. [OUT]
) {
  blic_scores(5, query_motif, target_motif, background, scores);
}

/**************************************************************************
 * Computes the llr score 
 **************************************************************************/
static double one_llr_score(
  int      num_dirichlets,
  ARRAY_T* query_counts, 
  ARRAY_T* target_counts,
  ARRAY_T* background
) {

  // Make a vector that sums the two sets of counts.
  int alphabet_size = get_array_length(query_counts);
  assert(alphabet_size == 4);  // BLiC only works for DNA.
  ARRAY_T* common_counts = allocate_array(alphabet_size);
  copy_array(query_counts, common_counts);
  sum_array(target_counts, common_counts);

  // Estimate the source distributions.
  ARRAY_T* query_source = estimate_source(num_dirichlets, query_counts);
  ARRAY_T* target_source = estimate_source(num_dirichlets, target_counts);
  ARRAY_T* common_source = estimate_source(num_dirichlets, common_counts);

  // First term.
  double first_term = 0.0;
  int i_alph;
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    first_term += get_array_item(i_alph, common_counts) 
      * log(get_array_item(i_alph, common_source));
  }

  // Second term.
  double second_term = 0.0;
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    second_term += get_array_item(i_alph, query_counts)
      * log(get_array_item(i_alph, query_source));
  }

  // Third term.
  double third_term = 0.0;
  for (i_alph = 0; i_alph < alphabet_size; i_alph++) {
    third_term += get_array_item(i_alph, target_counts)
      * log(get_array_item(i_alph, target_source));
  }

  return(first_term - (second_term + third_term));

}

/**************************************************************************
 * Computes log likelihood ratio scores for all pairs of positions in two motifs.
 **************************************************************************/
static void llr_scores(
  int num_dirichlets,     // Number of Dirichlets [IN]
  MOTIF_T* query_motif,   // Single query motif [IN].
  MOTIF_T* target_motif,  // Single target motif [IN].
  ARRAY_T* background,    // Background distribution [IN].
  MATRIX_T** scores       // BLiC score matrix, indexed by start. [OUT]
) {
  // Traverse the query and target motifs
  int index_query;
  for (index_query = 0; index_query < query_motif->length; index_query++) {
    ARRAY_T* query_counts = get_motif_counts(index_query, query_motif);

    int index_target;
    for (index_target = 0; index_target < target_motif->length; index_target++) {
      ARRAY_T* target_counts = get_motif_counts(index_target, target_motif);


      double my_score 
        = one_llr_score(
            num_dirichlets,
            query_counts,
            target_counts,
            background
          );
      free_array(target_counts);

      // Store the computed score.
      set_matrix_cell(index_query, index_target, my_score, *scores);
    }
    free_array(query_counts);
  }
}

/**************************************************************************
 * Computes log likelihood ratio scores for all pairs of positions in two motifs.
 * The number of Dirichlets is set to 1.
 **************************************************************************/
static void llr1_scores
  (MOTIF_T* query_motif,  // Single query motif [IN].
   MOTIF_T* target_motif, // Single target motif [IN].
   ARRAY_T* background,   // Background distribution [IN].
   MATRIX_T** scores)     // BLiC score matrix, indexed by start. [OUT]
{
  llr_scores(1, query_motif, target_motif, background, scores);
}

/**************************************************************************
 * Computes log likelihood ratio scores for all pairs of positions in two motifs.
 * The number of Dirichlets is set to 5.
 **************************************************************************/
static void llr5_scores
  (MOTIF_T* query_motif,  // Single query motif [IN].
   MOTIF_T* target_motif, // Single target motif [IN].
   ARRAY_T* background,   // Background distribution [IN].
   MATRIX_T** scores)     // BLiC score matrix, indexed by start. [OUT]
{
  llr_scores(5, query_motif, target_motif, background, scores);
}

/*****************************************************************************
 * prints the header of the html
 * 
 * DEPRECATED - the html is now generated by a stylesheet from the xml.
 *****************************************************************************/
void print_tomtom_html_header(
  FILE *output_path,
  char *query_filename,
  ARRAYLST_T *target_dbs,
  char *distance_measure,
  BOOLEAN_T sig_type_q,
  float sig_thresh
)
{
  int i;
  char *target_filename;
  fprintf(output_path,
    "<html>\n"
    "<head>\n"
    "  <title> Tomtom output</title>\n"
    "</head> <body>\n" 
    "<table border=1>\n"
    "  <tr> <th align=center colspan=2><h2>TOMTOM OUTPUT</h2></th></tr>\n"
    "  <tr> <td align=left colspan=2><b>Query File:</b> %s</td></tr>\n",
    query_filename);
  fprintf(output_path, "  <tr> <td align=left colspan=2><b>Target File:</b>");
  for (i = 0; i < arraylst_size(target_dbs); ++i) {
    target_filename = ((MOTIF_DB_T*)arraylst_get(i, target_dbs))->name;
    fprintf(output_path, " %s", target_filename);
  }
  fprintf(output_path, "</td></tr>\n");

  fprintf(output_path,
    "  <tr> <td align=left colspan=2><b>Distance Measure:</b> %s</td></tr>\n"
    "  <tr> <td align=left colspan=2><b>All Motif Matches with <i>%s</i>-value at most:</b> %.2g.",
    distance_measure,
    (sig_type_q ? "q" : "E"),
    sig_thresh
  );
  if (sig_type_q) {
    fprintf(output_path,
      " The <i>q</i>-value is the estimated false discovery rate if the occurrence"
      " is accepted as significant."
      " See Storey JD, Tibshirani R. \"Statistical significance for genome-wide studies\"."
      " <i>Proc. Natl Acad. Sci. USA (2003) 100:9440--9445</i>"
    );
  }
  fprintf(output_path, "</td></tr>\n");
} // print_tomtom_html_header

/*****************************************************************************
 * prints the html for a tomtom match
 * 
 * DEPRECATED - the html is now generated by a stylesheet from the xml.
 *****************************************************************************/
void print_tomtom_html_line(
  FILE *output_fp,              // path to output XML file
  TOMTOM_MATCH_T* match,        // stored results line
  char *logo_filename           // name of logo
  )
{
  MOTIF_T *target_motif, *query_motif;
  int query, target, offset, overlap;
  char query_id[MAX_MOTIF_ID_LENGTH+1], target_id[MAX_MOTIF_ID_LENGTH+1];
  double evalue, uncorrected_pvalue;
  char orientation;

  query_motif = match->query;
  target_motif = match->target;


  fprintf(output_fp,
    "  <tr>\n"
    "    <td>\n"
    "      <table>\n"
    "        <tr><th align=left nowrap>Target Motif:</th><td align=right>");
  // target id
  // skip over +/- at start of id
  if (target_motif->url && target_motif->url[0] != '\0') {
    fprintf(output_fp, "<a href=\"%s\">%s</a>", target_motif->url, get_motif_id(target_motif)+1);
  } else {
    fprintf(output_fp, "%s", get_motif_id(target_motif)+1);
  }
  fprintf(output_fp,
    "</td></tr>\n"
    "        <tr><th align=left nowrap>Target Description:</th><td align=right>%s</td></tr>\n"        // target description
    "        <tr><th align=left nowrap>Query Motif:</th><td align=right>",
    target_motif->id2);        // target description
    
  // query id
  if (query_motif->url && query_motif->url[0] != '\0') {
    fprintf(output_fp, "<a href=\"%s\">%s</a>", query_motif->url, get_motif_id(query_motif));
  } else {
    fprintf(output_fp, "%s", get_motif_id(query_motif));
  }
  fprintf(output_fp,
    "</td></tr>\n"
    "        <tr><th align=left nowrap>Query Description:</th><td align=right>%s</td></tr>\n"        // query description
    "        <tr><th align=left nowrap><i>p</i>-value:</th><td align=right><b>%.2g<b></td></tr>\n"        // p-value
    "        <tr><th align=left nowrap><i>E</i>-value:</th><td align=right><b>%.2g<b></td></tr>\n"        // E-value
    "        <tr><th align=left nowrap><i>q</i>-value:</th><td align=right><b>%.2g<b></td></tr>\n"        // q-value
    "        <tr><th align=left nowrap>Overlap:</th><td align=right>%d</td></tr>\n"                // overlap
    "        <tr><th align=left nowrap>Query Offset:</th><td align=right>%d</td></tr>\n"        // offset
    "        <tr><th align=left nowrap>Orientation:</th><td align=right>%c</td></tr>\n"                // orientation
    "        <tr><th align=left nowrap>Figures:</th>"                                                // links
    "<td align=right><a href='%s.eps'>[EPS]</a><a href='%s.png'>[PNG]</a></td></tr>\n"
    "      </table>\n"
    "    </td>\n"
    "    <td><img src='%s.png'></td>\n"                                                        // logo
    "  </tr>\n",
    query_motif->id2,        // query description
    match->pvalue,
    match->evalue,
    match->qvalue,
    match->overlap,
    match->offset,
    match->orientation,
    logo_filename, logo_filename,                // links
    logo_filename                                 // image
  );

} // print_tomtom_html_line


/*****************************************************************************
 * prints the summary for the tomtom html output
 * 
 * DEPRECATED - the html is now generated by a stylesheet from the xml.
 *****************************************************************************/
void print_tomtom_html_tailer(
  FILE *output_path,
  int argc,
  char **argv,
  int total_matches,
  BOOLEAN_T sig_type_q,
  float sig_thresh
) 
{
  int i;

  fprintf(output_path, 
    "  <tr>\n"
    "    <td> <b>Total matches with <i>%s</i>-value &le; %g:</b> %d\n"
    "    </td>\n"
    "  </tr>\n", 
    (sig_type_q ? "q" : "E"),
    sig_thresh,
    total_matches);
  fprintf(output_path, "<tr><th>Command line:</th><td>\n");
  for (i=0; i<argc; i++) {
    fprintf(output_path, "%s\n", argv[i]);
  }
  fprintf(output_path, "</td></tr>\n");
  fprintf(output_path,
    "</table>\n\n"
    "</body>\n"
    "</html>\n"
  );
} // print_tomtom_html_tailer


/*****************************************************************************
 * outputs the xml for a motif.
 *
 * Rather than using the motif id embeded in the motif object a passed string is
 * printed. This is because the motif may have been reverse complemented and the
 * strand printed on the front of the id (and the caller may not want that in 
 * the file). As query and target motifs may have conflicting ids a prefix
 * is appended to ensure uniqueness.
 *
 * The indent and tab strings specifiy respectivly, the text string appended to the
 * front of each line, and the text string appended in front of nested tags.
 * 
 * the buffer must be expandable using realloc and the buffer_len must be set
 * to reflect the currently allocated size or zero if the buffer is unallocated
 * and set to null.
 *****************************************************************************/
static void print_xml_motif(FILE *file, char *prefix, char *name, MOTIF_T *motif, char *indent, char *tab, char **buffer, int *buffer_len) {
  int i, len;
  char *alt, *url;
  MATRIX_T *freqs;
  double A, C, G, T;

  alt = get_motif_id2(motif);
  len = get_motif_length(motif);
  url = get_motif_url(motif);
  freqs = get_motif_freqs(motif);

  //put the name in the buffer (as we're about to use it twice)
  replace_xml_chars2(name, buffer, buffer_len, 0, TRUE);
  fprintf(file, "%s<motif ", indent);
  fprintf(file, "id=\"%s%s\" name=\"%s\" ", prefix, *buffer, *buffer);
  if (alt && alt[0] != '\0') fprintf(file, "alt=\"%s\" ", replace_xml_chars2(alt, buffer, buffer_len, 0,  TRUE));
  fprintf(file, "length=\"%d\" nsites=\"%g\" evalue=\"%g\" ", len, get_motif_nsites(motif), get_motif_evalue(motif));
  if (url && url[0] != '\0') fprintf(file, "url=\"%s\"", replace_xml_chars2(url, buffer, buffer_len, 0,  TRUE));
  fprintf(file, ">\n");
  for (i = 0; i < len; ++i) {
    //TODO check that the order is always ACGT
    A = get_matrix_cell(i, 0, freqs);
    C = get_matrix_cell(i, 1, freqs);
    G = get_matrix_cell(i, 2, freqs);
    T = get_matrix_cell(i, 3, freqs);
    
    fprintf(file, "%s%s<pos i=\"%d\" A=\"%g\" C=\"%g\" G=\"%g\" T=\"%g\"/>\n", indent, tab, (i+1), A, C, G, T);
  }
  fprintf(file, "%s</motif>\n", indent);
}

/*****************************************************************************
 * outputs the xml for the start of a query including the query motif.
 * 
 * the buffer must be expandable using realloc and the buffer_len must be set
 * to reflect the currently allocated size or zero if the buffer is unallocated
 * and set to null.
 *****************************************************************************/
static void print_xml_query_start(FILE *file, MOTIF_T *motif, char **buffer, int *buffer_len) {
  fprintf(file, "\t\t\t<query>\n");
  print_xml_motif(file, "q_", get_motif_id(motif), motif, "\t\t\t\t", "\t", buffer, buffer_len);
}

/*****************************************************************************
 * outputs a tomtom match for the xml
 * 
 * the buffer must be expandable using realloc and the buffer_len must be
 * set to the correct value (the buffer is allowed to be unallocated initially but
 * in that case the buffer_len must be zero and buffer must be NULL)
 *****************************************************************************/
static void print_xml_query_match(FILE *file, TOMTOM_MATCH_T *match, char **buffer, int *buffer_len) {
  char *target_id;
  int db_id;

  target_id = get_motif_id(match->target)+1;
  db_id = match->db->id;
  fprintf(file, "\t\t\t\t<match target=\"t_%d_%s\" orientation=\"%s\" offset=\"%d\" pvalue=\"%g\" evalue=\"%g\" qvalue=\"%g\"/>\n", 
      db_id, replace_xml_chars2(target_id, buffer, buffer_len, 0,  TRUE), (match->orientation == '+' ? "forward" : "reverse"), 
      match->offset, match->pvalue, match->evalue, match->qvalue);
}

/*****************************************************************************
 * outputs the xml tag designating the end of information on a query motif
 *****************************************************************************/
static void print_xml_query_end(FILE *file) {
  fprintf(file, "\t\t\t</query>\n");
}

/*****************************************************************************
 * outputs the xml results
 *
 * the buffer must be expandable using realloc and the buffer_len must be set
 * to reflect the currently allocated size or zero if the buffer is unallocated
 * and set to null.
 *****************************************************************************/
static void print_xml_results(FILE *xml_output, 
    int argc, char**argv, 
    char *distance_measure, 
    BOOLEAN_T sig_type_q, double sig_thresh, 
    char *bg_file, ARRAY_T *bg_freqs_target,
    time_t now,
    ARRAYLST_T *target_dbs,
    char *query_source, char *query_name, time_t query_last_mod, FILE *query_matches_temp,
    char **buffer, int *buffer_len) {

    int i;
    BOOLEAN_T bg_from_file;
    double A, C, G, T, cycles;

    //print dtd
    fprintf(xml_output, "%s", tomtom_dtd);
    //print xml body
    fprintf(xml_output, "<tomtom version=\"" VERSION "\" release=\"" ARCHIVE_DATE "\">\n");

    fprintf(xml_output, "\t<model>\n");
    fprintf(xml_output, "\t\t<command_line>tomtom");
    for (i = 1; i < argc; ++i) {
      char *arg;
      //note that this doesn't correctly handle the case of a " character in a filename
      //in which case the correct output would be an escaped quote or \"
      //but I don't think that really matters since you could guess the original command
      arg = replace_xml_chars2(argv[i], buffer, buffer_len, 0,  TRUE);
      if (strchr(argv[i], ' ')) {
        fprintf(xml_output, " &quot;%s&quot;", arg);
      } else {
        fprintf(xml_output, " %s", arg);
      }
    }
    fprintf(xml_output, "</command_line>\n");
    fprintf(xml_output, "\t\t<distance_measure value=\"%s\"/>\n", distance_measure);
    fprintf(xml_output, "\t\t<threshold type=\"%s\">%g</threshold>\n", (sig_type_q ? "qvalue" : "evalue"), sig_thresh);

    bg_from_file = (strcmp(bg_file, "motif-file") != 0);
    //TODO check ordering is always ACGT
    A = get_array_item(0, bg_freqs_target);
    C = get_array_item(1, bg_freqs_target);
    G = get_array_item(2, bg_freqs_target);
    T = get_array_item(3, bg_freqs_target);


    fprintf(xml_output, "\t\t<background from=\"%s\" A=\"%g\" C=\"%g\" G=\"%g\" T=\"%g\"", 
        (bg_from_file ? "file" : "first_target"), A, C, G, T);

    if (bg_from_file) {
      fprintf(xml_output, " file=\"%s\"/>\n", replace_xml_chars2(bg_file, buffer, buffer_len, 0,  TRUE));
    } else {
      fprintf(xml_output, "/>\n");
    }

    fprintf(xml_output, "\t\t<host>%s</host>\n", hostname());
    fprintf(xml_output, "\t\t<when>%s</when>\n", strtok(ctime(&now),"\n"));
    fprintf(xml_output, "\t</model>\n");
    fprintf(xml_output, "\t<targets>\n");
    for (i = 0; i < arraylst_size(target_dbs); ++i) {
      MOTIF_DB_T *db;
      RBNODE_T *node;
      int db_id_len, prefix_len, num;
      char *prefix;

      db = arraylst_get(i, target_dbs);

      db_id_len = 1;
      num = db->id;
      if (num < 0) {
        num = -num;
        db_id_len++;
      }
      while (num >= 10) {
        db_id_len++;
        num /= 10;
      }

      prefix_len = 2 + db_id_len + 2;
      prefix = (char*)mm_malloc(sizeof(char) * prefix_len);
      snprintf(prefix, prefix_len, "t_%d_", db->id);

      fprintf(xml_output, "\t\t<target_file index=\"%d\" ", db->id);
      fprintf(xml_output, "source=\"%s\" ", replace_xml_chars2(db->source, buffer, buffer_len, 0,  TRUE)); 
      fprintf(xml_output, "name=\"%s\" ", replace_xml_chars2(db->name, buffer, buffer_len, 0,  TRUE)); 
      fprintf(xml_output, "loaded=\"%d\" ", db->loaded);
      fprintf(xml_output, "excluded=\"%d\" ", db->excluded);
      fprintf(xml_output, "last_mod_date=\"%s\">\n", strtok(ctime(&(db->last_mod)),"\n"));
      for (node = rbtree_first(db->matched_motifs); node != NULL; node = rbtree_next(node)) {
        MOTIF_T *target;
        target = (MOTIF_T*)rbtree_value(node);
        print_xml_motif(xml_output, prefix, get_motif_id(target)+1, target, "\t\t\t", "\t", buffer, buffer_len);
      }
      fprintf(xml_output, "\t\t</target_file>\n");
    }
    fprintf(xml_output, "\t</targets>\n");
    fprintf(xml_output, "\t<queries>\n");
    fprintf(xml_output, "\t\t<query_file ");
    fprintf(xml_output, "source=\"%s\" ", replace_xml_chars2(query_source, buffer, buffer_len, 0,  TRUE)); 
    fprintf(xml_output, "name=\"%s\" ", replace_xml_chars2(query_name, buffer, buffer_len, 0,  TRUE)); 
    fprintf(xml_output, "last_mod_date=\"%s\">\n", strtok(ctime(&query_last_mod),"\n"));
    // Copy the matches to the output
    {
      int c;
      rewind(query_matches_temp);
      while ((c = fgetc(query_matches_temp)) != EOF) {
        c = fputc(c, xml_output);
        if (c == EOF) {
          break; // Error writing byte
        }
      }
    }
    fprintf(xml_output, "\t\t</query_file>\n");
    fprintf(xml_output, "\t</queries>\n");
    
    cycles = myclock(); //measures number of cycles since first call to myclock (unfortunately prone to overflow for large run times).
    fprintf(xml_output, "\t<runtime cycles=\"%.0f\" seconds=\"%.3f\"/>\n", cycles, cycles / CLOCKS_PER_SEC);
    fprintf(xml_output, "</tomtom>\n");

}

/*****************************************************************************
 * Parses the full string str to produce a float
 * if the complete string is not parsed or an error occurs it prints the 
 * appropriate error passing it the string (and remainder string for partial) and exits.
 *****************************************************************************/
static double parse_float(char *str, char *error_none, char *error_partial, char *error_range) {
  char *end_ptr;
  double result;
  result = strtod(str, &end_ptr);
  if (end_ptr == str) {
    if (error_none) die(error_none, str);
  }
  if (*end_ptr != '\0') {
    if (error_partial) die(error_partial, str, end_ptr);
  }
  if (errno == ERANGE) {
    if (error_range) die(error_range, str);
  }
  return result;
}

/*****************************************************************************
 * outputs the file name when given a path to that file.
 *
 * searches the path string for '/' chars and returns the pointer to the position
 * after the last '/'. Note that the path is not copied or modified in any way.
 * This is not intended to be perfect as I'm sure there's some valid way to have 
 * a '/' in a file name but I expect it will work most of the time.
 *****************************************************************************/
static char* file_name_from_path(char *path) {
  char *where, *result;
  result = path;
  for (where = path; *where != '\0'; where++) {
    if (*where == '/') result = where+1;
  }
  return result;
}

/*****************************************************************************
 * MAIN
 *****************************************************************************/
VERBOSE_T verbosity = NORMAL_VERBOSE;
#ifdef MAIN

int main(int argc, char *argv[])
{
  char *default_output_dirname = "tomtom_out"; // where to write HTML
  FILE *text_output;                        // where to write legacy output
  //deprecated - html output is now done with xml stylesheets
  //FILE *old_html_output = NULL;
  FILE *match_output = NULL;

  char *xml_convert_buffer = NULL;
  int xml_convert_buffer_len = 0;

  time_t now;

  /**********************************************
   * Initialize command line parameters
   **********************************************/
  char* bg_file = "motif-file"; // default to reading background freq. from first motif file 
  char* query_filename = NULL;
  char* query_outname = NULL;
  time_t query_last_mod;

  RBTREE_T *filter_ids = rbtree_create(compare_strings, NULL, NULL, NULL, NULL);
  RBTREE_T *filter_indexes = rbtree_create(compare_ints, copy_int, free, NULL, NULL);
  ARRAYLST_T *target_dbs = NULL;
  char *output_dirname = default_output_dirname;
  BOOLEAN_T clobber = TRUE;                    // clobber the default output dir
  BOOLEAN_T text_only = FALSE;                 // default is HTML
  BOOLEAN_T html = TRUE;                       // output html from xml (maybe later I'll add an option to turn it off)
  char* distance_measure = "pearson";
  void (*distance_func)(MOTIF_T *, MOTIF_T *, ARRAY_T *, MATRIX_T **) = pearson_scores;
  BOOLEAN_T sig_type_q = TRUE;
  float sig_thresh = 0.5;
  float query_pseudocount = DEFAULT_QUERY_PSEUDO; 
  float target_pseudocount = DEFAULT_TARGET_PSEUDO; 
  char* columnwise_scores_filename = NULL;
  BOOLEAN_T ssc = TRUE;                        // small sample correction in logos

  /**********************************************
   * Command line parsing.
   **********************************************/
  int option_count = 17;     // Must be updated if list below is
  cmdoption const motif_scan_options[] = {
    {"o", REQUIRED_VALUE},
    {"oc", REQUIRED_VALUE},
    {"bfile", REQUIRED_VALUE},
    {"m", REQUIRED_VALUE},
    {"mi", REQUIRED_VALUE},
    {"text", NO_VALUE},
    {"thresh", REQUIRED_VALUE},
    {"q-thresh", REQUIRED_VALUE},
    {"evalue", NO_VALUE},
    {"dist", REQUIRED_VALUE},
    {"internal", NO_VALUE},
    {"min-overlap", REQUIRED_VALUE},
    {"query-pseudo", REQUIRED_VALUE},
    {"target-pseudo", REQUIRED_VALUE},
    {"column-scores", REQUIRED_VALUE},//developer only
    {"no-ssc", NO_VALUE},
    {"verbosity", REQUIRED_VALUE},
  };

  int option_index = 0;

  /**********************************************
   * Define the usage message.
   **********************************************/
  char * usage = "\n   USAGE: tomtom [options] <query file> <target file>+\n"
  "\n"
  "   Options:\n"
  "     -o <output dir>\t\tname of directory for output files;\n\t\t\t\twill not replace existing directory\n"
  "     -oc <output dir>\t\tname of directory for output files;\n\t\t\t\twill replace existing directory\n"
  "     -bfile <background file>\tname of background file\n"
  "     -m <id>\tuse only query motifs with a specified id; may be repeated\n"
  "     -mi <index>\tuse only query motifs with a specifed index; may be repeated\n"
  "     -thresh <float>\t\tsignificance threshold; default: 0.5\n"
  "     -evalue\t\tuse E-value threshold; default: q-value\n"
  "     -dist allr|ed|kullback|pearson|sandelin|blic1|blic5|llr1|llr5\n"
  "\t\t\t\tdistance metric for scoring alignments;\n\t\t\t\tdefault: pearson\n"
  "     -internal\t\t\tonly allow internal alignments;\n\t\t\t\tdefault: allow overhangs\n"
  "     -min-overlap <int>\t\tminimum overlap between query and target;\n\t\t\t\tdefault: 1\n"
  "     -query-pseudo <float>\tdefault: %.1f\n"
  "     -target-pseudo <float>\tdefault: %.1f\n"
  "     -text\t\t\toutput in text format (default is HTML)\n"
  "     -no-ssc\t\t\tdon't apply small-sample correction to logos\n\t\t\t\tdefault: use small-sample correction\n"
  "     -verbosity [1|2|3|4]\tdefault: %d\n"
  // The "--column-scores <file>" option is undocumented -- developers only.
  "\n";

  (void) myclock();                     /* record CPU time */
  now = time(NULL);
  /**********************************************
   * Parse the command line.
   **********************************************/
  if (simple_setopt(argc, argv, option_count, motif_scan_options) != NO_ERROR) {
    die("Error: option name too long.\n");
  }
  BOOLEAN_T internal = FALSE;
  int min_overlap = 1;
  while (TRUE) {
    int c = 0;
    char* option_name = NULL;
    char* option_value = NULL;
    const char * message = NULL;

    /*  Read the next option, and break if we're done. */
    c = simple_getopt(&option_name, &option_value, &option_index);
    if (c == 0) {
      break;
    } else if (c < 0) {
      (void) simple_getopterror(&message);
      die("Error processing command line options (%s)\n", message);
    }

    if (strcmp(option_name, "thresh") == 0) {
      sig_thresh = parse_float(option_value, "Couldn't interpret value \"%s\" for -thresh argument.\n", 
          "Only partially interpreted value \"%s\" for -thresh argument as it contained uninterpretable text \"%s\".\n", 
          NULL);
    } else if (strcmp(option_name, "q-thresh") == 0) {        // kept for backward compatibility
      sig_thresh = parse_float(option_value, "Couldn't interpret value \"%s\" for -q-thresh argument.\n", 
          "Only partially interpreted value \"%s\" for -q-thresh argument as it contained uninterpretable text \"%s\".\n", 
          NULL);
      if (verbosity >= NORMAL_VERBOSE) {
        fprintf(stderr, "Warning: -thresh should be used in preference to -q-thresh as -q-thresh is deprecated.\n");
      }
    } else if (strcmp(option_name, "evalue") == 0) {
      sig_type_q = FALSE;
    } else if (strcmp(option_name, "bfile") == 0) {
      bg_file = option_value;
      if (!file_exists(bg_file)) {
        die("Can't find the specified background file \"%s\".", bg_file);
      }
    } else if (strcmp(option_name, "o") == 0) {
      output_dirname = option_value;
      clobber = FALSE;
    } else if (strcmp(option_name, "oc") == 0) {
      output_dirname = option_value;
      clobber = TRUE;
    } else if (strcmp(option_name, "text") == 0) {
      text_only = TRUE;
    } else if (strcmp(option_name, "no-ssc") == 0) {
      ssc = FALSE;
    } else if (strcmp(option_name, "m") == 0) {
      rbtree_put(filter_ids, option_value, NULL); 
    } else if (strcmp(option_name, "mi") == 0) {
      int num = atoi(option_value);
      rbtree_put(filter_indexes, &num, NULL); 
    } else if (strcmp(option_name, "dist") == 0) {
      distance_measure = option_value;
      //determine the selected distance function
      if (strcmp(distance_measure, "pearson") == 0) {
        distance_func = pearson_scores;
      } else if (strcmp(distance_measure, "allr") == 0) {
        distance_func = allr_scores;
      } else if (strcmp(distance_measure, "ed") == 0) {
        distance_func = ed_scores;
      } else if (strcmp(distance_measure, "kullback") == 0) {
        distance_func = kullback_scores;
      } else if (strcmp(distance_measure, "sandelin") == 0) {
        distance_func = sandelin_scores;
      } else if (strcmp(distance_measure, "blic1") == 0) {
        distance_func = blic1_scores;
      } else if (strcmp(distance_measure, "blic5") == 0) {
        distance_func = blic5_scores;
      } else if (strcmp(distance_measure,"llr1") == 0) {
        distance_func = llr1_scores;
      } else if (strcmp(distance_measure,"llr5") == 0) {
        distance_func = llr5_scores;
      } else {
        fprintf(stderr,
            "error: invalid distance measure \"%s\"\n-dist [allr|ed|kullback|pearson|sandelin|blic1|blic5|llr1|llr5]\n",
            distance_measure
            );
        exit(EXIT_FAILURE);
      }
    } else if (strcmp(option_name, "internal") == 0) {
      internal = TRUE;
    } else if (strcmp(option_name, "min-overlap") == 0) {
      //did some old version allow fractional overlaps? If not then why parse as a float and convert to integer?
      double temp_var;
      temp_var = atof(option_value);
      min_overlap = rint(temp_var);
    } else if (strcmp(option_name, "query-pseudo") == 0) {
      query_pseudocount = parse_float(option_value, 
          "Couldn't interpret value \"%s\" for -query-pseudo argument.\n", 
          "Only partially interpreted value \"%s\" for -query-pseudo argument as it contained uninterpretable text \"%s\".\n", 
          "The value \"%s\" for -query-pseudo is outside of the internally representable range.\n");
    } else if (strcmp(option_name, "target-pseudo") == 0) {
      target_pseudocount = parse_float(option_value, 
          "Couldn't interpret value \"%s\" for -target-pseudo argument.\n", 
          "Only partially interpreted value \"%s\" for -target-pseudo argument as it contained uninterpretable text \"%s\".\n", 
          "The value \"%s\" for -target-pseudo is outside of the internally representable range.\n");
    } else if (strcmp(option_name, "verbosity") == 0) {
      verbosity = atoi(option_value);
      if (verbosity < QUIET_VERBOSE || verbosity > DUMP_VERBOSE) {
        die("Unknown verbosity setting (%d).\n", verbosity);
      }
    } else if (strcmp(option_name, "column-scores") == 0) {
      columnwise_scores_filename = option_value;
    }
  }
  if (sig_type_q) {
    if (sig_thresh < 0 || sig_thresh > 1) die("Significance threshold (%g) must be within the inclusive range 0 to 1 if it is a q-value.\n", sig_thresh);
  } else {
    if (sig_thresh < 0) die("Significance threshold (%g) must be positive if it is an E-value.", sig_thresh);
  }
  if (argc < option_index + 2) { 
    fprintf(stderr, usage, DEFAULT_QUERY_PSEUDO, DEFAULT_QUERY_PSEUDO, NORMAL_VERBOSE);
    exit(EXIT_FAILURE);
  }
  // check the file names for the query and targets
  // allow one of them to use stdin
  {
    BOOLEAN_T stdin_used = FALSE;
    int target_db_id = 1;
    query_filename = argv[option_index++];
    query_last_mod = now;
    if (strcmp(query_filename, "-") == 0) {
      stdin_used = TRUE;
      query_outname = "standard input";
    } else if (file_exists(query_filename)) {
      struct stat stbuf;
      int result;
      query_outname = file_name_from_path(query_filename);
      //update the last mod date to the correct value
      result = stat(query_filename, &stbuf);
      if (result < 0) {
        //an error occured
        die("Failed to stat file \"%s\", error given as: %s\n", query_filename, strerror(errno));
      }
      query_last_mod = stbuf.st_mtime;
    } else {
      die("Query file \"%s\" does not exist.\n", query_filename);
    }
    target_dbs = arraylst_create();
    while (option_index < argc) { 
      MOTIF_DB_T *db;
      char *target_file, *target_name;
      target_file = argv[option_index++];
      target_name = NULL;
      if (strcmp(target_file, "-") == 0) {
        if (stdin_used) {
          die("Target file %d made use of standard input when it is already in use for another file.\n", arraylst_size(target_dbs)+1);
        }
        stdin_used = TRUE;
        target_name = "standard input";
      } else if (file_exists(target_file)) {
        target_name = file_name_from_path(target_file);
      } else {
        die("Target file \"%s\" does not exist.\n", target_file);
        target_name = NULL;//make the complier happy
      }
      db = create_motif_db(target_db_id++, target_file, target_name, now);
      arraylst_add(db, target_dbs);
    }
  }

  if (text_only == TRUE) {
    // Legacy: plain text output to standard out.
    text_output = stdout;
  } else {
    if (create_output_directory(output_dirname, clobber, (verbosity >= NORMAL_VERBOSE))) {
      // Failed to create output directory.
      die("Unable to create output directory %s.\n", output_dirname);
    }
    // Create the name of the output files (text, XML, HTML) and open them
    char *path;
    path = make_path_to_file(output_dirname, TXT_FILENAME);
    text_output = fopen(path, "w"); //FIXME CEG check for errors
    myfree(path);
    //deprecated - html output is now done with xml stylesheets
    //path = make_path_to_file(output_dirname, OLD_HTML_FILENAME);
    //old_html_output = fopen(path, "w"); //FIXME CEG check for errors
    //myfree(path);
    match_output = tmpfile();
    xml_convert_buffer_len = 100;//initial size, this may be expanded
    xml_convert_buffer = mm_malloc(sizeof(char) * xml_convert_buffer_len);
  }

  /**********************************************
   * Read the motifs.
   **********************************************/

  // General variables for reading the meme file
  BOOLEAN_T has_reverse_strand_query = FALSE;
  STRING_LIST_T* motif_occurrences = NULL;
  ARRAY_T* bg_freqs_query;

  //load all the query motifs
  ARRAYLST_T *query_motifs = arraylst_create();
  read_meme_file2(
      query_filename,
      bg_file,
      query_pseudocount,
      REQUIRE_PSPM,
      query_motifs,
      NULL, //motif occurrances (not used)
      &has_reverse_strand_query,
      &bg_freqs_query
      );
  //filter the query motifs
  if (rbtree_size(filter_indexes) > 0 || rbtree_size(filter_ids) > 0) {
    int last_index, SNUM = 0;//allows removing range down to first
    RBNODE_T *node;
    last_index = arraylst_size(query_motifs);
    rbtree_put(filter_indexes, &SNUM, NULL);
    node = rbtree_last(filter_indexes);
    for (; node != NULL; node = rbtree_prev(node)) {
      int num, i;
      //note num is protected location index plus 1
      num = *((int*)rbtree_key(node));
      if (num < 0) break;//check for bad numbers (note zero is accepted despite the first num being 1 
                            //because we want to remove the range down to index zero)
      for (i = last_index - 1; i >= num; --i) {
        MOTIF_T *motif;
        motif = (MOTIF_T*)arraylst_get(i, query_motifs);
        if (rbtree_lookup(filter_ids, get_motif_id(motif), FALSE, NULL) == NULL && 
            rbtree_lookup(filter_ids, get_motif_id2(motif), FALSE, NULL) == NULL) {
          arraylst_remove(i, query_motifs);
          free_motif(motif);
        }
      }
      last_index = i;
    }
  }
  //no need to maintain the filter sets now
  rbtree_destroy(filter_indexes);
  rbtree_destroy(filter_ids);

  ARRAY_T* bg_freqs_target = NULL;
  BOOLEAN_T has_reverse_strand_target = FALSE;

  //load all the target motifs
  int i, j;
  ARRAYLST_T *target_motifs = arraylst_create();
  {
    BOOLEAN_T created, printed_warning;
    MOTIF_T *motif;
    ARRAY_T *background;
    BOOLEAN_T has_reverse_strand;
    //for each target database
    for (i = 0, j = 0; i < arraylst_size(target_dbs); ++i) {
      RBTREE_T *ids;
      MOTIF_DB_T *db;
      //for duplicate checking within the database (other databases may share the same identifiers)
      ids = rbtree_create(compare_strings, NULL, NULL, NULL, NULL);
      
      db = (MOTIF_DB_T*)arraylst_get(i, target_dbs);
      //record the starting index of the motifs for this database
      db->list_index = j;
      background = NULL;
      has_reverse_strand = FALSE;
      printed_warning = FALSE;

      read_meme_file2(
          db->source,
          bg_file, // background freq. filename
          target_pseudocount,
          REQUIRE_PSPM,
          target_motifs,
          NULL, //motif occurrances (not used)
          &has_reverse_strand,
          &background
          );

      //record the number of motifs in this database
      db->loaded = arraylst_size(target_motifs) - db->list_index;

      //check that the background is the same
      if (bg_freqs_target) {
        //compare the background with the reference to ensure they are not different types, the bound is fairly lax
        if (has_reverse_strand != has_reverse_strand_target || !equal_arrays(0.25, bg_freqs_target, background)) {
          int k;
          //everything has been excluded from this database
          db->excluded = db->loaded;
          //remove the added motifs
          for (k = arraylst_size(target_motifs)-1; k >= j; --k) {
            MOTIF_T *a_motif = (MOTIF_T*)arraylst_remove(k, target_motifs);
            free_motif(a_motif);
            free(a_motif);
          }
          if (has_reverse_strand != has_reverse_strand_target) {
            if (has_reverse_strand) {
              fprintf(stderr, "Warning, the target file \"%s\" had a reverse strand unlike the reference and so was removed.\n", 
                  db->name);
            } else {
              fprintf(stderr, "Warning, the target file \"%s\" doesn't have a reverse strand unlike the reference and so was removed.\n", 
                  db->name);
            }
          } else {
            fprintf(stderr, "Warning, the target file \"%s\" had a different background and so was removed.\n", 
                db->name);
          }
        }
        free_array(background); 
      } else {
        bg_freqs_target = background;
        has_reverse_strand_target = has_reverse_strand;
      }
      //check for duplicate identifiers in the same file
      for (; j < arraylst_size(target_motifs); ++j) {
        motif = (MOTIF_T*)arraylst_get(j, target_motifs);
        rbtree_lookup(ids, get_motif_id(motif), TRUE, &created);
        if (!created) {
          MOTIF_T *a_motif = (MOTIF_T*)arraylst_remove(j--, target_motifs);
          if (!printed_warning) {
            if (verbosity >= NORMAL_VERBOSE) fprintf(stderr, "Warning, the following duplicate target motif ids in \"%s\" were removed:", 
                db->name);
            printed_warning = TRUE;
          }
          if (verbosity >= NORMAL_VERBOSE) fprintf(stderr, " %s", get_motif_id(motif));
          free_motif(a_motif);
          free(a_motif);
          db->excluded += 1;
        }
      }
      if (verbosity >= NORMAL_VERBOSE) {
        if (printed_warning) fprintf(stderr, "\n");
      } else if (verbosity <= QUIET_VERBOSE) {
        if (db->excluded) fprintf(stderr, "Warning, %d target motif(s) with duplicate ids in \"%s\" were removed.\n", db->excluded, db->name);
      }
      db->list_entries = j - db->list_index;
      //destroy the set of ids for this database
      rbtree_destroy(ids);
    }//for each target db
  }//end loading motifs

  /**********************************************
   * Reverse complement the motifs
   **********************************************/

  //put the existing motifs at the even indexes
  //and the reverse complement motifs at the odd indexes after the source motif
  add_reverse_complements2(target_motifs);

  //now we have double the motifs all the sizes in the motif dbs have to be doubled
  for (i = 0, j = 0; i < arraylst_size(target_dbs); ++i) {
    MOTIF_DB_T *db = (MOTIF_DB_T*)arraylst_get(i, target_dbs);
    db->list_index *= 2;
    db->list_entries *= 2;
  }


  /**********************************************
   * Check for a reasonable number of motifs to search
   **********************************************/

  int query_count = arraylst_size(query_motifs);
  int target_count = arraylst_size(target_motifs);

  // Check for small target databases (multiply by 2 for reverse complements).
  if (target_count < (MIN_TARGET_DATABASE_SIZE * 2)) {
    fprintf(stderr,
            "\nWarning: Target database size too small (%d). Provide at least %d motifs for accurate p-value computation\n\n",
            target_count/2, MIN_TARGET_DATABASE_SIZE);
  }


  /**********************************************
   * Strip ambiguous letters from motifs
   **********************************************/

  int query_i;
  int target_i;
  int total_target_length = 0;

  for(target_i = 0; target_i < target_count; ++target_i) {
    MOTIF_T* target = arraylst_get(target_i, target_motifs);
    total_target_length += target->length;
    // Remove the ambiguous letter columns
    while(get_num_cols(target->freqs) > target->alph_size) {
      remove_matrix_col(target->alph_size, target->freqs);
    }
  }

  // Remove the ambiguous letter columns for query
  for (query_i = 0; query_i < query_count; ++query_i) {
    MOTIF_T* query = arraylst_get(query_i, query_motifs);
    while (get_num_cols(query->freqs) > query->alph_size) {
      remove_matrix_col(query->alph_size, query->freqs);
    }
  }

  // Uniform background distribution array for the infinite alphabet
  ARRAY_T* bg_freqs_mod = allocate_array(total_target_length);
  init_array((ATYPE)(1.0 / total_target_length), bg_freqs_mod);

  // Recurse over the set of query motifs and target motifs
  int total_matches = 0;

  // Print output
  fprintf(text_output, "#Query ID\tTarget ID\tOptimal offset");
  fprintf(text_output, "\tp-value\tE-value\tq-value\tOverlap\t");
  fprintf(text_output, "Query consensus\tTarget consensus\tOrientation\n");
  //deprecated - html output is now done with xml stylesheets
  //if (! text_only) { 
  //  print_tomtom_html_header(old_html_output, query_filename, target_dbs, distance_measure, sig_type_q, sig_thresh); 
  //}

  /**********************************************
   * The main loop
   **********************************************/

  // Recurse over the set of queries. 
  for (query_i = 0; query_i < query_count; query_i++) {
    OBJECT_LIST_T* output_list = new_object_list(NULL, NULL, NULL, free);

    if (verbosity >= NORMAL_VERBOSE) fprintf(stderr, "Processing query %d out of %d \n", query_i + 1, query_count);
    int target_database_length = total_target_length;

    MOTIF_T *query_motif = (MOTIF_T*)arraylst_get(query_i, query_motifs); 
    MOTIF_T *target_motif, *original_motif;

    // Initialize an array to store the columnwise scores.
    MATRIX_T* all_columnwise_scores;
    all_columnwise_scores = 
      allocate_matrix(
        query_motif->length,
        target_database_length
        );

    int col_offset = 0;

    /**********************************************
     * Apply the distance function to each target
     **********************************************/

    // Recurse over the set of targets.
    for (target_i = 0; target_i < target_count; target_i++) {
      target_motif = (MOTIF_T*)arraylst_get(target_i, target_motifs);

        // It would be nice to allocate this outside the loop
        // but alas the target motifs may all be different sizes
        MATRIX_T* pairwise_column_scores;
        pairwise_column_scores = 
            allocate_matrix(
                query_motif->length,
                target_motif->length
                );

        // Compute score matrix based on the user choosen distance measure
        distance_func(query_motif, target_motif, bg_freqs_target, &pairwise_column_scores);

        // Update the background frequency array
        // tlb added documentation; 
        // the "frequency array" is [w_q, \sum w_t] and contains
        // scores of query motif columns (the rows) vs all target columns
        // with each new target's scores appended as a new set of columns
        // Each row of the frequency array contains the scores of a single
        // query column vs. all possible target columns.
        int row_i;
        for(row_i = 0; row_i < pairwise_column_scores->num_rows; row_i++) {
          int col_i;
          for(col_i = 0; col_i < pairwise_column_scores->num_cols; col_i++) {
            set_matrix_cell(
                row_i,
                col_i + col_offset,
                get_matrix_cell(row_i, col_i, pairwise_column_scores),
                all_columnwise_scores
                );
          }
        }

        col_offset += get_num_cols(pairwise_column_scores);

        free_matrix(pairwise_column_scores);

    }// Parsing all targets for the query done.

    // If requested, store the columnwise scores in an external file.
    if (columnwise_scores_filename != NULL) {
      FILE* columnwise_scores_file;
      fprintf(stderr, 
          "Storing %d by %d matrix of column scores in %s.\n", 
          get_num_rows(all_columnwise_scores),
          get_num_cols(all_columnwise_scores),
          columnwise_scores_filename);
      if (!open_file(columnwise_scores_filename, "w", FALSE, "column scores",
          "column scores", &columnwise_scores_file)) {
        exit(1);
      }
      print_matrix(all_columnwise_scores, 10, 7, FALSE, 
          columnwise_scores_file);
      fclose(columnwise_scores_file);
    }

    // Scale the pssm
    // this scales each score to a value in the range 0 to BINS
    PSSM_T* pssm = build_matrix_pssm(
        all_columnwise_scores,  // matrix
        NULL,                   // don't compute pv lookup table
        BINS);                  // range

    // TODO Bad--exposing internals of PSSM object.
    // Get the scaled pssm matrix and throw away the pssm object.
    SWAP (MATRIX_T *, pssm->matrix, all_columnwise_scores);
    free_pssm(pssm);

    // Compute the pvalue distributions
    int num_pv_lookup_array; // w * (w+1) / 2
    num_pv_lookup_array = (all_columnwise_scores->num_rows * 
        (all_columnwise_scores->num_rows + 1))/2;

    MATRIX_T* reference_matrix;
    reference_matrix = allocate_matrix(
        all_columnwise_scores->num_rows,
        all_columnwise_scores->num_rows);
    init_matrix(-1, reference_matrix);

    int num_array_lookups = (BINS * query_motif->length) + 1;
    MATRIX_T *pv_lookup_matrix = allocate_matrix(
        num_pv_lookup_array,
        num_array_lookups);
    init_matrix(0, pv_lookup_matrix);

    get_pv_lookup_new(
        all_columnwise_scores,  // scores of each query column
        target_database_length, // number of scores 
        BINS,                   // maximum score in a column
        bg_freqs_mod,           // Uniform
        reference_matrix,
        pv_lookup_matrix);

    /**********************************************
     * Extract the scaled scores for each target
     * and compute the optimal offset for the
     * best p-value
     **********************************************/

    // Compute the motif scores
    MOTIF_DB_T *target_db;
    int db_remain, db_i;
    db_i = 0;
    target_db = arraylst_get(db_i, target_dbs);
    db_remain = target_db->list_entries;
    for (target_i = 0, col_offset = 0; target_i < target_count; 
        col_offset += target_motif->length, ++target_i, --db_remain) {
      while (db_remain <= 0) {
        target_db = arraylst_get(++db_i, target_dbs);
        db_remain = target_db->list_entries;
      }

      target_motif = (MOTIF_T*)arraylst_get(target_i, target_motifs);

      //use the knowlege that the original motif is before its reverse complment
      //and hence is at all the even positions
      if (target_i % 2 == 0) {
        original_motif = target_motif;
      } else {
        original_motif = (MOTIF_T*)arraylst_get(target_i-1, target_motifs);
      }

      MATRIX_T* pairwise_subset;
      pairwise_subset = allocate_matrix(
          query_motif->length,
          target_motif->length);

      // Get a query-target pair from the overall score matrix.
      int col_i;
      for (col_i = 0; col_i < target_motif->length; col_i++) {
        ARRAY_T* temp_array = NULL;
        temp_array = get_matrix_column(col_i + col_offset,
                                       all_columnwise_scores);
        set_matrix_column(temp_array, col_i, pairwise_subset);
        free_array(temp_array);
      }

      // Compute the score and offset wrt the minimum pvalue CORRECTED.
      int optimal_offset = -1;
      double optimal_pvalue = 1;
      int overlapping_bases = -1;

      compare_motifs(
          pairwise_subset,
          &optimal_offset,
          reference_matrix,    // index into pv_lookup_matrix
          pv_lookup_matrix,    //   based on start, stop columns
          &optimal_pvalue,     //   of query motif
          &overlapping_bases,
          internal,
          min_overlap);

      // Override min-overlap if query is smaller than min-overlap
      int mo;
      mo = min_overlap;
      if (pairwise_subset->num_rows < mo) {
        mo = pairwise_subset->num_rows;
      }
      if (pairwise_subset->num_cols < mo) {
        mo = pairwise_subset->num_cols;
      }
      /* Correct the motif pvalue (minimization across several pvalues).
         Configurations mulitiplied by 2, because later the minimum of
         the two strands will be taken LATER. Both stands have the same
         number of total_configurations, thus muliplication by 2.*/
      int total_configurations = pairwise_subset->num_rows
        + pairwise_subset->num_cols - (2 * mo) + 1;

      /* Correct the total number of configurations for internal alignments */
      if (internal) {
        total_configurations = abs(pairwise_subset->num_rows
                                   - pairwise_subset->num_cols) + 1;
      }
      double motif_pvalue = EV(optimal_pvalue, total_configurations * 2);

      // Compute E-value.
      double motif_evalue = motif_pvalue * target_count / 2.0;

      /* Offset_by_length */
      int offset_by_length = optimal_offset - query_motif->length;

      // Store this match.
      TOMTOM_MATCH_T* my_match = new_tomtom_match(
          query_motif,
          target_motif,
          original_motif,
          target_db,
          offset_by_length,
          overlapping_bases,
          motif_pvalue,
          motif_evalue,
          offset_by_length);

      // Store this match.  Don't bother with the key.
      store_object((void*)my_match, NULL, motif_pvalue, output_list);

      free_matrix(pairwise_subset);
    } // Scores computed for all targets.

    // Free all the matrices/arrays related to the current query
    free_matrix(reference_matrix);
    free_matrix(pv_lookup_matrix);
    free_matrix(all_columnwise_scores);

    // Sort the stringlist.
    sort_objects(output_list);

    // Compute q-values.
    convert_tomtom_p_to_q(output_list);

    // Do some re-formatting of the output strings.
    OBJECT_LIST_T* output_final = new_object_list(NULL, NULL, NULL, free);
    select_strand(output_list, output_final);

    //output start of query
    if (!text_only) print_xml_query_start(match_output, query_motif, &xml_convert_buffer, &xml_convert_buffer_len);

    // Traverse the list of matches.
    char *query_consensus = get_cons(get_motif_freqs(query_motif));
    TOMTOM_MATCH_T* my_match = (TOMTOM_MATCH_T*)retrieve_next_object(output_final);
    while (my_match != NULL) {

      // Check significance threshold.
      if ((sig_type_q && my_match->qvalue <= sig_thresh) || 
          (!sig_type_q && my_match->evalue <= sig_thresh)) {
        MOTIF_DB_T *target_db;
        char *target_consensus;

        target_motif = my_match->target;
        original_motif = my_match->original;
        target_db = my_match->db;

        total_matches++;

        //put the motif into the tree so it can be printed later
        rbtree_put(target_db->matched_motifs, get_motif_id(original_motif)+1, original_motif);

        target_consensus = get_cons(get_motif_freqs(target_motif));

        // Print text output.
        fprintf(text_output,
            "%s\t%s\t%d\t%g\t%g\t%g\t%d\t%s\t%s\t%c\n",
            get_motif_id(query_motif),
            get_motif_id(target_motif)+1,//skip over the + or -
            my_match->offset,
            my_match->pvalue,
            my_match->evalue,
            my_match->qvalue,
            my_match->overlap,
            query_consensus,
            target_consensus,
            my_match->orientation
            );

        //free text only extra variables
        free(target_consensus);

        // Print XML output.
        if (!text_only) {
          // create alignment logo
          char *logo_filename 
            = create_logo(target_db, target_motif, query_motif, my_match->offset, ssc, output_dirname);

          print_xml_query_match(match_output, my_match, &xml_convert_buffer, &xml_convert_buffer_len);

          //deprecated - html output is now done with xml stylesheets
          //print_tomtom_html_line(old_html_output, my_match, logo_filename);

          myfree(logo_filename);
        }
      }

      // Get the next match.
      my_match = (TOMTOM_MATCH_T*)retrieve_next_object(output_final);
    } // loop over the list of matches

    if (!text_only) print_xml_query_end(match_output);

    free(query_consensus);
    free_object_list(output_final);
    free_object_list(output_list);

  }/*  The loop across the queries. */

  if (! text_only) {
    char *xml_path;
    FILE *xml_output;
    // print XML
    xml_path = make_path_to_file(output_dirname, XML_FILENAME);
    xml_output = fopen(xml_path, "w"); //FIXME CEG check for errors
    print_xml_results(
      xml_output,
      argc, argv,
      distance_measure, sig_type_q, sig_thresh,
      bg_file, bg_freqs_target,
      now,
      target_dbs,
      query_filename, query_outname, query_last_mod, match_output,
      &xml_convert_buffer,
      &xml_convert_buffer_len
    );
    free(xml_convert_buffer);
    // close matches temporary file
    fclose(match_output);
    // close xml output
    fclose(xml_output);
    //generate html from xml
    if (html) {
      char *stylesheet_path, *html_path;
      stylesheet_path = make_path_to_file(ETC_DIR, HTML_STYLESHEET);
      html_path = make_path_to_file(output_dirname, HTML_FILENAME);
      if (file_exists(stylesheet_path)) {
        print_xml_filename_to_filename_using_stylesheet(
            xml_path,         /* path to XML input file IN */
            stylesheet_path,  /* path to MAST XSL stylesheet IN */
            html_path       /* path to HTML output file IN */
            );
      } else {
        if (verbosity >= NORMAL_VERBOSE) 
          fprintf(stderr, 
            "Warning: could not find the stylesheet \"%s\" required for transformation of xml to html. Have you installed Tomtom correctly?\n", 
            stylesheet_path);
      }
      myfree(stylesheet_path);
      myfree(html_path);
    }
    myfree(xml_path);
    //deprecated - html output is now done with xmlstylesheets
    //print_tomtom_html_tailer(old_html_output, argc, argv, total_matches, sig_type_q, sig_thresh);
    //fclose(old_html_output);
  }
  // clean up motif dbs
  arraylst_destroy(destroy_motif_db, target_dbs);
  // clean up query and target motifs
  free_motifs(query_motifs);
  free_motifs(target_motifs);
  free_array(bg_freqs_query);
  free_array(bg_freqs_target);
  free_array(bg_freqs_mod);
  return (0);
}/* Main tomtom*/
#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
