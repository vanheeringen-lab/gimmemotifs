/*************************************************************************
 * FILE: scored-sites.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 15 November 2007
 * PROJECT: MEME
 * COPYRIGHT: 2007, UW
 * DESCRIPTION: A site and a list of scored sites, ready to be printed.
 *************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include "utils.h"
#include "array.h"
#include "qvalue.h"
#include "cisml.h"
#include "scored-sites.h"

// Maximum field width for any string in a scored site.
#define MAX_ENTRY 100

struct scored_site_t {
  char   seq_name[MAX_ENTRY];
  char   feature_type[MAX_ENTRY];
  int    start;
  int    end;
  int    seq_length;
  double logodds;
  double pvalue;
  double qvalue;
  char   strand;
  char   motif_id[MAX_ENTRY];
  char   window_seq[MAX_ENTRY];
};

struct scored_sites_t {
  int             num_sites; // Number of sites in this list.
  int             max_sites; // Allocated space in this list.
  int             total_sites; // Sites in list, plus those that weren't saved.
  BOOLEAN_T       has_qvalues; // Does the list contain qvalues?
  SCORED_SITE_T** sites;
};

/*************************************************************************
 * Allocate dynamic memory for an empty list of sites.
 *************************************************************************/
#define INITIAL_MAX_SITES 1000
SCORED_SITES_T* new_scored_sites
  ()
{
  SCORED_SITES_T* new_sites;  /* The list being created. */

  new_sites = (SCORED_SITES_T*)mm_calloc(1, sizeof(SCORED_SITES_T));
  new_sites->sites = (SCORED_SITE_T**)mm_calloc(INITIAL_MAX_SITES,
						sizeof(SCORED_SITE_T*));
  new_sites->num_sites = 0;
  new_sites->max_sites = INITIAL_MAX_SITES;
  new_sites->total_sites = 0;
  new_sites->has_qvalues = FALSE;
  return(new_sites);
}

/*************************************************************************
 * Allocate and initialize a single site.
 *************************************************************************/
static SCORED_SITE_T* new_scored_site
  (char*  seq_name,
   char*  feature_type,
   int    start,
   int    end,
   int    seq_length,
   double logodds,
   double pvalue,
   char   strand,
   char*  motif_id,
   char*  window_seq)
{
  SCORED_SITE_T* return_value = 
    (SCORED_SITE_T*)mm_malloc(sizeof(SCORED_SITE_T));

  strncpy(return_value->seq_name, seq_name, MAX_ENTRY);
  strncpy(return_value->feature_type, feature_type, MAX_ENTRY);
  return_value->start = start;
  return_value->end = end;
  return_value->seq_length = seq_length;
  return_value->logodds = logodds;
  return_value->pvalue = pvalue;
  return_value->qvalue = 0.0;
  return_value->strand = strand;
  strncpy(return_value->motif_id, motif_id, MAX_ENTRY);
  strncpy(return_value->window_seq, window_seq, MAX_ENTRY);

  return(return_value);
}

/*************************************************************************
 * Is the given list a null list?
 *************************************************************************/
static void check_null_sites
  (SCORED_SITES_T*  site_list)
{
  if (site_list == NULL) {
    die("Attempted to access null list of scored sites.\n");
  }
}

/*************************************************************************
 * Add a new site to a given list.
 *************************************************************************/
void add_scored_site
  (char*  seq_name,
   char*  feature_type,
   int    start,
   int    end,
   int    seq_length,
   double logodds,
   double pvalue,
   char   strand,
   char*  motif_id,
   char*  window_seq,
   SCORED_SITES_T* site_list)
{
  // Check for nulls.
  check_null_sites(site_list);

  // Create the new scored site object.
  SCORED_SITE_T* this_site = new_scored_site(
					     seq_name,
					     feature_type,
					     start,
					     end,
					     seq_length,
					     logodds,
					     pvalue,
					     strand,
					     motif_id,
					     window_seq
					     );

  // Reallocate space if there isn't any.
  if (site_list->num_sites >= site_list->max_sites) {
    site_list->sites = (SCORED_SITE_T**)mm_realloc(site_list->sites, 
						  (site_list->max_sites 
						   + INITIAL_MAX_SITES)
						   * sizeof(SCORED_SITE_T*));
    site_list->max_sites += INITIAL_MAX_SITES;
  }

  // Put the site in the list.
  site_list->sites[site_list->num_sites] = this_site;
  (site_list->num_sites)++;
  (site_list->total_sites)++;

}

/*************************************************************************
 * Compare two sites by pvalue for 'qsort'.
 *************************************************************************/
static int site_compare_by_pvalue
  (const void* elem1,
   const void* elem2)
{
  const double key1 = (*(SCORED_SITE_T **)elem1)->pvalue;
  const double key2 = (*(SCORED_SITE_T **)elem2)->pvalue;

  if (key1 < key2) {
    return(-1);
  } else if (key1 > key2) {
    return(1);
  }
  return(0);
}

/*************************************************************************
 * Compare two sites by sequence ID and position for 'qsort'.
 *************************************************************************/
static int site_compare_by_id
  (const void* elem1,
   const void* elem2)
{
  const char* key1 = (*(SCORED_SITE_T **)elem1)->seq_name;
  const char* key2 = (*(SCORED_SITE_T **)elem2)->seq_name;

  int comparison = strcmp(key1, key2);
  if (comparison < 0) {
    return(-1);
  } else if (comparison > 0) {
    return(1);
  } else {
    const int index1 = (*(SCORED_SITE_T **)elem1)->start;
    const int index2 = (*(SCORED_SITE_T **)elem2)->start;
    if (index1 < index2) {
      return(-1);
    } else if (index1 > index2) {
      return(1);
    }
  }
  return(0);
}

/*************************************************************************
 * Sort a given set of sites by pvalue or by sequence ID and position.
 *************************************************************************/
static void sort_sites
  (BOOLEAN_T sort_by_pvalue,
   SCORED_SITES_T* site_list)
{
  check_null_sites(site_list);

  // Tell the user what's up.
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Sorting %d sites ", site_list->num_sites);
    if (sort_by_pvalue) {
      fprintf(stderr, "by p-value.\n");
    } else {
      fprintf(stderr, "by sequence ID and start position.\n");
    }
  }

  if (sort_by_pvalue) {
    qsort(
	  (void *)(site_list->sites),
	  site_list->num_sites, 
	  sizeof(SCORED_SITE_T *),
	  site_compare_by_pvalue
	  );
  } else {
    qsort(
	  (void *)(site_list->sites),
	  site_list->num_sites, 
	  sizeof(SCORED_SITE_T *),
	  site_compare_by_id
	  );
  }
}

/*************************************************************************
 * Assuming that the list is sorted, replace p-values with q-values.
 *************************************************************************/
void convert_pvalues_to_qvalues
  (SCORED_SITES_T* site_list)
{

  // Tell the user what's up.
  int num_sites = site_list->num_sites;
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(
	    stderr,
	    "Computing q-values (%d stored p-values, %d total p-values).\n", 
	    num_sites,
	    site_list->total_sites
	    );
  }

  // Sort by p-value.
  sort_sites(TRUE, site_list);

  // Extract the p-values into an array.
  ARRAY_T* pvalues = allocate_array(num_sites);
  int i_site;
  for (i_site = 0; i_site < num_sites; i_site++) {
    set_array_item(i_site, site_list->sites[i_site]->pvalue, pvalues);
  }

  // Convert them to q-values.
  compute_qvalues(
    FALSE, // Don't stop with FDR.
		TRUE, // Don't estimate pi-zero.
		NULL, // Don't store pi-zero in a file.
		NUM_BOOTSTRAPS,
		NUM_BOOTSTRAP_SAMPLES,
		NUM_LAMBDA,
		MAX_LAMBDA,
		site_list->total_sites,
		pvalues,
    NULL // No sampled p-values
  );

  // Put them back into the list.
  for (i_site = 0; i_site < num_sites; i_site++) {
    site_list->sites[i_site]->qvalue = get_array_item(i_site, pvalues);
  }
  site_list->has_qvalues = TRUE;

  // Sort by sequence ID and position.
  sort_sites(FALSE, site_list);

}

/*************************************************************************
 * Increment total site count, without storing a new site.
 *************************************************************************/
void increment_total_sites
(SCORED_SITES_T* site_list)
{
  check_null_sites(site_list);
  (site_list->total_sites)++;
}
  
/*************************************************************************
 * Print a list of sites.
 *************************************************************************/
void print_scored_sites
  (char*           program_name,
   SCORED_SITES_T* site_list,
   FILE*           gff_file, // May be NULL.
   FILE*           cisml_file)
{
  int i_site;
  check_null_sites(site_list);

  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Printing %d out of %d sites.\n", site_list->num_sites, site_list->total_sites);
  }

  char* current_sequence = NULL;
  for (i_site = 0; i_site < site_list->num_sites; i_site++) {
    SCORED_SITE_T* this_site = site_list->sites[i_site];

    // Are we at the start of a new sequence?
    BOOLEAN_T new_sequence
      = ((current_sequence == NULL) ||
	 (strcmp(current_sequence, this_site->seq_name) != 0));

    // Only print GFF if the file handle is given.
    if (gff_file) {

      // If it's a new sequence, print the sequence line.
      if (new_sequence) {
	fprintf(gff_file, "%s", this_site->seq_name);
	fprintf(gff_file, "\t%s", program_name);
	fprintf(gff_file, "\tsequence");
	fprintf(gff_file, "\t1"); // Start
	fprintf(gff_file, "\t%d", this_site->seq_length); // End
	fprintf(gff_file, "\t."); // Score
	fprintf(gff_file, "\t."); // Strand
	fprintf(gff_file, "\t."); // Frame
	fprintf(gff_file, "\t.\n"); // Comment
      }
	
      // Print the motif line.
      fprintf(gff_file, "%s", this_site->seq_name);
      fprintf(gff_file, "\t%s", program_name);
      fprintf(gff_file, "\tmotif");
      fprintf(gff_file, "\t%d", this_site->start);
      fprintf(gff_file, "\t%d", this_site->end);
      fprintf(gff_file, "\t%6.4g", this_site->pvalue);
      fprintf(gff_file, "\t%c", this_site->strand);
      fprintf(gff_file, "\t.");  // Frame
      fprintf(gff_file, "\tmotif \"%s\"", this_site->motif_id);
      fprintf(gff_file, ";sequence \"%s\"", this_site->window_seq);
      fprintf(gff_file, ";log-odds \"%g\"", this_site->logodds);
      fprintf(gff_file, "\n");
    }

    current_sequence = this_site->seq_name;
  }
}

/*************************************************************************
 * Free memory associated with a list of sites.
 *************************************************************************/
void free_scored_sites
  (SCORED_SITES_T* site_list)
{
  int i_site;
  if (site_list == NULL) {
    return;
  }

  for (i_site = 0; i_site < site_list->num_sites; i_site++) {
    myfree(site_list->sites[i_site]);
  }

  myfree(site_list->sites);
  myfree(site_list);
}

