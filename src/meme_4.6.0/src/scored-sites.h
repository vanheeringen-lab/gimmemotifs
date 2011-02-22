/*************************************************************************
 * FILE: scored-sites.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 15 November 2007
 * PROJECT: MEME
 * COPYRIGHT: 2007, UW
 * DESCRIPTION: A site and a list of scored sites, ready to be printed.
 *************************************************************************/
#ifndef SCORED_SITES_H
#define SCORED_SITES_H

#include <stdio.h>

/****************************************************************************
 * Structured objects corresponding to one or multiple output lines.
 ****************************************************************************/
typedef struct scored_site_t SCORED_SITE_T;
typedef struct scored_sites_t SCORED_SITES_T;

/*************************************************************************
 * Allocate dynamic memory for a list of sites.
 *************************************************************************/
SCORED_SITES_T* new_scored_sites
  ();

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
   SCORED_SITES_T* site_list);

/*************************************************************************
 * Assuming that the list is sorted and that the scores in the list
 * are p-values, replace scores with q-values.
 *************************************************************************/
void convert_pvalues_to_qvalues
  (SCORED_SITES_T* site_list);

/*************************************************************************
 * Increment total site count, without storing a new site.
 *************************************************************************/
void increment_total_sites
(SCORED_SITES_T* site_list);

/*************************************************************************
 * Print a list of sites.
 *************************************************************************/
void print_scored_sites
  (char*           program_name,
   SCORED_SITES_T* site_list,
   FILE*           gff_file, // May be NULL.
   FILE*           cisml_file);

/*************************************************************************
 * Free memory associated with a list of sites.
 *************************************************************************/
void free_scored_sites
  (SCORED_SITES_T* site_list);

#endif
