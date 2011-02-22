#ifndef SPAMO_OUTPUT_H
#define SPAMO_OUTPUT_H

#include <stdio.h>

#include "motif.h"
#include "spamo-matches.h"

extern const char* xml_indicator;

extern const char* spamo_dtd;

/**************************************************************************
 * Outputs indent multiples of tab
 * Assumes indent is non-negative
 **************************************************************************/
void output_indent(FILE *file, char *tab, int indent);

/**************************************************************************
 * Outputs xml for a fasta database
 **************************************************************************/
void output_sequence_database(FILE *xml_output, SEQUENCE_DB_T *db, char *tab, int indent, char **buffer, int *buffer_len);

/**************************************************************************
 * Outputs xml for a motif database
 **************************************************************************/
void output_motif_database(FILE *xml_output, MOTIF_DB_T *db, char *tab, int indent, char **buffer, int *buffer_len);

/**************************************************************************
 * Outputs xml for a motif
 **************************************************************************/
void output_motif(FILE *xml_output, int db_id, MOTIF_T *motif, char *tab, int indent, char **buffer, int *buffer_len);

/**************************************************************************
 * Outputs xml for the histogram of the spacings
 **************************************************************************/
void output_secondary_motif(FILE *xml_output, SECONDARY_MOTIF_T *smotif, LINKLST_T *rmotifs, 
    int margin, int bin, char *tab, int indent, char **buffer, int *buffer_len);

/**************************************************************************
 * Creates histograms
 **************************************************************************/
void create_histograms
(int margin, int binsize, int binmax, double sigthresh, char *dir, 
 MOTIF_T *primary_motif, LINKLST_T *secondary_motifs, 
 BOOLEAN_T make_eps, BOOLEAN_T make_png, BOOLEAN_T make_for_redundant);

/**************************************************************************
 * Create an output file and dump the sequence matches to file.
 **************************************************************************/
void output_sequence_matches(char *dir, int margin, RBTREE_T *sequences, 
    MOTIF_T *primary_motif, SECONDARY_MOTIF_T *secondary_motif, int *matches);

#endif
