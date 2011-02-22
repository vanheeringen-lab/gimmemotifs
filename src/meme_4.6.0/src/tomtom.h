/****************************************************************************
 * FILE: tomtom.h
 * AUTHOR: Shobhit Gupta
 * CREATE DATE: 02/09/2006
 * PROJECT: TOMTOM
 * DESCRIPTION: Compare two motifs using various statistics.
 *              This file is only necessary for the check utility.
 * COPYRIGHT: 2006, UW
 ****************************************************************************/
#include "matrix.h"
#include "array.h"
#include "motif.h"
#include "string-list.h"
#include "utils.h"

/**************************************************************************
 * Computes sandelin scores for all columns of the query matrix against
 * target motif. This can then be used to by compute_overlap_score.
 **************************************************************************/
void sandelin_scores(
         MOTIF_T*,//query_motif
         MOTIF_T*,//target_motif
         ARRAY_T*,//background_frequencies
         MATRIX_T**//output score matrix
         );

/**************************************************************************
 * Computes Euclidean scores for all columns of the query matrix against target
 * motif. This can then be used to by compute_overlap_score.
 **************************************************************************/
void ed_scores(
    MOTIF_T*,//query_motif
    MOTIF_T*,//target_motif
    ARRAY_T*,//background_frequencies
    MATRIX_T**//output score matrix
    );

/**************************************************************************
 * Computes Kullback scores for all columns of the query matrix against target
 * motif. This can then be used to by compute_overlap_score.
 **************************************************************************/
void kullback_scores(
         MOTIF_T*,//query_motif
         MOTIF_T*,//target_motif
         ARRAY_T*,//background_frequencies
         MATRIX_T**//output score matrix
         );

/**************************************************************************
 * Computes ALLR scores for all columns of the query matrix against target
 * motif. This can then be used to by compute_overlap_score.
 **************************************************************************/
void allr_scores(
        MOTIF_T*,//query_motif
        MOTIF_T*,//target_motif
        ARRAY_T*,//background_frequencies
        MATRIX_T**//output score matrix
        );

/**************************************************************************
 * Computes Pearson scores for all columns of the query matrix against target
 * motif. This can then be used to by compute_overlap_score.
 **************************************************************************/
void pearson_scores(
        MOTIF_T*,//query_motif
        MOTIF_T*,//target_motif
        ARRAY_T*,//background_frequencies
        MATRIX_T**//output score matrix
        );

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
