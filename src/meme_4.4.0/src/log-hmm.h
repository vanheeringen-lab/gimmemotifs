/**************************************************************************
 * FILE: log-hmm.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 1-28-98
 * PROJECT: MHMM
 * DESCRIPTION: Convert between a normal HMM and a log HMM.
 **************************************************************************/
#ifndef LOG_HMM_H
#define LOG_HMM_H

#include "mhmm-state.h"
#include "utils.h"

#define DEFAULT_PROTEIN_PAM 250	// Default PAM distance for proteins.
#define DEFAULT_DNA_PAM 1      	// Default PAM distance for DNA.
#define DEFAULT_PROTEIN_BETA 10	// Default beta for proteins.
#define DEFAULT_DNA_BETA 1 	// Default beta for DNA.

/**************************************************************************
 * Convert an HMM to or from log form.
 *
 * Allocate space when converting to log form. FIXME
 **************************************************************************/
void convert_to_from_log_hmm
  (BOOLEAN_T to_log,         // Convert to log form?
   BOOLEAN_T zero_spacer_emit_lo,     // Set spacer emission log-odds = 0?
   double    gap_open,       // Cost to open a gap; ignored if < 0
   double    gap_extend,     // Cost to extend a gap; ignored if < 0
   ARRAY_T*  background,     // The background distribution.
   char *    sc_filename,    // Name of score file (ignored if background NULL).
   int 	     pam_dist,       // PAM distance (ignored if sc_filename not NULL).
   double    beta,           // Weight on pseudocounts.
   MHMM_T*   source_hmm,     // The HMM to convert.
   MHMM_T**  target_hmm);    // The same HMM in log/log-odds form

#endif
