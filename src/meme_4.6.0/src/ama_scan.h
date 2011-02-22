/**********************************************************************
 * FILE: ama_scan.h
 * AUTHOR: Fabian Buske / Robert McLeay for refactoring
 * PROJECT: MEME
 * COPYRIGHT: 2007-2008, UQ
 * VERSION: $Revision: 1.0$
 * DESCRIPTION: Routines to perform average motif affinity scans
 *
 **********************************************************************/

#ifndef AMA_SCAN_H_
#define AMA_SCAN_H_

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "meme-io.h"
#include "string-list.h"
#include "macros.h"
#include "motif.h"
#include "matrix.h"
#include "array.h"
#include "pssm.h"
#include "seq.h"
#include "alphabet.h"
#include "cisml.h"

extern char* motif_name; // Use this motif name in the output.

extern VERBOSE_T verbosity;

#define AVG_ODDS 0
#define MAX_ODDS 1
#define SUM_ODDS 2

#define min(a,b)      (a<b)?a:b
#define max(a,b)      (a>b)?a:b

/**********************************************************************
  ama_sequence_scan(
 **********************************************************************/
void ama_sequence_scan(
  SEQ_T* sequence,              // the sequence to scan (IN)  
  double *logcumback,           // cumulative bkg probability of sequence (IN)
  PSSM_PAIR_T* pssm_pair,       // the pos/neg pssms (IN)  
  int scoring,                  // AVG_ODDS or MAX_ODDS (IN)
  BOOLEAN_T pvalues,            // compute p-values
  int last,			// use only the last <n> positions,
				// specify 0 to use all positions (IN)
  SCANNED_SEQUENCE_T* scanned_seq, // the scanned sequence results (OUT)
  BOOLEAN_T* need_postprocessing // Flag indicating the need for postprocessing (OUT)
); 
#endif /* AMA_SCAN_H_ */
