/**********************************************************************
 * FILE: ama_scan.c
 * AUTHOR: Fabian Buske / Robert McLeay for refactoring / T L Bailey
 * PROJECT: MEME
 * COPYRIGHT: 2007-2008, UQ
 * VERSION: $Revision: 1.0$
 * DESCRIPTION: Routines to perform average motif affinity scans
 *
 **********************************************************************/

#include "ama_scan.h"
#include "background.h"

#define sqrt2 sqrt(2.0)

/*************************************************************************
 * Calculate the odds score for each motif-sized window at each
 * site in the sequence using the given nucleotide frequencies.
 *
 * This function is a lightweight version based on the one contained in
 * motiph-scoring. Several calculations that are unnecessary for gomo
 * have been removed in order to speed up the process.
 * Scores sequence with up to two motifs.
 *************************************************************************/
double score_sequence(
  SEQ_T*        seq,		// sequence to scan (IN)
  double *logcumback,		// cumulative bkg probability of sequence (IN)
  PSSM_PAIR_T*  pssm_pair,	// pos and neg pssms (IN)
  int method, 			// method used for scoring (IN)
  int last, 			//score only last <n> or
				//score all if <n> is zero (IN)
  BOOLEAN_T* isFeasible		// FLAG indicated if there is at least one position
				// where the motif could be matched against (OUT)
)
{
  assert(pssm_pair != NULL);
  assert(seq != NULL);

  PSSM_T* pos_pssm = pssm_pair->pos_pssm;
  assert(pos_pssm != NULL);
  PSSM_T* neg_pssm = pssm_pair->neg_pssm;
  int n_motifs = neg_pssm ? 2 : 1;

  char* raw_seq = get_raw_sequence(seq);
  int seq_length = get_seq_length(seq);
  int w = get_num_rows(pos_pssm->matrix);
  int n = seq_length - w + 1;

  if (verbosity >= DUMP_VERBOSE) {
    fprintf(stderr, "Debug n_motifs: %d seq_length: %d w: %d n: %d.\n", n_motifs, seq_length, w, n);
  }

  // Get alphabet;
  char* alphabet = get_alphabet(FALSE);
  int alph_size = get_alph_size(ALPH_SIZE);

  // Dependent on the "last" parameter, change the starting point
  int start;
  int N_scored;
  if (last > 0 && last < seq_length) {
    start = seq_length - last;
    N_scored  = n_motifs * (last - w + 1);	// number of sites scored
  } else {
    start = 0;
    N_scored  = n_motifs * n;			// number of sites scored
  }

  // For each motif (positive and reverse complement)
  double max_odds = 0.0;
  double sum_odds = 0.0;
  double requested_odds = 0.0;
  int i;

  if (verbosity >= HIGHER_VERBOSE) {
    fprintf(stderr, "Starting scan at position %d .\n", start);
  }

  for (i=0; i<n_motifs; i++) { 	// pos (and negative) motif
    PSSM_T* pssm = (i==0 ? pos_pssm : neg_pssm);	// choose +/- motif
    // For each site in the sequence
    int seq_index;
    for (seq_index = start; seq_index < n; seq_index++) {	// site 
      double odds = 1.0;
      // For each position in the motif window
      int motif_position;
      for (motif_position = 0; motif_position < w; motif_position++) { // column
        int i_site = seq_index + motif_position;
        char c = raw_seq[i_site];
        // Check for gaps at this site
        if (c == '-' || c == '.') { N_scored--; odds = 0; break; }
        // Check for ambiguity codes at this site
        int alph_index = alphabet_index(c, alphabet);
        if (alph_index >= alph_size || alph_index < 0) { N_scored--; odds = 0; break; }
        // multiple odds by value in appropriate motif cell
        odds *= get_matrix_cell(motif_position, alph_index, pssm->matrix);
      } // column
      //
      // Apply sequence-dependent background model.
      //
      if (logcumback) {
        int i_site = seq_index;
        double log_p = logcumback[i_site+w] - logcumback[i_site];	// log Pr(x | background)
        //printf("log_p:: %g motif_pos %d\n", log_p, motif_position);
        double adjust = exp(w*log(1/4.0) - log_p);	// Pr(x | uniform) / Pr(x | background)
        odds *= adjust;
      }
      // Add odds to growing sum.
      sum_odds += odds;				// sum of odds
      if (odds > max_odds) max_odds = odds;	// max of odds
    } // site
  } // motif

  if (verbosity >= HIGHER_VERBOSE) {
    fprintf(stderr, "Scored %d positions with the sum odds %f and the max odds %f.\n", N_scored, sum_odds, max_odds);
  }

  // has there been anything matched at all?
  if (N_scored == 0){
      if (verbosity >= NORMAL_VERBOSE) {
	    fprintf(stderr,"Sequence \'%s\' offers no location to match the motif against (sequence length too short?)\n",get_seq_name(seq));
      }
	  *isFeasible = FALSE;
	  return 0.0;
    // return odds as requested (MAX or AVG scoring)
  } else if (method == AVG_ODDS) {
    requested_odds = sum_odds / N_scored;	// mean
  } else if (method == MAX_ODDS) {
    requested_odds = max_odds;			// maximum
  } else if (method == SUM_ODDS) {
	requested_odds = sum_odds ;	// sum
  }

  return(requested_odds);
} // score_sequence

/**********************************************************************
  ama_sequence_scan()

  Scan a given sequence with a specified motif using either
  average motif affinity scoring or maximum one. In addition z-scores
  may be calculated. Also the scan can be limited to only the end of
  the passed sequences.

  The motif has to be converted to odds in advance (in order
  to speed up the scanning).

  The result will be stored in the scanned_sequence parameter.
 **********************************************************************/
void ama_sequence_scan(
  SEQ_T* sequence,		// the sequence to scan (IN)
  double *logcumback,		// cumulative bkg probability of sequence (IN)
  PSSM_PAIR_T* pssm_pair,	// the pos/neg pssms (IN)
  int scoring,			// AVG_ODDS or MAX_ODDS (IN)
  BOOLEAN_T pvalues,		// compute p-values (IN)
  int last,			// use only last <n> sequence positions
				// or 0 if all positions should be used
  SCANNED_SEQUENCE_T* scanned_seq, // the scanned sequence results (OUT)
  BOOLEAN_T* need_postprocessing // Flag indicating the need for postprocessing (OUT)
)
{
  assert(sequence != NULL);
  assert(pssm_pair != NULL);

  // FLAG indicates if sequence was suitable for motif matching
  BOOLEAN_T isFeasible = TRUE;

  // Score the sequence.
  double odds = score_sequence(sequence, logcumback, pssm_pair, scoring, last, &isFeasible);
	
  // Compute the p-value of the AVG_ODDS score.
  if (get_scanned_sequence_num_scanned_positions(scanned_seq) == 0L){
	  set_scanned_sequence_score(scanned_seq, odds);
	  // sequence has not been scanned before
	  if (!isFeasible){
          if (verbosity >= NORMAL_VERBOSE) {
          fprintf(stderr,"Sequence '%s' not suited for motif. P-value set to 1.0!\n",
            get_scanned_sequence_accession(scanned_seq));
          }
		set_scanned_sequence_pvalue(scanned_seq, 1.0);
	  } else if (odds < 0.0){
          if (verbosity >= NORMAL_VERBOSE) {
            fprintf(stderr,"Sequence '%s' got invalid (negative) odds score. P-value set to 1.0!\n",
		       get_scanned_sequence_accession(scanned_seq));
          }
		set_scanned_sequence_pvalue(scanned_seq, 1.0);
	  } else if (pvalues && scoring == AVG_ODDS) {
		double pvalue = get_ama_pv(odds, get_scanned_sequence_length(scanned_seq), 
		  get_total_gc_sequence(sequence), pssm_pair);
		set_scanned_sequence_pvalue(scanned_seq, pvalue);
	  }
	  // scanned_element is used to keep track how often a sequence has been scored
	  // this feature is used in downstream gomo where a one2many homolog relationship
	  // is encoded through the same sequence identifier
	  add_scanned_sequence_scanned_element(scanned_seq); 
  } else {
	  // sequence has been scored before
	  if(!has_scanned_sequence_score(scanned_seq))
		  // no score set yet, so do
		  set_scanned_sequence_score(scanned_seq, odds);
	  else {
		  // sum scores (take average later)
		  set_scanned_sequence_score(scanned_seq, odds + get_scanned_sequence_score(scanned_seq));
	  }
	  if (!isFeasible){
          if (verbosity >= NORMAL_VERBOSE) {
              fprintf(stderr,"Sequence '%s' not suited for motif. P-value set to 1.0!\n",
				  get_scanned_sequence_accession(scanned_seq));
          }
		  if (!has_scanned_sequence_pvalue(scanned_seq))
			  set_scanned_sequence_pvalue(scanned_seq, 1.0);
	  } else if (odds < 0.0){
          if (verbosity >= NORMAL_VERBOSE) {
              fprintf(stderr,"Sequence '%s' got invalid (negative) odds score. P-value set to 1.0!\n",
				  get_scanned_sequence_accession(scanned_seq));
          }
		  if (!has_scanned_sequence_pvalue(scanned_seq))
			  set_scanned_sequence_pvalue(scanned_seq, 1.0);
	  } else if (pvalues && scoring == AVG_ODDS) {
		  double pvalue = get_ama_pv(odds, get_scanned_sequence_length(scanned_seq), 
									 get_total_gc_sequence(sequence), pssm_pair);
		  if (!has_scanned_sequence_pvalue(scanned_seq))
			  set_scanned_sequence_pvalue(scanned_seq, pvalue);
		  else {
			  // keep minimum p-value only
			  set_scanned_sequence_pvalue(scanned_seq, min(pvalue, get_scanned_sequence_pvalue(scanned_seq)));
		  }
	  }
	  add_scanned_sequence_scanned_element(scanned_seq); 
	  *need_postprocessing = TRUE;
  }

} // ama_sequence_scan
