/********************************************************************
 * FILE: ama.c
 * AUTHOR: Fabian Buske and Timothy L. Bailey
 * CREATE DATE: 03/07/2008
 * PROJECT: MEME suite
 * COPYRIGHT: 2008, UQ
 *
 * AMA (average motif affinity) is part of an implementation of the
 * algorithm described in
 * "Associating transcription factor binding site motifs with target
 * Go terms and target genes"
 * authors: Mikael Boden and Timothy L. Bailey
 * published: Nucl. Acids Res (2008)
 *
 * The implementation is based on the fimo class in order to
 * use the scoring schemes described in the publication.
 *
 * AMA works on DNA data only.
 ********************************************************************/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <errno.h>
#include "matrix.h"
#include "alphabet.h"
#include "cisml.h"
#include "fasta-io.h"
#include "meme-io.h"
#include "string-list.h"
#include "simple-getopt.h"
#include "read_tamo.h"
#include "array.h"
#include "alphabet.h"
#include "macros.h"
#include "motif.h"
#include "matrix.h"
#include "gendb.h"
#include "pssm.h"

#include "ama_scan.h"

const int GFF_FORMAT = 0; // GFF file format
const int CISML_FORMAT = 1; // cisML file format

const int MEME_FORMAT = 1; // input format is meme
const int TAMO_FORMAT = 2; // input format is tamo

char* program_name = NULL;

// Nucleotide alphabet order as in motif.h
extern char alphabet[];

// As from our ama.h - this is done this way
// to avoid linking issues
VERBOSE_T verbosity = NORMAL_VERBOSE;

/**********************************************************************
 print_score()

 outputs the scores in a gff format
 **********************************************************************/
void print_score(char* program_name, CISML_T* cisml, FILE* gff_file) {

  // iterate over all patterns
  int num_patterns = get_cisml_num_patterns(cisml);
  if (num_patterns > 0) {
    PATTERN_T **patterns = get_cisml_patterns(cisml);
    int i = 0;
    for (i = 0; i < num_patterns; ++i) {
      // iterate over all sequences
      int num_seq = 0;
      num_seq = get_pattern_num_scanned_sequences(patterns[i]);
      SCANNED_SEQUENCE_T **sequences = get_pattern_scanned_sequences(
		      patterns[i]);
      int j = 0;
      for (j = 0; j < num_seq; ++j) {
	double score = 0.0;
        double pvalue = 1.0;

	if (has_scanned_sequence_score(sequences[j])) {
	  score = get_scanned_sequence_score(sequences[j]);
        }
	if (has_scanned_sequence_pvalue(sequences[j])) {
	  pvalue = get_scanned_sequence_pvalue(sequences[j]);
	}
	fprintf(gff_file, "%s", get_scanned_sequence_accession(
			sequences[j]));
	fprintf(gff_file, "\t%s", program_name);
	fprintf(gff_file, "\tsequence");
	fprintf(gff_file, "\t1"); // Start
	fprintf(gff_file, "\t%d", get_scanned_sequence_length(
			sequences[j])); // End
	fprintf(gff_file, "\t%g", score); // Score
	fprintf(gff_file, "\t%g", pvalue); // Score
	fprintf(gff_file, "\t."); // Strand
	fprintf(gff_file, "\t."); // Frame
	fprintf(gff_file, "\t.\n"); // Comment
      } // num_seq
    } // num_pattern
  } // pattern > 0
}

/*************************************************************************
 * Entry point for ama
 *************************************************************************/
int main(int argc, char *argv[]) {
  int max_seq_length = MAX_SEQ;
  STRING_LIST_T* selected_motifs = NULL;
  double pseudocount = 0.01;
  int output_format = CISML_FORMAT;
  program_name = "ama";
  int scoring = AVG_ODDS;
  BOOLEAN_T pvalues = FALSE;
  BOOLEAN_T normalize_scores = FALSE;
  int num_gc_bins = 1;
  int input_format = MEME_FORMAT;
  BOOLEAN_T scan_both_strands = TRUE;
  ARRAY_T* pos_bg_freqs = NULL;
  ARRAY_T* rev_bg_freqs = NULL;
  clock_t c0, c1; /* measuring cpu_time */
  CISML_T *cisml;
  int i;
  int last = 0;

  /**********************************************
   * COMMAND LINE PROCESSING
   **********************************************/

  const int num_options = 12;
  cmdoption const motif_scan_options[] = {
    { "i-format", REQUIRED_VALUE },
    { "max-seq-length", REQUIRED_VALUE },
    { "motif", REQUIRED_VALUE },
    { "motif-pseudo", REQUIRED_VALUE },
    { "rma", NO_VALUE },
    { "pvalues", NO_VALUE },
    { "norc", NO_VALUE },
    { "o-format", REQUIRED_VALUE },
    { "scoring", REQUIRED_VALUE },
    { "verbosity", REQUIRED_VALUE },
    { "gcbins", REQUIRED_VALUE },
    { "last", REQUIRED_VALUE }
  };

  int option_index = 0;

  // Define the usage message.
  char usage[] = "USAGE: ama [options] <motif file> <sequence file> <background file>\n"
    "\n"
    "   Options:\n"
    "     --motif <id>\t\t\tUse only the motif identified by <id>.\n"
    "       \t\t\t\t\tThis option may be repeated.\n"
    "     --motif-pseudo <float>\t\tThe value <float> times the background\n"
    "       \t\t\t\t\tfrequency is added to the count of each letter when creating \n"
    "       \t\t\t\t\tthe likelihood ratio matrix. (default: %g)\n"
    "     --norc\t\t\t\tDisables the scanning of the reverse complement strand.\n"
    "     --scoring [avg-odds|max-odds]\tIndicates whether the average or \n"
    "       \t\t\t\t\tthe maximum odds should be calculated (default: avg-odds)\n"
    "     --rma\t\t\t\tScale motif scores to the range 0-1. (Relative Motif Affinity).\n"
    "       \t\t\t\t\tMotif scores are scaled by the maximum score achievable by that PWM.\n"
    "       \t\t\t\t\t(default: motif scores are not normalized)\n"
    "     --pvalues\t\t\t\tPrint p-value of avg-odds score in cisml output.\n"
    "       \t\t\t\t\tIgnored for max-odds scoring.\n"
    "       \t\t\t\t\t(default: p-values are not printed)\n"
    "     --gcbins <bins>\t\t\tCompensate p-values for GC content of each sequence\n"
    "       \t\t\t\t\tusing given number of GC range bins. Recommended bins: 41.\n"
    "       \t\t\t\t\t(default: p-values are based on frequencies in background file)\n"
    "     --o-format [gff|cisml]\t\tOutput file format (default: cisml)\n"
    "     --verbosity [1|2|3|4]\t\tControls amount of screen output (default: %d)\n"
    "     --max-seq-length <int>\t\tSet the maximum length allowed for \n"
    "       \t\t\t\t\tinput sequences. (default: %d)\n"
    "     --last <int>\t\t\tUse only scores of (up to) last <n> sequence \n"
    "       \t\t\t\t\tpositions to compute AMA.\n"
    "\n";

  // Parse the command line.
  if (simple_setopt(argc, argv, num_options, motif_scan_options) != NO_ERROR) {
	  die("Error processing command line options: option name too long.\n");
  }

  while (TRUE) {
    int c = 0;
    char* option_name = NULL;
    char* option_value = NULL;
    const char * message = NULL;

    // Read the next option, and break if we're done.
    c = simple_getopt(&option_name, &option_value, &option_index);
    if (c == 0) {
	    break;
    } else if (c < 0) {
	    (void) simple_getopterror(&message);
	    die("Error processing command line options (%s).\n", message);
    } else if (strcmp(option_name, "max-seq-length") == 0) {
	    max_seq_length = atoi(option_value);
    } else if (strcmp(option_name, "norc") == 0) {
	    scan_both_strands = FALSE;
    } else if (strcmp(option_name, "motif") == 0) {
	    if (selected_motifs == NULL) {
		    selected_motifs = new_string_list();
	    }
	    add_string(option_value, selected_motifs);
    } else if (strcmp(option_name, "motif-pseudo") == 0) {
	    pseudocount = atof(option_value);
    } else if (strcmp(option_name, "o-format") == 0) {
	    if (strcmp(option_value, "gff") == 0)
		    output_format = GFF_FORMAT;
	    else if (strcmp(option_value, "cisml") == 0)
		    output_format = CISML_FORMAT;
	    else {
		    if (verbosity >= NORMAL_VERBOSE)
			    fprintf(stderr,
					    "Output format not known. Use standard instead (cisML).\n");
		    output_format = CISML_FORMAT;
	    }
    } else if (strcmp(option_name, "verbosity") == 0) {
	    verbosity = atoi(option_value);
    } else if (strcmp(option_name, "scoring") == 0) {
	    if (strcmp(option_value, "max-odds") == 0)
		    scoring = MAX_ODDS;
	    else if (strcmp(option_value, "avg-odds") == 0)
		    scoring = AVG_ODDS;
	    else
		    die("Specified scoring scheme not known.\n", message);
    } else if (strcmp(option_name, "pvalues") == 0) {
	    pvalues = TRUE;
    } else if (strcmp(option_name, "rma") == 0) {
	    normalize_scores = TRUE;
      fprintf(stderr, "Normalizing motif scores using RMA method.\n");
    } else if (strcmp(option_name, "gcbins") == 0) {
	    num_gc_bins = atoi(option_value);
	    pvalues = TRUE;
	    if (num_gc_bins <= 1) die("Number of bins in --gcbins must be greater than 1.\n", message);
    }
    // tamo format is not encouraged
    else if (strcmp(option_name, "i-format") == 0) {
	    if (strcmp(option_value, "tamo") == 0)
		    input_format = TAMO_FORMAT;
	    else if (strcmp(option_value, "meme") == 0)
		    input_format = MEME_FORMAT;
    }
    else if (strcmp(option_name, "last") == 0) {
      int i = 0;
      if (option_value[0] == '-') ++i;
      while (option_value[i] != '\0') {
        if (!isdigit(option_value[i])) {
          die("Specified parameter 'last' contains non-numeric characters.\n");
        }
        ++i;
      }
      last = atoi(option_value);
      if (errno != 0) {
        die("Specified parameter 'last' could not be parsed as a number as:\n%s\n",strerror(errno));
      }
      if (last < 0) {
        die("Specified parameter 'last' had negative value (%d) when only postive or zero values are allowed \n", last);
      }
    }
  }

  // Must have sequence and motif
  if (argc != option_index + 3) {
    fprintf(stderr, usage, pseudocount, verbosity, max_seq_length);
    exit(EXIT_FAILURE);
  }
  char* motif_filename = argv[option_index];
  option_index++;
  char* fasta_filename = argv[option_index];
  option_index++;
  char* bg_filename = argv[option_index];
  option_index++;

  // measure time
  c0 = clock();

  // Create cisml data structure for recording results
  cisml = allocate_cisml(program_name, motif_filename, fasta_filename);
  set_cisml_background_file(cisml, bg_filename);

  /**********************************************
   * Read the motifs and background model.
   **********************************************/
  int num_motifs = 0;
  int num_motifs_all = 0;
  MOTIF_T motifs[2 * MAX_MOTIFS];
  PSSM_PAIR_T* pssm_pairs[MAX_MOTIFS];	// note pssm_pairs is an array of pointers

  BOOLEAN_T has_reverse_strand;
  BOOLEAN_T read_file = FALSE;

  if (input_format == MEME_FORMAT) {
    //this reads any meme file, xml, txt and html
    read_meme_file(motif_filename, bg_filename, pseudocount,
        REQUIRE_PSPM, //read pspm, not pssm
        &num_motifs, motifs, NULL,//motif occurrences
        &has_reverse_strand, &pos_bg_freqs);
  } else if (input_format == TAMO_FORMAT) {
   set_alphabet(verbosity, "ACGT");
    // read in tamo formatted file
    MOTIF_T* m = NULL;
    read_tamo(motif_filename, pseudocount, &num_motifs, &m);
    int l;
    int max_motifs = min(num_motifs, MAX_MOTIFS);
    for (l = 0; l < max_motifs; ++l) {
	    copy_motif(&(m[l]), &(motifs[l]));
    }
    for (l = num_motifs; l > 0; --l) {
	    free_motif(&(m[l - 1]));
    }
    pos_bg_freqs = read_background_file(bg_filename);
  }

  if (verbosity >= NORMAL_VERBOSE) fprintf(stderr, "Number of motifs in file %d.\n", num_motifs);

  // make a CISML pattern to hold scores for each motif
  PATTERN_T** patterns = NULL;
  Resize(patterns, num_motifs, PATTERN_T*);
  int motif_index;
  for (motif_index = 0; motif_index < num_motifs; motif_index++) {
    MOTIF_T* motif = &(motifs[motif_index]);
    char* motif_id = get_motif_id(motif);
    // Drop the strand indicator from the motif id if doing both strands
    if (scan_both_strands && (*motif_id == '+' || *motif_id == '-')) motif_id++;
    patterns[motif_index] = allocate_pattern(motif_id, "");
    add_cisml_pattern(cisml, patterns[motif_index]);
  }

  // make reverse complement motifs and background frequencies.
  num_motifs_all = num_motifs;
  if (scan_both_strands == TRUE) {
    add_reverse_complements(&num_motifs_all, motifs);
    assert(num_motifs_all == (2 * num_motifs));
    rev_bg_freqs = allocate_array(get_array_length(pos_bg_freqs));
    complement_dna_freqs(pos_bg_freqs, rev_bg_freqs);
  }

  /**************************************************************
   * Convert motif matrices into log-odds matrices.
   * Scale them.
   * Compute the lookup tables for the PDF of scaled log-odds scores.
   **************************************************************/
  int ns = scan_both_strands ? 2 : 1;	// number of strands
  for (motif_index = 0; motif_index < num_motifs; motif_index++) {
    /*
     *  Note: If scanning both strands, we complement the motif frequencies
     *  but not the background frequencies so the motif looks the same.
     *  However, the given frequencies are used in computing the p-values
     *  since they represent the frequencies on the negative strands.
     *  (If we instead were to complement the input sequence, keeping the
     *  the motif fixed, we would need to use the complemented frequencies
     *  in computing the p-values.  Is that any clearer?)
    */
    double range = 300;		// 100 is not very good; 1000 is great but too slow
    PSSM_T* pos_pssm =
      build_motif_pssm(&(motifs[ns*motif_index]), pos_bg_freqs, pos_bg_freqs, range, num_gc_bins, TRUE);
    PSSM_T* neg_pssm = (scan_both_strands ?
      build_motif_pssm(&(motifs[ns*motif_index+1]), rev_bg_freqs, pos_bg_freqs, range, num_gc_bins, TRUE)
      : NULL
    );
    pssm_pairs[motif_index] = create_pssm_pair(pos_pssm, neg_pssm);
  }

  // Open the FASTA file for reading.
  FILE* fasta_file = NULL;
  if (open_file(fasta_filename, "r", FALSE, "FASTA", "sequences", &fasta_file) == 0) {
    die("Couldn't open the file %s.\n", fasta_filename);
  }
  if (verbosity >= NORMAL_VERBOSE) {
    if (last == 0) {
      fprintf(stderr, "Using entire sequence\n");
    } else {
      fprintf(stderr, "Limiting sequence to last %d positions.\n", last);
    }
  }

  /**************************************************************
   * Read in all sequences and score with all motifs
   **************************************************************/
  int seq_counter = 0;
  SEQ_T* sequence = NULL;
  while (read_one_fasta(fasta_file, max_seq_length, &sequence)) {
    seq_counter++;

    char* seq_name = get_seq_name(sequence);
    int seq_len = get_seq_length(sequence);
    int scan_len;
    if (last != 0) {
      scan_len = last;
    } else {
      scan_len = seq_len;
    }

    // Get the GC content of the sequence if binning p-values by GC
    // and store it in the sequence object.
    if (num_gc_bins > 1) {
      ARRAY_T *freqs = get_sequence_freqs(sequence);
      set_total_gc_sequence(sequence,
        get_array_item(1,freqs) + get_array_item(2,freqs));	// f(C) + f(G)
      free_array(freqs);			// clean up
    } else {
      set_total_gc_sequence(sequence, -1);	// flag ignore
    }

    /**************************************************************
     * Process all motifs.
     **************************************************************/
    int ns = scan_both_strands ? 2 : 1;
    for (motif_index = 0; motif_index < num_motifs; motif_index++) {
      PATTERN_T *pattern = patterns[motif_index];
      MOTIF_T* motif = &(motifs[ns*motif_index]);
      char* motif_id = get_motif_id(motif);
      char* bare_motif_id = motif_id;
      // Drop the strand indicator from the motif id
      if (*bare_motif_id == '+' || *bare_motif_id == '-') bare_motif_id++;
      if (verbosity >= HIGH_VERBOSE) {
        fprintf(stderr, "Using motif %s of width %d.\n", motif_id, motif->length);
      }
      if ((selected_motifs == NULL) || (have_string(bare_motif_id, selected_motifs) == TRUE)) {
        if (verbosity >= HIGHER_VERBOSE) {
          fprintf(stderr, "Scanning %s sequence with length %d "
              "abbreviated to %d with motif %s with length %d.\n",
              seq_name, seq_len, scan_len, motif_id, motif->length);
        }
        // check if scanned component of sequence has sufficient length for the motif
        if (scan_len < motif->length) {
          fprintf(stderr, "%s too short for motif %s. Skipping!\n", seq_name, motif_id);
          continue;
        }
        // Create a scanned_sequence record and save it in the pattern.
        SCANNED_SEQUENCE_T *scanned_seq = allocate_scanned_sequence(seq_name, seq_name, pattern);
        set_scanned_sequence_length(scanned_seq, scan_len);

        // scan the sequence using average/maximum motif affinity
        ama_sequence_scan(sequence, pssm_pairs[motif_index], scoring, pvalues, last, scanned_seq);
        // Normalize the scores to RMA format if necessary.
        if (TRUE == normalize_scores) {
                int k;
                double maxscore = 1;
                PSSM_T* pssm = pssm_pairs[motif_index]->pos_pssm;
                for (k=0;k<motif->length;k++) {
                        double maxprob = 0;
                        if (maxprob < get_matrix_cell(k, alphabet_index('A', alphabet), pssm->matrix))
                                maxprob = get_matrix_cell(k, alphabet_index('A', alphabet), pssm->matrix);
                        if (maxprob < get_matrix_cell(k, alphabet_index('C', alphabet), pssm->matrix))
                                maxprob = get_matrix_cell(k, alphabet_index('C', alphabet), pssm->matrix);
                        if (maxprob < get_matrix_cell(k, alphabet_index('G', alphabet), pssm->matrix))
                                maxprob = get_matrix_cell(k, alphabet_index('G', alphabet), pssm->matrix);
                        if (maxprob < get_matrix_cell(k, alphabet_index('T', alphabet), pssm->matrix))
                                maxprob = get_matrix_cell(k, alphabet_index('T', alphabet), pssm->matrix);
                        maxscore *= maxprob;
                }
                //fprintf(stderr, "Normalizing score for motif %s on sequence %s. Was %g, now %g\n", motif_id, seq_name,
                //    get_scanned_sequence_score(scanned_seq),get_scanned_sequence_score(scanned_seq)/maxscore );
                set_scanned_sequence_score(scanned_seq, get_scanned_sequence_score(scanned_seq)/maxscore);
        }
      } else {
        if (verbosity >= HIGH_VERBOSE) fprintf(stderr, "Skipping motif %s.\n", motif_id);
      }
    } // All motifs parsed

    free_seq(sequence);
  } // read sequences

  fclose(fasta_file);
  if (verbosity >= HIGH_VERBOSE) fprintf(stderr, "(%d) sequences read in.\n", seq_counter);
  if (verbosity >= NORMAL_VERBOSE) fprintf(stderr, "Finished          \n");

  // output results
  if (output_format == GFF_FORMAT) {
    print_score(program_name, cisml, stdout);
  } else if (output_format == CISML_FORMAT) {
    print_cisml(stdout, cisml, TRUE, NULL, FALSE);
  } else {
    die("Output format invalid!\n");
  }

  /**************************************
   * Clean up.
   **************************************/
  for (motif_index = 0; motif_index < num_motifs; motif_index++) {
    free_motif(&motifs[motif_index]);
    if (scan_both_strands) free_motif(&motifs[num_motifs + motif_index]);
    free_pssm_pair(pssm_pairs[motif_index]);
  }
  free_array(pos_bg_freqs);
  if (scan_both_strands) free_array(rev_bg_freqs);
  free_cisml(cisml);
  free_string_list(selected_motifs);

  // measure time
  if (verbosity >= NORMAL_VERBOSE) { // starting time
    c1 = clock();
    fprintf(stderr, "cycles (CPU);            %ld cycles\n", (long) c1);
    fprintf(stderr, "elapsed CPU time:        %f seconds\n", (float) (c1-c0) / CLOCKS_PER_SEC);
  }
  return (0);
}
