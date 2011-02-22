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

#define DEFINE_GLOBALS
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
#include "io.h"
#include "string-list.h"
#include "simple-getopt.h"
#include "array.h"
#include "macros.h"
#include "motif.h"
#include "matrix.h"
#include "gendb.h"
#include "pssm.h"
#include "background.h"
#include "hash_alph.h"
#include "red-black-tree.h"

#include "ama_scan.h"

const int GFF_FORMAT = 0; // GFF file format
const int CISML_FORMAT = 1; // cisML file format
const int DIRECTORY_FORMAT = 2; // all file formats available to a single directory

char* program_name = NULL;

#define DEFAULTDIR "ama_out"

static char *default_output_dirname = DEFAULTDIR;  // default output directory name
static const char *text_filename = "ama.txt";
static const char *cisml_filename = "ama.xml";

// Nucleotide alphabet order as in motif.h
extern char alphabet[];

// As from our ama.h - this is done this way
// to avoid linking issues
VERBOSE_T verbosity = NORMAL_VERBOSE;

void print_score(char* program_name, CISML_T* cisml, FILE* gff_file);

void post_process(CISML_T* cisml, PSSM_PAIR_T** pssm_pairs, BOOLEAN_T normalize_scores);

/*
 * Compare two strings
 */
int str_casecompar(void *v1, void *v2) {
	return strcasecmp((const char*)v1, (const char *)v2);
}

/*
 * Copy an integer
 */
void* int_copy(void *v1) {
	int *out;
	out = mm_malloc(sizeof(int));
	*out = *((int*)v1);
	return (void*)out;
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
  BOOLEAN_T combine_duplicates = FALSE;
  int num_gc_bins = 1;
  int sdbg_order = -1;				// don't use sequence background
  BOOLEAN_T scan_both_strands = TRUE;
  ARRAY_T* pos_bg_freqs = NULL;
  ARRAY_T* rev_bg_freqs = NULL;
  clock_t c0, c1; /* measuring cpu_time */
  CISML_T *cisml;
  char * out_dir = NULL;
  BOOLEAN_T clobber = FALSE;
  int i;
  int last = 0;

  /**********************************************
   * COMMAND LINE PROCESSING
   **********************************************/

  const int num_options = 15;
  cmdoption const motif_scan_options[] = {
    { "max-seq-length", REQUIRED_VALUE },
    { "motif", REQUIRED_VALUE },
    { "motif-pseudo", REQUIRED_VALUE },
    { "rma", NO_VALUE },
    { "pvalues", NO_VALUE },
    { "sdbg", REQUIRED_VALUE },
    { "norc", NO_VALUE },
    { "cs", NO_VALUE },
    { "o-format", REQUIRED_VALUE },
    { "o", REQUIRED_VALUE },
    { "oc", REQUIRED_VALUE },
    { "scoring", REQUIRED_VALUE },
    { "verbosity", REQUIRED_VALUE },
    { "gcbins", REQUIRED_VALUE },
    { "last", REQUIRED_VALUE }
  };

  int option_index = 0;

  // Define the usage message.
  char usage[] = "USAGE: ama [options] <motif file> <sequence file> [<background file>]\n"
    "\n"
    "   Options:\n"
    "     --sdbg <order>\t\t\tUse Markov background model of\n"
    "       \t\t\t\t\torder <order> derived from the sequence\n"
    "       \t\t\t\t\tto compute its likelihood ratios.\n"
    "       \t\t\t\t\tOverrides --pvalues, --gcbins and --rma;\n"
    "       \t\t\t\t\t<background file> is required unless\n"
    "       \t\t\t\t\t--sdbg is given.\n"
    "     --motif <id>\t\t\tUse only the motif identified by <id>.\n"
    "       \t\t\t\t\tThis option may be repeated.\n"
    "     --motif-pseudo <float>\t\tThe value <float> times the background\n"
    "       \t\t\t\t\tfrequency is added to the count of each\n"
    "       \t\t\t\t\tletter when creating the likelihood \n"
    "       \t\t\t\t\tratio matrix (default: %g).\n"
    "     --norc\t\t\t\tDisables the scanning of the reverse\n"
    "       \t\t\t\t\tcomplement strand.\n"
    "     --scoring [avg-odds|max-odds]\tIndicates whether the average or \n"
    "       \t\t\t\t\tthe maximum odds should be calculated\n"
    "       \t\t\t\t\t(default: avg-odds)\n"
    "     --rma\t\t\t\tScale motif scores to the range 0-1.\n"
    "       \t\t\t\t\t(Relative Motif Affinity).\n"
    "       \t\t\t\t\tMotif scores are scaled by the maximum\n"
    "       \t\t\t\t\tscore achievable by that PWM. (default:\n"
    "       \t\t\t\t\tmotif scores are not normalized)\n"
    "     --pvalues\t\t\t\tPrint p-value of avg-odds score in cisml\n"
    "       \t\t\t\t\toutput. Ignored for max-odds scoring.\n"
    "       \t\t\t\t\t(default: p-values are not printed)\n"
    "     --gcbins <bins>\t\t\tCompensate p-values for GC content of\n"
    "       \t\t\t\t\teach sequence using given number of \n"
    "       \t\t\t\t\tGC range bins. Recommended bins: 41.\n"
    "       \t\t\t\t\t(default: p-values are based on\n"
    "       \t\t\t\t\tfrequencies in background file)\n"
    "     --cs\t\t\t\tEnable combining sequences with same\n"
    "       \t\t\t\t\tidentifier by taking the average score\n"
    "       \t\t\t\t\tand the Sidac corrected p-value.\n"
    "     --o-format [gff|cisml]\t\tOutput file format (default: cisml)\n"
    "       \t\t\t\t\tignored if --o or --oc option used\n"
    "     --o <directory>\t\t\tOutput all available formats to\n"
    "       \t\t\t\t\t<directory>; give up if <directory>\n"
    "       \t\t\t\t\texists\n"
    "     --oc <directory>\t\t\tOutput all available formats to\n"
    "       \t\t\t\t\t<directory>; if <directory> exists\n"
    "       \t\t\t\t\toverwrite contents\n"
    "     --verbosity [1|2|3|4]\t\tControls amount of screen output\n"
    "       \t\t\t\t\t(default: %d)\n"
    "     --max-seq-length <int>\t\tSet the maximum length allowed for \n"
    "       \t\t\t\t\tinput sequences. (default: %d)\n"
    "     --last <int>\t\t\tUse only scores of (up to) last <n>\n"
    "       \t\t\t\t\tsequence positions to compute AMA.\n"
    "\n";

  // Parse the command line.
  if (simple_setopt(argc, argv, num_options, motif_scan_options) != NO_ERROR) {
    die("Error processing command line options: option name too long.\n");
  }
    
    BOOLEAN_T setoutputformat = FALSE;
    BOOLEAN_T setoutputdirectory = FALSE;

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
    } else if (strcmp(option_name, "cs") == 0) {
		combine_duplicates = TRUE;
    } else if (strcmp(option_name, "motif") == 0) {
	if (selected_motifs == NULL) {
	  selected_motifs = new_string_list();
	}
	add_string(option_value, selected_motifs);
    } else if (strcmp(option_name, "motif-pseudo") == 0) {
	pseudocount = atof(option_value);
    } else if (strcmp(option_name, "o-format") == 0) {
        if (setoutputdirectory) {
            if (verbosity >= NORMAL_VERBOSE)
                fprintf(stderr, "output directory specified, ignoring --o-format\n");
        } else {
            setoutputformat = TRUE;
            if (strcmp(option_value, "gff") == 0)
                output_format = GFF_FORMAT;
            else if (strcmp(option_value, "cisml") == 0)
                output_format = CISML_FORMAT;
            else {
                if (verbosity >= NORMAL_VERBOSE)
                  fprintf(stderr, "Output format not known. Using standard instead (cisML).\n");
                  output_format = CISML_FORMAT;
            }
        }
    } else if (strcmp(option_name, "o") == 0 || strcmp(option_name, "oc") == 0) {
        setoutputdirectory = TRUE;
        if (setoutputformat) {
            if (verbosity >= NORMAL_VERBOSE)
                fprintf(stderr, "output directory specified, ignoring --o-format\n");
        }
        clobber = strcmp(option_name, "oc") == 0;
        out_dir = (char*) malloc (sizeof(char)*(strlen(option_value)+1));
        strcpy(out_dir, option_value);
        output_format = DIRECTORY_FORMAT;
    } else if (strcmp(option_name, "verbosity") == 0) {
	verbosity = atoi(option_value);
    } else if (strcmp(option_name, "scoring") == 0) {
      if (strcmp(option_value, "max-odds") == 0)
	scoring = MAX_ODDS;
      else if (strcmp(option_value, "avg-odds") == 0)
	scoring = AVG_ODDS;
      else if (strcmp(option_value, "sum-odds") == 0)
	scoring = SUM_ODDS;
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
    } else if (strcmp(option_name, "sdbg") == 0) {
      sdbg_order = atoi(option_value);			// >=0 means use sequence bkg
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

  // --sdbg overrides --pvalues and --gcbins and --rma
  int req_args = 3;
  if (sdbg_order >= 0) {
    pvalues = FALSE;
    normalize_scores = FALSE;
    num_gc_bins = 1;
    req_args = 2;
  }

  // Check all required arguments given
  if (sdbg_order >= 0 && argc > option_index + req_args) {
    die("<background file> cannot be given together with --sdbg.\n");
  } else if (argc != option_index + req_args) {
    fprintf(stderr, usage, pseudocount, verbosity, max_seq_length);
    exit(EXIT_FAILURE);
  }

  // Get required arguments. 
  char* motif_filename = argv[option_index];
  option_index++;
  char* fasta_filename = argv[option_index];
  option_index++;
  char* bg_filename;
  if (req_args == 3) {			// required unless --sdbg given
    bg_filename = argv[option_index];
    option_index++;
  } else {
    bg_filename = "--uniform--";	// So PSSMs will use uniform background;
					// we can multiply them out later.
  }

  // measure time
  c0 = clock();

  // Set up hash tables for computing reverse complement if doing --sdbg
  if (sdbg_order >= 0) setup_hash_alph(DNAB);

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

  //this reads any meme file, xml, txt and html
  read_meme_file(motif_filename, bg_filename, pseudocount,
    REQUIRE_PSPM, //read pspm, not pssm
    &num_motifs, motifs, NULL,//motif occurrences
    &has_reverse_strand, &pos_bg_freqs);

  if (verbosity >= NORMAL_VERBOSE) 
    fprintf(stderr, "Number of motifs in file %d.\n", num_motifs);

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
      build_motif_pssm(
        &(motifs[ns*motif_index]), 
        pos_bg_freqs, 
        pos_bg_freqs, 
        NULL, // Priors not used
        0.0L, // alpha not used
        range, 
        num_gc_bins, 
        TRUE
      );
    PSSM_T* neg_pssm = (scan_both_strands ?
      build_motif_pssm(
        &(motifs[ns*motif_index+1]), 
        rev_bg_freqs, 
        pos_bg_freqs, 
        NULL, // Priors not used
        0.0L, // alpha not used
        range, 
        num_gc_bins, 
        TRUE
      )
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
  int seq_loading_num = 0;  // keeps track on the number of sequences read in total
  int seq_counter = 0;		// holds the index to the seq in the pattern
  int unique_seqs = 0;      // keeps track on the number of unique sequences
  BOOLEAN_T need_postprocessing = FALSE;
  SEQ_T* sequence = NULL;
  RBTREE_T* seq_ids = rbtree_create(str_casecompar,NULL,free,int_copy,free);
  RBNODE_T* seq_node;
  BOOLEAN_T created;
  while (read_one_fasta(fasta_file, max_seq_length, &sequence)) {
    ++seq_loading_num;
	created = FALSE;
    char* seq_name = get_seq_name(sequence);
    int seq_len = get_seq_length(sequence);
    int scan_len;
    if (last != 0) {
      scan_len = last;
    } else {
      scan_len = seq_len;
    }
	  
	// red-black trees are only required if duplicates should be combined
	if (combine_duplicates){
		//lookup seq id and create new entry if required, return sequence index
		char *tmp_id = mm_malloc(strlen(seq_name)+1); // required copy for rb-tree
		strncpy(tmp_id,seq_name,strlen(seq_name)+1);
		seq_node = rbtree_lookup(seq_ids, tmp_id, TRUE, &created);
		if (created) {// assign it a loading number
			rbtree_set(seq_ids, seq_node, &unique_seqs);
			seq_counter = unique_seqs;
			++unique_seqs;
		} else {
			seq_counter = *((int*)rbnode_get(seq_node));
		}
	}
	  
    //
    // Set up sequence-dependent background model and compute
    // log cumulative probability of sequence.
    //
    double *logcumback = NULL;                    // array of log cumulative probs.
    if (sdbg_order >= 0) {
      Resize(logcumback, seq_len+1, double);
      char* raw_seq = get_raw_sequence(sequence);
      char* alpha = get_alphabet(FALSE);
      BOOLEAN rc = FALSE;
      double *a_cp = get_markov_from_sequence(raw_seq, alpha, rc, sdbg_order, 0);
      log_cum_back(raw_seq, a_cp, sdbg_order, logcumback);
      myfree(a_cp);
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
		SCANNED_SEQUENCE_T* scanned_seq = NULL;

		
		if (!combine_duplicates || get_pattern_num_scanned_sequences(pattern) <= seq_counter){
			// Create a scanned_sequence record and save it in the pattern.
			scanned_seq = allocate_scanned_sequence(seq_name, seq_name, pattern);
			set_scanned_sequence_length(scanned_seq, scan_len);
		} else {
			// get existing sequence record
			scanned_seq = get_pattern_scanned_sequences(pattern)[seq_counter];
			set_scanned_sequence_length(scanned_seq, max(scan_len, get_scanned_sequence_length(scanned_seq)));
		}
		
		// check if scanned component of sequence has sufficient length for the motif
		if (scan_len < motif->length) {
			// set score to zero and p-value to 1 if not set yet
			if(!has_scanned_sequence_score(scanned_seq)){
				set_scanned_sequence_score(scanned_seq, 0.0);
			}
			if(pvalues && !has_scanned_sequence_pvalue(scanned_seq)){
				set_scanned_sequence_pvalue(scanned_seq, 1.0);
			} 
			add_scanned_sequence_scanned_element(scanned_seq); 
			if (get_scanned_sequence_num_scanned_positions(scanned_seq) > 0L) need_postprocessing = TRUE;
			if (verbosity >= HIGH_VERBOSE) fprintf(stderr, "%s too short for motif %s. Score set to 0!\n", seq_name, motif_id);
		} else {  
			// scan the sequence using average/maximum motif affinity
			ama_sequence_scan(sequence, logcumback, pssm_pairs[motif_index], scoring, 
							  pvalues, last, scanned_seq, &need_postprocessing);
		}

      } else {
        if (verbosity >= HIGH_VERBOSE) fprintf(stderr, "Skipping motif %s.\n", motif_id);
      }
    } // All motifs parsed

    free_seq(sequence);
    if (sdbg_order >= 0) myfree(logcumback);

  } // read sequences

  fclose(fasta_file);
  if (verbosity >= HIGH_VERBOSE) fprintf(stderr, "(%d) sequences read in.\n", seq_loading_num);
  if (verbosity >= NORMAL_VERBOSE) fprintf(stderr, "Finished          \n");

	
  // if any sequence identifier was multiple times in the sequence set  then
  // postprocess of the data is required
  if (need_postprocessing || normalize_scores){
	  post_process(cisml, pssm_pairs, normalize_scores);
  }
	
  // output results
  if (output_format == DIRECTORY_FORMAT) {
      if (!out_dir)
          die("no directory name specified");
      if (create_output_directory(out_dir, clobber, verbosity > QUIET_VERBOSE)) { // only warn in higher verbose modes
          fprintf(stderr, "failed to create output directory `%s' or already exists\n", out_dir);
          exit(1);
      }
      char *path = make_path_to_file(out_dir, text_filename);
      FILE * text_output, *cisml_output;
      text_output = fopen(path, "w"); //FIXME check for errors: MEME doesn't either and we at least know we have a good directory
      myfree(path);
      path = make_path_to_file(out_dir, cisml_filename);
      cisml_output = fopen(path, "w"); //FIXME check for errors
      myfree(path);
      print_cisml(cisml_output, cisml, TRUE, NULL, FALSE);
      print_score(program_name, cisml, text_output);
      fclose(cisml_output);
      fclose(text_output);
  } else if (output_format == GFF_FORMAT) {
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
  myfree(patterns);
  //destroy sequences
  rbtree_destroy(seq_ids);
	
  // measure time
  if (verbosity >= NORMAL_VERBOSE) { // starting time
    c1 = clock();
    fprintf(stderr, "cycles (CPU);            %ld cycles\n", (long) c1);
    fprintf(stderr, "elapsed CPU time:        %f seconds\n", (float) (c1-c0) / CLOCKS_PER_SEC);
  }
  return (0);
}

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
	  char* pattern_id = get_pattern_accession(patterns[i]);
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
	fprintf(gff_file, "\t%s\n",pattern_id); // Comment
      } // num_seq
    } // num_pattern
  } // pattern > 0
}


/**********************************************************************
 post_process()
 
 adjust/normalize scores and p-values
 **********************************************************************/
void post_process(CISML_T* cisml, PSSM_PAIR_T** pssm_pairs, BOOLEAN_T normalize_scores){
	int m_index, seq_index;
	for (m_index=0; m_index<get_cisml_num_patterns(cisml);++m_index){
		PATTERN_T* pattern = get_cisml_patterns(cisml)[m_index];
		double maxscore = 1;
		
		// FIXME: This should be done to the PSSM, not the individual scores!!!
		// Normalize the scores to RMA format if necessary.
		if (normalize_scores) {
			int k;
			PSSM_T* pssm = pssm_pairs[m_index]->pos_pssm;
			for (k=0;k<pssm->w;k++) {
				double maxprob = -BIG;		// These are scores, not probabilities!!!
				char *letters = "ACGT";		// DNA hard-coded.
				char c;
				while ((c = *(letters++))) {
					double prob = get_matrix_cell(k, alphabet_index(c, alphabet), pssm->matrix);
					if (maxprob < prob) maxprob = prob;
				}
				maxscore *= maxprob;
			}
		}
		
		// adjust each scanned sequence
		for (seq_index=0; seq_index<get_pattern_num_scanned_sequences(pattern);++seq_index){
			SCANNED_SEQUENCE_T* scanned_seq = get_pattern_scanned_sequences(pattern)[seq_index];
			// only adjust scores and p-values if more than one copy was scored
			// num_scanned_positions is (mis-)used in ama to indicate the number of times 
			// a sequence identifier 0occured in the set
			if (get_scanned_sequence_num_scanned_positions(scanned_seq) > 1L){
				// take average score
				if(has_scanned_sequence_score(scanned_seq)){
					double avg_odds = get_scanned_sequence_score(scanned_seq)/get_scanned_sequence_num_scanned_positions(scanned_seq);
					set_scanned_sequence_score(scanned_seq, avg_odds);
				}
				// adjust the minimum p-value for multiple hypothesis testing
				if(has_scanned_sequence_pvalue(scanned_seq)){
					double corr_pvalue = 1.0-pow(1.0-get_scanned_sequence_pvalue(scanned_seq),get_scanned_sequence_num_scanned_positions(scanned_seq));
					set_scanned_sequence_pvalue(scanned_seq, corr_pvalue);
				}
			}
			
			// normalize if requested
			if (normalize_scores) {
				set_scanned_sequence_score(scanned_seq, get_scanned_sequence_score(scanned_seq)/maxscore);
			}
		}
	}
}

