
#include <assert.h>
#include <errno.h>
#include <fnmatch.h>
#include <getopt.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

#include "matrix.h" // includes array.h with special defines set so must be first
#include "array-list.h"
#include "config.h"
#include "dir.h"
#include "fasta-io.h"
#include "io.h"
#include "linked-list.h"
#include "macros.h"
#include "meme-io.h"
#include "motif.h"
#include "projrel.h"
#include "spamo-cisml.h"
#include "spamo-matches.h"
#include "spamo-output.h"
#include "spamo-scan.h"
#include "red-black-tree.h"
#include "utils.h"
#include "xml-out.h"
#include "xml-util.h"

// default file names
#define SPAMO_XML_FILE "spamo.xml"
#define SPAMO_XSL_FILE "spamo-to-html.xsl"
#define SPAMO_HTML_FILE "spamo.html"

// default verbosity
VERBOSE_T verbosity = NORMAL_VERBOSE;

// output a message provided the verbosity is set appropriately
#define DEBUG_MSG(debug_level, debug_msg) { \
  if (verbosity >= debug_level) { \
    fprintf(stderr, debug_msg); \
  } \
}

#define DEBUG_FMT(debug_level, debug_msg_format, ...) { \
  if (verbosity >= debug_level) { \
    fprintf(stderr, debug_msg_format, __VA_ARGS__); \
  } \
}

typedef struct SPAMO_OPTIONS {
  // pseudo random number generator seed
  unsigned int prng_seed;
  // program start time
  time_t start;
  // output directory
  char *outdir;
  // should existing output be overwritten?
  BOOLEAN_T clobber;
  // should EPS files be created?
  BOOLEAN_T create_eps;
  // should PNG files be created?
  BOOLEAN_T create_png;
  // should images be made for redundant motifs?
  BOOLEAN_T create_for_redundant;
  // should the matches be output for motifs declared significant?
  BOOLEAN_T dump_seqs;
  // if the same file is specified for the primary motif and the secondary 
  // motif should the primary be included?
  BOOLEAN_T keep_primary;
  // should the hits be loaded from CISML files?
  BOOLEAN_T load_cismls;
  // pseudocount added to motifs loaded
  double pseudocount;
  // motif background for distributing pseudocounts
  char *bg_filename;
  // maximum distance from primary motif to secondary motif and minimum 
  // distance to edge for primary motif.
  int margin;
  // bin size for calculations (1 is highly recommended)
  int bin;
  // range for significance testing
  int test_range;
  // bit threshold for trimming motifs
  double trim_bit_threshold;
  // fraction of bases needed to be equal before they are considered redundant
  double sequence_similarity;
  // minimum score to be accepted by the scanning code
  double score_threshold;
  // pvalue threshold for spacings to be considered significant.
  double pvalue_cutoff;
  // the minimum overlap in the best site
  int overlap;
  // the fraction of the intersection as relative to the smaller set
  double joint;
  // name of primary motif
  char *primary_name;
  // index of primary motif
  int primary_index;
  // wildcard patterns for motif names to include
  ARRAYLST_T *include_patterns;
  // wildcard patterns for motif names to exclude
  ARRAYLST_T *exclude_patterns;
  // the sequences to be scanned
  char *sequences_file;
  // the primary motif file
  char *primary_file;
  // the files for the secondary motifs
  ARRAYLST_T *secondary_files;
  // a CISML file containing the results of pre-scanning the sequences
  // with the primary motif
  char *primary_cisml;
  // CISML files containing the results of pre-scanning the sequences
  // with the secondary meme files
  ARRAYLST_T *secondary_cismls;
  //DEVELOPER ONLY OPTIONS
  BOOLEAN_T output_peak_sets;
} SPAMO_OPTIONS_T;

/**************************************************************************
 * Frees memory allocated to store the options.
 **************************************************************************/
static void cleanup_options(SPAMO_OPTIONS_T *options) {
  arraylst_destroy(NULL, options->include_patterns);
  arraylst_destroy(NULL, options->exclude_patterns);
  arraylst_destroy(NULL, options->secondary_files);
  arraylst_destroy(NULL, options->secondary_cismls);
}

/**************************************************************************
 * Prints a usage message and exits. 
 * If given an error message it prints that first and will exit with
 * return code of EXIT_FAILURE.
 **************************************************************************/
static void usage(SPAMO_OPTIONS_T *options, char *format, ...) {
  va_list argp;

  char *usage = 
    "Usage:\n" 
    "    spamo [options] <sequences> <primary motif> <secondary motifs>+\n"
    "Options:\n"
    "    -o <directory>   create the directory and write output files in it;\n"
    "                     not compatible with -oc\n"
    "    -oc <directory>  create the directory but if it already exists allow\n"
    "                     overwriting the contents; default: spamo_out\n"
    "    -loadcismls      load CISML files instead of scanning; if this flag\n"
    "                     is specified then each motif must have a CISML file\n"
    "                     specified after it; not compatible with -trim;\n"
    "                     default: scan sequences\n"
    "    -eps             output histograms in eps format; usable with -png\n"
    "    -png             output histograms in png format; usable with -eps\n"
    "    -dumpseqs        dump matching trimmed sequences to output files;\n"
    "                     matches are initially in sequence name order;\n"
    "                     see documentation for column descriptions\n"
    "    -numgen <seed>   specify the random seed for initializing the pseudo-\n"
    "                     random number generator used in breaking ties;\n"
    "                     the seed is included in the output so experiments\n"
    "                     can be repeated; special value 'time' seeds to the\n"
    "                     system clock; default: 1\n"
    "    -margin <size>   edge margin excluded for primary motif matches and\n"
    "                     the maximum distance from the primary motif to the\n"
    "                     secondary motif; default: 150\n"
    "    -bin <size>      size of bins used for output; default: 1\n"
    "    -range <size>    the range from the primary to include in significance\n"
    "                     tests; default: 20\n"
    "    -shared <fract>  fraction of shared trimmed sequence content that\n"
    "                     is required to exclude the sequence as redundant;\n"
    "                     default: 0.5\n"
    "    -cutoff <pvalue> cutoff for spacings to be considered significant;\n"
    "                     default: 0.05\n"
    "    -overlap <size>  number of bases that the most significant spacing\n"
    "                     must overlap before further redundancy testing is\n"
    "                     done; default: 2\n"
    "    -joint <fract>   fraction of sites making up the most significant\n"
    "                     spacing that must be in both sets for the less\n"
    "                     significant motif to be considered redundant;\n "
    "                     default: 0.5\n"
    "    -pseudo <count>  pseudocount added to loaded motifs;\n"
    "                     default: 0.1\n"
    "    -bgfile <file>   file containing background frequency information\n"
    "                     used in applying pseudocounts\n"
    "    -trim <bits>     trim the edges of motifs based on the information\n"
    "                     content; positions on the edges with information\n"
    "                     content less than bits will not be used in\n"
    "                     scanning; ignored when loading CISML files;\n"
    "                     default: 0.25\n"
    "    -primary <name>  name of motif to select as the primary motif;\n"
    "                     overrides -primaryi\n"
    "    -primaryi <num>  index of motif to select as the primary motif\n"
    "                     counting from 1; default: 1\n"
    "    -keepprimary     if the same file is specified for the primary and\n"
    "                     secondary motifs then by default the primary motif\n"
    "                     is excluded but specifying this option keeps it\n"
    "    -inc <pattern>   name pattern to select as secondary motif; may be\n"
    "                     repeated; default: all motifs are used\n"
    "    -exc <pattern>   name pattern to exclude as secondary motif; may be\n"
    "                     repeated; default: all motifs are used\n"
    "    -help            print this usage message\n"
    "    -verbosity <v>   set the verbosity of the output: 1 (quiet) - 5 (dump);\n"
    "                     default: 2 (normal)\n"
    "Description:\n"
    "    SpaMo looks for significant spacings between a primary motif and a\n"
    "    library of secondary motifs.\n";

  if (format) {
    va_start(argp, format);
    vfprintf(stderr, format, argp);
    va_end(argp);
    fprintf(stderr, "\n");
    fputs(usage, stderr);
    fflush(stderr);
  } else {
    puts(usage);
  }
  cleanup_options(options);
  if (format) exit(EXIT_FAILURE);
  exit(EXIT_SUCCESS);
}

// define constants for the arguments that would clash with others if given
// a single character code
#define OUTPUT 500
#define CLOBBER 501
#define CREATE_EPS 502
#define CREATE_PNG 503
#define PRIMARY_INDEX 504
#define PSEUDOCOUNT 505
#define BACKGROUND 506
#define DEV_PEAKS 507

/**************************************************************************
 * Checks the consistency of arguments and loads them into the options 
 * structure for easy access.
 **************************************************************************/
static void process_arguments(int argc, char **argv, SPAMO_OPTIONS_T *options) {
  int option_index = 1;
  BOOLEAN_T bad_argument = FALSE;
  // Note options are specified as follows:
  // <name> <has argument> <(not used)> <int to return>
  struct option spamo_options[] = {
    {"o", required_argument, NULL, OUTPUT},
    {"oc", required_argument, NULL, CLOBBER},
    {"loadcismls", no_argument, NULL, 'l'},
    {"eps", no_argument, NULL, CREATE_EPS},
    {"png", no_argument, NULL, CREATE_PNG},
    {"dumpseqs", no_argument, NULL, 'd'},
    {"numgen", required_argument, NULL, 'n'},
    {"margin", required_argument, NULL, 'm'},
    {"bin", required_argument, NULL, 'b'},
    {"range", required_argument, NULL, 'r'},
    {"shared", required_argument, NULL, 's'},     //sequence similarity threshold
    {"cutoff", required_argument, NULL, 'c'},     //pvalue cutoff
    {"overlap", required_argument, NULL, 'o'},    //redundant site overlap
    {"joint", required_argument, NULL, 'j'},      //redundant joint fraction
    {"pseudocount", required_argument, NULL, PSEUDOCOUNT},
    {"bgfile", required_argument, NULL, BACKGROUND},
    {"trimmotifs", required_argument, NULL, 't'},
    {"primary", required_argument, NULL, 'p'},
    {"primaryi", required_argument, NULL, PRIMARY_INDEX},
    {"keepprimary", no_argument, NULL, 'k'},
    {"include", required_argument, NULL, 'i'},
    {"exclude", required_argument, NULL, 'x'},
    {"help", no_argument, NULL, 'h'},
    {"verbosity", required_argument, NULL, 'v'},
    {"zpeaks", no_argument, NULL, DEV_PEAKS},
    {NULL, 0, NULL, 0} //boundary indicator
  };

  // set option defaults
  options->start = time(NULL);
  options->prng_seed = 1;
  options->clobber = TRUE;
  options->outdir = "spamo_out";
  options->load_cismls = FALSE;
  options->create_eps = FALSE;
  options->create_png = FALSE;
  options->create_for_redundant = TRUE;
  options->dump_seqs = FALSE;
  options->margin = 150;
  options->bin = 1;
  options->test_range = 20;
  options->trim_bit_threshold = 0.25;
  options->sequence_similarity = 0.5;
  options->score_threshold = 7;
  options->pvalue_cutoff = 0.05;
  options->overlap = 2;
  options->joint = 0.5;
  options->primary_name = NULL;
  options->primary_index = 0;
  options->keep_primary = FALSE;
  options->include_patterns = arraylst_create();
  options->exclude_patterns = arraylst_create();
  options->sequences_file = NULL;
  options->primary_file = NULL;
  options->secondary_files = arraylst_create();
  options->primary_cisml = NULL;
  options->secondary_cismls = arraylst_create();
  //meme file loading defaults (same as FIMO)
  options->pseudocount = 0.1;
  options->bg_filename = NULL;
  
  // developer options
  options->output_peak_sets = FALSE;

  // parse optional arguments
  while (1) {
    int opt = getopt_long_only(argc, argv, "", spamo_options, NULL);
    if (opt == -1) break;
    switch (opt) {
      case OUTPUT:        //-o <dir>
        options->clobber = FALSE;
      case CLOBBER:       //-oc <dir>
        options->outdir = optarg;
        break;
      case 'l':           //-loadcismls
        options->load_cismls = TRUE;
        break;
      case CREATE_EPS:    //-eps
        options->create_eps = TRUE;
        break;
      case CREATE_PNG:    //-png
        options->create_png = TRUE;
        break;
      case 'd':           //-dumpseqs
        options->dump_seqs = TRUE;
        break;
      case 'n':           //-numgen <seed>|time
        if (strcmp(optarg, "time") == 0) {
          options->prng_seed = options->start;
        } else {
          options->prng_seed = strtoul(optarg, NULL, 10);
        }
        break;
      case 'm':           //-margin <num>
        options->margin = strtol(optarg, NULL, 10);
        break;
      case 'b':           //-bin <num>
        options->bin = strtol(optarg, NULL, 10);
        break;
      case 'r':           //-range <num>: pvalue testing range
        options->test_range = strtol(optarg, NULL, 10);
        break;
      case 's':           //-shared <fraction>: shared sequence
        options->sequence_similarity = strtod(optarg, NULL);
        break;
      case 'c':           //-cutoff <pvalue>: pvalue significance cutoff
        options->pvalue_cutoff = strtod(optarg, NULL);
        break;
      case 'o':           //-overlap <fraction>: minimum peak overlap
        options->overlap = (int)strtol(optarg, NULL, 10);
        break;
      case 'j':           //-joint <fraction>: minimum shared sites
        options->joint = strtod(optarg, NULL);
        break;
      case PSEUDOCOUNT:   //-pseudocount <fraction>: motif pseudocount
        options->pseudocount = strtod(optarg, NULL);
        if (options->pseudocount < 0) {
          usage(options, "Pseudocount must be positive but got \"%s\".", optarg);
        }
        break;
      case BACKGROUND:    //-bgfile <file>: background to distribute pseudocount
        options->bg_filename = optarg;
        break;
      case 't':           //-trimmotifs <bits>: bit threshold for trimming
        options->trim_bit_threshold = strtod(optarg, NULL);
        break;
      case 'p':           //-primary <name>: name of primary motif
        options->primary_name = optarg;
        break;
      case PRIMARY_INDEX: //-primaryi <index>: index of primary motif
        options->primary_index = (int)strtol(optarg, NULL, 10);
        if (options->primary_index < 1) {
          usage(options, "Primary index must be larger than zero but got \"%s\".", optarg);
        }
        options->primary_index -= 1;
        break;
      case 'k':           //-keepprimary: don't exclude the primary motif from the matched
        options->keep_primary = TRUE;
        break;
      case 'i':           //-include <pattern>
        arraylst_add(optarg, options->include_patterns);
        break;
      case 'x':           //-exclude <pattern>
        arraylst_add(optarg, options->exclude_patterns);
        break;
      case 'h':           //-help
        usage(options, NULL);
        break;
      case 'v':           //-verbosity <1-5>
        verbosity = (int)strtol(optarg, NULL, 10);
        if (verbosity < 1 || verbosity > 5) {
          usage(options, "Verbosity must be between 1 and 5 inclusive. Got \"%s\".", optarg);
        }
        break;
      case DEV_PEAKS:
        options->output_peak_sets = TRUE;
        break;
      case '?':           //unrecognised or ambiguous argument
        bad_argument = TRUE;
    }
  }
  if (bad_argument) usage(options, "One or more unknown or ambiguous options were supplied.");
  option_index = optind;
  // parse required arguments
  // get sequences file
  if (option_index >= argc) usage(options, "No sequences file!");
  options->sequences_file = argv[option_index++];
  if (!file_exists(options->sequences_file))
    usage(options, "Sequences file \"%s\" does not exist!", options->sequences_file);
  // get primary motif file
  if (option_index >= argc) usage(options, "No primary motif file!");
  options->primary_file = argv[option_index++];
  if (!file_exists(options->primary_file))
    usage(options, "Primary motif file \"%s\" does not exist!", options->primary_file);
  // get primary motif CISML file
  if (options->load_cismls) {
    if (option_index >= argc) usage(options, "Expected CISML file to be paired with primary motif file \"%s\".", argv[option_index-1]);
    options->primary_cisml = argv[option_index++];
    if (!file_exists(options->primary_cisml))
      usage(options, "Primary motif CISML file \"%s\" does not exist!", options->primary_cisml);
  }
  // get secondary motifs and CISML files
  if (option_index >= argc) usage(options, "No secondary motif files!");
  while (option_index < argc) {
    // get secondary motif file
    if (!file_exists(argv[option_index]))
      usage(options, "Secondary motif file \"%s\" does not exist!", argv[option_index]);
    arraylst_add(argv[option_index++], options->secondary_files);
    // get secondary motif CISML file
    if (options->load_cismls) {
      if (option_index >= argc) 
        usage(options, "Expected CISML file to be paired with secondary motif file \"%s\".", argv[option_index-1]);
      if (!file_exists(argv[option_index]))
        usage(options, "Secondary motif CISML file \"%s\" does not exist!", argv[option_index]);
      arraylst_add(argv[option_index++], options->secondary_cismls);
    }
  }
}

/**************************************************************************
 * Create the output directory.
 **************************************************************************/
static void create_spamo_output_directory(SPAMO_OPTIONS_T *options) {
  //make output directory
  if (create_output_directory(options->outdir, options->clobber, verbosity >= NORMAL_VERBOSE)) {
    die("Can not continue without an output directory.\n");
  }
}

/**************************************************************************
 * Loads the primary motif using the information specified in the options.
 **************************************************************************/
static void load_primary_motif(SPAMO_OPTIONS_T *options, MOTIF_DB_T **primary_db, MOTIF_T **primary_motif, ARRAY_T **background) {
  BOOLEAN_T has_reverse_strand = FALSE;
  int i;

  DEBUG_MSG(NORMAL_VERBOSE, "Loading Primary Motif\n");

  *primary_db = create_motif_db(0, options->primary_file, options->primary_cisml);

  //load the motifs
  ARRAYLST_T *motifs = arraylst_create();
  read_meme_file2(
      options->primary_file,
      options->bg_filename,
      options->pseudocount,
      REQUIRE_PSPM,
      motifs,
      NULL, //motif occurrances (not used)
      &has_reverse_strand,
      background
  );
  
  (*primary_db)->loaded = arraylst_size(motifs);

  *primary_motif = NULL;
  if (options->primary_name) {
    //select a motif with the requested name
    for (i = 0; i < arraylst_size(motifs); ++i) {
      MOTIF_T *motif = (MOTIF_T*)arraylst_get(i, motifs);
      if (strcmp(get_motif_id(motif), options->primary_name) == 0) {
        *primary_motif = motif;
        ++i;
        break;
      } else {
        //not the motif we want
        free_motif(motif);
        free(motif);
      }
    }
  } else {
    //if the index is impossible make sure it doesn't break the loop
    if (options->primary_index > arraylst_size(motifs)) {
      options->primary_index = arraylst_size(motifs);
    }
    //destroy all motifs before the selected index
    for (i = 0; i < options->primary_index; ++i) {
      MOTIF_T *motif = (MOTIF_T*)arraylst_get(i, motifs);
      free_motif(motif);
      free(motif);
    }
    //get the motif at the selected index
    if (i < arraylst_size(motifs)) *primary_motif = (MOTIF_T*)arraylst_get(i++, motifs);
  }
  //clean up remaining motifs
  for (; i < arraylst_size(motifs); ++i) {
    MOTIF_T *motif = (MOTIF_T*)arraylst_get(i, motifs);
    free_motif(motif);
    free(motif);
  }
  arraylst_destroy(NULL, motifs);
  if (*primary_motif == NULL) die("No acceptable primary motif found.");
  if (!options->load_cismls) trim_motif_by_bit_threshold(*primary_motif, options->trim_bit_threshold);
  (*primary_db)->excluded = (*primary_db)->loaded - 1;
}

/**************************************************************************
 * Loads the sequences, finds the primary motif matches and trims
 * the sequences to only the portion surrounding the best primary match.
 * Also filters unmatched and redundant sequences.
 **************************************************************************/
static void calculate_trimmed_filtered_sequences(
  SPAMO_OPTIONS_T *options, 
  MOTIF_T *primary_motif, 
  ARRAY_T *background, 
  char **trimmed_sequence_data, 
  SEQUENCE_DB_T **sequence_db,
  RBTREE_T **sequences
) {
  RBTREE_T *seqs;
  RBNODE_T *node, *node2, *next;
  SEQUENCE_T *sequence, *sequence2;
  SEQ_T **all_sequences, *current;
  int count, i, trim_len, match, match2, similarity, max_similarity, name_len, tooshort_count, nomatch_count, duplicate_count;
  char *seq_dest, *p1, *p2, *name;
  FILE *fasta_fp;
  BOOLEAN_T created;
  
  *sequence_db = create_sequence_db(options->sequences_file);

  // calculate the trimmed length
  trim_len = 2 * options->margin + get_motif_trimmed_length(primary_motif);
  max_similarity = trim_len * options->sequence_similarity;

  DEBUG_MSG(NORMAL_VERBOSE, "Loading All Sequences\n");
  // read the entire fasta file into memory...
  if (!open_file(options->sequences_file, "r", FALSE, "FASTA sequences", "sequences", &fasta_fp)) exit(EXIT_FAILURE);
  read_many_fastas(fasta_fp, MAX_SEQ, &count, &all_sequences);
  fclose(fasta_fp);
  DEBUG_FMT(NORMAL_VERBOSE, "Loaded %d Sequences\n", count);
  (*sequence_db)->loaded = count;

  tooshort_count = 0;

  // create our set, ensuring no duplicate ids
  seqs = rbtree_create(rbtree_strcmp, rbtree_strcpy, free, NULL, free);
  for (i = 0; i < count; ++i) {
    current = all_sequences[i];
    // check that the sequence is long enough to have at minimum a single site for the primary
    if (get_seq_length(current) < trim_len) {
      tooshort_count++;
      continue;
    }

    node = rbtree_lookup(seqs, get_seq_name(current), TRUE, &created);
    if (!created) die("Duplicate sequence identifier \"%s\" in FASTA file\n", get_seq_name(current));
    sequence = (SEQUENCE_T*)mm_malloc(sizeof(SEQUENCE_T));
    sequence->index = -1;
    sequence->length = get_seq_length(current);
    sequence->data = get_raw_sequence(current);
    sequence->name = get_seq_name(current);
    sequence->primary_match = 0;
    rbtree_set(seqs, node, sequence);
  }
  if (tooshort_count > 0) {
    DEBUG_FMT(QUIET_VERBOSE, "Warning: Eliminated %d sequences that were too "
        "short to be scanned with the current margin (%d) and primary motif "
        "length (%d). Sequences must be at minimum %d long.\n", tooshort_count, 
        options->margin, get_motif_trimmed_length(primary_motif), trim_len);
  }
  (*sequence_db)->excluded_tooshort = tooshort_count;

  DEBUG_MSG(NORMAL_VERBOSE, "Determining Best Primary Matches\n");
  // if the primary has already been scanned then use the CISML input, otherwise we scan it ourselves
  if (options->load_cismls) {
    load_spamo_primary(options->primary_cisml, options->margin, options->score_threshold, seqs, primary_motif);
  } else {
    scan_spamo_primary(options->margin, options->score_threshold, background, primary_motif, seqs);
  }
  DEBUG_MSG(NORMAL_VERBOSE, "Eliminating Sequences Without Primary Matches\n");
  //remove sequences that don't have a match
  nomatch_count = 0;
  node = rbtree_first(seqs);
  while (node != NULL) {
    next = rbtree_next(node);
    sequence = (SEQUENCE_T*)rbtree_value(node);
    if (!(sequence->primary_match)) {
      DEBUG_FMT(DUMP_VERBOSE, "Eliminating \"%s\": no primary match.\n", sequence->name);
      rbtree_delete(seqs, node, NULL, NULL);
      nomatch_count++;
    }
    node = next;
  }
  DEBUG_FMT(NORMAL_VERBOSE, "Eliminated %d Sequences\n", nomatch_count);
  (*sequence_db)->excluded_nomatch = nomatch_count;

  DEBUG_MSG(NORMAL_VERBOSE, "Eliminating Similar Sequences\n");
  duplicate_count = 0;
  //remove sequences that are too similar to each other in the margin around the match
  for (node = rbtree_first(seqs); node != NULL; node = rbtree_next(node)) {
    sequence = (SEQUENCE_T*)rbtree_value(node);
    match = sequence->primary_match;
    if (match < 0) match = -match;
    node2 = rbtree_next(node);
    while (node2 != NULL) {
      next = rbtree_next(node2);
      sequence2 = (SEQUENCE_T*)rbtree_value(node2);
      match2 = sequence->primary_match;
      if (match2 < 0) match2 = -match2;
      p1 = sequence->data+(match - options->margin - 1);
      p2 = sequence2->data+(match2 - options->margin - 1);
      similarity = 0;
      for (i = 0; i < trim_len; ++i, ++p1, ++p2) {
        if (*p1 == *p2) {
          ++similarity;
        } else if ((trim_len - i) < (max_similarity - similarity)) {
          // even if the rest of the sequence is identical it will not be similar enough
          break;
        }
      }
      if (similarity > max_similarity) {
        DEBUG_FMT(DUMP_VERBOSE, "Eliminating \"%s\": %d%% similarity to \"%s\".\n", sequence2->name, 
            (int)(((double)similarity / trim_len) * 100 + 0.5), sequence->name);
        rbtree_delete(seqs, node2, NULL, NULL);
        duplicate_count++;
      }
      node2 = next;
    }
  }
  DEBUG_FMT(NORMAL_VERBOSE, "Eliminated %d Sequences\n", duplicate_count);
  (*sequence_db)->excluded_similar = duplicate_count;

  DEBUG_MSG(NORMAL_VERBOSE, "Trimming Sequences\n");
  //allocate space for the trimmed sequences (in sequential memory)
  *trimmed_sequence_data = mm_malloc(sizeof(char) * rbtree_size(seqs) * (trim_len + 1));
  //trim the sequences down to the margin around the match and allocate an index to each node in the tree
  seq_dest = *trimmed_sequence_data;
  for (node = rbtree_first(seqs), i = 0; node != NULL; node = rbtree_next(node), ++i) {
    sequence = (SEQUENCE_T*)rbtree_value(node);
    //allocate an index
    sequence->index = i;
    //copy the sequence to trim
    match = sequence->primary_match;
    if (match < 0) match = -match;
    memcpy(seq_dest, sequence->data+(match - options->margin - 1), sizeof(char) * trim_len);
    sequence->data = seq_dest;
    sequence->length = trim_len;
    seq_dest+=trim_len;
    seq_dest[0] = '\0';
    seq_dest++;
    //copy the sequence name
    name = sequence->name;
    name_len = strlen(name);
    sequence->name = mm_malloc(sizeof(char) * (name_len + 1));
    memcpy(sequence->name, name, sizeof(char) * name_len);
    sequence->name[name_len] = '\0';
  }
  //now we have the name attribute allocated we must change the free function
  rbtree_alter_value_free(seqs, destroy_sequence);
  //deallocate the SEQ_T sequences
  DEBUG_MSG(NORMAL_VERBOSE, "Freeing Unneeded Sequences\n");
  for (i = 0; i < count; ++i) {
    current = all_sequences[i];
    free_seq(current); 
  }
  free(all_sequences);
  *sequences = seqs;
}

/**************************************************************************
 * Clumps together multiple motifs that have the same warning message.
 * Used by load_secondary_motifs.
 **************************************************************************/
static inline void clump_warning(int *current, int type, char *msg, char *file, char *motif_id) {
  if (*current != type) {//clump warnings of the same type
    if (verbosity >= NORMAL_VERBOSE) {
      //add a new line after a warning of a previous type
      if (*current)  fprintf(stderr, "\n");
      //output the message
      fprintf(stderr, msg, file);
    }
    *current = type;
  }
  if (verbosity >= NORMAL_VERBOSE) 
    fprintf(stderr, " %s", motif_id);
}

/**************************************************************************
 * Loads the secondary motifs filtering them with the information in the
 * options structure.
 **************************************************************************/
#define NO_WARNING 0
#define DUPLICATE_WARNING 1
#define MARGIN_WARNING 2
static void load_secondary_motifs(SPAMO_OPTIONS_T *options, ARRAYLST_T **secondary_dbs, RBTREE_T **secondary_motifs, ARRAY_T **secondary_background) {
  BOOLEAN_T has_reverse_strand, check_primary, keep, created;
  RBNODE_T *node;
  int i, j, k, warn_type; 
  ARRAYLST_T *motif_list;
  SECONDARY_KEY_T key;
  MOTIF_DB_T *db;
  MOTIF_T *motif;
  char *file, *cisml;
  
  DEBUG_MSG(NORMAL_VERBOSE, "Loading Secondary Motifs\n");
  cisml = NULL;
  //load all the motifs
  motif_list = arraylst_create();
  *secondary_dbs = arraylst_create_sized(arraylst_size(options->secondary_files));
  *secondary_motifs = rbtree_create(secondary_key_compare, secondary_key_copy, free, NULL, destroy_secondary_motif);
  *secondary_background = NULL;
  for (i = 0; i < arraylst_size(options->secondary_files); ++i) {//loop over the motif db files
    //deallocate memory for any existing background
    if (*secondary_background) {
      free_array(*secondary_background);
    }
    //reset the warning type for the warning clumping
    warn_type = NO_WARNING;
    //get the file name of the next motif db (and CISML if specified)
    file = (char*)arraylst_get(i, options->secondary_files);
    if(options->load_cismls) cisml = (char*)arraylst_get(i, options->secondary_cismls);
    //create a store of information about this motif database
    db = create_motif_db(i+1, file, cisml);//start the id at 1 because we use 0 for the primary
    arraylst_add(db, *secondary_dbs);
    //update the key to use this database
    key.db_id = db->id;
    //calculate if the primary motif needs to be checked for exclusion in this file
    check_primary = (options->keep_primary == FALSE && strcmp(options->primary_file, file) == 0);
    //load the next set of motifs
    read_meme_file2(
        file,
        options->bg_filename,
        options->pseudocount,
        REQUIRE_PSPM,
        motif_list,
        NULL, //motif occurrances (not used)
        &has_reverse_strand,
        secondary_background
    );
    //update the number of loaded motifs
    db->loaded = arraylst_size(motif_list);
    //as long as we're not loading from cisml then trim the motifs
    if (!options->load_cismls) trim_motifs_by_bit_threshold(motif_list, options->trim_bit_threshold);
    //scan the arraylist
    for (j = arraylst_size(motif_list)-1; j >= 0; j--) {//loop over the motifs
      //get the current motif
      motif = (MOTIF_T*)arraylst_get(j, motif_list);
      //update the key to use this motif
      key.motif_id = get_motif_id(motif);
      //check if it is the primary motif (and should be excluded)
      if (check_primary) {
        if (options->primary_name) {
          if (strcmp(get_motif_id(motif), options->primary_name) == 0) goto destroy_motif;
        } else {
          if (j == options->primary_index) goto destroy_motif;
        }
      }
      //check for inclusion
      keep = (arraylst_size(options->include_patterns) == 0);
      for (k = 0; k < arraylst_size(options->include_patterns); k++) {
        if (fnmatch((char*)arraylst_get(k, options->include_patterns), get_motif_id(motif), 0) == 0) {
          keep = TRUE;
          break;
        }
      }
      if (!keep) goto destroy_motif;
      //check for exclusion
      for (k = 0; k < arraylst_size(options->exclude_patterns); k++) {
        if (fnmatch((char*)arraylst_get(k, options->exclude_patterns), get_motif_id(motif), 0) == 0) {
          goto destroy_motif;
        }
      }
      //check that it fits in the margin
      if ((options->margin - get_motif_trimmed_length(motif) + 1) < 1) {
        clump_warning(&warn_type, MARGIN_WARNING, "Warning, the following "
            "secondary motifs in \"%s\" were excluded as they didn't fit in "
            "the margin:", file, get_motif_id(motif));
        goto destroy_motif;
      }
      //check for duplicate ids
      node = rbtree_lookup(*secondary_motifs, &key, TRUE, &created);
      if (!created) {
        clump_warning(&warn_type, DUPLICATE_WARNING, "Warning, the following "
            "duplicate secondary motifs in \"%s\" were excluded:", file, 
            get_motif_id(motif));
        goto destroy_motif;
      }
      //success!
      rbtree_set(*secondary_motifs, node, create_secondary_motif(options->margin, options->bin, db, motif));
      continue;
destroy_motif:
      //motif is not wanted
      free_motif(motif);
      free(motif);
      db->excluded++;
    }
    arraylst_clear(NULL, motif_list);
    if (warn_type && verbosity >= NORMAL_VERBOSE) fprintf(stderr, "\n");
  }
  arraylst_destroy(NULL, motif_list);
  //check that we found suitable motifs
  if (rbtree_size(*secondary_motifs) == 0) die("No acceptable secondary motifs found.");
}


/**************************************************************************
 * Loads the matches to the secondary motifs from either CISML files or by
 * directly interfacing with FIMO to scan the sequences.
 * Note that as the sequences stored in memory have been trimmed this
 * adjusts the positions in the CISML to be relative to the primary minus
 * the margin so as to match up with the trimmed variant.
 **************************************************************************/
static void load_secondary_matches(
  SPAMO_OPTIONS_T *options, 
  RBTREE_T *sequences, 
  MOTIF_T *primary_motif, 
  ARRAYLST_T *secondary_dbs,
  RBTREE_T *secondary_motifs, 
  ARRAY_T *secondary_background
) {
  RBNODE_T *node;
  MOTIF_DB_T *db;
  SECONDARY_MOTIF_T *smotif;
  int i, count, total, report_step, *matches, *hits, hits_size, test_max, total_tests;
  //calculate the total tests (used for p-value correction later)

  test_max = (int)(options->test_range / options->bin) + (options->test_range % options->bin ? 1 : 0);
  total_tests = calculate_test_count(options->margin, options->bin, test_max, secondary_motifs);

  if (options->load_cismls) {
    DEBUG_MSG(NORMAL_VERBOSE, "Loading Secondary Matches\n");
    //parse each cisml file
    for (i = 0; i < arraylst_size(secondary_dbs); ++i) {
      db = (MOTIF_DB_T*)arraylst_get(i, secondary_dbs);
      DEBUG_FMT(HIGH_VERBOSE, "Loading CISML \"%s\"\n", db->cisml);
      load_spamo_secondary(db->cisml, options->margin, 
          options->score_threshold, options->bin, options->pvalue_cutoff, 
          test_max, total_tests, options->dump_seqs, options->outdir,
          sequences, primary_motif, secondary_motifs, db->id);
    }
    //check everything loaded
    for (node = rbtree_first(secondary_motifs), count = 0; node != NULL; node = rbtree_next(node), count++) {
      smotif = (SECONDARY_MOTIF_T*)rbtree_value(node);
      if (!(smotif->loaded)) {
        die("Missing CISML match data for \"%s\"\n", get_motif_id(smotif->motif));
      }
    }
  } else {
    DEBUG_MSG(NORMAL_VERBOSE, "Scanning For Secondary Matches\n");
    //any motifs without data are scanned with FIMO
    matches = (int*)mm_malloc(sizeof(int) * rbtree_size(sequences));
    hits_size = 4 * options->margin;
    hits = (int*)mm_malloc(sizeof(int) * hits_size);
    total = rbtree_size(secondary_motifs);
    report_step = (int)(((double)total / 100) + 0.5);
    if (report_step == 0) report_step = 1; // mod 0 causes Arithmetic exception
    for (node = rbtree_first(secondary_motifs), count = 0; node != NULL; node = rbtree_next(node), count++) {
      smotif = (SECONDARY_MOTIF_T*)rbtree_value(node);
      if (!(smotif->loaded)) {
        //scan
        scan_spamo_secondary(options->margin, options->score_threshold, secondary_background,
            smotif->motif, sequences, matches, hits, hits_size);
        //process
        process_matches(options->margin, options->bin, options->pvalue_cutoff, 
            test_max, total_tests, primary_motif, sequences, smotif, matches);
        //dump sequences (if requested)
        if (options->dump_seqs && smotif->sig_count > 0) {
          output_sequence_matches(options->outdir, options->margin, sequences, 
              primary_motif, smotif, matches);
        }
      }
      if (count % report_step == 0) {
        DEBUG_FMT(HIGH_VERBOSE, "Scanned %d%% ", (int)(((double)(count + 1) / total) * 100 + 0.5));
      } else {
        DEBUG_MSG(HIGH_VERBOSE, ".");
      }
    }
    DEBUG_MSG(HIGH_VERBOSE, "\n");
    free(matches);
    free(hits);
  }
}

/**************************************************************************
 * Compares two secondary motifs by best pvalue
 **************************************************************************/
int secondary_motif_pvalue_comparator(void *v1, void *v2) {
  SECONDARY_MOTIF_T* smotif1 = (SECONDARY_MOTIF_T*)v1;
  SECONDARY_MOTIF_T* smotif2 = (SECONDARY_MOTIF_T*)v2;
  if (smotif1->sigs->pvalue < smotif2->sigs->pvalue) {
    return -1;
  } else if (smotif1->sigs->pvalue == smotif2->sigs->pvalue) {
    return 0;
  } else {
    return 1;
  }
}

/**************************************************************************
 * Returns a sorted array of the secondary motifs with a significant spacing, 
 * sorted by best spacing pvalue
 **************************************************************************/
LINKLST_T* sort_secondary_motifs(RBTREE_T *secondary_motifs) {
  RBNODE_T *node;
  LINKLST_T *list;
  SECONDARY_MOTIF_T *smotif;
  DEBUG_MSG(NORMAL_VERBOSE, "Sorting Significant Secondary Motifs\n");
  //copy entries with significant match into a list
  list = linklst_create();
  for (node = rbtree_first(secondary_motifs); node != NULL; node = rbtree_next(node)) {
    smotif = (SECONDARY_MOTIF_T*)rbtree_value(node);
    if (smotif->sig_count) linklst_add(smotif, list);
  }
  //sort by pvalue
  linklst_sort(secondary_motif_pvalue_comparator, list);
  return list;
}

/**************************************************************************
 * Calculates the maximum value in a bin to be used in output histograms
 **************************************************************************/
int calculate_bin_max(LINKLST_T* sorted_secondary_motifs) {
  LINK_T *node;
  SECONDARY_MOTIF_T *smotif;
  int binmax; 
  //calculate the maximum histogram count
  binmax = 0;
  for (node = linklst_first(sorted_secondary_motifs); node != NULL; node = linklst_next(node)) {
    smotif = (SECONDARY_MOTIF_T*)linklst_get(node);
    if (smotif->max_in_one_bin > binmax) {
      binmax = smotif->max_in_one_bin;
    }
  }
  return binmax;
}

/**************************************************************************
 * Calculate the maximum overlap of the two best peaks
 **************************************************************************/
static inline int peak_overlap(int bin_size, SECONDARY_MOTIF_T *smot1, SECONDARY_MOTIF_T *smot2) {
  int mot1_closest, mot1_furthest, mot2_closest, mot2_furthest, overlap;

  mot1_closest = smot1->sigs->bin * bin_size + 1;
  mot1_furthest = (smot1->sigs->bin + 1) * bin_size + get_motif_trimmed_length(smot1->motif) - 1;

  mot2_closest = smot2->sigs->bin * bin_size + 1;
  mot2_furthest = (smot2->sigs->bin + 1) * bin_size + get_motif_trimmed_length(smot2->motif) - 1;

  //check for overlap
  if (mot1_furthest < mot2_closest) return FALSE; //no overlap
  if (mot2_furthest < mot1_closest) return FALSE; //no overlap

  if (mot1_closest < mot2_closest) {
    overlap = mot1_furthest - mot2_closest + 1;
  } else {
    overlap = mot2_furthest - mot1_closest + 1;
  }
  return overlap;
}

/**************************************************************************
 * Count the intersection of two sorted sets of numbers.
 **************************************************************************/
static inline int set_intersect_count(int *list1, int count1, int *list2, int count2) {
  int same = 0;
  while (count1 > 0 && count2 > 0) {
    if (*list1 == *list2) {
      ++same;
      ++list1;
      --count1;
      ++list2;
      --count2;
    } else if (*list1 < *list2) {
      ++list1;
      --count1;
    } else {
      ++list2;
      --count2;
    }
  }
  return same;
}

/**************************************************************************
 * Determines if a motif has a result set so similar that it should be
 * declared redundant.
 **************************************************************************/
int is_redundant(SPAMO_OPTIONS_T *options, SECONDARY_MOTIF_T *best, SECONDARY_MOTIF_T *other) {
  int same, overlap, i;
  double joint;
  //calculate the overlap
  //check that they're on the same side
  if (SIDE(best->sigs->quad) != SIDE(other->sigs->quad)) return FALSE;//different side

  //check that the overlap is large enough
  overlap = peak_overlap(options->bin, best, other);
  if (overlap < options->overlap) return FALSE;

  //check that the sequence set is similar enough
  same = set_intersect_count(best->seqs, best->seq_count, other->seqs, other->seq_count);
  if (best->seq_count < other->seq_count) {
    joint = (double)same / best->seq_count;
  } else {
    joint = (double)same / other->seq_count;
  }
  return joint > options->joint;
}

/**************************************************************************
 * Groups secondary motifs
 **************************************************************************/
void group_secondary_motifs(SPAMO_OPTIONS_T *options, LINKLST_T* sorted_secondary_motifs) {
  LINK_T *best_node, *other_node, *next_node;
  SECONDARY_MOTIF_T *other;
  GROUPED_MOTIF_T *group;
  DEBUG_MSG(NORMAL_VERBOSE, "Grouping Redundant Secondary Motifs\n");
  for (best_node = linklst_first(sorted_secondary_motifs); 
      best_node != NULL; best_node = linklst_next(best_node)) {
    group = create_grouped_motif((SECONDARY_MOTIF_T*)linklst_get(best_node));
    other_node = linklst_next(best_node);
    while (other_node != NULL) {
      next_node = linklst_next(other_node);
      other = (SECONDARY_MOTIF_T*)linklst_get(other_node);
      if (is_redundant(options, group->best, other)) {
        linklst_remove(other_node, sorted_secondary_motifs);
        linklst_add(other, group->others);
      }
      other_node = next_node;
    }
    linklst_set(group, best_node);
  }
}


/**************************************************************************
 * Debugging
 * Outputs a grid of similarity information. Currently not finished...
 **************************************************************************/
void output_similarity_grid(SPAMO_OPTIONS_T *options, LINKLST_T* sorted_secondary_motifs) {
  char *sim_file;
  FILE *sim_output;
  LINK_T *node1, *node2;
  SECONDARY_MOTIF_T *smot1, *smot2;
  int i, j;

  //make xml output file
  sim_file = make_path_to_file(options->outdir, "similarity_grid.txt");
  if ((sim_output = fopen(sim_file, "w")) == NULL) {
    die("Could not open similarity file \"%s\" for writing.\n", sim_file);
  }

  //print the position, strand, and extent of the peak, print the set of sequences (by id)
  for (j = 0, node1 = linklst_first(sorted_secondary_motifs); 
      node1 != NULL; j++, node1 = linklst_next(node1)) {
    smot1 = (SECONDARY_MOTIF_T*)linklst_get(node1);
    fprintf(sim_output, "%3d %50s len=%2g side=%s bin=%3d set:", j, get_motif_id(smot1->motif), get_motif_length(smot1->motif), 
        (SIDE(smot1->sigs->quad) ? "d" : "u"), smot1->sigs->bin);
    for (i = 0; i < smot1->seq_count; i++) {
      fprintf(sim_output, " %6d", smot1->seqs[i]); 
    }
    fputs("\n", sim_output);
  }


  for (node1 = linklst_first(sorted_secondary_motifs); 
      node1 != NULL; node1 = linklst_next(node1)) {
    smot1 = (SECONDARY_MOTIF_T*)linklst_get(node1);
    for (node2 = linklst_first(sorted_secondary_motifs);
        node2 != NULL; node2 = linklst_next(node2)) {
      smot2 = (SECONDARY_MOTIF_T*)linklst_get(node2);

    }
  }


  fclose(sim_output);
  free(sim_file);
}

/**************************************************************************
 * Outputs xml for the program parameters/model
 **************************************************************************/
void output_model(FILE *xml_output, int argc, char **argv, SPAMO_OPTIONS_T *options, int bin_max, 
    char *tab, int indent, char **buffer, int *buffer_len) {
  int i;
  output_indent(xml_output, tab, indent);
  fprintf(xml_output, "<model>\n");
  output_indent(xml_output, tab, indent + 1);
  fprintf(xml_output, "<command_line>spamo ");
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
  output_indent(xml_output, tab, indent + 1);
  fprintf(xml_output, "<seed>%u</seed>\n", options->prng_seed);
  output_indent(xml_output, tab, indent + 1);
  fprintf(xml_output, "<margin>%d</margin>\n", options->margin);
  output_indent(xml_output, tab, indent + 1);
  fprintf(xml_output, "<bin_size>%d</bin_size>\n", options->bin);
  output_indent(xml_output, tab, indent + 1);
  fprintf(xml_output, "<bin_pvalue_calc_range>%d</bin_pvalue_calc_range>\n", options->test_range);
  output_indent(xml_output, tab, indent + 1);
  fprintf(xml_output, "<bin_pvalue_cutoff>%g</bin_pvalue_cutoff>\n", options->pvalue_cutoff);
  output_indent(xml_output, tab, indent + 1);
  fprintf(xml_output, "<seq_max_shared_fract>%g</seq_max_shared_fract>\n", options->sequence_similarity);
  output_indent(xml_output, tab, indent + 1);
  fprintf(xml_output, "<seq_min_hit_score>%g</seq_min_hit_score>\n", options->score_threshold);
  output_indent(xml_output, tab, indent + 1);
  fprintf(xml_output, "<redundant_overlap>%d</redundant_overlap>\n", options->overlap);
  output_indent(xml_output, tab, indent + 1);
  fprintf(xml_output, "<redundant_joint>%g</redundant_joint>\n", options->joint);
  output_indent(xml_output, tab, indent + 1);
  fprintf(xml_output, "<motif_pseudocount>%g</motif_pseudocount>\n", options->pseudocount);
  output_indent(xml_output, tab, indent + 1);
  fprintf(xml_output, "<motif_trim>%g</motif_trim>\n", options->trim_bit_threshold);
  output_indent(xml_output, tab, indent + 1);
  fprintf(xml_output, "<bin_max>%d</bin_max>\n", bin_max);
  output_indent(xml_output, tab, indent + 1);
  fprintf(xml_output, "<host>%s</host>\n", hostname());
  output_indent(xml_output, tab, indent + 1);
  fprintf(xml_output, "<when>%s</when>\n", strtok(ctime(&(options->start)),"\n"));
  output_indent(xml_output, tab, indent);
  fprintf(xml_output, "</model>\n");
}




/**************************************************************************
 * Outputs xml for the results
 **************************************************************************/
void output_results(
  int argc, 
  char **argv, 
  SPAMO_OPTIONS_T *options, 
  int bin_max, 
  SEQUENCE_DB_T *sequence_db,
  MOTIF_DB_T *primary_db,
  MOTIF_T *primary_motif, 
  ARRAYLST_T *secondary_dbs,
  LINKLST_T *secondary_motifs, 
  time_t start_time, 
  clock_t start_clock
) {
  FILE *xml_output;
  char *xml_file, *buffer, *tab;
  int buffer_len;
  LINK_T *node;
  GROUPED_MOTIF_T *gmotif;
  int i;
  time_t end_time;
  clock_t end_clock;
  tab = "  ";
  buffer = NULL;
  buffer_len = 0;
  DEBUG_MSG(NORMAL_VERBOSE, "Outputting Results\n");
  //create histograms
  create_histograms(options->margin, options->bin, bin_max, options->pvalue_cutoff, options->outdir, 
      primary_motif, secondary_motifs, options->create_eps, options->create_png, options->create_for_redundant);
  //make xml output file
  xml_file = make_path_to_file(options->outdir, SPAMO_XML_FILE);
  if ((xml_output = fopen(xml_file, "w")) == NULL) {
    die("Could not open output file \"%s\" for writing.\n", xml_file);
  }

  fputs(xml_indicator, xml_output);
  fputs(spamo_dtd, xml_output);
  fputs("<spamo version=\"" VERSION "\" release=\"" ARCHIVE_DATE "\">\n", xml_output);
  output_model(xml_output, argc, argv, options, bin_max, tab, 1, &buffer, &buffer_len);
  fprintf(xml_output, "%s<files>\n", tab);
  output_sequence_database(xml_output, sequence_db, tab, 2, &buffer, &buffer_len);
  output_motif_database(xml_output, primary_db, tab, 2, &buffer, &buffer_len);
  for (i = 0; i < arraylst_size(secondary_dbs); ++i) {
    output_motif_database(xml_output, (MOTIF_DB_T*)arraylst_get(i, secondary_dbs), tab, 2, &buffer, &buffer_len);
  }
  fprintf(xml_output, "%s</files>\n", tab);
  fprintf(xml_output, "%s<primary_motif>\n", tab);
  output_motif(xml_output, 0, primary_motif, tab, 2, &buffer, &buffer_len);
  for (node = linklst_first(secondary_motifs); node != NULL; node = linklst_next(node)) {
    gmotif = (GROUPED_MOTIF_T*)linklst_get(node);
    output_secondary_motif(xml_output, gmotif->best, gmotif->others, options->margin, options->bin, tab, 2, &buffer, &buffer_len);
  }
  fprintf(xml_output, "%s</primary_motif>\n", tab);
  time(&end_time);
  end_clock = clock();
  fprintf(xml_output, "%s<run_time cpu=\"%3.1f\" real=\"%3.1f\"/>\n", tab, ((double)(end_clock - start_clock)) / CLOCKS_PER_SEC, difftime(end_time, start_time));
  fputs("</spamo>\n", xml_output);
  fclose(xml_output);
  if (buffer != NULL) free(buffer);
  buffer = NULL;
  free(xml_file);
  
}

/**************************************************************************
 * Outputs html from xml
 **************************************************************************/
void output_html(SPAMO_OPTIONS_T *options) {
  char *xsl_path, *xml_path, *html_path;
  xml_path = make_path_to_file(options->outdir, SPAMO_XML_FILE);
  xsl_path = make_path_to_file(ETC_DIR, SPAMO_XSL_FILE);
  html_path = make_path_to_file(options->outdir, SPAMO_HTML_FILE);
  if (!file_exists(xsl_path)) {
    fprintf(stderr, "Warning: Can not find SpaMo xml to html stylesheet file \"%s\". Skiping creation of html.\n", xsl_path);
    goto cleanup;
  }
  if (!file_exists(xml_path)) {
    fprintf(stderr, "Warning: Can not find SpaMo xml file \"%s\". Skiping creation of html.\n", xsl_path);
    goto cleanup;
  }
  print_xml_filename_to_filename_using_stylesheet(xml_path, xsl_path, html_path);
cleanup:
  free(xsl_path);
  free(xml_path);
  free(html_path);
}


/**************************************************************************
 * entry point for spamo
 **************************************************************************/
int main(int argc, char **argv) {
  SPAMO_OPTIONS_T options;
  char *trimmed_sequence_data;
  SEQUENCE_DB_T *sequence_db;
  RBTREE_T *sequences;
  MOTIF_DB_T *primary_db;
  MOTIF_T *primary_motif;
  ARRAY_T *primary_background;
  ARRAYLST_T *secondary_dbs;
  RBTREE_T *secondary_motifs;
  ARRAY_T *secondary_background;
  LINKLST_T *sorted_secondary_motifs;
  int bin_max;
  time_t start_time;
  clock_t start_clock;

  time(&start_time);
  start_clock = clock();

  process_arguments(argc, argv, &options);
  
  srand(options.prng_seed);

  create_spamo_output_directory(&options);

  load_primary_motif(&options, &primary_db, &primary_motif, &primary_background);

  calculate_trimmed_filtered_sequences(&options, primary_motif, primary_background, &trimmed_sequence_data, &sequence_db, &sequences);

  load_secondary_motifs(&options, &secondary_dbs, &secondary_motifs, &secondary_background);

  load_secondary_matches(&options, sequences, primary_motif, secondary_dbs, secondary_motifs, secondary_background);

  sorted_secondary_motifs = sort_secondary_motifs(secondary_motifs);

  bin_max = calculate_bin_max(sorted_secondary_motifs);

  //debugging
  if (options.output_peak_sets) {
    output_similarity_grid(&options, sorted_secondary_motifs);
  }

  group_secondary_motifs(&options, sorted_secondary_motifs);

  output_results(argc, argv, &options, bin_max, sequence_db, primary_db, primary_motif, secondary_dbs, sorted_secondary_motifs, start_time, start_clock);
  
  //cleanup
  //free the groupings of secondary motifs (but not the secondary motifs)
  linklst_destroy_all(sorted_secondary_motifs, destroy_grouped_motif);
  //free secondary motifs
  rbtree_destroy(secondary_motifs);
  free_array(secondary_background);
  arraylst_destroy(destroy_motif_db, secondary_dbs);
  //free sequences information (not including the actual trimmed sequences)
  rbtree_destroy(sequences);
  destroy_sequence_db(sequence_db);
  //free trimmed sequences
  free(trimmed_sequence_data);
  //free primary motif
  free_motif(primary_motif);
  free(primary_motif);
  free_array(primary_background);
  destroy_motif_db(primary_db);

  //create the html
  output_html(&options);

  //free the command line options information
  cleanup_options(&options);

  return EXIT_SUCCESS;
}

