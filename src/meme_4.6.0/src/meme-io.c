/**************************************************************************
 * FILE: meme-io.c
 * CREATE DATE: 3/5/2001
 * AUTHOR: William Stafford Noble
 * PROJECT: MHMM
 * COPYRIGHT: 2001-2008, WSN
 * DESCRIPTION: Read a collection of motifs from a MEME 3.x or 4.x output file.
 **************************************************************************/
#ifdef MAIN
#define DEFINE_GLOBALS
#endif

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>
#include <limits.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "io.h"
#include "meme-io.h"
#include "metameme.h"
#include "mhmm-state.h"
#include "rdb-matrix.h"  // For reading background files
#include "string-list.h"
#include "xml-util.h"
#include "alphabet.h"
#include "hash_alph.h"
#include "buffer.h"

#define NEED_PSSM 2
#define NEED_PSPM 8

/* Maximum allowed width of one input line. */
#define MAX_LINE 10000

/* Initialize global variables associated with SPACER_T. */
char * SPACER_STRS[] = {"invalid", "nrdb", "data"};
int NUM_SPACER_T = 3;

// String that marks the letter frequency section.
char* LETTER_FREQUENCY =
"Background letter frequencies (from";
// tlb; letter frequencies might be from a -bfile
// "Background letter frequencies (from dataset with add-one prior applied):\n";

/* String that marks the beginning of the motif occurrences section. */
/* These have change slightly over time so we maintain               */
/* different version for the sake of compatibility                   */
char* START_OCCURRENCES_1 =
  "<INPUT TYPE = HIDDEN NAME = motif-summary VALUE = \"\n";
char* END_OCCURRENCES_1 = "\">\n";
char* START_OCCURRENCES_2 =
  "<INPUT TYPE = \"HIDDEN\" NAME = \"motif-summary\" VALUE = \"\n";
char* START_OCCURRENCES_3 =
  "<input type=\"hidden\" name=\"motif-summary\" value=\" \n";
char* END_OCCURRENCES_2 = " \">\n";

char* END_OF_TEXT_OCCURRENCES = 
  "----------------------------------------"
  "----------------------------------------";

/* String that marks the end of the motif. */
char* END_OF_MOTIF = "Combined";

#define MAX_XPATH_EXPRESSION 200

/*************************************************************************
 * Character matching functions
 ************************************************************************/
inline static int is_newline(void *config, int letter) {
  return (letter == '\n' || letter == '\r');
}

inline static int is_space_or_ab(void *config, int letter) {
  return (isspace(letter) || letter == '<');
}

inline static int is_blank(void *config, int letter) {
  return (isspace(letter) && letter != '\n' && letter != '\r');
}

inline static int is_quote(void *config, int letter) {
  return (letter == '"');
}

inline static int is_quote_or_nl(void *config, int letter) {
  return (letter == '"' || letter == '\n' || letter == '\r');
}

/*************************************************************************
 * Buffer functions
 ************************************************************************/

/*
 * find_bm_landmarks
 * searches for one of the boyer-moore strings and returns the index of the first one found.
 * if require is true then it must find one before the end of the file.
 */
static int find_bm_landmarks(FILE *fp, BUF_T *buf, BOOLEAN_T require, BOOLEAN_T autoskip, int landmarks, BMSTR_T **loc) {
  int found;
  if (landmarks <= 0) die("Bad number (%d) of landmarks specified.\n", landmarks);
  found = buf_fread_until(buf, fp, NULL, landmarks, loc);
  if (found == -1) {
    die("File error while finding %d landmark(s) \"%s\", error was %s\n", landmarks, bmstr_text(loc[0]), strerror(ferror(fp)));
  } else if (found == 0) {
    if (require)
      die("Error did not find any of the %d landmark(s) \"%s\" before reaching the end of the file\n", landmarks, bmstr_text(loc[0]));
    return 0;
  } else if (found > 0 && found <= landmarks) {
    if (autoskip) buf_skip(buf, bmstr_length(loc[found-1]));
    return found;
  } else {
    die("Unexpected return of %d from buf_fread_until.\n", found);
  }
  exit(1); //this is to make the compiler happy, it thinks that this is reachable, but it isn't as die doen't return.
}

/*
 * find_landmarks
 * searches for one of the strings and returns the index of the first one found.
 * if require is true then it must find one before the end of the file.
 */
static int find_landmarks(FILE *fp, BUF_T *buf, BOOLEAN_T require, BOOLEAN_T autoskip, int landmarks, ...) {
  BMSTR_T *loc[landmarks];
  int i, found;
  va_list landmark_lst;
  va_start(landmark_lst, landmarks);
  for (i = 0; i < landmarks; ++i) loc[i] = bmstr_create(va_arg(landmark_lst, char*));
  va_end(landmark_lst);
  found = find_bm_landmarks(fp, buf, require, autoskip, landmarks, loc);
  for (i = 0; i < landmarks; ++i) bmstr_destroy(loc[i]);
  return found;
}

/*
 * look_for
 * looks for an expected string and returns false if it's not found 
 */
static BOOLEAN_T look_for(FILE *fp, BUF_T *buf, char *expect) {
  int value;
  value = buf_fread_unexpected(buf, fp, expect, FALSE);
  if (value == -1) {
    die("File error while expecting \"%s\", error was %s\n", expect, strerror(ferror(fp)));
  }
  return (value != 1);
}

/*
 * expect
 * checks that the following characters match the expected value.
 */
static void expect(FILE *fp, BUF_T *buf, char *expect) {
  int value;
  value = buf_fread_unexpected(buf, fp, expect, FALSE);
  if (value == -1) {
    die("File error while expecting \"%s\", error was %s\n", expect, strerror(ferror(fp)));
  } else if (value == 1) {
    die("Error did not read the expected \"%s\"\n", expect);
  }
}

/*
 * consume
 * wrapper for buf_fread_consume that handles file errors and allows specification of the 
 * minimum expected delimiters to consume
 */
static int consume(FILE *fp, BUF_T *buf, int (*is_delim) (void*, int), void *config, int negate_is_delim, int minlen, char *desc) {
  int count;
  count = buf_fread_consume(buf, fp, is_delim, config, negate_is_delim);
  if (count == -1)
    die("File error reading %s, error was %s\n", desc, strerror(ferror(fp)));
  if (count < minlen)
    die("Error reading %s, too few delimiters to consume, expected %d, got %d\n", minlen, count);
  return count;
}

/*
 * parse_word
 * gets the next word as delimited by spaces. Optionally consumes the blanks before the first word.
 */
static char* parse_word(FILE *fp, BUF_T *buf, char *dest, int destlen, int minlen, BOOLEAN_T consume_space, char *desc) {
  char *word;
  int len;
  if (consume_space) {
    if (buf_fread_consume(buf, fp, is_blank, NULL, FALSE) == -1)
      die("File error reading %s, error was %s\n", desc, strerror(ferror(fp)));
  }
  word = buf_fread_token(buf, fp, buf_is_delim_space, NULL, FALSE, dest, destlen, &len);
  if (word == NULL) {
    if (len == 0) {
      die("File error reading %s, error was %s\n", desc, strerror(ferror(fp)));
    } else {
      die("Error reading %s, exceeded buffer\n", desc);
    }
  }
  if (len < minlen) {
    die("Error reading %s, word too short expected a minimum length of %d\n", desc, minlen);
  }
  return word;
}

/*
 * parse_pos_int
 * gets the next series of digits and converts them to a positive integer
 */
static int parse_pos_int(FILE *fp, BUF_T *buf, char *desc) {
  char tokbuf[20], *token, *end;
  int len;
  long int value;
  token = buf_fread_token(buf, fp, buf_is_delim_digit, NULL, TRUE, tokbuf, 20, &len);
  if (token == NULL) {
    if (len == 0) {
      die("File error reading %s, error was %s\n", desc, strerror(ferror(fp)));
    } else {
      die("Error reading %s, exceeded buffer\n", desc);
    }
  } else if(len == 0) {
    die("Error reading %s, no digits found\n", desc);
  }
  value = strtol(token, &end, 10);
  if (value >= INT_MAX) {
    die("Error reading %s, number too big to convert\n", desc);
  }
  return (int)value;
}

/*
 * parse_double
 * gets the next space delimited token and converts it to a double
 */
static double parse_double(FILE *fp, BUF_T *buf, char *desc) {
  char tokbuf[30], *token, *end;
  int len;
  double value;
  token = buf_fread_token(buf, fp, buf_is_delim_space, NULL, FALSE, tokbuf, 30, &len);
  if (token == NULL) {
    if (len == 0) {
      die("File error reading %s, error was %s\n", desc, strerror(ferror(fp)));
    } else {
      die("Error reading %s, exceeded buffer\n", desc);
    }
  } else if(len == 0) {
      die("Error reading %s, no number found\n", desc);
  }
  value = strtod(token, &end);
  if (*end != '\0') {
    die("Error reading %s, some of token \"%s\" wasn't number\n", desc, token);
  }
  if (errno == ERANGE) {
    //Error range is delibrately ignored as there are very small E and P values that can't be parsed
    errno = 0;
    //die("Error reading %s, number \"%s\" out of convertable range\n", desc, token);
  }
  return value;
}

/*
 * parse_prob
 * gets the next space delimited token skiping leading space and converts it to a probability
 */
static ATYPE parse_prob(FILE *fp, BUF_T *buf, char *desc) {
  char tokbuf[30], *token;
  int len;
  ATYPE value;
  //consume leading spaces
  if (buf_fread_consume(buf, fp, buf_is_delim_space, NULL, FALSE) == -1)
      die("File error reading %s, error was %s\n", desc, strerror(ferror(fp)));
  //get token
  token = buf_fread_token(buf, fp, buf_is_delim_space, NULL, FALSE, tokbuf, 30, &len);
  if (token == NULL) {
    if (len == 0) {
      die("File error reading %s, error was %s\n", desc, strerror(ferror(fp)));
    } else {
      die("Error reading %s, exceeded buffer\n", desc);
    }
  } else if(len == 0) {
      die("Error reading %s, no probability found\n", desc);
  }
  if (sscanf(token, PROB_SCAN, &value) != 1) {
    die("Error reading %s, could not convert probability\n", desc);
  }
  return value;
}

/**************************************************************************
 * Replace the elements an array of frequences with the average
 * over complementary bases.
 **************************************************************************/
void average_freq_with_complement(ARRAY_T *freqs) {

  assert(which_alphabet() == DNA_ALPH);

  char *alphabet = get_alphabet(FALSE);
  int a_index = alphabet_index('A', alphabet);
  int t_index = alphabet_index('T', alphabet);
  int g_index = alphabet_index('G', alphabet);
  int c_index = alphabet_index('C', alphabet);
  double at_freq =
    (get_array_item(a_index, freqs) + get_array_item(t_index, freqs)) / 2.0;
  double gc_freq =
    (get_array_item(g_index, freqs) + get_array_item(c_index, freqs)) / 2.0;
  set_array_item(a_index, at_freq, freqs);
  set_array_item(t_index, at_freq, freqs);
  set_array_item(g_index, gc_freq, freqs);
  set_array_item(c_index, gc_freq, freqs);

}

/**************************************************************************
 * Read an emission distribution from a file in MEME -bfile format.
 **************************************************************************/
ARRAY_T *read_background_file(
  const char* bg_filename      // Name of the file to read from.
)
{
  int i, j, alen1, alen2;
  char *alpha1;        // alphabet in file
  char *alpha2;        // target alphabet
  RDB_MATRIX_T *rdb_matrix;
  FILE *bg_file;       /* Pointer to the background file. */
  ARRAY_T *background;      // Array of probs to return.

  /* Open the file for reading. */
  if (open_file(bg_filename, "r", FALSE, "background", "frequencies", &bg_file)
      == 0)
    exit(1);

  // Read in the background file.
  rdb_matrix = read_rdb_matrix(" \t", FALSE, 1, FALSE, NULL, bg_file);

  /* get alphabet from the row names; discard tuples */
  alen1 = get_num_strings(rdb_matrix->row_names);
  alpha1 = (char *)mm_malloc(sizeof(char) * (alen1 + 1));
  for (i=j=0; i<alen1; i++) {
    char *wa = get_nth_string(i, rdb_matrix->row_names);
    if (strlen(wa) > 1) continue;  // ignore tuples
    alpha1[j++] = islower((int)wa[0]) ? toupper(wa[0]) : wa[0];
  }
  alpha1[j] = '\0';
  alen1 = strlen(alpha1);    // length of alphabet in file

  // Get the target alphabet without ambigs.
  alpha2 = get_alphabet(FALSE);
  alen2 = get_alph_size(ALPH_SIZE);  // length of target alphabet

  // Check that alphabets are the same length.
  if (alen1 != alen2)
    die("The -bg alphabet %s should be %s.\n", alpha1, alpha2);

  /* Allocate the background array. */
  background = allocate_array(get_alph_size(ALL_SIZE));

  // Reorder the probabilities in target alphabet order.
  for (i=0; i<alen2; i++) {
    int ii = strchr(alpha1, alpha2[i]) - alpha1;
    if (ii < 0 || ii >= alen1)
      die("The -bg alphabet %s ought to be %s.\n", alpha1, alpha2);
    set_array_item(i, get_matrix_cell(ii, 0, rdb_matrix->matrix), background);
  }

  /* Extend the distribution to account for ambiguous characters. */
  fill_in_ambiguous_chars(FALSE, background);

  /* Close the file. */
  fclose(bg_file);

  free_rdb_matrix(rdb_matrix);
  myfree(alpha1);

  return(background);
}

/***************************************************************************
 * Set a background distribution
 *  - by reading values from a file if filename is given, or
 *  - equal to the NRDB frequencies if filename is NULL.
 *
 ***************************************************************************/
ARRAY_T* get_background(char* bg_filename) {
  ARRAY_T* background;

  if ((bg_filename == NULL) || (strcmp(bg_filename, "nrdb") == 0)) {
    background = allocate_array(get_alph_size(ALL_SIZE));
    get_nrdb_frequencies(background);
    fill_in_ambiguous_chars(FALSE, background);
  } else {
    background = read_background_file(bg_filename);
  }

  if (verbosity > NORMAL_VERBOSE) {
    fprintf(stderr, "Background distribution: ");
    print_array(background, 5, 3, TRUE, stderr);
  }

  return(background);
}

/***********************************************************************
 * Read the version information from a MEME file.
 * Also looks for a <html> tag. If it finds the version before the html tag
 * then it assumes that it isn't html.
 * Accept only version 3 or 4.
 ***********************************************************************/
void read_version(FILE* motif_file, BUF_T *buffer, int* major, int* minor, int* bug, int* is_html)
{
  int found;
  BMSTR_T *vpats[3];
  //init output vars
  *major = -1;
  *minor = 0;
  *bug = 0;
  *is_html = 0;
  //setup search patterns
  vpats[0] = bmstr_create2("<html>", TRUE);//ignore case
  vpats[1] = bmstr_create("MEME version ");
  vpats[2] = bmstr_create("Meta-MEME version ");
  
  //look for the html or version strings
  found = buf_fread_until(buffer, motif_file, NULL, 3, vpats);
  if (found == 1) {//found html string
    *is_html = 1;
    //continue search for only version strings
    buf_skip(buffer, bmstr_length(vpats[0]));
    found = buf_fread_until(buffer, motif_file, NULL, 2, vpats+1);
    if (found > 0) ++found;
  }
  if (found >= 2) {//found a version string
    buf_skip(buffer, bmstr_length(vpats[found-1]));
    *major = parse_pos_int(motif_file, buffer, "major version number");
    if (look_for(motif_file, buffer, ".")) {
      *minor = parse_pos_int(motif_file, buffer, "minor version number");
      if (look_for(motif_file, buffer, ".")) {
        *bug = parse_pos_int(motif_file, buffer, "bugfix version number");
      }
    }
  } else {
    die("Can't find MEME version.\n");
  }
  if (*major == 1 || *major == 2) {
    die("Sorry, MEME version %d output files are no longer supported.\n",
      *major);
  } else if (*major != 3 && *major != 4) {
    die("Unknown MEME version: %d.%d.%d\n", *major, *minor, *bug);
  }
  bmstr_destroy(vpats[0]);
  bmstr_destroy(vpats[1]);
  bmstr_destroy(vpats[2]);
}

/***********************************************************************
 * Read the version information from MEME XML.
 * Accept only version 3.
 * Caller is responsible for freeing returned xmlChar
 ***********************************************************************/
static xmlChar* read_version_from_xml(xmlDocPtr meme_doc) {

  xmlXPathContextPtr xpathCtx = NULL;
  xmlXPathObjectPtr xpathObj = NULL;
  xmlChar* version = NULL;

  xpathCtx = xmlXPathNewContext(meme_doc);
  if(xpathCtx == NULL) {
    die("Error: unable to create new XPath context while reading version.\n");
  }
  xpathObj = xmlXPathEvalExpression(BAD_CAST "/MEME", xpathCtx);
  if(xpathObj == NULL) {
      die("Error: unable to evaluate xpath expression /MEME.\n");
  }
  version = xmlGetProp(xpathObj->nodesetval->nodeTab[0], BAD_CAST "version");
  if (version == NULL) {
    die("Error: missing version attribute in MEME.\n");
  }
  fprintf(stderr, "MEME version is %s\n", version);

  return(version);
}

/***********************************************************************
 * read_background
 * checks the source of the background frequencies and if it's the 
 * default values or from an external file then load the background 
 * frequencies and return 1. If the background frequencies must be 
 * loaded from the meme file then return 0.
 ***********************************************************************/
static int read_background(
  const char *bg_filename, 
  const char *meme_filename, 
  ARRAY_T** background) 
{
  // Establish the background frequencies.
  if (bg_filename == NULL) {
    // Default is to use pre-calculated bg freq. from NR database.
    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(
        stderr,
        "Using background frequencies from NR sequence database.\n"
      );
    }
    *background = allocate_array(get_alph_size(ALL_SIZE));
    get_nrdb_frequencies(*background);
  } else if (strcmp(bg_filename, "--uniform--") == 0) {
    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(
        stderr,
        "Using uniform background frequencies.\n"
      );
    }
    *background = allocate_array(get_alph_size(ALL_SIZE));
    get_uniform_frequencies(*background);
  } else {
    if (strcmp(bg_filename, "motif-file") == 0) {
      // If bg_filename matches "motif-file" read bg freq. from motif file.
      if (verbosity >= NORMAL_VERBOSE) {
        fprintf(
          stderr,
          "Using background frequencies from file %s.\n",
          meme_filename
        );
      }
      return 0;
    } else {
      // Otherwise read bg freqs. from external bg file.
      if (verbosity >= NORMAL_VERBOSE) {
        fprintf(
          stderr,
          "Using background frequencies from file %s.\n",
          bg_filename
        );
      }
      *background = read_background_file(bg_filename);;
    }
  }
  return 1;
}

/***********************************************************************
 * Read the alphabet from a file and allocate the given array.
 ***********************************************************************/
static void read_alphabet
  (FILE *motif_file, BUF_T *buffer)
{
  char *alphabet;
  find_landmarks(motif_file, buffer, TRUE, TRUE, 1, "ALPHABET= ");
  //consume space
  consume(motif_file, buffer, is_blank, NULL, FALSE, 0, "whitespace between alphabet landmark and alphabet"); 
  //read all a-z or A-Z letters
  alphabet = buf_fread_token(buffer, motif_file, buf_is_delim_alpha,
      NULL, TRUE, NULL, 100, NULL);
  if (alphabet == NULL) 
    die("Found alphabet but it's larger than 100 "
      "chars which can't match any supported alphabet\n");
  //as long as we find something smaller than 100 chars try to use it...
  set_alphabet(verbosity, alphabet);
  free(alphabet);
}

/***********************************************************************
 * Read a MEME file to find out whether both strands are included.
 ***********************************************************************/
static BOOLEAN_T read_strand
  (FILE *motif_file, BUF_T *buffer)
{
  find_landmarks(motif_file, buffer, TRUE, TRUE, 1, "strands: ");
  //remove any whitespace
  consume(motif_file, buffer, is_blank, NULL, FALSE, 0, "whitespace between strand landmark and strand specifier"); 
  //expect a plus
  expect(motif_file, buffer, "+");
  //remove any whitespace
  consume(motif_file, buffer, is_blank, NULL, FALSE, 0, "whitespace between strand specifier + and strand specifier -"); 
  //look for a minus
  return look_for(motif_file, buffer, "-");
}

/***********************************************************************
 * Read MEME XML to find out whether both strands are included.
 ***********************************************************************/
static BOOLEAN_T read_strand_from_xml(xmlXPathContextPtr xpath_ctxt) {

  xmlXPathObjectPtr xpath_obj = NULL;
  xmlChar* strand = NULL;
  BOOLEAN_T result = FALSE;

  xpath_obj = xpath_query(xpath_ctxt, "/MEME/model/strands");
  strand = xmlXPathCastNodeToString(xpath_obj->nodesetval->nodeTab[0]);
  if (strncmp("both", (const char *) strand, 5) == 0) {
    result = TRUE;
  }
  xmlFree(strand);
  xmlXPathFreeObject(xpath_obj);
  return result;
}

/***********************************************************************
 * Read a set of letter frequencies from MEME output.
 ***********************************************************************/
static ARRAY_T* read_bg_freqs
  (FILE *infile, BUF_T *buffer)      // An open file in MEME output format.
{
  int i, len, alph_size;
  char letter[2], number[30];
  ARRAY_T *frequencies;
  ATYPE value;

  //now we're at the start of the data
  frequencies = allocate_array(get_alph_size(ALL_SIZE));
  alph_size = get_alph_size(ALPH_SIZE);
  for (i = 0; i < alph_size; ++i) {
    //consume space
    consume(infile, buffer, buf_is_delim_space, NULL, FALSE, 0, "whitespace before letters in background frequency"); 
    //check the letter
    if (parse_word(infile, buffer, letter, 2, 1, FALSE, "background letter")[0] != get_alph_char(i)) {
      die("Error reading frequencies for %c.\n", get_alph_char(i));
    }
    //consume space
    consume(infile, buffer, is_blank, NULL, FALSE, 1, "whitespace between letter and background frequency"); 
    //set the value
    set_array_item(i, parse_prob(infile, buffer, "background letter frequency"), frequencies);
  }

  // Make sure the frequencies add up to 1.0. 
  normalize_subarray(0, alph_size, 0.0, frequencies);

  // Fill in ambiguous characters. 
  fill_in_ambiguous_chars(FALSE, frequencies);
  return frequencies;
}


/***********************************************************************
 * Read a motif's position specific probability matrix otherwise known
 * as the letter frequencies
 ***********************************************************************/
static void read_pspm(
  FILE *    infile,     /* An open file in MEME output format with the
                           pointer before the requested motif. IN */
  BUF_T *buffer,        /* buffer for reading the MEME file */
  ARRAY_T* background,  /* Pointer to array of bg freqs. */
  double pseudocount,   /* pseudocount IN */
  int loading_mode,     /* should the PSSM or PSPM be loaded? */
  MOTIF_T * a_motif     /* The motif. OUT */
) {
  int found, alph_size, is_oldmethod, i, j;
  ARRAY_T* these_freqs;
  char *token, tokbuf[30];
  
  //set defaults, allows detection of parsing problems
  alph_size = INT_MIN;

  //read each token on the line
  while ((token = parse_word(infile, buffer, tokbuf, 30, 0, TRUE, "pspm header"))[0] != '\0') {
    consume(infile, buffer, is_blank, NULL, FALSE, 1, "spacing"); 
    if (strcmp(token, "alength=") == 0) {
      alph_size = parse_pos_int(infile, buffer, "alphabet size");
    } else if (strcmp(token, "w=") == 0) {
      a_motif->length = parse_pos_int(infile, buffer, "motif width");
    } else if (strcmp(token, "E=") == 0) {
      a_motif->evalue = parse_double(infile, buffer, "motif E-value");
    } else if (strcmp(token, "nsites=") == 0) {
      a_motif->num_sites = parse_double(infile, buffer, "motif sites");
    } else if (strcmp(token, "n=") == 0) {
      // Issue a warning.
      static BOOLEAN_T first_time = TRUE;
      if (first_time) {
        fprintf(stderr, "\nWarning: This is an old MEME file that contains\n"
                        "posterior probabilities rather than frequencies.\n"
                        "Meta-MEME will still work, but it would be better\n"
                        "to run an updated version of MEME on your data.\n");
        first_time = FALSE;
      }
      a_motif->num_sites = parse_double(infile, buffer, "motif sites");
    }
  }
  //check that everything we need was loaded and the values are sane
  if (alph_size == INT_MIN) die("Alphabet size not parsed\n");
  if (a_motif->length == INT_MIN) die("Motif length not parsed\n");
  //test optional values
  if (a_motif->evalue == -1) {
    if (loading_mode & REQUIRE_EVALUE) 
      die("Motif %s E-value not parsed\n", a_motif->id);
  } else if (a_motif->evalue < 0) {
    die("\nThe record for motif %s indicated the Evalue was %lf.\n"
      "The E-value for a motif must be a positive number.\n",
      a_motif->id,
      a_motif->evalue);
  }

  if (a_motif->num_sites == -1) {
    if (loading_mode & REQUIRE_NSITES) 
      die("Motif %s number of sites not parsed\n", a_motif->id);
  } else if (a_motif->num_sites <= 0.0) {
    die("\nThe record for motif %s indicated the number of sites was %lf.\n"
      "The number of sites for a motif must be a positive number.\n",
      a_motif->id,
      a_motif->num_sites);
  }

  // Allocate memory for the matrix. 
  a_motif->freqs = allocate_matrix(a_motif->length, get_alph_size(ALL_SIZE));

  // Read the letter frequency matrix. 
  for (i = 0; i < a_motif->length; i++) {
    these_freqs = get_matrix_row(i, a_motif->freqs);

    double row_sum = 0.0;
    for (j = 0; j < alph_size; j++) {
      double letter_freq = 0.0;
      double bg_letter_freq = get_array_item(j, background);
      letter_freq = parse_prob(infile, buffer, "letter frequency");
      // Adjust freq using pseudocount
      letter_freq =
        ((a_motif->num_sites * letter_freq) + pseudocount * bg_letter_freq)
        / (a_motif->num_sites  + pseudocount);
      row_sum += letter_freq;
      set_array_item(j, letter_freq, these_freqs);
    }

    // Normalize the first alph_size positions. (MEME prints six digits of precision.) 
    normalize_subarray(0, alph_size, 0.00001, these_freqs);

    // Compute values for ambiguous characters. 
    fill_in_ambiguous_chars(FALSE, these_freqs);
  }

}

/***********************************************************************
 * Read a motif's position specific scoring matrix otherwise known
 * as the letter scores
 ***********************************************************************/
static void read_pssm(
  FILE *    infile,     /* An open file in MEME output format with the
                           pointer before the requested motif. IN */
  BUF_T *buffer,        /* buffer for reading the MEME file */
  int   loading_mode,   /* Should PSSM or PSPM be loaded? */
  MOTIF_T * a_motif     /* The motif. OUT */
) {
  int found, alph_size, is_pair, i, j;
  ARRAY_T* these_freqs;
  char *token, tokbuf[30];

  //set defaults, allows detection of parsing problems
  alph_size = INT_MIN;
  a_motif->length = INT_MIN;
  a_motif->evalue = -1;

  //read each token on the line
  while ((token = parse_word(infile, buffer, tokbuf, 30, 0, TRUE, "pssm header"))[0] != '\0') {
    consume(infile, buffer, is_blank, NULL, FALSE, 1, "spacing"); 
    if (strcmp(token, "alength=") == 0) {
      alph_size = parse_pos_int(infile, buffer, "alphabet size");
    } else if (strcmp(token, "w=") == 0) {
      a_motif->length = parse_pos_int(infile, buffer, "motif width");
    } else if (strcmp(token, "E=") == 0) {
      a_motif->evalue = parse_double(infile, buffer, "motif E-value");
    }
  }
  if (alph_size == INT_MIN) die("Alphabet size not parsed\n");
  if (a_motif->length == INT_MIN) die("Motif length not parsed\n");
  if (a_motif->evalue == -1) {
    if (loading_mode & REQUIRE_EVALUE) die("Motif evalue not parsed\n");
  } else if (a_motif->evalue < 0) {
    die("\nThe record for motif %s indicated the evalue of the motif was %f.\n"
      "The width of a motif must be a positive number.\n",
      a_motif->id,
      a_motif->evalue);
  }
  //ensure it is initialized
  a_motif->scores = NULL;
  //Does the caller actually want the pssm?
  if (!(loading_mode & WANT_PSSM)) return;

  // Allocate memory for the matrix. 
  a_motif->scores = allocate_matrix(a_motif->length, alph_size);

  // Read the letter frequency matrix. 
  for (i = 0; i < a_motif->length; i++) {
    these_freqs = get_matrix_row(i, a_motif->scores);

    for (j = 0; j < alph_size; j++) {
      set_array_item(j, parse_prob(infile, buffer, "letter score"), these_freqs);
    }
  }
}

/***********************************************************************
 * Read a motif frequency matrix from MEME XML into a motif.
 ***********************************************************************/
static void read_pspm_from_xml(
  xmlXPathContextPtr xpath_ctxt,  // XML document stream IN
  ARRAY_T* background, // Pointer to array of bg freqs. IN
  double pseudocount,  // Pseudoount to be applied to motif freq. IN
  char* motif_id,      // XML id attribute for motif IN
  MOTIF_T* motif       // Motif structure OUT
) {

  int      motif_col;
  int      alpha_col;
  ARRAY_T* these_freqs;
  ATYPE    value;
  xmlXPathObjectPtr xpath_obj;
  char xpath_expression[MAX_XPATH_EXPRESSION];

  motif->freqs = allocate_matrix(motif->length, get_alph_size(ALL_SIZE));

  // Read the letter frequency matrix,
  // one array per position in the motif
  for (motif_col = 0; motif_col < motif->length; motif_col++) {

    these_freqs = get_matrix_row(motif_col, motif->freqs);

    // Build the XPATH expression to get one row of the freq matrix.
    snprintf(
      xpath_expression,
      MAX_XPATH_EXPRESSION,
      "/MEME/motifs/motif[@id='%s']/probabilities/alphabet_matrix/"
      "alphabet_array[position()=%d]/value",
      motif_id,
      motif_col + 1
    );
    xpath_obj = xpath_query(xpath_ctxt, xpath_expression);
    // The number of columns in the freq matrix should should
    // match the size of the alphabet.
    assert(motif->alph_size == xpath_obj->nodesetval->nodeNr);

    // Fill in the values of the row from the XPATH result.
    for (alpha_col = 0; alpha_col < motif->alph_size; alpha_col++) {
      double bg_letter_freq = get_array_item(alpha_col, background);
      xmlNodePtr currValueNode = xpath_obj->nodesetval->nodeTab[alpha_col];
      value = xmlXPathCastNodeToNumber(currValueNode);
      value = ((motif->num_sites * value) + pseudocount * bg_letter_freq)
        / (motif->num_sites  + pseudocount);
      set_array_item(alpha_col, value, these_freqs);
    }

    // Normalize the first alph_size positions. (MEME prints six
    // digits of precision).
    normalize_subarray(0, motif->alph_size, 0.00001, these_freqs);

    /* Compute values for ambiguous characters. */
    fill_in_ambiguous_chars(FALSE, these_freqs);

    xmlXPathFreeObject(xpath_obj);

  }

}

/***********************************************************************
 * Read a motif score matrix from MEME XML into a motif.
 ***********************************************************************/
static void read_pssm_from_xml(
  xmlXPathContextPtr xpath_ctxt,  // XML document stream IN
  char* motif_id,      // XML id attribute for motif IN
  MOTIF_T* motif       // Motif structure OUT
) {

  int      motif_col;
  int      alpha_col;
  ARRAY_T* these_freqs;
  ATYPE    value;
  xmlXPathObjectPtr xpath_obj;
  char xpath_expression[MAX_XPATH_EXPRESSION];

  motif->scores = allocate_matrix(motif->length, motif->alph_size);

  // Read the letter frequency matrix,
  // one array per position in the motif
  for (motif_col = 0; motif_col < motif->length; motif_col++) {

    these_freqs = get_matrix_row(motif_col, motif->scores);

    // Build the XPATH expression to get one row of the freq matrix.
    snprintf(
      xpath_expression,
      MAX_XPATH_EXPRESSION,
      "/MEME/motifs/motif[@id='%s']/scores/alphabet_matrix/"
      "alphabet_array[position()=%d]/value",
      motif_id,
      motif_col + 1
    );
    xpath_obj = xpath_query(xpath_ctxt, xpath_expression);
    // The number of columns in the freq matrix should should
    // match the size of the alphabet.
    assert(motif->alph_size == xpath_obj->nodesetval->nodeNr);

    // Fill in the values of the row from the XPATH result.
    for (alpha_col = 0; alpha_col < motif->alph_size; alpha_col++) {
      xmlNodePtr currValueNode = xpath_obj->nodesetval->nodeTab[alpha_col];
      value = xmlXPathCastNodeToNumber(currValueNode);
      set_array_item(alpha_col, value, these_freqs);
    }
    xmlXPathFreeObject(xpath_obj);
  }

}


/***********************************************************************
 * Read the next motif from a file into a motif data structure.
 *
 * RETURN: Boolean - Was a motif found?
 ***********************************************************************/
static MOTIF_T* read_motif(
  FILE    *infile,      /* An open file in MEME output format. IN */
  BUF_T   *buffer,      /* An open buffer of the infile IN */
  BMSTR_T **landmarks,     /* Search strings, motif start marker and motifs end marker IN */
  int      is_html,     /* Is the file html? IN */
  ARRAY_T *background,  /* Pointer to array of bg freqs. IN */
  double   pseudocount, /* Pseudocount IN */
  int      loading_mode /* Should PSSM or PSPM be loaded IN */
) {
  MOTIF_T *a_motif;
  char id1buf[MAX_MOTIF_ID_LENGTH+1], id2buf[MAX_MOTIF_ID_LENGTH+1];
  char *id1, *id2;
  int found, len;

  found = buf_fread_until(buffer, infile, NULL, 2, landmarks+3);

  if (found == 2 || found == 0) {
    return NULL;
  } else if (found != 1) {
    die("Error reading motif, error was given as %s\n", strerror(ferror(infile))); 
  }

  // Parse the motif
  a_motif = (MOTIF_T*)mm_malloc(sizeof(MOTIF_T));
  // skip over the motif start string
  buf_skip(buffer, bmstr_length(landmarks[3]));
  // skip spaces or tabs (not newlines though)
  len = buf_fread_consume(buffer, infile, is_blank, NULL, FALSE);
  if (len == -1) die("Error reading motif, error was given as %s\n", strerror(ferror(infile))); 
  if (!is_html) {
    // read id
    id1 = buf_fread_token(buffer, infile, buf_is_delim_space, NULL, FALSE, id1buf, MAX_MOTIF_ID_LENGTH+1, &len);
    if (id1 == NULL || len == 0) die("Error reading motif ID.\n");
    // skip spaces or tabs (not newlines though)
    len = buf_fread_consume(buffer, infile, is_blank, NULL, FALSE);
    if (len == -1) die("Error reading motif, error was given as %s\n", strerror(ferror(infile))); 
    // read id2
    id2 = buf_fread_token(buffer, infile, buf_is_delim_space, NULL, FALSE, id2buf, MAX_MOTIF_ID_LENGTH+1, &len);
    if (id2 == NULL) die("Error attempting to read motif ID 2.\n");
    //check that it's actually an id (id2 is optional and only used by TOMTOM currently)
    if (strcmp(id2, "width") == 0) {
      //looks like we're reading an ordinary meme file which doesn't have a second id (yet)
      id2[0] = '\0'; 
    }
  } else {
    //html only has one id, additionally we need to look for a < as that may terminate the id instead of a space
    id1 = buf_fread_token(buffer, infile, is_space_or_ab, NULL, FALSE, id1buf, MAX_MOTIF_ID_LENGTH+1, &len);
    if (id1 == NULL || len == 0) die("Error reading motif ID.\n");
    id2 = id2buf;
    id2[0] = '\0';
  }

  // Record the ID of this motif. 
  set_motif_id(id1, a_motif);
  set_motif_id2(id2, a_motif);
  if (verbosity >= HIGH_VERBOSE) {
    fprintf(stderr, " %s", get_motif_id(a_motif));
  }

  //set some defaults
  a_motif->length = INT_MIN;
  a_motif->alph_size = -1;
  a_motif->ambigs = -1;
  a_motif->evalue = -1;
  a_motif->num_sites = -1;
  a_motif->complexity = -1;
  a_motif->freqs = NULL;
  a_motif->scores = NULL;
  a_motif->url = NULL;
  a_motif->trim_left = 0;
  a_motif->trim_right = 0;
  
  //search for pssm or pspm
  found = find_bm_landmarks(infile, buffer, TRUE, FALSE, 3, landmarks);
  do {
    switch (found) {
    case 1:
      if (a_motif->scores != NULL) 
        die("PSSM landmark repeats in motif %s.\n", id1);
      buf_skip(buffer, bmstr_length(landmarks[0]));
      read_pssm(infile, buffer, loading_mode, a_motif);
      break;
    case 2:
      if (a_motif->freqs != NULL) 
        die("PSPM landmark repeats in motif %s.\n", id1);
      buf_skip(buffer, bmstr_length(landmarks[1]));
      read_pspm(infile, buffer, background, pseudocount, loading_mode, a_motif);
      break;
    case 3:
      buf_skip(buffer, bmstr_length(landmarks[2]));
      consume(infile, buffer, buf_is_delim_space, NULL, FALSE, 0, "spacer");
      a_motif->url = buf_fread_token(buffer, infile, buf_is_delim_space, NULL, FALSE, NULL, MAX_MOTIF_URL_LENGTH+1, NULL);
      if (a_motif->url == NULL) die("Error reading motif URL.\n");
      break;
    } 
    //search for pssm, pspm, or next motif, or eof
    found = find_bm_landmarks(infile, buffer, FALSE, FALSE, 5, landmarks);
  } while (found != 0 && found <= 3);//if next motif or eof then stop searching

  if (loading_mode & NEED_PSPM && a_motif->freqs == NULL) {
    die("Motif %s missing PSPM.\n", id1);
  }

  if (loading_mode & NEED_PSSM && a_motif->scores == NULL) {
    die("Motif %s missing PSSM.\n", id1);
  }

  if (a_motif->url == NULL) {
    set_motif_url("", a_motif);
  }

  // Store the alphabet size in the motif. 
  a_motif->alph_size = get_alph_size(ALPH_SIZE);
  a_motif->ambigs = get_alph_size(AMBIG_SIZE);

  if (a_motif->freqs != NULL) {
    // Compute and store the motif complexity.
    a_motif->complexity = compute_motif_complexity(a_motif);
    if (!(loading_mode & WANT_PSPM)) {
      //the caller doesn't want the pspm, destroy it
      free_matrix(a_motif->freqs);
      a_motif->freqs = NULL;
    }
  }

  return a_motif;
}

/***********************************************************************
 * Read just the motifs from a MEME file.
 ***********************************************************************/
static void read_motifs(
  FILE    *motif_file,       // MEME file. IN
  BUF_T   *buffer,           // buffer IN
  int      is_html,          // Is the file html? IN
  ARRAY_T *background,       // Pointer to array of bg freqs. IN
  double   pseudocount,      // Pseudocount IN
  int      loading_mode,     // Should PSSM or PSPM be loaded IN
  ARRAYLST_T *motifs            // The retrieved motifs.  OUT
) {
  int initial_motif_count;
  BMSTR_T *landmarks[5];
  MOTIF_T *motif;
  initial_motif_count = arraylst_size(motifs);
  if (verbosity >= HIGH_VERBOSE) {
    fprintf(stderr, "Reading motif");
  }
  //setup search strings
  landmarks[0] = bmstr_create("log-odds matrix: ");
  landmarks[1] = bmstr_create("letter-probability matrix: ");
  landmarks[2] = bmstr_create("URL ");
  landmarks[3] = bmstr_create("MOTIF ");
  landmarks[4] = bmstr_create(END_OF_MOTIF);

  while ((motif = read_motif(motif_file, buffer, landmarks, is_html, background, pseudocount, loading_mode)) != NULL) {
    arraylst_add(motif, motifs);
  }
  arraylst_fit(motifs);

  bmstr_destroy(landmarks[0]);
  bmstr_destroy(landmarks[1]);
  bmstr_destroy(landmarks[2]);
  bmstr_destroy(landmarks[3]);
  bmstr_destroy(landmarks[4]);

  if (verbosity >= HIGH_VERBOSE) {
    fprintf(stderr, ".\n");
  }
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Read %d motifs.\n", (arraylst_size(motifs) - initial_motif_count));
  }
}

/***********************************************************************
 * Read the motifs from an XML document.
 * The caller is responsible for allocating the array of motifs.
 ***********************************************************************/
static void read_motifs_from_xml(
  xmlXPathContextPtr xpath_ctxt,    // MEME XPath context.
  ARRAY_T*  background,             // Pointer to array of bg freqs. IN
  double    pseudocount,            // pseudocount to be applied to motif IN
  int       loading_mode,          // Load PSSM or PSPM IN
  ARRAYLST_T* motifs               // The retrieved motifs.  OUT
) {
  int initial_motif_count;
  xmlXPathObjectPtr xpath_obj = NULL;
  xmlChar* property = NULL;
  char* path = "/MEME/motifs/motif"; // This path will select all motifs.

  initial_motif_count = arraylst_size(motifs);

  // Use XPATH to get the set of motifs.
  xpath_obj = xpath_query(xpath_ctxt, path);

  int num_motifs = (xpath_obj->nodesetval ? xpath_obj->nodesetval->nodeNr : 0);
  arraylst_preallocate(arraylst_size(motifs) + num_motifs, motifs);

  if (verbosity >= HIGH_VERBOSE) {
    fprintf(stderr, "Reading motif ");
  }

  xmlNodePtr currMotifNode = NULL;
  int i_motif = 0;
  for (i_motif = 0; i_motif < num_motifs; i_motif++) {
    MOTIF_T *motif = (MOTIF_T*)mm_malloc(sizeof(MOTIF_T));
    memset(motif, 0, sizeof(MOTIF_T));
    arraylst_add(motif, motifs);

    currMotifNode = xpath_obj->nodesetval->nodeTab[i_motif];
    if (currMotifNode == NULL) {
      die("Error: missing motif %d\n", i_motif);
    }

    // Get the motif name attribute
    property = read_xml_node_property(currMotifNode, "name");
    set_motif_id((char *) property, motif);
    set_motif_id2("", motif);
    xmlFree(property);

    // Get the motif length attribute
    property = read_xml_node_property(currMotifNode, "width");
    motif->length = atoi((char *) property);
    xmlFree(property);

    // Get the motif evalue attribute
    property = read_xml_node_property(currMotifNode, "e_value");
    motif->evalue = atof((char *) property);
    xmlFree(property);

    // Get the motif sites attribute
    property = read_xml_node_property(currMotifNode, "sites");
    motif->num_sites = atof((char *) property);
    xmlFree(property);

    // Get the optional motif url
    property = xmlGetProp(currMotifNode, BAD_CAST "url");
    if (property) {
      set_motif_url((char*)property, motif);
      xmlFree(property);
    } else {
      set_motif_url("", motif);
    }

    // Store the alphabet size in the motif.
    motif->alph_size = get_alph_size(ALPH_SIZE);
    motif->ambigs = get_alph_size(AMBIG_SIZE);

    // Get the motif id attribute
    xmlChar* motif_xml_id = read_xml_node_property(currMotifNode, "id");


    //set some defaults
    motif->freqs = NULL;
    motif->scores = NULL;
    motif->trim_left = 0;
    motif->trim_right = 0;

    // As no calculations depend on it, only get the PSSM
    if (loading_mode & WANT_PSSM) {
      // Get the score matrix for the motif.
      read_pssm_from_xml(
        xpath_ctxt,
        (char *) motif_xml_id,
        motif
      );
    }
    // Get the freq. matrix for the motif.
    read_pspm_from_xml(
      xpath_ctxt,
      background,
      pseudocount,
      (char *) motif_xml_id,
      motif
    );
    if (loading_mode & NEED_PSSM && motif->scores == NULL) {
      die("Motif %s missing PSSM.\n", motif->id);
    }
    if (loading_mode & NEED_PSPM && motif->freqs == NULL) {
      die("Motif %s missing PSPM.\n", motif->id);
    }
    xmlFree(motif_xml_id);

    // Compute and store the motif complexity.
    motif->complexity = compute_motif_complexity(motif);

    //if the caller doesn't want the pspm, don't store it longer than required
    if (!(loading_mode & WANT_PSPM)) {
      free_matrix(motif->freqs);
      motif->freqs = NULL;
    }

    if (verbosity >= HIGH_VERBOSE) {
      fprintf(stderr, " %s", get_motif_id(motif));
    }

  }

  if (verbosity >= HIGH_VERBOSE) {
    fprintf(stderr, ".\n");
  }

  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Read %d motifs.\n", (arraylst_size(motifs) - initial_motif_count));
  }

  /* Cleanup */
  xmlXPathFreeObject(xpath_obj);

}


/***********************************************************************
 * Create two copies of each motif.  The new IDs are preceded by "+"
 * and "-", and the "-" version is the reverse complement of the "+"
 * version.

void add_reverse_complements
  (int* num_motifs,
   MOTIF_T* motifs)
{
  int i_motif;                  // Index of the current motif.
  char motif_id[MAX_MOTIF_ID_LENGTH + 1]; // Name of the current motif;

  // Copy motifs.
  for (i_motif = 0; i_motif < *num_motifs; i_motif++) {
    copy_motif(&(motifs[i_motif]), &(motifs[*num_motifs + i_motif]));
    assert(motifs[i_motif].length != 0);
    assert(motifs[*num_motifs + i_motif].length != 0);
  }

  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Reverse complementing:");
  }

  // Add plusses to the first n motifs.
  motif_id[0] = '+';
  for (i_motif = 0; i_motif < *num_motifs; i_motif++) {
    strcpy(&(motif_id[1]), get_motif_id(&(motifs[i_motif])));
    set_motif_id(motif_id, &(motifs[i_motif]));

    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(stderr, " %s", motif_id);
    }
  }

  // Add minuses to the second n motifs.
  motif_id[0] = '-';
  for (i_motif = 0; i_motif < *num_motifs; i_motif++) {
    strcpy(&(motif_id[1]), get_motif_id(&(motifs[i_motif + *num_motifs])));
    set_motif_id(motif_id, &(motifs[i_motif + *num_motifs]));

    if (verbosity >= NORMAL_VERBOSE) {
      fprintf(stderr, " %s", motif_id);
    }
  }

  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "\n");
  }

  // Make the reverse complements.
  for (i_motif = 0; i_motif < *num_motifs; i_motif++) {
    reverse_complement_motif(&(motifs[i_motif + *num_motifs]));
  }

  // Double the motif counter.
  *num_motifs *= 2;
}
 ***********************************************************************/

/***********************************************************************
 * Create two copies of each motif.  The new IDs are preceded by "+"
 * and "-", and the "-" version is the reverse complement of the "+"
 * version.
 *
 * John Hawkins 2008 - I am changing this function so that the inverse
 * motifs are placed directly after the original in the array. This
 * helps with the implmentation of the BLS scan mode and has no apparent
 * effect on the other scan modes
 ***********************************************************************/
void add_reverse_complements
  (int* num_motifs,
   MOTIF_T* motifs)
{
  int i_motif;                  // Index of the current motif.
  char motif_id[MAX_MOTIF_ID_LENGTH + 1]; // Name of the current motif;

  // Copy motifs and change the IDs

  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Reverse complementing:");
  }

  //fprintf(stderr, "num motifs: %d\n", *num_motifs);

  // JH - Changed here so that copies end up directly after each original
  for (i_motif = (*num_motifs)-1; i_motif > -1; i_motif--) {

	// First copy the original to its new place
	if(i_motif > 0) {
	   copy_motif( &(motifs[i_motif]), &(motifs[ 2 * i_motif ]));
  }
	// Now make the second copy
	copy_motif( &(motifs[i_motif]), &(motifs[ 2 * i_motif + 1 ]));
  // Free freq. matrix from original.
	if(i_motif > 0) {
    free_matrix(motifs[i_motif].freqs);
  }

	// Now add plusses to the original motifs.
        motif_id[0] = '+';
	strcpy(&(motif_id[1]), get_motif_id(&(motifs[2 * i_motif])));
	set_motif_id(motif_id, &(motifs[2 * i_motif]));
	if (verbosity >= NORMAL_VERBOSE) {
           fprintf(stderr, " %s", motif_id);
        }

	// and add minuses to the copy motifs.
	motif_id[0] = '-';
	strcpy(&(motif_id[1]), get_motif_id(&(motifs[2 * i_motif + 1])));
    	set_motif_id(motif_id, &(motifs[2 * i_motif + 1]));
	if (verbosity >= NORMAL_VERBOSE) {
      		fprintf(stderr, " %s", motif_id);
    	}

	// Make the copy a reverse complement.
    	reverse_complement_motif(&(motifs[2 * i_motif + 1]));

  }
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "\n");
  }
  // Double the motif counter.
  *num_motifs *= 2;

}
 /***********************************************************************/

/***********************************************************************
 * Create two copies of each motif.  The new IDs are preceded by "+"
 * and "-", and the "-" version is the reverse complement of the "+"
 * version.
 *
 * The reverse complement is always listed directly after the original.
 ***********************************************************************/
void add_reverse_complements2
  (ARRAYLST_T* motifs)
{
  int i, count;
  char motif_id[MAX_MOTIF_ID_LENGTH + 1]; // Name of the current motif;

  count = arraylst_size(motifs);
  //double the array size
  arraylst_preallocate(count*2, motifs);
  for (i = 0; i < count; ++i) {
    MOTIF_T *new_motif = (MOTIF_T*)mm_malloc(sizeof(MOTIF_T));
    arraylst_add(new_motif, motifs);
  }
  //move and reverse complement the original motifs
  for (i = count-1; i >= 0; --i) {
    //move the motif to its new position
    arraylst_swap(i, 2*i, motifs);
    //get the motif and the one that will become the reverse complement
    MOTIF_T *orig = arraylst_get(2*i, motifs);
    MOTIF_T *rc = arraylst_get(2*i + 1, motifs);
    //copy and reverse complement the motif to the position after
    copy_motif(orig, rc);
    reverse_complement_motif(rc);//FIXME scores are not reverse complemented, only freqs
    //prepend a plus to the original motif name
    motif_id[0] = '+';
    strncpy(motif_id+1, get_motif_id(orig), MAX_MOTIF_ID_LENGTH);
    motif_id[MAX_MOTIF_ID_LENGTH] = '\0';//truncate if we must
    //set the motif ids to have the + or - in front
    set_motif_id(motif_id, orig);
    motif_id[0] = '-';
    set_motif_id(motif_id, rc);
  }
}
 /***********************************************************************/


/***********************************************************************
 * Sets the trim boundaries for all the motifs by scanning inward until 
 * a position with an information content larger or equal to the 
 * threshold is encountered.
 ***********************************************************************/
void trim_motifs_by_bit_threshold(
  ARRAYLST_T* motifs, 
  double threshold_bits
) {
  int i, len;
  MOTIF_T *motif;
  len = arraylst_size(motifs);
  for (i = 0; i < len; ++i) {
    motif = (MOTIF_T*)arraylst_get(i, motifs);
    trim_motif_by_bit_threshold(motif, threshold_bits);
  }
}
 /***********************************************************************/


/***********************************************************************
 * Read the motif occurrences info from a old MEME html file.
 *
 * Each line contains the following items
 *  - sequence id,
 *  - sequence p-value,
 *  - number n of motif occurrences, and
 *  - length of sequence.
 *
 * This is followed by n triples containing
 *  - motif id,
 *  - occurrence position, and
 *  - occurrence p-value.
 *
 * Return value: Does the given file contain a motif occurrence section?
 ***********************************************************************/
static BOOLEAN_T read_motif_occurrences_from_old_html(
  FILE          *motif_file,        /* MEME file. IN */
  BUF_T         *buffer,            /* buffer. IN */
  STRING_LIST_T **motif_occurrences  /* List of motif occurence strings. OUT */
)
{
  int length;
  char *line;

  assert(motif_file != NULL);
  assert(buffer != NULL);
  assert(motif_occurrences != NULL);

  if (!find_landmarks(motif_file, buffer, FALSE, TRUE, 3, START_OCCURRENCES_1, START_OCCURRENCES_2, START_OCCURRENCES_3)) {
    *motif_occurrences = NULL;
    return FALSE;
  }

  if (buf_fread_next_line(buffer, motif_file) < 0)
    die("File error reading motif occurrences, error was %s\n", strerror(ferror(motif_file)));
 
  //now read the motif occurrences
  *motif_occurrences = new_string_list();
  line = NULL;
  while (buf_remaining(buffer)) {
    line = buf_fread_token(buffer, motif_file, is_newline, NULL, FALSE, NULL, MAX_LINE, &length);
    if (line == NULL)
      die("File error reading motif occurrences, error was %s\n", strerror(ferror(motif_file)));
    if (buf_fread_consume(buffer, motif_file, is_newline, NULL, FALSE) < 0)
      die("File error reading motif occurrences, error was %s\n", strerror(ferror(motif_file)));
    // for compatibility I need to add on a new line char
    line = mm_realloc(line, sizeof(char)*(length+2));
    line[length] = '\n';
    line[length+1] = '\0';
    // test for end
    if (strcmp(END_OCCURRENCES_1, line) == 0) break;
    if (strcmp(END_OCCURRENCES_2, line) == 0) break;
    // not end so add line to list 
    // string list is really not a good datastructure choice in this case as it
    // should only be used with a lot of strings of exactly the same length.
    // however it would be a lot of bother to change it...
    add_string(line, *motif_occurrences);
    free(line);
    line = NULL;
  }
  if (line != NULL) free(line);

  return TRUE;

  /*
  char      line[MAX_LINE]; // Buffer for reading.

  assert(motif_occurences != NULL);

  // Look for the beginning of the motif occurrence section.
  while (TRUE) {
    if (fgets(line, MAX_LINE, motif_file) == NULL) {
      *motif_occurences = NULL;
      return(FALSE);
    }
    // The tag used to mark the begining of the
    // motif occurence section has changed slightly
    // over time. Check for both versions.
    if (strcmp(line, START_OCCURRENCES_1) == 0) {
      break;
    }
    if (strcmp(line, START_OCCURRENCES_2) == 0) {
      break;
    }
    if (strcmp(line, START_OCCURRENCES_3) == 0) {
      break;
    }
  }

  *motif_occurences = new_string_list();

  // Read each line, corresponding to each sequence.
  while (fgets(line, MAX_LINE, motif_file)) {
    // Check to see if we've reached the end of the section.
    // The tag used to mark the end of the
    // motif occurence section has changed slightly
    // over time. Check for both versions.
    if (strcmp(line, END_OCCURRENCES_1) == 0) {
      break;
    }
    if (strcmp(line, END_OCCURRENCES_2) == 0) {
      break;
    }
    add_string(line, *motif_occurences);
  }

  return(TRUE);
  */
}

/***********************************************************************
 * Read the motif occurrences info from a meme txt file
 *
 * Each line contains the following items
 *  - sequence id,
 *  - sequence p-value,
 *  - number n of motif occurrences, and
 *  - length of sequence.
 *
 * This is followed by n triples containing
 *  - motif id,
 *  - occurrence position, and
 *  - occurrence p-value.
 *
 * Return value: Does the given file contain a motif occurrence section?
 ***********************************************************************/
static BOOLEAN_T read_motif_occurrences_from_txt(
  FILE          *motif_file,        /* MEME file. IN */
  BUF_T         *buffer,            /* buffer. IN */
  ARRAYLST_T    *motifs,            /* previously parsed motifs (for the lengths) */
  STRING_LIST_T **motif_occurrences  /* List of motif occurence strings. OUT */
)
{
  char *seqid, *cpvalue, *diagram, *curr_id, *curr_pvalue, **occur_id, **occur_pvalue, *line;
  int i, unparsed, grow, seq_len, curr_pos, occurrences, occur_alloc, *occur_pos, line_len, line_pos;
  //find the heading 2 lines above the block diagram
  if (!find_landmarks(motif_file, buffer, FALSE, TRUE, 1, "COMBINED P-VALUE")) return FALSE;
  //skip end of current line
  if(buf_fread_next_line(buffer, motif_file) < 0) 
    die("File error while reading motif occurrence section, error was %s\n", strerror(ferror(motif_file)));
  //skip line containing underlines eg: ------
  if(buf_fread_next_line(buffer, motif_file) < 0) 
    die("File error while reading motif occurrence section, error was %s\n", strerror(ferror(motif_file)));
  //now we parse and construct lines
  occur_alloc = 2;
  occur_pos = (int*)mm_malloc(sizeof(int)*occur_alloc);
  occur_id = (char**)mm_malloc(sizeof(char*)*occur_alloc);
  occur_pvalue = (char**)mm_malloc(sizeof(char*)*occur_alloc);
  //allocate initial line
  line_len = 10;
  line = (char*)mm_malloc(sizeof(char)*line_len);
  *motif_occurrences = new_string_list();
  while (TRUE) {
    seqid = parse_word(motif_file, buffer, NULL, MAX_LINE+1, 0, FALSE, "sequence id");
    if (strcmp(seqid, END_OF_TEXT_OCCURRENCES) == 0) {
      // The line
      // -----------------------------------------------------------------------------
      // marks the end of the summary of motif occurrences.
      free(seqid);
      break;
    }
    cpvalue = parse_word(motif_file, buffer, NULL, MAX_LINE+1, 1, TRUE, "combined p-value");
    diagram = parse_word(motif_file, buffer, NULL, MAX_LINE+1, 1, TRUE, "motif diagram");
    //now process the diagram
    seq_len = 0;
    occurrences = 0;
    curr_pos = -1;
    curr_id = NULL;
    curr_pvalue = NULL;

    for (i = 0, unparsed = 0; diagram[i] != '\0'; ++i) {
      if (diagram[i] == '_') {
        //found an underline, it could be before or after the number
        if (unparsed < i) { //this is the underline after the number
          int num;
          diagram[i] = '\0';
          num = atoi(diagram+unparsed);
          if (num <= 0) die("Expected a positive number.");
          seq_len += num;
        }
        unparsed = i + 1;
      } else if (diagram[i] == '[') {//start of a motif occurrence
        unparsed = i + 1;
      } else if (diagram[i] == '(' && unparsed < i) {//work out which motif it is
        int m_index;
        //parse the motif id 
        diagram[i] = '\0';
        curr_id = diagram+unparsed;
        if (diagram[unparsed] == '-' || diagram[unparsed] == '+') {
          ++unparsed;
        } 
        for (m_index = 0; m_index < arraylst_size(motifs); ++m_index) {
          MOTIF_T *motif = (MOTIF_T*)arraylst_get(m_index, motifs);
          if (strcmp(motif->id, diagram+unparsed) == 0) {
            curr_pos = seq_len;
            seq_len += motif->length;
            break;
          }
        }
        if (m_index >= arraylst_size(motifs)) {
          die("Could not identify motif id in block diagram\n");
        }
        unparsed = i + 1;
      } else if (diagram[i] == ')' && unparsed < i) {
        diagram[i] = '\0';
        curr_pvalue = diagram+unparsed;
        unparsed = i + 1;
      } else if (diagram[i] == ']') {
        if (occurrences == occur_alloc) {
          occur_alloc += 3;
          occur_pos = (int*)mm_realloc(occur_pos, sizeof(int)*occur_alloc);
          occur_id = (char**)mm_realloc(occur_id, sizeof(char*)*occur_alloc);
          occur_pvalue = (char**)mm_realloc(occur_pvalue, sizeof(char*)*occur_alloc);
        }
        if (curr_pos == -1 || curr_id == NULL || curr_pvalue == NULL) {
          die("Motif not parsed correctly from combined block diagram\n");
        }
        occur_pos[occurrences] = curr_pos;
        occur_id[occurrences] = curr_id;
        occur_pvalue[occurrences] = curr_pvalue;
        curr_pos = -1;
        curr_id = NULL;
        curr_pvalue = NULL;
        ++occurrences;
        unparsed = i + 1;
      }
    }
    if (unparsed < i) {
      int num;
      diagram[i] = '\0';
      num = atoi(diagram+unparsed);
      if (num <= 0) die("Expected a positive number.");
      seq_len += num;
    }
    // create the line
    line_pos = 0;
    do {
      char *buf = line+line_pos;
      int   len  = line_len - line_pos;
      grow = snprintf(buf, sizeof(char)*len, "%s %s %d %d", seqid, cpvalue, occurrences, seq_len);
      if (grow >= len) {
        // Resize if we weren't able to print the full string
        // NB use >= test because the return value of snprint doesn't 
        // include the trailing '\0'
        line_len = line_pos + grow + 1;
        line = (char*)mm_realloc(line, line_len); 
      } else {
        line_pos += grow;
        grow = 0;
      }
    } while (grow);
    for (i = 0; i < occurrences; ++i) {
      do {
        char *buf = line+line_pos;
        int   len  = line_len - line_pos;
        grow = snprintf(buf, sizeof(char)*len, " %s %d %s", occur_id[i], occur_pos[i], occur_pvalue[i]);
        if (grow >= len) {
          // Resize if we weren't able to print the full string
          // NB use >= test because the return value of snprint doesn't 
          // include the trailing '\0'
          line_len = line_pos + grow + 1;
          line = (char*)mm_realloc(line, line_len); 
        } else {
          line_pos += grow;
          grow = 0;
        }
      } while (grow);
    }
    add_string(line, *motif_occurrences); 
    if (buf_fread_consume(buffer, motif_file, is_newline, NULL, FALSE) < 0)
      die("File error reading motif occurrences, error was %s\n", strerror(ferror(motif_file)));
    free(seqid);
    free(cpvalue);
    free(diagram);
  }
  free(line);
  free(occur_pos);
  free(occur_id);
  free(occur_pvalue);
  return TRUE;
}
/***********************************************************************
 * Look up the sequence name corresponding to a sequence id in
 * the XML file.
 *
 * Returns the squence name. Caller is responsible for freeing the name.
 ***********************************************************************/
static xmlChar* get_sequence_name_from_id(xmlXPathContextPtr ctxt, xmlChar* id) {

  xmlXPathObjectPtr sequence_obj;
  char xpath_expression[MAX_XPATH_EXPRESSION];
  xmlNodePtr node = NULL;

  int char_written = snprintf(
    xpath_expression,
    MAX_XPATH_EXPRESSION,
    "/MEME/training_set/sequence[@id='%s']",
    id
  );

  // Check that we didn't have to truncate the string.
  if (char_written >= MAX_XPATH_EXPRESSION) {
    die("Unable to get sequence name from sequence id. XPath expression exceeded maxiumum allowed size.");
  }
  sequence_obj = xpath_query(ctxt, xpath_expression);
  if (sequence_obj == NULL) {
    die("Unable to query XML file for name of sequence with id %s.", id);
  }
  node = sequence_obj->nodesetval->nodeTab[0];
  if (node == NULL) {
    die("Unable to find name for sequence with id %s.", id);
  }
  xmlChar* sequence_name = xmlGetProp(node, BAD_CAST "name");
  xmlXPathFreeObject(sequence_obj);

  if (sequence_name == NULL) {
    die("Unable to find name for sequence with id %s.", id);
  }

  return sequence_name;
}

/***********************************************************************
 * Look up the motif name corresponding to a motif id in
 * the XML file.
 *
 * Returns the squence name. Caller is responsible for freeing the name.
 ***********************************************************************/
static xmlChar* get_motif_name_from_id(xmlXPathContextPtr ctxt, xmlChar* id) {

  xmlXPathObjectPtr motif_obj;
  char xpath_expression[MAX_XPATH_EXPRESSION];
  xmlNodePtr node = NULL;

  int char_written = snprintf(
    xpath_expression,
    MAX_XPATH_EXPRESSION,
    "/MEME/motifs/motif[@id='%s']",
    id
  );

  // Check that we didn't have to truncate the string.
  if (char_written >= MAX_XPATH_EXPRESSION) {
    die("Unable to get motif name from motif id. XPath expression exceeded maxiumum allowed size.");
  }
  motif_obj = xpath_query(ctxt, xpath_expression);
  if (motif_obj == NULL) {
    die("Unable to query XML file for name of motif with id %s.", id);
  }
  node = motif_obj->nodesetval->nodeTab[0];
  if (node == NULL) {
    die("Unable to find name for motif with id %s.", id);
  }
  xmlChar* motif_name = xmlGetProp(node, BAD_CAST "name");
  xmlXPathFreeObject(motif_obj);

  if (motif_name == NULL) {
    die("Unable to find name for motif with id %s.", id);
  }

  return motif_name;
}
/***********************************************************************
 * Look up the sequence length corresponding to a sequence id
 * in the XML file.
 *
 * Returns the squence length as a string. Caller is responsible
 * for freeing the string.
 ***********************************************************************/
static xmlChar* get_sequence_length_from_id(xmlXPathContextPtr ctxt, xmlChar* id) {

  xmlXPathObjectPtr sequence_obj;
  char xpath_expression[MAX_XPATH_EXPRESSION];
  xmlNodePtr node = NULL;

  int char_written = snprintf(
    xpath_expression,
    MAX_XPATH_EXPRESSION,
    "/MEME/training_set/sequence[@id='%s']",
    id
  );

  // Check that we didn't have to truncate the string.
  if (char_written >= MAX_XPATH_EXPRESSION) {
    die("Unable to get sequence name from sequence id. XPath expression exceeded maxiumum allowed size.");
  }
  sequence_obj = xpath_query(ctxt, xpath_expression);
  if (sequence_obj == NULL) {
    die("Unable to query XML file for length of sequence with id %s.", id);
  }
  node = sequence_obj->nodesetval->nodeTab[0];
  if (node == NULL) {
    die("Unable to find length for sequence with id %s.", id);
  }
  xmlChar* sequence_length = xmlGetProp(node, BAD_CAST "length");
  xmlXPathFreeObject(sequence_obj);

  if (sequence_length == NULL) {
    die("Unable to find length for sequence with id %s.", id);
  }

  return sequence_length;
}

/***********************************************************************
 * Read the motif occurrences info from MEME XML.
 *
 * Each line contains the following items
 *  - sequence id,
 *  - sequence p-value,
 *  - number n of motif occurrences, and
 *  - length of sequence.
 *
 * This is followed by n triples containing
 *  - motif id,
 *  - occurrence position, and
 *  - occurrence p-value.
 *
 * Return value: Does the given file contain a motif occurrence section?
 ***********************************************************************/
static BOOLEAN_T read_motif_occurrences_from_xml(
  xmlXPathContextPtr xpath_ctxt,    // MEME XPath context.
  STRING_LIST_T** motif_occurrences  // List of motif occurrence strings. OUT
) {
  assert(motif_occurrences != NULL);

  *motif_occurrences = new_string_list();

  // Get all the scanned sites.
  xmlXPathObjectPtr scanned_sites_obj = NULL;
  scanned_sites_obj = xpath_query(xpath_ctxt, "/MEME/scanned_sites_summary/scanned_sites");

  // Get scanned sites for each sequence.
  int num_sequences =
    (scanned_sites_obj->nodesetval ? scanned_sites_obj->nodesetval->nodeNr : 0);
  xmlNodePtr currOccurenceNode = NULL;
  int i = 0;
  for (i = 0; i < num_sequences; i++) {

    currOccurenceNode = scanned_sites_obj->nodesetval->nodeTab[i];
    xmlChar* sequence_id = xmlGetProp(currOccurenceNode, BAD_CAST "sequence_id");
    xmlChar* pvalue = xmlGetProp(currOccurenceNode, BAD_CAST "pvalue");
    xmlChar* num_sites = xmlGetProp(currOccurenceNode, BAD_CAST "num_sites");
    xmlChar* sequence_name = get_sequence_name_from_id(xpath_ctxt, sequence_id);
    xmlChar* sequence_length = get_sequence_length_from_id(xpath_ctxt, sequence_id);
    char occurrence_string[255];
    strcpy(occurrence_string, (char *) sequence_name);
    strcat(occurrence_string, " ");
    strcat(occurrence_string, (char *) pvalue);
    strcat(occurrence_string, " ");
    strcat(occurrence_string, (char *) num_sites);
    strcat(occurrence_string, " ");
    strcat(occurrence_string, (char *) sequence_length);
    strcat(occurrence_string, " ");
    xmlFree(sequence_id);
    xmlFree(pvalue);
    xmlFree(num_sites);
    xmlFree(sequence_name);
    xmlFree(sequence_length);

    xmlNodePtr currChild = currOccurenceNode->children;
    while (currChild != NULL) {

      if (currChild->type == XML_ELEMENT_NODE) {
        xmlChar* motif_id = xmlGetProp(currChild, BAD_CAST "motif_id");
        xmlChar* motif_name = get_motif_name_from_id(xpath_ctxt, motif_id);
        xmlChar* strand = xmlGetProp(currChild, BAD_CAST "strand");
        xmlChar* position = xmlGetProp(currChild, BAD_CAST "position");
        xmlChar* site_pvalue = xmlGetProp(currChild, BAD_CAST "pvalue");
        if (strncmp("plus", (char *) strand, 4) == 0) {
          strcat(occurrence_string, "+");
        }
        else if (strncmp("minus", (char *) strand, 5) == 0) {
          strcat(occurrence_string, "-");
        }
        strcat(occurrence_string, (char *) motif_name);
        strcat(occurrence_string, " ");
        strcat(occurrence_string, (char *) position);
        strcat(occurrence_string, " ");
        strcat(occurrence_string, (char *) site_pvalue);
        strcat(occurrence_string, " ");
        xmlFree(motif_id);
        xmlFree(motif_name);
        xmlFree(strand);
        xmlFree(position);
        xmlFree(site_pvalue);
      }
      currChild = currChild->next;

    }

    add_string(occurrence_string, *motif_occurrences);

  }

  xmlXPathFreeObject(scanned_sites_obj);

  return(TRUE);

}

/*************************************************************************
 * Setup motif-to-motif occurrence and spacer length frequency
 * transition matrices in a naive fashion, without taking into account
 * any motif occurrence information.
 *************************************************************************/
#define SPACER_LENGTH 9.0 // Expected length of a spacer.
void compute_naive_transitions_and_spacers
  (const int  nmotifs,     // The number of motifs.
   MATRIX_T** transp_freq, // Motif-to-motif occurrence matrix.
   MATRIX_T** spacer_ave)  // Average spacer length matrix.
{
  int   i;
  int   j;
  PROB_T prob;

  // Allocate the matrices.
  *transp_freq = allocate_matrix(nmotifs + 2, nmotifs + 2);
  *spacer_ave = allocate_matrix(nmotifs + 2, nmotifs + 2);

  // Set up transition frequencies and spacer lengths from start state.
  prob = 1.0 / (PROB_T)nmotifs;
  for (j = 1; j <= nmotifs; j++) {
    set_matrix_cell(0, j, prob, *transp_freq);
    set_matrix_cell(0, j, SPACER_LENGTH, *spacer_ave);
  }

  /* Set up transition frequencies and spacer lengths between motifs
     and to the end state. */
  prob = 1.0 / ((PROB_T)(nmotifs + 1));
  for (i = 1; i <= nmotifs; i++) {
    for (j = 1; j <= nmotifs+1; j++) {
      set_matrix_cell(i, j, prob, *transp_freq);
      set_matrix_cell(i, j, SPACER_LENGTH, *spacer_ave);
    }
  }
}

/***********************************************************************
 * Read a MEME file in XML format using the libxml2 XML parsing tools.
 ***********************************************************************/
BOOLEAN_T read_meme_xml_file(
   char*      meme_filename,        // Name of MEME file  IN
   char*      bg_filename,            // Name of background freq. file IN
   double     pseudocount,          // Pseudocount to be applied motif freqs. IN
   int        loading_mode,           // Load PSSM or PSPM IN
   ARRAYLST_T* motifs,               // The retrieved motifs
   STRING_LIST_T** motif_occurrences, // Strings desc. motif occurrences  OUT
   BOOLEAN_T* has_reverse_strand,     // Does this file have both strands? OUT
   ARRAY_T**  background              // Background emission distribution  OUT
)
{
  xmlParserCtxtPtr ctxt = NULL;         // The parser context
  xmlDocPtr meme_doc = NULL;            // The resulting document tree
  xmlXPathContextPtr xpath_ctxt = NULL; // XPath context.

  ctxt = xmlNewParserCtxt();
  if (ctxt == NULL) {
    die("Failed to create XML parser.\n");
  }

  // Parse and validate the file.
  meme_doc = xmlCtxtReadFile(
    ctxt,
    meme_filename,
    NULL,  // Encoding
    XML_PARSE_RECOVER | XML_PARSE_DTDVALID | XML_PARSE_NOERROR | XML_PARSE_NOWARNING
  );

  // Did it parse?
  if (meme_doc == NULL || ctxt->valid == 0) {
   xmlFreeDoc(meme_doc);
   xmlFreeParserCtxt(ctxt);
   xmlCleanupParser();
   return FALSE;
  }
  if (verbosity > NORMAL_VERBOSE) {
    fprintf(stderr, "File %s is a valid MEME XML file.\n", meme_filename);
  }

  // Set up XPath context from parsed XML
  xpath_ctxt = xmlXPathNewContext(meme_doc);

  // Read the alphabet.
  read_alphabet_from_xml(xpath_ctxt);

  // Read the strandedness.
  if (which_alphabet() == DNA_ALPH) {
    *has_reverse_strand = read_strand_from_xml(xpath_ctxt);
  } else {
    *has_reverse_strand = FALSE;
  }

  // Establish the background frequencies.
  if (bg_filename == NULL) {
    // Default is to use pre-calculated bg freq. from NR database.
    if (verbosity >= HIGH_VERBOSE) {
      fprintf(
        stderr,
        "Using background frequencies from NR sequence database.\n"
      );
    }
    *background = allocate_array(get_alph_size(ALL_SIZE));
    get_nrdb_frequencies(*background);
  }
  else {
    if (strcmp(bg_filename, "motif-file") == 0) {
      // If bg_filename matches "motif-file" read bg freq. from motif file.
      if (verbosity >= NORMAL_VERBOSE) {
        fprintf(
          stderr,
          "Using background frequencies from file %s.\n",
          meme_filename
        );
      }
      *background = read_bg_freqs_from_xml(xpath_ctxt);
    }
    else {
      // Otherwise read bg freqs. from external bg file.
      if (verbosity >= NORMAL_VERBOSE) {
        fprintf(
          stderr,
          "Using background frequencies from file %s.\n",
          bg_filename
        );
      }
      *background = read_background_file(bg_filename);;
    }
  }

  // Read the motifs
  read_motifs_from_xml(
    xpath_ctxt,
    *background,
    pseudocount,
    loading_mode,
    motifs
  );

  // Read the motif occurences
  if (motif_occurrences) {
    read_motif_occurrences_from_xml(xpath_ctxt, motif_occurrences);
  }

  /* free up the resulting document */
  xmlXPathFreeContext(xpath_ctxt);
  xmlFreeDoc(meme_doc);
  xmlFreeParserCtxt(ctxt);
  xmlCleanupParser();

  return TRUE;

}

typedef struct html_reader_data HTMLRD_T;

/*
 * structure grouping data needed by
 * the html reading methods
 */
struct html_reader_data {
  FILE *fp;
  BUF_T *buffer;
  BMSTR_T **m_input; //matches <input type="hidden"
  BMSTR_T **m_name; //matches name="
  BMSTR_T **m_value; //matches value="
};

/*
 * creates the html reader data
 */
HTMLRD_T* create_html_reader_data(FILE *fp, BUF_T *buf) {
  HTMLRD_T *data = mm_malloc(sizeof(HTMLRD_T));
  data->fp = fp;
  data->buffer = buf;
  data->m_input = mm_malloc(sizeof(BMSTR_T*));
  data->m_input[0] = bmstr_create("<input type=\"hidden\"");
  data->m_name = mm_malloc(sizeof(BMSTR_T*));
  data->m_name[0] = bmstr_create("name=\"");
  data->m_value = mm_malloc(sizeof(BMSTR_T*));
  data->m_value[0] = bmstr_create("value=\"");
  return data;
}

/*
 * destroys the html reader data
 */
void destroy_html_reader_data(HTMLRD_T *data) {
  fclose(data->fp);
  buf_destroy(data->buffer);
  bmstr_destroy(data->m_input[0]);
  free(data->m_input);
  bmstr_destroy(data->m_name[0]);
  free(data->m_name);
  bmstr_destroy(data->m_value[0]);
  free(data->m_value);
}


/***********************************************************************
 * Function for reading a hidden input field from a html file.
 * assumes that the hidden input is layed out like:
 * <input type="hidden"{possible attributes}name="{field name}"{possible attributes}value="{field value}"{possible attributes}{optional /}>
 * so the attributes name and value must be in that order and not be repeated. If it is not as above then it probably won't work
 * This returns the name field and positions the buffer at the start of the value field
 ***********************************************************************/
char* read_hidden_input_name(HTMLRD_T* data) {
  char *name;
  int namelen;
  //find the hidden input field
  if (!find_bm_landmarks(data->fp, data->buffer, FALSE, TRUE, 1, data->m_input)) return NULL;
  //find the name
  if (!find_bm_landmarks(data->fp, data->buffer, FALSE, TRUE, 1, data->m_name)) return NULL;
  //read the name
  name = buf_fread_token(data->buffer, data->fp, is_quote, NULL, FALSE, NULL, 51, &namelen);
  if (namelen < 0)
    die("Meme html file input fields may not have name attributes larger than 50 characters.\n");
  if (name == NULL)
    die("File read error while attempting to read the name of a html input field. Error was given as %s.\n", strerror(ferror(data->fp)));
  //find the value
  if (!find_bm_landmarks(data->fp, data->buffer, FALSE, TRUE, 1, data->m_value)) {
    free(name);//can't find value so cleanup
    return NULL;
  }
  return name;
}


/***********************************************************************
 * Read a MEME html file.
 * version 4.3.2 and above
 ***********************************************************************/
static void read_meme_html_file(
   const char     *meme_filename,       // Name of MEME file  IN
   const char     *bg_filename,         // Name of background freq. file IN
   FILE           *motif_file,          // File ptr for file containing motifs IN
   BUF_T          *buffer,              // Buffer for motifs in
   double         pseudocount,          // Value of pseudocount IN
   int            loading_mode,         // should load PSSM or PSPM IN
   ARRAYLST_T*    motifs,               // The retrieved motifs
   STRING_LIST_T  **motif_occurrences,  // Strings desc. motif occurrences  OUT
   BOOLEAN_T      *has_reverse_strand,  // Does this file have both strands? OUT
   ARRAY_T        **background          // Background emission distribution  OUT
)
{
  //after the epoc of version 4.3.2 treat html differently to text
  int initial_motif_count;
  HTMLRD_T *htmld = create_html_reader_data(motif_file, buffer);
  BMSTR_T *freqstart[1];
  BMSTR_T *scorestart[1];
  int state = 0;
  int num_motifs = 0;
  char *field;
  int index, len;

  initial_motif_count = arraylst_size(motifs);

  freqstart[0] = bmstr_create("letter-probability matrix: ");
  scorestart[0] = bmstr_create("log-odds matrix: ");

  // read in each hidden field in the html file and parse depending on the state
  while ((field = read_hidden_input_name(htmld)) != NULL) {
    switch (state) {
      case 0:
        if (strcmp(field, "alphabet") == 0) {
          char *alphabet, alphabuf[100]; 
          int len;
          alphabet = buf_fread_token(buffer, motif_file, is_quote, NULL, FALSE, alphabuf, 100, &len);
          if (alphabet == NULL)
            die("Failed to read alphabet.\n");
          set_alphabet(verbosity, alphabet);
          if (which_alphabet() == DNA_ALPH) {
            state = 1;
          } else {
            *has_reverse_strand = FALSE;
            //see if we need to read the bg from the meme file
            if (!read_background(bg_filename, meme_filename, background)) {
              state = 2;//get bg from meme file
            } else {
              state = 3;//we got the bg elsewhere
            }
          }
        }
        break;
      case 1:
        if (strcmp(field, "strands") == 0) {
          char *strand, strbuf[10];
          int len;
          strand = buf_fread_token(buffer, motif_file, is_quote, NULL, FALSE, strbuf, 10, &len);
          if (strand == NULL)
            die("Failed to read strand.\n");
          *has_reverse_strand =  (strcmp(strand, "both") == 0);
          //see if we need to read the bg from the meme file
          if (!read_background(bg_filename, meme_filename, background)) {
            state = 2;//get bg from meme file
          } else {
            state = 3;//we got the bg elsewhere
          }
        }
        break;
      case 2:
        if (strcmp(field, "bgfreq") == 0) {
          *background = read_bg_freqs(motif_file, buffer);
          state = 3;
        }
        break;
      case 3:
        if (strcmp(field, "combinedblock") == 0) {
          if (motif_occurrences) {
            char *line;
            *motif_occurrences = new_string_list();
            //the line after the quote (") is the start
            buf_fread_next_line(buffer, motif_file);
            while (buf_remaining(buffer) || !feof(motif_file)) {
              line = buf_fread_token(buffer, motif_file, is_quote_or_nl, NULL, FALSE, NULL, MAX_LINE, &len); 
              if (line == NULL) die("Couldn't read combined block line.\n"); 
              if (buf_getc(buffer) == '"') {
                free(line);
                break;
              }
              add_string(line, *motif_occurrences);
              buf_fread_consume(buffer, motif_file, is_newline, NULL, FALSE);
              free(line);
            }
          }
          state = 4;
        }
        break;
      case 4:
        if (strcmp(field, "nmotifs") == 0) {
          int i;
          num_motifs = parse_pos_int(motif_file, buffer, "number of motifs");
          arraylst_preallocate(arraylst_size(motifs) + num_motifs, motifs);
          for (i = 0; i < num_motifs; i++) {
            MOTIF_T *a_motif = (MOTIF_T*)mm_malloc(sizeof(MOTIF_T));
            //set some defaults that won't be set in the meme file
            a_motif->trim_left = 0;
            a_motif->trim_right = 0;
            arraylst_add(a_motif, motifs);
          }
          state = 5;
        }
        break;
      case 8://look for the optional url
        if (strncmp(field, "url", 3) == 0 && (index = atoi(field+3)) > 0 && index <= num_motifs) {
          MOTIF_T *a_motif = arraylst_get(index-1, motifs);
          a_motif->url = buf_fread_token(buffer, motif_file, is_quote, NULL, FALSE, NULL, MAX_MOTIF_URL_LENGTH+1, NULL);
          if (a_motif->url == NULL) {
            die("File error while reading url for motif %d\n", index);
          }
          state = 5;
          break;
        }
      case 5://parse the block and get the id
        if (strncmp(field, "motifblock", 10) == 0 && (index = atoi(field+10)) > 0 && index <= num_motifs) {
          char *id1, idbuf[MAX_MOTIF_ID_LENGTH + 1];
          MOTIF_T *a_motif = arraylst_get(index-1, motifs);
          buf_fread_consume(buffer, motif_file, buf_is_delim_space, NULL, FALSE);
          expect(motif_file, buffer, "BL   MOTIF ");
          id1 = parse_word(motif_file, buffer, idbuf, MAX_MOTIF_ID_LENGTH + 1, 1, FALSE, "motif id");
          set_motif_id(id1, a_motif);
          set_motif_id2("", a_motif);
          a_motif->url = NULL;
          state = 6;
        }
        break;
      case 6://parse the position specific scoring matrix
        if (strncmp(field, "pssm", 4) == 0 && (index = atoi(field+4)) > 0 && index <= num_motifs) {
          MOTIF_T *a_motif = arraylst_get(index-1, motifs);
          find_bm_landmarks(motif_file, buffer, TRUE, TRUE, 1, scorestart);
          read_pssm(motif_file, buffer, loading_mode, a_motif);
          state = 7;
        }
        break;
      case 7://parse the position specific probability matrix
        if (strncmp(field, "pspm", 4) == 0 && (index = atoi(field+4)) > 0 && index <= num_motifs) {
          MOTIF_T *a_motif = arraylst_get(index-1, motifs);
          find_bm_landmarks(motif_file, buffer, TRUE, TRUE, 1, freqstart);
          read_pspm(motif_file, buffer, *background, pseudocount, loading_mode, a_motif);

          // Store the alphabet size in the motif. 
          a_motif->alph_size = get_alph_size(ALPH_SIZE);
          a_motif->ambigs = get_alph_size(AMBIG_SIZE);

          // Compute and store the motif complexity.
          a_motif->complexity = compute_motif_complexity(a_motif);

          if (!(loading_mode & WANT_PSPM)) {
            //the caller doesn't want the pspm, destroy it
            free_matrix(a_motif->freqs);
            a_motif->freqs = NULL;
          }
          state = 8;
        }
        break;
    }
    free(field);
  }
  for (index = arraylst_size(motifs) - num_motifs; index < arraylst_size(motifs); ++index) {
    MOTIF_T *a_motif = arraylst_get(index, motifs);
    if (a_motif->url == NULL) {
      set_motif_url("", a_motif);
    }
  }
  bmstr_destroy(freqstart[0]);
  bmstr_destroy(scorestart[0]);
  destroy_html_reader_data(htmld);
  if (state != 5 && state != 8) {
    die("Failed to parse properly\n");
  } 
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Read %d motifs.\n", (arraylst_size(motifs) - initial_motif_count));
  }
}

/***********************************************************************
 * Replicate motifs to an array
 ***********************************************************************/
static void replicate_motifs_to_array(ARRAYLST_T *list, MOTIF_T* array, int *num_motifs) {
  int i;
  *num_motifs = arraylst_size(list);
  for (i = 0; i < *num_motifs; ++i) {
    copy_motif((MOTIF_T*)arraylst_get(i, list), array+i);
  }
}

/***********************************************************************
 * Read a MEME file.
 ***********************************************************************/
void read_meme_file2(
   char*      meme_filename,          // Name of MEME file  IN
   char*      bg_filename,            // Name of background freq. file IN
   double     pseudocount,            // Value of pseudocount IN
   int        loading_mode,           // Load PSSM or PSPM IN
   ARRAYLST_T* motifs,                // The retrieved motifs IN OUT
   STRING_LIST_T** motif_occurrences, // Strings desc. motif occurrences  OUT
   BOOLEAN_T* has_reverse_strand,     // Does this file have both strands? OUT
   ARRAY_T**  background              // Background emission distribution  OUT
)
{
  // Open the given MEME file.
  FILE *motif_file; // MEME file containing the motifs.
  BUF_T *buffer;    // Buffer for the file.
  
  if (open_file(meme_filename, "r", TRUE, "motif", "motifs", &motif_file) == 0)
    exit(1);

  //initilize the buffer
  buffer = buf_create(100);
  //read only the first line, which will tell us if it's xml without making it impossible to parse the xml
  buf_fgets(buffer, motif_file);
  if (ferror(motif_file)) die("Error reading file \"%s\", error was given as %s\n", meme_filename, strerror(ferror(motif_file)));
  buf_flip(buffer);

  //see if the first line looks like xml
  if (!buf_unexpected(buffer, "<?xml", FALSE)) { 
    //looks like xml so try to read the MEME file as XML.
    BOOLEAN_T read_file = FALSE;
    read_file = read_meme_xml_file(
       meme_filename,
       bg_filename,
       pseudocount,
       loading_mode,
       motifs,
       motif_occurrences,
       has_reverse_strand,
       background
    );
    if (read_file) {
      buf_destroy(buffer);
      fclose(motif_file);
      return;
    }
    //uhoh, well this probably won't work but try reading as the other options
    //and hope that we're not reading from stdin...
  }

  // Check to be sure we can read this version.
  int vmajor, vminor, vbug, is_html;
  read_version(motif_file, buffer, &vmajor, &vminor, &vbug, &is_html);

  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Reading %s file vmajor: %d, vminor: %d, vbug: %d\n", (is_html ? "html" : "txt"), vmajor, vminor, vbug);
  }

  if (is_html && (vmajor > 4 || (vmajor == 4 && (vminor > 3 || (vminor == 3 && vbug >= 2))))) {
    read_meme_html_file(
       meme_filename,
       bg_filename,
       motif_file,
       buffer,
       pseudocount,
       loading_mode,
       motifs,
       motif_occurrences,
       has_reverse_strand,
       background
        );
    return;
  }

  // Read the alphabet.
  read_alphabet(motif_file, buffer);

  // Read the strandedness of the motifs.
  if (which_alphabet() == DNA_ALPH) {
    *has_reverse_strand = read_strand(motif_file, buffer);
  } else {
    *has_reverse_strand = FALSE;
  }

  if (!read_background(bg_filename, meme_filename, background)) {
    find_landmarks(motif_file, buffer, TRUE, TRUE, 1, LETTER_FREQUENCY);
    buf_fread_next_line(buffer, motif_file);
    *background = read_bg_freqs(motif_file, buffer);
  }

  // Read the specified motifs
  read_motifs(motif_file, buffer, is_html, *background, pseudocount, loading_mode, motifs);

  // Read the motif occurences
  if (motif_occurrences) {
    if (is_html) {
      read_motif_occurrences_from_old_html(motif_file, buffer, motif_occurrences);
    } else {
      read_motif_occurrences_from_txt(motif_file, buffer, motifs, motif_occurrences);
    }
  }
  /* Done with the buffer */
  buf_destroy(buffer);
  /* Close the MEME file. */
  fclose(motif_file);
}


/***********************************************************************
 * Read a MEME file.
 ***********************************************************************/
void read_meme_file(
   char*      meme_filename,          // Name of MEME file  IN
   char*      bg_filename,            // Name of background freq. file IN
   double     pseudocount,            // Value of pseudocount IN
   int        loading_mode,           // Load PSSM or PSPM IN
   int*       num_motifs,             // Number of motifs retrieved  OUT
   MOTIF_T*   motifs,                 // The retrieved motifs
   STRING_LIST_T** motif_occurrences, // Strings desc. motif occurrences  OUT
   BOOLEAN_T* has_reverse_strand,     // Does this file have both strands? OUT
   ARRAY_T**  background              // Background emission distribution  OUT
)
{
  int i;
  ARRAYLST_T *motif_list = arraylst_create();

  read_meme_file2(
      meme_filename,
      bg_filename,
      pseudocount,
      loading_mode,
      motif_list,
      motif_occurrences,
      has_reverse_strand,
      background
      );
  *num_motifs = arraylst_size(motif_list);
  if (*num_motifs > MAX_MOTIFS) {
    fprintf(stderr, "Warning: this method is limited to only loading %d motifs but the meme file contained %d motifs.", MAX_MOTIFS, *num_motifs);
    *num_motifs = MAX_MOTIFS;
  }
  //clear the destination or it might cause problems
  memset(motifs, 0, sizeof(MOTIF_T) * (*num_motifs));
  //copy over the motifs from the list
  for (i = 0; i < *num_motifs; ++i) {
    copy_motif((MOTIF_T*)arraylst_get(i, motif_list), motifs+i);
  }
  //free the list of motifs
  free_motifs(motif_list);
}

/***********************************************************************
 * free an array list and the contained motifs
 ***********************************************************************/
void free_motifs(
  ARRAYLST_T *motif_list
) {
  int i;
  MOTIF_T *tmp;
  for (i = 0; i < arraylst_size(motif_list); ++i) {
    tmp = (MOTIF_T*)arraylst_get(i, motif_list);
    free_motif(tmp);
    free(tmp);
    tmp = NULL;
  }
  arraylst_destroy(NULL, motif_list);
}

/***********************************************************************
 * Print out a given motif.
 ***********************************************************************/
static void write_motif(
  MOTIF_T* a_motif,     /* A data structure containing the motif. */
  FILE*    outfile      /* An ASCII version of the motif in MEME  */
) {                     /* output format.                         */
  int i;
  int j;

  fprintf(outfile, "MOTIF %s\n\n", a_motif->id);
  fprintf(outfile, "BL   MOTIF %s width=%d seqs=%6.3f\n", a_motif->id,
    a_motif->length, a_motif->num_sites);
  fprintf(outfile, "letter-probability matrix: alength= %d ",
    a_motif->alph_size);
  fprintf(outfile, "w= %d ", a_motif->length);
  fprintf(outfile, "nsites= %6.3f ", a_motif->num_sites);
  fprintf(outfile, "E= %g ", a_motif->evalue);
  fprintf(outfile, "complexity= %6.3f\n", a_motif->complexity);
  for (i = 0; i < a_motif->length; i++) {
    for (j = 0; j < a_motif->alph_size; j++) {
      fprintf(outfile, "%9.6f ", get_matrix_cell(i, j, a_motif->freqs));
    }
    fprintf(outfile, "\n");
  }
  fprintf(outfile, "\n\n");
}

#ifdef MAIN
VERBOSE_T verbosity = INVALID_VERBOSE;
#include "simple-getopt.h"

/*************************************************************************
 * int main
 *************************************************************************/
int main(int argc, char *argv[])
{
  // Data structures.
  int       num_motifs;         // The number of motifs in the model.
  MOTIF_T   motifs[2 * MAX_MOTIFS]; // The motifs.
  ARRAY_T*  background;         // Background probs for alphabet.
  ORDER_T*  order_spacing;      // Linear HMM order and spacing.

  // Command line parameters.
  char *    meme_filename;      // Input file containg motifs.
  BOOLEAN_T ids_only;           // Print only the motif ids?
  BOOLEAN_T reorder;            // Reorder the motifs?
  BOOLEAN_T has_reverse_strand; // MEME file contains both strands
  STRING_LIST_T* motif_occurrences; // Strings describing motif occurences.
  STRING_LIST_T* requested_motifs; // Indices of requested motifs.
  int       request_n;          // The user asked for the first n motifs.
  double    e_threshold;        // E-value threshold for motif inclusion.
  double    complexity_threshold; // For eliminating low complexity motifs.
  double    p_threshold;        // p-value threshold for motif occurences.
  BOOLEAN_T keep_unused;        // Keep unused transitions?
  char*     order_string;       // Motif order and spacing. *

  // Local variables.
  int       i_motif;
  int       alph_size;
  int       i_alph;

  /**********************************************
   * COMMAND LINE PROCESSING
   **********************************************/

  // Define command line options.
  cmdoption const options[] = {
    {"indices", NO_VALUE},
    {"reorder", NO_VALUE},
    {"motif", REQUIRED_VALUE},
    {"nmotifs", REQUIRED_VALUE},
    {"ethresh", REQUIRED_VALUE},
    {"lowcomp", REQUIRED_VALUE},
    {"pthresh", REQUIRED_VALUE},
    {"keep-unused", NO_VALUE},
    {"order", REQUIRED_VALUE},
    {"transpseudo", REQUIRED_VALUE},
    {"verbosity", REQUIRED_VALUE}
  };
  int option_count = 11;
  int option_index = 0;

  // Define the usage message.
  char      usage[1000] = "";
  strcpy(usage, "USAGE: meme-io [options] <motif file>\n");
  strcat(usage, "\n");
  strcat(usage, "   Options:\n");
  strcat(usage, "     --indices\n");
  strcat(usage, "     --reorder\n");
  strcat(usage, "     --motif <motif #> (may be repeated)\n");
  strcat(usage, "     --nmotifs <#>\n");
  strcat(usage, "     --ethresh <E-value>\n");
  strcat(usage, "     --lowcomp <value>\n");
  strcat(usage, "     --pthresh <p-value>\n");
  strcat(usage, "     --keep-unused\n");
  strcat(usage, "     --order <string>\n");
  strcat(usage, "     --verbosity 1|2|3|4|5 (default=2)\n");
  strcat(usage, "\n");

  // Initialize.
  order_spacing = NULL;

  // Make sure various options are set to NULL or defaults.
  meme_filename = NULL;
  ids_only = FALSE;
  reorder = FALSE;
  request_n = 0;
  requested_motifs = new_string_list();
  e_threshold = 0.0;
  complexity_threshold = 0.0;
  p_threshold = 0.0;
  keep_unused = FALSE;
  order_string = NULL;
  verbosity = NORMAL_VERBOSE;

  simple_setopt(argc, argv, option_count, options);

  // Parse the command line.
  while (1) {

    int c = 0;
    char* option_name = NULL;
    char* option_value = NULL;
    const char* message = NULL;

    // Read the next option, and break if we're done.
    c = simple_getopt(&option_name, &option_value, &option_index);
    if (c == 0) {
      break;
    } else if (c < 0) {
      simple_getopterror(&message);
      die("Error processing command line options (%s)\n", message);
    }

    if (strcmp(option_name, "indices") == 0) {
      ids_only = TRUE;
    } else if (strcmp(option_name, "reorder") == 0) {
      reorder = TRUE;
    } else if (strcmp(option_name, "motif") == 0) {
      add_string(option_value, requested_motifs);
    } else if (strcmp(option_name, "nmotifs") == 0) {
      request_n = atoi(option_value);
    } else if (strcmp(option_name, "ethresh") == 0) {
      e_threshold = atof(option_value);
    } else if (strcmp(option_name, "lowcomp") == 0) {
      complexity_threshold = atof(option_value);
    } else if (strcmp(option_name, "pthresh") == 0) {
      p_threshold = atof(option_value);
    } else if (strcmp(option_name, "keep-unused") == 0) {
      keep_unused = TRUE;
    } else if (strcmp(option_name, "order") == 0) {
      order_string = option_value;
    } else if (strcmp(option_name, "verbosity") == 0) {
      verbosity = (VERBOSE_T)atoi(option_value);
    }
  }

  // Read the single required argument.
  if (option_index + 1 != argc) {
    fprintf(stderr, "%s", usage);
    exit(1);
  }
  meme_filename = argv[option_index];

  // Set up motif requests.
  if (request_n != 0) {
    if (get_num_strings(requested_motifs) != 0) {
      die("Can't combine the -motif and -nmotifs options.\n");
    } else {
      for (i_motif = 0; i_motif < request_n; i_motif++) {
        char motif_id[MAX_MOTIF_ID_LENGTH + 1];
        sprintf(motif_id, "%d", i_motif + 1);
        add_string(motif_id, requested_motifs);
      }
    }
  }

  // Make sure motifs weren't selected redundantly.
  if ((get_num_strings(requested_motifs) != 0) && (e_threshold != 0.0)) {
    die("Can't use -motif or -nmotifs with -ethresh.");
  }
  if ((get_num_strings(requested_motifs) != 0) && (order_string != NULL)) {
    die("Can't use -motif or -nmotifs with -order.");
  }
  if ((order_string != NULL) && (e_threshold != 0.0)) {
    die("Can't use -ethresh and -order.");
  }

  // Parse the order string.
  order_spacing = create_order(order_string);

  /**********************************************
   * READING THE MOTIFS
   **********************************************/

  // Read all the motifs from the MEME file.
  read_meme_file(
     meme_filename,
     "motif-file", // Name of background freq. file
     0.0, // pseudocount
     WANT_PSSM | WANT_PSPM, //load the PSPMs
     &num_motifs,
     motifs,
     &motif_occurrences,
     &has_reverse_strand,
     &background
  );

  /**********************************************
   * WRITING THE MOTIFS
   **********************************************/
  if (ids_only) {
    for (i_motif = 0; i_motif < num_motifs; i_motif++) {
      printf("%s ", get_motif_id(&(motifs[i_motif])));
    }
    printf("\n");
  } else {

    printf("Meta-MEME version %s\n", VERSION);
    printf("ALPHABET= %s\n\n", get_alphabet(FALSE));
    printf("strands: +\n\n");

    // Print the background.
    printf("%s\n", LETTER_FREQUENCY);
    {
      int alph_size;
      int i;

      alph_size = get_alph_size(ALPH_SIZE);
      for (i = 0; i < alph_size; i++) {
        printf("%c %5.3f ", get_alph_char(i), get_array_item(i, background));
      }
      printf("\n\n");
    }

    // Print the requested motifs.
    for (i_motif = 0; i_motif < num_motifs; i_motif++) {
      write_motif(&(motifs[i_motif]), stdout);
    }

    // Print the letter frequencies.
    printf("Letter frequencies:\n");
    alph_size = get_alph_size(ALPH_SIZE);
    for (i_alph = 0; i_alph < alph_size; i_alph++) {
      printf("%c %5.3f ", get_alph_char(i_alph),
      get_array_item(i_alph, background));
    }
    printf("\n");
  }

  // Clean up and exit.
  free_string_list(requested_motifs);
  free_order(order_spacing);
  for (i_motif = 0; i_motif < num_motifs; i_motif++) {
    free_motif(&(motifs[i_motif]));
  }
  free_array(background);
  return(0);
}
#endif


/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */

