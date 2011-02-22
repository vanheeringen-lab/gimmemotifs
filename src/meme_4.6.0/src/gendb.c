/**************************************************************************
 * FILE: gendb.c
 * AUTHOR: Timothy L. Bailey
 * CREATE DATE: 10/29/02
 * PROJECT: MHMM
 * COPYRIGHT: 2002, TLB
 * DESCRIPTION: Generate sequences using a Markov model.
 **************************************************************************/

#ifdef MAIN
#define DEFINE_GLOBALS
#endif
#include "macros.h"
#include "background.h"
#include "fcodon.h"
#include "gendb.h"
#include "hash_alph.h"
#include "seq.h"

#define LOWEST 50			/* shortest sequence */
#define HIGHEST 2000			/* longest sequence */
#define SEQS 10				/* number of sequences */
#define MAX_ORDER 10			// largest Markov model order

/**************************************************************************/
/*
	get_letters

	Convert an alphabet into a set of letter strings.
	Get the size of alphabet, size of strings and default 0-order
	Markov frequencies.
*/
/**************************************************************************/
char **get_letters(
  int type,			// Type of alphabet.
  char **alph_p,		// Alphabet.
  int *r_p,			// Alphabet length.
  int *c_p,			// Length of letter string.
  double **f_p			// 0-order Markov frequencies.
)
{
  int i, j, k;
  char *alph, *old_alph, *ptr;
  int r, c, alen, old_alen;
  double *f = NULL;
  char **letters;			/* letters/codons to print */

  /* set up alphabet, frequency table, letters/codon array */
  switch (type) {
    case 0:
      alph = PROTEINB; f = nrfreq; alen = strlen(alph); r = alen; c = 2;
      create_2array(letters, char, r, c);
      for (i=0; i<r; i++) {letters[i][0] = alph[i]; letters[i][1] = '\0';}
      break;
    case 1:
      alph = DNAB; f = ntfreq; alen = strlen(alph); r = alen; c = 2;
      create_2array(letters, char, r, c);
      for (i=0; i<r; i++) {letters[i][0] = alph[i]; letters[i][1] = '\0';}
      break;
    case 2:
      alph = DNAB; f = fcodon; alen = strlen(alph); r = alen*alen*alen; c = 4;
      create_2array(letters, char, r, c);
      for (i=0; i<alen; i++) {			/* create codon string */
        for (j=0; j<alen; j++) {
          for (k=0; k<alen; k++) {
            int l = k+(alen * (j + (alen * i)));
            letters[l][0] = alph[i]; letters[l][1] = alph[j];
	    letters[l][2] = alph[k]; letters[l][3] = '\0';
          }
        }
      }
      break;
    case 3:
      alph = DNA0; alen = strlen(alph); r = alen; c = 2;
      create_2array(letters, char, r, c);
      for (i=0; i<r; i++) {letters[i][0] = alph[i]; letters[i][1] = '\0';}
      /* remove ambiguity characters from frequencies */
      mm_resize(f, r, double);
      old_alph = DNAB;
      old_alen = strlen(old_alph);
      for (i=0; i<alen; i++) f[i] = 0;
      for (i=0; i<old_alen; i++) {
        if ((ptr = strchr(alph, old_alph[i]))) {
          f[ptr-alph] += ntfreq[i];
        } else {
          /* divide ambiguous frequencies evenly among letters */
          for (j=0; j<alen; j++) { f[j] += ntfreq[i]/alen; }
        }
      }
      break;
    case 4:
      alph = PROTEIN0; alen = strlen(alph); r = alen; c = 2;
      create_2array(letters, char, r, c);
      for (i=0; i<r; i++) {letters[i][0] = alph[i]; letters[i][1] = '\0';}
      /* remove ambiguity characters from frequencies */
      mm_resize(f, r, double);
      old_alph = PROTEINB;
      old_alen = strlen(old_alph);
      for (i=0; i<alen; i++) f[i] = 0;
      for (i=0; i<old_alen; i++) {
        if ((ptr = strchr(alph, old_alph[i]))) {
          f[ptr-alph] += nrfreq[i];
        } else {
          /* divide ambiguous frequencies evenly among letters */
          for (j=0; j<alen; j++) { f[j] += nrfreq[i]/alen; }
        }
      }
      break;
    default:
      fprintf(stderr, "Unknown data type %d.\n", type);
      exit(1);
  }

  // Set up return.
  *alph_p = alph;
  *r_p = r;
  *c_p = c;
  *f_p = f;
  return(letters);
} // get_letters

/**************************************************************************/
/*
	get_cum_distr

	Get a cumulative distributions for each prefix "w".
	In other words
		sum_{x<a} Pr(wa)
	where "w" is a prefix and "x" and "a" are single characters.

	Reads a Markov model in bfile format from file bfile if given.
	Otherwise, uses the Markov model in array f().
*/
/**************************************************************************/
double *get_cum_distr(
  char *bfile,			// Markov model file.
  double *f,			// Markov model.
  char *alph,			// Alphabet of model.
  int  r,			// Alphabet length.
  int use_order,		// Order desired for model; use actual order if -1.
  int *order			// Order of returned model.
)
{
  int i, n;
  double *cum = NULL;		// Cumulative probability distribution.

  //
  // Get Markov model.
  //
  if (bfile == NULL) {			// Get from array.

    /* letter/codon distribution */
    mm_resize(cum, r, double);
    for (i=0; i<r; i++) cum[i] = f[i];
    *order = 0;

  } else {				// Get from -bg file.

    // Get the conditional probabilities.
    cum = read_markov_model(bfile, NULL, alph, FALSE, FALSE, order);

    if (*order > MAX_ORDER) {
      fprintf(stderr, "Order of Markov model %d exceeds MAX_ORDER %d.\n",
	*order, MAX_ORDER);
      exit(1);
    }
    if (*order < use_order) {
      fprintf(stderr, "The <bfile> is only order %d so you can't use order %d.\n",
	*order, use_order);
      exit(1);
    }

    *order = use_order < 0 ? *order : use_order;
  } // Cumulative distribution.


  // Get the number of probabilities.
  n = r;
  for (i=1; i<=*order; i++) n += pow(r, i+1);

  // Convert the conditional probabilites to cumulative probs.
  for (i=0; i<n; i++) {
    if (i%r == r-1) {
      cum[i] = 1.0;				// Last letter.
    } else if (i%r != 0) {
      cum[i] += cum[i-1];			// Not first letter.
    }
    if (0) {
      char *string;
      string = i2s(i);
      fprintf(stderr, "i %d word %s cum %f\n", i, string, cum[i]);
      free(string);
    }
  }

  return(cum);
} // get_cum_distr

/**************************************************************************/
/*
	print_random_seqs

	Print a number of random sequences.  Each sequence has
	random length and composition determined by the cumulative
	distribution input.

	Returns a single seq record (instead of printing) if
	out is NULL.
*/
/**************************************************************************/
SEQ_T *print_random_seqs (
  FILE *out,				// Stream to print on.
  int seed,				// Random number seed.
  int nseqs, 				// Number of sequences to print.
  int min,				// Minimum sequence length.
  int max,				// Maximum sequence length.
  char **letters,			// Array of letter strings.
  int r,				// Number of letter strings.
  int c,				// Length of letter strings.
  int order,				// Order of Markov model.
  double *cum				// Cumulative distribution(s) defining model.
)
{
  int i, j;
  int n;				// Length of sequence.
  char *buffer = NULL;			// Buffer for sequences.
  char *id = NULL;			// Sequence name.
  char first_letter = letters[0][0];	// First letter in alphabet.

  // Create the buffer for the string.
  mm_resize(buffer, max*(c-1)+1, char);

  /* set up random number generator */
  if (seed != 0) srand48(seed);

  /* print random sequences */
  for (i=0; i<nseqs; i++) {			/* sequence */

    // Decide length of sequence to print.
    n = (int) (min + drand48() * (max - min + 1));

    // Print FASTA ID line.
    if (out != NULL) {
      fprintf(out, ">SEQ_%-d %d\n", i+1, n);
    } else {
      mm_resize(id, 50, char);
      sprintf(id, ">SEQ_%-d %d\n", i+1, n);
    }

    /*
      Generate letters by
		1) random x ~ [0,1)
        2) binary search of cum for cum[i-1]<x<=cum[i]
		3) letter/codon is letters[i-1]
    */
    for (j=0; j<n; j++) {			/* generate letters/codons */
      double x = drand48();			/* random number */
      int lo = 0;
      int hi = r;
      int offset = 0;				// Offset into cum array.

      if (order >= 1) {				// Markov model.
        int start_ptr;
        // Find the offset into the cumulative prob array by looking
	// for the offset of the preceeding "order" characters.
        buffer[j] = first_letter;		// Now contains index into array.
        buffer[j+1] = '\0';
        start_ptr = j > order ? j-order : 0;	// Start of index string.
        offset = s2i(buffer + start_ptr);
        //fprintf(stderr, "b: %s\n offset: %d x: %f\n", buffer+start_ptr, offset, x);
      }

      while (hi-lo > 1) {			/* binary search */
        int mid = (lo+hi)/2;			/* midpoint */
        if (x > cum[mid+offset]) { lo = mid; } else { hi = mid; }
      }
      //fprintf(stderr, "%11.8f %s\n", x, letters[x<cum[lo+offset] ? lo : lo+1]);
      //fprintf(stderr, "%s", letters[x<cum[lo+offset] ? lo : lo+1]);
      buffer[j] = letters[x<cum[lo+offset] ? lo : lo+1][0];
    } /* generate letters/codons */
    buffer[j] = '\0';

    // Print the sequence.
    if (out != NULL) {
      for (j=0; j<n; j+=50) {
	fprintf(out, "%-50.50s\n", buffer+j);
      }
    } else {
      SEQ_T *seq = allocate_seq(id, "", 0, buffer);
      set_complete(TRUE, seq);
      return(seq);
    }

    if (i%100 == 0) fprintf(stderr, "%d\r", i);
  } /* sequence */

  return(NULL);
} // print_random_seqs

/**************************************************************************/
/*
	gendb

	Generate a file of synthetic sequences.
*/
/**************************************************************************/
SEQ_T *gendb(
  FILE *out,			// Output stream; return output if null.
  int type,			// Type of alphabet.
				// 0: protein w/ambigs
				// 1: dna w/ambigs
				// 2: codons
				// 3: dna w/o ambigs
				// 4: protein w/o ambigs
  char *bfile,      		// Name of Markov model file.
  int use_order,		// Order of Markov model to use.
  double *f,			// 0-order model; used if bfile==NULL.
  int nseqs,			// Number of sequences.
  int min,			// Shortest sequence.
  int max,			// Longest sequence.
  int seed 			// Random seed.
)
{
  char *alph;				/* alphabet */
  double *def_f;			/* letter or codon frequencies */
  double *cum;				/* letter/codon cumulative distr. */
  char **letters;			/* letters/codons to print */
  int r;				/* number of letters/codons */
  int c;				/* length of letter/codon strings */
  int order;

  // Get the letters and alphabet stuff.
  letters = get_letters(type, &alph, &r, &c, &def_f);
  if (f == NULL) f = def_f;

  // Get the cumulative distribution(s).
  cum = get_cum_distr(bfile, f, alph, r, use_order, &order);

  // Print the random sequences.
  return (print_random_seqs(out, seed, nseqs, min, max, letters, r, c, order, cum));
} // gendb

#ifdef MAIN
#include "simple-getopt.h"

int main(int argc, char **argv) {

  // Defaults for command line options.
  char* bfile     = NULL;               // Background file.
  int seed        = 0;                  // For random numbers.
  int min         = LOWEST;             // Minimum sequence length.
  int max         = HIGHEST;            // Maximum sequence length.
  int use_order   = -1;			// Use order defined in bfile.
  BOOLEAN_T type  = 0;			// (See usage statement.)
  BOOLEAN_T dummy = FALSE;		// Don't print dummy sequence.

  // Define command line options.
  const int num_options = 7;
  cmdoption const options[] = {
    {"bfile",  REQUIRED_VALUE},
    {"seed",   REQUIRED_VALUE},
    {"minseq", REQUIRED_VALUE},
    {"maxseq", REQUIRED_VALUE},
    {"order",  REQUIRED_VALUE},
    {"type",   REQUIRED_VALUE},
    {"dummy",  NO_VALUE}
  };
  int option_index = 0;
  simple_setopt(argc, argv, num_options, options);

  // Define the usage message.
  char      usage[1000] = "";
  strcpy(usage, "USAGE: gendb [options] <numseqs>\n");
  strcat(usage, "\n");
  strcat(usage, "   Options:\n");
  strcat(usage, "     --bfile <file>\n");
  strcat(usage, "     --seed <int> (default=0)\n");
  strcat(usage, "     --minseq <int> (default=50)\n");
  strcat(usage, "     --maxseq <int> (default=2000)\n");
  strcat(usage, "     --order <int> use Markov model of given order\n");
  strcat(usage, "     --type 0|1|2|3|4\n");
  strcat(usage, "            0=protein (default)\n");
  strcat(usage, "            1=dna with ambiguous characters\n");
  strcat(usage, "            2=codons\n");
  strcat(usage, "            3=dna w/o ambiguous characters\n");
  strcat(usage, "            4=protein w/o ambiguous characters\n");
  strcat(usage, "     --dummy\n");

  // Parse the command line.
  while (1) {
    int c = 0;
    char* option_name = NULL;
    char* option_value = NULL;
    const char* message = "";

    // Read the next option, and break if we're done.
    c = simple_getopt(&option_name, &option_value, &option_index);
    if (c == 0) {
      break;
    } else if (c < 0) {
      die("Error processing command line options (%s)\n", message);
    }

    if (strcmp(option_name, "bfile") == 0) {
      bfile = option_value;
    } else if (strcmp(option_name, "seed") == 0) {
      seed = atoi(option_value);
    } else if (strcmp(option_name, "minseq") == 0) {
      min = (int)atof(option_value);
    } else if (strcmp(option_name, "maxseq") == 0) {
      max = (int)atof(option_value);
    } else if (strcmp(option_name, "order") == 0) {
      use_order = atoi(option_value);
    } else if (strcmp(option_name, "type") == 0) {
      type = atoi(option_value);
    } else if (strcmp(option_name, "dummy") == 0) {
      dummy = TRUE;
    }
  }
  if (option_index + 1 != argc) {
    fprintf(stderr, "%s", usage);
    exit(EXIT_FAILURE);
  }
  int nseqs = atoi(argv[option_index]);
  fprintf(stderr, "Generating %d sequences.\n", nseqs);

  // Print dummy sequence if asked to.
  if (dummy) {
    printf(">SEQ_0 n= %d  s= %d  l= %d  seed= %d\n", nseqs, min, max, seed);
  }

  gendb(
	stdout,
	type,
	bfile,
	use_order,
	NULL, // 0-order model.  -- WSN 8/13/03
	nseqs,
	min,
	max,
	seed
	);

  exit(0);
} // main
#endif // MAIN
