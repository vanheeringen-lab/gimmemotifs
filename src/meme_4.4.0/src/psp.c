#include "psp.h"
#include <assert.h>
#include <math.h>

/* local functions */
static BOOLEAN read_psp(
  FILE *data_file,              /* file containing PSP data */
  char **sample_name,           /* unique identifier of saple == positional prior */
  int *w,			/* PSP motif width */
  double **pospriors,           /* positional prior probabilities */
  long int *length		/* number of priors == sequence length */
);
static long read_psp_data(
  FILE *data_file,		/* data file of sequences */
  double **pospriors,		/* posprior data */
  char *name 			/* name of sequence */
);
static char *numberstring(
  FILE *data_file,		/* data file */
  BOOLEAN point,		/* allow '.' at start? */
  char *name,			/* name of sequence */
  BOOLEAN skipToEOL             /* if TRUE skip rest of line */
);
static BOOLEAN all_zeros (
  BOOLEAN forward,
  double* data_end,
  int N
);

/* local macro */

/* assume ASCII or similar ordering of '0'..'9' */
#define numprefix(c) ((c >='0' && c <= '9') || (c == '.') || (c == '+') || (c == '-'))

/*
POSPRIORS format:
 [
 <header>
 <doubles>+
 ]*
 
 <header>       := ">"<id>" "<w>"
 Everything up to the first space is 
 assumed to be the unique name of the sequence.
 The w is the motif width for which the priors were scored.
 Priors contain P_{i,j} and are organized as follows: 
   -psp2: j =	..-2,-1,0,+1,+2,..	[position "0" in middle should be zero]
   -psp:  j = 	+1,+2..

 Only one instance of a unique name is permitted, and each unique name must
 match a previously created sequence name. Any previously created sequences
 or which a prior value is not read in has priors initialized to uniform
 probability (see complete_positional_priors() in init.c).
*/

/**********************************************************************/
/*
 read_psp_file
 
 Read a MEME PSP file.

 Creates all the priors and add them to dataset. Do this after reading
 in the sequences because the priors have to be matched to similarly 
 named sequences.

*/
/**********************************************************************/
extern void read_psp_file(
  char *psp_filename,		// PSP file name
  DATASET *dataset,		// the dataset
  BOOLEAN psp_revcomp, 		// PSP file contains both strands
  BOOLEAN meme_revcomp,		// MEME using both strands
  MOTYPE mtype			// model type
)
{
  int i, new_w;
  long int length, correct_length;
  SAMPLE *sample = NULL;	// look up in hash table based on ID
  char *sample_name = NULL;
  double *pospriors = NULL;
  int N = 0;
  BOOLEAN first_record = TRUE;

  FILE *prior_file = fopen(psp_filename,"r");
  if (!prior_file) {
    fprintf(stderr, "Failed to open PSP file: %s.\n", psp_filename);
    exit(1);
  }

  // Read in each PSP record and update the information in dataset struct.
  while (read_psp(prior_file, &sample_name, &new_w, &pospriors, &length)) {

    // Die if there was an error, which is signaled by length < 0.
    N++;
    if (length < 0) {
      fprintf(stderr, "Error occurred in PSP record number %d.\n", N);
      exit(1);
    }

    // Check that width hasn't changed.
    if (first_record) {
      dataset->psp_w = new_w;
      first_record = FALSE;
    } else if (new_w != dataset->psp_w) {
      fprintf(stderr, "PSP width W: `%d' in PSP record with ID: `%s'.\n"
        "doesn't match initial record's W: `%d'.\n", new_w, sample_name, dataset->psp_w);
      exit(1);
    }

    // Check that record ID matches a sequence ID.
    if ((sample=get_sample_by_name(sample_name)) == NULL) {
      fprintf(stderr, "Warning: PSP record ID: `%s' doesn't match any sequence's ID.\n",
        sample_name);
      //exit(1);
      continue;
    }

    // Check that the correct number of entries are in PSP record.
    correct_length = psp_revcomp ? 2*sample->length+1 : sample->length;
    if (length != correct_length) {
      fprintf(stderr, "Number of PSP values for sequence ID: `%s' is incorrect.\n"
	"Is: %ld but should be: %ld.\n", sample_name, length, correct_length);
      exit(1);
    }

    // shift -psp prior right and add zero at start
    if (!psp_revcomp) {
      double *psp = (double *) mymalloc((length+1) * (int) sizeof(double));
      psp[0] = 0;
      for (i=0; i<length; i++) psp[i+1] = pospriors[i];
      myfree(pospriors);
      pospriors = psp;
      length += 1;
    }

    // Check that all entries in prior are in [0,1].
    // Check that double-stranded prior is symmetrical.
    // Check that prior sums to 0<sum<=1.
    double psp_sum = 0;
    for (i=0; i<length; i++) {
      psp_sum += pospriors[i];
      if (pospriors[i] < 0 || pospriors[i]>1) {		// Illegal entry?
	fprintf(stderr, "Illegal value (%f) in PSP record %s at position %d.\n"
	  "All entries must be in range [0..1]\n", pospriors[i], sample_name, i);
        exit(1);
      }
      if (psp_revcomp && i < length/2) {		// Non-symmetrical?
        if (pospriors[i] != pospriors[length-1-i]) {
          double avg = (pospriors[i] + pospriors[length-1-i])/2.0;
	  fprintf(stderr, "Warning: Position %d in double-stranded PSP record %s\n"
            "not symmetrical: (+) %g (-) %g.\n"
            "Setting values on opposite strands to their average: %g.\n", 
            i, sample_name, pospriors[i], pospriors[length-1-i], avg);
          pospriors[i] = pospriors[length-1-i] = avg;
        }
      }
    }

    // check that 0 < sum <= 1
    double epsilon = 1e-6;
    if (psp_sum <= 0 || psp_sum>(1+epsilon)) {
      fprintf(stderr, "PSP record %s does not sum to 0<sum<=1 (%f).\n",
        sample_name, psp_sum);
      exit(1);
    }

    //
    // Set length to sample length
    //
    length = sample->length;

    // Add priors for negative strand if needed, scaled by 1/2.  
    if (meme_revcomp && !psp_revcomp) {
      double *psp = length + (double *) mymalloc((2*length+1) * (int) sizeof(double));
      psp[0] = pospriors[0]; 	// middle element is now probability of no site
      for (i=1; i<=length; i++) psp[i] = psp[-i] = pospriors[i]/2;
      myfree(pospriors);	// replace pospriors with psp
      pospriors = psp;
    }

    // Combine positive and negative strand priors if scanning only one strand.
    if (!meme_revcomp && psp_revcomp) {
      for (i=1; i<=length; i++) {
        pospriors[length+i] += pospriors[length-i];
        pospriors[length-i] = 0;
      }
    }

    // Set 0 index of psp_original to P_{i,0}
    if (meme_revcomp && psp_revcomp) { 	// index from middle
      sample->psp_original = &(pospriors[length]);
    } else {	// index from left end or middle set already
      sample->psp_original = pospriors;
    }

    // set the probability of no site and
    // make prior sum to 1 if Oops model
    double diff = 1 - psp_sum;
    if (mtype == Oops) {
      double factor = 1/(1 - diff);
      sample->psp_original[0] = 0;
      for (i=1; i<=length; i++) {
        sample->psp_original[i] *= factor;
        if (meme_revcomp || psp_revcomp) sample->psp_original[-i] *= factor;
      }
    } else if (mtype == Zoops) {
      sample->psp_original[0] = diff;
    }

  } // read next record

} /* read_psp_file */

/**********************************************************************/
/*
	psp_renormalize

	Initialize the log_psp array for each sequence from its
	psp_original array.

  	Method:

  	If the w_log_psp width is different from new_w,
  	set the log_psp array with (log P_ij)) position specific prior values.

  	If there is no PSP (s->psp_original==NULL), set up the uniform prior.

  	Otherwise, if new_w > psp_w, P_ij is the geometric mean of the original 
  	priors for width=psp_w sites it totally contains, 
  	normalized to sum to 1-P_i0.

  	If new_w < psp_w, P_ij is the original prior value for width psp_w,
	normalized to sum to 1-P_i0.


 */
/**********************************************************************/
extern void psp_renormalize(
  DATASET *dataset,		/* the dataset */
  int new_w,			/* new motif width */
  BOOLEAN revcomp,		/* reverse complement? */
  MOTYPE mtype			/* OOPS, ZOOPS or TCM? */
)
{
  // check for improper/unneccesary call
  assert(mtype == Oops || mtype == Zoops);	// PM FIXME need to add Tcm

  if (dataset->log_psp_w == new_w) return;	// nothing needed

  int n_samples = dataset->n_samples;
  int psp_w = dataset->psp_w;			// original psp width
  int i, j, k;
  for (i=0; i<n_samples; i++) { 		// loop over sequences
    SAMPLE *sample = dataset->samples[i];	// sequence record
    double *psp_original = sample->psp_original;// input PSP
    double *log_psp = sample->log_psp;		// renormalized log of PSP
    int length = sample->length;		// sequence length
    int new_m = length-new_w+1;			// number of possible new positions

    if (new_m < 1) continue;                    // skip sequence too short for motif

    // compute log_psp for this sequence
    if (psp_original == NULL) {			// Create a uniform prior.
      int ns = revcomp ? 2*new_m : new_m;	// possible sites on 1 or 2 strands
      double uniform = LOG(1.0/ns);		// (log) uniform prior
      log_psp[0] = LOG(0);			// for cleanliness
      sample->max_log_psp = uniform;		// always set
      if (revcomp) {
        for (j=1; j<=new_m; j++) log_psp[j] = log_psp[-j] = uniform;
      } else {
        for (j=1; j<=new_m; j++) log_psp[j] = uniform;
      }

    } else {					// Renormalize the input PSP.
      // set (log) P_i0 
      log_psp[0] = LOG(psp_original[0]);	

      // set rest of P_ij
      if (new_w <= psp_w) {			// width same or smaller
	//
        // use prior for original width
	//
        int old_m = length-new_w+1;		// number of possible old sites
        for (j=1; j<=old_m; j++) {		// copy original prior
          log_psp[j] = psp_original[j];
          if (revcomp) log_psp[-j] = psp_original[-j];
        }
        // set remaining legal sites priors to zero
        for (j=old_m; j<=new_m; j++) {
          log_psp[j] = 0;
          if (revcomp) log_psp[-j] = 0;
        }

      } else if (new_w > psp_w) {		// new w is larger
	//
        // compute geometric-meme-of-overlapped sites prior
	//
        int n_overlapped = MIN(psp_w, new_w - psp_w + 1);	// keep it real
        // compute product of overlapped sites priors
        for (j=1; j<=new_m; j++) {
          log_psp[j] = 1;
          if (revcomp) log_psp[-j] = 1;
          for (k=0; k<n_overlapped; k++) {
            log_psp[j] *= psp_original[j+k];
            if (revcomp) log_psp[-j] *= psp_original[-(j+k)];
          }
        }
        // compute geometric mean
        for (j=1; j<=new_m; j++) {
          log_psp[j] = pow(log_psp[j], 1.0/n_overlapped);
          if (revcomp) log_psp[-j] = pow(log_psp[-j], 1.0/n_overlapped);
        }
      } // new_w

      //
      // normalize so the total prior probability is 1-P_i0, and take log.
      //
      // get current sum
      double psp_sum = 0;
      for (j=1; j<=new_m; j++) {
        psp_sum += log_psp[j];
        if (revcomp) psp_sum += log_psp[-j];
      }
      // set scale factor
      double factor = (1 - psp_original[0]) / psp_sum;
      sample->max_log_psp = LOG(0);		
      // normalize and get maximum log_psp
      for (j=1; j<=new_m; j++) {
        log_psp[j] = LOG(factor * log_psp[j]);
        if (revcomp) log_psp[-j] = LOG(factor * log_psp[-j]);
        // record maximum log_psp in this sequence
        if (log_psp[j] > sample->max_log_psp) sample->max_log_psp = log_psp[j];
        if (revcomp && log_psp[-j] > sample->max_log_psp) sample->max_log_psp = log_psp[-j];
      }

    } // not uniform

    // set P_ij for illegal sites at end
    for (j=new_m+1; j<=length; j++) {
      log_psp[j] = LOG(0);
      if (revcomp) log_psp[-j] = LOG(0);
    }

  } // sequence loop

  // Renormalization complete
  dataset->log_psp_w = new_w;			// new log_psp width

}  // psp_renormalize

/**********************************************************************/
/*
	add_psp_to_log_not_o
*/
/**********************************************************************/
extern void add_psp_to_log_not_o(
  DATASET *dataset,		/* the dataset */
  int w,			/* motif width */
  BOOLEAN invcomp,		/* reverse complement? */
  MOTYPE mtype 			/* model type */
)
{
  int i,j;
  int n_samples = dataset->n_samples;
  SAMPLE **samples = dataset->samples;

  // Renormalize the PSP if width has changed.
  psp_renormalize(dataset, w, invcomp, mtype);

  for (i=0; i<n_samples; i++){		/* sequence */
    SAMPLE *s = samples[i];
    int lseq = s->length;
    double *weights = s->weights;       /* prb not in a previous site */
    double *not_o = s->not_o;		/* prb not overlapping a site */
    int *log_not_o = s->log_not_o;	/* log prb not overlapping a site */
    double max_log_psp = s->max_log_psp;	/* max pos-specific prior prob */
    double *log_psp = s->log_psp;	// log PSP

    for (j=0; j<=lseq-w; j++) { 	/* site start */

      // Sanity check: assuming that PSPs are always symmetrical
      if (invcomp &&  (log_psp[j+1] != log_psp[-(j+1)])) {
        fprintf(stderr, "add_psp_to_log_not_o: site %d differs: j %g -j %g\n", 
          j, log_psp[j+1], log_psp[-(j+1)]);
        exit(1);
      }

      // Scale PSP linearly so that maximum value is 1.0 
      // by dividing by the maximum psp value.
      int log_scaled_psp = (int) (SCALE_LOGS*(log_psp[j+1] - max_log_psp));

      // log_not_o: log(site not overlapped) + log(scaled_prior)
      log_not_o[j] = INT_LOG(not_o[j]) + log_scaled_psp;
     
    } /* for j */

  } /* for i */

} // add_psp_to_log_not_o

/**********************************************************************/
/*
 read_psp
 
 Read position-specific priors for a single sequence from the data file.
 Returns FALSE on EOF or bad sequence data.
 
 If PSP format, expect length to match equivalently named sequence length.
 If PSP2 format, expect length = 2*sequence length.
 
 Based on FASTA format.
 
 Note:
 checks after calling:
 1. sample name matches an existing sample
 2. as does the length
 3. all calls return same w
 4. depending on whether PSP or PSP2 format, check whether either
 w-1 zeros at end only or w-1 zeros at end and end of first length
 entries
 */
/**********************************************************************/

static BOOLEAN read_psp(
  FILE *psp_file,		/* file containing PSP */
  char **sample_name,           /* unique ID of sample == positional prior */
  int *w,                       /* width of motifs in prior */
  double **pospriors,           /* positional prior probabilities */
  long int *length		/* number of priors == sequence length */
)
{
  int i, c;
  char *name = NULL;

  *length = -1;			// In case of error, return TRUE with illegal length set.

  /* skip anything until first positional prior name */
  c = ' '; 
  while(c != EOF) { 
    if((c=fgetc(psp_file)) == '>') {  /* FASTA format */
      break;
    } 
    Skip_eol(c, psp_file);		/* go to end of line */
  }
  if (c==EOF) return FALSE;		/* no more sequences */
  
  /* get the sample name */
  /* read to first blank/tab/ or end of line/file */
  for (i=0; (c=fgetc(psp_file))!=EOF; ) {
    if ((c==' ') && i==0) {	/* skip blanks until name starts */
      continue;
    } else if (c==' ' || c=='\t' || c=='\n' || c=='\r') {
      break;				/* blank or nl ends name */
    } else {
      if ((i % RCHUNK) == 0) {
	Resize(name, i+RCHUNK, char);
      }
      name[i++] = c;			/* non-blank: add to name */
    }
  }
  Resize(name, i+1, char);
  name[i] = '\0';
  
  /* read in W (width of predicted motif; expect >= W-1 zeroes at end) */
  if (c != '\n' && c != '\r' && c != EOF) { 
    char *number;
    Skip_whi(c, psp_file);		/* skip whitespace to W */
    number = numberstring(psp_file, FALSE, name, TRUE); // ignore anything after W

    // Check that W is a number and is legal (>0)
    if (number==NULL || sscanf(number, "%d", w) != 1 || *w <= 0) {
      fprintf(stderr, "\nIllegal value for W: `%s' in PSP record with ID: `%s'.\n", 
        number, name);
      if (number != NULL) myfree(number);
      return(TRUE);			// bad W
    }
  }
  
  /* read in the actual prior data */
  *length = read_psp_data(psp_file, pospriors, name);
  
  /* sequence had bad data */
  if (*length < 0) {
    myfree(name);
    return FALSE;			// bad prior data format
  }
  
  *sample_name = name;
  
  return TRUE;
} /* read_psp */

/**********************************************************************/
/* 
 read_psp_data
 
 Read the positional prior data into a dynamic array.
 Assumes after the ">" line, and ends either at EOF or before next ">"
 Skips all white space.
 Checks for other characters disallowed in C floating-point numbers.
 
 Returns the number of positional priors or -1 on error.
 */
/**********************************************************************/
static long read_psp_data(
  FILE *psp_file,		/* data file of PSP sequences */
  double **pospriors, 		/* posprior data */
  char *name 			/* name of sequence */
)
{
  long length = 0, /* length of each number */
    N; /* seq length or 2x seq length if double-stranded */
  int c;
  char *number = NULL;
  double *newprior = NULL;
  double nextnumber;
  
  /* 
     read positional priors 
  */
  for(N=0; (c=fgetc(psp_file))!=EOF; ) {
    while (isspace(c)) c=fgetc(psp_file);  /* skip whitespace before number */
    ungetc(c,psp_file);
    if (c == '>') { 			/* end of "FASTA" sequence */
      break;
    }
    if (c==EOF)
      break;
    /* start of number; allow a '.' at start */		
    // want everything on the line this time
    number = numberstring (psp_file, TRUE, name, FALSE); 
    if (!number) {
      fprintf(stderr, "\nUnable to find valid number in PSP sequence %ld.\n", N+1);
      fprintf(stderr, "Fix PSP file.\n");
      return(-1); // bad data
    }
    
    if ((N % RCHUNK) == 0) {
      Resize(newprior, N+RCHUNK, double);
    }
    if (sscanf(number, "%lg", &nextnumber) != 1) {
      fprintf(stderr, "\nIllegal value for PSP `%s' in priors file `%s' at sequence %ld.\n", number, name, N+1);
      fprintf(stderr, "Fix PSP file.\n");
      if (number != NULL) myfree(number);
      return(-1);
    }
    newprior[N++] = nextnumber;
    myfree(number);
    length = 0;
  } /* read positional prior */
  *pospriors = newprior;
  myfree(number); // myfree checks if NULL and sets NULL after free

  return(N);
} /* read_psp_data */

/**********************************************************************/
/* 
 numberstring
 
 scan through the input up to next white space; assuming a correct initial
 character for a number, leave it up to sscanf to interpret the string
 
 Returns the next substring up to white space or EOF, NULL if error
*/
/**********************************************************************/
static char *numberstring(
  FILE *data_file,		/* data file */
  BOOLEAN point,		/* allow '.' at start? */
  char *name,			/* name of sequence */
  BOOLEAN skipToEOL             /* if TRUE skip rest of line */
)
{
  char *number,
    c;
  long length;
  number = NULL;
  length = 0;

  c=fgetc(data_file);
  if (!numprefix(c) || (!point && c == '.')) {		/* illegal character? */
    Resize(number, length+RCHUNK, char);
    number[length++] = c;
    return(number);
  }

  /* alocate space for ASCII version of number */
  while (!isspace(c) && c!= EOF) {	/* skip to whitespace or EOF */
    if ((length % RCHUNK) == 0) {
      Resize(number, length+RCHUNK, char);
    }
    number[length++] = c;
    c=fgetc(data_file);
  }

  /* put NULL at the end of the sequence */
  Resize(number, length+1, char);
  number[length] = '\0';
  if ((length % RCHUNK) == 0) {
    Resize(number, length+RCHUNK, char);
  }

  if (skipToEOL) {
    while (c!='\n' && c!='\r' && c!= EOF) {
      c=fgetc(data_file);
    }
  }

  return(number);
} /* numberstring */

/**********************************************************************/
/*
 all_zeros
 
 TRUE if the given range of the data is all zeros, false otherwise.
 
 */
/**********************************************************************/
static BOOLEAN all_zeros (
  BOOLEAN forward,
  double* data_end,
  int N
)
{
  int i, j;
  for (i = 0, j = 0; i < N; i++, j--)
    if (data_end[(forward? j : -j)] != 0.0) {
      int k, m;
      fprintf(stderr,"all_zeros: forward = %s j = %d, N = %d, data_end[(forward? j : -j)] = %lf\n",
	      (forward?"TRUE":"FALSE"), j, N, data_end[(forward? j : -j)]);
      for (k = 0, m = 0; k < N; k++, m--)
	fprintf(stderr,"data_end[%d] = %lf ",
		(forward? m : -m), data_end[(forward? m : -m)]);
      fprintf(stderr,"\n");
      return FALSE;
    }
  return TRUE;
} /* all_zeros */

