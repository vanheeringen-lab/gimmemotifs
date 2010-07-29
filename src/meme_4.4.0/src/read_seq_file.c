/*
 * $Id: read_seq_file.c 4046 2009-09-28 00:47:29Z james_johnson $
 * 
 */

/***********************************************************************
*								       *
*	MEME++							       *
*	Copyright 1994-2008, The Regents of the University of California    *
*	Author: Timothy L. Bailey				       *
*								       *
***********************************************************************/

/* 9-30-99 tlb; add sequence to SAMPLE */
/* 6-29-99 tlb; add resic */
/* 9-22-08 pm; add posprior */
/*
	Supports FASTA and SWISSPROT formats.

Note:
	Assumes SWISSPROT format if "ID" appears at beginning of first line.
	Ignores duplicate (same <id>) sequences.

All formats:
	<id>		:= <alpha>+
			   The *unique* name of the sequence.  May *not*
			   contain whitespace!

			   NOTE: if <id> is "WEIGHTS", the <description>
			   field is parsed as a string of sequence weights 
			   and any <sequence> field (if any) is ignored.
			   All weights are assigned in order to the 
			   sequences in the file. If there are more sequences
			   than weights, the remainder are given weight one.
			   Weights must be greater than zero and less than
			   or equal to one.  Weights may be specified by
			   more than one "WEIGHT" entry which may appear
			   anywhere in the file.

	<description>	:= <alpha>+
			   A description of the sequence;
		           may contain whitespace.

	<sequence>	:= <alpha>+
			   The DNA or protein sequence data; 
                           may contain whitespace.

	<text>		:= <alpha>*
			   Text which is ignored.


Pearson FASTA format:
	INPUT DATA FILE FORMAT:
		[
		<header>
		<sequence>+
		]*

	<header>	:= ">"<id>" "<description>
			    Everything up to the first space is 
			    assumed to be the unique name of the sequence.
	
SWISSPROT Format
	[
	<ID>
	<DE>+
	<SQ>
	<sequence>+
	//
	]*

	<ID>		:= "ID   "<id>" "<text>
	<DE>		:= "DE   "<description>
	<SQ>		:= "SQ   "<text>
*/

#include "meme.h"

static SAMPLE *create_sample(
  char *alpha,		/* the alphabet */
  long length,		/* length of sequence */
  char *name,		/* name of sample */
  char *sequence,  	/* the sequence */
  char *descript,	///< Description of the sample
  BOOLEAN revcomp	// make space for reverse complement
);

/* local variables */
#define DATA_HASH_SIZE 100000
/* 
   The hash table will be appended to if read_seq_file is called more
   than once.  This is so that positive sequences the negative sequence file
   will be ignored.
*/ 
static HASH_TABLE ht_seq_names = NULL;	/* hash of dataset seq names */

/**********************************************************************/
/*
        get_sample_by_name

        Check if sample name is defined in ht_seq_names hash table
*/
/**********************************************************************/
extern SAMPLE *get_sample_by_name(
  char *sample_name
)
{
  HASH_TABLE_ENTRY * hash_entry = hash_lookup_str(sample_name, ht_seq_names);
  return(hash_entry != NULL ? (SAMPLE *) hash_get_entry_value(hash_entry) : NULL);
} // get_sample_by_name

/**********************************************************************/
/*
	read_seq_file

	Open a sequence file and read in the sequences. 
	Returns a dataset->

	setup_hash_alphabet must have been called prior to first call to this.
*/
/**********************************************************************/
extern DATASET *read_seq_file(
  char *file_name,		/* name of file to open */
  char *alpha,			/* alphabet used in sequences */
  BOOLEAN use_comp,		/* use complementary strands, too */
  double seqfrac 		/* fraction of input sequences to use */
)
{
  int i, j;
  FILE *data_file;		/* file with samples to read */
  FILE *prior_file=NULL;	/* file with positional priors to read */
  char *sample_name;		/* name of sample read */
  char *sample_de;		/* descriptor text for sample */
  char *sequence;		/* sequence read */
  long length;			/* length of sequence */
  BOOLEAN error=FALSE;		/* none yet */
  SAMPLE *sample;		/* sample created */
  DATASET *dataset;		/* dataset created */
  int n_samples=0;		/* number of samples read */
  double *seq_weights=NULL;	/* sequence weights */
  int n_wgts=0;			/* number of sequence weights given */

  /* create a hash table of sequence names */
  if (!ht_seq_names) ht_seq_names = hash_create(DATA_HASH_SIZE);

  /* create a dataset */
  dataset = (DATASET *) mymalloc(sizeof(DATASET));
  dataset->alength = strlen(alpha);
  dataset->alphabet = alpha;
  dataset->psp_w = 0;			// indicates no PSP was read
  dataset->log_psp_w = 0;		// so log_psp will get initialized

  /* open data file */
  if (file_name == NULL) {
    fprintf(stderr, "You must specify a data file or `stdin'.\n");
    exit(1);
  } else if (strcmp(file_name, "stdin")) {
    data_file = fopen(file_name, "r"); 
    if (data_file == NULL) {
      fprintf(stderr, "Cannot open file `%s'.\n", file_name);
      exit(1);
    }
  } else {
    data_file = stdin;
  }

  /* initialize maximum length of sequences */
  dataset->max_slength = 0;
  dataset->min_slength = 10000000;

  dataset->n_samples = 0;	/* no samples yet */
  dataset->samples = NULL;	/* no samples */

  while (read_sequence(data_file, &sample_name, &sample_de, &sequence, 
    &length)) {

    /* skip sequence if an error occurred */
    if (length < 0) continue;

    /* parse weights if given; make (more than enough) room in array */
    if (strcmp(sample_name, "WEIGHTS")==0) {
      double wgt; 
      char *wgt_str = sample_de;
      Resize(seq_weights, n_wgts+(int)strlen(wgt_str), double);
      while (sscanf(wgt_str, "%lf", &wgt) == 1) {
        if (wgt <= 0 || wgt > 1) {
	  fprintf(stderr, 
            "Weights must be larger than zero and no greater than 1.\n");
	  exit(1);
        }
        seq_weights[n_wgts++] = wgt;			/* save weight */
        wgt_str += strspn(wgt_str, "      ");		/* skip white */
        wgt_str += strcspn(wgt_str, "     ");		/* skip token */
      }
      myfree(sample_name);
      myfree(sample_de);
      myfree(sequence);
      continue;
    }

    /* ignore duplicate (same sample name) sequences */ 
    if (hash_lookup_str(sample_name, ht_seq_names) != NULL) {
      fprintf(stderr, "Skipping sequence '%s'.\n", sample_name);
      myfree(sample_name);
      myfree(sample_de);
      myfree(sequence);
      continue;
    }
    hash_insert_str(sample_name, ht_seq_names);  /* put name in hash table */

    n_samples++;

    /* see if sequence will be used in random sample; store it if yes */
    if (drand48() >= 1 - seqfrac) {

      HASH_TABLE_ENTRY *hash_entry; // needed to add pointer to sample

      /* create the sample */
      sample = create_sample(alpha, length, sample_name, sequence, sample_de, use_comp);
      if (sample == NULL) {error = TRUE; continue;}

      /* record maximum length of actual sequences */
      dataset->max_slength = MAX(sample->length, dataset->max_slength);
      dataset->min_slength = MIN(sample->length, dataset->min_slength);

      /* put the sample in the array of samples */
      if ((dataset->n_samples % RCHUNK) == 0) {
        Resize(dataset->samples, dataset->n_samples + RCHUNK, SAMPLE *);
      }
      dataset->samples[dataset->n_samples++] = sample;
      hash_entry = hash_lookup_str (sample_name, ht_seq_names);
      if (!hash_entry) {
	fprintf(stderr, "hash error: added sample ID then failed to find it in read_seq_file ()\n");
	error = TRUE;
	break;
      }
      hash_set_entry_value(sample, hash_entry);
    }
    
  } /* sequences */
  if (length < 0) error = TRUE;			/* read_sequence error */
  
  /* resize the array of samples */
  if (dataset->n_samples) Resize(dataset->samples, dataset->n_samples, SAMPLE*);

  /* check that datafile contained at least one sample */
  if (!error) {
    if (n_samples == 0) {
      fprintf(stderr, "No sequences found in file `%s'.  Check file format.\n",
	      file_name);
      error = TRUE;
    } else if (dataset->n_samples == 0) {
      fprintf(stderr, 
        "No sequences sampled.  Use different seed or higher seqfrac.\n");
      error = TRUE;
    }
  }

  /* exit if there was an error */
  if (error) exit(1);

  /* calculate the prior residue frequencies and entropies 
     and |D|, size of the dataset */
  /* tlb; 5/9/97 wgt_total_res and weighted res_freq */
  dataset->res_freq = NULL;
  Resize(dataset->res_freq, dataset->alength, double);
  for (i=0; i<dataset->alength; i++) { dataset->res_freq[i] = 0; }
  dataset->total_res = 0;
  dataset->wgt_total_res = 0;
  for (i=0; i<dataset->n_samples; i++) {		/* sequence */
    long slen = dataset->samples[i]->length;
    double sw = dataset->samples[i]->sw = (n_wgts > i ? seq_weights[i] : 1);
    dataset->total_res += slen;
    dataset->wgt_total_res += slen*sw;
    for (j=0; j<dataset->alength; j++) {
      if (use_comp) { 	/* using complementary strand as well */
        dataset->res_freq[j] += sw * dataset->samples[i]->counts[j]/2.0;
        dataset->res_freq[hash(comp_dna(unhash(j)))] += 
          sw * dataset->samples[i]->counts[j]/2.0;
      } else {		/* not using complementary strands */
        dataset->res_freq[j] += sw * dataset->samples[i]->counts[j];
      }
    }
  }

  /* convert counts to frequencies */
  for (i=0; i<dataset->alength; i++)  
    dataset->res_freq[i] /= dataset->wgt_total_res;

  return dataset;
} /* read_seq_file */

/**********************************************************************/
/*
	create_sample

	Create a sample.

	Returns the sample or NULL on error.

*/
/**********************************************************************/
static SAMPLE *create_sample(
  char *alpha,		/* the alphabet */
  long length,		/* length of sequence */
  char *name,		/* name of sample */
  char *sequence,       /* the sequence */
  char *descript,	///< Description of the sample
  BOOLEAN revcomp	// make space for reverse complement
)
{
  long i, j; 
  SAMPLE *new1;
  int alength = strlen(alpha);			/* length of alphabet */

  /* disallow zero length samples */
  if (length == 0) {
    fprintf(stderr, "\nZero length sequences not allowed. (Sequence `%s').\n",
      name);
    return NULL;
  }

  /* create the record to hold the sample and its associated data */
  new1 = (SAMPLE *) mymalloc(sizeof(SAMPLE));

  /* assign the name and sequence data */
  new1->sample_name = name;
  new1->seq = strdup(sequence);
  new1->descript = descript;

  /* set up encoded version of sequence and weights */
  new1->res = (char *) mymalloc(length * (int) sizeof(char));
  new1->resic = (char *) mymalloc(length * (int) sizeof(char));
  new1->pYic = (char *) mymalloc(length * (int) sizeof(char));
  new1->weights = (double *) mymalloc(length * (int) sizeof(double));
  new1->not_o = (double *) mymalloc(length * (int) sizeof(double));
  new1->log_not_o = (int *) mymalloc(length * (int) sizeof(int));
  new1->logcumback = (double *) mymalloc((length+1) * (int) sizeof(double));

  new1->nsites = 0;
  new1->sites = NULL;
  new1->psp_original = NULL;
  if (revcomp) {
    // Offset pointer -lseq so that z[j], j in [-lseq...+lseq]
    new1->z = length + (double *) mymalloc((2*length+1) * (int) sizeof(double));
    new1->log_psp = length + (double *) mymalloc((2*length+1) * (int) sizeof(double));
  } else {
    new1->z = (double *) mymalloc((length+1) * (int) sizeof(double));
    new1->log_psp = (double *) mymalloc((length+1) * (int) sizeof(double));
  }
  for (i=0; i<length; i++) { 
    int c = (int) sequence[i];
    int e = hash(c);
    new1->res[i] = e;
    // FIXME: make resic obsolete
    new1->resic[length-i-1] = hash(comp_dna(c));
    new1->weights[i] = 1.0;
  }

  /* set up arrays to hold posterior probabilities */
  create_2array(new1->pY, int, 3, length);

  /* set up array to hold character counts) */
  new1->counts = (double *) calloc(alength, (int) sizeof(double));

  /* get character counts */
  for (i=0; i<length; i++) {
    int c = new1->res[i];
    if (c<alength) {				/* normal letter */
      new1->counts[c]++;
    } else {					/* 'X' */
      for (j=0; j<alength; j++) new1->counts[j] += 1.0/alength;
    }
  }

  /* record length of sequence */
  new1->length = length;

  return new1;
} /* create_sample */
