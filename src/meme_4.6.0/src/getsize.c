/*
 * $Id: getsize.c 4278 2009-12-23 09:58:37Z james_johnson $
 * 
 * $Log$
 *
 */

/***********************************************************************
*								       *
*	MEME
*	Copyright 1994--2009, The Regents of the University of California    *
*	Author: Timothy L. Bailey				       *
*								       *
***********************************************************************/
/* getsize.c */
/*	
	getsize <datafile> [options]

	Reads in a sequence file (Pearson/FASTA format).
	Prints the number of sequences, min L, max L, mean L, total residues
	and letters in alphabet used.

	Type "getsize" to see options.
*/
	
#define DEFINE_GLOBALS 
#include "read_sequence.h"
#include "general.h"
#include "hash_table.h"
#include "hash_alph.h"

#define DATA_HASH_SIZE 100000
#define LDNAB 30

/**********************************************************************/
/*
	main 
*/
/**********************************************************************/

extern int main(
  int argc,
  char *argv[]
)
{
  long i, j, k;
  char *datafile = NULL;
  FILE *data_file;
  int n_samples; 
  long max_slength, min_slength, length;
  char *sample_name, *id, *seq;
  int bfile = 0;			/* do not print bfile format */
  BOOLEAN l_only = FALSE;		/* do not print lengths only */
  BOOLEAN print_duplicates = TRUE;	/* print names of duplicates */
  BOOLEAN print_frequencies = FALSE;	/* do not print letter frequencies */
  BOOLEAN print_table = FALSE;		/* do not print frequencies table */
  BOOLEAN xlate_dna = FALSE;		/* do not translate DNA */
  BOOLEAN print_codons = FALSE;		/* do not print codon usage */
  HASH_TABLE ht_seq_names;		/* hash of dataset seq names */
  double codons[LDNAB][LDNAB][LDNAB];	/* counts of used codons */
  double alphabet[6][MAXASCII];		/* counts of used letters (per frame) */
  double total_res[6];			/* total letters per frame */
  char *letters = NULL;			/* string of used letters */

  /* set name of command */
  argv[0] = "getsize";

  /* get the command line arguments */
  i = 1;
  DO_STANDARD_COMMAND_LINE(1,
    USAGE(<datafile> [options]);
    NON_SWITCH(1, \n\t<datafile>\t\tfile containing sequences in FASTA format\n,
      switch (i++) {
        case 1: datafile = _OPTION_; break;
        default: COMMAND_LINE_ERROR;
      });
    FLAG_OPTN(1, l, \t\t\tjust print the length of each sequence, 
      l_only = TRUE);
    FLAG_OPTN(1, nd, \t\t\tdo not print warnings about duplicate sequences,
      print_duplicates = FALSE);
    FLAG_OPTN(1, dna, \t\t\tprint dna frequencies in bfile format, bfile = 1);
    FLAG_OPTN(1, prot, \t\t\tprint protein frequencies in bfile format, 
      bfile = 2);
    FLAG_OPTN(1, f, \t\t\tprint letter frequencies, print_frequencies = TRUE);
    FLAG_OPTN(1, ft, \t\t\tprint letter frequencies as latex table, 
      print_table = TRUE);
    FLAG_OPTN(1, x, \t\t\ttranslate DNA in six frames, xlate_dna = TRUE);
    FLAG_OPTN(1, codons, \t\tprint frame0 codon usage (implies -f xlate_dna), 
      print_codons = xlate_dna = print_frequencies = TRUE);
    USAGE(\n\tMeasure statistics of a FASTA file.);
  )

  /* 
    Setup hashing function for encoding strings as integers.
    If translating from DNA to protein, input alphabet must be DNAB;
    otherwise it may be all 26 letters plus asterisk.
  */
  if (xlate_dna) {
    setup_hash_alph(DNAB);			/* DNAB to position hashing */
    setup_hash_alph(PROTEINB);			/* PROTEINB to position hash */
    setup_hash_dnab2protb();			/* DNAB to PROTEINB hashing */
  } else if (bfile == 1) {			/* get DNA frequencies */
    setup_hash_alph(DNAB);			/* DNAB to position hashing */
  } else if (bfile == 2) {			/* get PROT frequencies */
    setup_hash_alph(PROTEINB);			/* PROTEINB to position hash */
  } else {
    setup_hash_alph("ABCDEFGHIJKLMNOPQRSTUVWXYZ*-");	/* full alphabet */
  }
  Resize(letters, MAXASCII, char);

  /* create a hash table of sequence names */
  ht_seq_names = hash_create(DATA_HASH_SIZE);

  /* open data file */
  if (datafile == NULL) {
    fprintf(stderr, "You must specify a data file or 'stdin'\n");
    exit(1);
  } else if ((strcmp(datafile, "stdin") == 0) || 
	     (strcmp(datafile, "-") == 0)) {
    data_file = stdin;
  } else {
    data_file = fopen(datafile, "r");
    if (data_file == NULL) {
      fprintf(stderr, "Cannot open file '%s'\n", datafile);
      exit(1);
    }
  }

  /* initialize counts of letters and codons used */
  for (i=0; i<6; i++) for (j=0; j<MAXASCII; j++) alphabet[i][j] = 0;
  for (i=0; i<LDNAB; i++) for (j=0; j<LDNAB; j++) for (k=0; k<LDNAB; k++)
    codons[i][j][k] = 0;

  /* initialize maximum length of sequences */
  max_slength = 0;
  min_slength = 10000000;

  n_samples = 0;			/* no samples yet */
  for (i=0; i<6; i++) total_res[i] = 0;	/* total residues (per frame) */

  while (read_sequence(data_file, &sample_name, &id, &seq, &length)) {
    /* Skip weights */
    if (strcmp(sample_name, "WEIGHTS")==0) continue;

    if (print_duplicates) { /* ignore duplicate (same sample name) sequences */
      if (hash_lookup_str(sample_name, ht_seq_names) != NULL) {
	fprintf(stderr, "Duplicate sequence: %s.\n", 
          sample_name);
	/* free up unneeded space */
	myfree(sample_name);
	myfree(id);
	myfree(seq);
	continue;
      }
      hash_insert_str(sample_name, ht_seq_names);	/* put name in table */
    } else {
      myfree(sample_name);
    }

    /*
       Count letters used in sequence.
    */
    if (!l_only) { 
      if (xlate_dna) {			/* translate DNA in six frames */
	for (i=0; i<2; i++) {		/* positive then negative strands */
	  if (i) invcomp_dna(seq, length);	/* negative strand */
	  for (j=0; j<3; j++) {		/* frame */
	    int f = j + 3 * i;
	    for (k=j; k<length-2; k+=3) {
	      /* hash DNAB codon to PROTEINB*/
	      int index = chash(TRUE, FALSE, seq+k);
	      alphabet[f][index]++;
	      total_res[f]++;
	      if (print_codons && f == 0) {		/* frame 0 */
		int i1 = dnabhash(seq[k]); 
		int i2 = dnabhash(seq[k+1]); 
		int i3 = dnabhash(seq[k+2]); 
		codons[i1][i2][i3]++;
	      }
	    }
	  } /* frame */
	} /* strand */
      } else {					/* don't translate */
	for (i=0; i<length; i++) alphabet[0][hash(seq[i])]++;
	total_res[0] += length;
      } /* xlate_dna */
    } /* !l_only */

    /* free up unneeded space */
    myfree(id);
    myfree(seq);

    /* record maximum length of actual sequences */
    max_slength = MAX(length, max_slength);
    min_slength = MIN(length, min_slength);

    if (l_only) printf("%ld\n", length);

    n_samples++;				/* number of non-duplicate sequences */
    if (print_frequencies && n_samples%1000 == 0) fprintf(stderr, "%d\r", n_samples);
  } /* read_sequence */

  /* print results */
  if (!l_only) {

    /* make string of letters used */
    for (i=0, j=0; i<MAXASCII; i++) 
      if (alphabet[0][i] || alphabet[1][i] || alphabet[2][i] || alphabet[3][i] 
        || alphabet[4][i] || alphabet[5][i]) letters[j++] = unhash(i);
    letters[j] = '\0';

    if (!bfile) {
      printf("%d %ld %ld %10.1f %.0f %s\n", n_samples, min_slength, max_slength,
        total_res[0]/n_samples, total_res[0], letters);
    }

    /*
      Print frequencies in bfile format
    */
    if (bfile) {
      char *alph_0 = bfile==1 ? DNA0 : PROTEIN0;	/* short alphabet */
      double tot = 0;					/* total legal let. */
      printf("# 0-order Markov frequencies from file %s\n", datafile);
      for (i=0; alph_0[i]; i++) tot += alphabet[0][hash(alph_0[i])];
      for (i=0; alph_0[i]; i++) {			/* letter */
        int k = hash(alph_0[i]);
        double freq = tot ? alphabet[0][k]/tot : 0;
	//printf("%c  %11.8f\n", alph_0[i], freq);
	printf("%c %.3e\n", alph_0[i], freq);
      } /* letters used */
    } /* bfile */

    /*
      Print letter frequencies as C arrays.
    */
    if (print_frequencies) {
      int nframes = xlate_dna ? 6 : 1;			/* number of frames */
      for (i=0; i<nframes; i++) {			/* frame */
        if (nframes) {
          printf("  double frame%ld[] = {\n", i);
        } else {
          printf("  double freq[] = {\n");
        }
	for (j=0; letters[j]; j++) {			/* letters used */
	  char c = letters[j];				/* letter */
	  int k = hash(c);				/* index */
	  printf("  %11.8f /* %c */,\n",alphabet[i][k]/total_res[i], c);
	} /* letters used */
        printf("  };\n");
      } /* frame */
      if (print_codons) {
        int alen = strlen(DNAB);			/* use DNAB alphabet */
        letters = DNAB;
        printf("  double fcodon[] = {\n");		/* start C array */
        for (i=0; i<alen; i++) {			/* pos 0 */
          for (j=0; j<alen; j++) { 			/* pos 1 */
            for (k=0; k<alen; k++) {			/* pos 2 */
              char c1 = letters[i]; char c2 = letters[j]; char c3 = letters[k]; 
              printf("  %11.8f /* %c%c%c */,\n", 
                codons[i][j][k]/total_res[0], c1, c2, c3);
            } /* pos 2 */
          } /* pos 1 */
        } /* pos 0 */
        printf("  };\n");
      } /* print_codons */
    } /* print_frequencies */

    /*
      Print letter frequencies as latex table.
    */
    if (print_table) {
      int nframes = xlate_dna ? 6 : 1;			/* number of frames */
      for (j=0; letters[j]; j++) {			/* letters used */
	char c = letters[j];				/* letter */
        int k = hash(c);				/* index */
        if (xlate_dna) {
          printf("%c & %7.4f", c, nrfreq[k]);
        } else {
          printf("%c", c);
        }
        for (i=0; i<nframes; i++) { 			/* frame */
	  printf(" & %7.4f", alphabet[i][k]/total_res[i]);
        } /* frame */
        printf(" \\\\\n");
      } /* letters used */
    } /* print_table */
  }

  return(0); 
} /* getsize */

