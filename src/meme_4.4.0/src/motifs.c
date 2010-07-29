/*
 * $Id: motifs.c 4452 2010-04-14 03:45:52Z james_johnson $
 * 
 * $Log$
 * Revision 1.2.2.1  2006/01/13 02:34:43  twhitington
 * Replace hash.(c,h) with new and improved hash_table.(c,h).
 *
 * Revision 1.2  2005/10/25 19:06:39  nadya
 * rm old macro for Header, all info is taken care of by Id and Log.
 *
 * Revision 1.1.1.1  2005/07/29 17:23:18  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/***********************************************************************
*                                                                      *
*       MEME++                                                         *
*       Copyright 1994, The Regents of the University of California    *
*       Author: Timothy L. Bailey                                      *
*                                                                      *
***********************************************************************/

/**********************************************************************/
/*
        read_motifs

        Read in a motif file.
        Put positives for motifs into hash tables.
        If requested, the dataset 
        is read in and sequences are encoded  as integers using hash().
        Only, n_samples, sample_name, length, and res are valid in dataset.
        
        Two formats are supported:

                .motif file (known sites=positions in sequence)
                MOTIF <motif_name> <motif_width>
                [<seq#> <pos> <seq_name>]+

                .tag file (known family members=sequences)
                MOTIF <motif_name> 0
                [<seq#> <seq#> <seq_name>]+

                <seq#> = 0;     "?" in Prosite; ie, possible positive   
                <seq#> = 1;     "T/N" in Prosite; ie, tp or fn  

        Returns number of motifs read in.

*/
/**********************************************************************/

#include "meme.h"
#define MOTIF_HASH_SIZE 100000
#define DATA_HASH_SIZE 100000
#define RCHUNK 100

extern int read_motifs (
  FILE *fdata,                          /* opened dataset file */
  char *filename,                       /* motif file */
  MOTIF motif[NMOTIFS],                 /* motif info */
  BOOLEAN save_dataset,                 /* return dataset in memory */
  DATASET *dataset                      /* integer-encoded dataset */
)
{
  int i, nmotifs;
  FILE *fptr;
  char seq_name[MSN+1];
  HASH_TABLE ht_seq_names;              /* hash of dataset seq names */
  char *sample_name;                    /* name of sample from dataset */
  char *id;                             /* id of sample from dataset */
  char *sequence;                       /* sample from dataset */
  long length;                          /* length of sample from dataset */
  int seq_no = 0;                       /* number of sequences in dataset */

  /* create a hash table of sequence names */
  ht_seq_names = hash_create(DATA_HASH_SIZE);

  /* hash the names in the dataset 
     so that motif files can be checked for bad ones; column is set to 0 
  */
  while (read_sequence(fdata, &sample_name, &id, &sequence, &length)) {//FIXME this is incompatible with stdin! This looks like a bug
    /* skip sequence if there was an error */
    if (length < 0) continue;

    if (hash_lookup_str(sample_name, ht_seq_names) != NULL) {
      /* printf("Duplicate sequence: %s\n", sample_name); */
      myfree(sample_name);
      myfree(id);
      myfree(sequence);
      continue; 
    }
    hash_insert_str(sample_name, ht_seq_names); /* put name in hash table */
    if (save_dataset) {
      /* create a sample and put the sample in the array of samples */
      if ((seq_no % RCHUNK) == 0) {
              Resize(dataset->samples,
                  (seq_no + RCHUNK) * (int) sizeof(SAMPLE *), SAMPLE *);
      }
      dataset->samples[seq_no] = (SAMPLE *) mymalloc(sizeof(SAMPLE));
      dataset->samples[seq_no]->sample_name = sample_name;
      dataset->samples[seq_no]->length = length;
      /* psp fields all to be set before 1st use */
      dataset->samples[seq_no]->psp_original = NULL;
      dataset->samples[seq_no]->max_log_psp = LOG(0);
      dataset->samples[seq_no]->log_psp = NULL;

      /* encode the sequence as integers */
      for (i=0; i < length; ++i) {
        int c = (int) sequence[i];
        int e = hash(c);
        if (e == -1) {
          printf("\nIllegal character %c in sequence %s.  ", c, sample_name);
          printf("Change alphabet or fix data file.\n");
          return 0;
        }
        sequence[i] = e;
      }
      dataset->samples[seq_no]->res = sequence;
      dataset->n_samples = ++seq_no;
      /*printf("\r%s", sample_name);*/
    } else {
      myfree(sample_name);
      myfree(id);
      myfree(sequence);
    }
  }
  rewind(fdata);

  /* open motif file */
  fptr = fopen(filename, "r");
  if (fptr == NULL) {
    fprintf(stderr, "Cannot open file '%s'\n", filename);
    exit(1);
  }

  /* clear number of positives for each motif and best ROC */
  for (i=0; i<NMOTIFS; i++) {
    motif[i].pos = 0;           /* known motif positives */
    motif[i].roc = 0;           /* this motif not found yet */
    motif[i].recall = 0;        /* this motif not found yet */
  }

  /* read in the motifs and store them */
  nmotifs = 0;
  while ( (i = fscanf(fptr, " MOTIF %s %d", motif[nmotifs].name,
    &motif[nmotifs].width) ) == 2) {
    int seq, off;                       /* positon of positive occurrence */
    motif[nmotifs].ht = hash_create(MOTIF_HASH_SIZE);

    if (nmotifs == NMOTIFS) {
      fprintf(stderr, 
        "Too many known motifs (> NMOTIFS) in file %s!\n", filename);
      fprintf(stderr, "Recompile with larger value of NMOTIFS.\n");
      exit(1);
    }
    while ((i = fgetc(fptr)) != '\n'); /* ignore anything else on line */

    while ( (i = fscanf(fptr, "%d %d ", &seq, &off) ) == 2) {
      int c;
      /* get the name of the sequence from the motif file */
      for (i=0; (c=fgetc(fptr)) != EOF; ) {
        if (c=='\n') {ungetc(c, fptr); break;}
        else if (i < MSN) seq_name[i++] = c;
      }
      if (i == 0) printf("motif file must include sequence names!\n");
      seq_name[i] = '\0';                       /* null terminate name */
     
      /* check that name is valid */
      //if (!hash_lookup(seq_name, 0, ht_seq_names) != NULL) { 
      if (hash_lookup_str(seq_name, ht_seq_names) == NULL) { 
        if (seq != 0) {         /* don't worry about possibles */
          printf("sample `%s' in motif %d not in dataset\n", seq_name, nmotifs);
        }
      } else  {                                 /* store name and offset/tag */
        /* put in table only if found in dataset */
        hash_insert(seq_name, off, motif[nmotifs].ht);  /* site in hash table */
        /* if seq == 0, it is a tag file and this is NOT a real positive 
           ie, it is marked as "?" in Prosite 
        */
        if (seq != 0) motif[nmotifs].pos++;             /* number of pos's */
      }
      while ((i = fgetc(fptr)) != '\n'); /* ignore anything else on line */
    }

    nmotifs++;                          /* number of motifs */
    if (i == EOF) break;
  }

  if (i != EOF) {
    fprintf(stderr, "Error reading motif file %s\n", filename);
    exit(1);
  }

  return nmotifs;
}
