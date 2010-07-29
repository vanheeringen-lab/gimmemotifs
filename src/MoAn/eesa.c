/* SESA is a wrapper for the ESA implementation */
/* Copyright (C) 2006 Bioinformatics Centre */
/* Author: Eivind Valen */

/* This program is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU General Public License */
/* as published by the Free Software Foundation; either version 2 */
/* of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program; if not, write to the Free Software */
/* Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA. */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#include "build.h"
#include "search.h"
#include "hittable.h"
#include "error.h"
#include "pssm.h"
#include "search.h"
#include "eesa.h"



/* INITIALIZATION */

inline EESA init_eesa() {
  EESA eesa = (EESA) malloc(sizeof(struct EESAStruct));
  eesa->seqnames = (char **) malloc(sizeof(char *) * DEFAULT_MAX_SEQ);
  eesa->seqborders = (unsigned int *) malloc(sizeof(unsigned int) * DEFAULT_MAX_SEQ);
  eesa->maxSeq = DEFAULT_MAX_SEQ;
  eesa->nSeq = 0;

  return eesa;
}

inline Hits init_hits(int nScores) {
  Hits hits = (Hits) malloc(sizeof(struct HitsStruct));
  hits->nScores = nScores;
  hits->pScores = (struct ExtHitEntry *) malloc(sizeof(struct ExtHitEntry) * nScores);
  return hits;
}

/* LOADING */

inline void add_sequence(EESA eesa, char *seq, unsigned int size, unsigned int pos) {
  unsigned int nSeq = eesa->nSeq;

  if (nSeq == eesa->maxSeq) {
    eesa->maxSeq *= 2;
    eesa->seqnames = (char **) realloc(eesa->seqnames, sizeof(char *) * eesa->maxSeq);
    eesa->seqborders = (unsigned int *) realloc(eesa->seqnames, sizeof(unsigned int) * eesa->maxSeq);
  }

  eesa->seqborders[nSeq] = pos;
  eesa->seqnames[nSeq] = (char *) malloc(sizeof(char) * size);
  strcpy(eesa->seqnames[nSeq], seq);
  eesa->nSeq++;
}

EESA load_fasta(const char *pFilename, char *pAlphabet, char *pIgnore, int free_pStr) {
  char *sequences;
  FILE *file;
  int size;
  char c;
  char namebuffer[MAX_NAME_SIZE+1];
  EESA eesa = init_eesa(); 

  if(!pFilename) {
    setError("Need to supply a filename");
    return NULL;
  }

  /* Open sequencefile */
  if(!(file = fopen(pFilename, "r"))) {
    char st[512];
    sprintf(st, "Could not open '%s' for reading.", pFilename);
    setError(st);
    return NULL;
  }

  /* Check file size. Will overestimate by size of labels, but who cares. */
  fseek(file, 0, SEEK_END);
  sequences = malloc(sizeof(char) * (ftell(file) + 1));
  rewind(file);

  /* File to string */

/*   fread(sequences, 1, size, file);  */

  size = 0;
  while ((c = fgetc(file)) != EOF) {     
    if (c == '\n') continue;

    /* Fasta header */
    if (c == '>') {
      if (eesa->nSeq) sequences[size++] = '|';

      /* Lasts till end of line (or til MAX_NAME_SIZE) */
      int n = 0;
      while ((c = fgetc(file)) != '\n') { 
	if (n < MAX_NAME_SIZE) 
	  namebuffer[n++] = c;
      }
      namebuffer[n++] = '\0';
      add_sequence(eesa, namebuffer, n, size);
      continue;
    }

    sequences[size++] = toupper(c);
  } 

  sequences[size] = '\0';
  fclose(file);

  char *newIgnore = malloc(sizeof(char) * strlen(pIgnore) + 2);
  strcpy(newIgnore, pIgnore);
  strcat(newIgnore, "|");

  eesa->sequences = sequences;

/*   fprintf(stderr, "FRE: %i\n", free_pStr);  */
/*   fprintf(stderr, "SIZ: %i\n", strlen(sequences));  */
/*   fprintf(stderr, "SIZ: %i\n", size);  */
/*   fprintf(stderr, "ALP: %s\n", pAlphabet);  */
/*   fprintf(stderr, "IGN: %s\n", newIgnore);  */
/*   fprintf(stderr, "SEQ: (%s)\n", sequences); */

  /* Build the ESA */
  eesa->esa = build_ESA(sequences, size, pAlphabet, newIgnore, free_pStr); 

/*   int pos; */
/*   for (pos = 0; pos < size; pos++) */
/*     fprintf(stderr, "%c", eesa->esa->pStr[pos]); */


  return (eesa->esa ? eesa : NULL);
}





/* SEARCHING */

Hits single_search(EESA eesa, PSSM pssm) {
  Hits results;

  struct HitTable *hittable = search(eesa->esa, pssm);
  
  results = convert_hittable(hittable, eesa);
  release_hits(hittable);
  
  return results;
}


int compare_hits(const void *pEntry1, const void *pEntry2) {
  return ((struct HitEntry *)pEntry1)->position - ((struct HitEntry *)pEntry2)->position;
}


Hits convert_hittable(struct HitTable *hittable, EESA eesa) {
  Hits hits = init_hits(hittable->nScores);
  unsigned int i, seq;

  qsort((void*)hittable->pScores, hittable->nScores, sizeof(struct HitEntry), compare_hits);

  i = hittable->nScores;
  seq = eesa->nSeq - 1;
  while (i--) {
    struct HitEntry old = hittable->pScores[i];

    while (old.position < eesa->seqborders[seq]) {
      seq--;
    }

    hits->pScores[i].seq = seq;
    hits->pScores[i].pos = old.position - eesa->seqborders[seq];
    hits->pScores[i].score = old.score;
  }


  return hits;
}


/* int count_hits(struct HitTable *hittable, EESA eesa) { */
/*   qsort((void*)hittable->pScores, hittable->nScores, sizeof(struct HitEntry), compare_hits); */
  

/* } */



