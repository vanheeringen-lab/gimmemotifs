/*---------------------------------------------------------------*
 * Wrapper for PSSM Searcher                                     *
 *                                                               *
 * By: Eivind Valen                                              *
 * University of Copenhagen, 2006                                *
 *---------------------------------------------------------------*/

#ifndef ESA_H_
#define ESA_H_

#include <stdio.h>
#include "shared.h"
#include "hittable.h"
#include "pssm.h"

/***************** DEFINITIONS *****************/
#define LOG_ZERO -10
#define MAX_NAME_SIZE 1024
#define DEFAULT_MAX_SEQ 10

/***************** DATA STRUCTURES *****************/

/**
 * EESAStruct
 *
 * seqnames   - array of the name of each sequence

 * seqborders - The normal library considers everything to be one
 *              sequence. This array contains the border so the sequence 
 *              and position within the sequence can be retrieved.
 *
 * The two arrays global2seq and global2pos is used if you need to
 * lookup position and sequence very frequently (and quickly) and are
 * not too concerned about space. This avoids the search through the
 * seqborder array, at the expense of one char and one int pr position
 * in the sequences!
 */
struct EESAStruct {
  unsigned int  nSeq;            /* Number of sequencees */
  unsigned int  maxSeq;          /* Current maximum number of sequencees */
  char **seqnames;               /* Name of sequences */
  unsigned int  *seqborders;     /* Borders between sequences */
  char *sequences;               /* The actual sequences */				  
  int seq_len;
  ESA esa;
};


/**
 * PSSMSetStruct
 *
 * Contains a set of PSSMs
 */
struct PSSMSetStruct {
  int   nmb;           /* Number of PSSMs */
  PSSM *pssms;         /* Array of PSSMs */
};


/**
 * ExtHitEntry
 *
 * Extended hit entry.
 * Contains the sequence and sequence position.
 */
typedef struct ExtHitEntry {  
  unsigned int seq;    /* Sequence */
  unsigned int pos;    /* Position in sequence */
  double score;        /* The attained score */
} Hit;


/**
 * HitsStruct
 *
 * Set of extended hit entries.
 */
struct HitsStruct {
  unsigned int nScores;        /* Number of hits registered */
  struct ExtHitEntry *pScores; /* Pointer to hits */
};


typedef struct EESAStruct *EESA;
typedef struct HitsStruct *Hits;


/***************** FUNCTIONS *****************/
EESA load_fasta(const char *pFilename, char *pAlphabet, char *pIgnore, int free_pStr);
Hits single_search(EESA eesa, PSSM pssm);
Hits convert_hittable(struct HitTable *hittable, EESA eesa);
int compare_hits(const void *pEntry1, const void *pEntry2);

/***************** MACROS *****************/
#define load_fasta_DNA(pFilename, free_pStr)  load_fasta(pFilename, "ACGT", "N", free_pStr)
#define load_fasta_protein(pFilename, free_pStr)  load_fasta(pFilename, "ABCDEFGHIKLMNPQRSTVWYZ", "", free_pStr)
#define get_seq_len(eesa, seq) (seq ? eesa->seqborders[seq] - eesa->seqborders[seq-1] : eesa->seqborders[seq])


#endif
