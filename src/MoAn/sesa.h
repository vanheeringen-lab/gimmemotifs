
#ifndef SESA_FILE
#define SESA_FILE

#include "pssm.h" 

#define DEFAULT_MAX_SEQ 10
#define MAX_NAME_SIZE 1024

/**
 * The following struct represents the static enhanced suffix array.
 *
 * It is a huge waste of space and is not intended for large sequence
 * sets.
 */
typedef struct SESAStruct
{
  unsigned int  nSeq;              /* Number of sequencees */
  unsigned int  nSet;              /* Number of sets */
  unsigned int  *nSeqSet;          /* Number of sequences in each sets */
  unsigned int  maxSeq;            /* Current maximum number of sequencees */
  char **seqnames;                 /* Name of sequences */
  char *sequences;                 /* The actual sequences */
  unsigned int *pos;               /* Position in sequence */
  unsigned int *seq;               /* Sequence containing suffix */
  unsigned int *set;               /* Set containing sequence */
  unsigned int bg[4];              /* FIX: SHOULD BE PSSM */
  ESA esa;                         /* ESA */
} *SESA;


typedef struct HitTableSetStruct {
  unsigned int nTables;
  unsigned int table_size;
  unsigned int *pssmNr;
  struct HitTable **tables;
} *HitTableSet;

typedef struct  PSSMScoreTable {
  int nPssm;  
  int size;
  int *pssmNr;
} *PSSMScoreTable;

/* SESA EESA2SESA(EESA eesa); */

/* FUNCTIONS */
void free_SESA(SESA sesa);
struct HitTable *SESA_search(SESA sesa, PSSM pssm);
SESA build_SESA(int files, char **filenames, char *pAlphabet, char *pIgnore, int ds);
PSSMScoreTable SESA_exhaustive_search(SESA sesa, int length, int mismatches, int alphabetSize, int (*acceptanceFunc)(SESA, struct HitTable *));
PSSMScoreTable SESA_exclusive_search(SESA sesa, int length, int mismatches, int alphabetSize, int (*acceptanceFunc)(SESA, struct HitTable *));
/* int *get_pssmnr_list(char *seq, int seq_length, int motif_length); */

int chr2ind(char i);
char ind2chr(int i);


/* ACCESSORS */
#define getSeq(array, index)                    (array)->seq[index]
#define getPos(array, index)                    (array)->pos[index]

#endif
