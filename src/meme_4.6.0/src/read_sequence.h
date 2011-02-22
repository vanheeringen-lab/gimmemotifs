#ifndef READ_SEQUENCE_H
#define READ_SEQUENCE_H

/* size of memory chunks to allocate for sequence data */
#define RCHUNK 100

/* maximum size of sequence description text string */
#define MAXDELEN 10000

/* local types */
typedef enum {FASTA, SWISSPROT} FORMAT_TYPE;

#include "macros.h"

extern BOOLEAN read_sequence(
  FILE *data_file,              /* file containing sequences */
  char **sample_name,           /* unique identifier of sequence */
  char **sample_de,             /* descriptor of sequence */
  char **sequence,              /* protein or DNA letters */
  long *length                   /* length of sequence */
);


#endif
