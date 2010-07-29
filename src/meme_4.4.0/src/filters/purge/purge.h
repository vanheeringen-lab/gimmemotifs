/* purge.h - codes and constants for purge program. */
#if !defined (PURGE)
#define PURGE
#include <time.h>
#include "afnio.h"
#include "block.h"
#include "residues.h"
#include "seqset.h"
#include "pairaln.h"
#include "dheap.h"
#include "gblast.h"

Boolean	RmHomologs(long cutoff, char method, long minimum, Boolean query, 
	ss_type P);

        /********* N   A   C   T   G *******/
#define DNA_MTRX "-4  -4  -4  -4  -4 \
                  -4   5  -4  -4  -4 \
                  -4  -4   5  -4  -4 \
                  -4  -4  -4   5  -4 \
                  -4  -4  -4  -4   5 "

#define PURGE_USAGE	"\nusage: purge file score <options>\n\
  options:\n\
     [-n]    - sequences are DNA (default: protein)\n\
     [-b]    - use blast heuristic method (default for protein)\n\
     [-e]    - use an exhaustive method (default for DNA)\n\
     [-q]    - keep first sequence in the set\n\
     [-x]    - use xnu to mask protein tandem repeats\n\
\n\
  Purge creates an output file from the input file such that\n\
  no two sequences have local alignment score greater than <score>.\n\
  The output file is named <file>.<score>.\n\
  Substitution matrices: BLOSUM62 (protein), +5/-1 (DNA).\n\
\n"

#endif

