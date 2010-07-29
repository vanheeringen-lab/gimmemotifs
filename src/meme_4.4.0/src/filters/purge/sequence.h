/* entity.h - entity data type. */
#if !defined(SEQUENCE)
#define SEQUENCE
#include "stdinc.h"
#include "alphabet.h"
#include "dheap.h"
#include "random.h"
/**************************** Sequence ADT ***************************
		E == <I,S>	2-tuple
		I == sequence identification key
		S == biological sequence 
**********************************************************************/
/*************************** Sequence type *****************************/
typedef struct {
	unsigned short		I;	/* identifier for entity */
	unsigned short		n;	/* length of entity sequence */
	Boolean			xnu;	/* sequences xnu'ed */
	char			*S;	/* sequence */
	char			*X;	/* if !xnu'ed == S; else X'ed seq */
	char			*info;	/* description of entity */
} sequence_type;
typedef sequence_type	*e_type;
/******************************** private **********************************/
long     get_diagonal_ends_seq(char *seq1, char *seq2, int **R, long n,
        long *begin, long *end);
e_type  EmptySeq(long I, long length);
long     seq_error(char *s);
/******************************** PUBLIC **********************************/
e_type  ReadSeq(FILE *fptr, long I, long size, a_type A);
e_type  ReadSeqFA(char *infile, long I, a_type A);
e_type  MergeSeqs(e_type *E);
Boolean IdentSeqs(e_type E1, e_type E2);
Boolean NonNullIdentSeqs(register e_type E1, register e_type E2);
long     CenterSeqHSP(long offset, e_type E1, e_type E2, a_type A);
long     GetFastaInfo(char *DBS_NAME, long max, long **pcounts,
        unsigned short **psize, a_type A);
e_type	NilSeq(e_type E);			/* undefine entity */
double  *FreqResSeq(e_type E, double *freq, a_type A);
long	*CountsSeq(e_type E, a_type A);
void    MaskSeq(long start, long end, e_type E);
long    *NumResSeq(e_type E,a_type A);
void    PutSeqID(FILE *fptr,e_type E);
void    PutSeqInfo(FILE *fptr,e_type E);
void    PutSeq(FILE *fptr,e_type E,a_type A);
void    PutXnuSeq(FILE *fptr,e_type E,a_type A);
void	ProcessSeq(e_type E,a_type A,double lambda,double K,double H);
e_type  CopySeq(e_type E);
e_type	ShuffleSeq(e_type E);
e_type  RandomSeq(long length, long I, double *freq, a_type A);
e_type	RtnShuffleSeq(e_type E);
e_type  RandomizeSeq(e_type E, double *freq, a_type A);
void    PutSeqRegion(FILE *fptr,long start, long length, e_type E, a_type A);
void    PutSeqRegion2(FILE *fptr,long start, long length, e_type E,
        long flank, a_type A);
long     PutDiagonalSeq(FILE *fptr, long offset, e_type E1, e_type E2, a_type A);
char    RandomResSeq(register e_type E);
/**************************** Macros ********************************/
#define	LenSeq(E)		((E)->n)
#define	ResSeq(i,E)		((E)->S[(i)])
#define	XnuSeq(i,E)		((E)->X[(i)])
#define EqSeq(i,x,E)		((E)->S[(i)]=(char) (x))
#define EqSeqKey(x,E)		((E)->info=(x))
#define RmSeq(E)		((E)->I=NULL)
#define EqSeqI(i,E)		((E)->I=(unsigned short) (i))
#define SeqI(E)			((E)->I)
#define SeqKey(E)		((E)->info)
#define SeqPtr(E)		((E)->S)
#define XnuSeqPtr(E)		((E)->X)
#endif
