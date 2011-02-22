/***************** seqset.h - sequence set data type. ********************/
/*************************** Sequence Set ADT ******************************

	Define: P = <S,A>
	  where
		S is a set of sequences
		A is the sequence alphabet 

**************************************************************************/
#if !defined (SEQSET)
#define SEQSET
#include "stdinc.h"
#include "afnio.h"
#include "alphabet.h"
#include "karlin.h"
#include "sequence.h"
#include "random.h"
/******************************* PRIVATE ***********************************/
/****************************** seqset type ****************************/
typedef struct {
	char		*name;		/* input filename for entities */
	e_type		*entity;	/* array of sequence entities */
	long		nent;		/* number of input entities */
	long		max_leng;	/* sequence entity maximumlength */
	long		min_leng;	/* sequence entity minimum length */
	Boolean		xnu;		/* has data been xnu'ed? */
	long		*counts;	/* number of b residues */
	long		total;		/* total # residues */
	double		*tfreq;		/* residue total freqs for seqs */ 
	a_type		A;		/* sequence alphabet */
} seqset_type;
typedef seqset_type	*ss_type;

/******************************* private *************************************/
void	seqset_error(char *s);
ss_type  calcseqsetfreq(ss_type P);
ss_type  xnu_seqset(ss_type P);
long	count_seqset_entities(FILE *fptr, ss_type P,long nsize[]);
FILE    *OpenSeqSetFile(ss_type P);
ss_type  seqset(char *filename,a_type A);
ss_type  fptr_seqset(FILE *fptr,a_type A);
ss_type  SeqSet_fptr(FILE *fptr,a_type A);
ss_type  MkXnuSeqSet_fptr(FILE *fptr,a_type A);

/******************************* PUBLIC *************************************/
/************************** seqset operations ***************************/
ss_type  SeqSet(char *name,a_type A);
ss_type  SeqSet1(char *filename,e_type E,a_type A);
ss_type  MkXnuSeqSet(char *filename,a_type A);
ss_type  MkXnuSeqSet1(char *filename,e_type E,a_type A);
ss_type RmSeqSet(e_type E, ss_type P);
ss_type  NilSeqSet(ss_type P);
double  SeqSetEntropy(ss_type P);
long     *LengthsSeqSet(ss_type P);
ss_type	PutSeqSet(FILE *fptr,ss_type P);		/* show seqset data */
ss_type  PutSeqSetPIDs(FILE *fptr, ss_type P);
ss_type	PutSeqSetEs(FILE *fptr,ss_type P);
ss_type	PutSeqSetE(FILE *fptr,long i, ss_type P);
ss_type	PutSeqSettFreqs(FILE *fptr,ss_type P);
ss_type	ShuffleSeqSet(ss_type P); 		/* shuffle seqset */
ss_type  ShuffleSeqSet2(ss_type P);
double  LogL0SeqSet(ss_type P);
/************************** sequence set defines ***************************/
#define MAX_NUMBER_SEQS	15000
#define SeqSetA(P)		(P->A)
#define nLetSeqSet(P)	(nAlpha(P->A))
#define CountsSeqSet(r,P)	((P)->counts[(r)])
#define TotalSeqSet(P)	((P)->total)
#define	tFreqSeqSet(P)	((P)->tfreq)
#define SeqP(n,i,P)	ResSeq(i,(P)->entity[(n)])
#define NSeqsSeqSet(P)	((P)->nent)
#define SqLenSeqSet(n,P)	LenSeq((P)->entity[(n)])
#define	SeqSeqSet2(n,P)	SeqPtr((P)->entity[(n)])
#define	SeqSeqSet(n,P)	XnuSeqPtr((P)->entity[(n)])
#define SeqSetE(n,P)	((P)->entity[(n)])
#define NameSeqSet(P)	((P)->name)
#define MinSeqSeqSet(P)	((P)->min_leng)
#define MaxSeqSeqSet(P)	((P)->max_leng)
#define	CntsSeqSet(P)	((P)->counts)
#define	XnuSeqSet(P)	xnu_seqset(P)

#endif

