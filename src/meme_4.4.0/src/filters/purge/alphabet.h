/****************** alphabet.h - alphabet abstract data type.***************/
/* ident: Alpha/ field: A/ File: alphabet.h/ FileId: ALPHA/ type a_type
*************************** alphabet datastructure **********************


  char: X  A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
  code:	0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20

*************************************************************************/
#if !defined(ALPHA)
#define ALPHA
#include "stdinc.h"
#include <ctype.h>
/************************** alphabet datatype *****************************/
/* 5-6-00 tlb; make R int */
typedef struct {
	long     n;			/* number of LETTERS */
	char    *alphabet;		/* ALPHABET */
	char    *code2let;		/* CODE2LETTER */
	char    *code2lower;		/* CODE2LETTER lower case */
	char    *let2code;		/* LETTER2CODE */
/* 5-6-00 tlb; make R int */
	int     **R;			/* relatedness scoring matrix */
	char    *C;			/* complementary bases */
	long	loR;			/* lowest value in R */
	long	hiR;			/* highest value in R */
	char	*prs;			/* pairs string */
	char	**pairs;		/* residue pairs */
	Boolean	*paired;		/* is residue r paired? */
	long	npairs;			/* number of pairs */
} alphabet_type;
typedef	alphabet_type *a_type;

/******************************** PRIVATE **********************************/
void	alpha_error(char *s,a_type A);

/******************************** PUBLIC ***********************************/
/**************************** alphabet operations **************************/
a_type	MkAlpha(char *map_S,char *R);	/* define alphabet */
a_type  MkAlphabet(char *map_s,char *R,char *prs);
a_type  MkAlphabetC(char *map_s,char *R,char *prs,char *comp);
a_type	DefAlphaR(char *R,a_type A);	/* redefine R */
a_type  DefAlphaP(char *prs,a_type A);
a_type	PutAlpha(FILE *fptr,a_type  A); /* output alphabet A to file */
void	PutAlphaR(FILE *fptr,a_type  A);/* output alphabet Pairs to file */
void    PutAlphaPairs(FILE *fptr, a_type A);
void	NilAlpha(a_type A);		/* make alphabet A undefined */

/**************************** alphabet defines ****************************/
#define nAlpha(A)		((A)->n)
#define UndefAlpha(A)		0
#define AlphaR(A)		((A)->R)
#define valAlphaR(c,d,A)	((A)->R[(int)(c)][(int)(d)])
#define lowAlphaR(A)		((A)->loR)
#define highAlphaR(A)		((A)->hiR)
#define valAlphaP(c,d,A)	((A)->pairs[(c)][(d)])
#define	AlphaChar(i,A)		((A)->code2let[(int)(i)])
#define	AlphaCharLow(i,A)	((A)->code2lower[(int)(i)])
#define	AlphaCode(c,A)		((A)->let2code[(int)(c)])
#define	MapAlphaCode(A)		((A)->let2code)
#define	MapCodeAlpha(A)		((A)->code2let)
#define MemAlpha(c,A)		((long) (A)->let2code[(int)(c)])
#define NPairsAlpha(A)		((A)->npairs)
#define PairedAlpha(c,A)	((A)->paired[(int)(c)])
#define	AlphaComplement(i,A)	(((A)->C != NULL)?(A)->C[(int)(i)]:0)

/************ CONSTANTS *****************/
#define	ALPHA_NUM_SYMBOLS		127

#endif

