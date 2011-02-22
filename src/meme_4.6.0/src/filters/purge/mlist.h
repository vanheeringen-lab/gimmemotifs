/************ mlist.h - multiple linked list abstract data type.***********/
#if !defined(MLIST)
#define MLIST
#include "stdinc.h"
/************************** alphabet datatype *****************************/
typedef struct {
	long     M;			/* number of elements for lists */
	long     N;			/* number of lists */
	long	*first;			/* pointers to beginning of lists */
	long	*item;			/* items on lists */
	long	*next;			/* pointer to next element */
	long	free;			/* pinter to free list */
} multi_list_type;
typedef	multi_list_type *ml_type;

/******************************** PRIVATE **********************************/
void	mlist_error(char *s);

/******************************** PUBLIC ***********************************/
/**************************** mlist operations **************************/
ml_type	MkMList(long M, long N);
long     Add2MList(register long i, register long n, register ml_type L);
long     GetListMList(register long *list, register long n, 
        register ml_type L);
void	NilMList(ml_type L);


/***************************************************************************/
#define EmptyMList(n,L)		((L)->first[(n)]==0)

#endif

