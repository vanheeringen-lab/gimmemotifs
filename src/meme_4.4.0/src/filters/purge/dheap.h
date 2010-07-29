#if !defined(DHEAP)
#define DHEAP
/* Header file for d-heap data structure. Maintains a subset
 * of items in {1,...,m}, where item has a key.*/
#include "stdinc.h"

typedef double keytyp;

typedef	struct {
	long	N;			/* max number of items in heap */
	long	n;			/* number of items in heap */
	long	d;			/* base of heap */
	long	*h;			/* {h[1],...,h[n]} is set of items */
	long	*pos;			/* pos[i] gives position of i in h */
	keytyp	*kvec;			/* kvec[i] is key of item i */
} dheap_type;

typedef dheap_type *dh_type;

/************************ private operations ***************************/
/* parent of item, leftmost and rightmost children */
#define pHeap(x,H)          (((x)+((H)->d-2))/(H)->d)
#define leftHeap(x,H)       ((H)->d*((x)-1)+2)
#define rightHeap(x,H)      ((H)->d*(x)+1)
#define	MAX_KEY		    DBL_MAX

long	minchildHeap(long i,dh_type H);	/* returm smallest child of item */
void	siftupHeap(long i ,long x,dh_type H);
				/* move item up to restore heap order */
void	siftdownHeap(long i,long x,dh_type H);	
				/* move item down to restore heap order */
void    dheap_error(char *s);

/************************ public operations ***************************/
dh_type	dheap(long N,long D);
void	Nildheap(dh_type H);
void	insrtHeap(long i,keytyp k, dh_type H);	
					/* insert item with specified key */
long	rmHeap(long i,dh_type H);	/* remove item from heap */
long	delminHeap(dh_type H);		/* delete and return smallest item */
void	chkeyHeap(long i,keytyp k, dh_type H);		
					/* change the key of an item */
void	PutHeap(FILE *fptr, dh_type H);		/* print the heap */

/************************* macros definitions **************************/
#define	minHeap(H)	((H)->n==0?NULL:(H)->h[1])
#define	ItemsInHeap(H)	((H)->n)
#define minkeyHeap(H)	((H)->kvec[(H)->h[1]])
#define minItemHeap(H)	((H)->h[1])
#define keyHeap(i,H)	((H)->kvec[(i)])
#define memHeap(i,H)	((H)->pos[(i)]!=0?TRUE:FALSE)
#define	emptyHeap(H)	((H)->n==0)
#define	fullHeap(H)	((H)->n >= (H)->N)
#endif



