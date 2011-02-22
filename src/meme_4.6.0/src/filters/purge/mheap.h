#if !defined (MMHEAP)
#define	MMHEAP
#include "afnio.h"
#include "dheap.h"
#include "stdinc.h"
/*************************** ADT MMHEAP ***************************
	
	min-max heap

**********************************************************************/

/*************************** MMHEAP type **************************/
typedef struct {
	dh_type		heap;		/* minheap */
	dh_type		maxheap;	/* maxheap */
	long		hpsz;		/* heap size */
	long		nfree;		/* #items available */
	long		*avail;		/* list of available item no */
} mheap_type;

typedef mheap_type *mh_type;

/******************************* private *****************************/

/******************************* Public ******************************/
/***************************** operations ****************************/
mh_type Mheap(long hpsz, long d);
long	DelMinMheap(mh_type H);
long	DelMaxMheap(mh_type H);
long	RmMheap(long i, mh_type H);
long	InsertMheap(keytyp key, mh_type H);
mh_type	NilMheap(mh_type H);

/**************************** macro operations **********************/
#define ItemsInMheap(H)		ItemsInHeap((H)->heap)
#define MinKeyMheap(H)		minkeyHeap((H)->heap)
#define MinItemMheap(H)		minItemHeap((H)->heap)
#define MaxKeyMheap(H)		(-minkeyHeap((H)->maxheap))
#define EmptyMheap(H)		emptyHeap((H)->heap)
#define SizeMheap(H)		((H)->hpsz)
#define FullMheap(H)		fullHeap((H)->heap)
#endif

