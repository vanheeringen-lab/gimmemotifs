#if !defined(BLOCKS)
#define BLOCKS
#include <stdio.h>
#include <math.h>
#include "stdinc.h"

/*************************** ADT Block *******************************

    B = subset of { 0,...,n-1 } where n is an element of Postive integers.

**********************************************************************/
typedef struct {
	unsigned short	*list;		/* list of non-zero words */
	unsigned short	*list2;		/* 2nd list of non-zero words */
	unsigned char		*b;	/* as bytes = 8 bits */
	unsigned long	*i;		/* as long ints = 32 bits */
	long		N;		/* maximum number in set */
	long		nbyte;		/* number of bytes in array */
	long		nint;		/* number of ints in array */
} block_type;

typedef	block_type	*b_type;

/*************************** Private *******************************/
extern unsigned char *CARDINALITY_BLOCKS;
extern unsigned long *SET_BIT_BLOCKS;
extern unsigned short *B1_list,*B2_list,*IB_list;
extern unsigned char *B1_byte,*B2_byte,*IB_byte;

void    initialize_blocks();
void    add_list_block(b_type B);
void    update_list_block(b_type B);
void    blockerror(char *s);
/*************************** Public ********************************/
b_type	Block(long N);			/* create Null block */
b_type	BlockL(long N);			/* create Null block with list */
void    CopyBlock(b_type B1,b_type B2);
void    CopyBlockL(b_type B1,b_type B2);
void    CopyBlockL2(b_type B1, b_type B2);
void    CopyBlockLL(register b_type B1,register b_type B2);
void    ClearBlock(b_type B);			/* zero out a set in B */
void    FillBlock(b_type B);
void	AddBlock(long element,b_type B); 	/* add an element to set B */
void	DeleteBlock(long element,b_type B);
unsigned long    MemberBlock(register long element, register b_type B);
b_type  UnionBlock(b_type B1, b_type B2);
void    UnionBlock3(b_type B1, b_type B2, b_type UB);
long	UnionBlockCL(b_type B1, b_type B2);
void	UnionBlockL(b_type B1, b_type B2);
void    UnionBlockLF(b_type B1, b_type B2, b_type UB);
void    IntersectNotBlock(b_type B1, b_type B2);
void    IntersectBlock3(b_type B1, b_type B2);
void	IntersectBlock1(b_type B1, b_type B2,b_type IB);
long     IntersectBlockCF(b_type B1,b_type B2,b_type IB);
void	IntersectBlockF(b_type B1,b_type B2,b_type IB);
void	IntersectBlockL(b_type B1,b_type B2);
long     IntersectBlockList(long  *L, b_type B1, b_type B2);
long     IntersectBlockLCXOR(b_type eB,b_type B, b_type IB);
long	CardBlock(b_type B);	/* return cardinality of a set in B */
long     CardBlockL(b_type B);
long	CardInterBlock(b_type B1,b_type B2);
long	CardInterBlockLF(b_type B1,b_type B2);
void	CardInterBlockLMult(register long *C, b_type *B1,b_type B2);
long     *ListBlock(b_type B);	/* return list of members in B */
long     *ListBlockL(b_type B);	/* return list of members in B */
b_type	PutBlock(FILE *fptr,b_type B);	/* print block B */
b_type	NilBlock(b_type B);		/* destroy ISets */
/*************************** Macros **********************/
#define MAX_BLOCK_SIZE	524000		/* = list of 65500 */
#define BLOCK_BUFFER_SIZE	65505	/* for temp list */
#define END_BLOCK_LIST	65530
#define BlockN(B)	((B)->N)	/* return max cardinality of sets */
#define EmptyBlockL(B)	((B)->list[0]==END_BLOCK_LIST)
#define ClearBlockL(B)	((B)->list[0]=END_BLOCK_LIST)

#endif
