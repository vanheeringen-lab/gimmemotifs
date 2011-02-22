#include "block.h"
unsigned long *SET_BIT_BLOCKS=NULL;
unsigned char *CARDINALITY_BLOCKS=NULL;
unsigned char *BIT_NUMBER_BLOCKS=NULL;
unsigned short *B1_list,*B2_list,*IB_list;
unsigned char *B1_byte,*B2_byte,*IB_byte;

void	initialize_blocks(void)
/* initialize lookup tables used for bitwise operations */
{
	long i,b,n;
	unsigned char bit;

	NEW(SET_BIT_BLOCKS,32,unsigned long);
	for(i=0,b=31;i<32;b--,i++) { SET_BIT_BLOCKS[i] = 1 << b; }
	NEW(CARDINALITY_BLOCKS,256,unsigned char);
	for(i=0;i<256;i++){ 
		for(b=n=0;b<8;b++){ if(i & (1<<b)) n++; }
		CARDINALITY_BLOCKS[i] = n;
	}
	NEW(BIT_NUMBER_BLOCKS,9,unsigned char);
	for(bit=1,b=7;b>=0;bit=(bit<<1),b--) BIT_NUMBER_BLOCKS[b]=bit;
}

/***************** List operations for blocks *********************/
void	add_list_block(b_type B) /* add a list to block B */
{ MEW(B->list,(B->nbyte+1),unsigned short); B->list[0] = END_BLOCK_LIST; }

void	update_list_block(b_type B)
/* set all bytes in block B not on the list to zero */
{
	long	i,a;
	for(i=0,a=0;i<B->nbyte;i++) {
		if(i==B->list[a]) a++;
		else B->b[i]=0;
	}
}

/*************** Create and Copy operations for blocks ********************/
b_type  Block(long N) /* create and return the null Block B */
{
	b_type	B;

	if(N > MAX_BLOCK_SIZE) blockerror("Too many segments.");
	if(CARDINALITY_BLOCKS == NULL) initialize_blocks();
	MEW(B,1,block_type);
	B->list = NULL; B->N = N; 
	B->list2 = NULL;
	B->nint = (long) ceil((double)N/(double)32.0);
	B->nbyte= 4*B->nint;
	NEW(B->i,B->nint,unsigned long);
	B->b = (unsigned char*) B->i;
	return (B);
}

b_type  BlockL(long N) /* create and return a null block B with a list */
{
	b_type	B;

	B=Block(N);
	MEW(B->list,(B->nbyte+1),unsigned short);
	B->list[0] = END_BLOCK_LIST;
	return (B);
}

void	CopyBlock(b_type B1,b_type B2)
/* copy block B2 into block B1 - block B2 is modified */
{
	long j; 
	if(B1->N != B2->N) { fprintf(stderr,"N1=%ld N2=%ld\t",B1->N,B2->N); 
		blockerror("incompatible sets!"); }
	if(B2->list != NULL) update_list_block(B2);
	for(j=0;j<B1->nint;j++) B1->i[j] = B2->i[j];
}

void	CopyBlockL(b_type B1,b_type B2)
/* copy block B2 into block B1 maintaining a list - B1 is modified */
{
	long j,k; 
	if(B2->N != B1->N) { fprintf(stderr,"N1=%ld N2=%ld\t",B1->N,B2->N); 
		blockerror("incompatible sets!"); }
	if(B2->list != NULL) update_list_block(B2);
	if(B1->list == NULL) add_list_block(B1);
	for(j=0,k=0;j<B2->nbyte;j++)
		if((B1->b[j] = B2->b[j])) B1->list[k++]=j;
	B1->list[k]=END_BLOCK_LIST;
}

void	CopyBlockL2(b_type B1, b_type B2)
/* copy block B2 into block B1 maintaining a list - B1 is modified */
{
	register unsigned char *b1=B1->b,*b2=B2->b;
	register unsigned short *list1=B1->list,*list2=B2->list;
	register long i;

	for( ;(i= *list1= *list2) != END_BLOCK_LIST; list1++,list2++){
		b1[i] = b2[i];
	} 
}

void	CopyBlockLL(register b_type B1,register b_type B2)
{
	register long i,e;

	for(i=0;(e=B1->list[i]=B2->list[i]) != END_BLOCK_LIST; i++){ 
		B1->b[e] = B2->b[e];
	}
}

/**************** Modify operations for blocks ********************/
void	ClearBlock(b_type B)
/* set block B equal to the empty set */
{ long j; for(j=0;j<B->nint;j++) B->i[j] = (unsigned long) 0; }

void	FillBlock(b_type B)
/* set block B equal to the set of integers from 0 to N-1 */
{	long j; 
	for(j=0;j<B->nbyte;j++) B->b[j] = (unsigned char) 0xff; 
	for(j=B->N;j<8*B->nbyte;j++) DeleteBlock(j,B);
}

void	DeleteBlock(long element,b_type B)
/* delete an element from block B */
{
	long	nbit,nint;
	if(element >= B->nbyte*8) blockerror("input numbers too large");
	nbit = element & 0x1f;		/* element mod 32 */
	nint = element >> 5;		/* element / 32 */
	if(MemberBlock(element,B)) B->i[nint] ^= SET_BIT_BLOCKS[nbit];
}

void	AddBlock(long element,b_type B)
/* add an element to block B */
{
	long	nbit,nint;
	if(element >= B->N) blockerror("input numbers too large");
	nbit = element & 0x1f;		/* element mod 32 */
	nint = element >> 5;		/* element / 32 */
	B->i[nint] |= SET_BIT_BLOCKS[nbit];
}

unsigned long	MemberBlock(register long element, register b_type B)
/* return 0 if element is not a member of B; else return non-zero */
{
	register long	nbit,nint;

	if(element >= B->nbyte*8) blockerror("input numbers too large");
	nbit = element & 0x1f;
	nint = element >> 5;
	return (B->i[nint] & SET_BIT_BLOCKS[nbit]);
}

/**************** Union operations for blocks ********************/
void	UnionBlock3(b_type B1, b_type B2, b_type UB)
/* modifies B1 to equal B1 U B2; B1 is returned */
{
	long	j;
	for(j=0;j<B1->nint;j++) 
		UB->i[j] = B1->i[j] | B2->i[j];
}

b_type	UnionBlock(b_type B1, b_type B2)
/* modifies B1 to equal B1 U B2; B1 is returned */
{
	long	j;
	if(B1->N != B2->N) {
		fprintf(stderr,"N1=%ld, N2=%ld\n",B1->N,B2->N);
		blockerror("N1!=N2 for Union operation");
	}
	if(B1->list != NULL) update_list_block(B1);
	if(B2->list != NULL) update_list_block(B2);
	for(j=0;j<B1->nint;j++) 
		B1->i[j] = B1->i[j] | B2->i[j];
	return B1;
}

long	UnionBlockCL(b_type B1, b_type B2)
/* modifies B1 to equal B1 U B2; B1 is returned; B1 uses a list */
{
	long	b,k,n=0;
	unsigned short *list;

	if(B1->list2==NULL) MEW(B1->list2,(B1->nbyte+1),unsigned short);
	list = B1->list2;
	B1_byte = B1->b; B2_byte = B2->b;
	B2_list = B2->list; B1_list = B1->list;
	for(k=0; TRUE ; ){
	   if(*B1_list < *B2_list){ b= *B1_list; B1_list++;
	   } else if(*B1_list > *B2_list){
		b= *B2_list; B2_list++; 
	   	B1_byte[b] = B2_byte[b];
	   } else if(*B1_list==END_BLOCK_LIST){ break;	/* b1 == b2 */
	   } else {
	  	b= *B1_list; B1_list++; B2_list++;
	   	B1_byte[b] = B1_byte[b] | B2_byte[b];
	   }
           n+=CARDINALITY_BLOCKS[(B1_byte[b])];
	   list[k++]=b; 
	}
	list[k]= END_BLOCK_LIST;
	B1->list2 = B1->list;
	B1->list = list;
	return n;
}

void	UnionBlockLF(b_type B1, b_type B2, b_type UB)
/* modifies UB to equal B1 U B2; B1 is returned; all use a list */
{
	register unsigned short *list=UB->list;
	register unsigned short *list1=B1->list;
	register unsigned short *list2=B2->list;
	register long b;
	register unsigned char *b1=B1->b,*b2=B2->b,*ub=UB->b;

	while(TRUE){
	   if(*list1 < *list2){
		b= *list1; list1++;
	   	ub[b] = b1[b];
	   } else if(*list1 > *list2){
		b= *list2; list2++; 
	   	ub[b] = b2[b];
	   } else {		/* b1 == b2 */ 
		if(*list1==END_BLOCK_LIST) break;
	  	b= *list1; list1++; list2++;
	   	ub[b] = b1[b] | b2[b];
	   }
	   *list=b; 
	   list++;
	}
	*list= END_BLOCK_LIST;
}

void	UnionBlockL(b_type B1, b_type B2)
/* modifies B1 to equal B1 U B2; B1 is returned; B1 uses a list */
{
	long	b,k;
	unsigned short *list;

	if(B1->list2==NULL) MEW(B1->list2,(B1->nbyte+1),unsigned short);
	list = B1->list2;
	B1_byte = B1->b; B2_byte = B2->b;
	B2_list = B2->list; B1_list = B1->list;
	for(k=0; TRUE ; ){
	   if(*B1_list < *B2_list){ b= *B1_list; B1_list++;
	   } else if(*B1_list > *B2_list){
		b= *B2_list; B2_list++; 
	   	B1_byte[b] = B2_byte[b];
	   } else if(*B1_list==END_BLOCK_LIST){ break;	/* b1 == b2 */
	   } else {
	  	b= *B1_list; B1_list++; B2_list++;
	   	B1_byte[b] = B1_byte[b] | B2_byte[b];
	   }
	   list[k++]=b; 
	}
	list[k]= END_BLOCK_LIST;
	B1->list2 = B1->list;
	B1->list = list;
}

/*************** Intersect operations **********************/
void	IntersectNotBlock(b_type B1, b_type B2)
/* modifies B1 to equal B1 intersect not B2 */
{
	long	j;
	if(B1->N != B2->N) blockerror("N1!=N2 for Intersection operation");
	for(j=0;j<B1->nint;j++) B1->i[j] = B1->i[j] & ~(B2->i[j]);
}

void	IntersectBlock1(b_type B1, b_type B2,b_type IB)
/* modifies IB to equal B1 intersect B2 */
{
	long	j;
	if(B1->N != B2->N || B1->N != IB->N) 
		blockerror("N1!=N2 for Intersection operation");
	for(j=0;j<B1->nint;j++) {
		IB->i[j] = B1->i[j] & B2->i[j];
	}
}

long	IntersectBlockList(long	*L, b_type B1, b_type B2)
/* returns a list of the elements in the intersect of B1 and B2 */
{
	long	C,j,b,bit;
	unsigned char	s;

	for(C=j=0; (b=B2->list[j]) != END_BLOCK_LIST; j++) {
	   if((s=B1->b[b] & B2->b[b])) {
	     for(bit=0;bit<8;bit++)  {
		if(s & BIT_NUMBER_BLOCKS[bit]) {
			L[C++] = (b*8 + bit);
		}
	     }
	   }
	}
	L[C]= -1;
	return C;
}

void	IntersectBlock3(b_type B1, b_type B2)
/* modifies B1 to equal B1 intersect B2 */
{
	long	j;
	if(B1->N != B2->N) blockerror("N1!=N2 for Intersection operation");
	for(j=0;j<B1->nint;j++) B1->i[j] = B1->i[j] & B2->i[j];
}

void	IntersectBlockL(b_type B1, b_type B2)
/* modifies B2 to equal B1 intersect B2 */
{
	long i,k,e;

	if(B1->list != NULL) update_list_block(B1);
	for(i=k=0;(e=B2->list[i]) != END_BLOCK_LIST; i++){ 
	   if((B2->b[e] = B1->b[e]) & B2->b[e]) {B2->list[k++]=e; } 
	} B2->list[k]=END_BLOCK_LIST;
}

long	IntersectBlockCF(b_type B1, b_type B2,b_type IB)
/* modifies IB to equal B1 intersect B2; returns cardinality of IB */
{
	register unsigned char *b1=B1->b,*b2=B2->b,*ib=IB->b;
	register unsigned short *list=B2->list;
	register unsigned short *newlist=IB->list;
	register unsigned char *card = CARDINALITY_BLOCKS;
	register long i,k,a,n=0;

	/**if(B1->N != B2->N || B1->N != IB->N) {
		fprintf(stderr,"N1=%d, N2=%d; NIB=%d\n",B1->N,B2->N,IB->N);
		blockerror("N1!=N2 for Intersection operation");
	} DEBUG*/
	for(n=a=k=0; (i=list[a]) != END_BLOCK_LIST; a++){
	   if((ib[i]=b1[i]) & b2[i]) {
		newlist[k++]=i; n+=card[(ib[i])];
	   }
	}
	newlist[k]= END_BLOCK_LIST;
	return n;
}

void	IntersectBlockF(b_type B1, b_type B2,b_type IB)
/* modifies IB to equal B1 intersect B2 */
{
	register unsigned char *b1=B1->b,*b2=B2->b,*ib=IB->b;
	register unsigned short *list=B2->list;
	register unsigned short *newlist=IB->list;
	register long e;

	while((e= *list) != END_BLOCK_LIST){
		if((ib[e] = b1[e]) & b2[e]) {*newlist=e; newlist++; } 
 		list++; 
	} *newlist=END_BLOCK_LIST;
}

/*************** Cardinality operations **********************/
long     CardInterBlockLF(b_type B1,b_type B2)
/* return cardinality of the intersection of B1 and B2 - using lists */
{
	register unsigned char *b1=B1->b,*b2=B2->b;
	register unsigned short *list=B2->list;
	register unsigned char *card = CARDINALITY_BLOCKS;
	register long	i,n=0;

	while((i= *list) != END_BLOCK_LIST) {
		n+=card[(b1[i] & b2[i])]; list++;
	} 
	return n;
}

long     CardInterBlock(b_type B1,b_type B2)
/* return cardinality of the intersection of B1 and B2 */
{
	long	i,b,byte,n=0;
	if(B1->N != B2->N) {
		fprintf(stderr,"N1=%ld; N2=%ld\n",B1->N,B2->N);
		blockerror("N1!=N2 for CardIntersect operation");
	}
	for(i=0,byte=0; i < B1->nint; i++,byte+=4) {
	   if(B1->i[i] & B2->i[i])
	     for(b=0;b<4;b++)  {
		n+=CARDINALITY_BLOCKS[(B1->b[byte+b])&
				     (B2->b[byte+b])];
	     }
	}
	return n;
}

long     CardBlockL(b_type B)
{
	register unsigned char *b=B->b, *card = CARDINALITY_BLOCKS;
	register unsigned short *list=B->list;
	register long	i,n=0;

	while((i= *list) != END_BLOCK_LIST) {
		n+=card[b[i]]; list++;
	} 
	return n;
}

long     CardBlock(b_type B)
/* return cardinality of block B */
{
	long	i,b,byte,n=0;
	if(B->list != NULL) update_list_block(B);
	for(i=0,byte=0; i < B->nint; i++,byte+=4) {
	   if(B->i[i])
		for(b=0;b<4;b++) 
			n+=CARDINALITY_BLOCKS[(B->b[byte+b])];
	}
	return n;
}

/******************* List operations **********************/
long	*ListBlock(b_type B)
/* return list of elements in block B; list is  terminated by -1 */
{
	long	C,i,j=0,*L;
	
	if(B->list != NULL) update_list_block(B);
	C=CardBlock(B);
	NEW(L,C+1,long);
	for(i=0;i<B->N;i++) if(MemberBlock(i,B)) { L[j]=i; j++; }
	L[j]= -1;
	return L;
}

long	*ListBlockL(b_type B)
/* return list of elements in block B; list is  terminated by -1 */
{
	long	C,p,b,bit,j=0,*L;
	
	C=CardBlockL(B);
	NEW(L,C+1,long);
	for(p=0,b=B->list[p]; b != END_BLOCK_LIST; p++,b=B->list[p]){
		for(bit=0;bit<8;bit++) {
			if(B->b[b] & BIT_NUMBER_BLOCKS[bit]){
				L[j]=b*8+bit;j++;
			}
		}
	}
	L[j]= -1;
	return L;
}

/******************* Output operations **********************/
b_type	PutBlock(FILE *fptr,b_type B)
/* print block B to file pointed to by fptr */
{
	long	i,j=0,*m;
	
	if(B->list != NULL) update_list_block(B);
	m = ListBlock(B);/**/
	/* m = ListBlockL(B);*/
	fprintf(fptr," {");
	for(i=0;m[i] != -1; i++) {
		if(m[i+1] != -1) fprintf(fptr,"%5ld,",m[i]);
		else fprintf(fptr,"%5ld",m[i]);
		if(i%10==9 && m[i+1] != -1) fprintf(fptr,"\n  ");
	}
	if(j>0) fprintf(fptr,"%5ld }\n",m[i]);
	else fprintf(fptr," }\n");
	free(m);
	return B;
}

b_type  NilBlock(b_type B)
/* destroy block B */
{
	if(B != NULL) {
		if(B->list !=NULL) free(B->list);
		if(B->list2 !=NULL) free(B->list2);
		free(B->i); free(B);
	}
	return (b_type) NULL;
}

void	blockerror(char *s) { fprintf(stderr,"Block Error: %s\n",s); exit(1); }

/*************** operations **********************/
long	IntersectBlockLCXOR(b_type eB,b_type B, b_type IB)
/* Sets IIB = eB intersect B and sets EOB = B xor IB */
{
	register unsigned char *eb=eB->b,*b=B->b,*ib=IB->b;
	register unsigned short *list=B->list,*eolist=B->list;
	register unsigned short *newlist=IB->list;
	register unsigned char *card = CARDINALITY_BLOCKS;
	register long e,k,ik,n;

	for(n=ik=k=0; (e= *list) != END_BLOCK_LIST; list++){
	   if((ib[e]=eb[e]) & b[e]) {
		newlist[ik++]=e; n+=card[(ib[e])];
		if(b[e]^ib[e]){ eolist[k++] = e; }
	   } else { eolist[k++] = e; }
	}
	eolist[k] = newlist[ik] = END_BLOCK_LIST;
	return n;
}

