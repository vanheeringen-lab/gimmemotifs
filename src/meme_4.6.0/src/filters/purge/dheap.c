#include "dheap.h"

dh_type	dheap(long N,long D)
/* Initialize a heap to store items in {1,...,N}.*/
{
	dh_type	H; long i;

	NEW(H,1,dheap_type);
	H->N = N; H->d = D; H->n = 0;
	NEW(H->h,N+1,long); NEW(H->pos,N+1,long); NEW(H->kvec,N+1,keytyp);
	for (i=1; i<= N; i++) H->pos[i] = 0;
	return H;
}

void	Nildheap(dh_type H)
{ free(H->h); free(H->pos); free(H->kvec); free(H); }

void	insrtHeap(long i,keytyp k, dh_type H)
/* insert item i with specified key */
{
	if(i<1) dheap_error("fatal! item to insert is < 1");
	if(H->pos[i]!=0) rmHeap(i,H);
	H->kvec[i] = k; H->n++; siftupHeap(i,H->n,H);
}

long	rmHeap(long i,dh_type H)
/* Remove item i from heap. */
{
	long j;
	if(H->pos[i]==0) return 0;
	j = H->h[H->n--];
	if (i != j && H->kvec[j] <= H->kvec[i]) siftupHeap(j,H->pos[i],H);
	else if (i != j && H->kvec[j]>H->kvec[i]) siftdownHeap(j,H->pos[i],H);
	H->pos[i] = 0;
	return i;
}

long	delminHeap(dh_type H)
/* delete and return item with smallest key */
{
	long i;
	if (H->n == 0) return 0;
	i = H->h[1];
	rmHeap(H->h[1],H);
	return i;
}


void	siftupHeap(long i ,long x,dh_type H)
/* Shift i up from position x to restore heap order.*/
{
	long px = pHeap(x,H);
	while (x > 1 && H->kvec[H->h[px]] > H->kvec[i]) {
		H->h[x] = H->h[px]; H->pos[H->h[x]] = x;
		x = px; px = pHeap(x,H);
	}
	H->h[x] = i; H->pos[i] = x;
}

void	siftdownHeap(long i,long x,dh_type H)
/* Shift i down from position x to restore heap order.*/
{
	long cx = minchildHeap(x,H);
	while (cx != 0 && H->kvec[H->h[cx]] < H->kvec[i]) {
		H->h[x] = H->h[cx]; H->pos[H->h[x]] = x;
		x = cx; cx = minchildHeap(x,H);
	}
	H->h[x] = i; H->pos[i] = x;
}

long	minchildHeap(long x,dh_type H)
/* Return the position of the child of the item at position x
   having minimum key. */
{
	long y, minc;
	if ((minc = leftHeap(x,H)) > H->n) return 0;
	for (y = minc + 1; y <= rightHeap(x,H) && y <= H->n; y++) {
		if (H->kvec[H->h[y]] < H->kvec[H->h[minc]]) minc = y;
	}
	return minc;
}

void	chkeyHeap(long i,keytyp k, dh_type H)
/* Change the key of i and restore heap order.*/
{
	keytyp ki;
	if(H->pos[i]==0) return;
	ki = H->kvec[i]; H->kvec[i] = k;
	     if (k < ki) siftupHeap(i,H->pos[i],H);
	else if (k > ki) siftdownHeap(i,H->pos[i],H);
}

void	PutHeap(FILE *fptr,dh_type H)
/* Print the contents of the heap. */
{
	long x;
	fprintf(fptr,"   h:");
	for (x = 1; x <= H->n; x++) fprintf(fptr," %2ld",H->h[x]);
	fprintf(fptr,"\nkvec:");
	for (x = 1; x <= H->n; x++) fprintf(fptr," %3f",H->kvec[H->h[x]]);
	fprintf(fptr,"\n pos:");
	for (x = 1; x <= H->n; x++) fprintf(fptr," %2ld",H->pos[H->h[x]]);
	fprintf(fptr,"\n");
}

void    dheap_error(char *s){fprintf(stderr,"dheap: %s\n",s);exit(1);}

