#include "mlist.h"

ml_type MkMList(long M, long N)
/******************************************************************
  N = number of lists; M = number of elements.
  sets all N lists to empty and creates the list of free cells.
 ******************************************************************/
{
	ml_type	L;
	long	i;

	MEW(L,1,multi_list_type);
	L->M = M; L->N = N;
	NEW(L->first,N+2,long);	/* all set to zero = empty lists */
	MEW(L->item,M+2,long); MEW(L->next,M+2,long);
	for(L->free=1, i=1; i < M; i++){
		L->next[i]=i+1;
	}
	L->next[i] = 0;
	return L;
}

long	Add2MList(register long i, register long n, register ml_type L)
/******************************************************************
 add item i to the nth list.
 ******************************************************************/
 {
	register long	new,t;

	if(n > L->N || n <= 0) mlist_error("list number out of range");
	if((new = L->free) == 0) mlist_error("out of memory");
	L->free = L->next[new]; 
	t = L->first[n];		/* first[n] = t[?:?]->  */
	L->item[new] = i;
	L->next[new] = t;		/** new[item:t]-> t[?:?]-> **/
	L->first[n] = new; 		/** first[n] = new[i:t]-> t[?:?]-> */
        return 0;
}

long	GetListMList(register long *list, register long n, 
	register ml_type L)
/******************************************************************
 convert list to contain the nth list in L; returns list length.
 WARNING: assumes list[ ] is long enough.
 ******************************************************************/
{
	register long	i,c;
	
	for(i=0, c = L->first[n] ; c != 0; c=L->next[c]){
		list[i++] = L->item[c];
	}
	return i;
}

void    NilMList(ml_type L)
{
	free(L->first);
	free(L->item);
	free(L->next);
	free(L);
}

void    mlist_error(char *s){ fprintf(stderr,"mlist: %s\n",s); exit(1); }
