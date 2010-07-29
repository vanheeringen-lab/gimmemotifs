/* purge.c - homology purge program */
#include "purge.h"

Boolean	RmHomologs(long cutoff, char method, long minimum, Boolean query, 
	ss_type P)
/** Use block bit array for edges to save space **/
{
	long maxleng,i,k,entry,*lengs,*list;
	long item,item1,*hits,N;
	b_type	*edges;
	FILE	*fptr;
	dh_type	H;
	gb_typ	B = NULL;
	keytyp	key;
	e_type	E,E1,E2;
	char	str[100];
	Boolean	made_file=FALSE,related;

	N = NSeqsSeqSet(P);
	H = dheap(N+2,3);
	NEW(edges,N+2,b_type); 
	NEW(list,N+2,long); NEW(lengs,N+2,long); NEW(hits,N+2,long);
	for(maxleng=0,entry =1 ;entry <= N;entry++) {
		E=SeqSetE(entry,P);
		lengs[entry]=LenSeq(E);
		edges[entry] = Block(N+1);
		maxleng = tMAX(long,maxleng,LenSeq(E));
	}
	maxleng *= 2;
	entry = 1;
	if(cutoff < 99999) {
	  for(item =1 ;item <= N;item++) {
	    E1=SeqSetE(item,P);
	    if(method == 'b') B=MakeGBlast(11, E1, SeqSetA(P));
	    for(item1 =item+1 ;item1 <= N;item1++) {
	        E2=SeqSetE(item1,P);
		if(method == 'h'){	/** fasta heuristic method **/
		     related = RelateSeqFastp(E1,E2,SeqSetA(P),cutoff);
		} else if(method == 'b'){ /** gblast heuristic method **/
			related = FastMatcherGBlast(E2, B, cutoff);
/*****
			if(item==3 && item1==40) MatcherGBlast(stderr,E2, B);
*****/
		} else related = RelatedSeqs(cutoff, E1, E2, SeqSetA(P));
		if(related){
			AddBlock(item,edges[item1]);
			AddBlock(item1,edges[item]);
			hits[item]++; hits[item1]++;
		}
	    }
	    if(method == 'b') NilGBlast(B);
#ifdef VERBOSE
	    fprintf(stderr,"\r%ld",entry++); /**/
#endif
	  }
	}
#ifdef VERBOSE
	fprintf(stderr,"\n%ld items compared\n", entry-1); 
#endif
 if(!query){
	for(entry=1; entry <= N; entry++){
		key = -(keytyp) hits[entry];
		key += (keytyp)lengs[entry]/(keytyp)maxleng;
		/*if(hits[entry] > 0)
		     fprintf(stdout,"#%d: %d hits; key = %f\n",
			entry,hits[entry],key);*/
		insrtHeap(entry,key,H);
	}
	while(TRUE){
	    if(emptyHeap(H)) print_error("error in RmHomologs( )");
	    else if(minkeyHeap(H) >= 0) break;
	    else {
	 	item=delminHeap(H);
		for(i=1;i <= N;i++){
		  if(MemberBlock(i,edges[item]) && memHeap(i,H)){
			hits[i]--;
			key = - (keytyp) hits[i];
			key += (keytyp)lengs[i]/(keytyp)maxleng;
			chkeyHeap(i,key,H);
		  }
		}
	    }
	}
  } else {	/** don't purge #1 from set; i.e., don't put in heap. **/
	for(entry=2; entry <= N; entry++){
		key = -(keytyp) hits[entry];
		key += (keytyp)lengs[entry]/(keytyp)maxleng;
		/*if(hits[entry] > 0)
		     fprintf(stdout,"#%d: %d hits; key = %f\n",
			entry,hits[entry],key);*/
		insrtHeap(entry,key,H);
	}
	while(TRUE){
	    if(emptyHeap(H) || minkeyHeap(H) >= 0) break;
	    else {
	 	item=delminHeap(H);
		for(i=2;i <= N;i++){
		  if(MemberBlock(i,edges[item]) && memHeap(i,H)){
			hits[i]--;
			key = - (keytyp) hits[i];
			key += (keytyp)lengs[i]/(keytyp)maxleng;
			chkeyHeap(i,key,H);
		  }
		}
	    }
	}
	insrtHeap(1,-1e99,H);
  }
	if(ItemsInHeap(H) < N+1){		/* TLB; allow all seqs thru */
	 if(ItemsInHeap(H) >= minimum){
	  sprintf(str,"%s.%c%ld",NameSeqSet(P),method,cutoff);
	  fptr = open_file(str,"","w");
	  for(entry=1,k=0; entry <= N; entry++){
	    E=SeqSetE(entry,P); 
	    if(memHeap(entry,H)) { PutSeqSetE(fptr,entry,P); k++; }
	  }
	  fclose(fptr);
#ifdef VERBOSE
	  fprintf(stderr,"%ld sequences put into file: '%s'\n",k,str);
#endif
	  made_file = TRUE;
	 } else fprintf(stderr,
		"less than %ld sequences left (%ld); file not created\n",
		minimum,ItemsInHeap(H));
	} else fprintf(stderr,"no sequences removed; file not created\n");
	for(entry=1; entry <= N; entry++) {
		NilBlock(edges[entry]);
	}
	free(lengs); free(hits); free(list); free(edges);
	Nildheap(H);
	return made_file;
}

