/* gblast.c - */
#include "gblast.h"

gb_typ MakeGBlast(long T, e_type E, a_type A)
{
	gb_typ B;
	long	n,q,q0,q1,q2,r0,r1,r2,s,t;
	long	c,d,e,i,x,y,z,hits,nlet;
	char	**best,*b0,*b1,*b2; /** letters ordered by decreasing score **/
/* 5-6-00 tlb; make R int */
	int 	**R;
	keytyp	key;
	dh_type	H;
	long        time1;

	NEW(B,1,gblast_type);
	B->A = A;
	B->E = E;
	/***** order letters by decreasing score *****/
	NEWP(best,nAlpha(A)+2,char);
	H = dheap(nAlpha(A)+3,3);
	for(i=0,c=0; c <= nAlpha(A); c++){
	    NEW(best[c],nAlpha(A)+2,char);
	    for(d=0; d <= nAlpha(A); d++){
		key = (keytyp) -valAlphaR(c,d,A);
		insrtHeap(d+1,key, H);
	    }
	    /*** fprintf(stderr,"%c: ", AlphaChar(c,A)); **/
	    for(i=0; (d = delminHeap(H)) != 0; i++){
		best[c][i] = d - 1;
	        /** fprintf(stderr,"%c ", AlphaChar(d-1,A)); ***/
	    }
	    /*** fprintf(stderr,"\n"); **/
	}
	Nildheap(H);
	/***** create finite automaton *****/
	/*** create transition function d(r,q). ***/
	n = nAlpha(A) +1; n = 1 + n + n*n;  
	NEWP(B->d,nAlpha(A)+2,long);
	for(c=0; c <= nAlpha(A); c++) NEW(B->d[c],n+2,long);
	for(q=0, c=0; c <= nAlpha(A); c++){
	    q++; B->d[c][0] = q; q0=q;		/** q0 = state "^c" **/
	    for(d=0; d <= nAlpha(A); d++){
	    	q++; B->d[d][q0] = q;        	/** q1 = state "cd" **/
	    }
	}
	B->nQ = q;				/** nQ = #states **/
	for(q=0, c=0; c <= nAlpha(A); c++){
	    q0 = B->d[c][0];			/** q = state "^c" **/
	    for(d=0; d <= nAlpha(A); d++){
	    	q1 = B->d[d][q0];		/** q = state "cd" **/
	    	q2 = B->d[d][0];		/** q2 = state "^d" **/
		for(e=0; e <= nAlpha(A); e++){
			q = B->d[e][q2];	/** q2 = state "de" **/
			B->d[e][q1] = q;	/** "cde" = "^de" = "de" **/
		}
	    }
	}
	/*** create acceptance states & positions pos[q][r][1..] ***/
	time1=time(NULL);
	nlet = nAlpha(A) +1;
	MEW(B->pos,nlet+2,ml_type);
	n = B->nQ;
	for(c=0; c <= nAlpha(A); c++) B->pos[c] = MkMList(30*LenSeq(E),n);
        /*** fprintf(stderr,"\nallocation time: %ld seconds\n",
		time(NULL)-time1); ****/
	time1=time(NULL);
	n = LenSeq(E) - 1;
	NEW(B->tmp,n+3,long); 
	c = XnuSeq(1,E);
	q = B->d[c][0];		/** q = state "^c" **/
	R = AlphaR(A);
	for(hits=0,i=2; i <= n; c=d, i++){
	   d = XnuSeq(i,E);
	   q = B->d[d][q];		/** q = state "xd" **/
	   e = XnuSeq(i+1,E);
/**** fprintf(stderr,"%c%c%c: (q=%d; d(%c,q)=%d)", 
		AlphaChar(c,A),AlphaChar(d,A),AlphaChar(e,A),
		q,AlphaChar(e,A),B->d[e][q]);
	   ** look through related 3-words for "matches". **/
	   b0 = best[c];
	   b1 = best[d];
	   b2 = best[e];
	   for(x=0; x <= nAlpha(A); x++){
	        r0 = b0[x];
		q0 = B->d[r0][0];
		s = T - (long) valAlphaR(c,r0,A);
		for(y=0; y <= nAlpha(A); y++){	
	           r1 = b1[y];
		   q1 = B->d[r1][q0];
		   t = s - (long) valAlphaR(d,r1,A);
		   for(z=0; z <= nAlpha(A); z++){	
			r2 = b2[z];
			if(t <= (long)valAlphaR(e,r2,A)){
				Add2MList(i-1, q1, B->pos[r2]);
				hits++;
/**** fprintf(stderr," %c%c%c(%d)", AlphaChar(r0,A), AlphaChar(r1,A),
					AlphaChar(r2,A),s+t+u);
				** add to acceptor states **/
			} else { break; } 
		   }
		   if(z==0) break;
		}
		if(y==0) break;
	   }
	   /*** fprintf(stderr,"\n");  ****/
	}
        /*** fprintf(stderr,"\nneighborhood time: %ld seconds\n",
		time(NULL)-time1); ****/
	for(c=0; c <= nAlpha(A); c++) free(best[c]); free(best);
	/** fprintf(stderr,"hits = %d\n",hits); 
        fprintf(stderr,"\tcreate time: %ld seconds\n", time(NULL)-time1);
	****/
	return B;
}

long	WinExtendGBlast(e_type E1, long i1, e_type E2, long i2, a_type A,
	long flank)
/** extend only flank residues in either direction **/
{
	long	n,n1,n2;
	register char	*p1,*p2;
	register long	k,score,max;

	n1 = LenSeq(E1); n2 = LenSeq(E2);
	p1 = XnuSeqPtr(E1); p2=XnuSeqPtr(E2);
	for(score=k=0, p1+=i1, p2+=i2; k<3; k++,p1++,p2++) {
		score += valAlphaR(*p1,*p2,A); 
	} 
	n = tMIN(long,n2-i2,n1-i1);
/** TEST **/
	n = tMIN(long,n,k+flank);
/** TEST **/
	for(max=score; k <= n; k++,p1++,p2++) {
		score += valAlphaR(*p1,*p2,A); 
		if(score > max) max = score;
		/** else if(score < (max-40)) break; ***/
		/**/ else if(score <= (max-40) || score < 0) break; /****/
	}
	n = tMIN(long,i1,i2);
/** TEST **/
	n = tMIN(long,n,flank);
/** TEST **/
	p1 = XnuSeqPtr(E1)+i1-1; p2=XnuSeqPtr(E2)+i2-1;
	for(score=max, k=1; k < n; k++,p1--,p2--){
		score += valAlphaR(*p1,*p2,A); 
		if(score > max) max = score;
		/** else if(score < (max-40)) return max; ***/
		/**/ else if(score < (max-40) || score < 0) return max; /****/
	}
	return max;
}

long	MatcherGBlastOffset(e_type E, gb_typ B, long *os)
{
	long	num,n,q,score,s,c,d,e,i,j,hits=0,max,*tmp=B->tmp,maxoff=0;
	a_type	A=B->A;
/** TEST **/
	long	flank = *os;
/** TEST **/

	n = LenSeq(E) - 1; c = XnuSeq(1,E);
	q = B->d[c][0];		/** q = state "^c" **/
	for(max=0, i=2; i <= n; c=d, i++){
	   d = XnuSeq(i,E);
	   q = B->d[d][q];		/** q = state "xd" **/
	   e = XnuSeq(i+1,E);		/** e = next token (mealy model) **/
	   if(!EmptyMList(q,B->pos[e])){ /** then q + e signals acceptance **/
		num=GetListMList(tmp, q, B->pos[e]);
		for(j=0;  j < num; j++){
			s = tmp[j];
/** TEST **/
			score=WinExtendGBlast(B->E, s, E, i-1, A, flank);
			if(score > max){ max = score; maxoff = s+1; }
/** TEST **/
/******
			score=ExtendGBlast(B->E, s, E, i-1, A);
			if(score > max){ max = score; maxoff = s-i+1; }
******/
		} hits+=j;
	   }
	}
	*os = maxoff;
	return max;
}

long	MatcherGBlast(FILE *fptr, e_type E, gb_typ B)
{
	long	num,n,q,score,s,c,d,e,i,j,hits=0,max,*tmp=B->tmp,maxoff=0;
	a_type	A=B->A;
	long        time1;

	time1=time(NULL);
	n = LenSeq(E) - 1; c = XnuSeq(1,E);
	q = B->d[c][0];		/** q = state "^c" **/
	for(max=0, i=2; i <= n; c=d, i++){
	   d = XnuSeq(i,E);
	   q = B->d[d][q];		/** q = state "xd" **/
	   e = XnuSeq(i+1,E);		/** e = next token (mealy model) **/
	   if(!EmptyMList(q,B->pos[e])){ /** then q + e signals acceptance **/
		num=GetListMList(tmp, q, B->pos[e]);
		for(j=0;  j < num; j++){
			s = tmp[j];
			score=ExtendGBlast(B->E, s, E, i-1, A);
			if(score > max){ max = score; maxoff = s-i+1; }
/*** DEBUG ***
			if(score == 47) 
			   fprintf(stderr,
				"%d - %c%c%c: s=%d\n",
				i-1,AlphaChar(c,A),
				AlphaChar(d,A),
				AlphaChar(e,A),
				s);
*** DEBUG ***/
		} hits+=j;
	   }
	}
	if(max > 0 && fptr != NULL) {
		fprintf(fptr,"hits = %ld; max score = %ld\n",hits,max); 
		fprintf(fptr,"score = %ld; offset = %ld\n", max,maxoff);
		PutDiagonalSeq(fptr, maxoff, B->E, E, A);
        	fprintf(fptr,"\tmatcher time: %ld seconds\n", 
			time(NULL)-time1);
	}
	return max;
}

long	ExtendGBlast(e_type E1, long i1, e_type E2, long i2, a_type A)
{
	long	n,n1,n2;
	register char	*p1,*p2;
	register long	k,score,max;

	n1 = LenSeq(E1); n2 = LenSeq(E2);
	p1 = XnuSeqPtr(E1); p2=XnuSeqPtr(E2);
	for(score=k=0, p1+=i1, p2+=i2; k<3; k++,p1++,p2++) {
		score += valAlphaR(*p1,*p2,A); 
	} 
	n = tMIN(long,n2-i2,n1-i1);
	for(max=score; k <= n; k++,p1++,p2++) {
		score += valAlphaR(*p1,*p2,A); 
		if(score > max) max = score;
		/** else if(score < (max-40)) break; ****/
		/**/ else if(score <= (max-40) || score < 0) break; /****/
	}
	n = tMIN(long,i1,i2);
	p1 = XnuSeqPtr(E1)+i1-1; p2=XnuSeqPtr(E2)+i2-1;
	for(score=max, k=1; k < n; k++,p1--,p2--){
		score += valAlphaR(*p1,*p2,A); 
		if(score > max) max = score;
		/** else if(score < (max-40)) return max; ****/
		/**/ else if(score < (max-40) || score < 0) return max; /****/
	}
	return max;
}

void	NilGBlast(gb_typ B)
{
	long	c;

	for(c=0; c <= nAlpha(B->A); c++) {
		NilMList(B->pos[c]);
		free(B->d[c]);
	}
	free(B->tmp); free(B->d); free(B->pos);
	free(B);
}

Boolean	FastMatcherGBlast(e_type E, gb_typ B, long score)
{
	long	num,n,q,s,c,d,e,i,j,*tmp=B->tmp;

	n = LenSeq(E) - 1; 
	c = XnuSeq(1,E);
	q = B->d[c][0];		/** q = state "^c" **/
	for(i=2; i <= n; c=d, i++){
	   d = XnuSeq(i,E);
	   q = B->d[d][q];		/** q = state "xd" **/
	   e = XnuSeq(i+1,E);		/** e = next token (mealy model) **/
	   if(!EmptyMList(q,B->pos[e])){ /** then q + e signals acceptance **/
		num=GetListMList(tmp, q, B->pos[e]);
		for(j=0;  j < num; j++){
		   s = tmp[j];
/***** if(score <= ExtendGBlast(B->E, s, E, i-1, B->A)) return TRUE; *****/
		   if(FastExtendGBlast(B->E, s, E, i-1, AlphaR(B->A),score))
				return TRUE; /*****/
		}
	   }
	}
	return FALSE;
}

/* 5-6-00 tlb; make R int */
Boolean	FastExtendGBlast(e_type E1, long i1, e_type E2, long i2,
	register int **R, long score)
/****************************************************************
 ****************************************************************/
{
	register long	s,max;
	register char	*p1,*p2,*end;

	p1 = XnuSeqPtr(E1)+i1; p2=XnuSeqPtr(E2)+i2;
	if((max=LenSeq(E2)-i2) > (s=LenSeq(E1)-i1)) end = p1 + s;
	else end = p1 + max;
	s = R[(int)p1[0]][(int)p2[0]] + R[(int)p1[1]][(int)p2[1]] + R[(int)p1[2]][(int)p2[2]];  
	for(max=s, p1+=3, p2+=3; p1 <= end; p1++,p2++) {
		if((s += R[(int)*p1][(int)*p2]) > max){
			if(s >= score) return TRUE;
			max = s;
		}
		/** else if(s <= (max-40)) break; ****/
		/**/ else if(s <= (max-40) || s < 0) break; /****/
	}
	p1 = XnuSeqPtr(E1)+i1; p2=XnuSeqPtr(E2)+i2;
	if(i1 <= i2) end = p1 - i1;
	else end = p1 - i2;
	for(s=max, p1--, p2--; p1 > end; p1--,p2--){
		if((s += R[(int)*p1][(int)*p2]) > max){
			if(s >= score) return TRUE;
			max = s;
		} 
		/** else if(s <= (max-40)) return FALSE; ****/
		/**/ else if(s <= (max-40) || s < 0) return FALSE; /****/
	}
	return FALSE;
}

/******************************* private *******************************/

long	gblast_error(char *s) { fprintf(stderr,"gblast: %s\n",s); exit(1); }
