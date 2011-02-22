#include "pairaln.h"

/****************************************************************
 first:  (seq2=1)
    from s1 = 1
               v
               DSIMVNPPY             s1 = s2 = 1;
               IATNPPYVGHIKG    	   s1++;
               ^
   to s1 = n1-3;
               v
         DSIMVNPPY                       s1 = n1 - 3 = 9 - 3 = 6
               IATNPPYVGHIKG             (s1...s1+end && s2...end) 
               ^				s2 = 1;
 second:  
              v                                   v
              IATNPPYVGHIKG    		IATNPPYVGHIKG
               DSIMVNPPY        to                DSIMVNPPY             
               ^                                  ^			
******************************************************************/
/* 5-6-00 tlb; make R int */
Boolean	relate_seq(register char *seq1, register char *seq2,
	register int **R, register long n, register long cutoff)
{
	register long min=0,sum=0;

	while(n > 0){
        	if(min < (sum += R[(int)seq1[n]][(int)seq2[n]])){
                       if(sum-min >= cutoff) return TRUE;
        	} else min = sum;
		n--;
        };
	return FALSE;
}

Boolean	RelatedSeqs(long cutoff, e_type E1, e_type E2, a_type A)
/*********************************************************************
Compare all diagonals of E1 and E2 to see if the two have an MSP
with score > cutoff.
 *********************************************************************/
{
	long	i,s,n1,n2,end,w;
	char	*seq1,*seq2;

	w = cutoff/highAlphaR(A);
	n1 = LenSeq(E1); seq1 = XnuSeqPtr(E1); 
	n2 = LenSeq(E2); seq2 = XnuSeqPtr(E2);
        for(end = n1 - w, s=0; s < end; s++,i--) {
	   	i = tMIN(long,n1-s,n2);
		if(relate_seq(seq1+s, seq2, AlphaR(A),i,cutoff)) return TRUE;
        }
        for(end = n2 - w, s=1; s <= end; s++) {
		i = tMIN(long,n2-s,n1);
		if(relate_seq(seq2+s, seq1, AlphaR(A),i,cutoff)) return TRUE;
	}
	return FALSE;
}

/* 5-6-00 tlb; make R int */
long	repeat_score_seq(register char *seq, register int **R, 
		register long o, register long n)
/********************************************************************
             /ptr=seq1+s
seq1: 	MQNKSQKETGDILGISQMHVSRL
seq2:	     MPPLFVMNNEILMHLRALKKTKKDVS
	     |...... n .......|
 ********************************************************************/
{
	register long min=0,sum=0,score=-9999;

	for(seq++; n > o; seq++){
        	if(min > (sum += R[(int)seq[0]][(int)seq[o]])){
			min = sum;
		} else if(score < (sum-min)) score = (sum-min);
		n--;
        };
	return score;
}

long	RepeatScoreSeq(e_type E, a_type A, long *offset)
/** look for internal repeats **/
{
	long	i,best,s,n,score;
	char	*seq;

	n = LenSeq(E); seq = XnuSeqPtr(E); 
        for(score= -9999,best=i=1; i < n; i++) {
		s = repeat_score_seq(seq,AlphaR(A),i,n);
		if(s > score) { best = i; score = s; }
        }
	*offset = best;
	return score;
}

/* 5-6-00 tlb; make R int */
long	diagonal_score_seq(register char *seq1, register char *seq2,
	register int **R, register long n)
/********************************************************************
             /ptr=seq1+s
seq1: 	MQNKSQKETGDILGISQMHVSRL
seq2:	     MPPLFVMNNEILMHLRALKKTKKDVS
	     |...... n .......|
 ********************************************************************/
{
	register long min=0,sum=0,score=-9999;

	while(n > 0){
        	if(min > (sum += R[(int)seq1[n]][(int)seq2[n]])){
			min = sum;
		} else if(score < (sum-min)) score = (sum-min);
		n--;
        };
	return score;
}

long	AlignSeqFastp(e_type E1, e_type E2, a_type A)
/*******************************************************************
 Align two sequences using the fastp algorithm.
 see W.R. Pearson & D.J. Lipman PNAS (1988) 85:2444-2448 and 
     Wilbur WJ; Lipman DJ (1983) Rapid similarity searches of nucleic acid 
     and protein data banks.  Proc Natl Acad Sci U S A 80: 726-730.
 *******************************************************************/
{
	long	i,v,x,s,*D,**Q,*nQ,nsave=10;
	char	*seq1,*seq2;
	long	n,r,n1,n2,score=0,item;
	mh_type	H;
	Boolean	hits=FALSE;
	long	*off,*off0;
	/** may want to remove spacer stuff - doesn't seem to help **/
	keytyp	*key,*key0;  /** freq. of 1to5-spaces between matches **/
	long	*lastpos,*lastpos0,maxspace=10;

	n1 = LenSeq(E1); n2 = LenSeq(E2);
	seq1 = XnuSeqPtr(E1); seq2 = XnuSeqPtr(E2);
	/** Construct a lookup table for the query sequence (n1) using k=1 
	    words with positions of all words in seq1. **/
	NEW(nQ,nAlpha(A)+3,long);
	MEW(Q,nAlpha(A)+3,long*);
	for(i=0; i<= nAlpha(A); i++) NEW(Q[i],n1+3,long);
	for(s=1; s <= n1; s++) {
		r = seq1[s];
		Q[r][nQ[r]] = s; nQ[r]++;
	}
        /** For each residue in a database sequence look up positions for 
	    all identical residues in the table & calculate the difference 
	    in position for all matches.  (Store all matches with the same 
            offset in a second table). ***/
	v = n1 + n2; 
	NEW(off0,v+4,long); off = off0 + n2 + 2;
	NEW(lastpos0,v+4,long); lastpos = lastpos0 + n2 + 2;
	NEW(key0,v+4,keytyp); key = key0 + n2 + 2;
	for(s=1; s <= n2; s++) {
		r = seq2[s];
		for(i=0; i< nQ[r]; i++){
			x = Q[r][i] - s; /** find offset **/
			off[x]++;
			if(lastpos[x]!=0){
				n = s - lastpos[x]; /** get spacer **/
				if(n <= maxspace) key[x]+=4.0/(keytyp)n;
			} else lastpos[x] = s;
		}
	}
        /** locate regions of similarity using offsets with high # matches **/
	H = Mheap(nsave,3);
	NEW(D,nsave+1,long);
	for(i = 1-n2; i < n1; i++){
		if(i < 0) v = n2 + i; /** get length of diagonal **/
		else v = n1 - i;
		key[i] += (keytyp)off[i] - (keytyp)v/(keytyp)nAlpha(A); /**/
		if(key[i] > 3.0) {
		   item=InsertMheap(-key[i],H);
		   D[item] = i;
		   /**** printf("key[%d]=%g (v=%d; percent = %g)\n",
			i,key[i],v, 100.0*(double)key[i]/(double)v); **/
		   /**** printf("off[%d]=%d (O-E=%f; v=%d; ave=%g)\n",
                        i,off[i],(keytyp)v/(keytyp)nAlpha(A),
			v, (double)off[i]/(double)v); **/
		}
	}
        /** compute the score for the ten offsets of highest similarity **/
	/*** printf("\nn1 = %d; n2 = %d\n",n1,n2);**/
	while((item=DelMinMheap(H))!=0){
		i = D[item];
		/*** fprintf(stderr,"key[%d] = %g; s[i] = ",i,key[i]); **/
		if(i > 0){	/** -> start at seq1[i] **/
		   n = tMIN(long,n1-i,n2);
		   s = diagonal_score_seq(seq1+i, seq2, AlphaR(A),n);
		} else {	/** -> start at seq2[i] **/
		   n = tMIN(long,n2+i,n1);
		   s= diagonal_score_seq(seq2-i, seq1, AlphaR(A),n);
		}
		if(s > score) { score = s; v = i; }
		/*** fprintf(stderr,"%d; score = %d\n",s,score); **/
		/*** PutDiagonalSeq(stdout, i, E1, E2, A); TEST **/
		hits=TRUE;
	}
	/** output high scoring segment **/
	if(hits) PutDiagonalSeq(stdout, v, E1, E2, A);
	else printf("no match\n");
	/** deallocate memory **/
	NilMheap(H); free(D);
	for(i=0; i<= nAlpha(A); i++) free(Q[i]);
	free(Q); free(nQ); 
	free(lastpos0); free(key0); free(off0);
	return score;
}

Boolean	RelateSeqFastp2(e_type E1, e_type E2, a_type A,long score)
/** use this to test fastp alignment, etc. **/
{
	if(score <= AlignSeqFastp(E1, E2, A)) return TRUE;
	else return FALSE;
}

Boolean	RelateSeqFastp(e_type E1, e_type E2, a_type A,long score)
/*******************************************************************
 See if two sequences are related with at or above a cutoff  score 
     using the fastp algorithm.
 *******************************************************************/
{
	long	i,x,*D,**Q,*nQ,r,s,n,n1,n2,item,nsave=10;
	char	*seq1,*seq2;
	mh_type	H;
	long	*off,*off0,*lastpos,*lastpos0; 
	keytyp	*key,*key0;
	long	maxspace=10;

	n1 = LenSeq(E1); n2 = LenSeq(E2);
	seq1 = XnuSeqPtr(E1); seq2 = XnuSeqPtr(E2);
	NEW(nQ,nAlpha(A)+3,long);
	MEW(Q,nAlpha(A)+3,long*);
	for(i=0; i<= nAlpha(A); i++) MEW(Q[i],n1+3,long);
	for(i=1; i <= n1; i++) {
		r = seq1[i]; Q[r][nQ[r]] = i; nQ[r]++;
	}
	n = n1 + n2; 
	NEW(off0,n+4,long); off = off0 + n2 + 2;
	NEW(lastpos0,n+4,long); lastpos = lastpos0 + n2 + 2;
	NEW(key0,n+4,keytyp); key = key0 + n2 + 2;
	for(i=1; i <= n2; i++) {
		r = seq2[i];
		for(n=0; n < nQ[r]; n++){
			x = Q[r][n] - i; /** find offset **/
			off[x]++;
			if(lastpos[x]!=0){
				s = i - lastpos[x]; /** get spacer **/
				if(s <= maxspace) key[x]+=4.0/(keytyp)s;
			} else lastpos[x] = i;
		}
	}
       /** locate regions of similarity using offsets with high # matches **/
	H = Mheap(nsave,3);
	MEW(D,nsave+2,long);
	for(i = 1-n2; i < n1; i++){
		if(i < 0) n = n2 + i; /** get length of diagonal **/
		else n = n1 - i;
		key[i] += (keytyp)off[i] - (keytyp)n/(keytyp)nAlpha(A);
		if(key[i] > 3.0) {
			item=InsertMheap(-key[i],H);
			D[item] = i;
		}
	}
        /** compute the score for the ten offsets of highest similarity **/
	for(i=0; i<= nAlpha(A); i++) free(Q[i]);
	free(Q); free(nQ); free(lastpos0); free(key0); free(off0);
	while((item=DelMinMheap(H))!=0){
		i = D[item];
		if(i > 0){	/** -> start at seq1[i] **/
		   n = tMIN(long,n1-i,n2); 
		   if(relate_seq(seq1+i,seq2,AlphaR(A),n,score)){
			NilMheap(H); free(D); return TRUE;
		   }
		} else {	/** -> start at seq2[i] **/
		   n = tMIN(long,n2+i,n1); 
		   if(relate_seq(seq2-i,seq1,AlphaR(A),n,score)){
			NilMheap(H); free(D); return TRUE;
		   }
		}
	}
	NilMheap(H); free(D); return FALSE;
}

