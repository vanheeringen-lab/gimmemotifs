#include "sequence.h"

e_type	ReadSeq(FILE *fptr, long I, long size, a_type A)
{
	char	*id,*seq;
	long	length,i,c='x';
	e_type	E;

        while(c != EOF){ if((c=fgetc(fptr)) == '>') break; }
	if(c==EOF) return NULL;
	NEW(id,80,char);
        for(i=0; (c=fgetc(fptr))!=EOF; i++){ 
		if(c=='\n') { if(i<80)id[i]= '\0'; break; }
		else if(i<79) id[i] = c;
		else if(i==79) id[i] = '\0';
	}
	NEW(seq,size+1, char);
        for(length=1; (c=fgetc(fptr))!=EOF; ){ 
		if(c == '>') { ungetc(c,fptr); break; }
		else if(isalpha(c)) {
			if(islower(c)) c = toupper(c);
			seq[length++] = AlphaCode(c,A);
		}
        }
	length--;
	if(length > size) seq_error("size and length inconsistency");
	NEW(E,1,sequence_type); 
	E->info=id;  E->n = (unsigned short) length;
	E->I = (unsigned short) I; 
	E->S = seq;  E->X = seq; E->xnu = FALSE;
	return E; 
}

Boolean	NonNullIdentSeqs(register e_type E1, register e_type E2)
/** If E1 has the same sequence as E2 return TRUE; else return FALSE **/
{
	register long	s,r1,r2;

	if(E1->n != E2->n) return FALSE;
	for(s=1; s<= (long) E1->n; s++){
		r1 = E1->X[s]; r2 = E2->X[s]; 
		if(r1 != 0 && r2 != 0 && r1 != r2) return FALSE;
	}
	return TRUE;
}

Boolean	IdentSeqs(e_type E1, e_type E2)
/** If E1 has the same sequence as E2 return TRUE; else return FALSE **/
{
	long	s;

	if(E1->n != E2->n) return FALSE;
	for(s=1; s<= (long) E1->n; s++){
		if(E1->X[s] != E2->X[s]) return FALSE;
	}
	return TRUE;
}

e_type	MergeSeqs(e_type *E)
/** merge array of sequences into one long sequence **/
{
	e_type mE;
	long  length,s,i,n;
	unsigned short I;

	if(E[0]==NULL) return NULL;
	for(I=length=n=0; E[n] != NULL; n++){
		length += E[n]->n;
		I = tMAX(unsigned short,I,E[n]->I);
	} I++;
	if(n==1) return CopySeq(E[0]);
	mE = EmptySeq(I,length);
	NEW(mE->info,80,char);
	for(i=0; i<80; i++) mE->info[i] = E[0]->info[i];
	for(s=1,n=0; E[n] != NULL; n++){
	   for(i=1; i<= (long) E[n]->n; i++,s++){
		mE->S[s] = E[n]->S[i];
	   }
	}
	return mE;
}

e_type	EmptySeq(long I, long length)
/*  create and return an empty sequence of length */
{
	char	*seq;
	e_type	E;

	NEW(seq,length+1, char);
	NEW(E,1,sequence_type);
	E->info=NULL;  E->n = (unsigned short) length;
        E->I = (unsigned short) I;
        E->S = seq;  E->X = seq; E->xnu = FALSE;
        return E;
}

long	GetFastaInfo(char *DBS_NAME, long max, long **pcounts, 
	unsigned short **psize, a_type A)
/*****************************************************************
 Get information on the number of sequences, lengths (=size) and residue 
    compositions for a fasta file 
 *****************************************************************/
{
	long	i,length;
	long	*counts;
	unsigned short	*size;
        /* 5-6-00 tlb; make c int */
	int c;
        FILE    *fptr;

        if((fptr = fopen(DBS_NAME,"r")) == NULL) {
                fprintf(stderr,"Could not open file \"%s\"\n",DBS_NAME);
                print_error("File does not exist.");
        }
	NEW(size,max+2,unsigned short);
	NEW(counts,nAlpha(A)+3,long);
        while((c=fgetc(fptr))!=EOF){ if(c=='>') break; }
	if(c==EOF) print_error("File not in fasta format!");
	/*** fprintf(stderr,"Determining database composition...");****/
        for(i=1,length=0;c!=EOF;length=0,i++) {
                if(c=='>') while((c=fgetc(fptr))!=EOF){ if(c=='\n') break; }
                while(c!='>') {
                   if(isalpha(c)) {
			length++; counts[(int)(AlphaCode(c,A))]++;
                   } else if(!isspace(c)) {
                        fprintf(stderr,"sequence %ld: %c undefined ",i,c);
			/****
                        fprintf(stderr," (ignored)\n");
			****/
                        fprintf(stderr," (replaced with null character: %c)\n",
				AlphaChar(UndefAlpha(A),A));
			length++; counts[UndefAlpha(A)]++;
                   }
                   if((c=fgetc(fptr))==EOF) break;
                }
                if(i > max) print_error("too many sequences; reset maximum");
                size[i] = length;
        } i--; 
        fclose(fptr);
	/*** fprintf(stderr," (%d sequences)\n",i); ***/
	*pcounts = counts; *psize = size;
	return i;
}

long    *CountsSeq(e_type E, a_type A)
{
	long j,r,*cnts;

	NEW(cnts,nAlpha(A)+3,long);
	for(j=1; j <= (long)E->n; j++){ r=E->S[j]; cnts[r]++; }
	return cnts;
}

double	*FreqResSeq(e_type E, double *freq, a_type A)
{
	long j,r;

	for(r = 0; r <= nAlpha(A); r++) freq[r] = 0.0;
	for(j=1; j <= (long)E->n; j++) { r=E->S[j]; freq[r]+=1.0; }
	for(r = 0; r <= nAlpha(A); r++) freq[r] /= (double) E->n;
	return freq;
}

long	*NumResSeq(e_type E,a_type A)
/** return residue frequencies **/
{
	long 	j,r;
	long	*num;

	NEW(num,nAlpha(A)+3,long);
	for(j=1; (long)j <= (long)E->n; j++) {
		r = E->S[j];
		num[r]++;
	}
	return num;
}

e_type	CopySeq(e_type E)
{
	e_type E2;
	long 	i,len;

	NEW(E2,1,sequence_type); 
	if(E==NULL) seq_error("CopySeq( ) can't copy null sequence.");
	if(E->info != NULL){
		len = strlen(E->info);
		NEW(E2->info,len + 2,char);
		for(i=0;i<=len; i++) E2->info[i]= E->info[i];
	}
	E2->n = E->n; E2->I = E->I; E2->xnu = E->xnu;
	NEW(E2->S, E->n + 2, char);
	for(i=1;i<=(long)E->n; i++) E2->S[i]= E->S[i];
	if(E->xnu){
		NEW(E2->X, E->n + 2, char);
		for(i=1;i<=(long)E->n; i++) E2->X[i]= E->X[i];
	} else { E2->X = E2->S; }
	return E2;
}

/************************* Randomization Routines ***********************/
e_type	RandomSeq(long length, long I, double *freq, a_type A)
/* return a randomly generated new sequence of length with residue 
   frequency " freq" */
{
	long 	i,c;
	double	r;
	e_type	E;

	E = EmptySeq(I, length);
	for(i=1;i<=(long)E->n; i++) {
		r = (double) Random()/(double) LONG_MAX; /* 0 <= r <= 1 */
		for(E->S[i]=0,c=1; c <= nAlpha(A); c++){
			r=r-freq[c];
			if(r <= 0.0){ E->S[i]=c; break; }
		}
	}
	return E;
}

e_type	RandomizeSeq(e_type E, double *freq, a_type A)
/* destroy E and return a randomly generated new sequence of same 
length with residue frequency " freq" */
{
	long 	i,c;
	double	r;

	if(E->xnu) {free(E->X); E->xnu=FALSE; E->X=E->S;}
	for(i=1;i<=(long)E->n; i++) {
		r = (double) Random()/(double) LONG_MAX; /* 0 <= r <= 1 */
		for(E->S[i]=0,c=1; c <= nAlpha(A); c++){
			r=r-freq[c];
			if(r <= 0.0){ E->S[i]=c; break; }
		}
	}
	if(E->info!=NULL)free(E->info); E->info = NULL; return E;
}

char	RandomResSeq(register e_type E)
/* return a randomly obtained residue in sequence E */
{
	register long	i;
	register double	r;

	do {
		r = (double) Random()/(double) LONG_MAX; /* 0 <= r <= 1 */
		i = (long) (r*E->n + 1); 
	} while(i < 1 && i > (long) E->n);
	return E->S[i];
}

e_type	RtnShuffleSeq(e_type E)
/* randomly permute the sequence E.S using a heap */
{
	e_type	rE;
	long 	i,n;
	dh_type	H;

	NEW(rE,1,sequence_type); 
	rE->xnu=FALSE; rE->I = 0; n = rE->n = E->n;
	NEW(rE->S,n+1,char);
	H = dheap(n+1,4);
	for(i=1;i<=n;i++){ insrtHeap(i,((keytyp)Random()),H);}
	for(i=1;i<=n;i++) rE->S[i] = E->S[delminHeap(H)];
	Nildheap(H); rE->info=NULL; rE->X = rE->S;
	return rE;
}

e_type	ShuffleSeq(e_type E)
/* randomly permute the sequence E.S using a heap; don't move 'X's. */
{
	char	*S;
	long 	i;
	dh_type	H;

	if(E->xnu) {free(E->X); E->xnu=FALSE; }
	NEW(S,E->n+1,char);
	H = dheap(E->n+1,4);
	for(i=1;i<=(long)E->n; i++){
		S[i]=E->S[i];
		if(S[i] != 0){
			insrtHeap(i,((keytyp)Random()),H);
		}
	}
	for(i=1;i<=(long)E->n; i++){
		if(E->S[i] != 0){
			E->S[i] = S[delminHeap(H)];
		}
	}
	Nildheap(H); free(S); 
	if(E->info!=NULL)free(E->info); E->info = NULL; 
	E->X = E->S;
	return E;
}

e_type	ShuffleSeq2(e_type E)
/* randomly permute the sequence E.S using a heap */
{
	char	*S;
	long 	i;
	dh_type	H;

	if(E->xnu) {free(E->X); E->xnu=FALSE; }
	NEW(S,E->n+1,char);
	H = dheap(E->n+1,4);
	for(i=1;i<=(long)E->n; i++) 
		{ S[i]=E->S[i];insrtHeap(i,((keytyp)Random()),H);}
	for(i=1;i<=(long)E->n; i++) E->S[i] = S[delminHeap(H)];
	Nildheap(H); free(S); 
	if(E->info!=NULL)free(E->info); E->info = NULL; 
	E->X = E->S;
	return E;
}

/************************* Print Routines ***********************/
void	PutSeqInfo(FILE *fptr,e_type E)
{	if(E->info !=NULL) fprintf(fptr,"%s\n", E->info); 
	else fprintf(fptr,"random%ld\n",(long)E->I);}

void	PutSeqID(FILE *fptr,e_type E)
{
	long	k;

	if(E->info !=NULL) {
	   for(k=0; !isspace(E->info[k]); k++){
		if(E->info[k]==0) break;
		fprintf(fptr,"%c", E->info[k]);
	   }
       	} else fprintf(fptr,"random%ld ",(long)E->I);
}

void    PutXnuSeq(FILE *fptr,e_type E,a_type A)
{
	long 	j;

	fprintf(fptr,">"); PutSeqInfo(fptr,E);
	for(j=1; (long)j <= (long)E->n; j++) {
	   fprintf(fptr,"%c",AlphaChar(E->X[j],A));
	   if(j%10 == 0) fprintf(fptr," ");
	   if(j%50 == 0) fprintf(fptr,"\n");
	}
	fprintf(fptr,"\n\n");
}

void    PutSeq(FILE *fptr,e_type E,a_type A)
{
	long 	j;

	fprintf(fptr,">"); PutSeqInfo(fptr,E);
	for(j=1; (long)j <= (long)E->n; j++) {
	   fprintf(fptr,"%c",AlphaChar(E->S[j],A));
	   if(j%10 == 0) fprintf(fptr," ");
	   if(j%50 == 0) fprintf(fptr,"\n");
	}
	fprintf(fptr,"\n\n");
}

void    MaskSeq(long start, long end, e_type E)
/** set the region from start to end to "null" residues (0) **/
{
	long	i;
	char	*s;
	
	if(start < 1 || start > end || end > (long) LenSeq(E)){
		seq_error("MaskSeq( ) - out of range.");
	} 
	if(E->xnu) s = E->X;
	else s = E->S;
	for(i=start; i <= end; i++) s[i] = 0;
}

void    PutSeqRegion(FILE *fptr,long start, long length, e_type E, a_type A)
{ PutSeqRegion2(fptr,start, length, E, 10, A); }

void    PutSeqRegion2(FILE *fptr,long start, long length, e_type E, 
	long flank, a_type A)
{
	long	e,end,i;

	fprintf(fptr,"%4ld  ",start);
	e = start + length - 1;
	end = e + flank;
	for(i=start-flank; i <= end; i++){
		if(i < 1 || i > (long) E->n) fprintf(fptr," ");
		else if(i == e) {
			fprintf(fptr,"%c ", AlphaChar(E->S[i],A));
		} else if(i == start) {
			fprintf(fptr," %c", AlphaChar(E->S[i],A));
		} else {
			fprintf(fptr,"%c", AlphaChar(E->S[i],A));
		}
	}
	e = tMIN(long,e,LenSeq(E));
	fprintf(fptr," %4ld",e);
}

/** a_type	ALPHABETA=NULL; ** DEBUG **/
/* 5-6-00 tlb; change R to int */
long	get_diagonal_ends_seq(char *seq1, char *seq2, int **R, long n,
	long *begin, long *end)
{
	register long min=LONG_MAX, sum=0,score=-9999;
	long	min_n=0;

	if(n <= 0) { *begin = n+1; *end = n; return 0; }
	score = min = sum = R[(int)seq1[n]][(int)seq2[n]];
	*begin = *end = min_n = n;
	n--;
	while(n > 0){
        	if(min > (sum += R[(int)seq1[n]][(int)seq2[n]])){
			min = sum; min_n = n-1;
		} 
		if(score < (sum-min)){
			score = (sum-min);
			*end = min_n; *begin = n;
		}
		/***** fprintf(stderr,"%4d: %c-%c (%d (sum = %d): %d..%d)\n",
			n,AlphaChar(seq2[n],ALPHABETA),
			AlphaChar(seq1[n],ALPHABETA),score, sum,
			*begin,*end); *****/
		n--;
        };
	/*** printf("\tEND\n\n"); ****/
	return score;
}

long	CenterSeqHSP(long offset, e_type E1, e_type E2, a_type A)
/* return the center of the high scoring pair with offset **/
{
	long	v,n,n1,n2,score=0,e,start,end;
	char	*seq1,*seq2;
	long	center;

	n1 = LenSeq(E1); n2 = LenSeq(E2);
	v = offset;
	if(v < (1-n2) || v >= n1) seq_error("offset out of range");
	seq1 = XnuSeqPtr(E1); seq2 = XnuSeqPtr(E2);
	if(v > 0){		/** case 1: offset **/
	   n = tMIN(long,n1-v,n2);
	   score=get_diagonal_ends_seq(seq1+v,seq2,AlphaR(A),n,&start,&e);
	   start = v+start; end = v+e-1;
	   
	} else {
	   n = tMIN(long,n2+v,n1);
	   score=get_diagonal_ends_seq(seq2-v,seq1,AlphaR(A),n,&start,&e);
	   end = e - 1;
	}
	n = (end - start + 1)/2;
	center = start + n;
	return center;
}

long	PutDiagonalSeq(FILE *fptr, long offset, e_type E1, e_type E2, a_type A)
/*******************************************************************
 print the MSP for diagonal at offset between seq1 and seq2.

 *******************************************************************/
{
	long	v,i,j,n,n1,n2,score=0,b,e,start,end,w=40;
	char	*seq1,*seq2;

/** ALPHABETA=A; *****/
	n1 = LenSeq(E1); n2 = LenSeq(E2);
	v = offset;
	if(v < (1-n2) || v >= n1) seq_error("offset out of range");
	seq1 = XnuSeqPtr(E1); seq2 = XnuSeqPtr(E2);
	if(v > 0){
	   n = tMIN(long,n1-v,n2);
	   score=get_diagonal_ends_seq(seq1+v,seq2,AlphaR(A),n,&start,&e);
	   fprintf(fptr,"\n\n");
	   for(b=start; b <= e; b+=w){
		end = tMIN(long,e,b+w-1);
		fprintf(fptr,"%4ld ",v+b);
		for(i=v+b,j=b; j <= end; i++,j++)
			fprintf(fptr,"%c",AlphaChar(seq1[i],A));
	  	fprintf(fptr," %4ld ",v+end);
		PutSeqID(fptr,E1);
	  	fprintf(fptr,"\n     ");
		for(i=v+b,j=b; j <= end; i++,j++)
			if(seq1[i]==seq2[j]) fprintf(fptr,":");
			else if(valAlphaR(seq1[i],seq2[j],A) >= 0)
							fprintf(fptr,".");
			else fprintf(fptr," ");
		fprintf(fptr,"\n%4ld ",b);
		for(j=b; j <= end; j++)
			fprintf(fptr,"%c", AlphaChar(seq2[j],A));
		fprintf(fptr," %4ld ",end);
		PutSeqID(fptr,E2);
		fprintf(fptr,"\n\n");
	   }
	   fprintf(fptr,"   score = %ld\n\n",score);
	} else {
	   n = tMIN(long,n2+v,n1);
	   score=get_diagonal_ends_seq(seq2-v,seq1,AlphaR(A),n,&start,&e);
	   fprintf(fptr,"\n\n");
	   for(b=start; b <= e; b+=w){
		end = tMIN(long,e,b+w-1);
		fprintf(fptr,"%4ld ",b);
		for(j=b; j <= end; j++)
			fprintf(fptr,"%c",AlphaChar(seq1[j],A));
		fprintf(fptr," %4ld ",end);
		PutSeqID(fptr,E1);
	  	fprintf(fptr,"\n     ");
		for(i=b-v,j=b; j <= end; i++,j++)
			if(seq2[i]==seq1[j]) fprintf(fptr,":");
			else if(valAlphaR(seq2[i],seq1[j],A) >= 0)
							fprintf(fptr,".");
			else fprintf(fptr," ");
		fprintf(fptr,"\n%4ld ",b-v);
		for(i=b-v,j=b; j <= end; i++,j++)
			fprintf(fptr,"%c", AlphaChar(seq2[i],A));
	  	fprintf(fptr," %4ld ",end-v);
		PutSeqID(fptr,E2);
		fprintf(fptr,"\n\n");
	   }
	   fprintf(fptr,"   score = %ld\n\n",score);
	}
	return score;
}

/********************************************************************
   ** C program to read in sequences in Fasta format
   ** test them for internal repeats
   ** and print them out with the internal repeats flagged
	Usage: xnu fasta-sequences-file\n
	         [-60] [-120] [-250]     PAM matrix to use\n");
	         [-p probability-cut]    default 0.01\n");
	         [-s score-cut]\n");
	         [-n search-width]       default 10, 0->all diagonals\n");
	         [-a] [-d]               ascending or descending\n");
	         [-.] [-x] [-o]          output options\n");
		(argv[i],"-p")==0) { pcut = atof(argv[++i]); }
		(argv[i],"-n")==0) { ncut = atoi(argv[++i]); }
		(argv[i],"-s")==0) { scut = atoi(argv[++i]); }

		(argv[i],"-60")==0) { m = pam60; lambda = lambda60;}
		(argv[i],"-250")==0) { m = pam250; lambda = lambda250;}

		(argv[i],"-a")==0) { ascend = 1; descend = 0;}
		(argv[i],"-d")==0) { ascend = 0; descend = 1;}
*********************************************************************/

void	ProcessSeq(e_type E,a_type A,double lambda,double K,double H)
{
	long	i,j,k,off,sum,beg,end,top,noff;
	long	topcut,fallcut,len;
	char 	*is,*hit;
	double	s0;
	long ascend=1;
	long descend=1;
	long ncut=4;
	long mcut=1;
	double pcut=0.001;
	long scut=0;

	if(E->xnu) return;	/* already done */
	if ((len=LenSeq(E)) == 0) return;
	NEW(E->X,len+2,char); E->xnu = TRUE;
	is = E->S +1; hit = E->X +1;
	for (i=0; i<len; i++) hit[i]=FALSE;
	noff = len-1;
	if(ncut>0) noff=ncut;
	if(scut!=0) topcut = scut;
	else { 
		s0 = - log( pcut*H / (noff*K) ) / lambda;
		if (s0>0) topcut = floor(s0 + log(s0)/lambda + 0.5);
		else topcut = 0;
	}
	fallcut = (long)log(K/0.001)/lambda;
	for (off=mcut; off<=noff; off++) {
		sum=top=0; beg=off; end=0;
		for(i=off; i<len; i++) {
			sum += valAlphaR(is[i],is[i-off],A);
			if (sum>top) { top=sum; end=i; }
			if (top>=topcut && top-sum>fallcut) {
				for (k=beg; k<=end; k++) {
					if (ascend) hit[k] = TRUE;
					if (descend) hit[k-off] = TRUE;
				}
				sum=top=0; beg=end=i+1;
			} else if (top-sum>fallcut) {
				sum=top=0;
				beg=end=i+1;
			}
			if (sum<0) { beg=end=i+1; sum=top=0; }
		}
		if (top>=topcut) {
			for (k=beg; k<=end; k++) {
				if (ascend) hit[k] = TRUE;
				if (descend) hit[k-off] = TRUE;
			}
		}
	}
	for (i=0; i<len; i+=60) {
		for (j=i; j<i+60 && j<len; j++) {
		   if (hit[j]) { hit[j] = UndefAlpha(A); }
		   else { hit[j] = is[j]; }
		}
	}
}

e_type  ReadSeqFA(char *infile, long I, a_type A)
/** read a sequence in from a fasta file ***/
{
        /* 5-6-00 tlb; make c int */
	int c;
        long     length;
	FILE	*fptr;
	e_type	E;

        if((fptr = fopen(infile,"r")) == NULL) {
                fprintf(stderr,"Could not open file \"%s\"\n",infile);
                seq_error("File does not exist!");
        }
        length=0;
        while((c=fgetc(fptr))!=EOF){ if(c=='>') break; }
        if(c=='>') while((c=fgetc(fptr))!=EOF){ if(c=='\n') break; }
        while(c!='>') {
             if(isalpha(c)) length++;
             else if(!isspace(c)) {
                  fprintf(stderr,"seq %ld: illegal character -> %c",I,c);
                  while((c=fgetc(fptr)) != EOF) {
                          fprintf(stderr,"%c",c);
                          if(c == '\n') break;
                  }
                  fprintf(stderr,"\n");
                  seq_error("fatal error.");
             }
             if((c=fgetc(fptr))==EOF) break;
        }
        fclose(fptr);
        if((fptr = fopen(infile,"r")) == NULL) {
                fprintf(stderr,"Could not open file \"%s\"\n",infile);
                seq_error("File does not exist!");
        }
	E = ReadSeq(fptr, I, length, A);
        fclose(fptr);
	return E;
}

e_type	NilSeq(e_type E)
{
	if(E!=NULL){
		free(E->S);
		if(E->info!=NULL)free(E->info);
		if(E->xnu) free(E->X);
		free(E);
	}
	return NULL;
}

long	seq_error(char *s){fprintf(stderr,"Seq: %s\n",s);exit(1);}

