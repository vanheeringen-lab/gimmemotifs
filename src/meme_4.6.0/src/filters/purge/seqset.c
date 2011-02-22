#include "seqset.h"

ss_type	MkXnuSeqSet1(char *filename, e_type E, a_type A) 
/** return NULL if can't xnu **/
{ return xnu_seqset(SeqSet1(filename, E, A)); }

ss_type	MkXnuSeqSet(char *filename,a_type A) 
/** return NULL if can't xnu **/
{ return xnu_seqset(seqset(filename,A)); }

ss_type	MkXnuSeqSet_fptr(FILE *fptr,a_type A) 
/** return NULL if can't xnu **/
{ return xnu_seqset(fptr_seqset(fptr,A)); }

ss_type	SeqSet1(char *filename, e_type E,a_type A)
/* create a seqset of one sequence E */
{
	ss_type	P;
	e_type	E2 = CopySeq(E);
	long	s,r;

	NEW(P,1,seqset_type);
	P->name = String(filename); P->A = A;
	P->nent = 1; 
	NEW(P->entity,2,e_type); 
	NEW(P->counts,nAlpha(A)+1,long);
	P->entity[1]=E2;
	for(s=1; s<= (long) LenSeq(E2); s++){
		r = ResSeq(s,E2);
		P->counts[r]++;
	}
	P->max_leng = P->min_leng = P->total = LenSeq(E2);
	P->tfreq = NULL;
	calcseqsetfreq(P);
	P->xnu = FALSE;
	return (P);
}

ss_type	RmSeqSet(e_type E, ss_type P)
/** remove all sequences in P that are identical to E **/
{
	e_type	E2;
	long	s,r,n,k;

        for(n=1;n<=P->nent;){
		E2 = P->entity[n];
		if(IdentSeqs(E2, E)){
		   for(s=1; s<= (long) LenSeq(E2); s++){
			r = ResSeq(s,E2); P->counts[r]--;
		   }
		   NilSeq(E2);
		   P->entity[n] = P->entity[P->nent];
		   P->nent--; 
		} else n++;
	}
	P->total = P->max_leng = 0; 
	P->min_leng = LONG_MAX;
        for(n=1;n<=P->nent;){
		E2 = P->entity[n];
		k = LenSeq(E2);
		if(P->max_leng < k) P->max_leng = k;
		if(P->min_leng > k) P->min_leng = k;
		P->total += LenSeq(E2);
	}
	calcseqsetfreq(P);
	return P;
}

ss_type	SeqSet(char *filename,a_type A) { return seqset(filename,A); }

ss_type	SeqSet_fptr(FILE *fptr,a_type A) 
{ return fptr_seqset(fptr,A); }

ss_type	fptr_seqset(FILE *fptr,a_type A) 
/* create a seqset from the input file with segment length k and */
/* alphabet A. */
{
	ss_type	P;
	long	i,s,r,nsize[MAX_NUMBER_SEQS+1];
	e_type	E;

	NEW(P,1,seqset_type);
	NEW(P->name,25,char);
	strcpy(P->name,"temp_file"); P->A = A;
	P->nent = count_seqset_entities(fptr,P,nsize); 
	NEW(P->entity,P->nent+1,e_type); 
	NEW(P->counts,nAlpha(A)+1,long); P->total = 0;
	rewind(fptr);
	for(i=1; i<=P->nent; i++){
	   E = ReadSeq(fptr,i,nsize[i],A);
	   P->entity[i]=E;
	   for(s=1;s<=(long)LenSeq(E);s++){r=ResSeq(s,E);P->counts[r]++;}
	   P->total += LenSeq(E);
	}
	P->tfreq = NULL; calcseqsetfreq(P); P->xnu = FALSE;
	return (P);
}

ss_type	seqset(char *filename,a_type A) 
/* create a seqset from the input file with segment length k and */
/* alphabet A. */
{
	ss_type	P;
	FILE	*fptr;
	long	i,s,r,nsize[MAX_NUMBER_SEQS+1];
	e_type	E;

	NEW(P,1,seqset_type);
	NEW(P->name,strlen(filename)+2,char);
	strcpy(P->name,filename); P->A = A;
	fptr = OpenSeqSetFile(P);
	P->nent = count_seqset_entities(fptr,P,nsize); 
        fclose(fptr);
	NEW(P->entity,P->nent+1,e_type); 
	NEW(P->counts,nAlpha(A)+1,long); P->total = 0;
	fptr = OpenSeqSetFile(P);
	for(i=1; i<=P->nent; i++){
	   E = ReadSeq(fptr,i,nsize[i],A);
	   P->entity[(int)i]=E;
	   for(s=1;s<=(long)LenSeq(E);s++){r=ResSeq(s,E);P->counts[r]++;}
	   P->total += LenSeq(E);
	}
	fclose(fptr); P->tfreq = NULL; calcseqsetfreq(P); P->xnu = FALSE;
	return (P);
}

FILE    *OpenSeqSetFile(ss_type P)
{
        FILE    *fptr;
        if((fptr = fopen(P->name,"r")) == NULL) {
                fprintf(stderr,"Could not open file \"%s\"\n",P->name);
                seqset_error("File does not exist!\n");
        }
        return fptr;
}

ss_type	NilSeqSet(ss_type P)
{
	long i;
	for(i=1; i<=P->nent; i++){
		if(P->entity[i] != NULL) NilSeq(P->entity[i]); 
	}
	free(P->name); free(P->entity); free(P->counts);
	if(P->tfreq != NULL) free(P->tfreq); 
	free(P);
	return (ss_type) NULL;
}

ss_type	xnu_seqset(ss_type P)
/*** return NULL if can't xnu SeqSet ***/
{
        long	s,i,j,low,high,len,v;
        double  *pr,lambda,K,H;
        char    r,r1,r2;
	e_type E;
	double	*freq;
	a_type A=P->A;
	char xnualpha[] = {"ARNDCQEGHILKMFPSTWYVBZX*-"};
	double xnufreq[20] = {
        0.081, 0.057, 0.045, 0.054, 0.015, 0.039, 0.061, 0.068, 0.022, 0.057,
        0.093, 0.056, 0.025, 0.040, 0.049, 0.068, 0.058, 0.013, 0.032, 0.067 };

	freq = tFreqSeqSet(P);
	if(P->xnu) return P;
	else P->xnu = TRUE;
        low = lowAlphaR(A);
        high = highAlphaR(A);
        len = high - low + 1;
        NEW(pr,len+1,double);
        for(i=0; i<=nAlpha(A); i++){
	   P->counts[i]=0;
           for(j=0; j<=nAlpha(A); j++){
                v = valAlphaR(i,j,A) - low;
                pr[v] += freq[i] * freq[j];
           }
        }
        if(!karlin(low,high,pr,&lambda,&K,&H)) {
		fprintf(stderr,"\nusing blast amino acid frequencies\n");
		if(nAlpha(A) == 20) {	
        	   for(i=0; i<= len; i++) pr[i] = 0.0;/****/
        	   for(i=0; i<20; i++){
			r1 = AlphaCode(xnualpha[i],A);
    			for(j=0; j<20; j++){
			   r2 = AlphaCode(xnualpha[j],A);
			   v = valAlphaR(r1,r2,A) - low;
                	   pr[v] += xnufreq[i] * xnufreq[j];
           		}
        	    }
        	    if(!karlin(low,high,pr,&lambda,&K,&H)) {
			seqset_error("this should not happen.");
		    }
		} else seqset_error("fatal error in xnu_seqset( ).");
	} else H=ExpectedInformation(A, lambda, freq);
        for(i=1;i<=NSeqsSeqSet(P);i++){
	   E = P->entity[i];
           if (LenSeq(E) != 0) {
        	ProcessSeq(E,A,lambda,K,H);
		for(s=1; s<= (long) LenSeq(E); s++){
		   r = XnuSeq(s,E); P->counts[(int)r]++;
		}
	   }

	}
        free(pr);
	calcseqsetfreq(P);
	return P;
}

double  LogL0SeqSet(ss_type P)
{
        double  *freq,L0,n,r;
        long     b;

        freq = tFreqSeqSet(P);
        for(L0=0.0, b=1; b<= nAlpha(P->A); b++){
                    if(CountsSeqSet(b,P) > 0){
                        n = (double) CountsSeqSet(b,P);
                        r = freq[b];
                        L0 += n * log(r);
                    }
        }
        return (1.4427*L0);
}

long	*LengthsSeqSet(ss_type P)
/* returns an array containing the sequence lengths */
{
	long	*len_seq,n;

	NEW(len_seq,NSeqsSeqSet(P) +1,long);
	for(n=1; n<= NSeqsSeqSet(P); n++)len_seq[n]=SqLenSeqSet(n,P);
	return len_seq;
}

/******************** Counting and Numbering Operations *******************/
long     count_seqset_entities(FILE *fptr,ss_type P, long nsize[])
{
        long i=0,j,length,c; 

	P->max_leng = 0; P->min_leng = 1000000; 
	while((c=fgetc(fptr))!=EOF){ if(c=='>') break; }
        for(i=1,length=0;c!=EOF;length=0,i++) { 
		if(c=='>') while((c=fgetc(fptr))!=EOF){ if(c=='\n') break; }
		while(c!='>') {
		   if(isalpha(c)) length++;
		   else if(!isspace(c)) {
			fprintf(stderr,"seq %ld: illegal character -> %c",i, (char)c);
			for(j=0; (c=fgetc(fptr)) != EOF; j++) {
				if(c == '\n') break;
				fprintf(stderr,"%c",(char)c);
				if(j > 10 || isspace(c)) break;
			}  
			fprintf(stderr,"\n");
			seqset_error("input file error - fatal.");
		   } 
           	   if((c=fgetc(fptr))==EOF) break; 
	     	}
		if(i >= MAX_NUMBER_SEQS) 
		   seqset_error("too many sequences; reset MAX_NUMBER_SEQS");
	   	P->max_leng = tMAX(long,P->max_leng,length);
		P->min_leng = tMIN(long,P->min_leng,length);
		nsize[i] = length;
	}
	i--;
        return i;
}

ss_type	PutSeqSet(FILE *fptr,ss_type P)
{
	fprintf(fptr,"\n  input file:\n");
	fprintf(fptr,"\tname: \"%s\"\n\ttotal sequences: %ld",
			P->name,P->nent);
	fprintf(fptr,"\n\tsequence lengths: %ld-%ld residues\n",
		       P->min_leng,P->max_leng);
	return P;
}

/*********************** Put SeqSet Entities Operations **********************/
ss_type	PutSeqSetEs(FILE *fptr,ss_type P)
/* print all sequence entities in seqset using fasta format */
{
	long     i;
	for(i=1;i<=P->nent; i++) {
 	   if(SeqI(P->entity[i]) != 0) PutSeqSetE(fptr, i,P);
	}
	fprintf(fptr,"\n\n");
	return P;
}

ss_type	PutSeqSetE(FILE *fptr, long i, ss_type P) /* print the ith entity */
{
	e_type	E;

	if(i <= P->nent && i > 0) { E = P->entity[i]; PutSeq(fptr,E,P->A); }
	return P;
}

ss_type	PutSeqSetPIDs(FILE *fptr, ss_type P)
/* print entity ids for selected entities */
{
	e_type	E;
	long     i;
	for(i=1;i<=P->nent;i++) {
	   E = P->entity[i];
 	   if(SeqI(E) != 0){
		fprintf(fptr,"#%-3ld ",(long)SeqI(E));
		PutSeqInfo(fptr,E);
	   }
	}
	fprintf(fptr,"\n\n");
	return P;
}

/********************** Frequency Operations ********************/

ss_type	calcseqsetfreq(ss_type P)
/* calculate the residue frequencies for seqset */
{
    long		s;
    a_type	A=P->A;

    if(P->tfreq==NULL) NEW(P->tfreq,nAlpha(A)+1,double);
    for(s=0;s<=(long) nAlpha(A);s++) {
	P->tfreq[s] = (double) P->counts[s]/(double) P->total;
    }
    return P;
}

double	SeqSetEntropy(ss_type P)
{
	long	i;
	double	*freq,H;

	freq = tFreqSeqSet(P);
	for(H=0.0,i = 1; i <= nAlpha(P->A); i++){
		if(freq[i] > 0.0) H += freq[i] * log(freq[i]);
	}
	return (-1.442695041*H);
}

ss_type	PutSeqSettFreqs(FILE *fptr,ss_type P)
{
	long i; double T=0.0;

	fprintf(fptr,"RES    %-6s %-s\n","NUM","FREQ");
	for(i=0;i<=nAlpha(P->A);T+=P->tfreq[i],i++)
	    fprintf(fptr,"%c(%2ld): %-6ld %-2.3f\n",
			AlphaChar(i,P->A),i,P->counts[i],P->tfreq[i]);
	fprintf(fptr,"TOTAL: %-6ld %-2.3f\n\n", P->total, T);
	return P;
}

/************************* Randomization Routines ***********************/

ss_type	ShuffleSeqSet2(ss_type P)
{ 
	long r,s,i,n,item;
	dh_type H;
	e_type	E;
	char	*S;

	for(n=0,i=1;i<=P->nent;i++) { E=P->entity[i]; n += LenSeq(E); }
	H = dheap(n+2,4);
	NEW(S,n+2,char);
	for(item=i=1;i<=P->nent;i++) {
		E=P->entity[i];
		for(s=1; s<= (long) LenSeq(E); s++){
			r = ResSeq(s,E);
			insrtHeap(item,((keytyp)Random()),H);
			S[item++]=r;
		}
	}
	for(i=1;i<=P->nent;i++) {
		E=P->entity[i];
		for(s=1; s<= (long) LenSeq(E); s++){
			item=delminHeap(H);
			if(item==0) seqset_error("shuffleSeqSet2 error");
			r=S[item];
			EqSeq(s,r,E);
		}
		EqSeqI(i,P->entity[i]);
	}
	Nildheap(H); free(S);
	return P;
}

ss_type	ShuffleSeqSet(ss_type P)
{ 
	long i;
	for(i=1;i<=P->nent;i++) {
		ShuffleSeq(P->entity[i]);
		EqSeqI(i,P->entity[i]);
	}
	return P;
}

void	seqset_error(char *s)
{fprintf(stderr,"Seq_Set: %s\n",s); exit(1);}

