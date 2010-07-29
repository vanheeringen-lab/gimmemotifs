#include "alphabet.h"

a_type	MkAlphabetC(char *map_s,char *R,char *prs,char *comp)
/* plus complementary bases */
{
	a_type	A;
	long	i;
	char	c,d;

	A = MkAlphabet(map_s,R,prs);
	NEW(A->C,nAlpha(A)+2,char); /** if not defined set to 0 **/
	for(i=0;isalpha(comp[i]);i+=2){
		if(!isalpha(comp[i+1])) break;
		else {
			c = AlphaCode(comp[i],A);
			d = AlphaCode(comp[i+1],A);
			A->C[(int)c]=d; A->C[(int)d]=c;
		}
	}
	return A;
}

a_type	MkAlphabet(char *map_s,char *R,char *prs)
/* plus relatedness matrix input as ordered pairs */
{
	a_type	A;
	long i;

	A=MkAlpha(map_s,R);
	NEWP(A->pairs,nAlpha(A)+2,char);
	for(i=0;i<=nAlpha(A);i++) NEW(A->pairs[i],nAlpha(A)+2,char);
	DefAlphaP(prs,A);
	return A;
}

a_type	DefAlphaP(char *prs,a_type A)
{
	long i,j,n;
	char c,d;

	for(i=0;i<=nAlpha(A);i++) {
	   A->paired[i]=FALSE;
	   for(j=0;j<=nAlpha(A);j++) { 
		if(i==j) A->pairs[i][j] = 0; 
		else A->pairs[i][j] = 99; 
	   }
	}
	for(n=1,A->npairs=i=0;isalpha(prs[i]);i+=2){
		if(!isalpha(prs[i+1])) break;
		else {
			c = AlphaCode(prs[i],A);
			d = AlphaCode(prs[i+1],A);
			A->pairs[(int)c][(int)d]=A->pairs[(int)d][(int)c]=n++;
			A->paired[(int)c]=A->paired[(int)d]=TRUE;
			A->npairs++;
		}
	}
	if(A->prs!=NULL) free(A->prs);
	NEW(A->prs,strlen(prs)+2,char);
	for(j=0,i=1;i<=A->npairs;i++) {
		A->prs[j]=prs[j];j++;
		A->prs[j]=prs[j];j++;
	}
	return A;
}

a_type	MkAlpha(char *map_s,char *R)
{
	a_type	A;
	long	n,i;
	char	c;

	NEW(A,1,alphabet_type);
	for(n=0;(c=map_s[n])!= '\0';n++) 		/* Get # letters */
		if(!isalpha(c)) alpha_error("Illegal alphabet string",A);
	A->n = n-1;
	NEW(A->alphabet,A->n+2,char);			/* ALPHABET */
	strncpy(A->alphabet,map_s,nAlpha(A)+1);	
	NEW(A->code2let,(strlen(map_s)+1),char);	/* CODE2LETTER */
	strcpy(A->code2let,map_s); 
	NEW(A->code2lower,(strlen(map_s)+1),char); 	/* LOWER*/
	strcpy(A->code2lower,map_s); 
	NEW(A->let2code,ALPHA_NUM_SYMBOLS,char);	      /* LETTER2CODE */
	for(i=0;i<ALPHA_NUM_SYMBOLS;i++) A->let2code[i]= 0;	/* =error */
	for(i=0;map_s[i]!=0;i++) {
		c = map_s[(int)i]; A->let2code[(int)c] = i; 
		if(isupper(c)) { c = tolower(c); A->let2code[(int)c] = i; }
		else if(islower(c)) { c = toupper(c); A->let2code[(int)c] = i; }
	}
	for(i=0;map_s[i]!=0;i++) 
		if(isupper(map_s[i])) A->code2lower[i]=tolower(map_s[i]); 
	if(R==NULL) { A->R = NULL; }			/* RELATION */
	else {
/* tlb 5-6-00; make R type int */
	   NEWP(A->R,A->n+2,int);			
	   for(i=0;i<=nAlpha(A);i++) NEW(A->R[i],A->n+2,int);
	   A = DefAlphaR(R,A);
	}
	NEW(A->paired,nAlpha(A)+2,Boolean);
	for(i=0;i<=nAlpha(A);i++) A->paired[i]=FALSE;
	A->pairs = NULL;
	A->prs = NULL;
	A->C = NULL;
	return (A);
}

a_type	DefAlphaR(char *R,a_type A)
/* relatedness matrix */
{
	long i,j,value,lo,hi;
	if(strlen(R)==0) { A->R = NULL; return A; }
	for(lo=9999,hi=(-9999),i=0;i<=nAlpha(A);i++) {
	   for(j=0;j<=nAlpha(A);j++) {
		while(R[0] != '-' && !isdigit(R[0])) {
		   if(R[0]== 0) alpha_error("Illegal input in DefAlphaR()",A);
		   else R++;
		}
		if(sscanf(R,"%ld",&value) != 1)
			alpha_error("Illegal input in DefAlphaR()",A);
		A->R[i][j] = value;
		lo = tMIN(long,value,lo);
		hi = tMAX(long,value,hi);
		while(R[0] == '-' || isdigit(R[0])) R++;
	   }
	}
	A->loR = lo; A->hiR = hi;
	return (A);
}

a_type	PutAlpha(FILE *fptr,a_type A)
{
	long	i;

	fprintf(fptr,"Alphabet: %s (n=%ld)\n", A->alphabet,nAlpha(A));
	fprintf(fptr," Code:\n");
	for(i=0;i<=nAlpha(A);i++) fprintf(fptr,"  %c",A->code2let[i]);
	fprintf(fptr,"\n");
	for(i=0;i<=nAlpha(A);i++) fprintf(fptr,"%3ld",i);
	fprintf(fptr,"\n\n");
	PutAlphaR(fptr,A);
	fprintf(fptr,"\n");
	PutAlphaPairs(fptr,A);
	fprintf(fptr,"\n");
	return A;
}

void	PutAlphaPairs(FILE *fptr, a_type A)
{
	long	i,j;
	if(A->pairs != NULL) {
		fprintf(fptr, "  ordered pairs:\n     ");
		for(j=0,i=1; i<=A->npairs; i++){
			fprintf(fptr,"%c", A->prs[j]); j++;
			fprintf(fptr,"%c ", A->prs[j]); j++;
			if(i%10==0 && i < A->npairs){
				fprintf(fptr, "\n     ");
			}
		}
		fprintf(fptr,"\n");
	} else fprintf(fptr,"  ordered pairs:   NULL\n");
	fprintf(fptr,"\n");
}

void	PutAlphaR(FILE *fptr, a_type A)
{
	long	i,j;
	if(A->R != NULL) {
		fprintf(fptr, "  relatedness matrix:\n     ");
		for(i=0; i<=nAlpha(A); i++)
			fprintf(fptr,"%3c", AlphaChar(i,A));
		for(i=0; i<=nAlpha(A); i++) {
			fprintf(fptr,"\n   %c:",AlphaChar(i,A));
			for(j=0; j<=nAlpha(A); j++)
				fprintf(fptr,"%3d", valAlphaR(i,j,A));
		}
		fprintf(fptr,"\n");
	} else fprintf(fptr,"  R:   NULL\n");
	fprintf(fptr,"\n");
}

void	NilAlpha(a_type A)
{
	long	i;
	if(A != NULL) {
	   free(A->alphabet);
	   free(A->code2let);
	   free(A->code2lower);
	   free(A->let2code);
	   if(A->C != NULL) free(A->C);
	   if(A->R != NULL) {
		for(i=0;i<=nAlpha(A);i++) free(A->R[i]);
		free(A->R);
	   }
	   if(A->pairs != NULL) {
		for(i=0;i<=nAlpha(A);i++) free(A->pairs[i]);
		free(A->pairs);
	   }
	   if(A->prs!=NULL) free(A->prs);
	   free(A->paired);
	   free(A);
	}
}

void	alpha_error(char *s, a_type A) 
{ fprintf(stderr,"%s\n",s); PutAlpha(stderr,A); exit(1); }

