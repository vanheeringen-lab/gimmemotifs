#include "purge.h"

int	main(int argc,char *argv[])
{ 
	int	i,cutoff_score = 200,minimum=0,arg;
	char	score[100],*DBS_NAME,c=' ';
	ss_type	P = NULL;
	a_type	A;
	Boolean	success=FALSE,DNA=FALSE,xnu=FALSE,query=FALSE;

	if(argc < 3) print_error(PURGE_USAGE);
	DBS_NAME = argv[1];
	if(sscanf(argv[2],"%d",&i)!=1) print_error(PURGE_USAGE);
	cutoff_score = tMAX(int,1,i);
	sprintf(score,".%c%d",c,cutoff_score);
        for(arg = 3; arg < argc; arg++){
           if(argv[arg][0] != '-') print_error(PURGE_USAGE);
           switch(argv[arg][1]) {
                case 'b': c = 'b'; break;
                case 'e': c = 'e'; break;
                case 'f': c = 'f'; break;
                case 'q': query=TRUE; break;
                case 'x': xnu=TRUE; break;
                case 'n': DNA=TRUE; xnu=FALSE; break;
                default: print_error(PURGE_USAGE);
            }
        }
	if(c=='b' && DNA == TRUE)
		print_error("-b option is not allowed with -n option");
	if(xnu==TRUE && DNA == TRUE)
		print_error("-x option is not allowed with -n option");
	if(DNA) A = MkAlpha("NACGT",DNA_MTRX);
	else A = MkAlpha(AMINO_ACIDS,PROT_BLOSUM62);
	if(xnu) P = MkXnuSeqSet(DBS_NAME,A);
	else P = SeqSet(DBS_NAME,A);
        if (c==' ') { if (DNA) c = 'e'; else c = 'b'; }
	success = RmHomologs(cutoff_score, c, minimum, query, P);
	NilAlpha(A); NilSeqSet(P);
	if(success) return 0;
	else return 1;
}

