#include "afnio.h"

FILE	*open_file(char *fstring,char *subfile,char *cmnd)
{
	FILE	*fptr;
	char	s[100];

	while(fstring[0] == ' ') fstring++;
	strcpy(s,fstring);
	strcat(s,subfile);
	if((fptr = fopen(s,cmnd)) == NULL) {
		fprintf(stderr,"Could not open file \"%s\"\n",s);
		print_error("File does not exist!\n");
	}
	return(fptr);
}

long     ParseReals(char *str, double *values, char *msg)
/** WARNING: index starts at 0 not 1 as for ParseIntegers() **/
{
        long	n;
	double	k;

	if(!isdigit(str[0])) print_error(msg);
        for(n=1; str[0] != 0; ){
           if(str[0] == ',') { n++; str++; }
           else if(isdigit(str[0])){
                if(sscanf(str,"%lf", &k) != 1){ print_error(msg); }
                else { 
		   values[n-1] = k; 
		   while(isdigit(str[0]) || str[0] == '.' ) str++; 
		}
           } else print_error(msg);
        }
        return n;
}

long     ParseIntegers(char *str, long *values, char *msg)
{
        long     n,k;

	if(!isdigit(str[0])) print_error(msg);
        for(n=1; str[0] != 0; ){
           if(str[0] == ',') { n++; str++; }
           else if(isdigit(str[0])){
                if(sscanf(str,"%ld", &k) != 1){ print_error(msg); }
                else { values[n] = k; while(isdigit(str[0])) str++; }
           } else print_error(msg);
        }
        return n;
}

char	*String(char *s)
{
	char *S;
	if((S=(char*) calloc((strlen(s)+1),sizeof(char)))==NULL)
		print_error("Out of Memory");
	strcpy(S,s); return S; 
}

