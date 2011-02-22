/* tlb 11/10/06; changes made to NEW functions to comply with ANSII standard;
  changed MIN, MAX to remove conflicts with macros of same names elsewhere
*/
/* defines.h - generic codes and constants for afn biosequence programs. */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include <math.h>

/* VALUES */
#define NIL    		-1
#define FALSE		0	
#define TRUE		1	

#define ILLEGAL		-1.0
#define BELL		((char) 7)	

#define FILE_BEGIN	0
#define FILE_CURRENT	1
#define FILE_END	2

#define Boolean		char

/* CONSTANTS */
#define MAX_INTEGER		LONG_MAX

/* MACROS - standard macro definitions and static types */
#define	MEW(x,n,t)	if ( (x=(t*) malloc(((n)*sizeof(t))))==NULL) { \
			  fprintf(stderr,"Out of Memory."); exit(1); }

#define	NEW(x,n,t)	if ( (x=(t*) calloc(n,sizeof(t)))==NULL ) { \
			  fprintf(stderr,"Out of Memory."); exit(1); }

#define	NEWP(x,n,t)	if ( (x=(t**) calloc(n,sizeof(t*)))==NULL) { \
			  fprintf(stderr,"Out of Memory."); exit(1); }

#define	NEWPP(x,n,t)	(( (x=(t***) calloc(n,sizeof(t**)))==NULL) ? \
			(t***) (fprintf(stderr,"Out of Memory."),exit(1),0):x)

#define	NEWP3(x,n,t)	(( (x=(t****) calloc(n,sizeof(t***)))==NULL) ? \
			(t****) (fprintf(stderr,"Out of Memory."),exit(1),0):x)

#define	GETCHAR(m,C)	do{ fprintf(stderr,"%s ", m); \
			  if(fscanf(stdin,"%c",(C)) == 1) { \
                	    while(getchar()!='\n') if(feof(stdin)) exit(1);\
			    break;\
 			  } while(getchar()!='\n') if(feof(stdin)) exit(1);\
			} while(TRUE);

#define	GETINT(m,i)	do{ fprintf(stderr,"%s ",m); \
			  if(fscanf(stdin,"%d",(i)) == 1) { \
                	    while(getchar()!='\n') if(feof(stdin)) exit(1);\
			    break;\
 			  } while(getchar()!='\n') if(feof(stdin)) exit(1);\
			} while(TRUE);

#define print_error(str) for(fprintf(stderr,"%s\n",str); TRUE; exit(1))
#define DIGIT2INT(c)    ((int)(c - 48))
#define tMIN(t,x,y)	(((t)(x) < (t)(y)) ? (t)(x) : (t)(y))
#define tMAX(t,x,y)	(((t)(x) > (t)(y)) ? (t)(x) : (t)(y))
#define SUM(x)		(((x) * (x+1)) / 2)

