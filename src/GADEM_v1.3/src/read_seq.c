#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include "defines.h"

char *alloc_char(int );
char **alloc_char_char(int, int);
double *alloc_double(int );

char **read_seq(int *numSeq,int *seqLen,char **geneID,int maxNumSeq,int maxSeqLen,double *chipScore,char *fileName) {

   FILE *fp;
   register int i,j,k;
   int len,num_seq,cn,maxBufferLen,pos,dotCn,numCn,yesSpace;
   char  **seq,*buffer,*tok,*sscore;

   fp=fopen(fileName,"r");
   if (!fp) { perror(fileName); exit(0); }

   maxBufferLen=MAX_BUFFER_LENGTH;
   seq=alloc_char_char(maxNumSeq,maxSeqLen+1);
   buffer=alloc_char(maxBufferLen);
   sscore=alloc_char(100);

   i=0;
   if (fgets(buffer,maxBufferLen,fp)>0) {
      while (!feof(fp)){
         len=strlen(buffer);
         buffer[len]='\0';
         if (buffer[0]=='>') {
            yesSpace=0;
            for (j=1; j<len; j++) {
               if (buffer[j]==' ') { yesSpace=1; break; } 
            }
            if (yesSpace) tok=strtok(buffer," ");
            else          tok=strtok(buffer,"\n");
            strcpy(geneID[i],tok);
            len=strlen(tok);
            geneID[i][len]='\0';
  
            pos=0;
            for (j=0; j<len; j++) {
               if (geneID[i][j]=='_') pos=j+1; 
            }
            for (j=0; j<len-7; j++) {
               if ( geneID[i][j]  =='_' && 
                   (geneID[i][j+1]=='S' || geneID[i][j+1]=='s') && 
                   (geneID[i][j+2]=='C' || geneID[i][j+2]=='c') &&
                   (geneID[i][j+3]=='O' || geneID[i][j+3]=='o') && 
                   (geneID[i][j+4]=='R' || geneID[i][j+4]=='r') && 
                   (geneID[i][j+5]=='E' || geneID[i][j+5]=='e') && 
                   (geneID[i][j+6]=='=') 
                  ) {
                  pos=j+1+6;
               } 
            }

            if (pos==0) chipScore[i]=0; 
            else {
               numCn=0;
               for (k=0,j=pos; j<len; j++,k++) {
                  sscore[k]=geneID[i][j]; 
                  if (isdigit(geneID[i][j])) numCn++; 
               }
               sscore[k]='\0';

               dotCn=0;
               for (j=pos; j<len; j++) {
                  if (geneID[i][j]=='.') dotCn++; 
               }
               if (numCn+dotCn==len-pos) chipScore[i]=strtod(sscore,NULL); 
               else chipScore[i]=0;
            }

            cn=0;
            do {
               if (fgets(buffer,maxBufferLen,fp)) {
                  len=strlen(buffer);
                  buffer[len-1]='\0';
                  if (buffer[0]!='>') {
                     for (j=0; j<len-1; j++) {
                        if (cn<maxSeqLen) { 
                           seq[i][cn]=buffer[j]; cn++; 
                        }
                     }
                  }
                  else break;
               }
               else break;
            } while (buffer[0]!='>');
            seq[i][cn]='\0'; seqLen[i]=cn;
            if (cn>1) {
               i++; 
               if (i>=maxNumSeq) { 
                  printf("\n\nErro: maximal number of seqences reached!\n"); 
                  printf("Please reset MAX_NUM_SEQ in gadem.h and rebuild (see installation)\n\n");
                  exit(0); 
               }
            }
         }
      };
   }
   fclose(fp);
   if (buffer) { free(buffer); buffer=NULL; }

   num_seq=i;
   for (i=0; i<num_seq; i++) {
      for (j=0; j<seqLen[i]; j++) {
         switch(seq[i][j]) {
            case 'A': seq[i][j]='a'; break;
            case 'C': seq[i][j]='c'; break;
            case 'G': seq[i][j]='g'; break;
            case 'T': seq[i][j]='t'; break;
            case 'N': seq[i][j]='n'; break;
            case 'a': seq[i][j]='a'; break;
            case 'c': seq[i][j]='c'; break;
            case 'g': seq[i][j]='g'; break;
            case 't': seq[i][j]='t'; break;
            case 'n': seq[i][j]='n'; break;
            default:  seq[i][j]='n'; break;
         }
      }
   }
   *numSeq=num_seq;
   //for (i=0; i<num_seq; i++) printf("seq[%4d] chipScore: %5.2f\n",i+1,chipScore[i]);
   //for (i=0; i<2; i++) printf("seq[%4d] length: %5d\n",i+1,seqLen[i]);
   /*------------------------------------------------------------------------
   fp=fopen("debug.seq","w");
   for (i=0; i<num_seq; i++) {
      fprintf(fp,"%s\n",geneID[i]);
      fprintf(fp,"%s\n",seq[i]);
   }
   fclose(fp);
   for (i=0; i<num_seq; i++) printf("seq[%4d] length: %5d\n",i+1,seqLen[i]);
   ------------------------------------------------------------------------*/
   if (sscore) { free(sscore); sscore=NULL; }

   return (seq);
}

