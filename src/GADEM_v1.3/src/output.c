#include <stdio.h>
#include <stdlib.h>
#include "stdlib.h"
#include <string.h>
#include <math.h>
#include "defines.h"
#include "gadem.h"

void print_motif(Sites *site,int nsites,char **seq,char **rseq,int *seqLen,int pwmLen,int id,double **opwm) {

   FILE *f1;
   char *fileName;
   register int i,j;

   fileName=alloc_char(500);
   sprintf(fileName,"%d.seq",id);
   f1=fopen(fileName,"w");
   for (i=0; i<nsites; i++) {
      if (site[i].rev=='0') {
         if (site[i].pos<0) {
            for (j=site[i].pos; j<0; j++) fprintf(f1,"x"); 
            for (j=0; j<pwmLen+site[i].pos; j++) {
               switch(seq[site[i].seq][j]) {
                  case 'a': fprintf(f1,"a"); break;
                  case 'c': fprintf(f1,"c"); break;
                  case 'g': fprintf(f1,"g"); break;
                  case 't': fprintf(f1,"t"); break;
                  case 'n': fprintf(f1,"n"); break;
                  default: break;
               }
            }
         }
         else {
            for (j=site[i].pos; j<min(seqLen[site[i].seq],site[i].pos+pwmLen); j++) {
               switch(seq[site[i].seq][j]) {
                  case 'a': fprintf(f1,"a"); break;
                  case 'c': fprintf(f1,"c"); break;
                  case 'g': fprintf(f1,"g"); break;
                  case 't': fprintf(f1,"t"); break;
                  case 'n': fprintf(f1,"n"); break;
                  default: break;
               }
            }
         }
         if (site[i].pos+pwmLen-seqLen[site[i].seq]>0) {
            for (j=seqLen[site[i].seq]; j<site[i].pos+pwmLen; j++) fprintf(f1,"x"); 
         }
         fprintf(f1,"\n");
      }
      else {
         if (site[i].pos<0) {
            for (j=site[i].pos; j<0; j++) fprintf(f1,"x"); 
            for (j=0; j<pwmLen+site[i].pos; j++) {
               switch(rseq[site[i].seq][j]) {
                  case 'a': fprintf(f1,"a"); break;
                  case 'c': fprintf(f1,"c"); break;
                  case 'g': fprintf(f1,"g"); break;
                  case 't': fprintf(f1,"t"); break;
                  case 'n': fprintf(f1,"n"); break;
                  default: break;
               }
            }
         }
         else {
            for (j=site[i].pos; j<min(seqLen[site[i].seq],site[i].pos+pwmLen); j++) {
               switch(rseq[site[i].seq][j]) {
                  case 'a': fprintf(f1,"a"); break;
                  case 'c': fprintf(f1,"c"); break;
                  case 'g': fprintf(f1,"g"); break;
                  case 't': fprintf(f1,"t"); break;
                  case 'n': fprintf(f1,"n"); break;
                  default: break;
               }
            }
         }
         if (site[i].pos+pwmLen-seqLen[site[i].seq]>0) {
            for (j=seqLen[site[i].seq]; j<site[i].pos+pwmLen; j++) fprintf(f1,"x"); 
         }
         fprintf(f1,"\n");
      }
   }
   fclose(f1);
   if (fileName) { free(fileName); fileName=NULL; }

   // print out individual observed PWM in gadem format
   /*-----------------------------------------------------------------------
      fileName=alloc_char(500);
      sprintf(fileName,"%d.mx",id);
      f1=fopen(fileName,"w");

      fprintf(f1,"4\t%d\n",pwmLen);
      for (i=0; i<4; i++) {
         for (j=0; j<pwmLen; j++) {
            if (j<pwmLen-1) fprintf(f1,"%5.4f\t",opwm[j][i]);
            else fprintf(f1,"%5.4f\n",opwm[j][i]);
         }
      }
      fclose(f1);
      if (fileName) { free(fileName); fileName=NULL; }
   }
   -----------------------------------------------------------------------*/
}

void print_result_2(Sites *site,int nsites,int numSeq,char **seq,char **rseq,int *seqLen,char **geneID,
   double logev,double **opwm,int pwmLen,int id,char *sdyad,char *pwmConsensus,int numCycle,
   double pvalueCutoff,double maxpFactor,FILE *fq,FILE *fpwm) {

   register int i,j;
   int cn[4],maxHeaderLen;
   int *seqCn;

   seqCn=alloc_int(numSeq);

   maxHeaderLen=strlen(geneID[0]);
   for (i=1; i<numSeq; i++) {
      if (strlen(geneID[i])>maxHeaderLen) maxHeaderLen=strlen(geneID[i]); 
   }
   maxHeaderLen=min(maxHeaderLen,MAX_SEQ_HEADER);

   for (i=0; i<numSeq; i++) seqCn[i]=0;
   for (i=0; i<nsites; i++) seqCn[site[i].seq]++; 
  
   for (i=0; i<4; i++) cn[i]=0; 
   for (i=0; i<numSeq; i++) {
      if (seqCn[i]==0) cn[0]++; 
      if (seqCn[i]==1) cn[1]++; 
      if (seqCn[i]==2) cn[2]++; 
      if (seqCn[i]>2)  cn[3]++; 
   }
   if (seqCn) { free(seqCn); seqCn=NULL; }

   fprintf(fq,"Cycle[%3d] motif[%d]:\n",numCycle,id);
   fprintf(fq,"   spaced dyad:\t\t\t%s\n",sdyad);
   fprintf(fq,"   motif consensus:\t\t%s\n",pwmConsensus);
   fprintf(fq,"\t\t\t\tm=ac r=ag w=at s=cg y=ct k=gt b=cgt d=agt h=act v=acg\n");
   fprintf(fq,"   motif length(w):\t\t%d\n",pwmLen);
   fprintf(fq,"   maxpfactor:\t\t\t%5.4f\n",maxpFactor);
   fprintf(fq,"   number of sites:\t\t%d\n",nsites);
   fprintf(fq,"   ln(E-value):\t\t\t%9.2f\n",logev);
   fprintf(fq,"   pwm p-value cutoff:\t\t%e\n\n",pvalueCutoff);
   fprintf(fq,"   Seqs with 0,1,2,>2 sites: ");
   fprintf(fq,"   %d,%d,%d,%d\n",cn[0],cn[1],cn[2],cn[3]);
   fprintf(fq,"   %d(%4.2f%) of %d seqs have >=1 site\n",
                  cn[1]+cn[2]+cn[3],100.0*(double)(cn[1]+cn[2]+cn[3])/(double)numSeq,numSeq);
   fprintf(fq,"\nIndividual sites: [%2dbp flanking--MOTIF--%2dbp flanking] [strand] [seq] [pos] [p-value]\n",FLANKING_BASES,FLANKING_BASES);
   for (i=0; i<nsites; i++) {
      for (j=0; j<min(maxHeaderLen,strlen(geneID[site[i].seq])); j++)  fprintf(fq,"%c",geneID[site[i].seq][j]); 
      for (j=strlen(geneID[site[i].seq]); j<maxHeaderLen; j++) fprintf(fq," ");
      fprintf(fq,"\t");
 
      if (site[i].rev=='0') {
         if (site[i].pos<FLANKING_BASES) {
            if (site[i].pos<0) {
               for (j=0; j<FLANKING_BASES; j++) fprintf(fq," "); 
            }
            else {
               for (j=0; j<FLANKING_BASES-site[i].pos; j++) fprintf(fq," "); 
               for (j=0; j<site[i].pos; j++) fprintf(fq,"%c",seq[site[i].seq][j]);
            }
         }
         else {
            for (j=site[i].pos-FLANKING_BASES; j<site[i].pos; j++) fprintf(fq,"%c",seq[site[i].seq][j]);
         }

         if (site[i].pos<0) {
            for (j=site[i].pos; j<0; j++) fprintf(fq,"X"); 
            for (j=0; j<pwmLen+site[i].pos; j++) {
               switch(seq[site[i].seq][j]) {
                  case 'a': fprintf(fq,"A"); break;
                  case 'c': fprintf(fq,"C"); break;
                  case 'g': fprintf(fq,"G"); break;
                  case 't': fprintf(fq,"T"); break;
                  case 'n': fprintf(fq,"N"); break;
                  default: break;
               }
            }
         }
         else {
            for (j=site[i].pos; j<min(seqLen[site[i].seq],site[i].pos+pwmLen); j++) {
               switch(seq[site[i].seq][j]) {
                  case 'a': fprintf(fq,"A"); break;
                  case 'c': fprintf(fq,"C"); break;
                  case 'g': fprintf(fq,"G"); break;
                  case 't': fprintf(fq,"T"); break;
                  case 'n': fprintf(fq,"N"); break;
                  default: break;
               }
            }
         }
         if (site[i].pos+pwmLen-seqLen[site[i].seq]>0) {
            for (j=seqLen[site[i].seq]; j<site[i].pos+pwmLen; j++) fprintf(fq,"X"); 
         }
         // print flanking region
         for (j=site[i].pos+pwmLen; j<min(site[i].pos+pwmLen+FLANKING_BASES,seqLen[site[i].seq]); j++) 
            fprintf(fq,"%c",seq[site[i].seq][j]);
         if (site[i].pos+pwmLen+FLANKING_BASES>seqLen[site[i].seq]) { 
            for (j=0; j<site[i].pos+pwmLen+FLANKING_BASES-seqLen[site[i].seq]; j++) fprintf(fq," "); 
         }
         fprintf(fq,"\t+\t%d\t%d\t%e\n",site[i].seq+1,site[i].pos+1,site[i].pvalue);
      }
      else {
         if (site[i].pos<FLANKING_BASES) {
            if (site[i].pos<0) {
               for (j=0; j<FLANKING_BASES; j++) fprintf(fq," "); 
            }
            else {
               for (j=0; j<FLANKING_BASES-site[i].pos; j++) fprintf(fq," "); 
               for (j=0; j<site[i].pos; j++) fprintf(fq,"%c",rseq[site[i].seq][j]);
            }
         }
         else {
            for (j=site[i].pos-FLANKING_BASES; j<site[i].pos; j++) fprintf(fq,"%c",rseq[site[i].seq][j]);
         }

         if (site[i].pos<0) {
            for (j=site[i].pos; j<0; j++) fprintf(fq,"X"); 
            for (j=0; j<pwmLen+site[i].pos; j++) {
               switch(rseq[site[i].seq][j]) {
                  case 'a': fprintf(fq,"A"); break;
                  case 'c': fprintf(fq,"C"); break;
                  case 'g': fprintf(fq,"G"); break;
                  case 't': fprintf(fq,"T"); break;
                  case 'n': fprintf(fq,"N"); break;
                  default: break;
               }
            }
         }
         else {
            for (j=site[i].pos; j<min(seqLen[site[i].seq],site[i].pos+pwmLen); j++) {
               switch(rseq[site[i].seq][j]) {
                  case 'a': fprintf(fq,"A"); break;
                  case 'c': fprintf(fq,"C"); break;
                  case 'g': fprintf(fq,"G"); break;
                  case 't': fprintf(fq,"T"); break;
                  case 'n': fprintf(fq,"N"); break;
                  default: break;
               }
            }
         }
         if (site[i].pos+pwmLen-seqLen[site[i].seq]>0) {
            for (j=seqLen[site[i].seq]; j<site[i].pos+pwmLen; j++) fprintf(fq,"X"); 
         }
         // print flanking region
         for (j=site[i].pos+pwmLen; j<min(site[i].pos+pwmLen+FLANKING_BASES,seqLen[site[i].seq]); j++) 
            fprintf(fq,"%c",rseq[site[i].seq][j]);
         if (site[i].pos+pwmLen+FLANKING_BASES>seqLen[site[i].seq]) { 
            for (j=0; j<site[i].pos+pwmLen+FLANKING_BASES-seqLen[site[i].seq]; j++) fprintf(fq," "); 
         }
         fprintf(fq,"\t-\t%d\t%d\t%e\n",site[i].seq+1,seqLen[site[i].seq]-site[i].pos,site[i].pvalue);
      }
   }
   fprintf(fq,"\n");
   fprintf(fpwm,">m%d_c%d_%s\n",id,numCycle,pwmConsensus);
   for (i=0; i<4; i++) {
      switch (i) {
         case 0: fprintf(fpwm,"A "); break; 
         case 1: fprintf(fpwm,"C "); break; 
         case 2: fprintf(fpwm,"G "); break; 
         case 3: fprintf(fpwm,"T "); break; 
         default: break; 
      }
      for (j=0; j<pwmLen; j++) {
         if (j<pwmLen-1) fprintf(fpwm,"%5.3f\t",opwm[j][i]);
         else fprintf(fpwm,"%5.3f\n",opwm[j][i]);
      }
   }
   fflush(fq);
   fflush(fpwm);
}

