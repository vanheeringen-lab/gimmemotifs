#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gadem.h"
#include "defines.h"

int scan_em_seq_ptable(Pgfs *scoreDist,int pDimension,Sites *site,int numSeq,char **seq,char **rseq,
   int *seqLen,int **ipwm,int pwmLen,int scoreCutoff,double *bfreq,char *Iseq) {

   register int ii,i,m,n;
   int plusScore,minusScore,posOnPlus;
   int siteCn,numUniq,found;
   int *colAve;
   Sites *uniq;

   colAve=alloc_int(pwmLen);
   for (m=0; m<pwmLen; m++) {
      colAve[m]=0; for (i=0; i<4; i++) { colAve[m] +=(int)(ipwm[m][i]*bfreq[i]); } 
   }

   siteCn=0;
   for (ii=0; ii<numSeq; ii++) {
      if (Iseq[ii]=='0') continue;
      /* scan plus strand */
      for (i=0; i<seqLen[ii]-pwmLen+1; i++) {
         posOnPlus=seqLen[ii]-i-pwmLen;
         plusScore=0;
         for (m=0; m<pwmLen; m++) {
            switch(seq[ii][i+m]) {
               case 'a': plusScore +=ipwm[m][0]; break;
               case 'c': plusScore +=ipwm[m][1]; break;
               case 'g': plusScore +=ipwm[m][2]; break;
               case 't': plusScore +=ipwm[m][3]; break;
               default:  plusScore +=colAve[m];  break;
            }
         }

        // scan minus strand
         minusScore=0;
         for (m=0; m<pwmLen; m++) {
            switch(rseq[ii][posOnPlus+m]) {
               case 'a': minusScore +=ipwm[m][0]; break;
               case 'c': minusScore +=ipwm[m][1]; break;
               case 'g': minusScore +=ipwm[m][2]; break;
               case 't': minusScore +=ipwm[m][3]; break;
               default:  minusScore +=colAve[m];  break;
            }
         }
         if (plusScore>=scoreCutoff) {
            if (siteCn>0 && plusScore==site[siteCn-1].score && i-site[siteCn-1].pos<pwmLen && site[siteCn-1].seq==ii) { }
            else {
               site[siteCn].pos=i;
               site[siteCn].seq=ii;
               site[siteCn].rev='0';
               site[siteCn].score=plusScore;
               //site[siteCn].pvalue=find_pvalue(plusScore,scoreDist,pDimension);
               siteCn++;
               if (siteCn==MAX_SITES) {
                  printf("Error: max sites %d reached!\n",MAX_SITES); 
                  printf(" set <MAX_SITES> in defines.h accordingly\n");
                  exit(0); 
               }
            }
         }
         if (minusScore>=scoreCutoff) {
            if (siteCn>0 && minusScore==site[siteCn-1].score && i-site[siteCn-1].pos<pwmLen && site[siteCn-1].seq==ii) { }
            else {
               site[siteCn].pos=posOnPlus;
               site[siteCn].seq=ii;
               site[siteCn].rev='1';
               site[siteCn].score=minusScore;
               //site[siteCn].pvalue=find_pvalue(minusScore,scoreDist,pDimension);
               siteCn++;
               if (siteCn==MAX_SITES) {
                  printf("Error: max sites %d reached!\n",MAX_SITES); 
                  printf(" set <MAX_SITES> in defines.h accordingly\n");
                  exit(0); 
               }
            }
         }
      }
   }
   if (colAve) { free(colAve); colAve=NULL; }
   // printf("number of sites:%d\n",siteCn);

   sort_sites_score(site,siteCn);

   uniq=alloc_site(siteCn);
   numUniq=0;
   for (m=0; m<siteCn; m++) {
      found=0;
      for (n=0; n<numUniq; n++) {
         //---------->------------> or <-----------<------------
         if (site[m].rev==uniq[n].rev && site[m].seq==uniq[n].seq) {
            if (abs(uniq[n].pos-site[m].pos)<pwmLen) { found=1; break; }
         }
         else if (site[m].seq==uniq[n].seq && site[m].rev!=uniq[n].rev) {
            //|----------->w w<-----------|
            if ((seqLen[site[m].seq]-site[m].pos-uniq[n].pos>0 && seqLen[site[m].seq]-site[m].pos-uniq[n].pos<2*pwmLen && site[m].rev=='1')|| 
                (seqLen[uniq[n].seq]-uniq[n].pos-site[m].pos>0 && seqLen[uniq[n].seq]-uniq[n].pos-site[m].pos<2*pwmLen && uniq[n].rev=='1')) {
               found=1; break; 
            }
         }
         else { }
      }
      if (!found) {
         uniq[numUniq].pos   =site[m].pos;
         uniq[numUniq].seq   =site[m].seq;
         uniq[numUniq].rev   =site[m].rev;
         uniq[numUniq].pvalue=find_pvalue(site[m].score,scoreDist,pDimension);;
         numUniq++;
      }
   }
   for (m=0; m<numUniq; m++) {
      site[m].pos   =uniq[m].pos;
      site[m].seq   =uniq[m].seq;
      site[m].rev   =uniq[m].rev;
      site[m].pvalue=uniq[m].pvalue;
   }
   // printf("uniq sites: %d\n",numUniq);

   if (uniq) { free(uniq); uniq=NULL; }
   return (numUniq);
}

int scan_llr_pgf(Pgfs *scoreDist,int pDimension,Sites *site,int numSeq,char **seq,char **rseq,
   int *seqLen,int **ipwm,int pwmLen,int scoreCutoff,double *bfreq) {

   register int ii,i,m,n;
   int plusScore,minusScore,posOnPlus;
   int siteCn,numUniq,found;
   int *colAve;
   Sites *uniq;

   colAve=alloc_int(pwmLen);
   for (m=0; m<pwmLen; m++) {
      colAve[m]=0; for (i=0; i<4; i++) { colAve[m] +=(int)(ipwm[m][i]*bfreq[i]); }
   }
   
   siteCn=0;
   for (ii=0; ii<numSeq; ii++) {
      /* scan plus strand */
      for (i=0; i<seqLen[ii]-pwmLen+1; i++) {
         posOnPlus=seqLen[ii]-i-pwmLen;
         plusScore=0;
         for (m=0; m<pwmLen; m++) {
            switch(seq[ii][i+m]) {
               case 'a': plusScore +=ipwm[m][0]; break;
               case 'c': plusScore +=ipwm[m][1]; break;
               case 'g': plusScore +=ipwm[m][2]; break;
               case 't': plusScore +=ipwm[m][3]; break;
               default:  plusScore +=colAve[m];  break;
            }
         }

         minusScore=0;
         for (m=0; m<pwmLen; m++) {
            switch(rseq[ii][posOnPlus+m]) {
               case 'a': minusScore +=ipwm[m][0]; break;
               case 'c': minusScore +=ipwm[m][1]; break;
               case 'g': minusScore +=ipwm[m][2]; break;
               case 't': minusScore +=ipwm[m][3]; break;
               default:  minusScore +=colAve[m];  break;
            }
         }
         if (plusScore>=scoreCutoff) {
            if (siteCn>0 && plusScore==site[siteCn-1].score && i-site[siteCn-1].pos<pwmLen && site[siteCn-1].seq==ii) { }
            else {
               site[siteCn].pos=i;
               site[siteCn].seq=ii;
               site[siteCn].rev='0';
               site[siteCn].score=plusScore;
               siteCn++;
               if (siteCn==MAX_SITES) {
                  printf("Error: max sites %d reached!\n",MAX_SITES);
                  printf(" set <MAX_SITES> in defines.h accordingly\n");
                  exit(0);
               }
            }
         }
         if (minusScore>=scoreCutoff) {
            if (siteCn>0 && minusScore==site[siteCn-1].score && i-site[siteCn-1].pos<pwmLen && site[siteCn-1].seq==ii) { }
            else {
               site[siteCn].pos=posOnPlus;
               site[siteCn].seq=ii;
               site[siteCn].rev='1';
               site[siteCn].score=minusScore;
               siteCn++;
               if (siteCn==MAX_SITES) {
                  printf("Error: max sites %d reached!\n",MAX_SITES);
                  printf(" set <MAX_SITES> in defines.h accordingly\n");
                  exit(0);
               }
            }
         }
      }
   }
   if (colAve) { free(colAve); colAve=NULL; }

   sort_sites_score(site,siteCn);

   uniq=alloc_site(siteCn);
   numUniq=0;
   for (m=0; m<siteCn; m++) {
      found=0;
      for (n=0; n<numUniq; n++) {
         //---------->------------> or <-----------<------------
         if (site[m].rev==uniq[n].rev && site[m].seq==uniq[n].seq) {
            if (abs(uniq[n].pos-site[m].pos)<pwmLen) { found=1; break; }
         }
         else if (site[m].seq==uniq[n].seq && site[m].rev!=uniq[n].rev) {
            //|----------->w w<-----------|
            if ((seqLen[site[m].seq]-site[m].pos-uniq[n].pos>0 && seqLen[site[m].seq]-site[m].pos-uniq[n].pos<2*pwmLen && site[m].rev=='1')|| 
                (seqLen[uniq[n].seq]-uniq[n].pos-site[m].pos>0 && seqLen[uniq[n].seq]-uniq[n].pos-site[m].pos<2*pwmLen && uniq[n].rev=='1')) {
               found=1; break; 
            }
         }
         else { }
      }
      if (!found) {
         uniq[numUniq].pos   =site[m].pos;
         uniq[numUniq].seq   =site[m].seq;
         uniq[numUniq].rev   =site[m].rev;
         uniq[numUniq].score =site[m].score;
         numUniq++;
      }
   }
   for (m=0; m<numUniq; m++) {
      site[m].pos   =uniq[m].pos;
      site[m].seq   =uniq[m].seq;
      site[m].rev   =uniq[m].rev;
      site[m].score =uniq[m].score;
      site[m].pvalue=find_pvalue(site[m].score,scoreDist,pDimension);
   }

   if (uniq) { free(uniq); uniq=NULL; }
   return (numUniq);
}

int scan_llr_empirical(Sites *site,int numSeq,char **seq,char **rseq,int *seqLen,double **logepwm,int pwmLen,
   double *bfreq,double **bscore,double **rbscore,double llrCutoff,double *nullScore,int numTopKmer,int totalKmer) {

   register int ii,i,m,n;
   int posOnPlus;
   int siteCn,numUniq,found;
   double plusScore,minusScore;
   double *colAve;
   Sites *uniq;

   colAve=alloc_double(pwmLen);
   for (m=0; m<pwmLen; m++) {
      colAve[m]=0; for (i=0; i<4; i++) { colAve[m] +=(logepwm[m][i]*bfreq[i]); } 
   }

   siteCn=0;
   for (ii=0; ii<numSeq; ii++) {
      /* scan plus strand */
      for (i=0; i<seqLen[ii]-pwmLen+1; i++) {
         posOnPlus=seqLen[ii]-i-pwmLen;
         plusScore=0;
         for (m=0; m<pwmLen; m++) {
            switch(seq[ii][i+m]) {
               case 'a': plusScore +=logepwm[m][0]; break;
               case 'c': plusScore +=logepwm[m][1]; break;
               case 'g': plusScore +=logepwm[m][2]; break;
               case 't': plusScore +=logepwm[m][3]; break;
               default:  plusScore +=colAve[m];  break;
            }
         }
         plusScore -=bscore[ii][i];

        // scan minus strand
         minusScore=0;
         for (m=0; m<pwmLen; m++) {
            switch(rseq[ii][posOnPlus+m]) {
               case 'a': minusScore +=logepwm[m][0]; break;
               case 'c': minusScore +=logepwm[m][1]; break;
               case 'g': minusScore +=logepwm[m][2]; break;
               case 't': minusScore +=logepwm[m][3]; break;
               default:  minusScore +=colAve[m];  break;
            }
         }
         minusScore -=rbscore[ii][posOnPlus];

         if (plusScore>=llrCutoff) {
            if (siteCn>0 && plusScore==site[siteCn-1].score && i-site[siteCn-1].pos<pwmLen && site[siteCn-1].seq==ii) { }
            else {
               site[siteCn].pos=i;
               site[siteCn].seq=ii;
               site[siteCn].llr=plusScore;
               site[siteCn].rev='0';
               siteCn++;
               if (siteCn==MAX_SITES) {
                  printf("Error: max sites %d reached!\n",MAX_SITES); 
                  printf(" set <MAX_SITES> in defines.h accordingly\n");
                  exit(0); 
               }
            }
         }

         if (minusScore>=llrCutoff) {
            if (siteCn>0 && minusScore==site[siteCn-1].score && i-site[siteCn-1].pos<pwmLen && site[siteCn-1].seq==ii) { }
            else {
               site[siteCn].pos=posOnPlus;
               site[siteCn].seq=ii;
               site[siteCn].llr=minusScore;
               site[siteCn].rev='1';
               siteCn++;
               if (siteCn==MAX_SITES) {
                  printf("Error: max sites %d reached!\n",MAX_SITES); 
                  printf(" set <MAX_SITES> in defines.h accordingly\n");
                  exit(0); 
               }
            }
         }
      }
   }
   if (colAve) { free(colAve); colAve=NULL; }
   printf("number of sites:%d\n",siteCn);

   sort_sites_llr(site,siteCn); // sort by llr in descending order

   uniq=alloc_site(siteCn);
   numUniq=0;
   for (m=0; m<siteCn; m++) {
      found=0;
      for (n=0; n<numUniq; n++) {
         //---------->------------> or <-----------<------------
         if (site[m].rev==uniq[n].rev && site[m].seq==uniq[n].seq) {
            if (abs(uniq[n].pos-site[m].pos)<pwmLen) { found=1; break; }
         }
         else if (site[m].seq==uniq[n].seq && site[m].rev !=uniq[n].rev) {
            //|----------->w w<-----------|
            if (((seqLen[site[m].seq]-site[m].pos-uniq[n].pos)>0 && (seqLen[site[m].seq]-site[m].pos-uniq[n].pos)<2*pwmLen && site[m].rev=='1')|| 
                ((seqLen[uniq[n].seq]-uniq[n].pos-site[m].pos)>0 && (seqLen[uniq[n].seq]-uniq[n].pos-site[m].pos)<2*pwmLen && uniq[n].rev=='1')) {
               found=1; break; 
            }
         }
         else { }
         // if (uniq[n].seq==site[m].seq && uniq[n].rev!=site[m].rev && abs(seqLen[site[m].seq]-site[m].pos-pwmLen-uniq[n].pos)<pwmLen) { 
      }
      if (!found) {
         uniq[numUniq].pos = site[m].pos;
         uniq[numUniq].seq = site[m].seq;
         uniq[numUniq].rev = site[m].rev;
         uniq[numUniq].llr = site[m].llr;
         numUniq++;
      }
   }
   for (m=0; m<numUniq; m++) {
      site[m].pos =uniq[m].pos;
      site[m].seq =uniq[m].seq;
      site[m].rev =uniq[m].rev;
      site[m].llr =uniq[m].llr;

      site[m].pvalue=0;
      for (n=0; n<numTopKmer; n++) {
         if (site[m].llr>nullScore[n]) {
            site[m].pvalue=(double)n/(double)totalKmer; 
            break; 
         }
      }
   }
   printf("uniq sites: %d\n",numUniq);

   if (uniq) { free(uniq); uniq=NULL; }
   return (numUniq);
}
