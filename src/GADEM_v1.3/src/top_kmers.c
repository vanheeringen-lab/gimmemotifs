#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include "gadem.h"
#include "defines.h"

int word_for_dyad(Words *word,char **seq,char **rseq,int numSeq,int *seqLen,double *baseFreq,
   int *numTop3mer,int *numTop4mer,int *numTop5mer) {

   int numTrimer,numTetramer,numPentamer;
   int numGroup,numAvailableTopKmer;
   int *trimerCn,*tetramerCn,*pentamerCn;
   char **trimer, **tetramer,**pentamer;
   Ktuples *kp3,*kp4,*kp5;

#ifdef DEBUG
   FILE *fp;
   int i,j;
#endif

   numTrimer=64; numTetramer=256; numPentamer=1024;

   // memory allocations
   trimer    =alloc_char_char(numTrimer,5);
   tetramer  =alloc_char_char(numTetramer,6);
   pentamer  =alloc_char_char(numPentamer,7);
   trimerCn  =alloc_int(numTrimer);
   tetramerCn=alloc_int(numTetramer);
   pentamerCn=alloc_int(numPentamer);

   // enumerate
   enumerate_kmers(trimer,tetramer,pentamer);

   count_k_tuples(seq,rseq,numSeq,seqLen,trimer,  numTrimer,  3,trimerCn);
   count_k_tuples(seq,rseq,numSeq,seqLen,tetramer,numTetramer,4,tetramerCn);
   count_k_tuples(seq,rseq,numSeq,seqLen,pentamer,numPentamer,5,pentamerCn);

   kp3=alloc_ktuples(numTrimer,  4);
   kp4=alloc_ktuples(numTetramer,5);
   kp5=alloc_ktuples(numPentamer,6);

   score_kmers(kp3,baseFreq,numTrimer,  trimerCn,  trimer,  3);
   score_kmers(kp4,baseFreq,numTetramer,tetramerCn,tetramer,4);
   score_kmers(kp5,baseFreq,numPentamer,pentamerCn,pentamer,5);

   //sort by z-score
   sort_kmer_z(kp3,numTrimer);
   sort_kmer_z(kp4,numTetramer);
   sort_kmer_z(kp5,numPentamer);

   numGroup=0;
   if (*numTop3mer!=0 ) {
      numAvailableTopKmer=top_kmer(&(word[numGroup]),kp3,numTop3mer,3);
      if (numAvailableTopKmer!=0) {
         numGroup++; *numTop3mer=numAvailableTopKmer; 
      }
      else *numTop3mer=0;
   }
   if (*numTop4mer!=0) {
      numAvailableTopKmer=top_kmer(&(word[numGroup]),kp4,numTop4mer,4);
      if (numAvailableTopKmer!=0) {
         numGroup++; *numTop4mer=numAvailableTopKmer; 
      }
      else *numTop4mer=0;
   }
   if (*numTop5mer!=0) {
      numAvailableTopKmer=top_kmer(&(word[numGroup]),kp5,numTop5mer,5);
      if (numAvailableTopKmer!=0) {
         numGroup++; *numTop5mer=numAvailableTopKmer; 
      }
      else *numTop5mer=0;
   }
   printf("top 3  4, 5-mers: %d %d %d\n",*numTop3mer,*numTop4mer,*numTop5mer);

#ifdef DEBUG
   fp=fopen("kmer.debug","w");
   for (i=0; i<numTrimer; i++) {
      fprintf(fp,"%s\t%d\t%5.1f\t%5.3f\n",kp3[i].seq,kp3[i].count,kp3[i].z,kp3[i].p);
   }
   for (i=0; i<numTetramer; i++) {
      fprintf(fp,"%s\t%d\t%5.1f\t%5.3f\n",kp4[i].seq,kp4[i].count,kp4[i].z,kp4[i].p);
   }
   for (i=0; i<numPentamer; i++) {
      fprintf(fp,"%s\t%d\t%5.1f\t%5.3f\n",kp5[i].seq,kp5[i].count,kp5[i].z,kp5[i].p);
   }
   fclose(fp);
#endif

#ifdef DEBUG
   fp=fopen("topWords.debug","w");
   for (i=0; i<numGroup; i++) {
      fprintf(fp,"%d\n",word[i].count);
      for (j=0; j<word[i].count; j++) fprintf(fp,"%s\t%8.6f\t%8.6f\n",word[i].s1[j],word[i].prob_sta[j],word[i].prob_end[j]);
   }
   fclose(fp);
#endif

   if (trimer[0])   { free(trimer[0]);   trimer[0]=NULL;   }
   if (trimer)      { free(trimer);      trimer=NULL;      }
   if (tetramer[0]) { free(tetramer[0]); tetramer[0]=NULL; }
   if (tetramer)    { free(tetramer);    tetramer=NULL;    }
   if (pentamer[0]) { free(pentamer[0]); pentamer[0]=NULL; }
   if (pentamer)    { free(pentamer);    pentamer=NULL;    }
   if (trimerCn)    { free(trimerCn);    trimerCn=NULL;    }
   if (tetramerCn)  { free(tetramerCn);  tetramerCn=NULL;  }
   if (pentamerCn)  { free(pentamerCn);  pentamerCn=NULL;  }

   if (kp3) destroy_ktuples(kp3,numTrimer);
   if (kp4) destroy_ktuples(kp4,numTetramer);
   if (kp5) destroy_ktuples(kp5,numPentamer);
 
   return (numGroup);
}

void score_kmers(Ktuples *kp,double *baseFreq,int numKmer,int *kmerCn,char **kmer,int kmerLen) {

   register int i,j;
   int totalCn;
   double p;

   totalCn=0; for (i=0; i<numKmer; i++)  totalCn +=kmerCn[i];

   for (i=0; i<numKmer; i++) {
      strcpy(kp[i].seq,kmer[i]); kp[i].seq[kmerLen]='\0';
      p=1.0;
      for (j=0; j<kmerLen; j++) {
         switch (kmer[i][j]) {
            case 'a': p *=baseFreq[0]; break; 
            case 'c': p *=baseFreq[1]; break; 
            case 'g': p *=baseFreq[2]; break; 
            case 't': p *=baseFreq[3]; break; 
            default: break; 
         }
      } 
      kp[i].count=kmerCn[i];
      if (p!=0.0) {
         kp[i].z=((double)kmerCn[i]-(double)totalCn*p)/sqrt((double)totalCn*p*(1.0-p));
      }
      else kp[i].z=0.0;
      kp[i].eCn=(double)totalCn*p;
   }
}

int top_kmer(Words *word,Ktuples *kp,int *numSpecifiedTopKmer,int kmerLen) {

   register int i;
   int numAvailableTopKmer;
   double sum;

   // standardize z scores for top-ranked tetramers
   sum=0.0; 
   for (i=0; i<*numSpecifiedTopKmer; i++) {
      if (kp[i].z>=Z_CUTOFF) sum +=kp[i].z; 
      else break;
   }
   // it is possible that none of the top-ranked k-mer has z>=cutoff
   // bug fix
   if (sum>0.001) {
      for (i=0; i<*numSpecifiedTopKmer; i++) {
         if (kp[i].z>=Z_CUTOFF) {
            kp[i].p=kp[i].z/sum; 
         }
         else break;
      }
      numAvailableTopKmer=i;
   }
   else {
      switch(kmerLen) {
         case 3: numAvailableTopKmer=10; break; 
         case 4: numAvailableTopKmer=20; break; 
         case 5: numAvailableTopKmer=40; break; 
         case 6: numAvailableTopKmer=60; break; 
         default: break; 
      }
   }

   word->count=numAvailableTopKmer;
   word->prob_sta[0]=0.0; 
   word->prob_end[0]=kp[0].p;
   strcpy(word->s1[0],kp[0].seq); 
   word->s1[0][kmerLen]='\0';

   for (i=1; i<numAvailableTopKmer; i++) {
      word->prob_sta[i]=word->prob_end[i-1]; 
      word->prob_end[i]=word->prob_sta[i]+kp[i].p;
      strcpy(word->s1[i],kp[i].seq); 
      word->s1[i][kmerLen]='\0';
   }
  
   // for (i=1; i<numAvailableTopKmer; i++) {
   //   printf("%s\t%5.3f\t%5.3f\n",word[0].s1[i],word[0].prob_sta[i],word[0].prob_end[i]); 
   //}
   return (numAvailableTopKmer); 
}

void enumerate_kmers(char **trimer,char **tetramer,char **pentamer) {

   register int i,j,k,l,m;
   int cn3,cn4,cn5;

   cn3=0; cn4=0; cn5=0;
   for (i=0; i<4; i++) {
      for (j=0; j<4; j++) {
         for (k=0; k<4; k++) {
            switch(i) {
               case 0:  trimer[cn3][0]='a'; break;
               case 1:  trimer[cn3][0]='c'; break;
               case 2:  trimer[cn3][0]='g'; break;
               case 3:  trimer[cn3][0]='t'; break;
               default: break;
            }
            switch(j) {
               case 0:  trimer[cn3][1]='a'; break;
               case 1:  trimer[cn3][1]='c'; break;
               case 2:  trimer[cn3][1]='g'; break;
               case 3:  trimer[cn3][1]='t'; break;
               default: break;
            }
            switch(k) {
               case 0:  trimer[cn3][2]='a'; break;
               case 1:  trimer[cn3][2]='c'; break;
               case 2:  trimer[cn3][2]='g'; break;
               case 3:  trimer[cn3][2]='t'; break;
               default: break;
            }
            trimer[cn3][3]='\0'; cn3++;

            for (l=0; l<4; l++) {
               switch(i) {
                  case 0:  tetramer[cn4][0]='a'; break;
                  case 1:  tetramer[cn4][0]='c'; break;
                  case 2:  tetramer[cn4][0]='g'; break;
                  case 3:  tetramer[cn4][0]='t'; break;
                  default: break;
               }
               switch(j) {
                  case 0:  tetramer[cn4][1]='a'; break;
                  case 1:  tetramer[cn4][1]='c'; break;
                  case 2:  tetramer[cn4][1]='g'; break;
                  case 3:  tetramer[cn4][1]='t'; break;
                  default: break;
               }
               switch(k) {
                  case 0:  tetramer[cn4][2]='a'; break;
                  case 1:  tetramer[cn4][2]='c'; break;
                  case 2:  tetramer[cn4][2]='g'; break;
                  case 3:  tetramer[cn4][2]='t'; break;
                  default: break;
               }
               switch(l) {
                  case 0:  tetramer[cn4][3]='a'; break;
                  case 1:  tetramer[cn4][3]='c'; break;
                  case 2:  tetramer[cn4][3]='g'; break;
                  case 3:  tetramer[cn4][3]='t'; break;
                  default: break;
               }
               tetramer[cn4][4]='\0'; cn4++;

               for (m=0; m<4; m++) {
                  switch(i) {
                     case 0:  pentamer[cn5][0]='a'; break;
                     case 1:  pentamer[cn5][0]='c'; break;
                     case 2:  pentamer[cn5][0]='g'; break;
                     case 3:  pentamer[cn5][0]='t'; break;
                     default: break;
                  }
                  switch(j) {
                     case 0:  pentamer[cn5][1]='a'; break;
                     case 1:  pentamer[cn5][1]='c'; break;
                     case 2:  pentamer[cn5][1]='g'; break;
                     case 3:  pentamer[cn5][1]='t'; break;
                     default: break;
                  }
                  switch(k) {
                     case 0:  pentamer[cn5][2]='a'; break;
                     case 1:  pentamer[cn5][2]='c'; break;
                     case 2:  pentamer[cn5][2]='g'; break;
                     case 3:  pentamer[cn5][2]='t'; break;
                     default: break;
                  }
                  switch(l) {
                     case 0:  pentamer[cn5][3]='a'; break;
                     case 1:  pentamer[cn5][3]='c'; break;
                     case 2:  pentamer[cn5][3]='g'; break;
                     case 3:  pentamer[cn5][3]='t'; break;
                     default: break;
                  }
                  switch(m) {
                     case 0:  pentamer[cn5][4]='a'; break;
                     case 1:  pentamer[cn5][4]='c'; break;
                     case 2:  pentamer[cn5][4]='g'; break;
                     case 3:  pentamer[cn5][4]='t'; break;
                     default: break;
                  }
                  pentamer[cn5][5]='\0'; cn5++;
               }
            }
         }
      }
   }
}

// find all k-tuples. overlapping ones within w discarded
void count_k_tuples(char **seq,char **rseq,int numSeq,int *seqLen,
   char **kmer,int numKer,int kmerLen,int *kmerCn) {

   register int i,j,k,m,l;
   int used,slide,numUniq;
   int *uniq,*id;
   char *s1,*s2;

   s1=alloc_char(kmerLen+1);
   s2=alloc_char(kmerLen+1);
   uniq=alloc_int(2*kmerLen);
   id=alloc_int(2*kmerLen);

   for (m=0; m<numKer; m++) kmerCn[m]=0;
   for (i=0; i<numSeq; i++) {

      for (m=0; m<2*kmerLen; m++) id[m]=-1;
      slide=0;
      for (j=0; j<seqLen[i]-kmerLen+1; j++) {
         for (k=0; k<kmerLen; k++)  s1[k]=seq[i][j+k];  s1[k]='\0';
         for (k=0; k<kmerLen; k++)  s2[k]=rseq[i][seqLen[i]-kmerLen-j+k];  s2[k]='\0';

         for (m=0; m<numKer; m++) {
            if (strncmp(s1,kmer[m],kmerLen)==0) { id[slide]=m; break; } 
         }
         slide++;

         for (m=0; m<numKer; m++) {
            if (strncmp(s2,kmer[m],kmerLen)==0) { id[slide]=m; break; } 
         }
         slide++;

         if ((slide==2*kmerLen) || (j==seqLen[i]-kmerLen)) {
            for (l=0; l<2*kmerLen; l++) uniq[l]=-2;

            numUniq=0;
            for (k=0; k<2*kmerLen; k++) {

               used=0;
               for (l=0; l<numUniq; l++) {
                  if (id[k]==uniq[l]) { used=1; break; } 
               }

               if (!used && id[k]!=-1) {
                  uniq[numUniq]=id[k]; numUniq++; 
               }
            }
            for (l=0; l<numUniq; l++)  (kmerCn[uniq[l]])++;
            slide=0;
            for (m=0; m<2*kmerLen; m++) id[m]=-1;
         }
      }
   }
   if (s1)   { free(s1);   s1=NULL;   }
   if (s2)   { free(s2);   s2=NULL;   }
   if (id)   { free(id);   id=NULL;   }
   if (uniq) { free(uniq); uniq=NULL; }
}

double *base_frequency(int numSeq,char **seq,int *seqLen) {

   register int i,j;
   int bcount[4];
   int sum;
   double *freq;

   freq=alloc_double(4);

   for (j=0; j<4; j++) bcount[j]=0;

   for (i=0; i<numSeq; i++) {
      for (j=0; j<seqLen[i]; j++) {
         switch (seq[i][j]) {
            case 'a': (bcount[0])++; break;
            case 'c': (bcount[1])++; break;
            case 'g': (bcount[2])++; break;
            case 't': (bcount[3])++; break;
            default: break;
         }
      }
   }

   sum=0; for (j=0; j<4; j++) sum +=bcount[j];

   if (bcount[0]<PSEUDO_COUNT || bcount[1]<PSEUDO_COUNT || bcount[2]<PSEUDO_COUNT || bcount[3]<PSEUDO_COUNT) {
      for (j=0; j<4; j++) {
         freq[j]=((double)bcount[j]+PSEUDO_COUNT)/((double)sum+PSEUDO_COUNT*4.0);   
      }   
   }
   else {
      for (j=0; j<4; j++) {
         freq[j]=(double)bcount[j]/(double)sum;   
      }   
   }

   freq[0]=(freq[0]+freq[3])/2.0; freq[3]=freq[0];
   freq[1]=(freq[1]+freq[2])/2.0; freq[2]=freq[1];
  
   return (freq);
}

