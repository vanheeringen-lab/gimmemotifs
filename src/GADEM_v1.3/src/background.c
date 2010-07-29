#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include <math.h>
#include <ctype.h>
#include <sys/types.h>
#include <unistd.h>
#include "gadem.h"
#include "random.h"
#include "defines.h"

void nonACGT(BACKGROUND_Model *back);

int read_userBackgModel(char *fileName,BACKGROUND_Model *back) {

   FILE *fp;
   char *buffer,*tok,*s1;
   int *cn,maxOligomer,len,len2,tabFound,maxAllowedOligomer,order;
   register int ii,i;
   double freq,*sum;

   fp=fopen(fileName,"r");
   if (!fp) { perror(fileName); exit(0); }

   maxAllowedOligomer=MAX_ORDER+1;
   buffer=alloc_char(256);
   s1=alloc_char(maxAllowedOligomer+1);
   sum=alloc_double(maxAllowedOligomer);
   cn=alloc_int(maxAllowedOligomer);
   for (i=0; i<maxAllowedOligomer; i++) { cn[i]=0; sum[i]=0; }

   maxOligomer=0;
   while (!feof(fp)) {
      if ((fgets(buffer,256,fp))>0)  {
         if (buffer[0]!='#') {
            len=strlen(buffer);
            buffer[len-1]='\0';
            tabFound=0;
            for (i=0; i<len; i++) { 
               if (buffer[i]=='\t') { tabFound=1; break; } 
            }
            if (tabFound) {
               tok=strtok(buffer,"\t");
               len2=strlen(tok);
               tok[len2]='\0';
               strcpy(s1,tok);

               if (len2>maxAllowedOligomer) continue;

               if (len2>maxOligomer) maxOligomer=len2;
               tok=strtok(0,"\t");
               freq=atof(tok);
               switch (len2) {
                  case 1:  { back->monomerFreq[cn[0]] =freq; strcpy(back->monomer[cn[0]],s1);  sum[0]+=freq; cn[0]++; break; }
                  case 2:  { back->dimerFreq[cn[1]]   =freq; strcpy(back->dimer[cn[1]],s1);    sum[1]+=freq; cn[1]++; break; }
                  case 3:  { back->trimerFreq[cn[2]]  =freq; strcpy(back->trimer[cn[2]],s1);   sum[2]+=freq; cn[2]++; break; }
                  case 4:  { back->tetramerFreq[cn[3]]=freq; strcpy(back->tetramer[cn[3]],s1); sum[3]+=freq; cn[3]++; break; }
                  case 5:  { back->pentamerFreq[cn[4]]=freq; strcpy(back->pentamer[cn[4]],s1); sum[4]+=freq; cn[4]++; break; }
                  case 6:  { back->hexamerFreq[cn[5]] =freq; strcpy(back->hexamer[cn[5]],s1);  sum[5]+=freq; cn[5]++; break; }
                  case 7:  { back->heptamerFreq[cn[6]]=freq; strcpy(back->heptamer[cn[6]],s1); sum[6]+=freq; cn[6]++; break; }
                  case 8:  { back->octamerFreq[cn[7]] =freq; strcpy(back->octamer[cn[7]],s1);  sum[7]+=freq; cn[7]++; break; }
                  case 9:  { back->nonamerFreq[cn[8]] =freq; strcpy(back->nonamer[cn[8]],s1);  sum[8]+=freq; cn[8]++; break; }
                  default: break; 
               } 
            }
            else {
               tok=strtok(buffer," ");
               len2=strlen(tok);
               if (len2>10) { printf("Error: up to 9th order is allowed!\n"); exit(0); }
               tok=strtok(0," ");
               freq=atof(tok);
               switch (len2) {
                  case 1:  { back->monomerFreq[cn[0]] =freq; sum[0]+=freq; cn[0]++; break; }
                  case 2:  { back->dimerFreq[cn[1]]   =freq; sum[1]+=freq; cn[1]++; break; }
                  case 3:  { back->trimerFreq[cn[2]]  =freq; sum[2]+=freq; cn[2]++; break; }
                  case 4:  { back->tetramerFreq[cn[3]]=freq; sum[3]+=freq; cn[3]++; break; }
                  case 5:  { back->pentamerFreq[cn[4]]=freq; sum[4]+=freq; cn[4]++; break; }
                  case 6:  { back->hexamerFreq[cn[5]] =freq; sum[5]+=freq; cn[5]++; break; }
                  case 7:  { back->heptamerFreq[cn[6]]=freq; sum[6]+=freq; cn[6]++; break; }
                  case 8:  { back->octamerFreq[cn[7]] =freq; sum[7]+=freq; cn[7]++; break; }
                  case 9:  { back->nonamerFreq[cn[8]] =freq; sum[8]+=freq; cn[8]++; break; }
                  default: break; 
               } 
            } 
         } 
      }
   }
   fclose(fp);

   /*---------------------------------------------------------------------------------
   for (i=0; i<4; i++)      printf("%15.14f\n",back->monomerFreq[i]);  printf("\n");
   for (i=0; i<16; i++)     printf("%15.14f\n",back->dimerFreq[i]);    printf("\n");
   for (i=0; i<64; i++)     printf("%15.14f\n",back->trimerFreq[i]);   printf("\n");
   for (i=0; i<256; i++)    printf("%15.14f\n",back->tetramerFreq[i]); printf("\n");
   for (i=0; i<1024; i++)   printf("%15.14f\n",back->pentamerFreq[i]); printf("\n");
   for (i=0; i<4096; i++)   printf("%15.14f\n",back->hexamerFreq[i]);  printf("\n\n");
   for (i=0; i<16384; i++)  printf("%15.14f\n",back->heptamerFreq[i]); printf("\n\n");
   for (i=0; i<65536; i++)  printf("%15.14f\n",back->octamerFreq[i]);  printf("\n\n");
   for (i=0; i<262144; i++) printf("%15.14f\n",back->nonamerFreq[i]);  printf("\n\n");
   ---------------------------------------------------------------------------------*/

   // check probability sum
   for (i=0; i<maxOligomer; i++) {
      if (fabs(sum[i]-1.0)>0.001) printf("sum of marginal %d: %8.6f\n",i+1,sum[i]); 
   }
   if (maxOligomer==0) { printf("Error: no frequencies in %s\n",fileName); exit(0); }

   order=0;
   for (ii=0; ii<maxOligomer; ii++) {
      transition_1st(back->dimerFreq,back->transition1);                   // compute transition from marginal
      order++;
      if (order==maxOligomer-1) break;
      transition_2nd(back->trimerFreq,back->transition2);
      order++;
      if (order==maxOligomer-1) break;
      transition_3rd(back->tetramerFreq,back->transition3);
      order++;
      if (order==maxOligomer-1) break;
      transition_4th(back->pentamerFreq,back->transition4);
      order++;
      if (order==maxOligomer-1) break;
      transition_5th(back->hexamerFreq,back->transition5);
      order++;
      if (order==maxOligomer-1) break;
      transition_6th(back->heptamerFreq,back->transition6);
      order++;
      if (order==maxOligomer-1) break;
      transition_7th(back->octamerFreq,back->transition7);
      order++;
      if (order==maxOligomer-1) break;
      transition_8th(back->nonamerFreq,back->transition8);
      order++;
      if (order==maxOligomer-1) break;
   }

   order=0;
   for (ii=0; ii<maxOligomer; ii++) {
      for (i=0; i<4; i++)  back->monomerFreq[i]=log(back->monomerFreq[i]); // log marginal
      for (i=0; i<16; i++) back->transition1[i]=log(back->transition1[i]); // log transition
      order++;
      if (order==maxOligomer-1) break;

      for (i=0; i<16; i++) back->dimerFreq[i]  =log(back->dimerFreq[i]);
      for (i=0; i<64; i++) back->transition2[i]=log(back->transition2[i]);
      order++;
      if (order==maxOligomer-1) break;

      for (i=0; i<64; i++)  back->trimerFreq[i] =log(back->trimerFreq[i]);
      for (i=0; i<256; i++) back->transition3[i]=log(back->transition3[i]);
      order++;
      if (order==maxOligomer-1) break;

      for (i=0; i<256; i++)  back->tetramerFreq[i]=log(back->tetramerFreq[i]);
      for (i=0; i<1024; i++) back->transition4[i] =log(back->transition4[i]);
      order++;
      if (order==maxOligomer-1) break;

      for (i=0; i<1024; i++) back->pentamerFreq[i]=log(back->pentamerFreq[i]);
      for (i=0; i<4096; i++) back->transition5[i] =log(back->transition5[i]);
      order++;
      if (order==maxOligomer-1) break;

      for (i=0; i<4096; i++)  back->hexamerFreq[i]=log(back->hexamerFreq[i]);
      for (i=0; i<16384; i++) back->transition6[i]=log(back->transition6[i]);
      order++;
      if (order==maxOligomer-1) break;

      for (i=0; i<16384; i++) back->heptamerFreq[i]=log(back->heptamerFreq[i]);
      for (i=0; i<65536; i++) back->transition7[i] =log(back->transition7[i]);
      order++;
      if (order==maxOligomer-1) break;

      for (i=0; i<65536; i++)  back->octamerFreq[i]=log(back->octamerFreq[i]);
      for (i=0; i<262144; i++) back->transition8[i]=log(back->transition8[i]);
      order++;
      if (order==maxOligomer-1) break;
   }
   nonACGT(back);

   /*  -------------------------debugging...          ------------------------- 
   for (i=0; i<4; i++)      printf("%8.5f\n",back->monomerFreq[i]);  printf("\n");
   for (i=0; i<16; i++)     printf("%8.5f\n",back->dimerFreq[i]);    printf("\n");
   for (i=0; i<64; i++)     printf("%8.5f\n",back->trimerFreq[i]);   printf("\n");
   for (i=0; i<256; i++)    printf("%8.5f\n",back->tetramerFreq[i]); printf("\n");
   for (i=0; i<1024; i++)   printf("%8.5f\n",back->pentamerFreq[i]); printf("\n");
   for (i=0; i<4096; i++)   printf("%8.5f\n",back->hexamerFreq[i]);  printf("\n\n");
   for (i=0; i<16384; i++)  printf("%8.5f\n",back->heptamerFreq[i]); printf("\n\n");
   for (i=0; i<65536; i++)  printf("%8.5f\n",back->octamerFreq[i]);  printf("\n\n");
   for (i=0; i<262144; i++) printf("%8.5f\n",back->nonamerFreq[i]);  printf("\n\n");
  
   for (i=0; i<16; i++)     printf("%8.5f\n",back->transition1[i]);  printf("\n");
   for (i=0; i<64; i++)     printf("%8.5f\n",back->transition2[i]);  printf("\n");
   for (i=0; i<256; i++)    printf("%8.5f\n",back->transition3[i]);  printf("\n");
   for (i=0; i<1024; i++)   printf("%8.5f\n",back->transition4[i]);  printf("\n");
   for (i=0; i<4096; i++)   printf("%8.5f\n",back->transition5[i]);  printf("\n\n");
   for (i=0; i<16384; i++)  printf("%8.5f\n",back->transition6[i]);  printf("\n\n");
   for (i=0; i<65536; i++)  printf("%8.5f\n",back->transition7[i]);  printf("\n\n");
   for (i=0; i<262144; i++) printf("%8.5f\n",back->transition8[i]);  printf("\n\n");
   exit(0); 
   ---------------------------debugging...          -----------------------------*/

   if (buffer) { free(buffer); buffer=NULL; }
   if (sum)    { free(sum);    sum=NULL;    }
   if (cn)     { free(cn);     cn=NULL;     }
   if (s1)     { free(s1);     s1=NULL;     }

   return (maxOligomer-1);
}

BACKGROUND_Model *alloc_background(void) {
 
   BACKGROUND_Model *back;
   
   back=NULL;

   back=(BACKGROUND_Model *)calloc(1,sizeof(BACKGROUND_Model));

   back->monomerFreq=alloc_double(4+1);       // marginal
   back->dimerFreq=alloc_double(16+1);        // marginal
   back->trimerFreq=alloc_double(64+1);       // marginal
   back->tetramerFreq=alloc_double(256+1);    // marginal
   back->pentamerFreq=alloc_double(1024+1);   // marginal
   back->hexamerFreq=alloc_double(4096+1);    // marginal
   back->heptamerFreq=alloc_double(16384+1);  // marginal
   back->octamerFreq=alloc_double(65536+1);   // marginal
   back->nonamerFreq=alloc_double(262144+1);  // marginal
   //back->decamerFreq=alloc_double(1048576); // marginal
   back->transition1=alloc_double(16+1);      // 1st order
   back->transition2=alloc_double(64+1);      // 2nd order
   back->transition3=alloc_double(256+1);     // 3rd order
   back->transition4=alloc_double(1024+1);    // 4th order
   back->transition5=alloc_double(4096+1);    // 5th order
   back->transition6=alloc_double(16384+1);   // 6th order
   back->transition7=alloc_double(65536+1);   // 7th order
   back->transition8=alloc_double(262144+1);  // 8th order
   //back->transition9=alloc_double(1048576); // 9th order
   back->monomer=alloc_char_char(4,2);       // k-mer
   back->dimer=alloc_char_char(16,3);        // k-mer
   back->trimer=alloc_char_char(64,4);       // k-mer
   back->tetramer=alloc_char_char(256,5);    // k-mer
   back->pentamer=alloc_char_char(1024,6);   // k-mer
   back->hexamer=alloc_char_char(4096,7);    // k-mer
   back->heptamer=alloc_char_char(16384,8);  // k-mer
   back->octamer=alloc_char_char(65536,9);   // k-mer
   back->nonamer=alloc_char_char(262144,10); // k-mer

   return (back);
}

void generate_background(int numSeq,char **seq,char **rseq,int *seqLen,BACKGROUND_Model *back,int MarkovOrder) {

   register int i;

   int numMonomer,numDimer,numTrimer,numTetramer,numPentamer,numHexamer,numHeptamer,numOctamer,numNonamer;
   int *monomerCn,*dimerCn,*trimerCn,*tetramerCn,*pentamerCn,*hexamerCn,*heptamerCn,*octamerCn,*nonamerCn;

   numMonomer=4; numDimer=16; numTrimer=64; numTetramer=256; numPentamer=1024; numHexamer=4096;
   numHeptamer=16384; numOctamer=65536; numNonamer=262144;

   monomerCn =alloc_int(4);     dimerCn   =alloc_int(16);      
   trimerCn  =alloc_int(64);    tetramerCn=alloc_int(256);  
   pentamerCn=alloc_int(1024);   hexamerCn=alloc_int(4096);  
   heptamerCn=alloc_int(16384);  octamerCn=alloc_int(65536);  
   nonamerCn =alloc_int(262144);  

   // enumerate all possible k-mers k=1,2,3,4,5,6,7,8,9
   numerate_monomer_to_pentamer(back);

   if (MarkovOrder==0) { 
      // count 1-mer
      monomerCn=count_nucleotides(seq,rseq,numSeq,seqLen,back->monomer,  numMonomer, 1); 
      marginal_prob(monomerCn,numMonomer,back->monomerFreq);

      for (i=0; i<4; i++)  back->monomerFreq[i]=log(back->monomerFreq[i]);
   }
   else if (MarkovOrder==1) { 
      // count 1-mer
      monomerCn=count_nucleotides(seq,rseq,numSeq,seqLen,back->monomer,  numMonomer, 1); 
      marginal_prob(monomerCn,numMonomer,back->monomerFreq);

      // count 2-mer
      dimerCn=count_nucleotides(seq,rseq,numSeq,seqLen,back->dimer,      numDimer,   2); 
      marginal_prob(dimerCn,numDimer,back->dimerFreq);
      transition_1st(back->dimerFreq,back->transition1);

      for (i=0; i<4; i++)  back->monomerFreq[i]=log(back->monomerFreq[i]);
      for (i=0; i<16; i++) back->transition1[i]=log(back->transition1[i]);
   }
   else if (MarkovOrder==2) { 
      // count 3-mer
      dimerCn=count_nucleotides(seq,rseq,numSeq,seqLen,back->dimer,      numDimer,   2); 
      marginal_prob(dimerCn,numDimer,back->dimerFreq);

      trimerCn  =count_nucleotides(seq,rseq,numSeq,seqLen,back->trimer,  numTrimer,  3); 
      marginal_prob(trimerCn,numTrimer,back->trimerFreq);
      transition_2nd(back->trimerFreq,back->transition2);

      for (i=0; i<16; i++) back->dimerFreq[i]  =log(back->dimerFreq[i]);
      for (i=0; i<64; i++) back->transition2[i]=log(back->transition2[i]);
   }
   else if (MarkovOrder==3) { 
      // count 4-mer
      trimerCn  =count_nucleotides(seq,rseq,numSeq,seqLen,back->trimer,  numTrimer,  3); 
      marginal_prob(trimerCn,numTrimer,back->trimerFreq);

      tetramerCn=count_nucleotides(seq,rseq,numSeq,seqLen,back->tetramer,numTetramer,4); 
      marginal_prob(tetramerCn,numTetramer,back->tetramerFreq);
      transition_3rd(back->tetramerFreq,back->transition3);

      for (i=0; i<64; i++)  back->trimerFreq[i] =log(back->trimerFreq[i]);
      for (i=0; i<256; i++) back->transition3[i]=log(back->transition3[i]);
   }
   else if (MarkovOrder==4) {
      // count 5-mer
      tetramerCn=count_nucleotides(seq,rseq,numSeq,seqLen,back->tetramer,numTetramer,4); 
      marginal_prob(tetramerCn,numTetramer,back->tetramerFreq);

      pentamerCn=count_nucleotides(seq,rseq,numSeq,seqLen,back->pentamer,numPentamer,5); 
      marginal_prob(pentamerCn,numPentamer,back->pentamerFreq);
      transition_4th(back->pentamerFreq,back->transition4);

      for (i=0; i<256; i++)  back->tetramerFreq[i]=log(back->tetramerFreq[i]);
      for (i=0; i<1024; i++) back->transition4[i] =log(back->transition4[i]);
   }
   else if (MarkovOrder==5) {
      // count 6-mer
      pentamerCn=count_nucleotides(seq,rseq,numSeq,seqLen,back->pentamer,numPentamer,5); 
      marginal_prob(pentamerCn,numPentamer,back->pentamerFreq);

      hexamerCn =count_nucleotides(seq,rseq,numSeq,seqLen,back->hexamer, numHexamer, 6);
      marginal_prob(hexamerCn,numHexamer,back->hexamerFreq);
      transition_5th(back->hexamerFreq,back->transition5);

      for (i=0; i<1024; i++) back->pentamerFreq[i]=log(back->pentamerFreq[i]);
      for (i=0; i<4096; i++) back->transition5[i] =log(back->transition5[i]);
   }
   else if (MarkovOrder==6) {
      // count 7-mer
      hexamerCn =count_nucleotides(seq,rseq,numSeq,seqLen,back->hexamer, numHexamer, 6);
      marginal_prob(hexamerCn,numHexamer,back->hexamerFreq);

      heptamerCn =count_nucleotides(seq,rseq,numSeq,seqLen,back->heptamer,numHeptamer, 7);
      marginal_prob(heptamerCn,numHeptamer,back->heptamerFreq);
      transition_6th(back->heptamerFreq,back->transition6);

      for (i=0; i<4096; i++)  back->hexamerFreq[i]=log(back->hexamerFreq[i]);
      for (i=0; i<16384; i++) back->transition6[i]=log(back->transition6[i]);
   }
   else if (MarkovOrder==7) {
      // count 8-mer
      heptamerCn =count_nucleotides(seq,rseq,numSeq,seqLen,back->heptamer,numHeptamer, 7);
      marginal_prob(heptamerCn,numHeptamer,back->heptamerFreq);

      octamerCn =count_nucleotides(seq,rseq,numSeq,seqLen,back->octamer, numOctamer, 8);
      marginal_prob(octamerCn,numOctamer,back->octamerFreq);
      transition_7th(back->octamerFreq,back->transition7);

      for (i=0; i<16384; i++) back->heptamerFreq[i]=log(back->heptamerFreq[i]);
      for (i=0; i<65536; i++) back->transition7[i] =log(back->transition7[i]);
   }
   else if (MarkovOrder==8) {
      // count 9-mer
      octamerCn =count_nucleotides(seq,rseq,numSeq,seqLen,back->octamer, numOctamer, 8);
      marginal_prob(octamerCn,numOctamer,back->octamerFreq);

      nonamerCn =count_nucleotides(seq,rseq,numSeq,seqLen,back->nonamer, numNonamer, 9);
      marginal_prob(nonamerCn,numNonamer,back->nonamerFreq);
      transition_8th(back->nonamerFreq,back->transition8);

      for (i=0; i<65536; i++)  back->octamerFreq[i]=log(back->octamerFreq[i]);
      for (i=0; i<262144; i++) back->transition8[i]=log(back->transition8[i]);
   }
   else if (MarkovOrder>8) { 
      printf("\nError: max Markov order: 8\n"); exit(0); 
   }
   else { }

   nonACGT(back);

   /*  ------------------------------debugging... -----------------------------------------------------
   if (MarkovOrder==1) {
      for (i=0; i<4; i++)      printf("%s\t%10.8f\n",back->monomer[i],back->monomerFreq[i]);  printf("\n");
      for (i=0; i<16; i++)     printf("%s\t%10.8f\n",back->dimer[i],back->transition1[i]);  printf("\n");
   }
   else if (MarkovOrder==2) {
      for (i=0; i<16; i++)     printf("%s\t%10.8f\n",back->dimer[i],back->dimerFreq[i]);    printf("\n");
      for (i=0; i<64; i++)     printf("%s\t%10.8f\n",back->trimer[i],back->transition2[i]);  printf("\n");
   }
   else if (MarkovOrder==3) {
      for (i=0; i<64; i++)     printf("%s\t%10.8f\n",back->trimer[i],back->trimerFreq[i]);   printf("\n");
      for (i=0; i<256; i++)    printf("%s\t%10.8f\n",back->tetramer[i],back->transition3[i]);  printf("\n");
   }
   else if (MarkovOrder==4) {
      for (i=0; i<256; i++)    printf("%s\t%10.8f\n",back->tetramer[i],back->tetramerFreq[i]); printf("\n");
      for (i=0; i<1024; i++)   printf("%s\t%10.8f\n",back->pentamer[i],back->transition4[i]);  printf("\n");
   }
   else if (MarkovOrder==5) {
      for (i=0; i<1024; i++)   printf("%s\t%10.8f\n",back->pentamer[i],back->pentamerFreq[i]); printf("\n");
      for (i=0; i<4096; i++)   printf("%s\t%10.8f\n",back->hexamer[i],back->transition5[i]);  printf("\n\n");
   }
   else if (MarkovOrder==6) {
      for (i=0; i<4096; i++)   printf("%s\t%10.8f\n",back->hexamer[i],back->hexamerFreq[i]);  printf("\n\n");
      for (i=0; i<16384; i++)  printf("%s\t%10.8f\n",back->heptamer[i],back->transition6[i]);  printf("\n\n");
   }
   else if (MarkovOrder==7) {
      for (i=0; i<16384; i++)  printf("%s\t%10.8f\n",back->heptamer[i],back->heptamerFreq[i]); printf("\n\n");
      for (i=0; i<65536; i++)  printf("%s\t%10.8f\n",back->octamer[i],back->transition7[i]);  printf("\n\n");
   }
   else if (MarkovOrder==8) {
      for (i=0; i<65536; i++)  printf("%s\t%10.8f\n",back->octamer[i],back->octamerFreq[i]);  
      for (i=0; i<262144; i++) printf("%s\t%10.8f\n",back->nonamer[i],back->transition8[i]);  printf("\n\n");
   } 
   exit(0); 
   ------------------------------debugging... -----------------------------------------------------*/

   if (monomerCn)    { free(monomerCn);     monomerCn=NULL;    }
   if (dimerCn)      { free(dimerCn);       dimerCn=NULL;      }
   if (trimerCn)     { free(trimerCn);      trimerCn=NULL;     }
   if (tetramerCn)   { free(tetramerCn);    tetramerCn=NULL;   }
   if (pentamerCn)   { free(pentamerCn);    pentamerCn=NULL;   }
   if (hexamerCn)    { free(hexamerCn);     hexamerCn=NULL;    }
   if (heptamerCn)   { free(heptamerCn);    heptamerCn=NULL;   }
   if (octamerCn)    { free(octamerCn);     octamerCn=NULL;    }
   if (nonamerCn)    { free(nonamerCn);     nonamerCn=NULL;    }
}

void numerate_monomer_to_pentamer(BACKGROUND_Model *back) {

   register int i,j,k,l,m,n,o,p,q;
   int cn2,cn3,cn4,cn5,cn6,cn7,cn8,cn9;

   cn2=0; cn3=0; cn4=0; cn5=0; cn6=0; cn7=0; cn8=0; cn9=0;

   for (i=0; i<4; i++) {
      switch(i) {
         case 0:  back->monomer[i][0]='a'; break; 
         case 1:  back->monomer[i][0]='c'; break; 
         case 2:  back->monomer[i][0]='g'; break; 
         case 3:  back->monomer[i][0]='t'; break; 
         default: break;
      } 
      back->monomer[i][1]='\0';

      for (j=0; j<4; j++) {
         switch(i) {
            case 0:  back->dimer[cn2][0]='a'; break; 
            case 1:  back->dimer[cn2][0]='c'; break; 
            case 2:  back->dimer[cn2][0]='g'; break; 
            case 3:  back->dimer[cn2][0]='t'; break; 
            default: break;
         } 
         switch(j) { 
            case 0:  back->dimer[cn2][1]='a'; break; 
            case 1:  back->dimer[cn2][1]='c'; break; 
            case 2:  back->dimer[cn2][1]='g'; break; 
            case 3:  back->dimer[cn2][1]='t'; break; 
            default: break;
         } 
         back->dimer[cn2][2]='\0'; cn2++;

         for (k=0; k<4; k++) {
            switch(i) {
               case 0:  back->trimer[cn3][0]='a'; break; 
               case 1:  back->trimer[cn3][0]='c'; break; 
               case 2:  back->trimer[cn3][0]='g'; break; 
               case 3:  back->trimer[cn3][0]='t'; break; 
               default: break;
            } 
            switch(j) { 
               case 0:  back->trimer[cn3][1]='a'; break; 
               case 1:  back->trimer[cn3][1]='c'; break; 
               case 2:  back->trimer[cn3][1]='g'; break; 
               case 3:  back->trimer[cn3][1]='t'; break; 
               default: break;
            } 
            switch(k) {
               case 0:  back->trimer[cn3][2]='a'; break; 
               case 1:  back->trimer[cn3][2]='c'; break; 
               case 2:  back->trimer[cn3][2]='g'; break; 
               case 3:  back->trimer[cn3][2]='t'; break; 
               default: break;
            }
            back->trimer[cn3][3]='\0';  cn3++;

            for (l=0; l<4; l++) {
               switch(i) {
                  case 0:  back->tetramer[cn4][0]='a'; break; 
                  case 1:  back->tetramer[cn4][0]='c'; break; 
                  case 2:  back->tetramer[cn4][0]='g'; break; 
                  case 3:  back->tetramer[cn4][0]='t'; break; 
                  default: break;
               } 
               switch(j) { 
                  case 0:  back->tetramer[cn4][1]='a'; break; 
                  case 1:  back->tetramer[cn4][1]='c'; break; 
                  case 2:  back->tetramer[cn4][1]='g'; break; 
                  case 3:  back->tetramer[cn4][1]='t'; break; 
                  default: break;
               } 
               switch(k) {
                  case 0:  back->tetramer[cn4][2]='a'; break; 
                  case 1:  back->tetramer[cn4][2]='c'; break; 
                  case 2:  back->tetramer[cn4][2]='g'; break; 
                  case 3:  back->tetramer[cn4][2]='t'; break; 
                  default: break;
               } 
               switch(l) {
                  case 0:  back->tetramer[cn4][3]='a'; break; 
                  case 1:  back->tetramer[cn4][3]='c'; break; 
                  case 2:  back->tetramer[cn4][3]='g'; break; 
                  case 3:  back->tetramer[cn4][3]='t'; break; 
                  default: break;
               } 
               back->tetramer[cn4][4]='\0'; cn4++; 
       
               for (m=0; m<4; m++) {
                  switch(i) {
                     case 0:  back->pentamer[cn5][0]='a'; break; 
                     case 1:  back->pentamer[cn5][0]='c'; break; 
                     case 2:  back->pentamer[cn5][0]='g'; break; 
                     case 3:  back->pentamer[cn5][0]='t'; break; 
                     default: break;
                  } 
                  switch(j) { 
                     case 0:  back->pentamer[cn5][1]='a'; break; 
                     case 1:  back->pentamer[cn5][1]='c'; break; 
                     case 2:  back->pentamer[cn5][1]='g'; break; 
                     case 3:  back->pentamer[cn5][1]='t'; break; 
                     default: break;
                  } 
                  switch(k) {
                     case 0:  back->pentamer[cn5][2]='a'; break; 
                     case 1:  back->pentamer[cn5][2]='c'; break; 
                     case 2:  back->pentamer[cn5][2]='g'; break; 
                     case 3:  back->pentamer[cn5][2]='t'; break; 
                     default: break;
                  } 
                  switch(l) {
                     case 0:  back->pentamer[cn5][3]='a'; break; 
                     case 1:  back->pentamer[cn5][3]='c'; break; 
                     case 2:  back->pentamer[cn5][3]='g'; break; 
                     case 3:  back->pentamer[cn5][3]='t'; break; 
                     default: break;
                  } 
                  switch(m) { 
                     case 0:  back->pentamer[cn5][4]='a'; break; 
                     case 1:  back->pentamer[cn5][4]='c'; break; 
                     case 2:  back->pentamer[cn5][4]='g'; break; 
                     case 3:  back->pentamer[cn5][4]='t'; break; 
                     default: break;
                  }
                  back->pentamer[cn5][5]='\0'; cn5++; 

                  for (n=0; n<4; n++) {
                     switch(i) {
                        case 0:  back->hexamer[cn6][0]='a'; break; 
                        case 1:  back->hexamer[cn6][0]='c'; break; 
                        case 2:  back->hexamer[cn6][0]='g'; break; 
                        case 3:  back->hexamer[cn6][0]='t'; break; 
                        default: break;
                     } 
                     switch(j) { 
                        case 0:  back->hexamer[cn6][1]='a'; break; 
                        case 1:  back->hexamer[cn6][1]='c'; break; 
                        case 2:  back->hexamer[cn6][1]='g'; break; 
                        case 3:  back->hexamer[cn6][1]='t'; break; 
                        default: break;
                     } 
                     switch(k) {
                        case 0:  back->hexamer[cn6][2]='a'; break; 
                        case 1:  back->hexamer[cn6][2]='c'; break; 
                        case 2:  back->hexamer[cn6][2]='g'; break; 
                        case 3:  back->hexamer[cn6][2]='t'; break; 
                        default: break;
                     } 
                     switch(l) {
                        case 0:  back->hexamer[cn6][3]='a'; break; 
                        case 1:  back->hexamer[cn6][3]='c'; break; 
                        case 2:  back->hexamer[cn6][3]='g'; break; 
                        case 3:  back->hexamer[cn6][3]='t'; break; 
                        default: break;
                     } 
                     switch(m) { 
                        case 0:  back->hexamer[cn6][4]='a'; break; 
                        case 1:  back->hexamer[cn6][4]='c'; break; 
                        case 2:  back->hexamer[cn6][4]='g'; break; 
                        case 3:  back->hexamer[cn6][4]='t'; break; 
                        default: break;
                     }
                     switch(n) { 
                        case 0:  back->hexamer[cn6][5]='a'; break; 
                        case 1:  back->hexamer[cn6][5]='c'; break; 
                        case 2:  back->hexamer[cn6][5]='g'; break; 
                        case 3:  back->hexamer[cn6][5]='t'; break; 
                        default: break;
                     }
                     back->hexamer[cn6][6]='\0'; cn6++; 

                     for (o=0; o<4; o++) {
                        switch(i) {
                           case 0:  back->heptamer[cn7][0]='a'; break; 
                           case 1:  back->heptamer[cn7][0]='c'; break; 
                           case 2:  back->heptamer[cn7][0]='g'; break; 
                           case 3:  back->heptamer[cn7][0]='t'; break; 
                           default: break;
                        } 
                        switch(j) { 
                           case 0:  back->heptamer[cn7][1]='a'; break; 
                           case 1:  back->heptamer[cn7][1]='c'; break; 
                           case 2:  back->heptamer[cn7][1]='g'; break; 
                           case 3:  back->heptamer[cn7][1]='t'; break; 
                           default: break;
                        } 
                        switch(k) {
                           case 0:  back->heptamer[cn7][2]='a'; break; 
                           case 1:  back->heptamer[cn7][2]='c'; break; 
                           case 2:  back->heptamer[cn7][2]='g'; break; 
                           case 3:  back->heptamer[cn7][2]='t'; break; 
                           default: break;
                        } 
                        switch(l) {
                           case 0:  back->heptamer[cn7][3]='a'; break; 
                           case 1:  back->heptamer[cn7][3]='c'; break; 
                           case 2:  back->heptamer[cn7][3]='g'; break; 
                           case 3:  back->heptamer[cn7][3]='t'; break; 
                           default: break;
                        } 
                        switch(m) { 
                           case 0:  back->heptamer[cn7][4]='a'; break; 
                           case 1:  back->heptamer[cn7][4]='c'; break; 
                           case 2:  back->heptamer[cn7][4]='g'; break; 
                           case 3:  back->heptamer[cn7][4]='t'; break; 
                           default: break;
                        }
                        switch(n) { 
                           case 0:  back->heptamer[cn7][5]='a'; break; 
                           case 1:  back->heptamer[cn7][5]='c'; break; 
                           case 2:  back->heptamer[cn7][5]='g'; break; 
                           case 3:  back->heptamer[cn7][5]='t'; break; 
                           default: break;
                        }
                        switch(o) { 
                           case 0:  back->heptamer[cn7][6]='a'; break; 
                           case 1:  back->heptamer[cn7][6]='c'; break; 
                           case 2:  back->heptamer[cn7][6]='g'; break; 
                           case 3:  back->heptamer[cn7][6]='t'; break; 
                           default: break;
                        }
                        back->heptamer[cn7][7]='\0'; cn7++; 

                        for (p=0; p<4; p++) {
                           switch(i) {
                              case 0:  back->octamer[cn8][0]='a'; break; 
                              case 1:  back->octamer[cn8][0]='c'; break; 
                              case 2:  back->octamer[cn8][0]='g'; break; 
                              case 3:  back->octamer[cn8][0]='t'; break; 
                              default: break;
                           } 
                           switch(j) { 
                              case 0:  back->octamer[cn8][1]='a'; break; 
                              case 1:  back->octamer[cn8][1]='c'; break; 
                              case 2:  back->octamer[cn8][1]='g'; break; 
                              case 3:  back->octamer[cn8][1]='t'; break; 
                              default: break;
                           } 
                           switch(k) {
                              case 0:  back->octamer[cn8][2]='a'; break; 
                              case 1:  back->octamer[cn8][2]='c'; break; 
                              case 2:  back->octamer[cn8][2]='g'; break; 
                              case 3:  back->octamer[cn8][2]='t'; break; 
                              default: break;
                           } 
                           switch(l) {
                              case 0:  back->octamer[cn8][3]='a'; break; 
                              case 1:  back->octamer[cn8][3]='c'; break; 
                              case 2:  back->octamer[cn8][3]='g'; break; 
                              case 3:  back->octamer[cn8][3]='t'; break; 
                              default: break;
                           } 
                           switch(m) { 
                              case 0:  back->octamer[cn8][4]='a'; break; 
                              case 1:  back->octamer[cn8][4]='c'; break; 
                              case 2:  back->octamer[cn8][4]='g'; break; 
                              case 3:  back->octamer[cn8][4]='t'; break; 
                              default: break;
                           }
                           switch(n) { 
                              case 0:  back->octamer[cn8][5]='a'; break; 
                              case 1:  back->octamer[cn8][5]='c'; break; 
                              case 2:  back->octamer[cn8][5]='g'; break; 
                              case 3:  back->octamer[cn8][5]='t'; break; 
                              default: break;
                           }
                           switch(o) { 
                              case 0:  back->octamer[cn8][6]='a'; break; 
                              case 1:  back->octamer[cn8][6]='c'; break; 
                              case 2:  back->octamer[cn8][6]='g'; break; 
                              case 3:  back->octamer[cn8][6]='t'; break; 
                              default: break;
                           }
                           switch(p) { 
                              case 0:  back->octamer[cn8][7]='a'; break; 
                              case 1:  back->octamer[cn8][7]='c'; break; 
                              case 2:  back->octamer[cn8][7]='g'; break; 
                              case 3:  back->octamer[cn8][7]='t'; break; 
                              default: break;
                           }
                           back->octamer[cn8][8]='\0'; cn8++; 

                           for (q=0; q<4; q++) {
                              switch(i) {
                                 case 0:  back->nonamer[cn9][0]='a'; break; 
                                 case 1:  back->nonamer[cn9][0]='c'; break; 
                                 case 2:  back->nonamer[cn9][0]='g'; break; 
                                 case 3:  back->nonamer[cn9][0]='t'; break; 
                                 default: break;
                              } 
                              switch(j) { 
                                 case 0:  back->nonamer[cn9][1]='a'; break; 
                                 case 1:  back->nonamer[cn9][1]='c'; break; 
                                 case 2:  back->nonamer[cn9][1]='g'; break; 
                                 case 3:  back->nonamer[cn9][1]='t'; break; 
                                 default: break;
                              } 
                              switch(k) {
                                 case 0:  back->nonamer[cn9][2]='a'; break; 
                                 case 1:  back->nonamer[cn9][2]='c'; break; 
                                 case 2:  back->nonamer[cn9][2]='g'; break; 
                                 case 3:  back->nonamer[cn9][2]='t'; break; 
                                 default: break;
                              } 
                              switch(l) {
                                 case 0:  back->nonamer[cn9][3]='a'; break; 
                                 case 1:  back->nonamer[cn9][3]='c'; break; 
                                 case 2:  back->nonamer[cn9][3]='g'; break; 
                                 case 3:  back->nonamer[cn9][3]='t'; break; 
                                 default: break;
                              } 
                              switch(m) { 
                                 case 0:  back->nonamer[cn9][4]='a'; break; 
                                 case 1:  back->nonamer[cn9][4]='c'; break; 
                                 case 2:  back->nonamer[cn9][4]='g'; break; 
                                 case 3:  back->nonamer[cn9][4]='t'; break; 
                                 default: break;
                              }
                              switch(n) { 
                                 case 0:  back->nonamer[cn9][5]='a'; break; 
                                 case 1:  back->nonamer[cn9][5]='c'; break; 
                                 case 2:  back->nonamer[cn9][5]='g'; break; 
                                 case 3:  back->nonamer[cn9][5]='t'; break; 
                                 default: break;
                              }
                              switch(o) { 
                                 case 0:  back->nonamer[cn9][6]='a'; break; 
                                 case 1:  back->nonamer[cn9][6]='c'; break; 
                                 case 2:  back->nonamer[cn9][6]='g'; break; 
                                 case 3:  back->nonamer[cn9][6]='t'; break; 
                                 default: break;
                              }
                              switch(p) { 
                                 case 0:  back->nonamer[cn9][7]='a'; break; 
                                 case 1:  back->nonamer[cn9][7]='c'; break; 
                                 case 2:  back->nonamer[cn9][7]='g'; break; 
                                 case 3:  back->nonamer[cn9][7]='t'; break; 
                                 default: break;
                              }
                              switch(q) { 
                                 case 0:  back->nonamer[cn9][8]='a'; break; 
                                 case 1:  back->nonamer[cn9][8]='c'; break; 
                                 case 2:  back->nonamer[cn9][8]='g'; break; 
                                 case 3:  back->nonamer[cn9][8]='t'; break; 
                                 default: break;
                              }
                              back->nonamer[cn9][9]='\0'; cn9++; 
                           }
                        }
                     }
                  } 
               } 
            } 
         } 
      } 
   }
}

int *count_nucleotides(char **seq,char **rseq,int numSeq,int *seqLen,char **word,int numWord,int wordLen) {

   register int i,j,k,m;
   int *wordCn;
   char *s1;

   s1=alloc_char(wordLen+1);
   wordCn=alloc_int(numWord);

   for (m=0; m<numWord; m++) wordCn[m]=0;

   for (i=0; i<numSeq; i++) {

      // plus strand
      for (j=0; j<seqLen[i]-wordLen+1; j++) {
         for (k=0; k<wordLen; k++) s1[k]=seq[i][j+k]; s1[k]='\0';
         for (m=0; m<numWord; m++) {
            if (strncmp(s1, word[m], (size_t) wordLen)==0) { (wordCn[m])++; break; } 
         }
      }

      // reverse complementary strand
      for (j=0; j<seqLen[i]-wordLen+1; j++) {
         for (k=0; k<wordLen; k++) s1[k]=rseq[i][j+k]; s1[k]='\0';
         for (m=0; m<numWord; m++) {
            if (strncmp(s1, word[m], (size_t) wordLen)==0) { (wordCn[m])++; break; } 
         }
      }
   }
   if (s1) { free(s1); s1=NULL; }

   return (wordCn);
}

void transition_1st(double *dimerFreq,double *transition1) {

   register int i,j;
   int cn1,cn2;
   double sum;

   /* dimer */
   cn1=0; cn2=0;
   for (i=0; i<4; i++) {
      sum=0; for (j=0; j<4; j++) { sum +=dimerFreq[cn1]; cn1++; }

      if (sum<=PSEUDO_FREQ) {
         for (j=0; j<4; j++) { transition1[cn2]=0.25; cn2++; }
      }
      else {
         for (j=0; j<4; j++)  { 
            transition1[cn2]=(PSEUDO_FREQ+dimerFreq[cn2])/(4*PSEUDO_FREQ+sum);
            cn2++; 
         }
      }
   }
}

void transition_2nd(double *trimerFreq,double *transition2) {

   register int i,k,j;
   int cn1,cn2;
   double sum;

   /* trimer */
   cn1=0; cn2=0;
   for (i=0; i<4; i++) {
      for (k=0; k<4; k++) {
         sum=0; for (j=0; j<4; j++) { sum +=trimerFreq[cn1]; cn1++; }

         if (sum<=PSEUDO_FREQ) {
            for (j=0; j<4; j++) { transition2[cn2]=0.25; cn2++; }
         }
         else {
            for (j=0; j<4; j++)  { 
               transition2[cn2]=(PSEUDO_FREQ+trimerFreq[cn2])/(4*PSEUDO_FREQ+sum);
               cn2++; 
            }
         }
      }
   }
}


void transition_3rd(double *tetramerFreq,double *transition3) {

   register int i,k,l,j;
   int cn1,cn2;
   double sum;

   /* tetramer */
   cn1=0; cn2=0;
   for (i=0; i<4; i++) {
      for (k=0; k<4; k++) {
         for (l=0; l<4; l++) {
            sum=0; for (j=0; j<4; j++) { sum +=tetramerFreq[cn1]; cn1++; }

            if (sum<=PSEUDO_FREQ) {
               for (j=0; j<4; j++) { transition3[cn2]=0.25; cn2++; }
            }
            else {
               for (j=0; j<4; j++)  { 
                  transition3[cn2]=(PSEUDO_FREQ+tetramerFreq[cn2])/(4*PSEUDO_FREQ+sum);
                  cn2++; 
               }
            }
         }
      }
   }
}

void transition_4th(double *pentamerFreq,double *transition4) {

   register int i,k,l,m,j;
   int cn1,cn2;
   double sum;

   /* pentamer */
   cn1=0; cn2=0;
   for (i=0; i<4; i++) {
      for (k=0; k<4; k++) {
         for (l=0; l<4; l++) {
            for (m=0; m<4; m++) {
               sum=0; for (j=0; j<4; j++) { sum +=pentamerFreq[cn1]; cn1++; }

               if (sum<=PSEUDO_FREQ) {
                  for (j=0; j<4; j++) { transition4[cn2]=0.25; cn2++; }
               }
               else {
                  for (j=0; j<4; j++)  { 
                     transition4[cn2]=(PSEUDO_FREQ+pentamerFreq[cn2])/(4*PSEUDO_FREQ+sum);
                     cn2++; 
                  }
               }
            }
         }
      }
   }
}

void transition_5th(double *hexamerFreq,double *transition5) {

   register int i,k,l,m,n,j;
   int cn1,cn2;
   double sum;

   /* hexamer */
   cn1=0; cn2=0;
   for (i=0; i<4; i++) {
      for (k=0; k<4; k++) {
         for (l=0; l<4; l++) {
            for (m=0; m<4; m++) {
               for (n=0; n<4; n++) {
                  sum=0; for (j=0; j<4; j++) { sum +=hexamerFreq[cn1]; cn1++; }
   
                  if (sum<=PSEUDO_FREQ) {
                     for (j=0; j<4; j++) { transition5[cn2]=0.25; cn2++; }
                  }
                  else {
                     for (j=0; j<4; j++)  { 
                        transition5[cn2]=(PSEUDO_FREQ+hexamerFreq[cn2])/(4*PSEUDO_FREQ+sum);
                        cn2++; 
                     }
                  }
               }
            }
         }
      }
   }
}

void transition_6th(double *heptamerFreq,double *transition6) {

   register int i,k,l,m,n,o,j;
   int cn1,cn2;
   double sum;

   /* heptamer */
   cn1=0; cn2=0;
   for (i=0; i<4; i++) {
      for (k=0; k<4; k++) {
         for (l=0; l<4; l++) {
            for (m=0; m<4; m++) {
               for (n=0; n<4; n++) {
                  for (o=0; o<4; o++) {
                     sum=0; for (j=0; j<4; j++) { sum +=heptamerFreq[cn1]; cn1++; }
   
                     if (sum<=PSEUDO_FREQ) {
                        for (j=0; j<4; j++) { transition6[cn2]=0.25; cn2++; }
                     }
                     else {
                        for (j=0; j<4; j++)  { 
                           transition6[cn2]=(PSEUDO_FREQ+heptamerFreq[cn2])/(4*PSEUDO_FREQ+sum);
                           cn2++; 
                        }
                     }
                  }
               }
            }
         }
      }
   }
}

void transition_7th(double *octamerFreq,double *transition7) {

   register int i,k,l,m,n,o,p,j;
   int cn1,cn2;
   double sum;

   /* octamerFreq */
   cn1=0; cn2=0;
   for (i=0; i<4; i++) {
      for (k=0; k<4; k++) {
         for (l=0; l<4; l++) {
            for (m=0; m<4; m++) {
               for (n=0; n<4; n++) {
                  for (o=0; o<4; o++) {
                     for (p=0; p<4; p++) {
                        sum=0; for (j=0; j<4; j++) { sum +=octamerFreq[cn1]; cn1++; }
   
                        if (sum<=PSEUDO_FREQ) {
                           for (j=0; j<4; j++) { transition7[cn2]=0.25; cn2++; }
                        }
                        else {
                           for (j=0; j<4; j++)  { 
                              transition7[cn2]=(PSEUDO_FREQ+octamerFreq[cn2])/(4*PSEUDO_FREQ+sum);
                              cn2++; 
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
}

void transition_8th(double *nonamerFreq,double *transition8) {

   register int i,k,l,m,n,o,p,q,j;
   int cn1,cn2;
   double sum;

   /* nonamerFreq */
   cn1=0; cn2=0;
   for (i=0; i<4; i++) {
      for (k=0; k<4; k++) {
         for (l=0; l<4; l++) {
            for (m=0; m<4; m++) {
               for (n=0; n<4; n++) {
                  for (o=0; o<4; o++) {
                     for (p=0; p<4; p++) {
                        for (q=0; q<4; q++) {
                           sum=0; for (j=0; j<4; j++) { sum +=nonamerFreq[cn1]; cn1++; }
   
                           if (sum<=PSEUDO_FREQ) {
                              for (j=0; j<4; j++) { transition8[cn2]=0.25; cn2++; }
                           }
                           else {
                              for (j=0; j<4; j++)  { 
                                 transition8[cn2]=(PSEUDO_FREQ+nonamerFreq[cn2])/(4*PSEUDO_FREQ+sum);
                                 cn2++; 
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
}

void marginal_prob(int *count,int numKmer,double *freq) {

   register int i;
   double sum;

   //printf("number of kmer=%d\n",numKmer);
   //for (i=0; i<numKmer; i++) printf("%d\n",count[i]);

   sum=0; for (i=0; i<numKmer; i++) sum +=(double)count[i];

   if (sum<=PSEUDO_FREQ) {
      printf("Error: data contains no [a,c,g,t].\n");  exit(0);
   }
   else {
      for (i=0; i<numKmer; i++) {
         freq[i]=(PSEUDO_FREQ+(double)count[i])/(sum+numKmer*PSEUDO_FREQ);
      }
   }
}

void ll_score_backg_model(int numSeq,char **seq,double **bscore,
   int *seqLen,int pwmLen,BACKGROUND_Model *back,int MarkovOrder) {

   register int i,j,k;
   int *base;
   char *s1; 

   s1=alloc_char(pwmLen+1);
   base=alloc_int(pwmLen);

   for (i=0; i<numSeq; i++) {
      for (j=0; j<seqLen[i]-pwmLen+1; j++) {

         for (k=0; k<pwmLen; k++) {
            switch(seq[i][j+k]) {
               case 'a': base[k]=0;  break; 
               case 'c': base[k]=1;  break; 
               case 'g': base[k]=2;  break; 
               case 't': base[k]=3;  break; 
               default:  base[k]=-1; break; 
            } 
         }

         if (MarkovOrder==1) {
            bscore[i][j]=0;
            if (base[0]==-1) bscore[i][j] +=back->monomerFreq[4];        // all marginal and transitional log-ed
            else             bscore[i][j] +=back->monomerFreq[base[0]];  // all marginal and transitional log-ed
      
            for (k=0; k<pwmLen-1; k++) {
               if (base[k]==-1 || base[k+1]==-1) bscore[i][j] +=back->transition1[16]; // transitional prob. all to 0.25
               else bscore[i][j] +=back->transition1[base[k]*4+base[k+1]];
            }
         }
         else if (MarkovOrder==2) {
            bscore[i][j]=0;
            if (base[0]==-1 || base[1]==-1) bscore[i][j] +=back->dimerFreq[16]; // marginal
            else bscore[i][j] +=back->dimerFreq[base[0]*4+base[1]];

            for (k=0; k<pwmLen-2; k++) {
               if (base[k]==-1||base[k+1]==-1||base[k+2]==-1) bscore[i][j] +=back->transition2[64]; //transitional
               else bscore[i][j] +=back->transition2[base[k]*16+base[k+1]*4+base[k+2]];
            }
         }
         else if (MarkovOrder==3) {
            bscore[i][j]=0;
            if (base[0]==-1||base[1]==-1||base[2]==-1) 
               bscore[i][j] +=back->trimerFreq[64]; //marginal
            else 
               bscore[i][j] +=back->trimerFreq[base[0]*16+base[1]*4+base[2]];

            for (k=0; k<pwmLen-3; k++) {
               if (base[k]==-1||base[k+1]==-1||base[k+2]==-1||base[k+3]==-1) 
                  bscore[i][j] +=back->transition3[256]; //transitional 
               else 
                  bscore[i][j] +=back->transition3[base[k]*64+base[k+1]*16+base[k+2]*4+base[k+3]];
            }
         } 
         else if (MarkovOrder==4) {
            bscore[i][j]=0;
            if (base[0]==-1||base[1]==-1||base[2]==-1||base[3]==-1) 
               bscore[i][j] +=back->tetramerFreq[256]; //marginal 
            else 
               bscore[i][j] +=back->tetramerFreq[base[0]*64+base[1]*16+base[2]*4+base[3]];

            for (k=0; k<pwmLen-4; k++) {
               if (base[k]==-1||base[k+1]==-1||base[k+2]==-1||base[k+3]==-1||base[k+4]==-1) 
                  bscore[i][j] +=back->transition4[1024]; // transitional
               else 
                  bscore[i][j] +=back->transition4[base[k]*256+base[k+1]*64+base[k+2]*16+base[k+3]*4+base[k+4]];
            }
         } 
         else if (MarkovOrder==5) {
            bscore[i][j]=0;
            if (base[0]==-1||base[1]==-1||base[2]==-1||base[3]==-1||base[4]==-1) 
               bscore[i][j] +=back->pentamerFreq[1024]; //marginal 
            else 
               bscore[i][j] +=back->pentamerFreq[base[0]*256+base[1]*64+base[2]*16+base[3]*4+base[4]];

            for (k=0; k<pwmLen-5; k++) {
               if (base[k]==-1||base[k+1]==-1||base[k+2]==-1||base[k+3]==-1||base[k+4]==-1||base[k+5]==-1) 
                  bscore[i][j] +=back->transition5[4096]; 
               else 
                  bscore[i][j] +=back->transition5[base[k]*1024+base[k+1]*256+base[k+2]*64+base[k+3]*16+base[k+4]*4+base[k+5]];
            }
         }
         else if (MarkovOrder==6) {
            bscore[i][j]=0;
            if (base[0]==-1||base[1]==-1||base[2]==-1||base[3]==-1||base[4]==-1||base[5]==-1) 
               bscore[i][j] +=back->hexamerFreq[4096]; 
            else 
               bscore[i][j] +=back->hexamerFreq[base[0]*1024+base[1]*256+base[2]*64+base[3]*16+base[4]*4+base[5]];

            for (k=0; k<pwmLen-6; k++) {
               if (base[k]==-1||base[k+1]==-1||base[k+2]==-1||base[k+3]==-1||base[k+4]==-1||base[k+5]==-1||base[k+6]==-1) 
                  bscore[i][j] +=back->transition6[16384]; 
               else 
                  bscore[i][j] +=back->transition6[base[k]*4096+base[k+1]*1024+base[k+2]*256+base[k+3]*64+base[k+4]*16+base[k+5]*4+base[k+6]];
            }
         }
         else if (MarkovOrder==7) {
            bscore[i][j]=0;
            if (base[0]==-1||base[1]==-1||base[2]==-1||base[3]==-1||base[4]==-1||base[5]==-1||base[6]==-1) 
               bscore[i][j] +=back->heptamerFreq[16384]; 
            else 
               bscore[i][j] +=back->heptamerFreq[base[0]*4096+base[1]*1024+base[2]*256+base[3]*64+base[4]*16+base[5]*4+base[6]];

            for (k=0; k<pwmLen-7; k++) {
               if (base[k]==-1||base[k+1]==-1||base[k+2]==-1||base[k+3]==-1||base[k+4]==-1||base[k+5]==-1||base[k+6]==-1||base[k+7]==-1) 
                  bscore[i][j] +=back->transition7[65536]; 
               else 
                  bscore[i][j] +=back->transition7[base[k]*16384+base[k+1]*4096+base[k+2]*1024+base[k+3]*256+base[k+4]*64+base[k+5]*16+base[k+6]*4+base[k+7]];
            }
         }
         else if (MarkovOrder==8) {
            bscore[i][j]=0;
            if (base[0]==-1||base[1]==-1||base[2]==-1||base[3]==-1||base[4]==-1||base[5]==-1||base[6]==-1||base[7]==-1) 
               bscore[i][j] +=back->octamerFreq[65536]; 
            else 
               bscore[i][j]+=back->octamerFreq[base[0]*16384+base[1]*4096+base[2]*1024+base[3]*256+base[4]*64+base[5]*16+base[6]*4+base[7]];

            for (k=0; k<pwmLen-8; k++) {
               if (base[k]==-1||base[k+1]==-1||base[k+2]==-1||base[k+3]==-1||base[k+4]==-1||base[k+5]==-1||base[k+6]==-1||base[k+7]==-1||base[k+8]==-1) 
                  bscore[i][j] +=back->transition8[262144]; 
               else 
                  bscore[i][j] +=back->transition8[base[k]*65536+base[k+1]*16384+base[k+2]*4096+base[k+3]*1024+base[k+4]*256+base[k+5]*64+base[k+6]*16+base[k+7]*4+base[k+8]];
            }
         }
         else {
            // default: 0th order
            bscore[i][j]=0;
            for (k=0; k<pwmLen; k++) {
               switch (seq[i][j+k]) {
                  case 'a': bscore[i][j] +=back->monomerFreq[0]; break; 
                  case 'c': bscore[i][j] +=back->monomerFreq[1]; break; 
                  case 'g': bscore[i][j] +=back->monomerFreq[2]; break; 
                  case 't': bscore[i][j] +=back->monomerFreq[3]; break; 
                  default:  bscore[i][j] +=back->monomerFreq[4]; break;
               } 
            } 
         } 
      } 
   }
   if (base) { free(base); base=NULL; }
}

void nonACGT(BACKGROUND_Model *back) {

   back->monomerFreq[4]     =-log(4.0);
   back->dimerFreq[16]      =-log(16.0);
   back->trimerFreq[64]     =-log(64.0);
   back->tetramerFreq[256]  =-log(256.0);
   back->pentamerFreq[1024] =-log(1024.0);
   back->hexamerFreq[4096]  =-log(4096.0);
   back->heptamerFreq[16384]=-log(16384.0);
   back->octamerFreq[65536] =-log(65536.0);
   back->transition1[16]    = log(0.25);
   back->transition2[64]    = log(0.25);
   back->transition3[256]   = log(0.25);
   back->transition4[1024]  = log(0.25);
   back->transition5[4096]  = log(0.25);
   back->transition6[16384] = log(0.25);
   back->transition7[65536] = log(0.25);
   back->transition8[262144]= log(0.25);
}

void simulate_background_seq(double *bfreq,int numSeq,int *seqLen,char **bseq) {

   register int i,j,k;
   double rand,sum;
   // FILE *fp;

   for (i=0; i<numSeq; i++) {
      for (j=0; j<seqLen[i]; j++) { 
         rand=genrand();
         sum=0;
         for (k=0; k<4; k++) {
            sum+=bfreq[k];
            if (rand>sum-bfreq[k] && rand<=sum) { 
               switch (k) {
                  case 0: bseq[i][j]='a'; break; 
                  case 1: bseq[i][j]='c'; break; 
                  case 2: bseq[i][j]='g'; break; 
                  case 3: bseq[i][j]='t'; break; 
                  default: break; 
               } 
            } 
         }
      } 
      bseq[i][j]='\0';
   }
   /*-------------------------------------------------
   fp=fopen("simulated_0th.seq","w");
   for (i=0; i<numSeq; i++) {
      fprintf(fp,">test\n");
      fprintf(fp,"%s\n",bseq[i]); 
   } fclose(fp);
   -------------------------------------------------*/
}

void simulate_backg_seq_bmodel(BACKGROUND_Model *back,int MarkovOrder,int numSeq,int *seqLen,char **bseq) {

   register int i,j,k;
   int id1,id2,id3,id4,id5,id6,id7,id8,ptr;
   double rand,sum,p;
   // FILE *fp;

   if (MarkovOrder==0) {
      for (i=0; i<numSeq; i++) {
         for (j=0; j<seqLen[i]; j++) { 
            rand=genrand();
            sum=0;
            for (k=0; k<4; k++) {
               p=exp(back->monomerFreq[k]);
               sum+=p;
               if (rand>(sum-p) && rand<=sum) { bseq[i][j]=back->monomer[k][0];  break; }
            }
         } 
         bseq[i][j]='\0';
      }
   }
   else if (MarkovOrder==1) {
      for (i=0; i<numSeq; i++) {
         rand=genrand();
         sum=0;
         for (k=0; k<4; k++) {
            p=exp(back->monomerFreq[k]);
            sum +=p;
            if (rand>(sum-p) && rand<=sum) { bseq[i][0]=back->monomer[k][0]; break; }
         }

         for (j=0; j<seqLen[i]-1; j++) {
            rand=genrand(); id1=0;
            switch (bseq[i][j]) {
               case 'a': id1=0; break; 
               case 'c': id1=1; break; 
               case 'g': id1=2; break; 
               case 't': id1=3; break; 
               default: break;
            }
            sum=0; ptr=id1*4;
            for (k=0; k<4; k++) {
               p=exp(back->transition1[ptr+k]);
               sum+=p;
               if (rand>(sum-p) && rand<=sum) {
                  bseq[i][j+1]=back->dimer[ptr+k][1]; break;
               }
            }
         } 
         bseq[i][seqLen[i]]='\0';
      }
   }
   else if (MarkovOrder==2) {
      for (i=0; i<numSeq; i++) {
         rand=genrand();
         sum=0;
         for (k=0; k<16; k++) {
            p=exp(back->dimerFreq[k]);
            sum+=p;
            if (rand>(sum-p) && rand<=sum) { 
               bseq[i][0]=back->dimer[k][0]; bseq[i][1]=back->dimer[k][1]; break;
            }
         }

         for (j=0; j<seqLen[i]-2; j++) {
            rand=genrand(); id1=-1; id2=-1;
            switch (bseq[i][j]) {
               case 'a': id1=0; break; 
               case 'c': id1=1; break; 
               case 'g': id1=2; break; 
               case 't': id1=3; break; 
               default: break;
            }
            switch (bseq[i][j+1]) {
               case 'a': id2=0; break; 
               case 'c': id2=1; break; 
               case 'g': id2=2; break; 
               case 't': id2=3; break; 
               default: break;
            }
            sum=0; ptr=id1*16+id2*4;
            for (k=0; k<4; k++) {
               p=exp(back->transition2[ptr+k]);
               sum+=p;
               if (rand>(sum-p) && rand<=sum) {
                  bseq[i][j+2]=back->trimer[ptr+k][2]; break;
               }
            }
         } 
         bseq[i][seqLen[i]]='\0';
      }
   }
   else if (MarkovOrder==3) {
      for (i=0; i<numSeq; i++) {
         rand=genrand();
         sum=0;
         for (k=0; k<64; k++) {
            p=exp(back->trimerFreq[k]);
            sum+=p;
            if (rand>(sum-p) && rand<=sum) { 
               bseq[i][0]=back->trimer[k][0]; bseq[i][1]=back->trimer[k][1]; bseq[i][2]=back->trimer[k][2]; break;
            }
         }
         for (j=0; j<seqLen[i]-3; j++) {
            rand=genrand(); id1=-1; id2=-1; id3=-1;
            switch (bseq[i][j]) {
               case 'a': id1=0; break; 
               case 'c': id1=1; break; 
               case 'g': id1=2; break; 
               case 't': id1=3; break; 
               default: break;
            }
            switch (bseq[i][j+1]) {
               case 'a': id2=0; break; 
               case 'c': id2=1; break; 
               case 'g': id2=2; break; 
               case 't': id2=3; break; 
               default: break;
            }
            switch (bseq[i][j+2]) {
               case 'a': id3=0; break; 
               case 'c': id3=1; break; 
               case 'g': id3=2; break; 
               case 't': id3=3; break; 
               default: break;
            }
            sum=0; ptr=id1*64+id2*16+id3*4;
            for (k=0; k<4; k++) {
               p=exp(back->transition3[ptr+k]);
               sum+=p;
               if (rand>(sum-p) && rand<=sum) {
                  bseq[i][j+3]=back->tetramer[ptr+k][3]; break;
               }
            }
         } bseq[i][seqLen[i]]='\0';
      }
   }
   else if (MarkovOrder==4) {
      for (i=0; i<numSeq; i++) {
         rand=genrand();
         sum=0;
         for (k=0; k<256; k++) {
            p=exp(back->tetramerFreq[k]);
            sum+=p;
            if (rand>(sum-p) && rand<=sum) { 
               bseq[i][0]=back->tetramer[k][0]; bseq[i][1]=back->tetramer[k][1]; 
               bseq[i][2]=back->tetramer[k][2]; bseq[i][3]=back->tetramer[k][3]; break;
            }
         }
         for (j=0; j<seqLen[i]-4; j++) {
            rand=genrand(); id1=-1; id2=-1; id3=-1; id4=-1;
            switch (bseq[i][j]) {
               case 'a': id1=0; break; 
               case 'c': id1=1; break; 
               case 'g': id1=2; break; 
               case 't': id1=3; break; 
               default: break;
            }
            switch (bseq[i][j+1]) {
               case 'a': id2=0; break; 
               case 'c': id2=1; break; 
               case 'g': id2=2; break; 
               case 't': id2=3; break; 
               default: break;
            }
            switch (bseq[i][j+2]) {
               case 'a': id3=0; break; 
               case 'c': id3=1; break; 
               case 'g': id3=2; break; 
               case 't': id3=3; break; 
               default: break;
            }
            switch (bseq[i][j+3]) {
               case 'a': id4=0; break; 
               case 'c': id4=1; break; 
               case 'g': id4=2; break; 
               case 't': id4=3; break; 
               default: break;
            }
            sum=0; ptr=id1*256+id2*64+id3*16+id4*4;
            for (k=0; k<4; k++) {
               p=exp(back->transition4[ptr+k]);
               sum+=p;
               if (rand>(sum-p) && rand<=sum) {
                  bseq[i][j+4]=back->pentamer[ptr+k][4]; break;
               }
            }
         } bseq[i][seqLen[i]]='\0';
      }
   }
   else if (MarkovOrder==5) {
      for (i=0; i<numSeq; i++) {
         rand=genrand();
         sum=0;
         for (k=0; k<1024; k++) {
            p=exp(back->pentamerFreq[k]);
            sum+=p;
            if (rand>(sum-p) && rand<=sum) { 
               bseq[i][0]=back->pentamer[k][0]; bseq[i][1]=back->pentamer[k][1]; bseq[i][2]=back->pentamer[k][2]; 
               bseq[i][3]=back->pentamer[k][3]; bseq[i][4]=back->pentamer[k][4]; break;
            }
         }
         for (j=0; j<seqLen[i]-5; j++) {
            rand=genrand(); id1=-1; id2=-1; id3=-1; id4=-1; id5=-1;
            switch (bseq[i][j]) {
               case 'a': id1=0; break; 
               case 'c': id1=1; break; 
               case 'g': id1=2; break; 
               case 't': id1=3; break; 
               default: break;
            }
            switch (bseq[i][j+1]) {
               case 'a': id2=0; break; 
               case 'c': id2=1; break; 
               case 'g': id2=2; break; 
               case 't': id2=3; break; 
               default: break;
            }
            switch (bseq[i][j+2]) {
               case 'a': id3=0; break; 
               case 'c': id3=1; break; 
               case 'g': id3=2; break; 
               case 't': id3=3; break; 
               default: break;
            }
            switch (bseq[i][j+3]) {
               case 'a': id4=0; break; 
               case 'c': id4=1; break; 
               case 'g': id4=2; break; 
               case 't': id4=3; break; 
               default: break;
            }
            switch (bseq[i][j+4]) {
               case 'a': id5=0; break; 
               case 'c': id5=1; break; 
               case 'g': id5=2; break; 
               case 't': id5=3; break; 
               default: break;
            }
            sum=0; ptr=id1*1024+id2*256+id3*64+id4*16+id5*4;
            for (k=0; k<4; k++) {
               p=exp(back->transition5[ptr+k]);
               sum+=p;
               if (rand>(sum-p) && rand<=sum) {
                  bseq[i][j+5]=back->hexamer[ptr+k][5]; break;
               }
            }
         } bseq[i][seqLen[i]]='\0';
      }
   }
   else if (MarkovOrder==6) {
      for (i=0; i<numSeq; i++) {
         rand=genrand();
         sum=0;
         for (k=0; k<4096; k++) {
            p=exp(back->hexamerFreq[k]);
            sum+=p;
            if (rand>(sum-p) && rand<=sum) { 
               bseq[i][0]=back->hexamer[k][0]; bseq[i][1]=back->hexamer[k][1]; bseq[i][2]=back->hexamer[k][2]; 
               bseq[i][3]=back->hexamer[k][3]; bseq[i][4]=back->hexamer[k][4]; bseq[i][5]=back->hexamer[k][5]; break;
            }
         }
         for (j=0; j<seqLen[i]-6; j++) {
            rand=genrand(); id1=-1; id2=-1; id3=-1; id4=-1; id5=-1; id6=-1;
            switch (bseq[i][j]) {
               case 'a': id1=0; break; 
               case 'c': id1=1; break; 
               case 'g': id1=2; break; 
               case 't': id1=3; break; 
               default: break;
            }
            switch (bseq[i][j+1]) {
               case 'a': id2=0; break; 
               case 'c': id2=1; break; 
               case 'g': id2=2; break; 
               case 't': id2=3; break; 
               default: break;
            }
            switch (bseq[i][j+2]) {
               case 'a': id3=0; break; 
               case 'c': id3=1; break; 
               case 'g': id3=2; break; 
               case 't': id3=3; break; 
               default: break;
            }
            switch (bseq[i][j+3]) {
               case 'a': id4=0; break; 
               case 'c': id4=1; break; 
               case 'g': id4=2; break; 
               case 't': id4=3; break; 
               default: break;
            }
            switch (bseq[i][j+4]) {
               case 'a': id5=0; break; 
               case 'c': id5=1; break; 
               case 'g': id5=2; break; 
               case 't': id5=3; break; 
               default: break;
            }
            switch (bseq[i][j+5]) {
               case 'a': id6=0; break; 
               case 'c': id6=1; break; 
               case 'g': id6=2; break; 
               case 't': id6=3; break; 
               default: break;
            }
            sum=0; ptr=id1*4096+id2*1024+id3*256+id4*64+id5*16+id6*4;
            for (k=0; k<4; k++) {
               p=exp(back->transition6[ptr+k]);
               sum+=p;
               if (rand>(sum-p) && rand<=sum) {
                  bseq[i][j+6]=back->heptamer[ptr+k][6]; break;
               }
            }
         } bseq[i][seqLen[i]]='\0';
      }
   }
   else if (MarkovOrder==7) {
      for (i=0; i<numSeq; i++) {
         rand=genrand();
         sum=0;
         for (k=0; k<16384; k++) {
            p=exp(back->heptamerFreq[k]);
            sum+=p;
            if (rand>(sum-p) && rand<=sum) { 
               bseq[i][0]=back->heptamer[k][0]; bseq[i][1]=back->heptamer[k][1]; bseq[i][2]=back->heptamer[k][2]; bseq[i][3]=back->heptamer[k][3]; 
               bseq[i][4]=back->heptamer[k][4]; bseq[i][5]=back->heptamer[k][5]; bseq[i][6]=back->heptamer[k][6]; break;
            }
         }
         for (j=0; j<seqLen[i]-7; j++) {
            rand=genrand(); id1=-1; id2=-1; id3=-1; id4=-1; id5=-1; id6=-1; id7=-1;
            switch (bseq[i][j]) {
               case 'a': id1=0; break; 
               case 'c': id1=1; break; 
               case 'g': id1=2; break; 
               case 't': id1=3; break; 
               default: break;
            }
            switch (bseq[i][j+1]) {
               case 'a': id2=0; break; 
               case 'c': id2=1; break; 
               case 'g': id2=2; break; 
               case 't': id2=3; break; 
               default: break;
            }
            switch (bseq[i][j+2]) {
               case 'a': id3=0; break; 
               case 'c': id3=1; break; 
               case 'g': id3=2; break; 
               case 't': id3=3; break; 
               default: break;
            }
            switch (bseq[i][j+3]) {
               case 'a': id4=0; break; 
               case 'c': id4=1; break; 
               case 'g': id4=2; break; 
               case 't': id4=3; break; 
               default: break;
            }
            switch (bseq[i][j+4]) {
               case 'a': id5=0; break; 
               case 'c': id5=1; break; 
               case 'g': id5=2; break; 
               case 't': id5=3; break; 
               default: break;
            }
            switch (bseq[i][j+5]) {
               case 'a': id6=0; break; 
               case 'c': id6=1; break; 
               case 'g': id6=2; break; 
               case 't': id6=3; break; 
               default: break;
            }
            switch (bseq[i][j+6]) {
               case 'a': id7=0; break; 
               case 'c': id7=1; break; 
               case 'g': id7=2; break; 
               case 't': id7=3; break; 
               default: break;
            }
            sum=0; ptr=id1*16384+id2*4096+id3*1024+id4*256+id5*64+id6*16+id7*4;
            for (k=0; k<4; k++) {
               p=exp(back->transition7[ptr+k]);
               sum+=p;
               if (rand>(sum-p) && rand<=sum) {
                  bseq[i][j+7]=back->octamer[ptr+k][7]; break;
               }
            }
         } bseq[i][seqLen[i]]='\0';
      }
   }
   else if (MarkovOrder==8) {
      for (i=0; i<numSeq; i++) {
         rand=genrand();
         sum=0;
         for (k=0; k<65536; k++) {
            p=exp(back->octamerFreq[k]);
            sum+=p;
            if (rand>(sum-p) && rand<=sum) { 
               bseq[i][0]=back->octamer[k][0]; bseq[i][1]=back->octamer[k][1]; bseq[i][2]=back->octamer[k][2]; bseq[i][3]=back->octamer[k][3]; 
               bseq[i][4]=back->octamer[k][4]; bseq[i][5]=back->octamer[k][5]; bseq[i][6]=back->octamer[k][6]; bseq[i][7]=back->octamer[k][7]; break;
            }
         }
         for (j=0; j<seqLen[i]-8; j++) {
            rand=genrand(); id1=-1; id2=-1; id3=-1; id4=-1; id5=-1; id6=-1; id7=-1; id8=-1;
            switch (bseq[i][j]) {
               case 'a': id1=0; break; 
               case 'c': id1=1; break; 
               case 'g': id1=2; break; 
               case 't': id1=3; break; 
               default: break;
            }
            switch (bseq[i][j+1]) {
               case 'a': id2=0; break; 
               case 'c': id2=1; break; 
               case 'g': id2=2; break; 
               case 't': id2=3; break; 
               default: break;
            }
            switch (bseq[i][j+2]) {
               case 'a': id3=0; break; 
               case 'c': id3=1; break; 
               case 'g': id3=2; break; 
               case 't': id3=3; break; 
               default: break;
            }
            switch (bseq[i][j+3]) {
               case 'a': id4=0; break; 
               case 'c': id4=1; break; 
               case 'g': id4=2; break; 
               case 't': id4=3; break; 
               default: break;
            }
            switch (bseq[i][j+4]) {
               case 'a': id5=0; break; 
               case 'c': id5=1; break; 
               case 'g': id5=2; break; 
               case 't': id5=3; break; 
               default: break;
            }
            switch (bseq[i][j+5]) {
               case 'a': id6=0; break; 
               case 'c': id6=1; break; 
               case 'g': id6=2; break; 
               case 't': id6=3; break; 
               default: break;
            }
            switch (bseq[i][j+6]) {
               case 'a': id7=0; break; 
               case 'c': id7=1; break; 
               case 'g': id7=2; break; 
               case 't': id7=3; break; 
               default: break;
            }
            switch (bseq[i][j+7]) {
               case 'a': id8=0; break; 
               case 'c': id8=1; break; 
               case 'g': id8=2; break; 
               case 't': id8=3; break; 
               default: break;
            }
            sum=0; ptr=id1*65536+id2*16384+id3*4096+id4*1024+id5*256+id6*64+id7*16+id8*4;
            for (k=0; k<4; k++) {
               p=exp(back->transition8[ptr+k]);
               sum+=p;
               if (rand>(sum-p) && rand<=sum) {
                  bseq[i][j+8]=back->nonamer[ptr+k][8]; break;
               }
            }
         } bseq[i][seqLen[i]]='\0';
      }
   }
   else { printf("Error: max order: 8\n"); exit(0); }

   /*-------------------------------------------------
   fp=fopen("simulated_0th.seq","w");
   for (i=0; i<numSeq; i++) {
      fprintf(fp,">test\n");
      fprintf(fp,"%s\n",bseq[i]); 
   } fclose(fp);
   -------------------------------------------------*/
}

