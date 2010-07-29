/*---------------------------------------------------------------*
 * PSSM Searcher                                                 *
 *                                                               *
 * By: Jes Frellsen, Ida Moltke and Martin Thiim                 *
 * University of Copenhagen, spring 2006                         *
 *---------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "timer.h"



void translate(unsigned char *pDest, char *pSrc, int strSize, char *pAlphabet) {
  unsigned char *dest;
  unsigned char *src;	
  unsigned char *alphabet;	
  int alpha[256];
  int i, asize;

  src = (unsigned char*)pSrc;
  dest = pDest;
  alphabet = (unsigned char*)pAlphabet;

  for(i = 0; i < 256; i++)
    alpha[i] = -1;
  
  asize = strlen(pAlphabet);		
  for(i = 0; i < asize; i++) {
    alpha[ alphabet[i] ] = i;
  }

  //alpha[ (unsigned char) STRINGSEP ] = STRINGSEPTRANS;

  for(i = 0; i < strSize; i++) {
    if( alpha[src[i]] < 0) {
      printf("Error: Letter '%c' found in text but not in alphabet\n", src[i]);
      exit(0);
    }
    else
      dest[i] = alpha[ (unsigned char) src[i]];
  }
}




void count_radix(unsigned char *pStr, int *resI, int n, int nAlph){
  const int max_depth = 10;

  int i, k;
  int *counts = malloc(sizeof(int)*nAlph+1);
  int *base = malloc(sizeof(int)*nAlph+1);
  
  unsigned char **data   = malloc(sizeof(char*)*n);
  unsigned char **result = malloc(sizeof(char*)*n);
  unsigned char **tmp;

  //  Reset counts
  for(i=0; i<nAlph; i++) {
    counts[i] = 0;
  }
    
  //Count occurrences of each letter
  for(i=max_depth+1; i<n; i++) {
    counts[pStr[i]]++;
  }

  counts[nAlph] = max_depth+1;

  // Initialize the data-set
  for(i=0; i<n; i++) {
    data[i] = pStr + i;
  }

  // Sort for each position from the back
  for(k=max_depth; k>=0; k--) {
    //Update the counts
    counts[pStr[k]]++;
    counts[nAlph]--;

    // Calculate the letter borders
    base[0] = 0;
    for(i=1; i<(nAlph+1); i++) {
      printf("Count[%i] %i\n", i-1, counts[i-1]);
      base[i] = base[i-1] + counts[i-1];
    }
      printf("Count[%i] %i\n", 5, counts[5]);
    
    // Insert data
    for(i=0; i<n; i++) {
      result[ base[*(data[i]+k)]++ ] = data[i];
    }

    // Swap
    tmp = data;
    data = result;
    result = tmp;
  }

  printf("Copying results\n");
  for(i=0; i<n; i++) 
    resI[i] = data[i] - pStr; 

  free(data);
  free(result);
}



// MULTIKEY QUICKSORT

#ifndef min
#define min(a, b) ((a)<=(b) ? (a) : (b))
#endif

#define swap(a, b) { char *t=x[a]; \
                     x[a]=x[b]; x[b]=t; }
#define i2c(i) x[i][depth]

inline void vecswap(int i, int j, int n, unsigned char **x)
{   while (n-- > 0) {
        swap(i, j);
        i++;
        j++;
    }
}

void ssort1(unsigned char **x, int n, int depth)
{   int    a, b, c, d, r, v;
    if (n <= 1)
        return;
    a = rand() % n;
    swap(0, a);
    v = i2c(0);
    a = b = 1;
    c = d = n-1;
    for (;;) {
        while (b <= c && (r = i2c(b)-v) <= 0) {
            if (r == 0) { swap(a, b); a++; }
            b++;
        }
        while (b <= c && (r = i2c(c)-v) >= 0) {
            if (r == 0) { swap(c, d); d--; }
            c--;
        }
        if (b > c) break;
        swap(b, c);
        b++;
        c--;
    }
    r = min(a, b-a);     vecswap(0, b-r, r, x);
    r = min(d-c, n-d-1); vecswap(b, n-r, r, x);
    r = b-a; ssort1(x, r, depth);
    if (i2c(r) != 6 && depth < 5)
        ssort1(x + r, a + n-d-1, depth+1);
    r = d-c; ssort1(x + n-r, r, depth);
}

void mkqs(unsigned char *pStr, int *resI, int n, int nAlph){
  int i;
  unsigned char **data = malloc(sizeof(char*)*n);

  // Initialize the data-set
  for(i=0; i<n; i++) {
    data[i] = pStr + i;
  }

  ssort1(data, n, 0);

  for(i=0; i<n; i++) 
    resI[i] = (data[i] - pStr); 

  free(data);
  
}

int main() {
  int i,k,n;
  char *text;
  unsigned char *trans;
  FILE *f;
  
  initTimer();
  
  // Open file
  f = fopen("/tmp/dbtss", "r");
  if(f == 0) {
    printf("File not found.\n");
    return 0;
  }
  // Obtain file size.
  fseek (f, 0 , SEEK_END);
  n = ftell (f);    
  rewind (f);
  printf("File size: %d\n", n);

  text = malloc(sizeof(char)*(n));
  k = fread(text, 1, n, f);
  printf("Read size: %d\n", k);
  fclose(f);
  
  addTimer("End reading");

  // Translate
  trans = malloc(sizeof(char)*(n+50));
  translate(trans, text, n, "ACGT|");
  addTimer("End translating");

  free(text);

  // Initialize tail
  for(i=n; i<n+50; i++)
    trans[i]=5;
  addTimer("End initialize tail");

  // Sort
  int *result = malloc(sizeof(int)*n);
  count_radix(trans, result, n, 5);

  addTimer("End sorting");

  printTimer();

  freeTimer();

  //  for(i=0; i<n; i++){
  //    for(k=0; k<3; k++)
      //      printf("%i", trans[result[i]+k]);
      //    printf("\n");
  //  }

  return 0;

}

