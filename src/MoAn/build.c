/*---------------------------------------------------------------*
 * PSSM Searcher                                                 *
 *                                                               *
 * By: Jes Frellsen, Ida Moltke and Martin Thiim                 *
 * University of Copenhagen, spring 2006                         *
 *---------------------------------------------------------------*/

#include "memcheck.h"
 
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "ds_ssort/ds_ssort.h"
#include "ds_ssort/common.h"
#include "shared.h"
#include "build.h"
#include "error.h"
#include "timer.h"

// Calculate skip for the suffix array
int calcSkip(ESA esa)
{
  int size = esa->size;
  int *latest = malloc(sizeof(*latest) * (MAXPSSMSIZE+1));
  if(latest == NULL)
  {
  	setError("Couldn't allocate temporary memory for skip calculation.");
  	return 0;
  }
  int i, j, minS, curLcp;
  
  // latest points to a table of where we last encountered the various possible
  // lcp values - we start by initializing it to -1 (means we never encountered
  // any of the lcp's)
  for(i = 0; i < MAXPSSMSIZE+1; i++)
    latest[i] = -1;
  
  // Last suffix has skip = size
  curLcp = getLcp(esa, size-1);
  latest[curLcp] = size-1;
  setSkip(esa, size-1,size);
  
  // Start from the second-last row and calculate skip's upwards  
  for(i = size - 2; i >= 0; i--){
    // Get lcp for this row
    curLcp = getLcp(esa, i);
    
    // Register it in table
    latest[curLcp] = i;
    
    // Find the row with the lowest index that is greater than ours and that has a 
    // strictly lower lcp 
    minS = -1;
    for(j = curLcp - 1; j >= 0; j--){			
       if(minS == -1) {
  	      minS = latest[j];
       }
       else {
  	      break;
       }
    }
 
    // then the smallest skip value is found
    for(; j >= 0; j--) {
       if(minS > latest[j]&&(latest[j]> -1)){
	      minS = latest[j];
       }
    }
    
    // If a value was found this value is registrered
    if(minS != -1) {
       setSkip(esa, i, minS);			
    } 
    // Else skip is set to the size of the text
    else {
       setSkip(esa, i, size);			
    }
  }
  free(latest);
  return 1;
}


// Calculate longest common prefix naively for the suffix array
void calcLcpNaiv(ESA esa) {
   int size = esa->size;
   unsigned char *pStr = esa->pStr;

   int i, lcp;
   unsigned char *prevSuf, *curSuf;
  
   prevSuf = pStr + getSuf(esa,0);
  
   getStr(esa)[getSize(esa)] = ENDSYM;


   for(i = 1; i < size; i++){
      curSuf = pStr + getSuf(esa,i);
      for(lcp = 0; lcp < MAXPSSMSIZE; lcp++) { 
         if(!(prevSuf[lcp] == curSuf[lcp])) {
	        break;
         } 
      }
    
      setLcp(esa, i, lcp);
      prevSuf = curSuf;
   }
}

int translate(unsigned char *pDest, char *pSrc, int strSize, char *pAlphabet, char *pIgnore)
{
  unsigned char *dest;
  unsigned char *src;	
  unsigned char *alphabet;	
  unsigned char *ignore;
  int alpha[256];
  int i, asize;

  src      = (unsigned char*)pSrc;
  dest     = pDest;
  alphabet = (unsigned char*)pAlphabet;
  ignore   = (unsigned char*)pIgnore;
  
  // Initialize the translation dictionary
  for(i = 0; i < 256; i++)
    alpha[i] = -1;
  
  // Inset the alphabet
  asize = strlen(pAlphabet);		
  for(i = 0; i < asize; i++) {
    alpha[ alphabet[i] ] = i;
  }

  // Inset the ignore alphabet
  asize = strlen(pIgnore);
  for(i = 0; i < asize; i++) {
     alpha[ ignore[i] ] = BREAKSYM;
  }

  // '\0' is ignored
  alpha[ 0 ] = BREAKSYM;

  // Translate the text
  for(i = 0; i <= strSize; i++) {
     if( alpha[src[i]] < 0) {
       char errormsg[80];
       sprintf(errormsg, "Letter '%c' found in text but not in alphabet\n", src[i]);
       setError(errormsg);
       return 0;
     }
     else
        dest[i] = alpha[ (unsigned char) src[i]];
  }
  
  return 1;
}


ESA build_ESA(char *pStr, int size, char *pAlphabet, char *pIgnore, int free_pStr) {

        // Check if the string includes a zero termination
        if(pStr[size] != '\0') {
	   setError("The string MUST include a zero termination within the size\n");
	   if(free_pStr)
	     free(pStr);
	   return NULL;
	}

	initTimer();

	int overshoot;
	ESA esa = malloc(sizeof(*esa));
	if(!esa)
	{
		setError("Couldn't allocate memory for ESA.\n");
		if(free_pStr)
		{
			free(pStr);
			freeTimer();
		}
	  	return NULL;
	}
	unsigned char *text;
	int n = size + 1; // Include the zeroterninatin in the string

	// Calculate the overshoot
	overshoot=init_ds_ssort(500,2000);	

	text = malloc((n + overshoot)*sizeof *text);	
	if(!text)
	{
		setError("Couldn't allocate memory for translated text.\n");
		free(esa);
		if(free_pStr)
		{
			free(pStr);
			freeTimer();
		}
		return NULL;
	}


	// Translate the text and stop if it fails
	if(! translate(text, pStr, n-1, pAlphabet, pIgnore) ) {
	  free(text);
	  free(esa);
	  if(free_pStr)
	    free(pStr);
	  freeTimer();
	  return NULL;
	}

	// Free pStr if possible
	if(free_pStr)
	  free(pStr);

	// Save the text, alphabet and size in the esa structure
	setStr(esa, text);
	setSize(esa, n);
	setAlphabetSize(esa, strlen(pAlphabet));
	setIgnoreAlphabetSize(esa, strlen(pIgnore));
	setAlphabet(esa, pAlphabet);
	setIgnoreAlphabet(esa, pIgnore);
	
	addTimer("Initializing");
	
	// Do the sorting, calc. lcp and calc. skip
	esa->suf = malloc(sizeof(int) * n);
	if(!esa->suf)
	{
		free(text);
		free(esa);
		freeTimer();
		setError("Couldn't allocate memory for suffix column in suffix array.\n");
		return NULL;
	}

	ds_ssort(esa->pStr, esa->suf, n, MAXPSSMSIZE);
	addTimer("DS-Sort");
	
	esa->lcp = malloc(sizeof(unsigned char) * n);	
	if(!esa->lcp)
	{
		setError("Couldn't allocate memory for LCP column in suffix array.\n");				
		free(esa->suf);
		free(text);
		free(esa);
		freeTimer();
		return NULL;
	}
	
	calcLcpNaiv(esa);
	addTimer("Calc Lcp");
	
	// The line below can be commented in to verify that there are "errors" in the suffix array
	// it will scan the array for errors and report the minimum depth at which an error was found
	// the last parameter specifies the max depth to search to).
	// As a side effect it calculates lcp (when used for this purpose the depth parameter should equa
	// that used when calling ds_ssort).
	// verifyNaively(esa, n, MAX_DEPTH);

	esa->skip = malloc(sizeof(int) * n);	
	if(!esa->skip)
	{
		setError("Couldn't allocate memory for SKIP column in suffix array.\n");						
		free(esa->lcp);
		free(esa->suf);
		free(text);
		free(esa);
		freeTimer();
		return NULL;		
	}
	
	if(calcSkip(esa) == 0)
	{
		free(esa->skip);
		free(esa->lcp);
		free(esa->suf);
		free(text);
		free(esa);
		freeTimer();
		return NULL;
	}
	addTimer("Calc Skip");
	printTimer();
	freeTimer();

	return esa;
}


ESA build_ESA_from_file(const char *pFilename, char *pAlphabet, char *pIgnore, int free_pStr) {
  char *sequences;
  FILE *file;
  int size;

  /* Open sequencefile */
  if(!(file = fopen(pFilename, "r"))) {
    char st[512];
    sprintf(st, "Could not open '%s' for reading.", pFilename);
    setError(st);
    return NULL;
  }

  /* Check file size. */
  fseek(file, 0, SEEK_END);
  size = ftell(file);    
  rewind(file);

  /* File to string */
  sequences = malloc(sizeof(char) * (size+1));
  fread(sequences, 1, size, file);
  sequences[size] = '\0';
  fclose(file);

  /* Build the ESA */
  return build_ESA(sequences, size, pAlphabet, pIgnore, free_pStr); 
}


ESA read_ESA_from_file(char *pFileName, unsigned char **ppExtraData, int *pDataRead)
{
	struct ESAFileFormat ff;
	int overshoot;	
	ESA esa = malloc(sizeof(*esa));
	if(!esa)
	{
		setError("Couldn't allocate memory for ESA.\n");		
		return NULL;
	}
	
	unsigned char *text;
	FILE *f;
	int n;
	
	f = fopen(pFileName, "r");
	if(!f)
	{
		char st[512];
		sprintf(st, "Could not open '%s' for reading.", pFileName);
		setError(st);
		free(esa);
		return NULL;
	}		
	
    if(fread(&ff, sizeof(struct ESAFileFormat), 1, f) != 1)
    {
    	setError("An error occurred reading the file.");
    	free(esa);
    	fclose(f);
    	return NULL;
    }
    
    if(strncmp(ff.ID, HEADERNAME, HEADERLENGTH) != 0)
    {
    	setError("Header name mismatch in ESA structure file.");
    	free(esa);
    	fclose(f);    	    	
		return NULL;
    }
    	
    if(ff.major != MAJOR_VERSION)
    {
    	setError("Incompatible version of the ESA structure file.");    	
    	free(esa);
    	fclose(f);
		return NULL;
    } 
    
    esa->alphabetSize = ff.alphabetSize;
    esa->alphabet = malloc( (esa->alphabetSize + 1) * sizeof(char));
    if(!esa->alphabet)
    {
    	setError("Couldn't allocate space for alphabet.");
    	free(esa);
    	fclose(f);
    	return NULL;
    }
    
    strncpy(esa->alphabet, ff.alphabet, ff.alphabetSize);
    esa->alphabet[esa->alphabetSize] = '\0';

    esa->ignoreAlphabetSize = ff.ignoreAlphabetSize;
    esa->ignoreAlphabet = malloc( (esa->ignoreAlphabetSize + 1) * sizeof(char));    
    if(!esa->ignoreAlphabet)
    {
    	setError("Couldn't allocate space for ignore alphabet.");
    	free(esa->alphabet);
    	free(esa);
    	fclose(f);
    	return NULL;   	
    }
    
    strncpy(esa->ignoreAlphabet, ff.ignoreAlphabet, ff.ignoreAlphabetSize);
    esa->ignoreAlphabet[esa->ignoreAlphabetSize] = '\0';
    
    n = esa->size = ff.size;
    
	overshoot=init_ds_ssort(500,2000);	
	text=malloc((n + overshoot)*sizeof *text);	
	
	if(!text)
	{
		setError("Couldn't allocate space for text.");
    	free(esa->alphabet);
    	free(esa->ignoreAlphabet);
    	free(esa);
    	fclose(f);
    	return NULL;   	
	}
	
    if(fread(text, sizeof(unsigned char), ff.size+1, f) != ff.size+1)
    {
    	setError("Couldn't read text from file.");
    	free(text);
    	free(esa->alphabet);
    	free(esa->ignoreAlphabet);
    	free(esa);
    	fclose(f);
    	return NULL;   	
    }
    
	esa->pStr = text;
    
	esa->suf = malloc(sizeof(int) * n);
	if(!esa->suf)
	{
		setError("Couldn't allocate memory for suf column in suffix array.");		
    	free(text);
    	free(esa->alphabet);
    	free(esa->ignoreAlphabet);
    	free(esa);
    	fclose(f);
    	return NULL;   	
	}
	
    if(fread(esa->suf, sizeof(int), ff.size, f) != ff.size)
    {
		setError("Couldn't read suffix column from file.");		
    	free(esa->suf);
    	free(text);
    	free(esa->alphabet);
    	free(esa->ignoreAlphabet);
    	free(esa);
    	fclose(f);
    	return NULL;  	
    }
    
	esa->lcp = malloc(sizeof(unsigned char) * n);	    
	if(!esa->lcp)
	{
		setError("Couldn't allocate memory for lcp column in suffix array.");		
    	free(esa->suf);
    	free(text);
    	free(esa->alphabet);
    	free(esa->ignoreAlphabet);
    	free(esa);
    	fclose(f);
    	return NULL;  	
	}
	
   if(fread(esa->lcp, sizeof(unsigned char), ff.size, f) != ff.size)
   {
   		setError("Couldn't read lcp from file.");		
   		free(esa->lcp);
    	free(esa->suf);
    	free(text);
    	free(esa->alphabet);
    	free(esa->ignoreAlphabet);
    	free(esa);
    	fclose(f);
    	return NULL;  	
   }

	esa->skip = malloc(sizeof(int) * n);	
	if(!esa->skip)
	{
		setError("Couldn't allocate space for skip column in suffix array.");		
   		free(esa->lcp);
    	free(esa->suf);
    	free(text);
    	free(esa->alphabet);
    	free(esa->ignoreAlphabet);
    	free(esa);
    	fclose(f);
    	return NULL;  	
	}
	
	if(calcSkip(esa) == 0)
	{
		free(esa->skip);
		free(esa->lcp);
    	free(esa->suf);
    	free(text);
    	free(esa->alphabet);
    	free(esa->ignoreAlphabet);
    	free(esa);
    	fclose(f);
    	return NULL;  	
	}
	
    // Read extra data if it was asked for and if it is present.
	if(ppExtraData != NULL) {
		*pDataRead = ff.nExtraData;
		if(*pDataRead != -1)
		{
			*ppExtraData = malloc(sizeof(unsigned char) * ff.nExtraData);
			if(*ppExtraData == NULL)
			{
				setError("Error - couldn't allocate space for extra data.\n");
				free(esa->skip);
				free(esa->lcp);
		    	free(esa->suf);
		    	free(text);
		    	free(esa->alphabet);
		    	free(esa->ignoreAlphabet);
		    	free(esa);
		    	fclose(f);
				return NULL;
			}
			if(fread(*ppExtraData, sizeof(unsigned char), ff.nExtraData, f) != (unsigned int)ff.nExtraData)
			{
				setError("Error reading extra data.\n");
				free(*ppExtraData);
				free(esa->skip);
				free(esa->lcp);
		    	free(esa->suf);
		    	free(text);
		    	free(esa->alphabet);
		    	free(esa->ignoreAlphabet);
		    	free(esa);
		    	fclose(f);
				return NULL;
			}
		} else {
			*ppExtraData = NULL;
		}
	}

	fclose(f);

	return esa;
}


void free_ESA(ESA esa)
{
  if(esa) {
    free(esa->suf);
    free(esa->lcp);	
    free(esa->skip);	
    free(esa->pStr);
    free(esa);
  }
}

int writeESAToFile(char *pFileName, ESA esa, unsigned char *pExtraData, int nExtra)
{
	FILE *f;	
	f = fopen(pFileName, "w");	
	if(f == 0)
	{
	        char st[512];
		sprintf(st, "Could not open '%s' for writing.", pFileName);
		setError(st);
		return 0;
	}
	
	// Initialize file format structure
	struct ESAFileFormat ff;
	memset(&ff, 0, sizeof(struct ESAFileFormat));
	strcpy(ff.ID, HEADERNAME); 
	ff.major = MAJOR_VERSION;
	ff.minor = MINOR_VERSION;
	ff.size = esa->size;	
	ff.alphabetSize = getAlphabetSize(esa);
	strncpy(ff.alphabet, getAlphabet(esa), getAlphabetSize(esa));
	ff.ignoreAlphabetSize = getIgnoreAlphabetSize(esa);
	strncpy(ff.ignoreAlphabet, getIgnoreAlphabet(esa), getIgnoreAlphabetSize(esa));
	if(pExtraData == NULL)
	{
		ff.nExtraData = -1;
	} else {
		ff.nExtraData = nExtra;
	}
	ff.reserved1 = 0; // Presently the reserved values are not used
	ff.reserved2 = 0;
	
	// Write header
	fwrite(&ff, sizeof(struct ESAFileFormat), 1, f);
	// Write translated string
	fwrite(esa->pStr, sizeof(unsigned char), ff.size+1, f);		
	// Write array
	fwrite(esa->suf, sizeof(int), ff.size, f);
	fwrite(esa->lcp, sizeof(unsigned char), ff.size, f);
	
	// Possibly write extra data
	if(pExtraData != NULL)
	{
		fwrite(pExtraData, sizeof(unsigned char), nExtra, f);
	} 	

	fclose(f);
	
	return 1;
}
