/*---------------------------------------------------------------*
 * PSSM Searcher                                                 *
 *                                                               *
 * By: Jes Frellsen, Ida Moltke and Martin Thiim                 *
 * University of Copenhagen, spring 2006                         *
 *---------------------------------------------------------------*/

#include "memcheck.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "stringcon.h"
#include "error.h"

// Concatenate an array of strings to one string
concatenated concat(char **pSrc, int *lengths, int numOfStr) {
  int total, offset, length;
  int i,j;
  
  concatenated result = malloc(sizeof(*result));
  if(!result){setError("Couldn't allocate memory for concat. of strings.");return NULL;}
  
  struct HashEntry *pHashTable = malloc(sizeof(struct HashEntry) * numOfStr);   
  if(!pHashTable) { free(result); setError("Couldn't allocate memory for hash table."); return NULL;}
  memset(pHashTable, '\0', sizeof(struct HashEntry) * numOfStr);
    
  int *offsets = malloc(sizeof(int)*(numOfStr+1));  
  if(!offsets){free(result); free(pHashTable); setError("Couldn't allocate memory for concat. of strings.");return NULL;}

  char *pDest;
  // Find the total length and calculates the start of each bucket
  total = 0;
  offsets[0] = 0;
  for(i=0; i<numOfStr; i++) {
    total += (lengths[i] + 1);
    offsets[i+1] = offsets[i] + (lengths[i] + 1);
  }
  
  // Now create the hash table
  // Each hash entry represents a ``segment'' of the sequence position space
  // Ie. if we have 1.000 keys and 40.000.000 positions then key 0
  // represents the range 0-39.999 (inclusive), key 1 represents 40.000-79.999
  // etc.
  // Each entry contains the first sequence which is at least partially contained
  // in the segment represented by the entry along with the number
  // of entries partially contained in the segment. Note that the same
  // sequence may be contained in multiple entries if it crosses segment boundaries.

  // Calculate segment size.
  // Make segments one longer to compensate for case where total length is not divisble
  // by HASHTABLE_SIZE
  int SEGMENT_LENGTH = total / numOfStr;  
  if(total % numOfStr != 0) { SEGMENT_LENGTH++; }

  int curHashEntry = 0;
  int curSequence = 0;    
  while(curSequence < numOfStr)
  {
    int lowSegment = offsets[curSequence] / SEGMENT_LENGTH;
    int highSegment = (offsets[curSequence] + lengths[curSequence]) / SEGMENT_LENGTH;
    
    if((lowSegment <= curHashEntry) && (curHashEntry <= highSegment))
    {
       pHashTable[curHashEntry].numSequences++;
    }
    if(curHashEntry < highSegment)
    {
       // Open next segment, note we already cleared all 
       // memory to zero after allocation so numSequences
       // of the new entry starts at zero.
       curHashEntry++;
       pHashTable[curHashEntry].firstSequence = curSequence;       
    } else {
       curSequence++;
    }
  }  
 
  pDest = malloc(sizeof(char)*total);
  if(!pDest){free(result); free(pHashTable); free(offsets);setError("Couldn't allocate memory for concat. of strings.");return NULL;}

  // Copy the strings
  for(i=0; i<numOfStr; i++) {
    offset = offsets[i];
    length = lengths[i];
    for(j=0; j<length; j++) {
      pDest[offset+j] = pSrc[i][j];
    }
    pDest[offset+length] = '\0';
  }

  // Return the results
  result->pStr = pDest;
  result->sizeStr = total - 1;
  result->offsets = offsets;
  result->numOffsets = numOfStr;
  result->pHashTable = pHashTable;
  result->hashSegmentLength = SEGMENT_LENGTH;

  return result;
}

// Reverse lookup the string number from the position in a concatenated string
inline int reverseLookup(concatenated data, int pos) {
  int i;
  int *offsets = data->offsets;
  int segment = pos / data->hashSegmentLength;
  int first = data->pHashTable[segment].firstSequence;
  int last = first + data->pHashTable[segment].numSequences - 1;

  for(i = first; i <= last; i++) {
    if((pos >= offsets[i]) && (pos <= offsets[i+1]))
      return i;
  }

  return -1;
}

// DEPRECATED! (requires an extra 'int sequence' field in the ScoreData)
// Do the reverse lookup for a whole ScoreData-set
/*void RLScoreData(struct ScoreData scores, concatenated data) {

  int *offsets = data->offsets;
  int nScores = scores.nScores;
  struct ScoreEntry *entries = scores.pScores;
  
  int i, pos, seq;

  for(i=0; i<nScores; i++) {
    pos = entries[i].position;
    seq = reverseLookup(data, pos);
    
    entries[i].sequence = seq;
    entries[i].position = pos - offsets[seq];
  }



}*/

// Free the concatenated structure
void freeConcatenated(concatenated data) {
  if(data) {
    if(data->pStr)
      free(data->pStr);
    if(data->offsets)
      free(data->offsets);
    if(data->pHashTable) {
       free(data->pHashTable);
    }
    free(data);
  }
}


void freeConcatenatedNotData(concatenated data) {
  if(data)
    free(data);
}
