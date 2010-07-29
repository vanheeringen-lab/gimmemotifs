/*---------------------------------------------------------------*
 * PSSM Searcher                                                 *
 *                                                               *
 * By: Jes Frellsen, Ida Moltke and Martin Thiim                 *
 * University of Copenhagen, spring 2006                         *
 *---------------------------------------------------------------*/

#include "memcheck.h"

#include <assert.h>
#include <stdio.h>
#include "error.h"
#include "hittable.h"
#include "shared.h"

int g_inited = 0;
int g_overload = 0;
int g_curSize;
struct HitTable *g_pHitTable;
struct HitEntry* g_NextHitEntry;
struct HitEntry* g_EndOfEntries;

int resize_table()
{

	// If we know we are overloaded, just finish
	if(g_overload == 1)
	{
		return 0;			
	}

	// Out of space - try to extend table
	int oldSize = g_curSize;
	int newSize = (int)(g_curSize * INCREASE_FACTOR);
	struct HitEntry *pNewPointer = realloc(g_pHitTable->pScores, newSize*sizeof(struct HitEntry));

	if(pNewPointer == 0)
	{
		// Darn, no more memory. Don't update anything and instead report
		// out of memory condition														
		g_overload = 1;
		setError("Failed to increase table size.\n");
		return 0;
	}
	// Update all pointers
	g_curSize = newSize;
	g_pHitTable->pScores = pNewPointer;
	g_NextHitEntry = pNewPointer + oldSize;
	g_EndOfEntries = pNewPointer + newSize;		

	return 1;
}
	
int init_hittable()
{
	// Validate preconditions
	if(g_inited == 1)
		return 0;

	g_pHitTable = malloc(sizeof(struct HitTable));	
	if(!g_pHitTable) {
		setError("Couldn't allocate memory for hit table.");
		return 0;
	}
	g_curSize = START_SCORES;
	g_pHitTable->pScores = malloc ( sizeof(struct HitEntry) * g_curSize);	
	if(!g_pHitTable->pScores) {
		setError("Couldn't allocate memory for hit table.");
		free(g_pHitTable);
		return 0;		
	}
	g_NextHitEntry = g_pHitTable->pScores;
	g_EndOfEntries = g_pHitTable->pScores + g_curSize;
	g_overload = 0;
	g_inited = 1;	
	return 1;
}


void reset_hits() {
  g_NextHitEntry = g_pHitTable->pScores;
  g_pHitTable->nScores = 0;
}

struct HitTable* finish()
{
	if(g_inited == 0)
		return 0;
	// If we suffered memory overload, clean up the system...
	if(g_overload == 1)
	{
		free(g_pHitTable->pScores);
		free(g_pHitTable);
		g_inited = 0;
		return NULL;
	}
	g_pHitTable->nScores = g_NextHitEntry - g_pHitTable->pScores;
	g_inited = 0;	

	return g_pHitTable;
}

int release_hits(struct HitTable *pSD)
{
  if(pSD) {
    if(pSD->pScores)
      free(pSD->pScores);
    free(pSD);
  }
  return 1;
}

// Function for comparing two HitEntry structs. Returns:
//  *  1 if score of HitEntry b is larger then a, 
//  *  0 if score of HitEntry b is equal a and
//  * -1 of score of HitEntry b is smaller than b
int compare_hitentries(const void *a, const void *b)
{
  const double aScore = ((struct HitEntry*) a)->score;
  const double bScore = ((struct HitEntry*) b)->score;

  return (aScore < bScore) - (aScore > bScore);
}

int sortHitTable(struct HitTable *pData)
{
  if(!pData || !pData->pScores) {
    setError("Failed sorting hits.\n");
    return 0;
  }
  else {
    qsort(pData->pScores, pData->nScores, sizeof(struct HitEntry), compare_hitentries);
    return 1;
  }
}
