/*---------------------------------------------------------------*
 * PSSM Searcher                                                 *
 *                                                               *
 * By: Jes Frellsen, Ida Moltke and Martin Thiim                 *
 * University of Copenhagen, spring 2006                         *
 *---------------------------------------------------------------*/

#ifndef HITTABLE_H_
#define HITTABLE_H_

#include <assert.h>
#include <stdio.h>

// Struct for entries in the array of hit entries
struct HitEntry
{
	int position; // Position in sequence
	double score; // The attained score
};

// Structure representing the result set
struct HitTable
{
	int nScores;  // Number of hits registered
	struct HitEntry* pScores; // Pointer to hits
};

extern int g_inited;
extern int g_overload;
extern int g_curSize;
extern struct HitTable *g_pHitTable;
extern struct HitEntry* g_NextHitEntry;
extern struct HitEntry* g_EndOfEntries;

int init_hittable(void);
struct HitTable* finish(void);
int release_hits(struct HitTable *pData);
int resize_table(void);
int sortHitTable(struct HitTable *pData);
void reset_hits(void);

inline extern int register_hit(int pos, double score)
{
	assert(g_inited == 1);	
	if(g_NextHitEntry == g_EndOfEntries) {		
		// Try to resize table - return 0 if it fails.
		if(resize_table() == 0)
		{
		  fprintf(stderr, "Resizing\n");
			return 0;
		}
	}

	g_NextHitEntry->position = pos;
	g_NextHitEntry->score = score;
	g_NextHitEntry++;

	return 1;
}

#endif /*SCORING_H_*/
