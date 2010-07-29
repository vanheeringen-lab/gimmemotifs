/*-----------------------------------------------------------------------------------*
 * PSSM Searcher                                                                     *
 *                                                                                   *
 * By: Jes Frellsen, Ida Moltke and Martin Thiim                                     *
 * University of Copenhagen, spring 2006                                             *
 *-----------------------------------------------------------------------------------*/
#include "scoring.h"
#include <stdio.h>

int main()
{
	struct HitTable *pScoreData;
	int i;
	printf("Hello, World!");
	init_hittable_system(4, 100);
	register_score(5, 1.0);
	register_score(8, 1.5);
	next_matrix();
	register_score(7, 2.0);
	next_matrix();
	next_matrix();
	register_score(3, 5.0);
	pScoreData = finish();	
	printf("nMatrices:   %d\n", pScoreData->nMatrices);	
	for(i = 0; i < 4; i++)
	{	
		int j, k;
		struct HitEntry *pEntry;
		pEntry = pScoreData->pScores[i];
		k = pScoreData->pnScores[i];
		printf("Matrix %d has %d entries\n", i, k);		
		for(j = 0; j < k; j++)
		{
			printf("Entry %d has pos %d and score %f\n", j, pEntry[j].position, (float)pEntry[j].score);
		}
	}

	return 0;
}

