

// Copyright 2009, 2010 Christopher Benner <cbenner@ucsd.edu>
// 
// This file is part of HOMER
//
// HOMER is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// HOMER is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.



#include "Motif.h"

Floattable* rankvalues = NULL;
int rankDescending(const void* a, const void* b) {
	float va = rankvalues->search(*(char**)a);
	float vb = rankvalues->search(*(char**)b);
	if (va < vb) return -1;
	if (va > vb) return 1;
	return 0;
}

Floattable* sorttable2 = NULL;
int compareMotifsFH2(const void* a, const void* b) {
	float val1 = sorttable2->search(*(char**)a);
	float val2 = sorttable2->search(*(char**)b);
	if (val1-10 < EMPTY_FLOAT) val1 = 100;
	if (val2-10 <  EMPTY_FLOAT) val2 = 100;
	if (val1 > val2) return 1;
	if (val1 < val2) return -1;
	return 0;
}
int compareMotifsFH2R(const void* a, const void* b) {
	float val1 = sorttable2->search(*(char**)a);
	float val2 = sorttable2->search(*(char**)b);
	if (val1-10 < EMPTY_FLOAT) val1 = -1000;
	if (val2-10 <  EMPTY_FLOAT) val2 = -1000;
	if (val1 > val2) return -1;
	if (val1 < val2) return 1;
	return 0;
}

int sortProteinStat(const void* a, const void* b) {
	float val1 = ((ProteinStat*)a)->pvalue;
	float val2 = ((ProteinStat*)b)->pvalue;
	if (val1 > val2) return 1;
	if (val1 < val2) return -1;
	return 0;
}
int sortConStat(const void* a, const void* b) {
	float val1 = (*(ConStat**)a)->pvalue;
	float val2 = (*(ConStat**)b)->pvalue;
	if (val1 > val2) return 1;
	if (val1 < val2) return -1;
	return 0;
}
int sortConStatDeg(const void* a, const void* b) {
	float val1 = (*(ConStat**)a)->dpvalue;
	float val2 = (*(ConStat**)b)->dpvalue;
	if (val1 > val2) return 1;
	if (val1 < val2) return -1;
	return 0;
}
	
int sortPSSMsPvalue(const void* a, const void* b) {
	float val1 = (*(PSSM**)a)->pvalue;
	float val2 = (*(PSSM**)b)->pvalue;
	if (val1 > val2) return 1;
	if (val1 < val2) return -1;
	return 0;
}

int sortFloat(const void*a, const void* b) {
	if (*(float*)a < *(float*)b) return -1;
	if (*(float*)b < *(float*)a) return 1;
	return 0;
}

char** getMersOfLength(CommandLine* cmd, Hashtable* mers, int length, int &numMers) {
	char** keys = mers->keys();
	int numkeys = mers->total;
	char** seq = new char*[numkeys];
	numMers = 0;
	for (int i=0;i<numkeys;i++) {
		int L = strlen(keys[i]);
		if (L == length) {
			seq[numMers++] = keys[i];
		} else {
			delete [](keys[i]);
		}
	}
	delete []keys;
	return seq;
}

void split(char* string, char** cols, int &numCols, char delim) {
    cols[0] = string;
    numCols=1;
    int len = strlen(string);
    for (int i=0;i<len;i++) {
        if (string[i] == delim) {
            string[i] = '\0';
            cols[numCols] = &(string[i+1]);
            numCols++;
        } else if (string[i] == '\n') {
            string[i] = '\0';
        } else if (string[i] == '\r') {
            string[i] = '\0';
	    }
    }
}



//code optimized for 1 seq length;
void findMotifs(CommandLine* cmd, XMerData* xmer) {



	char** dnaSeeds = NULL;
	ProteinStat* proteinSeeds = NULL;
	int numSeeds = 0;
	Inttable* seq2skip = new Inttable(MOTIF_HASH_SIZE);

	FILE* tmp = fopen(TMP_MOTIF_FILE,"w");
	fprintf(tmp, "Current Motifs:\n");
	fclose(tmp);

	//==============================================================================
	//Global Optimization
	//get seeds for local motif optimization
	if (cmd->seedfile != NULL) {
		//seeds = readMerFile(cmd, cmd->seedfile, numSeeds);
	} else {
		if (cmd->seqType == DNA_SEQ_TYPE) {
			dnaSeeds = findBestSeeds(cmd,xmer,NULL,0,numSeeds);
		} else {
			proteinSeeds = (ProteinStat*) findBestSeeds(cmd,xmer,NULL,0,numSeeds);
		}
	}
	int motifLength = cmd->merLengths[0];
	//===============================================================================
	//local optimization

	PSSM** pssm = new PSSM*[cmd->numSeeds];
	int numpssms = 0;
	int curNumMers = xmer->numXMers;
	fprintf(stderr, "------ If this is taking too long, try using the binomial instead of hypergeometric ------\n");
	for (int i=0;i<numSeeds;i++) {
		//check for redundancy
		//make sure seeds are in common with other seeds
		//
		PSSM* curPSSM = NULL;
		if (cmd->seqType == DNA_SEQ_TYPE) {
			if (checkforSimilarSequence(cmd,dnaSeeds[i],xmer,seq2skip)) {
				fprintf(stderr, "Skipping: %s\n", dnaSeeds[i]);
				continue;
			} else {
				fprintf(stderr, "Locally Optimizing Motif %d, Seed %d - %s\n", 
								numpssms+1,i+1, dnaSeeds[i]);
				curPSSM = new PSSM(dnaSeeds[i],cmd->seqType);
			}
		} else {
			char* str = getAAStringfromIndex(proteinSeeds[i].v);
			if (checkforSimilarSequence(cmd,proteinSeeds[i].v,xmer,seq2skip)) {
				fprintf(stderr, "Skipping: %s\n", str);
				delete []str;
				continue;
			} else {
				fprintf(stderr, "Locally Optimizing Motif %d, Seed %d - %s\n", 
								numpssms+1,i+1, str);
				delete []str;
				curPSSM = new PSSM(proteinSeeds[i].v,cmd->seqType);
			}
		}

		//fprintf(stderr, "curNumMers = %d of %d\n", curNumMers,numMers);
		//includedSeq = calculateMerPSSMscore(cmd, motif, xmer, mers, numMers, 1);
		int curIndex = localOptimizeMerPSSM(cmd,curPSSM,
					xmer,xmer->stats,curNumMers);
		pssm[numpssms] = curPSSM;

		for (int j=curIndex;j<curNumMers;j++) {
			seq2skip->insert(1,xmer->mers[xmer->stats[j]->id]);
		}
		curNumMers = curIndex;

		tmp = fopen(TMP_MOTIF_FILE, "a");
		fprintf(tmp, "Seed: %d\n", i+1);
		pssm[numpssms]->setGapInfo(cmd->gapLengths, cmd->numGapLengths);
		pssm[numpssms]->expXform();
		pssm[numpssms]->print(tmp);
		fclose(tmp);

		numpssms++;
		if (numpssms >= cmd->numSeeds) break;
	}
	qsort(pssm, numpssms, sizeof(PSSM*),sortPSSMsPvalue);

	char* filename = new char[1000];
	sprintf(filename, "%s.motifs%d", cmd->outputfile, motifLength);
	FILE* fp = fopen(filename, "w");
	for (int j=0;j<numpssms;j++) {
		pssm[j]->setGapInfo(cmd->gapLengths, cmd->numGapLengths);
		pssm[j]->expXform();
		pssm[j]->print(fp);
		delete pssm[j];		
	}
	delete []pssm;
	fclose(fp);
	fprintf(stderr, "Done finding motifs for mersize = %d (%s)\n",motifLength,filename);
}


//wrapper for DNA
int checkforSimilarSequence(CommandLine* cmd, char* seq, XMerData* xmer,
						Inttable* seq2skip) {
	int* values = NULL;
	values = seqIndexOfIUPAC(cmd->seqType,seq,NULL);
	int rv = checkforSimilarSequence(cmd,values,xmer,seq2skip);
	delete []values;
	return rv;
}

int checkforSimilarSequence(CommandLine* cmd, int* values, XMerData* xmer,
						Inttable* seq2skip) {

	static int* misMatchMers = NULL;
	if (misMatchMers == NULL) misMatchMers = new int[MAX_MISMATCHMERS];
	int numMisMatchMers = 0;
	xmer->seqTree->findMisMatchMers(values, 0, 0, misMatchMers, numMisMatchMers);

	int numFound = 0;
	for (int i=0;i<numMisMatchMers;i++) {
		if (EMPTY_INT != seq2skip->search(xmer->mers[misMatchMers[i]])) {
			numFound++;
		}
	}
	if (numFound / numMisMatchMers > SKIPMER_PERCENT_SKIP) {
		return 1;
	} else {
		return 0;
	}
}

XMerData* readMerFile(CommandLine* cmd) {

	Hashtable* ht = new Hashtable(MOTIF_HASH_SIZE);

	FILE* fp = fopen(cmd->merfile, "r");
	if (fp == NULL) {
		fprintf(stderr, "Could not open mer file\n");
		exit(0);
	}
	int *lengths = new int[100];
	char** cols = new char*[10000];
	int numCols = 0;

	int numLengths = 0;
	char buf[BSIZE];
	char *m = NULL;
	int lineNum = 0;
	float numGenes = 0;
	float numPos =0;
	while (fgets(buf,BSIZE,fp) != NULL) {
		lineNum++;
		if (lineNum % 100000 == 0) fprintf(stderr, "%d\n", lineNum);
		float score=0;
		float A = 0;
		float P = 0;

		split(buf, cols, numCols,'\t');
		m = cols[0];
		if (numCols > 1) sscanf(cols[1],"%f",&score);
		if (numCols > 2) sscanf(cols[2],"%f",&numGenes);
		if (numCols > 3) sscanf(cols[3],"%f",&numPos);
		if (numCols > 4) sscanf(cols[4],"%f",&A);
		if (numCols > 5) sscanf(cols[5],"%f",&P);

		int seqlen = strlen(m);
		if (cmd->minmer > seqlen || cmd->maxmer < seqlen) continue;
		int bad = 0;
		for (int k=0;k<numLengths;k++) {
			if (lengths[k] == seqlen) {
				bad = 1;
				break;
			}
		}
		if (!bad) {
			lengths[numLengths] = seqlen;
			numLengths++;
		}
		ConStatApprox* cs = new ConStatApprox();
		cs->fp = P;
		cs->fn = A-P;
		cs->pvalue = score;
		ht->insert(cs,m);
	}
	XMerData* xmer = new XMerData();
	xmer->loadMerData(ht);
	xmer->enumSites = numGenes;
	xmer->enumPosSites = numPos;
	xmer->enumNegSites = numGenes-numPos;
	xmer->numGenes = (int) numGenes;
	xmer->numPosGenes = (int) numPos;
	xmer->numNegGenes = (int) (numGenes-numPos);
	xmer->enumGenes = numGenes;
	xmer->enumPosGenes = numPos;
	xmer->enumNegGenes = numGenes-numPos;
	xmer->initZOOPSCache();

	if (cmd->algorithm == ALG_HYPERGEO) {
		cmd->statmemory = new StatMemory(MAX_STAT_CACHE, (int)round(xmer->enumGenes),
					(int)round(xmer->enumPosGenes));
	}


	fprintf(stderr, "\nCalculating statistics for %d mers\n", xmer->numXMers);
	for (int i=0;i<xmer->numXMers;i++) {
		//float pw = 0.0,nw=0.0;


		ConStatApprox* cs = (ConStatApprox*) xmer->stats[i];
		float ogp = 0, ogn = 0;
		if (cmd->zoopsApprox) {
			//This allows us to approximate additional motifs
			ogp = cs->fp;
			ogn = cs->fn;
			int ep = (int)round(cs->fp*WEIGHT_RES);
			int en = (int)round(cs->fn*WEIGHT_RES);
			cs->fp = xmer->papprox[ep] / ((float)WEIGHT_RES);
			cs->fn = xmer->napprox[en] / ((float)WEIGHT_RES);
		}
		if (cmd->algorithm == ALG_HYPERGEO) {
			cs->pvalue = cmd->statmemory->getStat(
						(int)round(cs->fn+cs->fp),(int)round(cs->fp));
			/*cs->pvalue = loghypergeo((int)round(xmer->enumGenes),(int)round(xmer->enumPosGenes),
							(int)round(cs->fn+cs->fp),(int)round(cs->fp));
			*/
		} else if (cmd->algorithm == ALG_APPROXBINOMIAL) {
			float r = ((float)cs->fn)/(xmer->enumNegSites);
			//fprintf(stderr, "%d\t%d\t%f\n", numPosSites, cs->p, r);
			cs->pvalue = logbinomial((int)round(xmer->enumPosSites), (int)round(cs->fp), r,
									(int)round(xmer->enumNegSites));
		} else if (cmd->algorithm == ALG_SITEBINOMIAL) {
			float r = ((float)cs->fn)/(xmer->enumNegSites);
			//fprintf(stderr, "%d\t%d\t%f\n", numPosSites, cs->p, r);
			cs->pvalue = logbinomial((int)round(xmer->enumPosSites), (int)round(cs->fp), r,
									(int)round(xmer->enumNegSites));
		} else if (cmd->algorithm == ALG_SITEHYPERGEO) {
			cs->pvalue = loghypergeo((int)round(xmer->enumSites),(int)round(xmer->enumPosSites),
								(int)round(cs->fn+cs->fp), (int)round(cs->fp));
		} else if (cmd->algorithm == ALG_BINOMIAL) {
			float r = cs->fn / xmer->enumNegGenes;
			cs->pvalue = logbinomial((int)round(xmer->enumPosGenes), 
										(int)round(cs->fp), r, 
										(int)round(xmer->enumNegGenes));
		} else if (cmd->algorithm == ALG_FREQDIFF) {
			cs->pvalue = (cs->fn/xmer->enumNegSites) - (cs->fp/xmer->enumPosSites);
			//cs->pvalue = (cs->fn/xmer->enumNegGenes) - (cs->fp/xmer->enumPosGenes);
		}
		if (cmd->zoopsApprox) {
			cs->fp = ogp;
			cs->fn = ogn;
		}
		xmer->stats[i]->dpvalue = xmer->stats[i]->pvalue;
	}
	fprintf(stderr, "done\n");

	//cmd->numMerLengths = numLengths;
	//cmd->merLengths = lengths;
	delete []cols;
	fclose(fp);
	delete ht;
	return xmer;
}


int localOptimizeMerPSSM(CommandLine* cmd, PSSM* pssm, XMerData* xmer,
				ConStat** mers, int numMers) {
	pssm->adjustvalues();
	pssm->logXform();
	float lastPvalue = 10000;
	pssm->threshold = 100;

	for (int i=0;i<cmd->maxOptIterations;i++) {
		char* conString = pssm->consensus();
		fprintf(stderr, "optimizing %s, iteration %d ...\n", conString, i+1);
		delete []conString;

		(void)calculateMerPSSMscore(cmd, pssm,xmer,
							mers, numMers, 2);
		pssm->adjustvalues();
		pssm->expXform();
	//	pssm->print(stderr);
		if (pssm->pvalue >= lastPvalue) break;
		lastPvalue = pssm->pvalue;
	}
	return calculateMerPSSMscore(cmd, pssm,xmer,mers, numMers,1);
}

char* seqIndexOf(int seqType, char* mer, char* values) {
	int len = strlen(mer);
	if (values == NULL) {
		values = new char[len+1];
	}
	for (int i=0;i<len;i++) {
		if (seqType == DNA_SEQ_TYPE) {
			switch(mer[i]) {
				case 'A':
					values[i] = 0;
					break;
				case 'C':
					values[i] = 1;
					break;
				case 'G':
					values[i] = 2;
					break;
				case 'T':
					values[i] = 3;
					break;
				default:
					values[i] = -1;
			}
		} else if (seqType == PROTEIN_SEQ_TYPE) {
			switch(mer[i]) {
				case 'A':
					values[i] = 0;
					break;
				case 'C':
					values[i] = 1;
					break;
				case 'D':
					values[i] = 2;
					break;
				case 'E':
					values[i] = 3;
					break;
				case 'F':
					values[i] = 4;
					break;
				case 'G':
					values[i] = 5;
					break;
				case 'H':
					values[i] = 6;
					break;
				case 'I':
					values[i] = 7;
					break;
				case 'K':
					values[i] = 8;
					break;
				case 'L':
					values[i] = 9;
					break;
				case 'M':
					values[i] = 10;
					break;
				case 'N':
					values[i] = 11;
					break;
				case 'P':
					values[i] = 12;
					break;
				case 'Q':
					values[i] = 13;
					break;
				case 'R':
					values[i] = 14;
					break;
				case 'S':
					values[i] = 15;
					break;
				case 'T':
					values[i] = 16;
					break;
				case 'V':
					values[i] = 17;
					break;
				case 'W':
					values[i] = 18;
					break;
				case 'Y':
					values[i] = 19;
					break;
				default:
					values[i] = -1;
					break;
			}
		}
	}
	values[len] = -1;
	return values;
}
char* getAAStringfromIndex(int* values) {
	if (values == NULL) return NULL;
	int len =0;
	int index =0;
	while (values[index] != 0) {
		index++;
	}
	len = index+1;
	char* str = new char[len];
	for (int i=0;i<index;i++) {
		switch(values[i]) {
			case 1:
				str[i] = 'A';
				break;
			case 2:
				str[i] = 'C';
				break;
			case 4:
				str[i] = 'D';
				break;
			case 8:
				str[i] = 'E';
				break;
			case 16:
				str[i] = 'F';
				break;
			case 32:
				str[i] = 'G';
				break;
			case 64:
				str[i] = 'H';
				break;
			case 128:
				str[i] = 'I';
				break;
			case 256:
				str[i] = 'K';
				break;
			case 512:
				str[i] = 'L';
				break;
			case 1024:
				str[i] = 'M';
				break;
			case 2048:
				str[i] = 'N';
				break;
			case 4096:
				str[i] = 'P';
				break;
			case 8192:
				str[i] = 'Q';
				break;
			case 16384:
				str[i] = 'R';
				break;
			case 32768:
				str[i] = 'S';
				break;
			case 65536:
				str[i] = 'T';
				break;
			case 131072:
				str[i] = 'V';
				break;
			case 262144:
				str[i] = 'W';
				break;
			case 524288:
				str[i] = 'Y';
				break;
			default:
				str[i] = 'X';
				break;
		}
				
	}
	str[index] = '\0';
	return str;
}
char* getStringFromIUPACIndex(int* mer) {
	if (mer == NULL) return NULL;
	int len =0;
	int index =0;
	while (mer[index] != 0) {
		index++;
	}
	len = index+1;
	char* str = new char[len];
	for (int i=0;i<index;i++) {
		switch(mer[i]) {
			case 1:
				str[i] = 'A';
				break;
			case 2:
				str[i] = 'C';
				break;
			case 3:
				str[i] = 'M';
				break;
			case 4:
				str[i] = 'G';
				break;
			case 5:
				str[i] = 'R';
				break;
			case 6:
				str[i] = 'S';
				break;
			case 7:
				str[i] = 'V';
				break;
			case 8:
				str[i] = 'T';
				break;
			case 9:
				str[i] = 'W';
				break;
			case 10:
				str[i] = 'Y';
				break;
			case 11:
				str[i] = 'H';
				break;
			case 12:
				str[i] = 'K';
				break;
			case 13:
				str[i] = 'D';
				break;
			case 14:
				str[i] = 'B';
				break;
			case 15:
				str[i] = 'N';
				break;
			default:
				str[i] = 'N';
		}
	}
	str[index] = '\0';
	return str;
}

// The advantage of this funcion is for representing IUPAC symbols...
int* seqIndexOfIUPAC(int seqType, char* mer, int* values) {
	int len = strlen(mer);
	if (values == NULL) {
		values = new int[len+1];
	}
	for (int i=0;i<len;i++) {
		if (seqType == DNA_SEQ_TYPE) {
			switch(mer[i]) {
				case 'A':
					values[i] = 1;
					break;
				case 'C':
					values[i] = 2;
					break;
				case 'G':
					values[i] = 4;
					break;
				case 'T':
					values[i] = 8;
					break;
				case 'M':
					values[i] = 3;
					break;
				case 'R':
					values[i] = 5;
					break;
				case 'W':
					values[i] = 9;
					break;
				case 'S':
					values[i] = 6;
					break;
				case 'Y':
					values[i] = 10;
					break;
				case 'K':
					values[i] = 12;
					break;
				case 'V':
					values[i] = 7;
					break;
				case 'H':
					values[i] = 11;
					break;
				case 'D':
					values[i] = 13;
					break;
				case 'B':
					values[i] = 14;
					break;
				case 'N':
					values[i] = 15;
					break;
				default:
					values[i] = 15;
			}
		} else if (seqType == PROTEIN_SEQ_TYPE) {
			switch(mer[i]) {
				case 'A':
					values[i] = 1;
					break;
				case 'C':
					values[i] = 2;
					break;
				case 'D':
					values[i] = 4;
					break;
				case 'E':
					values[i] = 8;
					break;
				case 'F':
					values[i] = 16;
					break;
				case 'G':
					values[i] = 32;
					break;
				case 'H':
					values[i] = 64;
					break;
				case 'I':
					values[i] = 128;
					break;
				case 'K':
					values[i] = 256;
					break;
				case 'L':
					values[i] = 512;
					break;
				case 'M':
					values[i] = 1024;
					break;
				case 'N':
					values[i] = 2048;
					break;
				case 'P':
					values[i] = 4096;
					break;
				case 'Q':
					values[i] = 8192;
					break;
				case 'R':
					values[i] = 16384;
					break;
				case 'S':
					values[i] = 32768;
					break;
				case 'T':
					values[i] = 65536;
					break;
				case 'V':
					values[i] = 131072;
					break;
				case 'W':
					values[i] = 262144;
					break;
				case 'Y':
					values[i] = 524288;
					break;
				case 'X':
					values[i] = 1048576;
					break;
				default:
					values[i] = 1048576;
					break;
			}
		}
	}
	values[len] = 0;
	return values;
}

// eact flag basically says that the objects held by the mers hashtable are FullConStat
// objects which hold information about which genes each mer is present in

//returns the new number of mers, were mers is sorted from lowest scoring
// to highest scoing with the mers included in the motif at the end
// the returned value marks the index were the motif starts
int calculateMerPSSMscore(CommandLine* cmd, PSSM* pssm, XMerData* xmer, 
				ConStat** mers, int numMers, int updateFlag) {


	//pssm->adjustvalues();
	pssm->logXform();
	char* conString = pssm->consensus();
	char* seqIndex = seqIndexOf(cmd->seqType,conString,NULL);
	delete []conString;

	int size = 4;
	if (cmd->seqType == PROTEIN_SEQ_TYPE) size = 20;
//fprintf(stderr, "Finding Similarity of PSSM to sequence\n");
	for (int i=0;i<numMers;i++) {
		seqIndex = seqIndexOf(cmd->seqType,xmer->mers[mers[i]->id],seqIndex);
		mers[i]->dpvalue = calculateProbabilityPSSM(cmd, pssm, seqIndex, mers[i]->dir, 1);

	}
//fprintf(stderr, "Sorting similarity\n");
	qsort(mers, numMers, sizeof(ConStat*), sortConStatDeg);
//fprintf(stderr, "Finding similarity threshold\n");

	//This will return a sorted list of ConStats and the index that separates those
	// bound by the motifs from the others
	if (updateFlag == 0) {
		delete []seqIndex;
		//find threshold index
		for (int i=numMers-1;i>=0;i--) {
			if (mers[i]->dpvalue < pssm->threshold) {
				return i;
			}
		}
		return 0;
	}

	float numGenes = xmer->enumGenes;
	float numPosGenes = xmer->enumPosGenes;
	float numNegGenes = xmer->enumNegGenes;

//fprintf(stderr, "numGenes = %f\t%f\n", numPosGenes, numNegGenes);
	
	if (cmd->algorithm == ALG_SITEHYPERGEO || cmd->algorithm == ALG_SITEBINOMIAL ||
			cmd->algorithm == ALG_FREQDIFF) {
		numGenes = xmer->enumSites;
		numPosGenes = xmer->enumPosSites;
		numNegGenes = xmer->enumNegSites;
	}
	int plen = pssm->length;
	float threshLowerBound = plen*4.0/3.0;
	float ngenes=0.0,pgenes=0.0,bestN=0.0,bestP=0.0,lastN=0.0,lastP=0.0;
	float bestPvalue = 1.0,lastPvalue = 1.0,pvalue = 1.0;
	float bestThreshold = 0.0,lastThreshold=0.0,threshold=0.0;
	int bestIndex = 0, lastIndex = 0;

	
	int* geneHash = new int[xmer->numGenes+1];
	for (int i=0;i<xmer->numGenes;i++) {
		geneHash[i] = 0;
	}

	//int totalPossibleRank = (numGenes*(numGenes+1))/2;
	for (int i=numMers-1;i>=0;i--) {
		if (cmd->exactFlag) {
			ConStatExact* cs = (ConStatExact*)mers[i];
			threshold = cs->dpvalue;
			for (int j=0;j<cs->p;j++) {
				gene_id_t gid = cs->pgenes[j];
				if (geneHash[gid] == 1) continue;
				geneHash[gid] = 1;
				pgenes += xmer->geneWeights[gid];
			}
			for (int j=0;j<cs->n;j++) {
				gene_id_t gid = cs->ngenes[j];
				if (geneHash[gid] == 1) continue;
				geneHash[gid] = 1;
				ngenes += xmer->geneWeights[gid];
			}
		} else {
			ConStatApprox* cs = (ConStatApprox*)mers[i];
			threshold = cs->dpvalue;
			if (cmd->zoopsApprox) {
				pgenes += cs->fp;
				ngenes += cs->fn;
			} else {
				pgenes += (numPosGenes-pgenes)*cs->fp/numPosGenes;
				ngenes += (numNegGenes-ngenes)*cs->fn/numNegGenes;
			}
		}

		float ogp = 0.0;
		float ogn = 0.0;
		if (cmd->zoopsApprox) {
			ogp = pgenes;
			ogn = ngenes;
			if (pgenes >= ZOOPS_APPROX_CACHE) {
				pgenes = ZOOPS_APPROX_CACHE-1;
			}
			if (ngenes >= ZOOPS_APPROX_CACHE) {
				ngenes = ZOOPS_APPROX_CACHE-1;
			}
			pgenes = xmer->papprox[(int)round(WEIGHT_RES*pgenes)]/WEIGHT_RES;
			//fprintf(stderr, "zoops = %f\n", pgenes);
			ngenes = xmer->napprox[(int)round(WEIGHT_RES*ngenes)]/WEIGHT_RES;
		}

		pvalue = 1.0;	
		if (cmd->algorithm == ALG_HYPERGEO || cmd->algorithm == ALG_SITEHYPERGEO) {
			pvalue = cmd->statmemory->getStat(
								(int)round(pgenes+ngenes),(int)round(pgenes));
			/*pvalue = loghypergeo((int)round(numGenes),(int)round(numPosGenes),
								(int)round(pgenes+ngenes),(int)round(pgenes));*/
		} else if (cmd->algorithm == ALG_SITEBINOMIAL 
										|| cmd->algorithm == ALG_BINOMIAL) {                    
			float r = ngenes / numNegGenes;
			pvalue = logbinomial((int)round(numPosGenes),(int)round(pgenes),
								r, (int)round(numNegGenes));
		} else if (cmd->algorithm == ALG_RANK) {
			//rank stuff - if it will ever be worth it.
		} else if (cmd->algorithm == ALG_FREQDIFF) {
			pvalue = ngenes/numNegGenes - pgenes/numPosGenes;
			//pvalue = ngenes/numNegGenes - pgenes/numPosGenes;
		}


		if (i < numMers-2 && ngenes/numNegGenes > cmd->maxNegGenePercentage) break;

		if (cmd->zoopsApprox) {
			pgenes = ogp;
			ngenes = ogn;
		}

		//if genes have the same score we need to make sure we 
		// count all of them together
		if (i != numMers-1 && (threshold < lastThreshold) 
								&& (lastPvalue < bestPvalue)) {
			bestPvalue = lastPvalue;
			bestN = lastN;
			bestP = lastP;
			bestIndex = lastIndex;
			bestThreshold = lastThreshold;
		}
		lastPvalue = pvalue;
		lastN = ngenes;
		lastP = pgenes;
		lastIndex = i;
		lastThreshold = threshold;
	}


//fprintf(stderr, "%d/%d\n", numMers-bestIndex, numMers);
	

//fprintf(stderr, "Regenerating PSSM...\n");
	char outstr[1000];
	if (cmd->algorithm != ALG_RANK) {
		sprintf(outstr, "%.1f,%.1f,%.1f,%.1f,%.2e", numGenes,numPosGenes,bestP+bestN,
									bestP,exp(bestPvalue));
	} else {
		//sprintf(outstr, "%d,%d,%f", numGenes, intP, bestPvalue);
	}
	pssm->setOutputString(outstr);

	int minIndex = numMers-1;
	for (int i=numMers-1;i>=0;i--) {
		float curthreshold = mers[i]->dpvalue;
		if (curthreshold > bestThreshold-threshLowerBound 
						&& (float)(numMers-minIndex)/(float)numMers < 0.15) {
			minIndex--;
		} else {
			break;
		}
	}

	//adjust threshold to ensure it gets that last xmer
	bestThreshold -= 0.0001;
	pssm->threshold = bestThreshold;
	pssm->pvalue = bestPvalue;

	if (updateFlag == 1) {
		delete []seqIndex;
		delete []geneHash;
		return bestIndex;
	}

	fprintf(stderr, ">> Number above threshold %d - %.2f\t%s\n", numMers-bestIndex,
					bestPvalue,pssm->outputString);


	//initialize matrix that we will use
	float** matrix = new float*[pssm->length];
	for (int i=0;i<pssm->length;i++) {
		matrix[i] = new float[size];
		for (int j=0;j<size;j++) {
			matrix[i][j] = 0.0;
		}
	}


	if (cmd->exactFlag) {
		for (int i=0;i<xmer->numGenes;i++) {
			geneHash[i] = 0;
		}
		for (int i=numMers-1;i>=bestIndex;i--) {
			ConStatExact* cs = (ConStatExact*)mers[i];
			for (int j=0;j<cs->p;j++) {
				geneHash[cs->pgenes[j]]++;
			}
			for (int j=0;j<cs->n;j++) {
				geneHash[cs->ngenes[j]]++;
			}
		}
	}

	float branchSize=cmd->branchSize; //normally 1
	PSSM** testPSSMs = new PSSM*[(int)(100.0/branchSize)]; //100 is the legitimet max in this case
	int numTestPSSMs = 0;
		

	int lessThanLastOne = 0;
	float lastCheckT = mers[numMers-1]->dpvalue;
	float lastT = lastCheckT+1.0001;
	char tempIndex[1000];


	ngenes=0.0,pgenes=0.0,lastN=0.0,lastP=0.0;
	//float bestPvalue = 1.0,lastPvalue = 1.0,pvalue = 1.0;
	//float bestThreshold = 0.0,lastThreshold=0.0,threshold=0.0;



	for (int i=numMers-1;i>=minIndex;i--) {
		float score = mers[i]->dpvalue;
		if (lastT > score+1e-5 && lastT < lastCheckT-branchSize && i != numMers-1) {
			PSSM* trial = pssm->copy();
			trial->expXform();
			for (int j=0;j<plen;j++) {
				for (int k=0;k<size;k++) {
					trial->matrix[j][k] = matrix[j][k];
				}
			}
			trial->adjustvalues();
			testPSSMs[numTestPSSMs++] = trial;

			lastCheckT = lastT;
		}
	
		//score based on deleting mer from motif
		float p = 0,lp = 0;
		float P=0,N=0;
		float mP=0.0,mN=0.0;
		if (cmd->exactFlag) {
			ConStatExact* cs = (ConStatExact*)mers[i];
			for (int j=0;j<cs->p;j++) {
				int gid = cs->pgenes[j];
				(void)geneHash[gid];
				//int a = geneHash[gid];
				P+= xmer->geneWeights[gid];
				/*if (a == 0) {
				} else {
					P+= 1.0/((float)geneHash[gid])
								*xmer->geneWeights[gid];
				}*/
			}
			for (int j=0;j<cs->n;j++) {
				int gid = cs->ngenes[j];
				(void)geneHash[gid];
				//int a = geneHash[gid];
				N+= xmer->geneWeights[gid];
				/*if (a == 0) {
				} else {
					N+= 1.0/((float)geneHash[gid])*
								xmer->geneWeights[gid];
				}*/
			}
			if (score < bestThreshold) {
				mP = bestP + P;
				mN = bestN + N;
			} else {
				mP = bestP - P;
				mN = bestN - N;
			}
		} else {
			ConStatApprox* cs = (ConStatApprox*)mers[i];
			if (score < bestThreshold) {
				if (cmd->zoopsApprox) {
					mP = bestP + cs->fp;
					mN = bestN + cs->fn;
				} else {
					mP = bestP + (numPosGenes-bestP)
								*(cs->fp/numPosGenes);
					mN = bestN + (numNegGenes-bestN)
								*(cs->fn/numNegGenes);
				}
			} else {
				if (cmd->zoopsApprox) {
					mP = bestP - cs->fp;
					mN = bestN - cs->fn;
					/*if (mP >= ZOOPS_APPROX_CACHE) {
						mP = ZOOPS_APPROX_CACHE-1;
					}
					if (mN >= ZOOPS_APPROX_CACHE) {
						mN = ZOOPS_APPROX_CACHE-1;
					}*/
					//mP = xmer->papprox[(int)round(mP*WEIGHT_RES)] / WEIGHT_RES;
					//mN = xmer->napprox[(int)round(mN*WEIGHT_RES)] / WEIGHT_RES;
				} else {
					mP = bestP*(1.0 - cs->fp/numPosGenes);
					mN = bestN*(1.0 - cs->fn/numNegGenes);
				}
			}
		}
		if (mN < 0) mN = 0.0;
		if (mP < 0) mP = 0.0;

		if (cmd->zoopsApprox) {
			mP = xmer->papprox[(int)round(mP*WEIGHT_RES)] / WEIGHT_RES;
			mN = xmer->napprox[(int)round(mN*WEIGHT_RES)] / WEIGHT_RES;
		}

		if (cmd->algorithm == ALG_HYPERGEO) {
			p = cmd->statmemory->getStat(
							(int)round(mP+mN), (int)round(mP));
			/*p = loghypergeo((int)round(numGenes), (int)round(numPosGenes),
							(int)round(mP+mN), (int)round(mP));*/
		} else if (cmd->algorithm == ALG_SITEHYPERGEO) {
			p = loghypergeo((int)round(numGenes), (int)round(numPosGenes),
							(int)round(mP+mN),(int)round(mP));
		} else if (cmd->algorithm == ALG_SITEBINOMIAL) {
			float r = mN/numNegGenes;
			p = logbinomial((int)round(numPosGenes), (int)round(mP), r, 
							(int)round(numNegGenes));
		} else if (cmd->algorithm == ALG_BINOMIAL) {
			float r = mN /numNegGenes;
			p = logbinomial((int)round(numPosGenes), (int)round(mP), r, 
							(int)round(numNegGenes));
		} else if (cmd->algorithm == ALG_FREQDIFF) {
			p = mN/numNegGenes - mP/numPosGenes;
		}

		lp = (bestPvalue-p);
		if (score > threshold) {
			lp *= -1;
		}

		seqIndex = seqIndexOf(cmd->seqType,xmer->mers[mers[i]->id],seqIndex);
		if (lp > 0) {
			if (cmd->dualMotifs) {
				int slen = plen/2;
				int d = mers[i]->dir;
				if (d & 2) {
					if (d & 4) {
						for (int j=0;j<slen;j++) {
							tempIndex[j] = size-1-seqIndex[slen+slen-1-j];
						}
					} else {
						for (int j=0;j<slen;j++) {
							tempIndex[j]=seqIndex[slen+j];
						}
					}
					if (d & 8) {
						for (int j=0;j<slen;j++) {
							tempIndex[slen+j] = size-1-seqIndex[slen-1-j];
						}
					} else {
						for (int j=0;j<slen;j++) {
							tempIndex[slen+j]=seqIndex[j];
						}
					}
				} else {
					if (d & 4) {
						for (int j=0;j<slen;j++) {
							tempIndex[j] = size-1-seqIndex[slen-1-j];
						}
					} else {
						for (int j=0;j<slen;j++) {
							tempIndex[j]=seqIndex[j];
						}
					}
					if (d & 8) {
						for (int j=0;j<slen;j++) {
							tempIndex[slen+j] = size-1-seqIndex[2*slen-1-j];
						}
					} else {
						for (int j=0;j<slen;j++) {
							tempIndex[slen+j]=seqIndex[slen+j];
						}
					}
				}
				for (int j=0;j<plen;j++) {
					matrix[j][(int)tempIndex[j]] += lp;
				}
				
			} else {
				if (mers[i]->dir) {
					for (int j=0;j<plen;j++) {
						matrix[plen-1-j][size-1-seqIndex[j]] += lp;
					}
				} else {
					for (int j=0;j<plen;j++) {
						matrix[j][(int)seqIndex[j]] += lp;
					}
				}
			}
		}
		
		lastT = score;
	}


	ConStat** shortMers = &(mers[minIndex]);
	int numShortMers = numMers - minIndex;
	PSSM* bestPSSM = pssm;
	float lastOne = pssm->pvalue;
	int bestShortIndex = bestIndex;
	for (int i=0;i<numTestPSSMs;i++) {

		//reset the mers variable to start part way into the array
		// that way we don't have to reallocate anything and
		// wee can still limit the number of mers used
		// to change the PSSM

		int curIndex = calculateMerPSSMscore(cmd, testPSSMs[i],
				xmer, shortMers, numShortMers, 1);
		testPSSMs[i]->expXform();

		fprintf(stderr, "%d\t%.2e\t%.2f\t%s\n",
				i+1,testPSSMs[i]->pvalue,exp(testPSSMs[i]->pvalue),
							testPSSMs[i]->outputString);

		if (testPSSMs[i]->pvalue > lastOne) {
			lessThanLastOne++;
		} else {
			lessThanLastOne = 0;
		}
		lastOne = testPSSMs[i]->pvalue;
		if (testPSSMs[i]->pvalue < bestPSSM->pvalue) {
			if (bestPSSM != pssm) {
				delete bestPSSM;
			}
			bestPSSM = testPSSMs[i];
			bestShortIndex = curIndex+minIndex;
		} else {
			delete testPSSMs[i];
		}

		if (lessThanLastOne > MAX_DECREASE) {
			fprintf(stderr, "less than last one!!!\n");
			break;
		}
	}

	if (bestPSSM != pssm) {
		for (int i=0;i<plen;i++) {
			for (int j=0;j<size;j++) {
				pssm->matrix[i][j] = bestPSSM->matrix[i][j];
			}
		}
		pssm->threshold = bestPSSM->threshold;
		pssm->state = bestPSSM->state;
		pssm->pvalue = bestPSSM->pvalue;
		pssm->score = bestPSSM->score;
		pssm->constatID = bestPSSM->constatID;
	}

	for (int i=0;i<pssm->length;i++) {
		delete [](matrix[i]);
		matrix[i] = NULL;
	}
	delete matrix;
	delete []seqIndex;
	delete []testPSSMs;
	delete []geneHash;
	return bestShortIndex;
}

void writeMerFile(CommandLine* cmd, XMerData* xmer, FILE* fp) {

	//int numGenes = (int)round(xmer->enumGenes);
	//int numPosGenes = (int)round(xmer->enumPosGenes);
	if (cmd->algorithm == ALG_SITEHYPERGEO || cmd->algorithm == ALG_SITEBINOMIAL) {
		//numGenes = (int)round(xmer->enumSites);
		//numPosGenes = (int)round(xmer->enumPosSites);
	}

	qsort(xmer->stats, xmer->numXMers, sizeof(ConStat*), sortConStatDeg);

	if (cmd->exactFlag) {
		for (int i=0;i<xmer->numXMers;i++) {
			ConStatExact* cs = (ConStatExact*) xmer->stats[i];
			if (cmd->algorithm == ALG_RANK) {
				fprintf(fp, "%s\t%e\t%d\t%d\t%d\n", xmer->mers[i], cs->pvalue,
							xmer->numGenes,cs->sum(), cs->p);
			} else {
				float wp = 0.0,wn=0.0;
				for (int j=0;j<cs->p;j++) wp += xmer->geneWeights[cs->pgenes[j]];
				for (int j=0;j<cs->n;j++) wn += xmer->geneWeights[cs->ngenes[j]];

				fprintf(fp, "%s\t%e\t%.1f\t%.1f\t%.1f\t%.1f\n", xmer->mers[i],cs->pvalue,
							xmer->enumGenes,xmer->enumPosGenes,wp+wn,wp);
			}
		}
	} else {
		for (int i=0;i<xmer->numXMers;i++) {
			ConStatApprox* cs = (ConStatApprox*) xmer->stats[i];
			float pv = cs->pvalue;
			if (cmd->action == ACTION_DMERS ){
				pv = cs->dpvalue;
			}
			fprintf(fp, "%s\t%e\t%.1f\t%.1f\t%.1f\t%.1f\n", xmer->mers[cs->id],pv,
								xmer->enumGenes,xmer->enumPosGenes,cs->fn+cs->fp, cs->fp);
		}
	}
}



/* If ALG_SITEFREQ is used as the scoring function, numGenes and numPosGenes are used
   to store the total number of sites and total number of positive sites
   */
XMerData* indexGeneSeq(CommandLine *cmd) {

	Inttable* geneID = new Inttable();
	Floattable* geneValue = new Floattable();
	Floattable* geneSValue = new Floattable();
	int group = 0;
	char* buf = new char[BUFFER];
	char** cols = new char*[10000];
	char* name=NULL;
	char* seq=NULL;
	int numCols = 0;

	//float* conservation = new float[BUFFER];
	//int consLength = 0;
	int numGenes = 0;
	int numPosGenes = 0;
	//int minmer = cmd->minmer;
	//int maxmer = cmd->maxmer;
	int gapsize = cmd->gapsize;
	//int gapoffset = cmd->gapoffset;
	if (gapsize > 0) {
		//maxmer = minmer;
	}
	int maxPossible = cmd->merLengths[cmd->numMerLengths-1]
						+ cmd->gapLengths[cmd->numGapLengths-1];
	fprintf(stderr, "Max possible= %d\n", maxPossible);
	char* mer = new char[maxPossible+1];

	cmd->totalSeqLength = 0;
	if (cmd->seqLength != NULL) delete cmd->seqLength;
	cmd->seqLength = new Inttable();

	Hashtable* mers = new Hashtable(MOTIF_HASH_SIZE);

	//load group file
	float rvalue = 0;
	float svalue = 0;
	FILE* fp = fopen(cmd->statfile,"r");
	if (fp == NULL) {
		fprintf(stderr, "Stat file is null - can't index mers\n");
		exit(0);
	}
	float checkFloat = 0.0;

	while (fgets(buf,BUFFER,fp) != NULL) {

		split(buf, cols, numCols, '\t');
		if (numCols < 2) continue;
		name = cols[0];
		checkFloat = geneValue->search(name);

		if (cmd->algorithm == ALG_RANK) {
			if (checkFloat > EMPTY_FLOAT_CHECK) {
				continue;
			}
			sscanf(cols[1], "%f", &rvalue);
			geneValue->insert(rvalue,name);
		} else {
			sscanf(cols[1], "%d", &group);
			if (checkFloat > EMPTY_FLOAT_CHECK) {
				if (group != 0) {
					fprintf(stderr, "!!! Duplicate sequence ID (%s)\n", name);
				}
				continue;
			}

			if (numCols < 3) {
				rvalue = 1.0;
			} else {
				sscanf(cols[2], "%f", &rvalue);
			} 
			if (numCols < 4) {
				svalue = 1.0;
			} else {
				sscanf(cols[3], "%f", &svalue);
			}

			geneID->insert(group,name);
			geneValue->insert(rvalue,name);
			geneSValue->insert(svalue,name);
		}
	}       
	fclose(fp);

	fp = fopen(cmd->seqfile,"r");
	if (fp == NULL) {
		fprintf(stderr, "Could not open sequence file - can't index oligos\n");
		exit(0);
	}

	int numGenesStart = 0;
	if (cmd->algorithm == ALG_RANK) {
		numGenesStart = geneValue->total;
	} else {
		numGenesStart = geneID->total;
	}
	fprintf(stderr, "Num of Sequences in Stat file: %d\n", numGenesStart);

	//initialize variables to hold gene specific data
	char** geneNames= new char*[numGenesStart];
	//geneWeights doubles as the rank variable if ranks are used.
	float* geneWeights = new float[numGenesStart];
	int* geneStats = new int[numGenesStart];
	for (int i=0;i<numGenesStart;i++){ 
		geneNames[i] = NULL;
	}
	float backTotal = 0;
	if (cmd->seqType == DNA_SEQ_TYPE) {
		cmd->backFreq = new float[4];
		for (int i=0;i<4;i++) {
			cmd->backFreq[i] = 0;
		}
	}


	double numSites = 0;
	double numPosSites = 0;
	int geneINDEX = 0;
	Inttable* finishedGenes = new Inttable();
	
	while (fgets(buf,BUFFER,fp) != NULL) {

		split(buf,cols, numCols,'\t');
		if (numCols < 2) continue;
		name = cols[0];
		seq = cols[1];

		//check if sequence is in the group file
		int stat = 0;
		rvalue = 0;
		if (cmd->algorithm == ALG_RANK) {
			rvalue = geneValue->search(name);
			if (EMPTY_FLOAT_CHECK >  rvalue) continue;
		} else {
			stat = geneID->search(name);
			if (EMPTY_INT == stat) continue;
		}
		if (EMPTY_INT != finishedGenes->search(name)) {
			continue;
		} else {
			finishedGenes->insert(1,name);
		}

		//check if too many Ns are present
		int Ncount= 0;
		int nindex = 0;
		int seqlen = strlen(seq);
		if (seqlen < cmd->minSeqLen) continue;
		if (seqlen > cmd->maxSeqLen) continue;
		while (nindex < seqlen) {
			if (cmd->seqType == DNA_SEQ_TYPE && seq[nindex] == 'N') Ncount++;
			if (cmd->seqType == PROTEIN_SEQ_TYPE && seq[nindex] == 'X') Ncount++;
			nindex++;
		}
		if (((float)Ncount)/(float)seqlen > cmd->nratiocutoff) continue;

		//This gene has now passed the quality controls needed to be
		//considered as part of the dataset.  We now assign it an integer
		//ID to maintain a minimal memmory profile.
		gene_id_t geneid = (gene_id_t) geneINDEX++;
		geneNames[geneid] = new char[strlen(name)+1];
		strcpy(geneNames[geneid], name);
		geneWeights[geneid] = geneValue->search(name);
		float seqWeight = geneSValue->search(name);
		seqWeight *= geneWeights[geneid];
		geneStats[geneid] = stat;
		
		//fprintf(stderr, "%d\t%s\t", i, gene[i].name);
		if ( numGenes % 100 == 0) fprintf(stderr, "%d\t", numGenes);
		cmd->seqLength->insert(seqlen, name);
		cmd->totalSeqLength += seqlen;

		Inttable* seqhash = new Inttable();

		if (stat > 0) numPosGenes++;


		if (cmd->freqAdjust) {
			float w = seqWeight;
			if (cmd->seqType == DNA_SEQ_TYPE) {
				if (cmd->revoppFlag) {
					for (int j=0;j<seqlen;j++) {
						switch (seq[j]) {
							case 'N':
								break;
							case 'A':
								cmd->backFreq[0]+=w;
								cmd->backFreq[3]+=w;
								backTotal += 2*w;
								break;
							case 'C':
								cmd->backFreq[1]+=w;
								cmd->backFreq[2]+=w;
								backTotal += 2*w;
								break;
							case 'G':
								cmd->backFreq[1]+=w;
								cmd->backFreq[2]+=w;
								backTotal += 2*w;
								break;
							case 'T':
								cmd->backFreq[0]+=w;
								cmd->backFreq[3]+=w;
								backTotal += 2*w;
								break;
						}
					}
				} else {
					for (int j=0;j<seqlen;j++) {
						switch (seq[j]) {
							case 'N':
								break;
							case 'A':
								cmd->backFreq[0]+=w;
								backTotal += w;
								break;
							case 'C':
								cmd->backFreq[1]+=w;
								backTotal += w;
								break;
							case 'G':
								cmd->backFreq[2]+=w;
								backTotal += w;
								break;
							case 'T':
								cmd->backFreq[3]+=w;
								backTotal += w;
								break;
						}
					}
				}
			}
		}

		for (int j=0;j<seqlen;j++) {
			for (int pp=0;pp<cmd->numMerLengths;pp++) {
				int len = cmd->merLengths[pp];
				int gapoffset = len/2;
				if (j+len > seqlen) break; //check if near end of seq
				for (int qq=0;qq<cmd->numGapLengths;qq++) {
					int gapsize = cmd->gapLengths[qq];
					int bad = 0;
					if (j+len+gapsize > seqlen) break; // check if near end of seq
					
					for (int k=0;k<len;k++) {
						if (gapsize > 0) {
							if (k < gapoffset) {
								mer[k] = seq[j+k];
							} else {
								mer[k] = seq[j+k+gapsize];
							}
						} else {
							mer[k] = seq[j+k];
						}
						if (cmd->seqType == DNA_SEQ_TYPE && mer[k] == 'N') {
							bad = 1;
							break;
						} else if (cmd->seqType == PROTEIN_SEQ_TYPE && mer[k] == 'X') {
							bad = 1;
							break;
						}
					}
					mer[len] = '\0';
					if (bad) break;


					int present = seqhash->search(mer);
					if (present != EMPTY_INT) {
						seqhash->insert(present+1, mer);
						continue;
					}
					if (cmd->dualMotifs || cmd->allowFlip) {
						char* repmer = NULL;
						repmer = getRepDualMotif(cmd, mer);
						present = seqhash->search(repmer);
						if (present != EMPTY_INT) {
							seqhash->insert(present+1, repmer);
						} else {
							seqhash->insert(1,repmer);
						}
						delete []repmer;
						continue;
					} else if (cmd->revoppFlag) {
						char* rmer = revopp(mer);
						present = seqhash->search(rmer);
						if (present != EMPTY_INT) {
							seqhash->insert(present+1, rmer);
							delete []rmer;
							continue;
						}
						delete []rmer;
					}
					seqhash->insert(1, mer);
				}
			}
		}
		char** keys = seqhash->keys();
		int numkeys = seqhash->total;
		for (int j=0;j<numkeys;j++) {
			char* rv = NULL;
			int nsites = 1;
			//nsites = seqhash->search(keys[j]);
			if (cmd->algorithm == ALG_SITEBINOMIAL
					|| cmd->algorithm == ALG_SITEHYPERGEO
					|| cmd->algorithm == ALG_FISHER
					|| cmd->algorithm == ALG_FREQDIFF
					|| cmd->zoopsApprox) {	
				nsites = seqhash->search(keys[j]);
				if (cmd->zoopsApprox) {
					if (nsites > cmd->maxPerSeq) {
						nsites = cmd->maxPerSeq;
					}
				}
			}
			ConStat* cs = (ConStat*) mers->search(keys[j]);
			if (cs == NULL && cmd->revoppFlag) {
				rv = revopp(keys[j]);
				cs = (ConStat*) mers->search(rv);
			}
			if (cs == NULL) {
				if (cmd->exactFlag) {	
					cs = (ConStat*) new ConStatExact();
				} else {
					cs = (ConStat*) new ConStatApprox();
				}
				mers->insert(cs, keys[j]);
			}
			if (cmd->exactFlag) {
				((ConStatExact*)cs)->addGene(geneid,stat);
			} else {

				float wvalue = ((float)nsites)*seqWeight;
				if (stat) {
					((ConStatApprox*)cs)->fp += wvalue;
					numPosSites += wvalue;
				} else {
					((ConStatApprox*)cs)->fn += wvalue;
				}
				numSites+=wvalue;
			}
			if (rv != NULL) delete []rv;
			delete [](keys[j]);
		}
		delete []keys;
		delete seqhash;
		numGenes++;
	}
	delete geneValue;
	delete geneID;
	delete []buf;
	delete []cols;

	if (cmd->seqType == DNA_SEQ_TYPE) {
		if (cmd->freqAdjust) {
			if (backTotal == 0 ){
				fprintf(stderr, "Something is wrong with the input!!!!\n");
				exit(0);
			}
			for (int i=0;i<4;i++) {
				cmd->backFreq[i] /= backTotal;
				fprintf(stderr, "%f\n", cmd->backFreq[i]);
			}
		} else {
			for (int i=0;i<4;i++) {
				cmd->backFreq[i] = 0.25;
			}
		}	
	}
			
	/*
	Need to fix later - used for finding genes by rank
	if (cmd->algorithm == ALG_RANK) {
		char** gname = geneID->keys();	
		int numGname = geneID->total;
		rankvalues = geneValue;
		qsort(gname, numGname, sizeof(char*), rankDescending);
		fprintf(stderr, "\nNumGenes %d\n", numGname);
		for (int i=0;i<numGname;i++) {
			geneID->insert(i+1, gname[i]);
			delete [](gname[i]);
		}
		delete []gname;
		//ranks now stored in geneID
	}*/
	//float numNegGenes = numGenes - numPosGenes;
	if (cmd->zoopsApprox) {
		//adjust the number of negative genes to 
		//float avgLen = numPosSites / numPosGenes;
		//numNegGenes = (numSites-numPosSites) / avgLen;
	}

	//squeeze data into more effiecent data structures
	mers->setLowMem(1);
	XMerData* xmer = new XMerData();
	xmer->loadMerData(mers);
	xmer->loadGeneData(geneNames, geneWeights, geneStats, numGenes);
	xmer->enumSites = numSites;
	xmer->enumPosSites = numPosSites;
	xmer->enumNegSites = numSites - numPosSites;
	//xmer->enumNegGenes = numNegGenes;
	//xmer->enumGenes = numPosGenes+numNegGenes;

fprintf(stderr, "\nNumber of Target Genes: %.1f\nNumber of Backgd Genes: %.1f\n", 
										xmer->enumPosGenes, xmer->enumNegGenes);
fprintf(stderr, "Number of Target Sites: %.1f\nNumber of Backgd Sites: %.1f\n", 
										xmer->enumPosSites, xmer->enumNegSites);
	if (cmd->autoScale) {
		float psitesPerGene = numPosSites / xmer->enumPosGenes;
		float nsitesPerGene = (numSites-numPosSites) / 
								(xmer->enumNegGenes);
		float w = psitesPerGene / nsitesPerGene;
		fprintf(stderr, "Target Ratio: %e\n", psitesPerGene);
		fprintf(stderr, "Backgd Ratio: %e\n", nsitesPerGene);
		fprintf(stderr, "Autoscale adjustment factor: %f\n", w);
		if (!cmd->exactFlag) {
			for (int i=0;i<xmer->numXMers;i++) {
				((ConStatApprox*)xmer->stats[i])->fn *= w;
			}
		}
		xmer->enumNegSites *= w;
		xmer->enumSites = xmer->enumPosSites+xmer->enumNegSites;
	}
	
	delete mers;
	delete finishedGenes;



	//Since this uses non exact scoring, we scale every thing by 10 so
	// that we can perform the approximaiton to 1/10 of a gene
	xmer->initZOOPSCache();

	//prepare datastructures for memorizing statistics
	if (cmd->algorithm == ALG_HYPERGEO) {
		cmd->statmemory = new StatMemory(MAX_STAT_CACHE, (int)round(xmer->enumGenes),
					(int)round(xmer->enumPosGenes));
	}



	fprintf(stderr, "\nCalculating statistics for %d mers\n", xmer->numXMers);
	for (int i=0;i<xmer->numXMers;i++) {
		float pw = 0.0,nw=0.0;
		if (cmd->exactFlag) {
			ConStatExact* cs = (ConStatExact*) xmer->stats[i];
			for (int j=0;j<cs->p;j++) pw += geneWeights[cs->pgenes[j]];
			for (int j=0;j<cs->n;j++) nw += geneWeights[cs->ngenes[j]];

			if (cmd->algorithm == ALG_HYPERGEO) {
				cs->pvalue = cmd->statmemory->getStat(
							(int)round((double)cs->n+cs->p),(int)round((double)cs->p));
				/*cs->pvalue = loghypergeo((int)round(xmer->enumGenes),
							(int)round(xmer->enumPosGenes),
							(int)round(cs->n+cs->p),(int)round(cs->p));
				*/

			} else if (cmd->algorithm == ALG_BINOMIAL) {
				float r = ((float)cs->n) / ((float)(xmer->numNegGenes));
				cs->pvalue = logbinomial(xmer->numPosGenes, cs->p, r,
									xmer->numNegGenes);
			} else if (cmd->algorithm == ALG_FREQDIFF) {
				cs->pvalue = cs->n/xmer->enumNegSites - cs->p/xmer->enumPosSites;
				//cs->pvalue = cs->n/xmer->numNegGenes - cs->p/xmer->numPosGenes;
			} else if (cmd->algorithm == ALG_RANK) {
				//first need to replace gene names with ranks
				/*
				for (int j=0;j<cs->n;j++) {
					int r = geneID->search(cs->ngenes[j]);
					cs->addRank(r);
				}
				int ranksum = cs->sum();
	//			fprintf(stderr, "%d\t%d\t%d\n", ranksum, cs->p, cs->n);
				cs->pvalue = rankSumStat(ranksum,cs->p,numGenes);
				*/
			}
		} else {
			ConStatApprox* cs = (ConStatApprox*) xmer->stats[i];
			float ogp = 0, ogn = 0;
			if (cmd->zoopsApprox) {
				//This allows us to approximate additional motifs
				ogp = cs->fp;
				ogn = cs->fn;
				int ep = (int)round(cs->fp*WEIGHT_RES);
				int en = (int)round(cs->fn*WEIGHT_RES);
				cs->fp = xmer->papprox[ep] / ((float)WEIGHT_RES);
				cs->fn = xmer->napprox[en] / ((float)WEIGHT_RES);
			}
			if (cmd->algorithm == ALG_HYPERGEO) {
				cs->pvalue = cmd->statmemory->getStat(
							(int)round(cs->fn+cs->fp),(int)round(cs->fp));
				/*cs->pvalue = loghypergeo((int)round(xmer->enumGenes),(int)round(xmer->enumPosGenes),
								(int)round(cs->fn+cs->fp),(int)round(cs->fp));
				*/
			} else if (cmd->algorithm == ALG_APPROXBINOMIAL) {
				float r = ((float)cs->fn)/(xmer->enumNegSites);
				//fprintf(stderr, "%d\t%d\t%f\n", numPosSites, cs->p, r);
				cs->pvalue = logbinomial((int)round(xmer->enumPosSites), (int)round(cs->fp), r,
										(int)round(xmer->enumNegSites));
			} else if (cmd->algorithm == ALG_SITEBINOMIAL) {
				float r = ((float)cs->fn)/(xmer->enumNegSites);
				//fprintf(stderr, "%d\t%d\t%f\n", numPosSites, cs->p, r);
				cs->pvalue = logbinomial((int)round(xmer->enumPosSites), (int)round(cs->fp), r,
										(int)round(xmer->enumNegSites));
			} else if (cmd->algorithm == ALG_SITEHYPERGEO) {
				cs->pvalue = loghypergeo((int)round(xmer->enumSites),(int)round(xmer->enumPosSites),
									(int)round(cs->fn+cs->fp), (int)round(cs->fp));
			} else if (cmd->algorithm == ALG_BINOMIAL) {
				float r = cs->fn / xmer->enumNegGenes;
				cs->pvalue = logbinomial((int)round(xmer->enumPosGenes), 
											(int)round(cs->fp), r, 
											(int)round(xmer->enumNegGenes));
			} else if (cmd->algorithm == ALG_FREQDIFF) {
				//cs->pvalue = cs->fn/xmer->enumNegGenes - cs->fp/xmer->enumPosGenes;
				cs->pvalue = cs->fn/xmer->enumNegSites - cs->fp/xmer->enumPosSites;
			}
			if (cmd->zoopsApprox) {
				cs->fp = ogp;
				cs->fn = ogn;
			}
		}
		xmer->stats[i]->dpvalue = xmer->stats[i]->pvalue;
	}
	fprintf(stderr, "done\n");

	return xmer;
}

char* getRepDualMotif(CommandLine* cmd, char* m) {
	int len = strlen(m);
	int slen = len/2;
	int flip1 = -1;
	int flip2 = -1;
	
	if (cmd->allowFlip && cmd->dualMotifs == 0) {
		//need to use revopp of the whole sequence
		char* rv = revopp(m);
		char* s = new char[len+1];
		s[len]='\0';
		char* rs = new char[len+1];
		rs[len]='\0';

		for (int i=0;i<slen;i++) {
			if (m[i] < m[i+slen]) {
				flip1 = 0;
				break;
			} else if (m[i] > m[i+slen]) {
				flip1 = 1;
				break;
			}
		}
		if (flip1 > 0) {
			for (int i=0;i<slen;i++) {
				s[i] = m[i+slen];
				s[i+slen] = m[i];
			}
		} else {
			for (int i=0;i<slen;i++) {
				s[i] = m[i];
				s[i+slen] = m[i+slen];
			}
		}

		for (int i=0;i<slen;i++) {
			if (rv[i] < rv[i+slen]) {
				flip2 = 0;
				break;
			} else if (rv[i] > rv[i+slen]) {
				flip2 = 1;
				break;
			}
		}
		if (flip2 > 0) {
			for (int i=0;i<slen;i++) {
				rs[i] = rv[i+slen];
				rs[i+slen] = rv[i];
			}
		} else {
			for (int i=0;i<slen;i++) {
				rs[i] = rv[i];
				rs[i+slen] = rv[i+slen];
			}
		}
		delete []rv;
		flip1 = 0;
		for (int i=0;i<len;i++) {
			if (s[i] < rs[i]) {
				flip1 = 0;
				break;
			} else if (s[i] > rs[i]) {
				flip1 = 1;
				break;
			}
		}
		if (flip1 ==0) {
			delete []rs;
			return s;
		} else {
			delete []s;
			return rs;
		}
	}
	if (cmd->dualMotifs) {

		char* m1 = new char[slen+1];
		char* m2 = new char[slen+1];
		char* s = new char[len+1];
		s[len] = '\0';
		m1[slen]='\0';
		m2[slen]='\0';

		for (int i=0;i<slen;i++) {
			m1[i] = m[i];
			m2[i] = m[i+slen];
		}
		char* rv1 = revopp(m1);
		char* rv2 = revopp(m2);
		for (int i=0;i<slen;i++) {
			if (m1[i] < rv1[i]) {
				flip1 = 0;
				break;
			} else if (m1[i] > rv1[i]) {
				flip1 = 1;
				break;
			}
		}
		if (flip1 > 0) {
			delete []m1;
			m1 = rv1;
		} else {
			delete []rv1;
		}
		
		for (int i=0;i<slen;i++) {
			if (m2[i] < rv2[i]) {
				flip2 = 0;
				break;
			} else if (m2[i] > rv2[i]) {
				flip2 = 1;
				break;
			}
		}
		if (flip2 > 0) {
			delete []m2;
			m2 = rv2;
		} else {
			delete []rv2;
		}
		if (cmd->allowFlip == 0) {
			for (int i=0;i<slen;i++) {
				s[i] = m1[i];
				s[i+slen] = m2[i];
			}
		} else {
			flip1 = 0;
			for (int i=0;i<slen;i++) {
				if (m1[i] < m2[i]) {
					flip1 = 0;
					break;
				} else if (m1[i] > m2[i]) {
					flip1 = 1;
					break;
				}
			}
			if (flip1 > 0) {
				for (int i=0;i<slen;i++) {
						s[i] = m2[i];
					s[i+slen] = m1[i];
				} 
			} else {
				for (int i=0;i<slen;i++) {
					s[i] = m1[i];
					s[i+slen] = m2[i];
				} 
			}
		}
		delete []m1;
		delete []m2;
		return s;
	}
	return NULL;

}

XMerData::XMerData() {
	init();
}
void XMerData::loadMerData(Hashtable* hash) {
	init();
	numXMers = hash->total;
	mers = hash->keys();
	stats = new ConStat*[numXMers];
	for (int i=0;i<numXMers;i++) {
		stats[i] = (ConStat*) hash->remove(mers[i]);
		stats[i]->id = i;
	}
}

void XMerData::loadGeneData(char** gnames,float* gweights, 
									int* gstats,int ngenes) {
	geneNames = gnames;
	geneWeights = gweights;
	geneStats = gstats;
	numGenes = ngenes;
	numPosGenes = 0;
	numNegGenes = 0;
	for (int i=0;i<numGenes;i++) {
		if (geneStats[i] == 1) {
			numPosGenes++;
			enumPosGenes+=geneWeights[i];
		} else {
			numNegGenes++;
			enumNegGenes+=geneWeights[i];
		}
	}
	enumGenes = enumPosGenes+enumNegGenes;
}

void XMerData::init() {
	numGenes = 0;
	numPosGenes = 0;
	numNegGenes = 0;
	enumGenes = 0.0;
	enumPosGenes = 0.0;
	enumNegGenes = 0.0;
	numXMers = 0;
	stats = NULL;
	mers = NULL;
	seqTree = NULL;
	papprox = NULL;
	napprox = NULL;
}
void XMerData::initZOOPSCache() {

	float * pmemory = NULL;
	float * nmemory = NULL;
	pmemory = new float[ZOOPS_APPROX_CACHE];
	nmemory = new float[ZOOPS_APPROX_CACHE];
	pmemory[1] = 1;
	nmemory[1] = 1;
	for (int i=2;i<ZOOPS_APPROX_CACHE;i++) {
		//pmemory[i] = i;
		//nmemory[i] = i;
		pmemory[i] = pmemory[i-1] + (enumPosGenes*WEIGHT_RES-pmemory[i-1])
								/(WEIGHT_RES*enumPosGenes);
		nmemory[i] = nmemory[i-1] + (enumNegGenes*WEIGHT_RES-nmemory[i-1])
								/(WEIGHT_RES*enumNegGenes);
		
	}
	fprintf(stderr, "\nzoopsApprox implemented\n");
	papprox = pmemory;
	napprox = nmemory;
}

XMerData::~XMerData() {
	for (int i=0;i<numXMers;i++) {
		if (stats[i] != NULL) delete [](stats[i]);
		if (mers[i] != NULL) delete [](mers[i]);
	}
	delete []stats;
	delete []mers;
	if (papprox != NULL) delete []papprox;
	if (napprox != NULL) delete []napprox;
}


ConStat::ConStat() {
	pvalue = 1.0;
	dpvalue = 1.0;
	id = 0;
	dir = 0;
}
ConStatApprox::ConStatApprox() {
	fp = 0.0;
	fn = 0.0;
}
ConStatExact::ConStatExact() {
	p = 0;
	n = 0;
	pgenes = NULL;
	ngenes = NULL;
}
ConStatExact::~ConStatExact() {
	if (pgenes != NULL) delete []pgenes;
	if (ngenes != NULL) delete []ngenes;
}
int ConStatExact::sum() {
	int sum = 0;
	if (pgenes == NULL) return 0;
	for (int i=0;i<p;i++) {
		sum+=pgenes[i];
	}
	return sum;
}

void ConStatExact::addRank(int r) {
	gene_id_t * rr = new gene_id_t[p+1];
	for (int i=0;i<p;i++) rr[i] = pgenes[i];
	rr[p]=r;
	if (pgenes != NULL) delete []pgenes;
	pgenes = rr;
	p++;
}

void ConStatExact::addGene(gene_id_t g, int stat) {
	if (stat) {
		gene_id_t* pg = new gene_id_t[p+1];
		for (int i=0;i<p;i++) pg[i] = pgenes[i];
		pg[p] = g;
		if (pgenes != NULL) delete []pgenes;
		pgenes = pg;
		p++;
	} else {
		gene_id_t* ng = new gene_id_t[n+1];
		for (int i=0;i<n;i++) ng[i] = ngenes[i];
		ng[n] = g;
		if (ngenes != NULL) delete []ngenes;
		ngenes = ng;
		n++;
	}
}


float alignPSSMs(PSSM* a, PSSM* b, int &bestOffset, int &bestDir) {
	a->expXform();
	b->expXform();
	int plen1 = a->length;
	int plen2 = b->length;
	int minOverlap = (int)round((double)plen1*0.66);
	int minOverlap2 = (int)round((double)plen1*0.66);
	if (minOverlap2 > minOverlap) minOverlap = minOverlap2;
	
	float bestScore = FLT_MIN;
	bestOffset = 0;
	bestDir = 0;
	for (int i=-plen2+1;i<plen1;i++) {
		int overlap = 0;
		int index1 = 0,index2=0;
		float rscore = 0,score = 0;
		int dir =0;
		
		if (i < 0) {
			overlap = i+plen2;
			index1 = 0;
			index2 = -1*i;
			if (overlap > plen1) overlap = plen1;

		} else {
			overlap = plen1-i;
			index1 = i;
			index2 = 0;
			if (overlap > plen2) overlap = plen2;
		}
		for (int j=0;j<overlap;j++) {
			for (int k=0;k<4;k++) {
				score += a->matrix[index1][k]*b->matrix[index2][k];
				rscore += a->matrix[index1][k]*b->matrix[plen2-1-index2][3-k];
			}
			//subtract the expected value from aligning N vs. N
		//	score -= 0.25;
		//	rscore -= 0.25;
			index1++;
			index2++;
		}
		
		if (score < rscore) {
			dir = 1;
			score = rscore;
		}
		//normalize score
		score -= (overlap*4)*0.25*0.25;
		score /= sqrt((double)overlap*4);

		if (overlap >= minOverlap) {
			if (score > bestScore) {
				bestScore = score;
				bestDir = dir;
				bestOffset = i;
			}
		}
	}
	return bestScore;
}



float calculateProbabilityPSSM(CommandLine* cmd, PSSM* pssm,char* mer,unsigned char &dir,
			int allowRevOpp) {
	float score = 0,rscore = 0;
	int len = pssm->length;
	static char *m1=NULL;
	static char *m2=NULL;
	static PSSM pssm1;
	static PSSM pssm2;
	if (cmd->dualMotifs == 0 && cmd->allowFlip == 1) {
		if (m1 == NULL) m1 = new char[100];
		if (m2 == NULL) m2 = new char[100];
		//construct flipped
		int slen = len/2;
		for (int i=0;i<slen;i++) {
			m1[i] = mer[i+slen];
		}
		for (int i=slen;i<len;i++) {
			m1[i] = mer[i-slen];
		}
		unsigned char r = 0;
		unsigned char r1 = 0;
		cmd->allowFlip = 0;
		score = calculateProbabilityPSSM(cmd,pssm,mer,r,1);
		rscore = calculateProbabilityPSSM(cmd,pssm,m1,r1,1);
		cmd->allowFlip = 1;
		dir = 0;
		if (score > rscore) {
			if (r) {
				dir = dir | 4;
				dir = dir | 8;
			}
		} else {
			score = rscore;
			dir = dir | 2; 
			if (r) {
				dir = dir | 4;
				dir = dir | 8;
			}
		}
		return score;
	}
	if (cmd->dualMotifs == 1 && cmd->allowFlip == 0) {
		if (m1 == NULL) m1 = new char[100];
		if (m2 == NULL) m2 = new char[100];
		int slen = len/2;
		char* rv = revopp(mer);
		cmd->dualMotifs = 0;
		unsigned char d1=0,d2=0,d3=0;
		float s1 = calculateProbabilityPSSM(cmd, pssm, mer,d1,1);
		for (int i=0;i<slen;i++) {
			m1[i] = mer[i];
			m1[i+slen] = rv[i];
		}
		float s2 = calculateProbabilityPSSM(cmd, pssm, m1,d2,0);
		for (int i=0;i<slen;i++) {
			m1[i] = rv[slen+i];
			m1[i+slen] = mer[slen+i];
		}
		float s3 = calculateProbabilityPSSM(cmd, pssm, m1,d3,0);
		if (s2 > s3) {
			dir = 0 | 0 | 0 | 8;
			score = s2;
		} else {
			dir = 0 | 0 | 4 | 0;
			score = s3;
		}
		if (score < s1) {
			score = s1;
			if (d1) {
				dir = 0 | 0 | 4 | 8;
			} else {
				dir = 0;
			}
		}
		cmd->dualMotifs = 1;
		delete []rv;
		return score;
		
	}
	if (cmd->dualMotifs == 1 && cmd->allowFlip == 1) {
		if (m1 == NULL) m1 = new char[100];
		if (m2 == NULL) m2 = new char[100];
		int slen = len/2;
		for (int i=0;i<slen;i++) {
			m1[i] = mer[i];
			m2[i] = mer[i+slen];
		}
		pssm1.length = slen;
		pssm1.matrix = pssm->matrix;
		pssm2.length = slen;
		pssm2.matrix = &(pssm->matrix[slen]);

		unsigned char r11 = 0,r12 = 0,r21 = 0,r22 = 0;
		cmd->dualMotifs = 0;
		cmd->allowFlip = 0;
		float s12=calculateProbabilityPSSM(cmd,&pssm2,m1,r12,1);
		float s11=calculateProbabilityPSSM(cmd,&pssm1,m1,r11,1);
		float s21=calculateProbabilityPSSM(cmd,&pssm1,m2,r21,1);
		float s22=calculateProbabilityPSSM(cmd,&pssm2,m2,r22,1);
		cmd->dualMotifs = 1;
		cmd->allowFlip = 1;
		if (s11 + s22 < s12 + s21) {
			score = s11+s22;
			dir = 0;
			if (r11) {
				dir = dir | 4;
			}
			if (r22) {
				dir = dir | 8;
			}
		} else {
			score = s12+s21;
			dir = 2;
			if (r12) {
				dir = dir | 4;
			} 
			if (r21) {
				dir = dir | 8;
			}
		}
		pssm2.matrix = NULL;
		pssm1.matrix = NULL;
		return score;
	}
	if (cmd->seqType == PROTEIN_SEQ_TYPE) {
		for (int i=0;i<len;i++) {
			score += pssm->matrix[i][(int)mer[i]];
			rscore += pssm->matrix[i][(int)(19-mer[len-i-1])];
		}
	} else if (cmd->seqType == DNA_SEQ_TYPE) {
		for (int i=0;i<len;i++) {
			score += pssm->matrix[i][(int)mer[i]];
			rscore += pssm->matrix[i][(int)(3-mer[len-i-1])];
		}
	}
	dir = 0;
	if (score < rscore && cmd->revoppFlag) {
		dir = 1;
		score = rscore;
	}
	return score;
}

SeqTree::SeqTree(char nsize) {
	branch = NULL;
	id = EMPTY_GENE;
	size = nsize;
	IUPAC = 1;
	IUPACIndex = 0;
}
SeqTree::~SeqTree() {
	if (branch != NULL) {
		for (int i=0;i<size;i++) {
			if (branch != NULL) {
				delete branch[i];
			}
			branch[i] = NULL;
		}
		delete []branch;
	}
}
int SeqTree::getNextIUPAC(int* value, int clevel,int numIUPAC, int IUPACtype) {
	//int index = value[clevel];
	if (id != EMPTY_GENE) {
		value[clevel]=0;
		return 0;
	}
	if (branch == NULL) {
		return -1;
	}
	if (size == 4) {
		while (IUPAC < 16) {
			if ((IUPACtype == IUPAC_2N || IUPACtype== IUPAC_N) && 
						(IUPAC == 7 || IUPAC == 11 || IUPAC == 13 || IUPAC==14))  {
				IUPAC++;
				continue;
			}
			if ((IUPACtype == IUPAC_N) && 
						(IUPAC == 3 || IUPAC == 5 || IUPAC == 6 || IUPAC== 9 || 
							IUPAC == 10 || IUPAC == 12))  {
				IUPAC++;
				continue;
			}
			int curNum = numIUPAC;
			if (IUPAC != 1 && IUPAC != 2 && IUPAC != 4 && IUPAC != 8) {
				curNum--;
			}
			if (curNum < 0) {
				IUPACIndex = 0;
				IUPAC++;
				continue;
			}

			value[clevel] = IUPAC;
			while (IUPACIndex < 4) {
				//fprintf(stderr, "=== %d\t%d\t%d\t%d\n", IUPAC, IUPACIndex, 1 << IUPACIndex, IUPAC & (1 << IUPACIndex));
				if ((IUPAC & (1 << IUPACIndex)) == 0) {
					IUPACIndex++;
					continue;
				}
				if (branch[IUPACIndex] == NULL) {
					IUPACIndex++;
					continue;
				}
				//fprintf(stderr, "%d\t%d\t%d\n", IUPAC, IUPACIndex, 1 << IUPACIndex);
				int check = branch[IUPACIndex]->getNextIUPAC(value,clevel+1,curNum, IUPACtype);
				if (check ==1) return 1;
				if (check ==0) {
					IUPACIndex=4;
					//IUPACIndex++;
					return 1;
				}
				if (check == -1) {
					branch[IUPACIndex]->IUPAC = 1;
					IUPACIndex++;
				}

			}
			IUPACIndex = 0;
			IUPAC++;
		}
	} else if (size == 20) {
		while (IUPAC < 0x00100000) {
	
			int curNum = numIUPAC;
			int test = 0;
			for (int j=0;j<19;j++) {
				if ((1 << j) & IUPAC) {
					test++;
					if (test > 2) {
						curNum -= 1;
						break;
					}
				}
			}
	
			value[clevel] = IUPAC;
			while (curNum >= 0 && IUPACIndex < 20) {
				if ((IUPAC & (1 << IUPACIndex)) == 0) {
					IUPACIndex++;
					continue;
				}
				if (branch[IUPACIndex] == NULL) {
					IUPACIndex++;
					continue;
				}
				//fprintf(stderr, "%d\t%d\t%d\n", IUPAC, IUPACIndex, 1 << IUPACIndex);
				int check = branch[IUPACIndex]->getNextIUPAC(value,clevel+1,curNum, IUPACtype);
				if (check ==1) return 1;
				if (check ==0) {
					//IUPACIndex++;
					IUPACIndex = 20;
					return 1;
				}
				if (check == -1) {
					branch[IUPACIndex]->IUPAC = 1;
					IUPACIndex++;
				}
			}
			IUPACIndex = 0;
	
			IUPAC = IUPAC << 1;
			if (IUPAC > 0x000fffff) {
				if (IUPAC == 0x00100000 && (IUPACtype >= IUPAC_N)) {
					IUPAC = 0x000fffff;
				} else {
					break;
				}
			}
		}	
	}
	return -1;
}


//Inserts sequence into sequnce tree - the primary input (mer), is a sequence
// converted using the seqIndexOf function, which assigns A to 0, C to 1, G to 3, etc.
void SeqTree::insert(char* mer, int clevel, int nid) {
	int index = mer[clevel];
	if (index == -1) {
		id = nid;
		//fprintf(stderr,"End of the string!! id=%d level=%d\n",id,clevel);
		return;
	}
//fprintf(stderr, "clevel=%d index=%d\n", clevel, index);
//	int INDEX = 0;
//	while ( ((1 << INDEX) & index) == 0) {
//		INDEX++;
//		if (INDEX > 30) {
//			fprintf(stderr, "Problem inserting sequence - INDEX=%d\n",INDEX);
//			exit(0);
//		}
//	}

	if (branch == NULL) {
		branch = new SeqTree*[size];
		for (int i=0;i<size;i++) {
			branch[i] = NULL;
		}
	}
	if (branch[index] == NULL) {
		branch[index] = new SeqTree(size);
	}
	branch[index]->insert(mer,clevel+1,nid);
}


void SeqTree::findMisMatchMers(int* mer, int clevel, int mis, 
						int* ids, int &numIDs) {
	int index = mer[clevel];
	if (index == 0) {
		ids[numIDs++] = id;
		return;
	}
	if (branch == NULL) {
		return;
	}
	int ref = 1;
	for (int i=0;i<size;i++) {
		if (branch[i] == NULL) {
			continue;
		}
		int comp = ref << i;
		int mismatch = 1;
		if (comp & index) {
			mismatch = 0;
		}
		if (mis == 0 && mismatch == 1) continue;

//fprintf(stderr, "branch = %d\tmis = %d\tindex = %d\n", i, mis, index);
		branch[i]->findMisMatchMers(mer,clevel+1,mis-mismatch,ids,numIDs);
	}
}

IndexTree::IndexTree() {
	values = NULL;
	tree = NULL;
	size = 0;
	recallIndex = 0;
}
IndexTree::~IndexTree() {
	if (tree != NULL) {
		for (int i=0;i<size;i++) {
			if (tree[i] != NULL){ 
				delete tree[i];
				tree[i] = NULL;
			}
		}
		delete []tree;
	}
	tree = NULL;
	if (values != NULL) {
		delete []values;
	}
	values = NULL;
	size =0;
}
int IndexTree::insert(int* mer, int clevel) {
	int v = mer[clevel];
	if (v == 0) {
		if (size == 0) {
			size = 1;
			return 1;
		} else {
			return 0;
		}
	}
	for (int i=0;i<size;i++) {
		if (values[i] == v) {
			return tree[i]->insert(mer, clevel+1);
		}
	}
	int* newValues = new int[size+1];
	IndexTree** newTree = new IndexTree*[size+1];
	//int found = 0;
	for (int i=0;i<size;i++) {
		newValues[i] = values[i];
		newTree[i] = tree[i];
	}
	newValues[size] = v;
	newTree[size] = new IndexTree();
	if (values != NULL) {
		delete []values;
	}
	if (tree != NULL) {
		delete []tree;
	}
	tree = newTree;
	values = newValues;
	size++;
	return tree[size-1]->insert(mer, clevel+1);
}
int IndexTree::getNextSeq(int* value, int clevel) {
	//int v = value[clevel];
/*
	if (index == 0) {
		value[clevel]=0;
		return 0;
	}
*/
	if (tree == NULL) {
		return -1;
	}
	while (recallIndex < size) {
		value[clevel] = values[recallIndex];
		int rv = tree[recallIndex]->getNextSeq(value, clevel+1);
		if (rv == -1) {
		}
	}
	return 0;
		
}


char** findBestSeeds(CommandLine* cmd, XMerData* xmer, int *merIDs, 
										int numMerIDs, int &numSeeds) {

	static int* values = NULL;
	if (values == NULL) values = new int[1000];
	int motifLength = cmd->merLengths[0];

	//Initialize SeqTree to store xmers in an easy to search data structure
	SeqTree* root = NULL;
	SeqTree* target = NULL;
	if (cmd->seqType == DNA_SEQ_TYPE) {
		root = new SeqTree(4);
		if (cmd->targetTreeFlag) {
			target = new SeqTree(4);
		}
	} else {
		root = new SeqTree(20);
		if (cmd->targetTreeFlag) {
			target = new SeqTree(20);
		}
	}


	//Insert xmers into the SeqTrees
	for (int i=0;i<xmer->numXMers;i++) {
		char* cvalues = new char[1000];
		cvalues = seqIndexOf(cmd->seqType,xmer->mers[i],cvalues);
		root->insert(cvalues, 0, i);
		int tflag = 0;
		if (cmd->targetTreeFlag) {
			if (cmd->exactFlag) {
				if (((ConStatExact*)(xmer->stats[i]))->p > 0) {
					tflag = 1;
				}
			} else {
				if (((ConStatApprox*)(xmer->stats[i]))->fp > 0) {
					tflag = 1;
				}
			}
		}
		if (tflag) {
			target->insert(cvalues,0,i);
		}
//fprintf(stderr, "GOT it\n");
		if (cmd->revoppFlag) {
			char* rv = revopp(xmer->mers[i]);
			cvalues = seqIndexOf(cmd->seqType,rv,cvalues);
			root->insert(cvalues, 0, i);
			if (tflag) {
				target->insert(cvalues,0,i);
			}
			delete []rv;
		}
		delete []cvalues;
	}


	//variables needed for keeping track of ridiculous amounts of data
	static int* misMatchMers = NULL;
	if (misMatchMers == NULL) misMatchMers = new int[MAX_MISMATCHMERS];
	int numMisMatchMers = 0;
	static char* geneHash = NULL;
	if (geneHash == NULL) geneHash = new char[xmer->numGenes+1];
	static gene_id_t* geneList = NULL;
	if (geneList == NULL) geneList = new gene_id_t[xmer->numGenes+1];
	int numGeneList = 0;
	float ngenes =0.0;
	float pgenes =0.0;
	for (int i=0;i<xmer->numGenes;i++) {
		geneHash[i] = 0;
	}
	static char* usedMers = NULL;
	if (usedMers == NULL) usedMers = new char[xmer->numXMers];
	for (int i=0;i<xmer->numXMers;i++) {
		usedMers[i]= 0;
		xmer->stats[i]->dir = 0;
	}
	int numMers = numMerIDs;
	if (numMers == 0) {
		numMers = xmer->numXMers;
	}

	int numMotifs = 0;
	char** motifsDNA = NULL;
	ProteinStat* motifsProtein = NULL;
	if (cmd->seqType == PROTEIN_SEQ_TYPE) {
		motifsProtein = new ProteinStat[MAX_GLOBAL_AA];
	}

	Floattable* possible = new Floattable(MOTIF_HASH_SIZE);
	fprintf(stderr, "Generating possible IUPAC motifs...");
	if (cmd->targetTreeFlag == 0) {
		target = root;
	}
	while (-1 != target->getNextIUPAC(values, 0, cmd->numIUPAC, cmd->IUPACtype)) {
		if (cmd->seqType == DNA_SEQ_TYPE) {
			char* str = getStringFromIUPACIndex(values);
			if (possible->search(str) < EMPTY_FLOAT_CHECK) {
				if (cmd->revoppFlag) {
					char *rm = revopp(str);
					if (possible->search(rm) < EMPTY_FLOAT_CHECK) {
						possible->insert(1.0, str);
					}
					delete []rm;
				} else {
					possible->insert(1.0, str);
				}
			}
			delete []str;
		} else {
			if (numMotifs == MAX_GLOBAL_AA) {
				fprintf(stderr, 
					"Reached maximum space statically alocated for peptide global search\n");
				break;
			}
			motifsProtein[numMotifs].v = new int[motifLength+1];
			for (int i=0;i<motifLength+1;i++){ 
				motifsProtein[numMotifs].v[i] = values[i];
			}
			numMotifs++;
		}
	}
	if (cmd->targetTreeFlag) {
		delete target;
	}

	if (cmd->seqType == DNA_SEQ_TYPE) {
		motifsDNA = possible->keys();
		numMotifs = possible->total;
	}

	fprintf(stderr, " %d total\nChecking Motif p-values..\n", numMotifs);

	for (int z=0;z<=cmd->mismatches;z++) {
		int totalChecked = 0;
		fprintf(stderr, "Checking %d mismatches...\n", z);
		for (int i=0;i<numMotifs;i++) {
			if (i % 10000 == 0) {
				fprintf(stderr, "\t%d",i); 
			}
			float currentPvalue = 0;
			if (cmd->seqType == DNA_SEQ_TYPE) {
				currentPvalue = possible->search(motifsDNA[i]);
			} else {
				currentPvalue = motifsProtein[i].pvalue;
			}
			if (cmd->speed == SPEED_FAST) {
				if (currentPvalue > 1.1) {
					continue;
				}
			}
			float pvalue = 1.0;

			numMisMatchMers = 0;
			ngenes = 0.0;
			pgenes = 0.0;
			numGeneList = 0;

			if (cmd->seqType == DNA_SEQ_TYPE) {
				values = seqIndexOfIUPAC(cmd->seqType,motifsDNA[i],values);
				root->findMisMatchMers(values, 0, z, misMatchMers, numMisMatchMers);
			} else {
				root->findMisMatchMers(motifsProtein[i].v, 0, z, misMatchMers, numMisMatchMers);
			}
			totalChecked++;

//fprintf(stderr, "%s\n", motifs[i]);
			for (int j=0;j<numMisMatchMers;j++){ 
				int merID = misMatchMers[j];
//				fprintf(stderr, "\t%s\n", xmer->mers[merID]);
				if (usedMers[merID] != 0) continue;
				usedMers[merID] = 1;
				if (cmd->exactFlag) {
					ConStatExact* cs = (ConStatExact*)xmer->stats[merID];
					for (int k=0;k<cs->p;k++) {
						gene_id_t gid = cs->pgenes[k];
						if (geneHash[gid] == 1) continue;
						geneHash[gid] = 1;
						geneList[numGeneList++] = gid;
						pgenes += xmer->geneWeights[gid];
					}
					for (int k=0;k<cs->n;k++) {
						gene_id_t gid = cs->ngenes[k];
						if (geneHash[gid] == 1) continue;
						geneHash[gid] = 1;
						geneList[numGeneList++] = gid;
						ngenes += xmer->geneWeights[gid];
					}
				} else {
					ConStatApprox* cs = (ConStatApprox*)xmer->stats[merID];
					if (cmd->zoopsApprox) {
						pgenes += cs->fp;
						ngenes += cs->fn;
					} else {
						pgenes += (xmer->enumPosGenes-pgenes)*cs->fp/xmer->enumPosGenes;
						ngenes += (xmer->enumNegGenes-ngenes)*cs->fn/xmer->enumNegGenes;
					}
				}
			}

			//float ppg = pgenes;
			//float npg = ngenes;
			if (cmd->zoopsApprox) {
				pgenes = xmer->papprox[(int)round(pgenes*WEIGHT_RES)]/WEIGHT_RES;
				ngenes = xmer->napprox[(int)round(ngenes*WEIGHT_RES)]/WEIGHT_RES;
			}

			if (cmd->algorithm == ALG_HYPERGEO) {
				pvalue = cmd->statmemory->getStat(
							(int)round(pgenes+ngenes),(int)round(pgenes));
				/*pvalue = loghypergeo((int)round(xmer->enumGenes),
							(int)round(xmer->enumPosGenes),
							(int)round(pgenes+ngenes),(int)round(pgenes));
				*/
			} else if (cmd->algorithm == ALG_SITEHYPERGEO) {
				pvalue = loghypergeo((int)round(xmer->enumSites),
							(int)round(xmer->enumPosSites), 
							(int)round(ngenes+pgenes),(int)round(pgenes));
			} else if (cmd->algorithm == ALG_SITEBINOMIAL) {
				float r = (ngenes) / (xmer->enumNegSites);
				pvalue = logbinomial((int)round(xmer->enumPosSites), 
							(int)round(pgenes), r, (int)round(xmer->enumNegSites));
			} else if (cmd->algorithm == ALG_BINOMIAL) {
				float r = ngenes / (xmer->enumNegGenes);
				pvalue = logbinomial((int)round(xmer->enumPosGenes),
							(int)round(pgenes),r,(int)round(xmer->enumNegGenes));
			} else if (cmd->algorithm == ALG_RANK) {
				//pvalue = rankSumStat(ranksum, P, numGenes);
			} else if (cmd->algorithm == ALG_FREQDIFF) {
				//pvalue = ngenes/xmer->enumNegGenes - pgenes/xmer->enumPosGenes;
				pvalue = ngenes/xmer->enumNegSites - pgenes/xmer->enumPosSites;
			}
			if (merIDs != NULL) {
				fprintf(stderr, "%d\t%f\t%f\t%f\n",z, pvalue,
						pgenes/xmer->enumPosGenes,ngenes/xmer->enumNegGenes);
			}
//fprintf(stderr, ":: pvalue = %e\n", pvalue);

			//char* str = getStringFromIUPACIndex(values);
			//printf("%s\t%e\t%d\n",str,pvalue,z);
			//delete []str;

			// record if it is a good value...

			//reset lists used to keep track of data
			if (cmd->exactFlag) {
				for (int j=0;j<numGeneList;j++) {
					geneHash[geneList[j]] = 0;
				}
				numGeneList = 0;
			}
			for (int j=0;j<numMisMatchMers;j++) {
				usedMers[misMatchMers[j]] = 0;
			}
			numMisMatchMers = 0;


			if (pvalue < currentPvalue) {
				if ((cmd->speed == SPEED_FAST) && (pvalue > DEG_PVALUE_THRESH)) {
					pvalue = 2.0;
				}
				if (cmd->seqType == DNA_SEQ_TYPE) {
					possible->insert(pvalue, motifsDNA[i]);
				} else {
					motifsProtein[i].pvalue = pvalue;
				}
			}
		}
		fprintf(stderr,"\n\t%d of %d\n",totalChecked,numMotifs);
	}

	//delete common variables
	delete []values;
	xmer->seqTree = root;
	if (usedMers != NULL) {
		delete []usedMers;
	}

	if (cmd->seqType == DNA_SEQ_TYPE) {
		sorttable2 = possible;
		qsort(motifsDNA, numMotifs,sizeof(char*),compareMotifsFH2);
		FILE* fp = fopen(MER_MOTIF_FILE, "w");
		for (int i=0;i<numMotifs;i++) {
			if (i>= 10000) {
				break;
			}
			float pvalue = possible->search(motifsDNA[i]);
			fprintf(fp,"%s\t%e\n",motifsDNA[i],pvalue);
		}
		fclose(fp);
		delete possible;
		numSeeds = numMotifs;
		return motifsDNA;
	} else {
		qsort(motifsProtein, numMotifs, sizeof(ProteinStat), sortProteinStat);
		FILE* fp = fopen(MER_MOTIF_FILE, "w");
		for (int i=0;i<numMotifs;i++) {
			if (i>= 10000) {
				break;
			}
			char* str = getAAStringfromIndex(motifsProtein[i].v);
			fprintf(fp,"%s\t%e\n",str,motifsProtein[i].pvalue);
			delete []str;
		}
		fclose(fp);
		numSeeds = numMotifs;
		return (char**) motifsProtein;
			
	}
}

ProteinStat::ProteinStat() {
	v = NULL;
	pvalue = 1.0;
}
ProteinStat::~ProteinStat() {
	if (v != NULL) {
		delete []v;
	}
	v= NULL;
}

char* revopp(char* motif) {
	int slen = strlen(motif);
	char* rev = new char[slen+1];
	if (motif == NULL || rev == NULL) {
		//fprintf(stderr, "%x = rev\n%x = motif\n",rev,motif);
		fprintf(stderr, "Got a bug somewhere :)\n");
		exit(0);
	}
	rev[slen] = '\0';
	for (int i=0;i<slen;i++) {
		switch(motif[slen-1-i]) {
			case 'A':
				rev[i] = 'T';
				break;
			case 'C':
				rev[i] = 'G';
				break;
			case 'G':
				rev[i] = 'C';
				break;
			case 'T':
				rev[i] = 'A';
				break;
			case 'N':
				rev[i] = 'N';
				break;
			case 'R':
				rev[i] = 'Y';
				break;
			case 'Y':
				rev[i] = 'R';
				break;
			case 'M':
				rev[i] = 'K';
				break;
			case 'K':
				rev[i] = 'M';
				break;
			case 'W':
				rev[i] = 'S';
				break;
			case 'S':
				rev[i] = 'W';
				break;
			case 'B':
				rev[i] = 'V';
				break;
			case 'V':
				rev[i] = 'B';
				break;
			case 'D':
				rev[i] = 'H';
				break;
			case 'H':
				rev[i] = 'D';
				break;
			default:
				fprintf(stderr, "Problem with getting the reverse opposite!!!\n");
				rev[i] = 'N';
		}
	}
	return rev;
}

PSSM::PSSM() {
	init();
}

void PSSM::init() {
	matrix = NULL;
	state= 0;
	length =0;
	threshold = 0;
	name = NULL;
	outputString = new char[1];
	outputString[0] = '\0';
	gapLengths = NULL;
	numGapLengths = 0;
	v = 1.0;
	pvalue = 1.0;
	size = 4;
	type = 0;
	backFreq = new float[size];
	for (int i=0;i<size;i++) {
		backFreq[i] = 1.0 / (float)size;
	}
}
void PSSM::setBackFreq(float *freq) {
	float def = 1.0 / ((float)size);
	if (backFreq != NULL) {
		delete []backFreq;
	}
	backFreq = new float[size];
	if (freq != NULL) {
		for (int i=0;i<size;i++) {
			backFreq[i] = freq[i];
		}
	} else {
		for (int i=0;i<size;i++) {
			backFreq[i] = def;
		}
	}
}
void PSSM::setOutputString(char* s) {
	if (outputString != NULL) {
		delete []outputString;
	}
	outputString = new char[strlen(s)+1];
	strcpy(outputString, s);
}
PSSM::PSSM(PSSM* pssm) {
	init();	
	length=pssm->length;
	v = pssm->v;
	size = pssm->size;
	type = pssm->type;
	numGapLengths = pssm->numGapLengths;
	gapLengths = new int[numGapLengths];
	for (int i=0;i<numGapLengths;i++){ 
		gapLengths[i] = pssm->gapLengths[i];
	}
	pvalue =pssm->pvalue;
	threshold=pssm->threshold;
	setBackFreq(pssm->backFreq);
	matrix = new float*[length];
	state = pssm->state;
	for (int i=0;i<length;i++) {
		matrix[i] = new float[size];
		for (int j=0;j<size;j++) {
			matrix[i][j] = pssm->matrix[i][j];
		}
	}
}
void PSSM::setGapInfo(int* newGapLengths, int newNumGapLengths) {
	if (gapLengths != NULL) delete []gapLengths;
	if (newNumGapLengths == 0) {
		gapLengths = new int[1];
		numGapLengths = 1;
		gapLengths[0] = 0;
		return;
	}
	numGapLengths = newNumGapLengths;
	gapLengths = new int[numGapLengths];
	for (int i=0;i<numGapLengths;i++) {
		gapLengths[i] = newGapLengths[i];
	}
}
void PSSM::setName(char* newname) {
	if (newname == NULL) return;
	name = new char[strlen(newname)+1];
	strcpy(name, newname);
}
void PSSM::setSize(int seqType) {
	type = seqType;
	if (seqType == DNA_SEQ_TYPE) {
		size = 4;
	} else if (seqType == PROTEIN_SEQ_TYPE) {
		size = 20;
	}
}
PSSM::PSSM(char* newname, float** profile, int plength, int seqType) {
	init();
	setSize(seqType);
	setName(newname);
	length = plength;
	//fprintf(stderr, "%d\n", length);
	matrix = new float*[length];
	for (int i=0;i<length;i++) {
		matrix[i] = new float[size];
		for (int j=0;j<size;j++) {
			matrix[i][j] = profile[i][j];
		}
	}
}

PSSM::~PSSM() {
	if (matrix != NULL) {
		for (int i=0;i<length;i++) {
			if (matrix[i] != NULL) {
				delete []matrix[i];
			}
		}
		delete []matrix;
	}
	if (outputString != NULL) {
		delete []outputString;
	}
	if (name != NULL) {
		delete []name;
	}
	if (gapLengths != NULL) {
		delete []gapLengths;
	}
	if (backFreq != NULL) {
		delete []backFreq;
	}
}
PSSM::PSSM(char* motif,int seqType) {
	init();
	set(motif,seqType);
}
PSSM::PSSM(int* motif,int seqType) {
	init();
	set(motif,seqType);
}
void PSSM::set(char* motif, int seqType) {
	int* values = seqIndexOfIUPAC(seqType, motif, NULL);
	set(values, seqType);
	delete []values;
}

void PSSM::set(int* values, int seqType) {
	/****************************
		settings
	*****************************/
	int buffer = PROFILE_EDGE_SPACE;
	float initvalue = MOTIF_INIT_VALUE;	
	setSize(seqType);
	setBackFreq(NULL);
	//fprintf(stderr, "Setting profile to %s\n", motif);

	length = 0;
	while (values[length] != 0) {
		length++;
	}
	length += 2*buffer;
	
	matrix = new float*[length];
	for (int i=0;i<length;i++) {
		matrix[i] = new float[size];	
		if (i<buffer || i>=length-buffer) {
			for (int j=0;j<size;j++) {
				matrix[i][j] = 1/size;
			}
		} else {
			for (int j=0;j<size;j++) {
				matrix[i][j] = (1.0-initvalue)/((float)size-1.0);
			}
			for (int j=0;j<size;j++) {
				int mask = 1 << j;
				if (values[i-buffer] & mask) {
					matrix[i][j] = initvalue;
				}
			}
		}
	}
	adjustvalues();
}

void PSSM::print(FILE* fp) {
	if (matrix == NULL) {
		fprintf(fp, "null matrix\n");
		return;
	}
	char* m = consensus();
	char* n = name;
	if (name == NULL) {
		n = m;
	}
	
	char* gapinfostr = makeVariableRangeString(gapLengths, numGapLengths);
	if (type == DNA_SEQ_TYPE) fprintf(fp, ">");
	if (type == PROTEIN_SEQ_TYPE) fprintf(fp, "<");
	
	fprintf(fp,"%s\t%s\t%f\t%e\t%s\t%s\n",m,n,threshold,pvalue,gapinfostr,outputString);
	delete []gapinfostr;
	delete []m;
	for (int i=0;i<length;i++) {
		for (int j=0;j<size;j++) {
			if (j>0) {
				fprintf(fp, "\t");
			}
			fprintf(fp, "%.3f", matrix[i][j]);
		}
		fprintf(fp, "\n");
	}
}
void PSSM::printMatrix(FILE* fp) {
	fprintf(stderr, ">>>>>\n");
	for (int i=0;i<length;i++) {
		for (int j=0;j<size;j++) {
			if (j>0) {
				fprintf(fp, "\t");
			}
			fprintf(fp, "%.3f", matrix[i][j]);
		}
		fprintf(fp, "\n");
	}
	
}
void PSSM::logXform() {
	if (state) return;
	adjustvalues();
	for (int i=0;i<length;i++) {
		for (int j=0;j<size;j++) {
			matrix[i][j] = (float)log((double)matrix[i][j]/ backFreq[j]);
		}
	}
	state = 1;
}
void PSSM::expXform() {
	if (!state) return;
	for (int i=0;i<length;i++) {
		for (int j=0;j<size;j++) {
			matrix[i][j] = (float)exp((double)matrix[i][j])*(backFreq[j]);
			//if (matrix[i][j] != matrix[i][j]) matrix[i][j] = 1.0;
		}
	}
	state = 0;
}
void PSSM::adjustPrint() {
	expXform();
	adjustvalues();
}

void PSSM::adjustvalues() {
	expXform();
	for (int i=0;i<length;i++) {
		float total = 0;
		float min = 0;
		//in case a value is less than zero...
		for (int j=0;j<size;j++) {
			if (matrix[i][j] < min) {
				min = matrix[i][j];
			}
		}
		for (int j=0;j<size;j++) {
			matrix[i][j] -= min;
			total += matrix[i][j];
		}
		int under = 0;
		float totalOver = 0;
		for (int j=0;j<size;j++) {
			matrix[i][j] /= total;
			if (matrix[i][j] < MIN_MATRIX_PROB) {
				under++;
			} else {
				totalOver+=matrix[i][j];
			}
		}
		//float newtotalOver = totalOver - under*MIN_MATRIX_PROB;
		for (int j=0;j<size;j++) {
			if (matrix[i][j] < MIN_MATRIX_PROB) {
				matrix[i][j] = MIN_MATRIX_PROB;
			} else {
				matrix[i][j] *= (totalOver-under*MIN_MATRIX_PROB)/(totalOver);
			}
		}
	}
}

char* PSSM::consensus() {
	if (matrix == NULL) return NULL;
	
	char* motif = new char[length+1];
	static char DNA[5] = "ACGT";
	static char PROTEIN[21] = "ACDEFGHIKLMNPQRSTVWY";
	motif[length] = '\0';
	for (int i=0;i<length;i++) {
		int best = 0;
		float bestval = matrix[i][0];
		for (int j=1;j<size;j++) {
			if (matrix[i][j] > bestval) {
				bestval = matrix[i][j];
				best = j;
			}
		}
		if (size == 4) {
			if (bestval <= .35) {
				motif[i] = 'N';
				continue;
			}
			motif[i]=DNA[best];
		} else if (size == 20) {
			if (bestval <= .1) {
				motif[i] = 'X';
				continue;
			}
			motif[i]=PROTEIN[best];
		}
	}
	return motif;
}
void PSSM::mildize(float factor) {
	logXform();
	for (int i=0;i<length;i++) {
		for (int j=0;j<size;j++) {
			matrix[i][j] /= factor;
		}
	}
	expXform();
	adjustvalues();
}

PSSM* PSSM::copy() {
	PSSM* pssm = new PSSM(name, matrix, length, type);
	pssm->v = v;
	pssm->type = type;
	pssm->size = size;
	pssm->threshold = threshold;
	pssm->pvalue = pvalue;
	pssm->score  = score;
	pssm->state = state;
	pssm->constatID = constatID;
	pssm->setBackFreq(backFreq);
	if (numGapLengths > 0) pssm->setGapInfo(gapLengths, numGapLengths);
	pssm->setOutputString(outputString);
	return pssm;
}


#define MAXPLEN 100

void readMatrix(CommandLine* cmd, PSSM** &allMatrix, int &numMatrix) {

	FILE* fp = fopen(cmd->motiffile,"r");
	if (fp == NULL) {
		fprintf(stderr, "!!! Error: Cannot open motif file - %s\n",cmd->motiffile);
		exit(0);
	}
	numMatrix = 0;	
	char* buf = new char[BSIZE];
	//char* line = NULL;
	char* name = new char[BSIZE];
	char* outstr = new char[BSIZE];
	
	char** cols = new char*[10000];
	int numCols = 0;

	while (fgets(buf,BSIZE,fp) != NULL) {
		if (buf[0] == 'T' || buf[0] == '>') {
			numMatrix++;
		}
	}
	rewind(fp);

	float** profile = new float*[BSIZE];
	for (int i=0;i<BSIZE;i++) profile[i] = new float[20];
	int plen = 0;
	allMatrix = new PSSM*[numMatrix];
	numMatrix = 0;
	float nThreshold=0.0;
	float nlogP=0.0;
	//float v=0.0;
	int size=4;
	int* gapL = NULL;
	int numGL = 0;

	while (fgets(buf,BSIZE,fp) != NULL) {
		buf[strlen(buf)-1] = '\0'; //get rid of newline character

		split(buf, cols, numCols, '\t');

		if (cols[0][0] == '<' || cols[0][0] == '>') {
			if (plen > 0) {
				int type = DNA_SEQ_TYPE;
				//if (size == 4) type = DNA_SEQ_TYPE;
				if (size == 20) type = PROTEIN_SEQ_TYPE;
				allMatrix[numMatrix] = new PSSM(name,profile,plen,type);
				allMatrix[numMatrix]->state = 0;
				allMatrix[numMatrix]->expXform();
				allMatrix[numMatrix]->adjustvalues();
				allMatrix[numMatrix]->threshold = nThreshold;
				allMatrix[numMatrix]->setOutputString(outstr);
				allMatrix[numMatrix]->pvalue = nlogP;
				allMatrix[numMatrix]->setGapInfo(gapL, numGL);
				numMatrix++;
			}
			if (cols[0][0] == '<') size = 20;
			if (cols[0][0] == '>') size = 4;
			plen = 0;
			gapL = NULL;
			numGL = 0;
			nThreshold=0.0;
			nlogP=0.0;
			//v=0.0;
			numGL = 0;
			name[0] = '\0';
			outstr[0] = '\0';

			if (numCols > 1) strcpy(name,cols[1]);
			if (numCols > 2) sscanf(cols[2],"%f", &nThreshold);
			if (numCols > 3) sscanf(cols[3],"%f", &nlogP);
			if (numCols > 4) gapL = parseRangeVariable(cols[4], numGL);
			if (numCols > 5) strcpy(outstr,cols[5]);


			if (cols[0][0] == 'T') {
				strcpy(name,cols[0]);
			}
		} else {
			//char* val = NULL;
			if ((size == 4 && numCols > 3) || (size==20 && numCols > 19)) {
				for (int i=0;i<size;i++) {
					sscanf(cols[i],"%g",&(profile[plen][i]));
				}
				plen++;
			}
		}
	}
	if (plen > 0) {
		int type = DNA_SEQ_TYPE;
		if (size == 4) type = DNA_SEQ_TYPE;
		if (size == 20) type = PROTEIN_SEQ_TYPE;
		allMatrix[numMatrix] = new PSSM(name,profile,plen,type);
		allMatrix[numMatrix]->state = 0;
		allMatrix[numMatrix]->expXform();
		allMatrix[numMatrix]->adjustvalues();
		allMatrix[numMatrix]->threshold = nThreshold;
		allMatrix[numMatrix]->setOutputString(outstr);
		allMatrix[numMatrix]->pvalue = nlogP;
		allMatrix[numMatrix]->setGapInfo(gapL, numGL);
		numMatrix++;
	}

	for (int i=0;i<BSIZE;i++) {
		delete [](profile[i]);
	}
	delete []profile;
	delete []buf;
	delete []name;
	delete []outstr;
	delete []cols;
	fclose(fp);
}


Site::Site() {
	init();
}
Site::Site(char* geneid, char* nseq, PSSM* nmotif, int ndir, int noffset, 
					int ngapsize, float ncons, float nratio, float nscore) {
	init();
	if (nseq != NULL){ 
		setSeq(nseq);
	}
	if (geneid != NULL) {
		if (geneID != NULL) delete []geneID;
		geneID = new char[strlen(geneid)+1];
		strcpy(geneID, geneid);
	}
	motif = nmotif;
	dir = ndir;
	offset = noffset;
	conservation = ncons;
	ratio = nratio;
	score = nscore;
	gapsize = ngapsize;
	gapoffset = motif->length / 2;

}
void Site::setSeq(char* nseq) {
	if (seq != NULL) delete []seq;
	seq = new char[strlen(nseq)+1];
	strcpy(seq, nseq);
}
void Site::init() {
	seq = NULL;
	motif = NULL;
	geneID = NULL;
	dir = 0;
	offset = 0;
	conservation = 0.0;
	ratio = 0.0;
	score = 0.0;
	gapsize = 0;
	gapoffset = 0;
}
Site::~Site() {
	if (seq != NULL) delete []seq;
	if (geneID != NULL) delete []geneID;
}
void Site::print(FILE* fp) {
	if (dir == 0) {
		fprintf(fp,"%s\t%d\t%s\t%f\t+\t%s\t%f\n", geneID, offset, seq, conservation, motif->name, score); 
	} else {
		fprintf(fp,"%s\t%d\t%s\t%f\t-\t%s\t%f\n", geneID, offset, seq, conservation, motif->name, score); 
	}
}

Sequence::Sequence() {
	stat=0;
	seq = NULL;
	w = 0;
	name=NULL;
}
Sequence::~Sequence() {
	if (seq != NULL) {
		delete []seq;
	}
	seq = NULL;
	if (name != NULL) {
		delete []name;
	}
	name = NULL;
}
void Sequence::setSeq(char* newSeq) {
	if (newSeq == NULL) return;
	if (seq != NULL) {
		delete []seq;
	}
	seq = new char[strlen(newSeq)+1];
	strcpy(seq, newSeq);
}
void Sequence::setName(char* newName) {
	if (newName == NULL) return;
	if (name != NULL) {
		delete []name;
	}
	name = new char[strlen(newName)+1];
	strcpy(name, newName);
}


Sequence* loadSequences(CommandLine* cmd, int &numSeqs) {

	int group = 0;
	char* name = NULL;
	char* buf = new char[BUFFER];

	//float* conservation = new float[BUFFER];
	Floattable* geneWeights = new Floattable();
	Inttable* duplicates = new Inttable();

	int numStatGenes = 0;

	char** cols = new char*[10000];
	int numCols = 0;
	int checkInt = 0;

	//load group file
	Inttable* geneID = NULL;
	if (cmd->statfile != NULL) {
		geneID = new Inttable();

		FILE* fp = fopen(cmd->statfile,"r");
		if (fp == NULL) {
			fprintf(stderr, "Can't open group file\n");
			exit(0);
		}
		while (fgets(buf,BUFFER,fp) != NULL) {

			float w = 1.0;
			split(buf, cols, numCols, '\t');
			name = cols[0];
			if (numCols < 2) {
				continue;
			}
			
			group = 0;
			sscanf(cols[1], "%d", &group);
			if (cmd->weightFlag && numCols > 2) {
				sscanf(cols[2], "%f", &w);
			}

			checkInt = duplicates->search(name);
			if (checkInt == EMPTY_INT) {
				duplicates->insert(1,name);
			} else {
				checkInt++;
				if (group != 0) {
					fprintf(stderr, "!!! Duplicate Identifier (%s) - found %d times !!! (dangerous - use unique IDs!!)\n",name,checkInt);
				}
				duplicates->insert(checkInt,name);
				continue;
			}
		

			geneID->insert(group,name);
			geneWeights->insert(w, name);
			numStatGenes++;
		}       
		fclose(fp);
	}
	delete duplicates;

	
	FILE* fp = fopen(cmd->seqfile,"r");
	if (fp == NULL) {
		fprintf(stderr, "Could not open sequence file - need sequence file\n");
		exit(0);
	}

	int maxSeq = numStatGenes;
	if (cmd->statfile == NULL) {
		maxSeq = MAX_NUM_SEQUENCES;
	}

	numSeqs = 0;
	Sequence* sequences = new Sequence[maxSeq];

	duplicates = new Inttable();

	while (fgets(buf,BUFFER,fp) != NULL) {
		split(buf, cols, numCols, '\t');
		if (numCols < 2) continue;

		name = cols[0];
		char* seq = cols[1];

		//check if sequence is in the group file
		int stat = 0;
		if (geneID != NULL) {
			stat = geneID->search(name);
			if (EMPTY_INT == stat) continue;
		}

		checkInt = duplicates->search(name);
		if (checkInt == EMPTY_INT) {
			duplicates->insert(1,name);
		} else {
			duplicates->insert(checkInt++,name);
			//fprintf(stderr, "!!! Duplicate Sequence ID (%s) Very bad!!!, Skipping...\n", name);
			continue;
		}

		//check if too many Ns are present
		int Ncount= 0;
		int nindex = 0;
		int seqlen = strlen(seq);
		while (nindex < seqlen) {
			if (cmd->seqType == DNA_SEQ_TYPE && seq[nindex] == 'N') Ncount++;
			if (cmd->seqType == PROTEIN_SEQ_TYPE && seq[nindex] == 'X') Ncount++;
			nindex++;
		}
		if (((float)Ncount)/(float)seqlen > cmd->nratiocutoff) continue;
		if (seqlen > cmd->maxSeqLen) continue;
		if (seqlen < cmd->minSeqLen) continue;
	
		float w = geneWeights->search(name);

		sequences[numSeqs].setSeq(seq);
		sequences[numSeqs].setName(name);
		sequences[numSeqs].w=w;
		sequences[numSeqs].stat=stat;
		numSeqs++;
	}
	delete duplicates;

	fprintf(stderr, "%d total sequences loaded\n", numSeqs);
	if (geneID != NULL) {
		delete geneID;
	}
	delete geneWeights;
	delete []buf;
	delete []cols;

	return sequences;
}


void findMotifSites(CommandLine* cmd, Sequence* seqs, int numSeqs, PSSM* motif, FILE* outfile) {


	Floattable* geneScores = new Floattable();
	Floattable* geneWeights = new Floattable();
	Inttable* geneID = new Inttable();

	int numGenes = 0;
	int numPosGenes = 0;
	float totalWGenes = 0.0;
	float totalWPosGenes = 0.0;
	float totalWNegGenes = 0.0;
	float totalSeqLen = 0.0;
	float totalPosSeqLen = 0.0;
	float totalNegSeqLen = 0.0;

	for (int j=0;j<numSeqs;j++) {

		float w = seqs[j].w;
		char* name = seqs[j].name;
		char* seq = seqs[j].seq;
		int stat = seqs[j].stat;
		int seqlen = strlen(seq);


		int testInt = geneID->search(name);
		if (testInt != EMPTY_INT) {
			if (stat != 0) {
				fprintf(stderr, "!!! Found dupliate sequence Identifier (%s) !!! Skipping...\n",name);
			}
			continue;
		}
		numGenes++;
		totalWGenes += w;
		if (stat == 1) {
			numPosGenes++;
			totalWPosGenes+=w;
			totalPosSeqLen += seqlen;
		} else {
			totalWNegGenes+=w;
			totalNegSeqLen += seqlen;
		}
		totalSeqLen  += seqlen;
			
		
		geneID->insert(stat,name);
		geneWeights->insert(w, name);


		//fprintf(stderr, "%d\t%s\t", i, gene[i].name);
		if ( numGenes % 1000 == 0) fprintf(stderr, "%d\t", numGenes);

		Site** sites = NULL;
		int numSites = 0;
		int eraseFlag = 0;
		if (cmd->action == ACTION_REMOVE) eraseFlag = 1;

		Site* bestSite = scoreMatrix(cmd, motif, name, seq,
							eraseFlag, sites, numSites);

		if (cmd->action == ACTION_GENESCORE) {
			//fprintf(stdout, "%s\t%s\t%f\t%s\n", name, motif->name, bestSite->score,bestSite->seq);
			bestSite->print(outfile);
		}
		if (cmd->action == ACTION_FIND) {
			for (int i=0;i<numSites;i++) {
				if (outfile != NULL) sites[i]->print(outfile);
			}
		} else if (cmd->action == ACTION_REMOVE) {
			if (outfile != NULL) fprintf(outfile, "%s\t%s\n", name, seq);
		}
		geneScores->insert(bestSite->score, name);
		if (sites != NULL) {
			for (int i=0;i<numSites;i++) {
				delete sites[i];
			}
			delete []sites;
		}
		delete bestSite;
	}
	fprintf(stderr, "\n");

	float autoscale = 1.0;
	float* approxIn = NULL;
	float* approxOut = NULL;
	int automax = (int)(totalWNegGenes*10);

	if (cmd->autoScale && cmd->action != ACTION_FIND && cmd->action != ACTION_REMOVE 
				&& cmd->action != ACTION_GENESCORE) {
		float psitesPerGene = totalPosSeqLen / totalWPosGenes;
		float nsitesPerGene = totalNegSeqLen / totalWNegGenes;
		autoscale = psitesPerGene / nsitesPerGene;
		fprintf(stderr, "Target Ratio: %e\n", psitesPerGene);
		fprintf(stderr, "Backgd Ratio: %e\n", nsitesPerGene);
		fprintf(stderr, "Autoscale adjustment factor: %f\n", autoscale);
		approxIn = new float[automax];
		approxOut = new float[automax];
		approxIn[0] = 0;
		approxIn[1] = 1;
		approxOut[0] = 0;
		approxOut[1] = 1;
		for (int i=2;i<automax;i++) {
			approxOut[i] = approxOut[i-1]+(totalWNegGenes-approxOut[i-1])/(totalWNegGenes);
			int index = (int) round(approxOut[i]);
			if (index >= automax) {
				approxIn[automax] = i;
			} else {
				approxIn[index] = i;
			}
		}
	}
	
	if (cmd->statmemory == NULL && cmd->action != ACTION_FIND && cmd->action != ACTION_REMOVE 
				&& cmd->action != ACTION_GENESCORE) {
		cmd->statmemory = new StatMemory(MAX_STAT_CACHE, (int)round(totalWGenes),
					(int)round(totalWPosGenes));
	}


	if (geneID != NULL && cmd->action != ACTION_FIND && cmd->action != ACTION_REMOVE
			&& cmd->action != ACTION_GENESCORE) {
		sorttable2 = geneScores;
		char** keys = geneScores->keys();
		int numkeys = geneScores->total;
		qsort(keys, numkeys,sizeof(char*),compareMotifsFH2R);
		float N = 0;
		float m = 0;
		float n = 0;
		float bestP = 1e100;
		float bestT = 0;
		int bestN = 0;
		int bestn =0;
		
		for (int i=0;i<numkeys;i++) {
			int stat = geneID->search(keys[i]);
			float w = geneWeights->search(keys[i]);
			N+=w;
			if (stat == 1) {
				n+=w;
			} else {
				m+=w;
			}
			float scaledNegGenes = m;
			if (cmd->autoScale) {
				int negcount = (int)m;
				if (negcount >= automax) {
					negcount = automax;
				}
				float estimateSites = approxIn[negcount];
				int intSites = (int)round(estimateSites *autoscale);
				if (intSites >= automax) {
					intSites = automax;
				}
				scaledNegGenes = (float)approxOut[intSites];
			}
	
			float logp = cmd->statmemory->getStat(
							(int)round(n+scaledNegGenes),(int)round(n));
			//fprintf(stderr, "%s\t%d\t%f\t%f\n", keys[i], stat, geneScores->search(keys[i]), logp);
			if (logp < bestP || (cmd->action==ACTION_GETPVALUE && 
									geneScores->search(keys[i]) > motif->threshold)) {
				bestP = logp;
				bestT = geneScores->search(keys[i]);
				bestN = (int)round(N);
				bestn = (int)round(n);
			}

			if (N-n > (totalWGenes-totalWPosGenes)/2 || (cmd->action==ACTION_GETPVALUE && 
									geneScores->search(keys[i]) < motif->threshold)) {
				//fprintf(stderr,"Present in more than half of background...\n");
				break;
			}
		}
		for (int i=0;i<numkeys;i++) {
			if (keys[i] != NULL) {
				delete [](keys[i]);
				keys[i] = NULL;
			}
		}
		delete []keys;
		motif->threshold = bestT-0.001;
		motif->pvalue= bestP;
		char* outstr = new char[1000];
		sprintf(outstr,"%d,%d,%d,%d,%e",(int)round(totalWGenes),
						(int)round(totalWPosGenes),bestN,bestn,exp(bestP));
		motif->setOutputString(outstr);
		fprintf(stderr, "%s\t%s\n", motif->name, motif->outputString);
	}
	delete geneScores;
	delete geneWeights;
	delete geneID;
	if (cmd->autoScale) {
		delete []approxIn;
		delete []approxOut;
	}
}


Site* scoreMatrix(CommandLine* cmd, PSSM* motif, char* id, char* seqStr, int eraseFlag,
								Site** &sites, int &numsites) {
	int seqLen = strlen(seqStr);
	char* seq = NULL;
	char* rseq = NULL;
	int plen = motif->length;
	float threshold = motif->threshold;
	char* merSeq = new char[motif->gapLengths[motif->numGapLengths-1]+plen+1];
	indexSeq(seqStr, seq, rseq, cmd->seqType);

	Site** leadSites = NULL;
	int numLeadSites = 0;
	Site* bestLeadSite = scoreSeqWithMatrix(cmd,motif,id,seq,
							seqLen, leadSites, numLeadSites);
	delete []seq;
	Site** lagSites = NULL;
	int numLagSites = 0;
	Site* bestLagSite = NULL;

	if (cmd->revoppFlag) {
		bestLagSite = scoreSeqWithMatrix(cmd,motif,id,rseq,
								seqLen, lagSites, numLagSites);
	}
	if (rseq != NULL) delete []rseq;
	numsites = numLeadSites + numLagSites;
	sites = new Site*[numsites];
	int indexSites = 0;
	
	for (int i=-1;i<numLeadSites;i++) {
		Site* current = NULL;
		if (i<0) current = bestLeadSite;
		else current = leadSites[i];

		current->dir = 0;
		int pos = current->offset;
		int gapsize= current->gapsize;
		//int gapoffset = current->gapoffset;
		for (int j=0;j<plen+gapsize;j++) {
			merSeq[j] = seqStr[pos+j];
		}
		merSeq[plen+gapsize] = '\0';
		current->offset += cmd->offset;
		current->setSeq(merSeq);			
		if (current->score > threshold && eraseFlag) {
			eraseSequence(seqStr,pos,plen+gapsize,cmd->seqType);
		}
		if (i >= 0) {
			sites[indexSites++] = current;
		}
	}
	delete []leadSites;

	if (cmd->revoppFlag) {
		for (int i=-1;i<numLagSites;i++) {
			Site* current = NULL;
			if (i<0) current = bestLagSite;
			else current = lagSites[i];
	
			current->dir = 1;
			int gapsize= current->gapsize;
			int pos = seqLen - 0 - current->offset - gapsize - plen;
			//int gapoffset = current->gapoffset;
			for (int j=0;j<plen+gapsize;j++) {
				merSeq[j] = seqStr[pos+j];
			}
			merSeq[plen+gapsize] = '\0';
			current->offset = cmd->offset + pos;
			current->setSeq(merSeq);			
			if (current->score > threshold && eraseFlag) {
				eraseSequence(seqStr,pos,plen+gapsize,cmd->seqType);
			}
			if (i >= 0) {
				sites[indexSites++] = current;
			}
		}
		delete []lagSites;
	}

	Site* bestSite = bestLeadSite;
	if (cmd->revoppFlag) {
		if (bestLeadSite->score > bestLagSite->score) {
			delete bestLagSite;
			bestSite = bestLeadSite;
		} else {
			delete bestLeadSite;
			bestSite= bestLagSite;
		}
	}
	return bestSite;

}
//seq for this function is in number form i.e 0=A,1=C,2=G,3=T
Site* scoreSeqWithMatrix(CommandLine* cmd, PSSM* motif, char* id, char* seq, 
								int seqLen,Site** &sites, int &numsites) {

	motif->logXform();
	int plen = motif->length;
	float** matrix = motif->matrix;

	float bestScore = -1e10;
	int bestPos = 0;
	int bestGap = 0;
	int gapoffset = plen/2;
	char* mer = new char[plen+1];
	mer[plen]='\0';
	char defSeq = '\0';
	float threshold = motif->threshold;
	float consScore = 0.0;
	float ratioScore = 1.0;

	numsites = 0;
	sites = new Site*[seqLen];

	for (int i=0;i<=seqLen-plen;i++) {
		for (int j=0;j<motif->numGapLengths;j++) {
			int gapsize = motif->gapLengths[j];
			if (i+plen+gapsize > seqLen) continue;

			if (gapsize ==0) {
				for (int j=0;j<plen;j++) {
					mer[j] = seq[i+j];
				}
			} else if (gapsize > 0) {
				for (int j=0;j<gapoffset;j++) {
					mer[j] = seq[i+j];
				}
				for (int j=gapoffset;j<plen;j++) {
					mer[j] = seq[i+j+gapsize];
				}
			}
			float score = 0.0;
			int bad = 0;
			for (int j=0;j<plen;j++) {
				if (cmd->seqType == DNA_SEQ_TYPE && mer[j] == 4) {
					bad=1;
					break;
				} else if (cmd->seqType == PROTEIN_SEQ_TYPE && mer[j] == 20) {
					bad = 1;
					break;
				} else {
					score += matrix[j][(int)mer[j]];
				}
			}
			if (bad) continue;

			if (score > threshold) {
				sites[numsites] = new Site(id, &defSeq, motif, 0, i,gapsize,
										consScore, ratioScore, score);
				numsites++;
			}
			if (score > bestScore) {
				bestScore = score;
				bestPos = i;
				bestGap = gapsize;
			}
		}
	}

	Site *bestSite = new Site(id,&defSeq, motif, 0, bestPos, bestGap, 
										consScore,ratioScore,bestScore);
	
	Site** nsites = NULL;
	if (numsites > 0) {
		nsites = new Site*[numsites];
		for (int i=0;i<numsites;i++) {
			nsites[i] = sites[i];
		}
	}
	delete []sites;
	sites = nsites;

	delete []mer;
	return bestSite;
}


void indexSeq(char* seqStr, char* &seq, char* &rseq, int seqType) {
	if (seqType == DNA_SEQ_TYPE) indexDNA(seqStr, seq, rseq);
	if (seqType == PROTEIN_SEQ_TYPE) indexProtein(seqStr, seq, rseq);
}
void indexProtein(char* seqStr, char* &seq, char* &rseq) {
	seq = seqIndexOf(PROTEIN_SEQ_TYPE, seqStr, NULL);
	rseq = NULL;
}

void indexDNA(char* seqStr, char* &seq, char* &rseq) {
	int seqLen = strlen(seqStr);
	seq = new char[seqLen+1];
	int len = 0;
	for (int i=0;i<seqLen;i++) {
		if (seqStr[i] == 'A') {
			seq[len++] = 0;
		} else if (seqStr[i] == 'C') {
			seq[len++] = 1;
		} else if (seqStr[i] == 'G') {
			seq[len++] = 2;
		} else if (seqStr[i] == 'T') {
			seq[len++] = 3;
		} else {
			seq[len++] = 4;
		}
	}
	seqLen = len;
	seq[seqLen] = 5;
	rseq = new char[seqLen+1];
	for (int i=seqLen-1;i>=0;i--) {
		if (seq[i] == 4) rseq[seqLen-i-1] = 4;
		else rseq[seqLen-i-1] = 3-seq[i];
	}
	rseq[seqLen] = '\0';
	seq[seqLen] = '\0';
}
void eraseSequence(char* seq, int offset, int length, int seqType) {
	int seqLen = strlen(seq);
	//fprintf(stderr, "%d\t%d\n", offset, length);
	for (int i=offset;i<offset+length;i++) {
	//fprintf(stderr, "%c", seq[i]);
		if (seqLen <= i) break;
		if (i < 0) continue;
		if (seqType == DNA_SEQ_TYPE) seq[i] = 'N';
		if (seqType == PROTEIN_SEQ_TYPE) seq[i] = 'X';
	}
}


//========================= End of file =============================
