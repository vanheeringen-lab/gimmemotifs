#include <stdio.h>
#include <math.h>
#include <string.h>
#include <strings.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>

#ifndef SEQUENCE_H
#define SEQUENCE_H

#include "Hashtable.h"
#include "statistics.h"
#include "CommandLine.h"
//#include <Motif.h>

#define BSIZE 10000
#define MAXPLEN 100
#define BAD_N_RATIO 0.75
#define DNA_SEQ 0
#define PROTEIN_SEQ 1


class Gene {
public:
	char* name;
	char* seqStr;
	int seqLen;
	int type;
	int elen;
	char* seq;
	char* rseq;
	float* expression;
	float* conservation;
	float* pssmCons;
	int numCons;
	float exp;
	int stat;
	int numExps;

	float* svmSeq;

	Gene();
	void init();
	~Gene();
	void print (FILE* fp);
	void indexSeq();
	void copy(Gene* g);
	void capitalize();
	int findMotif(char*);
	void randomizeSeq(float* freq);
	void eraseSequence(int offset, int length);
	void printSeq(FILE* fp);
	//float scoreMatrix(Matrix* m);
	static char* revopp(char*);
	void createSVMSeq();
	void fakeConservation();
};

void readSeq(char* file, Gene* &genes, int& numGenes,int seqType);
void readStat(char* file, Gene* &genes, int& numGenes,int prune);
void readExp(char* stat, Gene* &genes, int& numGenes,int prune);
void readCons(char* stat, Gene* &genes, int& numGenes,int prune);

void copystr(char* &dest, char* src);
float* getBackgroundFrequency(Gene* genes, int numGenes);
float alignMotifs(char* m1, char* m2, int &offset, int &edge);
float revAlignMotifs(char* m1, char* m2, int &offset, int &edge,int &orientation);

void createCorrelationDistribution(char* filename, Gene* genes, int numGenes,
            Gene*** regGenes, int* numRegGenes, int numMotifs);
Gene* createRandomizedGenes(Gene* genes, int numGenes);
void calcPSSMCons(Gene* genes, int numGenes, float factor);
float alignSeq(char* seq1, char* seq2,int localFlag,  float match, float missmatch,
                                        float gapStart, float gapContinue);
void filterDegenerateGenes(Gene* &genes, int &numGenes, int type, float nratiocutoff);

#endif
