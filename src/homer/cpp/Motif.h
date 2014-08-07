
// Copyright 2009, 2010, 2011, 2012 Christopher Benner <cbenner@gmail.com>
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


#include <stdio.h>
#include "Hashtable.h"
#include "CommandLine.h"
#include "statistics.h"
#include <limits.h>
#include <math.h>
#include <float.h>
#include <time.h>

#ifndef MOTIF_H
#define MOTIF_H

#define HYPERGEO 1
#define PROFILE 2
#define GENERATIVE 3

#define MIN_MOTIF_SIZE 6
#define MIN_ABSOLUTE_NUM_MOTIFS 3
#define MIN_RELATIVE_RATIO_MOTIFS 0

#define STAT_THRESHOLD .5
#define PROFILE_EDGE_SPACE 0
#define NUM_TOP_SEEDS 50
#define NUM_TOP_MOTIFS 800
#define MAX_DECREASE 4
#define ZOOPS_APPROX_CACHE 10000000
#define MAX_GLOBAL_AA 10000000

#define MAX_NUM_SEQUENCES 10000000

#define MIN_MATRIX_PROB 0.001
#define MOTIF_INIT_VALUE 0.8
#define SKIPMER_PERCENT_SKIP 0.1

#define MAX_STAT_CACHE 100000000

#define MOTIF_HASH_SIZE 10000000
#define BUFFER 1000000
#define BSIZE 1000000

#define TMP_MOTIF_FILE ".tmp.motifs"
#define MER_MOTIF_FILE ".mer.motifs"
#define WEIGHT_RES 1.0

#define MAX_MISMATCHMERS 1000000
#define DEG_PVALUE_THRESH -1.386

//might want to recompile with short for datasets with less than 60,000
//that are having trouble fitting in memory
//#define gene_id_t unsigned short
//#define EMPTY_GENE 65535
#define gene_id_t unsigned int
#define EMPTY_GENE 1194967295

#define FILL_FRACTION 0.50

class ConStat;
class ConStatExact;
class ConStatApprox;
class SeqTree;

class PSSM {
public:
	float** matrix;
	int size;
	char* name;
	char* outputString;
	int length;
	float v;
	float threshold;
	float pvalue;
	float state;
	float score;
	float *backFreq;
	int type;
	int constatID;
	int numGapLengths;
	int* gapLengths;
	PSSM();
	PSSM(PSSM*);
	void setSize(int seqType);
	PSSM(char* motif, int seqtype);
	PSSM(int* values, int seqtype);
	PSSM(char* newname,float**,int, int seqType);
	~PSSM();
	void init();
	void set(char*,int seqType);
	void set(int*,int seqType);
	void setBackFreq(float*);
	void print(FILE*);
	void printMatrix(FILE*);
	void logXform();
	void expXform();
	void setGapInfo(int* newGapLengths, int newNumGaps);
	void setOutputString(char* s);
	PSSM* copy();
	void adjustvalues();
	void adjustPrint();
	void mildize(float factor);
	void setName(char*);
	float optThreshold(float fdr,float* background);
	char* consensus();
};



class Site {
public:
	char* seq;
	char* geneID;
	int dir;
	int gapsize;
	int gapoffset;
	PSSM* motif;
	float conservation;
	float ratio;
	int offset;
	float score;

	Site();
	Site(char* geneid, char* nseq, PSSM* nmotif, int ndir, int noffset,
                    int ngapsize, float ncons, float nratio, float nscore);
	void setSeq(char* nseq);

	~Site();
	void init();
	void print(FILE*);

};


class XMerData {
public:
	int numXMers;
	ConStat** stats;
	char** mers;

	int numGenes;
	int numPosGenes;
	int numNegGenes;
	float enumGenes; //effective number of genes (after weighting)
	float enumPosGenes;  //effective number of genes (after weighting)
	float enumNegGenes;
	float enumSites; //effective number of genes (after weighting)
	float enumPosSites;  //effective number of genes (after weighting)
	float enumNegSites;
	float* papprox;
	float* napprox;
	int* geneStats;
	float* geneWeights;
	char** geneNames;
	SeqTree * seqTree;

	XMerData();
	~XMerData();

	void init();
	void loadMerData(Hashtable* m);
	void loadGeneData(char** names, float* weights, int* gstats, int ngenes);
	void autoScale();
	void initZOOPSCache();

};

class ConStat {
public: 
	int id;
	float pvalue;
	float dpvalue; // degenerate pvalue
	unsigned char dir; //hold information about alignment to motif

	ConStat();
};

class ConStatExact : public ConStat {
public:
    int p;
    int n;
	gene_id_t * pgenes;
	gene_id_t * ngenes;
	ConStatExact();
	~ConStatExact();
	void addGene(gene_id_t g,int stat);
	void addRank(int r);
	float fp();
	float fn();
	int sum();
};

class ConStatApprox : public ConStat {
public:
    float fp;
    float fn;
	ConStatApprox();
};


class SeqTree {

public:
	SeqTree** branch; //recusive nodes in the sequence tree
	int id;  //hold xmer index represented by a leaf
	char size; // size = 4 for DNA, 20 for AA - specifices number of branches

	//IUPAC and IUPACIndex are for finding all the IUPAC sequences in the tree
	int IUPAC; // current IUPAC symbol
	int IUPACIndex; // current branch for the current IUPAC symbol

	SeqTree(char newSize);
	~SeqTree();
	void insert(char* mer, int clevel, int  nid);
	void findMisMatchMers(int* mer, int clevel, int mis, int* ids, 
				int &numIDs);
	int getNextIUPAC(int* value, int clevel, int numIUPAC, int IUPACtype);
};

class ProteinStat {
public:
	int* v;
	float pvalue;
	ProteinStat();
	~ProteinStat();
};


class IndexTree {
public:
	int* values;
	IndexTree** tree;
	int size;
	int recallIndex;
	IndexTree();
	~IndexTree();
	int insert(int* mer, int clevel);
	int getNextSeq(int* value, int celvel);
};

class Sequence {
public:
	int stat;
	char* seq;
	float w;
	char* name;
	Sequence();
	~Sequence();
	void setSeq(char*);
	void setName(char*);
};


void split(char* string, char** cols, int &numCols, char delim);

Sequence* loadSequences(CommandLine *cmd, int &numSeqs);
void findMotifSites(CommandLine* cmd, Sequence* seqs, int numSeqs, PSSM* motif, FILE* outfile);

XMerData* indexGeneSeq(CommandLine *cmd);
void writeMerFile(CommandLine* cmd, XMerData* mers, FILE* fp);

void scoreDegenerateMers(CommandLine* cmd, XMerData* xmer,int* merIDs, int numMerIDs);
void findMotifs(CommandLine* cmd, XMerData* xmer);

int calculateMerPSSMscore(CommandLine* cmd, PSSM* pssm, XMerData* xmer,
                ConStat** mers, int numMers, int updateFlag);
int localOptimizeMerPSSM(CommandLine* cmd, PSSM* pssm, XMerData* xmer,
                ConStat** mers, int numMers);
PSSM** findMotifsOfLength(CommandLine* cmd, XMerData* xmer, ConStat** mers, int numMers,
                    int* seeds, int numSeeds, int length, int &numpssms);
float calculateProbabilityPSSM(CommandLine* cmd, PSSM* pssm,char* mer,unsigned char &dir,
            int allowRevOpp);

char** getMersOfLength(CommandLine* cmd, Hashtable* mers, int length, int &numMers);
int checkforSimilarSequence(CommandLine* cmd, char* seq,XMerData*, Inttable* seq2skip);
int checkforSimilarSequence(CommandLine* cmd, int* seq,XMerData*, Inttable* seq2skip);
char** getSimilarMers(CommandLine* cmd, char* mer, int &numMers);
Hashtable* readMerFile(CommandLine* cmd, int &numGenes, int &numPos);
char* seqIndexOf(int seqType,char* mer, char* values);
int* seqIndexOfIUPAC(int seqType,char* mer, int* values);
char* getStringFromIUPACIndex(int* mer);
float alignPSSMs(PSSM* a, PSSM* b, int &bestOffset, int &bestDir);

char** findBestSeeds(CommandLine* cmd, XMerData* xmer, int *merIDs, int numMerIDs, int &numSeeds);


XMerData* readMerFile(CommandLine* cmd);
char** mismatch3Seq(CommandLine* cmd, char* og,int &numKeys);
char** mismatch2Seq(CommandLine* cmd, char* og,int &numKeys);
char** mismatch1Seq(CommandLine* cmd, char* og,int &numKeys);
PSSM** extractPSSM(Floattable* hash, int &numPSSMs,int printFlag);
char* revopp(char* motif);
void readMatrix(CommandLine* cmd, PSSM** &allMatrix, int &numMatrix);
char* getAAStringfromIndex(int* values);
						
void indexSeq(char* seqStr, char* &seq, char* &rseq, int seqType);
void indexDNA(char* seqStr, char* &seq, char* &rseq);
void indexProtein(char* seqStr, char* &seq, char* &rseq);
void eraseSequence(char* seq, int offset, int length, int seqType);
char* getRepDualMotif(CommandLine* cmd, char* m);
char** dualmismatch1Seq(CommandLine* cmd,char* og,int &numKeys);
Site* scoreSeqWithMatrix(CommandLine* cmd, PSSM* motif, char* id, char* seq, 
                                int seqLen,Site** &sites, int &numsites);
Site* scoreMatrix(CommandLine* cmd, PSSM* motif, char* id, char* seqStr, int eraseFlag,
                                Site** &sites, int &numsites);



#endif
