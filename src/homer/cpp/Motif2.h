
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
#include "statistics.h"
#include <limits.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>

#ifndef MOTIF_H
#define MOTIF_H

#define MOTIF2_BUFFER 1000000

#define MOTIF_STRAND_POS 0
#define MOTIF_STRAND_NEG 1
#define MOTIF_STRAND_BOTH 2
#define MOTIF_STRAND_SEPARATE 3

#define MOTIF_SCORING_HYPERGEOMETRIC 0
#define MOTIF_SCORING_BINOMIAL 1
#define MOTIF_SCORING_ZSCORE 2
#define MOTIF_MISMATCH_PVALUE_MINIMUM -1.3863
//-0.6931 = 0.50
//-1.3863 = 0.25
//-2.3026 = 0.10
//-2.9957 = 0.05
//-4.6052 = 0.01
//-6.9078 = 0.001

#define MOTIF_MAX_BACKGROUND_PERCENT 0.5

#define MOTIF_SCORE_OLIGOS 0
#define MOTIF_REMOVE_OLIGOS 1

#define MOTIF_MISMATCH_MINIMUM_SKIP 1
#define MOTIF_MISMATCH_MINIMUM_COUNTS 1
#define MOTIF_MAX_SEED_OLIGOS_FRACTION 0.10
#define MOTIF_MAX_SEED_OLIGOS_TOTAL 1e5

#define MAX_OLIGO_NORMALIZATION_ITERATIONS 160
#define MIN_OLIGO_NORMALIZATION_IMPROVEMENT 0.01

#define MOTIF_SEED_MINIMUM_LOGP -10.0

#define OLIGO_SEARCH_KEEP 0
#define OLIGO_SEARCH_REMOVE 1

#define MOTIF_DEFAULT_NUMBER_OF_MOTIFS 25

#define MOTIF_INITIAL_DEGENERACY 0.1
#define MOTIF_PROBABILITY_MINIMUM 0.001
#define MOTIF_IUPAC_THRESHOLD 0.15

#define MOTIFS_NUM_POTENTIAL_MOTIFS 10
#define MOTIF_OPTIMZATION_DIFFERENT_THRESHOLDS 0
#define MOTIFS_MAX_OPT_ITERATIONS 10
#define MOTIF_MAX_HIT_RATE_TO_CONSIDER 0.005

#define MAXIMUM_MOTIF_LENGTH 10000
#define MAXIMUM_DENOVO_MOTIF_LENGTH 50

#define DENOVO_SPEED_AUTOADJUST 0
#define DENOVO_SPEED_CONSTANT 1
#define DENOVO_MASK_FULL 0
#define DENOVO_MASK_QUICK 1

#define SCORER_STANDARD_ZOOPS 0
#define SCORER_APPROXIMATE_ZOOPS 1

#define KNOWN_JUST_SCORE_ZOOPS 0
#define KNOWN_OPTIMIZE_THREASHOLD 1
#define KNOWN_MASK_MOTIFS 2
#define KNOWN_MSCORE_MOTIFS 3
#define KNOWN_FIND_ALL 4

#define MOTIF_OFFSET_NULL -12345678

#define MOTIF_OVERLAP_PERCENT 0.75

#define MOTIF_DEFAULT_REDUCE_PERCENTAGE 0.20

#define HOMER_MINIMUM_SEQ_WEIGHT 0.001

#define DNA_ALPHA 4

#define DOUBLE_EPS 1e-20
#define FLOAT_EPS 1e-10

#define MAX_STAT_CACHE_SIZE 500

#define HOMER_VERBOSE_LEVEL0 0
#define HOMER_VERBOSE_LEVEL1 1

class DeNovoMotifs;
class KnownMotifs;
class OligoArray;
class Oligo;
class Motif;
class EnrichmentScorer;
struct ThreadArgs_findMisMatchEnrichmentThread;
struct ThreadArgs_optimizeMotifThresholdThread;
struct ThreadArgs_scanMotifThread;
class SequenceStats;

class OligoArray {
public:

	int* tree; //if tree instances are 0 or negative, they reference oligos, 1=tree root
	int numTree;
	int maxTree;
	int alphaSize;
	int oligoLength;

	char* alpha;
	char* rvAlpha;
	char* rvMap;
	char* alphaIndex;
	char* rvAlphaIndex;
	
	double totalTarget;
	double totalBackground;

	Oligo** staticOligos;
	int numOligos;

	Oligo** activeOligos;
	int numActiveOligos;

	int maxOligos;
	int rvFlag;
	int maxAutoNormIters;

	pthread_mutex_t sortMutex;

	//stack variables for "recursive" looping - save's on new and delete time to keep them here

	int numCPUs;

	Oligo*** utilityOligos;
	int* numUtilityOligos;

	int** indexStack;
	int** seqIndex;
	int** alphaStack;
	int** positionStack;
	int** misMatchStack;
	int** fwdMisMatchStack;
	int** rvMisMatchStack;
	double** fwdScoreStack;
	double** rvScoreStack;
	


	OligoArray();
	~OligoArray();
	void init(int length, int totalDataSetBp, int alpha, int strandFlag,int numCPUs);
	void optimizeMotifThreshold(Motif* motif);
	void optimizeMemoryUsage();
	void printOligo(char* seq, FILE* fp);
	int addOligo(char* seq, double vt, double vbg, double addValue); //addValue = -1 for normal add
	int searchOligo(char* seq, int maxMisMatches, double& vt, double& vbg,int removeFlag,int cpu);
	void sortActiveOligos();
	void sortUtilityOligos(int cpu);
	void initializeActiveOligoArray();
	void initializeOligoPvalues(EnrichmentScorer* scorer);
	void normalizeOligos(int nlen, char* normFile);
	double normalizeOligosIteration(int nlen, char* normFile,int printFlag);
	void scoreOligosWithMotif(Motif* motif,double minimumScore, int removeFlag,int cpu);
	void calculateFoldEnrichment(double pseudoCount);
	void adjustOligoInstances(char* oligoSeq, double targetValue, double backgroundValue);
	void readOligoCountFile(char* file, int strand, int cpus);

private:
	int addOligoTree(char* seq, double vt, double vbg, double addvalue,int rvMode, int knownID);
	int searchOligoTree(char* seq,int maxMisMatches, double& vt, double& vbg);

};


class Oligo {
public:
	char* seq;
	double numTarget;
	double numBackground;
	float normWeight; // used if -olen is used to normalize short oligos
	float logp;  // used to figure out contribution to motifs when optimizing matrix
	double value; //generic value - used for different things
					// pvalue during findMisMatchEnrichment
	char flag; //flag used to make sure we don't doulbe count oligos when traversing tree
	char flag2; //generic flag - used to mark oligos to skip during mismatch analysis
	double* cpuData;
	Oligo();
	Oligo(char* seq);
	~Oligo();
	void init();
	void print(FILE*);
	static int oligoCmp(const void* a, const void* b);
	static int oligoCmpReverse(const void* a, const void* b);
	static int cpu;
	static int oligoCmpReverseParallel(const void* a, const void* b);
	static void revopp(char* seq);
};

class Sequence {
public:
	char* s;
	int g; // group, i.e. target (=1) or background (=0)
	double w; // i.e. 1.0
	char* name;
	int length;
	int id;

	Sequence();
	Sequence(char* name, char* seq, int group, double weight);
	~Sequence();
	void init();
	void setSeq(char* name, char* seq, int group, double weight);
	int motifScan(Motif* m,int strand, int offset,int mode, double &bestScore, 
					int* sitePositions, char* siteStrands, FILE* fp, pthread_mutex_t fileLock);
};
class LinkedSequence {
public:
	Sequence* s;
	LinkedSequence* link;
	
	LinkedSequence(Sequence* seq, LinkedSequence* link);
	~LinkedSequence();
};

class SequenceArray {
public:
	Sequence** seqs;
	int numSeqs;
	LinkedSequence* link;
	double numTarget;
	double numBackground;
	int bpTarget;
	int bpBackground;
	int longestSequence;
	EnrichmentScorer* scorer;	
	double maxBackgroundPercent;
	int maxAutoNormIters;

	SequenceArray();
	~SequenceArray();
	void deleteLinks();
	void addSequence(Sequence* s);
	void optimizeArray();
	void readSequenceFiles(char* seqFile, char* groupFile);
	void printGroupFile(FILE* fp);
	void printSequenceFile(FILE* fp);
	void motifScan(Motif* m,int strand, int offset,int mode,FILE* fp, pthread_mutex_t fileLock);
	void normalizeSequence(int nlen, int strand, char* normFile, int equalFlag);
	double normalizeSequenceIteration(int nlen, int strand, char* normFile, int printFlag,int neutralFlag);
	void maskMotifFromOligoArray(Motif* motif,OligoArray* oligoArray,char strand);
	void initializeScorer(int scoringMethod,int zoopsFlag,unsigned long long int maxCacheSize);
	void optimizeMotifThreshold(Motif* motif, char strand);
	void parseFasta2SeqAndGroupFiles(char* inputFASTAfile,char* bgFASTAfile,char* sFile, char* gFile);

};


class Motif {
public:
	double** m;
	int length;
	double threshold;
	int logFlag;
	char* name;
	char* consensus;
	double numTarget;
	double numBackground;
	double freqTarget;
	double freqBackground;
	double numSeqTarget;
	double numSeqBackground;
	double logp;
	double minimumScore;
	int badFlag;

	double similarPercent;
	double similarPercentOfOther;
	Motif* similarMotif;

	SequenceStats* seqStats;
	int numSeqs;

	double statAvgTarget;
	double statAvgBackground;
	double statStdTarget;
	double statStdBackground;
	double statStrandLogRatio;
	double siteMultiplicity;


	double* fwdMaxPossible;
	double* rvMaxPossible;

	static int alphaSize;
	static char* alphaIndex;

	Motif();
	Motif(char* seq);
	Motif(Motif*);
	Motif(char* name, char* consensus, double threshold, double logp, double** matrix, int mlen);
	~Motif();
	void init();
	void setName(char*);
	void logxform();
	void expxform();
	void normalize();
	void justNormalize();
	void print(FILE*);
	void getConsensus();
	char getSymbol(char*);
	void setMinimumScore(double percent);
	void setBlank(int len);
	void addOligo(Oligo* o);
	void initializeSeqStats(int);
	void optimizeThreshold(EnrichmentScorer* scorer, double maxBackPercent);
	void adjustLength(int);
	void calculateDistributionStats();

	static double scoreOverlap(Motif* m1, Motif* m2, double &percent1, double &percent2);
	static int motifCmp(const void* a, const void* b);
	static int nameID;
	static Motif** readMotifFile(char* file, int &numMotifs);
};


class DeNovoMotifs {
public:
	OligoArray* oligoArray;
	unsigned long long int totalOligos;
	unsigned long long int numOligoEstimate;
	unsigned int oligoLength;
	unsigned int normLength;
	int bufferLength;
	unsigned char alphaSize;
	unsigned char strand;
	double maximumExpectedPerBp;
	unsigned long long int maxCacheSize;
	int noDuplicateOligos;

	unsigned int numCPUs;
	int mutexIndex;
	int mutexTotal;
	int nextMeter;
	int incMeter;
	pthread_mutex_t mutex;

	double numTargets;
	double numBackground;
	Oligo** seedOligos;
	int numSeedOligos;

	Motif** motifs;
	int numMotifs;

	SequenceArray* fullSequences;

	char* seqFile;
	char* groupFile;
	char* inputFASTAfile;
	char* bgFASTAfile;
	int fastaFlag;
	char* normFile;
	char* tmpFile;
	FILE* sitefp;

	int randFlag;

	int speedFlag;
	double siteReduceThreshold;

	int optimizeGivenFlag;
	int scoringMethod;
	int finalZoopsFlag;
	double maxBackgroundPercent;

	int maxMisMatches;
	int misMatchTargetOnlyFlag;
	double misMatchMinPvalue;
	int minimumMisMatchToStartSkipping;
	double minimumCountsToStartSkipping;
	double minimumFractionToStartSkipping;

	double minimumSeedLogp;
	int totalMotifsToFind;
	int numTrialMotifs;
	int numConcurrentMotifs;
	int maxOptimizationIterations;
	int motifOptimizationMethod;
	int maskMethod;
	int maxAutoNormIters;

	int zoopsFlag;
	EnrichmentScorer* scorer;	
	//double (DeNovoMotifs::*enrichmentFunction)(double pv, double nv);

	DeNovoMotifs();
	~DeNovoMotifs();
	
	void init();
	void readSequence2Tree();
	void loadFullSequences();
	void checkParameters();
	void getOligoEnrichment();
	void normalizeOligos();


	void findMisMatchEnrichment(int maxMisMatches,FILE* fp,int verboseLevel);
	void findMisMatchEnrichmentThread(ThreadArgs_findMisMatchEnrichmentThread* threadArgs);
	void optimizeMotifMatrix(Motif* &motif);
	void optimizeMotifThreshold(Motif* m, int cpu);
	void optimizeMotifThresholdThread(ThreadArgs_optimizeMotifThresholdThread*);
	void initializeEnrichmentScoring();
	double scoreEnrichmentHypergeometric(double pv, double nv);
	double scoreEnrichmentBinomial(double pv, double nv);
	void getSeedOligos();
	void optimizeSeeds2Motifs();
	void optimizeGivenMotifs();
	void scoreWithKnownMotifs(int optFlag);
	void freeOligos();
	void parseFasta();
	void printMotifs(FILE* fp); // sorts motifs array
};

class EnrichmentScorer {
public:

	unsigned int numTargets;
	unsigned int numBackground;

	double** cache;
	int cacheLength;
	unsigned long long int maxCacheSize;

	int zoopsFlag;
	double* targetZoopsApprox;
	double* bgZoopsApprox;
	unsigned int maxZoopsApprox;

	double (EnrichmentScorer::*enrichmentFunction)(double pv, double nv);

	EnrichmentScorer();
	~EnrichmentScorer();
	double scoreEnrichment(double pv, double nv);
	void initializeEnrichmentScoring(int method, double numTarget, double numBack, int zoopsFlag, 
									unsigned long long int newMaxCacheSize);
	double scoreEnrichmentHypergeometric(double pv, double nv);
	double scoreEnrichmentBinomial(double pv, double nv);
	double scoreEnrichmentZscore(double pv, double nv);
	double getZoopsApproxTarget(double value);
	double getZoopsApproxBackground(double value);
};


class KnownMotifs {
public:
	SequenceArray* seqArray;
	unsigned long long int totalBp;
	unsigned int normLength;
	unsigned char alphaSize;
	unsigned char strand;
	unsigned long long int maxCacheSize;
	int offset;

	double numTargets;
	double numBackground;

	Motif** motifs;
	int numMotifs;

	unsigned int numCPUs;
	int mutexIndex;
	int mutexTotal;
	double nextMeter;
	double incMeter;
	pthread_mutex_t mutex;
	pthread_mutex_t fileLock;

	char* seqFile;
	char* groupFile;
	char* inputFASTAfile;
	char* bgFASTAfile;
	int fastaFlag;
	char* knownFile;
	char* normFile;
	FILE* sitefp;
	FILE* normfp;
	char* motifFilename;

	double maxBackgroundPercent;
	double siteReduceThreshold;

	int scoringMethod;
	EnrichmentScorer* scorer;	

	int mscoreFlag;
	int threasholdOptimizationFlag;
	int zoopsFlag;

	int knownFlag;
	int findFlag;
	int maskFlag;
	int normOnlyFlag;
	int maxAutoNormIters;
	int neutralFlag;

	KnownMotifs();
	~KnownMotifs();
	void checkParameters();
	void loadMotifs();
	void loadSequence();
	void parseFasta();
	void getMotifEnrichment();
	void printMotifEnrichment(FILE* fp);
	void printMotifs(FILE* fp);
	void scanMotifThread(ThreadArgs_scanMotifThread* args);
	void findOverlappingMotifs();
	
};

class GroupInfo {
public:
	double weight;
	unsigned int group;
};

class SequenceStats {
public:
	Sequence* seq;
	double score;
	int numSites;
	int* pos;
	char* strands;
	SequenceStats();
	~SequenceStats();
	void initialize();
	void setSites(int* positions, char* strands, int num);
	static int scoreCmp(const void* a, const void* b);
	static int defaultOrder(const void* a, const void* b);
};

void splitMotif2(char* string, char** cols, int &numCols, char delim);
void cleanSequence(char* seq);
int checkSequence(char* seq);
int sortDouble(const void* a, const void* b);

void printCMDhomer();
void printCMDdenovo();
void programDeNovo(int argc, char** argv);
void programKnown(int argc, char** argv);


struct ThreadArgs_findMisMatchEnrichmentThread {
public:
	DeNovoMotifs* denovo;
	int cpu;
	int maxMisMatches;
	int totalChecked;
	int numSkipped;
	int dueToTargetTags;
	int dueToLessThanLast;
	int dueToLowPvalue;
	int lastChecked;
	int lastCheckedTotal;
	int verboseLevel;
	ThreadArgs_findMisMatchEnrichmentThread(DeNovoMotifs* obj,int cpu,int mis, int total,int vlevel);
};
void* DeNovoMotifs_findMisMatchEnrichmentThread(void* threadArgs);

struct ThreadArgs_optimizeMotifThresholdThread {
public:
	DeNovoMotifs* denovo;
	int cpu;
	Motif** motifs;
	ThreadArgs_optimizeMotifThresholdThread(DeNovoMotifs* obj,int cpu,Motif** motifs);
};
void* DeNovoMotifs_optimizeMotifThresholdThread(void* threadArgs);

struct ThreadArgs_scanMotifThread {
public:
	KnownMotifs* known;
	int cpu;
	ThreadArgs_scanMotifThread(KnownMotifs* obj,int cpu);
};
void* KnownMotifs_scanMotifThread(void* threadArgs);



#endif
