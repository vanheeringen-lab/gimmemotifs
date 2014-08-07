
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
#include <strings.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "Hashtable.h"
#include "statistics.h"

#ifndef COMMANDLINE_H
#define COMMANDLINE_H

#define ACTION_MOTIFS 0
#define ACTION_MERS 1
#define ACTION_DMERS 2
#define ACTION_FIND 3
#define ACTION_OPTPVALUE 4
#define ACTION_REFINE 5
#define ACTION_CLUSTER 6
#define ACTION_SORTMERS 7
#define ACTION_REMOVE 8
#define ACTION_GENESCORE 9
#define ACTION_REFINETHRESH 10
#define ACTION_GETPVALUE 11
#define NO_ACTION -1

#define ALG_HYPERGEO 0 
#define ALG_BINOMIAL 1
#define ALG_APPROXBINOMIAL 6
#define ALG_SITEHYPERGEO 2
#define ALG_SITEBINOMIAL 3
#define ALG_FISHER 4
#define ALG_RANK 5
#define ALG_FREQDIFF 7 

#define IUPAC_N 1
#define IUPAC_2N 2
#define IUPAC_FULL 3

#define DNA_SEQ_TYPE 0
#define PROTEIN_SEQ_TYPE 1

#define SPEED_NORMAL 0
#define SPEED_FAST 1


class CommandLine {
public:
	int argc;
	char** argv;
		
	int algorithm;
	char* seqfile;
	char* consfile;
	char* statfile;
	char* merfile;
	char* motiffile;
	char* seedfile;
	char* outputfile;
	int* merLengths;
	int numMerLengths;
	int revoppFlag;
	int seqType;
	float consThresh;
	int offset;
	float nratiocutoff;
	int rmalignseeds;
	int numPreCheckSeeds;
	int minmer;
	int maxmer;
	int mismatches;
	int numIUPAC;
	int IUPACtype;
	int numSeeds;
	int exactFlag;
	int weightFlag;
	int maxOptIterations;
	int action;
	int minSeqLen;
	int maxSeqLen;
	int autoScale;
	int gapsize;
	int allowFlip;
	int targetTreeFlag;
	int gapoffset;
	int eraseFlag;
	int fisherSize;
	int zoopsApprox;
	int dualMotifs;
	int speed;
	int* gapLengths;
	int numGapLengths;
	int maxPerSeq;
	int freqAdjust;
	float *backFreq;
	float branchSize;
	StatMemory* statmemory;

	//Glbal Variables - fudge factors
	Inttable *seqLength;
	float maxNegGenePercentage;
	int totalSeqLength;
	int numSites;
	int numPosSites;
	
	CommandLine(int argc,char** argv);
	~CommandLine();
	void parseInput();
	void checkDependencies();
	void printHelp();

	
};


int * parseRangeVariable(char* input, int &numLengths);
char* makeVariableRangeString(int* lengths, int numLengths);

#endif

