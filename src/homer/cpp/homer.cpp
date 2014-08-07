
// Copyright 2009, 2010, 2011, 2012 Christopher Benner <cbenner@ucsd.edu>
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
#include "Motif.h"
#include "CommandLine.h"

int main (int argc, char** argv) {
	srand(time(NULL));

	CommandLine* cmd = new CommandLine(argc, argv);
	cmd->parseInput();
	cmd->checkDependencies();

	XMerData* xmer = NULL;
	PSSM** pssm = NULL;
	int numpssms = 0;
	int numSeqs = 0;
	Sequence* seqs = NULL;

	//fprintf(stderr, "Action = %d\n", cmd->action);

	if (cmd->action == ACTION_MERS) {

		xmer = indexGeneSeq(cmd);
		writeMerFile(cmd, xmer,stdout);
		
	} else if (cmd->action == ACTION_DMERS) {

		if (cmd->merfile == NULL) {
			//need to generate mers from sequence file
			xmer = indexGeneSeq(cmd);
		} else {
			//xmer = readMerFile(cmd);
			cmd->exactFlag = 0; //cant have exact flag if using merfile
		}
		//scoreDegenerateMers(cmd, xmer,NULL,0);
		writeMerFile(cmd,xmer,stdout);

	} else if (cmd->action == ACTION_SORTMERS) {
	
		xmer = readMerFile(cmd);
		writeMerFile(cmd, xmer,stdout);

	} else if (cmd->action == ACTION_FIND || cmd->action == ACTION_GENESCORE) {
		readMatrix(cmd, pssm, numpssms);
		seqs = loadSequences(cmd, numSeqs);
		for (int i=0;i<numpssms;i++) {
			pssm[i]->print(stderr);
			findMotifSites(cmd,seqs, numSeqs,pssm[i],stdout);
		}
	} else if (cmd->action==ACTION_REMOVE) {
		readMatrix(cmd, pssm, numpssms);
		char tempfile[100] = ".tmp.seq";
		char tempfile2[100] = ".tmp.seq2";
		char c[1000];
		sprintf(c,"cp %s %s", cmd->seqfile, tempfile);
		(void)system(c);
		cmd->seqfile = tempfile;
		
		for (int i=0;i<numpssms;i++) {
			FILE *fp = fopen(tempfile2, "w");
			seqs = loadSequences(cmd, numSeqs);
			findMotifSites(cmd,seqs, numSeqs, pssm[i],fp);
			fclose(fp);
			delete []seqs;
			numSeqs = 0;
			sprintf(c, "mv %s %s",tempfile2, tempfile);
			(void)system(c);
		}
		FILE *FP = fopen(tempfile, "r");
		char buf[BUFFER];
		while (fgets(buf,BUFFER,FP) != NULL) {
			printf("%s", buf);
		}
		fclose(FP);

	} else if (cmd->action == ACTION_OPTPVALUE) {
		
		readMatrix(cmd, pssm, numpssms);
		seqs = loadSequences(cmd, numSeqs);
		for (int i=0;i<numpssms;i++) {
			findMotifSites(cmd,seqs,numSeqs,pssm[i],NULL);
			pssm[i]->expXform();
			pssm[i]->print(stdout);
		}

	} else if (cmd->action == ACTION_GETPVALUE) {
		
		readMatrix(cmd, pssm, numpssms);
		seqs = loadSequences(cmd, numSeqs);
		for (int i=0;i<numpssms;i++) {
			findMotifSites(cmd,seqs,numSeqs,pssm[i],NULL);
			pssm[i]->expXform();
			pssm[i]->print(stdout);
		}

	} else if (cmd->action == ACTION_CLUSTER) {

		xmer = readMerFile(cmd);
		findMotifs(cmd, xmer);

	} else if (cmd->action == ACTION_MOTIFS) {
		
		if (cmd->seqfile != NULL) {
			xmer = indexGeneSeq(cmd);
		} else if (cmd->merfile != NULL) {
			//mers = readMerFile(cmd, numGenes, numPosGenes);
		}
		findMotifs(cmd, xmer);

	} else if (cmd->action == ACTION_REFINE || cmd->action == ACTION_REFINETHRESH) {
		
		readMatrix(cmd, pssm, numpssms);
		//int minlength = 1000;
		//int maxlength = 0;
		cmd->numMerLengths = 0;
		cmd->merLengths = new int[numpssms];
		fprintf(stderr, "Mer lengths of ");
		for (int i=0;i<numpssms;i++) {
			int bad = 0;
			for (int j=0;j<cmd->numMerLengths;j++){ 
				if (pssm[i]->length == cmd->merLengths[j]) {
					bad = 1;
					break;
				}
			}
			if (!bad) {
				cmd->merLengths[cmd->numMerLengths++] = pssm[i]->length;
				fprintf(stderr, " %d", pssm[i]->length);
			}
		}
		fprintf(stderr, "\n");
		
		if (cmd->seqfile != NULL) {
			xmer = indexGeneSeq(cmd);
		} else if (cmd->merfile != NULL) {
			//mers = readMerFile(cmd, numGenes, numPosGenes);
		}
		ConStat** cs = new ConStat*[xmer->numXMers];
		for (int i=0;i<numpssms;i++) {
			int plen = pssm[i]->length;
			int numMers = 0;
			for (int j=0;j<xmer->numXMers;j++) {
				if (strlen(xmer->mers[j]) == (unsigned int) plen) {
					cs[numMers++] = xmer->stats[j];
				}
			}

			//int curIndex = 0;
			if (cmd->action == ACTION_REFINE) {
				(void)localOptimizeMerPSSM(cmd,pssm[i],xmer,cs,numMers);
			} else if (cmd->action == ACTION_REFINETHRESH)  {
				(void)calculateMerPSSMscore(cmd,pssm[i],xmer,cs,numMers,1);
			}
				
			pssm[i]->expXform();
			pssm[i]->print(stdout);
		}
		
	}
}
