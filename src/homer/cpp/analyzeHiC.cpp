

// Copyright 2009 - 2014 Christopher Benner <cbenner@salk.edu>
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

#include "SeqTag.h"

void printCMD();
void printCircosFiles(char* prefix, int resolution, TagLibrary* tags, PeakLibrary* peaks1,
				PeakLibrary* peaks2,int,int,int,int,char**,int,char**,int,char**,int,TagLibrary*);

int main(int argc, char** argv) {

	int maxFiles = 1000;
	char* directory = NULL;
	char* bgDirectory = NULL;
	char* customBgFile = NULL;
	char* customBgBgFile = NULL;
	int defaultResolution = 1;
	int resolution = 10000000;
	int fullModelFlag = 1;
	char* chr1 = NULL;
	char* chr2 = NULL;
	int start1 = -1;
	int start2 = -1;
	int override = 0;
	int relativeFlag=0;
	
	int washuFlag = 0;
	int superRes = -1;
	int numCPUs = 1;
	int vsGenomeFlag = 0;
	//int boundaryFlag = 0;
	int boundaryScale = 0;
	int end1 = 2000000000;
	int end2 = 2000000000;
	int minDist = -1;
	int maxDist = -1;
	int nologFlag = 0;
	int chopifyFlag = 0;
	int fixedFlag = 0;
	int forceFlag = 0;
	int bgOnlyFlag = 0;
	int randomReadTotal = 100000000;
	char* randomizeModel = NULL;
	int siteFlag = 0;
	char* restrictionSite = NULL;
	int maxMisMatches = 0;
	char* poutFile = NULL;
	char* pout2File = NULL;
	char** geneFiles = new char*[maxFiles];
	char** bedFiles = new char*[maxFiles];
	char** tagDirs = new char*[maxFiles];
	for (int i=0;i<maxFiles;i++) {
		geneFiles[i] = NULL;
		tagDirs[i] = NULL;
		bedFiles[i] = NULL;
	}
	int numTagDirs = 0;
	int numBedFiles = 0;
	int numGeneFiles = 0;
	
	char* outputFile = NULL;
	int histFlag = 0;
	int size = 0;
	int centerFlag = 0;
	int actionFlag = HIC_MASK_NORM_DISTANCE | HIC_MASK_NORM_SEQDEPTH;
	char* fname = new char[10000];

	char* interactionFile = NULL;
	char* circosPrefix = NULL;
	char* lociFile = NULL;
	double pvalueThreshold = 1e-3;
	double zscoreThreshold = 1.0;
	//double clusterThreshold = 0.5;
	int clusteringIterations = 1;
	double stdFilter = 4.0;
	double minFilter = 0.20;
	double minExpReads = 0.0;

	char* peakFile1 = NULL;
	char* peakFile2 = NULL;
	int mask = 0;

	float maxTbp = 1.0;
	int maxCircosInteractions = 5000;
	
	char* inputInteractions  = NULL;
	char* rawbedFile = NULL;
	char* fourCbedFile = NULL;
	char* rawbedFilePrefix = NULL;
	char* expectMatrixOutputFile = NULL;
	int expectMatrixFlag=0;

	HiCparams* params = new HiCparams();


	//HomerConfig*  hc = new HomerConfig();

	if (argc < 2) {
		printCMD();
	}
	//strcpy(cmd,argv[0]);
	//for (int i=1;i<argc;i++) {
	//	strcat(cmd," ");
	//	strcat(cmd,argv[i]);
	//}

	for (int i=1;i<argc;i++) {
		
		if (i==1) {
			if (argv[i][0] == '-') {
				fprintf(stderr, "\n!!!!!!!!!!!!\n\tNEED to specify paired-end tag directory with first argument!!!\n");
				printCMD();
			}
			directory = argv[i];
			continue;
		}
		if (argv[i][0] == '-') {
			if (strcmp(argv[i],"-chr")==0) {
				chr1 = argv[++i];
				if (chr2 == NULL) chr2 = chr1;
			} else if (strcmp(argv[i],"-pos")==0) {
				parseUCSCpositionStr(argv[++i],chr1,start1,end1);
				if (chr2 == NULL) chr2 = chr1;
				if (start2 == -1) start2 = start1;
				if (end2 == 2000000000) end2 = end1;
			} else if (strcmp(argv[i],"-vsGenome")==0) {
				vsGenomeFlag = 1;
			} else if (strcmp(argv[i],"-pos2")==0) {
				parseUCSCpositionStr(argv[++i],chr2,start2,end2);
			} else if (strcmp(argv[i],"-center")==0) {
				centerFlag = 1;
			} else if (strcmp(argv[i],"-boundary")==0) {
				actionFlag = actionFlag | HIC_MASK_BOUNDARIES;
				//boundaryFlag = 1;
				sscanf(argv[++i],"%d",&boundaryScale);
			} else if (strcmp(argv[i],"-pout")==0) {
				poutFile = argv[++i];
			} else if (strcmp(argv[i],"-pout2")==0) {
				pout2File = argv[++i];
			} else if (strcmp(argv[i],"-randomize")==0) {
				randomizeModel = argv[++i];
				sscanf(argv[++i],"%d",&randomReadTotal);
			} else if (strcmp(argv[i],"-model")==0) {
				customBgFile = argv[++i];
			} else if (strcmp(argv[i],"-modelBg")==0) {
				customBgBgFile = argv[++i];
			} else if (strcmp(argv[i],"-createModel")==0) {
				customBgFile = argv[++i];
				forceFlag = 1;
			} else if (strcmp(argv[i],"-washu")==0) {
				washuFlag = 1;
			} else if (strcmp(argv[i],"-chopify")==0) {
				chopifyFlag = 1;
			} else if (strcmp(argv[i],"-relative")==0) {
				relativeFlag = 1;
			} else if (strcmp(argv[i],"-fixed")==0) {
				fixedFlag = 1;
			} else if (strcmp(argv[i],"-fullModel")==0) {
				fullModelFlag = 1;
			} else if (strcmp(argv[i],"-quickModel")==0) {
				fullModelFlag = 0;
			} else if (strcmp(argv[i],"-override")==0) {
				override= 1;
			} else if (strcmp(argv[i],"-chr2")==0) {
				chr2 = argv[++i];
			} else if (strcmp(argv[i],"-i")==0) {
				inputInteractions = argv[++i];
				actionFlag = actionFlag | HIC_MASK_NO_MATRIX;
			} else if (strcmp(argv[i],"-peakStats")==0) {
				rawbedFilePrefix = argv[++i];
				actionFlag = actionFlag | HIC_MASK_INTERACTIONSTATSFILES;
			} else if (strcmp(argv[i],"-4C")==0) {
				fourCbedFile = argv[++i];
				actionFlag = actionFlag | HIC_MASK_INTERACTION4CBEDFILE;
			} else if (strcmp(argv[i],"-iraw")==0) {
				rawbedFile = argv[++i];
				actionFlag = actionFlag | HIC_MASK_INTERACTIONBEDFILE;
			} else if (strcmp(argv[i],"-irawtags")==0) {
				rawbedFile = argv[++i];
				actionFlag = actionFlag | HIC_MASK_INTERACTIONBEDFILE;
				actionFlag = actionFlag | HIC_MASK_INTERACTIONTAGFILE;
			} else if (strcmp(argv[i],"-ped")==0) {
				bgDirectory = argv[++i];
			} else if (strcmp(argv[i],"-restrictionSite")==0) {
				restrictionSite = argv[++i];	
			} else if (strcmp(argv[i],"-mis")==0) {
				sscanf(argv[++i], "%d", &maxMisMatches);
			} else if (strcmp(argv[i],"-cpu")==0) {
				sscanf(argv[++i],"%d",&numCPUs);
			} else if (strcmp(argv[i],"-size")==0) {
				sscanf(argv[++i],"%d",&size);
			} else if (strcmp(argv[i],"-nolog")==0) {
				nologFlag = 1;
			} else if (strcmp(argv[i],"-minDist")==0) {
				sscanf(argv[++i],"%d",&minDist);
			} else if (strcmp(argv[i],"-maxDist")==0) {
				sscanf(argv[++i],"%d",&maxDist);
			} else if (strcmp(argv[i],"-std")==0) {
				sscanf(argv[++i],"%lf",&stdFilter);
			} else if (strcmp(argv[i],"-min")==0) {
				sscanf(argv[++i],"%lf",&minFilter);
			} else if (strcmp(argv[i],"-b")==0) {
				i++;
				while (argv[i][0] != '-') {
					bedFiles[numBedFiles++] = argv[i++];
					if (i >= argc) break;
				}
				i--;
			} else if (strcmp(argv[i],"-d")==0) {
				i++;
				while (argv[i][0] != '-') {
					tagDirs[numTagDirs++] = argv[i++];
					if (i >= argc) break;
				}
				i--;
			} else if (strcmp(argv[i],"-g")==0) {
				i++;
				while (argv[i][0] != '-') {
					geneFiles[numGeneFiles++] = argv[i++];
					if (i >= argc) break;
				}
				i--;
			} else if (strcmp(argv[i],"-res")==0) {
				i++;
				if (strcmp(argv[i],"site") == 0) {
					siteFlag = 1;
				} else {
					sscanf(argv[i],"%d",&resolution);
				}
				defaultResolution = 0;
			} else if (strcmp(argv[i],"-superRes")==0) {
				i++;
				if (strcmp(argv[i],"site") == 0) {
					siteFlag = 1;
				} else {
					sscanf(argv[i],"%d",&superRes);
					if (defaultResolution==1) {
						resolution = superRes;
					}
				}
				defaultResolution = 0;
			} else if (strcmp(argv[i],"-hist")==0) {
				sscanf(argv[++i],"%d",&resolution);
				defaultResolution = 0;
				expectMatrixFlag=1;
				histFlag = 1;
			} else if (strcmp(argv[i],"-start")==0) {
				sscanf(argv[++i],"%d",&start1);
				if (start2 == -1) start2 = start1;
			} else if (strcmp(argv[i],"-start2")==0) {
				sscanf(argv[++i],"%d",&start2);
			} else if (strcmp(argv[i],"-end")==0) {
				sscanf(argv[++i],"%d",&end1);
				if (end2 == 2000000000) end2 = end1;
			} else if (strcmp(argv[i],"-end2")==0) {
				sscanf(argv[++i],"%d",&end2);
			} else if (strcmp(argv[i],"-clusterFixed")==0) {
				//sscanf(argv[++i],"%lf",&clusterThreshold);
				actionFlag = actionFlag | HIC_MASK_CLUSTER;
				actionFlag = actionFlag | HIC_MASK_CLUSTERFIXED;
			} else if (strcmp(argv[i],"-cluster")==0) {
				//sscanf(argv[++i],"%lf",&clusterThreshold);
				actionFlag = actionFlag | HIC_MASK_CLUSTER;
			} else if (strcmp(argv[i],"-corrDepth")==0) {
				sscanf(argv[++i],"%lf",&minExpReads);
			} else if (strcmp(argv[i],"-corrIters")==0) {
				sscanf(argv[++i],"%d",&clusteringIterations);
			} else if (strcmp(argv[i],"-nomatrix")==0) {
				actionFlag = actionFlag | HIC_MASK_NO_MATRIX;
			} else if (strcmp(argv[i],"-pvalue")==0) {
				sscanf(argv[++i],"%lf",&pvalueThreshold);
			} else if (strcmp(argv[i],"-zscore")==0) {
				sscanf(argv[++i],"%lf",&zscoreThreshold);
			} else if (strcmp(argv[i],"-loci")==0) {
				lociFile = argv[++i];
				//actionFlag = actionFlag | HIC_MASK_LOGPVALUES;
				//actionFlag = actionFlag | HIC_MASK_LOGPVALUESEXACT;
				actionFlag = actionFlag | HIC_MASK_LOCISCORE;
			} else if (strcmp(argv[i],"-interactions")==0) {
				interactionFile = argv[++i];
				actionFlag = actionFlag | HIC_MASK_LOGPVALUES;
				actionFlag = actionFlag | HIC_MASK_LOGPVALUESEXACT;
				actionFlag = actionFlag | HIC_MASK_INTERACTIONS;
			} else if (strcmp(argv[i],"-circos")==0) {
				//interactionFile = argv[++i];
				i++;
				circosPrefix = argv[i];
				actionFlag = actionFlag | HIC_MASK_LOGPVALUES;
				actionFlag = actionFlag | HIC_MASK_LOGPVALUESEXACT;
				actionFlag = actionFlag | HIC_MASK_INTERACTIONS;
			} else if (strcmp(argv[i],"-logp")==0) {
				actionFlag = actionFlag | HIC_MASK_LOGPVALUES;
				actionFlag = actionFlag | HIC_MASK_LOGPVALUESEXACT;
			} else if (strcmp(argv[i],"-logpn")==0) {
				actionFlag = actionFlag | HIC_MASK_LOGPVALUES;
				mask = HIC_MASK_LOGPVALUESEXACT;
				actionFlag = actionFlag & (~mask);
			} else if (strcmp(argv[i],"-corr")==0) {
				actionFlag = actionFlag | HIC_MASK_CORRELATION;
				expectMatrixFlag=1;
				if (actionFlag & HIC_MASK_LOGPVALUES) {
				} else {
					nologFlag = 1;
				}
			} else if (strcmp(argv[i],"-force")==0) {
				forceFlag = 1;
			} else if (strcmp(argv[i],"-bgonly")==0) {
				forceFlag = 1;
				bgOnlyFlag = 1;
			} else if (strcmp(argv[i],"-norm")==0) {
				actionFlag = actionFlag | HIC_MASK_NORM_SEQDEPTH;
				actionFlag = actionFlag | HIC_MASK_NORM_DISTANCE;
			} else if (strcmp(argv[i],"-zscoreNorm")==0) {
				actionFlag = actionFlag | HIC_MASK_NORM_SEQDEPTH;
				actionFlag = actionFlag | HIC_MASK_NORM_DISTANCE;
				actionFlag = actionFlag | HIC_MASK_NORM_ZSCORE;
			} else if (strcmp(argv[i],"-rawAndExpected")==0) {
				actionFlag = actionFlag | HIC_MASK_NORM_SEQDEPTH;
				actionFlag = actionFlag | HIC_MASK_NORM_DISTANCE;
				actionFlag = actionFlag | HIC_MASK_RAWANDEXPECTED;
				expectMatrixFlag=1;
				expectMatrixOutputFile = argv[++i];
			} else if (strcmp(argv[i],"-raw")==0) {
				mask = HIC_MASK_NORM_SEQDEPTH | HIC_MASK_NORM_DISTANCE;
				actionFlag = actionFlag & (~mask);
			} else if (strcmp(argv[i],"-expected")==0) {
				actionFlag = actionFlag | HIC_MASK_EXPECTED;
			} else if (strcmp(argv[i],"-simpleNorm")==0) {
				mask = HIC_MASK_NORM_DISTANCE;
				actionFlag = actionFlag & (~mask);
				actionFlag = actionFlag | HIC_MASK_NORM_SEQDEPTH;
			} else if (strcmp(argv[i],"-p")==0) {
				peakFile1 = argv[++i];
			} else if (strcmp(argv[i],"-p2")==0) {
				peakFile2 = argv[++i];
			} else if (strcmp(argv[i],"-o")==0) {
				outputFile = argv[++i];
			} else {
				fprintf(stderr, "!!! Couldn't recognize \"%s\" !!!\n", argv[i]);
				printCMD();
			}
		}
	}
				
	if ((actionFlag & HIC_MASK_INTERACTION4CBEDFILE
			|| actionFlag & HIC_MASK_INTERACTIONSTATSFILES)
				&& numCPUs > 1) {
		fprintf(stderr, "!!! Warning, options not compatible with -cpu %d, switching to single CPU\n", numCPUs);
		numCPUs = 1;
	}

	if (peakFile1 == NULL && histFlag) {
		fprintf(stderr, "!!! Histogram mode only works with a peak file (\"-p <peakfile>\")!!!\n");
		exit(0);
	}

	if (siteFlag) {
		resolution = -1;
		fixedFlag = 1;
	} else {
		if (superRes < 0) {
			superRes = resolution;
		}
	}
	params->setRes(resolution, superRes);
	params->minDist = minDist;
	params->maxDist = maxDist;
	params->boundaryScale = boundaryScale;
	if (nologFlag == 1) params->logFlag = 0;
	params->minExpReads = minExpReads;


	//Will need to check this for site resolution...
	if (randomizeModel != NULL) {
		HiCBgModel* model = new HiCBgModel(superRes,randomizeModel,1);
		FILE* fp = stdout;
		if (outputFile != NULL) {
			fp = fopen(outputFile, "w");
			if (fp == NULL) {
				fprintf(stderr, "!!! Couldn't open %s for writing!!!\n", outputFile);
				exit(0);
			}
		}
		model->createRandomReads(fp,randomReadTotal);
		if (fp != stdout) fclose(fp);
		exit(0);
	}

	//fprintf(stderr, "chr1=%s\tstart1=%d\tend1=%d\n", chr1, start1, end1);
	//fprintf(stderr, "chr2=%s\tstart2=%d\tend2=%d\n", chr2, start2, end2);

	TagLibrary* tags = NULL;
	if (strcmp(directory,"none")==0) {
		//placeholder
		tags = new TagLibrary(directory);
	} else {
		tags = new TagLibrary(directory);
		tags->readTagDirectory();
		tags->setMaxTBP(maxTbp);
		//fprintf(stderr, "ASDF=%d\n", tags->pairedEndFlag);
		if (!tags->pairedEndFlag) {
			fprintf(stderr, "!!! Tag Directory %s doesn't appear to be a Paired-end Tag Directory!!!\n",
						directory);
			exit(0);
		}
		tags->numCPUs = numCPUs;
		tags->minDistHiC = minDist;
		tags->maxDistHiC = maxDist;
	}

	if (siteFlag) {
		resolution = -1;
		if (restrictionSite == NULL) {
			restrictionSite = tags->restrictionSite;	
			if (restrictionSite == NULL) {
				fprintf(stderr, "!!! Need to specify a restriction site for \"site based\" analysis!!!\n");
				exit(0);
			} else {
				fprintf(stderr, "\tUsing %s as the restriction site\n", restrictionSite);
			}
		}
	}

	TagLibrary* bgtags = NULL;
	if (bgDirectory != NULL) {
		bgtags = new TagLibrary(bgDirectory);
		bgtags->readTagDirectory();
		bgtags->setMaxTBP(maxTbp);
		if (!bgtags->pairedEndFlag) {
			fprintf(stderr, "!!! Background Tag Directory %s doesn't appear to be a Paired-end Tag Directory!!!\n",
						bgDirectory);
			exit(0);
		}
		bgtags->numCPUs = numCPUs;
		bgtags->minDistHiC = minDist;
		bgtags->maxDistHiC = maxDist;
	}

	GenomeInteractionLibrary* gil = new GenomeInteractionLibrary();
	GenomeInteractionLibrary* gilBg = new GenomeInteractionLibrary();

	PeakLibrary* peaks1 = NULL;
	PeakLibrary* peaks2 = NULL;

	if (inputInteractions != NULL) {
		int resEstimate = gil->read(inputInteractions);
		if (resEstimate == -1) {
			fprintf(stderr, "!!! Something is wrong with your input interaction file: %s!!!\n",
								inputInteractions);
			exit(0);
		}
		peaks1 = gil->peakList;

		if (defaultResolution == 1 && siteFlag == 0) {
			fprintf(stderr, "\tUsing Estimated Resolution Size: %d\n", resEstimate);
			resolution = resEstimate;
			superRes = resEstimate;
		} else {
			fprintf(stderr, "\tUsing Given Resolution Size: %d\n", superRes);
			gil->setInteractionPeakSize(superRes);
			if (fixedFlag==0) {
				peaks1->setPeakSize(superRes);
			}
		}
	} else {
		if (superRes > resolution) {
			gil->removeOverlap=1;
			gilBg->removeOverlap=1;
		}
	}
	if (strcmp(directory,"none")==0 && circosPrefix != NULL) {
		printCircosFiles(circosPrefix,resolution,tags,peaks1,peaks2,start1,end1,
					start2,end2,geneFiles,numGeneFiles,tagDirs,numTagDirs,bedFiles, numBedFiles,bgtags);

		//needs to happen after makeHiCMatrix
		sprintf(fname,"%s.circos.interactions.txt",circosPrefix);
		FILE* fp = fopen(fname,"w");
		gil->printCircos(fp,maxCircosInteractions);
		fclose(fp);

		sprintf(fname,"circos -conf %s.circos.conf",circosPrefix);
		int ok = system(fname);
		if (ok) {
			fprintf(stderr, "\t!!! Problem running circos.  Command:\n");
			fprintf(stderr, "\t\t%s\n",fname);
		} else {
			fprintf(stderr, "\tCircos output should be in %s.circos.png\n",circosPrefix);
		}
		exit(0);
	}

		

	// prepare file to output interesting interaction reads
	if (rawbedFile != NULL ) {
		gil->bedFile = fopen(rawbedFile,"w");
		if (gil->bedFile == NULL) {
			fprintf(stderr, "!!! Could not open file %s for writing!!!\n", rawbedFile);
			exit(0);
		}
		if (actionFlag & HIC_MASK_INTERACTIONBEDFILE && !(actionFlag & HIC_MASK_INTERACTIONTAGFILE)) {
			fprintf(gil->bedFile,"track name=\"Paired-end reads for %s\" itemRgb=On visibility=3 colorByStrand=\"255,0,0 255,0,0\"\n",directory);
		}
	}
	if (fourCbedFile != NULL ) {
		gil->fourCbedFile = fopen(fourCbedFile,"w");
		if (gil->fourCbedFile == NULL) {
			fprintf(stderr, "!!! Could not open file %s for writing!!!\n", fourCbedFile);
			exit(0);
		}
		if (actionFlag & HIC_MASK_INTERACTION4CBEDFILE) {
			fprintf(gil->fourCbedFile,"track name=\"In silico 4C for %s\" type=bedGraph yLineMark=\"0.0\" alwaysZeron=on maxHeightPixels=100:75:11 visibility=full autoScale=on\n",directory);
		}
	}

	// prepare files for regional stats
	if (actionFlag & HIC_MASK_INTERACTIONSTATSFILES && rawbedFilePrefix != NULL) {
		if (numCPUs > 1) {
			fprintf(stderr, "!!! Unfortunately, you cannot use -peakStats option with multiple CPUs for now\n");
			exit(0);
		}
		char* filename = new char[100000];
		sprintf(filename, "%s.coverage.bedGraph", rawbedFilePrefix);
		FILE* fp = fopen(filename,"w");
		if (fp == NULL) {
			fprintf(stderr, "!!! Could not open %s for writing!!!\n", filename);
			exit(0);
		}
		gil->bedFileCoverage = fp;
		sprintf(filename, "%s.interFrac.bedGraph", rawbedFilePrefix);
		fp = fopen(filename,"w");
		if (fp == NULL) {
			fprintf(stderr, "!!! Could not open %s for writing!!!\n", filename);
			exit(0);
		}
		gil->bedFileInterFrac = fp;
		sprintf(filename, "%s.localFrac.bedGraph", rawbedFilePrefix);
		fp = fopen(filename,"w");
		if (fp == NULL) {
			fprintf(stderr, "!!! Could not open %s for writing!!!\n", filename);
			exit(0);
		}
		gil->bedFileLocalFrac = fp;
		delete []filename;


		fprintf(gil->bedFileCoverage,"track name=\"Interaction Read Coverage for %s\" type=bedGraph autoScale=on\n",directory);
		fprintf(gil->bedFileInterFrac,"track name=\"Fraction of Interchromosomal interactions for %s\" type=bedGraph autoScale=off viewLimits=0:1\n",directory);
		fprintf(gil->bedFileLocalFrac,"track name=\"Fraction of local/noise/re-ligation interactions for %s\" type=bedGraph autoScale=off viewLimits=0:1\n",directory);
	}
		

	// prepare file to report significant region stats
	if (gil != NULL && lociFile != NULL) {
		gil->lociFile = fopen(lociFile,"w");
		if (gil->lociFile == NULL) {
			fprintf(stderr, "!!! Could not open file %s for writing!!!\n", lociFile);
			exit(0);
		}
		fprintf(gil->lociFile, "PositionID\tchr\tstart\tend\tstrand\tSum(logp)\n");
	}


	// prepare file to output significant interactions
	if (gil != NULL && (circosPrefix != NULL || interactionFile != NULL)) {
		gil->recordInteractions = 1;
		gilBg->recordInteractions = 1;
		if (interactionFile != NULL) {
			sprintf(fname,"%s",interactionFile);
			gil->interactionFile = fopen(fname,"w");
			sprintf(fname,"%s.bg.txt",interactionFile);
			//gilBg->interactionFile = fopen(fname,"w");
		}
	}

	if (minDist == -1) minDist = resolution/2;
	gil->resolution = superRes;
	gilBg->resolution = superRes;
	gil->threshold = log(pvalueThreshold);
	gilBg->threshold = log(pvalueThreshold);
	gil->zscoreThreshold = zscoreThreshold;
	gilBg->zscoreThreshold = zscoreThreshold;
	gil->minDist = minDist;
	gilBg->minDist = minDist;
	gil->maxDist = maxDist;
	gilBg->maxDist = maxDist;
	int peaks1IsGenomeFlag = 0;

	if (peaks1 == NULL) {
		if (peakFile1 == NULL) {
			if (siteFlag) {
				//stuff is not yet ready
				//peaks1 = tags->getCoverageRestrictionFragments(chr1,start1,end1,superRes,
				//						restrictionSite,maxMisMatches,genomeDirectory);
			} else {
				peaks1 = tags->getCoveragePeaks(chr1,start1,end1,resolution,superRes);
			}
			if (chr1 == NULL) {
				peaks1IsGenomeFlag = 1;
			}
		} else {
			peaks1 = new PeakLibrary(peakFile1,PEAK_READ_MODE_NORMAL);
			if (chopifyFlag) {
				int preN = peaks1->numPeaks;
				PeakLibrary* p = peaks1->getCoveragePeaks(resolution,superRes);
				delete peaks1;
				peaks1 = p;
				fprintf(stderr, "\tChopify: %d peaks chopped into %d peaks separated by %d of size %d\n", 
															preN,peaks1->numPeaks,resolution,superRes);
			} else if (fixedFlag == 0) {
				peaks1->setPeakSize(superRes);
			}
		}
		if (peakFile2 == NULL && vsGenomeFlag == 0) {
			if (chr1 != chr2 || start1 != start2 || end1 != end2) {
				peaks2 = tags->getCoveragePeaks(chr2,start2,end2,resolution,superRes);
			}
		} else {
			peaks2 = new PeakLibrary(peakFile2,PEAK_READ_MODE_NORMAL);
			if (chopifyFlag) {
				int preN = peaks2->numPeaks;
				PeakLibrary* p = peaks2->getCoveragePeaks(resolution,superRes);
				delete peaks2;
				peaks2 = p;
				fprintf(stderr, "\tChopify: %d peaks chopped into %d peaks separated by %d of size %d\n",
															preN,peaks2->numPeaks,resolution,superRes);
			} else if (fixedFlag == 0) {
				peaks2->setPeakSize(superRes);
			}
		}
	}
	if (vsGenomeFlag) {
		if (peaks2 != NULL) delete peaks2;
		peaks2 = tags->getCoveragePeaks(NULL, -1, 2000000000,resolution,superRes);
	}

	Doubletable* chrSizes = tags->getChrSizesFromTags();
	if (peaks1 != NULL) peaks1->setMaxChrPositions(chrSizes);
	if (peaks2 != NULL) peaks2->setMaxChrPositions(chrSizes);


	if (poutFile != NULL) {
		FILE* fp = fopen(poutFile, "w");
		if (fp==NULL) {
			fprintf(stderr, "!!! Couldn't open %s for writing peaks !!!\n", poutFile);
			exit(0);
		}
		peaks1->print(fp);
		fclose(fp);
	}
	if (pout2File != NULL && peaks2 != NULL) {
		FILE* fp = fopen(pout2File, "w");
		if (fp==NULL) {
			fprintf(stderr, "!!! Couldn't open %s for writing 2nd set of peaks !!!\n", pout2File);
			exit(0);
		}
		peaks2->print(fp);
		fclose(fp);
	}

	HiCBgModel* bgModel = NULL;
	PeakLibrary* modelPeaks1 = NULL;
	PeakLibrary* modelPeaks2 = NULL;

	if (customBgFile != NULL) {
		bgModel = new HiCBgModel(superRes,customBgFile,1+forceFlag);
		modelPeaks1 = peaks1;
		modelPeaks2 = peaks2;
	} else {	
		bgModel = new HiCBgModel(superRes,directory,0);
		peaks1IsGenomeFlag = 0;
	}
	if (actionFlag & HIC_MASK_NORM_DISTANCE) {
		if (bgModel->ready && forceFlag == 0) {
			fprintf(stderr, "\tFound HiC background model for %d bp resolution (-force to overwrite)\n", superRes);
		} else {
			fprintf(stderr, "\tNo HiC background model found for %d bp resolution. Creating...\n",superRes);

			if (superRes != resolution) {
				fprintf(stderr, "\n\t!!!!Warning... It is HIGHLY recommended that the -res and -superRes are the same\n");
				fprintf(stderr, "\twhen generating a background model!!! (consider building the bg model first...)\n\n");
			}
			//tags->getPETagDistribution(2000,300000000,superRes,0);
			tags->makeHiCBgModel(bgModel,modelPeaks1,modelPeaks2,peaks1IsGenomeFlag,fullModelFlag,params);
			fprintf(stderr, "\tFinished creating background model.\n\n");
			if (bgOnlyFlag) exit(0);
		}	
		bgModel->scale();
	}
	bgModel->stdFilter = stdFilter;
	bgModel->minFilter = minFilter;

	HiCBgModel* bgBgModel = NULL;
	if (bgDirectory != NULL) {
		if (customBgBgFile != NULL) {
			bgBgModel = new HiCBgModel(superRes,customBgBgFile,1);
		} else {	
			bgBgModel = new HiCBgModel(superRes,bgDirectory,0);
		}
		if (actionFlag & HIC_MASK_NORM_DISTANCE) {
			if (bgBgModel->ready && forceFlag == 0) {
				fprintf(stderr, "\tFound HiC background model (for bg directory) for %d bp resolution (-force to overwrite)\n", superRes);
			} else {
				fprintf(stderr, "\tNo HiC background model (for bg directory) found for %d bp resolution. Creating...\n",superRes);
				bgtags->makeHiCBgModel(bgBgModel,modelPeaks1,modelPeaks2,peaks1IsGenomeFlag,fullModelFlag,params);
				fprintf(stderr, "\tFinished creating background model.\n\n");
				if (bgOnlyFlag) exit(0);
			}	
			bgBgModel->scale();
		}
		bgBgModel->stdFilter = stdFilter;
		bgBgModel->minFilter = minFilter;
		bgBgModel->bgFlag = 1;
		//bgBgModel->normalizeVariance(bgModel);
	}




	if (histFlag == 0) {

		//Matrix mode
		double totalMatrix = 0;
		int reportX = peaks1->numPeaks;
		int reportY = 0;
		if (peaks2 != NULL) {
			totalMatrix = (double)peaks1->numPeaks*(double)peaks2->numPeaks;
			reportY = peaks2->numPeaks;
		} else {
			totalMatrix = (double)peaks1->numPeaks*(double)peaks1->numPeaks;
			reportY = peaks1->numPeaks;
		}
		double maxMatrixSize = ((double)HIC_MATRIX_LIMIT)*((double)HIC_MATRIX_LIMIT);
		if (!(actionFlag & HIC_MASK_NO_MATRIX) && relativeFlag==0 && totalMatrix > maxMatrixSize) {
			fprintf(stderr, "!!! Matrix way to big (%d x %d = %.2le elements) !!!\n", reportX, reportY,totalMatrix);
			fprintf(stderr, "\n\tTry using \"-nomatrix\" (max is %.2le)\n\n", maxMatrixSize);
			if (override != 1) {
				fprintf(stderr, "!!! If you're sure about this, add \"-override\"\n");
				exit(0);
			}
		} else if (!(actionFlag & HIC_MASK_NO_MATRIX) && relativeFlag != 0) {
			if (maxDist > 0) {
				relativeFlag = ((double)maxDist)/((double)resolution);
			} else {
				relativeFlag = peaks1->numPeaks;
			}
			if (maxMatrixSize < peaks1->numPeaks*relativeFlag) {
				fprintf(stderr, "!!! When making a relative matrix consider lowering the -maxDist <#> parameter.\n");
				fprintf(stderr, "!!! If you're sure about this, add \"-override\"\n");
				if (override != 1) {
					fprintf(stderr, "!!! If you're sure about this, add \"-override\"\n");
					exit(0);
				}
			}
		}
	
		if (expectMatrixFlag) {
			int p1 = peaks1->numPeaks;
			int p2 = p1;
			if (peaks2 != NULL) p2 = peaks2->numPeaks;
			if (relativeFlag) {
				p2 = relativeFlag*2+1;
			}
			if (bgModel!=NULL) bgModel->initializeExpectMatrix(p1,p2);
			if (bgBgModel != NULL) bgBgModel->initializeExpectMatrix(p1,p2);
		}
		params->relativeFlag = relativeFlag;
		double** matrix = tags->makeHiCMatrix(peaks1, peaks2, superRes,bgModel,actionFlag,gil,params);
		
		double** bgMatrix = NULL;
		if (bgtags != NULL) {
			if (gil->recordInteractions) {
				gil->recordInteractionBg=1;
				gil->sortInteractionIndexes();
			}
			if (gil->bedFile != NULL) {
				if (actionFlag & HIC_MASK_INTERACTIONBEDFILE && !(actionFlag & HIC_MASK_INTERACTIONTAGFILE)) {
					fprintf(gil->bedFile,"track name=\"Paired-end reads for %s\" itemRgb=On visibility=3 colorByStrand=\"0,0,0 0,0,0\"\n",bgDirectory);
				}
			}
			if (gil->fourCbedFile != NULL) {
				if (actionFlag & HIC_MASK_INTERACTION4CBEDFILE) {
					fprintf(gil->fourCbedFile,"track name=\"In silico 4C for %s\" type=bedGraph yLineMark=\"0.0\" alwaysZeron=on maxHeightPixels=100:75:11 visibility=full autoScale=on\n",bgDirectory);
				}
			}

			bgMatrix = bgtags->makeHiCMatrix(peaks1, peaks2, superRes,bgBgModel,actionFlag,gil,params);
		}

		if (gil != NULL && gil->recordInteractions && gil->interactionFile != NULL) {
			gil->print(gil->interactionFile,centerFlag);
		}
		int np2 = peaks1->numPeaks;
		if (peaks2 != NULL) np2 = peaks2->numPeaks;
		if (relativeFlag) {
			np2 = relativeFlag*2+1;
		}

		if (nologFlag && matrix != NULL) {
			
			fprintf(stderr, "\tReporting linear (not log2) ratios\n");
			for (int i=0;i<peaks1->numPeaks;i++) {
				for (int j=0;j<np2;j++) {
					matrix[i][j] = pow(2,matrix[i][j]);
				}
			}
			if (bgMatrix != NULL) {
				for (int i=0;i<peaks1->numPeaks;i++) {
					for (int j=0;j<np2;j++) {
						bgMatrix[i][j] = pow(2,bgMatrix[i][j]);
					}
				}
			}
		}


		if (actionFlag & HIC_MASK_NO_MATRIX || matrix == NULL) {
		} else {

			int symmetric = 0;
			if (peaks2 == NULL)  {
				symmetric = 1;
				peaks2 = peaks1;
			}


			//correlation
			if ((actionFlag & HIC_MASK_CORRELATION) && !(actionFlag & HIC_MASK_NO_MATRIX)) {
				if (relativeFlag) {
					fprintf(stderr, "\tWarning - correlation in relative mode might produce a problem...\n");
					exit(0);
				}
				if (peaks1 != peaks2) {
					fprintf(stderr, "\tCorrelation analysis can only be performed with a symetric matrix, skipping...\n");
				} else {

					int maxDiffIndex = -1;
					if (maxDist >= 0) {
						maxDiffIndex = maxDist/resolution;
					}

					double **expectedMatrix = NULL;
					if (minExpReads > 0.0) {
						if (bgModel != NULL && bgModel->expectMatrix != NULL) {
							expectedMatrix = bgModel->expectMatrix;
						}
						//actionFlag = actionFlag | HIC_MASK_EXPECTED;
						//expectedMatrix = tags->makeHiCMatrix(peaks1, peaks2, superRes,bgModel,actionFlag,gil,minExpReads,boundaryScale);
						//mask = HIC_MASK_EXPECTED;
						//actionFlag = actionFlag & (~mask);
					}
					for (int z=0;z<clusteringIterations;z++) {
						if (z > 0) fprintf(stderr, "\tClustering Iteration: %d\n", z+1);
						matrix = tags->calcCorrelationMatrix(matrix, peaks1->numPeaks,expectedMatrix,minExpReads,
																						maxDiffIndex);
					}
					if (0) { //experimental
						for (int i=0;i<peaks1->numPeaks-2;i++) {
							double corr = TagLibrary::calcCorrelation(matrix[i],matrix[i+2],peaks1->numPeaks,
							expectedMatrix[i],expectedMatrix[i+2],minExpReads,maxDiffIndex,i,i+2);
							fprintf(stdout, "%s\t%d\t%d\t%lf\n", peaks1->peakOrder[i+1]->chr,
											peaks1->peakOrder[i+1]->start,peaks1->peakOrder[i+1]->end,corr);
						}
						exit(0);
					}


					if (bgMatrix != NULL) {
						double **expectedMatrixBg = NULL;
						if (minExpReads > 0.0) {
							if (bgBgModel != NULL && bgBgModel->expectMatrix != NULL) {
								expectedMatrixBg = bgBgModel->expectMatrix;
							}
							//actionFlag = actionFlag | HIC_MASK_EXPECTED;
							//expectedMatrixBg = tags->makeHiCMatrix(peaks1, peaks2, superRes,bgModel,actionFlag,gil,minExpReads,boundaryScale);
							//mask = HIC_MASK_EXPECTED;
							//actionFlag = actionFlag & (~mask);
						}
						for (int z=0;z<clusteringIterations;z++) {
							if (z > 0) fprintf(stderr, "\tClustering Bg Iteration: %d\n", z+1);
							bgMatrix = tags->calcCorrelationMatrix(bgMatrix, peaks1->numPeaks,
																expectedMatrixBg,minExpReads,maxDiffIndex);
						}
						double* em = NULL;
						double* emBg = NULL;
						for (int i=0;i<peaks1->numPeaks;i++) {
							if (minExpReads > 0.0) {
								em = expectedMatrix[i];
								emBg = expectedMatrixBg[i];
							}
							double corr = TagLibrary::calcCorrelation(matrix[i],bgMatrix[i],peaks1->numPeaks,
														em,emBg,minExpReads,
														maxDiffIndex,i,i);
							fprintf(stdout, "%s\t%lf\n", peaks1->peakOrder[i]->name,corr);
						}
						exit(0);

						/*if (expectedMatrixBg != NULL) {
							for (int z=0;z<peaks1->numPeaks;z++) {
								delete [](expectedMatrixBg[z]);
							}
							delete []expectedMatrixBg;
						}*/
						expectedMatrixBg=NULL;
					}
					/*if (expectedMatrix != NULL) {
						for (int z=0;z<peaks1->numPeaks;z++) {
							delete [](expectedMatrix[z]);
						}
						delete []expectedMatrix;
					}*/
					expectedMatrix=NULL;
				}
			}


			if (bgMatrix != NULL) {
				for (int i=0;i<peaks1->numPeaks;i++) {
					for (int j=0;j<np2;j++) {
						matrix[i][j] = matrix[i][j]-bgMatrix[i][j];
					}
				}
			}

			fprintf(stderr, "\tPrinting Interaction Matrix: %d x %d", peaks1->numPeaks, np2);
			if (relativeFlag) fprintf(stderr, " (relative mode)");
			fprintf(stderr, "\n");
			FILE* fp = NULL;
			if (outputFile != NULL) {
				fp = fopen(outputFile, "w");
				if (fp == NULL) {
					fprintf(stderr, "!!! Could not open %s for writing !!!\n", outputFile);
					exit(0);
				}
			} else {
				fp = stdout;
			}
			if (washuFlag == 0) {
				fprintf(fp, "HiCMatrix (directory=%s)\tRegions",directory);
				if (relativeFlag) {
					for (int i=0;i<np2;i++) {
						int offset = (i-relativeFlag)*resolution;
						fprintf(fp, "\t%d",offset);
					}
				} else {
					for (int i=0;i<peaks2->numPeaks;i++) {
						fprintf(fp,"\t%s", peaks2->peakOrder[i]->name);
					}
				}
				fprintf(fp,"\n");
			}
			int washuID = 1;
			int minWashUOffset = -1;
			if (minDist > 0) {
				minWashUOffset = minDist/resolution;
			}
			for (int i=0;i<peaks1->numPeaks;i++) {
				if (washuFlag == 0) fprintf(fp,"%s\t%s",peaks1->peakOrder[i]->name,peaks1->peakOrder[i]->name);

				for (int j=0;j<np2;j++) {
					if (washuFlag == 0) {
						fprintf(fp,"\t%.3le",matrix[i][j]);
					} else {
						if (relativeFlag) {
							if (j==relativeFlag) continue;
							int jj = j+i-relativeFlag;
							if (jj < 0 || jj >= peaks1->numPeaks) continue;
							if (abs(j-relativeFlag) < minWashUOffset) continue;
							if (strcmp(peaks1->peakOrder[i]->chr,peaks1->peakOrder[i]->chr)!=0) continue;
							//if (peaks1->peakOrder[i]->chr != peaks1->peakOrder[jj]->chr) continue;
							fprintf(fp,"%s\t%d\t%d\t%s:%d-%d,%.3lf\t%d\t.\n",peaks1->peakOrder[i]->chr,
										peaks1->peakOrder[i]->start,peaks1->peakOrder[i]->end,
										peaks1->peakOrder[jj]->chr,peaks1->peakOrder[jj]->start,
										peaks1->peakOrder[jj]->end,matrix[i][j],washuID++);
							//fprintf(fp,"%s\t%d\t%d\t%s:%d-%d,%.3lf\t%d\t.\n",peaks1->peakOrder[jj]->chr,
							//			peaks1->peakOrder[jj]->start,peaks1->peakOrder[jj]->end,
							//			peaks1->peakOrder[i]->chr,peaks1->peakOrder[i]->start,
							//			peaks1->peakOrder[i]->end,matrix[i][j],washuID++);
						} else {
							if (peaks1 == peaks2 && j == i) continue;
							//if (peaks1 == peaks2 && j <= i) continue;
							fprintf(fp,"%s\t%d\t%d\t%s:%d-%d,%.3lf\t%d\t.\n",peaks1->peakOrder[i]->chr,
										peaks1->peakOrder[i]->start,peaks1->peakOrder[i]->end,
										peaks2->peakOrder[j]->chr,peaks2->peakOrder[j]->start,
										peaks2->peakOrder[j]->end,matrix[i][j],washuID++);
							//fprintf(fp,"%s\t%d\t%d\t%s:%d-%d,%.3lf\t%d\t.\n",peaks2->peakOrder[j]->chr,
							//			peaks2->peakOrder[j]->start,peaks2->peakOrder[j]->end,
							//			peaks1->peakOrder[i]->chr,peaks1->peakOrder[i]->start,
							//			peaks1->peakOrder[i]->end,matrix[i][j],washuID++);
							if (peaks1 != peaks2) {
								//fprintf(fp,"%s\t%d\t%d\t%s:%d-%d,%.3le\t%d\t.\n",peaks2->peakOrder[j]->chr,
								//		peaks2->peakOrder[j]->start,peaks2->peakOrder[j]->end,
								//		peaks1->peakOrder[i]->chr,peaks1->peakOrder[i]->start,
								//		peaks1->peakOrder[i]->end,matrix[i][j],washuID++);
							}
						}
					}
				}

				if (washuFlag == 0) fprintf(fp,"\n");
			}
			if (outputFile != NULL) {
				fclose(fp);
			}
			if (expectMatrixOutputFile != NULL && bgModel != NULL && bgModel->expectMatrix != NULL) {
				FILE* fexp = fopen(expectMatrixOutputFile,"w");
				if (fexp == NULL) {
					fprintf(stderr, "!!! Could not open %s for writing !!!\n", expectMatrixOutputFile);
					exit(0);
				}
				double** expectMatrix = bgModel->expectMatrix;

				fprintf(fexp, "Expected HiCMatrix (directory=%s)\tRegions",directory);
				if (relativeFlag) {
					for (int i=0;i<np2;i++) {
						int offset = (i-relativeFlag)*resolution;
						fprintf(fp, "\t%d",offset);
					}
				} else {
					for (int i=0;i<np2;i++) {
						fprintf(fexp,"\t%s", peaks2->peakOrder[i]->name);
					}
				}
				fprintf(fexp,"\n");
				for (int i=0;i<peaks1->numPeaks;i++) {
					fprintf(fexp,"%s\t%s",peaks1->peakOrder[i]->name,peaks1->peakOrder[i]->name);
					for (int j=0;j<np2;j++) {
						if (1) {//matrix[i][j] > -1.99) {
							//fprintf(fp,"\t%.3lf",matrix[i][j]);
							fprintf(fexp,"\t%.3le",expectMatrix[i][j]);
						} else {
							fprintf(fexp,"\tNA");
						}
					}
					fprintf(fexp,"\n");
				}
				fclose(fexp);
			}

			//clustering
			if ((actionFlag & HIC_MASK_CLUSTER) && !(actionFlag & HIC_MASK_NO_MATRIX)) {
				if (peaks1 != peaks2) {
					fprintf(stderr, "\tClustering can only be performed with a symetric matrix, skipping...\n");
				} else {
					fprintf(stderr, "\tClustering regions\n");
					if (actionFlag & HIC_MASK_LOGPVALUES) {
					} else {
						for (int i=0;i<peaks1->numPeaks;i++) {
							for (int j=0;j<peaks1->numPeaks;j++) {
								//matrix[i][j] = exp(matrix[i][j]);
								matrix[i][j] = -1*matrix[i][j];
							}
						}
					}
					int removeBadRegions=1;
					if (removeBadRegions && bgModel != NULL) {
						int newNumPeaks = 0;
						PeakLibrary* newpeaks = new PeakLibrary();
						int* map = new int[peaks1->numPeaks];
						for (int i=0;i<peaks1->numPeaks;i++) {
							if (bgModel->minCoverage <= peaks1->peakOrder[i]->v 
									&& bgModel->maxCoverage >= peaks1->peakOrder[i]->v) {
								map[newNumPeaks++] = i;
								newpeaks->addPeak(peaks1->peakOrder[i]);
							}
						}
	fprintf(stderr, "\tOmitting %d(of %d) peaks due to low/high read depth (>%.1lf std or <%.2lf of average depth)\n", newNumPeaks,peaks1->numPeaks,stdFilter,minFilter);
						newpeaks->setDefaultPeakOrder();
						double** nmatrix = new double*[newNumPeaks];
						for (int i=0;i<newNumPeaks;i++) {
							nmatrix[i] = matrix[map[i]];
							for (int j=0;j<newNumPeaks;j++) {
								nmatrix[i][j] = nmatrix[i][map[j]];
							}

						}
						delete []matrix;
						delete []map;
						delete peaks1;
						matrix = nmatrix;
						peaks1 = newpeaks;
					}
					//Clustering
					char** pnames = new char*[peaks1->numPeaks];
					for (int i=0;i<peaks1->numPeaks;i++) {
						pnames[i] = peaks1->peakOrder[i]->name;
					}
					TreeCluster* cluster = new TreeCluster(matrix,peaks1->numPeaks);
					if (actionFlag & HIC_MASK_CLUSTERFIXED) {
						cluster->fixedPositionFlag = 1;
					}
					cluster->fixedClusterScore = 0;
					cluster->cluster();
			
					if (actionFlag & HIC_MASK_LOGPVALUES) {
					} else {
						for (int i=0;i<peaks1->numPeaks;i++) {
							for (int j=0;j<peaks1->numPeaks;j++) {
								matrix[i][j] = -1*matrix[i][j];
							}
						}
					}
					cluster->printCDT(outputFile, pnames, matrix);
					//cluster->printClusters(outputFile, pnames, clusterThreshold);
					if (actionFlag & HIC_MASK_CLUSTERFIXED) {
						if (outputFile != NULL) {
							sprintf(fname, "%s.clusterSizes.txt", outputFile);
						} else {
							sprintf(fname, "out.clusterSizes.txt");
						}
						FILE* fpSizes = fopen(fname,"w");
						cluster->printClusterSizes(fpSizes,1000);
						fclose(fpSizes);
					}
					delete cluster;
					delete []pnames;
				}
			}

			if (matrix != NULL) {
				for (int i=0;i<peaks1->numPeaks;i++) {
					delete [](matrix[i]);
				}
				delete []matrix;
			}
			if (bgMatrix != NULL) {
				for (int i=0;i<peaks1->numPeaks;i++) {
					delete [](bgMatrix[i]);
				}
				delete []bgMatrix;
			}
			if (symmetric) peaks2 = NULL;
			
		}
		
			
		if (circosPrefix != NULL) {
			printCircosFiles(circosPrefix,resolution,tags,peaks1,peaks2,start1,end1,
						start2,end2,geneFiles,numGeneFiles,tagDirs,numTagDirs,bedFiles, numBedFiles,bgtags);

			//needs to happen after makeHiCMatrix
			sprintf(fname,"%s.circos.interactions.txt",circosPrefix);
			FILE* fp = fopen(fname,"w");
			gil->printCircos(fp,maxCircosInteractions);
			fclose(fp);

			sprintf(fname,"circos -conf %s.circos.conf",circosPrefix);
			int ok = system(fname);
			if (ok) {
				fprintf(stderr, "\t!!! Problem running circos.  Command:\n");
				fprintf(stderr, "\t\t%s\n",fname);
			} else {
				fprintf(stderr, "\tCircos output should be in %s.circos.png\n",circosPrefix);
			}
		}

	} else {

		//Histogram Mode
		if (size <= 0) size = resolution*100;
		//peaks1 = new PeakLibrary(peakFile1,PEAK_READ_MODE_NORMAL);
			
		tags->makeHiCHistogram(outputFile,peaks1, size, superRes,bgModel,actionFlag,params);

	}
	
	if (poutFile != NULL) {
		FILE* fp = fopen(poutFile, "w");
		if (fp==NULL) {
			fprintf(stderr, "!!! Couldn't open %s for writing peaks !!!\n", poutFile);
			exit(0);
		}
		peaks1->print(fp);
		fclose(fp);
	}
	if (pout2File != NULL && peaks2 != NULL) {
		FILE* fp = fopen(pout2File, "w");
		if (fp==NULL) {
			fprintf(stderr, "!!! Couldn't open %s for writing 2nd set of peaks !!!\n", pout2File);
			exit(0);
		}
		peaks2->print(fp);
		fclose(fp);
	}


	if (gil != NULL && gil->lociFile != NULL) fclose(gil->lociFile);
	if (gil != NULL && gil->interactionFile != NULL) fclose(gil->interactionFile);
	if (chrSizes != NULL) delete chrSizes;

}

void printCircosFiles(char* prefix, int resolution, TagLibrary* tags, PeakLibrary* peaks1, PeakLibrary* peaks2,
				int start1, int end1, int start2, int end2,char** geneFiles, int numGeneFiles,
				char** tagDirs, int numTagDirs,char** bedFiles, int numBedFiles, TagLibrary* bgtags) {


	Doubletable* chrSizes = tags->getChrSizesFromTags();
	Inttable* chrs = new Inttable(1000);
	char** keys = peaks1->chrs->keys();
	for (int i=0;i<peaks1->chrs->total;i++) {
		double size = chrSizes->search(keys[i]);
		if (size > EMPTY_DOUBLE_CHECK) {
			chrs->insert((int)size,keys[i]);
		}
		delete [](keys[i]);
	}
	delete []keys;
	if (peaks2 != NULL && peaks1 != peaks2) {
		char** keys = peaks2->chrs->keys();
		for (int i=0;i<peaks2->chrs->total;i++) {
			double size = chrSizes->search(keys[i]);
			if (size > EMPTY_DOUBLE_CHECK) {
				chrs->insert((int)size,keys[i]);
			}
			delete [](keys[i]);
		}
		delete []keys;
	}

	char** currentChrs = chrs->keys();
	qsort(currentChrs,chrs->total,sizeof(char*),&chrcmp);
	if (chrs->total < 1) {
		fprintf(stderr, "!!! Something seems to be wrong - no valid chromosomes are available...\n");
	}


	char* fname = new char[10000];
	char* fname2 = new char[10000];

	sprintf(fname,"%s.circos.conf",prefix);
	FILE* cfp = fopen(fname,"w");

	fprintf(cfp,"# automated circos file generated by HOMER (analyzeHiC)\n");
	fprintf(cfp,"\n<colors>\n\t<<include etc/colors.conf>>\n\t<<include etc/brewer.conf>>\n</colors>\n\n");
	fprintf(cfp,"\n<fonts>\n\t<<include etc/fonts.conf>>\n</fonts>\n");

	fprintf(cfp,"\n<ideogram>\n");
	fprintf(cfp,"\n<spacing>\n");
	fprintf(cfp,"\tdefault = 50u\n");
	fprintf(cfp,"\tbreak = 10u\n");
	fprintf(cfp,"\taxis_break_at_edge = yes\n");
	fprintf(cfp,"\taxis_break = yes\n");
	fprintf(cfp,"\taxis_break_style = 2\n");
	fprintf(cfp,"\n\t<break_style 1>\n");
	fprintf(cfp,"\t\tstroke_color = black\n");
	fprintf(cfp,"\t\tfill_color = blue\n");
	fprintf(cfp,"\t\tthickness = 0.25r\n");
	fprintf(cfp,"\t\tstroke_thickness = 2\n");
	fprintf(cfp,"\t</break>\n");
	fprintf(cfp,"\n\t<break_style 2>\n");
	fprintf(cfp,"\t\tstroke_color = black\n");
	fprintf(cfp,"\t\tthickness = 1.5r\n");
	fprintf(cfp,"\t\tstroke_thickness = 3\n");
	fprintf(cfp,"\t</break>\n");
	fprintf(cfp,"</spacing>\n");
	fprintf(cfp,"\nthickness = 10p\n");
	fprintf(cfp,"stroke_thickness = 2\n");
	fprintf(cfp,"stroke_color = black\n");
	fprintf(cfp,"fill = no\n");
	fprintf(cfp,"fill_color = black\n");
 	fprintf(cfp,"radius = 0.80r\nshow_label = yes\nlabel_with_tag = yes\nlabel_font = default\n");
	fprintf(cfp,"label_radius = dims(ideogram,radius) + 0.075r\nlabel_size = 60p\n");
	fprintf(cfp,"band_stroke_thinkness = 2\nshow_bands = yes\nfill_bands = yes\n");

	fprintf(cfp,"\n</ideogram>\n");


	fprintf(cfp,"\n<image>\n");
	fprintf(cfp,"\tdir = ./\n");
	fprintf(cfp,"\tfile = ./%s.circos.png\n",prefix);
	fprintf(cfp,"\t24bit = yes\n\tpng = yes\n\tsvg = yes\n\tradius = 1500p\n\tbackground = white\n");
	fprintf(cfp,"\tangle_offset = -90\n\tauto_alpha_colors = yes\n\tauto_alpha_steps = 20\n");
	fprintf(cfp,"</image>\n");

	fprintf(cfp,"\nkaryotype = %s.circos.karyotype.txt\n",prefix);
	fprintf(cfp,"\nchromosomes_display_default = no\n");
	fprintf(cfp,"\nchromosomes_units = %d\n", resolution);
	fprintf(cfp,"\nchromosomes = ");

	double totalBp = 0;	
	for (int i=0;i<chrs->total;i++) {
		int curLength = chrs->search(currentChrs[i]);

		if (i>0) fprintf(cfp,";");
		fprintf(cfp,"%s",currentChrs[i]);
		if (start1 > 1 || end1 < curLength) {
			int end = end1;
			if (end > chrs->search(currentChrs[i])) {
				end = chrs->search(currentChrs[i]);
			}
			int chrCurSize = end - start1 + 1;
			end = end/resolution+1;
			int start = start1/resolution-1;
			if (start < 0) start = 0;
			fprintf(cfp,":%d-%d",start,end);
			totalBp += (double)(chrCurSize);
			//fprintf(stderr, "adding %lf %d %d\n",(double)chrCurSize,start, end);
		} else {
			totalBp += (double)chrs->search(currentChrs[i]);
			//fprintf(stderr, "adding %lf\n",(double)chrs->search(currentChrs[i]));
		}
	}
	fprintf(stderr, "\tTotal bp to visualize in Circos: %.0lf\n", totalBp);
	double maxSizeForNames = 10e6;	

	double extrasRadius = numTagDirs*0.1 + numGeneFiles*0.15 + numBedFiles*0.025;
	if (totalBp > maxSizeForNames) {
		extrasRadius = numTagDirs*0.1 + numGeneFiles*0.10 + numBedFiles*0.025;
	}

	double defaultInfoSpace = 0.5;

	double extrasScaleFactor = 1.0;
	if (extrasRadius > defaultInfoSpace) {
		extrasScaleFactor = defaultInfoSpace/extrasRadius;
	}
	double interactionRadius = 1.0 - extrasRadius*extrasScaleFactor-0.005;


	double unadj = totalBp/(double)resolution/75.0;
	int ticRes = (int)(log(unadj)/log(10.0));
	int ticThickness = (int)pow(10.0,ticRes);
	int ticThickness2 = (int)pow(10.0,ticRes+1.0);
	if (ticThickness < 1) ticThickness = 1;
	if (ticThickness2 < 10) ticThickness2 = 10;
	//fprintf(stderr, "multiplier: %d %lf totalbp = %lf %d %d\n",ticRes,unadj,totalBp,ticThickness, ticThickness2);
	fprintf(cfp,"\nshow_ticks = yes\nshow_tick_labels=yes\n");
	fprintf(cfp,"\n<ticks>\n");
	fprintf(cfp,"\tradius = dims(ideogram,radius_outer)\n");
	//int multiplyer = (int)(log(ticRes)/log(10.0));
	fprintf(cfp,"\tmultiplier = 1e-6\n");
	fprintf(cfp,"\t<tick>\n");
	fprintf(cfp,"\t\tspacing = %du\n",ticThickness);
	fprintf(cfp,"\t\tsize = 8p\n\t\tthickness = 2p\n\t\tcolor = black\n\t\tshow_label = no\n");
	fprintf(cfp,"\t\tlabel_size = 20p\n\t\tformat = %%.3f\n");
	fprintf(cfp,"\t</tick>\n");
	fprintf(cfp,"\t<tick>\n");
	fprintf(cfp,"\t\tspacing = %du\n",ticThickness2);
	fprintf(cfp,"\t\tsize = 8p\n\t\tthickness = 2p\n\t\tcolor = black\n\t\tshow_label = yes\n");
	fprintf(cfp,"\t\tlabel_size = 40p\n\t\tlabel_offset = 5p\n\t\tformat = %%.3f\n");
	fprintf(cfp,"\t</tick>\n");
	fprintf(cfp,"</ticks>\n");


	
	fprintf(cfp,"\n<links>\n");	
	fprintf(cfp,"\n\tz = 5\n");	
	fprintf(cfp,"\n\tradius = %.2lfr\n",interactionRadius);	
	fprintf(cfp,"\tbezier_radius = 0.4r\n");	
	fprintf(cfp,"\tbezier_radius_purity = 0.45\n");	
	fprintf(cfp,"\tcrest = 0.4\n");	
	fprintf(cfp,"\tperturb = yes\n");	

	fprintf(cfp,"\n\t<link segdup>\n");	
	fprintf(cfp,"\t\tshow = yes\n\t\tcolor = black_a8\n\t\tthickness = 2\n\t\trecord_limit=5000\n");	
	fprintf(cfp,"\t\tfile = %s.circos.interactions.txt\n", prefix);

	fprintf(cfp,"\t\t<rules>\n");
	if (bgtags != NULL) {
		fprintf(cfp,"\t\t\t<rule>\n");
		fprintf(cfp,"\t\t\t\timportance = 210\n");
		fprintf(cfp,"\t\t\t\tcondition = _THICKNESS1_ < 0\n");
		fprintf(cfp,"\t\t\t\tthickness = eval(abs(_THICKNESS1_))\n");
		fprintf(cfp,"\t\t\t\tcolor = red_a3\n");
		fprintf(cfp,"\t\t\t\tz = 10\n");
		fprintf(cfp,"\t\t\t\tflow = continue\n");
		fprintf(cfp,"\t\t\t</rule>\n");
	}
	fprintf(cfp,"\t\t\t<rule>\n");
	fprintf(cfp,"\t\t\t\timportance = 205\n");
	fprintf(cfp,"\t\t\t\tcondition = _THICKNESS1_ > 20\n");
	fprintf(cfp,"\t\t\t\tthickness = 20\n");
	fprintf(cfp,"\t\t\t</rule>\n");
	fprintf(cfp,"\t\t\t<rule>\n");
	fprintf(cfp,"\t\t\t\timportance = 200\n");
	fprintf(cfp,"\t\t\t\tcondition = _THICKNESS1_ < -20\n");
	fprintf(cfp,"\t\t\t\tthickness = 20\n");
	fprintf(cfp,"\t\t\t</rule>\n");
	fprintf(cfp,"\t\t\t<rule>\n");
	fprintf(cfp,"\t\t\t\timportance = 195\n");
	fprintf(cfp,"\t\t\t\tcondition = _THICKNESS1_ < 0\n");
	fprintf(cfp,"\t\t\t\tthickness = eval(abs(_THICKNESS1_))\n");
	fprintf(cfp,"\t\t\t</rule>\n");
	fprintf(cfp,"\t\t\t<rule>\n");
	fprintf(cfp,"\t\t\t\timportance = 190\n");
	fprintf(cfp,"\t\t\t\tcondition = _THICKNESS1_ < 1\n");
	fprintf(cfp,"\t\t\t\tthickness = 1\n");
	fprintf(cfp,"\t\t\t</rule>\n");
	fprintf(cfp,"\t\t</rules>\n");

	fprintf(cfp,"\t</link>\n");
	fprintf(cfp,"\n</links>\n");


	// plots

	int maxColors = 6;
	char** colors = new char*[maxColors];
	colors[0] = (char*)"red";
	colors[1] = (char*)"green";
	colors[2] = (char*)"blue";
	colors[3] = (char*)"orange";
	colors[4] = (char*)"purple";
	colors[5] = (char*)"grey";
	int curColor = 0;

	fprintf(cfp,"<plots>\n");

	char* command = new char[10000];
	int ucscRes = (int) (totalBp / 3e3);
	if (ucscRes < 1) ucscRes = 1;
	fprintf(stderr, "\tSetting UCSC resolution at %d\n", ucscRes);
	char* ucscPos = new char[1000];
	strcpy(ucscPos,"genome");
	if (chrs->total == 1) {
		int start11 = 0;
		if (start1 > start11) start11 = start1;
		sprintf(ucscPos,"%s:%d-%d",currentChrs[0],start11,end1);
	}

	double curRadius = interactionRadius;

	for (int i=0;i<numBedFiles;i++) {

		sprintf(fname,"%s.circos.bed%d.txt", prefix,i+1);

		double radiusThickness = extrasScaleFactor * 0.025;

		fprintf(cfp,"\t<plot>\n");
		fprintf(cfp,"\t\tshow = yes\n\t\tz=5\n\t\ttype = tile\n");
		fprintf(cfp,"\t\tr0 = %.2lfr\n\t\tr1 = %.2lfr\n",curRadius,curRadius+radiusThickness);
		curRadius += radiusThickness;

		fprintf(cfp,"\t\tcolor = %s\n", colors[curColor]);
		if (++curColor >= maxColors) curColor = 0;

		fprintf(cfp,"\t\tlayers = 2\n\t\tmargin = 0.2u\n\t\tthickness=15\n\t\tpadding=2\n");
		fprintf(cfp,"\t\torientation = out\n\t\tlayers_overflow=hide\n");
		//fprintf(cfp,"\t\tbackground = no\n");
		fprintf(cfp,"\t\tfile = %s\n", fname);
		fprintf(cfp,"\t</plot>\n");

		//------ prep data --------
		FILE* fp = fopen(fname,"w");
		PeakLibrary* peaks = new PeakLibrary(bedFiles[i],PEAK_READ_MODE_NORMAL);
		peaks->setDefaultPeakOrder();
		for (int j=0;j<peaks->numPeaks;j++) {
			int good = 1;
			if (chrs->total == 1) {
				if (strcmp(currentChrs[0],peaks->peakOrder[j]->chr)!=0) good = 0;
				if (peaks->peakOrder[j]->end < start1) good = 0;
				if (peaks->peakOrder[j]->start > end1) good = 0;
			}
			if (good) {
				fprintf(fp, "%s %d %d\n", peaks->peakOrder[j]->chr,peaks->peakOrder[j]->start,
											peaks->peakOrder[j]->end);
			}
		}
		fclose(fp);
		delete peaks;
	}

	curColor = 0;
	for (int i=0;i<numTagDirs;i++) {

		double radiusThickness = extrasScaleFactor * 0.10;

		fprintf(cfp,"\t<plot>\n");
		fprintf(cfp,"\t\tshow = yes\n\t\tz=5\n\t\ttype = histogram\n");
		fprintf(cfp,"\t\tr0 = %.2lfr\n\t\tr1 = %.2lfr\n",curRadius,curRadius+radiusThickness);
		curRadius += radiusThickness;
		fprintf(cfp,"\t\tcolor = %s\n\t\tfill_color = %s\n", colors[curColor],colors[curColor]);
		if (++curColor >= maxColors) curColor = 0;
		fprintf(cfp,"\t\tfill_under = yes\n\t\tthickness = 1\n\t\textend_bin = no\n");
		//fprintf(cfp,"\t\tbackground = no\n\t\taxis = no\n");
		
		sprintf(fname,"%s.circos.histogram%d.txt", prefix,i+1);
		fprintf(cfp,"\t\tfile = %s\n\t\tmin = %d\n\t\tmax = %d\n", fname,0,75);
		fprintf(cfp,"\t</plot>\n");

		sprintf(command,"makeUCSCfile %s -circos %s -fsize 1e10 -res %d -o %s",tagDirs[i],ucscPos,
					ucscRes, fname);
		//fprintf(stderr, "cmd: %s\n", command);
		(void)system(command);
	}


	for (int i=0;i<numGeneFiles;i++) {

		sprintf(fname,"%s.circos.genes%d.txt", prefix,i+1);
		sprintf(fname2,"%s.circos.geneNames%d.txt", prefix,i+1);

		double radiusThickness = extrasScaleFactor * 0.025;
		if (totalBp > maxSizeForNames) {
			radiusThickness = extrasScaleFactor * 0.10;
		}

		fprintf(cfp,"\t<plot>\n");
		fprintf(cfp,"\t\tshow = yes\n\t\tz=5\n\t\ttype = tile\n");
		fprintf(cfp,"\t\tr0 = %.2lfr\n\t\tr1 = %.2lfr\n",curRadius,curRadius+radiusThickness);
		curRadius += radiusThickness;

		fprintf(cfp,"\t\tcolor = purple\n");
		fprintf(cfp,"\t\tlayers = 5\n\t\tmargin = 0.2u\n\t\tthickness=15\n\t\tpadding=8\n");
		fprintf(cfp,"\t\torientation = out\n\t\tlayers_overflow=hide\n");
		//fprintf(cfp,"\t\tbackground = no\n");

		fprintf(cfp,"\t\tfile = %s\n", fname);

		fprintf(cfp,"\t</plot>\n");

		//------ text section/gene names --------

		radiusThickness = extrasScaleFactor * 0.10;

		fprintf(cfp,"\t<plot>\n");
		if (totalBp <= maxSizeForNames) {
			fprintf(cfp,"\t\tshow = yes\n");
		} else {
			fprintf(cfp,"\t\tshow = no\n");
		}
		fprintf(cfp,"\t\tz=5\n\t\ttype = text\n");
		fprintf(cfp,"\t\tr0 = %.2lfr\n\t\tr1 = %.2lfr\n",curRadius,curRadius+radiusThickness);
		if (totalBp <= maxSizeForNames) {
			curRadius += radiusThickness;
		}

		fprintf(cfp,"\t\tcolor = black\n");
		fprintf(cfp,"\t\tshow_links = yes\n\t\tlink_dims=4p,4p,8p,4p,4p\n\t\tlink_thickness=4p\n");
		fprintf(cfp,"\t\tlink_color = red\n\t\tlabel_size=25p\n\t\tlabel_font=condensed\n");
		fprintf(cfp,"\t\tpadding = 0p\n\t\trpadding=0p\n");

		fprintf(cfp,"\t\tlabel_snuggle=yes\n\t\tmax_snuggle_distance = 1r\n");
		fprintf(cfp,"\t\tsnuggle_sampling = 2\n\t\tsnuggle_tolerance = 0.25r\n");
		fprintf(cfp,"\t\tsnuggle_link_overlap_test = yes\n\t\tsnuggle_link_overlap_tolerance = 2p\n");
		fprintf(cfp,"\t\tsnuggle_refine = yes\n");

		fprintf(cfp,"\t\tfile = %s\n", fname2);

		fprintf(cfp,"\t</plot>\n");
	
		FILE* fp = fopen(fname,"w");
		FILE* fp2 = fopen(fname2,"w");
		PeakLibrary* peaks = new PeakLibrary(geneFiles[i],PEAK_READ_MODE_NORMAL);
		peaks->setDefaultPeakOrder();
		for (int j=0;j<peaks->numPeaks;j++) {
			int good = 1;
			if (chrs->total == 1) {
				if (strcmp(currentChrs[0],peaks->peakOrder[j]->chr)!=0) good = 0;
				if (peaks->peakOrder[j]->end < start1) good = 0;
				if (peaks->peakOrder[j]->start > end1) good = 0;
			}
			if (good) {
				fprintf(fp, "%s %d %d\n", peaks->peakOrder[j]->chr,peaks->peakOrder[j]->start,
											peaks->peakOrder[j]->end);
				fprintf(fp2, "%s %d %d %s\n", peaks->peakOrder[j]->chr,peaks->peakOrder[j]->start,
											peaks->peakOrder[j]->end,peaks->peakOrder[j]->name);
			}
		}
		fclose(fp);
		fclose(fp2);
		delete peaks;
	}

	fprintf(cfp,"</plots>\n");

	fprintf(cfp,"\n\n<<include etc/housekeeping.conf>>\n");

	fclose(cfp);


	sprintf(fname,"%s.circos.karyotype.txt",prefix);
	cfp = fopen(fname,"w");
	for (int i=0;i<chrs->total;i++) {
		fprintf(cfp,"chr - %s %s 0 %d %s\n", currentChrs[i],currentChrs[i],chrs->search(currentChrs[i]),
												currentChrs[i]);
		delete [](currentChrs[i]);
	}
	fclose(cfp);
	delete []currentChrs;
	delete chrs;
	delete chrSizes;
	delete []fname;
	delete command;
	delete ucscPos;

}


void printCMD() {
	fprintf(stderr, "\n\tUsage: analyzeHiC <PE tag directory> [options]\n");
	fprintf(stderr, "\t\t...\n");
	fprintf(stderr, "\n\tResolution Options:\n");
	fprintf(stderr, "\t\t-res <#> (Resolution of matrix in bp or use \"-res site\" [see below], default: 10000000)\n");
	fprintf(stderr, "\t\t-superRes <#> (size of region to count tags, i.e. make twice the res, default: same as res)\n");
	//fprintf(stderr, "\t\tFull Resolution Options:\n");
	//fprintf(stderr, "\t\t\t-res site (to trigger full restriction site resolution analysis)\n");
	//fprintf(stderr, "\t\t\t-restrictionSite <ACGT...> (restriction site used, -mis <#> to allow mismatches)\n");

	fprintf(stderr, "\n\tRegion Options:\n");
	fprintf(stderr, "\t\t-chr <name> (create matrix on this chromosome, default: whole genome)\n");
	fprintf(stderr, "\t\t-start <#> (start matrix at this position, default:0)\n");
	fprintf(stderr, "\t\t-end <#> (end matrix at this position, default: no limit)\n");
	fprintf(stderr, "\t\t-pos chrN:xxxxxx-yyyyyy (UCSC formatted position instead of -chr/-start/-end)\n");
	fprintf(stderr, "\t\t-chr2 <name>, -start2 <#>, -end2 <#>, or -pos2 (Use these positions on the\n");
	fprintf(stderr, "\t\t\ty-axis of the matrix.  Otherwise the matrix will be sysmetric)\n");
	fprintf(stderr, "\t\t-p <peak file> (specify regions to make matrix, unbalanced, use -p2 <peak file>)\n");
	fprintf(stderr, "\t\t-vsGenome (normally makes a square matrix, this will force 2nd set of peaks to be the genome)\n");
	fprintf(stderr, "\t\t-chopify (divide up peaks into regions the size of the resolution, default: use peak midpoints)\n");
	fprintf(stderr, "\t\t-fixed (do not scale the size of the peaks to that of the resolution)\n");
	fprintf(stderr, "\t\t-relative (use with -maxDist <#>, outputs diagonal of matrix up to maxDistance)\n");
	fprintf(stderr, "\t\t-pout <filename> (output peaks used for analysis as a peak file, -pout2 <file> for 2nd set)\n");

	fprintf(stderr, "\n\tInteraction Matrix Options:\n");
	fprintf(stderr, "\t\t-norm (normalize by dividing each position by expected counts [log ratio], default)\n");
	fprintf(stderr, "\t\t-zscoreNorm (normalize log ratio by distance dependent stddev, add \"-nolog\" to return linear values)\n");
	fprintf(stderr, "\t\t-logp (output log p-values)\n");
	fprintf(stderr, "\t\t-simpleNorm (Only normalize based on total interactions per location [log ratio], not distance)\n");
	fprintf(stderr, "\t\t-raw (report raw interaction counts)\n");
	fprintf(stderr, "\t\t-expected (report expected interaction counts based on average profile)\n");
	fprintf(stderr, "\t\t-rawAndExpected <filename for expectd matrix> (raw counts sent to stdout)\n");
	fprintf(stderr, "\t\t-corr (report Pearson's correlation coeff, use \"-corrIters <#>\" to recursively calculate)\n");
	fprintf(stderr, "\t\t\t-corrDepth <#> (merge regions in correlation so that minimum # expected tags per data point)\n");
	fprintf(stderr, "\t\t-o <filename> (Output file name, default: sent to stdout)\n");
	fprintf(stderr, "\t\t-cluster (cluster regions, uses \"-o\" to name cdt/gtr files, default: out.cdt)\n");
	fprintf(stderr, "\t\t-clusterFixed (clusters adjacent regions, good for linear domains)\n");
	//fprintf(stderr, "\t\t\t(# defined the threshold at which to define clusters/domains, i.e. 0.5)\n");
	fprintf(stderr, "\t\t-std <#> (# of std deviations from mean interactions per region to exclude, default:4)\n");
	fprintf(stderr, "\t\t-min <#> (minimum fraction of average read depth to include in analysis, default: 0.2)\n");
	fprintf(stderr, "\t\t-override (Allow very large matrices to be created... be carful using this)\n");

	fprintf(stderr, "\n\tBackground Options:\n");
	fprintf(stderr, "\t\t-fullModel (perform exhaustive background model calculations, default: approxmiate)\n");
	fprintf(stderr, "\t\t-quickModel (perform approximate background model calculations, default: approxmiate)\n");
	fprintf(stderr, "\t\t-force (force the creation of a fresh genome-wide background model)\n");
	fprintf(stderr, "\t\t-bgonly (quit after creating the background model)\n");
	fprintf(stderr, "\t\t-createModel <custom bg model output file> (Create custome bg from regions specified, i.e. -p/-pos)\n");
	fprintf(stderr, "\t\t-model <custom bg model input file> (Use Custom background model, -modelBg for -ped)\n");
	fprintf(stderr, "\t\t-randomize <bgmodel> <# reads> (and the output is a PE tag file, initail PE tag directory not used\n");
	fprintf(stderr, "\t\t\tUse makeTagDirectory ... -t outputfile to create the new directory)\n");

	fprintf(stderr, "\tNon-matrix stuff:\n");
	fprintf(stderr, "\t\t-nomatrix (skip matrix creation - use if more than 100,000 loci)\n");
	//fprintf(stderr, "\t\t-loci <filename> (calculate collective p-value of each loci)\n");
	fprintf(stderr, "\t\t-interactions <filename> (output interactions, add \"-center\" to adjust pos to avg of reads)\n");
	fprintf(stderr, "\t\t-pvalue <#> (p-value cutoff for interactions, default 0.001)\n");
	fprintf(stderr, "\t\t-zscore <#> (z-score cutoff for interactions, default 1.0)\n");
	fprintf(stderr, "\t\t-minDist <#> (Minimum interaction distance, default: resolution/2)\n");
	fprintf(stderr, "\t\t-maxDist <#> (Maximum interaction distance, default: -1=none)\n");
	fprintf(stderr, "\t\t-boundary <#> (score boundaries at the given scale i.e. 100,000)\n");

	fprintf(stderr, "\n\tMaking Histograms:\n");
	fprintf(stderr, "\t\t-hist <#> (create a histogram matrix around peak positions, # is the resolution)\n");
	fprintf(stderr, "\t\t-size <#> (size of region in histogram, default = 100 * resolution)\n");

	fprintf(stderr, "\n\tComparing HiC experiments:\n");
	fprintf(stderr, "\t\t-ped <background PE tag directory>\n");

	fprintf(stderr, "\n\tCreating BED file to view with Wash U Epigenome Browser:\n");
	fprintf(stderr, "\t\t-washu (Both matrix and interaction outputs will be in WashH BED format)\n");

	fprintf(stderr, "\n\tCreating Circos Diagrams:\n");
	fprintf(stderr, "\t\t-circos <filename prefix> (creates circos files with the given prefix)\n");
	fprintf(stderr, "\t\t-d <tag directory 1> [tag directory 2] ... (will plot tag densities with circos)\n");
	fprintf(stderr, "\t\t-b <peak/BED file> (similar to visiualization of genes/-g, but no labels)\n");
	fprintf(stderr, "\t\t-g <gene location file> (shows gene locations)\n");

	fprintf(stderr, "\n\tGiven Interaction Analysis Mode (no matrix is produced):\n");
	fprintf(stderr, "\t\t-i <interaction input file> (for analyzing specific sets of interactions)\n");
	fprintf(stderr, "\t\t-iraw <output BED filename> (output raw reads from interactions, or -irawtags <file>)\n");
	fprintf(stderr, "\t\t-4C <output BED file> (outputs tags interacting with specified regions)\n");
	fprintf(stderr, "\t\t-peakStats <output BED file prefix> (outputs several UCSC bed/bedGraph files with stats)\n");
	fprintf(stderr, "\t\t-cpu <#> (number of CPUs to use, default: 1)\n");
	
	fprintf(stderr, "\n");
	exit(0);
}

