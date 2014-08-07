
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

#include "SeqTag.h"

#define BUFFER 10000

void printCMD();
char** addPeakFile2List(char** pfiles, int &numFiles, char* file);

int main(int argc, char** argv) {

	char** peakFiles = NULL;
	int numPeakFiles = 0;
	char* prefix = NULL;
	char* vennFile = NULL;
	char* matrixFile = NULL;
	int codeFlag = 0;

	char* chrFilter = NULL;
	int startFilter = -1000000000;
	int endFilter = 1000000000;
	char* coverageFile = NULL;

	char strand = STRAND_BOTH;
	int maxDistance = 100;
	int maxCoBound = 0;
	long int gsize = DEFAULT_GSIZE;
	//int bpFlag = 0;
	int givenFlag = 1;

	if (argc < 2) {
		printCMD();
	}
	for (int i=1;i<argc;i++) {
		if (argv[i][0] == '-') {
			if (strcmp(argv[i],"-matrix")==0) {
				matrixFile = argv[++i];
			} else if (strcmp(argv[i],"-prefix")==0) {
				prefix = argv[++i];
			} else if (strcmp(argv[i],"-code")==0) {
				codeFlag = 1;
			} else if (strcmp(argv[i],"-coverage")==0) {
				coverageFile = argv[++i];
			} else if (strcmp(argv[i],"-filter")==0) {
				parseUCSCpositionStr(argv[++i],chrFilter,startFilter,endFilter);
			} else if (strcmp(argv[i],"-venn")==0) {
				vennFile = argv[++i];
			} else if (strcmp(argv[i],"-strand")==0) {
				strand = STRAND_SEPARATE;
			} else if (strcmp(argv[i],"-gsize")==0) {
				sscanf(argv[++i], "%ld",&gsize);
			} else if (strcmp(argv[i],"-cobound")==0) {
				sscanf(argv[++i], "%d",&maxCoBound);
			} else if (strcmp(argv[i],"-bp")==0) {
				//bpFlag = 1;
			} else if (strcmp(argv[i],"-d")==0) {
				i++;
				if (strcmp(argv[i],"given")==0) {
					givenFlag = 1;
				} else {
					sscanf(argv[i], "%d",&maxDistance);
					if (maxDistance < 1) maxDistance = 1;
					givenFlag = 0;
				}
			} else if (strcmp(argv[i],"-file")==0) {
				i++;
				FILE* fp = fopen(argv[i],"r");
				if (fp==NULL) {
					fprintf(stderr, "!!! Can't open [-file] %s\n", argv[i]);
					exit(0);
				}
				char* buf = new char[BUFFER];
				while (fgets(buf, BUFFER, fp) != NULL) {
					if (buf[0] == '#') continue;
					char* filename = new char[strlen(buf)+1];
					strcpy(filename, buf);
					peakFiles = addPeakFile2List(peakFiles, numPeakFiles, filename);
				}
				delete []buf;
				fclose(fp);
			} else {
				fprintf(stderr, "!!! Can't recognize \"%s\"\n", argv[i]);
				printCMD();
			} 
		} else {
			peakFiles = addPeakFile2List(peakFiles, numPeakFiles, argv[i]);
		}
	}

	char* cmdline = getCMDstr(argc,argv);

	if (numPeakFiles < 1) {
		fprintf(stderr, "!!! Can't do much without a  peak file !!!\n");
		exit(0);
	}

	PeakLibrary* refPeaks = new PeakLibrary(10000000);
	fprintf(stderr, "\tMax distance to merge: ");
	if (givenFlag) {
		fprintf(stderr, "direct overlap required (-d given)\n");
	} else {
		fprintf(stderr, "%d bp\n", maxDistance);
	}
	if (givenFlag == 1) {
		maxDistance = -1;
	}

	if (maxCoBound == 0) {
		if (numPeakFiles < 2) {
			fprintf(stderr, "\tMerging single peak file... \n");
			refPeaks->readPeakFile(peakFiles[0],PEAK_READ_MODE_NORMAL);
			PeakLibrary* newpeaks = refPeaks->mergeNearbyPeaks(maxDistance,strand);
			if (chrFilter != NULL) {
				PeakLibrary* filteredPeaks = newpeaks->filterPeaksOutsideRange(chrFilter, 
															startFilter,endFilter);
				delete newpeaks;
				newpeaks = filteredPeaks;
			}
			if (coverageFile != NULL) {
				long long int totalCoverage = newpeaks->calculateCoverage();
				fprintf(stderr, "\tTotal coverage: %lld bp\n", totalCoverage);
				FILE* fp=fopen(coverageFile,"w");
				if (fp != NULL) {	
					fprintf(fp,"Peak file\tNumber of Peaks\tTotal Coverage(bp)\n");
					fprintf(fp,"%s\t%d\t%lld\n",peakFiles[0],newpeaks->numPeaks,totalCoverage);
					fclose(fp);
				} else {
					fprintf(stderr, "!!! Couldn't open coverage file (%s) !!!\n", coverageFile);
				}
			}
			fprintf(stdout, "# cmd = %s\n", cmdline);
			newpeaks->print(stdout);
		} else {
			fprintf(stderr, "\tMerging peaks... \n");
			refPeaks->mergePeaks(peakFiles, numPeakFiles, strand, maxDistance,
							gsize, prefix, vennFile, matrixFile,codeFlag,cmdline);
		}
	} else {
		if (numPeakFiles < 2) {
			fprintf(stderr, "!!! Can't look for co-bound peaks without 2 or more peak files !!!\n");
			exit(0);
		}
		fprintf(stderr, "\tCalculating co-bound peaks relative to reference: %s\n", peakFiles[0]);
		refPeaks->readPeakFile(peakFiles[0], PEAK_READ_MODE_NORMAL);
		char** pfiles = &(peakFiles[1]);
		refPeaks->getCoBoundPeaks(pfiles, numPeakFiles-1, strand, maxDistance, gsize, 
													prefix,maxCoBound,matrixFile,cmdline);
	}

	return 0;

}

char** addPeakFile2List(char** pfiles, int &numFiles, char* file) {
	int alreadyPresent = 0;	
	char** peakfiles = new char*[numFiles+1];
	if (numFiles > 0) {
		for (int i=0;i<numFiles;i++) {
			peakfiles[i] = pfiles[i];
			if (strcmp(peakfiles[i],file)==0) {
				alreadyPresent=1;
			}
		}
	}
	delete []pfiles;
	if (alreadyPresent == 0) {
		peakfiles[numFiles] = file;
		numFiles++;
	}
	return peakfiles;
}

void printCMD()  {
	fprintf(stderr, "\n\tUsage: mergePeaks [options] <primary peak file> [additional peak/annotation files...]\n");

	fprintf(stderr, "\n\tMerges and/or compares peak/position files (peak files listed twice are only considered once)\n"); 
	fprintf(stderr, "\n\tGeneral Options:\n");
	fprintf(stderr, "\t\t-strand (Only merge/consider peaks on the same strand, default: either strand)\n");
	fprintf(stderr, "\t\t-d <#|given> (Maximum distance between peak centers to merge, default: given)\n");
	fprintf(stderr, "\t\t\tUsing \"-d given\" looks for literal overlaps in peak regions - DEFAULT since v4.4\n");
	fprintf(stderr, "\t\t\tUse \"-d given\" when features have vastly different sizes (i.e. peaks vs. introns)\n");
	fprintf(stderr, "\t\t-file <filename> (file listing peak files to compare - for lots of peak files)\n");
	fprintf(stderr, "\t\t-gsize <#> (Genome size for significance calculations, default: 2e9)\n");

	fprintf(stderr, "\n\tMerging Peaks Options (default):\n");
	fprintf(stderr, "\t\t-prefix <filename> (Generates separate files for overlapping and unique peaks)\n");
	fprintf(stderr, "\t\t\tBy default all peaks are sent to stdout\n");
	fprintf(stderr, "\t\t-matrix <filename> (Generates files with pairwise comparison statistics)\n");
	fprintf(stderr, "\t\t\tfilename.logPvalue.matrix.txt - ln p-values for overlap, +values for divergence\n");
	fprintf(stderr, "\t\t\tfilename.logRatio.matrix.txt - ln ratio of observed/expected overlaps\n");
	fprintf(stderr, "\t\t\tfilename.count.matrix.txt - peak overlap counts\n");
	fprintf(stderr, "\t\t-venn <filename> (output venn diagram numbers to file, default: to stderr)\n");
	fprintf(stderr, "\t\t-code (report peak membership as binary instead of by file names)\n");

	fprintf(stderr, "\n\tClassify peaks by how many are co-bound by other peak files vs. reference(1st file)\n");
	fprintf(stderr, "\t\t-cobound <#> (Maximum number of co-bound peaks to consider)\n");
	fprintf(stderr, "\t\t\tWill output sets of peaks that are co-bound by various numbers of factors\n");
	fprintf(stderr, "\t\t\tto files coBoundBy0.txt, coBoundBy1.txt, coboundBy2.txt, ...\n");
	fprintf(stderr, "\t\t\tOr <prefix>.coBoundBy0.txt, <prefix>.coBoundBy1.txt, ...\n");
	fprintf(stderr, "\t\t-matrix <filename> (generates similar files to above with pairwise overlap statistics)\n");

	fprintf(stderr, "\n\tSingle peak file:\n");
	fprintf(stderr, "\t\t(If a single peak file is given, peaks within the maximum distance will be merged)\n");
	fprintf(stderr, "\t\t-filter chrN:XXX-YYY (only analyze peaks within range)\n");
	fprintf(stderr, "\t\t-coverage <output file> (returns the total bp covered by each peak file - use \"-d given\"\n");
	fprintf(stderr, "\n");
	exit(0);
			
}
