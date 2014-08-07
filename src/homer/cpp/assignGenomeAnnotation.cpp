
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

int main(int argc, char** argv) {

	char* peakfile = NULL;
	char* annfile = NULL;
	int directoryFlag = 0;
	char* annOutputFile = NULL;
	char* statsOutputFile = NULL;
	char* prioritizeFile = NULL;


	int tagAdjust = 0;
	int prioritizingIters = 1;

	
	if (argc < 3) {
		printCMD();
	}
	peakfile = argv[1];
	annfile = argv[2];
	for (int i=1;i<argc;i++) {
		if (i==1 || i==2) {
			if (argv[i][0] == '-') {
				fprintf(stderr, "!!! First and 2nd arguments are <peak file> and <annotation file>\n");
				printCMD();
			}
			continue;
		}
		if (argv[i][0] == '-') {
			if (strcmp(argv[i],"-tags")==0) {
				directoryFlag = 1;
			} else if (strcmp(argv[i],"-prioritize")==0) {
				prioritizeFile = argv[++i];
			} else if (strcmp(argv[i],"-ann")==0) {
				annOutputFile = argv[++i];
			} else if (strcmp(argv[i],"-np")==0) {
				sscanf(argv[++i],"%d",&prioritizingIters);
			} else if (strcmp(argv[i],"-stats")==0) {
				statsOutputFile = argv[++i];
			} else if (strcmp(argv[i],"-fragLength")==0) {
				i++;
				if (strcmp(argv[i],"auto")==0) {
					tagAdjust = TAGADJUST_AUTO;
				} else {
					sscanf(argv[i], "%d",&tagAdjust);
					tagAdjust = tagAdjust/2;
				}
			} else {
				fprintf(stderr, "!!! Can't recognize \"%s\"\n", argv[i]);
				printCMD();
			} 
		}
	}


	PeakLibrary* rawAnnotations = new PeakLibrary(30000000);
	rawAnnotations->readPeakFile(annfile,PEAK_READ_MODE_ANNOTATION);
	PeakLibrary* annotations = NULL;
	//system("sleep 100");

	if (prioritizeFile != NULL) {
		PeakLibrary* tmp = rawAnnotations;
		for (int i=0;i<prioritizingIters;i++) {
			annotations = tmp->prioritizeAnnotations();
			delete tmp;
			tmp = annotations;
		}
		//annotations->printSorted(stdout);
		FILE* fpann = fopen(prioritizeFile,"w");
		if (fpann == NULL) {
			fprintf(stderr, "!!! Could not open annotation file: %s for writing!!!\n", prioritizeFile);
		} else {
			annotations->printAnnotation(fpann);
			fclose(fpann);
		}
	} else {
		annotations = rawAnnotations;
	}

	FILE* annOutput = NULL;
	if (annOutputFile != NULL) {
		annOutput = fopen(annOutputFile, "w");
		if (annOutput == NULL) {
			fprintf(stderr, "!! Could not open %s for writing annotations!\n", annOutputFile);
		}
	}
	FILE* statOutput = NULL;
	if (statsOutputFile != NULL) {
		statOutput = fopen(statsOutputFile, "w");
		if (statOutput == NULL) {
			fprintf(stderr, "!! Could not open %s for annotation stats!\n", statsOutputFile);
		}
	}
	if (statOutput == NULL) {
		statOutput = stdout;
	}

	if (directoryFlag == 0) {

		PeakLibrary* peaks = new PeakLibrary();
		peaks->readPeakFile(peakfile, PEAK_READ_MODE_NORMAL);
		peaks->annotatePeakLocations(annotations, statOutput, annOutput);

	} else {

		TagLibrary* tags = new TagLibrary(peakfile);
		tags->readTagDirectory();
		tags->setSingleRead(1);
		tags->setTagAdjust(tagAdjust);
		tags->annotateTagLocations(annotations, statOutput, annOutput);

	}

	if (annOutput != NULL) {
		fclose(annOutput);
	}
	if (statOutput != NULL && statOutput != stdout) {
		fclose(statOutput);
	}

	return 0;

}


void printCMD() {
	fprintf(stderr, "\n\tUsage: assignGenomeAnnotation <peak file or tag directory> <annotation file> [options]\n");
	fprintf(stderr, "\n\tAssigns peaks or tags to genomic annotations based on their location\n"); 
	fprintf(stderr, "\n\tGeneral Options:\n");
	fprintf(stderr, "\t\t-tags (Use if first argument is a tag directory)\n");
	fprintf(stderr, "\t\t-fragLength <#|auto> (Approximate fragment length, default: auto)\n");
	fprintf(stderr, "\t\t-prioritize <outputfile> (annotation file is just a peak file - hasn't been prioritized yet)\n");
	//fprintf(stderr, "\t\t\t-np <#> (number prioritizing iteractions...)\n");
	fprintf(stderr, "\t\t\tA prioritized file will be created for future use with the program\n");
	fprintf(stderr, "\t\t\tThis option should be used if the annotation file isn't prioritized\n");
	fprintf(stderr, "\n\tOutput Options:\n");
	fprintf(stderr, "\t\t-ann <filename> (File to output annotations for each peak/tag, by default not created)\n");
	fprintf(stderr, "\t\t-stats <filename> (File to output annotation statistics, default to stdout)\n");
	fprintf(stderr, "\n");
	exit(0);
			

}

