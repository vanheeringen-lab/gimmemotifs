// Copyright 2009, 2010, 2011, 2012, 2013 Christopher Benner <cbenner@gmail.com>
// // 
// // This file is part of HOMER
// //
// // HOMER is free software: you can redistribute it and/or modify
// // it under the terms of the GNU General Public License as published by
// // the Free Software Foundation, either version 3 of the License, or
// // (at your option) any later version.
// //
// // HOMER is distributed in the hope that it will be useful,
// // but WITHOUT ANY WARRANTY; without even the implied warranty of
// // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// // GNU General Public License for more details.

#include "SeqTag.h"

void printCMD();

int main(int argc, char** argv) {


	char** tagDirs = new char*[10000];
	int numTagDirs=0;
	int res = 1000;
	int peakOutFlag = 1;

	if (argc < 2) {
		printCMD();
	}

	for (int i=1;i<argc;i++) {
		if (strcmp(argv[i],"-res")==0) {
			sscanf(argv[++i],"%d",&res);
		} else if (strcmp(argv[i],"-pout")==0) {
			peakOutFlag = 1;
		} else if (strcmp(argv[i],"-d")==0) {
			i++;
			for (;i<argc;i++) {
				if (strncmp(argv[i],"-",1)==0) {
					i--;
					break;
				} else {
					tagDirs[numTagDirs] = new char[strlen(argv[i])+1];
					strcpy(tagDirs[numTagDirs],argv[i]);
					numTagDirs++;
				}
			}
		}
	}

	int superRes = res;
	int start = -1;
	int end = 2000000000;
	char* chr = NULL;

	if (numTagDirs < 1) {
		fprintf(stderr, "!!! Need at least one tag directory!\n");
		printCMD();
	}
	char* tagDirName = tagDirs[0];
	TagLibrary* tags = new TagLibrary(tagDirName);
	tags->readTagDirectory();
	//fprintf(stderr, "%lf\n", tags->totalTags);
	PeakLibrary* peaks = tags->getCoveragePeaks(chr,start,end,res,superRes);
	if (peakOutFlag) {
		peaks->print(stdout);
	}
	peaks->setPeakTagSizeFixed(0,0);
	//int countMode = COUNT_MODE_TOTAL;
}


void printCMD() {
	fprintf(stderr, "\n\tUsage: getGenomeTilingPeaks <options>\n");
	fprintf(stderr, "\tOutput sent to stdout\n");
	fprintf(stderr, "\n\tOptions:\n"); 
	fprintf(stderr, "\t\t-d <tag directory1> [tag directory 2] ... (tag directories to analyze)\n"); 
	fprintf(stderr, "\t\t-res <#> (resolution of analysis)\n"); 
	//fprintf(stderr, "\n\tOther functions:\n");
	//fprintf(stderr, "\t\t-pout (print out peak file of all regions in the genome)\n");
	fprintf(stderr, "\n");
	exit(0);
}
