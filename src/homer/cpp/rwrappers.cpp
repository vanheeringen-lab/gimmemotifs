
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


//just a bunch of routines to help graph if R is available on the system...

#include "rwrappers.h"

RWrapper::RWrapper() {
	rstatus = RSTATUS_UNKNOWN;
	colors = NULL;
	checkForR();
	initColors(5);
}
RWrapper::~RWrapper() {
	deleteColors();
}

void RWrapper::deleteColors() {
	if (colors != NULL) {
		for (int i=0;i<numColors;i++) {
			if (colors[i] != NULL) delete [](colors[i]);
		}
		delete []colors;
	}
	colors = NULL;
	numColors = 0;
}

void RWrapper::initColors(int n) {
	if (colors != NULL) {
		deleteColors();
	}
	numColors = n;
	if (n == 0) return;
	srand(time(NULL));
	colors = new char*[numColors];
	for (int i=0;i<n;i++) {
		colors[i] = new char[20];
		if (i==0) {
			sprintf(colors[i],"\"blue\"");
		} else if (i==1) {
			sprintf(colors[i],"\"red\"");
		} else if (i==2) {
			sprintf(colors[i],"\"lightgreen\"");
		} else if (i==3) {
			sprintf(colors[i],"\"purple\"");
		} else if (i==4) {
			sprintf(colors[i],"\"black\"");
		} else {
			sprintf(colors[i],"c(%d,%d,%d)", rand() % 256, rand() % 256, rand() % 256);
		}
	}
}
int RWrapper::checkForR() {
    if (rstatus == RSTATUS_UNKNOWN) {
        char command[1000];
        sprintf(command, "R --version");
        int ok = system(command);
        if (ok == 0) {
            rstatus = RSTATUS_GOODTOGO;
        } else {
            rstatus = RSTATUS_NOTFOUND;
        }
    }
    fprintf(stderr, "rstatus=%d\n", rstatus);
    return rstatus;
}

void RWrapper::linePlot(char* inputFile, char* outputFile, int* cols, int numCols,
				char* title, char* xlabel, char* ylabel) {

	if (numCols > 5) initColors(numCols);
		
	char* filename = new char[1000];
	char* command = new char[1000];
	sprintf(filename, "%s.script.R",outputFile);
	FILE* fp = fopen(filename, "w");
	if (fp == NULL) {
		fprintf(stderr, "Can't open a temporary file %s\n", filename);
		return;
	}
	fprintf(fp, "png(%s)\n", outputFile);
	fprintf(fp, "x <- read.delim(\"%s\")\n", inputFile);
	fprintf(fp, "plot(x[,1],x[,2],type=\"n\")\n");
	for (int i=0;i<numCols;i++) {
		fprintf(fp, "lines(x[,1],x[,%d],col=%s)\n", cols[i],colors[i]);
	}
	fprintf(fp, "dev.off()\n");


	fclose(fp);
		

	sprintf(command, "R --no-save < \"%s\"", filename);
	int ok = system(command);
//	sprintf(command, "rm \"%s\"", filename);
//	ok = system(command);

	delete []filename;
}




