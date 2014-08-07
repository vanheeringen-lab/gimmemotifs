

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

#include "Motif2.h"

int main(int argc, char** argv) {
	if (argc < 2) {
		printCMDhomer();
	}
	if (strcmp(argv[1],"denovo")==0) {
		programDeNovo(argc,argv);
	} else if (strcmp(argv[1],"known")==0) {
		programKnown(argc,argv);
	} else if (strcmp(argv[1],"find")==0) {
		programKnown(argc,argv);
	} else if (strcmp(argv[1],"norm")==0) {
		programKnown(argc,argv);
	} else if (strcmp(argv[1],"mask")==0) {
		programKnown(argc,argv);
	} else {
		fprintf(stderr, "!!! homer2 command \"%s\" not recognized !!!\n",argv[1]);
		printCMDhomer();
	}
	return 0;
}
