
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
#include <stdlib.h>

#ifndef RWRAPPERS_H
#define RWRAPPERS_H

#define RSTATUS_UNKNOWN 0
#define RSTATUS_NOTFOUND 1
#define RSTATUS_GOODTOGO 2

class RWrapper {
	int rstatus;
	
	RWrapper();
	~RWrapper();
	int checkForR();
	void linePlot(char* inputFile, char* outputFile, int* cols, int numCols,
			char* title, char* xlabel, char* ylabel);

}




#endif
