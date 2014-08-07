
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

#ifndef CLUSTERING_H
#define CLUSTERING_H

class TreeCluster;
class TreeNode;

class TreeCluster {
public:

	double** matrix;
	int numMatrix;
	char** names;

	TreeNode* clusters;
	int* mask;

	int fixedClusterScore;
	int fixedPositionFlag;
	int reverseFlag;

	TreeCluster();
	TreeCluster(double ** m, int numMatrix);
	~TreeCluster();

	double* minDistCache;
	int* minIndexCache;
	
	void cluster();
	int* getOrder();
	int* getClusters(double threshold);
	void init();
	void init(double** m, int numMatrix);
	void setReverseFlag();
	void printCDT(char* prefix,char** names, double** m);
	void printClusters(FILE* fp,char**names,double thresh);
	double findMinimumPair(int &left, int &right);
	double initializeCache(int &left,int &right);
	double updateCache(int &left,int &right);
	void printClusterSizes(FILE* fp,int maxClusterSizeToReport);
	void scoreClusterFiles(char** files, int numFiles, FILE* fp);
	void loadGTR(char* gtrFile);
};



class TreeNode {
public:
	int left;
	int right;
	double distance;
	double weight;
	
	int matrixIndex;

	TreeNode();
	~TreeNode();
	void print(FILE*);
};



#endif
