#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <cmath>
#include <time.h>
#include <limits.h>

#include "Hashtable.h"
#include "statistics.h"
#include "SeqTag.h"
#include "Clustering.h"

#ifndef HOMERTOOLS_H
#define HOMERTOOLS_H

#define BARCODE_MINFREQ 0.02
#define TRIM_MAXREADLENGTH 100000

#define BUFFER 1000000



class StrSort {
public:
	char* str;
	double v;
};
int decendStrSort(const void* a, const void* b);

void split(char* string, char** cols, int &numCols, char delim);

#endif
