#include "mergesort.h"


void mergeSorter(void* input, int left, int right, int size,
				int (*compare)(const void*,const void*), void* scratch) {
	if (right == left+1) return;

	int length = right-left;
	int midpoint_distance = length/2;
	int l=left,r=left+midpoint_distance;

	mergeSorter(input,left,left+midpoint_distance,size,compare,scratch);
	mergeSorter(input,left+midpoint_distance,right,size,compare,scratch);
	int same = 0;
	int bad = 0;
	for (int i=0;i<length;i++) {
		if (l < left+midpoint_distance && 
				(r == right || 0 > compare((void*)((unsigned int)input+size*l),(void*)((unsigned int)input+size*r)))) {
			memcpy((void*)((unsigned int)scratch+size*i),(void*)((unsigned int)input+size*l),size);
			l++;
		} else {
			memcpy((void*)((unsigned int)scratch+size*i),(void*)((unsigned int)input+size*r),size);
			r++;
		}
	}
	for (int i=left;i<right;i++){ 
		memcpy((void*)((unsigned int)input+size*i),(void*)((unsigned int)scratch+size*(i-left)),size);
	}
}

void mergesort(void* input, int nelements, int size, 
				int (*compare)(const void*,const void*)) {
	void* scratch = (void*) malloc(size*nelements);
	mergeSorter(input,0,nelements,size,compare,scratch);
	free(scratch);
}

