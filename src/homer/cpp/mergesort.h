#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void mergesort(void* input, int nelements, int size,
                int (*compare)(const void*,const void*));
