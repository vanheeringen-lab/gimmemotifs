#ifndef __ALLOC_H__
#define __ALLOC_H__

char *alloc_char(int size);
char **alloc_char_char(int size1, int size2);

float *alloc_float(int size);
float **alloc_float_float(int size1, int size2);
float ***alloc_float_float_float(int size1, int size2, int size3);

int *alloc_int(int size);
int **alloc_int_int(int size1, int size2);

double *alloc_double(int size);
double **alloc_double_double(int size1, int size2);
double ***alloc_double_double_double(int size1, int size2, int size3);

#endif
