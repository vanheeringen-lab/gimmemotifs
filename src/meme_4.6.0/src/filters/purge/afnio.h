#include <stdio.h>
#include "stdinc.h"

FILE	*open_file(char *fstring,char *subfile,char *cmnd);
long     ParseIntegers(char *str, long *values, char *msg);
long     ParseReals(char *str, double *values, char *msg);
char	*String(char *s);

