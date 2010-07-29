#ifndef GLAM2_OUTPUT_H
#define GLAM2_OUTPUT_H

#include <stdio.h>
#include "glam2_glam2.h"

void print_aln(FILE *fp, glam2_aln *aln, data *d);

void print_aln_info(FILE *fp, glam2_aln *aln, data *d);

void print_alns(FILE *fp, glam2_aln *alns, data *d);

#endif
