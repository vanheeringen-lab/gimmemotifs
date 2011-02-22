#ifndef GLAM2_SITE_SAMPLE_H
#define GLAM2_SITE_SAMPLE_H

#include "glam2_glam2.h"

void unalign(glam2_aln *aln, const int seq_pick, const fasta *f);
void realign(glam2_aln *aln, const int seq_pick, const fasta *f);

void site_sample(glam2_aln *aln, int seq_pick, data *d, double temperature);

#endif
