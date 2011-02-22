#include <stdio.h>

#include <stdlib.h> /* added for exit */

#include "getfa.h"

int nword = 16;
int nk = 3;

void dust();

void dust_rect() {}

void set_dust_level(int value);
void getfafun(char *name, void (*fun)());

/*void main(argc, argv)
int argc;
char *argv[];
{
	FASTA *fa;
	int level = 11;

	if (argc < 2) {
		fprintf(stderr, "Usage: dust fasta-file [ cut-off ]\n");
		exit(1);
	}
	fa = getfa(argv[1], perform_wo64);
	if (argc > 2) {
		level = atoi(argv[2]);
	}
	printf("Sequence length: %d\n", fa->len);
	perform_wo64(fa->len, fa->seq, level);
}*/

int main(argc, argv)
int argc;
char *argv[];
{
	FASTA *fa;
	REGION *reg;
	int level = 20;
	int i;

	if (argc < 2) {
		fprintf(stderr, "Usage: dust fasta-file [ cut-off ]\n");
		exit(1);
	}
	if (argc > 2) {
		level = atoi(argv[2]);
	}
	set_dust_level(level);
	if (argc >= 4) {
		fa = getfa(argv[1]);
		reg = dust_segs(fa->len, fa->seq);
		for (i=0; reg[i].to != -1; i++) {
			printf("%6d..%d\n", reg[i].from, reg[i].to);
		}
	} else {
		getfafun(argv[1], dust);
	}
	return(0);
}
