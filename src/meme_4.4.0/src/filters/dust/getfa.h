typedef struct {
	int len;
	char *header;
	char *seq;
} FASTA;

FASTA *getfa();
void *mymalloc();

#define MAXREG    1001

typedef struct {
	int from;
	int to;
	int score;
} REGION;

REGION *dust_segs();



