
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
#include <stdlib.h>
#include <string.h>
#include "Hashtable.h"

#define MAXIMUM_ALLOWED 1
#define MAXLINE 100000
#define MAXTABS 10000
#define ARRAYSIZE 10000

class Peak {
public: 
	int start;
	int end;
	int dir;
	int finished;
	char* name;
	char* chr;
	Peak* next;
	Peak();
	Peak(char* n, char* c, int s, int e, int d);
	~Peak();
};

class ChrPeaks {
public:
	Peak** peaks;
	int numPeaks;
	int maxPeaks;
	void addPeak(Peak*);
	void sortPeaks();

	ChrPeaks();
	~ChrPeaks();
};

int checkInt(char* str);
int comparePeaks(const void* a, const void* b);
void cleanUpSeq(char* seq);
void revopp(char* seq);
void revopp2seq(char* revoppseq, char* seq);
void readSeq(int mode, Hashtable* oligos, int &fileIndex, char** argv, int argc, 
					int &curPosition, int &finished, int checkSize, int readSize);
void reportUniqueSequences(Hashtable* oligos, char** argv);

struct Oligo {
public:
	unsigned short c; //file index
	unsigned int p;
	unsigned int v;
	char d;//direction
};
	

#define READ_SEQ 0
#define CHECK_SEQ 1

int main (int argc, char** argv) {

	if (argc < 4) {
		fprintf(stderr, "\n\tUsage: getMappableRegions <# seq at a time|100000000> <fragment length> <FASTA file1> [FASTA file2...]\n");
		fprintf(stderr, "\n\tWill return tag file of uniquely mappable positions to STDOUT\n");
		fprintf(stderr, "\n\tThis program can take a long time to finish, and uses a bunch of memory.\n");
		fprintf(stderr, "\tIt works (naively) by loading # sequences at a time into a hashtable, then\n");
		fprintf(stderr, "\tscreening the rest of the genome for matches (both strands).  But it works :)\n");
		fprintf(stderr, "\tFor Homer to use the output, use \"homerTools special uniqmap\".\n");
		fprintf(stderr, "\n");
		return(0);
	}


	char* curfile = new char[1000];
	curfile[0] = '\0';
	int curPosition = 0;
	int curFileIndex = 3;

	int checkSize = 0;
	sscanf(argv[1],"%d",&checkSize);
	int readSize = 0;
	sscanf(argv[2],"%d",&readSize);

	int oligoHashSize = checkSize*2;

	Hashtable* oligos = new Hashtable(oligoHashSize);

	int finished = 0;
	int mode = READ_SEQ;

	while (!finished) {
		mode = READ_SEQ;
		readSeq(mode, oligos, curFileIndex, argv, argc, curPosition, finished, checkSize, readSize);

		mode = CHECK_SEQ;
		int checkFinished = 0;
		int checkPosition = 0;
		int checkFileIndex = 3;
		readSeq(mode, oligos, checkFileIndex, argv, argc, checkPosition, checkFinished, checkSize, readSize);
		reportUniqueSequences(oligos,argv);
		if (finished) {
			break;
		}
		delete oligos;
		oligos = new Hashtable(oligoHashSize);
	}
	return 0;
}


void readSeq(int mode, Hashtable* oligos, int &fileIndex, char** argv, int argc, 
					int &curPosition, int &finished, int checkSize, int readSize) {

	char* buffer = new char[MAXLINE+1];	

	char* seq1 = new char[100000];
	char* seq2 = new char[100000];
	char* revoppSeq = new char[1000000];
	char* curSeq = seq1;
	char* nextSeq = seq2;
	char* tmpSeq = NULL;
	char* seq = new char[100000];
	if (mode == READ_SEQ) {
		fprintf(stderr, "\tReading Sequences to the HASHTABLE\n");
	} else if (mode == CHECK_SEQ) {
		fprintf(stderr, "\tChecking for Redundant Sequences\n");
	}
	int totalOligos = 0;

	for (;fileIndex<argc;fileIndex++) {
		char* seqfile = argv[fileIndex];
		fprintf(stderr, "\tLooking through sequences from %s starting at %d\n", seqfile, curPosition);

		int reportSize = 10000000;
		int lastSizeReported = 0;


		FILE* fp = fopen(seqfile, "r");
		if (fp == NULL) {
			fprintf(stderr, "\t!!Could not open file for %s (.fa or .fa.masked)\n",seqfile);
			continue;
		}

		int totalLen=0;
		curSeq[0] = '\0';
		int curStart = 0;

		int curPos = 0;

		while (fgets(buffer,MAXLINE,fp)!=NULL) {
			int lineLen = strlen(buffer);
			if (lineLen > 0 && buffer[lineLen-1] == '\n') {
				buffer[lineLen-1] = '\0';
			}
			if (buffer[0] == '>') continue;
			cleanUpSeq(buffer);
			int L = strlen(buffer);
			totalLen += L;

			if (curSeq[0] == '\0') {
				curStart = totalLen-L;
			}
			strcat(curSeq, buffer);
//fprintf(stderr, "%d,",curPos);
			if (totalLen > lastSizeReported+reportSize) {
				lastSizeReported+=reportSize;
				fprintf(stderr , "\t\t%d\n",lastSizeReported);
			}
	
			revopp2seq(revoppSeq,curSeq);

			while (curPos+readSize < totalLen) {
				if (curPosition > curPos) {
					curPos++;
					continue;
				}
				int index = curPos-curStart;
				char* oligo = &(curSeq[index]);
				char nullCharacter = curSeq[index+readSize];
				curSeq[index+readSize] = '\0';
				if (mode == READ_SEQ) {
					Oligo* o = (Oligo*) oligos->search(oligo);
					if (o == NULL) {
						o = new Oligo();
						o->c=(unsigned short)fileIndex;
						o->p=curPos+1;
						o->v=0;
						o->d=0;
						oligos->insert(o,oligo);
					}
				} else if (mode == CHECK_SEQ) {
					Oligo* o = (Oligo*)oligos->search(oligo);
					if (o != NULL) {
						o->v++;
						if (o->v > MAXIMUM_ALLOWED) {
							o = (Oligo*)oligos->remove(oligo);
							delete o;
							o = NULL;
						}
					}
				}
				totalOligos++;
				curSeq[index+readSize] = nullCharacter;

				int revoppCurPos = (curPos+readSize)-1;
				index = totalLen-(curPos+readSize)-1;
				oligo = &(revoppSeq[index]);
				nullCharacter = revoppSeq[index+readSize];
				revoppSeq[index+readSize] = '\0';
				if (mode == READ_SEQ) {
					Oligo* o = (Oligo*) oligos->search(oligo);
					if (o == NULL) {
						o = new Oligo();
						o->c=(unsigned short)fileIndex;
						o->p=revoppCurPos+2;
						o->v=0;
						o->d=1;
						oligos->insert(o,oligo);
					}
				} else if (mode == CHECK_SEQ) {
					Oligo* o = (Oligo*)oligos->search(oligo);
					if (o != NULL) {
						o->v++;
						if (o->v > MAXIMUM_ALLOWED) {
							o = (Oligo*)oligos->remove(oligo);
							delete o;
							o = NULL;
						}
					}
				}
				revoppSeq[index+readSize] = nullCharacter;
				totalOligos++;
				
				curPos++;
			}

			if (curPosition < curPos) {
				curPosition = curPos;
			}
			int lastKeeper = curPos-curStart;
			strcpy(nextSeq, &(curSeq[lastKeeper]));
			tmpSeq = curSeq;
			curSeq = nextSeq;
			nextSeq = tmpSeq;
			curStart += lastKeeper;

			if (mode == READ_SEQ && totalOligos > checkSize) {
				break;
			}
		}
		if (feof(fp)) {
			curPosition = 0;
			if (fileIndex == argc-1) {
				finished=1;
			}
		}
		fclose(fp);
		if (mode == READ_SEQ && totalOligos > checkSize) {
			break;
		}

	}
	fprintf(stderr, "\tTotal of %d oligos",totalOligos);
	if (mode == READ_SEQ) {
		fprintf(stderr, " read\n");
	} else {
		fprintf(stderr, " checked\n");
	}

	delete []seq1;
	delete []seq2;
	delete []seq;
	delete []buffer;
}


void reportUniqueSequences(Hashtable* oligos, char** argv) {
	
	char** keys = oligos->keys();	
	int numkeys = oligos->total;
	for (int i=0;i<numkeys;i++) {
		Oligo* oligo = (Oligo*)oligos->search(keys[i]);
		if (oligo != NULL) {
			fprintf(stdout, "%s\t%d\t%d\n",argv[oligo->c],oligo->p,oligo->d);
			delete oligo;
		}	
		delete [](keys[i]);
	}
	delete []keys;
}
	

int checkInt(char* str) {
	int bad = 0;
	int len = strlen(str);
	for (int i=0;i<len;i++) {
		if (str[i] < 45) bad=1;
		if (str[i] == 46 || str[i] == 47) bad=1;
		if (str[i] > 57) bad=1;
	}
	return bad;
}
void cleanUpSeq(char* seq) {
	if (seq == NULL) return;
	int len = strlen(seq);
	for (int i=0;i<len;i++) {
		if (seq[i] == 'a') {
			seq[i] = 'A';
		} else if (seq[i] == 'c') {
			seq[i] = 'C';
		} else if (seq[i] == 'g') {
			seq[i] = 'G';
		} else if (seq[i] == 't') {
			seq[i] = 'T';
		} else if (seq[i] == 'n') {
			seq[i] = 'N';
		}
	}
}
void revopp2seq(char* revoppseq, char* seq) {
	if (seq == NULL) return;
	int len = strlen(seq);
	revoppseq[len]='\0';
	for (int i=0;i<len;i++) {
		if (seq[i] == 'A') {
			revoppseq[len-1-i] = 'T';
		} else if (seq[i] == 'C') {
			revoppseq[len-1-i] = 'G';
		} else if (seq[i] == 'G') {
			revoppseq[len-1-i] = 'C';
		} else if (seq[i] == 'T') {
			revoppseq[len-1-i] = 'A';
		} else if (seq[i] == 'N') {
			revoppseq[len-1-i] = 'N';
		} else {
			revoppseq[len-1-i] = seq[i];
			//fprintf(stderr, "Non-sequence character: %c\n",newseq[i]);
		}
	}
}
	
void revopp(char* seq) {
	if (seq == NULL) return;
	int len = strlen(seq);
	char* newseq = new char[len+1];
	strcpy(newseq, seq);
	for (int i=0;i<len;i++) {
		if (newseq[i] == 'A') {
			seq[len-1-i] = 'T';
		} else if (newseq[i] == 'C') {
			seq[len-1-i] = 'G';
		} else if (newseq[i] == 'G') {
			seq[len-1-i] = 'C';
		} else if (newseq[i] == 'T') {
			seq[len-1-i] = 'A';
		} else if (newseq[i] == 'N') {
			seq[len-1-i] = 'N';
		} else {
			seq[len-1-i] = newseq[i];
			//fprintf(stderr, "Non-sequence character: %c\n",newseq[i]);
		}
	}
	delete []newseq;
}
	

Peak::Peak() {
	start = -1;
	end = -1;
	dir = 0;
	finished = 0;
	name= NULL;
	chr = NULL;
}
Peak::Peak(char *n, char* c, int s, int e, int d) {
	start = s;
	end = e;
	dir = d;
	finished = 0;
	name = new char[strlen(n)+1];
	strcpy(name,n);
	chr = new char[strlen(c)+1];
	strcpy(chr,c);
}
Peak::~Peak() {
	if (name != NULL) {
		delete name;
	}
	if (chr != NULL) {
		delete chr;
	}
}

ChrPeaks::ChrPeaks() {
	peaks=NULL;
	numPeaks=0;
	maxPeaks=0;
}
ChrPeaks::~ChrPeaks() {
	if (peaks != NULL) {
		for (int i=0;i<numPeaks;i++) {
			delete peaks[i];
		}
		delete []peaks;
	}
	peaks=NULL;
	numPeaks=0;
}
void ChrPeaks::addPeak(Peak* peak) {
	if (peaks == NULL) {
		peaks = new Peak*[ARRAYSIZE];
		maxPeaks = ARRAYSIZE;
	} else if (numPeaks+1 >= maxPeaks) {
		Peak** newPeaks = new Peak*[ARRAYSIZE+maxPeaks];
		maxPeaks = ARRAYSIZE+maxPeaks;
		for (int i=0;i<numPeaks;i++){ 
			newPeaks[i] = peaks[i];
		}
		delete []peaks;
		peaks = newPeaks;
	}
	peaks[numPeaks] = peak;
	numPeaks++;
}
void ChrPeaks::sortPeaks() {
	qsort(peaks,numPeaks,sizeof(Peak*),comparePeaks);
}

int comparePeaks(const void* a, const void* b) {
	int s1 = (*(Peak**)a)->start;
	int s2 = (*(Peak**)b)->start;

	if (s1 < s2) {
		return -1;
	} else if (s1 > s2) {	
		return 1;
	} else {
		return 0;
	}
}
