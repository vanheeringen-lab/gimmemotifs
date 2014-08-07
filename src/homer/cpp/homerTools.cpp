
// Copyright 2009 - 2014 Christopher Benner <cbenner@salk.edu>
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

#include "homerTools.h"
#define _FILE_OFFSET_BITS 64

void printCMD();
void printCMDbarcodes();
void printCMDtrim();
void printCMDfreq();
void printCMDextract();
void printCMDspecial();
void barcodesProgram(int argc, char** argv);
void truSeqProgram(int argc, char** argv);
void trimProgram(int argc, char** argv);
void freqProgram(int argc, char** argv);
void extractProgram(int argc, char** argv);
void clusterProgram(int argc, char** argv);
void matrix3DProgram(int argc, char** argv);
void specialProgram(int argc, char** argv);
void decontaminateProgram(int argc, char** argv);

int main(int argc, char** argv) {

	if (argc < 2) {
		printCMD();
	}

	char* program = argv[1];
	if (strcmp(program,"barcodes")==0) {
		barcodesProgram(argc,argv);
	} else if (strcmp(program,"trim")==0) {
		trimProgram(argc,argv);
	} else if (strcmp(program,"truseq")==0) {
		truSeqProgram(argc,argv);
	} else if (strcmp(program,"freq")==0) {
		freqProgram(argc,argv);
	} else if (strcmp(program,"extract")==0) {
		extractProgram(argc,argv);
	} else if (strcmp(program,"decontaminate")==0) {
		decontaminateProgram(argc,argv);
	} else if (strcmp(program,"cluster")==0) {
		clusterProgram(argc,argv);
	} else if (strcmp(program,"matrix3D")==0) {
		matrix3DProgram(argc,argv);
	} else if (strcmp(program,"special")==0) {
		specialProgram(argc,argv);
	} else if (strcmp(program,"--help")==0) {
		printCMD();
	} else {
		fprintf(stderr, "!!! Could not recognize \"%s\" as a command !!!\n", argv[1]);
		printCMD();
	}
	return 0;
}

void printCMD() {
	fprintf(stderr, "\n\tUsage: homerTools <command> [--help | options]\n");
	fprintf(stderr, "\n\tCollection of tools for sequence manipulation\n");
	fprintf(stderr, "\n\tCommands: [type \"homerTools <command>\" to see individual command options]\n");
	fprintf(stderr, "\t\tbarcodes - separate FASTQ file by barcodes\n");
	fprintf(stderr, "\t\ttruseq - process truseq barcodes from unidentified indexes (illumina)\n");
	fprintf(stderr, "\t\ttrim - trim adapter sequences or fixed sizes from FASTQ files(also splits)\n");
	fprintf(stderr, "\t\tfreq - calculate position-dependent nucleotide/dinucleotide frequencies\n");
	fprintf(stderr, "\t\textract - extract specific sequences from FASTA file(s)\n");
	fprintf(stderr, "\t\tdecontaminate - remove bad tags from a contaminated tag directory\n");
	fprintf(stderr, "\t\tcluster - hierarchical clustering of a NxN distance matrix\n");
//	fprintf(stderr, "\t\tmatrix3D - convert distance matrix to 3D coordinates\n");
	fprintf(stderr, "\t\tspecial - specialized routines (i.e. only really useful for chuck)\n");
	fprintf(stderr, "\t\t\n"); 
	exit(0);
}


//------------------------------  trim -------------------------------------------------------------------


int matchAdapter5prime(char* seq, char* adapter,int maxMisMatches,int minMatchLength,int matchStart);
int matchAdapter3prime(char* seq, char* adapter,int maxMisMatches,int minMatchLength,int matchStart);
void printCMDtrim() {
	fprintf(stderr, "\n\tUsage: homerTools trim [options] <fastq file1> [file2] ...\n");
	fprintf(stderr, "\n\tProbably best to use only on option at a time...\n");
	fprintf(stderr, "\n\tOptions for command: trim\n");
	fprintf(stderr, "\t\t-3 <#|[ACGT]> (trim # bp or adapter sequence from 3' end of sequences)\n");
	fprintf(stderr, "\t\t-5 <#|[ACGT]> (trim # bp or adapter sequence from 5' end of sequences)\n");
	fprintf(stderr, "\t\t\t-mis <#> (Maximum allowed mismatches in adapter sequence, default: 0)\n");
	fprintf(stderr, "\t\t\t-minMatchLength <#> (minimum adapter sequence at edge to match, default: half adapter length)\n");
	fprintf(stderr, "\t\t\t-matchStart <#> (don't start searching for adapter until this position, default: 0)\n");
	fprintf(stderr, "\t\t-q <#> (Trim sequences once quality dips below threshold, default: none [range:0-40])\n");
	fprintf(stderr, "\t\t\t-qstart <#> (don't check quality until sequences are at least this long, default: 10)\n");
	fprintf(stderr, "\t\t\t-qwindow <#> (size of moving average to check for quality dropoff, default: 5)\n");
	fprintf(stderr, "\t\t-len <#> (Keep first # bp of sequence - i.e. make them the same length)\n");
	fprintf(stderr, "\t\t-stats <filename> (Output trimming statistics to filename, default: sent to stdout)\n");
	fprintf(stderr, "\t\t-min <#> (Minimum size of trimmed sequence to keep, default: 1)\n");
	fprintf(stderr, "\t\t-max <#> (Maximum read length, default: %d)\n", TRIM_MAXREADLENGTH);
	fprintf(stderr, "\t\t-suffix <filename suffix> (output is sent to InuptFileName.suffix, default: trimmed)\n");
	fprintf(stderr, "\t\t-lenSuffix <filename suffix> (length distribution is sent to InuptFileName.suffix, default: lengths)\n");
	//fprintf(stderr, "\t\t-split <#> (Split reads into two reads at bp #, output to trimmed1 and trimmed2)\n");
	//fprintf(stderr, "\t\t-revopp <#> (Return reverse opposite of read [if used with -split, only the 2nd\n");
	//fprintf(stderr, "\t\t\t\thalf of the read will be retuned as reverse opposite])\n");
	fprintf(stderr, "\n\tCommon Examples:\n");
	fprintf(stderr, "\t\tIllumina TruSeq:\n");
	fprintf(stderr, "\t\t\thomerTools trim -3 GATCGGAAGAGCACACGTCT -mis 1 -minMatchLength 4 -min 15\n");
	//fprintf(stderr, "\t\tIllumina TruSeq (PE - 2nd read):\n");
	//fprintf(stderr, "\t\t\thomerTools trim -3 AGATCGGAAGAGCGTCGTGT -mis 1 -minMatchLength 6\n");
	fprintf(stderr, "\t\tIllumina small RNA adapter:\n");
	fprintf(stderr, "\t\t\thomerTools trim -3 TCGTATGCCGTCTTCTGCTTGT -mis 1 -minMatchLength 4 -min 15\n");
	fprintf(stderr, "\t\t\n"); 
	exit(0);
}
void trimProgram(int argc, char** argv) {

	char* suffix = NULL;
	char* lenSuffix = NULL;
	char* statsFile = NULL;
	char** files = new char*[100000];	
	int numfiles = 0;
	if (argc < 3) {
		printCMDtrim();
	}
	int p3len = -1;
	int p5len = -1;
	int matchStart = 0;
	int maxMisMatches = 0;
	int minMatchLength = -1;
	int splitPos = -1;
	int minSizeToKeep = 1;
	int qthresh = -1;
	int qwindow = 5;
	int qstart = 10;
	
	int revoppFlag = 0;
	int maxReadLength = TRIM_MAXREADLENGTH;
	char* p3adapter = NULL;
	char* p5adapter = NULL;
	int fixedLength = -1;
	for (int i=2;i<argc;i++) {
		if (strcmp(argv[i],"-len")==0) {
			sscanf(argv[++i],"%d",&fixedLength);
			fprintf(stderr, "\tSequences will be output with a fixed length of %d\n",fixedLength);
		} else if (strcmp(argv[i],"-min")==0) {
			sscanf(argv[++i],"%d",&minSizeToKeep);
			fprintf(stderr, "\tMinimum length of sequence to keep %d\n",minSizeToKeep);
		} else if (strcmp(argv[i],"-revopp")==0) {
			revoppFlag = 1;
			fprintf(stderr, "\tReturning the reverse opposite of the read\n");
		} else if (strcmp(argv[i],"-q")==0) {
			sscanf(argv[++i],"%d",&qthresh);
			fprintf(stderr, "\tWill trim reads once quality scores dip below %d\n", qthresh);
		} else if (strcmp(argv[i],"-qwindow")==0) {
			sscanf(argv[++i],"%d",&qwindow);
			fprintf(stderr, "\tWill trim reads based on quality with a %d bp moving window\n", qwindow);
		} else if (strcmp(argv[i],"-qstart")==0) {
			sscanf(argv[++i],"%d",&qstart);
			fprintf(stderr, "\tWill trim reads based on quality starting at position %d\n", qstart);
		} else if (strcmp(argv[i],"-split")==0) {
			sscanf(argv[++i],"%d",&splitPos);
			fprintf(stderr, "\tSplitting reads at bp %d\n", splitPos);
		} else if (strcmp(argv[i],"-minMatchLength")==0) {
			sscanf(argv[++i],"%d",&minMatchLength);
			fprintf(stderr, "\tMinimum number of matches to adapter sequence at the edge: %d\n", minMatchLength);
		} else if (strcmp(argv[i],"-mis")==0) {
			sscanf(argv[++i],"%d",&maxMisMatches);
			fprintf(stderr, "\tMaximum number of mismatches in adapter sequence to consider %d\n", maxMisMatches);
		} else if (strcmp(argv[i],"-max")==0) {
			sscanf(argv[++i],"%d",&maxReadLength);
			fprintf(stderr, "\tMaximum length of sequence to keep %d\n", maxReadLength);
		} else if (strcmp(argv[i],"-matchStart")==0) {
			sscanf(argv[++i],"%d",&matchStart);
			fprintf(stderr, "\tWill only start matching adapter sequence after %d [for 3', before for 5']\n", matchStart);
		} else if (strcmp(argv[i],"-suffix")==0) {
			suffix=argv[++i];
			fprintf(stderr, "\tEach output file will be given the original filename with .%s as a suffix\n", suffix);
		} else if (strcmp(argv[i],"-lenSuffix")==0) {
			lenSuffix=argv[++i];
			fprintf(stderr, "\tRead length distribution files will be given the original filename with .%s as a suffix\n", lenSuffix);
		} else if (strcmp(argv[i],"-stats")==0) {
			statsFile=argv[++i];
			fprintf(stderr, "\tTrimming stats will be sent to file: %s\n",statsFile);
		} else if (strcmp(argv[i],"-3")==0) {
			fprintf(stderr, "\tSequences will be trimmed at 3' end");
			i++;
			if (argv[i][0] > 47 && argv[i][0] < 58) {
				sscanf(argv[i],"%d",&p3len);	
				fprintf(stderr, " by %d bp\n", p3len);
			} else {
				fprintf(stderr, " with adapter: %s\n", argv[i]);
				p3adapter = argv[i];
				if (minMatchLength < 0) {
					minMatchLength = strlen(p3adapter)/2;
				}
			}
		} else if (strcmp(argv[i],"-5")==0) {
			fprintf(stderr, "\tSequences will be trimmed at 5' end");
			i++;
			if (argv[i][0] > 47 && argv[i][0] < 58) {
				sscanf(argv[i],"%d",&p5len);	
				fprintf(stderr, " by %d bp\n", p5len);
			} else {
				fprintf(stderr, " with adapter: %s\n", argv[i]);
				p5adapter = argv[i];
				if (minMatchLength < 0) {
					minMatchLength = strlen(p5adapter)/2;
				}
			}
		} else if (argv[i][0] == '-') {
			printCMDtrim();
		} else {
			fprintf(stderr, "\tWill parse file %s\n", argv[i]);
			files[numfiles++] = argv[i];
		}
	}
	if (numfiles < 1) {
		fprintf(stderr, "!!! No input files... Not much to do!!!\n");
		exit(0);
	}

	char* buf = new char[BUFFER];
	char* current = new char[BUFFER];
	char* current2 = new char[BUFFER];
	char* curSeq = new char[BUFFER];
	char* curSeq2 = new char[BUFFER];
	char* curQual = new char[BUFFER];
	char* curQual2 = new char[BUFFER];
	char* curHeader = new char[BUFFER];
	char* curHeader2 = new char[BUFFER];
	char* newfilename = new char[10000];
	char* newfilename2 = new char[10000];
	int* readLengths = new int[maxReadLength+1];
	current[0]='\0';

	for (int i=0;i<numfiles;i++) {

		int zippedFlag = 0;
		int currentFormat = 0;
		char* workingFilename = unzipFileIfNeeded(files[i],zippedFlag,currentFormat);

		FILE* inputfp = fopen(workingFilename,"r");
		if (inputfp == NULL) {
			fprintf(stderr, "!!! Could not open input file: %s - skipping !!!\n", files[i]);
			continue;
		}
		FILE* outfp = stdout;
		FILE* outfp2 = NULL;

		strcpy(newfilename,workingFilename);
		strcat(newfilename,".");
		strcpy(newfilename2,workingFilename);
		strcat(newfilename2,".");
		if (suffix == NULL) {
			if (splitPos > 0) {
				strcat(newfilename,"trimmed1");
				strcat(newfilename2,"trimmed2");
			} else {
				strcat(newfilename,"trimmed");
			}
		} else {
			strcat(newfilename,suffix);
			if (splitPos > 0) {
				strcat(newfilename2,suffix);
				strcat(newfilename,"1");
				strcat(newfilename2,"2");
			}
		}
		outfp= fopen(newfilename, "w");
		if (outfp == NULL) {
			fprintf(stderr, "!!! Could not open output file: %s - skipping !!!\n", newfilename);
			continue;
		}
		if (splitPos > 0) {
			outfp2= fopen(newfilename2, "w");
			if (outfp2 == NULL) {
				fprintf(stderr, "!!! Could not open output file: %s - skipping !!!\n", newfilename2);
				continue;
			}
		}


		FILE* lenfp = stdout;
		strcpy(newfilename,workingFilename);
		strcat(newfilename,".");
		if (lenSuffix == NULL) {
			strcat(newfilename,"lengths");
		} else {
			strcat(newfilename,lenSuffix);
		}
		lenfp= fopen(newfilename, "w");
		if (lenfp == NULL) {
			fprintf(stderr, "!!! Could not open length distribution file: %s - skipping !!!\n", newfilename);
		}

		long long int lastStart = 0;
		long long int currentline = 0;
		int lastLinePlus=0;
		current[0]='\0';
		current2[0]='\0';
		curSeq[0]='\0';
		curSeq2[0]='\0';
		curHeader[0]='\0';
		curHeader2[0]='\0';
		curQual[0]='\0';
		curQual2[0]='\0';
		int start = 0;
		int end = 1000000000;
		int goodReads = 0;
		int totalReads = 0;
		int curLen = 0;
		int maxObservedLength = 0;
		int fastqFlag = 0;
		for (int i=0;i<maxReadLength+1;i++) readLengths[i]=0;

		while (fgets(buf, BUFFER, inputfp) != NULL) {
			currentline++;
			int lineLength = strlen(buf);
			if (lineLength > 1) {
				//trim off newline character
				if (lineLength > 0) buf[lineLength-1]='\0';

				if ((buf[0] == '@' && lastLinePlus==0) || (buf[0] == '>' && fastqFlag != 1)) {
					if ((buf[0] == '@' && lastLinePlus==0)) fastqFlag = 1;
					if (curLen >= minSizeToKeep && curLen < maxReadLength) {
						if (maxObservedLength < curLen) maxObservedLength = curLen;
						readLengths[curLen]++;
						fprintf(outfp, "%s\n%s\n+\n%s\n", curHeader,curSeq,curQual);
						if (splitPos > 0) {
							fprintf(outfp2, "%s\n%s\n+\n%s\n", curHeader2,curSeq2,curQual2);
						}
						goodReads++;
					} else {
						readLengths[0]++;
					}
					totalReads++;
					current[0]='\0';
					current2[0]='\0';
					curSeq[0]='\0';
					curSeq2[0]='\0';
					curHeader[0]='\0';
					curHeader2[0]='\0';
					curQual[0]='\0';
					curQual2[0]='\0';
					strcpy(curHeader,buf);
					if (splitPos > 0) {
						int x = strlen(curHeader);
						current[x-1] = '1';
						current[x] = '\0';
						strcpy(curHeader2,curHeader);
						current2[x-1] = '2';
					}
					lastStart = currentline;
					lastLinePlus=0;
					curLen = 0;
				} else if (buf[0] == '+' && currentline-2==lastStart) {
					//quality header - ignore
					lastLinePlus=1;
				} else {
					if (currentline-1==lastStart) {
						//sequence
						char* seq = buf;
						curLen=lineLength-1;
						start = 0;
						end = curLen;
						if (p5len > 0) {
							if (p5len >= curLen) {
								seq[0]='\0';
								curLen=0;
								start=0;
								end=0;
							} else {
								seq = &(seq[p5len]);
								start=p5len;
								curLen-=p5len;
							}
						}
						if (p3len > 0) {
							if (p3len >= curLen) {
								seq[0]='\0';
								curLen=0;
								start=0;
								end=0;
							} else {
								seq[curLen-p3len] = '\0';
								end=start+curLen-p3len;
								curLen-=p3len;
							}
						}
						if (p5adapter != NULL && curLen > 0) {
							int offset =  matchAdapter5prime(seq, p5adapter,maxMisMatches,minMatchLength, matchStart);
							if (offset >= 0) {
								seq = &(seq[offset+1]);
								curLen -= offset+1;
								start=offset+1;
								end = curLen;
							}
						}
						if (p3adapter != NULL && curLen > 0) {
							int offset =  matchAdapter3prime(seq, p3adapter,maxMisMatches,minMatchLength,matchStart);
							if (offset >= 0) {
								seq[offset]='\0';
								curLen = offset;
								end = offset;
							}
						}
						if (fixedLength > 0) {
							if (curLen > fixedLength) {
								seq[fixedLength] = '\0';
								curLen = fixedLength;
								end = start + fixedLength;
							}
						}
						if (splitPos > 0) {
							if (splitPos < (int)strlen(seq)) {
								char tmp = seq[splitPos];
								seq[splitPos] = '\0';
								strcat(curSeq,seq);
								seq[splitPos] = tmp;
								if (revoppFlag) {
									int L = strlen(seq);
									seq[L-1] = '\0';
									revopp(&(seq[splitPos]));
									seq[L-1] = '\0';
								}
								strcat(curSeq2,&(seq[splitPos]));
							} else {
								strcat(curSeq,seq);
							}
						} else {
							if (revoppFlag) {
								seq[strlen(seq)-1] = '\0';
								revopp(seq);
								seq[strlen(seq)-1] = '\0';
							}
							strcat(curSeq,seq);
						}
					} else {
						//quality scores
						char* qual = buf;
						if (end > lineLength-1) {
							end=0;
						} else {
							qual = &(buf[start]);
						} 
						qual[end]='\0';

						//trim read based on quality score
						if (qthresh > 0) {
							int L = strlen(qual);
							int total = 0;
							int N = 0;
							for (int i=qstart;i<L;i++) {
								if (i==qstart) {
									for (int j=i;j<i+qwindow;j++) {
										if (j>=L) break;
										N++;
										total+=qual[j]-33;
									}
								} else {
									total -= qual[i-1]-33;
									if (i+qwindow < L) {
										total += qual[i+qwindow]-33;
									} else {
										N--;
									}
								}
								if (N > 0 && (total/N < qthresh)) {	
									qual[i] = '\0';
									curSeq[i] = '\0';
									curLen = i-1;
									break;
								} else if (N < 1) {
								}
							}
						}

						if (splitPos > 0) {
							if (splitPos < (int)strlen(qual)) {
								char tmp = qual[splitPos];
								qual[splitPos] = '\0';
								strcat(curQual,qual);
								qual[splitPos] = tmp;
								if (revoppFlag) {
									int L = strlen(qual);
									qual[L-1] = '\0';
									reverse(&(qual[splitPos]));
									qual[L-1] = '\0';
								}
								strcat(curQual2,&(qual[splitPos]));
							} else {
								strcat(curQual,qual);
							}
						} else {
							if (revoppFlag) {
								qual[strlen(qual)-1] = '\0';
								reverse(qual);
								qual[strlen(qual)-1] = '\0';
							}
							strcat(curQual,qual);
						}
					}
					lastLinePlus=0;
				}
			} else {
				lastLinePlus=0;
				if (currentline-1==lastStart || (currentline-3==lastStart && lastLinePlus==1)) {
					//blank sequence or quality information
				} else {
					//whitespace
				}
			}
		}
		//take care of last one
		if (lastLinePlus==0 || fastqFlag != 1) {
			if (curLen >= minSizeToKeep && curLen < maxReadLength) {
				if (maxObservedLength < curLen) maxObservedLength = curLen;
				readLengths[curLen]++;
				fprintf(outfp, "%s\n%s\n+\n%s\n", curHeader,curSeq,curQual);
				if (splitPos > 0) {
					fprintf(outfp2, "%s\n%s\n+\n%s\n", curHeader2,curSeq2,curQual2);
				}
				goodReads++;
			} else {
				readLengths[0]++;
			}
		}

		fprintf(lenfp, "Length\t# reads\tFraction\n");
		for (int i=0;i<maxObservedLength+1;i++) {
			double ratio = ((double)(readLengths[i]))/((double)totalReads);
			fprintf(lenfp, "%d\t%d\t%lf%%\n", i, readLengths[i],ratio*100);
		}
		fclose(lenfp);
		fprintf(stderr, "\tTrimmed output: %d of %d reads\n", goodReads, totalReads);
		fclose(inputfp);
		fclose(outfp);
		if (outfp2 != NULL) fclose(outfp2);

		rezipFileIfNeeded(workingFilename,zippedFlag);

	}		

	delete []current;
	delete []newfilename;
	delete []newfilename2;
	delete []readLengths;
	delete []files;
	delete []buf;
	
	exit(0);
}

int matchAdapter5prime(char* seq, char* adapter, int maxMisMatches,int minMatchLength,int matchStart) {
	int seqLength = strlen(seq)-1; //seq is expected to have a \n character
	int adapterLength = strlen(adapter);
	for (int i=seqLength-1;i>=matchStart;i--) {
		int match=1;
		int numMis = 0;
		int numMatch = 0;
		for (int j=adapterLength-1;j>=0;j--) {
			int nindex = i-(adapterLength-1)+j;
			if (nindex<0) break;
			numMatch++;
			if (seq[nindex] == 'N') continue;
			if (adapter[j] == 'N') continue;
			if (seq[nindex] == adapter[j]) {
				continue;
			} else {
				numMis++;
				if (numMis > maxMisMatches) {
					match=0;
					break;
				}
			}
		}
		if (numMatch < minMatchLength) match = 0;
		if (match==1) {
			return i;
		}
	}
	return -1;
}

int matchAdapter3prime(char* seq, char* adapter, int maxMisMatches, int minMatchLength, int matchStart) {
	int seqLength = strlen(seq)-1; //seq is expected to have a \n character
	int adapterLength = strlen(adapter);
	for (int i=matchStart;i<seqLength;i++) {
		int match=1;
		int numMis = 0;
		int numMatch = 0;
		for (int j=0;j<adapterLength;j++) {
			if (i+j>=seqLength) break;
			numMatch++;
			if (seq[i+j] == 'N') continue;
			if (adapter[j] == 'N') continue;
			if (seq[i+j] == adapter[j]) {
				continue;
			} else {
				numMis++;
				if (numMis > maxMisMatches) {
					match=0;
					break;
				}
			}
		}
		if (numMatch < minMatchLength) match = 0;
		if (match==1) {
			return i;
		}
	}
	return -1;
}
//---------------------------- truseq -------------------------------------------------------------------

#define MAX_FREQ_TO_REPORT 30
void getTruSeqBarcode(char* buf,int lineLength,char* barcode);
int processBarcode(char* barcode, char* read, Hashtable* barcodeFPs, LongInttable* barcodeTotals,
							char* prefix,int printFlag);
void printCMDtruSeq() {
	fprintf(stderr, "\n\tUsage: homerTools truseq <fastq file1> [file2] ...\n");
	fprintf(stderr, "\n\tOptions for command: truseq\n");
	fprintf(stderr, "\t\t-min <#> (Minimum frequency of barcodes to keep: default=%.3f\n", BARCODE_MINFREQ);
	fprintf(stderr, "\t\t-freq <filename> (output file for barcode frequencies, default=file.freq.txt)\n");
	fprintf(stderr, "\t\t-qual <#> (Minimum quality score for barcode nucleotides, default=not used)\n");
	fprintf(stderr, "\t\t-qualBase <character> (Minimum quality character in FASTQ file, default=B)\n");
	fprintf(stderr, "\t\t\n"); 
	exit(0);
}
void truSeqProgram(int argc, char** argv) {

	char* freqFileName = NULL;
	char** files = new char*[100000];	
	int numfiles = 0;
	int qualThreshold = -120;
	char qualBase = 'B';
	int printFlag = 0;

	if (argc < 3) printCMDtruSeq();

	//int barcodeSize = 0;
	double minPercentage = BARCODE_MINFREQ;
	fprintf(stderr, "\tDefault minimum barcode percentage is %.0lf%%\n", minPercentage*100);

	for (int i=2;i<argc;i++) { 
		if (strcmp(argv[i],"-min")==0) {
			sscanf(argv[++i],"%lf",&minPercentage);
			fprintf(stderr, "\tMinimum percentage for barcodes set at %.2lf%%\n", minPercentage*100);
		} else if (strcmp(argv[i],"-freq")==0) {
			freqFileName = argv[++i];
			fprintf(stderr, "\tBarcode frequencies will be recorded in *.%s files\n", freqFileName);
		} else if (strcmp(argv[i],"-qualBase")==0) {
			sscanf(argv[++i],"%c",&qualBase);
		} else {
			files[numfiles++] = argv[i];
		}
	}

	if (numfiles < 1) {
		fprintf(stderr, "!!! No input files... Not much to do!!!\n");
		exit(0);
	}


	char* prefix = NULL;
	char* buf = new char[BUFFER];
	char* current = new char[BUFFER];
	char* barcode = new char[BUFFER];
	char* name= new char[10000];
	current[0]='\0';
	barcode[0]='\0';

	for (int i=0;i<numfiles;i++) {
		FILE* fp = fopen(files[i], "r");
		if (fp == NULL) {
			fprintf(stderr, "!!! Could not open %s\n", files[i]);
			continue;
		}
		prefix = files[i];
		Hashtable* barcodeFPs = new Hashtable();	
		LongInttable* barcodeTotals = new LongInttable();	
		long long int totalReads = 0;
		long long int lastStart = 0;
		long long int currentline = 0;
		int lastLinePlus=0;
		int format = SEQFILE_FORMAT_UNKNOWN ;
		char minQualityScore=0;
		current[0]='\0';
		barcode[0]='\0';
		while (fgets(buf, BUFFER, fp) != NULL) {
			currentline++;
			int lineLength = strlen(buf);
			if (lineLength > 1) {
				if (format == SEQFILE_FORMAT_UNKNOWN) {
					if (buf[0] == '@') {
						format = SEQFILE_FORMAT_FASTQ;
						fprintf(stderr, "\tFASTQ format detected\n");
					} else if (buf[0] == '>') {
						format = SEQFILE_FORMAT_FASTA;
						fprintf(stderr, "\tFASTA format detected\n");
					}
				}
				if ((format == SEQFILE_FORMAT_FASTQ && buf[0] == '@' && lastLinePlus==0)  
						|| (format == SEQFILE_FORMAT_FASTA && buf[0] == '>')) {
					if (minQualityScore >= qualThreshold) {
						totalReads += processBarcode(barcode,current,barcodeFPs,barcodeTotals,prefix,printFlag);
					}
					current[0]='\0';
					barcode[0]='\0';
					strcpy(current,buf);
					lastStart = currentline;
					lastLinePlus=0;
					minQualityScore=127;
					getTruSeqBarcode(buf,lineLength,barcode);
					//fprintf(stdout, "%s\n", barcode);

				} else if (buf[0] == '+' && currentline-2==lastStart) {
					//quality header
					lastLinePlus=1;
					strcat(current,buf);
				} else {
					lastLinePlus=0;
					strcpy(current,buf);
				}
			} else {
				lastLinePlus=0;
				if (currentline-1==lastStart || (currentline-3==lastStart && lastLinePlus==1)) {
					//blank sequence or quality information
					strcat(current,"\n");
				} else {
					//whitespace
				}
			}
		}
		if ((format == SEQFILE_FORMAT_FASTQ && lastLinePlus==0)  
								|| (format == SEQFILE_FORMAT_FASTA)) {
			if (minQualityScore >= qualThreshold) {
				totalReads += processBarcode(barcode,current,barcodeFPs,barcodeTotals,prefix,printFlag);
			}
		}
		fclose(fp);

		char** keys = barcodeTotals->keys();
		StrSort* data = new StrSort[barcodeTotals->total];
		for (int j=0;j<barcodeTotals->total;j++) {
			long long int barcodeTotal = barcodeTotals->search(keys[j]);
			if (printFlag) {
				FILE* bfp = (FILE*) barcodeFPs->search(keys[j]);
				fclose(bfp);
			}
			double ratio = ((double)barcodeTotal)/((double)totalReads);
			char* readInfo = new char[10000];
			sprintf(readInfo, "%s\t%lld\t%lf\n", keys[j], barcodeTotal,ratio);
			data[j].str = readInfo;
			data[j].v = ratio;
			if (ratio < minPercentage) {
				strcpy(name, prefix);
				strcat(name,".");
				strcat(name,keys[j]);
				remove(name);
			}
			delete [](keys[j]);
		}
		delete []keys;
	
		qsort(data,barcodeTotals->total,sizeof(StrSort),&decendStrSort);

		FILE* freqFP = NULL;
		strcpy(name,files[i]);
		strcat(name,".freq.txt");
		freqFP = fopen(name, "w");
		if (freqFP == NULL) {
			fprintf(stderr, "!!! Problem opening %s for writing (frequency file)!!!\n", name);
			exit(0);
		}
		fprintf(freqFP, "Barcode\tTotal Reads(of %lld)\tFrequency\n", totalReads);
		fprintf(stderr, "\tBarcode\tTotal Reads(of %lld)\tFrequency\n", totalReads);
		for (int i=0;i<barcodeTotals->total;i++) {
			fprintf(freqFP,"%s",data[i].str);
			if (i<MAX_FREQ_TO_REPORT) fprintf(stderr,"\t\t%s",data[i].str);
			delete [](data[i].str);
		}
		delete []data;
	
		if (freqFileName != NULL) {
			fclose(freqFP);
		}
		delete barcodeFPs;
		delete barcodeTotals;
	}

	delete []name;
	delete []files;
	delete []buf;
}
void getTruSeqBarcode(char* buf,int lineLength,char* barcode) {
	int lastColon = 0;
	int lastBp = 0;
	for (int i=0;i<lineLength;i++) {
		if (buf[i] == ':') {
			lastColon=i;
		} else if (buf[i] == 'A' || buf[i] == 'C' || buf[i]=='G' || buf[i]=='T' || buf[i] == 'N') {
			lastBp = i;
		}
	}
	if (lastBp-lastColon < 1) {
		barcode[0] = '\0';
	} else {
		strncpy(barcode,&(buf[lastColon+1]),lastBp-lastColon);
	}
}

//------------------------------  barcodes -------------------------------------------------------------------


void printCMDbarcodes() {
	fprintf(stderr, "\n\tUsage: homerTools barcodes <# bp in barcode> [options] <fastq file1> [file2] ...\n");
	fprintf(stderr, "\n\t3rd argument must be the number bp in the barcode\n");
	fprintf(stderr, "\n\tOptions for command: barcode\n");
	fprintf(stderr, "\t\t-min <#> (Minimum frequency of barcodes to keep: default=%.3f\n", BARCODE_MINFREQ);
	fprintf(stderr, "\t\t-freq <filename> (output file for barcode frequencies, default=file.freq.txt)\n");
	fprintf(stderr, "\t\t-qual <#> (Minimum quality score for barcode nucleotides, default=not used)\n");
	fprintf(stderr, "\t\t-qualBase <character> (Minimum quality character in FASTQ file, default=B)\n");
	fprintf(stderr, "\t\t\n"); 
	exit(0);
}
void barcodesProgram(int argc, char** argv) {

	char* freqFileName = NULL;
	char** files = new char*[100000];	
	int numfiles = 0;
	int qualThreshold = -120;
	char qualBase = 'B';
	int printFlag = 1;

	if (argc < 4) printCMDbarcodes();

	int barcodeSize = 0;
	double minPercentage = BARCODE_MINFREQ;
	sscanf(argv[2],"%d",&barcodeSize);
	if (barcodeSize < 1 || barcodeSize > 100) {
		fprintf(stderr, "!!! Barcode of size %d seems odd - make sure you use the program correctly!!!\n", barcodeSize);
	} else {
		fprintf(stderr, "\tUsing barcode size of %d bp\n", barcodeSize);
	}
	fprintf(stderr, "\tDefault minimum barcode percentage is %.0lf%%\n", minPercentage*100);

	for (int i=3;i<argc;i++) { 
		if (strcmp(argv[i],"-min")==0) {
			sscanf(argv[++i],"%lf",&minPercentage);
			fprintf(stderr, "\tMinimum percentage for barcodes set at %.2lf%%\n", minPercentage*100);
		} else if (strcmp(argv[i],"-freq")==0) {
			freqFileName = argv[++i];
			fprintf(stderr, "\tBarcode frequencies will be recorded in *.%s files\n", freqFileName);
		} else if (strcmp(argv[i],"-qualBase")==0) {
			sscanf(argv[++i],"%c",&qualBase);
			fprintf(stderr, "\tQuality base character set to %c (%d)\n", qualBase, qualBase);
		} else if (strcmp(argv[i],"-qual")==0) {
			sscanf(argv[++i],"%d",&qualThreshold);
			fprintf(stderr, "\tBarcode quality threshold set at %d\n", qualThreshold);
		} else if (argv[i][0] == '-') {
			printCMDbarcodes();
		} else {
			fprintf(stderr, "\tWill parse file %s\n", argv[i]);
			files[numfiles++] = argv[i];
		}
	}
	if (numfiles < 1) {
		fprintf(stderr, "!!! No input files... Not much to do!!!\n");
		exit(0);
	}


	char* prefix = NULL;
	char* buf = new char[BUFFER];
	char* current = new char[BUFFER];
	char* barcode = new char[BUFFER];
	char* name= new char[10000];
	current[0]='\0';
	barcode[0]='\0';

	for (int i=0;i<numfiles;i++) {
		FILE* fp = fopen(files[i], "r");
		if (fp == NULL) {
			fprintf(stderr, "!!! Could not open %s\n", files[i]);
			continue;
		}
		prefix = files[i];
		Hashtable* barcodeFPs = new Hashtable();	
		LongInttable* barcodeTotals = new LongInttable();	
		long long int totalReads = 0;
		long long int lastStart = 0;
		long long int currentline = 0;
		int lastLinePlus=0;
		int format = SEQFILE_FORMAT_UNKNOWN ;
		char minQualityScore=0;
		current[0]='\0';
		barcode[0]='\0';
		while (fgets(buf, BUFFER, fp) != NULL) {
			currentline++;
			int lineLength = strlen(buf);
			if (lineLength > 1) {
				if (format == SEQFILE_FORMAT_UNKNOWN) {
					if (buf[0] == '@') {
						format = SEQFILE_FORMAT_FASTQ;
						fprintf(stderr, "\tFASTQ format detected\n");
					} else if (buf[0] == '>') {
						format = SEQFILE_FORMAT_FASTA;
						fprintf(stderr, "\tFASTA format detected\n");
					}
				}
				if ((format == SEQFILE_FORMAT_FASTQ && buf[0] == '@' && lastLinePlus==0)  
						|| (format == SEQFILE_FORMAT_FASTA && buf[0] == '>')) {
					if (minQualityScore >= qualThreshold) {
						totalReads += processBarcode(barcode,current,barcodeFPs,barcodeTotals,prefix,printFlag);
					}
					current[0]='\0';
					barcode[0]='\0';
					strcpy(current,buf);
					lastStart = currentline;
					lastLinePlus=0;
					minQualityScore=127;
				} else if (buf[0] == '+' && currentline-2==lastStart) {
					//quality header
					lastLinePlus=1;
					strcat(current,buf);
				} else {
					if (lineLength-1 >= barcodeSize) {
						strcat(current,&(buf[barcodeSize]));
					} else {
						strcat(current,"\n");
					}
					if (currentline-1==lastStart) {
						//sequence - get barcode
						if (lineLength-1 >= barcodeSize) {
							strncpy(barcode,buf,barcodeSize);
						}
					} else {
						if (lineLength-1 < barcodeSize) {
							minQualityScore = -1;
						} else {
							minQualityScore = buf[0]-qualBase;
							for (int j=1;j<barcodeSize;j++) {
								if (buf[j]-qualBase < minQualityScore) minQualityScore = buf[j]-qualBase;
							}
						}
					}
					lastLinePlus=0;
				}
			} else {
				lastLinePlus=0;
				if (currentline-1==lastStart || (currentline-3==lastStart && lastLinePlus==1)) {
					//blank sequence or quality information
					strcat(current,"\n");
				} else {
					//whitespace
				}
			}
		}
		if ((format == SEQFILE_FORMAT_FASTQ && lastLinePlus==0)  
								|| (format == SEQFILE_FORMAT_FASTA)) {
			if (minQualityScore >= qualThreshold) {
				totalReads += processBarcode(barcode,current,barcodeFPs,barcodeTotals,prefix,printFlag);
			}
		}
		fclose(fp);

		char** keys = barcodeFPs->keys();
		StrSort* data = new StrSort[barcodeFPs->total];
		for (int j=0;j<barcodeFPs->total;j++) {
			long long int barcodeTotal = barcodeTotals->search(keys[j]);
			FILE* bfp = (FILE*) barcodeFPs->search(keys[j]);
			fclose(bfp);
			double ratio = ((double)barcodeTotal)/((double)totalReads);
			char* readInfo = new char[10000];
			sprintf(readInfo, "%s\t%lld\t%lf\n", keys[j], barcodeTotal,ratio);
			data[j].str = readInfo;
			data[j].v = ratio;
			if (ratio < minPercentage) {
				strcpy(name, prefix);
				strcat(name,".");
				strcat(name,keys[j]);
				remove(name);
			}
			delete [](keys[j]);
		}
		delete []keys;
	
		qsort(data,barcodeFPs->total,sizeof(StrSort),&decendStrSort);

		FILE* freqFP = NULL;
		strcpy(name,files[i]);
		strcat(name,".freq.txt");
		freqFP = fopen(name, "w");
		if (freqFP == NULL) {
			fprintf(stderr, "!!! Problem opening %s for writing (frequency file)!!!\n", name);
			exit(0);
		}
		fprintf(freqFP, "Barcode\tTotal Reads(of %lld)\tFrequency\n", totalReads);
		fprintf(stderr, "\tBarcode\tTotal Reads(of %lld)\tFrequency\n", totalReads);
		for (int i=0;i<barcodeFPs->total;i++) {
			fprintf(freqFP,"%s",data[i].str);
			fprintf(stderr,"\t\t%s",data[i].str);
			delete [](data[i].str);
		}
		delete []data;
	
		if (freqFileName != NULL) {
			fclose(freqFP);
		}
		delete barcodeFPs;
		delete barcodeTotals;
	}

	delete []name;
	delete []files;
	delete []buf;
}


int processBarcode(char* barcode, char* read, Hashtable* barcodeFPs, LongInttable* barcodeTotals,
												char* prefix,int printFlag) {
	if (barcode[0] == '\0') return 0;
	if (read[0] == '\0') return 0;
	long long int total = barcodeTotals->search(barcode);
	FILE* fp = NULL;
	if (total == EMPTY_INT) {
		barcodeTotals->insert(1,barcode);
	} else {
		barcodeTotals->insert(total+1,barcode);
	}	
	if (printFlag) {
		if (total == EMPTY_INT) {
			char* name= new char[10000];
			strcpy(name, prefix);
			strcat(name, ".");
			strcat(name, barcode);
			fp = fopen(name,"w");
			barcodeFPs->insert(fp,barcode);
			delete []name;
		} else {
			barcodeTotals->insert(total+1,barcode);
			fp = (FILE*) barcodeFPs->search(barcode);
		}
		fprintf(fp,"%s",read);
	}
	return 1;
}



int decendStrSort(const void* a, const void* b) {
	StrSort* aa = (StrSort*)a;
	StrSort* bb = (StrSort*)b;
	if (aa->v > bb->v) return -1;
	if (aa->v < bb->v) return 1;
	return 0;
}


//---------------------------- freq ---------------------------


void printCMDfreq() {
	fprintf(stderr, "\n\tUsage: homerTools freq [options] <sequence file>\n");
	fprintf(stderr, "\n\tOutputs nucleotide frequencies to stdout\n");
	fprintf(stderr, "\n\tOptions for command: freq\n");
	fprintf(stderr, "\t\t-format <tsv|fasta|fastq> (sequence file format, default: auto detect)\n");
	fprintf(stderr, "\t\t-offset <#> (offset of first base in output file, default: 0)\n");
	fprintf(stderr, "\t\t-maxlen <#> (Maximum length of sequences to consider, default: length of 1st seq)\n");
	fprintf(stderr, "\t\t-o <filename> (Output filename, default: output sent to stdout)\n");
	fprintf(stderr, "\t\t-gc <filename> (calculate CpG/GC content per sequence output to \"filename\")\n");
	fprintf(stderr, "\t\t\tOutputFormat: name<tab>CpG<tab>GC<tab>AG<tab>AC<tab>Length\n");
	fprintf(stderr, "\t\t\n"); 
	exit(0);
}
void freqProgram(int argc, char** argv) {

	char* outputFilename = NULL;
	char* gcFilename = NULL;
	char* seqFile = NULL;
	int fileFormat = SEQFILE_FORMAT_UNKNOWN;
	int maxLength = -1;
	int offset = 0;

	if (argc < 3) {
		printCMDfreq();
	}
	fprintf(stderr, "\n");
	for (int i=2;i<argc;i++) {
		if (strcmp(argv[i],"-offset")==0) {
			sscanf(argv[++i],"%d",&offset);
			//fprintf(stderr, "\tOffset set to %d\n",offset);
		} else if (strcmp(argv[i],"-maxlen")==0) {
			sscanf(argv[++i],"%d",&maxLength);
			//fprintf(stderr, "\tOffset set to %d\n",offset);
		} else if (strcmp(argv[i],"-format")==0) {
			i++;
			if (strcmp(argv[i],"tsv")==0) {
				fileFormat = SEQFILE_FORMAT_TSV;
			} else if (strcmp(argv[i],"fasta")==0) {
				fileFormat = SEQFILE_FORMAT_FASTA;
			} else if (strcmp(argv[i],"fastq")==0) {
				fileFormat = SEQFILE_FORMAT_FASTQ;
			} else if (strcmp(argv[i],"fa")==0) {
				fileFormat = SEQFILE_FORMAT_FASTA;
			} else if (strcmp(argv[i],"fq")==0) {
				fileFormat = SEQFILE_FORMAT_FASTQ;
			} else {
				fprintf(stderr, "!!! Could not recognize format: %s !!!\n", argv[i]);
				printCMDfreq();
			}
		} else if (strcmp(argv[i],"-o")==0) {
			outputFilename = argv[++i];
		} else if (strcmp(argv[i],"-gc")==0) {
			gcFilename = argv[++i];
		} else {
			if (seqFile == NULL) {
				seqFile = argv[i];
			} else {
				fprintf(stderr, "!!! Error- multiple sequence files specified: %s & %s\n", seqFile, argv[i]);
				printCMDfreq();
			}
		}
	}

	int autoFormat = determineSeqFileFormat(seqFile);
	if (fileFormat != SEQFILE_FORMAT_UNKNOWN) {
		if (autoFormat != fileFormat) {
			fprintf(stderr, "!!! Auto-check for your file format disagrees with your selected parameters\n");
			fprintf(stderr, "!!! Check your file and your options! Continuing anyway...\n");
		}
	} else {
		fileFormat = autoFormat;
	}

	NucleotideFreq* nf = new NucleotideFreq(seqFile, fileFormat, offset, maxLength,gcFilename);

	FILE* fp = stdout;
	if (outputFilename != NULL) {
		fp = fopen(outputFilename, "w");
		if (fp == NULL) {
			fprintf(stderr, "!!! Could not open %s for writing!!!\n", outputFilename);
			exit(0);
		}
	}

	nf->print(fp);

	if (outputFilename != NULL) {
		fclose(fp);
	}
	fprintf(stderr, "\n");
}

//---------------------------- extract ---------------------------


void printCMDextract() {
	fprintf(stderr, "\n\tUsage: homerTools extract <peak file/BED file> <Directory of FASTA files> [options]\n");
	fprintf(stderr, "\n\tThe <Directory of FASTA files> can be a single FASTA file instead.\n");
	fprintf(stderr, "\tIf using a Directory, files should be named with chromosomes,\n");
	fprintf(stderr, "\t\ti.e. chr1.fa or chr1.fa.masked or genome.fa/genome.fa.masked\n");
	fprintf(stderr, "\t\tIf having trouble, place all FASTA entries in single file instead of a directory\n");
	fprintf(stderr, "\t\tFASTA format: >chrname ... (anything after whitespace will be ignored)\n");
	fprintf(stderr, "\tThis program will output sequences to stdout in tab-delimited format\n");
	fprintf(stderr, "\n\tOptions for command: extract\n");
	fprintf(stderr, "\t\t-fa (output sequences in FASTA format - default is tab-delimited format)\n"); 
	fprintf(stderr, "\t\t-mask (mask out lower case sequence from genome)\n"); 
	fprintf(stderr, "\n\tAlternate Usage: homerTools extract stats <Directory of FASTA files>\n");
	fprintf(stderr, "\t\tDisplays stats about the genome files (such as length)\n");
	fprintf(stderr, "\t\t\n"); 
	exit(0);
}
void extractProgram(int argc, char** argv) {

	if (argc < 4) {
		printCMDextract();
	}

	char* posfile = argv[2];
	char* genomeDir = argv[3];

	int fastaFlag = 0;
	int statsFlag = 0;
	int maskFlag = 0;
	if (strcmp(posfile,"stats")==0) {
		statsFlag = 1;
	}

	fprintf(stderr, "\n");
	for (int i=4;i<argc;i++) {
		if (strcmp(argv[i],"-fa")==0) {
			fastaFlag = 1;
		} else if (strcmp(argv[i],"-mask")==0) {
			maskFlag = 1;
		} else {
			printCMDextract();
		}
	}

	FILE* fpout = stdout;

	if (statsFlag) {
		PeakLibrary::extractSequenceStats(genomeDir,fpout);
	} else {
		PeakLibrary* peaks = new PeakLibrary(posfile, PEAK_READ_MODE_NORMAL);
		peaks->extractSequence(genomeDir,fpout,fastaFlag,maskFlag);
	}

	exit(0);
}



//---------------------------- decontaminate ---------------------------


void printCMDdecontaminate() {
	fprintf(stderr, "\n\tUsage: homerTools decontaminate <Contaminated Tag Directory>\n");
	fprintf(stderr, "\t\t\t\t<Contaminant(background) Tag Directory> [options]\n");
	fprintf(stderr, "\n\tRemoves reads from an experiment that contains reads originating from another.\n");
	fprintf(stderr, "\tEither specify the fraction of contaminated reads (-frac <#>), or the program\n");
	fprintf(stderr, "\twill attempt to estimate (only works for quanitatively different experiments)\n");
	fprintf(stderr, "\tCreates contaminationHistogram.txt and contaminationScatterPlot.txt in the\n");
	fprintf(stderr, "\toutput directory to help guage the extent of contamination when estimating.\n");
	fprintf(stderr, "\n\tTwo tag directories are required after \"homerTools decontaminate\":\n");
	fprintf(stderr, "\t\t<Contaminated Tag Directory> <Contaminant Tag Directory>\n");
	fprintf(stderr, "\t\tThe Contaminated Tag Directory will be modified (recommend to copy the original)\n");
	fprintf(stderr, "\n\tOptions for command: decontaminate\n");
	fprintf(stderr, "\t\t-frac <#> (Estimate fraction of sample that is contaminated, default: auto)\n"); 
	fprintf(stderr, "\t\t-estimateOnly (Only estimate the contamination, do not decontaminate)\n"); 
	fprintf(stderr, "\t\t-o <output tag directory> (default: overrites contaminated tag directory)\n"); 
	fprintf(stderr, "\t\t-size <#> (Peak size for estimating contamination/Max distance from contaminant\n");
	fprintf(stderr, "\t\t\treads to remove contaminated reads, default: 250)\n"); 
	fprintf(stderr, "\t\t-min <#> (Minimum tag count to consider when estimating contamination, default: 20)\n"); 
	fprintf(stderr, "\t\t\n"); 
	exit(0);
}
void decontaminateProgram(int argc, char** argv) {

	if (argc < 4) {
		printCMDdecontaminate();
	}

	char* expDir = argv[2];
	char* inputDir = argv[3];
	char* outputDir = NULL;
	double fraction = -1.0;
	float minThreshold = 20.0;
	int estonly = 0;

	int size = 250;

	fprintf(stderr, "\n");
	for (int i=4;i<argc;i++) {
		if (strcmp(argv[i],"-frac")==0) {
			sscanf(argv[++i],"%lf",&fraction);
		} else if (strcmp(argv[i],"-o")==0) {
			outputDir= argv[++i];
		} else if (strcmp(argv[i],"-estimateOnly")==0) {
			estonly = 1;
		} else if (strcmp(argv[i],"-min")==0) {
			sscanf(argv[++i],"%f",&minThreshold);
		} else if (strcmp(argv[i],"-size")==0) {
			sscanf(argv[++i],"%d",&size);
		} else {
			printCMDdecontaminate();
		}
	}
	fprintf(stderr, "\n\tContaminated Directory: %s\n", expDir);
	fprintf(stderr, "\tContaminant Directory:  %s\n", inputDir);
	fprintf(stderr, "\n");
	if (outputDir == NULL) {
		if (estonly == 0) {
			fprintf(stderr, "\tWarning - this will overrite the data in %s\n", expDir);
			fprintf(stderr, "\tIt is recommended that you make a copy first\n");
			fprintf(stderr, "\tIf you wish to stop, hit CTRL+C NOW!!!\n");
			for (int i=5;i>0;i--) {
				fprintf(stderr, "\t\t%d\n",i);
				(void)system("sleep 1");
			}
		}
		outputDir = expDir;
	} else {
		fprintf(stderr, "\tOutput Directory:    %s\n", outputDir);
		char* command = new char[10000];
		sprintf(command, "mkdir -p \"%s\"",outputDir);
		(void)system(command);
		sprintf(command, "cp \"%s\"/* \"%s\"/",expDir, outputDir);
		(void)system(command);
	}

	TagLibrary* tags = new TagLibrary(outputDir);
	TagLibrary* input = new TagLibrary(inputDir);

	tags->readTagDirectory();
	input->readTagDirectory();
	if (fraction < 0.0) {
		fprintf(stderr, "\tEstimating fraction of reads in from contaminated experiement\n\n");
		fraction = tags->estimateContamination(input,size,minThreshold);
		if (estonly) exit(0);
	} else {
		fprintf(stderr, "\tFraction of reads in from contaminated experiement: %lf%% [user defined]\n",
											fraction*100.0);
	}
	tags->decontaminate(input,fraction,size);
	//should probably do more quality control...

	fprintf(stderr, "\n");
	exit(0);
}


//---------------------------- special ---------------------------

void specialUniqMapProgram(int argc, char** argv);

void printCMDspecial() {
	fprintf(stderr, "\n\tUsage: homerTools special <2nd command> ...\n");
	fprintf(stderr, "\n\tHighly specialized programs\n");
	fprintf(stderr, "\n\t2nd Commands:\n");
	fprintf(stderr, "\t\tuniqmap (parse Unique mappability files)\n"); 
	//fprintf(stderr, "\t\ttile (create peak file tiling the genome)\n"); 
	fprintf(stderr, "\t\t\n"); 
	exit(0);
}
void specialProgram(int argc, char** argv) {

	if (argc < 3) {
		printCMDspecial();
	}
	for (int i=2;i<argc;i++) {
		if (strcmp(argv[i],"uniqmap")==0) {
			specialUniqMapProgram(argc,argv);
		} else {
			printCMDspecial();
		}
	}
	exit(0);
}

//---------------------------- special uniqmap ---------------------------


void printCMDspecialUniqMap() {
	fprintf(stderr, "\n\tUsage: homerTools special uniqmap <uniqmap directory(output)> <uniq map file(input)>\n");
	fprintf(stderr, "\n\tOptions for command: special uniqmap\n");
	fprintf(stderr, "\t\tnone\n"); 
	fprintf(stderr, "\t\t\n"); 
	exit(0);
}
void specialUniqMapProgram(int argc, char** argv) {

	if (argc < 5) {
		printCMDspecialUniqMap();
	}

	char* directory = argv[3];
	char* inputFile = argv[4];

	char* command = new char[100000];
	strcpy(command, "mkdir -p \"");
	strcat(command, directory);
	strcat(command, "\"");
	(void)system(command);

	char* buf = new char[BUFFER];
	char** line = new char*[10000];
	int numCols = 0;

	int initializationSize = 400000000;

	FILE* fp = fopen(inputFile,"r");
	if (fp == NULL) {
		fprintf(stderr, "!! Could not open input file: %s\n", inputFile);
		exit(0);
	}

	Hashtable* uniqmapchr = new Hashtable();
	char* currentChr = new char[10000];
	currentChr[0] = '\0';


	UniqMapChrs* umc = NULL;
	int format = -1;

	while (fgets(buf, BUFFER, fp) != NULL) {
		split(buf, line, numCols, '\t');
		
		if (format == -1) {
			if (numCols < 3) continue;
			char* chr = line[0];
			int count = 0;
			while (chr[count] != '\0') {
				if (chr[count] == '.') {
					chr[count] = '\0';
					break;
				}
				count++;
			}
			if (strcmp(chr,currentChr) != 0) {
				strcpy(currentChr,chr);
				umc = (UniqMapChrs*) uniqmapchr->search(chr);
				if (umc == NULL) {
					umc = new UniqMapChrs(chr,directory,1); //1 to overwrite existing data
					umc->increaseMapSize(initializationSize);
					uniqmapchr->insert(umc,chr);
				}
			}
			int position = -1;
			int strand = -1;
			sscanf(line[1],"%d",&position);
			sscanf(line[2],"%d",&strand);
			umc->setMappable(position,strand);
		}
	}
	fclose(fp);

	sprintf(command, "%s/uniqMapStats.txt", directory);
	fp = fopen(command, "w");
	long long int gsize=0;
	long long int mappability = 0;
	fprintf(fp, "#Chr\tMappable bp (both strands)\tTotal bp (both strand)\n");

	char** keys = uniqmapchr->keys();
	for (int i=0;i<uniqmapchr->total;i++) {
		umc = (UniqMapChrs*) uniqmapchr->search(keys[i]);
		umc->saveFiles();
		gsize += umc->size;
		mappability += umc->numMappable;
		fprintf(fp, "%s\t%d\t%d\n",keys[i],umc->numMappable,umc->size);
		delete [](keys[i]);
	}
	delete []keys;
	fprintf(fp, "genome\t%lld\t%lld\n",mappability,gsize);
	fclose(fp);
	

	delete uniqmapchr;
	delete []buf;
	delete []currentChr;
	delete []line;
	delete []command;
	exit(0);
}


//---------------------------- special cluster ---------------------------


void printCMDcluster() {
	fprintf(stderr, "\n\tUsage: homerTools cluster [options]\n");
	fprintf(stderr, "\n\tOptions for command: cluster\n");
	fprintf(stderr, "\t\t-d <distance matrix file> (tab delimited)\n"); 
	fprintf(stderr, "\t\t\t-annCols <#> (number of annotation columns at beginning, def: 1)\n"); 
	fprintf(stderr, "\t\t-prefix <distance matrix file> (tab delimited)\n"); 
	fprintf(stderr, "\t\t-r (reverse distance metric = high values are \"closer\")\n"); 
	fprintf(stderr, "\n\tExtract Clusters from output gtr or atr format file:\n"); 
	fprintf(stderr, "\t\t-gtr <gtr file> (output from clustering - cdt must be present too)\n"); 
	fprintf(stderr, "\t\t-thresh <#> (threshold to define clusters)\n"); 
	fprintf(stderr, "\n\tScore regions against clusters (use either -gtr/-d to supply distance matrix):\n");
	fprintf(stderr, "\t\t-c <cluster file> [cluster file 2] ... (First column is region name)\n");
	fprintf(stderr, "\t\t\n"); 
	exit(0);
}


void clusterProgram(int argc, char** argv) {

	fprintf(stderr, "\n");
	if (argc < 3) {
		printCMDcluster();
	}

	char* distMatrixFile = NULL;
	char* prefix = NULL;
	char** clusterFiles = new char*[10000];
	int numClusterFiles = 0;
	int reverseFlag = 0;
	int numAnnCols = 1;

	char* gtrInputFile = NULL;
	double clusterThreshold = 0.0;

	for (int i=2;i<argc;i++) {
		if (strcmp(argv[i],"-d")==0) {
			distMatrixFile = argv[++i];
		} else if (strcmp(argv[i],"-prefix")==0) {
			prefix = argv[++i];
		} else if (strcmp(argv[i],"-gtr")==0) {
			gtrInputFile = argv[++i];
		} else if (strcmp(argv[i],"-thresh")==0) {
			sscanf(argv[++i],"%lf",&clusterThreshold);
		} else if (strcmp(argv[i],"-annCols")==0) {
			sscanf(argv[++i],"%d",&numAnnCols);
		} else if (strcmp(argv[i],"-c")==0) {
			i++;
			while (i<argc && argv[i][0] != '-') {
				clusterFiles[numClusterFiles++]=argv[i];
				i++;
			}
			if (i<argc && argv[i][0] == '-') i--;
		} else if (strcmp(argv[i],"-r")==0) {
			reverseFlag = 1;
		} else {
			fprintf(stderr, "!!! Invalid option: %s\n", argv[i]);
			printCMDcluster();
		}
	}

	FILE* clusterScoreFP = stdout;

	TreeCluster* cluster = NULL;

	if (gtrInputFile != NULL) {
		cluster = new TreeCluster();
		cluster->loadGTR(gtrInputFile);

		if (numClusterFiles < 1) {
			cluster->printClusters(stdout,NULL,clusterThreshold);
			exit(0);
		} else {
			cluster->scoreClusterFiles(clusterFiles, numClusterFiles,clusterScoreFP);
			exit(0);
		}
	}

	char* buf = new char[BUFFER];
	char** line = new char*[10000];
	int numCols = 0;

	FILE* fp = fopen(distMatrixFile,"r");
	if (fp == NULL) {
		fprintf(stderr, "!! Could not open input file: %s\n", distMatrixFile);
		exit(0);
	}


	int numMatrix = 0;
	double** distMatrix = NULL;
	char** names = NULL;
	
	int lineCount = 0;
	while (fgets(buf, BUFFER, fp) != NULL) {
		split(buf, line, numCols, '\t');
		lineCount++;
		if (lineCount == 1) {
			numMatrix = numCols-numAnnCols;
			distMatrix = new double*[numMatrix];
			names = new char*[numMatrix];
			for (int i=0;i<numMatrix;i++) {
				distMatrix[i] = new double[numMatrix];
				for (int j=0;j<numMatrix;j++) distMatrix[i][j]=FLT_MAX;
			}
			continue;
		}
		if (numCols < numMatrix+numAnnCols) {
			fprintf(stderr, "Something is wrong with your distance matrix!!!\n");
			exit(0);
		}
		int index = lineCount -2;
		names[index] = new char[strlen(line[0])+1];
		strcpy(names[index], line[0]);
		for (int i=numAnnCols;i<numCols;i++) {
			sscanf(line[i], "%lf", &(distMatrix[index][i-numAnnCols]));
		}
	}
	fclose(fp);

	if (reverseFlag) {
		for (int i=0;i<numMatrix;i++) {
			for (int j=0;j<numMatrix;j++) {	
				distMatrix[i][j] *= -1;
			}
		}
	}	


	cluster = new TreeCluster(distMatrix,numMatrix);
	if (numClusterFiles > 0) {
		cluster->names = names;
		cluster->scoreClusterFiles(clusterFiles, numClusterFiles,clusterScoreFP);
	} else {
		cluster->cluster();
		cluster->printCDT(prefix,names,distMatrix);
	}
	exit(0);
}

//----------------------- matrix3D program

void printCMDmatrix3D() {
	fprintf(stderr, "\n\tUsage: homerTools matrix3D <distance matrix file> [options]\n");
	fprintf(stderr, "\n\tOptions for command: cluster\n");
	fprintf(stderr, "\t\t-r (invert the scores i.e. high scores = close together)\n");
	fprintf(stderr, "\t\t-2 (input matrix has two columns of identifiers [to work well with treeview])\n");
	fprintf(stderr, "\t\t-c <cluster assignment file> (color clusters)\n");
	fprintf(stderr, "\t\t-iter <#> (number of iterations)\n");
	fprintf(stderr, "\t\t\n"); 
	exit(0);
}

void matrix3DProgram(int argc, char** argv) {

	fprintf(stderr, "\n");
	if (argc < 3) {
		printCMDmatrix3D();
	}

	srand(time(NULL));

	char* distMatrixFile = argv[2];
	//char* prefix = NULL;
	//int reverseFlag = 0;
	int col2Flag = 0;
	char* clusterFile = NULL;
	double iters = 1e6;

	for (int i=3;i<argc;i++) {
		if (strcmp(argv[i],"-r")==0) {
			//reverseFlag = 1;
		} else if (strcmp(argv[i],"-2")==0) {
			col2Flag = 1;
		} else if (strcmp(argv[i],"-c")==0) {
			clusterFile = argv[++i];
		} else if (strcmp(argv[i],"-iter")==0) {
			sscanf(argv[++i],"%lf",&iters);
		} else {
			fprintf(stderr, "!!! Invalid option: %s\n", argv[i]);
			printCMDcluster();
		}
	}

	int maxIterations = (int)iters; 
	char* buf = new char[BUFFER];
	char** line = new char*[10000];
	int numCols = 0;

	FILE* fp = fopen(distMatrixFile,"r");
	if (fp == NULL) {
		fprintf(stderr, "!! Could not open input file: %s\n", distMatrixFile);
		exit(0);
	}


	int numMatrix = 0;
	double** distMatrix = NULL;
	char** names = NULL;
	double fixedStep = -1;
	
	int lineCount = 0;
	while (fgets(buf, BUFFER, fp) != NULL) {
		split(buf, line, numCols, '\t');
		lineCount++;
		if (lineCount == 1) {
			numMatrix = numCols-1-col2Flag;
			distMatrix = new double*[numMatrix];
			names = new char*[numMatrix];
			for (int i=0;i<numMatrix;i++) {
				distMatrix[i] = new double[numMatrix];
				for (int j=0;j<numMatrix;j++) distMatrix[i][j]=FLT_MAX;
			}
			continue;
		}
		if (numCols < numMatrix+1) {
			fprintf(stderr, "Something is wrong with your distance matrix!!!\n");
			exit(0);
		}
		int index = lineCount -2;
		names[index] = new char[strlen(line[0])+1];
		strcpy(names[index], line[0]);
		for (int i=1+col2Flag;i<numCols;i++) {
			sscanf(line[i], "%lf", &(distMatrix[index][i-1-col2Flag]));
		}
	}




	Genome3D* genomeStruct = new Genome3D();
	genomeStruct->init(distMatrix,numMatrix,names,fixedStep);
	if (clusterFile != NULL)
		genomeStruct->addColorsFromClusters(clusterFile);

	genomeStruct->optimize(maxIterations);
	genomeStruct->print(stdout);


}
