
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


#include "CommandLine.h"

/*==========================================================
|  CommandLine::CommandLine() - Initializes all command line
|    parameters to defaults
\========================*/
CommandLine::CommandLine (int ARGC, char** ARGV) {
	argc = ARGC;
	argv = ARGV;
	algorithm=0;
	
	seqfile = NULL;
	consfile = NULL;
	statfile = NULL;
	merfile = NULL;
	motiffile = NULL;
	outputfile = NULL;
	seedfile = NULL;
	merLengths = NULL;
	numMerLengths =0;
	numGapLengths = 0;
	gapLengths = NULL;
	statmemory = NULL;
	autoScale = 1;
	targetTreeFlag = 1;

	action = NO_ACTION;
	weightFlag = 0;
	revoppFlag = 1;
	offset = 0;
	seqType = DNA_SEQ_TYPE;
	minmer = 10;
	maxmer = 10;
	mismatches = 2;
	numIUPAC = 0;
	dualMotifs = 0;
	allowFlip = 0;
	numSeeds = 50;
	IUPACtype = IUPAC_N;
	numPreCheckSeeds = 200;
	nratiocutoff = 0.9;
	maxOptIterations = 5;
	rmalignseeds = 1;
	consThresh = -1.0;
	exactFlag = 0;
	eraseFlag = 0;
	gapsize = 0;
	gapoffset = 0;
	fisherSize = 2;
	zoopsApprox = 1;
	speed = SPEED_FAST;
	maxPerSeq = 2;
	freqAdjust = 0;
	backFreq = NULL;
	branchSize = 0.5;

	maxNegGenePercentage = 0.5;
	minSeqLen = 0;
	maxSeqLen = 1000000;

	seqLength = NULL;
	totalSeqLength = 0;

	//global slush variables
}

CommandLine::~CommandLine ()
{
}


/*===============================================================
|  int CommandLine::parseInput(int,char**) - This function parses
|    the command line and updates the data members in the 
|    CommandLine object.  This is done so that additional parameters
|    can be added with minimul code changes to the main program
\============================================================*/
void
CommandLine::parseInput () {
	if (argc < 2) printHelp();
	for (int i = 1; i < argc; i++) {
		if (*argv[i] == '-') {
			if (argv[i][1] == '-') {
				printHelp();
			}
			if (strcmp(argv[i],"-s")==0) {
				seqfile = argv[i+1];
			} else if (strcmp(argv[i],"-c")==0) {
				consfile = argv[i+1];
			} else if (strcmp(argv[i],"-g")==0) {
				statfile = argv[i+1];
			} else if (strcmp(argv[i],"-m")==0) {
				motiffile = argv[i+1];
			} else if (strcmp(argv[i],"-mer")==0) {
				merfile = argv[i+1];
			} else if (strcmp(argv[i],"-dna")==0) {
				seqType = DNA_SEQ_TYPE;
			} else if (strcmp(argv[i],"-prot")==0) {
				seqType = PROTEIN_SEQ_TYPE;
			} else if (strcmp(argv[i],"-zoopsapprox")==0) {
				if (strcmp(argv[i+1],"OFF")==0 || strcmp(argv[i+1],"off")==0) {
					zoopsApprox = 0;
				} else if (argv[i+1][0] == '-') {
					zoopsApprox = 0;
				} else {
					zoopsApprox = 1;
					sscanf(argv[i+1],"%d", &maxPerSeq);
				}
			} else if (strcmp(argv[i],"-exact")==0) {
				exactFlag = 1;
			} else if (strcmp(argv[i],"-noautoscale")==0) {
				autoScale = 0;
			} else if (strcmp(argv[i],"-freqAdjust")==0) {
				freqAdjust = 1;
			} else if (strcmp(argv[i],"-dual")==0) {
				dualMotifs = 1;
			} else if (strcmp(argv[i],"-w")==0) {
				weightFlag = 1;
			} else if (strcmp(argv[i],"-T")==0) {
				targetTreeFlag = 0;
			} else if (strcmp(argv[i],"-flip")==0) {
				allowFlip = 1;
			} else if (strcmp(argv[i],"-offset")==0) {
				sscanf(argv[i+1],"%d", &offset);
			} else if (strcmp(argv[i],"-min")==0) {
				sscanf(argv[i+1],"%d", &minmer);
			} else if (strcmp(argv[i],"-max")==0) {
				sscanf(argv[i+1],"%d", &maxmer);
			} else if (strcmp(argv[i],"-mis")==0) {
				sscanf(argv[i+1],"%d", &mismatches);
			} else if (strcmp(argv[i],"-branch")==0) {
				sscanf(argv[i+1],"%f", &branchSize);
				if (branchSize < 1e-6) {
					fprintf(stderr, "Please set the -branch <#> paramater to a higher value (> 1e-6)\n");
					exit(0);
				}
			} else if (strcmp(argv[i],"-IUPAC")==0) {
				sscanf(argv[i+1],"%d", &numIUPAC);
			} else if (strcmp(argv[i],"-iupactype")==0) {
				sscanf(argv[i+1],"%d", &IUPACtype);
				if (IUPACtype < 1 || IUPACtype > 3) {
					fprintf(stderr, "Wrong input for IUPACtype!\n");
					exit(0);
				}
			} else if (strcmp(argv[i],"-N")==0) {
				sscanf(argv[i+1],"%f", &nratiocutoff);
			} else if (strcmp(argv[i],"-maxneg")==0) {
				sscanf(argv[i+1],"%f", &maxNegGenePercentage);
			} else if (strcmp(argv[i],"-S")==0) {
				sscanf(argv[i+1],"%d", &numSeeds);
			} else if (strcmp(argv[i],"-seqless")==0) {
				sscanf(argv[i+1],"%d", &minSeqLen);
			} else if (strcmp(argv[i],"-seqmore")==0) {
				sscanf(argv[i+1],"%d", &maxSeqLen);
			} else if (strcmp(argv[i],"-speed")==0) {
				if (strcmp(argv[i+1],"FAST")==0 || strcmp(argv[i+1],"fast")==0) {
					speed = SPEED_FAST;
				} else if (strcmp(argv[i+1],"NORMAL")==0 || strcmp(argv[i+1],"normal")==0) {
					speed = SPEED_NORMAL;
				}
			} else if (strcmp(argv[i],"-I")==0) {
				sscanf(argv[i+1],"%d", &maxOptIterations);
			} else if (strcmp(argv[i],"-o")==0) {
				outputfile = argv[i+1];
			} else if (strcmp(argv[i],"-rmalign")==0) {
				rmalignseeds = 0;
			} else if (strcmp(argv[i],"-preS")==0) {
				sscanf(argv[i+1],"%d",&numPreCheckSeeds);
			} else if (strcmp(argv[i],"-gap")==0) {
				gapLengths = parseRangeVariable(argv[i+1],numGapLengths);
			} else if (strcmp(argv[i],"-len")==0) {
				merLengths = parseRangeVariable(argv[i+1],numMerLengths);
			} else if (strcmp(argv[i],"-norevopp")==0) {
				revoppFlag = 0;
			} else if (strcmp(argv[i],"-alg")==0) {
				if (strcmp(argv[i+1],"hypergeo")==0) {
					algorithm = ALG_HYPERGEO;
				}
				if (strcmp(argv[i+1],"binomial")==0) {
					algorithm = ALG_BINOMIAL;
				}
				if (strcmp(argv[i+1],"sitehypergeo")==0) {
					algorithm = ALG_SITEHYPERGEO;
				}
				if (strcmp(argv[i+1],"sitebinomial")==0) {
					algorithm = ALG_SITEBINOMIAL;
				}
				if (strcmp(argv[i+1],"rank")==0) {
					algorithm = ALG_RANK;
				}
				if (strcmp(argv[i+1],"freqdiff")==0) {
					algorithm = ALG_FREQDIFF;
				}
				if (strcmp(argv[i+1],"fisher")==0) {
					algorithm = ALG_FISHER;
					sscanf(argv[i+2],"%d",&fisherSize);
				}
			} else if (strcmp(argv[i],"-a")==0) {
				if (strcmp(argv[i+1],"MOTIFS")==0) {
					action = ACTION_MOTIFS;
				} else if (strcmp(argv[i+1],"MERS")==0) {
					action = ACTION_MERS;
				} else if (strcmp(argv[i+1],"FIND")==0) {
					action = ACTION_FIND;
				} else if (strcmp(argv[i+1],"DMERS")==0) {
					action = ACTION_DMERS;
				} else if (strcmp(argv[i+1],"OPTPVALUE")==0) {
					action = ACTION_OPTPVALUE;
				} else if (strcmp(argv[i+1],"GETPVALUE")==0) {
					action = ACTION_GETPVALUE;
				} else if (strcmp(argv[i+1],"REFINE")==0) {
					action = ACTION_REFINE;
				} else if (strcmp(argv[i+1],"REFINETHRESH")==0) {
					action = ACTION_REFINETHRESH;
				} else if (strcmp(argv[i+1],"CLUSTER")==0) {
					action = ACTION_CLUSTER;
				} else if (strcmp(argv[i+1],"REMOVE")==0) {
					action = ACTION_REMOVE;
				} else if (strcmp(argv[i+1],"GENESCORE")==0) {
					action = ACTION_GENESCORE;
				} else if (strcmp(argv[i+1],"SORTMERS")==0) {
					action = ACTION_SORTMERS;
				}
			}
		}
	}
}

int* parseRangeVariable(char* input, int &numLengths) {
//fprintf(stderr, "%s\n", input);
	if (input == NULL) {
		numLengths = 0;
		return NULL;
	}
	numLengths = 1;
	int k=0;
	char* starts = new char[100];
	char* rangeFlags = new char[100];
	starts[0] = 0;
	rangeFlags[0] = 0;
	while (input[k] != '\0') {
		if (input[k] == ',') {
			input[k] = '\0';
			starts[numLengths] = k+1;
			rangeFlags[numLengths] = 0;
			numLengths++;
		} else if (input[k] == '-') {
			input[k] = '\0';
			starts[numLengths] = k+1;
			rangeFlags[numLengths] = 1;
			numLengths++;
		}
		k++;
	}
	int numNeeded = 0;
	int* lengths = new int[numLengths];
	for (k=0;k<numLengths;k++) {
		sscanf(&(input[(int)starts[k]]), "%d", &(lengths[k]));
		numNeeded++;
		if (rangeFlags[k] == 1) {
			numNeeded += lengths[k]-lengths[k-1];
		}
	}
	int * tmp = new int[numNeeded];
	int index = 0;
	for (k=0;k<numLengths;k++) {
		if (rangeFlags[k] == 0) {
			tmp[index++] = lengths[k];
		} else {
			for (int m=lengths[k-1]+1;m<=lengths[k];m++) {
				tmp[index++] = m;
			}
		}
	}
	delete []lengths;
	lengths = tmp;
	numLengths = index;
	delete []starts;
	delete []rangeFlags;
	return lengths;
}

int commandlinesortint(const void* a, const void* b) {
	if (*((int*)a) < *((int*)b)) {
		return -1;
	} else if (*((int*)a) > *((int*)b)) {
		return 1;
	} else {
		return 0;
	}
}

char* makeVariableRangeString(int* lengths, int numLengths) {
	char* str = NULL;
	if (numLengths == 0 || lengths == NULL) {
		str = new char[1];
		str[0] = '\0';
		return str;
	}
	str = new char[numLengths*10]; //make it plenty big
	char* curpos = str;
	qsort(lengths, numLengths, sizeof(int), commandlinesortint);
	int cflag = 0;
	int end = 0;
	//int start = 0;
	int z = 0;
	for (int i=0;i<numLengths;i++) {
		if (i==0) {
			sprintf(curpos, "%d", lengths[i]);
			z=0;
			while (curpos[z] != '\0') {z++;}
			curpos = &(curpos[z]);
			//start = lengths[i];
			continue;
		} else {
			if (lengths[i] == lengths[i-1]+1) {
				cflag = 1;
				end = lengths[i];
				if (i+1<numLengths) {
					continue;
				}
			}
			if (cflag && ((lengths[i] != lengths[i-1]+1) || i+1==numLengths)) {
				sprintf(curpos, "-%d",end);
				z=0;
				while (curpos[z] != '\0') {z++;}
				curpos = &(curpos[z]);
				if (i<numLengths-1) {
					sprintf(curpos, ",");
					z=0;
					while (curpos[z] != '\0') {z++;}
					curpos = &(curpos[z]);
				}
				cflag = 0;
			} else if (cflag == 0) { // && i+1!=numLengths) {
				sprintf(curpos, ",%d",lengths[i]);
				z=0;
				while (curpos[z] != '\0') {z++;}
				curpos = &(curpos[z]);
			}
				
		}
	}
	return str;
}


void CommandLine::checkDependencies() {

	if (action == NO_ACTION) {
		fprintf(stderr, "Need to specify an action -> -a [action]\n");
		exit(0);
	}
	if (action == ACTION_FIND || action == ACTION_OPTPVALUE 
				|| action == ACTION_GENESCORE) {
		if (motiffile == NULL) {
			fprintf(stderr, "Need to specify a motif file [-m]!\n");
			exit(0);
		}
	}
	if (outputfile == NULL) {
	
		outputfile = new char[4];
		strcpy(outputfile,"OUT");
	}
	if (seqType == PROTEIN_SEQ_TYPE) {
		revoppFlag = 0;
	}
	if (numGapLengths == 0) {
		gapLengths = new int[1];
		gapLengths[0] = 0;
		numGapLengths = 1;
	}
	//if (algorithm == ALG_BINOMIAL && exactFlag != 1) {
		//exactFlag = 1;
		//fprintf(stderr, "Need to use exact flag when using Bionomial algorithm\n");
	//}
	if (algorithm == ALG_RANK) {
		exactFlag = 1;
		//need to use exact flag for ranking alogirhtm
	}
	if (action == ACTION_MOTIFS) {	
		if (seqfile == NULL) {
			fprintf(stderr, "Need Sequence File\n");
			exit(0);
		}
		if (statfile == NULL) {
			fprintf(stderr, "Need Group File\n");
			exit(0);
		}
	}
	if (gapsize > 0) {
		maxmer = minmer;
		minmer =(int) (2*floor(((float)minmer))/2);
	}
	if (numMerLengths == 0) {
		numMerLengths = maxmer-minmer+1;
		merLengths = new int[numMerLengths];
		int j=0;
		for (int i=minmer;i<=maxmer;i++) {
			merLengths[j++] = i;
		}
	}

	if (action == ACTION_MOTIFS) {
		fprintf(stderr, "Finding motifs with lengths of %s\n", 
			makeVariableRangeString(merLengths,numMerLengths));
		fprintf(stderr, "Finding motifs with gaps of length %s\n", 
			makeVariableRangeString(gapLengths,numGapLengths));
		fprintf(stderr, "# of mismatches in glabal search phase: %d\n", mismatches); 
		fprintf(stderr, "Number of IUPAC symbols: %d at depth of %d\n", numIUPAC, IUPACtype);
	}
}

void CommandLine::printHelp () {
  fprintf (stderr, "%s : Empirical Motif Optimizer\n", argv[0]);
  fprintf (stderr, "usage: ./%s [data] [parameters] -a [action]\n", argv[0]);
  fprintf (stderr, "This program is meant to be called from other programs (i.e. findMotifsGenome.pl), and not used directly\n");
  fprintf (stderr, "Data options:\n");
  fprintf (stderr, "\t-dna|-prot : Sequence type (-dna)\n");
  fprintf (stderr, "\t-s <file> : Sequence File\n");
  fprintf (stderr, "\t-g <file> : Group/Stat File\n");
  fprintf (stderr, "\t-mer <file> : Mer File\n");
  fprintf (stderr, "\t-m <file> : PSSM Motif File\n");
  fprintf (stderr, "\t-o <file> : output file prefix\n");
  fprintf (stderr, "\t-seed <file> : seed file\n");
  fprintf (stderr, "\t-offset <#> : offset of sequence from TSS (-2000)\n");
  fprintf (stderr, "Parameter options:\n");
  fprintf (stderr, "\t-exact : remember mapping between mers and genes (default: approx)\n");
  fprintf (stderr, "\t-w : Weight sequences (according to addition columns in group file: 1st-gene 2nd-sequence)\n");
  fprintf (stderr, "\t-T : Test all sequences as candidate motifs (default: only test target sequences)\n");
  fprintf (stderr, "\t-noautoscale : Do not autoscale sequences to be equal in foreground and background\n");
  fprintf (stderr, "\t-freqAdjust : Compute log-odds using frequency, default (0.25)\n");
  fprintf (stderr, "\t-dual : find dual motifs in the form A<gap>B where A and B can be rev-opposites\n");
  fprintf (stderr, "\t-flip : find dual motifs in the form A<gap>B or B<gap>A\n");
  fprintf (stderr, "\t-zoopsapprox <OFF,#(max to count)> : (counts multiple motifs per sequence | default: 2)\n");
  fprintf (stderr, "\t-norevopp : don't search opposite strand (default->DNA:yes, Protein:no)\n");
  fprintf (stderr, "\t-min <#> : min mer size (6)\n");
  fprintf (stderr, "\t-max <#> : max mer size [also standard mer size] (10)\n");
//  fprintf (stderr, "\t-len <#,#,#-#> : check mers of size # (i.e. -len 6-8,10,12) (default=10\n");
  fprintf (stderr, "\t-len <#> : Find motifs of length # (default=10)\n");
  fprintf (stderr, "\t-gap <#,#,#-#> : Find motifs with gaps(0)(i.e. -gap 3 -gap 2,4,5 -gap 1-10\n");
  fprintf (stderr, "\t    Gaps will only be in the center of motif and will only use even lengthed motifs\n");
  fprintf (stderr, "\t-mis <#> : # of mismatches to check for degeneracy (1)\n");
  fprintf (stderr, "\t-IUPAC <#> : # of IUPAC codes per mer that can be used during global optimization (0)\n");
  fprintf (stderr, "\t-iupactype <1,2,or3> : Type of IUPAC symbols used\n");
  fprintf (stderr, "\t\t1: (default) Only N is used\n");
  fprintf (stderr, "\t\t2: Only N and 2 bp symbols are used (i.e. R = A or G\n");
  fprintf (stderr, "\t\t3: Full IUPAC code is used (includes 3-way symbols)\n");
  fprintf (stderr, "\t-S <#> : number of seeds to check during profile optimization(50)\n");
  fprintf (stderr, "\t-branch <#> : sets depth of optimization (closer to zero the more sensitive (0.5))\n");
  fprintf (stderr, "\t-I <#> : maximum number of iterations during optimization (5)\n");
  fprintf (stderr, "\t-rmalign : DO NOT remove aligned seeds\n");
  fprintf (stderr, "\t-maxneg <0 to 1> maximum percentage of negative genes that can contain the motif\n");
  fprintf (stderr, "\t-speed <NORMAL|FAST>: Program will heuristically avoid performing exhaustive\n");
  fprintf (stderr, "\t    calculations (default: FAST)\n");
  fprintf (stderr, "Scoring Functions:\n");
  fprintf (stderr, "\t-alg <method> : scoring algorithm (default: hypergeo)\n");
  fprintf (stderr, "\t    hypergeo - hypergeometric scoring (ZOOPS)\n");
  fprintf (stderr, "\t    binomial - binomical scoring [for variable length seq] (ZOOPS) (requires exact)\n");
  fprintf (stderr, "\t    approxbinomial - binomical scoring [for variable length seq] (ZOOPS) (requires exact)\n");
  fprintf (stderr, "\t    sitehypergeo - hypergeometric scoring across seq positions (very slow)\n");
  fprintf (stderr, "\t    sitebinomial - binomial scoring across seq positions\n");
  fprintf (stderr, "\t    fisher <#> - fisher exact test (slow, # scales exponentially)\n");
  fprintf (stderr, "\t      <# = largest repetition to consider [default=2]>\n");
  fprintf (stderr, "\t    rank - group file must have sortable numeric value\n");
  fprintf (stderr, "\t    freqdiff - used by most bayesian/nnet programs\n");
  fprintf (stderr, "\t    logit - used by most bayesian/nnet programs\n");
  fprintf (stderr, "Background Modeling options (this forces a binomial style scoring function):\n");
  fprintf (stderr, "\t-b <method> [method options...]\n");
  fprintf (stderr, "\t    markov <#> - generate hmm from target sequences using a hmm of order #\n");
  fprintf (stderr, "\t    bmarkov <#> - generate hmm from background sequences using a hmm of order #\n");
  fprintf (stderr, "\t    mosaic - generate mosaic hmm from background sequences **coming soon**\n");
  fprintf (stderr, "Filter Options:\n");
  fprintf (stderr, "\t-N <float> : filtering cutoff for ratio of N's in sequence (0.9)\n");
  fprintf (stderr, "\t-seqless <#> : filter sequences shorter than #\n");
  fprintf (stderr, "\t-seqmore <#> : filter sequences longer than #\n");
  fprintf (stderr, "Actions (-a):\n");
  fprintf (stderr, "\tMOTIFS - Find motifs <outfile>.motifs# where # = motif length\n");
  fprintf (stderr, "\tMERS - Create mer file (low memory) <stdout>\n");
  fprintf (stderr, "\tDMERS - Create degenerate mer file <stdout>\n");
  fprintf (stderr, "\tFIND - find motifs in sequence <stdout>\n");
  fprintf (stderr, "\tOPTPVALUE - optimize motif threshold and pvalue (exact)<stdout>\n");
  fprintf (stderr, "\tGETPVALUE - get the p-value enrichment for a given motif(exact)<stdout>\n");
  fprintf (stderr, "\tGENESCORE - returns highest motif score for each gene <stdout>\n");
  fprintf (stderr, "\tREFINE - optimize motif PSSM profile, threshold, and pvalue <stdout>\n");
  fprintf (stderr, "\tREFINETHRESH - optimize motif PSSM threshold and pvalue <stdout>\n");
  fprintf (stderr, "\tCLUSTER - cluster mers from seed file (can't use exact scoring) <outfile>\n");
  fprintf (stderr, "\tSORTMERS - sort a mer file according to pvalue <stdout>\n");
  fprintf (stderr, "\tREMOVE - removes motif from sequence (replaces with N's) <stdout>\n");
  fprintf (stderr, "\n");
  fprintf (stderr, "This program is meant to be called from other programs (i.e. findMotifsGenome.pl), and not used directly\n");
  fprintf (stderr, "\n");
  exit (0);
}
