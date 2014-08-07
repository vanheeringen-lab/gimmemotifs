
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

#include "Motif2.h"

int Motif::alphaSize = DNA_ALPHA;
char* Motif::alphaIndex = NULL;

void printCMDhomer() {
	fprintf(stderr, "\n\thomer2 <command> [options]\n");
	fprintf(stderr, "\n\tCommands:\n");
	fprintf(stderr, "\t\tdenovo (discover motifs using de novo motif analysis)\n");
	fprintf(stderr, "\t\tknown (score known motifs for enrichment)\n");
	fprintf(stderr, "\t\tfind (find specific instances motifs in the data)\n");
	fprintf(stderr, "\t\tnorm (normalize background sequences)\n");
	fprintf(stderr, "\t\tmask (mask instances of motifs)\n");
	fprintf(stderr, "\n");
	exit(0);
}

void printCMDdenovo() {
	fprintf(stderr, "\n\thomer2 denovo -s <seq file> -g <group file> [options]\n");
	fprintf(stderr, "\n\tDiscover motif de novo in a set of sequences.\n");
	fprintf(stderr, "\n\tPrimary Inputs:\n");
	fprintf(stderr, "\t\t-s <sequence file> (tab delimited sequence file)\n");
	fprintf(stderr, "\t\t-g <group file> (sequence group and weight assignments)\n");
	fprintf(stderr, "\t\t\t--- or ---\n");
	fprintf(stderr, "\t\t-i <input FASTA file>\n");
	fprintf(stderr, "\t\t-b <background FASTA file>\n");
	fprintf(stderr, "\n\tOptions:\n");
	fprintf(stderr, "\t\t-o <output motif file> (default: sent to stdout)\n");
	//fprintf(stderr, "\t\t-sites <predicted sites file> (default: none)\n");
	fprintf(stderr, "\t\t-len <#> (length of motif to search for, default: 10)\n");
	fprintf(stderr, "\t\t-mis <#> (maximum number of mismatches in global search phase, default: 2)\n");
	fprintf(stderr, "\t\t-strand <+|-|both> (search for motifs on specific strand, default: both)\n");
	fprintf(stderr, "\t\t-stat <hypergeo|binomial> (enrichment statistic, default: binomial)\n");
	//fprintf(stderr, "\t\t-stat <hypergeo|binomial|zscore> (enrichment statistic, default: binomial)\n");
	fprintf(stderr, "\t\t-S <#> (Total motifs to find, default: %d)\n",MOTIF_DEFAULT_NUMBER_OF_MOTIFS);
	fprintf(stderr, "\t\t-olen <#> (length of lower-order oligos to normalize in oligo table)\n");
	fprintf(stderr, "\t\t\t-oout <filename> (output normalization weights to file)\n");
	fprintf(stderr, "\t\t\t-omax <#> (max oligo normalization iterations, default: %d)\n", 
																MAX_OLIGO_NORMALIZATION_ITERATIONS);
	fprintf(stderr, "\t\t-p <#> (Number of processers to use, default: 1)\n");
	fprintf(stderr, "\n\tAdvanced Options:\n");
	fprintf(stderr, "\t\t-nozoops (skip proper zoops scoring at end, default: -zoops)\n");
	fprintf(stderr, "\t\t-tmp <filename> (temporary results file, default: \".tmp.motifs\")\n");
	fprintf(stderr, "\t\t-oligos <filename> (print enrichment of individual oligos)\n");
	fprintf(stderr, "\t\t-opt <motif file> (expand/futher optimize these motifs/skip global phase)\n");

	fprintf(stderr, "\n\t\tSpeed vs. Sensitivity:\n");
	fprintf(stderr, "\t\t-fullMask | -quickMask (choose one:)\n");
	fprintf(stderr, "\t\t\t-fullMask (as motifs are found, mask them from original sequences.\n");
	fprintf(stderr, "\t\t\t\tRequires more memory but yields \"cleaner\" results, default)\n");
	fprintf(stderr, "\t\t\t-quickMask (as motifs are found, mask bound \"oligos\" only, old way)\n");
	fprintf(stderr, "\t\t-e <#> (maximum expected motif instances per bp, default: %.3lf)\n",
																	MOTIF_MAX_HIT_RATE_TO_CONSIDER);
	fprintf(stderr, "\t\t-T <#> (number of trial matrices when optimizing, default: %d)\n",
																	MOTIFS_NUM_POTENTIAL_MOTIFS);
	fprintf(stderr, "\t\t-blen <#> (# of bp on either side to check for redundancy, default: 1)\n");
	fprintf(stderr, "\t\t-maxBack <#> (Max percentage of background that motifs may contain, default: 0.5)\n");
	fprintf(stderr, "\t\t-minlp <#> (minimum significance of seeds to optimize, default: %.3lf)\n",
																	MOTIF_SEED_MINIMUM_LOGP);

	fprintf(stderr, "\n\t\tMemory vs. Speed:\n");
	fprintf(stderr, "\t\t-cache <#> (size in MB of stat cache, helps for hypergeo, default: %d)\n",
																	MAX_STAT_CACHE_SIZE);
																	
	fprintf(stderr, "\n");
	exit(0);
}

void programDeNovo(int argc, char** argv) {

	DeNovoMotifs* denovo = new DeNovoMotifs();


	FILE* oligoOutputFile = NULL;
	FILE* outputFp = NULL;
	char* outputFileName = NULL;
	
	if (argc < 3) {
		printCMDdenovo();
	}	
	for (int i=2;i<argc;i++) {
		if (strcmp(argv[i],"-s")==0) {
			denovo->seqFile = argv[++i];
		} else if (strcmp(argv[i],"-g")==0) {
			denovo->groupFile = argv[++i];
		} else if (strcmp(argv[i],"-i")==0) {
			denovo->inputFASTAfile = argv[++i];
		} else if (strcmp(argv[i],"-b")==0) {
			denovo->bgFASTAfile = argv[++i];
		} else if (strcmp(argv[i],"-sites")==0) {
			denovo->sitefp = fopen(argv[++i],"w");
			if (denovo->sitefp == NULL) {
				fprintf(stderr, "Could not open %s for writing sites!!!\n",argv[i]);
				exit(0);
			}
		} else if (strcmp(argv[i],"-len")==0) {
			sscanf(argv[++i],"%d", &(denovo->oligoLength));
		} else if (strcmp(argv[i],"-olen")==0) {
			sscanf(argv[++i],"%d", &(denovo->normLength));
		} else if (strcmp(argv[i],"-omax")==0) {
			sscanf(argv[++i],"%d", &(denovo->maxAutoNormIters));
		} else if (strcmp(argv[i],"-T")==0) {
			sscanf(argv[++i],"%d", &(denovo->numTrialMotifs));
		} else if (strcmp(argv[i],"-minlp")==0) {
			sscanf(argv[++i],"%lf", &(denovo->minimumSeedLogp));
		} else if (strcmp(argv[i],"-p")==0) {
			sscanf(argv[++i],"%d", &(denovo->numCPUs));
		} else if (strcmp(argv[i],"-blen")==0) {
			sscanf(argv[++i],"%d", &(denovo->bufferLength));
		} else if (strcmp(argv[i],"-cache")==0) {
			unsigned long long int sizeMB = 0;
			sscanf(argv[++i],"%lld", &sizeMB);
			denovo->maxCacheSize = sizeMB*1000000/sizeof(double);
		} else if (strcmp(argv[i],"-zoops")==0) {
			denovo->finalZoopsFlag=1;
		} else if (strcmp(argv[i],"-nozoops")==0) {
			denovo->finalZoopsFlag=0;
		} else if (strcmp(argv[i],"-rand")==0) {
			denovo->randFlag=1;
		} else if (strcmp(argv[i],"-o")==0) {
			outputFileName=argv[++i];
			outputFp=fopen(outputFileName,"w");
			if (outputFp==NULL) {
				fprintf(stderr, "!!! Could not open %s for writing output!!!\n",outputFileName);
				exit(1);
			}
		} else if (strcmp(argv[i],"-oout")==0) {
			denovo->normFile = argv[++i];
		} else if (strcmp(argv[i],"-tmp")==0) {
			denovo->tmpFile = argv[++i];
		} else if (strcmp(argv[i],"-opt")==0) {
			denovo->optimizeGivenFlag=1;
			fprintf(stderr, "\tOptimizing given motifs...\n");
			denovo->motifs = Motif::readMotifFile(argv[++i],denovo->numMotifs);
		} else if (strcmp(argv[i],"-mis")==0) {
			sscanf(argv[++i],"%d", &(denovo->maxMisMatches));
		} else if (strcmp(argv[i],"-maxBack")==0) {
			sscanf(argv[++i],"%lf", &(denovo->maxBackgroundPercent));
		} else if (strcmp(argv[i],"-fullMask")==0) {
			denovo->maskMethod = DENOVO_MASK_FULL;
		} else if (strcmp(argv[i],"-quickMask")==0) {
			denovo->maskMethod = DENOVO_MASK_QUICK;
		} else if (strcmp(argv[i],"-e")==0) {
			sscanf(argv[++i],"%lf", &(denovo->maximumExpectedPerBp));
		} else if (strcmp(argv[i],"-S")==0) {
			sscanf(argv[++i],"%d", &(denovo->totalMotifsToFind));
		} else if (strcmp(argv[i],"-oligos")==0) {
			oligoOutputFile = fopen(argv[++i],"w");
			if (oligoOutputFile == NULL) {
				fprintf(stderr, "!!! Could not open oligo output file (%s)!!!\n", argv[i]);
				exit(0);
			}
		} else if (strcmp(argv[i],"-speed")==0) {
			i++;
			if (strcmp(argv[i],"auto")==0) {
				denovo->speedFlag = DENOVO_SPEED_AUTOADJUST;
			} else if (strcmp(argv[i],"constant")==0) {
				denovo->speedFlag = DENOVO_SPEED_CONSTANT;
			} else {
				fprintf(stderr, "!!! option \"-speed %s\" not valid !!!\n", argv[i]);
				exit(0);
			}
		} else if (strcmp(argv[i],"-stat")==0) {
			i++;
			if (strcmp(argv[i],"hypergeo")==0) {
				denovo->scoringMethod = MOTIF_SCORING_HYPERGEOMETRIC;
			} else if (strcmp(argv[i],"binomial")==0) {
				denovo->scoringMethod = MOTIF_SCORING_BINOMIAL;
			} else if (strcmp(argv[i],"zscore")==0) {
				denovo->scoringMethod = MOTIF_SCORING_ZSCORE;
			} else {
				fprintf(stderr, "!!! option \"-stat %s\" not valid !!!\n", argv[i]);
				exit(0);
			}
		} else if (strcmp(argv[i],"-strand")==0) {
			i++;
			if (strcmp(argv[i],"+")==0) {
				denovo->strand = MOTIF_STRAND_POS;
			} else if (strcmp(argv[i],"-")==0) {
				denovo->strand = MOTIF_STRAND_NEG;
			} else if (strcmp(argv[i],"both")==0) {
				denovo->strand = MOTIF_STRAND_BOTH;
			} else {
				fprintf(stderr, "!!! Could not recognize -strand \"%s\" !!!\n", argv[i]);
				printCMDdenovo();
			}
		} else {
			printCMDdenovo();
		}
	}
	fprintf(stderr, "\n");
	
	if (outputFileName == NULL) outputFp = stdout;

	denovo->checkParameters();
	if (denovo->fastaFlag) {
		denovo->parseFasta();
	}

	denovo->readSequence2Tree();
	denovo->normalizeOligos();

	//setup for motif finding
	denovo->initializeEnrichmentScoring();

	if (denovo->optimizeGivenFlag == 0) {
		//normal operation
		//Global optimization
		denovo->oligoArray->initializeActiveOligoArray();
		denovo->oligoArray->initializeOligoPvalues(denovo->scorer);
		denovo->findMisMatchEnrichment(denovo->maxMisMatches, oligoOutputFile,HOMER_VERBOSE_LEVEL1);	
		//local optimization
		if (denovo->maskMethod == DENOVO_MASK_FULL) {
			denovo->loadFullSequences();
		}
		denovo->optimizeSeeds2Motifs();
	} else {
		for (int i=0;i<denovo->numMotifs;i++) denovo->motifs[i]->adjustLength(denovo->oligoLength);
		denovo->optimizeGivenMotifs();
	}

	denovo->freeOligos();

	if (denovo->finalZoopsFlag && denovo->maskMethod == DENOVO_MASK_QUICK) {
		fprintf(stderr, "\n\tFinalizing enrichment with proper zoops scoring...\n");
		denovo->scoreWithKnownMotifs(KNOWN_OPTIMIZE_THREASHOLD);
	} else {
		fprintf(stderr, "\n\tFinalizing Enrichment Statistics (new in v3.4)\n");
		denovo->scoreWithKnownMotifs(KNOWN_JUST_SCORE_ZOOPS);
	}
	denovo->printMotifs(outputFp);

	if (denovo->tmpFile != NULL) {
		(void)unlink(denovo->tmpFile);
	} else {	
		(void)unlink(".tmp.motifs");
	}
	if (outputFileName != NULL) {
		fclose(outputFp);	
		fprintf(stderr, "\tOutput in file: %s\n", outputFileName);
	}
	fprintf(stderr, "\n");


	if (denovo->fastaFlag) {
		fprintf(stderr, "\tCleaning up temporary sequence and group files:\n");
		fprintf(stderr, "\t\t%s\n",denovo->seqFile);
		fprintf(stderr, "\t\t%s\n",denovo->groupFile);
		(void)unlink(denovo->seqFile);
		(void)unlink(denovo->groupFile);
	}

	delete denovo;

}

void printCMDknown() {
	fprintf(stderr, "\n\thomer2 known -s <seq file> -g <group file> -k <known motif file> [options]\n");
	fprintf(stderr, "\n\tFind the enrichment of known motifs in a set of sequences.\n");
	fprintf(stderr, "\n\tPrimary Inputs:\n");
	fprintf(stderr, "\t\t-s <sequence file> (tab delimited sequence file)\n");
	fprintf(stderr, "\t\t-g <group file> (sequence group and weight assignments)\n");
	fprintf(stderr, "\t\t\t--- or ---\n");
	fprintf(stderr, "\t\t-i <input FASTA file>\n");
	fprintf(stderr, "\t\t-b <background FASTA file>\n");
	fprintf(stderr, "\n\tOptions:\n");
	fprintf(stderr, "\t\t-o <output enrichment file> (default: sent to stdout)\n");
	fprintf(stderr, "\t\t-m <known motif file> (Known motif file to check enrichment for)\n");
	fprintf(stderr, "\t\t\t-mout <output motif file> (Output updated motifs with statistics)\n");
	fprintf(stderr, "\t\t-strand <+|-|both> (search for motifs on specific strand, default: both)\n");
	fprintf(stderr, "\t\t-stat <hypergeo|binomial> (enrichment statistic, default: binomial)\n");
	fprintf(stderr, "\t\t-nlen <#> (length of lower-order oligos to normalize, default: 0)\n");
	fprintf(stderr, "\t\t\t-nout <filename> (output normalization weights to file)\n");
	fprintf(stderr, "\t\t\t-nmax <#> (max normalization iterations, default: %d\n", 
																MAX_OLIGO_NORMALIZATION_ITERATIONS);
	fprintf(stderr, "\t\t-cache <#> (size in MB of stat cache: %d)\n",MAX_STAT_CACHE_SIZE);
	fprintf(stderr, "\t\t-p <#> (Number of processers to use, default: 1)\n");
	fprintf(stderr, "\t\t-opt (Optimize degeneracy threshold to get best enrichment, use -mout to get motifs)\n");
	fprintf(stderr, "\t\t-siteReduce <#> (Eliminate redundant motifs sharing > %% of sites)\n");
	fprintf(stderr, "\t\t-maxBack <#> (Max percentage of background that motifs may contain, default: 0.5)\n");
	fprintf(stderr, "\n");
	exit(0);
}
void printCMDfind() {
	fprintf(stderr, "\n\thomer2 find -s <seq file> -m <motif file> [options]\n");
	fprintf(stderr, "\n\tFind the instances of motif(s) in a set of sequences.\n");
	fprintf(stderr, "\n\tPrimary Inputs:\n");
	fprintf(stderr, "\t\t-s <sequence file> (tab delimited sequence file)\n");
	fprintf(stderr, "\t\t\t--- or ---\n");
	fprintf(stderr, "\t\t-i <input FASTA file>\n");
	fprintf(stderr, "\n\tOptions:\n");
	fprintf(stderr, "\t\t-o <output file> (default: sent to stdout)\n");
	fprintf(stderr, "\t\t-m <motif file> (Motif(s) to find instances of)\n");
	fprintf(stderr, "\t\t-offset <#> (offset to report motif instances from, default: midpoint)\n");
	fprintf(stderr, "\t\t-strand <+|-|both> (search for motifs on specific strand, default: both)\n");
	fprintf(stderr, "\t\t-p <#> (Number of processers to use, default: 1)\n");
	fprintf(stderr, "\t\t-mscore (instead of reporting sites, report best motif score per sequence)\n");
	fprintf(stderr, "\n");
	exit(0);
}
void printCMDmask() {
	fprintf(stderr, "\n\thomer2 mask -s <seq file> -m <motif file> [options]\n");
	fprintf(stderr, "\n\tRemove instances of motif(s) in a set of sequences.\n");
	fprintf(stderr, "\n\tPrimary Inputs:\n");
	fprintf(stderr, "\t\t-s <sequence file> (tab delimited sequence file)\n");
	fprintf(stderr, "\t\t\t--- or ---\n");
	fprintf(stderr, "\t\t-i <input FASTA file>\n");
	fprintf(stderr, "\n\tOptions:\n");
	fprintf(stderr, "\t\t-o <output tsv sequence file> (default: sent to stdout)\n");
	fprintf(stderr, "\t\t-m <motif file> (Motif(s) to find instances of)\n");
	fprintf(stderr, "\t\t-strand <+|-|both> (search for motifs on specific strand, default: both)\n");
	fprintf(stderr, "\t\t-p <#> (Number of processers to use, default: 1)\n");
	fprintf(stderr, "\n");
	exit(0);
}
void printCMDnorm() {
	fprintf(stderr, "\n\thomer2 norm -s <seq file> -g <group file> [options]\n");
	fprintf(stderr, "\n\tNormalize background sequences to remove short oligo enrichment\n");
	fprintf(stderr, "\n\tCreates new group file to use with motif finding\n");
	fprintf(stderr, "\n\tPrimary Inputs:\n");
	fprintf(stderr, "\t\t-s <sequence file> (tab delimited sequence file)\n");
	fprintf(stderr, "\t\t-g <group file> (sequence group and weight assignments)\n");
	fprintf(stderr, "\t\t\t--- or ---\n");
	fprintf(stderr, "\t\t-i <input FASTA file>\n");
	fprintf(stderr, "\t\t-b <background FASTA file>\n");
	fprintf(stderr, "\t\t\tOutput will be sequence and group file, even with FASTA input files\n");
	fprintf(stderr, "\n\tOptions:\n");
	fprintf(stderr, "\t\t-o <output weighted group file> (default: sent to stdout)\n");
	fprintf(stderr, "\t\t-s <sequence file> (tab delimited sequence file)\n");
	fprintf(stderr, "\t\t-g <group file> (sequence group and weight assignments)\n");
	fprintf(stderr, "\t\t-strand <+|-|both> (search for motifs on specific strand, default: both)\n");
	fprintf(stderr, "\t\t-nlen <#> (length of lower-order oligos to normalize, default: 0)\n");
	fprintf(stderr, "\t\t\t-nout <filename> (output normalization weights to file)\n");
	fprintf(stderr, "\t\t\t-nmax <#> (max normalization iterations, default: %d)\n", 
																MAX_OLIGO_NORMALIZATION_ITERATIONS);
	fprintf(stderr, "\t\t-neutral (set target/background to neutral i.e. 25%%, 6.25%%, etc. frequencies)\n");
	fprintf(stderr, "\t\t-p <#> (Number of processers to use, default: 1)\n");
	fprintf(stderr, "\n");
	exit(0);
}

void programKnown(int argc, char** argv) {

	KnownMotifs* known = new KnownMotifs();

	FILE* outputFp = NULL;
	char* outputFileName = NULL;

	known->findFlag = 0;
	if (strcmp(argv[1],"find")==0) {
		known->findFlag =1;
		known->knownFlag =0;
	}
	known->maskFlag = 0;
	if (strcmp(argv[1],"mask")==0) {
		known->maskFlag =1;
		known->knownFlag =0;
	}
	known->normOnlyFlag = 0;
	if (strcmp(argv[1],"norm")==0) {
		known->normOnlyFlag =1;
		known->knownFlag =0;
	}
	
	if (argc < 3) {
		if (known->maskFlag) printCMDmask();
		if (known->findFlag) printCMDfind();
		if (known->normOnlyFlag) printCMDnorm();
		printCMDknown();
	}	
	for (int i=2;i<argc;i++) {
		if (strcmp(argv[i],"-s")==0) {
			known->seqFile = argv[++i];
		} else if (strcmp(argv[i],"-g")==0) {
			known->groupFile = argv[++i];
		} else if (strcmp(argv[i],"-i")==0) {
			known->inputFASTAfile = argv[++i];
		} else if (strcmp(argv[i],"-b")==0) {
			known->bgFASTAfile = argv[++i];
		} else if (strcmp(argv[i],"-m")==0) {
			known->knownFile = argv[++i];
		} else if (strcmp(argv[i],"-p")==0) {
			sscanf(argv[++i],"%d", &(known->numCPUs));
		} else if (strcmp(argv[i],"-opt")==0) {
			known->threasholdOptimizationFlag = KNOWN_OPTIMIZE_THREASHOLD;
		} else if (strcmp(argv[i],"-mout")==0) {
			known->motifFilename = argv[++i];
		} else if (strcmp(argv[i],"-mscore")==0) {
			known->mscoreFlag = 1;
		} else if (strcmp(argv[i],"-neutral")==0) {
			known->neutralFlag = 1;
		} else if (strcmp(argv[i],"-siteReduce")==0) {
			sscanf(argv[++i],"%lf", &(known->siteReduceThreshold));
		} else if (strcmp(argv[i],"-maxBack")==0) {
			sscanf(argv[++i],"%lf", &(known->maxBackgroundPercent));
		} else if (strcmp(argv[i],"-offset")==0) {
			sscanf(argv[++i],"%d", &(known->offset));
		} else if (strcmp(argv[i],"-nlen")==0) {
			sscanf(argv[++i],"%d", &(known->normLength));
		} else if (strcmp(argv[i],"-nmax")==0) {
			sscanf(argv[++i],"%d", &(known->maxAutoNormIters));
		} else if (strcmp(argv[i],"-cache")==0) {
			unsigned long long int sizeMB = 0;
			sscanf(argv[++i],"%lld", &sizeMB);
			known->maxCacheSize = sizeMB*1000000/sizeof(double);
		} else if (strcmp(argv[i],"-o")==0) {
			outputFileName=argv[++i];
			outputFp=fopen(outputFileName,"w");
			if (outputFp==NULL) {
				fprintf(stderr, "!!! Could not open %s for writing output!!!\n",outputFileName);
				exit(1);
			}
		} else if (strcmp(argv[i],"-nout")==0) {
			known->normFile = argv[++i];
		} else if (strcmp(argv[i],"-stat")==0) {
			i++;
			if (strcmp(argv[i],"hypergeo")==0) {
				known->scoringMethod = MOTIF_SCORING_HYPERGEOMETRIC;
			} else if (strcmp(argv[i],"binomial")==0) {
				known->scoringMethod = MOTIF_SCORING_BINOMIAL;
			}
		} else if (strcmp(argv[i],"-strand")==0) {
			i++;
			if (strcmp(argv[i],"+")==0) {
				known->strand = MOTIF_STRAND_POS;
			} else if (strcmp(argv[i],"-")==0) {
				known->strand = MOTIF_STRAND_NEG;
			} else if (strcmp(argv[i],"both")==0) {
				known->strand = MOTIF_STRAND_BOTH;
			} else {
				fprintf(stderr, "!!! Could not recognize -strand \"%s\" !!!\n", argv[i]);
				if (known->maskFlag) printCMDmask();
				if (known->findFlag) printCMDfind();
				if (known->normOnlyFlag) printCMDnorm();
				printCMDknown();
			}
		} else {
			if (known->maskFlag) printCMDmask();
			if (known->findFlag) printCMDfind();
			if (known->normOnlyFlag) printCMDnorm();
			printCMDknown();
		}
	}
	fprintf(stderr, "\n");

	if (outputFileName == NULL) outputFp = stdout;
	if (known->findFlag) known->sitefp = outputFp;
	if (known->normOnlyFlag) known->normfp = outputFp;

	known->checkParameters();
	if (known->fastaFlag) {
		known->parseFasta();
	}
	known->loadSequence();

	if (known->normOnlyFlag == 0) {
		known->loadMotifs();
		known->getMotifEnrichment();
		if (known->knownFlag) {
			known->printMotifEnrichment(outputFp);
		}
		if (known->siteReduceThreshold > 0) known->findOverlappingMotifs();

		if (known->motifFilename != NULL) {
			FILE* mfp = fopen(known->motifFilename,"w");
			if (mfp != NULL) {
				known->printMotifs(mfp);
				fclose(mfp);
			}
		}
		if (known->maskFlag) {
			known->seqArray->printSequenceFile(outputFp);
		}
	}
	if (outputFileName != NULL) fclose(outputFp);

	if (known->fastaFlag) {
		fprintf(stderr, "\tCleaning up temporary sequence and group files:\n");
		fprintf(stderr, "\t\t%s\n",known->seqFile);
		fprintf(stderr, "\t\t%s\n",known->groupFile);
		(void)unlink(known->seqFile);
		(void)unlink(known->groupFile);
	}
}


void splitMotif2(char* string, char** cols, int &numCols, char delim) {
	cols[0] = string;
	numCols=1;
	int len = strlen(string);
	for (int i=0;i<len;i++) {
		if (string[i] == delim) {
			string[i] = '\0';
			cols[numCols] = &(string[i+1]);
			numCols++;
		} else if (string[i] == '\n') {
			string[i] = '\0';
		} else if (string[i] == '\r') {
			string[i] = '\0';
		}
	}
}


DeNovoMotifs::DeNovoMotifs() {
	init();
}
DeNovoMotifs::~DeNovoMotifs() {
	if (oligoArray != NULL) {
		delete oligoArray;
	}
	oligoArray = NULL;
	if (seedOligos != NULL) delete []seedOligos;
	if (scorer != NULL) delete scorer;
	if (fullSequences != NULL) delete fullSequences;
	// fill this in later
}
void DeNovoMotifs::freeOligos() {
	if (oligoArray != NULL) delete oligoArray;
	oligoArray = NULL;
}
void DeNovoMotifs::init() {
	seqFile = NULL;
	groupFile = NULL;
	inputFASTAfile = NULL;
	bgFASTAfile = NULL;
	numTargets = 0.0;
	numBackground = 0.0;
	maxMisMatches = 2;
	oligoArray = NULL;
	alphaSize = DNA_ALPHA;
	totalOligos = 0;
	strand = MOTIF_STRAND_BOTH;
	oligoLength = 10;
	normLength = 0;
	noDuplicateOligos=1;

	maskMethod = DENOVO_MASK_FULL;
	fullSequences = NULL;
	zoopsFlag = SCORER_APPROXIMATE_ZOOPS;
	seedOligos=NULL;
	normFile = NULL;
	tmpFile = NULL;
	sitefp = NULL;
	numSeedOligos = 0;
	bufferLength = -1;
	numOligoEstimate = 0;
	numCPUs=1;
	fastaFlag = 0;
	randFlag = 0;
	maximumExpectedPerBp = -1.0;
	maxCacheSize = MAX_STAT_CACHE_SIZE*1000000/sizeof(double);
	minimumSeedLogp = MOTIF_SEED_MINIMUM_LOGP;

	finalZoopsFlag = 1;
	scorer = NULL;
	siteReduceThreshold = MOTIF_DEFAULT_REDUCE_PERCENTAGE;

	optimizeGivenFlag = 0;
	maxBackgroundPercent = MOTIF_MAX_BACKGROUND_PERCENT;

	scoringMethod = MOTIF_SCORING_BINOMIAL;
	misMatchTargetOnlyFlag = 1;
	misMatchMinPvalue = MOTIF_MISMATCH_PVALUE_MINIMUM;
	minimumMisMatchToStartSkipping = MOTIF_MISMATCH_MINIMUM_SKIP;
	minimumCountsToStartSkipping = MOTIF_MISMATCH_MINIMUM_COUNTS;
	minimumFractionToStartSkipping = 0.000; // not used

	totalMotifsToFind = MOTIF_DEFAULT_NUMBER_OF_MOTIFS;
	numTrialMotifs = MOTIFS_NUM_POTENTIAL_MOTIFS;
	maxOptimizationIterations = MOTIFS_MAX_OPT_ITERATIONS;
	motifOptimizationMethod = MOTIF_OPTIMZATION_DIFFERENT_THRESHOLDS;
	
	maxAutoNormIters = MAX_OLIGO_NORMALIZATION_ITERATIONS;
	
}


void DeNovoMotifs::checkParameters() {
	if (numCPUs < 1) numCPUs = 1;
	if (numCPUs > 1) {
		int og = numTrialMotifs;
		int trialsPerCPU = ((int)(((((double)numTrialMotifs)/((double)numCPUs)))+0.5));
		if (trialsPerCPU < 1) trialsPerCPU = 1;
		numTrialMotifs = trialsPerCPU*numCPUs;
		if (og != numTrialMotifs) {
			fprintf(stderr, "\tNumber of Trial motifs (-T) set to %d (from %d) to work well with %d CPUs\n",
							numTrialMotifs, og, numCPUs);
		}
	}
	if (inputFASTAfile == NULL || bgFASTAfile == NULL) {
		if (seqFile == NULL ){
			fprintf(stderr, "!!! Need some type of sequence file for homer to do anything !!!\n");
			exit(0);
		}
		if (groupFile == NULL ){
			fprintf(stderr, "!!! Need a group file for homer to do anything !!!\n");
			exit(0);
		}
	} else {
		if (inputFASTAfile == NULL) {
			fprintf(stderr, "!!! Missing input FASTA file (-i) !!!\n");
			exit(0);
		}
		if (bgFASTAfile == NULL) {
			fprintf(stderr, "!!! Missing background FASTA file (-b) !!!\n");
			exit(0);
		}
		fastaFlag = 1;
		fprintf(stderr, "\tTreating input files as FASTA format\n");
	}
	if (bufferLength == -1) {
		if (oligoLength < 8) {
			bufferLength = 1;
		} else {
			bufferLength = 1 + (oligoLength-8)/4;
		}
		if (bufferLength > 1) {
			fprintf(stderr, "\t-blen automatically set to %d\n", bufferLength);
		}
		if (bufferLength > 1000) bufferLength = 0;
	}
	if (maximumExpectedPerBp < 0.0) {
		maximumExpectedPerBp = MOTIF_MAX_HIT_RATE_TO_CONSIDER;
		if (oligoLength == 10) {
			maximumExpectedPerBp = MOTIF_MAX_HIT_RATE_TO_CONSIDER*2;
		} else if (oligoLength == 9) {
			maximumExpectedPerBp = MOTIF_MAX_HIT_RATE_TO_CONSIDER*4;
		} else if (oligoLength < 9) {
			maximumExpectedPerBp = MOTIF_MAX_HIT_RATE_TO_CONSIDER*8;
		}
	}
}
void DeNovoMotifs::parseFasta() {
	SequenceArray* seqarray = new SequenceArray();
	seqFile = new char[10000];
	groupFile = new char[10000];
	srand(time(NULL));
	int r = rand();
	sprintf(seqFile, "0%d.seq", r);
	sprintf(groupFile, "0%d.group", r);
	seqarray->parseFasta2SeqAndGroupFiles(inputFASTAfile,bgFASTAfile,seqFile,groupFile);
	delete seqarray;
}


void DeNovoMotifs::normalizeOligos() {
	oligoArray->maxAutoNormIters = maxAutoNormIters;
	if (normLength > 0) {
		oligoArray->normalizeOligos(normLength,normFile);
	}
}

void DeNovoMotifs::readSequence2Tree() {

	Hashtable* groupdata = new Hashtable();

	fprintf(stderr, "\tScanning input files...\n");
	
	FILE* fp = fopen(groupFile, "r");
	if (fp == NULL) {
		fprintf(stderr, "!!! Could not open group file: %s !!!\n", groupFile);
		exit(0);
	}
	char* buf = new char[MOTIF2_BUFFER];
	char** cols = new char*[1000];
	int numCols;
	char* seq = NULL;
	char* name = NULL;

	int totalTargetsGroup = 0;
	while (fgets(buf, MOTIF2_BUFFER, fp) != NULL) {
		splitMotif2(buf,cols,numCols,'\t');
		if (numCols < 2) continue;
		GroupInfo* gi = new GroupInfo();
		gi->weight = 1.0;
		gi->group = 0;
		name = cols[0];
		if (numCols > 1) {
			sscanf(cols[1],"%d",&(gi->group));
			if (gi->group) totalTargetsGroup++;
		}
		if (numCols > 2) {
			sscanf(cols[2],"%lf",&(gi->weight));
		}
		groupdata->insert(gi, name);
	}
	fclose(fp);

	if (randFlag) {
		//randomize group labels
		double targetFrac = ((double)totalTargetsGroup)/((double)groupdata->total);
		srand(time(NULL));
		char** keys = groupdata->keys();
		for (int i=0;i<groupdata->total;i++) {
			GroupInfo* gi = (GroupInfo*)groupdata->search(keys[i]);
			if (((double)rand())/((double)RAND_MAX) < targetFrac) {
				gi->group = 1;
			} else {
				gi->group = 0;
			}
			delete [](keys[i]);
		}
		delete []keys;
	}


	fp = fopen(seqFile, "r");
	if (fp == NULL) {
		fprintf(stderr, "!!! Could not open seq file: %s !!!\n", seqFile);
		exit(0);
	}

	numOligoEstimate = 0;
	int totalSequences = 0;
	long long int totalLength = 0;
	while (fgets(buf, MOTIF2_BUFFER, fp) != NULL) {
		splitMotif2(buf,cols,numCols,'\t');
		if (numCols < 2) continue;
		name = cols[0];
		seq = cols[1];
//fprintf(stderr, " %s",name);
		GroupInfo* gi = (GroupInfo*)groupdata->search(name);
		if (gi == NULL) {
			continue;
		}
		totalSequences++;
		int slen = strlen(seq);
		int contribution = slen-oligoLength+1;
		if (contribution < 0) contribution = 0;
		numOligoEstimate += contribution;
		totalLength += slen;
	}
	if (numOligoEstimate <= 0) {
		fprintf(stderr, "!!! Something is wrong... are you sure you chose the right length for motif finding?\n");
		fprintf(stderr, "!!! i.e. also check your sequence file!!!\n");
		exit(0);
	}
	double avgSeqLength = totalLength/totalSequences;

	oligoArray = new OligoArray();
	oligoArray->init(oligoLength,numOligoEstimate, alphaSize, strand, numCPUs);
	double totalTargetOligos = 0.0;
	double totalBackgroundOligos = 0.0;

	rewind(fp);	

	//now read in for real
	int curNumSequences = 0;
	int barLength = totalSequences/80;
	int curThreshold = barLength;
	fprintf(stderr, "\tParsing sequences...\n");
	fprintf(stderr, "\t|0%%                                   50%%                                  100%%|\n");
	fprintf(stderr, "\t");
	while (fgets(buf, MOTIF2_BUFFER, fp) != NULL) {
		splitMotif2(buf,cols,numCols,'\t');
		if (numCols < 2) continue;
		GroupInfo* gi = (GroupInfo*)groupdata->search(cols[0]);
		if (gi == NULL) {
			continue;
		}
		name = cols[0];
		seq = cols[1];
		cleanSequence(seq);
		if (strand == MOTIF_STRAND_NEG) {
			Oligo::revopp(seq);
		}
			
		curNumSequences++;
		if (curNumSequences > curThreshold) {
			curThreshold += barLength;
			fprintf(stderr, "=");
		}
		double w = gi->weight;
		unsigned int g = gi->group;
		double pvalue = 0.0;
		double nvalue = 0.0;
		if (g==0) nvalue = w;
		if (g==1) pvalue = w;


		int length = strlen(seq);

		double numberOfGoodOligos = 0.0;
		for (int i=0;i<=length-(int)oligoLength;i++) {
			char tmp = seq[i+oligoLength];
			seq[i+oligoLength]='\0';
			char* oligo = &(seq[i]);
			int ok = checkSequence(oligo);
			if (ok) numberOfGoodOligos += 1.0;
			seq[i+oligoLength]=tmp;
		}
		if (numberOfGoodOligos < 1.0) {
			continue;
		}

		numTargets += pvalue;
		numBackground += nvalue;

		if (numberOfGoodOligos > avgSeqLength) {
			pvalue *= avgSeqLength/numberOfGoodOligos;
			nvalue *= avgSeqLength/numberOfGoodOligos;
		}
		for (int i=0;i<=length-(int)oligoLength;i++) {
			char tmp = seq[i+oligoLength];
			seq[i+oligoLength]='\0';
			char* oligo = &(seq[i]);
			int ok = checkSequence(oligo);
			if (ok) {
				totalOligos += oligoArray->addOligo(oligo,pvalue,nvalue,-1.0);
				//totalOligos += oligoArray->addOligo(oligo,pvalue,nvalue,(double)curNumSequences);
				totalTargetOligos += pvalue;
				totalBackgroundOligos += nvalue;
			}
			seq[i+oligoLength]=tmp;
		}
		
	}
	fprintf(stderr, "\n");
	fclose(fp);
	fprintf(stderr, "\tTotal number of Oligos: %lld\n", totalOligos);

	if (numTargets < 0.1 || numBackground < 0.1) {
		fprintf(stderr, "!!! Something is wrong - either no target sequences or no background sequences!!!\n");
		exit(0);
	}

	double oligoPerSeqTarget = totalTargetOligos/numTargets;
	//double oligoPerSeqBackround = totalBackgroundOligos/numBackground;
	double newNumBackground = totalBackgroundOligos/oligoPerSeqTarget;
	fprintf(stderr, "\tAutoadjustment for sequence coverage in background: %.2lfx\n",newNumBackground/numBackground);
	numBackground = newNumBackground;

	char** keys = groupdata->keys();
	for (int i=0;i<groupdata->total;i++) {
		GroupInfo* gi = (GroupInfo*)groupdata->search(keys[i]);
		delete gi;
		delete [](keys[i]);
	}
	delete []keys;
	delete groupdata;
	delete []buf;
	delete []cols;

	oligoArray->optimizeMemoryUsage();

}

void DeNovoMotifs::loadFullSequences() {

	if (fullSequences != NULL) delete fullSequences;
	fullSequences = new SequenceArray();
	fullSequences->readSequenceFiles(seqFile,groupFile);
	fullSequences->initializeScorer(scoringMethod,SCORER_STANDARD_ZOOPS,maxCacheSize);
	//if (normLength > 0 && !findFlag) {
	//	seqArray->normalizeSequence(normLength, strand, normFile);
	//}
}

int checkSequence(char* seq) {
	int i=0;
	while (seq[i] != '\0') {
		switch (seq[i]) {
			case 'a':
				seq[i] = 'A';
				break;
			case 'c':
				seq[i] = 'C';
				break;
			case 'g':
				seq[i] = 'G';
				break;
			case 't':
				seq[i] = 'T';
				break;
			case 'A':
			case 'C':
			case 'G':
			case 'T':
				break;
			default:
				return 0;
		}
		i++;
	}
	return 1;
}

ThreadArgs_findMisMatchEnrichmentThread::ThreadArgs_findMisMatchEnrichmentThread(
							DeNovoMotifs* obj,int cpuNum,int mis,int lasttotal,int vlevel) {
	denovo = obj;
	cpu = cpuNum;
	maxMisMatches = mis;
	totalChecked=0;
	numSkipped=0;
	dueToTargetTags=0;
    dueToLessThanLast=0;
   	dueToLowPvalue=0;
    lastChecked=0;
    lastCheckedTotal= lasttotal;
	verboseLevel = vlevel;
}


void DeNovoMotifs::findMisMatchEnrichment(int maxMisMatches,FILE* fp,int verboseLevel) {
	//oligoArray->resetOligoArray(); not needed with "activeOligos"
	int inc = (int)(((double)oligoArray->numActiveOligos)/80.0);
	int lastCheckedTotal = (int)oligoArray->numActiveOligos;
	if (verboseLevel != HOMER_VERBOSE_LEVEL0) {
		fprintf(stderr, "\n\tGlobal Optimization Phase: Looking for enriched oligos with up to %d mismatches...\n",maxMisMatches);
	}
	for (int j=0;j<=maxMisMatches;j++) {
		if (verboseLevel != HOMER_VERBOSE_LEVEL0) {
			fprintf(stderr, "\n\tScreening oligos %d (allowing %d mismatches):\n", oligoArray->numOligos, j);
			fprintf(stderr, "\t|0%%                                   50%%                                  100%%|\n");
			fprintf(stderr, "\t");
		}
		int totalChecked=0;
		int numSkipped=0;
		int dueToTargetTags=0;
		int dueToLessThanLast=0;
		int dueToLowPvalue=0;
		int lastChecked = 0;

		
		pthread_mutex_init(&mutex, NULL);
		mutexIndex = 0;
		mutexTotal = oligoArray->numActiveOligos;
		nextMeter = inc;
		incMeter = inc;
		ThreadArgs_findMisMatchEnrichmentThread** args = new ThreadArgs_findMisMatchEnrichmentThread*[numCPUs];

		pthread_t* thread = new pthread_t[numCPUs];
		pthread_attr_t attr;
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

		void* status;
		for (unsigned int i=0;i<numCPUs;i++) {
			args[i] = new ThreadArgs_findMisMatchEnrichmentThread(this,i,j,lastCheckedTotal,verboseLevel);
			(void)pthread_create(&(thread[i]), &attr, 
							&DeNovoMotifs_findMisMatchEnrichmentThread, args[i]);
		}
		pthread_attr_destroy(&attr);
		for (unsigned int i=0;i<numCPUs;i++) {
			(void)pthread_join(thread[i], &status);
			totalChecked += args[i]->totalChecked;
			numSkipped += args[i]->numSkipped;
			dueToTargetTags += args[i]->dueToTargetTags;
			dueToLessThanLast += args[i]->dueToLessThanLast;
			dueToLowPvalue += args[i]->dueToLowPvalue;
			if (args[i]->lastChecked > lastChecked) lastChecked = args[i]->lastChecked;
			delete args[i];
		}
		delete []thread;
		delete []args;
		pthread_mutex_destroy(&mutex);

		if (verboseLevel != HOMER_VERBOSE_LEVEL0) {
			fprintf(stderr, "\n\t\t%.2lf%% skipped, %.2lf%% checked (%d of %d), of those checked:\n", 
					(100.0*numSkipped)/oligoArray->numActiveOligos,
					(100.0*totalChecked)/oligoArray->numActiveOligos,
					totalChecked, oligoArray->numActiveOligos);
			fprintf(stderr, "\t\t%.2lf%% not in target, %.2lf%% increased p-value, %.2lf%% high p-value\n", 
					(100.0*dueToTargetTags)/oligoArray->numActiveOligos, 
					(100.0*dueToLessThanLast)/oligoArray->numActiveOligos,
					(dueToLowPvalue*100.0)/oligoArray->numActiveOligos);
		}
		inc = (int)(((double)lastChecked)/80.0);
		lastCheckedTotal = lastChecked;
	}
	oligoArray->sortActiveOligos();
	if (fp != NULL) {
		fprintf(fp, "OligoSequence\tTarget Occurences\tBackground Occurences\tEnrichment(ln p, up to %d mismatches)\n",
					maxMisMatches);
		for (int i=0;i<oligoArray->numActiveOligos;i++) {
			fprintf(fp, "%s\t%.1lf\t%.1lf\t%.3lf\n",oligoArray->activeOligos[i]->seq,oligoArray->activeOligos[i]->numTarget,
					oligoArray->activeOligos[i]->numBackground,oligoArray->activeOligos[i]->value);
		}
	}
}

void* DeNovoMotifs_findMisMatchEnrichmentThread(void* threadArgs) {

	ThreadArgs_findMisMatchEnrichmentThread* args = (ThreadArgs_findMisMatchEnrichmentThread*)threadArgs;
	args->denovo->findMisMatchEnrichmentThread(args);
	pthread_exit(threadArgs);
	return NULL;
}

void DeNovoMotifs::findMisMatchEnrichmentThread(ThreadArgs_findMisMatchEnrichmentThread* args) {

	int i=-1;

	while (i < mutexTotal) {

		pthread_mutex_lock(&mutex);
		i=mutexIndex++;
		if (i > nextMeter) {
			///////////////////Need to deal with "next"
			if (i <= args->lastCheckedTotal) {
				if (args->verboseLevel != HOMER_VERBOSE_LEVEL0) {
					fprintf(stderr, "=");
				}
			}
			//fprintf(stderr, "=%d %d",i,next);
			nextMeter += incMeter;
		}
		pthread_mutex_unlock(&mutex);
		if (i >= mutexTotal) break;

		double logp = 0.0;
		double vt = (double)oligoArray->activeOligos[i]->numTarget;
		double vbg = (double)oligoArray->activeOligos[i]->numBackground;
		int skipFlag = oligoArray->activeOligos[i]->flag2;
		if (!skipFlag && misMatchTargetOnlyFlag) {
			if (vt < 0.01) {
				skipFlag = 1;
				oligoArray->activeOligos[i]->flag2 = 1;
				args->dueToTargetTags++;
			}
		}
		if (!skipFlag) {
			if (args->maxMisMatches == 0) {
			} else {
				vt = 0.0;
				vbg = 0.0;
				oligoArray->searchOligo(oligoArray->activeOligos[i]->seq,
									args->maxMisMatches,vt,vbg,OLIGO_SEARCH_KEEP,args->cpu);
			}
	
			logp = scorer->scoreEnrichment(vt,vbg);
			if (oligoArray->activeOligos[i]->value > logp) {
				oligoArray->activeOligos[i]->value = logp;
			} else {
				if (minimumMisMatchToStartSkipping <= args->maxMisMatches 
							&& (vt+vbg > minimumCountsToStartSkipping)) {
					oligoArray->activeOligos[i]->flag2 = 1;
					args->dueToLessThanLast++;
				}
			}
			if ((minimumMisMatchToStartSkipping <= args->maxMisMatches)
						&& (oligoArray->activeOligos[i]->value > misMatchMinPvalue)
						&& (vt+vbg > minimumCountsToStartSkipping)) {
				oligoArray->activeOligos[i]->flag2 = 1;
				args->dueToLowPvalue++;
			}
			args->totalChecked++;
			args->lastChecked = i+1;
		} else {
			args->numSkipped++;
		}
	}
}


void DeNovoMotifs::getSeedOligos() {
	//assmes that oligoArray->oligos is sorted in terms of enrichment
	int maxSeedOligos = (int)MOTIF_MAX_SEED_OLIGOS_TOTAL;
	if (maxSeedOligos > oligoArray->numOligos) {
		maxSeedOligos = oligoArray->numOligos;
	}
	seedOligos = new Oligo*[maxSeedOligos];
	for (int i=0;i<maxSeedOligos;i++) {
		seedOligos[i] = oligoArray->activeOligos[i];
	}
	numSeedOligos = maxSeedOligos;
}

void DeNovoMotifs::initializeEnrichmentScoring() {
	scorer = new EnrichmentScorer();
	scorer->initializeEnrichmentScoring(scoringMethod,numTargets, numBackground,
															zoopsFlag,maxCacheSize);
	oligoArray->initializeOligoPvalues(scorer);
}

void DeNovoMotifs::optimizeMotifMatrix(Motif* &motif) {

	int cpu = 0;
	optimizeMotifThreshold(motif,cpu);
	if (motif->badFlag) return;

	if (motifOptimizationMethod == MOTIF_OPTIMZATION_DIFFERENT_THRESHOLDS) {

		Motif** newMotifs = new Motif*[numTrialMotifs];
		int newMotifIndex = 0;
		Motif* matrix = new Motif();
		double ogNumTarget = motif->numTarget;
		double ogNumBackground = motif->numBackground;
		double ogLogp = motif->logp;
		matrix->setBlank(oligoLength);

		Oligo** utilityOligos = oligoArray->utilityOligos[cpu];
		int numUtilityOligos = oligoArray->numUtilityOligos[cpu];

		int halfTrialMotifs = (int)((numTrialMotifs+1)/2.0);	
		if (halfTrialMotifs < 1) halfTrialMotifs=1;

		double inc = (utilityOligos[0]->value - motif->threshold)/(halfTrialMotifs);
		if (inc < (utilityOligos[0]->value)/((double)2*halfTrialMotifs)) {
			inc = (utilityOligos[0]->value)/((double)2*halfTrialMotifs);
		}

		double nextThreshold = utilityOligos[0]->value - inc;

		for (int i=0;i<numUtilityOligos;i++) {
			double currentThreshold = utilityOligos[i]->value;

			
			double nlogp = scorer->scoreEnrichment(ogNumTarget+utilityOligos[i]->numTarget,
													ogNumBackground+utilityOligos[i]->numBackground);
			utilityOligos[i]->value = 0;
			double diffPvalue = nlogp - ogLogp;
			//diffPvalue = -utilityOligos[i]->numTarget;
			if (diffPvalue > 0.0) {
				utilityOligos[i]->value = 0;
				//utilityOligos[i]->value = utilityOligos[i]->logp+diffPvalue;
				//if (utilityOligos[i]->value > 0) utilityOligos[i]->value = 0.0;
			} else {
				//utilityOligos[i]->value = diffPvalue;
				utilityOligos[i]->value = utilityOligos[i]->logp;
			}
				

			
			//set oligo contribution value
			matrix->addOligo(utilityOligos[i]);
			//fprintf(stderr, "======== %s %f\n", utilityOligos[i]->seq, utilityOligos[i]->logp);

			if (i+1<numUtilityOligos) {
				// need to keep going if a bunch of oligos have the same threshold
				double diff = currentThreshold - utilityOligos[i+1]->value;
				if (diff < FLOAT_EPS) {
					continue;
				}
				//fprintf(stderr, "diff=%le\n", diff);
			}

			if (currentThreshold < nextThreshold) {
				//matrix->print(stderr);
				newMotifs[newMotifIndex++] = new Motif(matrix);
				nextThreshold -= inc;
				if (newMotifIndex >= numTrialMotifs) break;
			}
		}


		//easily part to make parallel
		pthread_mutex_init(&mutex, NULL);
		mutexIndex = 0;
		mutexTotal = newMotifIndex;

		ThreadArgs_optimizeMotifThresholdThread** args = new ThreadArgs_optimizeMotifThresholdThread*[numCPUs];

		pthread_t* thread = new pthread_t[numCPUs];
		pthread_attr_t attr;
		pthread_attr_init(&attr);
		pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

		void* status;
		for (unsigned int i=0;i<numCPUs;i++) {
			args[i] = new ThreadArgs_optimizeMotifThresholdThread(this,i,newMotifs);
			(void)pthread_create(&(thread[i]), &attr, 
							&DeNovoMotifs_optimizeMotifThresholdThread, args[i]);
		}
		pthread_attr_destroy(&attr);
		for (unsigned int i=0;i<numCPUs;i++) {
			(void)pthread_join(thread[i], &status);
			delete args[i];
		}
		delete []thread;
		delete []args;
		pthread_mutex_destroy(&mutex);


		int bestIndex = -1;
		double bestIndexScore  = motif->logp; 
		for (int i=0;i<newMotifIndex;i++) {
			//newMotifs[i]->print(stderr);
			if (newMotifs[i]->logp < bestIndexScore) {
				bestIndex = i;
				bestIndexScore = newMotifs[i]->logp;
			}
		}
		//fprintf(stderr, "Best Index = %d of %d (new: %lf\told: %lf)\n",bestIndex, newMotifIndex,bestIndexScore,motif->logp);
		if (bestIndex > -1) {
			Motif* tmp = motif;
			motif = newMotifs[bestIndex];
			newMotifs[bestIndex] = tmp;
		}
		for (int i=0;i<newMotifIndex;i++) {
			delete newMotifs[i];
		}
		delete matrix;	
		delete []newMotifs;
	}	
}

void* DeNovoMotifs_optimizeMotifThresholdThread(void* threadArgs) {
	ThreadArgs_optimizeMotifThresholdThread* args = (ThreadArgs_optimizeMotifThresholdThread*)threadArgs;
	args->denovo->optimizeMotifThresholdThread(args);
	pthread_exit(threadArgs);
	return NULL;
}
ThreadArgs_optimizeMotifThresholdThread::ThreadArgs_optimizeMotifThresholdThread(
								DeNovoMotifs* obj,int newcpu,Motif** newmotifs) {
	denovo = obj;
	cpu = newcpu;
	motifs=newmotifs;
}
void DeNovoMotifs::optimizeMotifThresholdThread(ThreadArgs_optimizeMotifThresholdThread* args) {
	int i=-1;
	while (i<mutexTotal) {
		pthread_mutex_lock(&mutex);
		i=mutexIndex++;
		pthread_mutex_unlock(&mutex);
		if (i >= mutexTotal) break;

		args->motifs[i]->normalize();
		optimizeMotifThreshold(args->motifs[i],args->cpu);
	}
}

void DeNovoMotifs::optimizeMotifThreshold(Motif* motif, int cpu) {

	motif->setMinimumScore(maximumExpectedPerBp);
	oligoArray->scoreOligosWithMotif(motif,motif->minimumScore,MOTIF_SCORE_OLIGOS,cpu);
	oligoArray->sortUtilityOligos(cpu);

	double bestThreshold = 0;
	double bestScore = 0;
	double bestTvalue = 0;
	double bestBGvalue = 0;
	double currentTvalue = 0;
	double currentBGvalue = 0;
	double currentThreshold = 0;
	double maxBackground = maxBackgroundPercent*numBackground;
	Oligo** utilityOligos = oligoArray->utilityOligos[cpu];
	int numUtilityOligos = oligoArray->numUtilityOligos[cpu];
	if (numUtilityOligos < 1) {
		// basically, this is no longer a valid motif... need to bail gracefully...
		motif->badFlag = 1;
		return;
	}
	//int bestIndex = 0;
	//int lastChecked = 0;
	double bestPossibleThreshold = utilityOligos[0]->cpuData[cpu];
	double lowestPossibleThreshold = utilityOligos[numUtilityOligos-1]->cpuData[cpu];
	//fprintf(stderr, "\t%lf\t%lf\t%lf\n", motif->minimumScore, bestPossibleThreshold, lowestPossibleThreshold);
	for (int i=0;i<numUtilityOligos;i++) {
		currentTvalue += utilityOligos[i]->numTarget;
		currentBGvalue += utilityOligos[i]->numBackground;
		currentThreshold = utilityOligos[i]->cpuData[cpu];
		if (i+1<numUtilityOligos) {
			// need to keep going if a bunch of oligos have the same threshold
			double diff = currentThreshold - utilityOligos[i+1]->cpuData[cpu];
			if (diff < FLOAT_EPS) {
				continue;
			}
			//fprintf(stderr, "diff=%le\n", diff);
		}
		double logp = scorer->scoreEnrichment(currentTvalue,currentBGvalue);
		//fprintf(stderr, "%lf\t%lf\t%lf\t%lf\n", currentThreshold, logp,currentTvalue,numTargets);
		if (logp < bestScore) {
			//bestIndex = i;
			bestScore = logp;
			bestThreshold = currentThreshold;
			bestTvalue = currentTvalue;
			bestBGvalue = currentBGvalue;
		}
		if (currentThreshold < motif->minimumScore) break;
		if (logp > bestScore/3.0 && currentTvalue > numTargets) break;
		if (currentBGvalue > maxBackground) break;
		//lastChecked = i;
		/*fprintf(stderr, "%s\t%lf\t%lf\t%lf\t%lf\t%lf\n", utilityOligos[i]->seq, 
					currentTvalue, currentBGvalue, currentThreshold, logp, bestScore);*/
	}
	//fprintf(stderr, "%d of %d of %d\n", bestIndex, lastChecked,numUtilityOligos);
	double degfrac = (bestPossibleThreshold-bestThreshold)
							/(bestPossibleThreshold-lowestPossibleThreshold+0.01);
	if (degfrac > 0.9) {
		fprintf(stderr, "Optimal threshold is low - motif appears very degenerate!!!\n");
	}
	
	motif->threshold = bestThreshold-0.000001;
	motif->logp = bestScore;
	motif->numTarget = bestTvalue;
	motif->numBackground = bestBGvalue;
	motif->freqTarget = bestTvalue/numTargets;
	motif->freqBackground = bestBGvalue/numBackground;
	if (zoopsFlag == SCORER_APPROXIMATE_ZOOPS) {
		motif->freqTarget = scorer->getZoopsApproxTarget(bestTvalue)/numTargets;
		motif->freqBackground = scorer->getZoopsApproxBackground(bestBGvalue)/numBackground;
	}
//	fprintf(stderr, "\t%.3lf\t%.3lf\n", motif->minimumScore, motif->threshold);
}

void DeNovoMotifs::optimizeGivenMotifs() {
	int defaultCPU = 0;

	FILE* tmpfp = NULL;
	if (tmpFile == NULL) {
		tmpfp = fopen(".tmp.motifs","w");
		if (tmpfp == NULL) fprintf(stderr, "\t!! Could not open tmp results file (.tmp.motifs) !!\n");
	} else {
		tmpfp = fopen(tmpFile,"w");
		if (tmpfp == NULL) fprintf(stderr, "\t!! Could not open tmp results file (%s) !!\n", tmpFile);
	}

	fprintf(stderr, "\n\tLocal Optimization Phase:\n");
	for (int i=0;i<numMotifs;i++) {

		motifs[i]->getConsensus();
		fprintf(stderr, "\t\tInitial Consensus: %s\n", motifs[i]->consensus);
	
		double curLogp = 0.0;
		for (int j=0;j<maxOptimizationIterations;j++) {
			optimizeMotifMatrix(motifs[i]);
			motifs[i]->getConsensus();
			fprintf(stderr, "\t\tRound %d: %.2lf", j+1, motifs[i]->logp);
			fprintf(stderr, " %s ", motifs[i]->consensus);
			fprintf(stderr, "T:%.1lf(%.2f%%),B:%.1lf(%.2f%%),P:1e%d\n",motifs[i]->numTarget, 
									motifs[i]->freqTarget*100, motifs[i]->numBackground,
									motifs[i]->freqBackground*100,(int)(motifs[i]->logp/2.3026));
			if (curLogp-FLOAT_EPS < motifs[i]->logp) {
				break;
			}
			curLogp = motifs[i]->logp;
		}

		motifs[i]->getConsensus();
		fprintf(stderr, "\t\t=Final=: %.2lf", motifs[i]->logp);
		fprintf(stderr, " %s ", motifs[i]->consensus);
		fprintf(stderr, "T:%.1lf(%.2f%%),B:%.1lf(%.2f%%),P:1e%d\n",motifs[i]->numTarget, 
									motifs[i]->freqTarget*100, motifs[i]->numBackground,
									motifs[i]->freqBackground*100,(int)(motifs[i]->logp/2.3026));

		//do not remove oligos from the tree for this type of motif finding
		oligoArray->scoreOligosWithMotif(motifs[i],motifs[i]->minimumScore,MOTIF_SCORE_OLIGOS,defaultCPU);

		int totalBoundOligos=0;
		for (int j=0;j<oligoArray->numUtilityOligos[0];j++) {
			if (oligoArray->utilityOligos[0][j]->value > motifs[i]->threshold) {
				totalBoundOligos++;
			}
		}
		motifs[i]->getConsensus();
		fprintf(stderr, "\t\t%d of %d - %s Bound Oligos: %d\n",i+1,totalMotifsToFind,
											motifs[i]->consensus,totalBoundOligos);
		motifs[i]->print(tmpfp);
		if (tmpfp != NULL) fflush(tmpfp);
	}
	if (tmpfp != NULL) fclose(tmpfp);
}

void DeNovoMotifs::optimizeSeeds2Motifs() {

	//assmes that oligoArray->activeOligos is sorted in terms of enrichment
	if (oligoArray->numActiveOligos > oligoArray->numOligos*MOTIF_MAX_SEED_OLIGOS_FRACTION) {
		oligoArray->numActiveOligos = (int)(oligoArray->numOligos*MOTIF_MAX_SEED_OLIGOS_FRACTION);
	}
	if (oligoArray->numActiveOligos > MOTIF_MAX_SEED_OLIGOS_TOTAL) {
		oligoArray->numActiveOligos = (int)(MOTIF_MAX_SEED_OLIGOS_TOTAL);
	}
	if (minimumSeedLogp > MOTIF_SEED_MINIMUM_LOGP-0.01 
			&& minimumSeedLogp < MOTIF_SEED_MINIMUM_LOGP+0.01) {
		//default value - adjust for the target genes, the less of which lowers pvalues
		minimumSeedLogp = MOTIF_SEED_MINIMUM_LOGP+(4.0-log((double)numTargets)/log(10.0));
	}

	motifs = new Motif*[totalMotifsToFind];
	numMotifs = 0;
	int seedOligoIndex = 0;
	char* bufferSeq = new char[oligoLength+1];

	OligoArray* boundOligos = new OligoArray();
	boundOligos->init(oligoLength,numOligoEstimate/4,alphaSize,strand,1);
	int defaultCPU = 0;

	FILE* tmpfp = NULL;
	if (tmpFile == NULL) {
		tmpfp = fopen(".tmp.motifs","w");
		if (tmpfp == NULL) fprintf(stderr, "\t!! Could not open tmp results file (.tmp.motifs) !!\n");
	} else {
		tmpfp = fopen(tmpFile,"w");
		if (tmpfp == NULL) fprintf(stderr, "\t!! Could not open tmp results file (%s) !!\n", tmpFile);
	}


	fprintf(stderr, "\n\tLocal Optimization Phase:\n");
	for (int i=0;i<totalMotifsToFind;i++) {

		if (seedOligoIndex >= oligoArray->numActiveOligos) break;
		char* seedSeq = oligoArray->activeOligos[seedOligoIndex]->seq;
		double seedlogp = oligoArray->activeOligos[seedOligoIndex]->value;

		if (i > 4 && seedlogp > minimumSeedLogp && maskMethod == DENOVO_MASK_FULL) {
			fprintf(stderr, "\tRemaining seeds don't look promising (After initial 5 motifs, logp %.3lf > %.3lf)\n",
											seedlogp,minimumSeedLogp); 
			//fprintf(stderr, "\t\tTo analyze further, use -minlp <#> option\n"); 
			//fprintf(stderr, "\t\t%s\n", oligoArray->activeOligos[seedOligoIndex]->seq); 
			break;
		}
		oligoArray->activeOligos[seedOligoIndex++] = NULL;
		//fprintf(stderr, "\t%s",seedSeq);

		//first check if the oligo still exists and hasn't been masked out yet
		int there=0;
		for (int j=0;j<=bufferLength;j++) {
			double vt=0.0,vbg=0.0;
			for (int k=0;k<j;k++) bufferSeq[k]='N';
			for (int k=0;k<(int)oligoLength-j;k++) bufferSeq[k+j]=seedSeq[k];
			bufferSeq[oligoLength]='\0';
			there = boundOligos->searchOligo(bufferSeq,0,vt,vbg,OLIGO_SEARCH_KEEP,defaultCPU);
			if (there > 0) break;
			if (j==0) continue;
			vt=0.0,vbg=0.0;
			for (int k=0;k<j;k++) bufferSeq[oligoLength-k-1]='N';
			for (int k=0;k<(int)oligoLength-j;k++) bufferSeq[k]=seedSeq[k+j];
			bufferSeq[oligoLength]='\0';
			there = boundOligos->searchOligo(bufferSeq,0,vt,vbg,OLIGO_SEARCH_KEEP,defaultCPU);
			if (there > 0) break;
		}
		if (there > 0) { 
			i--;
			//fprintf(stderr, "\t%s... too similar, skipping\n",seedSeq);
			continue;
		}
		fprintf(stderr, "\t%d of %d Initial Sequence: %s... (%.3lf)\n",i+1,totalMotifsToFind,
															seedSeq,seedlogp);
		//fprintf(stderr, "...\n");


		motifs[i] = new Motif(seedSeq);
		motifs[i]->getConsensus();
		//fprintf(stderr, "\t\tInitial Consensus: %s\n", motifs[i]->consensus);
		numMotifs=i+1;
	
		double curLogp = 0.0;
		for (int j=0;j<maxOptimizationIterations;j++) {
			optimizeMotifMatrix(motifs[i]);
			if (motifs[i]->badFlag) {
				fprintf(stderr, "\t\t!! Motif no longer valid - masked out of existance...\n");
				break;
			}
			motifs[i]->getConsensus();
			fprintf(stderr, "\t\tRound %d: %.2lf", j+1, motifs[i]->logp);
			fprintf(stderr, " %s ", motifs[i]->consensus);
			fprintf(stderr, "T:%.1lf(%.2f%%),B:%.1lf(%.2f%%),P:1e%d\n",motifs[i]->numTarget, 
									motifs[i]->freqTarget*100, motifs[i]->numBackground,
									motifs[i]->freqBackground*100,(int)(motifs[i]->logp/2.3026));
			if (curLogp-FLOAT_EPS < motifs[i]->logp) {
				break;
			}
			curLogp = motifs[i]->logp;
		}
		if (motifs[i]->badFlag) {
			continue;
		}
		if (maskMethod == DENOVO_MASK_FULL) {
			//first optimize motif based on actual zoops scoring
			fullSequences->optimizeMotifThreshold(motifs[i],strand);
		}

		motifs[i]->getConsensus();
		fprintf(stderr, "\t\t=Final=: %.2lf", motifs[i]->logp);
		fprintf(stderr, " %s ", motifs[i]->consensus);
		fprintf(stderr, "T:%.1lf(%.2f%%),B:%.1lf(%.2f%%),P:1e%d\n",motifs[i]->numTarget, 
						motifs[i]->freqTarget*100, motifs[i]->numBackground,
						motifs[i]->freqBackground*100,(int)(motifs[i]->logp/2.3026));
		motifs[i]->print(tmpfp);
		if (tmpfp != NULL) fflush(tmpfp);


		if (maskMethod == DENOVO_MASK_FULL) {
			//remove sequences adjacent to motif from the original sequence that overlap
			fprintf(stderr, "\t\tPerforming exhaustive masking of motif...\n");
			fullSequences->maskMotifFromOligoArray(motifs[i],oligoArray,strand);

			//resort active Oligos to 
			fprintf(stderr, "\t\tReprioritizing potential motifs...\n");
			for (int i=0;i<oligoArray->numActiveOligos;i++) {
				if (oligoArray->activeOligos[i] == NULL) {
					oligoArray->activeOligos[i] = oligoArray->activeOligos[oligoArray->numActiveOligos-1];
					oligoArray->numActiveOligos--;
				} else {
					break;
				}
			}
			oligoArray->initializeOligoPvalues(scorer);
			findMisMatchEnrichment(maxMisMatches,NULL,HOMER_VERBOSE_LEVEL0);	
			seedOligoIndex=0;

		} else {

			//remove oligos from the tree to avoid finding the same motif again
			//remove bound oligos from data set
			oligoArray->scoreOligosWithMotif(motifs[i],motifs[i]->minimumScore,MOTIF_REMOVE_OLIGOS,defaultCPU);
	
			int totalBoundOligos=0;
			for (int j=0;j<oligoArray->numUtilityOligos[0];j++) {
				if (oligoArray->utilityOligos[0][j]->value > motifs[i]->threshold) {
					boundOligos->addOligo(oligoArray->utilityOligos[0][j]->seq,1.0,1.0,-1.0);
					totalBoundOligos++;
				}
			}
			fprintf(stderr, "\t\t%d of %d - %s Bound Oligos: %d\n",i+1,totalMotifsToFind,
											motifs[i]->consensus,totalBoundOligos);

		}


	}
	if (tmpfp != NULL) fclose(tmpfp);

	delete []bufferSeq;
	delete boundOligos;
}

void DeNovoMotifs::printMotifs(FILE* fp) {
	if (numMotifs < 1) return;

	char* name = new char[10000];	
	qsort(motifs, numMotifs, sizeof(Motif*), &(Motif::motifCmp));
	for (int i=0;i<numMotifs;i++) {
		motifs[i]->getConsensus();
		sprintf(name,"%d-%s", i+1, motifs[i]->consensus);
		motifs[i]->setName(name);
		motifs[i]->print(fp);
	}
	delete []name;
}
void DeNovoMotifs::scoreWithKnownMotifs(int optFlag) {

	KnownMotifs* known = new KnownMotifs();
	known->seqFile = seqFile;
	known->groupFile = groupFile;
	known->strand = strand;
	known->numCPUs = numCPUs;
	known->scoringMethod = scoringMethod;
	known->motifs = motifs;
	known->numMotifs = numMotifs;
	known->threasholdOptimizationFlag = optFlag;
	//known->siteReduceThreshold = siteReduceThreshold;

	known->loadSequence();
	known->getMotifEnrichment();
	if (known->siteReduceThreshold > 0) known->findOverlappingMotifs();

	motifs = known->motifs;
	numMotifs = known->numMotifs;
}

//============================ OligoArray ========================================

OligoArray::OligoArray() {
	tree = NULL;
	numTree = 0;
	maxTree = 0;
	alphaSize = DNA_ALPHA;
	oligoLength = 0;
	alpha = NULL;
	rvAlpha = NULL;
	rvMap = NULL;
	numCPUs = 1;

	pthread_mutex_init(&sortMutex, NULL);
	
	staticOligos = NULL;
	utilityOligos = NULL;
	numUtilityOligos = 0;
	numOligos = 0;
	activeOligos = NULL;
	numActiveOligos = 0;
	maxOligos = 0;

	rvFlag = 1;	

	indexStack = NULL;
	alphaStack = NULL;
	positionStack = NULL;
	fwdMisMatchStack = NULL;
	rvMisMatchStack = NULL;

	totalTarget = 0.0;	
	totalBackground = 0.0;	
	
	maxAutoNormIters = MAX_OLIGO_NORMALIZATION_ITERATIONS;
}

OligoArray::~OligoArray() {
	if (tree != NULL) delete []tree;
	if (staticOligos != NULL) {
		for (int i=0;i<numOligos;i++) {
			delete staticOligos[i];
		}
		delete []staticOligos;
	}
	if (alpha != NULL) delete []alpha;
	if (rvAlpha != NULL) delete []rvAlpha;
	if (rvMap != NULL) delete []rvMap;
	if (alphaIndex != NULL) delete []alphaIndex;
	if (rvAlphaIndex != NULL) delete []rvAlphaIndex;

	pthread_mutex_destroy(&sortMutex);

	if (numUtilityOligos != NULL) delete []numUtilityOligos;

	for (int i=0;i<numCPUs;i++) {	
		if (utilityOligos != NULL && utilityOligos[i] != NULL) delete [](utilityOligos[i]);
		if (indexStack != NULL && indexStack[i] != NULL) delete [](indexStack[i]);
		if (alphaStack != NULL && alphaStack[i] != NULL) delete [](alphaStack[i]);
		if (positionStack != NULL && positionStack[i] != NULL) delete [](positionStack[i]);
		if (fwdMisMatchStack != NULL && fwdMisMatchStack[i] != NULL) delete [](fwdMisMatchStack[i]);
		if (rvMisMatchStack != NULL && rvMisMatchStack[i] != NULL) delete [](rvMisMatchStack[i]);
		if (fwdScoreStack != NULL && fwdScoreStack[i] != NULL) delete [](fwdScoreStack[i]);
		if (rvScoreStack != NULL && rvScoreStack[i] != NULL) delete [](rvScoreStack[i]);
	}
	if (utilityOligos != NULL) delete []utilityOligos;
	if (indexStack != NULL) delete []indexStack;
	if (alphaStack != NULL) delete []alphaStack;
	if (positionStack != NULL) delete []positionStack;
	if (fwdMisMatchStack != NULL) delete []fwdMisMatchStack;
	if (rvMisMatchStack != NULL) delete []rvMisMatchStack;
	if (fwdScoreStack != NULL) delete []fwdScoreStack;
	if (rvScoreStack != NULL) delete []rvScoreStack;
}

void OligoArray::init(int length, int totalDataSetBp, int alphaSIZE, int strandFlag, int numcpu) {

	alphaSize = alphaSIZE;
	oligoLength = length;
	if (alphaSize == DNA_ALPHA) {
		alpha=new char[4];
		rvAlpha=new char[4];
		alpha[0] = 'A';
		alpha[1] = 'C';
		alpha[2] = 'G';
		alpha[3] = 'T';
		rvAlpha[0] = 'T';
		rvAlpha[1] = 'G';
		rvAlpha[2] = 'C';
		rvAlpha[3] = 'A';
		rvMap = new char[256];
		alphaIndex = new char[256];
		rvAlphaIndex = new char[256];
		rvMap[65]='T'; //rvMap['A']='T';
		rvMap[67]='G'; //rvMap['C']='G';
		rvMap[71]='C'; //rvMap['G']='C';
		rvMap[84]='A'; //rvMap['T']='A';
		alphaIndex[65]=0; //alphaIndex['A']=0;
		alphaIndex[67]=1; //alphaIndex['C']=1;
		alphaIndex[71]=2; //alphaIndex['G']=2;
		alphaIndex[84]=3; //alphaIndex['T']=3;
	}
	Motif::alphaIndex = alphaIndex;
	Motif::alphaSize = alphaSize;

	if (strandFlag == MOTIF_STRAND_BOTH) {
		rvFlag =1;
	} else {
		rvFlag = 0;
	}	
	
	maxOligos = INT_MAX;
	if (length < 20) {
		double dmax = pow(alphaSize,length)+10;
		if (rvFlag) {
			dmax = dmax/1.9;
		}
		if (dmax < (double)maxOligos) {
			maxOligos = (int) dmax;
		}
	}
	if (maxOligos > totalDataSetBp) {
		maxOligos = totalDataSetBp;
	}
	//int fullTreeEstimate = (int) (log(maxOligos)/log(4.0)+0.49);

	//maxTree = maxOligos*(3+oligoLength-fullTreeEstimate);
	maxTree = maxOligos*(5);
	numTree = 0;
	tree = new int[maxTree];
	for (int i=0;i<maxTree;i++) tree[i] = 0;

	staticOligos = new Oligo*[maxOligos+1];
	staticOligos[0] = NULL;
	for (int i=0;i<maxOligos;i++) {
		staticOligos[i] = NULL;
	}
	staticOligos[maxOligos] = NULL;
	numOligos = 0;

	numCPUs = numcpu;
	if (numCPUs < 1) numCPUs = 1;

	utilityOligos = new Oligo**[numCPUs];
	numUtilityOligos = new int[numCPUs];

	indexStack = new int*[numCPUs];
	alphaStack = new int*[numCPUs];
	positionStack = new int*[numCPUs];
	fwdMisMatchStack = new int*[numCPUs];
	rvMisMatchStack = new int*[numCPUs];
	fwdScoreStack = new double*[numCPUs];
	rvScoreStack = new double*[numCPUs];
	for (int i=0;i<numCPUs;i++) {
		if (i==0) {
			utilityOligos[i] = new Oligo*[maxOligos];
			for (int j=0;j<maxOligos;j++) utilityOligos[i][j] = NULL;
		} else {
			utilityOligos[i] = NULL;
		}

		numUtilityOligos[i]=0;
		indexStack[i] = new int[oligoLength];
		alphaStack[i] = new int[oligoLength];
		positionStack[i] = new int[oligoLength];
		fwdMisMatchStack[i] = new int[oligoLength];
		rvMisMatchStack[i] = new int[oligoLength];
		fwdScoreStack[i] = new double[oligoLength];
		rvScoreStack[i] = new double[oligoLength];
	}

	//fprintf(stderr, "\tEstimate of prefix tree depth = %d\n", fullTreeEstimate);
	//fprintf(stderr, "\tMax # of Oligos               = %d\n", maxOligos);
	//fprintf(stderr, "\tMax size of Tree              = %d\n", maxTree);

}
void OligoArray::optimizeMemoryUsage() {
	fprintf(stderr, "\n\tOligos: %d of %d max\n", numOligos, maxOligos);
	fprintf(stderr, "\tTree  : %d of %d max\n", numTree, maxTree);
	fprintf(stderr, "\tOptimizing memory usage...\n");

	maxTree = numTree+4;
	int* ntree = new int[maxTree];
	for (int i=0;i<maxTree;i++) {
		ntree[i] = tree[i];
	}
	delete []tree;
	tree = ntree;

	maxOligos = numOligos;
	Oligo** tmpOligos = new Oligo*[maxOligos];
	for (int i=0;i<maxOligos;i++) {
		tmpOligos[i] = staticOligos[i+1];
	}
	delete []staticOligos;
	staticOligos = new Oligo*[maxOligos+1];
	staticOligos[0]=0;
	for (int i=0;i<maxOligos;i++) {
		staticOligos[i+1] = tmpOligos[i];
		staticOligos[i+1]->cpuData = new double[numCPUs];;
	}

	for (int i=0;i<numCPUs;i++) {
		if (utilityOligos[i] != NULL) {
			delete [](utilityOligos[i]);
			utilityOligos[i] = NULL;
		}
		utilityOligos[i] = new Oligo*[maxOligos];
		for (int j=0;j<maxOligos;j++) {
			utilityOligos[i][j] = tmpOligos[j];
		}
	}
	delete []tmpOligos;
}


void OligoArray::readOligoCountFile(char* file,int strand, int ncpus) {
	fprintf(stderr, "\n\tReading oligos from file %s\n", file);
	FILE* fp = fopen(file,"r");
	if (fp == NULL) {
		fprintf(stderr, "!!! Couldn't open oligo read count file: %s !!!\n", file);
		return;
	}
	char* buf = new char[MOTIF2_BUFFER];
	char** cols = new char*[1000];
	int numCols;

	int firstOligoSize = -1;
	int numLines = 0;
	while (fgets(buf, MOTIF2_BUFFER, fp) != NULL) {
		splitMotif2(buf,cols,numCols,'\t');
		if (numCols < 3) continue;
		int ok = checkSequence(cols[0]);
		if (!ok) continue;
		numLines++;
		int curLength = strlen(cols[0]);
		if (firstOligoSize == -1) {	
			fprintf(stderr, "\t\tDetected Oligos of length %d from file %s\n", curLength,file);
			firstOligoSize = curLength;
		} else if (curLength != firstOligoSize) {
			fprintf(stderr, "\t\t!!! Error - different sized oligos detected (first one %d != %d) (%s)!!!\n",
									firstOligoSize,curLength,cols[0]);
		}
	}
	init(firstOligoSize,firstOligoSize*numLines,DNA_ALPHA,strand,ncpus);
	rewind(fp);
	double TV=0.0;
	double BV=0.0;
	while (fgets(buf, MOTIF2_BUFFER, fp) != NULL) {
		splitMotif2(buf,cols,numCols,'\t');
		if (numCols < 3) continue;
		int ok = checkSequence(cols[0]);
		if (!ok) continue;
		sscanf(cols[1],"%lf",&TV);
		sscanf(cols[2],"%lf",&BV);
		addOligo(cols[0],TV,BV,0.0);
	}
	fclose(fp);

	optimizeMemoryUsage();
	initializeActiveOligoArray();
	calculateFoldEnrichment(0.5);
}



void OligoArray::printOligo(char* seq, FILE* fp) {
	double vt=0.0,vbg=0.0;
	int cpu = 0;
	int found = searchOligo(seq,0,vt,vbg,0,cpu);
	if (found > 0) {
		utilityOligos[cpu][0]->print(fp);
		fprintf(stderr, "\t0\t%lf\t%lf\n", vt,vbg);
		searchOligo(seq,1,vt,vbg,0,cpu);
		fprintf(stderr, "\t1\t%lf\t%lf\n", vt,vbg);
		searchOligo(seq,2,vt,vbg,0,cpu);
		fprintf(stderr, "\t2\t%lf\t%lf\n", vt,vbg);
	} else {
		fprintf(stderr, "Can't find %s\n", seq);
	}
}

void OligoArray::sortActiveOligos() {
	if (activeOligos == NULL || numActiveOligos < 1) return;
	qsort(activeOligos, numActiveOligos, sizeof(Oligo*), &(Oligo::oligoCmp));
}
void OligoArray::sortUtilityOligos(int cpu) {
	pthread_mutex_lock(&sortMutex);
	Oligo::cpu = cpu;
	qsort(utilityOligos[cpu], numUtilityOligos[cpu], sizeof(Oligo*), &(Oligo::oligoCmpReverseParallel));
	pthread_mutex_unlock(&sortMutex);
}
void OligoArray::initializeActiveOligoArray() {
	if (activeOligos != NULL) delete []activeOligos;
	activeOligos = new Oligo*[numOligos];
	numActiveOligos = numOligos;
	for (int i=0;i<numOligos;i++) {
		activeOligos[i] = staticOligos[i+1];
		activeOligos[i]->flag=0;
		activeOligos[i]->flag2=0;
	}
}

int OligoArray::addOligo(char* seq, double vt, double vbg, double addValue) {
	totalTarget += vt;
	totalBackground += vbg;
	double tmpT=0.0,tmpBg=0.0;
	int cpu = 0;
	int found = searchOligo(seq,0,tmpT,tmpBg,OLIGO_SEARCH_KEEP,cpu);
	int rvMode = 0;
	if (found == 0) {
		rvMode = 0;
	} else {
		//if (found > 1) fprintf(stderr, "!!!!Problem %s\t%s\n", seq, utilityOligos[cpu][0]->seq);
		if (strcmp(seq,utilityOligos[cpu][0]->seq)==0) {
			rvMode=0;
		} else {
			rvMode=1;
		}
		if (rvMode == 1 && rvFlag == 0) {
			//fprintf(stderr, "!!!!WTF %s\t%s\n", seq, utilityOligos[0]->seq);
		}
	}
	int added = addOligoTree(seq,vt,vbg,addValue,rvMode,0);
	rvMode = 0;
	return added;
}

int OligoArray::addOligoTree(char* seq, double vt, double vbg, double addValue, int rvMode, int knownIndex) {

	int index = 0;
	int lastIndex = 0;
	int same = 0;
	int knownIndexToAdd=0;

	for (int i=0;i<oligoLength;i++) {
		char bp = '\0';
		if (rvMode) {
			bp = seq[oligoLength-i-1];
		} else {
			bp = seq[i];
		}

		if (index == numTree) {
			if (numTree+alphaSize > maxTree) {
				fprintf(stderr, "Out of room in oligo TREE! %d of %d\n", numTree+alphaSize, maxTree);
			}
			for (int j=0;j<alphaSize;j++) {
				tree[numTree]=0;
				numTree++;
			}
		} else if (index > numTree) {
			fprintf(stderr, "ERROR IN TREE GENERATION!!\n");
		}

		int indexDelta = alphaIndex[(int)bp];
		if (rvMode) indexDelta = 3-indexDelta;

		lastIndex = index+indexDelta;
		index = tree[lastIndex];

		if (index < 0 && i< oligoLength-1) {
			//this is a pointer to an oligo
			//check if equal
			same=1;
			char* otherSeq = staticOligos[-index]->seq;
			//fprintf(stderr, "\tChecking if %s = %s rvMode=%d starting at %d\n", seq, otherSeq,rvMode,i);
			for (int j=i+1;j<oligoLength;j++) {
				if (rvMode) {
					if (rvMap[(int)seq[oligoLength-j-1]] == otherSeq[j]) same++;
					else break;
				} else {
					if (seq[j] == otherSeq[j]) same++;
					else break;
				}
			}
			if (i+same < oligoLength) {
				//	fprintf(stderr, " !!!!! NOT THE SAME\n");
				//this means they aren't equal 
				knownIndexToAdd=-index;
				index = 0; //set index to zero so that a new tree node is created
			} else {
				//	fprintf(stderr, " ARE THE SAME\n");
				//sequences are the same - break and add the counts below
				//fprintf(stderr, "Sequences are the same!\n");
				break;
			}
		}
		if (index == 0) {
			if (i<oligoLength-1 && same > 0) {
				tree[lastIndex] = numTree;
				tree[numTree] = 0;
				index = numTree;
				same--;
			} else {
				//sequences haven't been added yet-
				break;
			}
		}
	}

	int oligoIndex = -tree[lastIndex];
	int added = 0;
	if (knownIndex != 0) {
		oligoIndex = knownIndex;
		tree[lastIndex]= -knownIndex;
	} else {
		if (oligoIndex == 0) {
			oligoIndex = ++numOligos;
			if (oligoIndex > maxOligos-1) {
				fprintf(stderr, "Ran out of room for more oligos %d of %d!!!!\n",oligoIndex,maxOligos);
			}
			tree[lastIndex] = -oligoIndex;
			if (rvMode) {
				fprintf(stderr, "THIS SHOULDN'T Happen!!!! (%s)\n",seq);
			}
			staticOligos[oligoIndex] = new Oligo(seq);
			//oligos[oligoIndex-1] = staticOligos[oligoIndex]; not needed
			added=1;
		} else if (oligoIndex < 0) {
			fprintf(stderr, "ERROR WITH OLIGO CREATION\n");
		} else {
		}
		if (addValue <= FLOAT_EPS || staticOligos[oligoIndex]->value + FLOAT_EPS < addValue) {
			staticOligos[oligoIndex]->numTarget += vt;
			staticOligos[oligoIndex]->numBackground += vbg;
		}
		staticOligos[oligoIndex]->value = addValue;
		
	}

	if (knownIndexToAdd != 0) {
		addOligoTree(staticOligos[knownIndexToAdd]->seq,staticOligos[knownIndexToAdd]->numTarget,
						staticOligos[knownIndexToAdd]->numBackground,addValue,rvMode,knownIndexToAdd);
	}
		

	return added;
}

void OligoArray::adjustOligoInstances(char* seq, double vt, double vbg) {
	double tmp1, tmp2;
	int cpu = 0;
	int maxMisMatches = 0;
	int found = searchOligo(seq, maxMisMatches, tmp1, tmp2, OLIGO_SEARCH_KEEP,cpu);
	if (found < 1) return;
	if (found > 1) {
		fprintf(stderr, "!!! Somthing is probably wrong... duplicates in the oligoarray\n");
	}
	double normWeight = utilityOligos[cpu][0]->normWeight;
	utilityOligos[cpu][0]->numTarget += vt;
	utilityOligos[cpu][0]->numBackground += vbg*normWeight;
	int remove = 1;
	if (utilityOligos[cpu][0]->numTarget < FLOAT_EPS) {
		utilityOligos[cpu][0]->numTarget = 0.0;
	} else {
		remove = 0;
	}
	if (utilityOligos[cpu][0]->numBackground < FLOAT_EPS) {
		utilityOligos[cpu][0]->numBackground = 0.0;
	} else {
		remove = 0;
	}
	if (remove) {
		//all instances of the oligo have been removed, so get it out of the tree
		found = searchOligo(seq, maxMisMatches, tmp1, tmp2, OLIGO_SEARCH_REMOVE,cpu);
	}
}

int OligoArray::searchOligo(char* seq,int maxMisMatches, double& vt, double& vbg, 
															int removeFlag,int cpu) {
	
	numUtilityOligos[cpu] = 0;

	indexStack[cpu][0] = 0;
	alphaStack[cpu][0] = 0;
	fwdMisMatchStack[cpu][0] = 0;
	rvMisMatchStack[cpu][0] = 0;
	int stackPosition = 0;
	int valIndex = 0;
	int oligoIndex = 0;
	
	int alphaPos = 0;
	int fwdMisMatch = 0;
	int rvMisMatch = 0;

	while (stackPosition > -1) {


		alphaPos = alphaStack[cpu][stackPosition];

		if (alphaPos >= alphaSize) {
			if (removeFlag == OLIGO_SEARCH_REMOVE) {
				int numEmpty = 0;
				for (int i=0;i<alphaSize;i++) {
					if (tree[indexStack[cpu][stackPosition]+i] == 0) numEmpty++;
				}
				if (numEmpty == alphaSize && stackPosition > 0) {
					tree[indexStack[cpu][stackPosition-1]+alphaStack[cpu][stackPosition-1]]=0;
				}
			}
			stackPosition--;
			continue;
		}


		fwdMisMatch = fwdMisMatchStack[cpu][stackPosition];
		rvMisMatch = rvMisMatchStack[cpu][stackPosition];

		if (rvFlag) {
			if ('N' != seq[oligoLength-stackPosition-1]) {
				if (rvAlpha[alphaPos] != seq[oligoLength-stackPosition-1])
					rvMisMatch++;
			}
		}
		if ('N' != seq[stackPosition]) {
			if (alpha[alphaPos] != seq[stackPosition]) fwdMisMatch++;
		}

		if (fwdMisMatch > maxMisMatches) {
			if (!rvFlag || rvMisMatch > maxMisMatches) {
				alphaStack[cpu][stackPosition]++;
				continue;
			}
		}

		valIndex = tree[indexStack[cpu][stackPosition]+alphaPos];
		if (valIndex < 0) {
			oligoIndex = -valIndex;
			//this is a pointer to an oligo
			//check if equal
			char* otherSeq = staticOligos[oligoIndex]->seq;
	
			for (int i=stackPosition+1;i<oligoLength;i++) {
				if (seq[i] != 'N' && seq[i] != otherSeq[i]) {
					fwdMisMatch++;
					if (fwdMisMatch > maxMisMatches) break;
				}
			}
			if (rvFlag) {
				for (int i=stackPosition+1;i<oligoLength;i++) {
					if (seq[oligoLength-i-1] != 'N' && rvMap[(int)seq[oligoLength-i-1]] != otherSeq[i]) {
						rvMisMatch++;
						if (rvMisMatch > maxMisMatches) break;
					}
				}
			}

			if (fwdMisMatch <= maxMisMatches || (rvFlag && rvMisMatch <= maxMisMatches)) {
				vt += staticOligos[oligoIndex]->numTarget;
				vbg += staticOligos[oligoIndex]->numBackground;
				staticOligos[oligoIndex]->flag = 1;
				utilityOligos[cpu][numUtilityOligos[cpu]++] = staticOligos[oligoIndex];
//fprintf(stderr, "\t%s\t%lf\t%lf\n", otherSeq,vt,vbg);
				if (removeFlag == OLIGO_SEARCH_REMOVE) {
					tree[indexStack[cpu][stackPosition]+alphaPos] = 0;
				}
			}
			alphaStack[cpu][stackPosition]++;
			continue;
		} else if (valIndex == 0) {
			//dead end
			alphaStack[cpu][stackPosition]++;
			continue;
		} else {
			//next node
			alphaStack[cpu][stackPosition+1]=0;
			indexStack[cpu][stackPosition+1]=valIndex;
			fwdMisMatchStack[cpu][stackPosition+1]=fwdMisMatch;
			rvMisMatchStack[cpu][stackPosition+1]=rvMisMatch;

			alphaStack[cpu][stackPosition]++;
			stackPosition++;

			continue;
		}
	}
	return numUtilityOligos[cpu];
} 

void OligoArray::scoreOligosWithMotif(Motif* motif, double minimumScore, int removeFlag, int cpu) {

	motif->logxform();
	numUtilityOligos[cpu]=0;
	
	indexStack[cpu][0] = 0;
	alphaStack[cpu][0] = 0;
	positionStack[cpu][0] = oligoLength-1;
	fwdScoreStack[cpu][0] = 0.0;
	rvScoreStack[cpu][0] = 0.0;

	int stackPosition = 0;
	int valIndex = 0;
	int oligoIndex = 0;
	
	int alphaPos = 0;

	int fwdPosition=0;
	int rvPosition=0;
	double fwdScore=0.0;
	double rvScore=0.0;

	while (stackPosition > -1) {

		alphaPos = alphaStack[cpu][stackPosition];
		if (alphaPos >= alphaSize) {
			if (removeFlag == MOTIF_REMOVE_OLIGOS) {
				int numEmpty = 0;
				for (int i=0;i<alphaSize;i++) {
					if (tree[indexStack[cpu][stackPosition]+i] == 0) numEmpty++;
				}
				if (numEmpty == alphaSize && stackPosition > 0) {
					tree[indexStack[cpu][stackPosition-1]+alphaStack[cpu][stackPosition-1]]=0;
				}
			}
			stackPosition--;
			continue;
		}

		fwdPosition = stackPosition;
		rvPosition = oligoLength-fwdPosition-1;

		fwdScore = fwdScoreStack[cpu][stackPosition];
		rvScore = rvScoreStack[cpu][stackPosition];

		fwdScore += motif->m[fwdPosition][alphaPos];
		rvScore += motif->m[rvPosition][3-alphaPos];

		if (fwdPosition < oligoLength-1) {
			if (fwdScore+motif->fwdMaxPossible[fwdPosition+1] < minimumScore) {
				if (!rvFlag || rvScore+motif->rvMaxPossible[rvPosition-1] < minimumScore) {
					alphaStack[cpu][stackPosition]++;
					continue;
				}
			}
		}

		valIndex = tree[indexStack[cpu][stackPosition]+alphaPos];
		if (valIndex < 0) {
			oligoIndex = -valIndex;
			//this is a pointer to an oligo
			char* seq = staticOligos[oligoIndex]->seq;
			for (int i=fwdPosition+1;i<oligoLength;i++) {
				char bp = seq[i];
				fwdScore += motif->m[i][(int)alphaIndex[(int)bp]];
				rvScore += motif->m[oligoLength-i-1][3-alphaIndex[(int)bp]];
			}
			//fprintf(stderr, "%s\t%lf\t%lf\n", seq, fwdScore, rvScore);
			double bestScore = fwdScore;
			char direction = 0;
			if (fwdScore < rvScore && rvFlag) {
				bestScore = rvScore;
				direction = 1;
			}
			staticOligos[oligoIndex]->cpuData[cpu] = bestScore;
			staticOligos[oligoIndex]->value = bestScore;
			staticOligos[oligoIndex]->flag = direction;
	
			if (removeFlag == MOTIF_REMOVE_OLIGOS) {
				if (bestScore >= motif->threshold) {
					//fprintf(stderr, "\t\t\t%s\n", staticOligos[oligoIndex]->seq);
					tree[indexStack[cpu][stackPosition]+alphaPos] = 0;
				}
			}
			utilityOligos[cpu][numUtilityOligos[cpu]] = staticOligos[oligoIndex];
			numUtilityOligos[cpu]++;
			alphaStack[cpu][stackPosition]++;
			continue;
		} else if (valIndex == 0) {
			//dead end
			alphaStack[cpu][stackPosition]++;
			continue;
		} else {
			//next node
			alphaStack[cpu][stackPosition+1]=0;
			indexStack[cpu][stackPosition+1]=valIndex;
			fwdScoreStack[cpu][stackPosition+1]=fwdScore;
			rvScoreStack[cpu][stackPosition+1]=rvScore;

			alphaStack[cpu][stackPosition]++;

			stackPosition++;

			continue;
		}
	}
}

void OligoArray::calculateFoldEnrichment(double pseudoCount) {

	double ratio = 1;
	if (totalBackground > 0) {
		ratio = totalTarget / totalBackground;
	}
	double AAA = 0.0;
	fprintf(stderr, "\tCalculating Oligo Enrichment\n");
	fprintf(stderr, "\t\tTotal Target Oligos = %lf\n\t\tTotal Background Oligos = %lf\n", totalTarget, totalBackground);
	for (int i=0;i<numActiveOligos;i++) {
		double nt = activeOligos[i]->numTarget;
		double nb = activeOligos[i]->numBackground*ratio;
		AAA += activeOligos[i]->numBackground;
		if (nt < pseudoCount) nt = pseudoCount;
		if (nb < pseudoCount) nb = pseudoCount;
		activeOligos[i]->value = nt/nb;
	}
}

void OligoArray::normalizeOligos(int nlen,char* normFile) {
	if (nlen < 1) return;
	double lastError = 1e100;
	double cur = 0.0;
	double next = 1.0;
	for (int i=0;i<maxAutoNormIters;i++) {
		int printFlag = 0;
		if (i==0) printFlag = 1;
		double curError = normalizeOligosIteration(nlen,normFile,printFlag);
		double improvementRate = (lastError-curError)/lastError;
		//fprintf(stderr, "Round %d:\t%le\t%le\t%le\n", i+1,lastError,curError,improvementRate);
		if (curError >= lastError) break;
		if (improvementRate < MIN_OLIGO_NORMALIZATION_IMPROVEMENT) break;
		lastError = curError;
		if (printFlag) {
			fprintf(stderr, "\tConverging on autonormalization solution:\n\t");
		} else {
			cur += 80.0/(double)maxAutoNormIters;
			if (cur > next) {
				fprintf(stderr, ".");
				next += 1.0;
			}
		}
	}
	fprintf(stderr, "\n");
	fprintf(stderr, "\n\tFinal normalization:\n");
	(void)normalizeOligosIteration(nlen,NULL,1);
}

double OligoArray::normalizeOligosIteration(int nlen,char* normFile, int printFlag) {

	FILE* fp = NULL;
	if (normFile != NULL && printFlag) {
		fp = fopen(normFile, "w");
		if (fp == NULL) {
			fprintf(stderr, "!!! Could not open %s for writing!!!\n",normFile);
		} else {
			fprintf(fp, "Oligo\tTargetCounts\tBackgroundCounts\tNormalizationFactor\n");
		}
	}

	double error = 0.0;

	char* mer = new char[1000];

	for (int curLen=1;curLen<=nlen;curLen++) {
		int total = (int) pow(4.0,(double)curLen);
		if (printFlag) fprintf(stderr, "\tAutonormalization: %d-mers (%d total)\n",curLen,total);
		double* targetLevels = new double[total];
		double* backgroundLevels = new double[total];
		double* normFactors = new double[total];
		int* pow4 = new int[curLen+1];
		pow4[0]=1;
		for (int i=1;i<=curLen;i++) pow4[i] = 4*pow4[i-1];
		int nbins = oligoLength-curLen+1;
		for (int i=0;i<total;i++) {
			targetLevels[i]=0.0;
			backgroundLevels[i]=0.0;
			normFactors[i] = 1.0;
		}

		double dnbins = (double)nbins;
		for (int i=1;i<numOligos+1;i++) {
			char* s = staticOligos[i]->seq;
			for (int j=0;j<nbins;j++) {
				int index = 0;
				for (int k=0;k<curLen;k++) {
					index += alphaIndex[(int)s[j+k]]*pow4[k];
				}
				targetLevels[index] += staticOligos[i]->numTarget/dnbins;
				backgroundLevels[index] += staticOligos[i]->numBackground/dnbins;
			}
		}

		double totalTarget = 0.0;
		double totalBackground = 0.0;
		for (int i=0;i<total;i++) {
			totalTarget += targetLevels[i];
			totalBackground += backgroundLevels[i];
		}
		double minimumBackgroundLevels = 0.5/totalBackground;
		double minimumTargetLevels = 0.5/totalTarget;
		for (int i=0;i<total;i++) {
			for (int j=0;j<curLen;j++) {
				int remainder = i % pow4[j+1];
				int asdf = (int)(remainder/pow4[j]);
				mer[j] = alpha[asdf];
			}
			mer[curLen]='\0';

			double tLevel = targetLevels[i];
			double bLevel = backgroundLevels[i];
			if (rvFlag) {
				int rindex = 0;
				for (int j=curLen-1;j>=0;j--) {
					rindex += (3-alphaIndex[(int)mer[j]])*pow4[curLen-j-1];
				}
				tLevel = (tLevel+targetLevels[rindex])/2.0;
				bLevel = (bLevel+backgroundLevels[rindex])/2.0;
			}

			if (bLevel < minimumBackgroundLevels) bLevel = minimumBackgroundLevels;
			if (tLevel < minimumTargetLevels) tLevel = minimumTargetLevels;

			double factor = totalBackground*tLevel/totalTarget/bLevel;
			if (printFlag) {
				if (curLen < 3) {
					fprintf(stderr, "\t\t%s\t%.2lf%%\t%.2lf%%\t%.3lf\n", mer,
									100.0*tLevel/totalTarget,
									100.0*bLevel/totalBackground,factor);
				}
				if (fp != NULL) {
					//fprintf(fp, "%s\t%.1lf\t%.1lf\t%.3lf\n", mer,tLevel,bLevel,factor);
					fprintf(fp, "\t\t%s\t%.2lf%%\t%.2lf%%\t%.3lf\n", mer,
									100.0*tLevel/totalTarget,
									100.0*bLevel/totalBackground,factor);
				}
			}
			normFactors[i] = factor;
			error += (factor-1.0)*(factor-1.0)/(double)total;
		}

		totalBackground = 0.0;
		for (int i=1;i<numOligos+1;i++) {
			char* s = staticOligos[i]->seq;
			double scoreUnit = (double)staticOligos[i]->numBackground;
			double newScore = 0;
			for (int j=0;j<nbins;j++) {
				int index = 0;
				for (int k=0;k<curLen;k++) {
					index += alphaIndex[(int)s[j+k]]*pow4[k];
				}
				newScore += normFactors[index];
			}
			newScore /= nbins;
			double score = (double)(newScore*scoreUnit);

			staticOligos[i]->numBackground = score;
			staticOligos[i]->normWeight *= newScore;
			totalBackground += staticOligos[i]->numBackground;
		}
		delete []targetLevels;
		delete []backgroundLevels;
		delete []normFactors;
	}
	if (fp != NULL) {
		fprintf(stderr, "\tNormalization weights can be found in file: %s\n", normFile);
		fclose(fp);
	}
	delete []mer;

	return error;
}
void OligoArray::initializeOligoPvalues(EnrichmentScorer* scorer) {
	//initialize individual oligo enrichments
	for (int i=1;i<numOligos+1;i++) {
		if (staticOligos[i] == NULL) {
			fprintf(stderr, "X");
			continue;
		}
		//fprintf(stderr, "- %x\n", activeOligos[i]);
		//activeOligos[i]->logp = -1*activeOligos[i]->numTarget;
		staticOligos[i]->value = 1.0;
		staticOligos[i]->flag2 = 0;
		staticOligos[i]->logp = (float)scorer->scoreEnrichment(staticOligos[i]->numTarget,
														staticOligos[i]->numBackground);
	}
}


// ----------------------- class PSSM -------------------------------

Motif::Motif() {
	init();
}
Motif::Motif(char* seq) {
	init();
	length = strlen(seq);
	consensus = new char[length+1];
	strcpy(consensus,seq);
	name = new char[length+1];
	strcpy(name,seq);
	m=new double*[length];
	for (int i=0;i<length;i++) {
		m[i]=new double[alphaSize];
		for (int j=0;j<alphaSize;j++) {
			m[i][j] = MOTIF_INITIAL_DEGENERACY;
			if (alphaIndex[(int)seq[i]] == j) {
				m[i][j] = 1.0 - MOTIF_INITIAL_DEGENERACY*3.0;
			}
		}
	}
	fwdMaxPossible = new double[length];
	rvMaxPossible = new double[length];
	for (int i=0;i<length;i++) {
		fwdMaxPossible[i]=0;
		rvMaxPossible[i]=0;
	}
	logFlag = 0;
}
Motif::Motif(Motif* motif) {
	init();
	motif->expxform();
	length = motif->length;
	m=new double*[length];
	for (int i=0;i<length;i++) {
		m[i]=new double[alphaSize];
		for (int j=0;j<alphaSize;j++) {
			m[i][j] = motif->m[i][j];
		}
	}
	fwdMaxPossible = new double[length];
	rvMaxPossible = new double[length];
	for (int i=0;i<length;i++) {
		fwdMaxPossible[i]=0;
		rvMaxPossible[i]=0;
	}
	logFlag = 0;
	numTarget = motif->numTarget;
	numBackground = motif->numBackground;
	freqTarget = motif->freqTarget;
	freqBackground = motif->freqBackground;
	normalize();
}
Motif::Motif(char* n, char* c, double t, double lp, double ** matrix, int len) {
	init();
	name = new char[strlen(n)+1];
	strcpy(name,n);
	consensus = new char[strlen(c)+1];
	strcpy(consensus,c);
	threshold = t; //-0.01;
	logp = lp;
	length = len;
	m = new double*[length];
	for (int i=0;i<length;i++) {
		m[i] = new double[alphaSize];
		for (int j=0;j<alphaSize;j++) {
			m[i][j] = matrix[i][j];
		}
	}
	logFlag = 0;
	//normalize();
}


Motif::~Motif() {
	if (m != NULL) {
		for (int i=0;i<length;i++) delete [](m[i]);
		delete []m;
	}
	if (seqStats != NULL) delete []seqStats;
	if (name != NULL) delete []name;
	if (consensus != NULL) delete []consensus;
	if (fwdMaxPossible != NULL) delete []fwdMaxPossible;
	if (rvMaxPossible != NULL) delete []rvMaxPossible;
}
void Motif::init() {
	m = NULL;
	name = NULL;
	consensus = NULL;
	length = 0;
	threshold = 0;
	logFlag = 0;
	numTarget = 0;
	numBackground = 0;
	freqTarget = 0;
	freqBackground = 0;
	seqStats = NULL;
	numSeqs = 0;
	logp=0;
	minimumScore = 0;
	badFlag =0;

	similarPercent = 0.0;
	similarPercentOfOther = 0.0;
	similarMotif = NULL;

	statAvgTarget = 0.0;
	statStdTarget = 0.0;
	statAvgBackground = 0.0;
	statStdBackground = 0.0;
	statStrandLogRatio = 0.0;
	siteMultiplicity = 0.0;
}
void Motif::setName(char* n) {
	if (name != NULL) delete []name;
	name = NULL;
	if (n == NULL) return;
	name = new char[strlen(n)+1];
	strcpy(name,n);
}
void Motif::addOligo(Oligo* oligo) {
	int rvFlag = oligo->flag;
	double value = -1*oligo->value;

	for (int i=0;i<length;i++) {
		int index = 0;
		if (rvFlag) index = 3-alphaIndex[(int)oligo->seq[length-i-1]];
		else index =alphaIndex[(int)oligo->seq[i]];
	
		m[i][index] += value;
	}
}
void Motif::initializeSeqStats(int num) {
	//if (seqStats != NULL) {
	if (seqStats != NULL && numSeqs != num) {
		delete []seqStats;
		seqStats = NULL;
	}
	numSeqs = num;
	if (seqStats == NULL) seqStats = new SequenceStats[numSeqs];
	for (int i=0;i<numSeqs;i++) {
		seqStats[i].initialize();
	}
}


void Motif::setBlank(int len) {
	length = len;
	m = new double*[length];
	for (int i=0;i<length;i++) {
		m[i] = new double[alphaSize];
		for (int j=0;j<alphaSize;j++) {
			m[i][j] = 0.0;
		}
	}
	fwdMaxPossible = new double[length];
	rvMaxPossible = new double[length];
	for (int i=0;i<length;i++) {
		fwdMaxPossible[i]=0;
		rvMaxPossible[i]=0;
	}
	logFlag = 0;
}

void Motif::logxform() {
	if (!logFlag) {
		for (int i=0;i<length;i++) {
			for (int j=0;j<alphaSize;j++) {
				if (m[i][j] < MOTIF_PROBABILITY_MINIMUM)
					m[i][j] = MOTIF_PROBABILITY_MINIMUM;
				m[i][j] = log(m[i][j]/0.25);
			}
		}
		logFlag =1;
	}
}
void Motif::expxform() {
	if (logFlag) {
		for (int i=0;i<length;i++) {
			for (int j=0;j<alphaSize;j++) {
				m[i][j] = exp(m[i][j])*0.25;
			}
		}
		logFlag = 0;
	}
	justNormalize();
}
void Motif::normalize() {
	justNormalize();
	logxform();
	expxform();
}
void Motif::justNormalize() {
	char tmp[1000];
	for (int i=0;i<length;i++) {
		double sum = 0;
		for (int j=0;j<alphaSize;j++) {
			sum += m[i][j];
		}
		for (int j=0;j<alphaSize;j++) {
			m[i][j] /= sum;
		}
		double totalerror = 0;
		double count = 0;
		for (int j=0;j<alphaSize;j++) {
			if (m[i][j] < MOTIF_PROBABILITY_MINIMUM) {
				totalerror += MOTIF_PROBABILITY_MINIMUM-m[i][j];
				count+=1.0;
			}
		}
		if (count > 0.0) {
			double adjust = totalerror/((double)alphaSize-count);
			for (int j=0;j<alphaSize;j++) {
				if (m[i][j] < MOTIF_PROBABILITY_MINIMUM) {
					m[i][j] = MOTIF_PROBABILITY_MINIMUM;
				} else {
					m[i][j] -= adjust;
				}
			}
		}
		for (int j=0;j<alphaSize;j++) {
			sprintf(tmp,"%.3lf",m[i][j]); // for rounding purposes
			sscanf(tmp,"%lf",&(m[i][j]));
		}
	}
}
void Motif::adjustLength(int newlen) {
	expxform();
	double** matrix = new double*[newlen];
	for (int i=0;i<newlen;i++) {
		matrix[i]=new double[alphaSize];
		for (int j=0;j<alphaSize;j++) matrix[i][j]=0.25;
	}
	if (length <= newlen) {
		int start = (newlen-length)/2;
		for (int i=start;i<length+start;i++) {
			for (int j=0;j<alphaSize;j++) matrix[i][j] = m[i-start][j];
		}
	} else {
		int start = (length-newlen)/2;
		for (int i=start;i<newlen+start;i++) {
			for (int j=0;j<alphaSize;j++) matrix[i-start][j] = m[i][j];
		}
	}
	for (int i=0;i<length;i++) delete [](m[i]);
	delete []m;
	m= matrix;
	length = newlen;

	if (fwdMaxPossible != NULL) delete []fwdMaxPossible;
	fwdMaxPossible = new double[length];
	for (int i=0;i<length;i++) fwdMaxPossible[i] = 0.0;

	if (rvMaxPossible != NULL) delete []rvMaxPossible;
	rvMaxPossible = new double[length];
	for (int i=0;i<length;i++) rvMaxPossible[i] = 0.0;
}



void Motif::print(FILE* fp) {
	if (fp == NULL) return;
	getConsensus();
	fprintf(fp, ">%s\t", consensus);
	if (name != NULL) {
		fprintf(fp, "%s", name);
	} else {
		fprintf(fp, "%s", consensus);
	}
	fprintf(fp, "\t%lf", threshold);
	fprintf(fp, "\t%lf", logp);
	fprintf(fp, "\t0");
	fprintf(fp, "\tT:%.1lf(%.2f%%),B:%.1lf(%.2f%%),P:1e%d",numTarget, freqTarget*100,  
									numBackground,freqBackground*100,(int)(logp/2.3026));
	if (siteMultiplicity > 0.5) {
		fprintf(fp,"\tTpos:%.1lf,Tstd:%.1lf,Bpos:%.1lf,Bstd:%.1lf,StrandBias:%.1lf,Multiplicity:%.2lf",
					statAvgTarget,statStdTarget,statAvgBackground,statStdBackground,
					statStrandLogRatio,siteMultiplicity);
	}
	if (similarMotif != NULL) {
		if (siteMultiplicity < 0.5) fprintf(fp, "\t");
		char* nn = similarMotif->name;
		if (nn == NULL) {
			if (similarMotif->consensus) similarMotif->getConsensus();
			nn = similarMotif->consensus;
		}
		fprintf(fp, "\tShares %.2lf%% of sites with %s(%.2lf%%)", similarPercent*100.0, nn,
														similarPercentOfOther*100.0);
	}
	fprintf(fp, "\n");
	for (int i=0;i<length;i++) {
		for (int j=0;j<alphaSize;j++) {
			if (j>0) fprintf(fp, "\t");
			fprintf(fp, "%.3lf",m[i][j]);
		}
		fprintf(fp, "\n");
	}
}

void Motif::getConsensus() {
	if (consensus != NULL) delete []consensus;
	consensus = new char[length+1];
	consensus[length] = '\0';
	expxform();
	char* alphaPresent = new char[alphaSize];
	for (int i=0;i<length;i++) {
		//int bestIndex = 0;
		double bestValue = 0;
		for (int j=0;j<alphaSize;j++) {
			alphaPresent[j]=0;
			if (m[i][j] > bestValue) {
				//bestIndex = j;
				bestValue = m[i][j];
			}
		}
		for (int j=0;j<alphaSize;j++) {
			double diff = bestValue - m[i][j];
			if (diff < MOTIF_IUPAC_THRESHOLD)
				alphaPresent[j]=1;
		}
		consensus[i] = getSymbol(alphaPresent);
	}
	delete []alphaPresent;
}
char Motif::getSymbol(char* p) {
	char bp = 'N';
	if (alphaSize == DNA_ALPHA) {
		if (p[0]) {
			if (p[1]) {
				if (p[2]) {
					if (p[3]) bp = 'N';
					else bp = 'V';
				} else {
					if (p[3]) bp = 'H';
					else bp = 'M';
				}
			} else {
				if (p[2]) {
					if (p[3]) bp = 'D';
					else bp = 'R';
				} else {
					if (p[3]) bp = 'W';
					else bp = 'A';
				}
			}
		} else {
			if (p[1]) {
				if (p[2]) {
					if (p[3]) bp = 'B';
					else bp = 'S';
				} else {
					if (p[3]) bp = 'Y';
					else bp = 'C';
				}
			} else {
				if (p[2]) {
					if (p[3]) bp = 'K';
					else bp = 'G';
				} else {
					if (p[3]) bp = 'T';
					else bp = 'N';
				}
			}
		}
	}
	return bp;
}

void Motif::setMinimumScore(double percent) {
	logxform();	
	int psize = (int) ceil(1.0/percent);
	if (psize < 2) {
		minimumScore = -1e60;
		return;
	}
	int psize4 = psize*alphaSize;
	double* lastResults=new double[psize];
	double* curResults=new double[psize4];
	for (int i=0;i<psize;i++) lastResults[i]=0.0;

	//set maximum scores, so you know when to quit
	for (int i=0;i<length;i++) {
		double max = 0;
		for (int j=0;j<alphaSize;j++) if (m[i][j] > max) max = m[i][j];
		rvMaxPossible[i]=max;
	}
	fwdMaxPossible[length-1]=rvMaxPossible[length-1];
	for (int i=length-2;i>=0;i--) fwdMaxPossible[i]=fwdMaxPossible[i+1]+rvMaxPossible[i];
	for (int i=1;i<length;i++) rvMaxPossible[i] += rvMaxPossible[i-1];

	//estimate percentile score
	for (int i=0;i<length;i++) {
		for (int j=0;j<psize;j++) {
			for (int k=0;k<alphaSize;k++) {
				curResults[j*alphaSize+k] = m[i][k]+lastResults[j];
			}
		}
		qsort(curResults, psize4, sizeof(double), &sortDouble);
		for (int j=0;j<psize;j++) {
			double sum = 0;
			for (int k=0;k<alphaSize;k++) {
				sum += curResults[j*alphaSize+k];
			}
			sum /= alphaSize;
			lastResults[j] = sum;
		}
	}
	minimumScore = lastResults[psize-2];
	//for (int i=0;i<psize;i++) {
	//	fprintf(stderr, "\t%lf\n", lastResults[i]);
	//}
	//fprintf(stderr, "\tMinimum Score set to %lf\n", minimumScore);
	delete []lastResults;
	delete []curResults;
}
 

int sortDouble(const void *a, const void* b) {
	if (*((double*)a) < *((double*)b)) return -1;
	if (*((double*)a) > *((double*)b)) return 1;
	return 0;
}
int Motif::nameID = 0;
int Motif::motifCmp(const void *a, const void *b) {
	double va=(*((Motif**)a))->logp;
	double vb=(*((Motif**)b))->logp;
	if (va < vb) return -1;
	if (va > vb) return 1;
	return 0;
}

Motif** Motif::readMotifFile(char* file, int &numMotifs) {

	numMotifs=0;
	if (file==NULL) return NULL;
	FILE* fp = fopen(file,"r");
	if (fp == NULL) {
		fprintf(stderr, "!!! Couldn't open motif file: %s !!!\n", file);
		return NULL;
	}
	char* buf = new char[MOTIF2_BUFFER];
	char** cols = new char*[1000];
	int numCols;

	int maxMotifs = 0;
	while (fgets(buf, MOTIF2_BUFFER, fp) != NULL) {
		if (buf[0] == '>') maxMotifs++;
	}

	int maxMotifLength = MAXIMUM_MOTIF_LENGTH;
	double** matrix = new double*[MAXIMUM_MOTIF_LENGTH];
	for (int i=0;i<maxMotifLength;i++) {
		matrix[i] = new double[DNA_ALPHA];
	}

	char* consensus = new char[MAXIMUM_MOTIF_LENGTH];
	char* name = new char[MAXIMUM_MOTIF_LENGTH];
	int curMotifLength = 0;
	consensus[0]='\0';
	name[0]='\0';
	double threshold = 0.0;
	double logp = 0.0;
	int started = 0;

	Motif** motifs = new Motif*[maxMotifs];

	rewind(fp);
	while (fgets(buf, MOTIF2_BUFFER, fp) != NULL) {
		splitMotif2(buf,cols,numCols,'\t');
		if (numCols < 3) continue;
		if (cols[0][0] == '>') {
			if (started) {
				motifs[numMotifs] = new Motif(name,consensus,threshold,
											logp,matrix,curMotifLength);
				numMotifs++;
			}
			curMotifLength = 0;
			if (numCols > 2) {
				started=1;
				strcpy(consensus,&(cols[0][1]));
				strcpy(name,cols[1]);
				sscanf(cols[2],"%lf",&threshold);
				logp = 0.0;
				if (numCols > 3) sscanf(cols[3],"%lf",&logp);
			} else {
				started=0;
			}
		} else if (numCols > 3) {
			for (int i=0;i<DNA_ALPHA;i++) {
				sscanf(cols[i],"%lf",&(matrix[curMotifLength][i]));
			}
			curMotifLength++;
		}
	}
	if (started) {
		motifs[numMotifs] = new Motif(name,consensus,threshold,
									logp,matrix,curMotifLength);
		numMotifs++;
	}
	
	fclose(fp);

	delete []buf;
	delete []cols;
	delete []consensus;
	delete []name;
 
	fprintf(stderr, "\t%d motifs loaded\n", numMotifs);

	//for (int i=0;i<numMotifs;i++) {
	//motifs[i]->print(stdout);
	//}
	//exit(0);


	return motifs;
}

void Motif::optimizeThreshold(EnrichmentScorer* scorer, double maxBackgroundPercent) {
	if (scorer == NULL) return;
	if (seqStats == NULL) return;
	qsort(seqStats, numSeqs, sizeof(SequenceStats),&(SequenceStats::scoreCmp)); 

	logp = 1.0;	
	double curNumTargets = 0;
	double curNumBackground = 0;
	double maxBackground =  scorer->numBackground*maxBackgroundPercent;
	for (int i=0;i<numSeqs-1;i++) {
		double curThresh = seqStats[i].score;
		double weight = seqStats[i].seq->w;
		if (seqStats[i].seq->g == 1) {
			curNumTargets+=weight;
		} else {
			curNumBackground+=weight;
		}
		if (curThresh < seqStats[i+1].score + 0.000001) {
			continue;
		}
		double curPvalue = scorer->scoreEnrichment(curNumTargets, curNumBackground);
		//fprintf(stderr, "%d\t%lf\t%lf\t%lf\t%lf\n", i,curThresh,curNumTargets,curNumBackground,curPvalue);
		if (curPvalue < logp) {
			logp = curPvalue;
			threshold = curThresh;
			numTarget = curNumTargets;
			numBackground = curNumBackground;
		}
		if (curNumBackground > maxBackground) break;
	}
	threshold -= 0.000001;
	freqTarget = numTarget/scorer->numTargets;
	freqBackground = numBackground/scorer->numBackground;
}

double Motif::scoreOverlap(Motif* m1, Motif* m2, double &percent1, double &percent2) {

	percent1 = 0;
	percent2 = 0;
	if (m1->numSeqs != m2->numSeqs) {
		fprintf(stderr, "!!! Different number of sequences for each motif !!!\n");
		exit(0);
	}
	qsort(m1->seqStats, m1->numSeqs, sizeof(SequenceStats),&(SequenceStats::defaultOrder)); 
	qsort(m2->seqStats, m2->numSeqs, sizeof(SequenceStats),&(SequenceStats::defaultOrder)); 

	int totalBp = 0;
	int totalSites1 = 0;
	int totalSites2 = 0;
	int overlap = 0;
	double expectedOverlap=0;

	int distance = 0;
	int maxMlen = m1->length;
	if (m2->length > maxMlen) maxMlen = m2->length;
	distance = (int) (((double)maxMlen)*MOTIF_OVERLAP_PERCENT);
	

	
	for (int i=0;i<m1->numSeqs;i++) {
		int seqlen = m1->seqStats[i].seq->length;
		int nsites1 = m1->seqStats[i].numSites;
		int nsites2 = m2->seqStats[i].numSites;
		totalBp += seqlen;
	
		int totalBp1 = nsites1*distance;
		int totalBp2 = nsites2*distance;
		int minSites = nsites1;
		if (nsites2 < minSites) minSites = nsites2;
		totalSites1 += nsites1;
		totalSites2 += nsites2;

		if (nsites1 == 0 || nsites2 == 0) continue;

		expectedOverlap += ((double)(totalBp1*totalBp2))/((double)seqlen*seqlen)*((double)(minSites));



		if (1) { //strand == MOTIF_STRAND_BOTH) {
			for (int j=0;j<nsites1;j++) {
				char d1 = m1->seqStats[i].strands[j];
				int p1 = m1->seqStats[i].pos[j];
				if (d1 == MOTIF_STRAND_POS) p1 += m1->length/2;
				else p1 -= m1->length/2;
				for (int k=0;k<nsites2;k++) {
					char d2 = m2->seqStats[i].strands[k];
					int p2 = m2->seqStats[i].pos[k];
					if (d2 == MOTIF_STRAND_NEG) p1 += m2->length/2;
					else p2 -= m2->length/2;

					if (p2-p1 > -distance && p2-p1 < distance) {
						overlap++;
					}
				}
			}
		} else { //if (strand != MOTIF_STRAND_BOTH) {
			//Strand specific code
			// The logic here is that sites are already ordered per sequence
			// first come + strand, in order i.e. 12+ 230+ 300+ then 309- 239- 34- 
			int startIndex2 = 0;
			for (int j=0;j<nsites1;j++) {
				char d1 = m1->seqStats[i].strands[j];
				int p1 = m1->seqStats[i].pos[j];
				for (int k=startIndex2;k<nsites2;k++) {
					char d2 = m2->seqStats[i].strands[k];
					int p2 = m2->seqStats[i].pos[k];
					if (d1 == MOTIF_STRAND_POS) {
						if (d2 == MOTIF_STRAND_NEG) {
							break;
						} else if (p2-p1 < -distance) {
							startIndex2 = k;
						} else if (p2-p1 < distance) {
							overlap++;
						} else {
							break;
						}
					} else { // d1 == MOTIF_STRAND_NEG
						if (d2 == MOTIF_STRAND_POS) {
							startIndex2 = k;
						} else if (p2-p1 > distance) {
							startIndex2 = k;
						} else if (p2-p1 > -distance) {
							overlap++;
						} else {
							break;
						}
					}
				}
			}
		}
	}
	int minNumSites = totalSites1;
	if (minNumSites > totalSites2) minNumSites = totalSites2;

	//fprintf(stderr, "%s = %d %s = %d Overlap = %d (%.2lf%%) (%.2lf%%) Expected = %lf\n",m1->name, totalSites1,m2->name, totalSites2, overlap, 100.0*overlap/totalSites1, 100.0*overlap/totalSites2,expectedOverlap);
	if (minNumSites == 0) return 0;

	percent1 = ((double)overlap)/(double)totalSites1;
	percent2 = ((double)overlap)/(double)totalSites2;

	double score = ((double)overlap)/(double)minNumSites;
	return score;
}

SequenceStats::SequenceStats() {
	initialize();
}
void SequenceStats::initialize() {
	seq = NULL;
	score = -1e100;
	numSites = 0;
	pos = NULL;
	strands = NULL;
}
SequenceStats::~SequenceStats() {
	if (pos != NULL) delete []pos;
	if (strands != NULL) delete []strands;
}
void SequenceStats::setSites(int* npos, char* nstrands, int nsites) {
	if (pos != NULL) {
		delete []pos;
		pos = NULL;
	}
	if (strands != NULL) {
		delete []strands;
		strands = NULL;
	}
	numSites = nsites;
	if (numSites < 1) return;
	pos = new int[numSites];
	strands = new char[numSites];
	for (int i=0;i<numSites;i++) {
		pos[i] = npos[i];
		strands[i] = nstrands[i];
	}
}

int SequenceStats::defaultOrder(const void *a, const void *b) {
	double va=((SequenceStats*)a)->seq->id;
	double vb=((SequenceStats*)b)->seq->id;
	if (va < vb) return 1;
	if (va > vb) return -1;
	return 0;
}
int SequenceStats::scoreCmp(const void *a, const void *b) {
	double va=((SequenceStats*)a)->score;
	double vb=((SequenceStats*)b)->score;
	if (va < vb) return 1;
	if (va > vb) return -1;
	return 0;
}



// ----------------------- class Oligo ------------------------------

int Oligo::cpu = 0;

Oligo::Oligo() {
	init();
}
void Oligo::init() {
	seq = NULL;
	numTarget = 0.0;
	numBackground = 0.0;
	value = 0.0;
	flag = 0;
	flag2 = 0;
	cpuData = NULL;
	normWeight = 1;
}
Oligo::~Oligo() {
	if (seq != NULL) {
		delete []seq;
	}
	if (cpuData != NULL) {
		delete []cpuData;
	}
}
Oligo::Oligo(char* s) {
	init();
	seq = new char[strlen(s)+1];
	strcpy(seq,s);
}
void Oligo::print(FILE* fp) {
	fprintf(fp, "%s\t%f\t%f\t%f\t%f\t%lf\t%d\t%d\n", seq, numTarget,numBackground,logp,normWeight,
				value,flag,flag2);
}
int Oligo::oligoCmp(const void *a, const void *b) {
	double va=(*((Oligo**)a))->value;
	double vb=(*((Oligo**)b))->value;
	if (va < vb) return -1;
	if (va > vb) return 1;
	return 0;
}
int Oligo::oligoCmpReverse(const void *a, const void *b) {
	double va=(*((Oligo**)a))->value;
	double vb=(*((Oligo**)b))->value;
	if (va < vb) return 1;
	if (va > vb) return -1;
	return 0;
}
int Oligo::oligoCmpReverseParallel(const void *a, const void *b) {
	double va=(*((Oligo**)a))->cpuData[cpu];
	double vb=(*((Oligo**)b))->cpuData[cpu];
	if (va < vb) return 1;
	if (va > vb) return -1;
	return 0;
}

void Oligo::revopp(char* seq) {
    if (seq == NULL) return;
    int len = strlen(seq);
    int mid = (int)((((double)len)/2.0)+0.5);
    for (int i=0;i<mid;i++) {
        int rIndex = len-1-i;
        char ogHere = seq[i];
        char ogThere = seq[rIndex];
        char there;
        char here;

        if (ogHere == 'A' || ogHere == 'a') {
            there = 'T';
        } else if (ogHere == 'C' || ogHere == 'c') {
            there = 'G';
        } else if (ogHere == 'G' || ogHere == 'g') {
            there = 'C';
        } else if (ogHere == 'T' || ogHere == 't') {
            there = 'A';
        } else {
            there = 'N';
        }
        if (ogThere == 'A' || ogThere == 'a') {
            here = 'T';
        } else if (ogThere == 'C' || ogThere == 'c') {
            here = 'G';
        } else if (ogThere == 'G' || ogThere == 'g') {
            here = 'C';
        } else if (ogThere == 'T' || ogThere == 't') {
            here = 'A';
        } else {
            here = 'N';
        }
		seq[i] = here;
        seq[rIndex] = there;
    }
}


//==================================================================================
Sequence::Sequence() {
	init();
}
Sequence::Sequence(char* n, char* seq, int group, double weight) {
	init();
	setSeq(n,seq,group,weight);
}
Sequence::~Sequence() {
	if (s != NULL) delete []s;
	if (name != NULL) delete []name;
}
void Sequence::init() {
	s=NULL;
	name=NULL;
	length = 0;
	w=1.0;
	g=0;
	id=0;
}
void Sequence::setSeq(char* n, char* seq, int group, double weight) {
	if (name != NULL) delete []name;
	if (n == NULL) {
		name = NULL;
	} else {
		name = new char[strlen(n)+1];
		strcpy(name,n);
	}
	w = weight;
	g = group;
	if (s != NULL) delete []s;
	if (seq == NULL) {
		s = NULL;
		length = 0;
	} else {
		length = strlen(seq);
		s = new char[length+1];
		strcpy(s,seq);
	}
}

int Sequence::motifScan(Motif* motif, int strand,int offset, int mode,double &bestScore,
										int* sitePositions, char* siteStrands,
										FILE* fp, pthread_mutex_t fileLock) {

	if (offset == MOTIF_OFFSET_NULL) {
		offset = length/-2;
	}

	bestScore = -1e10;
	int numSites = 0;
	if (strand == MOTIF_STRAND_POS || strand==MOTIF_STRAND_BOTH) {
		int totalLength = length-(motif->length)+1;
		for (int i=0;i<totalLength;i++) {
			double score = 0.0; 
			int bad = 0;
			for (int j=0;j<motif->length;j++) {
				switch(s[i+j]) {
					case 'A':
						score += motif->m[j][0];
						break;
					case 'C':
						score += motif->m[j][1];
						break;
					case 'G':
						score += motif->m[j][2];
						break;
					case 'T':
						score += motif->m[j][3];
						break;
					default:
						bad = 1;
						break;
				}
				if (bad) break;
			}
			if (bad) continue;
			if (score > bestScore) bestScore = score;
			if (score >= motif->threshold) {
				sitePositions[numSites] = i+offset;
				siteStrands[numSites] = MOTIF_STRAND_POS;
				numSites++;
				if (mode==KNOWN_JUST_SCORE_ZOOPS && fp == NULL) {
					break;
				}
				if (mode == KNOWN_MASK_MOTIFS) {
					for (int j=0;j<motif->length;j++) s[i+j] = 'N';
				} else if (mode == KNOWN_MSCORE_MOTIFS) {
				} else if (fp != NULL) {
					char tmp = s[i+motif->length];
					s[i+motif->length] = '\0';
					pthread_mutex_lock(&fileLock);
					fprintf(fp,"%s\t%d\t%s\t%s\t+\t%.6lf\n",name,i+offset,&(s[i]),motif->name,score);
					pthread_mutex_unlock(&fileLock);
					s[i+motif->length] = tmp;
				}
			}
		}
	}
	if (mode == KNOWN_JUST_SCORE_ZOOPS && numSites > 0) return 1;

	if (strand == MOTIF_STRAND_BOTH || strand==MOTIF_STRAND_NEG) {
		for (int i=length-1;i>=motif->length-1;i--) {
			double score = 0.0; 
			int bad = 0;
			for (int j=0;j<motif->length;j++) {
				switch(s[i-j]) {
					case 'A':
						score += motif->m[j][3];
						break;
					case 'C':
						score += motif->m[j][2];
						break;
					case 'G':
						score += motif->m[j][1];
						break;
					case 'T':
						score += motif->m[j][0];
						break;
					default:
						bad = 1;
						break;
				}
				if (bad) break;
			}
			if (bad) continue;
			if (score > bestScore) bestScore = score;
			if (score >= motif->threshold) {
				sitePositions[numSites] = i+offset;
				siteStrands[numSites] = MOTIF_STRAND_NEG;
				numSites++;
				if (mode==KNOWN_JUST_SCORE_ZOOPS && fp == NULL) {
					break;
				}
				if (mode == KNOWN_MASK_MOTIFS) {
					for (int j=0;j<motif->length;j++) s[i-j] = 'N';
				} else if (mode == KNOWN_MSCORE_MOTIFS) {
				} else if (fp != NULL) {
					char tmp = s[i+1];
					s[i+1] = '\0';
					pthread_mutex_lock(&fileLock);
					fprintf(fp,"%s\t%d\t%s\t%s\t-\t%.6lf\n",name,i+offset,&(s[i-motif->length+1]),motif->name,score);
					pthread_mutex_unlock(&fileLock);
					s[i+1] = tmp;
				}
			}
		}
	}
	return numSites;
}





//=================--------------------------------
LinkedSequence::LinkedSequence(Sequence* seq, LinkedSequence* parent) {
	s = seq;
	link = parent;
}
LinkedSequence::~LinkedSequence() {
	//if (link != NULL) delete link;
}
SequenceArray::SequenceArray() {
	seqs = NULL;
	numSeqs= 0;
	numTarget=0.0;
	numBackground=0.0;
	bpTarget=0;
	bpBackground=0;
	link=NULL;
	longestSequence = 0;
	scorer = NULL;
	maxBackgroundPercent = MOTIF_MAX_BACKGROUND_PERCENT;
	maxAutoNormIters = MAX_OLIGO_NORMALIZATION_ITERATIONS;
}
SequenceArray::~SequenceArray() {
	//if (link != NULL) delete link;
	deleteLinks();
	for (int i=0;i<numSeqs;i++) {
		delete seqs[i];
	}
	delete []seqs;
	if (scorer != NULL) delete scorer;
}
void SequenceArray::deleteLinks() {
	LinkedSequence* cur = link;
	while (cur != NULL) {
		LinkedSequence* next = cur->link;
		delete cur;
		cur  = next;
	}
}
void SequenceArray::initializeScorer(int scoringMethod,int zoopsFlag,unsigned long long int maxCacheSize) {
	scorer = new EnrichmentScorer();
	scorer->initializeEnrichmentScoring(scoringMethod,numTarget, numBackground,
															zoopsFlag,maxCacheSize);
}
void SequenceArray::addSequence(Sequence* s) {
	link = new LinkedSequence(s, link);
	numSeqs++;
	//fprintf(stderr, "%s\t%d\n", s->name, numSeqs);
}
void SequenceArray::optimizeArray() {
	seqs = new Sequence*[numSeqs];
	//fprintf(stderr, "numSeqs(optimize)=%d\n", numSeqs);
	numSeqs=0;
	numBackground = 0.0;
	numTarget = 0.0;
	bpTarget = 0;
	bpBackground = 0;
	longestSequence = 0;
	LinkedSequence* ls = link;
	while (ls != NULL) {
		seqs[numSeqs] = ls->s;
		seqs[numSeqs]->id = numSeqs;
		if (seqs[numSeqs]->length > longestSequence) {
			longestSequence = seqs[numSeqs]->length;
		}
		//fprintf(stderr, "New long: %d, %s %d\n", seqs[numSeqs]->length, seqs[numSeqs]->s,longestSequence);

		if (seqs[numSeqs]->g == 0) {
			numBackground += seqs[numSeqs]->w;
			bpBackground += seqs[numSeqs]->length;
		} else {
			numTarget += seqs[numSeqs]->w;
			bpTarget += seqs[numSeqs]->length;
		}
		ls = ls->link;
		numSeqs++;
	}
	longestSequence *= 2;
	//fprintf(stderr, "numSeqs(optimize)=%d\n", numSeqs);
}

void SequenceArray::parseFasta2SeqAndGroupFiles(char* inputFASTAfile,
							char* bgFASTAfile, char* seqFile, char* groupFile) {
	
	FILE* fpseq = fopen(seqFile,"w");
	if (fpseq == NULL) {
		fprintf(stderr, "!!! Could not open a temporary file (%s) for writing!!!\n", seqFile);
		exit(0);
	}
	FILE* fpgroup = fopen(groupFile,"w");
	if (fpgroup == NULL) {
		fprintf(stderr, "!!! Could not open a temporary file (%s) for writing!!!\n", groupFile);
		exit(0);
	}

	Inttable* names = new Inttable(1000000);	

	char* buf = new char[MOTIF2_BUFFER];

	char* curSeq = new char[1000000];
	char* curName = new char[MOTIF2_BUFFER];

	for (int i=0;i<2;i++) {

		if (i==0 && inputFASTAfile == NULL) continue;
		if (i==1 && bgFASTAfile == NULL) continue;
	
		FILE* fp = NULL;
		int g = 1;
		if (i==0) {		
			fp = fopen(inputFASTAfile,"r");
			if (fp == NULL) {
				fprintf(stderr, "Could not open FASTA file %s !!!\n", inputFASTAfile);
				exit(0);
			}
		} else {
			fp = fopen(bgFASTAfile,"r");
			if (fp == NULL) {
				fprintf(stderr, "Could not open FASTA file %s !!!\n", bgFASTAfile);
				exit(0);
			}
			g = 0;
		}
		curSeq[0] = '\0';
		curName[0] = '\0';
		while (fgets(buf, MOTIF2_BUFFER, fp) != NULL) {
			int lineLength = strlen(buf);
			if (lineLength > 0) buf[lineLength-1]='\0';
			//for carriage returns
			if (lineLength > 1 && buf[lineLength-2] == 13) buf[lineLength-2]='\0';
			

			if (buf[0] == '>') {
				if (curSeq[0] != '\0' && curName[0] != '\0') {
					fprintf(fpgroup,"%s\t%d\t1.0\n",curName,g);
					fprintf(fpseq,"%s\t%s\n",curName,curSeq);
				}
				curName[0] = '\0';
				curSeq[0] = '\0';
				if (lineLength > 1) {
					strcpy(curName,&(buf[1]));
					int c = names->search(curName);
					if (c != EMPTY_INT) {
						names->insert(++c,curName);
						sprintf(curName,"%s-%d",&(buf[1]),c);
					} else {
						names->insert(1,curName);
					}
				}
			} else {
				strcat(curSeq,buf);
			}
		}
		if (curSeq[0] != '\0' && curName[0] != '\0') {
			fprintf(fpgroup,"%s\t%d\t1.0\n",curName,g);
			fprintf(fpseq,"%s\t%s\n",curName,curSeq);
		}
		fclose(fp);

	}

	fclose(fpseq);
	fclose(fpgroup);

	delete []curSeq;
	delete []curName;
	delete names;
	delete []buf;
}


void SequenceArray::readSequenceFiles(char* seqFile, char* groupFile) {

	Hashtable* groupdata = new Hashtable();

	fprintf(stderr, "\tReading input files...\n");

	char* buf = new char[MOTIF2_BUFFER];
	char** cols = new char*[1000];
	int numCols;
	char* seq = NULL;
	char* name = NULL;


	FILE* fp = NULL;
	if (groupFile != NULL) {	
		fp = fopen(groupFile, "r");
		if (fp == NULL) {
			fprintf(stderr, "!!! Could not open group file: %s !!!\n", groupFile);
			exit(0);
		}
		while (fgets(buf, MOTIF2_BUFFER, fp) != NULL) {
			splitMotif2(buf,cols,numCols,'\t');
			if (numCols < 2) continue;
			GroupInfo* gi = new GroupInfo();
			gi->weight = 1.0;
			gi->group = 0;
			name = cols[0];
			if (numCols > 1) {
				sscanf(cols[1],"%d",&(gi->group));
			}
			if (numCols > 2) {
				sscanf(cols[2],"%lf",&(gi->weight));
			}
			groupdata->insert(gi, name);
		}
		fclose(fp);
	}

	fp = fopen(seqFile, "r");
	if (fp == NULL) {
		fprintf(stderr, "!!! Could not open seq file: %s !!!\n", seqFile);
		exit(0);
	}

	int totalSequences = 0;
	while (fgets(buf, MOTIF2_BUFFER, fp) != NULL) {
		splitMotif2(buf,cols,numCols,'\t');
		if (numCols < 2) continue;
		int group = 1;
		double weight = 1.0;
		if (groupFile != NULL) {
			GroupInfo* gi = (GroupInfo*)groupdata->search(cols[0]);
			if (gi == NULL) continue;
			group = gi->group;
			weight = gi->weight;
		}
		name = cols[0];
		seq = cols[1];
		cleanSequence(seq);
		//if (strlen(seq) < 1) continue;
		totalSequences++;
		Sequence* s = new Sequence(name,seq,group,weight);
		//fprintf(stderr, "%s\t%s\t%d\n", name,seq,s->length);
		//fprintf(stderr, "add seqeuene.\n");
		addSequence(s);
	}
	fprintf(stderr, "\t%d total sequences read\n", totalSequences);

	optimizeArray();
	/*
	double oligoPerSeqTarget = totalTargetOligos/numTargets;
	double oligoPerSeqBackround = totalBackgroundOligos/numBackground;
	double newNumBackground = totalBackgroundOligos/oligoPerSeqTarget;

	fprintf(stderr, "\tAutoadjustment for sequence coverage in background: %.2lfx\n",newNumBackground/numBackground);
	numBackground = newNumBackground;
	*/

	char** keys = groupdata->keys();
	for (int i=0;i<groupdata->total;i++) {
		GroupInfo* gi = (GroupInfo*)groupdata->search(keys[i]);
		delete gi;
		delete [](keys[i]);
	}
	delete []keys;
	delete groupdata;
	delete []buf;
	delete []cols;

}
void SequenceArray::printGroupFile(FILE* fp) {
	if (fp == NULL) {
		fprintf(stderr, "!!! Cannot output normalized group file !!!\n");
		return;
	}
	for (int i=0;i<numSeqs;i++) {
		fprintf(fp, "%s\t%d\t%lf\n",seqs[i]->name,seqs[i]->g,seqs[i]->w);
	}
}
void SequenceArray::printSequenceFile(FILE* fp) {
	if (fp == NULL) {
		fprintf(stderr, "!!! Cannot output sequence file !!!\n");
		return;
	}
	for (int i=0;i<numSeqs;i++) {
		fprintf(fp, "%s\t%s\n",seqs[i]->name,seqs[i]->s);
	}
}

void SequenceArray::maskMotifFromOligoArray(Motif* motif, OligoArray* oligos, char strand) {
	motif->logxform();
	motif->initializeSeqStats(numSeqs);
	double score = 0;
	if (longestSequence < 1) longestSequence = 2;
	int *sitePositions = new int[longestSequence];
	char *siteStrands = new char[longestSequence];
	int numPositions2Mask = 0;
	int *positions2Mask = new int[longestSequence*motif->length*2];
	char *mask = new char[longestSequence];
	for (int i=0;i<longestSequence;i++) {
		mask[i] = 0;
	}
	int mlen = motif->length;
	int offset = 0;
	double totalTsites = 0;
	double totalBsites = 0;
	pthread_mutex_t notused;
	pthread_mutex_init(&notused,NULL);
	pthread_mutex_unlock(&notused);
	for (int i=0;i<numSeqs;i++) {

		int nsites = seqs[i]->motifScan(motif,strand,offset,KNOWN_FIND_ALL,score,
								sitePositions,siteStrands,NULL,notused);
		motif->seqStats[i].seq = seqs[i];
		motif->seqStats[i].score = score;
		motif->seqStats[i].setSites(sitePositions, siteStrands,nsites);

		numPositions2Mask = 0;
		int slen = seqs[i]->length;
		double targetValue = 0;
		double backgroundValue = 0;
		if (seqs[i]->g) {
			targetValue = -1*seqs[i]->w;
			totalTsites += seqs[i]->w*nsites;
		} else {
			backgroundValue = -1*seqs[i]->w;
			totalBsites += seqs[i]->w*nsites;
		}

		for (int j=0;j<nsites;j++) {
			int p = sitePositions[j];
			//fprintf(stderr, "%s site %d of %d. %d, strand=%d\n", seqs[i]->name, j, nsites,p,siteStrands[j]);
			if (siteStrands[j] == MOTIF_STRAND_POS) {
				for (int k=-mlen+1;k<mlen;k++) {
					int pos = p+k;
					if (pos < 0) continue;
					if (pos > slen-mlen) break;
					if (mask[pos] != 0) continue; // don't double delete
					char* s = &(seqs[i]->s[pos]);
					char tmp = seqs[i]->s[pos+mlen];
					seqs[i]->s[pos+mlen] = '\0';
					//fprintf(stderr, "%d\t%s\n", k,s);
					int ok = checkSequence(s);
					if (ok) {
						oligos->adjustOligoInstances(s, targetValue, backgroundValue);
					}
					seqs[i]->s[pos+mlen] = tmp;

					mask[pos] = 1;
					positions2Mask[numPositions2Mask++]=pos;
				}
			} else {
				for (int k=-mlen+1;k<mlen;k++) {
					int pos = p-mlen+1+k;
					if (pos < 0) continue;
					if (pos > slen-mlen) break;
					if (mask[pos] != 0) continue; // don't double delete

					char* s = &(seqs[i]->s[pos]);
					char tmp = seqs[i]->s[pos+mlen];
					seqs[i]->s[pos+mlen] = '\0';
					//fprintf(stderr, "%d\t%s\n", k,s);
					int ok = checkSequence(s);
					if (ok) {
						oligos->adjustOligoInstances(s, targetValue, backgroundValue);
					}
					seqs[i]->s[pos+mlen] = tmp;
				
					mask[pos] = 1;
					positions2Mask[numPositions2Mask++]=pos;
				}
			}
		}
		//now literally mask the sequence
		for (int j=0;j<nsites;j++) {
			int p = sitePositions[j];
			if (siteStrands[j] == MOTIF_STRAND_POS) {
				for (int k=0;k<mlen;k++) seqs[i]->s[p+k] = 'N';
			} else {
				for (int k=0;k<mlen;k++) seqs[i]->s[p-k] = 'N';
			}
		}
		for (int j=0;j<numPositions2Mask;j++) {
			int p  = positions2Mask[j];
			mask[p] = 0;
		}
	}
	motif->calculateDistributionStats();

	delete []mask;
	delete []sitePositions;
	delete []siteStrands;
	delete []positions2Mask;
}
void SequenceArray::optimizeMotifThreshold(Motif* m, char strand) {
	m->logxform();
	m->numTarget = 0.0;
	m->numBackground = 0.0;
	m->initializeSeqStats(numSeqs);
	int *sitePositions = new int[longestSequence];
	char *siteStrands = new char[longestSequence];
	int mode = KNOWN_OPTIMIZE_THREASHOLD;
	int offset = 0;
	pthread_mutex_t notused;
	pthread_mutex_init(&notused,NULL);
	pthread_mutex_unlock(&notused);
	for (int i=0;i<numSeqs;i++) {
		double score = 0;
		//int nsites = 0;
		(void)seqs[i]->motifScan(m,strand,offset,mode,score,sitePositions,
											siteStrands,NULL,notused);
		if (mode == KNOWN_OPTIMIZE_THREASHOLD || mode==KNOWN_FIND_ALL) {
			m->seqStats[i].seq = seqs[i];
			m->seqStats[i].score = score;
		}
	}
	m->optimizeThreshold(scorer,maxBackgroundPercent);
	delete []sitePositions;
	delete []siteStrands;
}

void SequenceArray::motifScan(Motif* m, int strand,int offset, int mode,FILE* fp, pthread_mutex_t fileLock) {
	m->logxform();
	m->numTarget = 0.0;
	m->numBackground = 0.0;
	int *sitePositions = new int[longestSequence];
	char *siteStrands = new char[longestSequence];
	if (mode == KNOWN_OPTIMIZE_THREASHOLD || mode==KNOWN_FIND_ALL) {
		m->initializeSeqStats(numSeqs);
	}
	//fprintf(stderr, "numSeqs = %d, offset = %d\n", numSeqs,offset);
	for (int i=0;i<numSeqs;i++) {
		double score = 0;
		int nsites = 0;
		nsites = seqs[i]->motifScan(m,strand,offset,mode,score,sitePositions,
													siteStrands,fp,fileLock);
		//fprintf(stderr, "%d\n", nsites);
		if (nsites > 0) {
			if (seqs[i]->g == 1) {
				m->numTarget += seqs[i]->w;
			} else {
				m->numBackground += seqs[i]->w;
			}
		}
		if (mode == KNOWN_OPTIMIZE_THREASHOLD || mode==KNOWN_FIND_ALL) {
			m->seqStats[i].seq = seqs[i];
			m->seqStats[i].score = score;
			m->seqStats[i].setSites(sitePositions, siteStrands,nsites);
		}
		if (mode == KNOWN_MSCORE_MOTIFS) {
			pthread_mutex_lock(&fileLock);
			fprintf(fp,"%s\t%d\t%s\t%s\t+\t%.6lf\n",seqs[i]->name,0,m->consensus,m->name,score);
			pthread_mutex_unlock(&fileLock);
		}
	}
	if (numTarget > 0.0) m->freqTarget = m->numTarget / numTarget;
	if (numBackground > 0.0) m->freqBackground = m->numBackground / numBackground;

	if (mode == KNOWN_OPTIMIZE_THREASHOLD || mode==KNOWN_FIND_ALL) {
		m->calculateDistributionStats();
		if (mode == KNOWN_OPTIMIZE_THREASHOLD) {
			//m->deleteSites();
		}
	}

	delete []sitePositions;
	delete []siteStrands;
}
void Motif::calculateDistributionStats() {
	long double avgTarget= 0.0;
	double nTarget = 0;
	long double avgBack= 0.0;
	double nBack = 0;
	double nPos = 0.0;
	double nNeg = 0.0;
	double totalSites = 0.0;
	double totalSeqWithSites = 0.0;

	if (seqStats == NULL || numSeqs < 1) return;

	for (int i=0;i<numSeqs;i++) {
		double w = seqStats[i].seq->w;
		int g = seqStats[i].seq->g;
		if (seqStats[i].numSites < 1) continue;

		double uniqSites = 0;

		for (int j=0;j<seqStats[i].numSites;j++) {
			char d = seqStats[i].strands[j];
			int p = seqStats[i].pos[j];
			if (d == MOTIF_STRAND_POS) {
				p += length/2;
			} else {
				p -= length/2;
			}
			if (g == 1) {
				avgTarget += w*p;
				nTarget += w;
				if (d == MOTIF_STRAND_POS) nPos += w;
				else nNeg += w;
				int good = 1;
				for (int k=j+1;k<seqStats[i].numSites;k++) {
					char d2 = seqStats[i].strands[k];
					int p2 = seqStats[i].pos[k];
					if (d2 == MOTIF_STRAND_POS) {
						p2 += length/2;
					} else {
						p2 -= length/2;
					}
					if (p-p2>length/-2 && p-p2<length/2) {
						good = 0;
						break;
					}
				}
				if (good) uniqSites+=w;
			} else {
				avgBack += w*p;
				nBack += w;
			}
		}
		if (g == 1) {
			totalSeqWithSites += w;
			totalSites += uniqSites;
		}

	}

	if (nTarget > 0.0) avgTarget /= nTarget;
	if (nBack > 0.0) avgBack /= nBack;
	double stdTarget=0.0;
	double stdBack = 0.0;
	nTarget = 0.0;
	nBack = 0.0;

	for (int i=0;i<numSeqs;i++) {
		double w = seqStats[i].seq->w;
		int g = seqStats[i].seq->g;
		if (seqStats[i].numSites < 1) continue;

		for (int j=0;j<seqStats[i].numSites;j++) {
			char d = seqStats[i].strands[j];
			int p = seqStats[i].pos[j];
			if (d == MOTIF_STRAND_POS) {
				p += length/2;
			} else {
				p -= length/2;
			}
			if (g == 1) {
				double v = w*(p-avgTarget);
				stdTarget += v*v;
				nTarget += w;
			} else {
				double v = w*(p-avgBack);
				stdBack += v*v;
				nBack += w;
			}
		}
	}
	if (nTarget > 0.0) stdTarget /= nTarget;
	if (nBack > 0.0) stdBack /= nBack;
	stdTarget = sqrt(stdTarget);
	stdBack = sqrt(stdBack);

	statAvgTarget = avgTarget;
	statStdTarget = stdTarget;
	statAvgBackground = avgBack;
	statStdBackground = stdBack;

	statStrandLogRatio = 0.0;
	if (nNeg < 0.00001 && nPos > 0.00001) statStrandLogRatio = 10.0;
	else if (nNeg > 0.00001 && nPos < 0.00001) statStrandLogRatio = -10.0;
	else if (nNeg < 0.00001 && nPos < 0.00001) statStrandLogRatio = 0.0;
	else statStrandLogRatio = log(nPos/nNeg)/log(2);

	siteMultiplicity = 1.0;
	if (totalSeqWithSites > 0.01) siteMultiplicity = totalSites/totalSeqWithSites;

}


void SequenceArray::normalizeSequence(int nlen,int strand,char* normFile,int neutralFlag) {
	if (nlen < 1) return;

	double lastError = 1e100;
	double next = 1.0;
	double cur = 0.0;
	
	for (int i=0;i<maxAutoNormIters;i++) {
		int printFlag = 0;
		if (i==0) {
			printFlag = 1;
		}
		double curError = normalizeSequenceIteration(nlen,strand,normFile,printFlag,neutralFlag);
		//double improvementRate = (lastError-curError)/lastError;
		//fprintf(stderr, "Round %d:\tLastError: %le\tCurrentError:%le\tImprovment:%le\n", i+1,lastError,curError,improvementRate);
		if (curError >= lastError) break;
		//if (improvementRate < MIN_OLIGO_NORMALIZATION_IMPROVEMENT) break;
		lastError = curError;
		if (printFlag) {
			fprintf(stderr, "\tConverging on autonormalization solution:\n\t");
		} else {
			cur += 80.0/(double)maxAutoNormIters;
			if (cur > next) {
				fprintf(stderr, ".");
				next += 1.0;
			}
		}
	}
	fprintf(stderr, "\n");
	fprintf(stderr, "\tFinal normalization:");
	(void)normalizeSequenceIteration(nlen,strand,NULL,1,neutralFlag);
}

double SequenceArray::normalizeSequenceIteration(int nlen,int strand,char* normFile,int printFlag,int neutralFlag) {

	FILE* fp = NULL;
	if (normFile != NULL && printFlag) {
		fp = fopen(normFile, "w");
		if (fp == NULL) {
			fprintf(stderr, "!!! Could not open %s for writing!!!\n",normFile);
		} else {
			fprintf(fp, "Oligo\tTargetCounts\tBackgroundCounts\tNormalizationFactor\n");
		}
	}
	char* alphaIndex = new char[256];
	alphaIndex[65]=0; //alphaIndex['A']=0;
	alphaIndex[67]=1; //alphaIndex['C']=1;
	alphaIndex[71]=2; //alphaIndex['G']=2;
	alphaIndex[84]=3; //alphaIndex['T']=3;
	char* alpha = new char[256];
	alpha[0]='A';
	alpha[1]='C';
	alpha[2]='G';
	alpha[3]='T';
	int* totalGood = new int[numSeqs];
	double* neutralFreq = new double[nlen+1];
	neutralFreq[0] = 1.0;
	for (int i=0;i<nlen+1;i++) {
		neutralFreq[i]=pow(4.0,-1.0*i);
	}

	double error = 0.0;

	char* mer = new char[1000];
	for (int curLen=1;curLen<=nlen;curLen++) {
		int total = (int) pow(4.0,(double)curLen);
		if (printFlag) fprintf(stderr, "\tAutonormalization: %d-mers (%d total)\n",curLen,total);
		double* targetLevels = new double[total];
		double* backgroundLevels = new double[total];
		double* normFactors = new double[total];
		double* normFactorsT = new double[total];
		int* pow4 = new int[curLen+1];
		pow4[0]=1;
		for (int i=1;i<=curLen;i++) pow4[i] = 4*pow4[i-1];
		for (int i=0;i<total;i++) {
			targetLevels[i]=0.0;
			backgroundLevels[i]=0.0;
			normFactors[i] = 1.0;
			normFactorsT[i] = 1.0;
		}

		for (int i=0;i<numSeqs;i++) {
			char* s = seqs[i]->s;
			int goodOligos = 0;
			int nbins = seqs[i]->length-curLen;
			for (int j=0;j<nbins;j++) {
				int bad = 0;
				for (int k=0;k<curLen;k++) {
					if (s[j+k] == 'N') {
						bad = 1;
						break;
					}
				}
				if (bad) continue;
				goodOligos++;
			}
			totalGood[i] = goodOligos;
			double gOligos = (double)totalGood[i];
			for (int j=0;j<nbins;j++) {
				int index = 0;
				int bad = 0;
				for (int k=0;k<curLen;k++) {
					if (s[j+k] == 'N') {
						bad = 1;
						break;
					}
					index += alphaIndex[(int)s[j+k]]*pow4[k];
				}
				if (bad) continue;
				if (seqs[i]->g == 1) {
					targetLevels[index] += seqs[i]->w/gOligos;
				} else {
					backgroundLevels[index] += seqs[i]->w/gOligos;
				}
			}
		}

		double totalTarget = 0.0;
		double totalBackground = 0.0;
		for (int i=0;i<total;i++) {
			totalTarget += targetLevels[i];
			totalBackground += backgroundLevels[i];
		}
		double minimumTargetLevels = 0.5/totalTarget;
		double minimumBackgroundLevels = 0.5/totalBackground;
		for (int i=0;i<total;i++) {
			for (int j=0;j<curLen;j++) {
				int remainder = i % pow4[j+1];
				int asdf = (int)(remainder/pow4[j]);
				mer[j] = alpha[asdf];
			}
			mer[curLen]='\0';

			double tLevel = targetLevels[i];
			double bLevel = backgroundLevels[i];
			if (strand == MOTIF_STRAND_BOTH) {
				int rindex = 0;
				for (int j=curLen-1;j>=0;j--) {
					rindex += (3-alphaIndex[(int)mer[j]])*pow4[curLen-j-1];
				}
				tLevel = (tLevel+targetLevels[rindex])/2.0;
				bLevel = (bLevel+backgroundLevels[rindex])/2.0;
			}

			if (bLevel < minimumBackgroundLevels) bLevel = minimumBackgroundLevels;
			if (tLevel < minimumTargetLevels) tLevel = minimumTargetLevels;
			
			if (neutralFlag == 0) {
				double factor = totalBackground*tLevel/totalTarget/bLevel;
				if (printFlag) {
					if (curLen < 3) {
						fprintf(stderr, "\t\t%s\t%.2lf%%\t%.2lf%%\t%.3lf\n", mer,
											100.0*tLevel/totalTarget,
											100.0*bLevel/totalBackground,factor);
					}
					if (fp != NULL) {
						//fprintf(fp, "%s\t%.1lf\t%.1lf\t%.3lf\n", mer,tLevel,bLevel,factor);
						fprintf(fp, "%s\t%.2lf%%\t%.2lf%%\t%.3lf\n", mer, 
											100.0*tLevel/totalTarget,
											100.0*bLevel/totalBackground,factor);
					}
				}
				normFactors[i] = factor;
				error += (factor-1.0)*(factor-1.0)/(double)total;
			} else {
				double factorNeutralT = totalTarget/tLevel*neutralFreq[curLen];
				double factorNeutralB = totalBackground/bLevel*neutralFreq[curLen];
				if (printFlag) {
					if (curLen < 3) {
						fprintf(stderr, "\t\t%s\t%.2lf%%\t%.2lf%%\t%.3lf,%.3lf\n", mer,
											100.0*tLevel/totalTarget,
											100.0*bLevel/totalBackground,factorNeutralT,factorNeutralB);
					}
					if (fp != NULL) {
						//fprintf(fp, "%s\t%.1lf\t%.1lf\t%.3lf\n", mer,tLevel,bLevel,factor);
						fprintf(fp, "%s\t%.2lf%%\t%.2lf%%\t%.3lf,%.3lf\n", mer, 
											100.0*tLevel/totalTarget,
											100.0*bLevel/totalBackground,factorNeutralT,factorNeutralB);
					}
				}
				normFactorsT[i] = factorNeutralT;
				normFactors[i] = factorNeutralB;
				error += (factorNeutralT-1.0)*(factorNeutralT-1.0)/(double)total;
				error += (factorNeutralB-1.0)*(factorNeutralB-1.0)/(double)total;
			}
		}

		totalBackground = 0.0;
		totalTarget = 0.0;
		for (int i=0;i<numSeqs;i++) {
			char* s = seqs[i]->s;
			double scoreUnit = (double)seqs[i]->w;
			if (neutralFlag == 0 && seqs[i]->g == 1) {
				totalTarget += scoreUnit;
				continue;
			}
			double newScore = 0;
			int nbins = seqs[i]->length-curLen;
			double gOligos = (double)totalGood[i];
			for (int j=0;j<nbins;j++) {
				int index = 0;
				int bad = 0;
				for (int k=0;k<curLen;k++) {
					if (s[j+k] == 'N') {
						bad=1;
						break;
					}
					index += alphaIndex[(int)s[j+k]]*pow4[k];
				}
				if (bad) continue;
				if (seqs[i]->g == 1) {
					newScore += normFactorsT[index];
				} else {
					newScore += normFactors[index];
				}
			}
			if (gOligos > 0.5) {
				newScore /= gOligos;
			}
			double newWeight = newScore*scoreUnit;
			if (newWeight < HOMER_MINIMUM_SEQ_WEIGHT) newWeight = HOMER_MINIMUM_SEQ_WEIGHT;
			if (newWeight > 1.0/HOMER_MINIMUM_SEQ_WEIGHT) newWeight = 1.0/HOMER_MINIMUM_SEQ_WEIGHT;

			double penalty = newWeight;
			if (newWeight < 1) penalty = 1/newWeight;
			penalty*=penalty;
		
			double delta = newWeight-seqs[i]->w;
			double newWeight1 = seqs[i]->w+delta;
			//double newWeight2 = seqs[i]->w+delta;
			if (penalty > 1 && ((delta > 0 && newWeight > 1) || (delta < 0 && newWeight < 1))) {
				newWeight1 = seqs[i]->w+delta/penalty;
			} else {
			}
			//fprintf(stdout, "1\t%lf\t%lf\n",newWeight2,seqs[i]->w);
			//fprintf(stdout, "%lf\n",penalty);
			//fprintf(stdout, "%lf\t%lf\n",newWeight2,seqs[i]->w);

			//if (newWeight1 > 4) newWeight1 = 4;
			//if (newWeight1 < 0.25) newWeight1 = 0.25;
			
			seqs[i]->w = newWeight1;
			if (seqs[i]->g == 1) {
				totalTarget += seqs[i]->w;
			} else {
				totalBackground += seqs[i]->w;
			}
		}
		numTarget = totalTarget;
		numBackground = totalBackground;
		delete []targetLevels;
		delete []backgroundLevels;
		delete []normFactors;
		delete []pow4;
	}
	if (fp != NULL) {
		fprintf(stderr, "\tNormalization weights can be found in file: %s\n", normFile);
		fclose(fp);
	}
	delete []mer;
	delete []alphaIndex;
	delete []alpha;
	delete []totalGood;
	delete []neutralFreq;

	return error;

}

void cleanSequence(char* seq) {
	if (seq == NULL) return;
	int len = strlen(seq);
	for (int i=0;i<len;i++) {
		switch(seq[i]) {
			case 'A':
			case 'C':
			case 'G':
			case 'T':
			case 'N':
				break;
			case 'a':
				seq[i] = 'A';
				break;
			case 'c':
				seq[i] = 'C';
				break;
			case 'g':
				seq[i] = 'G';
				break;
			case 't':
				seq[i] = 'T';
				break;
			default:
				seq[i] = 'N';
		}
	}
}

//=====================================================================
KnownMotifs::KnownMotifs() {
	seqArray = new SequenceArray();
	totalBp = 0;
	normLength = 0;
	alphaSize = DNA_ALPHA;
	strand = MOTIF_STRAND_BOTH;
	maxCacheSize = MAX_STAT_CACHE_SIZE*1000000/sizeof(double);
	offset = MOTIF_OFFSET_NULL;

	numTargets = 0.0;
	numBackground= 0.0;

	motifs = NULL;
	numMotifs = 0;
	maxAutoNormIters = MAX_OLIGO_NORMALIZATION_ITERATIONS;

	numCPUs = 1;
	nextMeter = 0.0;
	incMeter = 0.0;
	
	seqFile = NULL;
	groupFile = NULL;
	inputFASTAfile = NULL;
	bgFASTAfile = NULL;
	fastaFlag = 0;
	knownFile = NULL;
	normFile = NULL;
	sitefp = NULL;
	normfp = NULL;
	motifFilename = NULL;
	neutralFlag = 0;

	maxBackgroundPercent = MOTIF_MAX_BACKGROUND_PERCENT;
	siteReduceThreshold = -1.0;

	//scoringMethod = MOTIF_SCORING_HYPERGEOMETRIC;
	scoringMethod = MOTIF_SCORING_BINOMIAL;
	scorer = NULL;

	mscoreFlag = 0;
	threasholdOptimizationFlag = KNOWN_JUST_SCORE_ZOOPS;
	zoopsFlag = SCORER_STANDARD_ZOOPS;
	
	knownFlag = 1;	
	findFlag = 0;
	maskFlag = 0;
	normOnlyFlag = 0;
}

KnownMotifs::~KnownMotifs() {
	delete seqArray;
	delete scorer;
}
void KnownMotifs::checkParameters() {
	if (numCPUs < 1) numCPUs = 1;
	if (seqFile == NULL){
		if (inputFASTAfile != NULL) {
			fastaFlag=1;
		} else {
			fprintf(stderr, "!!! Need a sequence file for homer to do anything !!!\n");
			exit(0);
		}
	}
	if (fastaFlag) {
		if (findFlag == 0 && maskFlag == 0 && bgFASTAfile == NULL){
			fprintf(stderr, "!!! Need a background FASTA file for homer to do anything !!!\n");
			exit(0);
		}
	} else {
		if (findFlag == 0 && maskFlag == 0 && groupFile == NULL){
			fprintf(stderr, "!!! Need a group file for homer to do anything !!!\n");
			exit(0);
		}
	}
	if (normOnlyFlag == 0 && knownFile == NULL ){
		fprintf(stderr, "!!! Need a file of known motifs for homer to do anything !!!\n");
		exit(0);
	}
}
void KnownMotifs::parseFasta() {
	SequenceArray* seqarray = new SequenceArray();
	seqFile = new char[10000];
	groupFile = new char[10000];
	srand(time(NULL));
	int r = rand();
	sprintf(seqFile, "0%d.seq", r);
	sprintf(groupFile, "0%d.group", r);
	seqarray->parseFasta2SeqAndGroupFiles(inputFASTAfile,bgFASTAfile,seqFile,groupFile);
	delete seqarray;
}


void KnownMotifs::loadMotifs() {
	motifs = Motif::readMotifFile(knownFile,numMotifs);
}

void KnownMotifs::loadSequence() {

	seqArray->maxAutoNormIters = maxAutoNormIters;
	seqArray->readSequenceFiles(seqFile,groupFile);
	if (normLength > 0 && !findFlag) {
		seqArray->normalizeSequence(normLength, strand, normFile,neutralFlag);
	}

	numTargets = seqArray->numTarget;
	numBackground = seqArray->numBackground;
	if (normOnlyFlag) {
		seqArray->printGroupFile(normfp);
	}
}


void KnownMotifs::getMotifEnrichment() {

	scorer = new EnrichmentScorer();
	if (!findFlag) {
		scorer->initializeEnrichmentScoring(scoringMethod,numTargets, numBackground,
												zoopsFlag,maxCacheSize);
	}

	if (maskFlag) {
		threasholdOptimizationFlag = KNOWN_MASK_MOTIFS;
	} else if (mscoreFlag) {
		threasholdOptimizationFlag = KNOWN_MSCORE_MOTIFS;
	} else if (siteReduceThreshold > 0) {
		threasholdOptimizationFlag = KNOWN_FIND_ALL;
	} else if (findFlag) {
		threasholdOptimizationFlag = KNOWN_FIND_ALL;
		fprintf(stderr, "\tFinding instances of %d motif(s)\n", numMotifs);
	} else {
		fprintf(stderr, "\tChecking enrichment of %d motif(s)\n", numMotifs);
	}

	fprintf(stderr, "\t|0%%                                    50%%                                  100%%|\n");
	fprintf(stderr, "\t");

	//parallel check
	pthread_mutex_init(&mutex, NULL);
	pthread_mutex_init(&fileLock, NULL);
	mutexIndex = 0;
	mutexTotal = numMotifs;
	incMeter = 1.0/((double)numMotifs);
	nextMeter = 0;

	ThreadArgs_scanMotifThread** args = new ThreadArgs_scanMotifThread*[numCPUs];

	pthread_t* thread = new pthread_t[numCPUs];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	void* status;
	for (unsigned int i=0;i<numCPUs;i++) {
		args[i] = new ThreadArgs_scanMotifThread(this,i);
		(void)pthread_create(&(thread[i]), &attr, 
						&KnownMotifs_scanMotifThread, args[i]);
	}
	pthread_attr_destroy(&attr);
	for (unsigned int i=0;i<numCPUs;i++) {
		(void)pthread_join(thread[i], &status);
		delete args[i];
	}
	delete []thread;
	delete []args;
	pthread_mutex_destroy(&mutex);
	pthread_mutex_destroy(&fileLock);

	fprintf(stderr, "\n");
}


void* KnownMotifs_scanMotifThread(void* threadArgs) {
	ThreadArgs_scanMotifThread* args = (ThreadArgs_scanMotifThread*)threadArgs;
	args->known->scanMotifThread(args);
	pthread_exit(threadArgs);
	return NULL;
}
ThreadArgs_scanMotifThread::ThreadArgs_scanMotifThread(KnownMotifs* km, int ncpu) {
	known = km;
	cpu = ncpu;
}

void KnownMotifs::scanMotifThread(ThreadArgs_scanMotifThread* args) {
	int i=-1;
	while (i<mutexTotal) {
		pthread_mutex_lock(&mutex);
		i=mutexIndex++;
		if (i < mutexTotal) {
			while (((double)i+1.0)/((double)numMotifs) >= nextMeter) {
				fprintf(stderr, "=");
				nextMeter += 1.0/80.0;
			}
		}
		pthread_mutex_unlock(&mutex);
		if (i >= mutexTotal) break;


		seqArray->motifScan(motifs[i],strand,offset,threasholdOptimizationFlag,sitefp,fileLock);
		if (zoopsFlag == SCORER_APPROXIMATE_ZOOPS) {
			motifs[i]->freqTarget = scorer->getZoopsApproxTarget(motifs[i]->numTarget)/seqArray->numTarget;
			motifs[i]->freqBackground = scorer->getZoopsApproxBackground(motifs[i]->numBackground)
																			/seqArray->numBackground;
		}
		if (threasholdOptimizationFlag == KNOWN_OPTIMIZE_THREASHOLD) {
			motifs[i]->optimizeThreshold(scorer,maxBackgroundPercent);
		} else if (threasholdOptimizationFlag == KNOWN_JUST_SCORE_ZOOPS) {
			if (findFlag == 0 && numBackground > 0.5 && numTargets > 0.5) {
				double logp = scorer->scoreEnrichment(motifs[i]->numTarget, motifs[i]->numBackground);
				motifs[i]->logp = logp;
			} else {
				motifs[i]->logp = 0.0;
			}
		}
	}

}



void KnownMotifs::printMotifEnrichment(FILE* fp) {
	if (fp == NULL) return;
	if (numMotifs < 1) {
		fprintf(stderr, "!!! No motifs to analyze !!!\n");
		return;
	}
	qsort(motifs, numMotifs, sizeof(Motif*),&(Motif::motifCmp));

	fprintf(fp, "Motif Name\tConsensus\tp-value\tLog p-value\tLog2 enrichment ratio");
	fprintf(fp, "\t# target sequences with Motif(total=%d)", (int) numTargets);
	fprintf(fp, "\t%% of target sequences with Motif");
	fprintf(fp, "\t# background sequences with Motif(total=%d)", (int) numBackground);
	fprintf(fp, "\t%% of background sequences with Motif");
	fprintf(fp, "\n");
	for (int i=0;i<numMotifs;i++) {
		fprintf(fp,"%s\t%s\t1e%d\t%.3lf",motifs[i]->name,motifs[i]->consensus,
					(int)(motifs[i]->logp/2.3026),motifs[i]->logp);
		double tfreq = motifs[i]->numTarget/numTargets;
		double bgfreq = motifs[i]->numBackground/numBackground;
		double ratio = tfreq;
		if (ratio < 1.0/numTargets) {
			ratio = 1.0/numTargets;
		}
		if (bgfreq < 1.0/numBackground) {
			ratio = ratio / (1.0/numBackground);
		} else {
			ratio = ratio / bgfreq;
		}
		ratio = log(ratio)/log(2.0);
		fprintf(fp, "\t%.3lf",ratio);
		fprintf(fp, "\t%d\t%.2lf%%",(int)motifs[i]->numTarget,tfreq*100);
		fprintf(fp, "\t%d\t%.2lf%%",(int)motifs[i]->numBackground,bgfreq*100);
		fprintf(fp, "\n");
	
	}
}
void KnownMotifs::printMotifs(FILE* fp) {
	if (numMotifs < 1) return;
	qsort(motifs, numMotifs, sizeof(Motif*), &(Motif::motifCmp));
	for (int i=0;i<numMotifs;i++) {
		motifs[i]->print(fp);
	}
}

void KnownMotifs::findOverlappingMotifs() {

	qsort(motifs, numMotifs, sizeof(Motif*), &(Motif::motifCmp));
	
	char* mask = new char[numMotifs];
	for (int i=0;i<numMotifs;i++) {
		mask[i] = 0;
	}
	for (int i=0;i<numMotifs;i++) {
		if (mask[i]) continue;
		for (int j=i+1;j<numMotifs;j++) {
			if (mask[j]) continue;
			double percentI = 0.0, percentJ=0.0;
			double score = Motif::scoreOverlap(motifs[i],motifs[j],percentI, percentJ);
			if (score > siteReduceThreshold) {
				mask[j] = 1;
				motifs[j]->similarPercent = percentJ;
				motifs[j]->similarPercentOfOther = percentI;
				motifs[j]->similarMotif = motifs[i];
			}
		}
	}
}

//=================================================================================

EnrichmentScorer::EnrichmentScorer() {
	numTargets = 0;
	numBackground = 0;
	zoopsFlag = 1;
	targetZoopsApprox = NULL;
	bgZoopsApprox = NULL;
	maxZoopsApprox = 0;
	enrichmentFunction = NULL;

	cache = NULL;
	maxCacheSize = MAX_STAT_CACHE_SIZE;
}
EnrichmentScorer::~EnrichmentScorer() {
	if (targetZoopsApprox != NULL) delete []targetZoopsApprox;
	if (bgZoopsApprox != NULL) delete []bgZoopsApprox;
	if (cache != NULL) {
		for (int i=0;i<cacheLength;i++) {
			if (cache[i] != NULL) delete [](cache[i]);
		}
		delete []cache;
	}
}

double EnrichmentScorer::scoreEnrichmentHypergeometric(double pv, double nv) {
	unsigned int targetValue = (unsigned int)pv;
	unsigned int bgValue = (unsigned int)nv;
	if (zoopsFlag == SCORER_APPROXIMATE_ZOOPS) {
		if (targetValue >= maxZoopsApprox) targetValue = maxZoopsApprox-1;
		if (bgValue >= maxZoopsApprox) bgValue = maxZoopsApprox-1;
		targetValue = (unsigned int)targetZoopsApprox[targetValue];
		bgValue = (unsigned int)bgZoopsApprox[bgValue];
	}
	unsigned int mtotal = targetValue+bgValue;
	double logp = 0.0;

	if (mtotal < (unsigned int)cacheLength) {
		if (cache[mtotal] == NULL) {
			cache[mtotal]=new double[mtotal+1];
			for (unsigned int i=0;i<mtotal+1;i++) cache[mtotal][i]=10.0;
		}
		if (cache[mtotal][targetValue] > 1.0) {
			cache[mtotal][targetValue] = loghypergeoD((unsigned int)(numTargets+numBackground),
					(unsigned int)numTargets, mtotal, targetValue);
		}
		logp = cache[mtotal][targetValue];
	} else {
		logp = loghypergeoD((unsigned int)(numTargets+numBackground),
				(unsigned int)numTargets, mtotal, targetValue);
	}
	return logp;
}


double EnrichmentScorer::scoreEnrichmentZscore(double pv, double nv) {
	double targetValue = pv;
	double bgValue = nv;
	if (zoopsFlag == SCORER_APPROXIMATE_ZOOPS) {
		unsigned int pvInt = (unsigned int)pv;
		unsigned int nvInt = (unsigned int)nv;
		if (pvInt >= maxZoopsApprox) pvInt = maxZoopsApprox-1;
		if (nvInt >= maxZoopsApprox) nvInt = maxZoopsApprox-1;
		targetValue = targetZoopsApprox[pvInt];
		bgValue = bgZoopsApprox[nvInt];
	}
	if (bgValue < 1.0) bgValue = 1.0;

	double bgRatio = bgValue/numBackground;
	double std = sqrt(numBackground*bgRatio*(1-bgRatio));
	if (std < 1.0) std = 1.0;
	double expected = numTargets*bgRatio;
	return -1* (targetValue-expected)/std;

}
double EnrichmentScorer::scoreEnrichmentBinomial(double pv, double nv) {
	unsigned int targetValue = (unsigned int)pv;
	unsigned int bgValue = (unsigned int)nv;
	if (zoopsFlag == SCORER_APPROXIMATE_ZOOPS) {
		if (targetValue >= maxZoopsApprox) targetValue = maxZoopsApprox-1;
		if (bgValue >= maxZoopsApprox) bgValue = maxZoopsApprox-1;
		targetValue = (unsigned int)targetZoopsApprox[targetValue];
		bgValue = (unsigned int)bgZoopsApprox[bgValue];
	}
	unsigned mtotal = targetValue+bgValue;
	double ratio = ((double)bgValue)/numBackground;
	double logp = 0.0;
	if (mtotal < (unsigned int)cacheLength) {
		if (cache[mtotal] == NULL) {
			cache[mtotal]=new double[mtotal+1];
			for (unsigned int i=0;i<mtotal+1;i++) cache[mtotal][i]=10.0;
		}
		if (cache[mtotal][targetValue] > 1.0) {
			cache[mtotal][targetValue] = logbinomialD((unsigned int)(numTargets),
						(unsigned int)targetValue, ratio, (unsigned int)numBackground);
		} else {
			//fprintf(stderr, ".");
		}
		logp = cache[mtotal][targetValue];
	} else {
		logp = logbinomialD((unsigned int)(numTargets),
						(unsigned int)targetValue, ratio, (unsigned int)numBackground);
	}
	return logp;
}
void EnrichmentScorer::initializeEnrichmentScoring(int scoringMethod,double numT, 
									double numB,int zFlag, unsigned long long int newMaxCacheSize) {

	numTargets = (unsigned int)numT;
	numBackground = (unsigned int)numB;
	zoopsFlag = zFlag;
	maxCacheSize = newMaxCacheSize;
	cacheLength = (int)sqrt((double)maxCacheSize*2.0);
	fprintf(stderr, "\tCache length = %d\n", cacheLength);
	cache = new double*[cacheLength];
	for (int i=0;i<cacheLength;i++) {
		cache[i] = NULL;
	}

	if (scoringMethod == MOTIF_SCORING_HYPERGEOMETRIC) {
		enrichmentFunction = &EnrichmentScorer::scoreEnrichmentHypergeometric;
		fprintf(stderr, "\tUsing hypergeometric scoring\n");
	} else if (scoringMethod == MOTIF_SCORING_BINOMIAL) {
		enrichmentFunction = &EnrichmentScorer::scoreEnrichmentBinomial;
		fprintf(stderr, "\tUsing binomial scoring\n");
	} else if (scoringMethod == MOTIF_SCORING_ZSCORE) {
		enrichmentFunction = &EnrichmentScorer::scoreEnrichmentZscore;
		fprintf(stderr, "\tUsing z-score scoring\n");
	}

	if (zoopsFlag == SCORER_APPROXIMATE_ZOOPS) {
		if (targetZoopsApprox != NULL) delete []targetZoopsApprox;
		if (bgZoopsApprox != NULL) delete []bgZoopsApprox;
			
		maxZoopsApprox = (int) (5.0*(numT+numB));

		//fprintf(stderr, "MaxZoopsApprox=%d\n", maxZoopsApprox);
		targetZoopsApprox = new double[maxZoopsApprox];
		bgZoopsApprox = new double[maxZoopsApprox];
		targetZoopsApprox[0]=0.0;
		bgZoopsApprox[0]=0.0;
		targetZoopsApprox[1]=1.0;
		bgZoopsApprox[1]=1.0;
		for (unsigned int i=2;i<maxZoopsApprox;i++) {
			targetZoopsApprox[i]=targetZoopsApprox[i-1]+(numT-targetZoopsApprox[i-1])/numT;
			bgZoopsApprox[i]=bgZoopsApprox[i-1]+(numB-bgZoopsApprox[i-1])/numB;
		}
	}
}

double EnrichmentScorer::scoreEnrichment(double pv, double nv) {
	
	return (this->*enrichmentFunction)(pv, nv);
}
double EnrichmentScorer::getZoopsApproxTarget(double v) {
	if (v >= maxZoopsApprox) v = maxZoopsApprox-1;
	return targetZoopsApprox[(int)v];
}
double EnrichmentScorer::getZoopsApproxBackground(double v) {
	if (v >= maxZoopsApprox) v = maxZoopsApprox-1;
	return bgZoopsApprox[(int)v];
}





