
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

#include "SeqTag.h"

void printCMD();

int main(int argc, char** argv) {

	int updateFlag = 0;
	char* directory = NULL;
	char** files = NULL;
	char** tagDirs = NULL;
	char** tagFiles = NULL;
	char* genome = NULL;
	char* genomeDirectory = NULL;
	char* name = NULL;
	char* normGCfile = NULL;
	char* normOligofile = NULL;
	char* restrictionSite = NULL;
	int numDirs = 0;
	int maxMisMatches = -1;
	int flipFlag = 0;
	double minmapq = 10.0;
	int rsmaxMisMatches = 0;
	int restrictionMode = 4;
	int midpointFlag = 0;
	int numFiles = 0;
	int sspeFlag = 0;
	int numTagFiles = 0;
	//int maskFlag = 0;
	int cytosineContext = CYTOSINE_CONTEXT_CG;
	int format = FORMAT_UNKNOWN;
	int mode = MODE_UNIQUE;
	char* cmd = new char[10000];
	int autoCorrelateFlag = 1;
	int tagCountDistributionFlag = 1;
	int tagLengthDistributionFlag = 1;
	int checkGCFlag = 0;
	//int chrOnlyFlag = 0;
	int singleFileFlag = 0;
	int freqStart = -50;
	int freqEnd = 50;
	int oligoStart = 100;
	int oligoEnd = -100;
	int restrictionSiteLength = -10000000;
	int peBackgroundLength = -10000000;
	int oligoLength = 0;
	int pairedEndFlag = 0;
	double maxNormRatio = 2.0;
	double minNormRatio = 0.25;
	double gcWindow = 0.025;
	double tbp = 0.0;
	int autoCorrRange = 4000;
	double autoCorrMaxTags = 1e9;
	float maxNormTBP = 5.0;
	int rmPETagBg = 0;
	int peStyleFlag = 0;
	int removeSpikeSize = -1;
	double removeSpikeFold = -1.0;
	int removeSelfLigationFlag = 0;
	int removeRestrictionEnds = 0;
	double parseAlignmentCpGMinValue = 10.0;

	int peLocalWindowSize = 20000;
	int peLargeWindowSize = 300000000;
	int peLargeResolution = 1000;
	int minReadLength = INT_MIN;
	int maxReadLength = INT_MAX;

	double maxPerror = 10.0;

	char defaultPENames[] = "petag";
	
	int fragLength = FRAGMENT_LEN_AUTO;

	HomerConfig*  hc = new HomerConfig();

	if (argc < 2) {
		printCMD();
	}
	strcpy(cmd,argv[0]);
	for (int i=1;i<argc;i++) {
		strcat(cmd," ");
		strcat(cmd,argv[i]);
	}
	for (int i=1;i<argc;i++) {
		
		if (i==1) {
			if (argv[i][0] == '-') {
				fprintf(stderr, "\n!!!!!!!!!!!!\n\tNEED to specify directory with first argument!!!\n");
				printCMD();
			}
			directory = argv[i];
			continue;
		}
		if (argv[i][0] == '-') {
			if (strcmp(argv[i],"-keep")==0) {
				mode = MODE_KEEPONE;
			} else if (strcmp(argv[i],"-keepOne")==0) {
				mode = MODE_KEEPONE;
			} else if (strcmp(argv[i],"-keepAll")==0) {
				mode = MODE_KEEPALL;
			} else if (strcmp(argv[i],"-unique")==0) {
				mode = MODE_UNIQUE;
			} else if (strcmp(argv[i],"-mask")==0) {
				//maskFlag = 1;
			} else if (strcmp(argv[i],"-update")==0) {
				updateFlag =1;
			} else if (strcmp(argv[i],"-flip")==0) {
				flipFlag =1;
			} else if (strcmp(argv[i],"-C")==0) {
				format = FORMAT_BOWTIE_COLOR;
			} else if (strcmp(argv[i],"-forceBED")==0) {
				mode = MODE_UNIQUE;
			} else if (strcmp(argv[i],"-assignMidPoint")==0) {
				midpointFlag = 1;;
			} else if (strcmp(argv[i],"-force5th")==0) {
				mode = MODE_BED_FORCE5TH;
			} else if (strcmp(argv[i],"-single")==0) {
				singleFileFlag = 1;
			} else if (strcmp(argv[i],"-chrOnly")==0) {
				//chrOnlyFlag = 1;
			} else if (strcmp(argv[i],"-illuminaPE")==0) {
				peStyleFlag = SEQFILE_FORMAT_ILLUMINAPE;
			} else if (strcmp(argv[i],"-genome")==0) {
				genome = argv[++i];
				int slen = strlen(genome);
				if (genome[slen-1] == 'r') {
					genome[slen-1]='\0';
					//maskFlag = 1;
				}
				genomeDirectory = hc->getGenomeDirectory(genome);
				if (genomeDirectory == NULL) exit(0);
			} else if (strcmp(argv[i],"-checkGC")==0) {
				checkGCFlag = 1;
			} else if (strcmp(argv[i],"-normGC")==0) {
				normGCfile = argv[++i];
				checkGCFlag = 1;
				Tag::precision = 3;
				PETag::precision = 3;
				fprintf(stderr, "\tTag value precision set to 3 for normalized tag data\n");
			} else if (strcmp(argv[i],"-normOligo")==0) {
				sscanf(argv[++i], "%d", &oligoLength);
				if (oligoLength < 1 || oligoLength > 20) {
					fprintf(stderr, " !!! Need to enter a valid oligo length: -normOligo <#> !!!\n");
					exit(0);
				}
				Tag::precision = 3;
				PETag::precision = 3;
				fprintf(stderr, "\tTag value precision set to 3 for normalized tag data\n");
			} else if (strcmp(argv[i],"-mCcontext")==0) {
				i++;
				if (strcmp(argv[i],"CG")) {
					cytosineContext = CYTOSINE_CONTEXT_CG;
				} else if (strcmp(argv[i],"CHG")) {
					cytosineContext = CYTOSINE_CONTEXT_CHG;
				} else if (strcmp(argv[i],"CHH")) {
					cytosineContext = CYTOSINE_CONTEXT_CHH;
				} else if (strcmp(argv[i],"all")) {
					cytosineContext = CYTOSINE_CONTEXT_ALL;
				} else {
					fprintf(stderr, "!!! %s not recognized for -mCcontext option!!!\n", argv[i]);
					exit(0);
				}
			} else if (strcmp(argv[i],"-oligoStart")==0) {
				sscanf(argv[++i], "%d", &oligoStart);
			} else if (strcmp(argv[i],"-oligoEnd")==0) {
				sscanf(argv[++i], "%d", &oligoEnd);
			} else if (strcmp(argv[i],"-restrictionSiteLength")==0) {
				sscanf(argv[++i], "%d", &restrictionSiteLength);
			} else if (strcmp(argv[i],"-minCounts")==0) {
				sscanf(argv[++i], "%lf", &parseAlignmentCpGMinValue);
			} else if (strcmp(argv[i],"-PEbgLength")==0) {
				sscanf(argv[++i], "%d", &peBackgroundLength);
			} else if (strcmp(argv[i],"-normFixedOligo")==0) {
				normOligofile = argv[++i];
			} else if (strcmp(argv[i],"-removePEbg")==0) {
				rmPETagBg = 1;
			} else if (strcmp(argv[i],"-removeSpikes")==0) {
				sscanf(argv[++i], "%d", &removeSpikeSize);
				sscanf(argv[++i], "%lf", &removeSpikeFold);
			} else if (strcmp(argv[i],"-mis")==0) {
				sscanf(argv[++i], "%d", &maxMisMatches);
			} else if (strcmp(argv[i],"-sspe")==0) {
				sspeFlag = 1;
			} else if (strcmp(argv[i],"-mapq")==0) {
				sscanf(argv[++i], "%lf", &minmapq);
			} else if (strcmp(argv[i],"-rsmis")==0) {
				sscanf(argv[++i], "%d", &rsmaxMisMatches);
			} else if (strcmp(argv[i],"-freqStart")==0) {
				sscanf(argv[++i], "%d", &freqStart);
			} else if (strcmp(argv[i],"-freqEnd")==0) {
				sscanf(argv[++i], "%d", &freqEnd);
			} else if (strcmp(argv[i],"-minlen")==0) {
				sscanf(argv[++i], "%d", &minReadLength);
			} else if (strcmp(argv[i],"-maxlen")==0) {
				sscanf(argv[++i], "%d", &maxReadLength);
			} else if (strcmp(argv[i],"-removeSelfLigation")==0) {
				removeSelfLigationFlag = 1;
			} else if (strcmp(argv[i],"-removeRestrictionEnds")==0) {
				removeRestrictionEnds = 1;
			} else if (strcmp(argv[i],"-iterNorm")==0) {
				sscanf(argv[++i], "%lf", &maxPerror);
				minNormRatio = 0.0;
				maxNormRatio = 1.0;
			} else if (strcmp(argv[i],"-minNormRatio")==0) {
				sscanf(argv[++i], "%lf", &minNormRatio);
			} else if (strcmp(argv[i],"-maxNormRatio")==0) {
				sscanf(argv[++i], "%lf", &maxNormRatio);
			} else if (strcmp(argv[i],"-restrictionSite")==0) {
				restrictionSite = argv[++i];
				fprintf(stderr, "\tRestriction site set to %s\n",restrictionSite);
			} else if (strcmp(argv[i],"-both")==0) {
				restrictionMode = 0;
				fprintf(stderr, "\tWill only keep reads where both ends are near restriction sites\n");
			} else if (strcmp(argv[i],"-one")==0) {
				restrictionMode = 3;
				fprintf(stderr, "\tWill keep reads with at least one end near a restriction site\n");
			} else if (strcmp(argv[i],"-onlyOne")==0) {
				restrictionMode = 2;
				fprintf(stderr, "\tWill keep reads with at only one end near a restriction site\n");
			} else if (strcmp(argv[i],"-none")==0) {
				restrictionMode = 1;
				fprintf(stderr, "\tWill keep reads where both ends are far from restriction site\n");
			} else if (strcmp(argv[i],"-tbp")==0) {
				sscanf(argv[++i], "%lf", &tbp);
			} else if (strcmp(argv[i],"-name")==0) {
				name = argv[++i];
			} else if (strcmp(argv[i],"-precision")==0) {
				int p = TAG_VALUE_RESOLUTION;
				sscanf(argv[++i],"%d",&p);
				Tag::precision = p;
				PETag::precision = p;
				fprintf(stderr, "\tTag value precision set to %d\n", p);
			} else if (strcmp(argv[i],"-len")==0 || strcmp(argv[i],"-fragLength")==0) {
				if (i+1 >= argc ) {
					fprintf(stderr, "Error specifying -len\n");
					printCMD();
				}
				i++;
				if (strcmp(argv[i],"given") == 0) {
					fragLength = FRAGMENT_LEN_GIVEN;
				} else {
					sscanf(argv[i],"%d", &fragLength);
				}
			} else if (strcmp(argv[i],"-format")==0) {
				if (i+1 >= argc ) {
					fprintf(stderr, "Error specifying -format\n");
					printCMD();
				}
				i++;
				if (strcmp(argv[i],"bowtie") == 0) {
					if (format != FORMAT_BOWTIE_COLOR) format = FORMAT_BOWTIE;
				} else if (strcmp(argv[i],"bed") == 0) {
					format = FORMAT_BED;
				} else if (strcmp(argv[i],"eland_result") == 0) {
					format = FORMAT_ELANDRESULT;
				} else if (strcmp(argv[i],"eland_export") == 0) {
					format = FORMAT_ELANDEXPORT;
				} else if (strcmp(argv[i],"eland_extended") == 0) {
					format = FORMAT_ELANDEXTENDED;
				} else if (strcmp(argv[i],"mCpGbed") == 0 || strcmp(argv[i],"mcpgbed") == 0) {
					format = FORMAT_MCPGBED;
					Tag::precision = 2;
					PETag::precision = 2;
					fprintf(stderr, "\tTag value precision set to 2 for mC data\n");
				} else if (strcmp(argv[i],"allC") == 0 || strcmp(argv[i],"allc") == 0) {
					format = FORMAT_LISTER_ALLC;
					Tag::precision = 2;
					PETag::precision = 2;
					fprintf(stderr, "\tTag value precision set to 2 for mC data\n");
				} else if (strcmp(argv[i],"bismark") == 0) {
					format = FORMAT_BISMARK;
					Tag::precision = 2;
					PETag::precision = 2;
					fprintf(stderr, "\tTag value precision set to 2 for mC data\n");
				} else if (strcmp(argv[i],"sam") == 0) {
					format = FORMAT_SAM;
				} else if (strcmp(argv[i],"HiCsummary") == 0) {
					format = FORMAT_HICSUMMARY;
				} else {
					fprintf(stderr, "Error specifying -format\n");
					printCMD();
				}
			} else if (strcmp(argv[i],"-d")==0) {
				i++;
				for (;i<argc;i++) {
					if (argv[i][0] == '-') {
						i--;
						break;
					}
					char** newdirs = new char*[numDirs+1];
					for (int j=0;j<numDirs;j++) newdirs[j] = tagDirs[j];
					if (tagDirs != NULL) delete []tagDirs;
					tagDirs = newdirs;
					tagDirs[numDirs] = argv[i];
					numDirs++;
					fprintf(stderr, "\tWill add tag directory: %s\n", tagDirs[numDirs-1]);
				}
			} else if (strcmp(argv[i],"-t")==0) {
				i++;
				for (;i<argc;i++) {
					if (argv[i][0] == '-') {
						i--;
						break;
					}
					char** newtagfiles = new char*[numTagFiles+1];
					for (int j=0;j<numTagFiles;j++) newtagfiles[j] = tagFiles[j];
					if (tagFiles != NULL) delete []tagFiles;
					tagFiles = newtagfiles;
					tagFiles[numTagFiles] = argv[i];
					numTagFiles++;
					fprintf(stderr, "\tWill parse tag file: %s\n", tagFiles[numTagFiles-1]);
				}
			} else {
				fprintf(stderr, "!!! Couldn't recognize: %s !!!\n", argv[i]);
				printCMD();
			}
		} else {
			char** newfiles = new char*[numFiles+1];
			for (int j=0;j<numFiles;j++) newfiles[j] = files[j];
			if (files != NULL) delete []files;
			files = newfiles;
			files[numFiles++] = argv[i];
			if (pairedEndFlag == 0) {
				int index = 0;
				while (argv[i][index] != '\0') {
					if (argv[i][index] == ',') {
						pairedEndFlag = 1;
						fprintf(stderr, "\tMaking paired end tag directory\n");
						break;
					}
					index++;
				}
			}
			fprintf(stderr,"\tWill parse file: %s\n", files[numFiles-1]);
		}
	}

	if (normOligofile != NULL && (oligoStart > oligoEnd)) {
		fprintf(stderr, "!!! Need to properly set -oligoStart <#> and -oligoEnd <#> to use -normFixedOligo <...>!!!\n");
		exit(0);
	}
	if (updateFlag == 0 && numFiles == 0 && numDirs == 0 && numTagFiles == 0) {
		fprintf(stderr, "!!! No input files specified!!!\n");
		printCMD();
	}
	if (restrictionSite != NULL && genome == NULL) {
		fprintf(stderr, "!!! Must specify genome to check for restriction sites\n");
		exit(1);
	}
	fprintf(stderr, "\n");

	TagLibrary* tags = new TagLibrary(directory);
	if (format == FORMAT_MCPGBED || format == FORMAT_LISTER_ALLC || format == FORMAT_BISMARK) {
		tags->parseAlignmentCpGMinValue = parseAlignmentCpGMinValue;
		tags->parseAlignmentCcontext = cytosineContext;
		tags->mCflag=1;
	}
	if (updateFlag) {
		fprintf(stderr, "\tUsing existing Tag Directory (%s)\n",directory);
		tags->readTagDirectory();
	}
	tags->maxmismatches = maxMisMatches;
	tags->maxReadLength = maxReadLength;
	tags->minReadLength = minReadLength;
	tags->minmapq = minmapq;
	tags->peStyleFlag = peStyleFlag;
	tags->sspeFlag = sspeFlag;
	tags->flipFlag = flipFlag;
	tags->setFragLength(fragLength);
	tags->setMaxTBP(tbp);
	if (genome != NULL) tags->setGenome(genome);
	if (name != NULL) tags->setName(name);
	if (updateFlag == 0) {
		if (cmd != NULL) tags->setCMD(cmd);
		tags->setSingleFile(singleFileFlag);
		tags->pairedEndFlag = pairedEndFlag;
		tags->parseAlignmentFiles(files,numFiles,format,mode,tagDirs,numDirs,tagFiles,numTagFiles);
	} else {
		if (tbp > 0.0) {
			//need to save tbp modification...
			fprintf(stderr, "\tRestricting tags per bp...\n");
			tags->readAndSave();
		}
	}
	pairedEndFlag = tags->pairedEndFlag;
	//Tag directory stats

	if (tagLengthDistributionFlag == 1) {
		int max = MAX_READ_LENGTH;
		double *d = tags->getTagLengthDistribution(NULL, max);
		if (d != NULL) delete []d;
	}
	if (fragLength != FRAGMENT_LEN_AUTO && fragLength != FRAGMENT_LEN_GIVEN) {
		//fprintf(stderr, "\tForcing fragment length = %d\n", fragLength);
		tags->setFragLength(fragLength);
	}
	if (tagCountDistributionFlag) {
		int max = MAX_TAGS_PER_BP;
		double* d = tags->getTagCountDistribution(NULL, max);
		if (d != NULL) delete []d;
	}

	int peUpdateFlag = 0;
	if (pairedEndFlag) {
		int distLength = 0;
		double *dist = tags->getPETagDistribution(peLocalWindowSize,peLargeWindowSize,
										peLargeResolution,defaultPENames, distLength);
		if (dist != NULL) delete []dist;
	}

	if (autoCorrelateFlag) {
		tags->autoCorrelateTags(NULL, autoCorrRange, autoCorrMaxTags);
	} 
	if (fragLength != FRAGMENT_LEN_AUTO && fragLength != FRAGMENT_LEN_GIVEN) {
		fprintf(stderr, "\tForcing fragment length = %d\n", fragLength);
		tags->setFragLength(fragLength);
	}


	if (genomeDirectory != NULL && (checkGCFlag || oligoEnd >= oligoStart)) {
		OligoArray* oa = NULL;
		fprintf(stderr, "\tChecking GC bias...\n");
		int currEst = tags->fragmentLengthEstimate;
		if (pairedEndFlag) tags->setFragLength(currEst/2); // we want to only check half of the fragment
		NucleotideFreq* nf = tags->checkTagSeqBias(genomeDirectory,freqStart,
										freqEnd,oa, oligoStart, oligoEnd);
									//	freqEnd,oa, 1, 0); //oligoStart, oligoEnd);
		if (nf != NULL && normGCfile != NULL && pairedEndFlag == 0) {
			fprintf(stderr, "\tCorrecting GC bias...\n");
			char* nfile = normGCfile;
			if (strcmp(nfile,"default")==0) {
				nfile = new char[10000];
				sprintf(nfile, "%s/genomeGCcontent.txt",tags->directory);
			}
			tags->normalizeTagCountsGC(nfile, nf, minNormRatio, maxNormRatio, gcWindow,maxPerror);
			if (nfile != normGCfile) {
				delete []nfile;
			}
		}
		if (nf != NULL && normOligofile != NULL) {
			if (strcmp(normOligofile,"default") !=0) {
				int ncpus = 1;
				if (oa != NULL) 
					delete oa;
				oa = new OligoArray();
				oa->readOligoCountFile(normOligofile,MOTIF_STRAND_POS,ncpus);
			}
			tags->normalizeTagCountsFixedOligo(genomeDirectory, oa, oligoStart, oligoEnd, 
														minNormRatio, maxNormRatio);
		}
		if (nf != NULL) delete nf;

		if (pairedEndFlag) tags->setFragLength(currEst); // reset
	}
	if (genomeDirectory != NULL && oligoLength > 0) {
		Hashtable* oligos = new Hashtable();
		tags->normalizeTagCountsOligos(genomeDirectory, oligos, oligoLength,
										oligoStart, oligoEnd,minNormRatio, maxNormRatio,maxNormTBP);
	}

	if (pairedEndFlag && rmPETagBg) {
		if (peBackgroundLength < 100000) {
			peBackgroundLength = (int)(1.5*(double)tags->fragmentLengthEstimate);
		}
		tags->removePETagBackground(peBackgroundLength);
		peUpdateFlag = 1;
	}
	if (pairedEndFlag && restrictionSite != NULL) {
		if (genomeDirectory == NULL) {
			fprintf(stderr, "!!! Must specify valid genome !!!\n");
		} else {
			if (restrictionSiteLength < -1000000) {
				restrictionSiteLength = (int)(1.5*(double)tags->fragmentLengthEstimate);
			}
			tags->assignPETagsToRestrictionSites(restrictionSite,rsmaxMisMatches,genomeDirectory,
						restrictionMode,midpointFlag,removeSelfLigationFlag,removeRestrictionEnds,
						restrictionSiteLength);
			peUpdateFlag = 1;
		}
	}
	if (removeSpikeSize > 0) {
		tags->removeTagSpikes(removeSpikeSize,removeSpikeFold);
		peUpdateFlag = 1;
	}
	if (pairedEndFlag && peUpdateFlag) {
		int numDist = 0;
		double* dist = tags->getPETagDistribution(peLocalWindowSize,peLargeWindowSize,peLargeResolution,
								defaultPENames,numDist);
		if (dist != NULL) {
			delete []dist;
		}
	}

	tags->printTagInfo();

	delete hc;
}


void printCMD() {
	fprintf(stderr, "\n\tUsage: makeTagDirectory <directory> <alignment file 1> [file 2] ... [options]\n");
	fprintf(stderr, "\n\tCreates a platform-independent 'tag directory' for later analysis.\n"); 
	fprintf(stderr, "\tCurrently BED, eland, bowtie, and sam files are accepted. The program will try to\n");
	fprintf(stderr, "\tautomatically detect the alignment format if not specified.  Program will also\n");
	fprintf(stderr, "\tunzip *.gz, *.bz2, and *.zip files and convert *.bam to sam files on the fly\n");
	fprintf(stderr, "\tExisting tag directories can be added or combined to make a new one using -d/-t\n");
	fprintf(stderr, "\tIf more than one format is needed and the program cannot auto-detect it properly,\n");
	fprintf(stderr, "\tmake separate tag directories by running the program separately, then combine them.\n");
	fprintf(stderr, "\tTo perform QC/manipulations on an existing tag directory, add \"-update\"\n");
	fprintf(stderr, "\n\tOptions:\n");
	//fprintf(stderr, "\t\t-name <experiment name> (optional, names the experiment)\n");
	fprintf(stderr, "\t\t-fragLength <# | given> (Set estimated fragment length - given: use read lengths)\n");
	fprintf(stderr, "\t\t\tBy default treats the sample as a single read ChIP-Seq experiment\n");
	fprintf(stderr, "\t\t-format <X> where X can be: (with column specifications underneath)\n");
	fprintf(stderr, "\t\t\tbed - BED format files:\n");
	fprintf(stderr, "\t\t\t\t(1:chr,2:start,3:end,4:+/- or read name,5:# tags,6:+/-)\n");
	fprintf(stderr, "\t\t\t\t-force5th (5th column of BED file contains # of reads mapping to position)\n");
	fprintf(stderr, "\t\t\tsam - SAM formatted files (use samTools to covert BAMs into SAM if you have BAM)\n");
	fprintf(stderr, "\t\t\t\t-unique (keep if there is a single best alignment based on mapq)\n");
	fprintf(stderr, "\t\t\t\t\t-mapq <#> (Minimum mapq for -unique, default: 10, set negative to use AS:i:/XS:i:)\n");
	fprintf(stderr, "\t\t\t\t-keepOne (keep one of the best alignments even if others exist)\n");
	fprintf(stderr, "\t\t\t\t-keepAll (include all alignments in SAM file)\n");
	fprintf(stderr, "\t\t\t\t-mis (Maximum allowed mismatches, default: no limit, uses MD:Z: tag)\n");
	fprintf(stderr, "\t\t\t\t-sspe (strand specific, paired-end reads[flips strand of 2nd read to match])\n");
	fprintf(stderr, "\t\t\tbowtie - output from bowtie (run with --best -k 2 options)\n");
	fprintf(stderr, "\t\t\t\t(1:read name,2:+/-,3:chr,4:position,5:seq,6:quality,7:NA,8:misInfo)\n");
	fprintf(stderr, "\t\t\teland_result - output from basic eland\n");
	fprintf(stderr, "\t\t\t\t(1:read name,2:seq,3:code,4:#zeroMM,5:#oneMM,6:#twoMM,7:chr,\n");
	fprintf(stderr, "\t\t\t\t\t\t\t8:position,9:F/R,10-:mismatches\n");
	fprintf(stderr, "\t\t\teland_export - output from illumina pipeline (22 columns total)\n");
	fprintf(stderr, "\t\t\t\t(1-5:read name info,9:sequence,10:quality,11:chr,13:position,14:strand)\n");
	fprintf(stderr, "\t\t\teland_extended - output from illumina pipeline (4 columns total)\n");
	fprintf(stderr, "\t\t\t\t(1:read name,2:sequence,3:match stats,4:positions[,])\n");
	fprintf(stderr, "\t\t\tmCpGbed - encode style mCpG reporting in extended BED format, no auto-detect\n");
	fprintf(stderr, "\t\t\t\t(1:chr,2:start,3:end,4:name,5:,6:+/-,7:,8:,9:,10:#C,11:#mC)\n");
	fprintf(stderr, "\t\t\tallC - Lister style output files detailing the read information about all cytosines\n");
	fprintf(stderr, "\t\t\t\t(1:chr,2:pos,3:strand,4:context,#mC,#totalC,#unmC\n");
	fprintf(stderr, "\t\t\tbismark - Bismark style output files detailing the read information about all cytosines\n");
	fprintf(stderr, "\t\t\t\t(1:chr,2:pos,3:strand,4:#mC,5:#unmC,6:context,7:triseq\n");
	fprintf(stderr, "\t\t\t\t-minCounts <#> (minimum number of reads to report mC/C ratios, default: 10)\n");
	fprintf(stderr, "\t\t\t\t-mCcontext <CG|CHG|CHH|all> (only use C's in this context, default: CG)\n");
	fprintf(stderr, "\t\t\tHiCsummary - minimal paired-end read mapping information\n");
	fprintf(stderr, "\t\t\t\t(1:readname,2:chr1,3:5'pos1,4:strand1,5:chr2,6:5'pos2,7:strand2)\n");
	fprintf(stderr, "\t\t-flip (flip strand of each read, i.e. might want to use with some RNA-seq)\n");
	//fprintf(stderr, "\t\t\tmap - map formatted files (i.e. output MAQ)\n");
	//fprintf(stderr, "\t\t\teland_multi - output from illumina pipeline\n");
	//fprintf(stderr, "\t\t\t\t(1:read name,2:seq,3:code,4:#MM,5:#oneMM,6:#twoMM,7:chr,\n");
	//fprintf(stderr, "\t\t-C (color space mapping with bowtie)\n");
	//fprintf(stderr, "\t\t-keep (keep one mapping of each read regardless if multiple equal mappings exist)\n");
	//fprintf(stderr, "\t\t-chrOnly (Only keep tags mapping to chromosomes whose names start with \"chr\")\n");
	fprintf(stderr, "\t\t-force5th (5th column of BED file contains # of reads mapping to position)\n");
	fprintf(stderr, "\t\t-d <tag directory> [tag directory 2] ... (add Tag directory to new tag directory)\n");
	fprintf(stderr, "\t\t-t <tag file> [tag file 2] ... (add tag file i.e. *.tags.tsv to new tag directory)\n");
	fprintf(stderr, "\t\t-single (Create a single tags.tsv file for all \"chromosomes\" - i.e. if >100 chromosomes)\n");
	fprintf(stderr, "\t\t-update (Use current tag directory for QC/processing, do not parse new alignment files)\n");
	fprintf(stderr, "\t\t-tbp <#> (Maximum tags per bp, default: no maximum)\n");
	fprintf(stderr, "\t\t-precision <1|2|3> (number of decimal places to use for tag totals, default: 1)\n");
	fprintf(stderr, "\t\t-minlen <#> and -maxlen <#> (Filter reads with lengths outside this range)\n");
	fprintf(stderr, "\n\t\tGC-bias options:\n");
	fprintf(stderr, "\t\t-genome <genome version> (To see available genomes, use \"-genome list\")\n");
	fprintf(stderr, "\t\t\t-or- (for custom genomes):\n");
	fprintf(stderr, "\t\t-genome <path-to-FASTA file or directory of FASTA files>\n");
	fprintf(stderr, "\n\t\t-checkGC (check Sequence bias, requires \"-genome\")\n");
	fprintf(stderr, "\t\t\t-freqStart <#> (offset to start calculating frequency, default: -50)\n");
	fprintf(stderr, "\t\t\t-freqEnd <#> (distance past fragment length to calculate frequency, default: +50)\n");
	fprintf(stderr, "\t\t\t-oligoStart <#> (oligo bias start)\n");
	fprintf(stderr, "\t\t\t-oligoEnd <#> (oligo bias end)\n");
	fprintf(stderr, "\t\t-normGC <target GC profile file> (i.e. tagGCcontent.txt file from control experiment)\n");
	fprintf(stderr, "\t\t\tUse \"-normGC default\" to match the genomic GC distribution\n");
	fprintf(stderr, "\t\t-normFixedOligo <oligoFreqFile> (normalize 5' end bias, \"-normFixedOligo default\" ok)\n");
//	fprintf(stderr, "\t\t-normOligo <#> (length of oligos to normalize, i.e. -normOligo 10)\n");
	//fprintf(stderr, "\t\t\tUse \"-normOligo default\" to match the genomic oligo distribution\n");
	fprintf(stderr, "\t\t-minNormRatio <#> (Minimum deflation ratio of tag counts, default: 0.25)\n");
	fprintf(stderr, "\t\t-maxNormRatio <#> (Maximum inflation ratio of tag counts, default: 2.0)\n");
	fprintf(stderr, "\t\t-iterNorm <#> (Sets -max/minNormRatio to 1 and 0, iteratively normalizes such that the\n");
	fprintf(stderr, "\t\t\tresulting distrubtion is no more than #%% different than target, i.e. 0.1,default: off)\n");

	//HiC stuff
	fprintf(stderr, "\n\tPaired-end/HiC options\n");
	
	fprintf(stderr, "\t\t-illuminaPE (when matching PE reads, assumes last character of read name is 0 or 1)\n");
	fprintf(stderr, "\t\t-removePEbg (remove paired end tags within 1.5x fragment length on same chr)\n");
	fprintf(stderr, "\t\t\t-PEbgLength <#> (remove PE  reads facing on another within this distance, default: 1.5x fragLen)\n");
	fprintf(stderr, "\t\t-restrictionSite <seq> (i.e. AAGCTT for HindIII, assign data < 1.5x fragment length to sites)\n");
	fprintf(stderr, "\t\t\tMust specify genome sequence directory too. (-rsmis <#> to specify mismatches, def: 0)\n");
	fprintf(stderr, "\t\t\t-both, -one, -onlyOne, -none (Keeps reads near restriction sites, default: keep all)\n");
	fprintf(stderr, "\t\t\t-removeSelfLigation (removes reads linking same restriction fragment)\n");
	fprintf(stderr, "\t\t\t-removeRestrictionEnds (removes reads starting on a restriction fragment)\n");
	fprintf(stderr, "\t\t\t-assignMidPoint (will place reads in the middle of HindIII fragments)\n");
	fprintf(stderr, "\t\t\t-restrictionSiteLength <#> (maximum distance from restriction site, default: 1.5x fragLen)\n");
	fprintf(stderr, "\t\t-removeSpikes <size bp> <#> (remove tags from regions with > than # times\n");
	fprintf(stderr, "\t\t\tthe average tags per size bp, suggest \"-removeSpikes 10000 5\")\n");

	fprintf(stderr, "\n");
	exit(0);
			

}
