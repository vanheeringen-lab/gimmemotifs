

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

	char* directory = NULL;
	char* inputDirectory = NULL;

	int peakSize = 200;
	int fragmentLength = TAGADJUST_AUTO;
	int fragmentLengthInput = TAGADJUST_AUTO;
	double maxtbp = -1.0;
	double maxtbpInput = -1.0;
	int centerFlag = 0;
	int nfrFlag = 0;
	int nfrSize = 100;
	int superFlag = 0;
	int revFlag = 0;
	int normalizeTagThresh = 0;
	int tssAutoCorrelationMax = 50;

	long int genomeSize = DEFAULT_GSIZE;
	char* outputfile = NULL;
	char* gtfFile = NULL;

	
	double normTotal = DEFAULT_NORM_TOTAL;

	char strand = STRAND_BOTH;

	PeakFinder* pf = new PeakFinder();

	if (argc < 2) {
		printCMD();
	}

	char* cmd = new char[10000];
	strcpy(cmd, argv[0]);
	for (int i=1;i<argc;i++) {
		strcat(cmd," ");
		strcat(cmd,argv[i]);
	}
	pf->setCMD(cmd);

	directory = argv[1];
	for (int i=1;i<argc;i++) {
		if (i==1) {
			if (argv[i][0] == '-') {
				fprintf(stderr, "!!! First argument needs to be a <tag directory>\n");
				printCMD();
			}
			pf->setDirectory(argv[i]);
			continue;
		}
		if (argv[i][0] == '-') {
			if (strcmp(argv[i],"-i")==0) {
				inputDirectory = argv[++i];
			} else if (strcmp(argv[i],"-o")==0) {
				outputfile = argv[++i];
			} else if (strcmp(argv[i],"-fragLength")==0 || strcmp(argv[i],"-len")==0) {
				i++;
				if (strcmp(argv[i],"auto")==0) {
					fragmentLength = TAGADJUST_AUTO;
				} else {
					sscanf(argv[i], "%d",&fragmentLength);
				}
			} else if (strcmp(argv[i],"-inputFragLength")==0) {
				i++;
				if (strcmp(argv[i],"auto")==0) {
					fragmentLengthInput = TAGADJUST_AUTO;
				} else {
					sscanf(argv[i], "%d",&fragmentLengthInput);
				}
			} else if (strcmp(argv[i],"-gtf")==0) {
				gtfFile = argv[++i];
			} else if (strcmp(argv[i],"-rev")==0) {
				revFlag = 1;
			} else if (strcmp(argv[i],"-center")==0) {
				centerFlag = 1;
				pf->centerFlag = 1;
			} else if (strcmp(argv[i],"-nfr")==0) {
				nfrFlag = 1;
				pf->nfrFlag = 1;
			} else if (strcmp(argv[i],"-uniqmap")==0) {
				pf->uniqMapDirectory = argv[++i];
			} else if (strcmp(argv[i],"-region")==0) {
				pf->regionFlag = 1;
			} else if (strcmp(argv[i],"-gsize")==0) {
				fprintf(stderr, "\tReminder, this recently changed: put actualy genome size, not 2x like before...\n");
				double dsize = 0.0;
				sscanf(argv[++i], "%lf",&dsize);
				genomeSize = (long long int) dsize;
				pf->setGenomeSize(genomeSize);
			} else if (strcmp(argv[i],"-size")==0) {
				sscanf(argv[++i], "%d",&peakSize);
				pf->setPeakSize(peakSize);
			} else if (strcmp(argv[i],"-expectmC")==0) {
				sscanf(argv[++i], "%lf",&(pf->expectedMethylC));
			} else if (strcmp(argv[i],"-mCthresh")==0) {
				sscanf(argv[++i], "%lf",&(pf->methylCthreshold));
			} else if (strcmp(argv[i],"-minNumC")==0) {
				sscanf(argv[++i], "%d",&(pf->minNumC));
			} else if (strcmp(argv[i],"-pseudoCount")==0) {
				sscanf(argv[++i], "%lf",&(pf->pseudoTags));
			} else if (strcmp(argv[i],"-superWindow")==0) {
				sscanf(argv[++i], "%d",&(pf->superWindow));
			} else if (strcmp(argv[i],"-superSlope")==0) {
				sscanf(argv[++i], "%lf",&(pf->superSlope));
			} else if (strcmp(argv[i],"-typical")==0) {
				pf->typicalFile = argv[++i];
			} else if (strcmp(argv[i],"-endFold")==0) {
				double f = 10.0;
				sscanf(argv[++i], "%lf",&f);
				pf->endFold = f;
			} else if (strcmp(argv[i],"-tssFold")==0) {
				double f = 3.0;
				sscanf(argv[++i], "%lf",&f);
				pf->foldTranscriptStart = f;
			} else if (strcmp(argv[i],"-bodyFold")==0) {
				double f = 3.0;
				sscanf(argv[++i], "%lf",&f);
				pf->foldTranscriptBody = f;
			} else if (strcmp(argv[i],"-norm")==0) {
				sscanf(argv[++i], "%lf",&normTotal);
				pf->normTotal = normTotal;
			} else if (strcmp(argv[i],"-confPvalue")==0) {
				sscanf(argv[++i], "%lf",&(pf->transcriptConfidenceThreshold));
			} else if (strcmp(argv[i],"-tssSize")==0) {
				sscanf(argv[++i], "%d",&(pf->tssSize));
			} else if (strcmp(argv[i],"-minBodySize")==0) {
				sscanf(argv[++i], "%d",&(pf->minBodySize));
			} else if (strcmp(argv[i],"-maxBodySize")==0) {
				sscanf(argv[++i], "%d",&(pf->maxBodySize));
			} else if (strcmp(argv[i],"-minDist")==0) {
				int mindist = 0;
				sscanf(argv[++i], "%d",&mindist);
				if (mindist < 1) {
					fprintf(stderr, "!!! -minDist set to a value less than 1 - must be 1 or greater.  Using default...\n");
					mindist = 0;
				}
				pf->minDist = mindist;
			} else if (strcmp(argv[i],"-tbp")==0) {
				sscanf(argv[++i], "%lf",&maxtbp);
			} else if (strcmp(argv[i],"-inputtbp")==0) {
				sscanf(argv[++i], "%lf",&maxtbpInput);
			} else if (strcmp(argv[i],"-method")==0) {
				i++;
				if (strcmp(argv[i],"fold")==0) {
					pf->groseqMethod = GROSEQ_METHOD_FOLD;
				} else if (strcmp(argv[i],"level")==0) {
					pf->groseqMethod = GROSEQ_METHOD_LEVEL;
				} else {
					fprintf(stderr, "Didn't recognize -method %s !!!\n", argv[i]);
					printCMD();
				}
			} else if (strcmp(argv[i],"-style")==0) {
				i++;
				if (strcmp(argv[i],"factor")==0) {
					pf->style = PEAK_STYLE_CHIPSEQ;
					centerFlag = 1;
					pf->centerFlag = 1;
				} else if (strcmp(argv[i],"super")==0) {
					pf->style = PEAK_STYLE_SUPERENHANCERS;
					pf->stitchMode = REGION_MODE_HISTONE;
					pf->regionFlag = 1;
					//pf->setPeakSize(300);
					pf->minDist = 12500;
					//pf->poisson = 0.001;
					//pf->clonalFold = 0.0;
					superFlag = 1;
					pf->filterMode = PEAKFINDER_FILTER_MODE_FDR;
				} else if (strcmp(argv[i],"histone")==0) {
					pf->style = PEAK_STYLE_HISTONE;
					pf->stitchMode = REGION_MODE_HISTONE;
					pf->regionFlag = 1;
					pf->setPeakSize(500);
					pf->minDist = 1000;
					pf->poisson = 0.001;
					pf->filterMode = PEAKFINDER_FILTER_MODE_POISSON;
				} else if (strcmp(argv[i],"dnase")==0) {
					pf->style = PEAK_STYLE_DNASE;
					centerFlag = 0;
					pf->centerFlag = 0;
					pf->setPeakSize(150);
					fragmentLength=1;
					fragmentLengthInput=1;
				} else if (strcmp(argv[i],"tss")==0) {
					pf->style = PEAK_STYLE_TSS;
					centerFlag = 1;
					pf->centerFlag = 1;
					pf->setPeakSize(150);
					pf->strand = STRAND_SEPARATE;
					fragmentLength=1;
					fragmentLengthInput=1;
					maxtbp = 0;
					maxtbpInput = 0;
					pf->clonalFold = 0.0;
					//pf->minDist = 300;
					//pf->filterMode = PEAKFINDER_FILTER_MODE_POISSON;
				} else if (strcmp(argv[i],"groseq")==0) {
					pf->style = PEAK_STYLE_GROSEQ;
					pf->strand = STRAND_SEPARATE;
					pf->stitchMode = REGION_MODE_GROSEQ;
					pf->tagThresh = 25;
					pf->regionFlag = 1;
					pf->clonalFold = 0.0;
					fprintf(stderr, "\tWill find transcripts from GR0-Seq data\n");
					pf->filterMode = PEAKFINDER_FILTER_MODE_POISSON;
				} else if (strcmp(argv[i],"mC")==0) {
					pf->style = PEAK_STYLE_METHYLC;
					if (pf->mCflag < 1) pf->mCflag = PEAKS_FIND_UNMETHYLC;
					fragmentLength=1;
					fragmentLengthInput=1;
					maxtbp = 0;
					maxtbpInput = 0;
					pf->setPeakSize(500);
					pf->clonalFold = 0.0;
				} else {
					fprintf(stderr, "Didn't recognize -style %s !!!\n", argv[i]);
					printCMD();
				}
			} else if (strcmp(argv[i],"-strand")==0) {
				i++;
				if (strcmp(argv[i],"both")==0) {
					strand = STRAND_BOTH;
				} else if (strcmp(argv[i],"separate")==0) {
					strand = STRAND_SEPARATE;
				}
				pf->strand = strand;
			} else if (strcmp(argv[i],"-methylC")==0) {
				pf->mCflag = PEAKS_FIND_METHYLC;
			} else if (strcmp(argv[i],"-mC")==0 || strcmp(argv[i],"-unmethylC")==0) {
				if (pf->mCflag < 1) pf->mCflag = PEAKS_FIND_UNMETHYLC;
			} else if (strcmp(argv[i],"-fdr")==0) {
				double f = 0.001;;
				sscanf(argv[++i], "%lf",&f);
				pf->fdr = f;
				pf->filterMode = PEAKFINDER_FILTER_MODE_FDR;
			} else if (strcmp(argv[i],"-poisson")==0) {
				double f = 0.001;;
				sscanf(argv[++i], "%lf",&f);
				pf->poisson = f;
				pf->filterMode = PEAKFINDER_FILTER_MODE_POISSON;
			} else if (strcmp(argv[i],"-tagThreshold")==0 || strcmp(argv[i],"-minReadDepth")==0) {
				double f = 0.0;
				sscanf(argv[++i], "%lf",&f);
				pf->tagThresh = f;
				pf->filterMode = PEAKFINDER_FILTER_MODE_THRESH;
			} else if (strcmp(argv[i],"-minTagThreshold")==0) {
				double f = 0.0;
				sscanf(argv[++i], "%lf",&f);
				pf->minTagThresh = f;
			} else if (strcmp(argv[i],"-ntagThreshold")==0) {
				double f = 0.0;
				sscanf(argv[++i], "%lf",&f);
				pf->tagThresh = f;
				pf->minTagThresh = f;
				pf->filterMode = PEAKFINDER_FILTER_MODE_THRESH;
				normalizeTagThresh = 1;
			} else if (strcmp(argv[i],"-P")==0) {
				double f = 1.0;
				sscanf(argv[++i], "%lf",&f);
				pf->poissonInput = f;
			} else if (strcmp(argv[i],"-F")==0) {
				double f = 4.0;
				sscanf(argv[++i], "%lf",&f);
				pf->inputFold = f;
			} else if (strcmp(argv[i],"-L")==0) {
				double f = 4.0;
				sscanf(argv[++i], "%lf",&f);
				pf->localFold = f;
			} else if (strcmp(argv[i],"-LP")==0) {
				double f = 1.0;
				sscanf(argv[++i], "%lf",&f);
				pf->poissonLocal = f;
			} else if (strcmp(argv[i],"-C")==0) {
				double f = 4.0;
				sscanf(argv[++i], "%lf",&f);
				pf->clonalFold = f;
			} else if (strcmp(argv[i],"-inputSize")==0) {
				int s = 0;
				sscanf(argv[++i], "%d",&s);
				pf->inputSize = s;
			} else if (strcmp(argv[i],"-localSize")==0) {
				int s = 0;
				sscanf(argv[++i], "%d",&s);
				pf->localSize = s;
			} else {
				printCMD();
			}
		}
	}

	if (maxtbp >= -0.1 || maxtbpInput >= -0.1) {
		pf->setMaxTBP(maxtbp,maxtbpInput);
	}

	if (outputfile != NULL) {
		if (strcmp(outputfile,"auto") == 0) {
			outputfile = new char[10000];
			if (pf->style == PEAK_STYLE_HISTONE) {
				sprintf(outputfile,"%s/regions.txt",directory);
			} else if (pf->style == PEAK_STYLE_METHYLC) {
				sprintf(outputfile,"%s/regions.txt",directory);
			} else if (pf->style == PEAK_STYLE_SUPERENHANCERS) {
				sprintf(outputfile,"%s/superEnhancers.txt",directory);
			} else if (pf->style == PEAK_STYLE_TSS) {
				sprintf(outputfile,"%s/tss.txt",directory);
			} else if (pf->style == PEAK_STYLE_GROSEQ) {
				sprintf(outputfile,"%s/transcripts.txt",directory);
				if (gtfFile == NULL) {
					gtfFile = new char[10000];
					sprintf(gtfFile,"%s/transcripts.gtf",directory);
				}
			} else {
				sprintf(outputfile,"%s/peaks.txt",directory);
			}
		}
		pf->setOutputFile(outputfile);
	}


	TagLibrary* tags = new TagLibrary(directory);
	tags->readTagDirectory();
	tags->revStrand = revFlag;
	if (normalizeTagThresh) pf->tagThresh *= tags->totalTags/pf->normTotal;

	tags->setSingleRead(1);
	if (0) {
		int tmpInt = 0;
		double* tmp = tags->getTagCountDistribution(NULL,tmpInt);
		delete []tmp;
	}
	int tagAdjust;
	if (fragmentLength == TAGADJUST_AUTO) {
		tagAdjust=TAGADJUST_AUTO;
	} else {
		tagAdjust = (int) floor(((double)fragmentLength)/2.0);
		tags->fragmentLengthEstimate = fragmentLength;
	}
	if (pf->style == PEAK_STYLE_GROSEQ) {
		tagAdjust=0;
		if (fragmentLength != TAGADJUST_AUTO) {
			tags->fragmentLengthEstimate = fragmentLength;
		} else {
			tags->fragmentLengthEstimate = 150;
		}
	}
	tags->setTagAdjust(tagAdjust);

	TagLibrary* inputTags = NULL;
	if (inputDirectory != NULL) {
		inputTags = new TagLibrary(inputDirectory);
		inputTags->readTagDirectory();
		inputTags->revStrand = revFlag;
		inputTags->setSingleRead(1);
		int tagAdjustInput;
		if (fragmentLengthInput == TAGADJUST_AUTO) {
			tagAdjustInput=TAGADJUST_AUTO;
		} else {
			tagAdjustInput = (int) floor(((double)fragmentLengthInput)/2.0);
			inputTags->fragmentLengthEstimate = fragmentLengthInput;
		}
		inputTags->setTagAdjust(tagAdjustInput);
	}

	pf->setTagLibraries(tags,inputTags);
		
	fprintf(stderr, "\tFragment Length = %d\n", tags->fragmentLengthEstimate);


	FILE* fp = stdout;
	if (pf->outputFileName != NULL) {
		fp = fopen(pf->outputFileName, "w");
		if (fp == NULL) {
			fprintf(stderr, "Could not open %s for writing!!!\n", pf->outputFileName);
			exit(1);
		}
	}

	PeakLibrary* peaks = NULL;
	if (pf->style == PEAK_STYLE_GROSEQ) {
		//peaks = pf->findGroSeqRegions();
		peaks = pf->findGroSeqTranscripts();
		peaks->checkForOutOfBoundsCoordinates(tags);
		if (gtfFile != NULL) {
			FILE* gtfFp = fopen(gtfFile,"w");
			if (gtfFp == NULL) {
				fprintf(stderr, "!!! Couldn't open %s for writing the GTF file...!!!!\n", gtfFile);
			} else {
				peaks->printGTF(gtfFp);
				fclose(gtfFp);
			}
		}
		pf->printGroSeq(fp);
		//exit(0);
	} else if (pf->style == PEAK_STYLE_METHYLC) {
		peaks = pf->findmCPeaks();
		peaks->checkForOutOfBoundsCoordinates(tags);
		pf->print(fp);
	} else {
		peaks = pf->findPeaks();
		if (centerFlag) {
			peaks->centerPeaks(tags,pf->peakSize,pf->strand);
		} else if (nfrFlag) {
			peaks->centerNFR(tags,pf->peakSize,pf->strand,nfrSize);
		}
		if (pf->style == PEAK_STYLE_TSS) {
			pf->addHeader((char*)"Dispersion Ratio");
			pf->addHeader((char*)"Periodic Ratio");
			peaks->analyzeTSSpattern(tags,pf->peakSize,pf->strand,0);

			if (0) {
				char* infoStr = new char[1000];
				for (int i=0;i<=tssAutoCorrelationMax;i++) {
					sprintf(infoStr, "%d",i);
					pf->addHeader(infoStr);
				}
				delete []infoStr;
				peaks->analyzeReadAutocorrelation(tags,pf->peakSize,pf->strand, tssAutoCorrelationMax,0);
			}			
		}
		peaks->checkForOutOfBoundsCoordinates(tags);

		if (superFlag) {
			char* notes = NULL;
			PeakLibrary* typical = NULL;
			PeakLibrary* npeaks = peaks->getSuperEnhancers(pf->superSlope,pf->superWindow,notes,typical);
			delete peaks;
			peaks = npeaks;

			if (pf->typicalFile != NULL && typical->numPeaks > 0) {
				FILE* fptypical = fopen(pf->typicalFile,"w");
				if (fptypical != NULL) {
					pf->numSuper = -1*typical->numPeaks;
					pf->print(fptypical);
					typical->print(fptypical);
					fclose(fptypical);
				}
			}
			delete typical;
			pf->numSuper = peaks->numPeaks;
			if (notes != NULL) {
				//fprintf(fp, "%s",notes);
				delete notes;
			}
		}

		pf->print(fp);
	}


	peaks->print(fp);

	if (pf->outputFileName != NULL) {
		fclose(fp);
	}
}


void printCMD() {
	fprintf(stderr, "\n\tUsage: findPeaks <tag directory> [options]\n");
	fprintf(stderr, "\n\tFinds peaks in the provided tag directory.  By default, peak list printed to stdout\n"); 
	fprintf(stderr, "\n\tGeneral analysis options:\n");
	fprintf(stderr, "\t\t-o <filename|auto> (file name for to output peaks, default: stdout)\n"); 
	fprintf(stderr, "\t\t\t\"-o auto\" will send output to \"<tag directory>/peaks.txt\", \".../regions.txt\",\n");
	fprintf(stderr, "\t\t\tor \".../transcripts.txt\" depending on the \"-style\" option\n");
	fprintf(stderr, "\t\t-style <option> (Specialized options for specific analysis strategies)\n"); 
	fprintf(stderr, "\t\t\tfactor (transcription factor ChIP-Seq, uses -center, output: peaks.txt,  default)\n"); 
	fprintf(stderr, "\t\t\thistone (histone modification ChIP-Seq, region based, uses -region -size 500 -L 0, regions.txt)\n"); 
	fprintf(stderr, "\t\t\tgroseq (de novo transcript identification from GroSeq data, transcripts.txt)\n"); 
	fprintf(stderr, "\t\t\ttss (TSS identification from 5' RNA sequencing, tss.txt)\n"); 
	fprintf(stderr, "\t\t\tdnase (Hypersensitivity [crawford style (nicking)], peaks.txt)\n"); 
	fprintf(stderr, "\t\t\tsuper (Super Enhancers, superEnhancers.txt)\n"); 
	fprintf(stderr, "\t\t\tmC (Cytosine methylation (BS-seq/methylC-seq), regions.txt)\n"); 

	fprintf(stderr, "\n\tchipseq/histone options:\n");
	fprintf(stderr, "\t\t-i <input tag directory> (Experiment to use as IgG/Input/Control)\n"); 
	fprintf(stderr, "\t\t-size <#> (Peak size, default: auto)\n"); 
	fprintf(stderr, "\t\t-minDist <#> (minimum distance between peaks, default: peak size x2)\n"); 
	fprintf(stderr, "\t\t-gsize <#> (Set effective mappable genome size, default: 2e9)\n"); 
	fprintf(stderr, "\t\t-fragLength <#|auto> (Approximate fragment length, default: auto)\n"); 
	fprintf(stderr, "\t\t-inputFragLength <#|auto> (Approximate fragment length of input tags, default: auto)\n"); 
	fprintf(stderr, "\t\t-tbp <#> (Maximum tags per bp to count, 0 = no limit, default: auto)\n");
	fprintf(stderr, "\t\t-inputtbp <#> (Maximum tags per bp to count in input, 0 = no limit, default: auto)\n");
	fprintf(stderr, "\t\t-strand <both|separate> (find peaks using tags on both strands or separate, default:both)\n");
	fprintf(stderr, "\t\t-norm # (Tag count to normalize to, default 10000000)\n");

	fprintf(stderr, "\t\t-region (extends start/stop coordinates to cover full region considered \"enriched\")\n");
	fprintf(stderr, "\t\t-center (Centers peaks on maximum tag overlap and calculates focus ratios)\n");
	fprintf(stderr, "\t\t-nfr (Centers peaks on most likely nucleosome free region [works best with mnase data])\n");
	fprintf(stderr, "\t\t\t(-center and -nfr can be performed later with \"getPeakTags\"\n");

	fprintf(stderr, "\n\tPeak Filtering options: (set -F/-L/-C to 0 to skip)\n");
	fprintf(stderr, "\t\t-F <#> (fold enrichment over input tag count, default: 4.0)\n");
	fprintf(stderr, "\t\t  -P <#> (poisson p-value threshold relative to input tag count, default: 0.0001)\n");
	fprintf(stderr, "\t\t-L <#> (fold enrichment over local tag count, default: 4.0)\n");
	fprintf(stderr, "\t\t  -LP <#> (poisson p-value threshold relative to local tag count, default: 0.0001)\n");
	fprintf(stderr, "\t\t-C <#> (fold enrichment limit of expected unique tag positions, default: 2.0)\n");
	fprintf(stderr, "\t\t-localSize <#> (region to check for local tag enrichment, default: 10000)\n");
	fprintf(stderr, "\t\t-inputSize <#> (Size of region to search for control tags, default: 2x peak size)\n"); 
	fprintf(stderr, "\t\t-fdr <#> (False discovery rate, default = 0.001)\n"); 
	fprintf(stderr, "\t\t-poisson <#> (Set poisson p-value cutoff, default: uses fdr)\n"); 
	fprintf(stderr, "\t\t-tagThreshold <#> (Set # of tags to define a peak, default: 25)\n"); 
	fprintf(stderr, "\t\t-ntagThreshold <#> (Set # of normalized tags to define a peak, by default uses 1e7 for norm)\n"); 
	fprintf(stderr, "\t\t-minTagThreshold <#> (Absolute minimum tags per peak, default: expected tags per peak)\n"); 

	fprintf(stderr, "\n\tSuperEnhancer Options: (Need to specify \"-style super\"):\n"); 
	fprintf(stderr, "\t\t-superSlope <#> (Slope threshold to identify super vs. typical enh., default: %.2lf)\n",DEFAULT_SUPERENHANCER_SLOPE_THRESHOLD); 
	fprintf(stderr, "\t\t-superWindow <#> (moving window/number of peaks to use to calculate slope, default: %d)\n",DEFAULT_SUPERENHANCER_SLOPE_WINDOW); 
	fprintf(stderr, "\t\t-typical <filename> (Output typical enhancers to this file, default: not used)\n"); 
	fprintf(stderr, "\n\tMethylC-Seq/BS-Seq options (Need to specify \"-style mC\"):\n"); 
	//fprintf(stderr, "\t\t-mC (Cytosine methylation (BS-seq))\n"); 
	fprintf(stderr, "\t\t-unmethylC / -methylC (find unmethylated/methylated regions, default: -unmethyC)\n"); 
	fprintf(stderr, "\t\t-mCthresh <#> (methylation threshold of regions, default: avg methylation/2)\n"); 
	//fprintf(stderr, "\t\t-expectmC <#> (Expected methylation level, default: average methylation)\n"); 
	fprintf(stderr, "\t\t-minNumC <#> (Minimum number of cytosines per methylation peak, default: %d\n", MINIMUM_C_PER_PEAK_METHYLC); 
	fprintf(stderr, "\n\tGroSeq Options (Need to specify \"-style groseq\"):\n"); 
	fprintf(stderr, "\t\t-tssSize <#> (size of region for initiation detection/artifact size, default: %d)\n",GROSEQ_TSSSIZE); 
	fprintf(stderr, "\t\t-minBodySize <#> (size of regoin for transcript body detection, default: %d)\n",GROSEQ_MINBODYSIZE); 
	//fprintf(stderr, "\t\t-maxBodySize <#> (size of regoin for transcript body detection, default: %d)\n",GROSEQ_MAXBODYSIZE); 
	fprintf(stderr, "\t\t-tssFold <#> (fold enrichment for new initiation dectection, default: %.1lf)\n",GROSEQ_TSSFOLDCHANGE); 
	fprintf(stderr, "\t\t-bodyFold <#> (fold enrichment for new transcript dectection, default: %.1lf)\n",GROSEQ_BODYFOLDCHANGE); 
	fprintf(stderr, "\t\t-endFold <#> (end transcript when levels are this much less than the start, default: %.1lf)\n",GROSEQ_ENDFOLDCHANGE); 
	fprintf(stderr, "\t\t-method <fold|level> (method used for identifying new transcripts, default: fold)\n");
	fprintf(stderr, "\t\t-fragLength <#> (Approximate fragment length, default: 150)\n"); 
	fprintf(stderr, "\t\t-uniqmap <directory> (directory of binary files specifying uniquely mappable locations)\n"); 
	fprintf(stderr, "\t\t\tDownload from http://biowhat.ucsd.edu/homer/groseq/\n"); 
	fprintf(stderr, "\t\t-confPvalue <#> (confidence p-value: %.2le)\n",GROSEQ_CONFPVALUE); 
	fprintf(stderr, "\t\t-minReadDepth <#> (Minimum initial read depth for transcripts, default: auto)\n"); 
	fprintf(stderr, "\t\t-pseudoCount <#> (Pseudo tag count, default: %.1lf)\n", GROSEQ_DENOVO_PSEUDOCOUNT); 
	fprintf(stderr, "\t\t-rev (reverse strand of reads - for first-strand rna-seq/gro-seq)\n");
	fprintf(stderr, "\t\t-gtf <filename> (Output de novo transcripts in GTF format)\n"); 
	fprintf(stderr, "\t\t\t\"-o auto\" will produce <dir>/transcripts.txt and <dir>/transcripts.gtf\n");
	fprintf(stderr, "\t\t\n"); 
	exit(0);

}
