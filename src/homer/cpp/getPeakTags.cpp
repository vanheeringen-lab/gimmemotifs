
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

#include "SeqTag.h"

void printCMD();

int main(int argc, char** argv) {

	char* directory = NULL;
	char* peakfile = NULL;
	int start = -1000;
	int end = 1000;
	int fixedFlag = 0;
	int offset = NULL_OFFSET;
	int tagAdjust = 0;
	char strand = BOTH_STRANDS;
	int ratioFlag = 0;
	float tbpValue = 0;
	int outputMode = OUTPUT_MODE_COUNT;
	char* vcfFile = NULL;
	char* wigFile = NULL;
	char* bedGraphFile = NULL;
	int maxIndividuals = 1000;
	int numIndividuals = 0;
	int tagAutoCorrelation = 0;
	int allSNPsFlag = 0;
	int centerFlag = 0;
	int nfrFlag = 0;
	int nfrSize = 100;
	int fragmentLength = TAGADJUST_AUTO;
	char** individuals = new char*[maxIndividuals];

	FILE* fpout = stdout;
	
	if (argc < 3) {
		printCMD();
	}
	peakfile = argv[1];
	directory = argv[2];
	int startIndex = 3;

	if (strcmp(directory,(char*)"-vcf")==0) {
		if (argc < 4) printCMD();
		vcfFile = argv[3];
		startIndex = 4;
	} else if (strcmp(directory,(char*)"-wig")==0) {
		if (argc < 4) printCMD();
		wigFile = argv[3];
		startIndex = 4;
	} else if (strcmp(directory,(char*)"-bedGraph")==0) {
		if (argc < 4) printCMD();
		bedGraphFile = argv[3];
		startIndex = 4;
	}
	for (int i=startIndex;i<argc;i++) {
		if (argv[i][0] == '-') {
			if (strcmp(argv[i],"-start")==0) {
				sscanf(argv[++i], "%d",&start);
			} else if (strcmp(argv[i],"-end")==0) {
				sscanf(argv[++i], "%d",&end);
			} else if (strcmp(argv[i],"-offset")==0) {
				sscanf(argv[++i], "%d",&offset);
			} else if (strcmp(argv[i],"-tagAdjust")==0) {
				i++;
				if (strcmp(argv[i],"auto")==0) {
					tagAdjust = TAGADJUST_AUTO;
				} else {
					sscanf(argv[i], "%d",&tagAdjust);
				}
			} else if (strcmp(argv[i],"-fragLength")==0 || strcmp(argv[i],"-len")==0) {
				i++;
				if (strcmp(argv[i],"auto")==0) {
					fragmentLength = TAGADJUST_AUTO;
				} else {
					sscanf(argv[i], "%d",&fragmentLength);
				}
			} else if (strcmp(argv[i],"-tbp")==0) {
				sscanf(argv[++i], "%f",&tbpValue);
			} else if (strcmp(argv[i],"-all")==0) {
				allSNPsFlag = 1;
			} else if (strcmp(argv[i],"-strand")==0) {
				i++;
				if (strcmp(argv[i],"+")==0) {
					strand = POSITIVE_STRAND;
				} else if (strcmp(argv[i],"-")==0) {
					strand = NEGATIVE_STRAND;
				} else if (strcmp(argv[i],"both")==0) {
					strand = BOTH_STRANDS;
				}
			} else if (strcmp(argv[i],"-center")==0) {
				centerFlag = 1;
			} else if (strcmp(argv[i],"-tagAutocorrelation")==0) {
				sscanf(argv[++i], "%d",&tagAutoCorrelation);
			} else if (strcmp(argv[i],"-nfr")==0) {
				nfrFlag = 1;
			} else if (strcmp(argv[i],"-nfrSize")==0) {
				sscanf(argv[++i], "%d",&nfrSize);
			} else if (strcmp(argv[i],"-ratio")==0) {
				ratioFlag = 1;
			} else if (strcmp(argv[i],"-fixed")==0) {
				start = 0;
				end = 0;
				fixedFlag = 1;
			} else if (strcmp(argv[i],"-peaktags")==0) {
				outputMode = OUTPUT_MODE_PEAKTAGS;
			} else if (strcmp(argv[i],"-peaksnps")==0) {
				outputMode = OUTPUT_MODE_PEAKTAGS;
			} else if (strcmp(argv[i],"-tags")==0) {
				outputMode = OUTPUT_MODE_TAGS;
			} else if (strcmp(argv[i],"-count")==0) {
				outputMode = OUTPUT_MODE_COUNT;
			} else if (strcmp(argv[i],"-individuals")==0) {
				for (i++;i<argc;i++) {
					if (argv[i][0] == '-') {
						i--;
						break;
					} else {
						individuals[numIndividuals++]=argv[i];
					}
				}
			} else {
				printCMD();
			}
		}
	}


	PeakLibrary* peaks = new PeakLibrary(peakfile,PEAK_READ_MODE_NORMAL);

	if (vcfFile != NULL) {
		// vcf mode
		int tagStart = start;
		int tagEnd = end;
		if (fixedFlag == 0) {
			peaks->setPeakTagSizeRefPos(offset,tagStart,tagEnd);
		} else {
			peaks->setPeakTagSizeFixed(tagStart,tagEnd);
			//countMode = COUNT_MODE_TBP;
		}

		SNP* baseSNP = peaks->addSNPsfromVCF(vcfFile,individuals,numIndividuals,allSNPsFlag);

		fprintf(fpout, "#PeakID\t");
		for (int i=0;i<numIndividuals;i++) {
			if (i>0) fprintf(fpout,",");
			fprintf(fpout,"%s",individuals[i]);
		}
		fprintf(fpout,"\n");

		//peaks->printSNPs(stdout,start,end,outputMode);
		if (outputMode == OUTPUT_MODE_COUNT) {
			peaks->printSNPtotals(fpout);
		} else if (outputMode == OUTPUT_MODE_PEAKTAGS) {
			peaks->printSNPs(fpout);
		}
		SNP::deleteBaseSNP(baseSNP);
	} else if (wigFile != NULL || bedGraphFile != NULL) {

		int extraSpace = (int)(2.0*fabs((double)tagAdjust));
		int tagStart = start-extraSpace;
		int tagEnd = end+extraSpace;
		if (fixedFlag == 0) {
			peaks->setPeakTagSizeRefPos(offset,tagStart,tagEnd);
		} else {
			peaks->setPeakTagSizeFixed(tagStart,tagEnd);
			//countMode = COUNT_MODE_TBP;
		}
		if (bedGraphFile != NULL) {
			peaks->processBedGraphFile(bedGraphFile, start, end, strand, fpout, outputMode);
		} else if (wigFile != NULL) {
			peaks->processWiggleFile(wigFile, start, end, strand, fpout, outputMode);
		}

	} else {

		if (centerFlag || nfrFlag || tagAutoCorrelation) {
			tagAdjust = 0;
		}

		TagLibrary* tags = new TagLibrary(directory);
		tags->readTagDirectory();
		tags->setMaxTBP(tbpValue);
		tags->setTagAdjust(tagAdjust);
		tags->setSingleRead(1);

		tagAdjust = tags->fragmentLengthEstimate;
		if (fragmentLength != TAGADJUST_AUTO) {
			tags->fragmentLengthEstimate = fragmentLength;
		}
		int peakSize = end - start;
		if (fixedFlag == 1) {
			peakSize = 0;
		}
	
		if (tagAutoCorrelation) {
			fprintf(fpout, "PeakID");
			for (int i=0;i<=tagAutoCorrelation;i++){ 
				fprintf(fpout, "\t%d",i);
			}
			fprintf(fpout,"\n");
			peaks->setDefaultPeakOrder();
			peaks->analyzeReadAutocorrelation(tags,peakSize,strand,tagAutoCorrelation,0);
			for (int i=0;i<peaks->numPeaks;i++) {
				fprintf(fpout,"%s",peaks->peakOrder[i]->name);
				fprintf(fpout,"\t%s\n",peaks->peakOrder[i]->data);
			}
		} else if (centerFlag) {
			if (nfrFlag) {
				peaks->centerNFR(tags,peakSize,strand,nfrSize);
			} else {
				peaks->centerPeaks(tags,peakSize,strand);
			}
			peaks->checkForOutOfBoundsCoordinates(tags);
			peaks->print(fpout);
		} else {	

			int countMode = COUNT_MODE_TOTAL;
			int extraSpace = (int)(2.0*fabs((double)tagAdjust));
			int tagStart = start-extraSpace;
			int tagEnd = end+extraSpace;
			if (fixedFlag == 0) {
				peaks->setPeakTagSizeRefPos(offset,tagStart,tagEnd);
			} else {
				peaks->setPeakTagSizeFixed(tagStart,tagEnd);
				//countMode = COUNT_MODE_TBP;
			}
			if (ratioFlag) {
				countMode = COUNT_MODE_RATIO;
			}

	
			if (outputMode == OUTPUT_MODE_COUNT) {
				Doubletable* rv = NULL;
				if (nfrFlag) {
					peaks->addTagLibrary(tags);
					rv = peaks->scoreNFR(tags,0,nfrSize,strand);
				} else {
					if (fixedFlag == 0) {
						peaks->setPeakTagSizeRefPos(offset,start,end);
					} else {
						peaks->setPeakTagSizeFixed(start,end);
					}
					rv = peaks->countPeakTagsLowMemory(tags,strand,countMode);
				}
				char** keys = rv->keys();
				for (int i=0;i<rv->total;i++) {
					double v = rv->search(keys[i]);
					if (v < EMPTY_DOUBLE_CHECK) {
						fprintf(fpout, "%s\tNA\n", keys[i]);
					} else {
						fprintf(fpout, "%s\t%.2f\n", keys[i],v);
					}
				}	
			} else if (outputMode == OUTPUT_MODE_PEAKTAGS || outputMode==OUTPUT_MODE_TAGS) {
				peaks->addTagLibrary(tags);
				peaks->printRelativeTags(fpout,0,start,end,outputMode);
			}
		}
	}

	return 0;

}


void printCMD() {
	fprintf(stderr, "\n\tUsage: getPeakTags <peak file> <tag directory> [options]\n");
	fprintf(stderr, "\n\tExtracts tags near each peak from the tag directory and either counts them,\n"); 
	fprintf(stderr, "\tor outputs the individual tag information for each peak\n"); 
	fprintf(stderr, "\n\tGeneral Options:\n");
	fprintf(stderr, "\t\t-start <#> (position to start counting tags relative to peak center, default: -1000)\n");
	fprintf(stderr, "\t\t-end <#> (position to stop counting tags relative to peak center, default: 1000)\n");
	fprintf(stderr, "\t\t-offset <#> (relative position of start to peak center, default: midpoint)\n");
	fprintf(stderr, "\t\t-fixed (Count tags relative to actual peak start and stop - parameters -start and\n");
	fprintf(stderr, "\t\t\t\t-end will extend from start and end of peaks, not reference position)\n");
	fprintf(stderr, "\n\tOutput Options:\n");
	fprintf(stderr, "\t\t-count (DEFAULT, Will output total/ratio tag counts to stdout)\n");
	fprintf(stderr, "\t\t\t-strand <both|+|-> (Strand [relative to peak] to count tags from, default:both)\n");
	fprintf(stderr, "\t\t\t-tagAdjust <#> (bp to shift tag positions to estimate fragment centers, default: 0)\n");
	fprintf(stderr, "\t\t\t\t'-tagAdjust auto' uses half of the estimated tag fragment length\n");
	fprintf(stderr, "\t\t\t-tbp <#> (Maximum tags per bp to count, 0 = no limit, default: 0)\n");
	fprintf(stderr, "\t\t\t-ratio (report tags per bp)\n");
	fprintf(stderr, "\t\t\t-nfr (Return nucleosome free region score over the peak center)\n");
	fprintf(stderr, "\t\t\t-nfrSize <#> (nucleosome free region size, default 100)\n");
	fprintf(stderr, "\t\t-peaktags (output locations of tags relative to peak reference position to stdout)\n");
	fprintf(stderr, "\t\t\t\tColumns: 1:peakname,2:tags separated by \",\"\n");
	fprintf(stderr, "\t\t-tags (output genomic locations of tags i.e. tags.tsv file to stdout)\n");
	fprintf(stderr, "\t\t\t\tColumns: 1:none,2:chr,3:position,4:strand,5:value,6:length/-1\n");
	fprintf(stderr, "\t\t-tagAutocorrelation <#> (output autocorrelation between read in locus, #=max distance)\n");
	fprintf(stderr, "\n\tPeak Centering: find position with specific tag features, output is a peak file\n");
	fprintf(stderr, "\t\t-center (Center peaks on summit of peak)\n");
	fprintf(stderr, "\t\t\t-nfr (Center on nucleosome free region instead of maximum tag pile-up)\n");
	fprintf(stderr, "\t\t\t-nfrSize <#> (nucleosome free region size, default 100)\n");
	fprintf(stderr, "\t\t\t-fragLength <#> (sequencing fragment length estimate, default: auto)\n");
	fprintf(stderr, "\n\tWiggle/bedGraph mode: To get tag density from these instead of tag directories\n");
	fprintf(stderr, "\n\t\tusage: getPeakTags <peak file> -wig <wiggle file> [options]\n");
	fprintf(stderr, "\t\t -or-  getPeakTags <peak file> -bedGraph <bedGraph file> [options]\n");
	fprintf(stderr, "\n\tVCF mode: For extracing SNP information from specific regions of the genome\n");
	fprintf(stderr, "\n\t\tUsage: getPeakTags <peak file> -vcf <vcf file> [options]\n");
	fprintf(stderr, "\n\t\tWhere the 2nd argument MUST be \"-vcf\" and 3rd argument MUST be a vcf file\n");
	fprintf(stderr, "\n\t\tVCF options:\n");
	fprintf(stderr, "\t\t\t-start <#> | -end <#> | -offset <#> | -fixed (same as above)\n");
	fprintf(stderr, "\t\t\t-individuals <name 1> [name 2] ... (only report genotypes for these column headers)\n");
	fprintf(stderr, "\t\t\t-all (By default, only SNPs with non-reference alleles in at least on individual are\n");
	fprintf(stderr, "\t\t\t\treported, \"-all\" will report all SNPs regardless of genotype)\n");
	fprintf(stderr, "\t\t\tOutput options:\n");
	fprintf(stderr, "\t\t\t\t-count (DEFAULT, Will output total edit distance from reference)\n");
	fprintf(stderr, "\t\t\t\t-peaksnps (will individual SNPs and genotypes associated with each peak)\n");
	fprintf(stderr, "\n");
	exit(0);
			

}

