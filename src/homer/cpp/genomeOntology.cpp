
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

#define BUFFER 10000

void printCMD();
char** addPeakFile2List(char** pfiles, int &numFiles, char* file);
void analyzePooledPeaks(TagLibrary* tags, PeakLibrary** p, GenomeOntologyResults* gr, 
											char strand, int startIndex, int stopIndex);

int main(int argc, char** argv) {

	char** peakFiles = NULL;
	int numPeakFiles = 0;
	char* controlPeakFile = NULL;

	char strand = STRAND_BOTH;
	int maxDistance = 100;
	long int gsize = DEFAULT_GSIZE;
	//int bpFlag = 0;

	if (argc < 2) {
		printCMD();
	}
	for (int i=1;i<argc;i++) {
		if (argv[i][0] == '-') {
			if (strcmp(argv[i],"-strand")==0) {
				strand = STRAND_SEPARATE;
			} else if (strcmp(argv[i],"-gsize")==0) {
				sscanf(argv[++i], "%ld",&gsize);
			} else if (strcmp(argv[i],"-bg")==0) {
				controlPeakFile = argv[++i];
			} else if (strcmp(argv[i],"-bp")==0) {
				//bpFlag = 1;
			} else if (strcmp(argv[i],"-d")==0) {
				i++;
				if (strcmp(argv[i],"given")==0) {
					maxDistance = -1;
				} else {
					sscanf(argv[i], "%d",&maxDistance);
				}
			} else if (strcmp(argv[i],"-file")==0) {
				i++;
				FILE* fp = fopen(argv[i],"r");
				if (fp==NULL) {
					fprintf(stderr, "!!! Can't open [-file] %s\n", argv[i]);
					exit(0);
				}
				char* buf = new char[BUFFER];
				while (fgets(buf, BUFFER, fp) != NULL) {
					if (buf[0] == '#') continue;
					buf[strlen(buf)-1] = '\0';
					char* filename = new char[strlen(buf)+1];
					strcpy(filename, buf);
					peakFiles = addPeakFile2List(peakFiles, numPeakFiles, filename);
				}
				delete []buf;
				fclose(fp);
			} else {
				fprintf(stderr, "!!! Can't recognize \"%s\"\n", argv[i]);
				printCMD();
			} 
		} else {
			peakFiles = addPeakFile2List(peakFiles, numPeakFiles, argv[i]);
		}
	}

	if (numPeakFiles < 2) {
		fprintf(stderr, "!!! Can't do much with less than 2 peak files\n");
		exit(0);
	}


	int directoryFlag = 1;
	struct stat st_buf;
	stat(peakFiles[0],&st_buf);
	if (S_ISREG(st_buf.st_mode)) {
		directoryFlag = 0;
	}

	char** pfiles = &(peakFiles[1]);

	if (directoryFlag == 0) {
		// Normal version - analyze overlap of reference peak file with other stuff
		PeakLibrary* refPeaks = new PeakLibrary(10000000);
		
		refPeaks->readPeakFile(peakFiles[0], PEAK_READ_MODE_NORMAL);
		char** pfiles = &(peakFiles[1]);
		refPeaks->genomeOntology(pfiles, numPeakFiles-1, strand, maxDistance, gsize, controlPeakFile);
	
		delete refPeaks;

	} else {
		// Tag Directory Version - analyze tag counts in each set relative to expected
		// To keep things undercontrol, limit ourselves to MAX_PEAKS_AT_A_TIME


		//fprintf(stderr, "\tPerforming Tag Directory Analysis (%s)\n",peakFiles[0]);
		TagLibrary* tags = new TagLibrary(peakFiles[0]);
		tags->readTagDirectory();
		tags->setSingleRead(1);
		TagLibrary* input = NULL;
		if (controlPeakFile != NULL) {
			input = new TagLibrary(controlPeakFile);
			input->readTagDirectory();
			input->setSingleRead(1);
		}

		FILE* fp = stdout;

		int currentPeakFile = 0;
		int currentTotal = 0;
		int totalPeaks = 0;
		PeakLibrary**p = new PeakLibrary*[numPeakFiles-1];
		GenomeOntologyResults* gr = new GenomeOntologyResults[numPeakFiles-1];
		GenomeOntologyResults* grInput = new GenomeOntologyResults[numPeakFiles-1];

		for (int i=0;i<numPeakFiles-1;i++) {
			p[i] = new PeakLibrary(pfiles[i], PEAK_READ_MODE_COUNTING+i);
			totalPeaks += p[i]->numPeaks;

			gr[i].name = pfiles[i];
			gr[i].coverage = p[i]->calculateCoverage();
			gr[i].numPeaks = p[i]->numPeaks;
			gr[i].tagCount=0.0;
			gr[i].numPeaksReporting=0;
			if (input != NULL) grInput[i].tagCount = 0.0;
			if (input != NULL) grInput[i].numPeaksReporting = 0;

			if (currentTotal > 0 && currentTotal+p[i]->numPeaks > MAX_PEAKS_AT_A_TIME) {
				
				fprintf(stderr, "\t\tAnalyzing annotations %d to %d of %d (%d total peaks)\n",
							currentPeakFile+1,i,numPeakFiles-1,currentTotal);
				analyzePooledPeaks(tags,p,gr,strand,currentPeakFile,i-1);
				if (input != NULL) {
					fprintf(stderr, "\t\tAnalyzing annotations %d to %d of %d for Input (%d total peaks)\n",
							currentPeakFile+1,i,numPeakFiles-1,currentTotal);
					analyzePooledPeaks(input,p,grInput,strand,currentPeakFile,i-1);
				}
				for (int j=currentPeakFile;j<i;j++) {
					delete p[j];
				}
				currentPeakFile = i;
				currentTotal = p[i]->numPeaks;
			} else {
				currentTotal += p[i]->numPeaks;
			}
		}

		fprintf(stderr, "\t\tAnalyzing annotations %d to %d of %d (%d total peaks)\n",
							currentPeakFile+1,numPeakFiles-1,numPeakFiles-1,currentTotal);
		analyzePooledPeaks(tags,p,gr,strand,currentPeakFile,numPeakFiles-2);
		if (input != NULL) {
			fprintf(stderr, "\t\tAnalyzing annotations %d to %d of %d for Input (%d total peaks)\n",
							currentPeakFile+1,numPeakFiles-1,numPeakFiles-1,currentTotal);
			analyzePooledPeaks(input,p,grInput,strand,currentPeakFile,numPeakFiles-2);
		}

		for (int i=currentPeakFile;i<numPeakFiles-1;i++) {
			delete p[i];
		}

		fprintf(fp, "Annotation\tNumRegions\tCoverage(bp)\tTag Counts [Total=%le]",tags->totalTags );
		fprintf(fp, "\tExpected Tag Counts(gSize=%.1le)\tLog Fold Enrichment", (double)gsize);
		fprintf(fp, "\tLog P-value (+ for underrepresented)\tP-value");
		if (input != NULL) {
			fprintf(fp, "\tControl Tag Counts [Total=%le]\tExpected Input Tag Counts",input->totalTags);
			fprintf(fp, "\tLog Fold Enrichment");
			fprintf(fp, "\tLog P-value (+ underrepresented)\tP-value");
			fprintf(fp, "\tLog Fold Enrichment (Control vs. Input)");
			fprintf(fp, "\tLog P-value (Control vs. Input, + for Input enrichment)");
			fprintf(fp, "\tP-value (Control vs. Input)");
		}
		fprintf(fp, "\t%% of total");
		fprintf(fp, "\n");

		double tbp = tags->totalTags / (double)gsize;
		double tbpInput = 1;
		if (input != NULL) {
			tbpInput = input->totalTags / (double)gsize;
		}
		double totalTagsD = tags->totalTags;
		if (totalTagsD < 1.0) totalTagsD = 1.0;

		for (int i=0;i<numPeakFiles-1;i++) {
			fprintf(fp,"%s\t%d\t%lld\t%.1lf",gr[i].name,gr[i].numPeaks,gr[i].coverage,gr[i].tagCount);
	

			double expected = tbp*gr[i].coverage;
			double ee = expected;
			double vv = gr[i].tagCount;
			if (ee < 1) ee = 1.0;
			if (vv < 1) vv = 1.0;
			double lratio = log(vv/ee);
			double logp = 0.0;
			double pvalue = 1.0;
			if (vv > ee) {
				logp = ilogCumulativePoisson((int)vv,ee);
				pvalue = exp(logp);
			} else {
				logp = -1*logCumulativePoisson((int)vv,ee);
				pvalue = exp(-1*logp);
			}
			fprintf(fp, "\t%.1lf\t%.2lf\t%.2lf\t%.2le",expected,lratio,logp,pvalue);

			if (input != NULL) {

				expected = tbpInput * gr[i].coverage;
				ee = expected;
				vv = grInput[i].tagCount;
				if (ee < 1) ee = 1.0;
				if (vv < 1) vv = 1.0;
				lratio = log(vv/ee);
				logp = 0.0;
				pvalue = 1.0;
				if (vv > ee) {
					logp = ilogCumulativePoisson((int)vv,ee);
					pvalue = exp(logp);
				} else {
					logp = -1*logCumulativePoisson((int)vv,ee);
					pvalue = exp(-1*logp);
				}
				fprintf(fp, "\t%.1lf\t%.1lf\t%.2lf\t%.2lf\t%.2le",grInput[i].tagCount,expected,lratio,logp,pvalue);
		
				//comparison
				float tC = gr[i].tagCount;
				float iC = grInput[i].tagCount;

				if (tags->totalTags > input->totalTags) {
					tC /= tags->totalTags;
					tC *= input->totalTags;
				} else {
					iC /= input->totalTags;
					iC *= tags->totalTags;
				}
				if (tC < 1) tC = 1;
				if (iC < 1) iC = 1;
				lratio = log(tC/iC);
				logp = 0.0;
				pvalue = 1.0;
				if (tC > iC) {
					logp = ilogCumulativePoisson((int)tC,iC);
					pvalue = exp(logp);
				} else {
					logp = -1*logCumulativePoisson((int)tC,iC);
					pvalue = exp(-1*logp);
				}
				fprintf(fp, "\t%.2lf\t%.2lf\t%.2le",lratio,logp,pvalue);

	
			}
			double total = gr[i].tagCount/totalTagsD;
			fprintf(fp, "\t%lf%%", total*100.0);

			fprintf(fp,"\n");

		}
		delete []p;
		delete []gr;
		delete tags;
		if (input != NULL) {
			delete []grInput;
			delete input;
		}

	}

	delete []peakFiles;
	return 0;
}

void analyzePooledPeaks(TagLibrary* tags, PeakLibrary** p, GenomeOntologyResults* gr, 
											char strand, int startIndex, int stopIndex) {
	PeakLibrary* pooled = new PeakLibrary(MAX_PEAKS_AT_A_TIME);
	for (int j=startIndex;j<=stopIndex;j++) {
		pooled->addPeakLibrary(p[j]);
	}

	Doubletable* counts = tags->getPeakTagCounts(pooled,strand);
	char** keys = counts->keys();
	for (int i=0;i<counts->total;i++) {
		double v = counts->search(keys[i]);
		int index = 0;
		int slen = strlen(keys[i]);

		for (int j=0;j<slen;j++) {
			if (keys[i][j] == '-') {
				keys[i][j] = '\0';
				break;
			}
		}
		sscanf(keys[i],"%d",&index);
		index -= PEAK_READ_MODE_COUNTING;
		gr[index].tagCount += v;
		gr[index].numPeaksReporting++;
		delete [](keys[i]);
	}
	delete counts;
	delete []keys;
	delete pooled;
}

char** addPeakFile2List(char** pfiles, int &numFiles, char* file) {
	int alreadyPresent = 0;	
	char** peakfiles = new char*[numFiles+1];
	if (numFiles > 0) {
		for (int i=0;i<numFiles;i++) {
			peakfiles[i] = pfiles[i];
			if (strcmp(peakfiles[i],file)==0) {
				alreadyPresent=1;
			}
		}
	}
	delete []pfiles;
	if (alreadyPresent == 0) {
		peakfiles[numFiles] = file;
		numFiles++;
	}
	return peakfiles;
}

void printCMD()  {
	fprintf(stderr, "\n\tThis program is not normally used directly - please try running GenomeOntology.pl\n");
	fprintf(stderr, "\n\tUsage: genomeOntology [options] <primary peak file> [additional peak/ann files...]\n");

	fprintf(stderr, "\n\tCalculates significance of overlap between primary peak file and all others\n"); 
	fprintf(stderr, "\tA tag directory can be given in stead of a peak file to perform an unbiased\n"); 
	fprintf(stderr, "\tanalysis of reads, although for this tags will be counted in annotations instead\n"); 
	fprintf(stderr, "\tof calculating overlaps (different, but still useful)\n"); 
	fprintf(stderr, "\n\t\tOutput is a table of p-values/stats is sent to stdout\n"); 
	fprintf(stderr, "\n\tGeneral Options:\n");
	//fprintf(stderr, "\t\t-strand (Only merge/consider peaks on the same strand, default: either strand)\n");
	fprintf(stderr, "\t\t-d <#|given> (Maximum distance between peak centers to consider, default: 100)\n");
	fprintf(stderr, "\t\t\tUsing \"-d given\" looks for literal overlaps in peak regions, and calculates\n");
	fprintf(stderr, "\t\t\tsignificance based on the total overlap in bp between peaks/annotations\n");
	fprintf(stderr, "\t\t\tUse \"-d given\" when features have vastly different sizes (i.e. introns vs. peaks)\n");
	fprintf(stderr, "\t\t-file <filename> (file listing peak files to compare - for lots of peak files)\n");
	fprintf(stderr, "\t\t-gsize <#> (Genome size for significance calculations, default: 2e9)\n");
	//fprintf(stderr, "\t\t-bp (Overlaps calculated as coverage in bp, otherwise by default overlaps are binary)\n");
	fprintf(stderr, "\n");
	exit(0);
			
}
