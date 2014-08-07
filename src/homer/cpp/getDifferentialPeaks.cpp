
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
	char* backgroundDirectory = NULL;
	char* peakfile = NULL;
	double foldChange=4.0;
	double poissonChange=0.0001;
	int fixedFlag = 1;
	int tagAdjust = TAGADJUST_AUTO;
	int tagAdjustBg = TAGADJUST_AUTO;
	int size = 0;
	int mode = DIFFPEAK_MODE_DIFF;
	char strand = 2;
	float tbpValue = 0;
	float tbpValueBg = 0;
	
	if (argc < 4) {
		printCMD();
	}
	peakfile = argv[1];
	directory = argv[2];
	backgroundDirectory = argv[3];
	for (int i=1;i<argc;i++) {
		if (i==1 || i==2 || i==3) {
			if (argv[i][0] == '-') {
				fprintf(stderr, "!!! First and 2nd arguments are <peak file> and <tag directory>\n");
				printCMD();
			}
			continue;
		}
		if (argv[i][0] == '-') {
			if (strcmp(argv[i],"-size")==0) {
				i++;
				if (strcmp(argv[i],"given")==0) {
					fixedFlag = 1;
				} else {
					fixedFlag = 0;
					sscanf(argv[i], "%d",&size);
				}
			} else if (strcmp(argv[i],"-same")==0) {
				mode = DIFFPEAK_MODE_SAME;
			} else if (strcmp(argv[i],"-rev")==0) {
				mode = DIFFPEAK_MODE_REV;
			} else if (strcmp(argv[i],"-tagAdjust")==0) {
				i++;
				if (strcmp(argv[i],"auto")==0) {
					tagAdjust = TAGADJUST_AUTO;
				} else {
					sscanf(argv[i], "%d",&tagAdjust);
				}
			} else if (strcmp(argv[i],"-tagAdjustBg")==0) {
				i++;
				if (strcmp(argv[i],"auto")==0) {
					tagAdjustBg = TAGADJUST_AUTO;
				} else {
					sscanf(argv[i], "%d",&tagAdjustBg);
				}
			} else if (strcmp(argv[i],"-tbp")==0) {
				sscanf(argv[++i], "%f",&tbpValue);
			} else if (strcmp(argv[i],"-tbpBg")==0) {
				sscanf(argv[++i], "%f",&tbpValueBg);
			} else if (strcmp(argv[i],"-strand")==0) {
				i++;
				if (strcmp(argv[i],"+")==0) {
					strand = 0;
				} else if (strcmp(argv[i],"-")==0) {
					strand = 1;
				} else if (strcmp(argv[i],"both")==0) {
					strand = 2;
				}
			} else if (strcmp(argv[i],"-fixed")==0) {
				fixedFlag = 1;
			} else if (strcmp(argv[i],"-F")==0) {
				sscanf(argv[++i], "%lf",&foldChange);
			} else if (strcmp(argv[i],"-P")==0) {
				sscanf(argv[++i], "%lf",&poissonChange);
			} else {
				printCMD();
			}
		}
	}


	TagLibrary* tags = new TagLibrary(directory);
	tags->readTagDirectory();
	tags->setSingleRead(1);
	tags->setMaxTBP(tbpValue);
	tags->setTagAdjust(tagAdjust);

	TagLibrary* background = new TagLibrary(backgroundDirectory);
	background->readTagDirectory();
	background->setSingleRead(1);
	background->setMaxTBP(tbpValueBg);
	background->setTagAdjust(tagAdjustBg);

	PeakLibrary* peaks = new PeakLibrary();
	peaks->readPeakFile(peakfile,PEAK_READ_MODE_NORMAL);

	int maxadj = abs(tags->tagAdjust);
	int maxadjBg = abs(background->tagAdjust);
	if (maxadjBg > maxadj) maxadj = maxadjBg;
	int halfSize = size/2;

	if (fixedFlag) {
		peaks->setPeakTagSizeFixed(-maxadj,maxadj);
		size = 0;
		fprintf(stderr, "\tUsing fixed size peaks\n");
	} else {
		peaks->setPeakTagSizeRefPos(NULL_OFFSET,-halfSize-maxadj,halfSize+maxadj);
	}

	int strStatsFlag = 1;	
	PeakLibrary* diffPeaks = peaks->getDifferentialPeaks(tags, background, foldChange, poissonChange,
					mode, -1*halfSize, halfSize, strand,strStatsFlag);


	FILE* fp = stdout;
	fprintf(fp,"# HOMER Differential Peaks\n");
	fprintf(fp,"# Original Peak File = %s\n", peakfile);
	fprintf(fp,"# Peak File Tag Directory = %s\n", directory);
	fprintf(fp,"# Background Tag Directory = %s\n", backgroundDirectory);
	fprintf(fp,"#\n");
	fprintf(fp,"# Mode = ");
	if (mode == DIFFPEAK_MODE_SAME) {
		fprintf(fp, "Similar Peaks (-same)\n");
	} else if (mode == DIFFPEAK_MODE_REV) {
		fprintf(fp, "Reverse - Peaks higher in background (-rev)\n");
	} else {
		fprintf(fp, "Peaks enriched over background\n");
	}
	fprintf(fp,"# Differential peaks = %d\n", diffPeaks->numPeaks);
	fprintf(fp,"# Fraction of original peaks = %.2lf%%\n", 100.0*((double)diffPeaks->numPeaks)/((double)peaks->numPeaks));
	fprintf(fp,"# Original number of peaks = %d\n", peaks->numPeaks);
	fprintf(fp,"# Size of region used to count tags = ");
	if (fixedFlag) {
		fprintf(fp,"given peak sizes (-fixed)\n");
	} else {
		fprintf(fp,"%d\n", size);
	}
	if (foldChange > 0) {
		fprintf(fp, "# Fold change cutoff = %.2lf\n", foldChange);
	}
	if (poissonChange < 1.0) {
		fprintf(fp, "# Poisson p-value cutoff = %.2le\n", poissonChange);
	}
	fprintf(fp, "# strand = %d (tags counted on 0=+strand, 1=-strand, 2=both strands)\n", strand);
	if (tagAdjust == TAGADJUST_AUTO) {
		fprintf(fp, "# peak directory tagAdjust = %d\n", (tags->fragmentLengthEstimate)/2);
	} else {
		fprintf(fp, "# peak directory tagAdjust = %d\n", tagAdjust);
	}
	if (tbpValue > 0.0) fprintf(fp, "# peak directory max tags per bp = %.1lf\n", tbpValue);

	if (tagAdjustBg == TAGADJUST_AUTO) {
		fprintf(fp, "# background directory tagAdjust = %d\n", (background->fragmentLengthEstimate)/2);
	} else {
		fprintf(fp, "# background directory tagAdjust = %d\n", tagAdjustBg);
	}
	if (tbpValueBg > 0) fprintf(fp, "# background directory max tags per bp = %.1lf\n", tbpValueBg);

	fprintf(fp, "#\n");
	char* cmdline = getCMDstr(argc,argv);
	fprintf(fp, "# cmd = %s\n",cmdline);
	fprintf(fp, "#\n");

	fprintf(fp, "#PeakID\tchr\tstart\tend\tstrand\tscore\tfocus ratio/other\tTotal Tags\tBackground Tags\tFold Change vs. Background\tp-value\n");

	diffPeaks->print(fp);	

	return 0;

}


void printCMD() {
	fprintf(stderr, "\n\tUsage: getDifferentialPeaks <peak file> <target tag directory> <background tag directory>  [options]\n");
	fprintf(stderr, "\n\tExtracts tags near each peak from the tag directories and counts them,\n"); 
	fprintf(stderr, "\toutputting peaks with significantly different tag densities\n"); 
	fprintf(stderr, "\n\tGeneral Options:\n");
	fprintf(stderr, "\t\t-F <#> (fold enrichment over background tag count, default: 4.0)\n");
	fprintf(stderr, "\t\t-P <#> (poisson enrichment p-value over background tag count, default: 0.0001)\n");
	fprintf(stderr, "\t\t-same (return similar peaks instead of different peaks)\n");
	fprintf(stderr, "\t\t-rev (return peaks with higher tag counts in background instead of target library)\n");
	fprintf(stderr, "\t\t-size <#> (size of region around peak to count tags, default: -fixed)\n");
	fprintf(stderr, "\t\t-fixed (Count tags relative to actual peak start and stop, default)\n");
	fprintf(stderr, "\t\t\t\"-size given\" is the same as \"-fixed\"\n");
	fprintf(stderr, "\n\tAdvanced Options:\n");
	fprintf(stderr, "\t\t-strand <both|+|-> (Strand [relative to peak] to count tags from, default:both)\n");
	fprintf(stderr, "\t\t-tagAdjust <#> (bp to shift tag positions to estimate fragment centers, default: auto)\n");
	fprintf(stderr, "\t\t\t'-tagAdjust auto' uses half of the estimated tag fragment length\n");
	fprintf(stderr, "\t\t-tagAdjustBg <#> (bp to shift background tag positions to estimate fragment centers, default: auto)\n");
	fprintf(stderr, "\t\t\t'-tagAdjustBg auto' uses half of the estimated tag fragment length\n");
	fprintf(stderr, "\t\t-tbp <#> (Maximum tags per bp to count, 0 = no limit, default: 0)\n");
	fprintf(stderr, "\t\t-tbpBg <#> (Maximum background tags per bp to count, 0 = no limit, default: 0)\n");
	fprintf(stderr, "\n");
	exit(0);
			

}
#include "SeqTag.h"

