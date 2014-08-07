
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
	char* inputDirectory = NULL;
	char* name = NULL;
	char* uniqMapDirectory = NULL;
	char strand = STRAND_BOTH;
	double fileSize = 5e7;
	double normLength = 100.0;
	int resolution = 1;
	double normTagCount = 1e7;
	char* color = NULL;
	int logFlag = 0;
	double pseudoCounts = 5.0;
	int method = UCSC_BEDGRAPH;
	int fragLength = TAGADJUST_AUTO;
	int inputFragLength = TAGADJUST_AUTO;
	char* outputfile = NULL;
	int lastTagFlag = 1;
	int tagAdjust = 0;
	int inputTagAdjust = 0;
	int negFlag = 0;
	double tbp = 0;
	double mintbp = 0;
	double inputtbp = 0;
	double mininputtbp = 0;
	int condenseFlag = UCSC_RES_MAX;
	int style = UCSC_CHIPSEQ;
	int circosFlag = 0;
	Peak* circosPeak = NULL;
	char* chromSizeFile = NULL;
	
	
	if (argc < 2) {
		printCMD();
	}

	directory = argv[1];

	for (int i=1;i<argc;i++) {
		if (i==1) {
			if (argv[i][0] == '-') {
				fprintf(stderr, "!!! First argument should be a <tag directory>\n");
				printCMD();
			}
			continue;
		}
		if (argv[i][0] == '-') {
			if (strcmp(argv[i],"-fsize")==0) {
				sscanf(argv[++i], "%lf",&fileSize);
			} else if (strcmp(argv[i],"-fragLength")==0) {
				i++;
				if (strcmp(argv[i],"auto")==0) {
					fragLength = TAGADJUST_AUTO;
				} else if (strcmp(argv[i],"given")==0) {
					fragLength = FRAGMENT_LEN_GIVEN;
				} else {
					sscanf(argv[i], "%d",&fragLength);
				}
			} else if (strcmp(argv[i],"-inputFragLength")==0) {
				i++;
				if (strcmp(argv[i],"auto")==0) {
					inputFragLength = TAGADJUST_AUTO;
				} else if (strcmp(argv[i],"given")==0) {
					inputFragLength = FRAGMENT_LEN_GIVEN;
				} else {
					sscanf(argv[i], "%d",&inputFragLength);
				}
			} else if (strcmp(argv[i],"-strand")==0) {
				i++;
				if (strcmp(argv[i],"both")==0) {
					strand = STRAND_BOTH;
				} else if (strcmp(argv[i],"separate")==0) {
					strand = STRAND_SEPARATE;
				} else if (strcmp(argv[i],"+")==0) {
					strand = STRAND_POSITIVE;
				} else if (strcmp(argv[i],"-")==0) {
					strand = STRAND_NEGATIVE;
				} else {
					fprintf(stderr, "\t!!!! -strand %s is not an option\n", argv[i]);
					printCMD();
				}
			} else if (strcmp(argv[i],"-log")==0 || strcmp(argv[i],"-log2")==0) {
				logFlag = 1;
			} else if (strcmp(argv[i],"-i")==0) {
				inputDirectory = argv[++i];
				fprintf(stderr, "\tUsing %s as an input experiment\n", inputDirectory);
			} else if (strcmp(argv[i],"-o")==0) {
				outputfile = argv[++i];
			} else if (strcmp(argv[i],"-normLength")==0) {
				sscanf(argv[++i], "%lf",&normLength);
			} else if (strcmp(argv[i],"-mintbp")==0) {
				sscanf(argv[++i], "%lf",&mintbp);
			} else if (strcmp(argv[i],"-mininputtbp")==0) {
				sscanf(argv[++i], "%lf",&mininputtbp);
			} else if (strcmp(argv[i],"-avg")==0) {
				condenseFlag = UCSC_RES_AVG;
			} else if (strcmp(argv[i],"-neg")==0) {
				negFlag = 1;
			} else if (strcmp(argv[i],"-uniqmap")==0) {
				uniqMapDirectory = argv[++i];
			} else if (strcmp(argv[i],"-circos")==0) {
				circosFlag = 1;
				i++;
				if (strcmp(argv[i],"genome")==0) {
					circosPeak = NULL;
				} else {
					circosPeak = new Peak();
					parseUCSCpositionStr(argv[i],circosPeak->chr,circosPeak->start,circosPeak->end);
					fprintf(stderr, "\tExtracting data for %s:%d-%d\n", circosPeak->chr, circosPeak->start, circosPeak->end);
				}
			} else if (strcmp(argv[i],"-noadj")==0) {
				normTagCount = 0.0;
			} else if (strcmp(argv[i],"-style")==0) {
				i++;
				if (strcmp(argv[i],"unmethylated")==0) {
					style = UCSC_UNMETHYLATED_CpG;
					fragLength = 1;
					normTagCount = 0.0;
				} else if (strcmp(argv[i],"methylated")==0) {
					style = UCSC_METHYLATED_CpG;
					fragLength = 1;
					normTagCount = 0.0;
				} else if (strcmp(argv[i],"tss")==0) {
					style = UCSC_TSS;
					strand = STRAND_SEPARATE;
					normLength=0;
					fragLength = 1;
				} else if (strcmp(argv[i],"dnase")==0) {
					style = UCSC_DNASE;
				} else if (strcmp(argv[i],"chipseq")==0) {
					style = UCSC_CHIPSEQ;
				} else if (strcmp(argv[i],"rnaseq")==0) {
					style = UCSC_RNASEQ;
					strand = STRAND_SEPARATE;
					fragLength = FRAGMENT_LEN_GIVEN;
				} else {
					fprintf(stderr, "\t!!!! -style %s is not an option\n", argv[i]);
					printCMD();
				}
			} else if (strcmp(argv[i],"-tbp")==0) {
				sscanf(argv[++i],"%lf",&tbp);
			} else if (strcmp(argv[i],"-inputtbp")==0) {
				sscanf(argv[++i],"%lf",&inputtbp);
			} else if (strcmp(argv[i],"-pseudo")==0) {
				sscanf(argv[++i],"%lf",&pseudoCounts);
			} else if (strcmp(argv[i],"-lastTag")==0) {
				lastTagFlag = 0;
			} else if (strcmp(argv[i],"-res")==0) {
				sscanf(argv[++i], "%d",&resolution);
			} else if (strcmp(argv[i],"-name")==0) {
				name = argv[++i];
			} else if (strcmp(argv[i],"-adjust")==0) {
				sscanf(argv[++i], "%d",&tagAdjust);
			} else if (strcmp(argv[i],"-inputAdjust")==0) {
				sscanf(argv[++i], "%d",&inputTagAdjust);
			} else if (strcmp(argv[i],"-norm")==0) {
				sscanf(argv[++i], "%lf",&normTagCount);
			} else if (strcmp(argv[i],"-bigWig")==0) {
				chromSizeFile = argv[++i];
				FILE* ff = fopen(chromSizeFile,"r");
				if (ff == NULL) {
					fprintf(stderr, "\t!!!! Must specify a valid chromosme size file (see UCSC documentation) !!!!\n");
					exit(0);
				}
				fclose(ff);
				method = UCSC_BIGWIG;
				if (fileSize < 1e8) {
					fprintf(stderr, "\tWhen using bigWig it's recommended that you set \"-fsize 1e20\"\n");
				}
			} else if (strcmp(argv[i],"-color")==0) {
				color = argv[++i];
			} else {
				printCMD();
			}
		}
	}
	if (method == UCSC_BIGWIG && strand == STRAND_SEPARATE) {
		fprintf(stderr, "!!! Error: Cannot use \"-bigWig\" and \"-strand separate\" in same command !!!\n");
		fprintf(stderr, "!!! Must do each strand separately: i.e. \"-strand +\" and then \"-strand -\" !!!\n");
		exit(0);
	}

	TagLibrary* tags = new TagLibrary(directory);
	tags->readTagDirectory();
	tags->setSingleRead(1);
	if (fragLength != TAGADJUST_AUTO) {
		tags->setFragLength(fragLength);
	}
	if (style == UCSC_DNASE) tagAdjust = -1*(tags->fragmentLengthEstimate/2);

	if (tagAdjust != 0) {
		tags->setTagAdjust(tagAdjust);
	}
	if (tbp > 0.0) tags->setMaxTBP(tbp);
	if (mintbp > 0.0) tags->setMinTBP(mintbp);
	if (tags->fragmentLengthEstimate == FRAGMENT_LEN_GIVEN) {
		fprintf(stderr, "\n\tVisualization fragment length = given\n");
	} else {
		fprintf(stderr, "\n\tVisualization fragment length = %d\n", tags->fragmentLengthEstimate);
	}


	TagLibrary* input = NULL;
	if (inputDirectory != NULL) {
		input = new TagLibrary(inputDirectory);
		input->readTagDirectory();
		input->setSingleRead(1);
		if (inputFragLength != TAGADJUST_AUTO) {
			input->setFragLength(inputFragLength);
		}
		if (style == UCSC_DNASE) inputTagAdjust = -1*(input->fragmentLengthEstimate/2);
	
		if (inputTagAdjust != 0) {
			input->setTagAdjust(inputTagAdjust);
		}
		if (inputtbp > 0) input->setMaxTBP(inputtbp);
		if (mininputtbp > 0) input->setMinTBP(mininputtbp);
		if (input->fragmentLengthEstimate == FRAGMENT_LEN_GIVEN) {
			fprintf(stderr, "\n\tVisualization of input fragment length = given\n");
		} else {
			fprintf(stderr, "\n\tVisualization of input fragment length = %d\n", input->fragmentLengthEstimate);
		}
	}	



	FILE* OUT = stdout;
	if (outputfile != NULL) {
		if (strcmp(outputfile,"auto")==0) {
			char* autoname = new char[10000];
			if (tags->name == NULL) tags->makeName();
			if (method == UCSC_BEDGRAPH) {
				sprintf(autoname,"%s/%s.ucsc.bedGraph",directory,tags->name);
			} else if (method == UCSC_BIGWIG) {
				sprintf(autoname,"%s/%s.ucsc.bigWig",directory,tags->name);
			}
			OUT = fopen(autoname, "w");
			if (OUT == NULL) {
				fprintf(stderr, "!!! Could not open autofile %s for writing!!!\n",autoname);
				exit(1);
			}
			outputfile = autoname;
			fprintf(stderr, "\tOutput file: %s\n", outputfile);
		} else {
			OUT = fopen(outputfile, "w");
			if (OUT == NULL) {
				fprintf(stderr, "!!! Could not open %s for writing!!!\n",outputfile);
				exit(1);
			}
		}

	}

	tags->printBedGraph(OUT, normTagCount, strand, resolution,negFlag, fileSize,color,name,method,
						lastTagFlag,uniqMapDirectory,style,condenseFlag,circosFlag,circosPeak,
						input, pseudoCounts, logFlag,normLength);

	if (outputfile != NULL) {
		fclose(OUT);
		if (method == UCSC_BEDGRAPH && !circosFlag) {
			char* str = new char[10000];
			fprintf(stderr, "\tGzipping file %s\n", outputfile);
			sprintf(str,"gzip -f \"%s\"",outputfile);
			(void)system(str);
		} else if (method == UCSC_BIGWIG) {
			char* str = new char[10000];
			char* tmpFile = new char[10000];
			srand(time(NULL));
			sprintf(tmpFile, "%d.tmp", rand());
			fprintf(stderr, "\tCreating bigWig from bedGraph %s\n", outputfile);
			sprintf(str,"bedGraphToBigWig \"%s\" \"%s\" \"%s\"",outputfile, chromSizeFile,tmpFile);
			(void)system(str);
			sprintf(str,"mv \"%s\" \"%s\"",tmpFile,outputfile);
			(void)system(str);
			delete []str;
			delete []tmpFile;
		}
	}

	return 0;

}


void printCMD() {
	fprintf(stderr, "\n\tUsage: makeUCSCfile <tag directory> [options]\n");
	fprintf(stderr, "\n\tCreates a bedgraph file for visualization using the UCSC Genome Browser\n"); 
	fprintf(stderr, "\n\tGeneral Options:\n");
	fprintf(stderr, "\t\t-fsize <#> (Size of file, when gzipped, default: 5e7)\n");
	fprintf(stderr, "\t\t-strand <both|separate|+|-> (control if reads are separated by strand, default: both)\n");
	fprintf(stderr, "\t\t-fragLength <# | auto | given> (Approximate fragment length, default: auto)\n");
	fprintf(stderr, "\t\t-adjust <#> (Adjust edge of tag 3' by # bp, negative for 5', default: none[good for dnase])\n");
	fprintf(stderr, "\t\t-tbp <#> (Maximum tags per bp to count, default: no limit)\n");
	fprintf(stderr, "\t\t-mintbp <#> (Minimum tags per bp to count, default: no limit)\n");
	fprintf(stderr, "\t\t-res <#> (Resolution, in bp, of file, default: 1)\n");
	fprintf(stderr, "\t\t\t-avg (report average coverage if resolution is larger than 1bp, default: max is reported)\n");
	//fprintf(stderr, "\t\t-uniqmap <directory> (Directory with uniqmap files - account for mappability of sequence)\n");
	fprintf(stderr, "\t\t-lastTag (To keep ucsc happy, last mapped tag is NOT extended by default\n");
	fprintf(stderr, "\t\t\tUsing this option will allow extending of data past the last tag position)\n");
	fprintf(stderr, "\t\t-norm <#> (Total number of tags to normalize experiment to, default: 1e7)\n");
	fprintf(stderr, "\t\t-normLength <#> (Expected length of fragment to normalize to [0=off], default: 100)\n");
	fprintf(stderr, "\t\t-noadj (Do not normalize tag counts)\n");
	fprintf(stderr, "\t\t-neg (plot negative values, i.e. for - strand transcription)\n");
	fprintf(stderr, "\t\t-CpG (Show unmethylated CpG ratios)\n");
	fprintf(stderr, "\t\t-color <(0-255),(0-255),(0-255)> (no spaces, rgb color for UCSC track, default: random)\n");
	fprintf(stderr, "\t\t-i <input tag directory> (normalize bedGraph to input data)\n");
	fprintf(stderr, "\t\t\t-pseudo <#> (Number of pseudo counts used to smooth out low coverage areas, default: 5)\n");
	fprintf(stderr, "\t\t\t-log (report log2 ratio instead of linear ratio)\n");
	fprintf(stderr, "\t\t\t-inputtbp <#>, -inputFragLength <#>, -inputAdjust <#> can also be set\n");
	//fprintf(stderr, "\t\t\t\tNOTE: For now, -inputFragLength is pegged to the -fragLength parameter)\n");
	fprintf(stderr, "\t\t-bigWig <chrom.size file> (creates a full resolution bigWig file and track line file)\n");
	fprintf(stderr, "\t\t\tThis requires bedGraphToBigWig to be available in your executable path\n");
	fprintf(stderr, "\t\t\tAlso, because how how bigWig files work, use \"-strand -\" and \"-strand +\"\n");
	fprintf(stderr, "\t\t\tin separate runs to make strand specific files: \"-strand separate\" will not work\n");
	fprintf(stderr, "\t\t\tConsider using makeBigWig.pl and makeMultiWigHub.pl if interested in bigWigs\n");
	fprintf(stderr, "\t\t-o <filename|auto> (send output to this file - will be gzipped, default: prints to stdout)\n");
	fprintf(stderr, "\t\t\tauto: this will place an appropriately named file in the tag directory\n");
	fprintf(stderr, "\t\t-name <...> (Name of UCSC track, default: auto generated)\n");
	fprintf(stderr, "\t\t-style <option> (See options below:)\n");
	fprintf(stderr, "\t\t\tchipseq (standard, default)\n");
	fprintf(stderr, "\t\t\trnaseq (strand specific, if unstranded add '-strand both' to end of command)\n");
	fprintf(stderr, "\t\t\ttss (strand specific, single bp fragment length)\n");
	fprintf(stderr, "\t\t\tdnase (fragments centered on tag position instead of downstream)\n");
	fprintf(stderr, "\t\t\tmethylated (single bp resolution of cytosine methylation)\n");
	fprintf(stderr, "\t\t\tunmethylated (single bp resolution of unmethylated cytosines)\n");
	fprintf(stderr, "\t\t-circos <chrN:XXX-YYY|genome> (output only a specific region for circos[no header])\n");
	fprintf(stderr, "\n");
	exit(0);
			

}
#include "SeqTag.h"

