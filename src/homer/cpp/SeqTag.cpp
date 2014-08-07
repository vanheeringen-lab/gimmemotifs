#include "SeqTag.h"
//Hard code homer install directory...
const char* HomerConfig::homeDirectory = "/genomics/homer/./";

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

#define BUFFER 100000

void split(char* string, char** cols, int &numCols, char delim) {
	cols[0] = string;
	numCols=1;
	char delim2 = 0;
	if (delim == WHITE_SPACE) {
		delim = '\t';
		delim2 = 32;
	}
	int len = strlen(string);
	for (int i=0;i<len;i++) {
		if (string[i] == delim || string[i] == delim2) {
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
int countChar(char* string, char delim) {
	int count = 0;
	int i=0;
	if (string == NULL) return count;
	while (string[i] != '\0') {
		if (string[i++] == delim) count++;
	}
	return count;
}
int cmpInt(const void* a, const void* b);
char* getCMDstr(int argc, char** argv) {
	int len = 0;
	for (int i=0;i<argc;i++) {
		len+=1;
		if (argv[i] != NULL) len += strlen(argv[i]);
	}
	char* str = new char[len+100];
	int pos = 0;
	for (int i=0;i<argc;i++) {
		if (i>0) {
			strcpy(&(str[pos])," ");
			pos++;
		}
		strcpy(&(str[pos]),argv[i]);
		pos += strlen(argv[i]);
	}
	return str;
}
	



HomerConfig::HomerConfig() {
	genomes = new Hashtable();
	promoters = new Hashtable();
	software = new Hashtable();
	readConfigFile();
}
HomerConfig::~HomerConfig() {
	char** keys = genomes->keys();
	for (int i=0;i<genomes->total;i++) {
		HomerConfigInfo* cfi = (HomerConfigInfo*)genomes->search(keys[i]);
		delete cfi;
		delete [](keys[i]);
	}
	delete []keys;
	delete genomes;

	keys = promoters->keys();
	for (int i=0;i<promoters->total;i++) {
		HomerConfigInfo* cfi = (HomerConfigInfo*)promoters->search(keys[i]);
		delete cfi;
		delete [](keys[i]);
	}
	delete []keys;
	delete promoters;

	keys = software->keys();
	for (int i=0;i<software->total;i++) {
		HomerConfigInfo* cfi = (HomerConfigInfo*)software->search(keys[i]);
		delete cfi;
		delete [](keys[i]);
	}
	delete []keys;
	delete software;

}
void HomerConfig::readConfigFile() {
	char* fname = new char[10000];
	sprintf(fname, "%s/config.txt",homeDirectory);
	FILE* fp = fopen(fname,"r");
	delete []fname;
	if (fp == NULL) {
		fprintf(stderr, "!!! HOMER not configured properly!!!!\n");
		fprintf(stderr, "!!! Could not open %s !!!\n",fname); 
		exit(0);
	}
	char* buffer = new char[BUFFER];
	char** cols = new char*[100];
	int numCols = 0;
	int section = HOMERCONFIG_SECTION_NONE;
	while (fgets(buffer, BUFFER, fp) != NULL) {
		split(buffer, cols, numCols,'\t');
		if (cols[0][0] == '#') continue;
		if (strcmp(cols[0],"SOFTWARE")==0) {
			section = HOMERCONFIG_SECTION_SOFTWARE;
			continue;
		} else if (strcmp(cols[0],"GENOMES")==0) {
			section = HOMERCONFIG_SECTION_GENOMES;
			continue;
		} else if (strcmp(cols[0],"PROMOTERS")==0) {
			section = HOMERCONFIG_SECTION_PROMOTERS;
			continue;
		}
		if (section == HOMERCONFIG_SECTION_SOFTWARE) {
			if (numCols < 5) continue;
			HomerConfigInfo *cfi = new HomerConfigInfo(cols[0],cols[1],cols[2],cols[3],cols[4],NULL);
			software->insert(cfi,cols[0]);
			continue;
		} else if (section == HOMERCONFIG_SECTION_GENOMES) {
			if (numCols < 6) continue;
			HomerConfigInfo *cfi = new HomerConfigInfo(cols[0],cols[1],cols[2],cols[3],cols[4],cols[5]);
			genomes->insert(cfi,cols[0]);
			continue;
		} else if (section == HOMERCONFIG_SECTION_PROMOTERS) {
			if (numCols < 6) continue;
			HomerConfigInfo *cfi = new HomerConfigInfo(cols[0],cols[1],cols[2],cols[3],cols[4],cols[5]);
			promoters->insert(cfi,cols[0]);
			continue;
		}
	}
	fclose(fp);
	delete []buffer;
	delete []cols;
}
char* HomerConfig::getGenomeDirectory(char* gname) {
	if (gname == NULL || strcmp(gname,"list")==0) listAvailableGenomes();

	HomerConfigInfo* cfi = (HomerConfigInfo*)genomes->search(gname);
	if (cfi == NULL || cfi->directory == NULL) {

	    struct stat st_buf;
    	stat(gname,&st_buf);
    	if (S_ISREG(st_buf.st_mode)) {
        	fprintf(stderr, "\tCustom genome sequence file: %s\n", gname);
    	} else if (S_ISDIR(st_buf.st_mode)) {
        	fprintf(stderr, "\tCustom genome sequence directory: %s\n", gname);
		} else {
			fprintf(stderr, "!!! Could not find genome \"%s\" !!!\n",gname);
			listAvailableGenomes();
			return NULL;
		}
		return gname;
	} else {
		char* dir = cfi->directory;
		char* dname = new char[2+strlen(dir)+strlen(homeDirectory)];
		sprintf(dname,"%s/%s",homeDirectory,dir);
		return dname;
	}
}
void HomerConfig::listAvailableGenomes() {
	fprintf(stderr, "\tAvailable Genomes:\n");
	fprintf(stderr, "\t\tName\tDescription\n");
	if (genomes->total == 0) {
		fprintf(stderr, "\n\t\t!!! No available genomes (run configureHomer.pl to load some) !!!\n");
	} else {
		char** keys = genomes->keys();
		for (int i=0;i<genomes->total;i++) {
			HomerConfigInfo* cfi = (HomerConfigInfo*)genomes->search(keys[i]);
			fprintf(stderr, "\t\t%s\t%s\n", cfi->name, cfi->description);
			delete [](keys[i]);
		}
		delete []keys;
	}
	fprintf(stderr, "\n");
}

HomerConfigInfo::HomerConfigInfo(char* nname,char* nver,char* ndesc,char* nurl,char* ndir,char* nparam) {
	name = new char[strlen(nname)+1];
	strcpy(name,nname);
	version = new char[strlen(nver)+1];
	strcpy(version,nver);
	description = new char[strlen(ndesc)+1];
	strcpy(description,ndesc);
	url = new char[strlen(nurl)+1];
	strcpy(url,nurl);
	directory = new char[strlen(ndir)+1];
	strcpy(directory,ndir);
	if (nparam == NULL) {
		parameters = NULL;
		numParameters = 0;
	} else {
		char** line = new char*[100];
		split(nparam, line, numParameters,'\t');
		parameters = new char*[numParameters];
		for (int i=0;i<numParameters;i++) {
			parameters[i] = new char[strlen(line[i])+1];
			strcpy(parameters[i],line[i]);
		}
		delete []line;
	}
}
HomerConfigInfo::~HomerConfigInfo() {
	if (name != NULL) delete []name;
	if (version != NULL) delete []version;
	if (directory != NULL) delete []directory;
	if (description != NULL) delete []description;
	if (url != NULL) delete []url;
	for (int i=0;i<numParameters;i++) {
		delete [](parameters[i]);
	}
	if (parameters != NULL) delete []parameters;
}

Genome::Genome() {
	name = NULL;
	directory = NULL;
	totalbp = 2000000000;
	totalMappablebp = 2000000000;
}
Genome::~Genome() {
	if (name != NULL) delete []name;
	if (directory != NULL) delete []directory;
}
void Genome::setName(char* newname) {
	if (newname == NULL) return;
	if (name != NULL) delete []name;
	name = new char[strlen(newname)+1];
	strcpy(name, newname);
}
void Genome::setDirectory(char* newname) {
	if (newname == NULL) return;
	if (directory != NULL) delete []name;
	directory = new char[strlen(newname)+1];
	strcpy(directory, newname);
}

// class PeakFinder ---------------------------------------------------------------------

PeakFinder::PeakFinder() {
	name = NULL;
	directory = NULL;
	uniqMapDirectory = NULL;
	outputFileName = NULL;
	genome = NULL;
	peakSize = 0;
	localSize = 10000;
	inputSize = 0;
	tagThresh = 0;
	mCflag = 0;
	minTagThresh = -1.0;
	minDist = 0;
	groseqMethod = GROSEQ_METHOD_FOLD;
	totalTags = 0.0;
	totalInputTags = 0.0;
	tagsUsedForClustering = 0.0;
	maxtbp = 0.0;
	maxtbpInput = 0.0;
	mintbp = 0.0;
	mintbpInput = 0.0;
	tagAdjust=0;
	tagAdjustInput=0;
	strand=BOTH_STRANDS;
	gsize = DEFAULT_GSIZE;
	expectedMethylC = AVERAGE_METHYLC;
	minNumC = MINIMUM_C_PER_PEAK_METHYLC;
	//gsize *= 2;
	fdr = 0.001;
	fdrThresh = 0.0;
	normTotal = DEFAULT_NORM_TOTAL;
	poisson = 0.0;
	poissonInput = 0.0001;
	poissonLocal = 0.0001;
	poissonThresh = 0.0;
	tbpAuto=1;
	pvalue = 1.0;
	method = 0;
	tbpThreshold = 0.01;
	tbp = 0.0;
	tbpInput = 0.0;
	numSuper = 0;
	superWindow = DEFAULT_SUPERENHANCER_SLOPE_WINDOW;
	superSlope = DEFAULT_SUPERENHANCER_SLOPE_THRESHOLD;
	typicalFile = NULL;
	diffMode = DIFFPEAK_MODE_DIFF;
	filterMode = PEAKFINDER_FILTER_MODE_FDR;
	fdrTable = NULL;
	poissonTable = NULL;
	tagsInPeaks = 0.0;
	fdrSize = PEAKFINDER_FDRSIZE;
	extraheader = NULL;
	stitchMode = REGION_MODE_HISTONE;
	hysteresis = 1.0;
	foldTranscriptStart = GROSEQ_TSSFOLDCHANGE;
	foldTranscriptBody = GROSEQ_BODYFOLDCHANGE;
	transcriptConfidenceThreshold = GROSEQ_CONFPVALUE;
	endFold = GROSEQ_ENDFOLDCHANGE;
	tssSize = GROSEQ_TSSSIZE;
	minBodySize = GROSEQ_MINBODYSIZE;
	maxBodySize = GROSEQ_MAXBODYSIZE;
	pseudoTags = GROSEQ_DENOVO_PSEUDOCOUNT;
	methylCthreshold = -10000.0;

	centerFlag=0;
	nfrFlag=0;
	nfrSize=100;
	style = PEAK_STYLE_CHIPSEQ;
	
	inputFold = 4.0;
	localFold = 4.0;
	clonalFold = 2.0;
	numPeaks = 0;
	numFilterPeaks = 0;
	inputPeaks =0;
	clonalPeaks = 0;
	localPeaks = 0;

	cmd = NULL;

	tags=NULL;
	input=NULL;
}
PeakFinder::~PeakFinder() {
	if (genome != NULL) delete []genome;
	if (name != NULL) delete []name;
	if (extraheader != NULL) delete []extraheader;
	if (directory != NULL) delete []directory;
}
void PeakFinder::addHeader(char* str) {
	if (str == NULL) return;
	char* header = extraheader;
	int len = 0;
	if (header != NULL) len = strlen(header)+1;
	extraheader = new char[len + strlen(str)+2];
	extraheader[0] = '\0';
	if (header != NULL) {
		strcpy(extraheader, header);
		delete []header;
	}
	strcat(extraheader, "\t");
	strcat(extraheader, str);
}
void PeakFinder::print(FILE* fp) {
	if (fp == NULL) return;

	fprintf(fp, "# HOMER Peaks\n");
	fprintf(fp, "# Peak finding parameters:\n");
	fprintf(fp, "# tag directory = %s\n",directory);
	fprintf(fp, "#\n");

	if (numSuper > 0) {
		fprintf(fp, "# total super enhancers = %d\n", numSuper);
		fprintf(fp, "# super enhancer stitching window = %d\n", minDist);
		fprintf(fp, "# super enhancer slope threshold = %.2lf\n", superSlope);
		fprintf(fp, "# super enhancer slope moving window = %d\n", superWindow);
		fprintf(fp, "#\n");
	} else if (numSuper < 0) {
		fprintf(fp, "# total typical(non-super) enhancers = %d\n", -1*numSuper);
		fprintf(fp, "# super enhancer stitching window = %d\n", minDist);
		fprintf(fp, "# super enhancer slope threshold = %.2lf\n", superSlope);
		fprintf(fp, "# super enhancer slope moving window = %d\n", superWindow);
		fprintf(fp, "#\n");
	}
	fprintf(fp, "# total peaks = %d\n",numPeaks);
	fprintf(fp, "# peak size = %d\n",peakSize);
	if (strand == STRAND_BOTH) {
		fprintf(fp, "# peaks found using tags on both strands\n");
	} else {
		fprintf(fp, "# peaks found separately for each strand\n");
	}
	fprintf(fp, "# minimum distance between peaks = %d\n",minDist);
	fprintf(fp, "# fragment length = %d\n", tagAdjust*2);
	//fprintf(fp, "# genome = %s\n",genome);
	fprintf(fp, "# genome size = %.0lld\n",gsize);
	fprintf(fp, "# Total tags = %.1lf\n", totalTags);
	fprintf(fp, "# Total tags in peaks = %.1lf\n", tagsInPeaks);
	double ratio = tagsInPeaks/totalTags*100.0;
	fprintf(fp, "# Approximate IP efficiency = %.2lf%%\n", ratio);
	fprintf(fp, "# tags per bp = %.6f\n", tbp);
	fprintf(fp, "# expected tags per peak = %.3f\n", tpp);
	fprintf(fp, "# maximum tags considered per bp = %.1f\n", maxtbp);
	fprintf(fp, "# effective number of tags used for normalization = %.1f\n", normTotal);
	if (regionFlag) {
		fprintf(fp, "# Individual peaks have been stitched together into variable length regions\n");
	}
	if (centerFlag) {
		fprintf(fp, "# Peaks have been centered at maximum tag pile-up\n");
	}
	if (nfrFlag) {
		fprintf(fp, "# Peaks/Regions have been centered on the most likely Nucleosome Free Region\n");
	}
	if (filterMode == PEAKFINDER_FILTER_MODE_FDR) {
		fprintf(fp, "# FDR rate threshold = %.9lf\n", fdr);
		fprintf(fp, "# FDR effective poisson threshold = %le\n", poisson);
		fprintf(fp, "# FDR tag threshold = %.1f\n", fdrThresh);
	} else if (filterMode == PEAKFINDER_FILTER_MODE_POISSON) {
		fprintf(fp, "# Poisson p-value threshold = %le\n", poisson);
		fprintf(fp, "# Poisson tag threshold = %.1f\n", poissonThresh);
	} else if (filterMode == PEAKFINDER_FILTER_MODE_POISSON) {
		fprintf(fp, "# Manual tag threshold = %.1f\n", tagThresh);
	}
	fprintf(fp, "# number of putative peaks = %d\n",numFilterPeaks);
	if (input != NULL && (inputFold > 0 || poissonInput < 1.0)) {
		fprintf(fp, "#\n");
		fprintf(fp, "# input tag directory = %s\n",input->directory);
		if (inputFold > 0) fprintf(fp, "# Fold over input required = %.2lf\n",inputFold);
		if (poissonInput < 1.0) fprintf(fp, "# Poisson p-value over input required = %.2le\n",poissonInput);
		fprintf(fp, "# Putative peaks filtered by input = %d\n",inputPeaks);
	}
	if (localFold > 0 || poissonLocal < 1.0) {
		fprintf(fp, "#\n");
		fprintf(fp, "# size of region used for local filtering = %d\n", localSize);
		if (localFold > 0) fprintf(fp, "# Fold over local region required = %.2lf\n",localFold);
		if (poissonLocal < 1.0) fprintf(fp, "# Poisson p-value over local region required = %.2le\n",poissonLocal);
		fprintf(fp, "# Putative peaks filtered by local signal = %d\n",localPeaks);
	}
	if (clonalFold > 0) {
		fprintf(fp, "#\n");
		fprintf(fp, "# Maximum fold under expected unique positions for tags = %.2lf\n",clonalFold);
		fprintf(fp, "# Putative peaks filtered for being too clonal = %d\n",clonalPeaks);
	}

	if (uniqMapDirectory != NULL) fprintf(fp, "#\n# uniqMapDirectory = %s\n",uniqMapDirectory);
	fprintf(fp, "#\n# cmd = %s\n",cmd);
	fprintf(fp, "#\n# Column Headers:\n");
	fprintf(fp, "#PeakID\tchr\tstart\tend\tstrand");
	if (style == PEAK_STYLE_METHYLC) {
		fprintf(fp,"\tAvg. mC%%");
	} else {
		fprintf(fp,"\tNormalized Tag Count");
	}
	if (numSuper != 0) {
		fprintf(fp, "\tsuperEnhancer slope");
	} else if (style == PEAK_STYLE_METHYLC) {
		fprintf(fp, "\ttotal cytosines");
	} else if (regionFlag) {
		fprintf(fp, "\tregion size");
	} else if (centerFlag) {
		fprintf(fp, "\tfocus ratio");
	} else {
		fprintf(fp, "\tNot used");
	}
	//fprintf(fp, "\tUniquely mappable (bp)");

	if (extraheader != NULL) {
		fprintf(fp, "%s", extraheader);
	}
		


	fprintf(fp, "\n");
}

void PeakFinder::printGroSeq(FILE* fp) {
	if (fp == NULL) return;

	fprintf(fp, "# HOMER Nascent RNA (GroSeq) transcripts\n");
	fprintf(fp, "# Transcript finding parameters:\n");
	fprintf(fp, "# tag directory = %s\n",directory);
	fprintf(fp, "#\n");

	fprintf(fp, "# total transcripts = %d\n",numPeaks);
	if (strand == STRAND_BOTH) {
		fprintf(fp, "# transcripts found using tags on both strands\n");
	} else {
		fprintf(fp, "# transcripts found separately for each strand\n");
	}
	fprintf(fp, "# fragment length = %d\n", tagAdjust*2);
	//fprintf(fp, "# genome = %s\n",genome);
	//fprintf(fp, "# genome size = %.0lld\n",gsize);
	fprintf(fp, "# Total tags = %.1lf\n", totalTags);
	//fprintf(fp, "# Total tags in transcripts = %.1lf\n", tagsInPeaks);
	//double ratio = tagsInPeaks/totalTags*100.0;
	//fprintf(fp, "# Approximate IP effeciency = %.2lf%%\n", ratio);
	//fprintf(fp, "# tags per bp = %.6f\n", tbp);
	//fprintf(fp, "# expected tags per peak = %.3f\n", tpp);
	if (maxtbp > 0) {
		fprintf(fp, "# maximum tags considered per bp = %.1f\n", maxtbp);
	} else {
		fprintf(fp, "# maximum tags considered per bp = no limit\n");
	}
	//fprintf(fp, "# effective number of tags used for normaliztion = %.1f\n", normTotal);
	fprintf(fp, "# TSS region size = %d\n",tssSize);
	fprintf(fp, "# minimum transcript body resolution = %d\n",minBodySize);
	fprintf(fp, "# maximum transcript body resolution = %d\n",maxBodySize);
	fprintf(fp, "# transcript body size tag threshold = %.1lf\n",tagThresh);
	fprintf(fp, "# transcript body size tag threshold confidence p-value = %.8lf\n",
											transcriptConfidenceThreshold);
	fprintf(fp, "# initiation fold change threshold =  %.3f\n", foldTranscriptStart);
	fprintf(fp, "# transcript body fold change threshold =  %.3f\n", foldTranscriptBody);
	fprintf(fp, "# transcript termination fold change threshold =  %.3f\n", endFold);

	if (uniqMapDirectory != NULL) fprintf(fp, "#\n# uniqMapDirectory = %s\n",uniqMapDirectory);

	fprintf(fp, "#\n# cmd = %s\n",cmd);
	fprintf(fp, "#\n# Column Headers:\n");
	fprintf(fp, "#PeakID\tchr\tstart\tend\tstrand\tInitial read depth\tlength (bp)\t");

	if (extraheader != NULL) {
		fprintf(fp, "%s", extraheader);
	}

	fprintf(fp, "\n");
}

double PeakFinder::getGroSeqThreshold(double fold, double confidenceThreshold) {

	double t= 0.0;
	FILE* tablefp = NULL;	
	if (tablefp != NULL) fprintf(tablefp, "Count\tFraction\ttotal\n");
	for (int i=1;i<200;i++) {
		double sumExtreme=0.0;
		double total=0.0;
		double I = (double)i;
		for (int j=0;j<2000;j++) {
			double J = (double)j;
			double p = exp(logPoisson(j,I));
			total+=p;
			if (j==0) {
				sumExtreme += p;
			} else {
				double f = J/I;
				if (J<I) f = I/J;
				if (f >= fold) sumExtreme += p;
			}
		}
		if (sumExtreme < confidenceThreshold) {
			t = I;
			break;
		}
		if (tablefp != NULL) fprintf(tablefp, "%d\t%le\t%lf\n", i, sumExtreme, total);
	}
	fprintf(stderr, "\tTag threshold for %.1lf-fold set at %.1lf (%.2le)\n", fold,t,confidenceThreshold);
	return t;

}


void PeakFinder::setTagLibraries(TagLibrary* exp, TagLibrary* i) {
	tags = exp;
	input = i;
}
void PeakFinder::setMaxTBP(double maxT, double maxI) {
	tbpAuto = 0;
	maxtbp = maxT;
	maxtbpInput = maxI;
	if (maxtbp < 0) maxtbp = maxtbpInput;
	if (maxtbpInput < 0) maxtbpInput = maxtbp;
}
void PeakFinder::setGenomeSize(long long int genomeSize) {
	gsize = genomeSize;
}
void PeakFinder::setPeakSize(int size) {
	peakSize = size;
}
void PeakFinder::setCMD(char* str) {
	if (str == NULL) return;
	cmd=new char[strlen(str)+1];
	strcpy(cmd,str);
}
void PeakFinder::setOutputFile(char* fname) {
	if (fname == NULL) return;
	outputFileName=new char[strlen(fname)+1];
	strcpy(outputFileName,fname);
}
void PeakFinder::setName(char* newname) {
	if (newname == NULL) return;
	name=new char[strlen(newname)+1];
	strcpy(name,newname);
}
void PeakFinder::setDirectory(char* newname) {
	if (newname == NULL) return;
	directory=new char[strlen(newname)+1];
	strcpy(directory,newname);
}
void PeakFinder::setGenome(Genome* g) {
	genome = g;
	gsize = genome->totalMappablebp;
}
void PeakFinder::determineMaxTBP() {
	maxtbp = 1.0;
	maxtbpInput = 1.0;
	fprintf(stderr, "\tTotal Tags = %.1lf\n", totalTags);
	fprintf(stderr, "\tTags per bp = %.6lf\n", tbp);
	if (tbp > tbpThreshold) maxtbp = floor(tbp/tbpThreshold);
	if (tbpInput > tbpThreshold) maxtbpInput = floor(tbp/tbpThreshold);
	fprintf(stderr, "\tMax tags per bp set automatically to %.1f\n",maxtbp); 
}
void PeakFinder::checkParameters() {
	if (tags != NULL) {

		if (tags->gsizeEstimate > 10 && tags->gsizeEstimate < gsize &&
						gsize  ==  ((long long int)DEFAULT_GSIZE)) {
			gsize = tags->gsizeEstimate;
			fprintf(stderr, "\t!!! Estimated genome size (from tag directory) is smaller than default\n");
			fprintf(stderr, "\t    genome size.  Using estimate (%lld) [to change specify -gsize]\n",gsize);
			//gsize *= 2;
		}
		tbp = tags->totalTags/(double)gsize;
		//fprintf(stderr, "%lf\t%lf\t%lf\n", tbp,tags->totalTags,(double)gsize);
		if (strand != STRAND_BOTH) tbp /= 2.0;
		tags->tbp = tbp;
		totalTags = tags->totalTags;
		tagAdjust = (int)(tags->fragmentLengthEstimate/2);
	}
	if (input != NULL) {
		tbpInput = input->totalTags/(double)gsize;
		if (strand != STRAND_BOTH) tbpInput /= 2.0;
		input->tbp = tbpInput;
		totalInputTags = input->totalTags;
		//tagAdjustInput = input->fragmentLengthEstimate;
		tagAdjustInput = (int)(input->fragmentLengthEstimate/2);
	}
	if (peakSize == 0) {
		if (tags != NULL) {
			peakSize = tags->peakSizeEstimate;
		} else {
			peakSize = 200;
		}
	}
	if (minDist == 0) {
		minDist = (int) floor(peakSize * FINDPEAKS_MINDIST_DEFAULT);
	}
	if (tbpAuto) determineMaxTBP();
	tags->setMaxTBP(maxtbp);
	if (input != NULL) {
		tags->setMaxTBP(maxtbp);
	}

	if (strand == STRAND_SEPARATE) {
		fprintf(stderr, "\tFinding tags on separate strands: doubling effective genome size\n");
	}

}
PeakLibrary* PeakFinder::findGroSeqTranscripts() {
	fprintf(stderr, "\tGRO-Seq Transcript Analysis (new version June 2012)...\n");
	//fprintf(stderr, "\tGRO-Seq Transcript Analysis...\n");
	if (tagThresh < 1.0) {
		tagThresh = getGroSeqThreshold(foldTranscriptBody,transcriptConfidenceThreshold);
	}
	PeakLibrary* regions = tags->findGroSeqTranscripts(style,input,uniqMapDirectory,strand,
					tssSize,minBodySize,maxBodySize,tagThresh,foldTranscriptStart,
					foldTranscriptBody,endFold,pseudoTags,inputFold,groseqMethod);
	tagAdjust = tags->fragmentLengthEstimate/2;
	totalTags = tags->totalTags;
	numPeaks = regions->numPeaks;
	fprintf(stderr, "\t%d total transcripts discovered\n", numPeaks);
	//PeakLibrary* TagLibrary::findRegions(int mode, TagLibrary* input, char* uniqMapDirectory, char strand,
				//int minSize, int maxGapSize, double hysteresis, double threshold, double inputFold)
	return regions;
}


PeakLibrary* PeakFinder::findGroSeqRegions() {
	fprintf(stderr, "\tGRO-Seq Transcript Analysis...\n");
	//fprintf(stderr, "\tGRO-Seq Transcript Analysis...\n");
	if (tagThresh < 1.0) {
		tagThresh = getGroSeqThreshold(foldTranscriptBody,transcriptConfidenceThreshold);
	}
	PeakLibrary* regions = tags->findGroSeqRegions(style,input,uniqMapDirectory,strand,
					tssSize,minBodySize,maxBodySize,tagThresh,foldTranscriptStart,
					foldTranscriptBody,endFold,pseudoTags,inputFold);
	tagAdjust = tags->fragmentLengthEstimate/2;
	totalTags = tags->totalTags;
	numPeaks = regions->numPeaks;
	fprintf(stderr, "\t%d total transcripts discovered\n", numPeaks);
	//PeakLibrary* TagLibrary::findRegions(int mode, TagLibrary* input, char* uniqMapDirectory, char strand,
				//int minSize, int maxGapSize, double hysteresis, double threshold, double inputFold)
	return regions;
}

PeakLibrary* PeakFinder::findmCPeaks() {

	if (expectedMethylC < AVERAGE_METHYLC+1) {
		expectedMethylC = tags->totalTags/((double)tags->totalPositions);
	}
	fprintf(stderr, "\tAverage(or expected) mC = %.2lf%%\n", expectedMethylC*100.0);
	fprintf(stderr, "\tMinimum size of peaks = %d\n", peakSize);
	if (methylCthreshold < -1000.0) {
		if (mCflag == PEAKS_FIND_UNMETHYLC) {
			methylCthreshold = expectedMethylC/2.0;
		} else {
			methylCthreshold = expectedMethylC + (1.0-expectedMethylC)*1.0/2.0;
		}
		fprintf(stderr, "\tmC%% threshold automatically set:\n");
	} 
	fprintf(stderr, "\tmC%% threshold = %.2lf%%\n", methylCthreshold*100.0);
	fprintf(stderr, "\tMinimum number of cytosines per peak = %d\n", minNumC);
	if (mCflag == PEAKS_FIND_UNMETHYLC) {
		fprintf(stderr, "\tWill find unmethylated regions...\n");
	} else if (mCflag == PEAKS_FIND_METHYLC) {
		fprintf(stderr, "\tWill find methylated regions...\n");
	}
	fprintf(stderr, "\n");

	PeakLibrary* peaks = tags->findmCPeaks(peakSize, strand, mCflag, methylCthreshold, minNumC, input);
	numPeaks = peaks->numPeaks;

	return peaks;

}

PeakLibrary* PeakFinder::findPeaks() {

	if (uniqMapDirectory != NULL) {
		realGsize = 0;
		UniqMapChrs::getUniqMapStats(uniqMapDirectory,realGsize, gsize, NULL, NULL);
		gsize /= 2;
		realGsize /= 2;
		fprintf(stderr, "\tMappable genome size set to %lld (out of %lld)\n", gsize,realGsize);
	}
	checkParameters();

	int curMinDist = minDist;
	if (regionFlag == 0) {
		fprintf(stderr, "\tFinding peaks of size %d, no closer than %d\n", peakSize, minDist);
	} else {
		fprintf(stderr, "\tInitially finding peaks of size %d bp for stitching into regions", peakSize);
		fprintf(stderr, " with %d bp stitching size)\n", minDist);
		curMinDist = (int)(-1*(peakSize/PEAKFRACTION_REGION_MINDIST));
		// negative minDist will signal "findPutativePeaks" to find peaks with purpose of stitching
		// them together later
	}
	if (minTagThresh < 0) {
		minTagThresh = tbp*peakSize-1.0;
	}
	PeakLibrary* putativePeaks = tags->findPutativePeaks(peakSize,curMinDist,strand,minTagThresh);

	tagsUsedForClustering = tags->getAdjustedTagTotal();
	fprintf(stderr, "\t\tTags Used for cluster (less clonal tags) = %.1lf / %.1lf\n",
							tagsUsedForClustering, totalTags);

	if (uniqMapDirectory != NULL) {
		char s = STRAND_BOTH;
		if (strand == STRAND_SEPARATE) s = STRAND_POSITIVE;
		long long int tmpgsize = 0;
		long long int actualGsize = 0;
		putativePeaks->getMappabilityCount(uniqMapDirectory,tags->tagAdjust,s,actualGsize,tmpgsize);
		actualGsize /= 2;
		fprintf(stderr, "\tMappable genome size set at %lld (of %lld)\n", gsize, actualGsize);
	}

	// remove peaks not meeting fdr levels or poisson p-value or absolute tag thresh
	approxFdrTable();
	addHeader((char*)"findPeaks Score");
	PeakLibrary* filteredPeaks = filterPeaks(putativePeaks);
	numFilterPeaks = filteredPeaks->numPeaks;
	delete putativePeaks;

	int halfPeakSize = (peakSize+1)/2;
	filteredPeaks->setPeakTagSizeRefPos(NULL_OFFSET,-1*halfPeakSize,halfPeakSize);	

	if (input != NULL && inputFold > 0) {
		int before = filteredPeaks->numPeaks;
		tags->setMaxTBP(0);
		input->setMaxTBP(0);
		if (tags->totalTags > input->totalTags) {
			addHeader((char*)"Total Tags (normalized to Control Experiment)");
			addHeader((char*)"Control Tags");
		} else {
			addHeader((char*)"Total Tags");
			addHeader((char*)"Control Tags (normalized to IP Experiment)");
		}
		addHeader((char*)"Fold Change vs Control");
		addHeader((char*)"p-value vs Control");
		int halfInputSize = halfPeakSize;
		if (inputSize == 0) {
			halfInputSize = halfPeakSize*2;
		} else {
			halfInputSize = inputSize/2;
		}
		int strOutputFlag =1;
		if (regionFlag) strOutputFlag = 0;
		PeakLibrary* inputFiltered = filteredPeaks->getDifferentialPeaks(tags,input,
				inputFold,poissonInput,diffMode,-1*halfInputSize, halfInputSize, strand,strOutputFlag);
		delete filteredPeaks;
		filteredPeaks = inputFiltered;
		inputPeaks = before-filteredPeaks->numPeaks;
	}

	if (regionFlag && style != PEAK_STYLE_SUPERENHANCERS) {
		fprintf(stderr, "\tStitching together putative peaks into regions\n");
		PeakLibrary* regions = filteredPeaks->stitchRegions(minDist,stitchMode);
		delete filteredPeaks;
		filteredPeaks = regions;

		if (input != NULL && inputFold > 0) {
			if (style == PEAK_STYLE_SUPERENHANCERS) {
				//don't filter for superenhancers
			} else {
				fprintf(stderr, "\tChecking regions against input...\n");
				int before = filteredPeaks->numPeaks;
				tags->setMaxTBP(0);
				input->setMaxTBP(0);
				filteredPeaks->setPeakTagSizeFixed(0,0);
				int strOutputFlag = 1;
				PeakLibrary* inputFiltered = filteredPeaks->getDifferentialPeaks(tags,input,
					inputFold,poissonInput,diffMode,0,0, strand,strOutputFlag);
				delete filteredPeaks;
				filteredPeaks = inputFiltered;
				inputPeaks = before-filteredPeaks->numPeaks;
			}
		}

	} else {
		if (localFold > 0) {
			int before = filteredPeaks->numPeaks;
			tags->setMaxTBP(0);
			addHeader((char*)"Fold Change vs Local");
			addHeader((char*)"p-value vs Local");
			PeakLibrary* localFiltered = filteredPeaks->filterLocalPeaks(tags,
						peakSize,localSize,localFold,poissonLocal,diffMode, strand);
			delete filteredPeaks;
			filteredPeaks = localFiltered;
			localPeaks = before-filteredPeaks->numPeaks;
		}
	}

	if (clonalFold > 0) {
		int before = filteredPeaks->numPeaks;
		tags->setMaxTBP(0);
		addHeader((char*)"Clonal Fold Change");
		PeakLibrary* clonalFiltered = filteredPeaks->filterClonalPeaks(tags,
					peakSize,clonalFold,diffMode, strand);
		delete filteredPeaks;
		filteredPeaks = clonalFiltered;
		clonalPeaks = before-filteredPeaks->numPeaks;
	}

	if (style == PEAK_STYLE_SUPERENHANCERS) {
		fprintf(stderr, "\tStitching together putative peaks into regions\n");
		PeakLibrary* regions = filteredPeaks->stitchRegions(minDist,stitchMode);
		delete filteredPeaks;
		filteredPeaks = regions;
	}


	//normalize tag counts
	tags->setMaxTBP(0);

	filteredPeaks->setPeakTagSizeFixed(0,0);
	//filteredPeaks->setPeakTagSizeRefPos(NULL_OFFSET,-1*halfPeakSize,halfPeakSize);	

	filteredPeaks->setDefaultPeakOrder();
	filteredPeaks->tagsInPeaks = 0.0;

	//int tagsIndex = filteredPeaks->addTagLibrary(tags);
	//Doubletable* expTags = filteredPeaks->countPeakTags(tagsIndex,0,0,strand,COUNT_MODE_TOTAL);
	Doubletable* expTags = filteredPeaks->countPeakTagsLowMemory(tags,strand,COUNT_MODE_TOTAL);
	double normFactor = normTotal/tags->totalTags;

	if (style == PEAK_STYLE_SUPERENHANCERS && input != NULL) {
		//int inputTagsIndex = filteredPeaks->addTagLibrary(input);
		//Doubletable* inputTags = filteredPeaks->countPeakTags(inputTagsIndex,0,0,strand,COUNT_MODE_TOTAL);
		Doubletable* inputTags = filteredPeaks->countPeakTagsLowMemory(input,strand,COUNT_MODE_TOTAL);
		double inputNormFactor = normTotal/input->totalTags;

		for (int i=0;i<filteredPeaks->numPeaks;i++) {
			filteredPeaks->peakOrder[i]->v = expTags->search(filteredPeaks->peakOrder[i]->name);
			filteredPeaks->tagsInPeaks += filteredPeaks->peakOrder[i]->v;
			filteredPeaks->peakOrder[i]->v *= normFactor;
			filteredPeaks->peakOrder[i]->v -= inputTags->search(filteredPeaks->peakOrder[i]->name)*inputNormFactor;
		}
		delete inputTags;
	} else {
		for (int i=0;i<filteredPeaks->numPeaks;i++) {
			filteredPeaks->peakOrder[i]->v = expTags->search(filteredPeaks->peakOrder[i]->name);
			filteredPeaks->tagsInPeaks += filteredPeaks->peakOrder[i]->v;
		}
		double normFactor = normTotal/tags->totalTags;
		filteredPeaks->normalizePeakScore(normFactor);
	}

	tagsInPeaks = filteredPeaks->tagsInPeaks;
	delete expTags;

	numPeaks = filteredPeaks->numPeaks;
	fprintf(stderr, "\tTotal Peaks identified = %d\n", numPeaks);
	filteredPeaks->setDefaultPeakOrder();
	
	return filteredPeaks;
}

void PeakLibrary::normalizePeakScore(float normFactor) {
	char** keys = peaks->keys();
	for (int i=0;i<peaks->total;i++) {
		Peak* p = (Peak*) peaks->search(keys[i]);
		p->v *= normFactor;
		delete [](keys[i]);
	}
	delete []keys;
}
void PeakLibrary::checkForOutOfBoundsCoordinates(TagLibrary* tags) {

	if (tags == NULL) return;	
	if (tags->chrs == NULL) return;	

	char** keys = chrs->keys();
	for (int i=0;i<chrs->total;i++) {
		ChrPeaks* cp = (ChrPeaks*)chrs->search(keys[i]);
		if (cp == NULL) continue;
		ChrTags* ct = (ChrTags*)tags->chrs->search(keys[i]);
		if (ct == NULL) continue;
		cp->checkForOutOfBoundsCoordinates(ct);

		delete [](keys[i]);	
	}
	delete []keys;
}

PeakLibrary* PeakFinder::filterPeaks(PeakLibrary* putativePeaks) {

	double *peakHeights = new double[fdrSize];
	for (int i=0;i<fdrSize;i++) {
		peakHeights[i] = 0;
	}

	char** keys  = putativePeaks->peaks->keys();
	for (int i=0;i<putativePeaks->peaks->total;i++) {	
		Peak* p = (Peak*) putativePeaks->peaks->search(keys[i]);
		int index = (int)p->v;
		if (index >= fdrSize) {
			index = fdrSize -1;
		}
		peakHeights[index]+=1.0;
	}
	fprintf(stderr, "\t\tThreshold\tPeak Count\tExpected Peak Count\tFDR\n");
	for (int i=fdrSize-2;i>=0;i--) {
		peakHeights[i] += peakHeights[i+1];
		double cfdr = 0.0;
		if (peakHeights[i] > 0) {
			cfdr = fdrTable[i]/peakHeights[i];
		}
		if (cfdr < fdr) {
			fdrThresh = (float)i;
		}
		if (cfdr > fdr/10.0 && i <= 50) {
			fprintf(stderr, "\t\t%d\t%.3lf\t%.3lf\t%lf\n", i,peakHeights[i],fdrTable[i],cfdr);
		}
	}
	if (filterMode == PEAKFINDER_FILTER_MODE_FDR) {
		poisson = poissonTable[(int)fdrThresh];
		fprintf(stderr, "\t%.2lf%% FDR Threshold set at %.1f (poisson pvalue ~ %.2le)\n", 
															100.0*fdr,fdrThresh,poisson);
	}
	
	PeakLibrary* goodPeaks = new PeakLibrary();	
	int numGood=0;
	char* strScore = new char[10000];
	tagsInPeaks = 0.0;

	double threshold2Use = fdrThresh;
	if (filterMode == PEAKFINDER_FILTER_MODE_POISSON) {
		threshold2Use = poissonThresh;
	} else if (filterMode == PEAKFINDER_FILTER_MODE_THRESH) {
		threshold2Use = tagThresh;
	}

	for (int i=0;i<putativePeaks->peaks->total;i++) {
		Peak* p = (Peak*) putativePeaks->peaks->search(keys[i]);
		if (p->v >= threshold2Use-000000.1) {
			tagsInPeaks += (double)p->v;
			sprintf(strScore,"%f",p->v);
			p->addData(strScore);
			goodPeaks->addPeak(p);
			numGood++;
		}
		delete [](keys[i]);
	}
	delete []keys;
	delete []peakHeights;
	delete []strScore;
	fprintf(stderr, "\t%d peaks passed threshold\n", numGood);
	goodPeaks->sortChr();

	return goodPeaks;
}


void PeakFinder::approxFdrTable() {

	tbp = tagsUsedForClustering/((double)gsize);

	//tpp = tags per peak	
	tpp = tbp*((double)peakSize);
	if (strand != STRAND_BOTH) {
		tpp = tbp*((double)peakSize/2.0);
	}
	fprintf(stderr, "\tExpected tags per peak = %lf (tbp = %lf)\n",tpp,tbp);	

	//this is the toughest part - and is a little empirical based on randomized tag positions
	double numTests = 2.0* gsize / ((double)peakSize);
	if (regionFlag) {
		numTests = 2.0* gsize / ((double)(peakSize/PEAKFRACTION_REGION_MINDIST/2.0));
	} else {
		numTests = 2.0* gsize / ((double)(minDist/2.0));
	}
	if (strand == STRAND_SEPARATE) {
		numTests *= 2;
		fprintf(stderr, "\tFinding tags on separate strands: doubling effective genome size\n");
	}


	fdrSize = PEAKFINDER_FDRSIZE;
	fdrTable = new double[fdrSize];
	poissonTable = new double[fdrSize];

	double cum = 0;
	for (int i=0;i<fdrSize;i++) {
		double lp = logPoisson(i,tpp);
		lp = exp(lp);
		//fprintf(stderr, "\t%d\t%.9lf\t%.9lf\t%.3f\n",i, lp, 1-cum,numTests*(1-cum));
		double poissonCum = 1-cum;
		fdrTable[i]=numTests*poissonCum;
		poissonTable[i]=poissonCum;
		if (poissonThresh < 0.01 && poissonCum < poisson) {
			poissonThresh = (float)i;
			if (filterMode == PEAKFINDER_FILTER_MODE_POISSON) {
				fprintf(stderr, "\tPoisson Threshold set at %d tags\n", i);
			}
		}
		cum += lp;
	}

}


// class PeakLibrary --------------------------------------------------------------------

PeakLibrary::PeakLibrary() {
	initialize(PEAK_LIBRARY_DEFAULT_SIZE);
}
PeakLibrary::PeakLibrary(int expectedNumPeaks) {
	initialize(expectedNumPeaks);
}
PeakLibrary::PeakLibrary(char* fname, int mode) {
	if (fname == NULL) {
		initialize(PEAK_LIBRARY_DEFAULT_SIZE);
		return;
	}
	FILE* fp = fopen(fname,"r");
	if (fp == NULL) {
		fprintf(stderr, "!!! Could not open peak file: %s !!!\n", fname); 
		initialize(PEAK_LIBRARY_DEFAULT_SIZE);
		return;
	}
	int numLines = 0;
	char* buf = new char[BUFFER];
	while (fgets(buf, BUFFER, fp) != NULL) {
		numLines++;
	}
	fclose(fp);
	delete []buf;
	initialize(numLines+1000);
	readPeakFile(fname,mode);
}
void PeakLibrary::initialize(int expectedNumPeaks) {
	chrs = new Hashtable(10000);
	peaks = new Hashtable(expectedNumPeaks);
	duplicates = new Inttable((int)(expectedNumPeaks/5));
	exps = NULL;
	name = NULL;
	genome = NULL;
	numPeaks = 0;
	numExps = 0;
	avgPeakSize = 0.0;
	tagsInPeaks = 0.0;
	fixedFlag = 0;
	peakOrder = NULL;
	chrSizes = NULL;
	duplicateWarningFlag = 1;
}
PeakLibrary* PeakLibrary::copy() {
	PeakLibrary* cp = new PeakLibrary(numPeaks);
	if (peakOrder == NULL) setDefaultPeakOrder();
	for (int i=0;i<numPeaks;i++) {
		cp->addPeak(peakOrder[i]);
	}
	return cp;
}
PeakLibrary::~PeakLibrary() {
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		ChrPeaks* ct = (ChrPeaks*) chrs->search(keys[i]);
		delete ct;
		delete [](keys[i]);
	}
	delete []keys;
	delete chrs;
	keys = peaks->keys();
	for (int i=0;i<peaks->total;i++) {
		Peak* p = (Peak*) peaks->search(keys[i]);
		delete p;
		delete [](keys[i]);
	}
	delete []keys;
	delete peaks;
	delete duplicates;
	if (name != NULL) delete []name;
	if (genome != NULL) delete []genome;
	if (exps != NULL) {
		delete []exps;
	}
	if (peakOrder != NULL) delete []peakOrder;
	
}

PeakLibrary* PeakLibrary::prioritizeAnnotations() {


	char* annname = new char[1000];	
	char** chr = chrs->keys();
	qsort(chr,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		//add peak corresponding to leftover space
		strcpy(annname,"Intergenic");
		addPeak((char*)"Intergenic",chr[i],1,1500000000,0,0,0,0.0,(char*)"N",-1,1900000000);
		//strcpy(annname,"IntergenicEnd");
		//addPeak((char*)"Intergenic",chr[i],1500000001,1500000009,0,0,0,0.0,(char*)"N",-1,1950000000);
	}
	sortChr();
	
	PeakLibrary* annotations = new PeakLibrary(10000000);
	annotations->duplicateWarningFlag = 0;

	fprintf(stderr, "\tPrioritizing Annotations: ");
	for (int i=0;i<chrs->total;i++) {
		ChrPeaks* cp = (ChrPeaks*) chrs->search(chr[i]);
		//fprintf(stderr, "%s,", chr[i]);
		fprintf(stderr, "."); //, chr[i]);
		cp->prioritizeAnnotations(annotations);
		delete [](chr[i]);
	}
	delete []chr;
	fprintf(stderr, "\n");

	annotations->sortChr();
	//annotations->printSorted(stdout);
	//exit(0);

	return annotations;
}
int PeakLibrary::getAveragePeakSize() {

	double total = 0;
	double N = 0;
	char** keys = peaks->keys();
	for (int i=0;i<numPeaks;i++) {
		Peak* p = (Peak*) peaks->search(keys[i]);
		int d = p->end - p->start;
		total += (double)d;
		N += 1.0;
		delete [](keys[i]);
	}
	//fprintf(stderr, "total=%lf\tN=%lf\n", total, N);
	delete []keys;
	if (N > 0.5) {
		total /= N;
	}
	return (int)(total+0.5);
}
void PeakLibrary::processBedGraphFile(char* bedGraphFile, int relativeStart, int relativeEnd,
											char strand, FILE* outputFile, int outputMode) {
	char* buf = new char[BUFFER];
	char** line = new char*[BUFFER];
	int numCols = 0;
	char* lastChr = new char[BUFFER];
	lastChr[0]='\0';

	Hashtable* linkedTags = new Hashtable(numPeaks*2);
	Hashtable* linkedFirstTags = new Hashtable(numPeaks*2);
	ChrPeaks* cp = NULL;
	int cpIndex = 0;

	int start = 0;
	int end = 0;
	float v = 0;

	FILE* fp = fopen(bedGraphFile, "r");
	while (fgets(buf, BUFFER, fp) != NULL) {
		split(buf, line, numCols,WHITE_SPACE);
		if (numCols < 4) continue;
		if (line[0][0] == '#') continue;
		if (strcmp(line[0],"track")==0) continue;
		
		if (strcmp(line[0],lastChr) != 0) {
			//new chromosome
			cp = (ChrPeaks*) chrs->search(line[0]);
			cpIndex=0;

			strcpy(lastChr,line[0]);
		}
		if (cp == NULL) continue;
		sscanf(line[1],"%d",&start);
		sscanf(line[2],"%d",&end);
		sscanf(line[3],"%f",&v);

		for (int i=cpIndex;i<cp->numPeaks;i++) {
			int pstart = cp->peaks[i]->refPos+relativeStart;
			int pend = cp->peaks[i]->refPos+relativeEnd;
			//int refPos = cp->peaks[i]->refPos;
			if (cp->peaks[i]->strand == STRAND_NEGATIVE) {
				pstart = cp->peaks[i]->refPos-relativeEnd;
				pend = cp->peaks[i]->refPos-relativeStart;
			}
			if (cp->peaks[i]->fixedFlag) {
				pstart = cp->peaks[i]->start+relativeStart;
				pend = cp->peaks[i]->end+relativeEnd;
				//refPos = cp->peaks[i]->start;
				if (cp->peaks[i]->strand == STRAND_NEGATIVE) {
					pstart = cp->peaks[i]->start-relativeEnd;
					pend = cp->peaks[i]->end-relativeStart;
					//refPos = cp->peaks[i]->end;
				}
			}
			if (pstart > end) {
				break;
			}
			if (i == cpIndex && pend < start) {
				cpIndex++;
			}

			if (end > pstart && start < pend) {
				int tstart = start;
				if (tstart < pstart) tstart = pstart;
				int tend = end;
				if (tend > pend) tend = pend;
				int length = tend-tstart;
	
				int peakPos = cp->peaks[i]->refPos;
				if (cp->peaks[i]->fixedFlag) {
					peakPos = cp->peaks[i]->start;
					if (cp->peaks[i]->strand == STRAND_NEGATIVE) {
						peakPos = cp->peaks[i]->end;
					}
				}
				int pos = tstart-peakPos;
				if (cp->peaks[i]->strand == STRAND_NEGATIVE) {
					pos = peakPos-tend;
				}

				LinkedTag* link = NULL;
				link = new LinkedTag(pos,STRAND_POSITIVE,length,v,NULL);

				LinkedTag* firstLink = (LinkedTag*)linkedFirstTags->search(cp->peaks[i]->name);
				if (firstLink == NULL) linkedFirstTags->insert(link,cp->peaks[i]->name);

				LinkedTag* lastLink = (LinkedTag*)linkedTags->search(cp->peaks[i]->name);
				if (lastLink != NULL) {
					lastLink->tag = link;
				}
				linkedTags->insert(link,cp->peaks[i]->name);
			}
		}
	} 
	fclose(fp);

	char** keys = linkedFirstTags->keys();
	for (int i=0;i<linkedFirstTags->total;i++) {  
		LinkedTag* firstLink = (LinkedTag*) linkedFirstTags->search(keys[i]);
		if (firstLink == NULL) continue;
		fprintf(outputFile, "%s\t",keys[i]);
		int numLinks = 0;
		LinkedTag* curTag = firstLink;
		double total = 0.0;
		double totalBp = 0.0;
		while (curTag != NULL) {
			total += ((double)curTag->len)*curTag->v;
			totalBp += ((double)curTag->len);
			if (outputMode == OUTPUT_MODE_PEAKTAGS) {
				if (numLinks > 0) fprintf(outputFile, ",");
				fprintf(outputFile,"%d|%d=%.3lf",curTag->p,curTag->p+curTag->len,curTag->v);
			}
			numLinks++;
			LinkedTag* next = curTag->tag;
			delete curTag;
			curTag = next;
		}
		if (outputMode == OUTPUT_MODE_COUNT) {
			if (totalBp > 0.0) total /= totalBp;
			fprintf(outputFile, "%.3lf",total);
		}
		fprintf(outputFile, "\n");
		delete [](keys[i]);
	}
	delete []keys;
	delete linkedFirstTags;
	delete linkedTags;

	delete []buf;
	delete []line;

}
void PeakLibrary::processWiggleFile(char* wiggleFile, int relativeStart, int relativeEnd,
											char strand, FILE* outputFile, int outputMode) {
	char* buf = new char[BUFFER];
	char** line = new char*[BUFFER];
	int numCols = 0;
	char* lastChr = new char[BUFFER];
	lastChr[0]='\0';

	Hashtable* linkedTags = new Hashtable(numPeaks*2);
	Hashtable* linkedFirstTags = new Hashtable(numPeaks*2);
	ChrPeaks* cp = NULL;
	int cpIndex = 0;

	int start = 0;
	int end = 0;
	float v = 0;


	int step = 1;
	int span = 1;
	int position = 0;

	int fixedFlag=0;

	FILE* fp = fopen(wiggleFile, "r");
	while (fgets(buf, BUFFER, fp) != NULL) {
		split(buf, line, numCols, WHITE_SPACE);
		if (line[0][0] == '#') continue;
		if (strcmp(line[0],"track")==0) continue;
		if (strcmp(line[0],"variableStep")==0) {
			span = 1;
			fixedFlag=0;
			for (int i=1;i<numCols;i++) {
				if (strlen(line[i]) > 6) {
					if (strncmp(line[i],"chrom=",6)==0) {
						char* c = &(line[i][6]);
						if (strcmp(c,lastChr) != 0) {
							//new chromosome
							cp = (ChrPeaks*) chrs->search(c);
							cpIndex=0;
							strcpy(lastChr,c);
						}
					}
				} 
				if (strlen(line[i]) > 5) {
					if (strncmp(line[i],"span=",5)==0) {
						char* s = &(line[i][5]);
						sscanf(s,"%d",&span);
					}
				}
			}
			continue;
		}
		if (strcmp(line[0],"fixedStep")==0) {
			span = 1;
			step = 1;
			fixedFlag=1;
			for (int i=1;i<numCols;i++) {
				if (strlen(line[i]) > 6) {
					if (strncmp(line[i],"chrom=",6)==0) {
						char* c = &(line[i][6]);
						if (strcmp(c,lastChr) != 0) {
							//new chromosome
							cp = (ChrPeaks*) chrs->search(c);
							cpIndex=0;
							strcpy(lastChr,c);
						}
					}
				} 
				if (strlen(line[i]) > 5) {
					if (strncmp(line[i],"span=",5)==0) {
						char* s = &(line[i][5]);
						sscanf(s,"%d",&span);
					}
				}
				if (strlen(line[i]) > 5) {
					if (strncmp(line[i],"step=",5)==0) {
						char* s = &(line[i][5]);
						sscanf(s,"%d",&step);
					}
				}
				if (strlen(line[i]) > 6) {
					if (strncmp(line[i],"start=",6)==0) {
						char* s = &(line[i][6]);
						sscanf(s,"%d",&position);
						start = position;
					}
				}
			}
			continue;
		}
		if (cp == NULL) continue;
					
		v= 0.0;
		if (fixedFlag==0) {
			if (numCols < 2) continue;
			sscanf(line[0],"%d",&start);
			sscanf(line[1],"%f",&v);
		} else if (fixedFlag==1) {
			if (numCols < 1) continue;
			sscanf(line[0],"%f",&v);
		}
	
		end = start+span;
	
		for (int i=cpIndex;i<cp->numPeaks;i++) {
			int pstart = cp->peaks[i]->refPos+relativeStart;
			int pend = cp->peaks[i]->refPos+relativeEnd;
			//int refPos = cp->peaks[i]->refPos;
			if (cp->peaks[i]->strand == STRAND_NEGATIVE) {
				pstart = cp->peaks[i]->refPos-relativeEnd;
				pend = cp->peaks[i]->refPos-relativeStart;
			}
			if (cp->peaks[i]->fixedFlag) {
				pstart = cp->peaks[i]->start+relativeStart;
				pend = cp->peaks[i]->end+relativeEnd;
				//refPos = cp->peaks[i]->start;
				if (cp->peaks[i]->strand == STRAND_NEGATIVE) {
					pstart = cp->peaks[i]->start-relativeEnd;
					pend = cp->peaks[i]->end-relativeStart;
					//refPos = cp->peaks[i]->end;
				}
			}
			if (pstart > end) {
				break;
			}
			if (i == cpIndex && pend < start) {
				cpIndex++;
			}

			if (end > pstart && start < pend) {
				//fprintf(stderr, "\t%d\t%d\t%f\n", start,end,v);
				int tstart = start;
				if (tstart < pstart) tstart = pstart;
				int tend = end;
				if (tend > pend) tend = pend;
				int length = tend-tstart;
	
				int peakPos = cp->peaks[i]->refPos;
				if (cp->peaks[i]->fixedFlag) {
					peakPos = cp->peaks[i]->start;
					if (cp->peaks[i]->strand == STRAND_NEGATIVE) {
						peakPos = cp->peaks[i]->end;
					}
				}
				int pos = tstart-peakPos;
				if (cp->peaks[i]->strand == STRAND_NEGATIVE) {
					pos = peakPos-tend;
				}

				LinkedTag* link = NULL;
				link = new LinkedTag(pos,STRAND_POSITIVE,length,v,NULL);

				LinkedTag* firstLink = (LinkedTag*)linkedFirstTags->search(cp->peaks[i]->name);
				if (firstLink == NULL) linkedFirstTags->insert(link,cp->peaks[i]->name);

				LinkedTag* lastLink = (LinkedTag*)linkedTags->search(cp->peaks[i]->name);
				if (lastLink != NULL) {
					lastLink->tag = link;
				}
				linkedTags->insert(link,cp->peaks[i]->name);
			}
		}
		if (fixedFlag == 1) {
			start += step;
		}
	} 
	fclose(fp);

	char** keys = linkedFirstTags->keys();
	for (int i=0;i<linkedFirstTags->total;i++) {  
		LinkedTag* firstLink = (LinkedTag*) linkedFirstTags->search(keys[i]);
		if (firstLink == NULL) continue;
		fprintf(outputFile, "%s\t",keys[i]);
		int numLinks = 0;
		LinkedTag* curTag = firstLink;
		double total = 0.0;
		double totalBp = 0.0;
		while (curTag != NULL) {
			total += ((double)curTag->len)*curTag->v;
			totalBp += ((double)curTag->len);
			if (outputMode == OUTPUT_MODE_PEAKTAGS) {
				if (numLinks > 0) fprintf(outputFile, ",");
				fprintf(outputFile,"%d|%d=%.3lf",curTag->p,curTag->p+curTag->len,curTag->v);
			}
			numLinks++;
			LinkedTag* next = curTag->tag;
			delete curTag;
			curTag = next;
		}
		if (outputMode == OUTPUT_MODE_COUNT) {
			if (totalBp > 0.0) total /= totalBp;
			fprintf(outputFile, "%.3lf",total);
		}
		fprintf(outputFile, "\n");
		delete [](keys[i]);
	}
	delete []keys;
	delete linkedFirstTags;
	delete linkedTags;

	delete []buf;
	delete []line;

}

void PeakLibrary::readPeakFile(char* filename, int mode) {

	FILE* fp = fopen(filename,"r");
	if (fp == NULL) {
		fprintf(stderr, "Could not open peak file (%s)\n", filename);
		return;
	}

	char* buf = new char[BUFFER];
	char* countingName = new char[BUFFER];
	char** line = new char*[BUFFER];
	char* chr = NULL;
	char* name = NULL;
	char* ann = NULL;
	int numCols = 0;


	avgPeakSize = 0.0;
	int NN = 0;
	int start = 0;
	int end = 0;
	int midpoint = NULL_REF;
	float tagCount = 0;
	float focusRatio = 0;
	int dir = 0;
	unsigned int currentPriority = 0;
	unsigned int priority = 0;
	unsigned int intergenicPriority = 1000000000;
	char defaultPrefix[8] = "default";
	char defaultName[50] = "";
	int defaultCount = 0;
	int mappability = -1;

	while (fgets(buf, BUFFER, fp) != NULL) {
		split(buf, line, numCols,'\t');
		if (numCols < 3) continue;
		if (line[0][0] == '#') continue;
	
		//check if there is a header (i.e. 3rd column should be integer
		if (line[2][0] != '-' && (line[2][0] < 48 || line[2][0] > 57)) {
			continue;
		}

		int format = FORMAT_UNKNOWN;

		name = NULL;
		chr = NULL;
		start = -1;
		end = -1;
		dir = -1;
		

		int col1notNum = 0;
		int col2notNum = 0;
		int col3notNum = 0;
		col1notNum = checkInt(line[1]);
		col2notNum = checkInt(line[2]);
		if (col2notNum) {
			continue;
		}
		int col3strand = -1;
		int col4strand = -1;
		int col5strand = -1;

		if (numCols == 3) {	
			if (col1notNum) {
				continue;
			}
			format = FORMAT_BED;
		} else {
			col3notNum = checkInt(line[3]);
			if (col1notNum && col3notNum) {
				continue;
			}
			if (numCols == 4) {
				if (col3notNum) {
					format = FORMAT_BED;
					col3strand = checkStrand(line[3]);
					if (col3strand == -1) {
						name = line[3];
					} else {
						dir = col3strand;
					}
				} else {
					format = FORMAT_PEAK;
				}
			} else {
				col4strand = checkStrand(line[4]);
				col3strand = checkStrand(line[3]);
				if (col4strand != -1 && col3notNum==0) {
					format = FORMAT_PEAK;
					dir = col4strand;
				} else if (col3notNum) {
					format = FORMAT_BED;
					if (col3strand != -1) {
						dir = col3strand;
					} else {
						name = line[3];
					}
					if (numCols > 5) {
						col5strand = checkStrand(line[5]);
						if (col5strand != -1) {
							dir = col5strand;
						}
					}
				} else {
					continue;
				}
			}
		}

		if (format == FORMAT_BED) {
			chr = line[0];
			sscanf(line[1],"%d",&start);
			sscanf(line[2],"%d",&end);
			start++;
		} else if (format == FORMAT_PEAK) {
			name = line[0];
			chr = line[1];
			sscanf(line[2],"%d",&start);
			sscanf(line[3],"%d",&end);
		} else {
			continue;
		}

		if (dir == -1) {
			dir = 0;
		}
		if (name == NULL) {
			defaultCount++;
			sprintf(defaultName,"%s-%d", defaultPrefix, defaultCount);
			name = defaultName;
		}
			
		//initialize optional items
		midpoint = NULL_REF;
		tagCount = 0.0;
		focusRatio = 0.0;
		ann = NULL;

		avgPeakSize += (end-start);
		NN++;

		if (mode == PEAK_READ_MODE_NORMAL && format== FORMAT_PEAK) {
			if (numCols > 5) sscanf(line[5],"%f",&tagCount);
			if (numCols > 6) sscanf(line[6],"%f",&focusRatio);
		} else if (mode == PEAK_READ_MODE_ANNOTATION && format == FORMAT_PEAK) {
			if (numCols > 5) {
				ann = line[5];
				if (strncmp(ann,"I",1)==0) {
					currentPriority = intergenicPriority++;
				} else {
					currentPriority = priority++;
				}
			} else {
				currentPriority = priority++;
			}
		} else if (mode >= PEAK_READ_MODE_COUNTING) {
			sprintf(countingName, "%d-%d",mode,NN);
			name = countingName;
		}

		addPeak(name,chr,start,end,midpoint,(char)dir,tagCount,focusRatio,ann,mappability,priority);
		//if (numPeaks % 1000000==0) fprintf(stderr, "\t\tloaded %d peaks\n", numPeaks);
	}
	fclose(fp);

	avgPeakSize /= (double)NN;

	delete []buf;
	delete []line;
	delete []countingName;
	sortChr();
}

void PeakLibrary::sortChr() {

	numPeaks = 0;

	char** chrKeys = chrs->keys();
	qsort(chrKeys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		ChrPeaks* ct = (ChrPeaks*)chrs->search(chrKeys[i]);
		ct->sort();
		numPeaks += ct->numPeaks;
		delete [](chrKeys[i]);
	}
	delete []chrKeys;
}
void PeakLibrary::setMaxChrPositions(Doubletable* nchrSizes) {
	if (nchrSizes != NULL) {
		chrSizes = nchrSizes;
	}

	if (chrSizes != NULL) {
		char** chr = chrs->keys();
		for (int i=0;i<chrs->total;i++) {
			int max = (int)chrSizes->search(chr[i]);
			if (max == EMPTY_INT) max = 0;
			ChrPeaks* cp = (ChrPeaks*)chrs->search(chr[i]);
			for (int j=0;j<cp->numPeaks;j++) {
				if (cp->peaks[j]->start < 0) cp->peaks[j]->start = 0;
				if (cp->peaks[j]->end < 0) cp->peaks[j]->end = 1;
				if (cp->peaks[j]->start > max) cp->peaks[j]->start = max;
				if (cp->peaks[j]->end > max) cp->peaks[j]->end = max;
			}
			delete [](chr[i]);
		}
		delete []chr;
	}
}
PeakLibrary* PeakLibrary::getCoveragePeaks(int resolution,int superRes) {

	PeakLibrary* cpeaks = new PeakLibrary();

	char** keys = chrs->keys();
	char* pname = new char[500];

	if (superRes < 0) superRes = resolution;
	int halfSuperRes = superRes/2;
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		ChrPeaks* cp = (ChrPeaks*)chrs->search(keys[i]);
		for (int j=0;j<cp->numPeaks;j++) {
			int cstart = cp->peaks[j]->start;
			int cend = cp->peaks[j]->end;
			if (cstart < 0) cstart = 0;

			for (int j=cstart;j<cend;j+=resolution) {
				sprintf(pname,"%s-%d",keys[i],j); 
				int midpoint = j+resolution/2;
				int s = midpoint - halfSuperRes;
				int e = midpoint + halfSuperRes;
				cpeaks->addPeak(pname,keys[i],s,e,midpoint,STRAND_POSITIVE,
								0.0, 0.0, NULL, 0,0);
			}
		}

		delete [](keys[i]);
	}
	delete []keys;
	delete []pname;
	cpeaks->sortChr();
	return cpeaks;
}

Doubletable* PeakLibrary::reportTotalValuesPerChr() {

	Doubletable* chrTotals = new Doubletable(1000);
	char** keys = peaks->keys();
	double gTotal = 0.0;
	for (int i=0;i<peaks->total;i++) {
		Peak* p = (Peak*)peaks->search(keys[i]);
		double curTotal = chrTotals->search(p->chr);
		if (curTotal < EMPTY_DOUBLE_CHECK) {
			curTotal = 0.0;
		}
		curTotal += p->v;
		gTotal += p->v;
		chrTotals->insert(curTotal, p->chr);

		delete [](keys[i]);
	}
	delete []keys;
	chrTotals->insert(gTotal, (char*)"genome");

	if (1) {
		char** keys = chrTotals->keys();
		qsort(keys,chrTotals->total,sizeof(char*),&chrcmp);
		fprintf(stderr, "\tChr\tTotal\n");
		for (int i=0;i<chrTotals->total;i++) {
			fprintf(stderr, "\t%s\t%lf\n", keys[i],chrTotals->search(keys[i]));
			delete [](keys[i]);
		}
		delete []keys;
	}

	return chrTotals;
}
	
Peak* PeakLibrary::checkForPeak(char* name, char* chr, int start, int end, char strand, int &info) {
	info = PEAK_LIBRARY_NOVEL_PEAK;
	Peak* p = NULL;
	if (name != NULL) {
		p = (Peak*) peaks->search(name);
		if (p != NULL) {
			info = PEAK_LIBRARY_SAME_NAME;
			return p;
		}
	}
	if (chr != NULL) {
		ChrPeaks* cp = (ChrPeaks*) chrs->search(chr);
		if (cp != NULL) {
			for (int i=0;i<cp->numPeaks;i++) {
				if (cp->peaks[i]->start == start && cp->peaks[i]->end==end && cp->peaks[i]->strand == strand) {
					p = cp->peaks[i];
					info = PEAK_LIBRARY_SAME_POSITION;
					return p;
				}
			}	
		}
	}
	return p;
}
Peak* PeakLibrary::addPeak(char* name, char* chr,int start, int end, int midpoint, char dir, 
					float value, float ratio, char* extraData,int mappability, unsigned int priority) {
	Peak* oldPeak = (Peak*)peaks->search(name);
	Peak* p=NULL;
	static int warningIssued = 0;
	if (oldPeak != NULL) {
		int dupCount = duplicates->search(name);
		if (dupCount == EMPTY_INT) {
			dupCount = 1;
		}
		dupCount++;
		//fprintf(stderr, "dupCount=%d \t%s\n",dupCount,name);
		duplicates->insert(dupCount, name);

		char* newname = new char[strlen(name)+20];
		sprintf(newname,"%s--%d", name,dupCount);
		if (priority < 1) {
			if (warningIssued == 0 && duplicateWarningFlag) {
				fprintf(stderr,"!!! Duplicate peak name (%s) !!! This will can/will cause problems !!!\n", name);
				fprintf(stderr,"\tNew name for this peak is %s\n", newname);
			} else if (warningIssued % 1000 == 0 && duplicateWarningFlag) {
				fprintf(stderr, "!!! warning over %d peaks with duplicate names\n", warningIssued);
			}
			warningIssued++;
		}
		oldPeak = (Peak*)peaks->search(newname);
		if (oldPeak != NULL && duplicateWarningFlag) {
			fprintf(stderr, "There's a problem!!! %s - %d %s\n", newname,dupCount,name);
		}
		p = new Peak(newname,name,chr,start,end,midpoint,dir,value,ratio,extraData,mappability,priority);
		delete []newname;
	} else {
		p = new Peak(name,NULL,chr,start,end,midpoint,dir,value,ratio,extraData,mappability,priority);
	}
	peaks->insert(p, p->name);
	ChrPeaks* cp = (ChrPeaks*)chrs->search(chr);
	if (cp == NULL) {
		cp = new ChrPeaks();
		chrs->insert(cp,chr);
	}
	cp->addPeak(p);
	numPeaks++;
	return p;
}
void PeakLibrary::addPeak(Peak* p) {
	char* n = p->name;
	if (p->ogname != NULL) {
		n = p->ogname;
	}
	addPeak(n,p->chr,p->start,p->end,p->refPos,p->strand,p->v,p->focusRatio,p->data,
											p->uniqMap,p->priority);
}
void PeakLibrary::print(FILE* fp) {
	char** keys = peaks->keys();
	sortKeys(keys);
	for (int i=0;i<peaks->total;i++) {
		Peak* p = (Peak*) peaks->search(keys[i]);
		p->print(fp);
		delete [](keys[i]);
	}	
	delete []keys;
}
void PeakLibrary::printSorted(FILE* fp) {
	sortChr();
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		ChrPeaks* cp = (ChrPeaks*) chrs->search(keys[i]);
		cp->print(fp);
		delete [](keys[i]);
	}	
	delete []keys;
}
void PeakLibrary::printGTF(FILE* fp) {
	sortChr();
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		ChrPeaks* cp = (ChrPeaks*) chrs->search(keys[i]);
		for (int j=0;j<cp->numPeaks;j++) {
			cp->peaks[j]->printGTF(fp);
		}
		delete [](keys[i]);
	}	
	delete []keys;
}
void PeakLibrary::printAnnotation(FILE* fp) {
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		ChrPeaks* cp = (ChrPeaks*) chrs->search(keys[i]);
		for (int j=0;j<cp->numPeaks;j++) {
			cp->peaks[j]->printAnnotation(fp);
		}
		delete [](keys[i]);
	}	
	delete []keys;
}
Hashtable* sortingPeaks=NULL;
int sortPeakNames(const void* a, const void* b) {
	Peak* p1 = (Peak*)sortingPeaks->search(*((char**)a));
	Peak* p2 = (Peak*)sortingPeaks->search(*((char**)b));
	float v1 = p1->v;
	float v2 = p2->v;
	if (v1 < v2) return 1;
	if (v2 < v1) return -1;
	char* c1 = p1->chr;
	char* c2 = p2->chr;
	int c = chrcmp((void*)&c1,(void*)&c2);
	if (c != 0) return c;
	int pos1 = p1->start;
	int pos2 = p2->start;
	if (pos1 > pos2) return 1;
	if (pos1 < pos2) return -1;
	return 0;
}
void PeakLibrary::sortKeys(char** keys) {
	sortingPeaks = peaks;
	qsort(keys, numPeaks,sizeof(char*),&sortPeakNames);
}

PeakLibrary* PeakLibrary::getSuperEnhancers(double slopeThresh, int window, char* &notes, PeakLibrary* &typical) {

	PeakLibrary* npeaks = new PeakLibrary();
	typical = new PeakLibrary();
	if (numPeaks < 1) { 
		fprintf(stderr, "!!Warning: zero peaks to find superenhancers from...\n");
		return npeaks;
	}

	int halfWindow = window/2;


	char** keys = peaks->keys();
	sortKeys(keys);
	double stepSize = ((double)window)/(double)(peaks->total);
	Peak* firstPeak = (Peak*)peaks->search(keys[0]);
	double max = (double)firstPeak->v;

	int numSuperEnhancers = -1;

	for (int i=0;i<peaks->total;i++) {
		if (i < halfWindow || i >= peaks->total-halfWindow) {
			if (numSuperEnhancers == -1) {
				npeaks->addPeak((Peak*)peaks->search(keys[i]));
			} else {
				typical->addPeak((Peak*)peaks->search(keys[i]));
			}
			continue;
		}

		Peak* p1 = (Peak*) peaks->search(keys[i-halfWindow]);
		Peak* p2 = (Peak*) peaks->search(keys[i+halfWindow]);
		double slope = (p1->v-p2->v)/max/stepSize;
		//if (i % window == 0) fprintf(stderr, "\tNumber of Enhancers: %d\tSlope: %.2lf\n", i, slope);
		if (slope < slopeThresh && numSuperEnhancers == -1) {
			numSuperEnhancers = i;
		}
		Peak* p = (Peak*)peaks->search(keys[i]);
		p->focusRatio = slope;

		if (numSuperEnhancers == -1) {
			npeaks->addPeak((Peak*)peaks->search(keys[i]));
		} else {
			typical->addPeak((Peak*)peaks->search(keys[i]));
		}

		if (i == halfWindow) {
			for (int j=0;j<halfWindow;j++) {
				Peak* pp = (Peak*)peaks->search(keys[j]);
				pp->focusRatio = slope;
			}
		}
		if (i == peaks->total-halfWindow-1) {
			for (int j=i+1;j<peaks->total;j++) {
				Peak* pp = (Peak*)peaks->search(keys[j]);
				pp->focusRatio = slope;
			}
		}
	}	
	if (numSuperEnhancers == -1) {
		numSuperEnhancers = npeaks->numPeaks;
	}
	for (int i=0;i<peaks->total;i++) {
		delete [](keys[i]);
	}
	delete []keys;

	notes = new char[2000];
	fprintf(stderr, "\t%d super enhancers found (of %d total)\n", numSuperEnhancers, peaks->total);
	sprintf(notes, "# Super Enhancer classification is in the focus ratio/length column\n# This value is the slope on the super-enhancer plot (window size of %d).\n# Given a slope cutoff of 1, the top %d peaks in this file are Super Enhancers\n", window, numSuperEnhancers);

	return npeaks;
}


SNP* PeakLibrary::addSNPsfromVCF(char* vcfFile, char** individuals, int &numIndividuals, 
									int allFlag) {
	

	if (vcfFile == NULL) return NULL;
	FILE* fp = fopen(vcfFile,"r");
	if (fp == NULL) {
		fprintf(stderr, "!!! Could not open vcf file %s !!!\n",vcfFile);
		return NULL;
	}


	int columnOffset = 9;
	SNP* baseSNP=NULL;
	SNP* lastSNP=NULL;

	char** individualNames = NULL;
	int* individualIndex = NULL;
	int numIndividualIndex = 0;

	ChrPeaks* cp = NULL;
	int index = 0;

	char* lastChr = new char[10000];
	char* curChr = new char[10000];

	char** altSeq = new char*[1000];
	int numAltSeq=0;
	char** data = new char*[1000];
	int numData=0;
	curChr[0] = '\0';
	lastChr[0] = '\0';

	char* buf = new char[BUFFER];
	char** line = new char*[BUFFER];
	int numCols = 0;

	while (fgets(buf, BUFFER, fp) != NULL) {
		split(buf,line,numCols,'\t');
		if (line[0][0] == '#') {
			if (strlen(line[0]) < 2) continue;
			if (line[0][1] == '#') continue;
			if (strncmp(line[0],"#CHR",4)==0) {
				
				if (numIndividuals == 0) {
					numIndividualIndex = numCols-columnOffset;
					individualNames = new char*[numIndividualIndex];
					individualIndex = new int[numIndividualIndex];
					for (int i=columnOffset;i<numCols;i++) {
						individualIndex[i-columnOffset] = i;
						individualNames[i-columnOffset] = new char[strlen(line[i])+1];
						individuals[i-columnOffset] = new char[strlen(line[i])+1];
						strcpy(individualNames[i-columnOffset],line[i]);
						strcpy(individuals[i-columnOffset],line[i]);
					}
					numIndividuals = numIndividualIndex;
				} else {
					numIndividualIndex = numIndividuals;
					individualNames = new char*[numIndividuals];
					individualIndex = new int[numIndividuals];
					for (int j=0;j<numIndividuals;j++) {
						individualNames[j] = new char[strlen(individuals[j])+1];
						strcpy(individualNames[j],individuals[j]);
						individualIndex[j]=-1;
					}
					for (int i=columnOffset;i<numCols;i++) {
						for (int j=0;j<numIndividuals;j++) {
							if (strcmp(line[i],individuals[j])==0) {
								individualIndex[j]=i;
							}
						}
					}
					for (int j=0;j<numIndividuals;j++) {
						if (individualIndex[j]==-1) {
							fprintf(stderr, "Warning: Didn't find individual %s in VCF file %s\n", 
											individuals[j],vcfFile);
							fprintf(stderr, "         Maybe you didn't spell it right??? (case-sensitive)\n");
						}
					}
				}
				for (int j=0;j<numIndividualIndex;j++) {
					//fprintf(stderr, "%s in %d\n", individualNames[j], individualIndex[j]); 
				}
								
			}
			continue;
		}
		if (line[0][0] != 'c') {
			sprintf(curChr,"chr%s",line[0]);
		} else {
			strcpy(curChr,line[0]);
		}
		int curPos = 0;
		sscanf(line[1],"%d",&curPos);
		
		char* refSeq = line[3];	
		split(line[4],altSeq,numAltSeq,',');

		if (strcmp(curChr,lastChr) != 0) {
			cp = (ChrPeaks*) chrs->search(curChr);
			index = 0;
			fprintf(stderr, "\t\t%s\n", curChr);
			strcpy(lastChr,curChr);
		}

		if (cp == NULL) continue;

		while (index < cp->numPeaks) {
			if (cp->peaks[index]->tagEnd < curPos) {
				index++;
			} else {
				break;
			}
		}
		if (index >= cp->numPeaks) continue;

		SNP* snp = NULL;

		for (int i=index;i<cp->numPeaks;i++) {
			if (curPos < cp->peaks[i]->tagStart) break;
			if (curPos >= cp->peaks[i]->tagStart &&
					curPos <= cp->peaks[i]->tagEnd) {

				
				if (snp == NULL) {
					//i.e. don't do this unless we absolutely have too
					int interesting = 1;
	
					snp = new SNP(refSeq, altSeq, numAltSeq, numIndividualIndex);

					
					for (int j=0;j<numIndividualIndex;j++) {
						int gindex1 = 0;
						int gindex2 = 0;
						int colIndex = individualIndex[j];
						if (colIndex < 0) {


						} else if (strcmp(line[colIndex],".")==0) {

						} else {
							split(line[colIndex],data,numData,':');
							if (strcmp(data[0],".")==0) {
						
							} else {
								split(data[0],data,numData,'/');
								if (numData < 2) {
									split(data[0],data,numData,'|');
								}
								if (numData < 2) {
									fprintf(stderr, "Having trouble parsing...\n");
								} else {
									if (strcmp(data[0],".")!=0) {
										sscanf(data[0],"%d",&gindex1);
									}
									if (strcmp(data[1],".")!=0) {
										sscanf(data[1],"%d",&gindex2);
									}
								}
							}
						}
						if (gindex1 > 0 || gindex2 > 0) {
							interesting = 1;
						}
						snp->allele1[j] = gindex1;
						snp->allele2[j] = gindex2;
					}
					if (interesting==0) {
						delete snp;
						snp = NULL;
						break;
					}
					if (baseSNP == NULL) {
						baseSNP = snp;
					} else {
						lastSNP->nextSNP = snp;
						snp->count = lastSNP->count+1;
					}
					lastSNP = snp;
				}
				cp->peaks[i]->addSNP(snp,curPos);
			}
		}
	}
	fclose(fp);

	delete []lastChr;
	delete []curChr;
	delete []altSeq;
	delete []line;
	delete []buf;

	return baseSNP;

}


Doubletable* PeakLibrary::countPeakTagsLowMemory(TagLibrary* tags, char direction,int mode) {

	Doubletable* results = new Doubletable();
	char** keys = peaks->keys();
	for (int i=0;i<peaks->total;i++) {
		results->insert(0.0,keys[i]);
		delete [](keys[i]);
	}
	delete []keys;

	keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	if (tags->singleFile) tags->readSingleTagFile();

	for (int i=0;i<chrs->total;i++) {
		ChrPeaks* cp = (ChrPeaks*) chrs->search(keys[i]);
		ChrTags* ct = (ChrTags*) tags->chrs->search(keys[i]);
		delete [](keys[i]);
		if (ct == NULL || cp == NULL) {
			continue;
		}
		cp->countPeakTagsLowMemory(results,ct,direction,mode);
	}
	delete []keys;

	return results;
}

int PeakLibrary::addTagLibrary(TagLibrary* t) {
	TagLibrary** newexps = new TagLibrary*[numExps+1];
	if (exps != NULL) {
		for (int i=0;i<numExps;i++) {
			newexps[i] = exps[i];
		}
		delete []exps;
	}
	newexps[numExps] = t;
	int rv = numExps;
	numExps++;
	exps = newexps;

	char** keys = peaks->keys();
	for (int i=0;i<peaks->total;i++) {
		Peak* p = (Peak*) peaks->search(keys[i]);
		p->addExp();  // this forces a blank space to be added (linkedlist)
		delete [](keys[i]);
	}
	delete []keys;
	keys = NULL;

	keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	if (t->singleFile) t->readSingleTagFile();
	for (int i=0;i<chrs->total;i++) {
		ChrPeaks* cp = (ChrPeaks*) chrs->search(keys[i]);
		ChrTags* ct = (ChrTags*) t->chrs->search(keys[i]);
		delete [](keys[i]);
		if (ct == NULL || cp == NULL) {
			continue;
		}
		cp->addTagLibrary(ct,numExps-1);
	}
	delete []keys;

	sortPeakTags(numExps-1); //here we optimize the linkedlist for access
	return rv;
}

void PeakLibrary::sortPeakTags(int expIndex) {
	char** keys = peaks->keys();
	for (int i=0;i<peaks->total;i++) {
		Peak* p = (Peak*)peaks->search(keys[i]);
		p->sortTags(expIndex);
		delete [](keys[i]);
	}
	delete []keys;
}

Doubletable* PeakLibrary::scoreNFR(TagLibrary* tags, int expIndex, int nfrSize, char direction) {
	Doubletable* results = new Doubletable();
	int nucSize = 150;
	char** keys = peaks->keys();
	for (int i=0;i<peaks->total;i++) {
		Peak* p = (Peak*) peaks->search(keys[i]);
		double value = p->scoreNFR(expIndex,tags->fragmentLengthEstimate,nfrSize,nucSize,direction);
		results->insert(value,keys[i]);
		delete [](keys[i]);
	}
	delete []keys;
	return results;
}

Doubletable* PeakLibrary::countPeakTags(int expIndex, int start, int end, char direction,int mode) {
	Doubletable* results = new Doubletable();
	char** keys = peaks->keys();
	for (int i=0;i<peaks->total;i++) {
		Peak* p = (Peak*) peaks->search(keys[i]);
		double value = p->countTags(expIndex,start,end, direction,mode);
		results->insert(value,keys[i]);
		delete [](keys[i]);
	}
	delete []keys;
	return results;
}

void PeakLibrary::printSNPtotals(FILE* fp) {
	char** keys = peaks->keys();
	for (int i=0;i<peaks->total;i++) {
		Peak* p = (Peak*) peaks->search(keys[i]);
		p->printSNPtotals(fp);;
		delete [](keys[i]);
	}
	delete []keys;
}
void PeakLibrary::printSNPs(FILE* fp) {
	char** keys = peaks->keys();
	for (int i=0;i<peaks->total;i++) {
		Peak* p = (Peak*) peaks->search(keys[i]);
		p->printSNPs(fp);;
		delete [](keys[i]);
	}
	delete []keys;
}
void PeakLibrary::printRelativeTags(FILE* fp, int expIndex, int start, int end,int mode) {
	char** keys = peaks->keys();
	for (int i=0;i<peaks->total;i++) {
		Peak* p = (Peak*) peaks->search(keys[i]);
		p->printRelativeTags(fp,expIndex,start,end,mode);
		delete [](keys[i]);
	}
	delete []keys;
}

void PeakLibrary::setPeakTagSizeRefPos(int newOffset, int startOffset, int endOffset) {
//fprintf(stderr, "%d\t%d\t%d\n", newOffset, startOffset, endOffset);
	char** keys = peaks->keys();
	for (int i=0;i<peaks->total;i++) {
		Peak* p = (Peak*) peaks->search(keys[i]);
		p->setPeakTagSizeRefPos(newOffset,startOffset,endOffset);
		delete [](keys[i]);
	}
	delete []keys;
	sortChr();
}
void PeakLibrary::setPeakTagSizeFixed(int startOffset, int endOffset) {
	fixedFlag = 1;
	char** keys = peaks->keys();
	for (int i=0;i<peaks->total;i++) {
		Peak* p = (Peak*) peaks->search(keys[i]);
		p->setPeakTagSizeFixed(startOffset,endOffset);
		delete [](keys[i]);
	}
	delete []keys;
	sortChr();
}

void PeakLibrary::addPeakLibrary(PeakLibrary* p) {
	char** keys = p->peaks->keys();
	for (int i=0;i<p->peaks->total;i++) {
		Peak* peak = (Peak*) p->peaks->search(keys[i]);
		addPeak(peak);
		delete [](keys[i]);
	}
	delete []keys;
	sortChr();
}

void PeakLibrary::setDefaultPeakOrder() {
	if (peakOrder != NULL) {
		delete []peakOrder;
	}
	if (numPeaks < 1) return;
	peakOrder = new Peak*[numPeaks];
	char** keys = peaks->keys();
	for (int i=0;i<peaks->total;i++) {
		peakOrder[i] = (Peak*) peaks->search(keys[i]);
		delete [](keys[i]);
	}
	delete []keys;

	qsort(peakOrder,numPeaks,sizeof(Peak*),&cmpPeaks);
	for (int i=0;i<numPeaks;i++) {
		peakOrder[i]->index = i;
	}

}


PeakLibrary* PeakLibrary::filterClonalPeaks(TagLibrary* tags, int peakSize,
						double threshold, int mode, char strand) {

	char* outputstr = new char[10000];
	int halfPeakSize = (peakSize)/2;

	float currentMaxTBP = tags->maxtbp;
	setPeakTagSizeRefPos(NULL_OFFSET,-1*halfPeakSize,halfPeakSize);	

	tags->setMaxTBP(0);
	//int tagsIndex = addTagLibrary(tags);
	Doubletable* expTags = countPeakTagsLowMemory(tags,strand,COUNT_MODE_TOTAL);
	//Doubletable* expTags = countPeakTags(tagsIndex,-1*halfPeakSize,halfPeakSize,strand,COUNT_MODE_TOTAL);

	tags->setMaxTBP(1);
	//int posIndex = addTagLibrary(tags);
	Doubletable* posTags = countPeakTagsLowMemory(tags,strand,COUNT_MODE_TOTAL);

	tags->setMaxTBP(currentMaxTBP);

	double avgTagsPerPosition = tags->averageTagsPerPosition;
	if (avgTagsPerPosition < 1) avgTagsPerPosition = 1.0;

	//Doubletable* expTags = countPeakTags(tagsIndex,-1*halfPeakSize,halfPeakSize,strand,COUNT_MODE_TOTAL);
	//Doubletable* posTags = countPeakTags(posIndex,-1*halfPeakSize,halfPeakSize,strand,COUNT_MODE_TOTAL);

	int expectedSize = peakSize*10;
	double *expected = new double[expectedSize];
	for (int i=0;i<expectedSize;i++) {
		if (i==0) {
			expected[i]=0;
		} else if (i==1) {
			expected[i] = 1;
		} else {
			expected[i] = expected[i-1]+((double)peakSize-expected[i-1])/((double)peakSize);
		}
	}

	PeakLibrary* goodPeaks = new PeakLibrary();
	char** keys = expTags->keys();
	int numGood = 0;
	int totalChecked = 0;
	goodPeaks->tagsInPeaks = 0.0;
	for (int i=0;i<peaks->total;i++) {
		double pt = expTags->search(keys[i]);
		double ct = posTags->search(keys[i]);
		if (pt < EMPTY_DOUBLE_CHECK || ct < EMPTY_DOUBLE_CHECK) {
			delete [](keys[i]);
			continue;
		}
		totalChecked++;
	
		int index = (int)(pt/avgTagsPerPosition);
		if (index>expectedSize-1) index = expectedSize-1;
		double et = expected[index];
		double fold = et/ct;
//fprintf(stderr, "%s\t%lf\t%lf\t%d\t%lf\n", keys[i], pt, ct, index,et);
		if (ct>0 && fold < threshold) {
			goodPeaks->tagsInPeaks += pt;
			numGood++;	
			Peak* p = (Peak*) peaks->search(keys[i]);
			sprintf(outputstr, "%.2lf",fold);
			p->addData(outputstr);
			goodPeaks->addPeak(p);
		}
		delete [](keys[i]);
		
	}
	delete []keys;
	delete []outputstr;
	delete posTags;
	delete expTags;
	if (totalChecked == 0) {
		fprintf(stderr, "\t!! Something is wrong - no peaks were checked!\n");
	} else {
		double ratio = ((double)numGood)/((double)totalChecked);
		fprintf(stderr, "\tClonal filtering: %d of %d (%.2lf%% passed)\n",numGood,totalChecked,ratio*100.0);
	}
	goodPeaks->sortChr();
	return goodPeaks;
}

PeakLibrary* PeakLibrary::filterLocalPeaks(TagLibrary* tags, int peakSize, int localSize,
						double foldThresh, double poissonThresh, int mode, char strand) {

	char* outputstr = new char[10000];

	int halfPeakSize = (peakSize)/2;
	int halfLocalSize = (localSize)/2;

	//int tagsIndex = addTagLibrary(tags);

	setPeakTagSizeRefPos(NULL_OFFSET,-1*halfPeakSize,halfPeakSize);
	Doubletable* expTags = countPeakTagsLowMemory(tags,strand,COUNT_MODE_TOTAL);
	setPeakTagSizeRefPos(NULL_OFFSET,-1*halfLocalSize,halfLocalSize);
	Doubletable* localTags = countPeakTagsLowMemory(tags,strand,COUNT_MODE_TOTAL);

	//Doubletable* expTags = countPeakTags(tagsIndex,-1*halfPeakSize,halfPeakSize,strand,COUNT_MODE_TOTAL);
	//Doubletable* localTags = countPeakTags(tagsIndex,-1*halfLocalSize,halfLocalSize,strand,COUNT_MODE_TOTAL);

	double peakLength = (double)peakSize;
	double localLength = (double)(localSize-peakSize);
	double localAdjustFactor = peakLength/localLength;

	double tagPseudoCount = 0.0;
	if (tags->tbp > 0.0) tagPseudoCount = 0.5;
	//if (tags->tbp > 0.0) tagPseudoCount = localLength*tags->tbp;
	
	PeakLibrary* goodPeaks = new PeakLibrary();
	char** keys = peaks->keys();
	int numGood = 0;
	goodPeaks->tagsInPeaks = 0.0;
	int totalChecked = 0;
	for (int i=0;i<peaks->total;i++) {
		double pt = expTags->search(keys[i]);
		double lt = localTags->search(keys[i]);
		if (pt < EMPTY_DOUBLE_CHECK || lt < EMPTY_DOUBLE_CHECK) {
			delete [](keys[i]);
			continue;
		}
		totalChecked++;
		lt -= pt;
		if (lt < tagPseudoCount) lt = tagPseudoCount;

		double ltadjusted = lt*localAdjustFactor;
		double fold = pt/ltadjusted;
		int foldGood = 0;

		if (foldThresh > 0.000001) {
			if (fold > foldThresh) {
				foldGood = 1;
			}
		} else {
			foldGood = 1;
		}

		int ptInt = (int) pt;
		//since the input is our "expected" meausrement, that becomes lambda
		//double cumPvalue = cumulativePoisson(ptInt, ltadjusted);
		//double pp = 1-cumPvalue;
		double pp = exp(ilogCumulativePoisson(ptInt, ltadjusted));
		//fprintf(stderr, "%lf\t%lf\n", pp, threshold);
		int poissonGood = 0;
		if (poissonThresh < 1.0) {
			if (pp < poissonThresh) {
				poissonGood = 1;
			}
		} else {
			poissonGood = 1;
		}

		if (foldGood == 1 && poissonGood == 1) {
			goodPeaks->tagsInPeaks += pt;
			numGood++;	
			Peak* p = (Peak*) peaks->search(keys[i]);
			sprintf(outputstr, "%.2lf\t%.2le",fold, pp);
			p->addData(outputstr);
			goodPeaks->addPeak(p);
		}
		delete [](keys[i]);
	}
	delete []keys;
	delete []outputstr;
	delete localTags;
	delete expTags;
	if (totalChecked == 0) {
		fprintf(stderr, "\t!! Something is wrong - no peaks were checked!\n");
	} else {
		double ratio = ((double)numGood)/((double)totalChecked);
		fprintf(stderr, "\tLocal Background Filtering: %d of %d (%.2lf%% passed)\n",numGood,totalChecked,ratio*100.0);
	}
	goodPeaks->sortChr();
	return goodPeaks;
}


PeakLibrary* PeakLibrary::getDifferentialPeaks(TagLibrary* tags, TagLibrary* input, 
				double foldThresh,double poissonThresh, int mode, int start, int end, char strand,int strFlag) {


	char* outputstr = new char[10000];
	//int tagsIndex = addTagLibrary(tags);
	//int inputIndex = addTagLibrary(input);
	//Doubletable* expTags = countPeakTags(tagsIndex,start,end,strand,COUNT_MODE_TOTAL);
	//Doubletable* inputTags = countPeakTags(inputIndex,start,end,strand,COUNT_MODE_TOTAL);
	Doubletable* expTags = countPeakTagsLowMemory(tags,strand,COUNT_MODE_TOTAL);
	Doubletable* inputTags = countPeakTagsLowMemory(input,strand,COUNT_MODE_TOTAL);

	//double peakSize = ((double)(end-start));
	if (fixedFlag) {
		//peakSize = avgPeakSize;
	}
	//fprintf(stderr, "peakSize = %lf\n", peakSize);

	double inputPseudoCount = 0.0;
	double tagPseudoCount = 0.0;
	if (input->tbp > 0.0) inputPseudoCount = 0.5;
	if (tags->tbp > 0.0) tagPseudoCount = 0.5;
	//if (input->tbp > 0.0) inputPseudoCount = ((double)(peakSize))*input->tbp;
	//if (tags->tbp > 0.0) tagPseudoCount = ((double)(peakSize))*tags->tbp;
	//fprintf(stderr, "\tpeudo= %lf\n",  tagPseudoCount);

	double minTotal = tags->totalTags;
	if (input->totalTags < tags->totalTags) minTotal = input->totalTags;

	double tagsRatio = minTotal / tags->totalTags;
	double inputRatio = minTotal / input->totalTags;

	int numGood = 0;
	int totalChecked = 0;

	PeakLibrary *goodPeaks = new PeakLibrary();
	goodPeaks->tagsInPeaks = 0.0;
	char** keys = expTags->keys();
	for (int i=0;i<expTags->total;i++) {
		double tp = expTags->search(keys[i]);
		double ip = inputTags->search(keys[i]);
		if (tp < EMPTY_DOUBLE_CHECK || ip < EMPTY_DOUBLE_CHECK) {
			delete [](keys[i]);
			continue;
		}

		totalChecked++; 
		if (ip < inputPseudoCount) ip = inputPseudoCount;
		if (tp < tagPseudoCount) tp = tagPseudoCount;
		double ipn = ip * inputRatio;	
		double tpn = tp * tagsRatio;
		double foldchange = tpn;
		if (ipn > 0) {
			foldchange = tpn / ipn;
		}
		int foldGood = 0;
		if (foldThresh > 0.000001) {
			if (mode == DIFFPEAK_MODE_DIFF) {
				if (foldchange > foldThresh) {
					foldGood = 1;
				}
			} else if (mode == DIFFPEAK_MODE_SAME) {
				if (foldchange < foldThresh && foldchange > 1.0/foldThresh) {
					foldGood = 1;
				}
			} else if (mode == DIFFPEAK_MODE_REV) {
				if (foldchange < 1.0/foldThresh) {
					foldGood = 1;
				}
			}
		} else {
			foldGood = 1;
		}

		int tpInt = (int) tpn;
		//since the input is our "expected" meausrement, that becomes lambda
		//double cumPvalue = cumulativePoisson(tpInt, ipn);
		//double pp = 1-cumPvalue;
		double pp = exp(ilogCumulativePoisson(tpInt,ipn));
		//fprintf(stderr, "%.1lf\t%.1lf\t%lf\n", tpn, ipn, pp);
		int poissonGood = 0;
		if (poissonThresh < 1.0) {
			if (mode == DIFFPEAK_MODE_DIFF) {
				if (pp < poissonThresh && foldchange > 1) {
					poissonGood = 1;
				}
			} else if (mode == DIFFPEAK_MODE_SAME) {
				int ipInt = (int) ipn;
				//double cumPvalueR = cumulativePoisson(ipInt, tpn);
				//double ppR = 1-cumPvalue;
				double ppR = exp(ilogCumulativePoisson(ipInt, tpn));
				if (ppR > poissonThresh && pp > poissonThresh) {
					poissonGood = 1;
				}
			} else if (mode == DIFFPEAK_MODE_REV) {
				int ipInt = (int) ipn;
				//double cumPvalueR = cumulativePoisson(ipInt, tpn);
				//double ppR = 1-cumPvalue;
				double ppR = exp(ilogCumulativePoisson(ipInt, tpn));
				if (ppR < poissonThresh && foldchange < 1) {
					poissonGood = 1;
					pp = ppR;
				}
			}
		} else {
			poissonGood = 1;
		}

		if (poissonGood == 1 && foldGood == 1) {
			goodPeaks->tagsInPeaks += tp;
			Peak* p = (Peak*) peaks->search(keys[i]);
			if (strFlag) {
				sprintf(outputstr, "%.1lf\t%.1lf\t%.2lf\t%.2le", tpn, ipn, foldchange,pp);
				p->addData(outputstr);
			}
			goodPeaks->addPeak(p);
			numGood++;
		}
		delete [](keys[i]);
	}
	delete []keys;
	delete []outputstr;
	delete expTags;
	delete inputTags;

	if (totalChecked == 0) {
		fprintf(stderr, "\t!! Something is wrong - no peaks were checked!\n");
	} else {
		double ratio = ((double)numGood)/((double)totalChecked);
		fprintf(stderr, "\tDifferential Peaks: %d of %d (%.2lf%% passed)\n",numGood,totalChecked,ratio*100.0);
	}
	goodPeaks->sortChr();
	//fprintf(stderr, "totalTagsInPeaks=%lf\n", goodPeaks->tagsInPeaks);
	return goodPeaks;

}

PeakLibrary* PeakLibrary::stitchRegions(int maxDistance,int mode) {
	PeakLibrary* regions = new PeakLibrary();
	//regions->tagsInPeaks = 0.0;
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		ChrPeaks* cp = (ChrPeaks*)chrs->search(keys[i]);
		cp->stitchRegions(regions, maxDistance,mode);
		delete [](keys[i]);
	}
	delete []keys;
	regions->sortChr();
	return regions;
}

void PeakLibrary::analyzeReadAutocorrelation(TagLibrary* tags, int peakSize, char strand, int maxSize, int mode) {

	fprintf(stderr, "\tAnalysing TSS autocorrelation pattern (max length: %d)\n",maxSize);
	int fragLength = 1;
	int halfPeakSize = (int)(peakSize/2.0);
	tags->setTagAdjust(0);

	if (peakSize == 0) {
		setPeakTagSizeFixed(-halfPeakSize-fragLength,halfPeakSize+fragLength);
	} else {
		setPeakTagSizeRefPos(NULL_OFFSET,-halfPeakSize-fragLength,halfPeakSize+fragLength);
	}
	int expIndex = addTagLibrary(tags);

	char** keys = peaks->keys();
	for (int i=0;i<peaks->total;i++) {
		Peak* p = (Peak*)peaks->search(keys[i]);
		p->analyzeReadAutocorrelation(expIndex,maxSize,strand,mode);
		delete [](keys[i]);
	}
	delete []keys;
}
void PeakLibrary::analyzeTSSpattern(TagLibrary* tags, int peakSize, char strand, int mode) {

	fprintf(stderr, "\tAnalysing TSS initiation patterns %d\n",peakSize);
	int fragLength = 1;
	int halfPeakSize = (int)(peakSize/2.0);
	tags->setTagAdjust(0);

	if (peakSize == 0) {
		setPeakTagSizeFixed(-halfPeakSize-fragLength,halfPeakSize+fragLength);
	} else {
		setPeakTagSizeRefPos(NULL_OFFSET,-halfPeakSize-fragLength,halfPeakSize+fragLength);
	}
	int expIndex = addTagLibrary(tags);

	char** keys = peaks->keys();
	for (int i=0;i<peaks->total;i++) {
		Peak* p = (Peak*)peaks->search(keys[i]);
		p->analyzeTSSpattern(expIndex,fragLength,strand,mode);
		delete [](keys[i]);
	}
	delete []keys;
}

void PeakLibrary::centerPeaks(TagLibrary* tags, int peakSize, char strand) {

	fprintf(stderr, "\tCentering peaks of size %d using a fragment length of %d\n",
				peakSize, tags->fragmentLengthEstimate);
	int fragLength = tags->fragmentLengthEstimate;
	int halfPeakSize = (int)(peakSize/2.0);
	
	tags->setTagAdjust(0);
	if (peakSize == 0) {
		setPeakTagSizeFixed(-halfPeakSize-fragLength,halfPeakSize+fragLength);
	} else {
		setPeakTagSizeRefPos(NULL_OFFSET,-halfPeakSize-fragLength,halfPeakSize+fragLength);
	}
	int expIndex = addTagLibrary(tags);

	char** keys = peaks->keys();
	for (int i=0;i<peaks->total;i++) {
		Peak* p = (Peak*)peaks->search(keys[i]);
		p->centerPeak(expIndex,fragLength,strand);
		delete [](keys[i]);
	}
	delete []keys;
}

void PeakLibrary::centerNFR(TagLibrary* tags, int peakSize, char strand, int nfrSize) {

	fprintf(stderr, "\tCentering peaks on Nucleosome Free Region of size %d using a fragment length of %d\n",
				nfrSize, tags->fragmentLengthEstimate);
	int fragLength = tags->fragmentLengthEstimate;
	int halfPeakSize = (int)((nfrSize+peakSize)/2.0);
	int nucSize = 150;	

	tags->setTagAdjust(0);
	if (peakSize == 0) {
		setPeakTagSizeFixed(-1*halfPeakSize-fragLength-nucSize,halfPeakSize+fragLength+nucSize);
	} else {
		setPeakTagSizeRefPos(NULL_OFFSET,-halfPeakSize-fragLength-nucSize,halfPeakSize+fragLength+nucSize);
	}
	int expIndex = addTagLibrary(tags);

	char** keys = peaks->keys();
	for (int i=0;i<peaks->total;i++) {
		Peak* p = (Peak*)peaks->search(keys[i]);
		p->centerNFR(expIndex,fragLength,strand,nfrSize,nucSize);
		delete [](keys[i]);
	}
	delete []keys;
}


void PeakLibrary::annotatePeakLocations(PeakLibrary* annotations, FILE* statsFile, FILE* annFile) {

	Doubletable* stats = new Doubletable();
	Doubletable* sizeTotals = new Doubletable();

	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	fprintf(stderr, "\tAnnotating:");
	for (int i=0;i<chrs->total;i++) {
		ChrPeaks* cp = (ChrPeaks*) chrs->search(keys[i]);
		ChrPeaks* cpAnn = (ChrPeaks*) annotations->chrs->search(keys[i]);
		if (cpAnn == NULL) {
			//fprintf(stderr, "\tNo annotations for on chromosome %s\n", keys[i]);
			continue;
		}
		fprintf(stderr, "."); //, keys[i]);
		cp->annotatePeakLocations(cpAnn, annFile, stats);
		cpAnn->getAnnotationSizeTotals(sizeTotals);
		delete [](keys[i]);
	}
	delete []keys;
	fprintf(stderr, "\n");

	char intron[20] = "Intron";
	char intergenic[20] = "Intergenic";
	char exon[20] = "Exon";
	char promoter[20] = "Promoter";

	if (statsFile != NULL) {
		fprintf(statsFile,"Annotation\tNumber of peaks\tTotal size (bp)\tLog2 Enrichment\n");
	}
	fprintf(stderr,"\t\tAnnotation\tNumber of peaks\tTotal size (bp)\tLog2 Enrichment\n");

	double totalregions = 0.0;
	double totalsize = 0.0;
	keys = sizeTotals->keys();
	for (int i=0;i<sizeTotals->total;i++) {
		double s = sizeTotals->search(keys[i]);
		if (s < EMPTY_DOUBLE_CHECK) s = 0.0;
		double n = stats->search(keys[i]);
		if (n < EMPTY_DOUBLE_CHECK) n = 0.0;
		totalsize += s;
		totalregions += n;
	}
	if (totalsize < 0.001) totalsize = 1.0;
	if (totalregions < 0.001) totalregions = 1.0;

	keys = sizeTotals->keys();
	for (int i=0;i<sizeTotals->total;i++) {
		char* annname = keys[i];
		if (strcmp(keys[i],"I")==0) {
			annname = intron;
		} else if (strcmp(keys[i],"N")==0) {
			annname = intergenic;
		} else if (strcmp(keys[i],"P")==0) {
			annname = promoter;
		} else if (strcmp(keys[i],"E")==0) {
			annname = exon;
		}

		double s = sizeTotals->search(keys[i]);
		if (s < EMPTY_DOUBLE_CHECK) s = 0.0;
		double n = stats->search(keys[i]);
		if (n < EMPTY_DOUBLE_CHECK) n = 0.0;

		double ss = s;
		if (ss < 1) ss = 0.5;
		double expectedFraction = (ss)/totalsize;
		double nn = n;
		if (nn < 1) nn = expectedFraction;
		double obsFraction = nn/totalregions;
		double ratio = obsFraction/expectedFraction;
		if (ratio > 1e-200) {
			ratio = log(ratio)/log(2.0);
		}
		

		if (statsFile != NULL) {
			fprintf(statsFile,"%s\t%.1lf\t%.0lf\t%.3lf\n",annname,n,s,ratio);
		}
		fprintf(stderr,"\t\t%s\t%.1lf\t%.0lf\t%.3lf\n",annname,n,s,ratio);
		delete [](keys[i]);
	}
	delete []keys;
	delete stats;

}

void PeakLibrary::genomeOntology(char** peakFiles, int numPeakFiles, char strand, int maxDistance, 
							long long int gsize, char* controlFileName) {

	FILE* fp = stdout;
		
	//setup statistics
	long long int genomicBins = (int)(gsize);
	long long int refCoverage = calculateCoverage();
	int refNumPeaks = numPeaks;
	int refAvgSize = refCoverage/numPeaks;

	long long int controlCoverage = 0;
	int controlNumPeaks = 0;
	//int controlAvgSize = 0;

	PeakLibrary* control = NULL;

	if (controlFileName != NULL) {
		control = new PeakLibrary(controlFileName,PEAK_READ_MODE_NORMAL);
		controlCoverage = control->calculateCoverage();
		controlNumPeaks = control->numPeaks;
		//controlAvgSize = controlCoverage/controlNumPeaks;
	}
	
	fprintf(fp, "PeakFile/Annotation\t#features[ref=%d]\tCoverage(bp)[ref=%lld]\tAvgFeatureSize[ref=%d]",
						refNumPeaks,refCoverage,refAvgSize);
	fprintf(fp, "\tOverlap(#peaks)\tOverlap(bp)");
	if (maxDistance > 1) {
		fprintf(fp, "\tExpected Overlap(#peaks, gsize=%.2e)", (double)gsize);
	} else {
		fprintf(fp, "\tExpected Overlap(bp, gsize=%.2e)", (double)gsize);
	}
	fprintf(fp, "\tLog Ratio Enrichment\tLog P-value(+ underrepresented)\tP-value");
	if (control != NULL) {
		fprintf(fp, "\tControl Overlap(#peaks,[total=%d])\tControl Overlap(bp,[total=%lld])",
																controlNumPeaks,controlCoverage);
		if (maxDistance > 1) {
			fprintf(fp, "\tControl Expected Overlap(#peaks)");
		} else {
			fprintf(fp, "\tControl Expected Overlap(bp)");
		}
		fprintf(fp, "\tControl Log Ratio Enrichment\tControl Log P-value(+ underrepresented)\tControl P-value");
	
		fprintf(fp, "\tLog Ratio Enrichment(vs. Control)\tLog P-value(+ underrepresented, vs. Control)\tP-value(vs. Control)");
	}
	fprintf(fp, "\n");


	for (int i=0;i<numPeakFiles;i++) {
		//fprintf(stderr, "\t\tProcessing peaks in %s\n", peakFiles[i]);
		PeakLibrary* p = new PeakLibrary(peakFiles[i],  PEAK_READ_MODE_NORMAL);	
		int queryNumPeaks = p->numPeaks;
		long long int queryCoverage = p->calculateCoverage();


		long long int totalOverlap = 0;
		Hashtable* mapping = getOverlappingPeaks(p,strand,maxDistance,totalOverlap);
		//Inttable* overlappedPeaks = new Inttable;
		int totalCommon = mapping->total;
		char** keys = mapping->keys();
		for (int j=0;j<mapping->total;j++) {
			PeakMapping* pm = (PeakMapping*) mapping->search(keys[j]);
			//totalOverlap += pm->bp;
			//for (int k=0;k<pm->numPeaks;k++) {
			//	overlappedPeaks->insert(1,pm->peaks[k]->name);
			//}
			delete pm;
			delete [](keys[j]);
		}
		//int totalCommonRev = overlappedPeaks->total;
		//delete overlappedPeaks;
		delete []keys;
		delete mapping;


		int totalCommonControl = 0;
		long long int totalOverlapControl = 0;

		if (control != NULL) {
			totalOverlapControl = 0;
			Hashtable* mapping = control->getOverlappingPeaks(p,strand,maxDistance,totalOverlapControl);
			totalCommonControl = mapping->total;
			char** keys = mapping->keys();
			for (int j=0;j<mapping->total;j++) {
				PeakMapping* pm = (PeakMapping*) mapping->search(keys[j]);
				//totalOverlapControl += pm->bp;
				delete pm;
				delete [](keys[j]);
			}
			delete []keys;
			delete mapping;
		}
		delete p;
	
		if (queryNumPeaks < 1) continue;
		int queryAvgSize = queryCoverage/queryNumPeaks;
		fprintf(fp,"%s\t%d\t%lld\t%d",peakFiles[i],queryNumPeaks,queryCoverage,queryAvgSize);

		//int binSize = maxDistance * CORRECTION_FACTOR;

		if (maxDistance > 1) {
			genomicBins = gsize/(maxDistance*CORRECTION_FACTOR);
			unsigned int numExpected = (unsigned int) ((((double)refNumPeaks)*((double)queryNumPeaks))
																		/(double)genomicBins +0.5);
			double logp = 0.0;
			if ((unsigned int)totalCommon>=numExpected) {
				logp=loghypergeoD((unsigned int)genomicBins,refNumPeaks, queryNumPeaks,totalCommon);
			} else {
				logp=-1*iloghypergeoD((int)genomicBins,refNumPeaks, queryNumPeaks,totalCommon);
			}
			int top = 1;
			int bottom = 1;
			if (numExpected > 1) bottom = numExpected;
			if (totalCommon > 1) top = totalCommon;
			double logr = log(((double)top)/((double)bottom));
			double ppvalue = 1.0;
			if (logp < 0) {
				ppvalue = exp(logp);
			} else {
				ppvalue = exp(logp*-1);
			}
			fprintf(fp,"\t%d\t%lld\t%d\t%.2lf\t%.2lf\t%.2le",totalCommon,
											totalOverlap,numExpected,logr,logp,ppvalue);

			if (control != NULL) {
				int numExpected =(int) ((((double)controlNumPeaks)*
													((double)queryNumPeaks))/(double)genomicBins +0.5);
				int top = totalCommonControl;
				int bottom = numExpected;
				if (top < 1) top = 1;
				if (bottom < 1) bottom = 1;
				double logr = log(((double)top)/((double)bottom));
				double logp = 0.0;
				double ppvalue = 1.0;
				if (totalCommonControl>=numExpected) {
					logp=loghypergeoD((int)genomicBins,controlNumPeaks, queryNumPeaks,totalCommonControl);
					ppvalue = exp(logp);
				} else {
					logp=-1*iloghypergeoD((int)genomicBins,controlNumPeaks, queryNumPeaks,totalCommonControl);
					ppvalue = exp(logp*-1);
				}
				fprintf(fp,"\t%d\t%lld\t%d\t%.2lf\t%.2lf\t%.2le",totalCommonControl,totalOverlapControl,
																	numExpected,logr,logp,ppvalue);

				bottom = controlNumPeaks;
				if (bottom < 1) bottom = 1;
				int expectControl = (int)(((double)totalCommonControl)*((double)refNumPeaks)
																	/((double)bottom));
				top = totalCommon;
				bottom = expectControl;
				if (top < 1) top = 1;
				if (bottom < 1) bottom = 1;
				logr = log(((double)top)/((double)bottom));
				logp = 0.0;
				ppvalue = 1.0;
				double adjBins = (((double)refNumPeaks)*((double)queryNumPeaks)/((double)bottom));
				if (adjBins > INT_MAX) {
					//probably very very divergent, like just got an inf back. ok screw it
					if (adjBins > 1e30) adjBins = 1e30;
					double rate = ((double)queryNumPeaks)/adjBins;
					if (logr >= 0) {
						logp=logbinomialD(refNumPeaks,totalCommon,rate,INT_MAX);
						ppvalue = exp(logp);
					} else {
						logp=-1*ilogbinomialD(refNumPeaks,totalCommon,rate,INT_MAX);
						ppvalue = exp(-1*logp);
					}
				} else {
					if (logr >= 0.0) {
						logp=loghypergeoD((int)adjBins,refNumPeaks, queryNumPeaks,totalCommon);
						ppvalue = exp(logp);
					} else {
						logp=-1*iloghypergeoD((int)adjBins,refNumPeaks, queryNumPeaks,totalCommon);
						ppvalue = exp(logp*-1);
					}
				}
				fprintf(fp,"\t%.2lf\t%.2lf\t%.2le",logr,logp,ppvalue);
			}

		} else {

			int binSize= refAvgSize;
			genomicBins = gsize/binSize;

			long long int expOverlap = (long long int) (((double)refCoverage)*
													((double)queryCoverage)/((double)gsize));
			long long int top = totalOverlap;
			long long int bottom = expOverlap;
			if (top < 1) top = 1;
			if (bottom < 1) bottom = 1;
			double logr = (double) log(((double)top)/((double)bottom));
			double logp = 0.0;
			double ppvalue = 1.0;
			if (logr > 0.0) {
				logp=loghypergeoD((int)genomicBins,(int)(refCoverage/binSize),
									(int)(queryCoverage/binSize),(int)(totalOverlap/binSize));
				ppvalue = exp(logp);
			} else {
				logp=-1*iloghypergeoD((int)genomicBins,(int)(refCoverage/binSize), 
									(int)(queryCoverage/binSize),(int)(totalOverlap/binSize));
				ppvalue = exp(logp*-1);
			}
			fprintf(fp,"\t%d\t%lld\t%lld\t%.2lf\t%.2lf\t%.2le",totalCommon,
											totalOverlap,expOverlap,logr,logp,ppvalue);

			if (control != NULL) {

				long long int expOverlap = (long long int) (((double)controlCoverage)*
												((double)queryCoverage)/((double)gsize));
				long long int top = totalOverlapControl;
				long long int bottom = expOverlap;
				if (top < 1) top = 1;
				if (bottom < 1) bottom = 1;
				double logr = log(((double)top)/((double)bottom));
				double logp = 0.0;
				double ppvalue = 1.0;
				if (logr > 0.0) {
					logp=loghypergeoD((int)genomicBins,(int)(controlCoverage/binSize),
										(int)(queryCoverage/binSize),(int)(totalOverlapControl/binSize));
					ppvalue = exp(logp);
				} else {
					logp=-1*iloghypergeoD((int)genomicBins,(int)(controlCoverage/binSize),
										(int)(queryCoverage/binSize),(int)(totalOverlapControl/binSize));
					ppvalue = exp(logp*-1);
				}

				fprintf(fp,"\t%d\t%lld\t%lld\t%.2lf\t%.2lf\t%.2le",totalCommonControl,totalOverlapControl,
																	expOverlap,logr,logp,ppvalue);

				bottom = controlCoverage;
				if (bottom < 1) bottom = 1;
				long long int expectControl = (int)(((double)totalOverlapControl)*((double)refCoverage)
															/((double)bottom));
				top = totalOverlap;
				bottom = expectControl;
				if (top < 1) top = 1;
				if (bottom < 1) bottom = 1;
				logr = log(((double)top)/((double)bottom));
				logp = 0.0;
				ppvalue = 1.0;

				bottom = expectControl;
				if (bottom < 1) bottom = 1;
				double adjBins = (((double)refCoverage)*((double)queryCoverage)/((double)expectControl))
																		/((double)binSize);
				if (adjBins > INT_MAX) {
					//probably very very divergent, like just got an inf back. ok screw it
					if (adjBins > 1e30) adjBins = 1e30;
					double rate = ((double)queryCoverage/binSize)/adjBins;
					if (logr >= 0) {
						logp=logbinomialD((int)(refCoverage/binSize),
											(int)(totalOverlap/binSize),rate,INT_MAX);
						ppvalue = exp(logp);
					} else {
						logp=-1*ilogbinomialD((int)(refCoverage/binSize),
											(int)(totalOverlap/binSize),rate,INT_MAX);
						ppvalue = exp(-1*logp);
					}
				} else {
					if (logr >= 0.0) {
						logp=loghypergeoD((int)(adjBins),(int)(refCoverage/binSize), 
								(int)(queryCoverage/binSize),(int)(totalOverlap/binSize));
						ppvalue = exp(logp);
					} else {
						logp=-1*iloghypergeoD((int)(adjBins),(int)(refCoverage/binSize), 
								(int)(queryCoverage/binSize),(int)(totalOverlap/binSize));
						ppvalue = exp(logp*-1);
					}
				}
				fprintf(fp,"\t%.2lf\t%.2lf\t%.2le",logr,logp,ppvalue);
			}
		}
		fprintf(fp,"\n");

	}

	if (control != NULL) {
		delete control;
	}

}

void PeakLibrary::getCoBoundPeaks(char** peakFiles, int numPeakFiles, char strand, int maxDistance, 
						long long int gsize, char* prefix, int maxCoBound, char* matrixFile,char* cmdline) {
	
	Inttable* cobound = new Inttable(10000000);
	char** keys = peaks->keys();
	setDefaultPeakOrder();
	int* ogStarts = new int[peaks->total];
	int* ogEnds = new int[peaks->total];
	for (int i=0;i<peaks->total;i++) {
		cobound->insert(0,keys[i]);
		ogStarts[i] = peakOrder[i]->start;
		ogEnds[i] = peakOrder[i]->end;
	}
	//keep keys around for later

	if (maxDistance > 0) {
		setPeakSize(maxDistance);
	}
	int avgPeakSize = getAveragePeakSize();
	//fprintf(stderr, "avgPeakSize= %d\n", avgPeakSize);
	long long int totalCoverage = calculateCoverage();

	FILE* pfp = NULL;
	FILE* rfp = NULL;
	FILE* nfp = NULL;
	if (matrixFile != NULL) {
		char* fname = new char[10000];
		strcpy(fname,matrixFile);
		strcat(fname,".logPvalue.txt");
		pfp = fopen(fname, "w");
		if (pfp == NULL) fprintf(stderr, "!!! Could not open %s for writing!!!\n", fname);
		strcpy(fname,matrixFile);
		strcat(fname,".logRatio.txt");
		rfp = fopen(fname, "w");
		if (rfp == NULL) fprintf(stderr, "!!! Could not open %s for writing!!!\n", fname);
		strcpy(fname,matrixFile);
		strcat(fname,".count.txt");
		nfp = fopen(fname, "w");
		if (nfp == NULL) fprintf(stderr, "!!! Could not open %s for writing!!!\n", fname);
		delete []fname;
		fprintf(pfp, "natural log p-values for peak overlaps (+values for divergence) (cmd = %s)",cmdline);	
		fprintf(rfp, "natural log ratios of observed/expected peak overlaps (cmd = %s)", cmdline);	
		fprintf(nfp, "peak overlap totals (cmd = %s)", cmdline);	
		fprintf(pfp, "\tvs Reference Peaks\n");
		fprintf(rfp, "\tvs Reference Peaks\n");
		fprintf(nfp, "\tvs Reference Peaks\n");
	}

	int adjustedTotalPeaks = numPeaks;
	if (maxDistance > 0) {
		long long int totalOverlap=0;
		Hashtable* mapping = getOverlappingPeaks(this,strand,-1,totalOverlap);
		adjustedTotalPeaks = mapping->total;
		char** mkeys = mapping->keys();
		for (int i=0;i<mapping->total;i++) {
			PeakMapping* pm = (PeakMapping*) mapping->search(mkeys[i]);
			delete pm;
			delete [](mkeys[i]);
		}
		delete mapping;
		delete []mkeys;
	}


	fprintf(stderr, "\n\tComparing peaks: (peakfile, overlapping peaks, logRatio(obs/expected), logP)\n");
	for (int i=0;i<numPeakFiles;i++) {
		if (matrixFile != NULL) {
			fprintf(pfp, "%s", peakFiles[i]);
			fprintf(rfp, "%s", peakFiles[i]);
			fprintf(nfp, "%s", peakFiles[i]);
		}
	
		PeakLibrary* p = new PeakLibrary(peakFiles[i],  PEAK_READ_MODE_NORMAL);	

		if (maxDistance > 0) p->setPeakSize(maxDistance);
		int curAvgPeakSize = p->getAveragePeakSize();
		long long int curCoverage = p->calculateCoverage();
		int adjustedTotalPeaksI = p->numPeaks;

		if (maxDistance > 0) {
			long long int totalOverlap=0;
			Hashtable* mapping = p->getOverlappingPeaks(p,strand,-1,totalOverlap);
			adjustedTotalPeaksI = mapping->total;
			char** mkeys = mapping->keys();
			for (int j=0;j<mapping->total;j++) {
				PeakMapping* pm = (PeakMapping*) mapping->search(mkeys[j]);
				if (pm != NULL) delete pm;
				delete [](mkeys[j]);
			}
			delete mapping;
			delete []mkeys;
		}

		long long int totalOverlap=0;
		Hashtable* mapping = getOverlappingPeaks(p,strand,-1,totalOverlap);
		//fprintf(stderr, "totalOverlap=%lld\n", totalOverlap);
		int totalCommon = mapping->total;
		char** mkeys = mapping->keys();
		for (int j=0;j<mapping->total;j++) {
			int co = cobound->search(mkeys[j]);
			if (co == EMPTY_INT) {
				fprintf(stderr, "\tSomething might be wrong....\n");
			}	
			cobound->insert(co+1,mkeys[j]);
			PeakMapping* pm = (PeakMapping*) mapping->search(mkeys[j]);
			delete pm;
			delete [](mkeys[j]);
		}

		double logp = 0.0;
		double logr = 0.0;

		if (maxDistance > 0) {
			int binSize = maxDistance;
			long long int genomicBins = gsize/(binSize*2);
			long long int expOverlap = (long long int) (((double)adjustedTotalPeaks)*((double)adjustedTotalPeaksI)/
														((double)genomicBins));
			long long int top = totalCommon;
			long long int bottom = expOverlap;
			if (top < 1) top = 1;
			if (bottom < 1) bottom = 1;
			logr = log(((double)top)/((double)bottom));
			logp = 0.0;
			//fprintf(stderr, "%d\t%lld\t%lld\n", binSize,totalCoverage,curCoverage);
			if (logr > 0.0) {
				logp=loghypergeoD((int)genomicBins,adjustedTotalPeaks,adjustedTotalPeaksI,totalCommon);
				//fprintf(stderr, "%d\t%d\t%d\t%d\n",(int)genomicBins,adjustedTotalPeaks,adjustedTotalPeaksI,totalCommon);
			} else {
				logp=-1*iloghypergeoD((int)genomicBins,adjustedTotalPeaks,adjustedTotalPeaksI,totalCommon);
			}
			fprintf(stderr, "\t**\t%s\t%d\t%.2lf\t%.2lf\n", peakFiles[i],totalCommon,logr,logp);
			//fprintf(stderr, "\t\tStats disabled with fixed distance\n");
		} else {
			int binSize = curAvgPeakSize;
			int binSize2 = avgPeakSize;
			if (binSize2 < binSize) {
				binSize = binSize2;
			}
			long long int genomicBins = gsize/binSize;
	
			long long int expOverlap = (long long int) (((double)totalCoverage)*((double)curCoverage)/((double)gsize));
			long long int top = totalOverlap;
			long long int bottom = expOverlap;
			if (top < 1) top = 1;
			if (bottom < 1) bottom = 1;
			logr = log(((double)top)/((double)bottom));
			logp = 0.0;
			//fprintf(stderr, "%d\t%lld\t%lld\n", binSize,totalCoverage,curCoverage);
			if (logr > 0.0) {
				logp=loghypergeoD((int)genomicBins,(int)(totalCoverage/binSize),
									(int)(curCoverage/binSize),(int)(totalOverlap/binSize));
			} else {
				logp=-1*iloghypergeoD((int)genomicBins,(int)(totalCoverage/binSize),
									(int)(curCoverage/binSize),(int)(totalOverlap/binSize));
			}
			fprintf(stderr, "\t\t%s\t%d\t%.2lf\t%.2lf\n", peakFiles[i],totalCommon,logr,logp);
		}
		if (matrixFile != NULL) {
			fprintf(pfp, "\t%lf\n", logp);
			fprintf(rfp, "\t%lf\n", logr);
			fprintf(nfp, "\t%d\n", totalCommon);
		}
		delete []mkeys;
		delete mapping;
		delete p;
	}
	if (maxDistance > 0) {
		fprintf(stderr, "\n\t**\tPairwise stats are approx with fixed distance (-d %d) and -cobound #\n", maxDistance);
		fprintf(stderr, "\t\tThey get worse as the size increases and peaks from a single file start overlapping\n");
		fprintf(stderr, "\t\tTo get accurate ones, adjust peak sizes first with adjustPeakFile.pl and\n");
		fprintf(stderr, "\t\tthen rerun mergePeaks with the \"-d given\" option (only applies to -cobound #)\n");
	}
	if (matrixFile != NULL) {
		fclose(pfp);
		fclose(rfp);
		fclose(nfp);
	}

	fprintf(stderr, "\n");

	//if adjusted peak sizes, return to the original for output
	for (int i=0;i<numPeaks;i++) {
		peakOrder[i]->start = ogStarts[i];
		peakOrder[i]->end = ogEnds[i];
	}

	char* filename = new char[100000];
	FILE** fps = new FILE*[maxCoBound+1];
	for (int i=0;i<maxCoBound+1;i++) {
		if (prefix != NULL) {
			sprintf(filename, "%s.coBoundBy%d.txt",prefix, i);
		} else {
			sprintf(filename, "coBoundBy%d.txt", i);
		}
		fps[i] = fopen(filename,"w");
		if (fps[i] == NULL) {
			fprintf(stderr, "!!! Couldn't open %s for writing!!!!\n", filename);
		}
		fprintf(fps[i],"#Name (cmd = %s)\tchr\tstart\tend\tstrand\tstat\tstat2\n",cmdline);
	}
	delete []filename;

	int* totals = new int[numPeakFiles+1];
	for (int i=0;i<numPeakFiles+1;i++) totals[i]=0;

	for (int i=0;i<peaks->total;i++) {
		int co = cobound->search(keys[i]);
		if (co == EMPTY_INT) {
			fprintf(stderr, "!!! Something is wrong!!!\n");
		}
		if (co>=0 && co<numPeakFiles+1) {
			totals[co]++;
			if (co > maxCoBound) co = maxCoBound;
			Peak* p = (Peak*) peaks->search(keys[i]);
			p->print(fps[co]);
		} else {
			fprintf(stderr, "!!! number of co-bound peaks is impossible!!!\n");
		}
		delete [](keys[i]);
	}	
	delete []keys;

	for (int i=0;i<numPeakFiles+1;i++) {
		fprintf(stderr, "\tCo-bound by %d peaks: %d", i, totals[i]);
		if (i==maxCoBound) {
			int t = 0;
			for (int j=i;j<numPeakFiles+1;j++) {
				t+=totals[j];
			}
			fprintf(stderr, " (max: %d effective total)", t);
		}
		fprintf(stderr, "\n");
		if (i<=maxCoBound) {
			fclose(fps[i]);
		}
	}
	delete []totals;
	delete []fps;
	delete []ogStarts;
	delete []ogEnds;
	delete cobound;

}

void PeakLibrary::setPeakSize(int size) {
	char** keys = peaks->keys();
	for (int i=0;i<peaks->total;i++) {
		Peak* p = (Peak*)peaks->search(keys[i]);
		p->setPeakSize(size);
		delete [](keys[i]);
	}
	delete []keys;
}

PeakLibrary* PeakLibrary::mergeNearbyPeaks(int maxDistance, char strand) {

	PeakLibrary* merged = new PeakLibrary(numPeaks);
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {

		ChrPeaks* cp = (ChrPeaks*)chrs->search(keys[i]);
		cp->mergeNearbyPeaks(merged,maxDistance,strand);
		delete [](keys[i]);
	}
	delete []keys;
	merged->sortChr();
	fprintf(stderr, "\t%.2lf%% of %d peaks were within %d bp and merged\n",
					((double)(numPeaks-merged->numPeaks))/((double)numPeaks)*100.0,numPeaks,maxDistance);
	return merged;
}

PeakLibrary* PeakLibrary::filterPeaksOutsideRange(char* chr, int start, int end) {

	if (chr == NULL) return NULL;

	PeakLibrary* newpeaks = new PeakLibrary(numPeaks);
	ChrPeaks* cp = (ChrPeaks*) chrs->search(chr);
	for (int i=0;i<cp->numPeaks;i++) {
		if (cp->peaks[i]->start < end && cp->peaks[i]->end > start) {
			newpeaks->addPeak(cp->peaks[i]);
		}
	}
	newpeaks->sortChr();
	fprintf(stderr, "\t%.2lf%% of %d peaks were within %s:%d-%d\n",
					((double)(newpeaks->numPeaks))/((double)numPeaks)*100.0,numPeaks,chr,start,end);
	return newpeaks;
}

void PeakLibrary::mergePeaks(char** peakFiles, int numFiles, char strand, int maxDistance,
						long long int gsize, char* prefix, char* vennFile, 
						char* matrixFile,int codeFlag, char* cmdline) {

	PeakLibrary** peaksets = new PeakLibrary*[numFiles];
	char** filenamemod = new char*[numFiles];

	Inttable** mergemap = new Inttable*[numFiles];
	int* avgPeakSize = new int[numFiles];
	long long int* totalCoverage = new long long int[numFiles];

	int fixedDistanceFlag = 0;
	if (maxDistance >= 0) {
		fixedDistanceFlag = maxDistance;
	}

	int maxNumPeaks = 0;
	int** commonPeaks = new int*[numFiles];
	long long int** commonOverlap = new long long int*[numFiles];

	for (int i=0;i<numFiles;i++) {
		commonPeaks[i] = new int[numFiles];
		commonOverlap[i] = new long long int[numFiles];

		peaksets[i] = new PeakLibrary(peakFiles[i],PEAK_READ_MODE_NORMAL);
		if (peaksets[i]->numPeaks > maxNumPeaks) {
			maxNumPeaks = peaksets[i]->numPeaks;
		}

		if (maxDistance >= 0) {
			peaksets[i]->setPeakSize(maxDistance);
			avgPeakSize[i] = maxDistance;
		} else {
			avgPeakSize[i] = peaksets[i]->getAveragePeakSize();
		}
		totalCoverage[i] = peaksets[i]->calculateCoverage();
		//fprintf(stderr, "%s: %ld bp\n", peakFiles[i], cov);

		mergemap[i] = new Inttable();
		filenamemod[i] = new char[strlen(peakFiles[i])+1];
		strcpy(filenamemod[i], peakFiles[i]);
		for (unsigned int j=0;j<strlen(peakFiles[i]);j++) {
			if (filenamemod[i][j] == '/') {
				filenamemod[i][j] = '_';
			}
		}
	}
	maxDistance = -1;

	int totalPossibleMerges = maxNumPeaks*numFiles*(numFiles);
	//fprintf(stderr, "total possible %d\n", totalPossibleMerges);
	int* possibleMerges = new int[10000];
	//int numPossibleMerges = 0;
	MergePeak** mergedPeaks = new MergePeak*[totalPossibleMerges];
	int numMerged = 0;

	long long int genomicBins = (int)(gsize);

	FILE* pfp = NULL;
	FILE* rfp = NULL;
	FILE* nfp = NULL;
	if (matrixFile != NULL) {
		char* fname = new char[10000];
		strcpy(fname,matrixFile);
		strcat(fname,".logPvalue.matrix.txt");
		pfp = fopen(fname, "w");
		if (pfp == NULL) fprintf(stderr, "!!! Could not open %s for writing!!!\n", fname);
		strcpy(fname,matrixFile);
		strcat(fname,".logRatio.matrix.txt");
		rfp = fopen(fname, "w");
		if (rfp == NULL) fprintf(stderr, "!!! Could not open %s for writing!!!\n", fname);
		strcpy(fname,matrixFile);
		strcat(fname,".count.matrix.txt");
		nfp = fopen(fname, "w");
		if (nfp == NULL) fprintf(stderr, "!!! Could not open %s for writing!!!\n", fname);

		fprintf(pfp, "natural log p-values for peak overlaps (+values for divergence), cmd = %s",cmdline);	
		fprintf(rfp, "natural log ratios of observed/expected peak overlaps, cmd = %s",cmdline);	
		fprintf(nfp, "peak overlap totals, cmd = %s",cmdline);	
		for (int i=0;i<numFiles;i++) {
			fprintf(pfp, "\t%s", peakFiles[i]);
			fprintf(rfp, "\t%s", peakFiles[i]);
			fprintf(nfp, "\t%s", peakFiles[i]);
		}
		fprintf(pfp, "\n");
		fprintf(rfp, "\n");
		fprintf(nfp, "\n");
			
	}

	for (int i=0;i<numFiles;i++) {
		int totalPeaksI = peaksets[i]->numPeaks;

		for (int j=0;j<numFiles;j++) {

			int totalPeaksJ = peaksets[j]->numPeaks;
			fprintf(stderr, "\tComparing %s (%d total) and %s (%d total)\n", peakFiles[i], 
								totalPeaksI,peakFiles[j],totalPeaksJ);

			long long int totalOverlap = 0;
			Hashtable* mapping = peaksets[i]->getOverlappingPeaks(peaksets[j],strand,maxDistance,totalOverlap);
			//fprintf(stderr, "totalOverlap=%lld\n", totalOverlap);
			int totalCommon = mapping->total;
			commonPeaks[i][j]=totalCommon;
			commonOverlap[i][j] = totalOverlap;

			char** keys = mapping->keys();
			//fprintf(stderr, "%d\n",mapping->total);
			for (int k=0;k<mapping->total;k++) {
				PeakMapping* pm = (PeakMapping*) mapping->search(keys[k]);
				//totalOverlap += pm->bp;
				int numPossibleMerges = 0;
				int iIndex = mergemap[i]->search(keys[k]);
				if (iIndex != EMPTY_INT) {
					possibleMerges[numPossibleMerges] = iIndex;
					numPossibleMerges++;
				}
				for (int a=0;a<pm->numPeaks;a++) {
					int index = mergemap[j]->search(pm->peaks[a]->name);
					if (index != EMPTY_INT) {
						int found = 0;
						for (int b=0;b<numPossibleMerges;b++) {
							if (possibleMerges[b] == index) {
								found = 1;
								break;
							}
						}
						if (found == 0) {
							possibleMerges[numPossibleMerges] = index;
							numPossibleMerges++;
						}
					}
				}
				
				int mergerIndex = 0;
				if (numPossibleMerges == 1) {
					mergerIndex = possibleMerges[0];
				} else {
					mergerIndex = numMerged;
					numMerged++;
					mergedPeaks[mergerIndex] = new MergePeak(numFiles);
					for (int a=0;a<numPossibleMerges;a++) {
						int mpIndex = possibleMerges[a];
						MergePeak* mp = mergedPeaks[mpIndex];
						mergedPeaks[mpIndex] = NULL;
						for (int b=0;b<numFiles;b++) {
							for (int c=0;c<mp->numPeaks[b];c++) {
								mergemap[b]->insert(mergerIndex,mp->exps[b][c]->name);
							}
						}
						mergedPeaks[mergerIndex]->add(mp);
						delete mp;
					}
				}
				mergedPeaks[mergerIndex]->add(i,pm->source);
				mergemap[i]->insert(mergerIndex,keys[k]);
				for (int a=0;a<pm->numPeaks;a++) {
					mergedPeaks[mergerIndex]->add(j,pm->peaks[a]);
					mergemap[j]->insert(mergerIndex,pm->peaks[a]->name);
				}
				delete [](keys[k]);
				delete pm;
			}

			//fprintf(stderr, "numMerged = %d\n\n\n", numMerged);

			delete []keys;
			delete mapping;

		}
	}

	//comparison
	for (int i=0;i<numFiles;i++) {
		int totalPeaksI = commonPeaks[i][i];
		long long int totalCoverageI = totalCoverage[i];
		if (matrixFile != NULL) {
			fprintf(pfp, "%s", peakFiles[i]);
			fprintf(rfp, "%s", peakFiles[i]);
			fprintf(nfp, "%s", peakFiles[i]);
		}

		for (int j=0;j<numFiles;j++) {
			int totalPeaksJ = commonPeaks[j][j];
			long long int totalCoverageJ = totalCoverage[j];

			int totalCommon = commonPeaks[i][j];
			long long int totalOverlap = commonOverlap[i][j];
			double logr;
			double logp;

			if (fixedDistanceFlag != 0) {
				int binSize = fixedDistanceFlag;
				genomicBins = gsize/binSize;
				if (i==j) {
					int newCommon = totalPeaksI-totalCommon;
					totalCommon = newCommon;
				}
				int numExpected =(int) ((((double)totalPeaksI)*((double)totalPeaksJ))/(double)genomicBins + 0.5);
				//fprintf(stderr, "%d\t%lld\n", numExpected,genomicBins);
				if (totalCommon>=numExpected) {
					logp=loghypergeoD((int)genomicBins,totalPeaksI, totalPeaksJ,totalCommon);
				} else {
					logp=-1*iloghypergeoD((int)genomicBins,totalPeaksI, totalPeaksJ,totalCommon);
				}
				int top = totalCommon;
				int bottom = numExpected;
				if (top < 1) top = 1;
				if (bottom < 1) bottom = 1;
				logr = log(((double)top)/((double)bottom));
			} else {
				int binSize = avgPeakSize[i];
				int binSize2 = avgPeakSize[j];
				if (binSize2 < binSize) {
					binSize = binSize2;
				}
				if (binSize < 1) binSize = 1;
				genomicBins = gsize/binSize;

				long long int expOverlap = (long long int) (((double)totalCoverageI)*((double)totalCoverageJ)/((double)gsize));
				long long int top = totalOverlap;
				long long int bottom = expOverlap;
				if (top < 1) top = 1;
				if (bottom < 1) bottom = 1;
				logr = log(((double)top)/((double)bottom));
				logp = 0.0;
			//	fprintf(stderr, "%d\t%lld\t%lld\n", binSize,totalCoverageI,totalCoverageJ);
				if (logr > 0.0) {
					logp=loghypergeoD((int)genomicBins,(int)(totalCoverageI/binSize),
									(int)(totalCoverageJ/binSize),(int)(totalOverlap/binSize));
				} else {
					logp=-1*iloghypergeoD((int)genomicBins,(int)(totalCoverageI/binSize),
									(int)(totalCoverageJ/binSize),(int)(totalOverlap/binSize));
				}
				//fprintf(stderr, "\t%d\t%lld / %lld\t%d\t%lf\n", binSize, totalOverlap,expOverlap, totalCommon,logp);
			}
			//	if (matrixFile != NULL) {
			if (matrixFile != NULL) {
				fprintf(pfp, "\t%lf", logp);
				fprintf(rfp, "\t%lf", logr);
				fprintf(nfp, "\t%d", totalCommon);
			}
		}
		//if (matrixFile != NULL) {
		if (matrixFile != NULL) {
			fprintf(pfp, "\n");
			fprintf(rfp, "\n");
			fprintf(nfp, "\n");
		}
	}
	if (matrixFile != NULL) {
		fclose(pfp);
		fclose(rfp);
		fclose(nfp);
	}

	fprintf(stderr, "\n");

	char* codestring = new char[numFiles+1];
	codestring[numFiles] = '\0';
	Hashtable* combinations = new Hashtable();
	for (int i=0;i<numMerged;i++) {
		if (mergedPeaks[i] == NULL) continue;
		mergedPeaks[i]->combine();
		for (int j=0;j<numFiles;j++) {
			if (mergedPeaks[i]->numPeaks[j] > 0) {
				codestring[j] = '1';
			} else {
				codestring[j] = '0';
			}
		}
		MergePeaksArray* subgroup = (MergePeaksArray*)combinations->search(codestring);
		if (subgroup == NULL) {
			subgroup = new MergePeaksArray();
			combinations->insert(subgroup,codestring);
			subgroup->p = new MergePeak*[numMerged]; // numMerged is max
			subgroup->n = 0;
		}
		subgroup->p[subgroup->n]=mergedPeaks[i];
		subgroup->n++;
	}

	FILE* vennfp = NULL;
	if (vennFile == NULL) {
		vennfp = stderr;
	} else {
		vennfp = fopen(vennFile, "w");
		if (vennfp == NULL) {
			fprintf(stderr, "!!! Couldn't open venn diagram output file: %s !!!\n", vennFile);
			vennfp = stderr;
		}
	}
	if (codeFlag==0) {
		fprintf(vennfp, "%s", peakFiles[0]);
		for (int i=1;i<numFiles;i++) {
			fprintf(vennfp, "\t%s", peakFiles[i]);
		}
	} else {
		fprintf(vennfp, "Code");
	}
	fprintf(vennfp, "\tTotal\tName\n");


	if (prefix == NULL) {
		fprintf(stdout, "#name (cmd = %s)\tchr\tstart\tend\tstrand\tStat\tParent files\tTotal subpeaks",cmdline);
		for (int i=0;i<numFiles;i++){ 
			fprintf(stdout, "\t%s",peakFiles[i]);
		}
		fprintf(stdout, "\n");
	}

	char** keys = combinations->keys();
	for (int i=0;i<combinations->total;i++) {
		char* filename = new char[100000];
		char* filestr = new char[100000];
		filename[0] = '\0';
		filestr[0] = '\0';
		if (prefix != NULL) {
			strcpy(filename,prefix);
			strcat(filename, "_");
		}
		int totalNum=0;
		if (codeFlag==0) {
			for (int j=0;j<numFiles;j++) {
				if (j>0) {
					fprintf(vennfp, "\t");
				}
				if (keys[i][j] == '1') {
					totalNum++;
					fprintf(vennfp, "X");
					if (totalNum > 1) {
						strcat(filestr, "|");
						strcat(filename, "_");
					}
					strcat(filename, filenamemod[j]);
					strcat(filestr, peakFiles[j]);
				}
			}
		} else {
			fprintf(vennfp, "%s", keys[i]);
			strcat(filename, keys[i]);
			strcat(filestr, keys[i]);
		}
		//fprintf(stderr, "\nfilename=%s\nfilestr=%s\n\n", filename,filestr);
	
		FILE *fp = NULL;
		if (prefix != NULL) {
			fp = fopen(filename,"w");
			if (fp == NULL) {
				fprintf(stderr, "!!! Could not open file %s for writing\n", filename);
				exit(0);
			}
			fprintf(fp, "#name (cmd = %s)\tchr\tstart\tend\tstrand\tStat\tParent files\tTotal subpeaks",cmdline);
			for (int k=0;k<numFiles;k++){ 
				fprintf(fp, "\t%s",peakFiles[k]);
			}
			fprintf(fp, "\n");
		} else {
			fp = stdout;
		}
	
		MergePeaksArray* mpa = (MergePeaksArray*)combinations->search(keys[i]);
		//fprintf(stderr, "===== %s\n", keys[i]);
		for (int j=0;j<mpa->n;j++) {
			fprintf(fp,"%s\t%s\t%d\t%d\t", mpa->p[j]->merged->name, mpa->p[j]->merged->chr,
									mpa->p[j]->merged->start, mpa->p[j]->merged->end);
			if (mpa->p[j]->merged->strand == 0) {
				fprintf(fp, "+");
			} else {
				fprintf(fp, "-");
			}
			fprintf(fp, "\t%f\t%s\t%.0lf", mpa->p[j]->merged->v, filestr,mpa->p[j]->merged->focusRatio);
			for (int k=0;k<numFiles;k++) {
				if (mpa->p[j]->exps[k] != NULL && mpa->p[j]->numPeaks[k] > 0) {
					fprintf(fp, "\t%s", mpa->p[j]->exps[k][0]->name);
					for (int kk=1;kk<mpa->p[j]->numPeaks[k];kk++){ 
						fprintf(fp, ",%s", mpa->p[j]->exps[k][kk]->name);
					}
				} else {
					fprintf(fp, "\t");
				}
			}
			fprintf(fp, "\n");
		}
		int numCommon = mpa->n;
		delete mpa;

		if (prefix != NULL && fp != stdout) {
			fclose(fp);
		}
		fprintf(vennfp, "\t%d\t%s\n",numCommon, filestr);
		delete []filename;
		delete []filestr;
	}
	
	//print out leftovers
	for (int i=0;i<numFiles;i++) {
		break;
		char* filename = new char[100000];
		char* filestr = new char[100000];
		filename[0] = '\0';
		filestr[0] = '\0';
		if (prefix != NULL) {
			strcpy(filename,prefix);
			strcat(filename, "_");
		}
		if (codeFlag==0) {
			strcat(filename, filenamemod[i]);
			strcat(filestr, peakFiles[i]);
			for (int j=0;j<numFiles;j++) {
				if (j>0) {
					fprintf(vennfp, "\t");
				}
				if (j==i) {
					fprintf(vennfp, "X");
				}
			}
		} else {
			int startindex = strlen(filename);
			for (int j=0;j<numFiles;j++) {
				if (j == i) {
					fprintf(vennfp, "1");
					filename[startindex+j] = '1';
					filestr[j] = '1';
				} else {
					fprintf(vennfp, "0");
					filename[startindex+j] = '0';
					filestr[j] = '0';
				}
			}
			filename[startindex+numFiles] = '\0';
			filestr[numFiles] = '\0';
		}
	
		FILE *fp = NULL;
		if (prefix != NULL) {
			fp = fopen(filename,"w");
			if (fp == NULL) {
				fprintf(stderr, "!!! Could not open file %s for writing\n", filename);
				exit(0);
			}
		} else {
			fp = stdout;
		}

		int totalLeft = 0;
		char** keys = peaksets[i]->peaks->keys();
		for (int j=0;j<peaksets[i]->peaks->total;j++) {
			if (EMPTY_INT == mergemap[i]->search(keys[j]) ) {
				totalLeft++;
				Peak* p = (Peak*)peaksets[i]->peaks->search(keys[j]);
				//p->print(fp);
				fprintf(fp,"%d-%s\t%s\t%d\t%d\t", i+1,p->name, p->chr,p->start,p->end);
				if (p->strand == 0) {
					fprintf(fp, "+");
				} else {
					fprintf(fp, "-");
				}
				fprintf(fp, "\t%f\t%s\t1\n", p->v, filestr);
			}
			delete [](keys[j]);
		}
		delete []keys;
		delete []filename;
	
		fprintf(vennfp, "\t%d\t%s\n", totalLeft, filestr);
		if (prefix != NULL && fp != stdout) {
			fclose(fp);
		}
	}

	if (vennFile != NULL && vennfp != stderr) {
		fclose(vennfp);
	}

	for (int i=0;i<numFiles;i++) {
		delete peaksets[i];
		delete [](filenamemod[i]);
		delete mergemap[i];
		delete [](commonPeaks[i]);
		delete [](commonOverlap[i]);
	}
	delete []commonPeaks;
	delete []commonOverlap;
	delete []peaksets;
	delete []filenamemod;
	delete []mergemap;

	for (int i=0;i<numMerged;i++) {
		delete mergedPeaks[i];
	}
	delete []mergedPeaks;
}

long long int PeakLibrary::calculateCoverage() {

	long long int bp=0;
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		ChrPeaks* cp = (ChrPeaks*) chrs->search(keys[i]);
		bp += cp->calculateCoverage();
		delete [](keys[i]);
	}
	delete []keys;

	return bp;
}



MergePeaksArray::MergePeaksArray() {
	p = NULL;
	n=0;
}
MergePeaksArray::~MergePeaksArray() {
	if (p != NULL) delete []p;
}

MergePeak::MergePeak(int n) {
	numExps = n;
	exps = new Peak**[numExps];
	numPeaks = new int[numExps];
	for (int i=0;i<numExps;i++) {
		exps[i] = NULL;
		numPeaks[i] = 0;
	}
	merged = NULL;
} 
MergePeak::~MergePeak() {
	if (exps != NULL) {
		for (int i=0;i<numExps;i++) {
			if (exps[i] != NULL) {
				delete [](exps[i]);
			}
		}
		delete []exps;
	}
	if (numPeaks != NULL) {
		delete []numPeaks;
	}
}
void MergePeak::add(MergePeak* mp) {
	for (int i=0;i<mp->numExps;i++) {
		Peak** newpeaks = new Peak*[numPeaks[i]+mp->numPeaks[i]];
		for (int j=0;j<numPeaks[i];j++) {
			newpeaks[j] = exps[i][j];
		}
		int curNumPeaks = numPeaks[i];
		for (int j=0;j<mp->numPeaks[i];j++) {
			int dontAdd = 0;
			for (int k=0;k<curNumPeaks;k++) {
				if (strcmp(mp->exps[i][j]->name,newpeaks[k]->name)==0) {
					dontAdd=1;
					break;
				}
			}
			if (dontAdd) break;
			newpeaks[curNumPeaks++]=mp->exps[i][j];
		}
		delete [](exps[i]);
		exps[i] = newpeaks;
		numPeaks[i] = curNumPeaks;
	}
}
void MergePeak::add(int index, Peak* p) {
	int addFlag = 1;
	for (int i=0;i<numPeaks[index];i++) {
		if (strcmp(p->name, exps[index][i]->name)==0) {
			addFlag = 0;
			break;
		}
	}
	if (addFlag == 0) {
		return;
	}
	Peak** peaks = new Peak*[numPeaks[index]+1];
	for (int i=0;i<numPeaks[index];i++) {
		peaks[i] = exps[index][i];
	}
	if (exps[index] != NULL) {
		delete [](exps[index]);
	}
	exps[index] = peaks;
	exps[index][numPeaks[index]] = p;
	numPeaks[index]++;
}
void MergePeak::combine() {
	long long int position = 0;
	int n = 0;
	double vTotal = 0.0;
	int N = 0;
	for (int i=0;i<numExps;i++) {
		for (int j=0;j<numPeaks[i];j++) {
			if (n == 0) {
				merged = exps[i][j]->copy();
			}
			vTotal += (double)(exps[i][j]->v);
			int v = 1;
			position += ((long long int)exps[i][j]->start)*v;
			position += ((long long int)exps[i][j]->end)*v;
			if (exps[i][j]->start < merged->start) merged->start = exps[i][j]->start;
			if (exps[i][j]->end > merged->end) merged->end = exps[i][j]->end;
			n+=2*v;
			N++;
		}
	}
	if (n > 0) {
		position /= n;
	}

	//int length = (merged->end-merged->start)/2;
	//merged->start = ((int)position)-length;
	//merged->end = ((int)position)+length;
	merged->v = vTotal/((double)N);
	merged->focusRatio = ((double)N);
	delete [](merged->name);
	merged->name = new char[100];
	sprintf(merged->name, "Merged-%s-%lld-%d",merged->chr,position,N);
}

Hashtable* PeakLibrary::getOverlappingPeaks(PeakLibrary* peakset, char strand, int maxDistance,
													long long int &overlap) {

	overlap = 0;
	Hashtable* mapping = new Hashtable();
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		ChrPeaks* cp = (ChrPeaks*) chrs->search(keys[i]);
		ChrPeaks* cpSet = (ChrPeaks*) peakset->chrs->search(keys[i]);
		if (cpSet == NULL) {
			continue;
		}
		cp->getOverlappingPeaks(mapping, cpSet, strand,maxDistance,overlap);
		delete [](keys[i]);
	}
	delete []keys;

	return mapping;
}

void PeakLibrary::getMappabilityCount(char* uniqMapDirectory, int tagAdjust, char strand,
							long long int &totalPossible, long long int &totalMappable) {
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	fprintf(stderr, "\tCalculating mappability (tag adjustment at %d bp)...\n", tagAdjust);
	for (int i=0;i<chrs->total;i++) {
		ChrPeaks* cp = (ChrPeaks*) chrs->search(keys[i]);
		UniqMapChrs* umc = new UniqMapChrs(keys[i],uniqMapDirectory,0);
		totalPossible += umc->size*2;
		totalMappable += umc->numMappable;
		cp->getMappabilityCount(umc,tagAdjust,strand);
		delete umc;
		delete [](keys[i]);
	}
	delete []keys;
}



// class ChrPeaks -----------------------------------------------


ChrPeaks::ChrPeaks() {
	//peaks = new Peak*[PEAK_INC];
	peaks = NULL;
	peakList = NULL;
	maxPossible = PEAK_INC;
	numPeaks = 0;
}
ChrPeaks::~ChrPeaks() {
	if (peaks != NULL) delete []peaks;
	if (peakList != NULL) delete peakList;
	maxPossible = 0;
}
void ChrPeaks::addPeak(Peak* p) {
	if (peakList == NULL) {
		peakList = new LinkedList();
	}
	peakList->add(p);
}
void ChrPeaks::sort() {

	int numListPeaks = 0;
	Peak** array = NULL;
	if (peakList != NULL) {
		array = (Peak**) peakList->toArray(numListPeaks);
	}

	int newNumPeaks = numPeaks+numListPeaks;
	Peak** newpeaks = new Peak*[newNumPeaks];
	for (int i=0;i<numPeaks;i++) {
		newpeaks[i] = peaks[i];
	}
	for (int i=0;i<numListPeaks;i++) {
		newpeaks[i+numPeaks]=array[i];
	}
	delete []peaks;
	delete []array;
	delete peakList;
	peakList = NULL;

	peaks = newpeaks;
	numPeaks = newNumPeaks;
	if (numPeaks < 2) return;
	qsort(peaks,numPeaks,sizeof(Peak*),&cmpPeaks);
}
void ChrPeaks::checkForOutOfBoundsCoordinates(ChrTags* ct) {
	if (ct==NULL) return;
	//ct->getMaxPosition();
	//xxxxxxxxxxxxxxxxxxxxxxxx
	int minP = 1;
	int maxP = (int)ct->appearentSize;
	for (int i=0;i<numPeaks;i++) {
		if (peaks[i]->start < minP) peaks[i]->start = minP;
		if (peaks[i]->end < minP) peaks[i]->end = minP;
		if (peaks[i]->start < minP) peaks[i]->start = minP;
		if (peaks[i]->end < minP) peaks[i]->end = minP;

		if (peaks[i]->start > maxP) peaks[i]->start = maxP;
		if (peaks[i]->end > maxP) peaks[i]->end = maxP;
		if (peaks[i]->start > maxP) peaks[i]->start = maxP;
		if (peaks[i]->end > maxP) peaks[i]->end = maxP;
	}
}

//global clusterF**k variable
static Inttable* chrIndex = NULL;
int chrcmp(const void* chr1, const void* chr2) {
	if (chr1 == chr2) return 0;
	if (chr1 == NULL) return -1;
	if (chr2 == NULL) return 1;
	char* c1 = *((char**) chr1);	
	char* c2 = *((char**) chr2);	
	int sc = strcmp(c1,c2);
	if (sc == 0) return 0;

	if (chrIndex == NULL) {
		chrIndex = new Inttable(10000);
		char* tmp = new char[1000];
		int index = 1;
		for (int i=0;i<=100;i++) {
			sprintf(tmp,"chr%d",i);
			chrIndex->insert(index++,tmp);
			sprintf(tmp,"chr%d_random",i);
			chrIndex->insert(index++,tmp);
			sprintf(tmp,"chr%dL",i);
			chrIndex->insert(index++,tmp);
			sprintf(tmp,"chr%dL_random",i);
			chrIndex->insert(index++,tmp);
			sprintf(tmp,"chr%dR",i);
			chrIndex->insert(index++,tmp);
			sprintf(tmp,"chr%dR_random",i);
			chrIndex->insert(index++,tmp);
		}
		chrIndex->insert(index++,(char*)"chrI");
		chrIndex->insert(index++,(char*)"chrII");
		chrIndex->insert(index++,(char*)"chrIII");
		chrIndex->insert(index++,(char*)"chrIV");
		chrIndex->insert(index++,(char*)"chrV");
		chrIndex->insert(index++,(char*)"chrVI");
		chrIndex->insert(index++,(char*)"chrVII");
		chrIndex->insert(index++,(char*)"chrVIII");
		chrIndex->insert(index++,(char*)"chrIX");
		chrIndex->insert(index++,(char*)"chrU");
		chrIndex->insert(index++,(char*)"chrU_random");
		chrIndex->insert(index++,(char*)"chrX");
		chrIndex->insert(index++,(char*)"chrX_random");
		chrIndex->insert(index++,(char*)"chrXI");
		chrIndex->insert(index++,(char*)"chrXII");
		chrIndex->insert(index++,(char*)"chrXIII");
		chrIndex->insert(index++,(char*)"chrXIV");
		chrIndex->insert(index++,(char*)"chrXV");
		chrIndex->insert(index++,(char*)"chrXVI");
		chrIndex->insert(index++,(char*)"chrXVII");
		chrIndex->insert(index++,(char*)"chrXVIII");
		chrIndex->insert(index++,(char*)"chrXIX");
		chrIndex->insert(index++,(char*)"chrY");
		chrIndex->insert(index++,(char*)"chrY_random");
		chrIndex->insert(index++,(char*)"chrZ");
		chrIndex->insert(index++,(char*)"chrZ_random");
		chrIndex->insert(index++,(char*)"chrM");
		chrIndex->insert(index++,(char*)"chrM_random");
		chrIndex->insert(index++,(char*)"chrMT");
		chrIndex->insert(index++,(char*)"chrMT_random");
		chrIndex->insert(index++,(char*)"chrUn");
		chrIndex->insert(index++,(char*)"chrUn_random");
		chrIndex->insert(index++,(char*)"genome");
		chrIndex->insert(index++,(char*)"null");
		delete []tmp;
	}
	int i1 = chrIndex->search(c1);
	int i2 = chrIndex->search(c2);
	//if (i1 == EMPTY_INT) fprintf(stderr, "Couldn't find anything for %s\n", c1);
	//if (i2 == EMPTY_INT) fprintf(stderr, "Couldn't find anything for %s\n", c2);
	if (i1 != EMPTY_INT && i2 != EMPTY_INT) {
		if (i1 < i2) return -1;
		if (i1 > i2) return 1;
		return 0;
	}
	return sc;
}

int cmpPeaks(const void* a, const void* b) {
	char* ac = (*(Peak**)a)->chr;
	char* bc = (*(Peak**)b)->chr;
	int cc = chrcmp(&ac,&bc);
	if (cc != 0) return cc;
	int ap = (*(Peak**)a)->tagStart;
	int bp = (*(Peak**)b)->tagStart;
	if (ap < bp) return -1;
	if (ap > bp) return 1;
	int ad = (*(Peak**)a)->tagEnd;
	int bd = (*(Peak**)b)->tagEnd;
	if (ad < bd) return -1;
	if (ad > bd) return 1;
	char al = (*(Peak**)a)->strand;
	char bl = (*(Peak**)b)->strand;
	if (al < bl) return -1;
	if (al > bl) return 1;
	return 0;
}

void ChrPeaks::countPeakTagsLowMemory(Doubletable* results, ChrTags* ct, char direction, int mode) {

	ct->loadTags();

	int peakIndex = 0;
	//int tagIndex = 0;
	//establish when we can start forgetting about peaks
	// peaks are already sorted by their starting positions;
	int* finishedLength = new int[numPeaks];
	double* peakTotals = new double[numPeaks];
	double* peakPositions = new double[numPeaks];
	for (int i=0;i<numPeaks;i++) {
		peakTotals[i] = 0.0;
		peakPositions[i] = 0.0;
		if (i==0) {
			finishedLength[i] = peaks[i]->tagEnd;
		} else {
			if (finishedLength[i-1] > peaks[i]->tagEnd) {
				finishedLength[i] = finishedLength[i-1];
			} else {
				finishedLength[i] = peaks[i]->tagEnd;
			}
		}
	}

	//fprintf(stderr, "TotalPeaks=%d, totalTags=%d\n",numPeaks,ct->totalPositions);

	for (int i=0;i<ct->totalPositions;i++) {
		int p = ct->tags[i].p;
		int d = ct->tags[i].d;
		float v = ct->tags[i].v;

		if (p < peaks[peakIndex]->tagStart) continue;
		int stop=0;
		while (p > finishedLength[peakIndex]) {
			peakIndex++;
			if (peakIndex >= numPeaks) {
				stop = 1;
				break;
			}
		}
		if (stop) break;
		for (int j=peakIndex;j<numPeaks;j++) {
			if (direction == STRAND_SEPARATE && peaks[j]->strand != d) {
				continue;
			}
			if (p >= peaks[j]->tagStart && p <= peaks[j]->tagEnd) {
				peakTotals[j]+=v;
				peakPositions[j]+=1.0;
			}
			if (p <= peaks[j]->tagStart) {
				break;
			}
		}
	}

	for (int i=0;i<numPeaks;i++) {

		double total = peakTotals[i];
		double totalPosUsed = peakPositions[i];
		double totalBP = (double) peaks[i]->end-peaks[i]->start+1;
		
		if (mode == COUNT_MODE_TBP) {
			if (totalBP > 0) total /= totalBP;
		}
		if (mode == COUNT_MODE_RATIO) {
			if (totalPosUsed > 0.0) {
				total /= totalPosUsed;
			} else {
				total = EMPTY_DOUBLE;
			}
		}
		results->insert(total,peaks[i]->name);
	}
	delete []finishedLength;
	delete []peakTotals;
	delete []peakPositions;
	ct->freeTags();

}


void ChrPeaks::addTagLibrary(ChrTags* ct, int expIndex) {
	ct->loadTags();

	int peakIndex = 0;
	//int tagIndex = 0;
	//establish when we can start forgetting about peaks
	// peaks are already sorted by their starting positions;
	int* finishedLength = new int[numPeaks];
	for (int i=0;i<numPeaks;i++) {
		if (i==0) {
			finishedLength[i] = peaks[i]->tagEnd;
		} else {
			if (finishedLength[i-1] > peaks[i]->tagEnd) {
				finishedLength[i] = finishedLength[i-1];
			} else {
				finishedLength[i] = peaks[i]->tagEnd;
			}
		}
	}

	//fprintf(stderr, "TotalPeaks=%d, totalTags=%d\n",numPeaks,ct->totalPositions);

	for (int i=0;i<ct->totalPositions;i++) {
		int p = ct->tags[i].p;

		if (p < peaks[peakIndex]->tagStart) continue;
		int stop=0;
		while (p > finishedLength[peakIndex]) {
			peakIndex++;
			if (peakIndex >= numPeaks) {
				stop = 1;
				break;
			}
		}
		if (stop) break;
		for (int j=peakIndex;j<numPeaks;j++) {
			if (p >= peaks[j]->tagStart && p <= peaks[j]->tagEnd) {
				peaks[j]->addTag(&(ct->tags[i]),expIndex);
			}
			if (p <= peaks[j]->tagStart) {
				break;
			}
		}
	}
	delete []finishedLength;

	ct->freeTags();
}
void ChrPeaks::getMappabilityCount(UniqMapChrs* umc, int tagAdjust, char strand) {
	char s = 0;
	for (int i=0;i<numPeaks;i++) {
		if (peaks[i]->strand == STRAND_POSITIVE) {
			if (strand == STRAND_POSITIVE) s = STRAND_POSITIVE;
			else if (strand == STRAND_NEGATIVE) s = STRAND_NEGATIVE;
			else s = STRAND_BOTH;
		} else {
			if (strand == STRAND_POSITIVE) s = STRAND_NEGATIVE;
			else if (strand == STRAND_NEGATIVE) s = STRAND_POSITIVE;
			else s = STRAND_BOTH;
		}
		if (tagAdjust == 0) {
			peaks[i]->uniqMap = umc->countRegion(peaks[i]->start,peaks[i]->end,s);
		} else {
			int count = 0;
			if (s ==  STRAND_POSITIVE || s == STRAND_BOTH)
				count += umc->countRegion(peaks[i]->start-tagAdjust,peaks[i]->end-tagAdjust,STRAND_POSITIVE);
			if (s ==  STRAND_NEGATIVE || s == STRAND_BOTH)
				count += umc->countRegion(peaks[i]->start+tagAdjust,peaks[i]->end+tagAdjust,STRAND_NEGATIVE);
			
			peaks[i]->uniqMap = count;
		}
			
		//fprintf(stderr, "%s %d\n", peaks[i]->name, peaks[i]->uniqMap);
	}
}

void ChrPeaks::prioritizeAnnotations(PeakLibrary* annotations) {

	Peak** active = new Peak*[numPeaks*2];
	Peak** stillActive = new Peak*[numPeaks*2];
	Peak** tmpcopy = NULL;
	int numActive = 0;

	int numStillActive = 0;
	int skipFlag = 0;
	//int resortFlag = 0;
	for (int i=0;i<numPeaks;i++) {
		//int curStart = peaks[i]->start;
		//int curEnd = peaks[i]->end;
		//int curPriority = peaks[i]->priority;
//fprintf(stderr, "%s\t%s\t%d\n", peaks[i]->name, peaks[i]->data, peaks[i]->priority);

		numStillActive = 0;
		skipFlag = 0;
		int originalStart= peaks[i]->start;
		for (int j=0;j<numActive;j++) {
			if (j >= numPeaks*2) {
				fprintf(stderr, "Error - to many annotations - algorithm error\n");
			}
			//int start = active[j]->start;
			//int end = active[j]->end;
			//int priority = active[j]->priority;

			if (active[j]->end < originalStart) {
				if (active[j]->end-active[j]->start>=0) {
					annotations->addPeak(active[j]);
				}
				//stillActive[numStillActive++] = peaks[i];
			} else if (active[j]->end < peaks[i]->start) {
				stillActive[numStillActive++] = active[j];
			} else if (active[j]->end <= peaks[i]->end) {
				if (active[j]->start >= peaks[i]->start) {
					//contained within peaks[i]
					if (active[j]->priority < peaks[i]->priority) {
						//split it up
						Peak* sub1 = peaks[i]->copy();
						sub1->end = active[j]->start-1;
						stillActive[numStillActive++] = sub1;

						peaks[i]->start = active[j]->end+1;
						stillActive[numStillActive++] = active[j];
					} else {

						// active[j] peak is totally masked by new peak
					}
				} else {
					if (active[j]->priority < peaks[i]->priority) {
						peaks[i]->start = active[j]->end+1;
						stillActive[numStillActive++] = active[j];
					} else {
						active[j]->end = peaks[i]->start-1;
						stillActive[numStillActive++] = active[j];
					}
				}
			} else {
				if (active[j]->start < peaks[i]->start) {
					if (active[j]->priority < peaks[i]->priority) {
						//don't keep peaks[i]
						skipFlag = 1;
						stillActive[numStillActive++] = active[j];
					} else {
						//split it up active[j]
						Peak* sub1 = active[j]->copy();
						sub1->end = peaks[i]->start-1;
						stillActive[numStillActive++] = sub1;
						active[j]->start = peaks[i]->end+1;
						stillActive[numStillActive++] = active[j];
						//resortFlag = 1;
					}
				} else if (active[j]->start <= peaks[i]->end) {
					if (active[j]->priority < peaks[i]->priority) {
						peaks[i]->end = active[j]->start-1;
						stillActive[numStillActive++] = active[j];
					} else {
						active[j]->start = peaks[i]->end+1;
						stillActive[numStillActive++] = active[j];
					}
				} else {
					stillActive[numStillActive++] = active[j];
				}
			}
		}
		if (skipFlag==0) {
			stillActive[numStillActive++] = peaks[i];
		}
		tmpcopy = active;
		active = stillActive;
		stillActive = tmpcopy;
		numActive = numStillActive;
		qsort(active,numActive,sizeof(Peak*),&cmpPeaks);
	}
	for (int i=0;i<numActive;i++) {
		if (active[i]->end-active[i]->start>=0) {
			annotations->addPeak(active[i]);
		}
	}
	
	delete []active;
	delete []stillActive;
	
}

void ChrPeaks::mergeNearbyPeaks(PeakLibrary* merged, int maxDistance, char strand) {
	int count = 1;	
	int* mask = new int[numPeaks];
	for (int i=0;i<numPeaks;i++) mask[i]=0;
	
	for (int i=0;i<numPeaks;i++) {
		if (mask[i]) continue;
		int ptotal = 1;
		for (int j=1;j<numPeaks;j++) {
			if (strand == STRAND_SEPARATE && peaks[i]->strand != peaks[j]->strand) continue;
			if (maxDistance < 0) { //given overlap
				if (peaks[j]->start <= peaks[i]->end) {
					if (peaks[j]->end > peaks[i]->end) {
						peaks[i]->end = peaks[j]->end;
					}
					mask[j]=1;
					ptotal++;
				} else {
					break;
				}
			} else {
				if (fabs(peaks[j]->refPos-peaks[i]->refPos) < maxDistance) {
					if (peaks[j]->end > peaks[i]->end) {
						peaks[i]->end = peaks[j]->end;
					}
					mask[j]=1;
					ptotal++;
				} else {
					if (peaks[j]->start-peaks[i]->end > maxDistance) {
						break;
					}
				}
			}
		}
		peaks[i]->refPos = (peaks[i]->start+peaks[i]->end)/2;
		peaks[i]->focusRatio = (double)ptotal;
		if (ptotal>1) {
			sprintf(peaks[i]->name,"Merged-%s-%d",peaks[0]->chr,count++);
		}
		merged->addPeak(peaks[i]);
	}

	delete []mask;
}

void ChrPeaks::getAnnotationSizeTotals(Doubletable* sizeTotals) {
	for (int i=0;i<numPeaks;i++) {
		int s = peaks[i]->end-peaks[i]->start+1;
		if (s < 1) continue;
		if (peaks[i]->start < 0) continue;
		if (peaks[i]->end > 1400000000) continue;
		char* str = peaks[i]->data;
		double tsize = sizeTotals->search(str);
		if (tsize < EMPTY_DOUBLE_CHECK) {
			tsize = 0.0;
		}
		tsize += s;
		sizeTotals->insert(tsize, str);
	}
}

void ChrPeaks::annotatePeakLocations(ChrPeaks* annotations, FILE* annFile,Doubletable* stats) {


	int annIndex = 0;
	//int curPos = 0;
	for (int i=0;i<numPeaks;i++) {
		int p = peaks[i]->refPos;
		int lastEnd = -1;
		for (int j=annIndex;j<annotations->numPeaks;j++) {
			int s = annotations->peaks[j]->start;
			int e = annotations->peaks[j]->end;
			if (e > lastEnd) lastEnd = e;
			if (lastEnd < peaks[i]->start) annIndex++;
			if (p <= e && p >= s) {
				char* str = annotations->peaks[j]->data;
				double statCount = stats->search(str);
				if (statCount < EMPTY_DOUBLE_CHECK) {
					statCount = 0.0;
				}
				statCount += 1.0;
				stats->insert(statCount, str);
				if (annFile != NULL) {
					char* nn = annotations->peaks[j]->name;
					if (annotations->peaks[j]->ogname != NULL)
						nn = annotations->peaks[j]->ogname;
					
					fprintf(annFile,"%s\t%s\t%s\n", peaks[i]->name, str, nn);
				}
				break;
			}
		}
	}

}
void ChrPeaks::getOverlappingPeaks(Hashtable* mapping, ChrPeaks* peakset, char strand,int maxDistance, 
																				long long int &overlap) {

	int doneThroughPos = 0;
	int doneThroughNeg = 0;

	int sameFlag = 0;
	if (this == peakset) {
		sameFlag = 1;
	}
	int setIndex = 0;
	for (int i=0;i<numPeaks;i++) {
		//int done = 0;
		int safeDistance = maxDistance;
		if (maxDistance < 1) {
			safeDistance = 1;
		}
		while (setIndex < peakset->numPeaks &&
						peakset->peaks[setIndex]->end+safeDistance < peaks[i]->start) {
			setIndex++;
		}
		if (setIndex == peakset->numPeaks) {
			//done=1;
			break;
		}

		PeakMapping *pm = NULL;

		int doneThrough = doneThroughPos;
		if (strand != STRAND_BOTH && peaks[i]->strand == STRAND_NEGATIVE) {
			doneThrough = doneThroughNeg;
		}

		for (int j=setIndex;j<peakset->numPeaks;j++) {
			if (peakset->peaks[j]->start-safeDistance > peaks[i]->end) {
				break;
			}

			if (strand != STRAND_BOTH && peaks[i]->strand != peakset->peaks[j]->strand) {
				continue;
			}

			int diff = peakset->peaks[j]->refPos - peaks[i]->refPos;
			if (diff < 0) diff *= -1;
			int samePeak = 0;
			if (peakset->peaks[j] == peaks[i]) samePeak = 1;

			int add=0;
			if (maxDistance > -1) {
				if (diff < maxDistance) {
					add=1;
				}
			} else {
/*fprintf(stderr, "Comparing:. %d %s vs. %d %s\n", i,peaks[i]->name,j,peakset->peaks[j]->name);
peaks[i]->print(stderr);
peakset->peaks[j]->print(stderr);
fprintf(stderr, "Donethrough=%d\n", doneThrough);*/
				int cend = peakset->peaks[j]->end;
				if (peaks[i]->end < cend) cend = peaks[i]->end;
				int cstart = peakset->peaks[j]->start;
				if (peaks[i]->start > cstart) cstart = peaks[i]->start;
				if (doneThrough > cstart) cstart = doneThrough;

				if (peakset->peaks[j]->start >= peaks[i]->start &&
							peakset->peaks[j]->start <= peaks[i]->end) {
					if (!samePeak && cend > doneThrough) {
						overlap += cend-cstart+1;
						doneThrough = cend;
					}
					add=1;
				} else if (peakset->peaks[j]->end >= peaks[i]->start &&
							peakset->peaks[j]->end <= peaks[i]->end) {
					if (!samePeak && cend > doneThrough) {
						overlap += cend-cstart+1;
						doneThrough = cend;
					}
					add=1;
				} else if (peakset->peaks[j]->start <= peaks[i]->start && 
							peakset->peaks[j]->end >= peaks[i]->end) {
					if (!samePeak && cend > doneThrough) {
						overlap += cend-cstart+1;
						doneThrough = cend;
					}
					add=1;
				}
//fprintf(stderr, "afterhrough=%d\n\n", doneThrough);
			}
			if (add ==1) {
				if (sameFlag && peaks[i] == peakset->peaks[j]) {
			//		continue;
				}
				if (pm == NULL) {
					pm = new PeakMapping(peaks[i]);
				}
				pm->add(peakset->peaks[j],diff);
			}
		}

		doneThrough = peaks[i]->end;
		if (strand != STRAND_BOTH) {
			if (peaks[i]->strand == STRAND_NEGATIVE) {
				doneThroughNeg = doneThrough;
			} else {
				doneThroughPos = doneThrough;
			}
		} else {
			doneThroughPos = doneThrough;
		}

		if (pm != NULL) {
			pm->calculateCoverage();
			mapping->insert(pm, peaks[i]->name);
		}
	}
}

PeakMapping::PeakMapping(Peak* p) {
	source = p;
	peaks = NULL;
	distance = NULL;
	numPeaks = 0;
	bp = 0;
}
PeakMapping::~PeakMapping() {
	//fprintf(stderr, "numPeaks=%d\n", numPeaks);
	if (peaks != NULL) {
		delete []peaks;
	}
	if (distance != NULL) {
		delete []distance;
	}
}
void PeakMapping::add(Peak* p , int dist) {
	Peak** pp = new Peak*[1+numPeaks];
	int* d = new int[1+numPeaks];
	for (int i=0;i<numPeaks;i++) {
		pp[i] = peaks[i];
		d[i] = distance[i];
	}
	if (peaks != NULL) {
		delete []peaks;
		delete []distance;
	}
	pp[numPeaks] = p;
	d[numPeaks] = dist;
	peaks = pp;
	distance = d;
	numPeaks++;
}
void PeakMapping::calculateCoverage() {
	bp = 0;
	if (numPeaks == 0) return;
	qsort(peaks, numPeaks, sizeof(Peak*), &cmpPeaks);
	//sorted by starting positions
	int* pos = new int[2*numPeaks+2];
	for (int i=0;i<2*numPeaks+2;i++) pos[i]=-1;
	pos[0] = peaks[0]->start;
	pos[1] = peaks[0]->end;
	int currentIndex = 0;
	for (int i=1;i<numPeaks;i++) {	
		if (pos[currentIndex+1] < peaks[i]->start) {
			currentIndex+=2;
			pos[currentIndex] = peaks[i]->start;
			pos[currentIndex+1] = peaks[i]->end;
		} else {
			if (pos[currentIndex+1] < peaks[i]->end) {
				pos[currentIndex+1] = peaks[i]->end;
			}
		}
	}
	int start = source->start;
	int end = source->end;
	for (int i=0;i<numPeaks*2;i++) {
		if (pos[i] == -1) break;
		if (pos[i] > end) break;
		if (pos[i+1] < start) continue;
		if (pos[i] < start && pos[i+1] > end) {
			bp = end-start;
			break;
		}
		if (pos[i] < start && pos[i+1] < end) {
			bp += pos[i+1]-start;
		} else if (pos[i] > start && pos[i+1] < end) {
			bp += pos[i+1]-pos[i];
		} else if (pos[i] > start && pos[i+1] > end) {
			bp += end-pos[i];
			break;
		}
	}
	delete []pos;
	//fprintf(stderr, "bp=%d (%d)\n",bp,numPeaks);
}


long long int ChrPeaks::calculateCoverage() {

	long long int bp = 0;	

	int* pos = new int[2*numPeaks+2];
	for (int i=0;i<2*numPeaks+2;i++) pos[i]=-1;
	pos[0] = peaks[0]->start;
	pos[1] = peaks[0]->end;
	int currentIndex = 0;
	for (int i=1;i<numPeaks;i++) {	
		if (pos[currentIndex+1] < peaks[i]->start) {
			currentIndex+=2;
			pos[currentIndex] = peaks[i]->start;
			pos[currentIndex+1] = peaks[i]->end;
		} else {
			if (pos[currentIndex+1] < peaks[i]->end) {
				pos[currentIndex+1] = peaks[i]->end;
			}
		}
	}
	for (int i=0;i<currentIndex+2;i+=2) {
		if (pos[i] == -1) break;
		bp += pos[i+1]-pos[i];
	}
	delete []pos;

	return bp;
}

void ChrPeaks::print(FILE* fp) {
	for (int i=0;i<numPeaks;i++) {
		peaks[i]->printAnnotation(fp);
	}
}

void ChrPeaks::stitchRegions(PeakLibrary* regions, int maxDistance,int mode) {

	double* totals = new double[numPeaks];
	int* totalDist = new int[numPeaks];
	int* uniqmap = new int[numPeaks];
	for (int i=0;i<numPeaks;i++) totals[i]=0.0;
	for (int i=0;i<numPeaks;i++) totalDist[i]=0;
	for (int i=0;i<numPeaks;i++) uniqmap[i]=0;
	int totalIndex = 0;
	int maxNumEstDist = 500;

	for (int k=0;k<2;k++) {

		for (int index=0;index<numPeaks;index++) {
			int i=index;
			int peakLimit = numPeaks;
			if (k==1) {
				i = numPeaks-index-1;
				peakLimit = i+i+1;
			}

			if (peaks[i]->strand != k) continue;
			int curStart = peaks[i]->start;
			int curEnd = peaks[i]->end;
			//char curStrand = peaks[i]->strand;
			float v = peaks[i]->v;
			int lastAdded = i;
			totalIndex = 0;
			totals[totalIndex] = v;
			uniqmap[totalIndex] = peaks[i]->uniqMap;
			if (uniqmap[totalIndex] < 0) uniqmap[totalIndex] = peaks[i]->end - peaks[i]->start;
			if (uniqmap[totalIndex] < 10) uniqmap[totalIndex] = 10;
			totalDist[totalIndex] = peaks[i]->refPos;

			for (int indexj=i+1;indexj<peakLimit;indexj++) {
				int j=indexj;
				if (k==1) {
					j=i-(indexj-i);
				}
				if (peaks[j]->strand != k) continue;
				if (peaks[j]->refPos - peaks[lastAdded]->refPos < maxDistance
				 		&& peaks[j]->refPos - peaks[lastAdded]->refPos > -1*maxDistance) {

					int add = 0;
					if (mode == REGION_MODE_GROSEQ) {
						double score = 0.0;
						double NN = 0.0;
						//int diff = 0;
						int dindex = totalIndex;
						while (dindex >= 0 && abs(totalDist[totalIndex]-totalDist[dindex])
													< maxNumEstDist) {
							//score+=totals[dindex];
							double curScore = totals[dindex]/((double)uniqmap[dindex]);
//							if (curScore > score) {
//								score = curScore;
//							}
							score+=curScore;
							NN+=1.0;
							//diff = totalDist[totalIndex]-totalDist[dindex];
							dindex--;
						}
						if (NN > 0.0) {
							//score /= NN;
							score/=NN;
						} else {
							score = 1;
						}
						
						double pmap = (double)peaks[j]->uniqMap;	
						if (pmap < 0) pmap= (double)(peaks[j]->end-peaks[j]->start);
						if (uniqmap[totalIndex] < 10) uniqmap[totalIndex] = 10;

						double fold = (peaks[j]->v/pmap)/(score);

						//fprintf(stderr, "%d\t%lf\t%lf\t%lf\t%lf\n", peaks[j]->refPos,fold,peaks[j]->v,pmap,score);
						if (fold > REGION_GROSEQ_FOLDCHANGE) {
							break;
						} else {
							add=1;
						}
					} else if (mode == REGION_MODE_HISTONE) {
						add=1;
					}

					if (add) {
						totalIndex++;
						if (k==1) {
							curStart = peaks[j]->start;
						} else {
							curEnd = peaks[j]->end;
						}
						totalDist[totalIndex] = peaks[j]->refPos;
						totals[totalIndex] = peaks[j]->v;
						uniqmap[totalIndex] = peaks[j]->uniqMap;
						if (uniqmap[totalIndex] < 0) uniqmap[totalIndex] = peaks[j]->end - peaks[j]->start;
						if (uniqmap[totalIndex] < 10) uniqmap[totalIndex] = 10;
						v += peaks[j]->v;
						lastAdded = j;
					}
				} else {
					break;
				}
			}
			peaks[i]->start = curStart;
			peaks[i]->end = curEnd;
			//fprintf(stderr, "done = %d\t%d\n", curStart, curEnd);
			peaks[i]->focusRatio = curEnd-curStart;
			peaks[i]->v = v;
			regions->addPeak(peaks[i]);
			index = lastAdded;
			if (k==1) {
				index = numPeaks-lastAdded-1;
			}
		}
	}
	delete []totals;
	delete []totalDist;
}

// class Peaks -----------------------------------------------


Peak::Peak() {
	init();
}
void Peak::init() {
	name = NULL;
	chr = NULL;
	refPos = 0;
	start = 0;
	end = 0;
	tagStart = 0;
	tagEnd = 0;
	fixedFlag = 0;
	priority = 0;
	strand = STRAND_POSITIVE;
	seq = NULL;
	uniqMap = 0;
	v = 0.0;
	index = 0;
	focusRatio = 0.0;
	data = NULL;
	ogname = NULL;
	snps = NULL;
	exps = NULL;
	numTags = NULL;
	maxTags = NULL;
	numExps=0;
}
Peak::Peak(char* newname, char* originalName, char* newchr, int newstart, int newend, int newRef, char dir, 
					float value, float ratio, char* otherdata, int mappability,unsigned int newpriority) {
	if (newname == NULL) {
		int L = strlen(newchr)+13+1+6;
		name = new char[L];
		sprintf(name,"%s-%d-%d",newchr,newstart,dir);
		//fprintf(stderr, "name=%s \t %d\n",name,L);
	} else {
		name = new char[strlen(newname)+1];
		strcpy(name,newname);
	}
	ogname = NULL;
	if (originalName != NULL) {
		ogname = new char[strlen(originalName)+1];
		strcpy(ogname,originalName);
	}
	chr = new char[strlen(newchr)+1];
	strcpy(chr,newchr);
	start = newstart;
	end = newend;
	tagStart = newstart;
	tagEnd = newend;
	refPos = newRef;
	strand = dir;
	v = value;
	seq = NULL;
	priority = newpriority;
	focusRatio = ratio;
	fixedFlag =0;
	uniqMap = mappability;
	exps = NULL;
	numTags = NULL;
	maxTags = NULL;
	numExps = 0;
	snps = NULL;
	if (otherdata != NULL) {
		data = new char[strlen(otherdata)+1];
		strcpy(data,otherdata);
	} else {
		data = NULL;
	}
	if (refPos == NULL_REF) {
		refPos = (start+end)/2;
	}
}
Peak::~Peak() {
	if (name != NULL) delete []name;
	if (ogname != NULL) delete []ogname;
	if (chr != NULL) delete []chr;
	if (data != NULL) delete []data;
	if (exps != NULL) {
		for (int i=0;i<numExps;i++) {
			delete exps[i];
		}
		delete []exps;
	}
	if (numTags != NULL) delete []numTags;
	if (maxTags != NULL) delete []maxTags;
	if (seq != NULL) delete []seq;
	if (snps != NULL) delete snps;
}
Peak* Peak::copy() {
	Peak* np = new Peak(name,ogname,chr,start,end,refPos,strand,v,focusRatio,data,uniqMap,priority);
	np->uniqMap = uniqMap;
	/*if (exps != NULL & numExps>0) {
		np->numExps = numExps;
		np->exps= new Tag*[numExps];
		np->numTags = new int[numExps];
		np->maxTags = new int[numExps];
		for (int i=0;i<numExps;i++) {
			np->numTags[i] = numTags[i];
			np->maxTags[i] = maxTags[i];
			np->exps[i] = new Tag[maxTags[i]];
			for (int j=0;j<numTags[i];j++) {
				np->exps[i][j].copy(&(exps[i][j]));
			}
		}
	}*/
	return np;
}
void Peak::printSNPs(FILE* fp) {
	fprintf(fp, "%s\t",name);
	if (snps != NULL) snps->print(fp,strand);
	fprintf(fp, "\n");
}
void Peak::printSNPtotals(FILE* fp) {
	fprintf(fp, "%s\t",name);
	if (snps != NULL) snps->printTotals(fp);
	fprintf(fp, "\n");
}
void Peak::print(FILE* fp) {
	char dir = '+';
	if (strand == 1) dir = '-';
	if (v < 1.0) {
		fprintf(fp, "%s\t%s\t%d\t%d\t%c\t%.3f\t%.3f",name,chr,start,end,dir,v,focusRatio);
	} else if (v < 10.0) {
		fprintf(fp, "%s\t%s\t%d\t%d\t%c\t%.2f\t%.3f",name,chr,start,end,dir,v,focusRatio);
	} else {
		fprintf(fp, "%s\t%s\t%d\t%d\t%c\t%.1f\t%.3f",name,chr,start,end,dir,v,focusRatio);
	}
	//fprintf(fp, "\t%d", uniqMap);
	if (data != NULL) {
		fprintf(fp, "\t%s", data);
	}
	fprintf(fp, "\n");
}
void Peak::printGTF(FILE* fp) {
	char dir = '+';
	if (strand == 1) dir = '-';
	fprintf(fp, "%s\thomer\texon\t%d\t%d\t%lf\t%c\t.\tgene_id \"%s\"; transcript_id \"%s\"\n",
																chr,start,end,v,dir,name,name);
}
void Peak::printAnnotation(FILE* fp) {
	char dir = '+';
	if (strand == 1) dir = '-';
	fprintf(fp, "%s\t%s\t%d\t%d\t%c",name,chr,start,end,dir);
	if (data != NULL) {
		fprintf(fp, "\t%s", data);
	}
	fprintf(fp, "\t%d", priority);
	fprintf(fp, "\n");
}
void Peak::setSeq(char* str) {
	if (seq != NULL) {
		delete []seq;
	}
	seq = new char[strlen(str)+1];
	strcpy(seq,str);
}
void Peak::setPeakSize(int size) {
	int mid = (start+end)/2;
	int halfSize = size/2;
	start = mid-halfSize;
	end = mid+halfSize;
}

void Peak::addExp() {
	Tag** newexps = new Tag*[numExps+1];
	int* newnumTags = new int[numExps+1];
	if (exps != NULL) {
		for (int i=0;i<numExps;i++) {
			newexps[i] = exps[i];
			newnumTags[i] = numTags[i];
		}
		delete []exps;
		delete []numTags;
	}

	//we're going to hide this in the variable for now...
	LinkedList* linkedlist = new LinkedList();
	newexps[numExps]= (Tag*)linkedlist;
	newnumTags[numExps]=0;

	numTags = newnumTags;
	exps = newexps;

	numExps++;
}

void Peak::addData(char* str) {
	if (str == NULL) return;
	char* olddata = data;
	int len = 0;
	if (olddata != NULL) len = strlen(olddata)+1;	
	data = new char[len + strlen(str)+2];
	data[0] = '\0';
	if (olddata != NULL) {
		strcpy(data, olddata);
		strcat(data, "\t");
		delete []olddata;
	}
	strcat(data, str);
}

void Peak::addTag(Tag* t, int expIndex) {

	Tag* nt = new Tag();
	nt->copy(t);
	if (strand == 0) {
		nt->p = nt->p - refPos;
	} else {
		nt->p = refPos - nt->p;
		if (nt->d == 0) {
			nt->d = 1;
		} else {
			nt->d = 0;
		}
	}
	LinkedList* linkedlist = (LinkedList*) exps[expIndex];
	linkedlist->add(nt);
	//fprintf(stderr, "+ %.1f ",nt->v);
	//numTags[expIndex]++;
}
void Peak::setPeakTagSizeFixed(int startOffset,int endOffset) {
	setOffset(0);
	fixedFlag = 1;
	if (strand == 0) {
		tagStart = start + startOffset;
		tagEnd = end + endOffset;
		refPos = start;
	} else {
		tagStart = start - endOffset;
		tagEnd = end - startOffset;
		refPos = end;
	}
}
void Peak::setPeakTagSizeRefPos(int newOffset,int startOffset,int endOffset) {
	setOffset(newOffset);
	if (strand ==0) {
		tagStart = refPos + startOffset;
		tagEnd = refPos + endOffset;
	} else {
		tagStart = refPos - endOffset;
		tagEnd = refPos - startOffset;
	}
}

void Peak::setOffset(int newoffset) {
	if (newoffset == NULL_OFFSET) {
		//refPos = (int)floor(((float)(start+end))/2.0);
		refPos = (start+end)/2;
		return;
	}
	if (strand == 0) {
		refPos = start - newoffset;
	} else {
		refPos = end + newoffset;
	}
}
void Peak::sortTags(int expIndex) {

	LinkedList* list = (LinkedList*)exps[expIndex];
	Tag** array = (Tag**)list->toArray(numTags[expIndex]);
	Tag* tagset = new Tag[numTags[expIndex]];
	delete list;
	exps[expIndex] = tagset;
	for (int i=0;i<numTags[expIndex];i++) {
		exps[expIndex][i].copy(array[i]);
		delete array[i];
	}
	delete []array;

	for (int i=0;i<numExps;i++) {
		if (i==expIndex || expIndex == ALL_PEAK_EXPS) {
			qsort(exps[i],numTags[i],sizeof(Tag),&cmpTags);
		}
	}
}

double Peak::countTags(int expIndex, int regionStart, int regionEnd, char direction, int mode) {

	int startPos = regionStart;
	int endPos = regionEnd;
	if (fixedFlag) endPos = (end-start)+endPos;

	if (direction == STRAND_SEPARATE) {
		direction = 0;
	}

	double total = 0;
	int totalBP = 0;
	double totalPosUsed = 0.0;

	for (int i=0;i<numExps;i++) {
		if (i == expIndex || expIndex == ALL_PEAK_EXPS) {
			totalBP += endPos-startPos;
			for (int j=0;j<numTags[i];j++) {
				if (exps[i][j].p > endPos) break;
				if (exps[i][j].p >= startPos) {
					if (direction == BOTH_STRANDS || exps[i][j].d == direction) {
						total += exps[i][j].v;
						totalPosUsed += 1.0;
					}
				}
			}
		}
	}

	if (mode == COUNT_MODE_TBP) {
		if (totalBP > 0) total /= (float)totalBP;
	}
	if (mode == COUNT_MODE_RATIO) {
		if (totalPosUsed > 0.0) {
			total /= totalPosUsed;
		} else {
			total = EMPTY_DOUBLE;
		}
	}

	return total;
}

void Peak::printRelativeTags(FILE* fp, int expIndex, int regionStart, int regionEnd,int mode) {

	int startPos = regionStart;
	int endPos = regionEnd;
	if (fixedFlag) endPos = (end-start)+endPos;

	int totalBP = 0;
	int started = 0;
	if (mode == OUTPUT_MODE_PEAKTAGS) {
		fprintf(fp,"%s\t",name);
	}

	for (int i=0;i<numExps;i++) {
		if (i == expIndex || expIndex == ALL_PEAK_EXPS) {
			totalBP += endPos-startPos;
			for (int j=0;j<numTags[i];j++) {
				if (exps[i][j].p > endPos) break;
				if (exps[i][j].p >= startPos) {
					//fprintf(stderr, "%d\t%d\t%f\n", j, exps[i][j].p, exps[i][j].v);

					if (mode == OUTPUT_MODE_PEAKTAGS) {
						float pstrand = -1e20;
						float nstrand = -1e20;
						if (exps[i][j].d == 0) {
							pstrand = exps[i][j].v;
						} else {
							nstrand = exps[i][j].v;
						}
						while (j+1<numTags[i]) {
							if (exps[i][j].p == exps[i][j+1].p) {
								if (exps[i][j+1].d == 0) {
									if (pstrand > -1e19) pstrand += exps[i][j+1].v;
									else pstrand = exps[i][j+1].v;
								} else {
									if (nstrand > -1e19) nstrand += exps[i][j+1].v;
									else nstrand = exps[i][j+1].v;
								}
								j++;
							} else {
								break;
							}
						}
						if (started) fprintf(fp,",");
						started=1;
						fprintf(fp,"%d=",exps[i][j].p);
						if (pstrand > -1e19) {
							fprintf(fp,"%.1f",pstrand);
						} else {
							fprintf(fp,"NA");
						}
						fprintf(fp,"|");
							if (nstrand > -1e19) {
							fprintf(fp,"%.1f",nstrand);
						} else {
							fprintf(fp,"NA");
						}
					} else {
						char dd = 0;
						int pos = 0;
						if (strand==0) {
							pos = exps[i][j].p+refPos;
							if (exps[i][j].d == 1) {
								dd=1;
							}
						} else {
							pos = -1*exps[i][j].p+refPos;
							if (exps[i][j].d == 0) {
								dd=1;
							} else {
								dd=0;
							}
						}
						fprintf(fp,"\t%s\t%d\t%d\t%.1f\t%d\n",chr,pos,dd,exps[i][j].v,exps[i][j].len);	
					}
				}
			}
		}
	}
	if (mode == OUTPUT_MODE_PEAKTAGS) {
		fprintf(fp,"\n");
	}
}

void Peak::analyzeReadAutocorrelation(int expIndex,int maxSize,char strand,int mode) {

	double *ac = new double[maxSize+1];
	for (int i=0;i<=maxSize;i++) ac[i] = 0.0;
	double acTotal = 0.0;
	for (int i=0;i<numTags[expIndex];i++) {
		Tag* t = &(exps[expIndex][i]);
		if (strand == STRAND_SEPARATE && t->d == 1) continue;

		for (int j=i;j<numTags[expIndex];j++) {
			Tag* t2 = &(exps[expIndex][j]);
			if (strand == STRAND_SEPARATE && t2->d == 1) continue;
			int diff = abs(t->p - t2->p);
			double v = t->v*t2->v;
			acTotal+=v;
			if (diff > maxSize) diff = maxSize;
			ac[diff] += v;
		}
	}

	if (acTotal < 1.0) acTotal = 1.0;

	char* infoStr = new char[1000];
	infoStr[0]='\0';
	for (int i=0;i<=maxSize;i++) {
		sprintf(infoStr,"%le",ac[i]/acTotal);
		addData(infoStr);
	}
	delete []infoStr;

}
void Peak::analyzeTSSpattern(int expIndex,int fragLength,char strandInfo,int mode) {

	double dispersion = calculateTSSdispersion(expIndex,strandInfo, DEFAULT_TSS_DISPERSION_SIZE);
	double periodic  = calculateTSSperiodic(expIndex,strandInfo, DEFAULT_TSS_DISPERSION_SIZE);
	char* infoStr = new char[1000];
	sprintf(infoStr,"%lf\t%lf",dispersion,periodic);
	addData(infoStr);
	delete []infoStr;
}
double Peak::calculateTSSperiodic(int expIndex,int strandInfo, int focusRadius) {

	double totalReads = 0;
	double periodReads = 0;
	for (int i=0;i<numTags[expIndex];i++) {
		Tag* t = &(exps[expIndex][i]);
		if (strandInfo == STRAND_SEPARATE && t->d == 1) continue;
		if (t->p <= focusRadius && t->p >= -1*focusRadius) continue;
		int mod = abs(t->p) % 10;
		if (mod <= 2 || mod >= 8) periodReads += t->v;
		totalReads += t->v;
	}

	double ratio = 0.0;
	if (totalReads > 0.0) {
		ratio = periodReads / totalReads;
	}
	return ratio;
}
double Peak::calculateTSSdispersion(int expIndex,int strandInfo, int focusRadius) {

	double totalReads = 0;
	double focusReads = 0;
	for (int i=0;i<numTags[expIndex];i++) {
		Tag* t = &(exps[expIndex][i]);
		if (strandInfo == STRAND_SEPARATE && t->d == 1) continue;
		if (t->p <= focusRadius && t->p >= -1*focusRadius) {
			focusReads += t->v;
		}
		totalReads += t->v;
	}

	double ratio = 0.0;
	if (totalReads > 0.0) {
		ratio = focusReads / totalReads;
	}
	return ratio;
}

void Peak::centerPeak(int expIndex,int fragLength,char strandInfo) {

	int numCovTags=0;
	Tag* coverageTags = getCoverageTags(expIndex,fragLength,numCovTags,strandInfo);
	if (numCovTags < 1) {
		focusRatio = 0.0;
		return;
	}


	float centerSum=0;
	float centerN = 0;

	float max = 0;	
	for (int i=0;i<numCovTags;i++) {
		//if (coverageTags[i].v >= max) {
		if (coverageTags[i].v > max) {
			int diff = 1;
			int p = coverageTags[i].p;
			if (i < numCovTags-1) {
				diff = coverageTags[i+1].p - coverageTags[i].p;
				p += diff/2;
			}
			if (max < coverageTags[i].v) {
				max = coverageTags[i].v;
				centerSum = diff*p; 
				centerN = diff;
			} else {
				centerSum += diff*p; 
				centerN += diff;
			}
		}
	}
	delete []coverageTags;

	int centerP = refPos;
	int offset = 0;
	if (centerN ==0) {
		//fprintf(stderr, "Prolem\n");
	} else {
		centerP = (int)(((double)centerSum)/((double)centerN));
		offset = centerP;
		if (strand == STRAND_POSITIVE) {
			refPos += offset;
		} else {
			refPos -= offset;
		}
		int halfPeakSize = (int)((end-start)/2.0);
		start = refPos - halfPeakSize;
		end = refPos + halfPeakSize;
	}
	
	double totalTags=1.0;
	double goodTags=0.0;
	for (int i=0;i<numTags[expIndex];i++) {
		if (exps[expIndex][i].d == 0) {
			if (exps[expIndex][i].p-PEAKRATIO_TOLERANCE_BP < offset) {
				goodTags += exps[expIndex][i].v;
			}
		} else {
			if (exps[expIndex][i].p+PEAKRATIO_TOLERANCE_BP > offset) {
				goodTags += exps[expIndex][i].v;
			}
		}
		totalTags += exps[expIndex][i].v;
	}	

	if (totalTags > 0) {
		focusRatio = goodTags/totalTags;
	} else {
		focusRatio = 0.0;
	}

}


double Peak::scoreNFR(int expIndex, int fragLength, int nfrSize, int nucSize, char strandInfo) {
//void Peak::scoreNFR(int expIndex,int fragLength,char strandInfo,int nfrSize, int nucSize) {
	int numCovTags=0;
	Tag* coverageTags = getCoverageTags(expIndex,fragLength,numCovTags,strandInfo);

	if (numCovTags < 1) return 0;

	int poffset = 0;

	int nfrSizeHalf = nfrSize/2;
	int covOffset = -1*(nucSize+nfrSizeHalf)+coverageTags[0].p;
	int covEnd = coverageTags[numCovTags-1].p - covOffset + (nucSize+nfrSizeHalf);

	double* cov = new double[covEnd+1];
	int index = 0;
	for (int i=0;i<=covEnd;i++) {
		int p = i+covOffset;
		if (index > numCovTags-1) {
			cov[i] = 0;
			continue;
		}
		if (p < coverageTags[0].p) {
			cov[i] = 0.0;
		} else if (p < coverageTags[index+1].p) {
			cov[i] = coverageTags[index].v;
		} else {
			index++;
			if (index > numCovTags-1) {
				cov[i] = 0.0;
			} else {
				cov[i] = coverageTags[index].v;
			}
		}
	}
		

	//initialize Scores
	double nuc1Score = 0.0;
	double nuc2Score = 0.0;
	double nfrScore = 0.0;

	for (int i=0-poffset-nucSize-nfrSizeHalf;i<0-poffset-nfrSizeHalf;i++) nuc1Score += cov[i-covOffset];
	for (int i=0-poffset-nfrSizeHalf;i<0-poffset+nfrSizeHalf;i++) nfrScore += cov[i-covOffset];
	for (int i=0-poffset+nfrSizeHalf;i<0-poffset+nfrSizeHalf+nucSize;i++) nuc2Score += cov[i-covOffset];

	double maxNucScore = nuc1Score;
	if (nuc2Score > nuc1Score) maxNucScore = nuc2Score;
	maxNucScore /= nucSize;
	if (maxNucScore < 1.0) maxNucScore = 1.0;
	
	
	double scoreValue = (nuc1Score+nuc2Score)/(2*(double)nucSize)-nfrScore/((double)nfrSize);
	scoreValue /= maxNucScore;

	delete []coverageTags;
	delete []cov;

	return scoreValue;
}



void Peak::centerNFR(int expIndex,int fragLength,char strandInfo,int nfrSize, int nucSize) {
	int numCovTags=0;
	Tag* coverageTags = getCoverageTags(expIndex,fragLength,numCovTags,strandInfo);

	if (numCovTags < 1) {
		focusRatio = 0.0;
		return;
	}

	int nfrSizeHalf = nfrSize/2;
	int covOffset = -1*(nucSize+nfrSizeHalf)+coverageTags[0].p;
	int covEnd = coverageTags[numCovTags-1].p - covOffset + (nucSize+nfrSizeHalf);

	double* cov = new double[covEnd+1];
	int index = 0;
	for (int i=0;i<=covEnd;i++) {
		int p = i+covOffset;
		if (index > numCovTags-1) {
			cov[i] = 0;
			continue;
		}
		if (p < coverageTags[0].p) {
			cov[i] = 0.0;
		} else if (p < coverageTags[index+1].p) {
			cov[i] = coverageTags[index].v;
		} else {
			index++;
			if (index > numCovTags-1) {
				cov[i] = 0.0;
			} else {
				cov[i] = coverageTags[index].v;
			}
		}
	}
		

	//initialize Scores
	double nuc1Score = 0;
	double nuc2Score = 0;
	double nfrScore = 0;
	for (int i=0;i<nucSize;i++) nuc1Score += cov[i];
	for (int i=nucSize;i<nucSize+nfrSize;i++) nfrScore += cov[i];
	for (int i=nucSize+nfrSize;i<2*nucSize+nfrSize;i++) nuc2Score += cov[i];
	
	int nuc1edge1Offset = -nucSize-nfrSizeHalf;
	int nuc1edge2Offset = -nfrSizeHalf;
	int nuc2edge1Offset = nfrSizeHalf;
	int nuc2edge2Offset = nucSize+nfrSizeHalf;

	double centerSum=0;
	double centerN = 0;
	double max = -1e10;
	for (int i=nucSize+nfrSizeHalf;i<covEnd-(nucSize+nfrSizeHalf);i++) {
		nuc1Score -= cov[i+nuc1edge1Offset];
		nuc1Score += cov[i+nuc1edge2Offset];
		nuc2Score -= cov[i+nuc2edge1Offset];
		nuc2Score += cov[i+nuc2edge2Offset];
		nfrScore -= cov[i+nuc1edge2Offset];
		nfrScore += cov[i+nuc2edge1Offset];
		double ss = (nuc1Score+nuc2Score)/(2*(double)nucSize)-nfrScore/((double)nfrSize);

		if (ss >= max) {
			int p = i+covOffset;
			if (max < ss) {
				max = ss;
				centerSum = (double)p;
				centerN = 1.0;
			} else {
				centerSum += (double)p; 
				centerN += 1.0;
			}
		}
//		fprintf(stderr, "%d\t%lf\t%lf\n", i+covOffset, cov[i],ss);
	}
	focusRatio = max;
	delete []coverageTags;
	delete []cov;

	int centerP = refPos;
	int offset = 0;
	if (centerN < 0.5) {
		//fprintf(stderr, "Prolem\n");
	} else {
		centerP = (int)(((double)centerSum)/((double)centerN));
		offset = centerP;
		if (strand == STRAND_POSITIVE) {
			refPos += offset;
		} else {
			refPos -= offset;
		}
		int halfPeakSize = (int)((end-start)/2.0);
		start = refPos - halfPeakSize;
		end = refPos + halfPeakSize;
	}
//print(stderr);
//exit(0);
}





Tag* Peak::getCoverageTags(int expIndex, int fragLength, int& coveragePositions, char strandInfo) {

	coveragePositions = numTags[expIndex]*2;
	Tag* coverageTags = new Tag[coveragePositions];

	int cIndex = 0;
	for (int i=0;i<numTags[expIndex];i++) {
		Tag* t = &(exps[expIndex][i]);
		if (strandInfo == STRAND_SEPARATE && t->d == 1) {
			//fprintf(stderr, "skipping...\n");
			continue;
		}
		if (strandInfo == POSITIVE_STRAND && t->d == 1) {
			//fprintf(stderr, "skipping...\n");
			continue;
		}
		if (strandInfo == NEGATIVE_STRAND && t->d == 0) {
			//fprintf(stderr, "skipping...\n");
			continue;
		}
		coverageTags[cIndex].v=t->v;
		coverageTags[cIndex].p=t->p;
		coverageTags[cIndex].d=0;
		coverageTags[cIndex].len=t->len;
		cIndex++;
		coverageTags[cIndex].d = 0;
		coverageTags[cIndex].len = t->len;
		if (t->d == 0) {
			coverageTags[cIndex].p = t->p+fragLength;
			coverageTags[cIndex].v = -1*t->v;
		} else {
			coverageTags[cIndex].p = t->p-fragLength;
			coverageTags[cIndex-1].v = -1*t->v;
			coverageTags[cIndex].v = t->v;
		}
	//fprintf(stderr, "%d\n", coverageTags[cIndex].p+refPos);
		cIndex++;
	}
	coveragePositions = cIndex;
	if (coveragePositions < 1) {
		delete []coverageTags;
		return NULL;
	}

	qsort(coverageTags,coveragePositions,sizeof(Tag),&cmpTags);

	int last = 0;
	for (int i=1;i<coveragePositions;i++) {
		if (coverageTags[i].p == coverageTags[last].p) {
			coverageTags[last].v += coverageTags[i].v;
		} else {
			last++;
			if (last == i) continue;
			coverageTags[last].copy(&(coverageTags[i]));
		}
	}
	coveragePositions = last+1;

	float total = 0.0;
	for (int i=0;i<coveragePositions;i++) {
		total += coverageTags[i].v;
		coverageTags[i].v = total;
	}

	return coverageTags;
}


SNP::SNP() {
	init();
}
SNP::SNP(char* ref, char** alts, int numAlts, int numPeeps) {
	init();
	genotypes = new char*[numAlts+1];
	genotypes[0] = new char[strlen(ref)+1];
	strcpy(genotypes[0],ref);
	for (int i=1;i<=numAlts;i++) {
		genotypes[i] = new char[strlen(alts[i-1])+1];
		strcpy(genotypes[i],alts[i-1]);
	}
	numGenotypes = numAlts+1;
	allele1 = new char[numPeeps];
	allele2 = new char[numPeeps];
	numIndividuals = numPeeps;
	for (int i=0;i<numIndividuals;i++) {
		allele1[i] = 0;
		allele2[i] = 0;
	}
}
SNP::~SNP() {
	if (genotypes != NULL) {
		for (int i=0;i<numGenotypes;i++) {
			if (genotypes[i] != NULL) delete [](genotypes[i]);
		}
		delete []genotypes;
	}
	if (rvGenotypes != NULL) {
		for (int i=0;i<numGenotypes;i++) {
			if (rvGenotypes[i] != NULL) delete [](rvGenotypes[i]);
		}
		delete []rvGenotypes;
	}
	if (allele1 != NULL) delete []allele1;
	if (allele2 != NULL) delete []allele2;
	if (editDistances != NULL) delete []editDistances;
	allele1=NULL;
	allele2=NULL;
	editDistances=NULL;
}
void SNP::deleteBaseSNP(SNP* snp) {
	SNP* current = snp;
	while (current != NULL) {
		SNP* next = current->nextSNP;
		delete current;
		current = next;
	}
}
void SNP::init() {
	genotypes = NULL;
	rvGenotypes = NULL;
	numGenotypes = 0;
	allele1 = NULL;
	allele2 = NULL;
	numIndividuals=0;
	nextSNP = NULL;
	editDistances = NULL;
	count=0;
}
void SNP::addEditDistances(double* totals) {
	if (editDistances == NULL) calculateEditDistances();
	for (int i=0;i<numIndividuals;i++) {
		totals[i] += editDistances[i];
	}
}
void SNP::calculateEditDistances() {

	double* compare = new double[numGenotypes];
	compare[0] = 0.0;
	int L1 = strlen(genotypes[0]);
	for (int i=1;i<numGenotypes;i++) {
		int L2 = strlen(genotypes[i]);
		char* shorter = NULL;
		char* longer = NULL;
		int shorterL = L1;
		int longerL = L1;
		if (L1 < L2) {
			shorter = genotypes[0];
			shorterL = L1;
			longer = genotypes[i];
			longerL = L2;
		} else if (L1 > L2) {
			shorter = genotypes[i];
			shorterL = L2;
			longer = genotypes[0];
			longerL = L1;
		} else {
			//same length = edit distance set to length of SNP
			compare[i] = (double) L1;
			continue;
		}
		int maxMatch = 0;
		int match = 0;
		for (int j=0;j<shorterL;j++) {
			if (shorter[j] == longer[j]) {
				match++;
			} else {
				break;
			}
		}
		if (match > maxMatch) maxMatch = match;
		match = 0;
		for (int j=0;j<shorterL;j++) {
			if (shorter[shorterL-1-j] == longer[longerL-1-j]) {
				match++;
			} else {
				break;
			}
		}
		if (match > maxMatch) maxMatch = match;
		compare[i] = (double)maxMatch;
	}
	
	editDistances = new double[numIndividuals];
	for (int i=0;i<numIndividuals;i++) {
		editDistances[i] = (compare[(int)allele1[i]] + compare[(int)allele2[i]])/2;
	}
	delete []compare;
}
void SNP::printReference(FILE* fp,char strand) {
	if (strand == STRAND_NEGATIVE) {
		if (rvGenotypes == NULL) {
			rvGenotypes = new char*[numGenotypes];
			for (int i=0;i<numGenotypes;i++) {
				rvGenotypes[i]=new char[strlen(genotypes[i])+1];
				strcpy(rvGenotypes[i],genotypes[i]);
				revopp(rvGenotypes[i]);
			}
		}
		fprintf(fp, "%s", rvGenotypes[0]);
	} else {
		fprintf(fp, "%s", genotypes[0]);
	}
}
void SNP::print(FILE* fp,char strand) {
	char** geno = genotypes;
	if (strand == STRAND_NEGATIVE) {
		if (rvGenotypes == NULL) {
			rvGenotypes = new char*[numGenotypes];
			for (int i=0;i<numGenotypes;i++) {
				rvGenotypes[i]=new char[strlen(genotypes[i])+1];
				strcpy(rvGenotypes[i],genotypes[i]);
				revopp(rvGenotypes[i]);
			}
		}
		geno = rvGenotypes;
	}
	for (int i=0;i<numIndividuals;i++) {
		if (i>0) fprintf(fp,";");
		if (geno == NULL) {
			fprintf(fp,"./.");
		} else {
			if (geno[(int)allele1[i]] == NULL) {
				fprintf(fp,".");
			} else {
				fprintf(fp, "%s", geno[(int)allele1[i]]);
			}
			fprintf(fp, "/");
			if (geno[(int)allele2[i]] == NULL) {
				fprintf(fp,".");
			} else {
				fprintf(fp, "%s", geno[(int)allele2[i]]);
			}
		}
	}
}


PeakSNPs::PeakSNPs() {
	init();
}
PeakSNPs::~PeakSNPs() {
	if (snps != NULL) delete []snps;
	if (offsets != NULL) delete []offsets;
}
void PeakSNPs::init() {
	snps = NULL;
	numSNPs = 0;
	offsets = NULL;
}
void PeakSNPs::addSNP(SNP* snp, int pos) {
	SNP** oldsnps = snps;
	int* oldOffsets = offsets;
	snps = new SNP*[numSNPs+1];
	offsets = new int[numSNPs+1];
	for (int i=0;i<numSNPs;i++) {
		snps[i] = oldsnps[i];
		offsets[i] = oldOffsets[i];
	}
	if (oldsnps != NULL) delete []oldsnps;
	if (oldOffsets != NULL) delete []oldOffsets;
	snps[numSNPs] = snp;
	offsets[numSNPs] = pos;
	numSNPs++;
}
void PeakSNPs::printTotals(FILE* fp) {
	if (snps == NULL) return;

	int n = snps[0]->numIndividuals;
	double* totals = new double[n];
	for (int i=0;i<n;i++) totals[i] = 0.0;
	
	for (int i=0;i<numSNPs;i++) {
		snps[i]->addEditDistances(totals);
	}
	for (int i=0;i<n;i++) {
		if (i>0) fprintf(fp,",");
		fprintf(fp, "%.1lf", totals[i]);
	}
}
void PeakSNPs::print(FILE* fp,char strand) {
	for (int i=0;i<numSNPs;i++) {
		if (i>0) fprintf(fp,",");
		fprintf(fp,"%d(",offsets[i]);
		snps[i]->printReference(fp,strand);
		fprintf(fp,")=");
		snps[i]->print(fp,strand);
	}
}


void Peak::addSNP(SNP* snp, int pos) {
	if (snps == NULL) {
		snps = new PeakSNPs();
	}
	int diff = pos - refPos;
	if (strand == NEGATIVE_STRAND) diff *=-1;
	snps->addSNP(snp,diff);
}


//#########################################################################
//#########################################################################
//#########################################################################
//#########################################################################
//#########################################################################
//#########################################################################
//#########################################################################
//#########################################################################





TagLibrary::TagLibrary(char* dir) {
	directory = NULL;
	if (dir != NULL) {
		directory = new char[strlen(dir)+1];
		strcpy(directory,dir);
	}
	chrs= new Hashtable(100000);
	chrNames = NULL;
	uniqmapchrs = NULL;
	uniqmapDirectory = NULL;
	name = NULL;
	cmd = NULL;
	genome = NULL;
	memType = TAGLIBRARY_LOWMEM;
	totalTags = 0;
	parseAlignmentCpGMinValue = 0.0;
	parseAlignmentCcontext = CYTOSINE_CONTEXT_CG;
	mCflag = 0;
	sspeFlag = 0;
	flipFlag = 0;
	totalPositions = 0;
	medianTagCount= 0;
	averageTagsPerPosition=0.0;
	averageTagLength=0.0;
	tagAdjust = 0;
	maxmismatches = -1;
	minmapq = 10.0;
	gcNormalized = 0;
	gcAverage = -1.0;
	oligoNormalized = 0;
	tbp = 0.0;
	peStyleFlag = SEQFILE_FORMAT_UNKNOWN;
	maxtbp = 0.0;
	mintbp = 0.0;
	singleFile = 0;
	singleFilename = NULL;
	pairedEndFlag = 0;
	forceSingleReadFlag = 0;
	restrictionSite = NULL;
	tagFile = NULL;
	fragmentLengthSetting = FRAGMENT_LEN_AUTO;
	fragmentLengthEstimate = NULL_INT_VALUE;
	peakSizeEstimate = NULL_INT_VALUE;
	gsizeEstimate = 0;
	revStrand = 0;
	numCPUs = 1;
}
TagLibrary::~TagLibrary() {
	if (chrs != NULL) {
		char** keys = chrs->keys();
		qsort(keys,chrs->total,sizeof(char*),&chrcmp);
		for (int i=0;i<chrs->total;i++) {
			ChrTags* ct = (ChrTags*) chrs->search(keys[i]);
			delete ct;
			delete [](keys[i]);
		}
		delete []keys;
		delete chrs;
	}
	if (uniqmapchrs != NULL) {
		char** keys = uniqmapchrs->keys();
		qsort(keys,chrs->total,sizeof(char*),&chrcmp);
		for (int i=0;i<uniqmapchrs->total;i++) {
			UniqMapChrs* ct = (UniqMapChrs*) uniqmapchrs->search(keys[i]);
			delete ct;
			delete [](keys[i]);
		}
		delete []keys;
		delete uniqmapchrs;
	}
	if (chrNames != NULL) delete chrNames;
	if (name != NULL) delete []name;
	if (genome != NULL) delete []genome;
	if (cmd != NULL) delete []cmd;
	if (uniqmapDirectory != NULL) delete []uniqmapDirectory;

}
void TagLibrary::setName(char* newname) {
	if (newname == NULL) return;
	if (name != NULL) delete []name;
	name = new char[strlen(newname)+1];
	strcpy(name, newname);
}
void TagLibrary::setGenome(char* newgenome) {
	if (newgenome == NULL) return;
	if (genome != NULL) delete []genome;
	genome = new char[strlen(newgenome)+1];
	strcpy(genome, newgenome);
}
void TagLibrary::setSingleFile(int singleFileFlag) {
	singleFile= singleFileFlag;
	singleFilename = new char[10000];
	sprintf(singleFilename,"%s/genome.tags.tsv",directory);
}
void TagLibrary::setFragLength(int fraglength) {
	fragmentLengthSetting = fraglength;
	fragmentLengthEstimate = fraglength;
}
void TagLibrary::setCMD(char* newcmd) {
	if (newcmd == NULL) return;
	if (cmd != NULL) delete []cmd;
	cmd = new char[strlen(newcmd)+1];
	strcpy(cmd,newcmd);
}
void TagLibrary::addUniqMap(char* dir) {
	
	uniqmapDirectory = new char[strlen(dir)+1];
	strcpy(uniqmapDirectory,dir);
//WORK ON THIS NEXT

}

void TagLibrary::makeDirectory() {
	fprintf(stderr, "\tCreating directory: %s and removing existing *.tags.tsv\n",directory);
	char* filename = new char[10000];
	strcpy(filename, "mkdir -p \"");
	strcat(filename, directory);
	strcat(filename, "\"");
	(void)system(filename);
	strcpy(filename, "rm -f \"");
	strcat(filename, directory);
	strcat(filename, "\"/*.tags.tsv");
	(void)system(filename);
	delete []filename;
}	

void TagLibrary::setSingleRead(int flag) {
	forceSingleReadFlag = flag;
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*) chrs->search(keys[i]);
		ct->forceSingleReadFlag = flag;
		delete [](keys[i]);
	}
	delete []keys;
}

char* unzipFileIfNeeded(char* file, int& zipFlag, int &curFormat) {
	zipFlag = 0;
	if (file == NULL) {
		return NULL;
	}
	int len = strlen(file);
	char* newName = new char[len+1];
	char* command = new char[100+2*len+1];
	strcpy(newName,file);
	if (len > 4) {
		if (file[len-3] == '.' && file[len-2] == 'g' && file[len-1] == 'z') {
			//gzipped file
			fprintf(stderr, "\tTreating %s as a GNU zip file\n", file); 
			newName[len-3] = '\0';
			zipFlag = ZIPPED_FLAG_GZ;
			sprintf(command, "gunzip -c \"%s\" > \"%s\"",file, newName);
			(void)system(command);
		} else if (file[len-4] == '.' && file[len-3] == 'b' && file[len-2] == 'z' && file[len-1] == '2') {
			//bz2 zipped file
			fprintf(stderr, "\tTreating %s as a bz2 zip file\n", file); 
			newName[len-4] = '\0';
			zipFlag = ZIPPED_FLAG_BZ2;
			sprintf(command, "bunzip2 -c \"%s\" > \"%s\"",file,newName);
			(void)system(command);
		} else if (file[len-4] == '.' && file[len-3] == 'z' && file[len-2] == 'i' && file[len-1] == 'p') {
			//zip file
			fprintf(stderr, "\tTreating %s as a zip file\n", file); 
			newName[len-4] = '\0';
			zipFlag = ZIPPED_FLAG_ZIP;
			sprintf(command, "unzip -p \"%s\" > \"%s\"",file,newName);
			(void)system(command);
		} else if (file[len-4] == '.' && file[len-3] == 'b' && file[len-2] == 'a' && file[len-1] == 'm') {
			//zip file
			fprintf(stderr, "\tTreating %s as a bam file\n", file); 
			newName[len-4] = '\0';
			zipFlag = ZIPPED_FLAG_BAM;
			curFormat = FORMAT_SAM;
			sprintf(command, "samtools view -h \"%s\" > \"%s\"",file,newName);
			(void)system(command);
		} else {
			//not recognized, or not zipped...
		}

	}
	delete []command;
	return newName;
}
void rezipFileIfNeeded(char* file, int zipFlag) {
	int len = strlen(file);
	char* command = new char[100+2*len+1];
	if (zipFlag) {
		sprintf(command, "rm \"%s\"", file);
		//fprintf(stderr, "Ending: %s\n", command);
		(void)system(command);
	}
	delete []file;
	delete []command;
}
		


void TagLibrary::parseAlignmentFiles(char** files, int numFiles, int format, int mode,
							char** tagDirs, int numDirs,char** tagFiles,int numTagFiles) {
		
	makeDirectory();


	for (int i=0;i<numFiles;i++) {
		fprintf(stderr, "\n");
		if (pairedEndFlag) {
			fprintf(stderr, "\tReading paired end alignment files %s\n", files[i]);
			readPEAlignment(files[i],format,mode);
		} else {
			int zippedFlag = 0;
			int currentFormat = format;
			char* workingFilename = unzipFileIfNeeded(files[i],zippedFlag,currentFormat);
			fprintf(stderr, "\tReading alignment file %s\n", files[i]);
			readAlignment(workingFilename,currentFormat,mode,TAGLIBRARY_SINGLE_FLAG);
			rezipFileIfNeeded(workingFilename,zippedFlag);
		}
	}
	for (int i=0;i<numDirs;i++) {
		fprintf(stderr, "\tAdding tag directory %s\n", tagDirs[i]);
		addTagDirectory(tagDirs[i]);
	}
	for (int i=0;i<numTagFiles;i++) {
		fprintf(stderr, "\tAdding tag file %s\n", tagFiles[i]);
		readNewTagFile(tagFiles[i]);
	}
	if (totalTags < 1) {
		fprintf(stderr, "\t!!! Something is wrong - no reads were added to tag directory !!!\n");
		fprintf(stderr, "\t!!! Check your input files or the makeTagDirectory command options... !!!\n");
		exit(0);
	}
	//fprintf(stderr, "Troubleshooting...\n");	
	//int ok = system("sleep 4");
	//ok = 1;
	optimizeTagFiles();
	fprintf(stderr, "\tTotal Tags = %.1f\n", totalTags);
	fprintf(stderr, "\tTotal Positions = %lld\n", totalPositions);
	printTagInfo();
}

void TagLibrary::addTagDirectory(char* tagDir) {
	
	DIR *dir = opendir(tagDir);
	char * file = NULL;
	struct dirent* entry = NULL;
	char* fullfile = new char[10000];
	while ((entry = readdir(dir))) {
		file = entry->d_name;		
		int len = strlen(file);
		if (len < 9) continue;	
		if (strcmp(".tags.tsv",&(file[len-9])) == 0) {
			strcpy(fullfile,tagDir);
			strcat(fullfile,"/");
			strcat(fullfile,file);

			readNewTagFile(fullfile);

		}
	}
	closedir(dir);
	delete []fullfile;
	

}

void TagLibrary::readTagDirectory() {

	readTagInfo();


	if (singleFile) {
		char* file = new char[10000];
		strcpy(file, directory);
		strcat(file, "/genome.tags.tsv");
		tagFile = fopen(file,"r");
		if (tagFile == NULL) {
			fprintf(stderr, "!!! For some reason, could not open file: %s\n", file);
			exit(0);
		}
		readSingleTagFile();
		//////////////////////////////////


	} else {
	
		DIR *dir = opendir(directory);
		char * file = NULL;
		struct dirent* entry = NULL;
		char* fullfile = new char[10000];
		char* chr = new char[10000];
		while ((entry = readdir(dir))) {
			file = entry->d_name;		
			int len = strlen(file);
			
			if (len < 9) continue;	
			if (strcmp(".tags.tsv",&(file[len-9])) == 0) {
				strcpy(chr,file);
				chr[len-9]='\0';
				strcpy(fullfile,directory);
				strcat(fullfile,"/");
				strcat(fullfile,file);
				//fprintf(stderr, "\tReading %s\n", fullfile);
	
				ChrTags* chrtags = (ChrTags*) chrs->search(chr);
				if (chrtags == NULL) {
					chrtags = new ChrTags(chr);
					chrs->insert(chrtags,chr);
					chrtags->mCflag = mCflag;
				}
	
				if (singleFile) {
					chrtags->singleFile = 1;
				} else {
					char* file2 = getTagFileName(chr);
					chrtags->setTagFile(file2);
					delete []file2;
				}
				if (pairedEndFlag) {
					chrtags->pairedEndFlag = pairedEndFlag;
				}
			}
		}
		closedir(dir);
		delete []fullfile;
	}
}


void TagLibrary::readPEAlignment(char* file, int format, int mode) {
	char** files = new char*[100];
	int numFiles = 0;
	split(file, files, numFiles, ',');
	if (numFiles < 2) {
		fprintf(stderr, "!!! Missing comma between file names for Paired End data !!!\n");
		fprintf(stderr, "!!! %s !!!\n",file);
		exit(0);
	}

	int zippedFlag = 0;
	int currentFormat = format;
	char* workingFilename = unzipFileIfNeeded(files[0],zippedFlag,currentFormat);
	//fprintf(stderr, "\tReading alignment file %s\n", files[i]);
	readAlignment(workingFilename,currentFormat,mode,TAGLIBRARY_PE_READ1_FLAG);
	rezipFileIfNeeded(workingFilename,zippedFlag);
	zippedFlag = 0;
	currentFormat = format;
	workingFilename = unzipFileIfNeeded(files[1],zippedFlag,currentFormat);
	//fprintf(stderr, "\tReading alignment file %s\n", files[i]);
	readAlignment(workingFilename,currentFormat,mode,TAGLIBRARY_PE_READ2_FLAG);
	rezipFileIfNeeded(workingFilename,zippedFlag);

	delete []files;

	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	if (chrNames == NULL) {
		chrNames = new Hashtable(1000);
	}
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*) chrs->search(keys[i]);
		if (ct->tagfpR1 != NULL) fclose(ct->tagfpR1);
		if (ct->tagfpR2 != NULL) fclose(ct->tagfpR2);
		ct->tagfpR1 = NULL;
		ct->tagfpR2 = NULL;
		if (NULL == chrNames->search(keys[i])) {
			chrNames->insert(keys[i],keys[i]);
		}
	}

	//FILE* badfp = fopen("nomatch.tsv","w");

	//joinPEreads
	char* fname = new char[10000];
	char* buf = new char[BUFFER];
	char** line = new char*[BUFFER];
	int numCols = 0;
	fprintf(stderr, "\t\tMatching paired reads...\n");
	for (int i=0;i<chrs->total;i++) {
		fprintf(stderr, "\t\t\t%s\n", keys[i]);
		ChrTags* ct = (ChrTags*) chrs->search(keys[i]);
		sprintf(fname,"%s.R1",ct->tagFile);
		FILE* fp = fopen(fname,"r");

		Hashtable* read1tags = new Hashtable(10000000);	
		//Inttable* found = new Inttable(1000000);	
		while (fgets(buf, BUFFER, fp) != NULL) {
			split(buf,line,numCols,'\t');
			if (peStyleFlag == SEQFILE_FORMAT_ILLUMINAPE) {
				trimIlluminaReadName(line[0]);
			}
			int p = 0;
			int len = 0;
			int d = 0;
			float v = 0.0;
			char* c = (char*) chrNames->search(line[1]);
			//fprintf(stderr, "|%s|\t|%s|\n", c,line[1]);
			sscanf(line[2],"%d",&p);
			sscanf(line[3],"%d",&d);
			sscanf(line[4],"%f",&v);
			sscanf(line[5],"%d",&len);
			PETag* petag = new PETag(line[0],c,p,(char)d,v,len);
			read1tags->insert(petag,line[0]);
			//found->insert(0,line[0]);
		}
		fclose(fp);

		for (int j=0;j<chrs->total;j++) {
			ChrTags* ct2 = (ChrTags*) chrs->search(keys[j]);
			sprintf(fname,"%s.R2",ct2->tagFile);
			fp = fopen(fname,"r");

			while (fgets(buf, BUFFER, fp) != NULL) {
				split(buf,line,numCols,'\t');
				if (peStyleFlag == SEQFILE_FORMAT_ILLUMINAPE) {
					trimIlluminaReadName(line[0]);
				}
				PETag* petag = (PETag*)read1tags->search(line[0]);
				if (petag != NULL) {
					int p = 0;
					int d = 0;
					float v = 0.0;
					int len = 0;
					char* c = (char*) chrNames->search(line[1]);
					sscanf(line[2],"%d",&p);
					sscanf(line[3],"%d",&d);
					sscanf(line[4],"%f",&v);
					sscanf(line[5],"%d",&len);
					petag->p2 = p;
					petag->d2 = (char)d;
					petag->len2 = len;
					petag->chr2 = c;
					ct->printAlignedPETag(petag,0);
					ct2->printAlignedPETag(petag,1);
					//int fnum = found->search(line[0]);
					//found->insert(fnum+1,line[0]);
				} else {
				}
			}
			fclose(fp);
		}
		char** hkeys = read1tags->keys();
		for (int j=0;j<read1tags->total;j++) {
			PETag* petag = (PETag*) read1tags->search(hkeys[j]);
			//int fnum = found->search(hkeys[j]);
			//if (fnum == 0) petag->print(badfp);
			//if (fnum > 1) {
			//	fprintf(stderr, "Multiple matches: %d %s ",fnum , hkeys[j]);
			//	petag->print(stderr);
			//}
			delete petag;
			delete [](hkeys[j]);
		}
		delete []hkeys;
		delete read1tags;
		//delete found;
	}

	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*) chrs->search(keys[i]);
		sprintf(fname,"rm -f \"%s.R1\"",ct->tagFile);
		(void)system(fname);
		sprintf(fname,"rm -f \"%s.R2\"",ct->tagFile);
		(void)system(fname);
	}	
	//joinPEreads
	delete []fname;
	delete []buf;
	delete []line;
	delete []keys;


}

void TagLibrary::trimIlluminaReadName(char* name) {
	int slen = strlen(name);
	if (slen > 1) {
		name[slen-1] = '\0';
	}
}


void TagLibrary::readAlignment(char* file, int format, int mode, int PEflag) {
	
	FILE* fp = fopen(file,"r");
	if (fp == NULL) {
		fprintf(stderr, "Could not open %s\n", file);
		return;
	}

	char* buf = new char[BUFFER];
	char* name = new char[BUFFER];
	char* lastname = new char[BUFFER];
	char* chr = new char[BUFFER];	
	char* lastchr = new char[BUFFER];	
	char** line = new char*[BUFFER];
	char** line2 = new char*[BUFFER];
	int numCols = 0;
	int numCols2 = 0;
	char cigarCodes[1000];
	int cigarLens[1000];
	int numCodes;
	name[0]='\0';
	lastname[0]='\0';
	chr[0] = '\0';
	lastchr[0] = '\0';


	PETag* petag = new PETag();
	int pos = 0;
	int secondFlag = 0;
	int mismatches = 0;
	int mismatchQuality = 0;
	int dir = 0;
	int lastpos = 0;
	int lastmismatches = 0;
	int lastmismatchQuality = 0;
	int lastdir = 0;
	float numTags = 0;
	int len = -1;

	double avgReadsPerPos = 0;
	int totalReads = 0;

	while (fgets(buf, BUFFER, fp) != NULL) {
		//sscanf(buf, "%s\t%s\t%d\t%d\t%d",name,chr,&pos,&dir,&tagCount);
		split(buf, line, numCols,'\t');

		if (numCols < 3) continue;


		if (format == FORMAT_UNKNOWN) {
			if (numCols > 6 && format == FORMAT_UNKNOWN) {
				if (line[1][0] == '-' || line[1][0] == '+') {
					fprintf(stderr, "\tGuessing that your alignment file is bowtie format\n");
					format = FORMAT_BOWTIE;
				}
			}
			if (numCols > 2 && format == FORMAT_UNKNOWN) {
				if (line[0][0] == 'c' && line[0][1] == 'h' && line[0][2] == 'r') {
					fprintf(stderr, "\tGuessing that your alignment file is BED format\n");
					format = FORMAT_BED;
				}
			}
			if (numCols > 20 && format == FORMAT_UNKNOWN) {
				fprintf(stderr, "\tGuessing that your alignment file is eland_export format\n");
				format = FORMAT_ELANDEXPORT;
			}
			if (numCols == 4 && format == FORMAT_UNKNOWN) {
				if (strcmp(line[2],"NM")==0 || strcmp(line[2],"QC")==0) {
					format = FORMAT_ELANDEXTENDED;
				} else {
					int colonCount = 0;
					for (unsigned int i=0;i<strlen(line[2]);i++) {
						if (line[2][i] == ':') colonCount++;
					}
					if (colonCount == 2) format = FORMAT_ELANDEXTENDED;
				}
				fprintf(stderr, "\tGuessing that your alignment file is eland_extended format\n");
			}
			if (line[0][0] == '@' && format == FORMAT_UNKNOWN) {
				if (strcmp(line[0],"@SQ")==0 || strcmp(line[0],"@HD")==0 || strcmp(line[0],"@RG")==0) {
					fprintf(stderr, "\tGuessing that your alignment file is SAM format\n");
					format = FORMAT_SAM;
				}
			}
			if (format == FORMAT_UNKNOWN) {
				format = FORMAT_ELANDRESULT;
				fprintf(stderr, "\tGuessing (by default) that your alignment file is eland_result format\n");
			}
		}

		if (format == FORMAT_HICSUMMARY) {

			if (numCols < 7) continue;

			pairedEndFlag = 1;

			petag->name = line[0];
			petag->chr1 = line[1];
			sscanf(line[2],"%d",&pos);
			petag->p1 = pos;
			char dir = STRAND_POSITIVE;
			if (line[3][0] == '-') dir = STRAND_NEGATIVE;
			if (flipFlag) {
				if (dir == STRAND_POSITIVE) dir = STRAND_NEGATIVE;
				else dir = STRAND_POSITIVE;
			}
			petag->d1 = dir;
			petag->len1 = 1;
			petag->v = 1.0;
			

			petag->chr2 = line[4];
			int pos2 = 0;
			sscanf(line[5],"%d",&pos2);
			petag->p2 = pos2;
			char dir2 = STRAND_POSITIVE;
			if (line[6][0] == '-') dir2 = STRAND_NEGATIVE;
			if (flipFlag) {
				if (dir2 == STRAND_POSITIVE) dir2 = STRAND_NEGATIVE;
				else dir2 = STRAND_POSITIVE;
			}
			petag->d2 = dir2;
			petag->len2 = 1;
			addAlignedPETag(petag);

			petag->chr1=line[4];
			petag->chr2=line[1];
			petag->p2 = petag->p1;
			petag->p1 = pos2;
			petag->d1 = dir2;
			petag->d2 = dir;
			addAlignedPETag(petag);

		} else if (format == FORMAT_SAM) {
			if (line[0][0] == '@') continue;
			if (line[2][0] == '*') continue;
			strcpy(name,line[0]);
			strcpy(chr,line[2]);
			int samFlag = 0;
			int initLen = 0;
			len = strlen(line[9]);
			if (len < minReadLength || len > maxReadLength) continue;
			sscanf(line[3],"%d",&pos);
			if (pos == 0) continue;
			if (line[5][0] == '*') continue;

			sscanf(line[1],"%d",&samFlag);
			dir = STRAND_POSITIVE;
			int peFlag=0;
			if (samFlag & 0x10) {
				dir = STRAND_NEGATIVE;
			}
			if (flipFlag) {
				if (dir==STRAND_POSITIVE) dir =STRAND_NEGATIVE;
				else dir=STRAND_POSITIVE;
			}


			if (sspeFlag && (samFlag & 0x80)) {
				if (dir == STRAND_POSITIVE) {
					dir = STRAND_NEGATIVE;
				} else {
					dir = STRAND_POSITIVE;
				}
			}
			if (samFlag & 0x1) {
				peFlag = 1;
			}
			if (maxmismatches > -1) {
				int numMis = 0;
				for (int i=11;i<numCols;i++) {
					if (strlen(line[i]) < 6) continue;
					if (strncmp(line[i],"MD:Z:",5)==0) {
						char* misStr = &(line[i][5]);
						numMis = parseMDZstr(misStr);
						break;
					}
				}
				if (numMis > maxmismatches) continue;
			}

			if (mode == MODE_UNIQUE) {
				if (samFlag & 0x100) continue;
				if (minmapq < 0) {
					int good = 1;
					double alnScore = -1e8;
					double secScore = -1e9;
					for (int i=11;i<numCols;i++) {
						if (strlen(line[i]) < 6) continue;
						if (strncmp(line[i],"X0:i:",5)==0) {
							if (strcmp(&(line[i][5]),"1")!=0) {
								good = 0;
							}
						} else if (strncmp(line[i],"AS:i:",5)==0) {
							sscanf(&(line[i][5]),"%lf",&alnScore);
						} else if (strncmp(line[i],"XS:i:",5)==0) {
							sscanf(&(line[i][5]),"%lf",&secScore);
						}
					}
					if (alnScore > secScore) {
					} else {
						good = 0;
					}
					if (good == 0) continue;
				} else {
					double mapqScore = 0.0;
					sscanf(line[4], "%lf",&mapqScore);
					if (mapqScore < minmapq) continue;
				}
			} else if (mode == MODE_KEEPONE) {
				if (samFlag & 0x100) continue;
			}
			int startPos = getRightCoordFromCIGAR(line[5],dir,cigarCodes,cigarLens,numCodes,initLen);
			if (startPos == CIGAR_ERROR) {
				fprintf(stderr, "\tError in read: %s\n", name);
				startPos = 0;
			}
			pos += startPos;
			len = initLen;
			float v=1.0;
			if (peFlag) v=0.5;
			addAlignedTag(name,chr,pos,(char)dir,len,v,PEflag);
		}

		if (format == FORMAT_BOWTIE || format == FORMAT_BOWTIE_COLOR) {
			if (numCols < 6) {
				//fprintf(stderr, "problem\n");
				continue;
			}
			
			strcpy(name,line[0]);
			strcpy(chr,line[2]);
			dir = STRAND_POSITIVE;
			if (line[1][0] == '-') dir = STRAND_NEGATIVE;
			if (flipFlag) {
				if (dir == STRAND_POSITIVE) dir = STRAND_NEGATIVE;
				else dir = STRAND_POSITIVE;
			}
			sscanf(line[3],"%d",&pos);
			pos++;
			int taglen = strlen(line[4]);
			if (dir == 1) {
				pos += taglen-1;
			}
			len = taglen;
			if (len < minReadLength || len > maxReadLength) continue;
			mismatches = 0;

			int quality = getQualityTotal(line[5]);
			mismatchQuality = -1*quality;

			if (format == FORMAT_BOWTIE_COLOR) {
				//split(line[5],line2,numCols2,'!');
				//mismatches = numCols2;
			} else {
				if (numCols > 7) {
					split(line[7],line2,numCols2,',');
					if (line2[0][0] != '\0') {
						mismatches = numCols2;
					}
				}	
			}
			//fprintf(stderr, "%s\t%d\t%d\n", line[0],mismatchQuality, mismatches);
			int good = 0;

			if (mode == MODE_UNIQUE) {
				if (strcmp(name,lastname)!=0) {
					if (secondFlag == 0) {
						good = 1;
					}
					secondFlag = 0;
				} else {
					secondFlag = 1;
					if (lastmismatchQuality < mismatchQuality 
									|| lastmismatches < mismatches) {
						good = 1;
					}
				}
			} else if (mode == MODE_KEEPONE) {
				if (secondFlag == 0) good = 1;
				if (strcmp(name,lastname)!=0) {
					secondFlag = 0;
				} else {
					secondFlag = 1;
				}
			}

			if (good) {
				if (lastchr[0] != '\0') {
					addAlignedTag(lastname,lastchr,lastpos,(char)lastdir,len,1.0,PEflag);
				}
			}
			good = 0;
			char* tmp = chr;
			chr = lastchr;
			lastchr = tmp;
			tmp = name;
			name = lastname;
			lastname = tmp;
			lastmismatches = mismatches;
			lastdir = dir;
			lastpos = pos;
		} else if (format == FORMAT_BED) {
			strcpy(chr,line[0]);
			name[0] = '\0';

			dir = STRAND_POSITIVE;
			numTags = 1.0;
			if (numCols > 3) {
				if (line[3][0] == '+' || line[3][0] =='-') {
					if (line[3][0] == '-') dir = STRAND_NEGATIVE;
				} else if (numCols > 5) {
					strcpy(name,line[3]);
					if (line[5][0] =='-') {
						dir = STRAND_NEGATIVE;
					}
				}
			}
			if (flipFlag) {
				if (dir == STRAND_POSITIVE) dir = STRAND_NEGATIVE;
				else dir = STRAND_POSITIVE;
			}
			if (numCols > 4 && mode == MODE_BED_FORCE5TH) {
				sscanf(line[4],"%f",&numTags);
			}
			int start = 0;
			int end = 0;
			sscanf(line[1], "%d", &start);
			sscanf(line[2], "%d", &end);
			if (dir == 0) {
				pos = start+1; //0 base index
			} else {
				pos = end;
			}
			len = end-start;
			if (len < minReadLength || len > maxReadLength) continue;
			addAlignedTag(name,chr,pos,(char)dir,len,numTags,PEflag);
			avgReadsPerPos += numTags;
			totalReads++;
		} else if (format == FORMAT_MCPGBED) {
			strcpy(chr,line[0]);
			name[0] = '\0';

			dir = 0;
			numTags = 1.0;
			if (numCols < 11) {
				continue;
			}
			double mC = 0.0;
			double unmC = 0.0;
			sscanf(line[9],"%lf",&unmC);
			sscanf(line[10],"%lf",&mC);

			double totalC = mC+unmC;
			if (totalC < parseAlignmentCpGMinValue) continue;
			double ratio = 0.0;
			if (totalC > 0.0) {
				ratio = mC/totalC;
			}

			int start = 0;
			int end = 0;
			sscanf(line[1], "%d", &start);
			sscanf(line[2], "%d", &end);
			if (line[5][0] =='-') {
				dir = 1;
			}
			if (dir == 0) {
				pos = start+1; //0 base index
			} else {
				pos = end;
			}
			len = end-start;
			if (len < minReadLength || len > maxReadLength) continue;
			addAlignedTag(name,chr,pos,(char)dir,(int)totalC,ratio,PEflag);
			avgReadsPerPos += ratio;
			totalReads++;
		} else if (format == FORMAT_LISTER_ALLC) {
			strcpy(chr,line[0]);
			name[0] = '\0';
			dir = 0;
			numTags = 1.0;
			if (numCols < 7) {
		        continue;
			}
			
			if (parseAlignmentCcontext == CYTOSINE_CONTEXT_CG) {
				if (strcmp(line[3],"CG")!=0) continue;
			} else if (parseAlignmentCcontext == CYTOSINE_CONTEXT_CHH) {
				if (strcmp(line[3],"CHH")!=0) continue;
			} else if (parseAlignmentCcontext == CYTOSINE_CONTEXT_CHG) {
				if (strcmp(line[3],"CHG")!=0) continue;
			} else if (parseAlignmentCcontext == CYTOSINE_CONTEXT_ALL) {
			}

			if (chr[0] != 'c' && chr[0] != 'C') {
				addChrToName(chr);
			}
			double mC = 0.0;
			double unmC = 0.0;
			sscanf(line[6],"%lf",&unmC);
			sscanf(line[4],"%lf",&mC);

			double totalC = mC+unmC;
			if (totalC < parseAlignmentCpGMinValue) continue;
			double ratio = 0.0;
			if (totalC > 0.0) {
				ratio = mC/totalC;
			}
			sscanf(line[1], "%d", &pos);
			pos++; //files are 0-based index
			if (line[2][0] == '-') {
				dir = 1;
			}
			len = 1;
			addAlignedTag(name,chr,pos,(char)dir,(int)totalC,ratio,PEflag);
			avgReadsPerPos += ratio;
		} else if (format == FORMAT_BISMARK) {
			strcpy(chr,line[0]);
			name[0] = '\0';
			dir = 0;
			numTags = 1.0;
			if (numCols < 7) {
		        continue;
			}
			
			if (parseAlignmentCcontext == CYTOSINE_CONTEXT_CG) {
				if (strcmp(line[5],"CG")!=0) continue;
			} else if (parseAlignmentCcontext == CYTOSINE_CONTEXT_CHH) {
				if (strcmp(line[5],"CHH")!=0) continue;
			} else if (parseAlignmentCcontext == CYTOSINE_CONTEXT_CHG) {
				if (strcmp(line[5],"CHG")!=0) continue;
			} else if (parseAlignmentCcontext == CYTOSINE_CONTEXT_ALL) {
			}

			if (chr[0] != 'c' && chr[0] != 'C') {
				addChrToName(chr);
			}
			double mC = 0.0;
			double unmC = 0.0;
			sscanf(line[4],"%lf",&unmC);
			sscanf(line[3],"%lf",&mC);

			double totalC = mC+unmC;
			if (totalC < parseAlignmentCpGMinValue) continue;
			double ratio = 0.0;
			if (totalC > 0.0) {
				ratio = mC/totalC;
			}
			sscanf(line[1], "%d", &pos);
			if (line[2][0] == '-') {
				dir = 1;
			}
			len = 1;
			addAlignedTag(name,chr,pos,(char)dir,(int)totalC,ratio,PEflag);
			avgReadsPerPos += ratio;
		} else if (format == FORMAT_ELANDRESULT) {
			if (numCols < 9) continue;
			if (strcmp(line[6],"") == 0) continue;
			if (strcmp(line[7],"") == 0) continue;
			strcpy(name,line[0]);
			strcpy(chr,line[6]);
			cleanUpChrName(chr);
			sscanf(line[7],"%d",&pos);
			dir = STRAND_POSITIVE;
			len = strlen(line[1]);
			if (len < minReadLength || len > maxReadLength) continue;
			if (line[8][0] == 'R') {
				dir = STRAND_NEGATIVE;
				pos += len-1;
			}
			if (flipFlag) {
				if (dir == STRAND_POSITIVE) dir = STRAND_NEGATIVE;
				else dir = STRAND_POSITIVE;
			}
			addAlignedTag(name,chr,pos,(char)dir,len,1.0,PEflag);
		} else if (format == FORMAT_ELANDEXTENDED) {
			if (numCols < 4) continue;
			split(line[2],line2,numCols2,':');
			if (numCols2 != 3) continue;
			strcpy(name,line[0]);
			int mis = -1;
			for (int i=0;i<numCols;i++) {
				int numMatches;
				sscanf(line2[i],"%d",&numMatches);
				if (numMatches == 0) continue;
				if (numMatches > 1 && mode != MODE_KEEPONE) break;
				mis = i;
				break;
			}
			if (mis < 0) continue;
			split(line[3],line2,numCols2,',');
			char* curChr=NULL;
			for (int i=0;i<numCols2;i++) {
				char* posStr = line2[i];
				char curDir = -1;
				int curMis = 0;
				int posFlag = 0;
				int slen = strlen(line2[i]);
				for (int j=0;j<slen;j++) {
					if (line2[i][j] == ':') {
						line2[i][j] = '\0';
						curChr = line2[i];
						posStr = &(line2[i][j+1]);
					} else if (line2[i][j] == 'F') {
						curDir = 0;
						posFlag = 1;
						line2[i][j] = '\0';
					} else if (line2[i][j] == 'R') {
						curDir = 1;
						posFlag = 1;
						line2[i][j] = '\0';
					} else if (posFlag) {
						switch (line2[i][j]) {
							case 'A':
							case 'C':
							case 'G':
							case 'T':
							case 'N':
							case 'a':
							case 'c':
							case 'g':
							case 't':
							case 'n':
								curMis++;
							default:
								break;
						}
					}
				}
				if (curMis == mis) {
					strcpy(chr,curChr);
					cleanUpChrName(chr);
					sscanf(posStr, "%d", &pos);
					len = strlen(line[1]);
					if (curDir == 1) {
						pos += len-1;
					}
					if (flipFlag) {
						if (curDir == STRAND_POSITIVE) curDir = STRAND_NEGATIVE;
						else curDir = STRAND_POSITIVE;
					}
					if (len < minReadLength || len > maxReadLength) {
					} else {
						addAlignedTag(name,chr,pos,curDir,len,1.0,PEflag);
					}
					break;
				}
			}
		} else if (format == FORMAT_ELANDEXPORT) {
			if (numCols < 20) continue;
			strcpy(chr,line[10]);
			if (strcmp(chr,"") == 0) continue;
			if (strcmp(line[12],"") == 0) continue;
			strcpy(name,line[0]);
			cleanUpChrName(chr);
			sscanf(line[12],"%d",&pos);
			dir = STRAND_POSITIVE;
			len = strlen(line[8]);
			if (line[13][0] == 'R') {
				dir = STRAND_NEGATIVE;
				pos += len-1;
			}
			if (flipFlag) {
				if (dir == STRAND_POSITIVE) dir = STRAND_NEGATIVE;
				else dir = STRAND_POSITIVE;
			}
			if (len < minReadLength || len > maxReadLength) continue;
			addAlignedTag(name,chr,pos,(char)dir,len,1.0,PEflag);
		}
	}
	fclose(fp);

	if ((format == FORMAT_BOWTIE || format == FORMAT_BOWTIE_COLOR) && secondFlag == 0) {
		if (lastchr[0] != '\0') {
			if (len < minReadLength || len > maxReadLength) {
			} else {
				addAlignedTag(lastname,lastchr,lastpos,(char)lastdir,len,1.0,PEflag);
			}
		}
	}

	if (format == FORMAT_BED) {
		if (totalReads > 1) {
			avgReadsPerPos /= (double)totalReads;
			if (avgReadsPerPos > 10.0) {
				fprintf(stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
				fprintf(stderr, "Average reads per BED file: %.1lf\n", avgReadsPerPos);
				fprintf(stderr, "Good chance that the 5th column of your BED file has a weird value in it!\n");
				fprintf(stderr, "By default, this is read as the number of reads for that position\n");
				fprintf(stderr, "To count each entry as only one read (ignore the 5th column) use -forceBED\n");
				fprintf(stderr, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
			}
		}
	}
			

	delete petag;
	delete []buf;
	delete []chr;
	delete []lastchr;
	delete []name;
	delete []lastname;
	delete []line;
	delete []line2;
}
void TagLibrary::addChrToName(char* str) {
	int len = strlen(str);
	str[len+3]='\0';
	for (int i=len-1;i>=0;i--) {
		str[i+3] = str[i];
	}
	str[0] = 'c';
	str[1] = 'h';
	str[2] = 'r';
}
int TagLibrary::parseMDZstr(char* misStr) {
	int numMis = 0;
	int i=0;
	while (misStr[i] != '\0') {
		if ((misStr[i] > 64 && misStr[i] < 91) || (misStr[i] > 96 && misStr[i] < 123)) {
			numMis++;
		}
		i++;
	}
	return numMis;
}
int TagLibrary::getRightCoordFromCIGAR(char* str, int dir, char* cigarCodes, int* cigarLens, int &numCodes, int &initLen) {
	numCodes = 0;
	int i = 0;
	char* start = str;
	int lastStartI = 0;
	int totalLen= 0;

	int firstStart = -1;
	int firstEnd = -1;
	int lastStart = -1;
	int lastEnd = -1;
	int firstActive = 0;
	int lastActive = 0;
	//fprintf(stderr, "%s\n",str);
	while (str[i] != '\0') {
		if (str[i] < 48 || str[i] > 57) {
			if (i-lastStartI == 0) {
				fprintf(stderr, "Parsing problem with CIGAR string...\n");
				return CIGAR_ERROR;
			}
			cigarCodes[numCodes] = str[i];
			str[i] = '\0';
			sscanf(start, "%d", &(cigarLens[numCodes]));
			//fprintf(stderr, " %c|%s|%d ", cigarCodes[numCodes],start,cigarLens[numCodes]);
			if (cigarCodes[numCodes] == 'M' || cigarCodes[numCodes] == '=') {
				if (firstStart < 0) {
					firstStart = totalLen;
					firstActive=1;
				}
				if (firstActive) {
					firstEnd = totalLen+cigarLens[numCodes];
				}
				if (lastActive==0) {
					lastStart = totalLen;
					lastActive=1;
				}
				if (lastActive) {
					lastEnd = totalLen+cigarLens[numCodes];
				}
				totalLen += cigarLens[numCodes];
			} else if (cigarCodes[numCodes] == 'D' || cigarCodes[numCodes] == 'N' || cigarCodes[numCodes] == 'X'
						|| cigarCodes[numCodes] == 'H' || cigarCodes[numCodes] == 'P') {
				totalLen += cigarLens[numCodes];
				firstActive = 0;
				lastActive = 0;
			} else if (cigarCodes[numCodes] == 'I' || cigarCodes[numCodes] == 'S') {
				totalLen -= cigarLens[numCodes];
			}
			start = &(str[i+1]);
			lastStartI = i+1;
			numCodes++;
		}
		i++;
	}
	int rv = totalLen;
	if (dir == STRAND_POSITIVE) {
		rv = firstStart;
		initLen = firstEnd-firstStart;
	} else {
		rv = lastEnd-1;
		initLen = lastEnd-lastStart;
	}
	//fprintf(stderr, "\t%d\t%d\t%d\n", rv, initLen,dir);
	return rv;
}

int TagLibrary::getQualityTotal(char* quality) {
	int sum = 0;
	if (quality == NULL) return sum;
	int i=0;
	while (quality[i] != '\0') {
		sum += quality[i++];
	}
	return sum;
}

void TagLibrary::cleanUpChrName(char* chr) {
	int len = strlen(chr);
	for (int i=0;i<len-2;i++) {
		if (chr[i] == '.') {
			if (chr[i+1] == 'f') {
				if (chr[i+2] == 'a') {
					chr[i] = '\0';
					return;
				}
			}
		}
	}
}



void TagLibrary::readSingleTagFile() {

	if (chrs != NULL && chrs->total > 0) {
		char** keys = chrs->keys();
		qsort(keys,chrs->total,sizeof(char*),&chrcmp);
		for (int i=0;i<chrs->total;i++) {
			ChrTags* ct = (ChrTags*)chrs->search(keys[i]);
			ct->freeTags();
			delete [](keys[i]);
		}
		delete []keys;
	}
		
	if (tagFile != NULL) {
		fclose(tagFile);
	}
	if (singleFilename == NULL) {
		fprintf(stderr, "Single tag file is not active\n");
		exit(0);
	}
	tagFile = fopen(singleFilename,"r");
	if (tagFile == NULL) {
		fprintf(stderr, "!!! Could not open %s for reading !!!\n",singleFilename);
		exit(0);
	}

	char* buf = new char[BUFFER];
	char* chr = NULL;
	char* name = new char[BUFFER];	
	char** line = new char*[BUFFER];
	int numCols = 0;
	int pos = 0;
	float tagCount= 0;
	int dir = 0;
	int len = -1;

	int optimizedFlag = 1;
	int numChangedChr = 0;
	char* lastChr = new char[10000];
	int lastPos = 0;


	while (fgets(buf, BUFFER, tagFile) != NULL) {
		//sscanf(buf, "%s\t%s\t%d\t%d\t%d",name,chr,&pos,&dir,&tagCount);
		split(buf, line, numCols,'\t');
		if (numCols < 5) continue;

		chr = line[1];
		sscanf(line[2],"%d",&pos);
		sscanf(line[3],"%d",&dir);
		sscanf(line[4],"%f",&tagCount);
		len = -1;
		if (numCols > 5) {
			sscanf(line[5],"%d",&len);
		}


		ChrTags* ct = (ChrTags*)chrs->search(chr);
		if (ct == NULL) {
			ct = new ChrTags(chr);
			ct->singleFile = 1;
			ct->mCflag = mCflag;
			ct->pairedEndFlag = pairedEndFlag;
			ct->tagfp = tagFile;
			chrs->insert(ct,chr);
		}
		if (pos > ct->appearentSize) {
			ct->appearentSize = pos;
		}
		ct->addTag(pos,(char)dir,len,tagCount);
		if (strcmp(chr,lastChr) != 0) {
			numChangedChr++;
			strcpy(lastChr,chr);
		} else {
			if (lastPos > pos) optimizedFlag = 0;
		}
		lastPos = pos;

	}
	fclose(tagFile);
	tagFile = NULL;

	if (numChangedChr > chrs->total) {
		optimizedFlag = 0;
	}
	if (optimizedFlag == 0) {
		fprintf(stderr, "\tOptimizing single genome.tags.tsv file...\n");
		tagFile = fopen(singleFilename,"w");
	}


	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);

	gsizeEstimate = 0;
	totalTags = 0.0;
	totalPositions = 0;
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*)chrs->search(keys[i]);
		ct->optimizedFlag=optimizedFlag;
		ct->optimizeTags();
		if (optimizedFlag == 0) {
			ct->tagfp = tagFile;
			ct->print();
		}
		ct->loaded =1;
		totalTags += ct->totalTags;
		ct->ogTotalTags = ct->totalTags;
		totalPositions += ct->totalPositions;
		gsizeEstimate += ct->appearentSize;
		//fprintf(stderr, "\t%s\t%lld\n", keys[i],ct->appearentSize);
	}
	if (optimizedFlag == 0) {
		fprintf(stderr, "\tEstimated genome size = %lld\n", gsizeEstimate);
		fclose(tagFile);
		tagFile = NULL;
	}

	delete []lastChr;
	delete []buf;
	delete []name;
	delete []line;
}



void TagLibrary::readNewTagFile(char* filename) {


	FILE* fp = fopen(filename,"r");
	if (fp == NULL) {
		fprintf(stderr, "Could not open %s\n", filename);
		return;
	}

	char* buf = new char[BUFFER];
	char* chr = NULL;
	char* name = NULL;
	char** line = new char*[BUFFER];
	int numCols = 0;
	int pos = 0;
	int pos2 = 0;
	float tagCount= 0;
	int dir = 0;
	int dir2 = 0;
	int len = -1;
	int len2 = -1;

	PETag* petag = new PETag();

	//pairedEndFlag = -1;

	while (fgets(buf, BUFFER, fp) != NULL) {
		//sscanf(buf, "%s\t%s\t%d\t%d\t%d",name,chr,&pos,&dir,&tagCount);
		split(buf, line, numCols,'\t');
		if (numCols < 5) continue;

		name = line[0];
		chr = line[1];
		sscanf(line[2],"%d",&pos);
		sscanf(line[3],"%d",&dir);
		sscanf(line[4],"%f",&tagCount);
		len = -1;
		if (numCols > 5) {
			sscanf(line[5],"%d",&len);
		}
		if (numCols < 9) {
			if (len < minReadLength || len > maxReadLength) continue;
			addAlignedTag(line[0],line[1],pos,(char)dir,len,tagCount,TAGLIBRARY_SINGLE_FLAG);
		} else {
			petag->name = name;
			petag->chr1 = chr;
			petag->p1 = pos;
			petag->d1 = dir;
			petag->v = tagCount;
			petag->len1 = len;
			sscanf(line[7],"%d",&pos2);
			sscanf(line[8],"%d",&dir2);
			if (numCols > 9) {
				sscanf(line[9],"%d",&len2);
			}
			petag->chr2 = line[6];
			petag->p2 = pos2;
			petag->d2 = dir2;
			petag->len2 = len2;
			if (len2 < minReadLength || len2 > maxReadLength) continue;
			addAlignedPETag(petag);
			pairedEndFlag = 1;
		}
	}
	fclose(fp);

	delete petag;
	delete []buf;
	delete []line;
}

void TagLibrary::readTagInfo() {

	if (directory == NULL) {
		fprintf(stderr, "No directory specified!\n");
		return;
	}
	char* buf = new char[BUFFER];
	char** line = new char*[BUFFER];
	int numCols=0;

	strcpy(buf, directory);
	strcat(buf, "/tagInfo.txt");
	FILE* fp = fopen(buf, "r");
	if (fp == NULL) {
		fprintf(stderr, "!!! Could not open tag information file: %s\n", buf);
		fprintf(stderr, "Probably not a valid tag directory.  Quitting...\n");
		delete []buf;
		exit(0);
		return;
	} 
	
	while (fgets(buf, BUFFER, fp) != NULL) {
		split(buf,line,numCols,'\t');

		if (strncmp(line[0],"genome=",7) == 0) {
			setGenome(&(line[0][7]));
			if (numCols > 2) {
				sscanf(line[1],"%lld",&totalPositions);
				sscanf(line[2],"%lf",&totalTags);
			}
		} else if (strncmp(line[0],"pairedEnd=true",14) == 0) {
			pairedEndFlag = 1;
		} else if (strncmp(line[0],"mC=true",7) == 0) {
			mCflag = 1;
		} else if (strncmp(line[0],"restrictionSite=",16) == 0) {
			char* rs = &(line[0][16]);
			setRestrictionSite(rs);
		} else if (strncmp(line[0],"name=",5) == 0) {
			setName(&(line[0][5]));
		} else if (strncmp(line[0],"cmd=",4) == 0) {
			setCMD(&(line[0][4]));
		} else if (strncmp(line[0],"singleTagFile=true",18) == 0) {
			singleFile = 1;
			singleFilename = new char[10000];
			sprintf(singleFilename,"%s/genome.tags.tsv",directory);
		} else if (strncmp(line[0],"fragmentLengthEstimate=",23) == 0) {
			if (strncmp(&(line[0][23]), "given",5)==0) {
				fragmentLengthSetting = FRAGMENT_LEN_GIVEN;
				fragmentLengthEstimate = FRAGMENT_LEN_GIVEN;
			} else {
				sscanf(&(line[0][23]),"%d",&fragmentLengthEstimate);
				fragmentLengthSetting = fragmentLengthEstimate;
			}
		} else if (strncmp(line[0],"gsizeEstimate=",14) == 0) {
			sscanf(&(line[0][14]),"%lld",&gsizeEstimate);
		} else if (strncmp(line[0],"peakSizeEstimate=",17) == 0) {
			sscanf(&(line[0][17]),"%d",&peakSizeEstimate);
		} else if (strncmp(line[0],"averageTagLength=",17) == 0) {
			sscanf(&(line[0][17]),"%lf",&averageTagLength);
		} else if (strncmp(line[0],"averageTagsPerPosition=",23) == 0) {
			sscanf(&(line[0][23]),"%lf",&averageTagsPerPosition);
		} else if (strncmp(line[0],"Estimated Fragment Length=",26) == 0) {
			if (strncmp(&(line[0][26]), "given",5)==0) {
				fragmentLengthSetting = FRAGMENT_LEN_GIVEN;
				fragmentLengthEstimate = FRAGMENT_LEN_GIVEN;
			} else {
				sscanf(&(line[0][26]),"%d",&fragmentLengthEstimate);
				fragmentLengthSetting = fragmentLengthEstimate;
			}
			char** og = new char*[10];
			int numOG = 0;
			split(line[0],og, numOG, ',');
			if (numOG > 1) {
				if (strncmp(og[1]," Estimated Peak Width=",22) == 0) {
					sscanf(&(og[1][22]),"%d",&peakSizeEstimate);
				}
			}
			delete []og;
		} else if (strncmp(line[0],"genome",6) == 0 && numCols > 2) {
			sscanf(line[1],"%lld",&totalPositions);
			sscanf(line[2],"%lf",&totalTags);
		}
	}

	if (fragmentLengthEstimate < 0 || fragmentLengthEstimate > 10000)  {
		if (fragmentLengthEstimate == TAGADJUST_AUTO) {
			fprintf(stderr, "\tCouldn't figure fragment length (tagInfo.txt file parsing error) setting to %d\n",
												 TAGADJUST_DEFAULT*2);
			fragmentLengthEstimate = TAGADJUST_DEFAULT*2;
		} else if (fragmentLengthSetting == FRAGMENT_LEN_GIVEN) {
			fprintf(stderr, "\tUsing actual fragment/read lengths\n");
		} else {
			fprintf(stderr, "!!! Do you mean for the fragment length estimate to be %d?? If not, set manually!!\n",
					fragmentLengthEstimate);
		}
	}
	//fprintf(stderr, "%d\t%lf\n", totalPositions,totalTags);
	//fprintf(stderr, "%d\t%d\n", fragmentLengthEstimate, peakSizeEstimate);
	tbp = totalTags/((double)2e9);

	if (singleFile) {
		//preadSingleTagFile();
	}

	delete []buf;
	delete []line;

}


void TagLibrary::optimizeTagFiles() {

	fprintf(stderr, "\n\tOptimizing tag files...\n");
	totalTags = 0;
	totalPositions = 0;

	setMaxTBP(maxtbp);
	setMinTBP(mintbp);

	if (singleFile) {
		readSingleTagFile();
		return;
	}

	char** chrKeys = chrs->keys();
	qsort(chrKeys,chrs->total,sizeof(char*),&chrcmp);
	gsizeEstimate = 0;
	for (int i=0;i<chrs->total;i++) {
		//fprintf(stderr, "Optimizing %s\n", chrKeys[i]);
		ChrTags* ct = (ChrTags*)chrs->search(chrKeys[i]);
		ct->optimizeTagFile();
		totalTags += ct->totalTags;
		totalPositions += ct->totalPositions;
		gsizeEstimate += ct->appearentSize;
		//fprintf(stderr, "\t\t%s: %lld\n", chrKeys[i], ct->appearentSize);
		delete [](chrKeys[i]);
	}
	fprintf(stderr, "\tEstimated genome size = %lld\n", gsizeEstimate);

	if (gsizeEstimate > 1) {
		tbp = totalTags/gsizeEstimate;	
		fprintf(stderr, "\tEstimated average read density = %.6lf per bp\n", tbp);
	}

	delete []chrKeys;

}

void TagLibrary::addAlignedPETag(PETag* petag) {
	char* chr = petag->chr1;
	ChrTags* chrtags = (ChrTags*) chrs->search(chr);
	if (chrtags == NULL) {
		chrtags = new ChrTags(chr);
		chrtags->pairedEndFlag = 1;
		chrs->insert(chrtags,chr);

		if (singleFile == 0) {
			char* file = getTagFileName(chr);
			chrtags->setTagFile(file);
			chrtags->openTagFile((char*)"w");
			delete []file;
		} else {
			if (tagFile == NULL) {
				singleFilename = new char[10000];
				sprintf(singleFilename,"%s/genome.tags.tsv",directory);
				tagFile = fopen(singleFilename, "w");
			}
			chrtags->tagfp = tagFile;
			chrtags->singleFile = 1;
		}
	}
	chrtags->printAlignedPETag(petag, 0);
	totalTags += petag->v;
}

void TagLibrary::addAlignedTag(char* name,char* chr,int pos,char dir, int length, 
														float value,int PEflag) {

	static int warning = 0;

	if (pos < 1) pos = 1;
	ChrTags* chrtags = (ChrTags*) chrs->search(chr);
	if (chrtags == NULL) {
		chrtags = new ChrTags(chr);
		chrtags->pairedEndFlag = pairedEndFlag;
		chrtags->mCflag = mCflag;
		chrs->insert(chrtags,chr);

		if (singleFile == 0) {
			if (chrs->total > 1000 && warning == 0) {
				warning = 1;
				fprintf(stderr, "!!! Warning: more than 1000 chromosomes detected.  This can casue file I/O problems\n");
				fprintf(stderr, "!!! If the command fails, consider reruning makeTagDirectories with \"-single\"\n");
				fprintf(stderr, "!!! so that it uses a single tag file instead of one per chromosome\n");
			}

			char* file = getTagFileName(chr);
//fprintf(stderr, "||%s||\n", file);
			chrtags->setTagFile(file);
			chrtags->openTagFile((char*)"w");
			delete []file;
		} else {

			if (tagFile == NULL) {
				char* file = new char[10000];
				sprintf(file, "%s/genome.tags.tsv", directory);
				tagFile = fopen(file, "w");
				delete []file;
			}
			chrtags->tagfp = tagFile;
			chrtags->singleFile = 1;
		}
	}
	if (pairedEndFlag && chrtags->tagfpR1 == NULL) {
		char* readFile = new char[10000];
		sprintf(readFile, "%s.R1", chrtags->tagFile);
		chrtags->tagfpR1 = fopen(readFile, "w");
		sprintf(readFile, "%s.R2", chrtags->tagFile);
		chrtags->tagfpR2 = fopen(readFile, "w");
		delete []readFile;
	}
	chrtags->printAlignedTag(name,pos,dir,length,value,PEflag);
	totalTags += value;
}

char* TagLibrary::getTagFileName(char* chrname) {
	char* file = new char[10000];
	sprintf(file, "%s/%s.tags.tsv",directory,chrname);
	return file;
}
char* TagLibrary::getDirFileName(char* filename) {
	char* file = new char[10000];
	strcpy(file,directory);
	strcat(file, "/");
	strcat(file, filename);
	return file;
}


void TagLibrary::addTag(char* chr,int pos,char dir, int length, float value) {
	ChrTags* chrtags = (ChrTags*) chrs->search(chr);
	if (chrtags == NULL) {
		chrtags = new ChrTags(chr);
		char* file = getTagFileName(chr);
		chrtags->setTagFile(file);
		chrs->insert(chrtags,chr);
	}
	chrtags->addTag(pos,dir,length,value);
}


void TagLibrary::printTagInfo() {
	printTagInfo(NULL);
}
void TagLibrary::printTagInfo(FILE* fp) {

	char* filename = new char[10000];
	strcpy(filename, directory);
	strcat(filename, "/tagInfo.txt");

	FILE* info = NULL;
	if (fp == NULL) {
		info = fopen(filename,"w");
	} else {
		info = fp;
	}
	if (info == NULL) {
		fprintf(stderr, "Could not open file for output in directory %s\n", directory);
		return;
	}
	
	if (name != NULL) {
		fprintf(info, "name=%s\tUnique Positions\tTotal Tags\n",name);
	} else {
		fprintf(info, "name\tUnique Positions\tTotal Tags\n");
	}

	if (genome != NULL) {
		fprintf(info, "genome=%s\t%lld\t%.1lf\n", genome, totalPositions, totalTags);
	} else {
		fprintf(info, "genome\t%lld\t%.1lf\n", totalPositions, totalTags);
	}
	if (pairedEndFlag) {
		fprintf(info, "pairedEnd=true\n");
	}
	if (mCflag) {
		fprintf(info, "mC=true\n");
	}
	if (restrictionSite != NULL) {
		fprintf(info, "restrictionSite=%s\n", restrictionSite);
	}
	fprintf(info, "fragmentLengthEstimate=%d\t\t\n",fragmentLengthEstimate);
	fprintf(info, "peakSizeEstimate=%d\t\t\n",peakSizeEstimate);
	fprintf(info, "tagsPerBP=%lf\t\t\n",tbp);
	fprintf(info, "averageTagsPerPosition=%.3lf\t\t\n",averageTagsPerPosition);
	fprintf(info, "averageTagLength=%.3lf\t\t\n",averageTagLength);
	fprintf(info, "gsizeEstimate=%lld\t\t\n",gsizeEstimate);
	fprintf(info, "averageFragmentGCcontent=%.3lf\t\t\n",gcAverage);
	if (singleFile) {
		fprintf(info, "singleTagFile=true\t\t\n");
	}
	if (gcNormalized) {
		fprintf(info, "Tag Directory has been GC-normalized\t\t\n");
	}
	if (oligoNormalized) {
		fprintf(info, "Tag Directory has been oligo-normalized\t\t\n");
	}


	char** chr = chrs->keys();
	qsort(chr,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*)chrs->search(chr[i]);
		fprintf(info,"%s\t%d\t%.1f\n",chr[i],ct->totalPositions,ct->totalTags);
		delete [](chr[i]);
	}

	if (cmd != NULL) {
		fprintf(info, "cmd=%s\t\t\n", cmd);
	}
	if (fp == NULL) fclose(info);
	delete []chr;
	delete []filename;

}


double* TagLibrary::getPETagDistribution(int windowSize, int largeWindowSize, int largeResolution, 
							char* outputPrefix, int &distLength) {


	FILE *fp = NULL;
	FILE *localfp = NULL;

	if (outputPrefix != NULL) {	
		char* fname = new char[100000];
		sprintf(fname, "%s.FreqDistribution_%d.txt",outputPrefix,largeResolution);
		char* file = getDirFileName(fname);
		fp = fopen(file, "w");
		if (fp == NULL) {
			fprintf(stderr, "Cannot open %s for writing\n", file);
			return NULL;
		}
		delete []file;
	
		sprintf(fname, "%s.LocalDistribution.txt",outputPrefix);
		file = getDirFileName(fname);
		localfp = fopen(file, "w");
		if (localfp == NULL) {
			fprintf(stderr, "Cannot open %s for writing\n", file);
			return NULL;
		}
		delete []fname;
	}

	double* sameStrand = new double[windowSize];
	double* diffStrand = new double[windowSize];
	double* smoothed = new double[windowSize];
	for (int i=0;i<windowSize;i++) {
		sameStrand[i] = 0.0;
		diffStrand[i] = 0.0;
		smoothed[i] = 0.0;
	}


	int largeLength = (largeWindowSize/largeResolution)+2;
	double* largeWindow = new double[largeLength];
	double* largeWindowN = new double[largeLength];
	for (int i=0;i<largeLength;i++) {
		largeWindow[i] = 0.0;
		largeWindowN[i] = 0.0;
	}
	distLength = largeLength;


	//int genomeIndexSize = gsizeEstimate/largeResolution+1;

	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	if (singleFile) readSingleTagFile();
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*) chrs->search(keys[i]);
		int chrLength  = ct->getPETagDistribution(sameStrand, diffStrand,windowSize,
											largeWindow,largeResolution,largeLength);
		delete [](keys[i]);
		int indexLen = chrLength/largeResolution+1;
		for (int j=0;j<indexLen;j++) {
			int z = j;
			if (j > largeLength-2) z = largeLength-2;
			largeWindowN[z] += (double)(indexLen-j);
		}
		largeWindowN[largeLength-1] += (double)indexLen;
	}
	delete []keys;

	double sum = 0.0;
	for (int i=0;i<largeLength;i++) {
		if (largeWindowN[i] > 0.0) {
			largeWindow[i] /= largeWindowN[i];
		}
		sum += largeWindow[i];
	}
	for (int i=0;i<largeLength;i++) {
		largeWindow[i]/=sum;
	}
	

	//int minFragLen = AUTOCORRELATION_OFFSETMIN;
	//if (averageTagLength > 0.0) minFragLen = (int) averageTagLength+3;

	fragmentLengthEstimate = 0;
	peakSizeEstimate = 0;
	
	if (localfp != NULL) {
		int halfWindow = windowSize/2;
		fprintf(localfp, "Local Distance in bp between PE tags\tSame Strand\tOpposite Strands\n");
		for (int i=0;i<windowSize;i++) {
			int offset = i-halfWindow;
			fprintf(localfp, "%d\t%lf\t%lf\n",offset,sameStrand[i],diffStrand[i]);
		}
		fclose(localfp);
	}
	if (fp != NULL) {
		fprintf(fp, "Distance between PE tags\tFraction of total PE tags");
		//fprintf(fp, "(Interchromosomal:%le)\n", largeWindow[largeLength-1]/totalTags);
		fprintf(fp, "(Interchromosomal:%le)\n", largeWindow[largeLength-1]);
		for (int i=0;i<largeLength-1;i++) {
			if (i==largeLength-2) fprintf(fp, "More than ");
			//fprintf(fp, "%d\t%le\n", i*largeResolution,largeWindow[i]/totalTags);
			fprintf(fp, "%d\t%le\n", i*largeResolution,largeWindow[i]);
		}
		fclose(fp);
	}

	delete []sameStrand;
	delete []diffStrand;
	delete []smoothed;

	return largeWindow;
}

void TagLibrary::autoCorrelateTags(FILE* nfp, int windowSize, double maxTags) {
	FILE *fp = nfp;
	if (fp == NULL) {
		char* file = getDirFileName((char*)"tagAutocorrelation.txt");
		fp = fopen(file, "w");
		if (fp == NULL) {
			fprintf(stderr, "Cannot open %s for writing\n", file);
			return;
		}
	}
	if (mCflag) {
		windowSize = 20000;
	}

	int avgWindow = AUTOCORRELATION_HALFWINDOW;
	double totalCount = 0;
	double* sameStrand = new double[windowSize+1];
	double* diffStrand = new double[windowSize+1];
	double* sameStrandN = new double[windowSize+1];
	double* diffStrandN = new double[windowSize+1];
	double* smoothed = new double[windowSize+1];
	for (int i=0;i<windowSize+1;i++) {
		sameStrand[i] = 0.0;
		diffStrand[i] = 0.0;
		sameStrandN[i] = 0.0;
		diffStrandN[i] = 0.0;
		smoothed[i] = 0.0;
	}

	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	if (singleFile) readSingleTagFile();
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*) chrs->search(keys[i]);
		ct->autoCorrelateTags(sameStrand, diffStrand,windowSize,maxTags,totalCount,
									sameStrandN, diffStrandN, mCflag);
		delete [](keys[i]);
		if (totalCount > maxTags) break;
	}
	delete []keys;

	if (mCflag) {
		for (int i=0;i<windowSize+1;i++) {
			if (diffStrandN[i] > 0) diffStrand[i] /= diffStrandN[i];
			if (sameStrandN[i] > 0) sameStrand[i] /= sameStrandN[i];
		}
	}
				


	int halfWindow = windowSize/2;
	double backEst = 0;
	double maximum = 0;
	int lastEstimate = 0;

	int minFragLen = AUTOCORRELATION_OFFSETMIN;
	if (averageTagLength > 0.0) minFragLen = (int) averageTagLength+3;
	int bottomFlag = 0;

	fragmentLengthEstimate = 0;
	peakSizeEstimate = 0;
	
	for (int i=0;i<windowSize;i++) {
		int offset = i-halfWindow;
		double avg = 0;
		double n= 0;
		for (int j=-avgWindow;j<avgWindow;j++) {
			if (j+i < 0) continue;
			if (j+i >= windowSize) break;
			avg += diffStrand[j+i];
			n+=1.0;
		}
		if (n>0.5) avg /= n;
		smoothed[i] = avg;

		if (offset == AUTOCORRELATION_BACKOFFSET-AUTOCORRELATION_HALFWINDOW) {
			backEst = smoothed[i];
		}
		if (offset > minFragLen) {
			if (smoothed[i] > maximum && bottomFlag) {
				fragmentLengthEstimate = offset;
				maximum = smoothed[i];
			} else if (bottomFlag == 0) {
				if (smoothed[i] > smoothed[i-1]) bottomFlag =1;
			}
			if (smoothed[i] < backEst) {
				if (lastEstimate != fragmentLengthEstimate) {
					peakSizeEstimate = offset-fragmentLengthEstimate;
					lastEstimate = fragmentLengthEstimate;
				}
			}
		}		
	}

	if (mCflag) {
		fragmentLengthEstimate = 1;
		peakSizeEstimate = 1000;
	} else {
		fprintf(stderr, "\tFragment Length Estimate: %d\n", fragmentLengthEstimate);

		if (fragmentLengthEstimate <= minFragLen) {
			fprintf(stderr, "\t\t!!! No reliable estimate for fragment length\n"); 
			fprintf(stderr, "\t\t!!! PLEASE SET MANUALLY (currently set to 150, edit tagInfo.txt to change)\n"); 
			fragmentLengthEstimate = 150;
		}
		fprintf(stderr, "\tPeak Width Estimate: %d\n", peakSizeEstimate);
		if (peakSizeEstimate < fragmentLengthEstimate) {
			peakSizeEstimate = fragmentLengthEstimate;
			fprintf(stderr, "\t\t!!! No reliable estimate for peak size\n");
			fprintf(stderr, "\t\tSetting Peak width estimate to be equal to fragment length estimate\n");
		}
	}


	double posBackSignal = 0.0;
	double negBackSignal = 0.0;
	double posBackN=0.0;
	double negBackN=0.0;
	double posSignal = 0.0;
	double negSignal = 0.0;
	double posN=0.0;
	double negN=0.0;

	for (int i=0;i<windowSize;i++) {
		int offset = i-halfWindow;
		if (offset < fragmentLengthEstimate*-1 || offset > fragmentLengthEstimate*2) {
			posBackSignal += sameStrand[i];
			posBackN += 1.0;
			negBackSignal += diffStrand[i];
			negBackN += 1.0;
		} else {
			posSignal += sameStrand[i];
			posN += 1.0;
			negSignal += diffStrand[i];
			negN += 1.0;
		}
	}
		
	/*for (int i=0;i<windowSize;i++) {
		int offset = i-halfWindow;
		if (i < AUTOCORRELATION_ENRICHMENT || i > windowSize-AUTOCORRELATION_ENRICHMENT) {
			posBackSignal += sameStrand[i];
			posBackN += 1.0;
			negBackSignal += diffStrand[i];
			negBackN += 1.0;
		} else if ((offset >= 0 && offset*2 < AUTOCORRELATION_ENRICHMENT)
					||  (offset < 0 && offset*-2 < AUTOCORRELATION_ENRICHMENT)) {
			posSignal += sameStrand[i];
			posN += 1.0;
		} else if ((offset >= fragmentLengthEstimate 
						&& (offset-fragmentLengthEstimate)*2 < AUTOCORRELATION_ENRICHMENT)
					||  (offset < fragmentLengthEstimate 
						&& (offset-fragmentLengthEstimate)*-2 < AUTOCORRELATION_ENRICHMENT)) {
			negSignal += diffStrand[i];
			negN += 1.0;
		}
	}*/

	if (posN > 0.0) posSignal /= posN;
	if (negN > 0.0) negSignal /= negN;
	if (posBackN > 0.0) posBackSignal /= posBackN;
	if (negBackN > 0.0) negBackSignal /= negBackN;
	double posEnrichment = 1.0;
	double negEnrichment = 1.0;
	double posVsNegEnrichment = 1.0;
	if (posBackSignal > 0.0) posEnrichment = posSignal/posBackSignal;
	if (negBackSignal > 0.0) negEnrichment = negSignal/negBackSignal;
	if (negSignal > 0.0) posVsNegEnrichment = posSignal/negSignal;

	if (mCflag==0) {
		fprintf(stderr, "\tAutocorrelation quality control metrics:\n");
		fprintf(stderr, "\t\tSame strand fold enrichment: %.1lf\n", posEnrichment);
		fprintf(stderr, "\t\tDiff strand fold enrichment: %.1lf\n", negEnrichment);
		fprintf(stderr, "\t\tSame / Diff fold enrichment: %.1lf\n", posVsNegEnrichment);
		fprintf(stderr, "\n");
		if (posVsNegEnrichment > AUTOCORRELATION_RNA_THRESHOLD) {
			fprintf(stderr, "\t\tGuessing sample is strand specific RNA-Seq\n");
			fprintf(stderr, "\t\tSetting fragment length estimate to 75, edit tagInfo.txt to change\n");
			fragmentLengthEstimate = 75;
			peakSizeEstimate = fragmentLengthEstimate;
		} else if (posEnrichment > AUTOCORRELATION_CHIP_THRESHOLD
						|| negEnrichment > AUTOCORRELATION_CHIP_THRESHOLD) {
			if (negEnrichment <= 0.01 || fabs(log(posEnrichment/negEnrichment))>0.50) {
				fprintf(stderr, "\t\tGuessing sample is ChIP-Seq - uneven enrichment between same strand and\n");
				fprintf(stderr, "\t\tdifferent strands - may have problems such as clonal amplification.\n");
				if (fragmentLengthEstimate > 400) {
					fprintf(stderr, "\t\tMay be strand-unspecific RNA-Seq (i.e. fragmented cDNA lib).  If so,\n");
					fprintf(stderr, "\t\tmanually set the fragment length estimate in tagInfo.txt\n");
				}
			} else {
				fprintf(stderr, "\t\tGuessing sample is ChIP-Seq - autocorrelation looks good.\n");
			}
		} else {
			fprintf(stderr, "\t\tGuessing sample is ChIP-Seq - may have low enrichment with lots of background\n");
		}
		fprintf(stderr, "\n");
	
		fprintf(fp, "Distance in bp(Fragment Length Estimate: %d", fragmentLengthEstimate);
		fprintf(fp, ")(Peak Width Estimate: %d)", peakSizeEstimate);
		fprintf(fp, "\tSame Strand (+ for Watson strand, - for Crick)\tOpposite Strand\n");
	} else {
		fprintf(fp, "Distance between nucleotides");
		fprintf(fp, "\tAvg abs change in mC%%, Same Strand (+ for Watson strand, - for Crick)\tAvg abs change in mC%%, Opposite Strand\n");
	}

	

	for (int i=0;i<windowSize;i++) {
		int offset = i-halfWindow;
		if (mCflag) {
			fprintf(fp, "%d\t%lf\t%lf\n", offset,sameStrand[i],diffStrand[i]);
		} else {
			fprintf(fp, "%d\t%.1f\t%.1f\n", offset,sameStrand[i],diffStrand[i]);
		}
	}
		
	if (nfp == NULL) {
		fclose(fp);
	}
	
	delete []sameStrand;
	delete []diffStrand;
	delete []sameStrandN;
	delete []diffStrandN;
	delete []smoothed;

}

void TagLibrary::setMaxTBP(float max) {

	maxtbp = max;
	
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*) chrs->search(keys[i]);
		ct->setMaxTBP(max);
		delete [](keys[i]);
	}
	if (keys != NULL) delete []keys;
}

void TagLibrary::setMinTBP(float min) {

	mintbp = min;
	
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*) chrs->search(keys[i]);
		ct->setMinTBP(min);
		delete [](keys[i]);
	}
	if (keys != NULL) delete []keys;
}

void TagLibrary::setTagAdjust(int dist) {
	if (dist == TAGADJUST_AUTO) {
		if (0) { //fragmentLengthEstimate <= AUTOCORRELATION_OFFSETMIN) {
			dist = TAGADJUST_DEFAULT;
			fragmentLengthEstimate = dist*2;
		} else {
			dist = (fragmentLengthEstimate)/2;
		}
	}
	tagAdjust = dist;
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*) chrs->search(keys[i]);
		ct->setTagAdjust(dist);
		ct->revStrand = revStrand;
		delete [](keys[i]);
	}
	delete []keys;
}

PeakLibrary* TagLibrary::findPutativePeaks(int peakSize, int minDist,char strand,float minCount) {

	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	PeakLibrary* putativePeaks = new PeakLibrary(10000000);
	if (singleFile) readSingleTagFile();
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*) chrs->search(keys[i]);
		ct->loadTags();
		ct->findPutativePeaks(putativePeaks,peakSize,minDist,strand,minCount);
		ct->freeTags();
		delete [](keys[i]);
	}
	delete []keys;

	putativePeaks->sortChr();
	return putativePeaks;
}

PeakLibrary* TagLibrary::findmCPeaks(int peakSize, char strand,int mCflag, double mCthresh, 
															int minNumC, TagLibrary* input) {


	PeakLibrary* peaks = new PeakLibrary();

	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);

	if (singleFile) readSingleTagFile();
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*) chrs->search(keys[i]);
		ChrTags* ctinput = NULL;
		if (input != NULL) {
			ctinput = (ChrTags*) input->chrs->search(keys[i]);
			ctinput->loadTags();
		}
		ct->findmCPeaks(peaks,peakSize,strand,mCflag,mCthresh, minNumC,ctinput);
	}
	fprintf(stderr, "\n\tTotal Regions = %d\n\n", peaks->numPeaks);
	return peaks;
}

//GRO-Seq transcript identification


PeakLibrary* TagLibrary::findGroSeqTranscripts(int mode, TagLibrary* input, char* uniqMapDirectory, char strand,
				int tssSize, int minBodySize, int maxBodySize, double threshold, double foldTranscriptStart,
				double foldTranscriptBody, double endFold, double pseudoTags, double inputFold,int groseqMethod) {

	int inputFragLen = 0;
	double inputNorm = 1.0;
	if (input != NULL) {
		inputFragLen = input->fragmentLengthSetting;
		inputNorm = totalTags/input->totalTags;
	}
	
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	PeakLibrary* putativePeaks = new PeakLibrary();

	if (singleFile) readSingleTagFile();
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*) chrs->search(keys[i]);
		ct->loadTags();
		ChrTags* ctinput = NULL;
		if (input != NULL) {
			ctinput = (ChrTags*) input->chrs->search(keys[i]);
			ctinput->loadTags();
		}
		UniqMapChrs * umc = NULL;
		if (uniqMapDirectory != NULL) {
			umc = new UniqMapChrs(keys[i],uniqMapDirectory,0);
		}

		int curStrand = strand;
		if (strand == STRAND_SEPARATE) {
			curStrand = STRAND_POSITIVE;
		}
		fprintf(stderr, "\t\t%s ...", keys[i]);
	
		
		PeakLibrary* ctpeaks = ct->findGroSeqTranscripts(mode,ctinput,umc,curStrand,tssSize,
									minBodySize,maxBodySize,
									threshold,foldTranscriptStart, foldTranscriptBody,endFold,inputFold,
									inputNorm,fragmentLengthEstimate,inputFragLen, pseudoTags,groseqMethod);
		if (ctpeaks != NULL) {
			putativePeaks->addPeakLibrary(ctpeaks);
			fprintf(stderr, " %d", ctpeaks->numPeaks);
			delete ctpeaks;
		}

		if (strand == STRAND_SEPARATE) {
			curStrand = STRAND_NEGATIVE;
			PeakLibrary* ctpeaks = ct->findGroSeqTranscripts(mode,ctinput,umc,curStrand,tssSize,
									minBodySize,maxBodySize,
									threshold,foldTranscriptStart, foldTranscriptBody,endFold,inputFold,
									inputNorm,fragmentLengthEstimate,inputFragLen, pseudoTags,groseqMethod);
			if (ctpeaks != NULL) {
				fprintf(stderr, "+ %d-", ctpeaks->numPeaks);
				putativePeaks->addPeakLibrary(ctpeaks);
				delete ctpeaks;
			}
		}
		fprintf(stderr, "\n");


		ct->freeTags();
		if (ctinput != NULL) {
			ctinput->freeTags();
			delete ctinput;
		}
		if (umc != NULL) {
			delete umc;
		}
		delete [](keys[i]);
	}
	delete []keys;

	putativePeaks->sortChr();
	return putativePeaks;
}



PeakLibrary* TagLibrary::findGroSeqRegions(int mode, TagLibrary* input, char* uniqMapDirectory, char strand,
				int tssSize, int minBodySize, int maxBodySize, double threshold, double foldTranscriptStart,
				double foldTranscriptBody, double endFold, double pseudoTags, double inputFold) {

	int inputFragLen = 0;
	double inputNorm = 1.0;
	if (input != NULL) {
		inputFragLen = input->fragmentLengthSetting;
		inputNorm = totalTags/input->totalTags;
	}
	
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	PeakLibrary* putativePeaks = new PeakLibrary();

	if (singleFile) readSingleTagFile();
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*) chrs->search(keys[i]);
		ct->loadTags();
		ChrTags* ctinput = NULL;
		if (input != NULL) {
			ctinput = (ChrTags*) input->chrs->search(keys[i]);
			ctinput->loadTags();
		}
		UniqMapChrs * umc = NULL;
		if (uniqMapDirectory != NULL) {
			umc = new UniqMapChrs(keys[i],uniqMapDirectory,0);
		}

		int curStrand = strand;
		if (strand == STRAND_SEPARATE) {
			curStrand = STRAND_POSITIVE;
		}
		fprintf(stderr, "\t\t%s ...", keys[i]);
	
		
		PeakLibrary* ctpeaks = ct->findGroSeqRegions(mode,ctinput,umc,curStrand,tssSize,
									minBodySize,maxBodySize,
									threshold,foldTranscriptStart, foldTranscriptBody,endFold,inputFold,
									inputNorm,fragmentLengthEstimate,inputFragLen, pseudoTags);
		if (ctpeaks != NULL) {
			putativePeaks->addPeakLibrary(ctpeaks);
			fprintf(stderr, " %d", ctpeaks->numPeaks);
			delete ctpeaks;
		}

		if (strand == STRAND_SEPARATE) {
			curStrand = STRAND_NEGATIVE;
			PeakLibrary* ctpeaks = ct->findGroSeqRegions(mode,ctinput,umc,curStrand,tssSize,
									minBodySize,maxBodySize,
									threshold,foldTranscriptStart, foldTranscriptBody,endFold,inputFold,
									inputNorm,fragmentLengthEstimate,inputFragLen, pseudoTags);
			if (ctpeaks != NULL) {
				fprintf(stderr, "+ %d-", ctpeaks->numPeaks);
				putativePeaks->addPeakLibrary(ctpeaks);
				delete ctpeaks;
			}
		}
		fprintf(stderr, "\n");


		ct->freeTags();
		if (ctinput != NULL) {
			ctinput->freeTags();
			delete ctinput;
		}
		if (umc != NULL) {
			delete umc;
		}
		delete [](keys[i]);
	}
	delete []keys;

	putativePeaks->sortChr();
	return putativePeaks;
}

double* TagLibrary::getTagLengthDistribution(FILE* nfp, int &max) {
	FILE *fp = nfp;
	if (fp == NULL) {
		char* file = NULL;
		if (mCflag) {
			file = getDirFileName((char*)"mCreadCoverageDistribution.txt");
		} else {
			file = getDirFileName((char*)"tagLengthDistribution.txt");
		}
			
		fp = fopen(file, "w");
		if (fp == NULL) {
			fprintf(stderr, "Cannot open %s for writing\n", file);
			return NULL;
		}
	}

	if (max == 0) {
		max = MAX_READ_LENGTH;
	}
	double* dist = new double[max];
	for (int i=0;i<max;i++) {
		dist[i] = 0;
	}


	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	if (singleFile) readSingleTagFile();
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*) chrs->search(keys[i]);
		ct->getTagLengthDistribution(dist,max);
		delete [](keys[i]);
	}
	delete []keys;


	float totalDist = 0;
	averageTagLength = 0.0;
	int maxValue = 0;
	for (int i=0;i<max;i++) {
		dist[i] /= (double)totalTags;
		if (dist[i] > 0) maxValue = i;
		averageTagLength += ((double)i)*dist[i];
		totalDist += dist[i];
	}

	if (mCflag) {
		fprintf(stderr, "\tAverage read depth (from read depth of data exceeding %.1lf threshold) = %.1lf\n", 
							parseAlignmentCpGMinValue,averageTagLength);
		fprintf(fp, "Read depth (Average read depth = %lf)", averageTagLength);
		fprintf(fp, "\tFraction of Positions\n");
	} else {
		fprintf(stderr, "\tAverage tag length = %.1lf\n", averageTagLength);
		fprintf(fp, "Tag Length (Average tag length = %lf)", averageTagLength);
		fprintf(fp, "\tFraction of Tags\n");
	}
	for (int i=0;i<max;i++) {
		fprintf(fp,"%d\t%lf\n",i,dist[i]);
		if (i>=maxValue) break;
	}
	if (nfp == NULL) {
		fclose(fp);
	}
	return dist;

}




double* TagLibrary::getTagCountDistribution(FILE* nfp, int &max) {
	FILE *fp = nfp;
	if (fp == NULL) {
		char* file = NULL;
		if (mCflag) {
			file = getDirFileName((char*)"mCratioDistribution.txt");
		} else {
			file = getDirFileName((char*)"tagCountDistribution.txt");
		}
		fp = fopen(file, "w");
		if (fp == NULL) {
			fprintf(stderr, "Cannot open %s for writing\n", file);
			return NULL;
		}
	}

	if (max == 0) {
		max = MAX_TAGS_PER_BP;
	}
	double* dist = new double[max];
	for (int i=0;i<max;i++) {
		dist[i] = 0;
	}

	int scaleFactor = 1;
	if (mCflag) {
		scaleFactor = 20;
	}

	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	if (singleFile) readSingleTagFile();
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*) chrs->search(keys[i]);
		ct->getTagCountDistribution(dist,max,scaleFactor);
		delete [](keys[i]);
	}
	delete []keys;


	float totalDist = 0;
	averageTagsPerPosition = 0;
	int maxValue = 0;
	for (int i=0;i<max;i++) {
		dist[i] /= (double)totalPositions;
		if (dist[i] > 0) maxValue = i;
		averageTagsPerPosition += ((double)i)/((double)scaleFactor)*dist[i];
		if (totalDist < 0.5 && totalDist+dist[i] > 0.5) {
			medianTagsPerPosition = i;
		}
		totalDist += dist[i];
	}

	if (mCflag) {
		double median = ((double)medianTagsPerPosition)/((double)scaleFactor);
		fprintf(stderr, "\tMedian mC/C = %.2lf\n",median);
		fprintf(stderr, "\tAverage mC/C = %.3lf\n", averageTagsPerPosition);
		fprintf(fp, "Tags per tag position (Median = %.2lf, tags per genomic bp = %.3lf)", 
					median,tbp);
		fprintf(fp, "\tFraction of Positions\n");
		for (int i=0;i<max;i++) {
			fprintf(fp,"%.2lf\t%lf\n",((double)i)/((double)scaleFactor),dist[i]);
			if (i>=maxValue) break;
		}
	} else {

		int idealTbp = (int)ceil(tbp*2+0.001);

		fprintf(stderr, "\tMedian tags per position = %d (ideal: %d)\n",medianTagsPerPosition,idealTbp);
		fprintf(stderr, "\tAverage tags per position = %.3lf\n", averageTagsPerPosition);
		if (medianTagsPerPosition > idealTbp*3) {
			fprintf(stderr, "\t\t!! Might have some clonal amplification in this sample if sonication was used\n");
			fprintf(stderr, "\t\tIf this is ChIP-Seq using sonicated fragments, consider adding the option \"-tbp %d\"\n",
								idealTbp);
			fprintf(stderr, "\t\tIgnore if analyzing RNA, MNase, etc. data\n");
		} else if (medianTagsPerPosition > idealTbp) {
			fprintf(stderr, "\t\t!! Might have some clonal amplification in this sample if sonication was used\n");
			fprintf(stderr, "\t\tIf using BED files, try using -forceBED to ignore the 5th column\n");
			fprintf(stderr, "\t\tIgnore if analyzing RNA, MNase, etc. data\n");
		}
		fprintf(fp, "Tags per tag position (Median = %d, tags per genomic bp = %.3lf)", 
					medianTagsPerPosition,tbp);
		fprintf(fp, "\tFraction of Positions\n");
		for (int i=0;i<max;i++) {
			fprintf(fp,"%d\t%lf\n",i,dist[i]);
			if (i>=maxValue) break;
		}
	}

	if (nfp == NULL) {
		fclose(fp);
	}

	return dist;
}
double TagLibrary::getAdjustedTagTotal() {
	double total = 0;
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*) chrs->search(keys[i]);
		total += ct->totalTags;
		delete [](keys[i]);
	}
	delete []keys;
	return total;
}
void TagLibrary::makeName() {
	if (directory == NULL) {
		setName((char*)"Unknown");
		fprintf(stderr, "For some reason this tag library doesn't have a directory?!?\n");
		return;
	}
	char* str = new char[strlen(directory)+1];
	char** cols = new char*[1000];
	int numCols = 0;
	strcpy(str,directory);
	split(str,cols,numCols,'/');
	char* newname = cols[numCols-1];
	if (cols[numCols-1][0] == '\0' && numCols > 1) {
		newname= cols[numCols-2];
	}
	setName(newname);
	delete []str;
	delete []cols;

}
void TagLibrary::printBedGraph(FILE* fp, double normTotalTags, char strand, int resolution, int negFlag,
			double fileSize, char* color, char* bedname, int method, int lastTagFlag,
			char* uniqMapDirectory,int style,int condenseFlag, int circosFlag, Peak* circosPeak,
			TagLibrary* input, double pseudoCounts, int logFlag, double normLength) {
	char* cc = color;
	int c1 = 0;
	int c2 = 0;
	int c3 = 0;
	if (cc == NULL) {
		cc = new char[100];
		srand ( time(NULL) );
		c1 = (int)(256*(((double)rand())/RAND_MAX));
		c2 = (int)(256*(((double)rand())/RAND_MAX));
		c3 = (int)(256*(((double)rand())/RAND_MAX));
		sprintf(cc,"%d,%d,%d",c1,c2,c3);
	}
	char* ucscname= NULL;
	if (bedname == NULL) {
		if (name == NULL) {
			makeName();
		}
		char* rationame = new char[1000];
		rationame[0] = '\0';
		if (input != NULL) {
			if (logFlag) {
				sprintf(rationame, " normalized to %s (log2)", input->directory);
			} else {
				sprintf(rationame, " normalized to %s", input->directory);
			}
		}
		ucscname = new char[10000];
		if (style == UCSC_METHYLATED_CpG) {
			sprintf(ucscname,"%s Methylated Cytosines Total Positions = %.2e%s",name, (double)totalPositions,rationame);
		} else if (style == UCSC_UNMETHYLATED_CpG) {
			sprintf(ucscname,"%s Unmethylated Cytosines Total Positions = %.2e%s",name, (double)totalPositions,rationame);
		} else {
			sprintf(ucscname,"%s Total Tags = %.2e, normalized to %.2e%s",name, totalTags,normTotalTags,rationame);
		}
			
	} else {
		ucscname = bedname;
	}

	char* fwdname = new char[10000];
	char curStrand = strand;
	if (strand == STRAND_SEPARATE || strand == STRAND_POSITIVE) {
		sprintf(fwdname,"%s + strand",ucscname);
		curStrand = STRAND_POSITIVE;
	} else {
		strcpy(fwdname,ucscname);
	}


	double reductionRatio = 1.0;
	if (fileSize > 1.0) {
		double datapoints = fileSize / 2.0 / BYTES_PER_LINE;
		if (strand == STRAND_POSITIVE || strand == STRAND_NEGATIVE) {
			datapoints /= 2.0;
		}
		reductionRatio = datapoints / ((double)totalPositions);
		//fprintf(stderr, "%lf\n%lf\n%lf\n", datapoints, totalTags, reductionRatio);
		if (reductionRatio < 1) {
			fprintf(stderr, "\tReduction Ratio set at %lf\n",reductionRatio);
		} else {
			fprintf(stderr, "\tNo need to remove tags to get desired file size\n");
		}
	}
	double normFactor = normTotalTags/totalTags;
	if (normTotalTags < 0.00001) normFactor = 1.0;
	double inputNormFactor = 1.0;
	if (input != NULL) {
		if (input->totalTags > 0.001) inputNormFactor =normTotalTags/input->totalTags;
		if (normTotalTags < 0.00001) inputNormFactor = 1.0;
		if (normFactor > inputNormFactor) {
			pseudoCounts*=normFactor;
		} else {
			pseudoCounts*=inputNormFactor;
		}
	}
	//might need to change this...
	int inputFragLength = fragmentLengthEstimate;

	
	
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);


	if (strand == STRAND_SEPARATE || strand == STRAND_POSITIVE || strand==STRAND_BOTH) {
		if (method == UCSC_BEDGRAPH && !circosFlag) {
			fprintf(fp, "track type=bedGraph name=\"%s\" description=\"%s\" ",fwdname, fwdname);
			fprintf(fp, "color=%s ",cc);
			fprintf(fp, "visibility=full yLineOnOff=on autoScale=on yLineMark=\"0.0\" ");
			fprintf(fp, "alwaysZero=on graphType=bar maxHeightPixels=128:75:11 ");
			if (style == UCSC_METHYLATED_CpG || style == UCSC_UNMETHYLATED_CpG) {
				fprintf(fp, "windowingFunction=mean smoothingWindow=off\n");
			} else {
				fprintf(fp, "windowingFunction=maximum smoothingWindow=off\n");
			}
		} else if (method == UCSC_BIGWIG && !circosFlag) {
			if (fp != stdout) {
				fprintf(stdout, "track type=bigWig name=\"%s\" description=\"%s\" ",fwdname, fwdname);
				fprintf(stdout, "bigDataUrl=http://server/path/file(Needs to be filled in) ");
				fprintf(stdout, "color=%s ",cc);
				fprintf(stdout, "visibility=full yLineOnOff=on autoScale=on yLineMark=\"0.0\" ");
				fprintf(stdout, "alwaysZero=on graphType=bar maxHeightPixels=128:75:11 ");
				if (style == UCSC_METHYLATED_CpG || style == UCSC_UNMETHYLATED_CpG) {
					fprintf(stdout, "windowingFunction=mean smoothingWindow=off\n");
				} else {
					fprintf(stdout, "windowingFunction=maximum smoothingWindow=off\n");
				}
			}
		}
	

		if (singleFile) readSingleTagFile();
		if (input != NULL && input->singleFile) input->readSingleTagFile();
		for (int i=0;i<chrs->total;i++) {
			if (circosFlag && circosPeak != NULL) {
				if (strcmp(keys[i],circosPeak->chr)!=0) continue;
			}

			fprintf(stderr, "\tGenerating bedGraph for %s\n", keys[i]);
			ChrTags* ct = (ChrTags*)chrs->search(keys[i]);
			ChrTags* inputct = NULL;
			if (input != NULL) {
				inputct = (ChrTags*)input->chrs->search(keys[i]);
			}

			UniqMapChrs* umc = NULL;
			if (uniqMapDirectory != NULL) umc = new UniqMapChrs(keys[i],uniqMapDirectory,0);

			ct->printBedGraph(fp,normFactor,curStrand,resolution,negFlag,fragmentLengthEstimate,
										reductionRatio,style,condenseFlag,lastTagFlag,umc,circosPeak,
										inputct, pseudoCounts, logFlag,inputNormFactor,inputFragLength,
										normLength);
			if (umc != NULL) delete umc;
			if (method == UCSC_BIGWIG && !circosFlag) {
				//long long end = 1;
				if (ct->totalPositions > 0) {
					if (fragmentLengthEstimate > 0) {
						//end = ct->appearentSize+fragmentLengthEstimate;
					} else {
						//end = ct->appearentSize+200;
					}
					//end = ct->appearentSize;
				}
			}
		}
	}


	if (strand == STRAND_SEPARATE || strand == STRAND_NEGATIVE) {
		curStrand = STRAND_NEGATIVE;
		sprintf(fwdname,"%s - strand",ucscname);
		if (color == NULL) {
			sprintf(cc,"%d,%d,%d",c1/2,c2/2,c3/2);
		}
		if (method == UCSC_BEDGRAPH && !circosFlag) {
			fprintf(fp, "track type=bedGraph name=\"%s\" description=\"%s\" ",fwdname, fwdname);
			fprintf(fp, "color=%s ",cc);
			fprintf(fp, "visibility=full yLineOnOff=on autoScale=on yLineMark=\"0.0\" ");
			fprintf(fp, "alwaysZero=on graphType=bar maxHeightPixels=128:75:11 ");
			if (style == UCSC_METHYLATED_CpG || style == UCSC_UNMETHYLATED_CpG) {
				fprintf(fp, "windowingFunction=mean smoothingWindow=off\n");
			} else {
				fprintf(fp, "windowingFunction=maximum smoothingWindow=off\n");
			}
		} else if (method == UCSC_BIGWIG && !circosFlag) {
			if (fp != stdout) {
				fprintf(stdout, "track type=bigWig name=\"%s\" description=\"%s\" ",fwdname, fwdname);
				fprintf(stdout, "bigDataUrl=http://server/path/file(Needs to be filled in) ");
				fprintf(stdout, "color=%s ",cc);
				fprintf(stdout, "visibility=full yLineOnOff=on autoScale=on yLineMark=\"0.0\" ");
				fprintf(stdout, "alwaysZero=on graphType=bar maxHeightPixels=128:75:11 ");
				if (style == UCSC_METHYLATED_CpG || style == UCSC_UNMETHYLATED_CpG) {
					fprintf(stdout, "windowingFunction=mean smoothingWindow=off\n");
				} else {
					fprintf(stdout, "windowingFunction=maximum smoothingWindow=off\n");
				}
			}
		}

		if (singleFile) readSingleTagFile();
		if (input != NULL && input->singleFile) input->readSingleTagFile();
		for (int i=0;i<chrs->total;i++) {
			if (circosFlag && circosPeak != NULL) {
				if (strcmp(keys[i],circosPeak->chr)!=0) continue;
			}

			fprintf(stderr, "\tGenerating reverse strand bedGraph for %s\n", keys[i]);
			ChrTags* ct = (ChrTags*)chrs->search(keys[i]);
			ChrTags* inputct = NULL;
			if (input != NULL) {
				inputct = (ChrTags*)input->chrs->search(keys[i]);
			}

			UniqMapChrs* umc = NULL;
			if (uniqMapDirectory != NULL) umc = new UniqMapChrs(keys[i],uniqMapDirectory,0);

			ct->printBedGraph(fp,normFactor,curStrand,resolution,negFlag,fragmentLengthEstimate,
										reductionRatio,style,condenseFlag,lastTagFlag,umc,circosPeak,
										inputct, pseudoCounts, logFlag,inputNormFactor, inputFragLength,
										normLength);
			if (umc != NULL) delete umc;

			if (method == UCSC_BIGWIG) {
				//long long end = 1;
				if (ct->totalPositions > 0) {
					if (fragmentLengthEstimate > 0) {
						//end = ct->appearentSize+fragmentLengthEstimate;
					} else {
						//end = ct->appearentSize+200;
					}
				}
			}
		}
	}
	

	if (color == NULL) {
		delete []cc;
	}
	if (bedname == NULL) {
		delete []ucscname;
	}
	delete []fwdname;
	for (int i=0;i<chrs->total;i++) {
		delete [](keys[i]);
	}
	delete []keys;
}


void TagLibrary::annotateTagLocations(PeakLibrary* annotations, FILE* statsFile, FILE* annFile) {

	Doubletable* stats = new Doubletable();

	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	if (singleFile) readSingleTagFile();
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*) chrs->search(keys[i]);
		ChrPeaks* cpAnn = (ChrPeaks*) annotations->chrs->search(keys[i]);
		if (cpAnn == NULL) {
			//fprintf(stderr, "\tNo annotations for on chromosome %s\n", keys[i]);
			continue;
		}
		fprintf(stderr, "\tAnnotating %s\n", keys[i]);
		ct->annotateTagLocations(cpAnn, annFile, stats);
		delete [](keys[i]);
	}
	delete []keys;

	char intron[20] = "Intron";
	char intergenic[20] = "Intergenic";
	char exon[20] = "Exon";
	char promoter[20] = "Promoter";

	keys = stats->keys();
	for (int i=0;i<stats->total;i++) {
		char* annname = keys[i];
		if (strcmp(keys[i],"I")==0) {
			annname = intron;
		} else if (strcmp(keys[i],"N")==0) {
			annname = intergenic;
		} else if (strcmp(keys[i],"P")==0) {
			annname = promoter;
		} else if (strcmp(keys[i],"E")==0) {
			annname = exon;
		}

		if (statsFile != NULL) {
			fprintf(statsFile,"%s\t%.1lf\n",annname, stats->search(keys[i]));
		}
		fprintf(stderr,"\t\t%s\t%.1lf\n",annname, stats->search(keys[i]));
		delete [](keys[i]);
	}
	delete []keys;
	delete stats;

}

Doubletable* TagLibrary::getPeakTagCounts(PeakLibrary* p, char strand) {
	Doubletable* counts = new Doubletable(MAX_PEAKS_AT_A_TIME);
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	if (singleFile) readSingleTagFile();
	//fprintf(stderr, "\tCounting:\n");
	for (int i=0;i<chrs->total;i++) {
		//fprintf(stderr, "\t\t%s\n",keys[i]);
		ChrTags* ct = (ChrTags*) chrs->search(keys[i]);
		ChrPeaks* cp = (ChrPeaks*) p->chrs->search(keys[i]);
		if (ct != NULL && cp != NULL) {
			ct->getPeakTagCounts(counts,cp,strand);
		}
		delete [](keys[i]);
	}
	delete []keys;
	//fprintf(stderr, "\n");
	return counts;
}
void TagLibrary::readAndSave() {
	char** keys = chrs->keys();
	totalTags = 0.0;
	totalPositions = 0;
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	if (singleFile) readSingleTagFile();
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*)chrs->search(keys[i]);
		ct->readAndSave();
		totalTags += ct->totalTags;
		totalPositions += ct->totalPositions;
		delete [](keys[i]);
	}
	delete []keys;
}
void TagLibrary::normalizeTagCountsGC(char* gcCtrlFile, NucleotideFreq* nf,
					double minNormRatio, double maxNormRatio, double gcWindow, double maxPercentError) {
	//checkTagSeqBias must have been performed before

	NucleotideFreq* ctrl = new NucleotideFreq();
	ctrl->readGCcontentFile(gcCtrlFile);

	nf->gcNorm = new double[nf->gcDistTotal];
	nf->gcOG = new double[nf->gcDistTotal];
	for (int i=0;i<nf->gcDistTotal;i++) {
		nf->gcNorm[i]=1.0;
		nf->gcOG[i]=nf->gcDist[i];
	}
	NucleotideFreq* nfcopy = nf->copy();
	int maxIterations = 50;
	char* normFileName = new char[10000];
	double lastError = FLT_MAX;
	for (int i=0;i<maxIterations;i++) {
		//sprintf(normFileName,"%s/tagGCnormalization.step%d.txt",directory,i+1);
		double perror = nfcopy->calculateGCNormalization(ctrl,minNormRatio, maxNormRatio,gcWindow,NULL);
		if (perror > lastError) {
			fprintf(stderr, "\tError increasing... stopping optimization.\n");
			break;
		}
		lastError = perror;
		fprintf(stderr, "\t\tNormalizing GC (Unaccounted for reads = %.2lf%%)\n",perror*100.0);
		for (int i=0;i<nf->gcDistTotal;i++) {
			nf->gcNorm[i] *= nfcopy->gcNorm[i];
			nfcopy->gcDist[i] *= nfcopy->gcNorm[i];
		}
		if (perror < maxPercentError) break;
	}
	sprintf(normFileName,"%s/tagGCnormalization.txt",directory);
	nf->printGCnormFile(normFileName);
	delete []normFileName;

	totalTags = 0.0;
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	if (singleFile) readSingleTagFile();
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*)chrs->search(keys[i]);
		ct->normalizeTagCountsGC(nf);
		delete [](keys[i]);
		totalTags += ct->totalTags;
	}
	delete []keys;

	gcNormalized=1;

}
void TagLibrary::normalizeTagCountsOligos(char* genomeDirectory, Hashtable* oligos,
				int oligoLength,int regionStart, int regionEnd,double minFold, double maxFold, float maxPerBp) {
	//checkTagSeqBias must have been performed before
	char* fname = new char[10000];

	fprintf(stderr, "\tChecking Oligo Enrichment Profiles...\n");
	int normalizeFlag = 0;
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	if (singleFile) readSingleTagFile();
	for (int i=0;i<chrs->total;i++) {
		fprintf(stderr, "\t\t%s\n", keys[i]);
		ChrTags* ct = (ChrTags*)chrs->search(keys[i]);
		ct->normalizeTagCountsOligos(genomeDirectory,oligos,oligoLength,regionStart,
										regionEnd,minFold,maxFold, maxPerBp, normalizeFlag);
		delete [](keys[i]);
	}
	delete []keys;

	char* rOligo = new char[oligoLength+1];
	keys =  oligos->keys();
	for (int i=0;i<oligos->total;i++) {
		OligoProfile* op = (OligoProfile*) oligos->search(keys[i]);
		if (op->revoppFlag) continue;
		strcpy(rOligo,keys[i]);
		revopp(rOligo);
		OligoProfile* rop = (OligoProfile*) oligos->search(rOligo);
		if (rop == NULL) {
			op->normalize();
		} else {
			OligoProfile::mergeRevopps(op,rop);
			op->normalize();
			rop->normalize();
//sprintf(fname,"%s/%s.profile.txt",directory,keys[i]);
//fp = fopen(fname, "w");
//op->print(fp);
//fclose(fp);
		}
		delete [](keys[i]);
	}
	delete []keys;
	delete []rOligo;


	totalTags = 0.0;
	keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	normalizeFlag = 1;
	fprintf(stderr, "\tNormalizing signal based on Oligo enrichment...\n");
	if (singleFile) readSingleTagFile();
	for (int i=0;i<chrs->total;i++) {
		fprintf(stderr, "\t\t%s\n", keys[i]);
		ChrTags* ct = (ChrTags*)chrs->search(keys[i]);
		ct->normalizeTagCountsOligos(genomeDirectory,oligos,oligoLength,regionStart,
										regionEnd,minFold,maxFold, maxPerBp,normalizeFlag);
		delete [](keys[i]);
		totalTags += ct->totalTags;
	}
	delete []keys;
	delete []fname;
	oligoNormalized = 1;
	
}

void TagLibrary::normalizeTagCountsFixedOligo(char* genomeDirectory, OligoArray* oligos,
									int oligoStart,int oligoEnd,double minFold, double maxFold) {
	//checkTagSeqBias must have been performed before

	totalTags = 0.0;
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	if (singleFile) readSingleTagFile();
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*)chrs->search(keys[i]);
		ct->normalizeTagCountsFixedOligo(genomeDirectory,oligos,oligoStart,oligoEnd,minFold,maxFold);
		delete [](keys[i]);
		totalTags += ct->totalTags;
	}
	if (ChrTags::allSeqs != NULL) {
		delete ChrTags::allSeqs;
		ChrTags::allSeqs = NULL;
	}
	delete []keys;
	oligoNormalized = 1;
}

void TagLibrary::removeTagSpikes(int spikeSize, double spikeFold) {

	fprintf(stderr, "\tRemoving tags spikes (%d bp regions with >%.1lfx tags than average)...\n",
							spikeSize,spikeFold);

	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	int oldTotalPositions = totalPositions;
	totalPositions = 0;
	double oldTotalTags = totalTags;
	totalTags = 0.0;

	PeakLibrary* regions = new PeakLibrary(1000000);

	double avg = 0.0;
	int N = 0;
	if (singleFile) readSingleTagFile();
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*)chrs->search(keys[i]);
		double chrTotal = 0.0;
		int chrN = 0;
		ct->getAverageCoverage(spikeSize, chrTotal, chrN, regions);
		avg += chrTotal;
		N += chrN;
	}
	regions->sortChr();
	regions->setDefaultPeakOrder();
	if (N > 0) {
		avg /= (double)N;
	}
	fprintf(stderr, "\t\tAverage tags per %d bp is %.1lf\n", spikeSize,avg);

	PeakLibrary* badRegions = new PeakLibrary(1000000);
	for (int i=0;i<regions->numPeaks;i++) {
		Peak* p = (Peak*)regions->peakOrder[i];
		if (p->v > avg*spikeFold) {
			badRegions->addPeak(p);
		}
	}
	badRegions->sortChr();
	badRegions->setDefaultPeakOrder();
	fprintf(stderr, "\t\t%d of %d regions exceeded tag coverage limit (%.1lf%%)\n", 
					badRegions->numPeaks, regions->numPeaks, 100.0*badRegions->numPeaks/regions->numPeaks);
	delete regions;

	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*)chrs->search(keys[i]);
		ct->removeTagsInPeaks(badRegions);
		totalPositions += ct->totalPositions;
		totalTags += ct->totalTags;
	}
	fprintf(stderr, "\t\tNew Total Positions: %lld (%.2lf%% kept)\n", totalPositions, 
												((double)totalPositions)/((double)oldTotalPositions)*100.0);
	fprintf(stderr, "\t\tNew Total Tags: %.1lf (%.2lf%% kept)\n", totalTags, totalTags/oldTotalTags*100.0);
	delete badRegions;

}
void TagLibrary::removePETagBackground(int fragLength) {

	fprintf(stderr, "\tRemoving PE tags separated by less than %d bp...\n", fragLength);

	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	int oldTotalPositions = totalPositions;;
	totalPositions = 0;
	double oldTotalTags = totalTags;
	totalTags = 0.0;
	if (singleFile) readSingleTagFile();
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*)chrs->search(keys[i]);
		ct->removePETagBackground(fragLength);
		totalPositions += ct->totalPositions;
		totalTags += ct->totalTags;
	}
	fprintf(stderr, "\tNew Total Positions: %lld (%.2lf%% kept)\n", totalPositions, 
												((double)totalPositions)/((double)oldTotalPositions)*100.0);
	fprintf(stderr, "\tNew Total Tags: %.1lf (%.2lf%% kept)\n", totalTags, totalTags/oldTotalTags*100.0);

}
void TagLibrary::assignPETagsToRestrictionSites(char* site,int maxMisMatches,char* genomeDirectory,
							int mode,int midFlag,int selfLigationFlag,int removeRestrictionEnds, int fragLength) {

	setRestrictionSite(site);

	int distLength = fragmentLengthEstimate*RESTRICTION_SITE_DISTRIBUTION_SCALE;
	double* posStrand = new double[distLength];
	double* negStrand = new double[distLength];
	for (int i=0;i<distLength;i++) {
		posStrand[i] =0.0;
		negStrand[i] =0.0;
	}

	//find sites in the genome first
	fprintf(stderr, "\tDetermining %s locations genome wide...\n", site);
	Hashtable* sites = new Hashtable(10000000);
	int totalSites = 0;
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	if (singleFile) readSingleTagFile();
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*)chrs->search(keys[i]);
		ct->findRestrictionSites(site,maxMisMatches,sites,genomeDirectory);
		int* s = (int*)sites->search(keys[i]);
		totalSites+=s[0];
		fprintf(stderr, "\t\t%s = %d\n", keys[i],s[0]);
	}
	fprintf(stderr, "\t\tTotal Genome-wide: %d\n",totalSites);

	//assign tags to sites.
	int oldTotalPositions = totalPositions;;
	totalPositions = 0;
	double oldTotalTags = totalTags;
	totalTags = 0.0;
	if (singleFile) readSingleTagFile();
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*)chrs->search(keys[i]);
		ct->assignPETagsToRestrictionSites(site,sites,fragLength,mode,midFlag,selfLigationFlag,
							removeRestrictionEnds,posStrand, negStrand, distLength);
		totalPositions += ct->totalPositions;
		totalTags += ct->totalTags;
	}
	fprintf(stderr, "\tNew Total Positions: %lld (%.2lf%% of original)\n", totalPositions, 
												((double)totalPositions)/((double)oldTotalPositions)*100.0);
	fprintf(stderr, "\tNew Total Tags: %.1lf (%.2lf%% of original)\n", totalTags, totalTags/oldTotalTags*100.0);
	

	for (int i=0;i<chrs->total;i++) {
		int* sitePos = (int*) sites->search(keys[i]);
		if (sitePos != NULL) delete []sitePos;
		delete [](keys[i]);
	}
	delete []keys;

	char* filename = new char[100000];
	sprintf(filename,"%s/petagRestrictionDistribution.%s.mis%d.txt",directory,site,maxMisMatches);
	FILE* fp = fopen(filename, "w");
	if (fp != NULL) {
		fprintf(fp, "Distance from %s site\tReads + strand\tReads - strand\n", site);
		for (int i=0;i<distLength;i++) {
			int pos = i-distLength/2;
			fprintf(fp, "%d\t%lf\t%lf\n",pos,posStrand[i],negStrand[i]);
		}
		fclose(fp);
	} else {
		fprintf(stderr, "!!! couldn't open %s for writing...\n", filename);
	}
	delete []negStrand;
	delete []posStrand;
	delete []filename;
	delete sites;

}

void TagLibrary::setRestrictionSite(char* site) {
	if (site == NULL) return;
	if (restrictionSite != NULL) delete []restrictionSite;
	int siteLength = strlen(site);
	restrictionSite = new char[siteLength+1];
	strcpy(restrictionSite,site);
}

NucleotideFreq* TagLibrary::checkTagSeqBias(char* genomeDirectory,int freqStart,int freqEnd,
											OligoArray* &oligos, int oligoStart, int oligoEnd) {

	fprintf(stderr, "\t\tCurrent Fragment length estimate: %d\n", fragmentLengthEstimate);
	if (fragmentLengthEstimate < 1) {
		fprintf(stderr, "\t\tSomething is wrong - this is too short setting to 41\n");
		fragmentLengthEstimate = 41;
	}
	int seqlen = freqEnd+fragmentLengthEstimate-freqStart;
	NucleotideFreq* nf = new NucleotideFreq(freqStart, seqlen, fragmentLengthEstimate);
	NucleotideFreq* nfctrl = new NucleotideFreq(freqStart, seqlen, fragmentLengthEstimate);
	NucleotideFreq* nfuniq = new NucleotideFreq(freqStart, seqlen, fragmentLengthEstimate);

	int numCPUs = 0;

	oligos = NULL;
	if (oligoStart <= oligoEnd) {
		fprintf(stderr, "\t\tChecking 5' Bias\n");
		oligos = new OligoArray();
		oligos->init(oligoEnd-oligoStart+1, 2000000000, DNA_ALPHA, MOTIF_STRAND_POS,numCPUs);
	}


	fprintf(stderr, "\tChecking Tag/Fragment sequence for bias...\n");
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	if (singleFile) readSingleTagFile();
	for (int i=0;i<chrs->total;i++) {
		fprintf(stderr, "\t\t%s\n", keys[i]);
		ChrTags* ct = (ChrTags*) chrs->search(keys[i]);
		ct->checkTagSeqBias(genomeDirectory,nf,nfuniq,nfctrl,freqStart,freqEnd,fragmentLengthEstimate,
								oligos,oligoStart, oligoEnd);
		delete [](keys[i]);
	}
	delete []keys;

	//clean up for genomes that need to load all at once
	if (ChrTags::allSeqs != NULL) {
		delete ChrTags::allSeqs;
		ChrTags::allSeqs = NULL;
	}

	char* freqFileName= new char[100000];
	//print nucleotide frequency
	sprintf(freqFileName,"%s/tagFreq.txt",directory);
	FILE* fp = fopen(freqFileName, "w");
	if (fp == NULL) {
		fprintf(stderr, "!!!Could not open %s for writing frequencies!!!\n",freqFileName);
		return NULL;
	}
	nf->print(fp);
	fclose(fp);	

	sprintf(freqFileName,"%s/tagFreqUniq.txt",directory);
	fp = fopen(freqFileName, "w");
	if (fp == NULL) {
		fprintf(stderr, "!!!Could not open %s for writing frequencies!!!\n",freqFileName);
		return NULL;
	}
	nfuniq->print(fp);
	fclose(fp);	

	//print GCcontent
	sprintf(freqFileName,"%s/tagGCcontent.txt",directory);
	fp = fopen(freqFileName, "w");
	if (fp == NULL) {
		fprintf(stderr, "!!!Could not open %s for writing GC content!!!\n",freqFileName);
		return NULL;
	}
	double avgtagGC = nf->printGC(fp);
	fclose(fp);

	sprintf(freqFileName,"%s/genomeGCcontent.txt",directory);
	fp = fopen(freqFileName, "w");
	if (fp == NULL) {
		fprintf(stderr, "!!!Could not open %s for writing GC content!!!\n",freqFileName);
		return NULL;
	}
	double avgGenomeGC = nfctrl->printGC(fp);
	fclose(fp);

	fprintf(stderr, "\tAvg Fragment GC%% = %.2lf%%\n", avgtagGC*100.0);
	fprintf(stderr, "\tAvg Expected GC%% = %.2lf%%\n", avgGenomeGC*100.0);
	gcAverage = avgtagGC;


	if (oligos != NULL) {	
		sprintf(freqFileName,"%s/oligoEnrichment%dto%d.txt",directory,oligoStart,oligoEnd);
		fp = fopen(freqFileName, "w");
		if (fp == NULL) {
			fprintf(stderr, "!!!Could not open %s for writing oligo enrichment!!!\n",freqFileName);
			return NULL;
		}
		fprintf(fp,"Oligo\tRead Count(%d to %d)\tGenome Count\tEnrichment\n",oligoStart, oligoEnd);
		oligos->initializeActiveOligoArray();
		oligos->calculateFoldEnrichment(0.5);
		oligos->sortActiveOligos();
		for (int i=oligos->numActiveOligos-1;i>=0;i--) {
			fprintf(fp, "%s\t%.1f\t%1.f\t%.4lf\n", oligos->activeOligos[i]->seq,
						oligos->activeOligos[i]->numTarget, oligos->activeOligos[i]->numBackground,
						oligos->activeOligos[i]->value);
		}
	}


	delete nfuniq;
	delete nfctrl;
	delete []freqFileName;

	return nf;
}

Doubletable* TagLibrary::getChrSizesFromTags() {

	double genomeTotal = 0.0;
	Doubletable* chrSizes = new Doubletable(1000);
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		ChrTags* ct = (ChrTags*)chrs->search(keys[i]);
		double maxPosition = (double)ct->getMaxPosition();
		genomeTotal += maxPosition;
		chrSizes->insert(maxPosition,keys[i]);
		delete [](keys[i]);
	}
	fprintf(stderr, "\tGenome Size=%.1lf\n",genomeTotal);
	chrSizes->insert(genomeTotal,(char*)"genome");
	
	return chrSizes;

}


//NOT FINISHED
PeakLibrary* TagLibrary::getCoverageRestrictionFragments(char* chr, int start, int end, int superRes,
													char* siteSeq, int maxMisMatches,char* genomeDirectory) {
/*
	setRestrictionSite(siteSeq);

	PeakLibrary* peaks = new PeakLibrary();

	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);

	if (superRes < 0) superRes = resolution;
	int halfSuperRes = superRes /2;
	for (int i=0;i<chrs->total;i++) {
		if (chr != NULL && strcmp(chr,keys[i])!=0) {
			delete [](keys[i]);
			continue;
		}
		ChrTags* ct = (ChrTags*)chrs->search(keys[i]);
		int maxPosition = ct->getMaxPosition();
		//fprintf(stderr, "%s\t%d\n", keys[i],maxPosition);
		int cstart = start;
		int cend = end;
		if (cstart < 0) cstart = 0;
		if (cend > maxPosition) cend = maxPosition;


		Hashtable* sites = new Hashtable(1000000);
		ct->findRestrictionSites(site,maxMisMatches,sites,genomeDirectory);
		int* sitePos = (int*)sites->search(keys[i]);
		int totalSites = sitePos[0];



		for (int j=cstart;j<cend;j+=resolution) {
			char* pname = new char[50];
			sprintf(pname,"%s-%d",keys[i],j); 
			int midpoint = j+resolution/2;
			int s = midpoint - halfSuperRes;
			int e = midpoint + halfSuperRes;
			peaks->addPeak(pname,keys[i],s,e,midpoint,STRAND_POSITIVE,
							0.0, 0.0, NULL, 0,0);
			delete []pname;
		}

		delete [](keys[i]);
	}
	delete []keys;
	peaks->sortChr();
	return peaks;
*/
	return NULL;
}



PeakLibrary* TagLibrary::getCoveragePeaks(char* chr, int start, int end, int resolution,int superRes) {
	PeakLibrary* peaks = new PeakLibrary();

	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);

	if (superRes < 0) superRes = resolution;
	int halfSuperRes = superRes /2;
	for (int i=0;i<chrs->total;i++) {
		if (chr != NULL && strcmp(chr,keys[i])!=0) {
			delete [](keys[i]);
			continue;
		}
		ChrTags* ct = (ChrTags*)chrs->search(keys[i]);
		int maxPosition = ct->getMaxPosition();
		//fprintf(stderr, "%s\t%d\n", keys[i],maxPosition);
		int cstart = start;
		int cend = end;
		if (cstart < 0) cstart = 0;
		if (cend > maxPosition) cend = maxPosition;

		for (int j=cstart;j<cend;j+=resolution) {
			char* pname = new char[50];
			sprintf(pname,"%s-%d",keys[i],j); 
			int midpoint = j+resolution/2;
			int s = midpoint - halfSuperRes;
			int e = midpoint + halfSuperRes;
			peaks->addPeak(pname,keys[i],s,e,midpoint,STRAND_POSITIVE,
							0.0, 0.0, NULL, 0,0);
			delete []pname;
		}

		delete [](keys[i]);
	}
	delete []keys;
	peaks->sortChr();
	return peaks;
}


double TagLibrary::estimateContamination(TagLibrary* input, int maxDistance, float minThreshold) {

	setTagAdjust(0);
	input->setTagAdjust(0);

	PeakFinder* pf = new PeakFinder();
	pf->setDirectory(directory);
	pf->regionFlag = 1;
	pf->setPeakSize(maxDistance);
	pf->tagThresh = minThreshold;
	pf->filterMode = PEAKFINDER_FILTER_MODE_THRESH;
	pf->strand = STRAND_SEPARATE;
	pf->clonalFold = 0.0;
	pf->setMaxTBP(0.0,0.0);
	pf->setTagLibraries(input,NULL);
	PeakLibrary* peaks = pf->findPeaks();

	peaks->setPeakTagSizeFixed(0,0);
	//peaks->addTagLibrary(this);
	//peaks->addTagLibrary(input);

	Doubletable* targetCounts = peaks->countPeakTagsLowMemory(this,POSITIVE_STRAND,COUNT_MODE_TOTAL);
	Doubletable* inputCounts = peaks->countPeakTagsLowMemory(input,POSITIVE_STRAND,COUNT_MODE_TOTAL);
	//Doubletable* targetCounts = peaks->countPeakTags(0,0,0,POSITIVE_STRAND,COUNT_MODE_TOTAL);
	//Doubletable* inputCounts = peaks->countPeakTags(1,0,0,POSITIVE_STRAND,COUNT_MODE_TOTAL);

	char* scatterPlotName = new char[10000];
	sprintf(scatterPlotName, "%s/contaminationScatterPlot.txt", directory);
	FILE* fp = fopen(scatterPlotName, "w");
	fprintf(fp, "RegionID\tchr\tstart\tend\tstrand\tContaminated (log2)\tContaminant (log2)\n");

	double** values = new double*[2];	
	values[0] = new double[targetCounts->total];
	values[1] = new double[targetCounts->total];
	int numValues = 0;

	double histInc = 0.1;
	double histMax = 20.0;
	int histLimit = (int)(histMax/histInc)*2+1;
	double* hist = new double[histLimit];	
	double* shist = new double[histLimit];	
	int halfPoint = (int)(histMax/histInc);
	for (int i=0;i<histLimit;i++) {
		hist[i] = 0.0;
		shist[i] = 0.0;
	}

	char** keys = targetCounts->keys();
	for (int i=0;i<targetCounts->total;i++) {
		double tc = targetCounts->search(keys[i]);
		double ic = inputCounts->search(keys[i]);
		Peak* p = (Peak*)peaks->peaks->search(keys[i]);
		if (tc < EMPTY_DOUBLE_CHECK) continue;
		if (ic < EMPTY_DOUBLE_CHECK) continue;
		if (p == NULL) continue;

		double ntc = tc/totalTags*1e7;
		double nic = ic/input->totalTags*1e7;
		values[0][numValues] = ntc;
		values[1][numValues] = nic;
		numValues++;
		double logT = log(1.0+ntc)/log(2.0);
		double logI = log(1.0+nic)/log(2.0);
		if (tc > minThreshold/1.0) {
			double diff = logI-logT;
			int index = (int)(diff/histInc)+halfPoint;
			if (index < 0) index = 0;
			if (index >= histLimit) index = histLimit-1;
			hist[index]+=1.0;
		}

		if (p->strand == STRAND_POSITIVE) {
			fprintf(fp, "%s\t%s\t%d\t%d\t+\t%lf\t%lf\n", p->name,p->chr,p->start,p->end,logT,logI);
		} else {
			fprintf(fp, "%s\t%s\t%d\t%d\t-\t%lf\t%lf\n", p->name,p->chr,p->start,p->end,logT,logI);
		}
		delete [](keys[i]);
	}
	fclose(fp);


	sprintf(scatterPlotName, "%s/contaminationHistogram.txt", directory);
	fp = fopen(scatterPlotName, "w");
	fprintf(fp, "Log Difference (contaminant/contaminated)\tContamination Estimate\tCount(min %f)\tSmoothed\n",minThreshold);
	int good = 0;
	for (int i=0;i<histLimit;i++) {
		if (hist[i] < 0.01 && good == 0) continue;
		good = 1;
		double value = ((double)(i-halfPoint))*histInc;
		
		int n = 1;
		shist[i] = hist[i];
		if (i>0) {
			shist[i]+=hist[i-1];
			n++;
		}
		if (i+1<histLimit) {
			shist[i]+=hist[i+1];
			n++;
		}
		shist[i] /= (double)n;
	

		double f = 1.0/(pow(2.0,value));
		if (f > 1.0) f=1.0;
		fprintf(fp, "%.2lf\t%.3lf\t%.1lf\t%.1lf\n", value, f,hist[i],shist[i]);
	
	}
	fclose(fp);

	double est = 0.0;
	double lastValue = 0.0;
	double total=0.0;
	for (int i=histLimit-1;i>=0;i--) {
		total += shist[i];
		if (shist[i] < lastValue && total > 10.0) {
			double value = ((double)(i-halfPoint))*histInc;
			double f = 1.0/(pow(2.0,value));
			est = f;
			break;
		}
		lastValue = shist[i];
	}
	fprintf(stderr, "\tContamination estimate: %.3lf (%.1lf%%)\n", est, est*100.0);

	delete peaks;
	delete pf;
	delete []hist;
	delete []shist;
	delete []scatterPlotName;
	delete [](values[0]);
	delete [](values[1]);
	delete []values;
	delete []keys;

	return est;
}
void TagLibrary::decontaminate(TagLibrary* input, double fraction, int maxDistance) {


	double normFactor = totalTags/input->totalTags;

	fprintf(stderr, "\tDecontaminating...\n");

	double startingCount = totalTags;

	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		fprintf(stderr, "\t\t%s...\n",keys[i]);

		ChrTags* ct = (ChrTags*)chrs->search(keys[i]);
		ct->loadTags();
			
		ChrTags* bgct = (ChrTags*)input->chrs->search(keys[i]);
		if (bgct != NULL) {
			bgct->loadTags();

			//fprintf(stderr, "decontaminating...");
			ct->decontaminate(bgct,maxDistance,fraction,normFactor);

			bgct->freeTags();
		}
		ct->dontSAVE = 0;	
		ct->print();
		ct->freeTags();
		delete [](keys[i]);
	}
	optimizeTagFiles();
	printTagInfo();

	double removed = startingCount - totalTags;
	fprintf(stderr, "\tRemoved %.1lf tags (%.2lf%%)\n", removed, removed/startingCount*100.0);

}



// --------------- PETag Stuff -----------------------------------
double PeakLibrary::adjustPETagTotalsWithModel(HiCBgModel* model, int useTotalsTooFlag) {

	if (model->refPeaks == NULL) {
		fprintf(stderr, "\tUsing old background model - update with \"-force\" to continue\n");
		exit(0);
	}
	int resolution = model->res;
	double dresolution = (double)model->res;

	Peak** refOrder = model->refPeaks->peakOrder;
	int numRefPeaks = model->refPeaks->numPeaks;

	double accuracy = 0.0;
	double accN = 0.0;

	double avgFactor = 0.0;
	//double avgV = 0.0;
	double NNN = 0.0;

	int refIndex = 0;
	for (int i=0;i<numPeaks;i++) {
		//peakOrder[i]->print(stderr);
		while (refIndex < numRefPeaks &&
				(chrcmp(&(refOrder[refIndex]->chr),&(peakOrder[i]->chr)) < 0)) {
			refIndex++;
		}
		while (refIndex < numRefPeaks && 
					chrcmp(&(refOrder[refIndex]->chr),&(peakOrder[i]->chr)) == 0 &&
					refOrder[refIndex]->refPos+resolution < peakOrder[i]->refPos) {
			refIndex++;
		}
		if (refIndex >= numRefPeaks) break;
		if (chrcmp(&(refOrder[refIndex]->chr),&(peakOrder[i]->chr)) != 0) continue;

		double factor = 0.0;
		double totalV = 0.0;
		double w = 0.0;

		double max = 0;
		if (refOrder[refIndex]->refPos >= peakOrder[i]->refPos-resolution 
		 				&& refOrder[refIndex]->refPos <= peakOrder[i]->refPos+resolution 
						&& refOrder[refIndex]->v > 0.5) {

			double f = log(refOrder[refIndex]->focusRatio);
			int dist = refOrder[refIndex]->refPos - peakOrder[i]->refPos;
			double w1 = 1.0-fabs(((double)dist)/dresolution);
			if (w1 > max) max = w1;
			factor += f*w1;
			w += w1;
			totalV += refOrder[refIndex]->v*w1;
		}
		if (refIndex+1 < numRefPeaks
						&& chrcmp(&(refOrder[refIndex+1]->chr),&(peakOrder[i]->chr)) == 0
						&& refOrder[refIndex+1]->refPos >= peakOrder[i]->refPos-resolution 
		 				&& refOrder[refIndex+1]->refPos <= peakOrder[i]->refPos+resolution 
						&& refOrder[refIndex+1]->v > 0.5) {
			double f = log(refOrder[refIndex+1]->focusRatio);
			int dist = refOrder[refIndex+1]->refPos - peakOrder[i]->refPos;
			double w1 = 1.0-fabs(((double)dist)/dresolution);
			if (w1 > max) max = w1;
			factor += f*w1;
			w += w1;
			totalV += refOrder[refIndex+1]->v*w1;
		}
		if (w > 0.00001) {
			//fprintf(stderr, "%s\t%lf\t%d\t%lf\t%lf\n", peakOrder[i]->name, factor,peakOrder[i]->refPos,w,dresolution);
			factor /= w;
			factor = exp(factor);
			peakOrder[i]->focusRatio = factor;
			if (useTotalsTooFlag) {
				peakOrder[i]->v = totalV/w;
			}
			if (peakOrder[i]->v > 0) {
				avgFactor += factor/peakOrder[i]->v;
				NNN+= 1.0;
			}
		} else {
			peakOrder[i]->focusRatio = -1.0;
			if (useTotalsTooFlag) {
				//peakOrder[i]->v = 0.0;
			}
		}
		if (max > 0.05) {
			accuracy += max;
			accN += 1.0;
		}
	}
	if (NNN > 0) {
		avgFactor /= NNN;
	}
	for (int i=0;i<numPeaks;i++) {
		if (peakOrder[i]->focusRatio < 0.0) {
			peakOrder[i]->focusRatio = peakOrder[i]->v*avgFactor;
		}
	}
	//fprintf(stderr, "numPeaks = %d\trefPeaks = %d\n", numPeaks,numRefPeaks);
	accuracy /= (double) accN;

	return accuracy;
}


// need to make parallel xxxxxxxxxxxxxxxxxxxxxxx
double PeakLibrary::adjustPETagTotalsForModel(HiCBgModel* model, int maxCPUs) {

	PeakLibrary* refPeaks = model->refPeaks;
	if (model->refPeaks == NULL) {
		fprintf(stderr, "\tUsing old background model - update with \"-force\" to continue\n");
		exit(0);
	}
	int resolution = model->res;

	Peak** refOrder = refPeaks->peakOrder;
	int numRefPeaks = refPeaks->numPeaks;
	double* tmpTotals = new double[numPeaks];
	double* errors = new double[numPeaks];
	for (int i=0;i<numPeaks;i++) {
		//fprintf(stderr, "%d\t%f\t%f\n", i, peakOrder[i]->v,peakOrder[i]->focusRatio);
		peakOrder[i]->focusRatio = 1.0;
		tmpTotals[i] = 1.0;
		errors[0] = 0.0;
	}
	double totalReads = 0.0;
	for (int i=0;i<numRefPeaks;i++) totalReads += refOrder[i]->v*refOrder[i]->focusRatio;
	model->totalModelReads = totalReads;

	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	int maxIterations = 1000;
	double totalChange = 0.0;
	double *changeHistory = new double[maxIterations];
	for (int i=0;i<maxIterations;i++) changeHistory[i] = 0.0;
	int numNoImprove = 0;
	srand ( time(NULL) );
	int numRestarts = 0;
	int largeWindow = 40;
	double lastChange = FLT_MAX;
	for (int j=0;j<maxIterations;j++) {
		totalChange = 0.0;

		for (int i=0;i<chrs->total;i++) {
			ChrPeaks* cp = (ChrPeaks*) chrs->search(keys[i]);
			double dd = cp->adjustTotalsBasedOnModel(refPeaks,resolution,model,tmpTotals,errors,maxCPUs);
			totalChange += dd;
			//fprintf(stdout, "%d\t%d\t%s\t%lf\t%lf\n",j,i,keys[i],dd,totalChange);
		}

		changeHistory[j] = totalChange;
		fprintf(stderr, "\tIteration %d, delta=%.1lf\n", j+1, totalChange);
		if (j>largeWindow) {
			if (totalChange > changeHistory[j-largeWindow]/4.0) {
				fprintf(stderr, "\tProbably not going to converge... hopefully good enough!!\n");
				break;
			}
		}


		if (totalChange > lastChange-0.5) {
			numNoImprove++;
		} else {
			numNoImprove = 0;
		}
		for (int i=0;i<numPeaks;i++) {
			//focus ratio used to hold new "total" value
			peakOrder[i]->focusRatio = tmpTotals[i];
		}
		double totalReads = 0.0;
		for (int i=0;i<numRefPeaks;i++) totalReads += refOrder[i]->v*refOrder[i]->focusRatio;
		model->totalModelReads = totalReads;

		if (totalChange < 1.0 || numRestarts > 3) break;


		if (numNoImprove >= 3) {
			break;
			fprintf(stderr, "\t\tRandomizing values, contining...\n");
			numRestarts++;
			double totalReads = 0.0;
			for (int i=0;i<numRefPeaks;i++) {
				refOrder[i]->focusRatio *= 0.9+0.20*((double)rand()/((double)RAND_MAX));
				totalReads += refOrder[i]->v*refOrder[i]->focusRatio;
			}
			model->totalModelReads = totalReads;
		}

		lastChange = totalChange;
	}
	model->modelError = totalChange;

	int badCount = 0;
	for (int i=0;i<numPeaks;i++) {
		if (errors[i] > 1.0) {
			if (badCount == 0) {
				fprintf(stderr, "\tPeaks with scaling problems:\n");
			}
			fprintf(stderr, "\t\t%s\t%d\tError:%lf\n",refOrder[i]->chr,refOrder[i]->start, errors[i]);
			badCount++;
		}
	}
	for (int i=0;i<chrs->total;i++) {
		delete [](keys[i]);
	}
	delete []keys;
	delete []tmpTotals;
	delete []errors;
	return totalChange;
}

double ChrPeaks::adjustTotalsBasedOnModel(PeakLibrary* refPeaks, int resolution, HiCBgModel* model,
									double* tmpTotals,double* errors, int numCPUs) {
	double totalChange = 0.0;
	if (numPeaks < 1) return totalChange;

	Peak** refOrder = refPeaks->peakOrder;
	int numRefPeaks = refPeaks->numPeaks;

	char* chr = peaks[0]->chr;
	HiCBgModelChr* chrModel = (HiCBgModelChr*) model->chrs->search(peaks[0]->chr);
	if (chrModel == NULL) {
		chrModel = model;
	}

	double totalModelReads = model->totalModelReads;
	if (totalModelReads < 1.0) totalModelReads = 1.0;

	//add up interchromosomal expected interacitons first - same for all peaks
	double expectInter = 0.0;
	int minIndex = INT_MAX;
	int maxIndex = 0;
	for (int i=0;i<numRefPeaks;i++) {
		if (strcmp(refOrder[i]->chr,chr) != 0) {
			double expected = refOrder[i]->v*refOrder[i]->focusRatio;
			//expected *= chrModel->interChrScaleFactor;
			HiCBgModelChr* otherChrModel = (HiCBgModelChr*) model->chrs->search(refOrder[i]->chr);
			if (otherChrModel == NULL) otherChrModel = chrModel;

			expected *= sqrt(chrModel->interChrScaleFactor*otherChrModel->interChrScaleFactor);
			expectInter += expected;
		} else {
			if (i < minIndex) minIndex = i;
			if (i > maxIndex) maxIndex = i+1;
		}
	}
	expectInter /= totalModelReads;


	//now do intra chromosomal interactions
	//for (int i=0;i<numPeaks;i++) {
	//
	pthread_mutex_init(&mutex, NULL);
	pthread_mutex_init(&mutex2, NULL);
	mutexIndex = 0;
	mutexTotal = numPeaks;

	ThreadArgs_adjustTotals **args = new ThreadArgs_adjustTotals*[numCPUs];
	
	pthread_t* thread = new pthread_t[numCPUs];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	void* status;
	for (int i=0;i<numCPUs;i++) {
		args[i] = new ThreadArgs_adjustTotals(i,this,errors,tmpTotals,expectInter,minIndex,maxIndex,
						totalModelReads,resolution, chrModel,refPeaks);
		(void)pthread_create(&(thread[i]),&attr,&ChrPeaks_adjustTotalsThread,args[i]);
	}
	pthread_attr_destroy(&attr);
	for (int i=0;i<numCPUs;i++) {
		(void)pthread_join(thread[i],&status);
		delete args[i];
	}
	delete []thread;
	delete []args;
	pthread_mutex_destroy(&mutex);
	pthread_mutex_destroy(&mutex2);

	for (int i=0;i<numPeaks;i++) {
		totalChange += errors[peaks[i]->index];
	}
	return totalChange;

}
ThreadArgs_adjustTotals::ThreadArgs_adjustTotals(int ncpu, ChrPeaks* ccpeaks, double *e, double* t,
					double ei, int mi, int maxi, double tmr, int res, HiCBgModelChr* cm,
					PeakLibrary *rp) {
	cpu = ncpu;
	cpeaks = ccpeaks;
	errors = e;
	tmpTotals = t;
	expectInter = ei;
	minIndex = mi;
	maxIndex = maxi;
	totalModelReads = tmr;
	resolution = res;
	chrModel = cm;
	refPeaks = rp;
}
void* ChrPeaks_adjustTotalsThread(void* threadArgs) {
	ThreadArgs_adjustTotals* args = (ThreadArgs_adjustTotals*)threadArgs;
	(args->cpeaks)->adjustTotalsThread(args->errors, args->tmpTotals, args->expectInter,args->minIndex,
									args->maxIndex, args->totalModelReads, args->resolution,
									args->chrModel,args->refPeaks);
	return threadArgs;
}

void ChrPeaks::adjustTotalsThread(double* errors, double* tmpTotals, double expectInter,
					int minIndex,int maxIndex, double totalModelReads, 
					int resolution,HiCBgModelChr* chrModel,PeakLibrary* refPeaks) {
	int i=-1;
	Peak** refOrder = refPeaks->peakOrder;
	while (i<mutexTotal) {
		pthread_mutex_lock(&mutex);
		i = mutexIndex;
		mutexIndex++;
		pthread_mutex_unlock(&mutex);
		if (i>=mutexTotal) break;


		double expTotal = expectInter;
		for (int j=minIndex;j<maxIndex;j++) {
			double expected = (refOrder[j]->v*refOrder[j]->focusRatio) / totalModelReads;
			int distance = abs(peaks[i]->refPos - refOrder[j]->refPos);
			int bin = (distance+resolution/2)/resolution;
			expected *= chrModel->scaleFactor[bin];
			expTotal += expected;
		}

		double adjTotal = peaks[i]->v*peaks[i]->focusRatio;
		expTotal *= adjTotal;

		double error = fabs(expTotal-peaks[i]->v);	
		errors[peaks[i]->index] = error;

		double ratio = 1.0;
		if (expTotal > 0) {
			ratio = peaks[i]->v/expTotal;
		}
		tmpTotals[peaks[i]->index] = peaks[i]->focusRatio * ratio;
	}
}



//
// Start of getPETagTotals parallel stuff
//
void TagLibrary::getPETagTotals(PeakLibrary* peaks) {

	char** keys = chrs->keys();
	if (singleFile) readSingleTagFile();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	fprintf(stderr, "\tCalculating PE Tag Coverage:");

	pthread_mutex_init(&mutex, NULL);
	mutexIndex = 0;
	mutexPeakIndex = 0;
	mutexTotal = chrs->total;
	mutexChrKeys = keys;
	ThreadArgs_getPETagTotals **args = new ThreadArgs_getPETagTotals*[numCPUs];

	pthread_t* thread = new pthread_t[numCPUs];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	void* status;
	for (int i=0;i<numCPUs;i++) {
		args[i] = new ThreadArgs_getPETagTotals(i,this,peaks);
		(void)pthread_create(&(thread[i]),&attr,&TagLibrary_getPETagTotalsThread,args[i]);
	}
	pthread_attr_destroy(&attr);
	for (int i=0;i<numCPUs;i++) {
		(void)pthread_join(thread[i],&status);
		delete args[i];
	}
	delete []thread;
	delete []args;
	pthread_mutex_destroy(&mutex);
	delete []keys;
	fprintf(stderr, "\n");
}
ThreadArgs_getPETagTotals::ThreadArgs_getPETagTotals(int ncpu, TagLibrary* ntags,
																PeakLibrary* npeaks) {
	cpu = ncpu;
	tags = ntags;
	peaks = npeaks;
}
void* TagLibrary_getPETagTotalsThread(void* threadArgs) {
	ThreadArgs_getPETagTotals* args = (ThreadArgs_getPETagTotals*)threadArgs;
	(args->tags)->getPETagTotalsThread(args->cpu, args->peaks);
	return threadArgs;
}

void TagLibrary::getPETagTotalsThread(int cpu, PeakLibrary* peaks) {
	int i=-1;
	int peakIndex = 0;
	while (i<mutexTotal) {
		pthread_mutex_lock(&mutex);
		if (mutexPeakIndex < peakIndex) mutexPeakIndex = peakIndex;
		i=mutexIndex++;
		peakIndex = mutexPeakIndex;
		pthread_mutex_unlock(&mutex);
		if (i>=mutexTotal) break;
		
		char* c = mutexChrKeys[i];
		fprintf(stderr, ".");
		ChrTags* ct = (ChrTags*) chrs->search(c);
		ct->getPETagTotals(peaks,peakIndex);
		delete []c;
	}

}


//
// Start of getPETagTotals parallel stuff
//
void TagLibrary::calcHiCMatrix(PeakLibrary* peaks1,PeakLibrary* peaks2,double** matrix,
					int resolution, HiCBgModel* model, int revFlag, int numHistBins, 
					int actionFlag, double totalInteractions, GenomeInteractionLibrary* gil,
					HiCparams* params) {

	char** keys = chrs->keys();
	if (singleFile) readSingleTagFile();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);


	pthread_mutex_init(&mutex, NULL);
	pthread_mutex_init(&mutexMatrix, NULL);
	if (gil != NULL) pthread_mutex_init(&(gil->mutex), NULL);
	if (model != NULL) pthread_mutex_init(&(model->mutex), NULL);
	mutexIndex = 0;
	mutexPeakIndex = 0;
	mutexTotal = chrs->total;
	mutexChrKeys = keys;
	ThreadArgs_makeHiCMatrix **args = new ThreadArgs_makeHiCMatrix*[numCPUs];

	pthread_t* thread = new pthread_t[numCPUs];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	void* status;
	for (int i=0;i<numCPUs;i++) {
		args[i] = new ThreadArgs_makeHiCMatrix(i,this,peaks1,peaks2,matrix,resolution,model,
										revFlag,numHistBins,actionFlag,totalInteractions,gil,
										params);
		(void)pthread_create(&(thread[i]),&attr,&TagLibrary_makeHiCMatrixThread,args[i]);
	}
	pthread_attr_destroy(&attr);
	for (int i=0;i<numCPUs;i++) {
		(void)pthread_join(thread[i],&status);
		delete args[i];
	}
	delete []thread;
	delete []args;
	pthread_mutex_destroy(&mutex);
	pthread_mutex_destroy(&mutexMatrix);
	if (gil != NULL) pthread_mutex_destroy(&(gil->mutex));
	if (model != NULL) pthread_mutex_destroy(&(model->mutex));
	fprintf(stderr, "\n");
	delete []keys;
}
ThreadArgs_makeHiCMatrix::ThreadArgs_makeHiCMatrix(int ncpu, TagLibrary* ntags,
				PeakLibrary* npeaks1, PeakLibrary* npeaks2,
                double** nmatrix, int nresolution, HiCBgModel* nmodel, int nrevFlag,
                int nnumHistBins, int nactionFlag , double ntotalInteractions,
                GenomeInteractionLibrary* ngil,HiCparams* nparams) {
	cpu = ncpu;
	tags = ntags;
	peaks1 = npeaks1;
	peaks2 = npeaks2;
	matrix = nmatrix;
	resolution = nresolution;
	model = nmodel;
	revFlag = nrevFlag;
	numHistBins = nnumHistBins;
	actionFlag = nactionFlag;
	totalInteractions = ntotalInteractions;
	gil = ngil;
	params = nparams;
}
void* TagLibrary_makeHiCMatrixThread(void* threadArgs) {
	ThreadArgs_makeHiCMatrix* args = (ThreadArgs_makeHiCMatrix*)threadArgs;
	(args->tags)->makeHiCMatrixThread(args->cpu, args->peaks1,args->peaks2,
				args->matrix,args->resolution, args->model,args->revFlag,
				args->numHistBins, args->actionFlag, args->totalInteractions,
				args->gil,args->params);
	return threadArgs;
}

void TagLibrary::makeHiCMatrixThread(int cpu,
				PeakLibrary* peaks1, PeakLibrary* peaks2,
                double** matrix, int resolution, HiCBgModel* model, int revFlag,
                int numHistBins, int actionFlag , double totalInteractions,
                GenomeInteractionLibrary* gil,HiCparams* params) {
	int i=-1;
	int peakIndex = 0;
	while (i<mutexTotal) {
		pthread_mutex_lock(&mutex);
		if (mutexPeakIndex < peakIndex) mutexPeakIndex = peakIndex;
		i=mutexIndex++;
		peakIndex = mutexPeakIndex;
		pthread_mutex_unlock(&mutex);
		if (i>=mutexTotal) break;
		
		char* c = mutexChrKeys[i];
		ChrTags* ct = (ChrTags*) chrs->search(c);
		params->fragLengthEstimate = fragmentLengthEstimate;
		ct->makeHiCMatrix(peaks1,peakIndex,peaks2,matrix,resolution,model,
					revFlag,numHistBins,actionFlag,totalInteractions,gil,&mutexMatrix,
					params);

		delete []c;
	}
}


void TagLibrary::makeHiCBgModel(HiCBgModel* model,PeakLibrary* &peaks1, PeakLibrary* peaks2, 
								int peaks1IsGenomeFlag, int fullModelFlag, HiCparams* params) {

	int actionFlag = HIC_MASK_NORM_SEQDEPTH | HIC_MASK_CREATEMODEL;
	double* dist = NULL;
	int distLength = 0;


	if (fullModelFlag) {
		fprintf(stderr, "\tGenerating Background using -fullModel\n");
	} else {
		fprintf(stderr, "\tGenerating Background using -quickModel\n");
	}

	PeakLibrary* refPeaks = NULL;
	if (peaks1 != NULL && peaks1IsGenomeFlag) {
		refPeaks = peaks1;
	} else {
		refPeaks = getCoveragePeaks(NULL,-1,1000000000,model->res,model->res);
		if (peaks1 == NULL) {
			peaks1 = refPeaks;
			peaks1IsGenomeFlag=1;
		}
	}
	if (peaks1IsGenomeFlag==0) {
		fprintf(stderr, "\tGenerating background model from a user specified subset of the genome\n");
	}
	model->totalRegions = refPeaks->numPeaks;
	model->totalModelReads = totalTags;
	peaks1->setDefaultPeakOrder();

	if (peaks2 != NULL) {
		fprintf(stderr, "\tNon-symetrical background...\n");
		peaks2->setDefaultPeakOrder();
	} else {
		peaks2 = peaks1;
	}
	if (refPeaks != peaks1) {
		refPeaks->setDefaultPeakOrder();
	}

	double totalInteractions = totalTags;
	//fprintf(stderr, "actionflag = %d\n", actionFlag);
	//delete model->chrSizes;
	model->chrSizes = getChrSizesFromTags();
	model->initialize();

	getPETagTotals(peaks1);
	if (peaks1 != peaks2) {
		fprintf(stderr, "\tCalculating read coverage of 2nd set of peaks\n");
		getPETagTotals(peaks2);
	}
	if (peaks1IsGenomeFlag == 0) {
		fprintf(stderr, "\tCalculating read coverage of genomic reference (differs from input files)\n");
		getPETagTotals(refPeaks);
	}

	double avg = 0;
	double std = 0;
	double N = 0;
	peaks1->setDefaultPeakOrder();
	for (int i=0;i<peaks1->numPeaks;i++) {
		avg += peaks1->peakOrder[i]->v;
		N+=1.0;
	}
	if (peaks1 != peaks2) {
		peaks1->setDefaultPeakOrder();
		for (int i=0;i<peaks2->numPeaks;i++) {
			avg += peaks2->peakOrder[i]->v;
			N+=1.0;
		}
	}
	avg /= N;
	for (int i=0;i<peaks1->numPeaks;i++) {
		std += (peaks1->peakOrder[i]->v - avg)*(peaks1->peakOrder[i]->v - avg);
	}
	if (peaks1 != peaks2) {
		for (int i=0;i<peaks2->numPeaks;i++) {
			std += (peaks2->peakOrder[i]->v - avg)*(peaks2->peakOrder[i]->v - avg);
		}
	}
	std /= N;
	std = sqrt(std);
	model->setCoverageLimits(avg,std);
	fprintf(stderr, "\tAvg interactions per peak = %.1lf +/- %.1lf\n", avg,std);
	params->relativeFlag = 0;
	params->boundaryScale = 0;

	if (fullModelFlag) {
		model->fullModelFlag = 1;
		fprintf(stderr, "\n\tFinding Interactions to average into expected profile (-fullModel)...\n");
		calcHiCMatrix(peaks1, peaks2,NULL, model->res, model,0,0,actionFlag,totalInteractions,NULL,params);
		model->normalize();

	} else {

        dist = getPETagDistribution(10,1000000000,
                                        model->res,NULL, distLength);	
		if (dist == NULL) {
			fprintf(stderr, "!!! Error- something went wrong with quick calculation of contact frequncies...\n");
			exit(0);
		}
		dist[0] /= totalTags;
		for (int i=1;i<distLength;i++) {
			dist[i] /= 2*totalTags;
		}
		model->useApproximation(dist,distLength);
	}

	model->scale();
	model->refPeaks = refPeaks;
	model->save();

	refPeaks->adjustPETagTotalsForModel(model,numCPUs);
	model->save();
	if (peaks1 != refPeaks) peaks1->adjustPETagTotalsWithModel(model,0);
	if (peaks1 != peaks2) peaks2->adjustPETagTotalsWithModel(model,0);

	model->stdFlag = 1;
	actionFlag = actionFlag | HIC_MASK_NORM_SEQDEPTH;
	actionFlag = actionFlag | HIC_MASK_NORM_DISTANCE;

	fprintf(stderr, "\n\tCalculating Variation:\n");
	if (fullModelFlag) {
		calcHiCMatrix(peaks1, peaks2,NULL, model->res, model,0,0,actionFlag,totalInteractions,NULL,params);
		model->normalizeStd();
	} else {
		model->setDefaultVariation();
	}

	//delete peaks1;

	model->stdFlag = 0;
	model->save();

	if (peaks1IsGenomeFlag==0) {
		delete refPeaks;
	}

	fprintf(stderr, "\tFinished Generating Background Model\n");
}

void TagLibrary::makeHiCHistogram(char* outputFile, PeakLibrary* peaks, int size, 
				int resolution, HiCBgModel* model, int actionFlag,HiCparams* params) {

	peaks->setPeakSize(resolution);
	int halfBins = size/resolution/2;
	int matrixLength = halfBins*2+1;

	double totalInteractions = totalTags;
	fprintf(stderr, "\ttotalInteractions: %lf\n", totalInteractions);
	PeakLibrary* histPeaks = new PeakLibrary(10*halfBins*peaks->numPeaks);

	char** keys = peaks->peaks->keys();
	for (int i=0;i<peaks->peaks->total;i++) {
		Peak* p = (Peak*)peaks->peaks->search(keys[i]);
		delete [](keys[i]);
		if (p == NULL) continue;
		Peak* np = p->copy();
		char* newPeakNameStr = new char[1000];
		sprintf(newPeakNameStr, "%d:0",i);
		delete [](np->name);
		np->name = newPeakNameStr;
		np->refPos = (np->start + np->end)/2;
		histPeaks->addPeak(np);

		for (int j=1;j<=halfBins;j++) {
			int x = j;	
			if (np->strand == STRAND_NEGATIVE) x *= -1;
	
			sprintf(newPeakNameStr, "%d:%d",i,j);
			np->name = newPeakNameStr;
			np->start = p->start + x*resolution;
			np->end = p->end + x*resolution;
			np->refPos = (np->start + np->end)/2;
			histPeaks->addPeak(np);

			sprintf(newPeakNameStr, "%d:-%d",i,j);
			np->name = newPeakNameStr;
			np->start = p->start + -1*x*resolution;
			np->end = p->end + -1*x*resolution;
			np->refPos = (np->start + np->end)/2;
			histPeaks->addPeak(np);
		}
		delete np;
	}
	delete []keys;
	histPeaks->sortChr();
	histPeaks->setDefaultPeakOrder();

	double **matrix = new double*[matrixLength*2];
	for (int i=0;i<matrixLength*2;i++) {
		matrix[i] = new double[matrixLength];
		for (int j=0;j<matrixLength;j++) {
			matrix[i][j]=0.0;
		}
	}

	// pre calculate the total number of interactions per location
	double accuracy = 0;
	if (actionFlag & HIC_MASK_NORM_SEQDEPTH && model != NULL) {
		accuracy = histPeaks->adjustPETagTotalsWithModel(model,1);
	}
	if (accuracy < 0.98) {
		fprintf(stderr, "\tCalculating PE Tag densities (poor accuracy: %.2lf%%)\n",accuracy*100.0);
		getPETagTotals(histPeaks);
		if (actionFlag & HIC_MASK_NORM_DISTANCE && model != NULL) {
			accuracy = histPeaks->adjustPETagTotalsWithModel(model,0);
		}
	} else {
		fprintf(stderr, "\tUsing model to derive PE Tag densities (accuracy: %.2lf%%)\n",accuracy*100.0);
	}
	if (model != NULL && !(actionFlag & HIC_MASK_LOGPVALUES)) {
		model->initializeExpectMatrix(matrixLength,matrixLength);
		actionFlag = actionFlag | HIC_MASK_RAWANDEXPECTED;
	}

	double avg = 0.0;
	double std = 0.0;
	for (int i=0;i<histPeaks->numPeaks;i++) {
		avg += histPeaks->peakOrder[i]->v;
	}
	avg /= histPeaks->numPeaks;
	for (int i=0;i<histPeaks->numPeaks;i++) {
		std += (histPeaks->peakOrder[i]->v - avg)*(histPeaks->peakOrder[i]->v - avg);
	}
	std /= histPeaks->numPeaks;
	std = sqrt(std);
	model->setCoverageLimits(avg,std);
	fprintf(stderr, "\tAvg interactions per peak = %.1lf +/- %.1lf\n", avg,std);



	calcHiCMatrix(histPeaks,histPeaks,matrix,resolution,model,0,matrixLength,
									actionFlag,totalInteractions,NULL,params);

	fprintf(stderr, "\tMatrix: %d x %d\n", matrixLength,matrixLength);


	if (model != NULL && !(actionFlag & HIC_MASK_LOGPVALUES) && actionFlag & HIC_MASK_NORM_SEQDEPTH) {
		double **expectMatrix = model->expectMatrix;
		for (int i=0;i<matrixLength;i++) {
			for (int j=0;j<matrixLength;j++) {
				matrix[i][j] = log((matrix[i][j]+HIC_RATIO_PSEUDO_COUNT)/
									(expectMatrix[i][j]+HIC_RATIO_PSEUDO_COUNT))/log(2.0);
			}
		}
	}
	for (int i=0;i<matrixLength;i++) {
		for (int j=0;j<matrixLength;j++) {
			if (matrix[i+matrixLength][j] > 0) {
				//matrix[i][j] /= matrix[i+matrixLength][j];
			}
		}
	}

	fprintf(stderr, "\tPreparing matrix...\n");

	FILE* fp = NULL;
	if (outputFile != NULL) {
		fp = fopen(outputFile, "w");
		if (fp == NULL) {
			fprintf(stderr, "!!! Could not open %s for writing !!!\n", outputFile);
			exit(0);
		}
	} else {
		fp = stdout;
	}
	fprintf(fp, "HiCMatrix (directory=%s)\t2D Histogramn (%d total regions)",directory,
										peaks->numPeaks);
	for (int i=-halfBins;i<=halfBins;i++) {
		fprintf(fp,"\t%d", i*resolution);
	}
	fprintf(fp,"\n");
	for (int i=0;i<matrixLength;i++) {
		fprintf(fp,"%d\t%d",(-halfBins+i)*resolution,(-halfBins+i)*resolution);
		for (int j=0;j<matrixLength;j++) {
			fprintf(fp,"\t%le",matrix[i][j]);
		}
		fprintf(fp,"\n");
	}
	if (outputFile != NULL) {
		fclose(fp);
	}


	for (int i=0;i<matrixLength;i++) {
		delete [](matrix[i]);
	}
	delete []matrix;
}

void TagLibrary::scoreInteractionBoundaries(PeakLibrary* peaks1,
				int resolution, HiCBgModel* model, int actionFlag) {
	int peakIndex1 = 0;
	peaks1->setDefaultPeakOrder();

	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) { 
		ChrTags* ct = (ChrTags*)chrs->search(keys[i]);
		ct->scoreInteractionBoundaries(peaks1, peakIndex1,resolution, model,actionFlag);
	}
}

double** TagLibrary::makeHiCMatrix(PeakLibrary* peaks1, PeakLibrary* peaks2, 
				int resolution, HiCBgModel* model, int actionFlag,
				GenomeInteractionLibrary* gil, HiCparams* params) {


	int relativeFlag = params->relativeFlag;
	peaks1->setDefaultPeakOrder();
	int numPeaks1 = peaks1->numPeaks;
	int numPeaks2 = numPeaks1;
	if (relativeFlag != 0) {
		numPeaks2 = relativeFlag*2+1;
		peaks2= peaks1;
	} else if (peaks2!=NULL) {
		peaks2->setDefaultPeakOrder();
		numPeaks2 = peaks2->numPeaks;
	} else {
		peaks2 = peaks1;
	}
	double totalInteractions = totalTags;
	fprintf(stderr, "\tTotal Interactions: %.1lf\n", totalInteractions);

	double **matrix = NULL;
	if (actionFlag & HIC_MASK_NO_MATRIX) {
	} else {
		matrix = new double*[numPeaks1];
		for (int i=0;i<numPeaks1;i++) {
			matrix[i] = new double[numPeaks2];
			for (int j=0;j<numPeaks2;j++) {
				matrix[i][j]=0.0;
			}
		}
	}

	int revFlag = 0;
	Peak** peakArray1 = peaks1->peakOrder;
	Peak** peakArray2 = peaks2->peakOrder;

	if (actionFlag & HIC_MASK_NORM_SEQDEPTH) {

		double accuracy = 0.0;
		if (actionFlag & HIC_MASK_NORM_DISTANCE && model != NULL) {
			accuracy = peaks1->adjustPETagTotalsWithModel(model,1);
		} else if (actionFlag & HIC_MASK_NORM_SEQDEPTH && model != NULL) {
			accuracy = peaks1->adjustPETagTotalsWithModel(model,2);
		}
		fprintf(stderr, "\tModel accuracy: %lf\n", accuracy);
		if (accuracy < 0.95) {
			getPETagTotals(peaks1);
			if (actionFlag & HIC_MASK_NORM_DISTANCE && model != NULL) {
				accuracy = peaks1->adjustPETagTotalsWithModel(model,0);
			}
		} else {
			fprintf(stderr, "\tUsing model to derive PE Tag densities\n");
		}

		if (peaks1 != peaks2) {
			fprintf(stderr, "\n\tNon-symmetric matrix, Calculating Read Coverage for 2nd Set of Regions...\n");
			accuracy = 0.0;
			if (actionFlag & HIC_MASK_NORM_DISTANCE && model != NULL) {
				accuracy = peaks2->adjustPETagTotalsWithModel(model,1);
			} else if (actionFlag & HIC_MASK_NORM_SEQDEPTH && model != NULL) {
				accuracy = peaks2->adjustPETagTotalsWithModel(model,2);
			}
			fprintf(stderr, "\tModel accuracy: %lf\n", accuracy);
			if (accuracy < 0.95) {
				getPETagTotals(peaks2);
				if (actionFlag & HIC_MASK_NORM_DISTANCE && model != NULL) {
					accuracy = peaks2->adjustPETagTotalsWithModel(model,0);
				}
			} else {
				fprintf(stderr, "\tUsing model to derive PE Tag densities for 2nd set of peaks\n");
			}
		}	
		fprintf(stderr, "\n");

		double avg = 0.0;
		double std = 0.0;
		for (int i=0;i<numPeaks1;i++) {
			avg += peakArray1[i]->v;
		}
		int localNumPeaks2 = numPeaks2;
		if (relativeFlag) {
			localNumPeaks2 = numPeaks1;
		}
		for (int i=0;i<localNumPeaks2;i++) {
			avg += peakArray2[i]->v;
		}
		avg /= (numPeaks1+localNumPeaks2);
		for (int i=0;i<numPeaks1;i++) {
			std += (peakArray1[i]->v-avg)*(peakArray1[i]->v-avg);
		}
		for (int i=0;i<localNumPeaks2;i++) {
			std += (peakArray2[i]->v-avg)*(peakArray2[i]->v-avg);
		}
		std /= (numPeaks1+localNumPeaks2);
		std = sqrt(std);
		fprintf(stderr, "\tAverage interaction count in regions = %.1lf +/- %.1lf\n", avg,std);
		model->setCoverageLimits(avg,std);
	}


	//calculate values for the matrix
	revFlag = 0;
	calcHiCMatrix(peaks1,peaks2,matrix,resolution,model,revFlag,0,actionFlag,
							totalInteractions,gil,params);


	if (0 && peaks1 != peaks2) {
		fprintf(stderr, "\n\tNon-symmetric matrix, reversing the calculations...\n");
		revFlag = 1;
		calcHiCMatrix(peaks2,peaks1,matrix,resolution,model,revFlag,0,actionFlag,
								totalInteractions,gil,params);

	}
	if (gil != NULL && gil->recordInteractions) {
		if (gil->recordInteractionBg == 0) {
			gil->totalTests = (double)numPeaks1*(double)numPeaks2/2.0;
			gil->optimizeInteractionArray();
		}
		fprintf(stderr, "\tRegions with too many or too few reads: %.2lf%%\n",
					100.0*((double)model->badRegions)/(double)(model->goodRegions+model->badRegions));
		fprintf(stderr, "\tTotal Significant Interactions: %d\n", gil->numInteractions);
	}

	return matrix;
}


double** TagLibrary::calcCorrelationMatrix(double **matrix, int numPeaks1, double **expectedMatrix,
											double minReads,int maxIndexDiff) {
	double **corrMatrix = new double*[numPeaks1];

	//expectedMatrix = NULL;
	if (expectedMatrix != NULL) {
		fprintf(stderr, "\tPooling regions during correlation calculation (need >%.1lf expected reads)\n", minReads);
	}

	int limit = maxIndexDiff/2;
	if (maxIndexDiff < 0) {
		limit = INT_MAX;
	}

	double barSize = 80.0;
	double progressCount=0.0;
	double nextThresh = 1.0/barSize;
	fprintf(stderr, "\tCalculating Correlation Matrix:\n");
	fprintf(stderr, "\t|0%%                                    |50%%                                100%%|\n\t");
	for (int i=0;i<numPeaks1;i++) {
		progressCount += 1.0/((double)numPeaks1);
		while (progressCount >= nextThresh) {
			fprintf(stderr, "=");
			nextThresh += 1.0/barSize;
		}

		corrMatrix[i] = new double[numPeaks1];
		for (int j=0;j<numPeaks1;j++) {

			if (abs(i-j) > limit) {
				corrMatrix[i][j] = 0.0;
				continue;
			}

			double* expectI = NULL;
			double* expectJ = NULL;
			if (expectedMatrix != NULL) {
				expectI = expectedMatrix[i];
				expectJ = expectedMatrix[j];
			}
			corrMatrix[i][j] = calcCorrelation(matrix[i],matrix[j],numPeaks1,expectI,expectJ,minReads,maxIndexDiff,i,j);
			if (corrMatrix[i][j] < -1.9) corrMatrix[i][j] = 0.0;
		}
	}
	fprintf(stderr, "\n");
	for (int i=0;i<numPeaks1;i++) {
		delete [](matrix[i]);
	}
	delete []matrix;
	return corrMatrix;
}
double TagLibrary::calcCorrelation(double* a1, double* a2, int n, double* e1, double* e2, double min,int maxIndexDiff,
										int refIndexI, int refIndexJ) {


	if (maxIndexDiff >=0 && abs(refIndexI-refIndexJ) > maxIndexDiff/2) return 0.0;
	double expectI = 0.0;
	double expectJ = 0.0;
	double curValueI = 0.0;
	double curValueJ = 0.0;
	double curN = 0.0;
	static double* v1 = NULL;
	static double* v2 = NULL;
	int N = 0;
	static int vLength = 0;
	if (v1 == NULL || vLength < n) {
		if (v1 != NULL) delete []v1;
		if (v2 != NULL) delete []v2;
		v1 = new double[n];
		v2 = new double[n];
		vLength = n;
	}
	double avgSize = 0.0;

	int cstart = 0;
	int cend = n;
	if (maxIndexDiff >= 0) {
		cstart = refIndexJ-maxIndexDiff;
		cend = refIndexI+maxIndexDiff;
		if (refIndexJ < refIndexI) {
			cstart = refIndexI-maxIndexDiff;
			cend = refIndexJ+maxIndexDiff;
		}
		if (cstart < 0) cstart = 0;
		if (cend > n) cend = n;
	}

	for (int k=cstart;k<cend;k++) {
		if (1) { //if (i != k && j != k && i != j) {
			
			if ((fabs(a1[k]) > 1e-10 || fabs(a2[k]) > 1e-10)) {
				curValueI += a1[k];
				curValueJ += a2[k];
				curN++;
			}

			int resetFlag = 0;
			int addFlag = 0;
			if (e1 != NULL) {
				expectI += e1[k];
				expectJ += e2[k];
				if (expectI >= min && expectJ >= min) {
					if (curN > 0) addFlag = 1;
					resetFlag = 1;
				}
			} else {
				if (curN > 0) addFlag = 1;
				resetFlag = 1;
			}
			
			if (addFlag) {
				v1[N] = curValueI/curN;
				v2[N] = curValueJ/curN;
				N++;
				avgSize += curN;
				//fprintf(stderr, "%d %d %d %lf %lf %lf\n", i,j,k,curN,expectI,expectJ);
			}
			if (resetFlag) {
				curValueI=0.0;
				curValueJ=0.0;
				curN=0.0;
				expectI=0.0;
				expectJ=0.0;
			}
		}
	}
	double rv = -2.0;
	if (N > 1) {
		rv = (double)correlation(v1,v2,N);
		//fprintf(stderr, "Avg Size = %.1lf\n", ((double)avgSize)/((double)N));
	}
	return rv;
}



// class ChrTags ------------------------------------------------------------


ChrTags::ChrTags(char* newchr) {
	chr = NULL;
	if (newchr != NULL) {
		chr = new char[strlen(newchr)+1];
		strcpy(chr, newchr);
	}
	linkTag = NULL;
	firstTag = NULL;
	tags = NULL;
	petags = NULL;
	mCflag = 0;
	totalTags = 0.0;
	ogTotalTags = 0.0;
	totalPositions = 0;
	adjusted = 0;
	optimizeOverride = 0;
	maxtbp = 0.0;
	mintbp = 0.0;
	tagAdjust = 0;
	numLinks = 0;
	loaded = 0;
	memType = TAGLIBRARY_LOWMEM;
	tagFile = NULL;
	tagfp = NULL;
	tagfpR1 = NULL;
	tagfpR2 = NULL;
	dontSAVE = 0;
	singleFile = 0;
	gcFreq = NULL;
	seq = NULL;
	chrLen = 0;
	chrNames = new Hashtable(1000);
	firstPETag = NULL;
	linkPETag = NULL;
	pairedEndFlag = 0;
	forceSingleReadFlag = 0;
	maxPosition = -1;
	appearentSize = 0;
	revStrand = 0;
}

ChrTags::~ChrTags() {
	if (tags != NULL) delete []tags;
	if (chr != NULL) delete []chr;
	if (tagFile != NULL) delete []tagFile;
	if (tagfp != NULL) fclose(tagfp);
	if (gcFreq != NULL) delete []gcFreq;
	if (seq != NULL) delete []seq;
	linkTag = NULL;
	firstTag = NULL;
	firstPETag = NULL;
	linkPETag = NULL;
	tags = NULL;
	if (chrNames != NULL) {
		char** keys = chrNames->keys();
		for (int i=0;i<chrNames->total;i++) {
			char* cname = (char*) chrNames->search(keys[i]);
			if (cname!=NULL) delete []cname;
			delete [](keys[i]);
		}
		delete []keys;
		delete chrNames;
	}
}

void ChrTags::freeTags() {
	if (tags != NULL) delete []tags;
	if (petags != NULL) delete []petags;
	tags = NULL;
	petags = NULL;
	loaded = 0;
	adjusted = 0;
}

void ChrTags::setTagFile(char* file) {
	if (file == NULL) return;
	tagFile = new char[strlen(file)+1];
	strcpy(tagFile, file);
}
void ChrTags::openTagFile(char* mode) {
	if (tagFile == NULL) return;
	closeTagFile();
	tagfp = fopen(tagFile, mode);
	if (tagfp == NULL) {
		fprintf(stderr, "!!!!! Could not open file %s for printing tags!!!!!\n", tagFile);
		fprintf(stderr, "!!!!! Is this a valid file name?  May need to:\n");
		fprintf(stderr, "\tTry a different name for the tag directory (is the same name as an existing file?)\n");
		fprintf(stderr, "\t-or- you may need to rename your chromosomes if they have weird characters!\n");
		fprintf(stderr, "\t-or- try using the \"-single\" flag.\n");
		exit(0);
	}
}

void ChrTags::closeTagFile() {
	if (tagfp != NULL) fclose(tagfp);
	if (tagfpR1 != NULL) fclose(tagfp);
	if (tagfpR2 != NULL) fclose(tagfp);
	tagfp = NULL;
	tagfpR1 = NULL;
	tagfpR2 = NULL;
}
void ChrTags::printAlignedTag(char* rname, int pos, char dir, int length, float value,int PEflag) {
	FILE *fp = tagfp;
	if (PEflag == TAGLIBRARY_PE_READ1_FLAG) {
		fp = tagfpR1;
	} else if (PEflag == TAGLIBRARY_PE_READ2_FLAG) {
		fp = tagfpR2;
	}
	if (Tag::precision == 1) {
		fprintf(fp,"%s\t%s\t%d\t%d\t%.1f\t%d\n",rname,chr,pos,dir,value,length);
	} else if (Tag::precision == 2) {
		fprintf(fp,"%s\t%s\t%d\t%d\t%.2f\t%d\n",rname,chr,pos,dir,value,length);
	} else {
		fprintf(fp,"%s\t%s\t%d\t%d\t%.3f\t%d\n",rname,chr,pos,dir,value,length);
	}
}
void ChrTags::printAlignedPETag(PETag* petag, int revFlag) {
	petag->print(tagfp,revFlag);
}
int ChrTags::getMaxPosition() {
	openTagFile((char*)"r");
	if (tagfp == NULL) {
		fprintf(stderr, "Could not open %s\n",tagFile);
		return -1;
	}
	int ok = fseek(tagfp,-10000,SEEK_END); //should be enough to see several lines at the end
	if (ok) {
		ok = fseek(tagfp,0,SEEK_SET);
	}
	char* buf = new char[BUFFER];
	char** line = new char*[BUFFER];
	int numCols=0;
	int pos = 0;

	maxPosition = -1;

	int lineCount = 0;
	while (fgets(buf, BUFFER, tagfp) != NULL) {
		//fprintf(stderr, "%s\n", buf);
		split(buf, line, numCols,'\t');
		if (numCols < 5) continue;
		lineCount++;
		if (lineCount < 2) continue; // don't want partial lines

		sscanf(line[2],"%d",&pos);
		if (pos > maxPosition) maxPosition = pos;
	}

	delete []buf;
	delete []line;
	closeTagFile();	
	return maxPosition;
}



void ChrTags::readTagFile(char* file) {
	setTagFile(file);
	readTagFile();
}
void ChrTags::loadTags() {
	if (loaded && adjusted) return;
	if (pairedEndFlag && forceSingleReadFlag == 1) {
		optimizeOverride=1;
		//fprintf(stderr, "Reading single...\n");
	}
	if (singleFile == 0 && loaded == 0) {
		readTagFile();
		ogTotalTags = totalTags;
	}
	adjustTags();
	if (pairedEndFlag && forceSingleReadFlag == 1) optimizeOverride=0;
}
void ChrTags::readTagFile() {
	openTagFile((char*)"r");
	if (tagfp == NULL) {
		fprintf(stderr, "Could not open %s\n",tagFile);
		return;
	}
	char* buf = new char[BUFFER];
	char* name = new char[BUFFER];	
	char** line = new char*[BUFFER];
	int numCols = 0;
	int pos = 0;
	float tagCount= 0;
	int dir = 0;
	int len = -1;
	char* chr1 = NULL;
	char* chr2 = NULL;
	int pos2 = 0;
	int dir2 = 0;
	int len2 = -1;
	appearentSize = 0;

	optimizedFlag=1;
	int lastPos = -1000000000;

	while (fgets(buf, BUFFER, tagfp) != NULL) {
		//sscanf(buf, "%s\t%s\t%d\t%d\t%d",name,chr,&pos,&dir,&tagCount);
		split(buf, line, numCols,'\t');
		if (numCols < 5) continue;

		sscanf(line[2],"%d",&pos);
		sscanf(line[3],"%d",&dir);
		sscanf(line[4],"%f",&tagCount);
		//if (tagCount < 0.00001) continue;
		len = -1;
		if (numCols > 5) {
			sscanf(line[5],"%d",&len);
		}
		if (pos > appearentSize) appearentSize = (long long int) pos;

		if (pairedEndFlag && forceSingleReadFlag == 0) {
			if (numCols < 10) {
				fprintf(stderr, "!!! Expecting paired end tags (%s) !!!\n",tagFile);
				exit(0);
			}
			chr1 = (char*)chrNames->search(line[1]);	
			if (chr1 == NULL) {
				char* cc = new char[strlen(line[1])+1];
				strcpy(cc, line[1]);
				chrNames->insert(cc,line[1]);
				chr1 = cc;
			}
			chr2 = (char*)chrNames->search(line[6]);	
			if (chr2 == NULL) {
				char* cc = new char[strlen(line[6])+1];
				strcpy(cc, line[6]);
				chrNames->insert(cc,line[6]);
				chr2 = cc;
			}
			sscanf(line[7],"%d",&pos2);
			sscanf(line[8],"%d",&dir2);
			sscanf(line[9],"%d",&len2);
			addPETag(chr1,pos,dir,len,chr2,pos2,dir2,len2,tagCount);
		
		} else {
			if (revStrand) {
				if (dir == STRAND_POSITIVE) {
					pos += len;
					dir = STRAND_NEGATIVE;
				} else {
					pos -= len;
					dir = STRAND_POSITIVE;
				}
				if (pos < 0) pos = 0;
				if (pos > appearentSize) appearentSize = (long long int) pos;
			}
			addTag(pos,(char)dir,len,tagCount);
		}

		if (pos < lastPos) optimizedFlag = 0;
		lastPos = pos;
	}


	//fprintf(stderr, "appearentSize = %lld\n", appearentSize);

	delete []buf;
	delete []name;
	delete []line;

	closeTagFile();

	optimizeTags();
	loaded =1;
}

void ChrTags::deleteSequence() {
	if (allSeqs != NULL) {
		char* storedSeq = (char*)allSeqs->search(chr);
		if (storedSeq != NULL && storedSeq != seq) {
			delete []storedSeq;
		}
		allSeqs->insert(NULL, chr);
	}
	if (seq != NULL) delete []seq;
	seq = NULL;
	chrLen = 0;
}
void ChrTags::getSequence(char* dest, int start, int end, char strand) {

	char* s = NULL;
	int len = 0;
	if (start >= 1 && end < chrLen) {
		s = &(seq[start-1]);
		len = end-start+1;
	}
	//fprintf(stderr, "%d\t%d\t%d\t%d\n", start, end, strand,chrLen);
	if (s != NULL && len > 0) {
		strncpy(dest, s, len);
		dest[len]='\0';
		if (strand == STRAND_POSITIVE) {
		} else if (strand == STRAND_NEGATIVE) {
			revopp(dest);
		}
		//fprintf(stderr, "%s\t%d\t%d\t%d\t%d\n%s\n", chr,start,end,strand,len,dest);
	} else {
		dest[0]='\0';
	}
}

//static variable
Hashtable* ChrTags::allSeqs = NULL;

void ChrTags::loadSequence(char* genomeDir) {


	if (allSeqs != NULL) {
		seq = (char*) allSeqs->search(chr);	
		chrLen = 0;
		if (seq != NULL) chrLen = strlen(seq);
		return;
	}

    int dirMode = 1;
    struct stat st_buf;
    stat(genomeDir,&st_buf);
    if (S_ISREG(st_buf.st_mode)) {
        dirMode = 0;
        //fprintf(stderr, "\tExtracting sequences from file: %s\n", genomeDir);
    } else {
        //fprintf(stderr, "\tExtracting sequences from directory: %s\n", genomeDir);
    }
	FILE* fp = NULL;
	int singleGenomeFile = 0;
	if (dirMode) {
		char* fname = new char[100000];
		sprintf(fname,"%s/%s.fa",genomeDir,chr);
		fp = fopen(fname,"r");
		if (fp == NULL) {
			sprintf(fname,"%s/%s.fa.masked",genomeDir,chr);
			fp = fopen(fname,"r");
		}
		if (fp == NULL) {
			sprintf(fname,"%s/genome.fa",genomeDir);
			fp = fopen(fname,"r");
			singleGenomeFile = 1;
		}
		if (fp == NULL) {
			sprintf(fname,"%s/genome.fa.masked",genomeDir);
			fp = fopen(fname,"r");
			singleGenomeFile = 1;
		}
		if (fp == NULL) {
			fprintf(stderr, "Couldn't open a sequence file for \"%s\" (%s)\n", chr,fname);
			delete []fname;
			return;
		}
		delete []fname;
	} else {
		fp = fopen(genomeDir,"r");
		if (fp == NULL) {
			fprintf(stderr, "Couldn't open the sequence file %s\n", genomeDir);
			return;
		}
		singleGenomeFile = 1;
	}

	if (singleGenomeFile) {
		allSeqs = new Hashtable();
		fprintf(stderr, "\tReading whole genome at once...\n");
	}	

	chrLen=0;
	if (seq != NULL) delete []seq;
	unsigned int likelyMaximum = 400000000;
	if (totalPositions > 0  && tags != NULL) {
		likelyMaximum = tags[totalPositions-1].p + 10000000;
	}
	if (singleGenomeFile) {
		likelyMaximum = 1000000000;
	}
	seq = new char[likelyMaximum];
	seq[0] = '\0';

	char* buffer = new char[BUFFER];
	int goodSeq = 0;
	int found = 0;
	char* curName= new char[BUFFER];

	while (fgets(buffer,BUFFER,fp)!=NULL) {
		int lineLen = strlen(buffer);
		if (lineLen > 0 && buffer[lineLen-1] == '\n') {
			buffer[lineLen-1] = '\0';
			lineLen--;
		}
		if (buffer[0] == '>') {
			if (singleGenomeFile && seq[0] != '\0') {
				char* seq2 = new char[strlen(seq)+1];
				strcpy(seq2,seq);
				allSeqs->insert(seq2, curName);
			}
			int count = 1;
			chrLen = 0;
			while (buffer[count] != '\0') {
				if (buffer[count] < 33) {
					buffer[count] = '\0';
					break;
				}
				count++;
			}
			strcpy(curName,&(buffer[1]));
			if (singleGenomeFile || strcmp(chr, curName)==0) {
				goodSeq = 1;
				seq[0] = '\0';
			} else {
				if (goodSeq == 1) {
					//done
					break;
				}
				goodSeq = 0;
			}
			continue;
		}
		if (goodSeq) {
			found=1;
			if (lineLen+chrLen > (int)likelyMaximum) break;

			for (int i=0;i<lineLen;i++) {
				if (buffer[i] == 'a') buffer[i]='A';
				else if (buffer[i] == 'c') buffer[i]='C';
				else if (buffer[i] == 'g') buffer[i]='G';
				else if (buffer[i] == 't') buffer[i]='T';
				else if (buffer[i] == 'n') buffer[i]='N';
			}

			memcpy(&(seq[chrLen]),buffer,lineLen);
			chrLen+=lineLen;
			seq[chrLen]='\0';
			//fprintf(stderr, "chrlen=%d, chrLen=%d\n",strlen(seq),chrLen);
		}
	}
	fclose(fp);
	if (singleGenomeFile && seq[0] != '\0') {
		char* seq2 = new char[strlen(seq)+1];
		strcpy(seq2,seq);
		allSeqs->insert(seq2, curName);
	}
	
	if (singleGenomeFile) {
		chrLen = 0;
		seq = (char*) allSeqs->search(chr);
		if (seq != NULL) chrLen = strlen(seq);
		return;
	}

	delete []buffer;
	delete []curName;
	if (found==0) {
		delete []seq;
		seq = NULL;
	}

}

void ChrTags::addTag(int pos, char dir, int length, float value) {
	if (linkTag == NULL) {
		firstTag = new LinkedTag(pos,dir,length,value,NULL);
		linkTag = firstTag;
	} else {
		LinkedTag* next = new LinkedTag(pos,dir,length,value, NULL);
		linkTag->tag = next;
		linkTag = next;
	}
	numLinks++;
}
void ChrTags::addPETag(char* c1, int p1, char d1, int len1, char* c2, int p2, 
														char d2, int len2,float v) {
	if (linkPETag == NULL) {
		firstPETag = new LinkedPETag(c1,p1,d1,len1,c2,p2,d2,len2,v,NULL);
		linkPETag = firstPETag;
	} else {
		LinkedPETag* next = new LinkedPETag(c1,p1,d1,len1,c2,p2,d2,len2,v,NULL);
		linkPETag->tag = next;
		linkPETag = next;
	}
	numLinks++;
}
void ChrTags::print() {
	if (dontSAVE) {
		fprintf(stderr, "\t!!!! Tags were adjusted - shouldn't save them!!! Error in code\n");
		return;
	}
	if (singleFile) {
		for (int i=0;i<totalPositions;i++) {
			tags[i].print(tagfp,chr);
		}
	} else {
		openTagFile((char*)"w");
		if (pairedEndFlag) {
			for (int i=0;i<totalPositions;i++) {
				petags[i].print(tagfp);
			}
		} else {
			for (int i=0;i<totalPositions;i++) {
				tags[i].print(tagfp,chr);
			}
		}
	}
}

void ChrTags::optimizeTagFile() {
	optimizeOverride = 1;
	readTagFile();
	adjustTags();
	dontSAVE=0;
	print();
	freeTags();
	optimizeOverride = 0;
}

void ChrTags::optimizeTags() {

	if (pairedEndFlag && forceSingleReadFlag==0) {
		optimizePETags();
		return;
	}

	if (numLinks == 0 || linkTag == NULL) return;
	
	totalPositions = 0;
	totalTags = 0;
	if (tags != NULL) {
		delete []tags;
	}
	tags = new Tag[numLinks];

	LinkedTag* current = firstTag;
	LinkedTag* next = NULL;
	while (current != NULL) {
		tags[totalPositions].copy(current);
		totalTags += tags[totalPositions].v;
		next = current->tag;
		delete current;
		totalPositions++;
		current = next;
	}
	linkTag = NULL;
	firstTag = NULL;
	numLinks = 0;

	//fprintf(stderr, "%s  - optimizedFlag = %d\n", chr, optimizedFlag);
	if (optimizedFlag && !optimizeOverride) return;
	
	qsort(tags,totalPositions,sizeof(Tag),&cmpTags);
	int last = 0;
	totalTags = 0;
	if (totalPositions > 0) {
		totalTags = tags[0].v;
	} else {
		fprintf(stderr, "No tags in tag file!!!\n");
		return;
	}	
	for (int i=1;i<totalPositions;i++) {
		totalTags += tags[i].v;
		int j=i-1;
		if (mCflag) {
			if (tags[i].p == tags[j].p && tags[i].d == tags[j].d) {
				double sum = tags[last].v*tags[last].len + tags[i].v*tags[i].len;
				double d = (double)(tags[last].len+tags[i].len);
				if (d > 0.0) sum /= d;
				tags[last].v = sum;
				tags[last].len = (int)d;
			} else {
				last++;
				if (last == i) continue;
				tags[last].copy(&(tags[i]));
			}
		} else {
			if (tags[i].p == tags[j].p && tags[i].d == tags[j].d && tags[i].len == tags[j].len) {
				tags[last].v += tags[i].v;
			} else {
				last++;
				if (last == i) continue;
				tags[last].copy(&(tags[i]));
			}
		}
	}
	totalPositions = last+1;

	if (maxtbp > 0.0) {
		totalTags = 0.0;
		for (int i=0;i<totalPositions;i++) {
			if (tags[i].v > maxtbp) tags[i].v = maxtbp;
			totalTags += tags[i].v;
		}
	}
	if (mintbp > 0.0) {
		totalTags = 0.0;
		for (int i=0;i<totalPositions;i++) {
			if (tags[i].v < mintbp) tags[i].v = mintbp;
			totalTags += tags[i].v;
		}
	}
}
void ChrTags::optimizePETags() {

	if (numLinks == 0 || linkPETag == NULL) return;

	totalPositions = 0;
	totalTags = 0;
	if (petags != NULL) {
		delete []petags;
	}
	petags = new PETag[numLinks];

	LinkedPETag* current = firstPETag;
	LinkedPETag* next = NULL;
	while (current != NULL) {
		petags[totalPositions].copy(current);
		totalTags += petags[totalPositions].v;
		next = current->tag;
		delete current;
		totalPositions++;
		current = next;
	}
	linkPETag = NULL;
	firstPETag = NULL;
	numLinks = 0;

	//fprintf(stderr, "%s  - optimizedFlag = %d\n", chr, optimizedFlag);
	if (optimizedFlag && !optimizeOverride) return;
	//fprintf(stderr, "Optimizing\n");
	qsort(petags,totalPositions,sizeof(PETag),&cmpPETags);

	int last = 0;
	totalTags = 0;
	if (totalPositions > 0) {
		totalTags = petags[0].v;
	}	
	for (int i=1;i<totalPositions;i++) {
		totalTags += petags[i].v;
		int j=i-1;
		if (j>=0 && petags[i].p1 == petags[j].p1 && petags[i].d1 == petags[j].d1 &&
					strcmp(petags[i].chr2,petags[j].chr2)==0 && petags[i].p2 == petags[j].p2 && 
					petags[i].d2 == petags[j].d2) {
			petags[last].v += petags[i].v;
		} else {
			last++;
			if (last == i) continue;
			petags[last].copy(&(petags[i]));
		}
	}
	totalPositions = last+1;

	if (maxtbp > 0.0) {
		totalTags = 0.0;
		for (int i=0;i<totalPositions;i++) {
			if (petags[i].v > maxtbp) petags[i].v = maxtbp;
			totalTags += petags[i].v;
		}
	}
	if (mintbp > 0.0) {
		totalTags = 0.0;
		for (int i=0;i<totalPositions;i++) {
			if (petags[i].v < mintbp) petags[i].v = mintbp;
			totalTags += petags[i].v;
		}
	}
}

int cmpTags(const void* a, const void* b) {
	int ap = ((Tag*)a)->p;
	int bp = ((Tag*)b)->p;
	if (ap < bp) return -1;
	if (ap > bp) return 1;
	char ad = ((Tag*)a)->d;
	char bd = ((Tag*)b)->d;
	if (ad < bd) return -1;
	if (ad > bd) return 1;
	char al = ((Tag*)a)->len;
	char bl = ((Tag*)b)->len;
	if (al < bl) return -1;
	if (al > bl) return 1;
	return 0;
}
int cmpPETags(const void* a, const void* b) {
	int ap = ((PETag*)a)->p1;
	int bp = ((PETag*)b)->p1;
	if (ap < bp) return -1;
	if (ap > bp) return 1;
	char ad = ((PETag*)a)->d1;
	char bd = ((PETag*)b)->d1;
	if (ad < bd) return -1;
	if (ad > bd) return 1;
	char* ac = ((PETag*)a)->chr2;
	char* bc = ((PETag*)b)->chr2;
	int x = strcmp(ac,bc);
	if (x < 0) return -1;
	if (x > 0) return 1;
	ap = ((PETag*)a)->p2;
	bp = ((PETag*)b)->p2;
	if (ap < bp) return -1;
	if (ap > bp) return 1;
	ap = ((PETag*)a)->d2;
	bp = ((PETag*)b)->d2;
	if (ap < bp) return -1;
	if (ap > bp) return 1;
	return 0;
}

int cmp2ndPETagPointers(const void* a, const void* b) {
	char* ac = (*(PETag**)a)->chr2;
	char* bc = (*(PETag**)b)->chr2;
	int x = chrcmp(&ac,&bc);
	if (x != 0) return x;
	int ap = (*(PETag**)a)->p2;
	int bp = (*(PETag**)b)->p2;
	if (ap < bp) return -1;
	if (ap > bp) return 1;
	ap = (*(PETag**)a)->d2;
	bp = (*(PETag**)b)->d2;
	if (ap < bp) return -1;
	if (ap > bp) return 1;
	return 0;
}
int cmp2ndPETags(const void* a, const void* b) {
	char* ac = ((PETag*)a)->chr2;
	char* bc = ((PETag*)b)->chr2;
	int x = chrcmp(&ac,&bc);
	if (x != 0) return x;
	int ap = ((PETag*)a)->p2;
	int bp = ((PETag*)b)->p2;
	if (ap < bp) return -1;
	if (ap > bp) return 1;
	ap = ((PETag*)a)->d2;
	bp = ((PETag*)b)->d2;
	if (ap < bp) return -1;
	if (ap > bp) return 1;
	return 0;
}
int ChrTags::getPETagDistribution(double* sameStrand, double* diffStrand, int windowSize,
						double * largeWindow, int resolution, int largeLength) {

	loadTags();

	int rvLen = 0;
	if (totalPositions > 0) {
		rvLen = petags[totalPositions-1].p1 - petags[0].p1;
	}
	
	int halfWindow = (windowSize)/2;
	for (int i=0;i<totalPositions;i++) {

		float v = petags[i].v;
		
		if (strcmp(petags[i].chr1,petags[i].chr2) != 0) {
			largeWindow[largeLength-1]+=v;
		} else {
			int diff = petags[i].p2 - petags[i].p1;
			int bin = (int)(fabs(((double)diff)/((double)resolution))+0.5);
	//fprintf(stderr, "%d\t%d\t%d\n", diff, bin, resolution);
			if (bin > largeLength-2) {
				largeWindow[largeLength-2]+=v;
			} else {
				largeWindow[bin]+=v;
			}

			if (diff > halfWindow || diff < -1*halfWindow) continue;
			
			if (petags[i].d1 == petags[i].d2) {
				if (petags[i].d1 == 0) {
					sameStrand[halfWindow+diff] += v;
				} else {
					sameStrand[halfWindow-diff] += v;
				}
			} else {
				if (petags[i].d1 == 0) {
					diffStrand[halfWindow+diff] += v;
				} else {
					diffStrand[halfWindow-diff] += v;
				}
			}
		}
	}
	freeTags();
	forceSingleReadFlag = 0;
	return rvLen;
}
void ChrTags::removeTagsInPeaks(PeakLibrary* regions) {
	loadTags();

	int newPositions = 0;
	double newTotal = 0.0;

	int peakIndex = 0;
	Peak** peaks = regions->peakOrder;
	int numPeaks = regions->numPeaks;
	int startPeakIndex = -1;
	int maxPeakIndex = 0;
	for (;peakIndex<numPeaks;peakIndex++) {
		if (strcmp(peaks[peakIndex]->chr,chr) == 0) {
			if (startPeakIndex == -1) startPeakIndex = peakIndex;
			maxPeakIndex = peakIndex+1;
		}
	}
	peakIndex = startPeakIndex;
	int skip = 0;
	if (peakIndex == -1) skip=1;

	if (skip) {
	} else if (pairedEndFlag) {
		PETag* newPETags = NULL;
		newPETags = new PETag[totalPositions];
		for (int i=0;i<totalPositions;i++) {
			int remove = 0;
			if (peakIndex < maxPeakIndex) {
				while (petags[i].p1 > peaks[peakIndex]->end) {
					peakIndex++;
					if (peakIndex >= maxPeakIndex) {
						peakIndex = maxPeakIndex-1;
						break;
					}
				}
				for (int j=peakIndex;j<maxPeakIndex;j++) {
					if (petags[i].p1 < peaks[j]->start) break;
					if (petags[i].p1 >= peaks[j]->start && petags[i].p1 <= peaks[j]->end) {
						remove = 1;
						break;
					}
				}
			}
			if (remove == 0) {
				newPETags[newPositions].copy(&(petags[i]));
				newTotal += newPETags[newPositions].v;
				newPositions++;
			}
		}
		delete []petags;
		petags = newPETags;
		totalPositions = newPositions;
		totalTags = newTotal;
		qsort(petags,totalPositions,sizeof(PETag),&cmp2ndPETags);
		
		newPETags = new PETag[totalPositions];
		newPositions = 0;
		newTotal = 0.0;
		peakIndex = 0;
		maxPeakIndex = numPeaks;
		for (int i=0;i<totalPositions;i++) {
			int remove = 0;
			if (peakIndex < numPeaks) {
				while (chrcmp(&(petags[i].chr2),&(peaks[peakIndex]->chr))<0) {
					peakIndex++;
					if (peakIndex >= maxPeakIndex) {
						peakIndex = maxPeakIndex-1;
						break;
					}
				}
				if (chrcmp(&(petags[i].chr2),&(peaks[peakIndex]->chr))==0) {
					while (petags[i].p2 > peaks[peakIndex]->end) {
						peakIndex++;
						if (peakIndex >= maxPeakIndex) {
							peakIndex = maxPeakIndex-1;
							break;
						}
					}
					for (int j=peakIndex;j<maxPeakIndex;j++) {
						if (petags[i].p2 < peaks[j]->start) break;
						if (petags[i].p2 >= peaks[j]->start && petags[i].p2 <= peaks[j]->end) {
							remove = 1;
							break;
						}
					}
				}
			}
			if (remove == 0) {
				newPETags[newPositions].copy(&(petags[i]));
				newTotal += newPETags[newPositions].v;
				newPositions++;
			}
		}
		delete []petags;
		petags = newPETags;
		totalPositions = newPositions;
		totalTags = newTotal;

	} else {
		Tag* newTags = NULL;
		newTags = new Tag[totalPositions];
		for (int i=0;i<totalPositions;i++) {
			int remove = 0;
			if (peakIndex < maxPeakIndex) {
				while (tags[i].p > peaks[peakIndex]->end) {
					peakIndex++;
					if (peakIndex >= maxPeakIndex) {
						peakIndex = maxPeakIndex-1;
						break;
					}
				}
				for (int j=peakIndex;j<maxPeakIndex;j++) {
					if (tags[i].p < peaks[j]->start) break;
					if (tags[i].p >= peaks[j]->start && tags[i].p <= peaks[j]->end) {
						remove = 1;
						break;
					}
				}
			}
			if (remove == 0) {
				newTags[newPositions].copy(&(tags[i]));
				newTotal += newTags[newPositions].v;
				newPositions++;
			}
		}
		delete []tags;
		tags = newTags;
		totalPositions = newPositions;
		totalTags = newTotal;
	
	}

	freeTags();
}
void ChrTags::getAverageCoverage(int spikeSize, double &chrTotal, int &chrN, PeakLibrary* regions) {

	loadTags();

	char* pname = new char[1000];

	int index = 0;	
	for (int i=0;i<appearentSize;i+=spikeSize) {
		double localTotal = 0;
		if (pairedEndFlag) {
			while (index < totalPositions && petags[index].p1 < i+spikeSize) {
				if (petags[index].p1 >= i) {
					localTotal += petags[index].v;
				}
				index++;
			}
		} else {
			while (index < totalPositions && tags[index].p < i+spikeSize) {
				if (tags[index].p >= i) {
					localTotal += tags[index].v;
				}
				index++;
			}
		}
		chrN++;
		chrTotal += localTotal;

		sprintf(pname,"%s:%d-%d", chr, i, i+spikeSize);

		regions->addPeak(pname,chr,i, i+spikeSize,i+spikeSize/2,STRAND_POSITIVE,localTotal,
						0.0, NULL, -1, 0);
	}
	freeTags();
	delete []pname;
}

void ChrTags::findRestrictionSites(char* site, int maxMisMatches, Hashtable* sites, char* genomeDirectory) {
	loadSequence(genomeDirectory);

	int siteLength = strlen(site);
	char* rvSite = new char[siteLength+1];
	strcpy(rvSite,site);
	revopp(rvSite);
	//int rvFlag = 1;
	if (strcmp(rvSite,site)==0) {
		//rvFlag = 0;
	}
	int halfSite = siteLength/2;

	int threshold = siteLength-maxMisMatches;

	int* sitePos = new int[chrLen/10];
	int numSites = 0;
	for (int i=0;i<chrLen-siteLength;i++) {

		int score=siteLength;
		for (int j=0;j<siteLength;j++) {
			if (seq[i+j]!=site[j] && site[j] != 'N') {
				score--;
				if (score < threshold) break;
			}
		}
		if (score >= threshold) {
			sitePos[numSites++]=i+halfSite;
			continue;
		}
		score=siteLength;
		for (int j=0;j<siteLength;j++) {
			if (seq[i+j]!=rvSite[j] && site[j] != 'N') {
				score--;
				if (score < threshold) break;
			}
		}
		if (score >= threshold) sitePos[numSites++]=i+halfSite;
	}
	int* allSites = new int[numSites+1];
	allSites[0]=numSites;
	for (int i=0;i<numSites;i++) {
		allSites[i+1] = sitePos[i];
	}
	delete []sitePos;
	sites->insert(allSites,chr);

	deleteSequence();
}

void ChrTags::removePETagBackground(int fragLength) {

	loadTags();
	for (int i=0;i<totalPositions;i++) {
		if (petags[i].chr1 == NULL || petags[i].chr2 == NULL) continue;

		if (strcmp(petags[i].chr1,petags[i].chr2)==0) {
			int diff = petags[i].p2-petags[i].p1;
			if (diff > -2*fragLength && diff < fragLength*2) {
				petags[i].chr1 = NULL;
				continue;
			}
		}
	}
	qsort(petags,totalPositions,sizeof(PETag),&cmpPETags);

	int last = -1;
	totalTags = 0;
	for (int i=0;i<totalPositions;i++) {
		if (petags[i].chr1 == NULL || petags[i].chr2 == NULL) continue;
		totalTags += petags[i].v;
		int j=i-1;
		if (last > -1 && petags[i].p1 == petags[j].p1 && petags[i].d1 == petags[j].d1 &&
					strcmp(petags[i].chr2,petags[j].chr2)==0 && petags[i].p2 == petags[j].p2 && 
					petags[i].d2 == petags[j].d2) {
			petags[last].v += petags[i].v;
		} else {
			last++;
			if (last == i) continue;
			petags[last].copy(&(petags[i]));
		}
	}
	totalPositions = last+1;

	dontSAVE = 0;
	print();
	freeTags();
}


void ChrTags::assignPETagsToRestrictionSites(char* site, Hashtable* siteHash,int fragLength,int mode,int midFlag,
				int removeSelfLigationFlag,int removeRestrictionEnds, double* posStrand, double* negStrand, int distLength) {
				//int removeSelfLigationFlag, double* posStrand, double* negStrand, int distLength) {
	
	loadTags();

	qsort(petags,totalPositions,sizeof(PETag),&cmpPETags);

	int* sites = (int*) siteHash->search(chr);
	int numSites = sites[0];
	int siteIndex = 1;
	double numSelfLigation = 0;
	double numEnds = 0;
	//char* nullChr = "null";
	int siteLength = strlen(site);
	int posEndDiff = siteLength/2-1;
	int negEndDiff = posEndDiff;
	if (siteLength % 2 == 0) {
		negEndDiff++;
	}
	//fprintf(stderr, "\tposEndDiff=%d\tnegEndDiff=%d\n",posEndDiff,negEndDiff);


	int distOffset = distLength/2;
	//fprintf(stderr, "distOffset=%d distLength=%d fragLength=%d\n", distOffset,distLength,fragLength);
	//This part only for tag distribution around site
	if (posStrand != NULL && negStrand != NULL && distLength > 0) {
		int minIndex = 0;
		for (int i=1;i<numSites;i++) {
			int pos = sites[i]-distOffset;
			for (int j=minIndex;j<totalPositions;j++) {
				int p1 = petags[j].p1;
				if (p1 < pos) {
					minIndex++;
				} else if (p1 >= pos+distLength) {
					break;
				} else {
					int diff = p1 - pos;
					if (petags[j].d1 == STRAND_POSITIVE) {
						posStrand[diff]+=petags[j].v;
					} else {
						negStrand[diff]+=petags[j].v;
					}
				}
			}
		}
	}

	for (int i=0;i<totalPositions;i++) {
		if (petags[i].chr1 == NULL || petags[i].chr2 == NULL) continue;

		int p1 = petags[i].p1;
		while (siteIndex <= numSites && sites[siteIndex] < p1-fragLength) {
			siteIndex++;
		}
		if (siteIndex > numSites) {
			while (i<totalPositions) {
				i++;
			}
			break;
		}
		if (p1 < sites[siteIndex]-fragLength) {
			continue;
		}
		int d1 = petags[i].d1;
		int best = fragLength+10;
		int bestIndex = -1;
		for (int j=siteIndex;j<=numSites;j++) {
			int diff = sites[j]-p1;
			if (diff > fragLength) break;
			if (d1==1) diff *= -1;
			if (removeRestrictionEnds) {
				if ((d1==0 && diff == posEndDiff) ||
							(d1==1 && diff == negEndDiff)) {
					petags[i].chr1 = NULL;
					bestIndex = -1;
					numEnds += petags[i].v;
					break;
				}
			}
			if (diff < 0) continue;
			if (diff > fragLength) continue;
			if (diff < RESTRICTION_SITE_MIN_DISTANCE) continue;
			if (diff < best) {
				best = diff;
				bestIndex = j;
			}
		}
		if (bestIndex != -1) {
			if (midFlag) {
				if (d1 == STRAND_POSITIVE) {
					if (bestIndex > 1) {
						petags[i].p1 = -1*(sites[bestIndex]+sites[bestIndex-1])/2;
					} else {
						petags[i].p1 *= -1;
					}
				} else {
					if (bestIndex < numSites-1) {
						petags[i].p1 = -1*(sites[bestIndex]+sites[bestIndex+1])/2;
					} else {
						petags[i].p1 *= -1;
					}
				}
			} else {
				petags[i].p1 *= -1;
				//petags[i].p1 = -1*sites[bestIndex];
				if (removeSelfLigationFlag && petags[i].chr1 == petags[i].chr2) {
					if (petags[i].d1 == STRAND_POSITIVE && petags[i].d2 == STRAND_NEGATIVE && bestIndex > 0) {
						if (petags[i].p2 > sites[bestIndex-1] && (petags[i].p2-sites[bestIndex-1] < fragLength)) {
							numSelfLigation += petags[i].v;
							petags[i].chr1 =NULL;
						}
					} else if (petags[i].d1 == STRAND_NEGATIVE && petags[i].d2 == STRAND_POSITIVE 
																				&& bestIndex < numSites-1) {
						if (petags[i].p2 < sites[bestIndex+1] && (sites[bestIndex+1]-petags[i].p2 < fragLength)) {
							numSelfLigation += petags[i].v;
							petags[i].chr1 =NULL;
						}
					}
							
				}
			}
		}
	}


	qsort(petags,totalPositions,sizeof(PETag),&cmp2ndPETags);
	sites = NULL;
	numSites = 0;
	siteIndex = 0;
	char* curChr= NULL;

	if (0 && posStrand != NULL && negStrand != NULL && distLength > 0) {
		int minIndex = 0;
		for (int i=1;i<numSites;i++) {
			int pos = sites[i];
			for (int j=minIndex;j<totalPositions;j++) {
				int p1 = petags[j].p2;
				if (p1 < pos-distOffset) {
					minIndex++;
				} else if (p1 >= pos+distOffset) {
					break;
				} else {
					int diff = p1 - (pos-distOffset);
					if (petags[j].d2 == STRAND_POSITIVE) {
						posStrand[diff]+=petags[j].v;
					} else {
						negStrand[diff]+=petags[j].v;
					}
				}
			}
		}
	}

	for (int i=0;i<totalPositions;i++) {

		if (petags[i].chr1 == NULL) continue;
		if (petags[i].chr2 == NULL) continue;
		if (curChr == NULL || strcmp(curChr,petags[i].chr2)!=0) {
			curChr = petags[i].chr2;
			sites = (int*) siteHash->search(petags[i].chr2);
			if (sites == NULL) {
				//fprintf(stderr, "Couldn't find sites for %s\n", petags[i].chr2);
				numSites = 0;
			} else {
				//fprintf(stderr, "\t\t\tfiltering interchromosomal interactions from %s to %s\n", chr, curChr);
				numSites = sites[0];
			}
			siteIndex = 1;
		}

		int p2 = petags[i].p2;
		while (siteIndex <= numSites && sites[siteIndex] < p2-fragLength) {
			siteIndex++;
		}
		if (siteIndex > numSites) {
			while (i<totalPositions && strcmp(petags[i].chr2,curChr)==0) {
				i++;
			}
			continue;
		}
		if (p2 < sites[siteIndex]-fragLength) {
			continue;
		}
		int d2 = petags[i].d2;
		int best = fragLength+10;
		int bestIndex = -1;
		for (int j=siteIndex;j<=numSites;j++) {
			int diff = sites[j]-p2;
			if (diff > fragLength) break;
			if (d2==1) diff *= -1;
			if (removeRestrictionEnds) {
				if ((d2==0 && diff == posEndDiff) ||
							(d2==1 && diff == negEndDiff)) {
					petags[i].chr1 = NULL;
					bestIndex = -1;
					numEnds += petags[i].v;
					break;
				}
			}
			if (diff < 0) continue;
			if (diff > fragLength) continue;
			if (diff < RESTRICTION_SITE_MIN_DISTANCE) continue;
			if (diff < best) {
				best = diff;
				bestIndex = j;
			}
		}
		if (bestIndex != -1) {
			if (midFlag) {
				if (d2 == STRAND_POSITIVE) {
					if (bestIndex > 1) {
						petags[i].p2 = -1*(sites[bestIndex]+sites[bestIndex-1])/2;
					} else {
						petags[i].p2 *= -1;
					}
				} else {
					if (bestIndex < numSites-1) {
						petags[i].p2 = -1*(sites[bestIndex]+sites[bestIndex+1])/2;
					} else {
						petags[i].p2 *= -1;
					}
				}
			} else {
				petags[i].p2 *= -1;
				//petags[i].p2 = -1*sites[bestIndex];
			}
		}
	}


	for (int i=0;i<totalPositions;i++) {
		if (petags[i].chr1 == NULL) continue;
		int site1 = 0;
		int site2 = 0;
		if (petags[i].p1 < 0) {
			site1 = 1;
			petags[i].p1 *= -1;
		}
		if (petags[i].p2 < 0) {
			site2 = 1;
			petags[i].p2 *= -1;
		}
		if (midFlag && removeSelfLigationFlag) {
			//catch religation sites and remove them
			if (petags[i].p1 == petags[i].p2 && strcmp(petags[i].chr1,petags[i].chr2)==0) {
				petags[i].chr1= NULL;
				numSelfLigation = 0;
				continue;
			}
		}
		int s = site1+site2;
		if (mode == 0) {
			if (s < 2) petags[i].chr1 = NULL;
		} else if (mode == 1) {
			if (s > 0) petags[i].chr1 = NULL;
		} else if (mode == 2) {
			if (s != 1) petags[i].chr1 = NULL;
		} else if (mode == 3) {
			if (s < 1) petags[i].chr1 = NULL;
		}
	}


	qsort(petags,totalPositions,sizeof(PETag),&cmpPETags);

	double oldNumberTags = totalTags;

	int last = -1;
	totalTags = 0;
	for (int i=0;i<totalPositions;i++) {
		if (petags[i].chr1 == NULL || petags[i].chr2 == NULL) continue;
		totalTags += petags[i].v;
		if (last > -1 && petags[i].p1 == petags[last].p1 && petags[i].d1 == petags[last].d1 &&
					strcmp(petags[i].chr2,petags[last].chr2)==0 && petags[i].p2 == petags[last].p2 && 
					petags[i].d2 == petags[last].d2) {
			petags[last].v += petags[i].v;
		} else {
			last++;
			if (last == i) continue;
			petags[last].copy(&(petags[i]));
		}
	}
	totalPositions = last+1;

	double fracRemoved =totalTags/oldNumberTags;
	double fracSelf = numSelfLigation/oldNumberTags*100.0;
	double fracEnds = numEnds/oldNumberTags*100.0;
	fprintf(stderr, "\t\t\tKept %.2lf%% of interactions from %s", 100.0*fracRemoved,chr);
	if (removeSelfLigationFlag) fprintf(stderr, " (%.2lf%% self-ligations)", fracSelf);
	if (removeRestrictionEnds) fprintf(stderr, " (%.2lf%% fragment ends)", fracEnds);
	fprintf(stderr, "\n");

	dontSAVE = 0;
	print();
	freeTags();

}

void ChrTags::autoCorrelateTags(double* sameStrand, double* diffStrand,int windowSize,
					double maxTags,double &totalCount, double* sameStrandN, double* diffStrandN, int mCflag) {
	forceSingleReadFlag = 1;
	loadTags();

	int halfWindow = (windowSize)/2;

	for (int i=0;i<totalPositions;i++) {
		double refV = (double)tags[i].v;
		for (int j=i+1;j<totalPositions;j++) {
			int diff = tags[j].p - tags[i].p;
			if (diff > halfWindow) break;
			double v = 0.0; // (double)(tags[i].v*tags[i].v);
			v = 1.0;
			if (mCflag) {
				v = fabs(tags[j].v-refV);
				//v = tags[j].v-refV;
			}
			if (tags[i].d == tags[j].d) {
				if (tags[i].d == 0) {
					sameStrand[halfWindow+diff] += v;
					sameStrandN[halfWindow+diff] += 1.0;
				} else {
					sameStrand[halfWindow-diff] += v;
					sameStrandN[halfWindow-diff] += 1.0;
				}
			} else {
				if (tags[i].d == 0) {
					diffStrand[halfWindow+diff] += v;
					diffStrandN[halfWindow+diff] += 1.0;
				} else {
					diffStrand[halfWindow-diff] += v;
					diffStrandN[halfWindow-diff] += 1.0;
				}
			}
			totalCount += v;
		}
		if (totalCount > maxTags) break;
	}
	freeTags();
	forceSingleReadFlag = 0;
}
void ChrTags::setMaxTBP(float max) {
	maxtbp = max;
}
void ChrTags::setMinTBP(float min) {
	mintbp = min;
}
void ChrTags::setTagAdjust(int dist) {
	tagAdjust = dist;
}
void ChrTags::adjustTags() {
	if (totalPositions < 1) return;
	if (adjusted) return;
//fprintf(stderr, "maxtbp=%f\tadj=%d\n",maxtbp,tagAdjust);
	if (pairedEndFlag && forceSingleReadFlag == 0) {
		if (maxtbp > FLOAT_ZERO || mintbp > FLOAT_ZERO) {
			totalTags = 0.0;
			for (int i=0;i<totalPositions;i++) {
				if (petags[i].v > maxtbp) {
					petags[i].v = maxtbp;
				}
				if (petags[i].v < mintbp) {
					petags[i].v = mintbp;
				}
				totalTags += petags[i].v;
			}
		}
		return;
	} else {
		int minPosition = tags[0].p;
		int maxPosition = tags[totalPositions-1].p;
		if (maxtbp > FLOAT_ZERO || mintbp > FLOAT_ZERO) {
			totalTags = 0;
			for (int i=0;i<totalPositions;i++) {
				if (tags[i].v > maxtbp) {
					tags[i].v = maxtbp;
				}
				if (tags[i].v < mintbp) {
					tags[i].v = mintbp;
				}
				totalTags += tags[i].v;
			}
		}
		dontSAVE = 1;
		if (tagAdjust != 0) {
			for (int i=0;i<totalPositions;i++) {
				if (tags[i].d == 0) {
					tags[i].p += tagAdjust;
				} else {
					tags[i].p -= tagAdjust;
				}
				if (tags[i].p > maxPosition) {
					tags[i].p = maxPosition;
				}
				if (tags[i].p < minPosition) {
					tags[i].p = minPosition;
				}
			}
			qsort(tags,totalPositions,sizeof(Tag),&cmpTags);
		}
	}
	adjusted = 1;
}

PeakLibrary* ChrTags::findGroSeqTranscripts(int mode, ChrTags* input, UniqMapChrs* umc, char strand,
					int tssSize, int minBodySize, int maxBodySize, 
					double threshold, double foldTranscriptStart, 
					double foldTranscriptBody, double endFold, double inputFold, double inputNorm, 
					int fragLength, int inputFragLength, double pseudoTags,int groseqMethod) {
	int maxPosition = 0;
	int initialFragLength = minBodySize/2;
	int tssFragLength = 1;
	//fprintf(stderr, "\tinitialFragLenth=%d\n", initialFragLength);
	//fprintf(stderr, "\ttssFragLenth=%d\n", tssFragLength);
	//fprintf(stderr, "\tBuilding profile...\n");
	float *profile = buildProfile(strand, 1.0,initialFragLength,maxPosition);
	//fprintf(stderr, "\tBuilding profileS...\n");
	float *profileS = buildProfile(strand, 1.0,tssFragLength,maxPosition);
	if (umc != NULL) {
		//fprintf(stderr, "\tCorrecting profile for mappability...\n");
		float * umcProfile = umc->buildUniqMapProfile(strand, initialFragLength, maxPosition);
		umc->adjustProfile(profile,umcProfile,strand,initialFragLength,maxPosition);
		delete []umcProfile;
	}


	PeakLibrary* transcripts = new PeakLibrary(1000000);

	//fprintf(stderr, "\tminBodySize = %d\n", minBodySize);
	//fprintf(stderr, "\ttssSize = %d\n", tssSize);
	double forSum = pseudoTags;
	double backSum = pseudoTags;
	double forSumS = pseudoTags;
	double backSumS = pseudoTags;
	double minLevel = pseudoTags;
	double minLevelStart = minLevel;
	if (groseqMethod == GROSEQ_METHOD_LEVEL) {
		minLevelStart = threshold*((double)fragLength)/((double)minBodySize);
		//fprintf(stderr, "minLevel=%lf\n", minLevelStart);
	}
	double minReads = threshold;
	//minReads = 20;
	Peak* trans = NULL;
	char* transName = new char[10000];
	int index = 1;


	//int inc = 1;
	if (strand == STRAND_NEGATIVE) {
		//inc = -1;
	}

	if (strand == STRAND_POSITIVE) {
		for (int i=0;i<minBodySize;i++) {
			if (i>=maxPosition) break;
			forSum+=profile[i];
		}
		for (int i=0;i<tssSize;i++) {
			if (i>=maxPosition) break;
			forSumS+=profileS[i];
		}
	} else {
		for (int i=maxPosition-1;i>maxPosition-1-minBodySize;i--) {
			if (i<0) break;
			forSum+=profile[i];
		}
		for (int i=maxPosition-1;i>maxPosition-1-tssSize;i--) {
			if (i<0) break;
			forSumS+=profileS[i];
		}
	}
	for (int x=1;x<maxPosition;x++) {
		int i= x;
		if (strand == STRAND_NEGATIVE) {
			i = maxPosition-x;
		}
			
		double bodyFold = 1.0;
		double tssFold = 0.0;

		
		if (strand == STRAND_POSITIVE) {
			if (i < maxPosition-minBodySize) {
				forSum+=profile[i+minBodySize-1];
				forSum-=profile[i-1];
			}
			if (i < minBodySize) {
				backSum += profile[i];
			} else {
				backSum += profile[i];
				backSum -= profile[i-minBodySize];
	
				if (backSum < pseudoTags) backSum= pseudoTags;
				bodyFold = forSum/backSum;
			}
			if (i < maxPosition-tssSize) {
				forSumS+=profileS[i+tssSize-1];
				forSumS-=profileS[i-1];
			}
			if (i < tssSize) {
				backSumS += profileS[i];
			} else {
				backSumS += profileS[i];
				backSumS -= profileS[i-tssSize];
	
				if (backSumS < pseudoTags) backSumS= pseudoTags;
				tssFold = forSumS-backSumS;
			}
		} else {
			if (i > minBodySize) {
				forSum+=profile[i-minBodySize+1];
				forSum-=profile[i+1];
			}
			if (i > maxPosition-minBodySize) {
				backSum += profile[i];
			} else {
				backSum += profile[i];
				backSum -= profile[i+minBodySize];
	
				if (backSum < pseudoTags) backSum= pseudoTags;
				bodyFold = forSum/backSum;
			}
			if (i >tssSize) {
				forSumS+=profileS[i-tssSize+1];
				forSumS-=profileS[i+1];
			}
			if (i > maxPosition-tssSize) {
				backSumS += profileS[i];
			} else {
				backSumS += profileS[i];
				backSumS -= profileS[i+tssSize];
	
				if (backSumS < pseudoTags) backSumS= pseudoTags;
				tssFold = forSumS-backSumS;
			}
		}
	
	
		if (profile[i] > minLevel) {
			if ((groseqMethod == GROSEQ_METHOD_LEVEL && profile[i] > minLevelStart) ||
						(groseqMethod == GROSEQ_METHOD_FOLD && 
					  	 (bodyFold > foldTranscriptBody && tssFold > foldTranscriptStart))) {
				if (trans == NULL) {
					sprintf(transName, "%s-%d-%d",chr,index++,strand);
					trans = new Peak(transName,transName, chr, i,i,i,strand,tssFold,0.0,NULL,0,0);
				} else if ((strand == STRAND_POSITIVE && i-trans->start < minBodySize)
							|| (strand == STRAND_NEGATIVE && trans->end-i < minBodySize)) {
					if (tssFold > trans->v) {
					//if (groseqMethod != GROSEQ_METHOD_LEVEL && tssFold > trans->v) {
						if (strand == STRAND_POSITIVE) trans->start = i;
						if (strand == STRAND_NEGATIVE) trans->end = i;
						trans->v = tssFold;
					}
				} else if (groseqMethod == GROSEQ_METHOD_LEVEL) {
		
				} else {
					if (strand == STRAND_POSITIVE) trans->end = i;
					if (strand == STRAND_NEGATIVE) trans->start = i;
					double density = trimGroSeqTranscript(profileS,trans,strand,fragLength,endFold,minBodySize);
			//	if (density*minBodySize < minReads) {
					if (density*(trans->end-trans->start) < minReads) {
						delete trans;
					} else {
						transcripts->addPeak(trans);
					}
					//trans->print(stdout);
					//delete trans;
					trans = NULL;
					sprintf(transName, "%s-%d-%d",chr,index++,strand);
					trans = new Peak(transName,transName, chr, i,i,i,strand,tssFold,0.0,NULL,0,0);
				}
			}
		} else {
			if (trans != NULL) {
				if (strand == STRAND_POSITIVE) trans->end = i;
				if (strand == STRAND_NEGATIVE) trans->start = i;
				double density = trimGroSeqTranscript(profileS,trans,strand,fragLength,endFold,minBodySize);
			//	if (density*minBodySize < minReads) {
				if (density*(trans->end-trans->start) < minReads) {
					delete trans;
				} else {
					transcripts->addPeak(trans);
				}
				//trans->print(stdout);
				//delete trans;
				trans = NULL;
			}
		}
	}
	//transcripts->print(stdout);
	//if (strand == STRAND_NEGATIVE) exit(0);

	delete []profile;
	delete []profileS;
	delete []transName;

	return transcripts;
}

double ChrTags::trimGroSeqTranscript(float* profile, Peak* trans,int strand, int fragLength,double endFold, int minBodySize) {

	double sum = 0.0;	
	double N = 0.0;
	double initValue = 0.0;
	if (strand == STRAND_POSITIVE) {
		int newEnd = trans->start;
		for (int i=trans->start;i<trans->end;i++) {
			if (i>=trans->start+minBodySize) break;
			sum+=profile[i];
			if (profile[i] > 0.0) newEnd = i+fragLength;
			N += 1.0;
		}
		initValue = sum;
		if (N > 1.0) {
			initValue/=N;
		} else {
			trans->end = newEnd;
			trans->focusRatio = 0;
			return 0.0;
		}
		for (int i=trans->start+minBodySize;i<trans->end;i++) {
			sum+=profile[i];
			sum-=profile[i-minBodySize];
			if (sum/N > initValue/endFold) {
				if (profile[i] > 0.0) newEnd = i+fragLength;
			}
		}
		trans->end = newEnd;
	} else {
		int newStart = trans->end;
		for (int i=trans->end;i>trans->start;i--) {
			if (i<=trans->end-minBodySize) break;
			sum+=profile[i];
			if (profile[i] > 0.0) newStart = i-fragLength;
			N += 1.0;
		}
		initValue = sum;
		if (N > 1.0) {
			initValue/=N;
		} else {
			trans->start = newStart;
			trans->focusRatio = 0;
			return 0.0;
		}
		for (int i=trans->end-minBodySize;i>trans->start;i--) {
			sum+=profile[i];
			sum-=profile[i+minBodySize];
			if (sum/N > initValue/endFold) {
				if (profile[i] > 0.0) newStart = i-fragLength;
			}
		}
		trans->start = newStart;
	}

	trans->focusRatio = (double)(trans->end-trans->start);
	trans->v = initValue;
	return initValue;
}


/*
	int numDepths = 1;
	int maxNumRes = numDepths*2;
	int* depths = new int[numDepths];
	depths[0] = 50;
	//depths[1] = 250;
	//depths[2] = 2500;
	float **smooth = new float*[maxNumRes+1];
	for (int i=0;i<maxNumRes;i++) {
		smooth[i] = new float[maxPosition+1];
	}
	smooth[2] = new float[maxPosition+1];
	fprintf(stderr, "\tBuilding tag coverages at different resolutions\n");

	for (int k=0;k<numDepths;k++) {	
		float thresh = depths[k];
		fprintf(stderr, "\tPreping %f\n", thresh);
		double lastSum = -1;
		double lastSumN = -1;
		int lastPos = 0;
		int baseIndex = k*2;
		for (int i=0;i<maxPosition;i++) {
			//if (i % 10000000 == 0) fprintf(stderr, "%d\n",i); 
			double sum = lastSum;
			double sumN = lastSumN;
			if (lastSum < 0.0) {
				sum=0.0;
				sumN=0.0;
			} else {
				sum -= profile[i-1];
				sumN -= umcProfile[i-1];
			}
			
			if (sum < thresh) {
				for (int j=lastPos;j<maxPosition;j++) {
					sum += profile[j];
					sumN += umcProfile[j];
					lastPos = j+1;
					if (sum > thresh) break;
				}
			}
				
			//int len = (lastPos-i);
			double len = sumN;
			double d = sum;
			if (sum < 0.001) d = 1.0;
			if (len < 1.5) {
				len = thresh/d;
			}
			//float v = (float)len;
			float v = log(1.0/len)/log(2.0);
			smooth[baseIndex][i] = v;
			lastSum = sum;
			lastSumN = sumN;
		}
		lastSum = -1;
		lastPos = maxPosition-1;
		lastSumN = -1;
		for (int i=maxPosition-1;i>=0;i--) {
			//if (i % 10000000 == 0) fprintf(stderr, "%d\n",i); 
			double sum = lastSum;
			double sumN = lastSumN;
			if (lastSum < 0.0) {
				sum=0.0;
				sumN=0.0;
			} else {
				sum -= profile[i+1];
				sumN -= umcProfile[i+1];
			}
			
			if (sum < thresh) {
				for (int j=lastPos;j>=0;j--) {
					sum += profile[j];
					sumN += umcProfile[j];
					lastPos = j-1;
					if (sum > thresh) break;
				}
			}
					
				
			//int len = (i-lastPos);
			double len = sumN;
			double d = sum;
			if (sum < 0.001) {
				d = 1.0;
			}
			if (len < 1.5) {
				len = thresh/d;
			}
			d = len;
			float v = log(1.0/(float)d)/log(2.0);
			smooth[baseIndex+1][i] = v;
			lastSum = sum;
			lastSumN = sumN;
		}
	}



	if (0) {

		int gstart = 74269921;
		int gend = 74294520;

		fprintf(stdout, "Coordinate\tprofile");
		fprintf(stdout , "\n");
		for (int i=gstart;i<=gend;i++) {
			fprintf(stdout, "%d\t%.3f", i, profile[i]);
			for (int j=0;j<maxNumRes;j++){ 
				fprintf(stdout, "\t%.3f", smooth[j][i]);
			}
			fprintf(stdout, "\n");

		}
		exit(0);
	}

	fprintf(stderr, "\n\tFinding Transcripts\n");
	char* transName = new char[10000];


	int windowSize = 1000;
	double forWindowSum = 1.0;	
	double lagWindowSum = 1.0;	
	for (int i=minBodySize;i<minBodySize+windowSize;i++) {
		if (i >= maxPosition) break;
		forWindowSum += smooth[0][i];
	}
	for (int i=0;i<maxPosition;i++) {
		smooth[2][i] = (smooth[0][i]+smooth[1][i])/2.0;
	}


	double logFoldTranscriptBody = log(foldTranscriptBody)/log(2.0);
	double logEndFold = log(endFold)/log(2.0);
	int index = 1;
	int minGiveUp = -1;
	Peak* trans = NULL;
	double cutoff = log(1.0/(double)maxBodySize)/log(2.0);
	double transLevel = 0.0;
	double transLevelN = 0.0;
	fprintf(stderr, "\tCutoff = %lf\n", cutoff);
	fprintf(stderr, "\tMinBodySize = %d\n", minBodySize);
	fprintf(stderr, "\tMaxBodySize = %d\n", maxBodySize);
	fprintf(stderr, "\tWindowSize = %d\n", windowSize);
	for (int i=0;i<maxPosition;i++) {
		//if (i % 100000 == 0) fprintf(stderr, "%d\n",i); 

		if (i < maxPosition-minBodySize-windowSize) {
			forWindowSum += smooth[2][i+windowSize+minBodySize];
			forWindowSum -= smooth[2][i+minBodySize];
		}
		if (i < windowSize+1) {
			lagWindowSum += smooth[2][i-1];
		} else {
			lagWindowSum += smooth[2][i-1];
			lagWindowSum -= smooth[2][i-windowSize-1];
		}
		double ratio = (forWindowSum-lagWindowSum)/(double)windowSize;

	
		int base = 0;	
		double level = smooth[base*2][i];
		double level2 = smooth[2][i];
		int dist = (int)(1.0/pow(2,level));
		while (dist < minBodySize && base < numDepths-1) {
			base++;
			level = smooth[base*2][i];
			dist = (int)(1.0/pow(2,level));
		}
		if (dist > minBodySize && base > 0) {
			base--;
			level = smooth[base*2][i];
			dist = (int)(1.0/pow(2,level));
		}
			
		double initScore = smooth[base*2][i]-smooth[base*2+1][i];

			

		if (level > cutoff) {
			if (ratio > logFoldTranscriptBody) {
				if (trans == NULL) {
					sprintf(transName, "%s-%d-%d-abovethresh",chr,index++,strand);
					trans = new Peak(transName,transName, chr, i,i+1,i,strand,initScore,initScore,NULL,0,0);
					transLevel = 0.0;
					transLevelN = 0.0;
					minGiveUp = (int)(1.0/pow(2,level))+i;
					if (minGiveUp < minBodySize) minGiveUp = minBodySize;
				} else {
					if (trans->v < initScore) {
						if (i-trans->start < minBodySize) {
							trans->v = initScore;
							trans->focusRatio = initScore;
							trans->start = i;
							trans->end = i+1;
							minGiveUp = (int)(1.0/pow(2,level))+i;
							if (minGiveUp < minBodySize) minGiveUp = minBodySize;
						} else {
							//trans->end = i;
							trans->print(stdout);
							delete trans;
							sprintf(transName, "%s-%d-%d-InitScore-old%lf-vs-%lf-%f-%f-%d",chr,index++,strand,trans->v,initScore,smooth[0][i],smooth[1][i],i);
							trans = new Peak(transName,transName, chr, i,i+1,i,strand,initScore,0.0,NULL,0,0);
							transLevel = 0.0;
							transLevelN = 0.0;
							minGiveUp = (int)(1.0/pow(2,level))+i;
							if (minGiveUp < minBodySize) minGiveUp = minBodySize;
						}
					} else if (i-trans->start > windowSize*2
									&& ratio > logFoldTranscriptBody) {
//fprintf(st	derr, "%d %lf\t%lf\n", i, ratio, logFoldTranscriptBody);
						//trans->end = i;
						trans->print(stdout);
						delete trans;
						sprintf(transName, "%s-%d-%d-BodyFold",chr,index++,strand);
						trans = new Peak(transName,transName, chr, i,i+1,i,strand,initScore,0.0,NULL,0,0);
						transLevel = 0.0;
						transLevelN = 0.0;
						minGiveUp = (int)(1.0/pow(2,level))+i;
						if (minGiveUp < minBodySize) minGiveUp = minBodySize;
					}
				}
			}
			if (trans != NULL && transLevelN > 0.5) {
				double avgLevel = transLevel/transLevelN;
				if (level > avgLevel) {
					trans->end=i;
				} else if (level < avgLevel-logEndFold) {
					trans->print(stdout);
					delete trans;
					trans = NULL;
				}
			}
		} else {
			if (level2 < cutoff && trans != NULL && i>minGiveUp) {
				trans->print(stdout);
				delete trans;
				trans = NULL;
			}
		}
		if (trans != NULL && trans->focusRatio > initScore) {
			trans->focusRatio = initScore;
			trans->end = i;
		}
		if (trans != NULL) {
			transLevel += smooth[2][i];
			transLevelN += 1.0;
		}
	}	

	exit(0);
	return NULL;
*/

PeakLibrary* ChrTags::findGroSeqRegions(int mode, ChrTags* input, UniqMapChrs* umc, char strand,
					int tssSize, int minBodySize, int maxBodySize, 
					double threshold, double foldTranscriptStart, 
					double foldTranscriptBody, double endFold, double inputFold, double inputNorm, 
					int fragLength, int inputFragLength, double pseudoTags) {

	int maxPosition = 0;
	float *profile = buildProfile(strand, 1.0,fragLength,maxPosition);

	double avgCoverage = 0;
	for (int i=0;i<maxPosition;i++) {
		avgCoverage+=profile[i];
	}
	avgCoverage /=(double)maxPosition;
	//fprintf(stderr, "Avg Coverage: %lf\n", avgCoverage);


	int maxInputPosition;
	float *inputProfile = NULL;
	if (input != NULL) {
		inputProfile = buildProfile(strand, inputNorm, inputFragLength, maxInputPosition);	

		double avgInputCoverage = 0;
		for (int i=0;i<maxInputPosition;i++) {
			avgInputCoverage+=inputProfile[i];
		}
		avgInputCoverage /=(double)maxInputPosition;
		fprintf(stderr, "Avg Input Coverage: %lf\n", avgInputCoverage);

		float pseudo = (float) (avgCoverage/2.0);
		float pseudoInput = (float) (avgInputCoverage/2.0);
		for (int i=0;i<maxPosition;i++) {
			float v = profile[i];
			if (v < pseudo) v = pseudo;
			float vi = pseudoInput;
			if (i<maxInputPosition) vi = inputProfile[i];	
			if (vi < pseudoInput) vi = pseudoInput;
			float ratio = v/vi;
			profile[i] = ratio;
		}
		delete []inputProfile;
	}

	//fprintf(stdout, "track type=bedGraph alwaysZero=on name=\"transcripts\"\n");
	if (umc != NULL) {

		float minMapLevel = 1.0/((float)fragLength);
		//fprintf(stderr, "%d\n", fragLength);
		float * umcProfile = umc->buildUniqMapProfile(strand, fragLength, maxPosition);

		float maplevel = 0;
		float historicAvg = 0;
		float historicN = 0;

		for (int i=1;i<maxPosition;i++) {
			int index = i;
			if (strand == STRAND_NEGATIVE) {
				index = maxPosition-i-1;
			}

			maplevel = umcProfile[index];
			if (maplevel < 0) {
				//maplevel = minMapLevel;
			}
			if (historicN > 0.5) {
				if (maplevel < 1.0-minMapLevel) {
					profile[index] = profile[index] + (historicAvg/historicN)*(1-maplevel);
				}
			}
			if (i<=fragLength) {
				historicN++;
				historicAvg += profile[index];
			} else {
				historicAvg += profile[index];
				if (strand == STRAND_NEGATIVE) {
					historicAvg -= profile[index+fragLength];
				} else {
					historicAvg -= profile[index-fragLength];
				}
			}
			if (historicAvg < 0) historicAvg = 0;
		}

		delete []umcProfile;
	}

	//fprintf(stderr, "Done with uniqmap\n");
	
	//fprintf(stdout, "track name=\"transcripts\"\n");
	PeakLibrary* transcripts = new PeakLibrary();

	int numBodySizes = (maxBodySize-minBodySize)/tssSize+1;
	int* bodySizes = new int[numBodySizes];
	int* dirBodySizes = new int[numBodySizes];
	bodySizes[0] = minBodySize;
	for (int i=1;i<numBodySizes;i++) {
		bodySizes[i] = (int) (bodySizes[i-1]*1.5);
		//fprintf(stderr, "\t\t%d\t%d\n", i, bodySizes[i]);
		//bodySizes[i] = minBodySize+i*tssSize;
		dirBodySizes[i] = bodySizes[i];
		if (bodySizes[i] > maxBodySize) {
			maxBodySize = bodySizes[i];
			numBodySizes = i+1;
			break;
		}
	}
	//fprintf(stderr, "\t# bodySizes = %d\n",numBodySizes);
	for (int i=0;i<numBodySizes;i++) {
		dirBodySizes[i] = bodySizes[i];
	}

	double* transLevel = new double[numBodySizes];
	double* transN = new double[numBodySizes];
	double* preLevel = new double[numBodySizes];
	double* preN = new double[numBodySizes];
	double* tssLevel = new double[numBodySizes];
	int* tssDeadZones = new int[numBodySizes];
	double* pseudoCounts = new double[numBodySizes];
	for (int i=0;i<numBodySizes;i++) {
		transLevel[i]=0.0;
		transN[i]=0.0;
		preLevel[i]=0.0;
		preN[i]=0.0;
		tssDeadZones[i] = -1000;
		pseudoCounts[i] = pseudoTags*(double)fragLength/(double)bodySizes[i];
	}

	int lastStart = 0;
	double lastStartLevel = 0;
	double transcriptFlag = 0;
	//int tssDeadZone = -100;
	int curBodySize = -1;

	int dirMaxBodySize = maxBodySize;
	int dirTssSize = tssSize;

	//int dirFactor = 1;
	if (strand== STRAND_NEGATIVE) {
		lastStart = 2000000000;
		//tssDeadZone = 2000000000;
		dirTssSize = -1*tssSize;
		//dirFactor = -1;
		dirMaxBodySize = -1*maxBodySize;
		for (int i=0;i<numBodySizes;i++) {
			dirBodySizes[i] *= -1;
			tssDeadZones[i] = 2000000000;
		}
	}

	char* transcriptName = new char[10000];
	double curTransTotal = 0.0;

	for (int i=0;i<maxPosition;i++) {

		//if (i%10000000 == 0) fprintf(stderr, "\t\t%d\n", i);

		for (int j=0;j<numBodySizes;j++) {
			tssLevel[j] = 0.0;
		}

		int gIndex = i;
		if (strand == STRAND_NEGATIVE) {
			gIndex = maxPosition-1-i;
		}
		int tss = gIndex-dirTssSize-dirMaxBodySize;

		for (int j=0;j<numBodySizes;j++) {
			int p = i-maxBodySize+bodySizes[j];
			int index = p;
			if (strand == STRAND_NEGATIVE) {
				index = maxPosition-1-p;
			}

			if (p < 0) {
			} else if (p < bodySizes[j]) {
				transN[j] += 1.0;
				transLevel[j] += profile[index];
			} else {
				transLevel[j] += profile[index];
				transLevel[j] -= profile[index-dirBodySizes[j]];
			}

			if (p > bodySizes[j]+tssSize  && p <= bodySizes[j]*2+tssSize) {
				preN[j] += 1.0;
				preLevel[j] += profile[index-dirTssSize-dirBodySizes[j]];
			} else if (p > bodySizes[j]*2+tssSize) {
				preLevel[j] += profile[index-dirTssSize-dirBodySizes[j]];
				preLevel[j] -= profile[index-dirTssSize-dirBodySizes[j]*2];
			}

			if (p > maxBodySize*2+tssSize) {

				//only in this case can we consider transcripts

				double trans = transLevel[j]/transN[j];
				double pre = preLevel[j]/preN[j];

				double inductionFold = profile[tss]/(pre+pseudoCounts[j]);
				double trailingFold = (trans+pseudoCounts[j])/(pre+pseudoCounts[j]);
				tssLevel[j]=-1*trailingFold;

				double tagCounts = transLevel[j]/fragLength;



					//high enough to be considered at this bodySize...

				if (inductionFold > foldTranscriptStart
								 && inductionFold > trailingFold) {
						//	) {

					if ((strand == STRAND_NEGATIVE && tss < tssDeadZones[j])
										|| (strand != STRAND_NEGATIVE && tss > tssDeadZones[j])) {
						if (trailingFold > foldTranscriptBody) {
							if (tagCounts >= threshold || j+1==numBodySizes) {
								tssLevel[j]=trans;
							} else {
							}
						} else {
							//spike/artifact...
						}
					} else {
						//To close to previous TSS
						if (transcriptFlag == 1 && (tagCounts >= threshold || j+1==numBodySizes)) {
							tssDeadZones[j] = tss+dirTssSize+dirBodySizes[j];
						}
					}
				}

				if (transcriptFlag == 1 && curBodySize == j) {
												
					if ((trans < lastStartLevel/endFold)) {
					//} else if (transcriptFlag == 1 && (curLevel < lastStartLevel/endFold)) {
						transcriptFlag = 0;
						//transcripts->addPeak(NULL,chr,lastStart,tss,(tss-lastStart)/2,
						//				strand,lastLevel,0,NULL,-1,0);
						double vv = curTransTotal / (double)(tss-lastStart);
						double transLength = tss-lastStart;
						if (strand==STRAND_NEGATIVE) {
							sprintf(transcriptName, "%s:%d-%d,-",chr,tss,lastStart);
							transcripts->addPeak(transcriptName,chr,tss,lastStart,lastStart,strand,
													-1*vv,-1*transLength,NULL,bodySizes[curBodySize],0);
							//fprintf(stdout, "%s\t%d\t%d\t-\n", chr, tss,lastStart);
						} else {
							sprintf(transcriptName, "%s:%d-%d,+",chr,lastStart,tss);
							transcripts->addPeak(transcriptName,chr,lastStart,tss,lastStart,strand,
													vv,transLength,NULL,bodySizes[curBodySize],0);
							//fprintf(stdout, "%s\t%d\t%d\t+\n", chr,lastStart,tss);
						}
					} else {
						curTransTotal += trans;
					}
				}
			}
		}
		int ok2check = 1;
/*		if (transcriptFlag) {
			ok2check = 0;
			for (int j=0;j<numBodySizes;j++) {
				if (curBodySize == j && tssLevel[j] > 0.001) {
					ok2check=1;
				}
			}
		}*/
		for (int j=0;j<numBodySizes;j++) {
			if (tssLevel[j] > 0.001) break;
			//if (tssLevel[j] > -1*((foldTranscriptBody-1)/2+1)) {
			if (tssLevel[j] > -1*(foldTranscriptBody)) {
				ok2check=0;
				break;
			}
		}

		if (ok2check) {
			for (int j=0;j<numBodySizes;j++) {
				if (tssLevel[j] > 0) {
					if (transcriptFlag == 1) {
						double vv = curTransTotal / (double)(tss-lastStart);
						double transLength = tss-lastStart;
						if (strand==STRAND_NEGATIVE) {
							sprintf(transcriptName, "%s:%d-%d,-",chr,tss,lastStart);
							transcripts->addPeak(transcriptName,chr,tss,lastStart,lastStart,strand,
												-1*vv,-1*transLength,NULL,bodySizes[curBodySize],0);
							//fprintf(stdout, "%s\t%d\t%d\t-\n", chr, tss,lastStart);
						} else {
							sprintf(transcriptName, "%s:%d-%d,+",chr,lastStart,tss);
							transcripts->addPeak(transcriptName,chr,lastStart,tss,lastStart,strand,
												vv,transLength,NULL,bodySizes[curBodySize],0);
							//fprintf(stdout, "%s\t%d\t%d\t+\n", chr, lastStart,tss);
						}
					}
					lastStart = tss;
					lastStartLevel = tssLevel[j];
					transcriptFlag = 1;
					//tssDeadZone = tss+dirTssSize/2;
					for (int k=0;k<numBodySizes;k++) {
						tssDeadZones[k] = tss+dirTssSize+dirBodySizes[k];
					}
					curTransTotal = 0.0;
					curBodySize = j;
					//fprintf(stdout, "%s\t%d\t%d\t%.1f\n", chr, tss,tss+100,inductionFold);
					//fprintf(stdout, "%s\t%d\t%d\t%.1f\n", chr, tss,tss+100,trailingFold);
					break;
				}
			}
		}
	}
	delete []profile;
	transcripts->sortChr();
	return transcripts;
}

float* ChrTags::buildProfile(char strand, double norm, int fraglen, int &maxPosition) {

	maxPosition = tags[totalPositions-1].p + fraglen+1000;
	float *levels = new float[maxPosition];
	for (int i=0;i<maxPosition;i++) levels[i] = 0;

	for (int i=0;i<totalPositions;i++) {
		if (strand != STRAND_BOTH && strand != tags[i].d) continue;
		if (tags[i].d == STRAND_POSITIVE) {
			for (int j=tags[i].p;j<tags[i].p+fraglen;j++) {
				if (j >= maxPosition) break;
				levels[j] += tags[i].v*norm;
			}
		} else {
			for (int j=tags[i].p;j>tags[i].p-fraglen;j--) {
				if (j < 1) break;
				levels[j] += tags[i].v*norm;
			}
		}
	}
	return levels;
}

void ChrTags::findmCPeaks(PeakLibrary* peaks, int peakSize, char strand,int mCflag,
									double mCthresh,int minNumC, ChrTags* input) {
	loadTags();
	if (input != NULL) {
		input->loadTags();
	}

	int minC = minNumC;

	int resFactor = 20;


	//Initially we are going to tile the chromosome and look for regions of interesting methylation
	//Regions that have enough cytosines are kept for futher analysis
	int numMaxTestPeaks = tags[totalPositions-1].p/peakSize*resFactor+1;
	PeakmC** testPeaks = new PeakmC*[numMaxTestPeaks];
	int numTestPeaks=0;
	
	int index = 0;
	for (int i=0;i<numMaxTestPeaks;i++) {
		PeakmC* tp = new PeakmC();
		int start = i*peakSize/((double)resFactor);
		int end = start+peakSize;
		tp->start = start;
		tp->end = end;
		tp->gstart = end;
		tp->gend = start;
		int sIndex =-1;
		int eIndex =-1;
		for (int j=index;j<totalPositions;j++) {
			if (tags[j].p < start) {	
				if (j==index) index++;
				continue;	
			} else if (tags[j].p >= end) {
				break;
			} else {
				if (sIndex ==-1) sIndex = j;
				eIndex = j;
				tp->numC +=1;
				tp->mCavg += tags[j].v;
				if (tags[j].v < mCthresh) {
					if (tags[j].p < tp->gstart) {
						tp->gstart=tags[j].p;
					}
					if (tags[j].p > tp->gend) {
						tp->gend = tags[j].p;
					}
				}
			}
		}
		tp->calcAvg();
		//fprintf(stderr, "tp->numC=%d\n", tp->numC);
		if (eIndex != -1 && tp->numC >= minC) {
			for (int j=sIndex;j<=eIndex;j++) {
				tp->mCstd += (tp->mCavg-tags[j].v)*(tp->mCavg-tags[j].v);
			}
			tp->mCstd /= (double)tp->numC;
			tp->mCstd = sqrt(tp->mCstd);
			testPeaks[numTestPeaks++] = tp;
			//fprintf(stdout, "%s\t%d\t%d\tRegion-%lf-%lf-%d\t1\t+\n", chr,start, end, tp->mCavg, tp->mCstd, tp->numC);
		} else {
			delete tp;
		}
	}


	//We then sort the regions by their degree of methylation
	PeakmC** sortPeaks = new PeakmC*[numMaxTestPeaks];
	for (int i=0;i<numTestPeaks;i++) sortPeaks[i] = testPeaks[i];

	qsort(sortPeaks,numTestPeaks,sizeof(PeakmC*),&cmpPeakmC);
	int midPoint = numTestPeaks/2;
	int mid4 = (int)numTestPeaks*(4.0/10.0);
	int mid6 = (int)numTestPeaks*(6.0/10.0);
	double avgStd = 0.0;
	double mCmedian = 0;
	if (numTestPeaks > 0) {
		mCmedian = sortPeaks[midPoint]->mCavg;
	}
	for (int i=mid4;i<mid6;i++) {
		avgStd += sortPeaks[i]->mCstd;
	}
	if (mid6-mid4 > 0) {
		avgStd /= (double)(mid6-mid4);
	}
	fprintf(stderr, "\t\t%s\tmedian=%lf\t+/- %lf (avg stdev of median quintile)", chr, mCmedian,avgStd);



	//stitch together regions that are below the threshold
	int lastStart = -1;
	int lastEnd = -1;
	int glastStart = -1;
	int glastEnd = -1;
	int lastN = 0;
	double lastmC = 0;
	int chrNumPeaks = 0;
	int badPeak = 1;
	for (int i=0;i<numTestPeaks;i++) {
		
		if ((mCflag == PEAKS_FIND_UNMETHYLC && testPeaks[i]->mCavg < mCthresh)
				|| (mCflag == PEAKS_FIND_METHYLC && testPeaks[i]->mCavg > mCthresh)) {
		//if (testPeaks[i]->mCavg+testPeaks[i]->mCstd < mCmedian-avgStd) {
			if (lastN > 0) {
				if (lastEnd+1>testPeaks[i]->start || badPeak == 0) {
					lastEnd = testPeaks[i]->end;
					glastEnd = testPeaks[i]->gend;
					lastN++;
					lastmC += testPeaks[i]->mCavg;
				} else {
					if (lastStart < 0) lastStart = 0;
					if (glastStart < 0) glastStart = 0;
					if (lastEnd > tags[totalPositions-1].p) lastEnd = tags[totalPositions-1].p;
					if (glastEnd > tags[totalPositions-1].p) glastEnd = tags[totalPositions-1].p;

					peaks->addPeak(NULL,chr,glastStart,glastEnd,(glastStart+glastEnd)/2,STRAND_POSITIVE,lastmC/lastN,lastN,NULL,0,0);
					chrNumPeaks++;
					//fprintf(stdout, "%s\t%d\t%d\t%s-Region%d\t%lf\t+\n",chr,lastStart,lastEnd,chr,regionID++,lastmC/lastN);
					lastN = 1;
					lastStart = testPeaks[i]->start;
					lastEnd = testPeaks[i]->end;
					glastStart = testPeaks[i]->gstart;
					glastEnd = testPeaks[i]->gend;
					lastmC = testPeaks[i]->mCavg;
					badPeak = 0;
				}
			} else {
				lastN = 1;
				lastStart = testPeaks[i]->start;
				lastEnd = testPeaks[i]->end;
				glastStart = testPeaks[i]->gstart;
				glastEnd = testPeaks[i]->gend;
				lastmC = testPeaks[i]->mCavg;
				badPeak = 0;
			}
		} else {
			badPeak =1;
		}
	}
	if (lastN != -1 && lastStart != -1) {
		if (lastStart < 0) lastStart = 0;
		if (glastStart < 0) glastStart = 0;
		if (lastEnd > tags[totalPositions-1].p) lastEnd = tags[totalPositions-1].p;
		if (glastEnd > tags[totalPositions-1].p) glastEnd = tags[totalPositions-1].p;

		peaks->addPeak(NULL,chr,glastStart,glastEnd,(glastStart+glastEnd)/2,STRAND_POSITIVE,lastmC/lastN,lastN,NULL,0,0);
		chrNumPeaks++;
		//fprintf(stdout, "%s\t%d\t%d\t%s-Region%d\t%lf\t+\n",chr,lastStart,lastEnd,chr,regionID++,lastmC/lastN);
	}
	fprintf(stderr, "\t%d regions\n", chrNumPeaks);
			

	for (int i=0;i<numTestPeaks;i++) {
		if (testPeaks[i] != NULL) delete testPeaks[i];
		testPeaks[i]=NULL;
	}
	delete []testPeaks;
	delete []sortPeaks;


	if (input != NULL) {
		input->freeTags();
	}
	freeTags();
}
PeakmC::PeakmC() {
	numC = 0;
	mCavg = 0.0;
	mCstd = 0.0;
	start = -1;
	end = -1;
}
void PeakmC::calcAvg() {
	if (numC > 0) {
		mCavg /= numC;
	}
}
PeakmC::~PeakmC() {
}

int cmpPeakmC(const void* a, const void* b) {
	double mCa = (*(PeakmC**)a)->mCavg;
	double mCb = (*(PeakmC**)b)->mCavg;
	if (mCa < mCb) return -1;
	if (mCa > mCb) return 1;
	return 0;
}

// negative value for minDist is used to indicate region based peak finding
void ChrTags::findPutativePeaks(PeakLibrary* putativePeaks, int peakSize, int minDist,char strand, double minCount) {

	doubleIndex* tagCount = new doubleIndex[totalPositions];
	int* centers = new int[totalPositions];
	char* mask = new char[totalPositions];
	char* pname = new char[100];
	int pid = 1;

	int regionFlag = 0;
	if (minDist < 0) {
		minDist *= -1;
		regionFlag = 1;
	}

	int lastPosIndex = -1;
	int lastNegIndex = -1;
	int lastPosFinish = -1;
	int lastNegFinish = -1;
	int lastPosPosition = -1;
	int maxEND = 0;
	if (totalPositions > 0) {
		maxEND = tags[totalPositions-1].p;
	}

	fprintf(stderr, "\t\tFinding peaks on %s (minCount=%.1lf), total tags positions = %d\n",
																chr,minCount,totalPositions);

	for (int i=0;i<totalPositions;i++) {
		float v = tags[i].v;
		char d = tags[i].d;
		tagCount[i].index = i;
		int p = tags[i].p;

		if (strand==STRAND_BOTH) {
			int start = i+1;
			if (lastPosFinish >= i) {
				if (lastPosPosition < p) {
					tagCount[i].v = tagCount[i-1].v-tags[i-1].v;
					tagCount[i].vp = tagCount[i-1].vp-tags[i-1].v*(double)((tags[i-1].p)-p);
					start = lastPosIndex+1;
				} else {
					//some tags could be at the same location...
					tagCount[i].v = tagCount[i-1].v;
					tagCount[i].vp = tagCount[i-1].vp;
					continue;
				}
			} else {
				tagCount[i].v = v;
				//tagCount[i].vp = v*(double)(p-p); i.e. 0
				tagCount[i].vp = 0;
			}
			for (int j=start;j<totalPositions;j++) {
				int diff = tags[j].p-p;
				if (diff > peakSize) break;
				tagCount[i].v += tags[j].v;
				tagCount[i].vp += tags[j].v*(double)((tags[j].p)-p);
				lastPosFinish = i;
			}
			lastPosPosition = p;
		} else if (strand==STRAND_SEPARATE) {
			int start = i+1;
			if (d==0) {
				if (lastPosFinish >= i) {
					if (tags[lastPosFinish].p == p) {
						tagCount[i].v = tagCount[lastPosIndex].v;
						tagCount[i].vp = tagCount[lastPosIndex].vp;
						lastPosIndex = i;
						continue;
					} else {
						tagCount[i].v = tagCount[lastPosIndex].v-tags[lastPosIndex].v;
						tagCount[i].vp = tagCount[lastPosIndex].vp 
										- tags[lastPosIndex].v*(double)(tags[lastPosIndex].p - p);
						start = lastPosFinish+1;
					}
				} else {
					tagCount[i].v = v;
					tagCount[i].vp = 0;
					//tagCount[i].vp = v*(double)(p-p);
				}
				lastPosIndex = i;
			} else {
				if (lastNegFinish >= i) {
					if (tags[lastNegFinish].p == p) {
						tagCount[i].v = tagCount[lastNegIndex].v;
						tagCount[i].vp = tagCount[lastNegIndex].vp;
						lastNegIndex = i;
						continue;
					} else {
						tagCount[i].v = tagCount[lastNegIndex].v-tags[lastNegIndex].v;
						tagCount[i].vp = tagCount[lastNegIndex].vp 
										- tags[lastNegIndex].v*(double)(tags[lastNegIndex].p-p);
						start = lastNegFinish+1;
					}
				} else {
					tagCount[i].v = v;
					tagCount[i].vp = 0;
					//tagCount[i].vp = v*(double)(p-p); i.e. 0
				}
				lastNegIndex = i;
			}
			for (int j=start;j<totalPositions;j++) {
				int diff = tags[j].p-p;
				if (diff > peakSize) break;
				if (d != tags[j].d) continue;
				tagCount[i].v += tags[j].v;
				tagCount[i].vp += tags[j].v*(double)(tags[j].p-p);
				if (d==0) {
					lastPosFinish = i;
				} else {
					lastNegFinish = i;
				}
			}
		}
	}

	for (int i=0;i<totalPositions;i++) {
		mask[i] = 0;
		if (tagCount[i].v < 0.000001) tagCount[i].v = 1.0;
		tagCount[i].position = tags[i].p+(int)floor(tagCount[i].vp/tagCount[i].v);
		if (tagCount[i].position < 0) {
			fprintf(stderr, "%d\t%d\t%d\t%lf\t%lf\n", i,tags[i].p, tagCount[i].position,tagCount[i].vp,tagCount[i].v);
		}
		centers[i] = tagCount[i].position;
	}

	qsort(tagCount,totalPositions,sizeof(doubleIndex),&cmpDoubleIndex);


	int d = minDist;
	if (d < 1) d =1;
	int expectedNumPeaks = appearentSize/d*2;
	if (expectedNumPeaks < 10) expectedNumPeaks = 10;

	int halfSize = (peakSize)/2;

	for (int i=0;i<totalPositions;i++) {

		unsigned int index = tagCount[i].index;
		double value = tagCount[i].v;
		int center = tagCount[i].position;

		if (value < minCount) break;

		char d = tags[index].d;
//fprintf(stderr, "%d\n", d);
		if (strand==STRAND_BOTH) d=0;

		if (!mask[index]) {
			int start = center-halfSize;
			int end = start+peakSize;
			if (start < 1) {
				start = 1;
			}
			if (end > maxEND) {	
				end = maxEND;
			}
		
			sprintf(pname, "%s-%d",chr,pid++);	
			//sprintf(pname, "%d",rand());	
			//fprintf(stderr, "%s\t%s\t%d\t%d\t%d\t%d\t%f\n",pname,chr,start,end,center,d,value);
			putativePeaks->addPeak(pname,chr,start,end,center,d,value,0,NULL,-1,0);
		}
		if (regionFlag == 0 || !mask[index]) {
			for (int j=index;j<totalPositions;j++) {
				if (centers[j]-center < minDist) {
					if (strand==STRAND_BOTH || d==tags[j].d) {
						mask[j]=1;
					}
				} else {
					break;
				}
			}
			for	(int j=index-1;j>=0;j--) {
				if (center-centers[j] < minDist) {
					if (strand==STRAND_BOTH || d==tags[j].d) {
						mask[j]=1;
					}
				} else {
					break;
				}
			}
		}
	}

	delete []tagCount;
	delete []centers;
	delete []mask;
}
int cmpDoubleIndex(const void* a, const void* b) {
	double av = ((doubleIndex*)a)->v;
	double bv = ((doubleIndex*)b)->v;
	if (av < bv) return 1;
	if (av > bv) return -1;
	return 0;
}

void ChrTags::getTagLengthDistribution(double* dist,int max) {
	forceSingleReadFlag = 1;
	loadTags();
	for (int i=0;i<totalPositions;i++) {
		int len = tags[i].len;
		if (len < 0) len = 0;
		if (len >= max) {
			len = max-1;
		}
		dist[len]+=tags[i].v;
	}
	freeTags();
	forceSingleReadFlag = 0;
}

void ChrTags::getTagCountDistribution(double* dist,int max,int scaleFactor) {
	//forceSingleReadFlag = 1;
	loadTags();
	for (int i=0;i<totalPositions;i++) {
		int v = 0;
		if (pairedEndFlag) {
			v = (int)floor((double)petags[i].v*scaleFactor);
		} else {
			v = (int)floor((double)tags[i].v*scaleFactor);
		}
		if (v >= max) {
			v = max-1;
		}
		dist[v]++;
	}
	freeTags();
	//forceSingleReadFlag = 0;
}


Tag* ChrTags::getCoverageTags(int &coveragePositions,double normFactor,char strand,int resolution,
								int setFragLength,int style,int lastTagFlag,UniqMapChrs* umc,
								double normLength) {

	loadTags();

	resolution = 1;
	coveragePositions = totalPositions*2;
	Tag* coverageTags = new Tag[coveragePositions];

	int fragLength = setFragLength;
	int halfFragLength = fragLength/2;
	double totalMappable = (double)fragLength;

	//int firstPosition = tags[0].p;
	int lastPosition = tags[totalPositions-1].p;
	double totalUsed = 0.0;

	int cIndex = 0;
	for (int i=0;i<totalPositions;i++) {
		if (strand != STRAND_BOTH && tags[i].d != strand) continue;

		int padj = 0;
		if (tags[i].d == STRAND_POSITIVE) padj = -1;

		if (style == UCSC_UNMETHYLATED_CpG) tags[i].v = -1*(1.0-tags[i].v);

		totalUsed += tags[i].v;
		tags[i].v *= (float)normFactor;
		if (setFragLength == FRAGMENT_LEN_GIVEN) {
			fragLength = tags[i].len;
		}
		if (normLength > 1e-10) {
			int L = fragLength;
			if (L < COVERAGE_MIN_NORMLENGTH) L = COVERAGE_MIN_NORMLENGTH;
			tags[i].v *= normLength/((double)L);
		}

		if (umc != NULL) {
			int start = tags[i].p+padj-halfFragLength;
			int end = tags[i].p+padj+halfFragLength;
			double map = (double)umc->countRegion(start,end,tags[i].d);
			double frac = map/totalMappable;
			if (frac < 0.10) frac = 0.10;
			tags[i].v /= frac;
			//fprintf(stderr, "%s:%d-%d\t%d\t%lf\n", chr, start, end, tags[i].d, frac);
		}

		coverageTags[cIndex].v=tags[i].v;
		coverageTags[cIndex].p=tags[i].p+padj;
		if (resolution > 1)
			coverageTags[cIndex].p=(int) (resolution*floor(((double)(tags[i].p+padj))
														/((double)resolution)+0.5));

		if (lastTagFlag && coverageTags[cIndex].p>lastPosition)
			coverageTags[cIndex].p=lastPosition;


		coverageTags[cIndex].d=0;
		coverageTags[cIndex].len=tags[i].len;
		cIndex++;
		coverageTags[cIndex].d = 0;
		coverageTags[cIndex].len = tags[i].len;

		if (setFragLength == FRAGMENT_LEN_GIVEN) {
			fragLength = coverageTags[cIndex].len;
		}


		if (tags[i].d == 0) {
			coverageTags[cIndex].p = tags[i].p+padj+fragLength;
			if (resolution > 1)
				coverageTags[cIndex].p=(int) (resolution*floor(((double)(tags[i].p+padj+fragLength))
																	/((double)resolution)+0.5));
			coverageTags[cIndex].v = -1*tags[i].v;
		} else {
			coverageTags[cIndex].p = tags[i].p+padj-fragLength;
			if (resolution > 1)
				coverageTags[cIndex].p=(int) (resolution*floor(((double)(tags[i].p+padj-fragLength))
																	/((double)resolution)+0.5));
			if (coverageTags[cIndex].p < 0)
				coverageTags[cIndex].p = 1;
			coverageTags[cIndex-1].v = -1*tags[i].v;
			coverageTags[cIndex].v = tags[i].v;
		}
		if (lastTagFlag && coverageTags[cIndex].p>lastPosition)
			coverageTags[cIndex].p=lastPosition;
		cIndex++;
	}

	qsort(coverageTags,cIndex,sizeof(Tag),&cmpTags);

	int last = 0;
	for (int i=1;i<cIndex;i++) {
		if (coverageTags[i].p == coverageTags[i-1].p) {
			coverageTags[last].v += coverageTags[i].v;
		} else {
			last++;
			if (last == i) continue;
			coverageTags[last].copy(&(coverageTags[i]));
		}
	}
	coveragePositions = last+1;

	float total = 0.0;
	for (int i=0;i<coveragePositions;i++) {
		total += coverageTags[i].v;
		coverageTags[i].v = total;
	}

	freeTags();

	return coverageTags;
}

Tag* ChrTags::getRatioTags(int &ratioPositions, Tag* coverageTags, int coveragePositions, Tag* inputTags, 
				int inputPositions, double pseudoCounts) {

	ratioPositions =coveragePositions+inputPositions;
	Tag* ratioTags = new Tag[ratioPositions];

	int inputIndex = 0;
	int ratioIndex = 0;
	double curValue = 0.0;
	double curInputValue = 0.0;
	for (int i=0;i<coveragePositions;i++) {
		int cp = coverageTags[i].p;
		while (inputIndex < inputPositions && inputTags[inputIndex].p < cp) {
			curInputValue = inputTags[inputIndex].v;
			ratioTags[ratioIndex].p = inputTags[inputIndex].p;
			ratioTags[ratioIndex].d = 0;
			ratioTags[ratioIndex].v = (curValue+pseudoCounts)/(curInputValue+pseudoCounts);
			ratioIndex++;
			inputIndex++;
		}
		curValue = coverageTags[i].v;
		if (inputIndex < inputPositions) {
			if (cp == inputTags[inputIndex].p) {
				curInputValue = inputTags[inputIndex].v;
				inputIndex++;
			}
		}
		ratioTags[ratioIndex].p = cp;
		ratioTags[ratioIndex].d = 0;
		ratioTags[ratioIndex].v = (curValue+pseudoCounts)/(curInputValue+pseudoCounts);
		ratioIndex++;
	}
	if (inputIndex < inputPositions) {
		for (int i=inputIndex;i<inputPositions;i++) {
			curInputValue = inputTags[inputIndex].v;
			ratioTags[ratioIndex].p = inputTags[inputIndex].p;
			ratioTags[ratioIndex].d = 0;
			ratioTags[ratioIndex].v = (curValue+pseudoCounts)/(curInputValue+pseudoCounts);
			ratioIndex++;

		}
	}

	return ratioTags;
}

void ChrTags::printBedGraph(FILE* fp,double normFactor,char strand,int resolution,int negFlag,
				int setFragLength, double reductionRatio,int style,int condenseFlag, int lastTagFlag,
				UniqMapChrs* umc,Peak* circosPeak,ChrTags* input, double pseudoCounts, int logFlag,
				double inputNormFactor, int inputFragLength, double normLength) {
	//fprintf(stderr, "resolution = %d condenseFlag = %d\n", resolution, condenseFlag);


	int coveragePositions = totalPositions*2;
	Tag* coverageTags = getCoverageTags(coveragePositions,normFactor,strand,resolution,setFragLength,
									style,lastTagFlag,umc,normLength);

	if (input != NULL) {
		int coveragePositionsInput = input->totalPositions*2;
		Tag* coverageTagsInput = input->getCoverageTags(coveragePositionsInput,inputNormFactor,strand,resolution,
									inputFragLength,style,lastTagFlag,umc,normLength);

		
		int ratioPositions =coveragePositions+coveragePositionsInput;
		Tag* ratioTags = getRatioTags(ratioPositions,coverageTags, coveragePositions, coverageTagsInput, 
								coveragePositionsInput,pseudoCounts);		
/*
for (int i =0;i<coveragePositions;i++) {
	fprintf(stderr, "Exp: %d\t%f\n", coverageTags[i].p,coverageTags[i].v);
	if (coverageTags[i].p > 3010000) break;
}
for (int i =0;i<coveragePositionsInput;i++) {
	fprintf(stderr, "Inp: %d\t%f\n", coverageTagsInput[i].p,coverageTagsInput[i].v);
	if (coverageTagsInput[i].p > 3010000) break;
}
for (int i =0;i<ratioPositions;i++) {
	fprintf(stderr, "Rat: %d\t%f\n", ratioTags[i].p,ratioTags[i].v);
	if (ratioTags[i].p > 3010000) break;
}*/

		delete []coverageTagsInput;
		delete []coverageTags;
		coverageTags = ratioTags;
		coveragePositions = ratioPositions;
	}

	if (reductionRatio < 1.0) {
		double adjust = ((double)totalTags)/((double)(ogTotalTags));
		reductionRatio /= adjust;
	}


	appearentSize = 0;
	if (reductionRatio >= 1.0) {
		if (resolution > 1) {
			//fprintf(stderr, "here... %d\n", condenseFlag);

			int curIndex = 0;
			for (int i=0;i<coverageTags[coveragePositions-1].p;i+=resolution) {
				while (curIndex < coveragePositions && coverageTags[curIndex].p < i) {
					curIndex++;
				}
				if (curIndex >= coveragePositions) break;
				int nextIndex = curIndex;
				while (nextIndex < coveragePositions && coverageTags[nextIndex].p < i+resolution) {
					nextIndex++;
				}
				if (nextIndex == coveragePositions) nextIndex = coveragePositions-1;
				double value = 0.0;
				double N = 0.0;
				if (coverageTags[curIndex].p > i+resolution) {
					if (curIndex ==0) continue;
					if (condenseFlag == UCSC_RES_MAX) {
						if (value < coverageTags[curIndex-1].v) value = coverageTags[curIndex-1].v;
					} else if (condenseFlag == UCSC_RES_AVG) {
						N += (double)resolution;
						value += ((double)resolution)*coverageTags[curIndex-1].v;
					}
				} else {
					if (curIndex != 0) {
						if (condenseFlag == UCSC_RES_MAX) {
							if (value < coverageTags[curIndex-1].v) value = coverageTags[curIndex-1].v;
						} else if (condenseFlag == UCSC_RES_AVG) {
							double d = (double)(coverageTags[curIndex].p-i);
							value += coverageTags[curIndex-1].v*d;
							N += d;
						} 
					} else {
						if (condenseFlag == UCSC_RES_AVG) {
							double d = (double)(coverageTags[curIndex].p-i);
							value += 0.0;
							N += d;
						}
					}
					for (int j=curIndex;j<nextIndex;j++) {
						if (curIndex == coveragePositions-1) break;
						if (condenseFlag == UCSC_RES_MAX) {
							if (value < coverageTags[j].v) value = coverageTags[j].v;
						} else if (condenseFlag == UCSC_RES_AVG) {
							double d = (double)(coverageTags[j+1].p-coverageTags[j].p);
							if (coverageTags[j+1].p > i+resolution) {
								d = (double) ((i+resolution)-coverageTags[j].p);
							}
							value += coverageTags[j].v*d;
							N += d;
						}
					}
				}
				if (condenseFlag == UCSC_RES_AVG) {
					if (N > 0.5) value /= N;
				}
				if (value > 0.05 || value < -0.05 || input!=NULL) {
					int good2print = 1;
					if (circosPeak != NULL) {
						if (circosPeak->end < i || circosPeak->start > i+resolution) {
							good2print = 0;
						}
					}
					if (good2print) {
						if (logFlag) {
							value = log(value)/log(2.0);
						}
						if (negFlag) {
							fprintf(fp,"%s\t%d\t%d\t%.1f\n",chr,i,i+resolution,-1*value);
						} else {
							fprintf(fp,"%s\t%d\t%d\t%.1f\n",chr,i,i+resolution,value);
						}
					}
				}
				appearentSize = i+resolution;
				
			}
		} else {
			for (int i=0;i<coveragePositions-1;i++) {
				if (coverageTags[i].v > 0.05 || coverageTags[i].v < -0.05) {
					int pp= coverageTags[i].p;
					if (pp < 1) pp= 1;
					if (pp < coverageTags[i+1].p) {
						int good2print = 1;
						if (circosPeak != NULL) {
							if (circosPeak->end < pp || circosPeak->start > coverageTags[i+1].p) {
								good2print = 0;
							}
						}
						if (good2print) {
							double value = coverageTags[i].v;
							if (logFlag) {
								value = log(value)/log(2.0);
							}
							if (negFlag) {
								fprintf(fp,"%s\t%d\t%d\t%.1f\n",chr,pp,coverageTags[i+1].p,-1*value);
							} else {
								fprintf(fp,"%s\t%d\t%d\t%.1f\n",chr,pp,coverageTags[i+1].p,value);
							}
						}
					}
				}
				appearentSize = coverageTags[i+1].p;
			}
		}

	} else {
		//remove tags to be able to show them all
		double normAdjust = normFactor/reductionRatio;
		double runningTotal = fabs(coverageTags[0].v/normAdjust);
		int lastValue = 0;
		int lastPosition = coverageTags[0].p;
		int lastIndex = 0;
		int ZZ=0;
		for (int i=1;i<coveragePositions;i++) {
			runningTotal += fabs((coverageTags[i].v-coverageTags[i-1].v)/normFactor*reductionRatio);
			if ((int)runningTotal > lastValue) {
				if (lastPosition != coverageTags[i].p) {
					int start = coverageTags[lastIndex].p;
					int end = coverageTags[i].p;
					float sum= 0.0;
					float range = 0.0;
					for (int j=lastIndex;j<i;j++) {
						float diff = (float)(coverageTags[j+1].p - coverageTags[j].p);
						sum += diff*coverageTags[j].v;
						range += diff;
					}
					float avgValue =0.0;
					if (range > 0) avgValue = sum/range;
					if (start < 1) start = 1;
					if ((avgValue > 0.05 || avgValue < -0.05 || input != NULL) && (start < end)) {
						int good2print = 1;
						if (circosPeak != NULL) {
							if (circosPeak->end < start || circosPeak->start > end) {
								good2print = 0;
							}
						}
						if (good2print) {
							if (logFlag) {
								avgValue = log(avgValue)/log(2.0);
							}
							if (negFlag) {
								fprintf(fp, "%s\t%d\t%d\t%.1f\n",chr,start,end,-1*avgValue);
							} else {
								fprintf(fp, "%s\t%d\t%d\t%.1f\n",chr,start,end,avgValue);
							}
						}
					}
					appearentSize = end;
	
					lastValue=(int)runningTotal;
					lastPosition = coverageTags[i].p;
					lastIndex = i;
					ZZ++;
				} else {
					lastIndex = i;
				}
			}
		}
	}

	delete []coverageTags;

}

void ChrTags::annotateTagLocations(ChrPeaks* annotations, FILE* annFile, Doubletable* stats) {
	
	loadTags();
	
	int annIndex = 0;
	for (int i=0;i<totalPositions;i++) {
		int p = tags[i].p;
		char d = tags[i].d;
		float v = tags[i].v;
		int len = tags[i].len;
		//int mid = p;
		if (len > 0) {
			if (d == 0) {
				//mid = p+len/2;
			} else {
				//mid = p-len/2;
			}
		}
		for (int j=annIndex;j<annotations->numPeaks;j++) {
			int s = annotations->peaks[j]->start;
			int e = annotations->peaks[j]->end;
			if (e < p) annIndex++;
			if (p <= e && p >= s) {
				char* str = annotations->peaks[j]->data;
				double statCount = stats->search(str);
				if (statCount < EMPTY_DOUBLE_CHECK) {
					statCount = 0.0;
				}
				statCount += v;
				stats->insert(statCount, str);
				if (annFile != NULL) {
					char* nn = annotations->peaks[j]->name;
					if (annotations->peaks[j]->ogname != NULL)
						nn = annotations->peaks[j]->ogname;
					
					fprintf(annFile,"\t%s\t%d\t%d\t%.1f\t%d\t%s\t%s\n", chr,p,d,v,len, nn,str);
				}
				break;
			}
		}
	}
	freeTags();
	
}

void ChrTags::getPeakTagCounts(Doubletable* counts, ChrPeaks* cp, char strand) {
	
	loadTags();

	double* totals = new double[cp->numPeaks];
	for (int i=0;i<cp->numPeaks;i++) {
		totals[i]=0.0;
	}

	int cpIndex = 0;
	for (int i=0;i<totalPositions;i++) {
		int p = tags[i].p;
		char d = tags[i].d;
		float v = tags[i].v;
		int len = tags[i].len;
		int mid = p;
		if (len > 0) {
			if (d == 0) {
				mid = p+len/2;
			} else {
				mid = p-len/2;
			}
		}


		while (cpIndex < cp->numPeaks && cp->peaks[cpIndex]->end < mid) {
			cpIndex++;
		}
		if (cpIndex == cp->numPeaks) break;

		for (int j=cpIndex;j<cp->numPeaks;j++) {
			if (mid < cp->peaks[j]->start) break;
			if (strand == STRAND_SEPARATE && d != cp->peaks[j]->strand) continue;
			if (cp->peaks[j]->start <= mid && mid <= cp->peaks[j]->end) {
				totals[j] += v;
			}
		}
	}
	freeTags();

	for (int i=0;i<cp->numPeaks;i++) {
		counts->insert(totals[i],cp->peaks[i]->name);
	}
	delete []totals;
}

void ChrTags::checkTagSeqBias(char* genomeDirectory,NucleotideFreq* nf,NucleotideFreq* nfuniq,
				NucleotideFreq* nfctrl,	int freqStart,int freqEnd, int fragLen,
				OligoArray* oligos, int oligoStart, int oligoEnd) {

	forceSingleReadFlag=1;
	loadTags();
	if (gcFreq == NULL) {
		gcFreq = new float[totalPositions];
	}
	loadSequence(genomeDirectory);
	if (seq == NULL) {
		freeTags();
		forceSingleReadFlag=0;
		return;
	}

	SeqFreqStats* sfs = new SeqFreqStats();
	int maxLength = -1*freqStart+freqEnd+fragLen+1000;
	char* fragSeq = new char[maxLength];
	int olen = oligoEnd-oligoStart+1;

	int start = 0;
	int end = 0;
	int gcStart = -1*freqStart;
	int gcEnd = gcStart+fragLen;
	for (int i=0;i<totalPositions;i++) {
		if (tags[i].d == STRAND_POSITIVE) {
			start = tags[i].p+freqStart;
			end = tags[i].p+fragLen+freqEnd;
		} else {
			start = tags[i].p-fragLen-freqEnd;
			end = tags[i].p-freqStart;
		}
		getSequence(fragSeq,start,end,tags[i].d);
		//fprintf(stderr, "%s\n", fragSeq);
		gcFreq[i]=-1.0;
		if (fragSeq[0] != '\0') {
			nf->addSequence(NULL, fragSeq,tags[i].v,gcStart,gcEnd,sfs);
			if (sfs->N > 0) {
				//fprintf(stderr, "gcFreq[%d]\t%lf\n", i, sfs->GC);
				gcFreq[i]=sfs->GC;
			}
			nfuniq->addSequence(NULL, fragSeq,1.0,gcStart,gcEnd,sfs);
		} else {
			gcFreq[i]=-1.0;
		}
		if (oligos != NULL) {
			if (tags[i].d == STRAND_POSITIVE) {
				start = tags[i].p+oligoStart;
				end = tags[i].p+oligoEnd;
			} else {
				start = tags[i].p-oligoEnd;
				end = tags[i].p-oligoStart;
			}
			getSequence(fragSeq,start,end,tags[i].d);
			//fprintf(stderr, "%s\n", fragSeq);
			if (fragSeq[0] != '\0') {
//fprintf(stderr, "adding from tags %s\n", fragSeq);
				int bad = 0;
				for (int j=0;j<olen;j++) {
					if (fragSeq[j] == 'N') {
						bad = 1;
						break;
					}
				}
				//if (bad == 0) oligos->addOligo(fragSeq,tags[i].v,0.0);
				if (bad == 0) oligos->addOligo(fragSeq,1.0,0.0,0.0);
			}
		}
	}
	if (fragSeq != NULL) delete []fragSeq;
	delete sfs;
	//delete []gcFreq;
	//gcFreq = NULL;

	nfctrl->addChr(seq,gcStart,gcEnd);
	if (oligos != NULL) {
		fprintf(stderr, "\t\t\tAnalyzing Oligos:");
		for (int i=0;i<chrLen-olen;i++) {
			if (i % 10000000 == 0) fprintf(stderr, ".");
			int bad = 0;
			for (int j=i;j<i+olen;j++) {
				if (seq[j] == 'N') {
					bad = j;
					break;
				}
			}
			if (bad) {
				i=bad;
				continue;
			}
			char* s = &(seq[i]);
			char tmp = seq[i+olen];
			seq[i+olen] = '\0';
			//fprintf(stderr, "adding %s\n", s);
			oligos->addOligo(s,0,1.0,0.0);
			seq[i+olen] = tmp;
		}
		fprintf(stderr, "\n");
	}

	deleteSequence();
	freeTags();
	forceSingleReadFlag=0;
}
void ChrTags::readAndSave() {
	loadTags();
	dontSAVE=0;
	print();
	freeTags();
}
void ChrTags::normalizeTagCountsGC(NucleotideFreq* nf) {

	loadTags();

	totalTags = 0.0;
	for (int i=0;i<totalPositions;i++) {
		if (gcFreq[i] > -1.0) {
			int index = (int)floor(gcFreq[i]/nf->gcInc+0.0001);
			double adjRatio = nf->gcNorm[index];
			tags[i].v *= adjRatio;
		}
		totalTags += tags[i].v;
	}

	dontSAVE = 0; // this normally protects us from screwing with the tag file
	print();
	freeTags();
}

OligoProfile::OligoProfile(int nlen, int noffset) {
	offset = noffset;
	length = nlen;
	revoppFlag = 0;
	normFlag = 0;
	pProfile = new double[length];
	nProfile = new double[length];
	for (int i=0;i<length;i++) nProfile[i]=0.0;
	for (int i=0;i<length;i++) pProfile[i]=0.0;
	N =0.0;
}
OligoProfile::~OligoProfile() {
	if (pProfile != NULL) delete []pProfile;
	if (nProfile != NULL) delete []nProfile;
}
void OligoProfile::mergeRevopps(OligoProfile* op1, OligoProfile* op2) {
	int start = op1->offset;
	int end = op1->offset+op1->length;
	int end2 = -1*op2->offset;
	int start2 = -1*(op2->offset+op2->length);
	if (start2 < start) start = start2;
	if (end2 > end) end = end2;

	int length = end-start;
	int offset = start;
	double *pProfile = new double[length];
	double *nProfile = new double[length];
	double *counts = new double[length];
	for (int i=0;i<length;i++) {
		pProfile[i] = 0.0;
		nProfile[i] = 0.0;
		counts[i] = 0.0;
	}
	double max = 0;
	int rOffset = op1->offset-start;
	for (int i=rOffset;i<rOffset+op1->length;i++) {
		pProfile[i] += op1->pProfile[i-rOffset];
		nProfile[i] += op1->nProfile[i-rOffset];
		counts[i] += op1->N;
		if (counts[i] > max) max = counts[i];
	}
	for (int i=0;i<op2->length;i++) {
		int noffset = -1*op2->offset-offset;
		pProfile[noffset-i] += op2->nProfile[i];
		nProfile[noffset-i] += op2->pProfile[i];
		counts[noffset-i] += op2->N;
		if (counts[i] > max) max = counts[i];
	}
	if (max > 1.0) {
		for (int i=0;i<length;i++) {
			if (counts[i] > 1.0 && counts[i] < max) {
				double ratio = max/counts[i];
				pProfile[i] *= ratio;
				nProfile[i] *= ratio;
			}
		}
	}

	delete [](op1->pProfile);
	delete [](op1->nProfile);
	op1->pProfile = new double[length];
	op1->nProfile = new double[length];
	op1->length = length;
	op1->offset = offset;
	op1->revoppFlag = 1;
	op1->N = max;
	for (int i=0;i<length;i++) {
		op1->pProfile[i] = pProfile[i];
		op1->nProfile[i] = nProfile[i];
	}
	if (op1 != op2) {
		delete [](op2->pProfile);
		delete [](op2->nProfile);
		op2->pProfile = new double[length];
		op2->nProfile = new double[length];
		op2->length = length;
		op2->offset = offset;
		op2->revoppFlag = 1;
		op2->N = max;
		for (int i=0;i<length;i++) {
			op2->pProfile[i] = pProfile[length-i-1];
			op2->nProfile[i] = nProfile[length-i-1];
		}
	}

	delete []counts;
	delete []pProfile;
	delete []nProfile;
}
void OligoProfile::normalize() {
	double pLevel = 0.0;
	double pN = 0.0;
	for (int i=OLIGO_NORMALIZATION_BUFFER-offset;i<length;i++) {
		pLevel += pProfile[i];
		pN+=1.0;
	}
	if (pN < 1.0 || pLevel < 1.0) {
		pLevel = 1.0;
		pN = 1.0;
	}
	pLevel /= pN;
	for (int i=0;i<length;i++) {
		pProfile[i] /= pLevel;
	}
	double nLevel = 0.0;
	double nN = 0.0;
	for (int i=0;i<-1*OLIGO_NORMALIZATION_BUFFER-offset;i++) {
		nLevel += nProfile[i];
		nN+=1.0;
	}
	if (nN < 1.0 || nLevel < 1.0) {
		nLevel = 1.0;
		nN = 1.0;
	}
	nLevel /= nN;
	for (int i=0;i<length;i++) {
		nProfile[i] /= nLevel;
	}
	normFlag = 1;
}
void OligoProfile::print(FILE* fp) {
	if (fp == NULL) return;
	fprintf(fp, "Position(%.0lf total)\tstrand +\tstrand -\n", N);
	for (int i=0;i<length;i++) {
		fprintf(fp, "%d\t%.2lf\t%.2lf\n", i+offset,pProfile[i],nProfile[i]);
	}
}

void ChrTags::normalizeTagCountsOligos(char* genomeDirectory, Hashtable* oligos, int oligoLength,
				int regionStart, int regionEnd, double minFold, double maxFold, float maxPerBp,int normFlag) {

	loadTags();
	loadSequence(genomeDirectory);
	if (seq == NULL) {
		freeTags();
		return;
	}

	int regionLength = regionEnd-regionStart;
	int maxRegion = regionStart;
	int halfOligoLength = oligoLength/2;

	float* ratios = NULL;
	if (normFlag) {
		ratios = new float[totalPositions];
		for (int i=0;i<totalPositions;i++) ratios[i]=1.0;
	}



	if (maxRegion < 0) maxRegion *= -1;
	if (regionEnd > 0 && regionEnd > maxRegion) maxRegion = regionEnd;
	if (regionEnd < 0 && regionEnd < -1*maxRegion) maxRegion = -1*regionEnd;

	int tagIndex = 0;
	for (int i=0;i<chrLen-oligoLength;i++) {

		int rp = i+halfOligoLength; // recenter to the middle of the oligo

		while (tagIndex < totalPositions && (tags[tagIndex].p+maxRegion < rp)) {
			tagIndex++;
		}
		if (tagIndex == totalPositions) break;

		if (rp < tags[tagIndex].p-maxRegion) {
			i = tags[tagIndex].p-maxRegion-halfOligoLength;
			if (i < 0) i=0;
			rp = i+halfOligoLength;
		}


		char* curOligo = &(seq[i]);
		char tmp = seq[i+oligoLength];
		seq[i+oligoLength] = '\0';

		OligoProfile*  op =  (OligoProfile*) oligos->search(curOligo);
		if (op == NULL) {
			op = new OligoProfile(regionLength,regionStart);
			oligos->insert(op,curOligo);
		}
		op->N++;
		for (int j=tagIndex;j<totalPositions;j++) {
			int p = tags[j].p;
			if (p > rp+maxRegion-1) break;
			if (p < rp-maxRegion) continue;
			char d = tags[j].d;
			int k = p-(rp+regionStart);
			if (normFlag) {
				if (d==STRAND_POSITIVE) {
					//tags[j].v /= (op->pProfile[k] + (oligoLength-1.0))/oligoLength;
					ratios[j] /= op->pProfile[k];
				} else {
					//tags[j].v /= (op->nProfile[k] + (oligoLength-1.0))/oligoLength;
					ratios[j] /= op->nProfile[k];
				}
			} else {
				float v = tags[j].v;
				if (v > maxPerBp) v = maxPerBp;
				if (d==STRAND_POSITIVE) {
					op->pProfile[k] += v;
				} else {
					op->nProfile[k] += v;
				}
			}
		}
		seq[i+oligoLength] = tmp;
	}

	if (normFlag) {
		for (int i=0;i<totalPositions;i++) {
			if (ratios[i] < (float)minFold) ratios[i] = (float)minFold;
			if (ratios[i] > (float)maxFold) ratios[i] = (float)maxFold;
			tags[i].v *= ratios[i];
		}
		delete []ratios;
		dontSAVE = 0; // this normally protects us from screwing with the tag file
		print();
	}
	freeTags();
	deleteSequence();
}


void ChrTags::normalizeTagCountsFixedOligo(char* genomeDirectory, OligoArray* oligos,
							int oligoStart, int oligoEnd, double minFold, double maxFold) {

	loadTags();
	loadSequence(genomeDirectory);
	if (seq == NULL) {
		freeTags();
		return;
	}

	int olen = oligoEnd-oligoStart+1;
	char* fragSeq = new char[olen*2];

	int start = 0;
	int end = 0;
	double vt=0,vbg=0;
	for (int i=0;i<totalPositions;i++) {
		if (tags[i].d == STRAND_POSITIVE) {
			start = tags[i].p+oligoStart;
			end = tags[i].p+oligoEnd;
		} else {
			start = tags[i].p-oligoEnd;
			end = tags[i].p-oligoStart;
		}
		getSequence(fragSeq,start,end,tags[i].d);
		if (fragSeq[0] != '\0') {
			int bad = 0;
			for (int j=0;j<olen;j++) {
				if (fragSeq[j] == 'N') {
					bad = 1;
					break;
				}
			}
			if (bad == 1) continue;
			int cpu=0;
			int found = oligos->searchOligo(fragSeq,0,vt,vbg,OLIGO_SEARCH_KEEP,cpu);
			if (found > 0) {
				double ratio = oligos->utilityOligos[0][0]->value;
				if (ratio > 0.0000001) {
					ratio = 1.0/ratio;
				}
				if (ratio < minFold) ratio = minFold;
				if (ratio > maxFold) ratio = maxFold;
				tags[i].v *= ratio;
				fprintf(stderr, "%s\t%d\t%d\t%s\t%f\t%lf\t%lf\n",chr,tags[i].p,tags[i].d,fragSeq,ratio,minFold,maxFold);
				//tags[i].v=ratio;
			}
			//if (bad == 0) oligos->addOligo(fragSeq,tags[i].v,0.0);
		}
		totalTags += tags[i].v;
	}

	dontSAVE = 0; // this normally protects us from screwing with the tag file
	print();
	freeTags();
	deleteSequence();
	delete []fragSeq;
}




void ChrTags::getPETagTotals(PeakLibrary* peaks1, int &peakIndex1) {

	// initialize peakIndex for chromosome
	int numPeaks1 = peaks1->numPeaks;
	Peak** peakArray1 = peaks1->peakOrder;
	

	while (peakIndex1 < numPeaks1 && 
		chrcmp(&(peakArray1[peakIndex1]->chr),&chr) < 0) peakIndex1++;

	if (peakIndex1 >= peaks1->numPeaks) return;

	if (strcmp(peaks1->peakOrder[peakIndex1]->chr,chr) != 0) return;

	loadTags();

	int lastIndex1 = peakIndex1;
	while (lastIndex1 < numPeaks1 && 
		chrcmp(&(peakArray1[lastIndex1]->chr),&chr) <= 0) lastIndex1++;

	int startIndex1 = peakIndex1;
	for (int i=startIndex1;i<lastIndex1;i++) {
		peakArray1[i]->v = 0.0;
		peakArray1[i]->focusRatio = 1.0;
	}

	int curMinPos = peakArray1[peakIndex1]->start;
	int curEndPos = peakArray1[peakIndex1]->end;
	int quit = 0;
	for (int i=0;i<totalPositions;i++) {
		int p = petags[i].p1;
		if (p < curMinPos) continue;
		quit = 0;
		while (p > curEndPos) {
			peakIndex1++;
			if (peakIndex1 >= lastIndex1) {
				quit=1;
				break;
			}
			curMinPos = peakArray1[peakIndex1]->start;
			curEndPos = peakArray1[peakIndex1]->end;
		}
		if (quit) break;
		float v = petags[i].v;
		for (int j=peakIndex1;j<lastIndex1;j++) {
			int cstart = peakArray1[j]->start; 
			int cend = peakArray1[j]->end; 
			if (cstart > p) break;
			if (cend < p) continue;
			peakArray1[j]->v += v;
		}
	}
	freeTags();
}


// if model != NULL, then the model is "built"
void ChrTags::scoreInteractionBoundaries(PeakLibrary* peaks1, int &peakIndex1,
				int resolution, HiCBgModel* model,int actionFlag) {

	// initialize peakIndex for chromosome
	int numPeaks1 = peaks1->numPeaks;
	Peak** peakArray1 = peaks1->peakOrder;

	while (peakIndex1 < numPeaks1 && 
		chrcmp(&(peakArray1[peakIndex1]->chr),&chr) < 0) peakIndex1++;

	if (peakIndex1 >= peaks1->numPeaks) return;
	if (strcmp(peaks1->peakOrder[peakIndex1]->chr,chr) != 0) return;

	fprintf(stderr, "\t\t%s\n", chr);
	loadTags();


	double* scores  = new double[totalPositions-1];
	double* totals  = new double[totalPositions-1];
	
	int curStartIndex = 0;
	int curEndIndex = 0;
	double avgScore = 0.0;
	double validScores = 0.0;
	int last = -1000;
	for (int i=0;i<totalPositions-1;i++) {
		if (i% 1000000==0) fprintf(stderr, "%d of %d\n", i,totalPositions);
		scores[i] = -1e10;
		int p1 = petags[i].p1;
		int pnext1 = petags[i+1].p1;
		int midpoint = (p1+pnext1)/2;
		if (midpoint < last+50) continue;
		last = midpoint;
		while (petags[curStartIndex].p1 < midpoint - resolution) {
			curStartIndex++;
		}
		while (curEndIndex < totalPositions-1 && petags[curEndIndex].p1 < midpoint + resolution) {
			if (strcmp(petags[curEndIndex].chr2,chr)!=0) petags[curEndIndex].p2 = -1;
			curEndIndex++;
		}
		double c11=0.0;
		double c12=0.0;
		double c21=0.0;
		double c22=0.0;
		for (int j=curStartIndex;j<curEndIndex;j++) {
			if (petags[j].p2 == -1) continue;
			if (petags[j].p1 < midpoint) {
				if (petags[j].p2 < midpoint) c11+=1.0; 
				else c12+=1.0; 
			} else {
				if (petags[j].p2 < midpoint) c21+=1.0; 
				else c22+=1.0; 
			}
		}
		double d1 = c11+c12;
		double d2 = c21+c22;
		if (d1 < 1.0) d1 = 1.0;
		if (d2 < 1.0) d2 = 1.0;
		scores[i] = (c11/d1+c22/d2)/2.0;
		totals[i] = d1+d2;
		validScores += 1.0;
		avgScore += scores[i];
		//fprintf(stderr, "%d\t%d\t%lf\t%.1lf\n",midpoint, curEndIndex-curStartIndex,scores[i],d1+d2);
	}
	avgScore /= validScores;
	double stdScore = 0.0;
	for (int i=0;i<totalPositions-1;i++) {
		if (scores[i] < -1e9) continue;
		stdScore += (avgScore-scores[i])*(avgScore-scores[i]);
	}
	stdScore = sqrt(stdScore/validScores);
	fprintf(stderr, "\t%s score=%lf+/-%lf\n", chr, avgScore, stdScore);
printf("track name=\"boundaries\" type=bedGraph\n");
	for (int i=0;i<totalPositions-1;i++) {
		if (scores[i] < -1e9) continue;
		int p1 = petags[i].p1;
		int pnext1 = petags[i+1].p1;
		if (p1 == pnext1) continue;
		scores[i] = logbinomial((int)totals[i],(int)(totals[i]*scores[i]),avgScore,(int)totals[i]);
		//double zscore = (scores[i]-avgScore)/stdScore;
		printf("%s\t%d\t%d\t%lf\n",chr,p1,pnext1,scores[i]);
	}
	delete []scores;
	delete []totals;

	exit(0);
}




// if model != NULL, then the model is "built"
void ChrTags::makeHiCMatrix(PeakLibrary* peaks1, int &peakIndex1, PeakLibrary* peaks2,
				double **matrix, int resolution, HiCBgModel* model,int revFlag, 
				int numHistBins, int actionFlag, double totalInteractions,
				GenomeInteractionLibrary* gil, pthread_mutex_t* mutexMatrix,
				HiCparams* params) {



	int fragmentLengthEstimate = params->fragLengthEstimate;
	int boundaryScale = params->boundaryScale;
	int relativeFlag = params->relativeFlag;

	int halfBins = numHistBins/2;
	double** expectMatrix = NULL;
	if (model != NULL && model->expectMatrix != NULL) {
		expectMatrix = model->expectMatrix;
	}
			

	int maxDist = -1;
	if (relativeFlag != 0) {
		maxDist = relativeFlag*resolution;
	}

	// initialize peakIndex for chromosome
	int numPeaks1 = peaks1->numPeaks;
	Peak** peakArray1 = peaks1->peakOrder;
	int numPeaks2 = peaks2->numPeaks;
	Peak** peakArray2 = peaks2->peakOrder;

	double peak2TotalAvg = 0.0;
	for (int i=0;i<peaks2->numPeaks;i++) {
		peak2TotalAvg += peakArray2[i]->v;
	}
	peak2TotalAvg /= (double)peaks2->numPeaks;
	

	while (peakIndex1 < numPeaks1 && 
		chrcmp(&(peakArray1[peakIndex1]->chr),&chr) < 0) peakIndex1++;

	if (peakIndex1 >= peaks1->numPeaks) return;
	if (strcmp(peaks1->peakOrder[peakIndex1]->chr,chr) != 0) return;

	fprintf(stderr, "\t\t%s\n", chr);
	loadTags();


	//----- First part identifies PE tags that should be considered in more detail	
	//----- basically tags associated with peaksArray1

	int fineRes = 10;
	int maxFineRes = resolution/fineRes*5;
	double* fineResProximity = NULL;
	if (actionFlag & HIC_MASK_BOUNDARIES) {
		fineResProximity = new double[maxFineRes];
		for (int i=0;i<maxFineRes;i++) { 
			fineResProximity[i] = 0.0;
		}
	}
	

	int lastIndex1 = peakIndex1;
	while (lastIndex1 < numPeaks1 && 
		chrcmp(&(peakArray1[lastIndex1]->chr),&chr) <= 0) lastIndex1++;

	int startIndex1 = peakIndex1;
	int curNumPeaks1 = lastIndex1-startIndex1;
	int* endTagIndex1 = new int[curNumPeaks1];
	int* startTagIndex1 = new int[curNumPeaks1];
	double* tagTotals = new double[numPeaks1];
	for (int i=0;i<curNumPeaks1;i++) {
		startTagIndex1[i]=-1;
		endTagIndex1[i]=-1;
		tagTotals[i] = 0.0;
	}

	int curMinPos = peakArray1[peakIndex1]->start;
	int curEndPos = peakArray1[peakIndex1]->end;
	int quit = 0;
	for (int i=0;i<totalPositions;i++) {
		int p = petags[i].p1;
		if (p < curMinPos) continue;
		quit = 0;
		while (p > curEndPos) {
			peakIndex1++;
			if (peakIndex1 >= lastIndex1) {
				quit=1;
				break;
			}
			curMinPos = peakArray1[peakIndex1]->start;
			curEndPos = peakArray1[peakIndex1]->end;
		}
		if (quit) break;
		float v = petags[i].v;
		if (actionFlag & HIC_MASK_BOUNDARIES) {
			int dist = maxFineRes-1;
			if (strcmp(petags[i].chr1,petags[i].chr2)==0) {
				dist = abs(petags[i].p1-petags[i].p2);
				dist /= fineRes;
				if (dist >= maxFineRes-1) dist=maxFineRes-1;
			}
			fineResProximity[dist]+=v;
		}
		for (int j=peakIndex1;j<lastIndex1;j++) {
			int cstart = peakArray1[j]->start; 
			int cend = peakArray1[j]->end; 
			if (cstart > p) break;
			if (cend < p) continue;
			int indexIndex = j-startIndex1;
			if (startTagIndex1[indexIndex] == -1) startTagIndex1[indexIndex] = i;
			endTagIndex1[indexIndex] = i;
			tagTotals[indexIndex] += v;
		}
	}

	if (actionFlag & HIC_MASK_BOUNDARIES) {
		double t= 0.0;
		for (int i=0;i<maxFineRes;i++) t+=fineResProximity[i];
		if (t < 1.0) t= 1.0;
		for (int i=0;i<maxFineRes;i++) {
			fineResProximity[i]/=t;
			//fprintf(stdout, "%d\t%le\n",i*fineRes,fineResProximity[i]);
		}
	}
		


	//get peak subsets for histogram mode
	int** histPeakIndexes = NULL;
	double* tagTotals2 = NULL;
	int* histNumPeaks = NULL;
	Peak** peakArray2Backup = NULL;
	int histIndex1=-1;
	int histPosition1=-1;
	double avgTagTotal = 0.0;
	if (numHistBins > 0) {
		histPeakIndexes = new int*[curNumPeaks1];
		histNumPeaks = new int[curNumPeaks1];
		int* mapIndex = new int[curNumPeaks1];
		for (int i=0;i<curNumPeaks1;i++) histNumPeaks[i]=0;
		int* map = new int[numPeaks1];
		int** map2 = new int*[numPeaks2];
		for (int i=0;i<numPeaks1;i++) map[i]=-1;
		for (int i=0;i<numPeaks2;i++) map2[i]=NULL;

		
		for (int i=0;i<curNumPeaks1;i++) {
			avgTagTotal+=tagTotals[i];
			int pIndex = startIndex1+i;
			sscanf(peakArray1[pIndex]->name, "%d:%d",&histIndex1,&histPosition1);
			map[histIndex1]=i;
			mapIndex[i]=histIndex1;
		}
		avgTagTotal/=curNumPeaks1;
		fprintf(stderr, "\tAverage Tags per bin = %lf\n", avgTagTotal);
		for (int i=0;i<numPeaks2;i++) {
			sscanf(peakArray2[i]->name, "%d:%d",&histIndex1,&histPosition1);
			if (map[histIndex1]==-1) continue;

			if (map2[histIndex1]==NULL) {
				map2[histIndex1]=new int[numHistBins+1];
				map2[histIndex1][0]=-1;
			}
			int j=0;
			while (map2[histIndex1][j]!=-1) j++;
			map2[histIndex1][j]=i;
			map2[histIndex1][j+1]=-1;
		}
		for (int i=0;i<curNumPeaks1;i++) {
			int mIndex = mapIndex[i];
			histPeakIndexes[i] = NULL;
			histNumPeaks[i] = 0;
			if (map2[mIndex]==NULL) {
				//nothing...
				continue;
			}
			histPeakIndexes[i] = new int[1+numHistBins];
			for (int j=0;j<numHistBins+1;j++) {
				if (map2[mIndex][j] == -1) break;
				histPeakIndexes[i][j] = map2[mIndex][j];
				histNumPeaks[i]++;
			}
		}
		for (int i=0;i<numPeaks2;i++) {
			if (map2[i] != NULL) delete [](map2[i]);
		}
		delete []map2;
		delete []mapIndex;
		delete []map;
		peakArray2Backup = peakArray2;
		peakArray2 = new Peak*[numHistBins*2];
		tagTotals2 = new double[numHistBins*2];
	}


	//----- This next part iterates through each peak from peakArray1,
	//----- basically tags associated with peaksArray1 and looks for ones
	//----- that line up with peakArray2
	int maxTags = 0;
	for (int i=0;i<curNumPeaks1;i++) {
		int count = endTagIndex1[i]-startTagIndex1[i]+1;
		if (count > maxTags) maxTags = count;
	}
	PETag** curPETags = new PETag*[maxTags];
	PETag** iTags = new PETag*[maxTags];

	PETag** curPETagsB1 = NULL;
	PETag** curPETagsB2 = NULL;
	//int countB1 = 0;
	//int countB2 = 0;
	if (actionFlag & HIC_MASK_BOUNDARIES) {
		curPETagsB1 = new PETag*[maxTags];
		curPETagsB2 = new PETag*[maxTags];
	}
		

	HiCBgModelChr* chrModel = NULL;
	if (actionFlag & HIC_MASK_CREATEMODEL) {
		if (model == NULL) {
			fprintf(stderr, "!!! Something is wrong... contact the idiot that wrote this software !!!!\n");
			exit(0);
		}
		chrModel = (HiCBgModelChr*)model->chrs->search(chr);
		if (chrModel == NULL) {
			pthread_mutex_lock(&(model->mutex));
			chrModel = model->addChrModel(chr);
			pthread_mutex_unlock(&(model->mutex));
		}
		if (chrModel != NULL && model->stdFlag) {
			chrModel->stdFlag = 1;
			//fprintf(stderr, "Setting STD flag!\n");
		}
		if (chrModel != NULL && model->bgFlag) {
			chrModel->bgFlag = 1;
		}
	} else if (model != NULL) {
		chrModel = (HiCBgModelChr*)model->chrs->search(chr);
		if (chrModel == NULL && actionFlag & HIC_MASK_NORM_DISTANCE) {
			fprintf(stderr, "!!! Warning: Could not find background model for %s !!!\n", chr);
			chrModel = model;
		}
	}
		

	double p2Total=0.0;
	double* p2Stats = new double[numPeaks2];
	double* p2Totals = new double[numPeaks2];
	double* p2StatsB1 = NULL;
	double* p2StatsB2 = NULL;
	double* expectedArray = new double[numPeaks2];
	double* expectedArrayB1 = NULL;
	double* expectedArrayB2 = NULL;
	if (actionFlag & HIC_MASK_BOUNDARIES) {
		p2StatsB1 = new double[numPeaks2];
		p2StatsB2 = new double[numPeaks2];
		expectedArrayB1 = new double[numPeaks2];
		expectedArrayB2 = new double[numPeaks2];
	}
		


	int gilInteractionIndex = 0;


	//fprintf(stderr, "\t\tSarting Loops\n");
	for (int i=0;i<curNumPeaks1;i++) {

		if (endTagIndex1[i] < 0) continue;
		int count = endTagIndex1[i]-startTagIndex1[i]+1;
		if (count < 1) continue;

		int pIndex = startIndex1+i;

		int peakMidPoint = 0;
		peakMidPoint = peakArray1[pIndex]->refPos;
		double peak1total = peakArray1[pIndex]->v;
		if (model != NULL) {
			pthread_mutex_lock(&(model->mutex));
			if (peak1total < model->minCoverage || peak1total > model->maxCoverage) {
				model->badRegions++;
				pthread_mutex_unlock(&(model->mutex));
				if (actionFlag & HIC_MASK_NORM_SEQDEPTH) {
					continue;
				}
			} else {
				model->goodRegions++;
				pthread_mutex_unlock(&(model->mutex));
			}
		} else if (peak1total < 1.0 && actionFlag & HIC_MASK_NORM_SEQDEPTH) {
			continue;
		}
		char* chr1 = peakArray1[pIndex]->chr;
		int start1 = peakArray1[pIndex]->start;
		int end1 = peakArray1[pIndex]->end;

		int index = 0;
		for (int j=startTagIndex1[i];j<=endTagIndex1[i];j++) {
			iTags[index] = NULL;
			curPETags[index++] = &(petags[j]);
		}
		qsort(curPETags, count, sizeof(PETag*),&cmp2ndPETagPointers);

		//if background, update interaction index to find previous interactions
		if (gil != NULL) {
			while (gilInteractionIndex < gil->numInteractions && gil->interactions != NULL
						&& gil->interactions[gilInteractionIndex]->peakIndex1 < pIndex) {
				gilInteractionIndex++;
			}
			if (actionFlag & HIC_MASK_INTERACTION4CBEDFILE && gil != NULL && gil->bedFile != NULL
					&& actionFlag & HIC_MASK_INTERACTION4CBEDFILE) {
				pthread_mutex_lock(&(gil->mutex));
				for (int k=0;k<count;k++) {
					curPETags[k]->print2ndHalf(gil->bedFile);
				}
				pthread_mutex_unlock(&(gil->mutex));
			}
		}


		//double totalPvalue = 0.0;
		double preScore  = 0.0;
		double postScore  = 0.0;

		int histIndex1 = -1;
		int histPosition1 = 0;
		if (numHistBins != 0) {
			sscanf(peakArray1[pIndex]->name, "%d:%d",&histIndex1,&histPosition1);
			if (histIndex1 < 0) {
				fprintf(stderr, "Problem parsing tmp historgram peak name %s\n",
						peakArray1[pIndex]->name);
				continue;
			}
			numPeaks2 = 0;
			for (int j=0;j<histNumPeaks[i];j++) {
				peakArray2[j] = peakArray2Backup[histPeakIndexes[i][j]];
				tagTotals2[j] = tagTotals[histPeakIndexes[i][j]-startIndex1];
				numPeaks2++;
			}
		}

		if (gil != NULL && actionFlag & HIC_MASK_INTERACTIONSTATSFILES) {
						
			if (i==curNumPeaks1-1) {
				continue;
			}
			double totalInterChrInteractions = 0.0;
			double totalAllInteractions = 0.0;
			double totalLocalInteractions = 0.0;
			for (int i=0;i<count;i++) {
				if (strcmp(curPETags[i]->chr2,chr1)==0) {
					if (fabs(curPETags[i]->p1-curPETags[i]->p2) < (double)(2*fragmentLengthEstimate)) {
						totalLocalInteractions += curPETags[i]->v;
					}	
				} else {
					totalInterChrInteractions += curPETags[i]->v;
				}
				totalAllInteractions += curPETags[i]->v;
			}
			//totalAllInteractions += 10.0;
			//totalLocalInteractions += 10.0;
			//totalInterChrInteractions += 10.0;
			double r = 0.0;
			double r2 = 0.0;
			if (totalAllInteractions > 0.5) r = totalInterChrInteractions / totalAllInteractions;
			if (totalAllInteractions > 0.5) r2 = totalLocalInteractions / totalAllInteractions;
			pthread_mutex_lock(&(gil->mutex));
			fprintf(gil->bedFileInterFrac,"%s\t%d\t%d\t%lf\n",chr1,start1,end1,r);
			fprintf(gil->bedFileLocalFrac,"%s\t%d\t%d\t%lf\n",chr1,start1,end1,r2);
			fprintf(gil->bedFileCoverage,"%s\t%d\t%d\t%lf\n",chr1,start1,end1,totalAllInteractions);
			pthread_mutex_unlock(&(gil->mutex));
			continue;
		}

		p2Total=0.0;
		double p2TotalInteractions = 0.0;
		for (int j=0;j<numPeaks2;j++) {
			p2Stats[j] = 0.0;
			expectedArray[j] = 0.0;
			if (actionFlag & HIC_MASK_BOUNDARIES) {
				p2StatsB1[j] = 0.0;
				p2StatsB2[j] = 0.0;
			}
		}

		int tagIndex = 0;
		char* curTagChr2 = curPETags[tagIndex]->chr2;
		int curTagPos2 = curPETags[tagIndex]->p2;
		int curITagIndex = 0;

		for (int j=0;j<numPeaks2;j++) {		

			char* chr2 = peakArray2[j]->chr;
			int curStart = peakArray2[j]->start;
			int curEnd = peakArray2[j]->end;

			if (relativeFlag) {
				if (strcmp(chr,chr2) == 0) {
					int diff = abs(peakArray2[j]->refPos-peakMidPoint);
					if (maxDist > 0 && diff > maxDist) continue;
				} else {
					continue;
				}	
			}
					

			if (gil != NULL) {
				if (strcmp(chr,chr2) == 0) {
					int diff = abs(peakArray2[j]->refPos-peakMidPoint);
					if (gil->maxDist > 0) {
						if (diff > gil->maxDist) continue;
					} else {
						//if (diff < gil->minDist) continue;
					}
				} else {
					if (gil->maxDist > 0) {
						continue;
					}
				}
			}

			double peak2total = peakArray2[j]->v;
			if (actionFlag & HIC_MASK_NORM_SEQDEPTH) {
				if (peak2total < 0.5) {
					//fprintf(stderr, "skipping...A %d\n", curStart);
					continue;
				}
				if (model != NULL && (peak2total < model->minCoverage || peak2total > model->maxCoverage)) {
					//fprintf(stderr, "skipping...B %d\n", curStart);
					continue;
				}
			}
		
			double interactions = 0.0;

			double midpoint1 = 0.0;
			double midpoint2 = 0.0;

			for (int k=0;k<curITagIndex;k++) iTags[k] = NULL;
			curITagIndex = 0;

			int skip = 0;
			int cmp = chrcmp(&chr2,&curTagChr2);
			if (cmp < 0) {
				skip = 1;
			}
			if (cmp == 0 && curTagPos2 > curEnd) skip = 1;
			while (!skip && (cmp > 0 || curTagPos2 < curStart)) {
				tagIndex++;
				if (tagIndex >= count) {
					skip = 1;
					break;
				}
				curTagChr2 = curPETags[tagIndex]->chr2;
				curTagPos2 = curPETags[tagIndex]->p2;
				cmp = chrcmp(&chr2,&curTagChr2);
				if (cmp < 0) {
					skip = 1;
					break;
				}
				if (cmp == 0 && curTagPos2 > curEnd) {
					skip = 1;
					break;
				}
			}
			if (!skip) {
				for (int k=tagIndex;k<count;k++) {
					char* c = curPETags[k]->chr2;
					int p = curPETags[k]->p2;
					float v = curPETags[k]->v;
					//maybe speed up with: if (curTagChr2 != curPETags[tagIndex].chr2) break;
					if (strcmp(c,chr2)!=0) break;
					if (p > curEnd) break;
					iTags[curITagIndex++] = curPETags[k];
					interactions += v;
					midpoint1 += ((double)curPETags[k]->p1)*v;
					midpoint2 += ((double)curPETags[k]->p2)*v;
				}
			}
	
			p2Stats[j] = interactions;
			p2TotalInteractions += interactions;
			p2Totals[j] = peak2total;
			if (interactions > 0.0) {
				midpoint1 /= interactions;
				midpoint2 /= interactions;
			}
			//if (interactions < 1e-10) continue;
			int preFlag = 0;
			int postFlag = 0;
			int distance = INT_MAX;
			if (chrcmp(&chr2,&chr1)==0) {
				distance = peakMidPoint - peakArray2[j]->refPos;
				if (distance < -1) {
					preFlag = 1;
					distance *= -1;
				} else {
					postFlag = 1;
				}
			}
			if (actionFlag & HIC_MASK_BOUNDARIES) {
				p2StatsB1[j] = interactions;
			}

			if (actionFlag & HIC_MASK_INTERACTIONBEDFILE && gil != NULL && gil->bedFile != NULL
						&& gilInteractionIndex < gil->numInteractions 
						&& gil->interactions[gilInteractionIndex]->peakIndex1 == pIndex
						&& gil->interactions[gilInteractionIndex]->peakIndex2 ==  j) {
				int format = 0;
				if (actionFlag & HIC_MASK_INTERACTIONTAGFILE) {
					format = 1;
				}
				pthread_mutex_lock(&(gil->mutex));
				gil->interactions[gilInteractionIndex]->printPETags(gil->bedFile,iTags, curITagIndex,format);
				pthread_mutex_unlock(&(gil->mutex));
			}

			if (actionFlag & HIC_MASK_NORM_SEQDEPTH) {
				double logp=0.0;
				double simpleExpect = peak2total * peak1total / (totalInteractions);
				double expect = simpleExpect;
				double stdFactor = 1.0;
				double stdScaleFactor = 1.0;
				double zscore = 0.0;
				if (actionFlag & HIC_MASK_NORM_DISTANCE) {
					// need to adjust values so that the distance normalization is consistent with the
					// total reads per location
					expect = (peak1total*peakArray1[pIndex]->focusRatio) * (peak2total*peakArray2[j]->focusRatio);
					expect /= model->totalModelReads;
					if (distance < INT_MAX) {
						if (1) {
							double distBin = ((double)distance)/((double)resolution);
							int bin1 = (int)floor(distBin);
							int bin2 = (int)ceil(distBin);
							if (bin1 == bin2) bin2++;
							double bin1weight = 1.0-(distBin-(double)bin1);
							double bin2weight = 1.0-((double)bin2-distBin);
							expect *= (chrModel->scaleFactor[bin1]*bin1weight+chrModel->scaleFactor[bin2]*bin2weight);
							stdFactor = (chrModel->stdScaleFactor[bin1]*bin1weight+chrModel->stdScaleFactor[bin2]*bin2weight);
							stdScaleFactor = stdFactor;
						} else {
							int bin = (distance+resolution/2)/resolution;
							expect *= chrModel->scaleFactor[bin];
							stdFactor = chrModel->stdScaleFactor[bin];
							stdScaleFactor = chrModel->stdScaleFactor[bin];
						}
					} else {
						HiCBgModelChr* otherChrModel = (HiCBgModelChr*)model->chrs->search(chr2);
						if (otherChrModel == NULL) {
							expect *= chrModel->interChrScaleFactor;
						} else {
							expect *= sqrt(chrModel->interChrScaleFactor*otherChrModel->interChrScaleFactor);
						}
						stdFactor = chrModel->interStdScale;
						stdScaleFactor = chrModel->interStdScale;
					}
				}
				expectedArray[j] = expect;
				p2Stats[j] = interactions/expect;

				if (!(actionFlag & HIC_MASK_CREATEMODEL)) {
					if (actionFlag & HIC_MASK_RAWANDEXPECTED) {
						// just report raw reads
						p2Stats[j] = interactions;
					} else {
						//report ratio
						p2Stats[j] = log((HIC_RATIO_PSEUDO_COUNT+interactions)
												/(HIC_RATIO_PSEUDO_COUNT+expect))/log(2.0);
						if (stdScaleFactor > 0.0001) {
							zscore = p2Stats[j]/stdScaleFactor;
						}
						if (chrModel != NULL && actionFlag & HIC_MASK_NORM_ZSCORE) {
							p2Stats[j] = zscore;
						}
					}
				} else if (chrModel != NULL && chrModel->stdFlag) {
					//estimating variance while creating model
					p2Stats[j] = log((HIC_RATIO_PSEUDO_COUNT+interactions)
											/(HIC_RATIO_PSEUDO_COUNT+expect))/log(2.0);
				} else {
					//creating model - 
					if (peak2total > 0) {
						//p2Stats[j] = interactions/peak2total*peak1total;
						p2Stats[j] = interactions;
					} else {
						p2Stats[j] = 0;
					}
				}

				if (actionFlag & HIC_MASK_EXPECTED) {
					p2Stats[j]=expect;
				}
				if (actionFlag & HIC_MASK_LOGPVALUES) {
					double minTotal = peak1total;
					if (peak2total < minTotal) minTotal = peak2total;
					if (minTotal > 2.0) {
						double x = interactions;
						if (x > minTotal) x = minTotal;
						double ratio = 0;
						if (actionFlag & HIC_MASK_LOGPVALUESEXACT) {
							ratio = expect/minTotal;
						} else {
							// scale interactions to simple expect values - this should NOT be used EVER
							fprintf(stderr, "Don't use this option (normalizeing logp)- something is wrong...\n");
							x = interactions/expect*simpleExpect;
							//fprintf(stderr, "scaling: %lf\t%lf\t%lf\n", interactions,expect,x);
							ratio = simpleExpect/minTotal;
						}
						if (ratio > 1.0) ratio = 1.0;
						if (ratio < 0) ratio = 0.0;
						if (x > expect) {
							logp = logbinomialD((unsigned int)minTotal,(unsigned int)x,
												ratio,(unsigned int)minTotal);
						} else {
							logp = ilogbinomialD((unsigned int)minTotal,(unsigned int)x,
												ratio,(unsigned int)minTotal);
							logp *= -1;
						}
					}
					p2Stats[j]=logp;
				}
				if (actionFlag & HIC_MASK_LOCISCORE && (gil == NULL || distance >= gil->minDist)) {
					if (preFlag) preScore += p2Stats[j];
					if (postFlag) postScore += p2Stats[j];
					//if (logp < gil->threshold
                  	//	              && zscore > gil->zscoreThreshold) {
					//	totalPvalue += logp;
					//}
				}

				//if background, update interaction index to find previous interactions
				if (gil != NULL && gil->interactions != NULL) {

					while (gilInteractionIndex < gil->numInteractions && 
							gil->interactions[gilInteractionIndex]->peakIndex1 <= pIndex &&
							gil->interactions[gilInteractionIndex]->peakIndex2 < j) {
						gilInteractionIndex++;
					}
					if (gilInteractionIndex < gil->numInteractions && 
							gil->interactions[gilInteractionIndex]->peakIndex1 == pIndex &&
							gil->interactions[gilInteractionIndex]->peakIndex2 ==  j) {
						//Found an interaction!!!

						if (gil->recordInteractionBg) {
							gil->interactions[gilInteractionIndex]->totalPeak1Bg = peak1total;
							gil->interactions[gilInteractionIndex]->totalPeak2Bg = peak2total;
							gil->interactions[gilInteractionIndex]->interactionsBg = interactions;
							gil->interactions[gilInteractionIndex]->expectedBg = expect;
							gil->interactions[gilInteractionIndex]->logpBg = logp;
							gil->interactions[gilInteractionIndex]->zscoreBg = zscore;

							double logpdiffTg = 0.0;
							double logpdiffBg = 0.0;
							double logpdiff=0.0;


							// from the point of view of the primary experiment	
							double minTotalTg = gil->interactions[gilInteractionIndex]->totalPeak1;
							if (minTotalTg > gil->interactions[gilInteractionIndex]->totalPeak2) {
								minTotalTg = gil->interactions[gilInteractionIndex]->totalPeak2;
							}
	
							double interactionsTg = gil->interactions[gilInteractionIndex]->interactions;
							double expectTg = gil->interactions[gilInteractionIndex]->expected;
							double zscoreTg = gil->interactions[gilInteractionIndex]->zscore;
							double bgExpectInteractions = (interactionsTg-expectTg)/zscoreTg*zscore+expectTg;
							
							if (minTotalTg > 2.0) {
								double x = interactionsTg;
								if (x > minTotalTg) x = minTotalTg;
								double ratio = 0;
								ratio = bgExpectInteractions/minTotalTg;
								if (ratio > 1.0) ratio = 1.0;
								if (ratio < 0) ratio = 0.0;
								//fprintf(stderr, "%d\t%d\t%lf\t%lf\n", (unsigned int) minTotalTg, (unsigned int)x, ratio, (double)minTotalTg);
								if (x > bgExpectInteractions) {
									logpdiffTg = logbinomialD((unsigned int)minTotalTg,(unsigned int)x,
														ratio,(unsigned int)minTotalTg);
								} else {
									logpdiffTg = ilogbinomialD((unsigned int)minTotalTg,(unsigned int)x,
														ratio,(unsigned int)minTotalTg);
									logpdiffTg *= -1;
								}
							}

							// from the point of view of the background experiment	
							double minTotalBg = peak1total;
							if (minTotalBg > peak2total) {
								minTotalBg = peak2total;
							}
	
							double testInteractions = (interactions-expect)/zscore*zscoreTg+expect;
							
							if (minTotalBg > 2.0) {
								double x = testInteractions;
								if (x > minTotalBg) x = minTotalBg;
								double ratio = 0;
								ratio = interactions/minTotalBg;
								if (ratio > 1.0) ratio = 1.0;
								if (ratio < 0) ratio = 0.0;
								if (x > interactions) {
									logpdiffBg = logbinomialD((unsigned int)minTotalBg,(unsigned int)x,
														ratio,(unsigned int)minTotalBg);
								} else {
									logpdiffBg = ilogbinomialD((unsigned int)minTotalBg,(unsigned int)x,
														ratio,(unsigned int)minTotalBg);
									logpdiffBg *= -1;
								}
							}

							if (fabs(logpdiffBg) < fabs(logpdiffTg)) {
								logpdiff = logpdiffBg;
							} else {
								logpdiff = logpdiffTg;
							}
							//logpdiff = logpdiffTg;
							//gil->interactions[gilInteractionIndex]->expectedBg = logpdiffBg;

							gil->interactions[gilInteractionIndex]->logpDiff = logpdiff;
							if (logpdiff < gil->threshold 
									&& gil->interactions[gilInteractionIndex]->zscore-zscore > gil->zscoreThreshold) {
								gil->interactions[gilInteractionIndex]->thickness *= -1;
							}
	
						} else {
							int thickness = (int)((((double)logp)*-1.0)/5.0 * (10000.0/(10000.0+(double)resolution/3.0)));
							if (thickness < 1) thickness = 1;
							thickness*=2;

							int istart1 = (int)midpoint1-resolution/2;
							int iend1 = (int)midpoint1+resolution/2;
							int istart2 = (int)midpoint2-resolution/2;
							int iend2 = (int)midpoint2+resolution/2;
							gil->interactions[gilInteractionIndex]->start1 = istart1;
							gil->interactions[gilInteractionIndex]->end1 = iend1;
							gil->interactions[gilInteractionIndex]->start2 = istart2;
							gil->interactions[gilInteractionIndex]->end2 = iend2;
							gil->interactions[gilInteractionIndex]->totalPeak1 = peak1total;
							gil->interactions[gilInteractionIndex]->totalPeak2 = peak2total;
							gil->interactions[gilInteractionIndex]->interactions= interactions;
							gil->interactions[gilInteractionIndex]->expected = expect;
							gil->interactions[gilInteractionIndex]->zscore = zscore;
							gil->interactions[gilInteractionIndex]->logp = logp;
							gil->interactions[gilInteractionIndex]->thickness = thickness;

						}
					}
				} else if (actionFlag & HIC_MASK_INTERACTIONS && gil != NULL && gil->recordInteractions 
								&& logp < gil->threshold && distance >= gil->minDist 
								&& zscore > gil->zscoreThreshold) {
					if (peaks1 != peaks2 || chrcmp(&chr1,&chr2) < 0 
								|| (strcmp(chr1,chr2) != 0 || curStart < start1)) {
						//int thickness = (int)(((double)logp)*-1.0)/5.0;
						int thickness = (int)((((double)logp)*-1.0)/5.0 * (10000.0/(10000.0+(double)resolution/3.0)));
						if (thickness < 1) thickness = 1;
						thickness*=2;

						int istart1 = (int)midpoint1-resolution/2;
						int iend1 = (int)midpoint1+resolution/2;
						int istart2 = (int)midpoint2-resolution/2;
						int iend2 = (int)midpoint2+resolution/2;
						GenomeInteraction* gi = new GenomeInteraction(NULL,chr1,istart1,iend1,peak1total,
														pIndex,peakArray1[pIndex],
														chr2,istart2,iend2,peak2total,
														j,peakArray2[j],
														interactions,expect,logp,zscore,thickness);
						pthread_mutex_lock(&(gil->mutex));
						gil->addInteraction(gi);
						pthread_mutex_unlock(&(gil->mutex));
					}
				}
			}
//fprintf(stderr, "%s\t%d\t%s\t%d\t%lf\n", chr, start1,  chr2, curStart, p2Stats[j]);
			p2Total += p2Stats[j];

			if (gil != NULL && gil->fourCbedFile != NULL) {
				int ee= curEnd;
				if (curStart < ee) {
					pthread_mutex_lock(&(gil->mutex));
					fprintf(gil->fourCbedFile, "%s\t%d\t%d\t%lf\n", chr2, curStart, ee, p2Stats[j]);
					pthread_mutex_unlock(&(gil->mutex));
				}
			}

		}
		//if (p2TotalInteractions < 1.0) continue;

		if (actionFlag & HIC_MASK_LOCISCORE && gil->lociFile != NULL) {
			pthread_mutex_lock(&(gil->mutex));
			fprintf(gil->lociFile, "%s\t%s\t%d\t%d\t%d\t%lf\t%lf\n", 
							peakArray1[pIndex]->name,chr1,start1,end1,peakArray1[pIndex]->strand,preScore,postScore);
			pthread_mutex_unlock(&(gil->mutex));
		}
		//fprintf(stderr, "p2Total = %.1lf vs. %.1lf\n", p2Total, peak1total);

		if (actionFlag & HIC_MASK_CREATEMODEL && chrModel != NULL) {
			double w = peak1total;
			chrModel->wtotal+=w;
		}

		double upBscore = 0.0;
		double downBscore = 0.0;
		double upEscore = 0.0;
		double downEscore = 0.0;
		double upBN = 0.0;
		double downBN = 0.0;

		for (int j=0;j<numPeaks2;j++) {		

			char* chr2 = peakArray2[j]->chr;
			int curStart = peakArray2[j]->start;
			int curEnd = peakArray2[j]->end;

			double peak2total = peakArray2[j]->v;
			if (actionFlag & HIC_MASK_NORM_SEQDEPTH && peak2total < 0.5) continue;
			if (actionFlag & HIC_MASK_NORM_SEQDEPTH && model != NULL && (peakArray2[j]->v < model->minCoverage 
							|| peakArray2[j]->v > model->maxCoverage)) continue;
			//fprintf(stderr, "sub arrayHERE peaks %d %s %d %lf %lf %lf\n",i,peakArray1[i]->name,peakArray1[i]->refPos,p2Total,model->minCoverage, model->maxCoverage);

			int histIndex2 = -1;
			int histPosition2 = -1;
			if (numHistBins != 0) {
				sscanf(peakArray2[j]->name, "%d:%d",&histIndex2,&histPosition2);
				if (histIndex2 != histIndex1) continue;
			}

			if (actionFlag & HIC_MASK_BOUNDARIES) {
				if (strcmp(chr1,chr2)==0) {
					int diff = (curEnd+curStart)/2 - peakMidPoint;
					int minDistCheck = resolution;
					if (minDistCheck < params->minDist) {
						minDistCheck = params->minDist;
					}
					if (0) {
						if (diff < -1*resolution && diff > -1*boundaryScale) {
							downBscore += p2StatsB1[j];
							downEscore += expectedArray[j];
							downBN+=1.0;
						} else if (diff > resolution && diff < boundaryScale) {
							upBscore += p2StatsB1[j];
							upEscore += expectedArray[j];
							upBN+=1.0;
						} 
					} else {
						if (diff < -1*minDistCheck && diff > -1*boundaryScale) {
							if (params->logFlag) {
								downBscore += p2Stats[j];
							} else {
								downBscore += pow(2.0,p2Stats[j]);
							}
							//downEscore += expectedArray[j];
							downBN+=1.0;
						} else if (diff > minDistCheck && diff < boundaryScale) {
							if (params->logFlag) {
								upBscore += p2Stats[j];
							} else {
								upBscore += pow(2.0,p2Stats[j]);
							}
							//upEscore += expectedArray[j];
							upBN+=1.0;
						} 
					}
				}
			}

			if (actionFlag & HIC_MASK_CREATEMODEL && chrModel != NULL) {
				double w = peak1total*peak2total/totalInteractions;
				//double normV = p2Stats[j]/(peak2total*((double)numPeaks2)/peak1total);
				double normV = p2Stats[j]/(w*(double)numPeaks2);
				//double normV = p2Stats[j]/(peak1total*peak2total)*totalInteractions/(double)numPeaks2;
				//double normV = p2Stats[j]/peak1total;
				if (chrModel->stdFlag) {
					normV = p2Stats[j];
				}
				if (strcmp(chr1,chr2)==0) {
					int diff = (curEnd+curStart)/2 - peakMidPoint;
					if (diff < 0) diff *= -1;

					int x = (int)floor(((double)diff)/((double)resolution)+0.5);
					if (x < chrModel->maxIndex) {
						//fprintf(stderr, "%d\t%lf\t%lf\n", x ,normV, p2Total);
						if (chrModel->stdFlag) {
//							fprintf(stderr, "std %lf\t%lf\n",normV,1.0*w);
							chrModel->std[x] += normV*normV*w;
							chrModel->expectedN[x] += 1.0*w;
						} else {
							chrModel->expected[x] += normV*w;
							chrModel->expectedN[x] += 1.0*w;
						}
					} else {
						fprintf(stderr, "!! Potential problem, x is greater than max index\n");
					}
				} else {
					if (chrModel->stdFlag) {
						chrModel->interStd += normV*normV*w;
						chrModel->interChrN += 1.0*w;
					} else {
						chrModel->interChr += normV*w;
						chrModel->interChrN += 1.0*w;
					}
				}	
			}
			if (matrix != NULL) {
				double normV = p2Stats[j];
				if (numHistBins != 0) {
					int xx = histPosition1+halfBins;
					int yy = histPosition2+halfBins;

					//if (peak1total >= avgPeak1*0.5 && p2Totals[j] >= avgPeak1*0.5) {
	
					pthread_mutex_lock(mutexMatrix);
					if (revFlag) {
						matrix[yy][xx] += normV;
						matrix[yy+numHistBins][xx]+=1.0;
						if (expectMatrix != NULL) {
							expectMatrix[yy][xx] += expectedArray[j];
							//expectMatrix[yy+numHistBins][xx]+=1.0;
						}
					} else {
						matrix[xx][yy] += normV;
						matrix[xx+numHistBins][yy] += 1.0;
						if (expectMatrix != NULL) {
							expectMatrix[xx][yy] += expectedArray[j];
							//expectMatrix[xx+numHistBins][yy] += 1.0;
						}
					}
					pthread_mutex_unlock(mutexMatrix);
					//}
				} else {
					int jj = j;
					if (relativeFlag) {
						jj = j-pIndex+relativeFlag;
						if (jj < 0 || jj > relativeFlag*2) continue;
					}
					if (revFlag) matrix[jj][pIndex] += normV;
					else matrix[pIndex][jj] += normV;
					//fprintf(stderr, "pIndex=%d jj=%d normV=%.3lf\n",pIndex,jj,normV);
					if (expectMatrix != NULL) {
						if (revFlag) expectMatrix[jj][pIndex] += expectedArray[j];
						else expectMatrix[pIndex][jj] += expectedArray[j];
					}
				}
			}
		}

		if (actionFlag & HIC_MASK_BOUNDARIES) {
			// cool - lets nail down some boundaries...
			//
			//if (upBN > 0) upBscore /= upBN;
			//if (downBN > 0) downBscore /= downBN;
			if (upEscore > 0) upBscore /= upEscore;
			if (downEscore > 0) downBscore /= downEscore;
			//double bscore = downBscore;
			//double bscore = upBscore;
			double bscore = upBscore-downBscore;

			fprintf(stdout,"%s\t%d\t%d\t%lf\n", chr,start1,end1,bscore);
			
			//double minCountsForCorrelation = 5.0;
			//if (i>2) {
			//	double corr	= TagLibrary::calcCorrelation(p2StatsB1,p2Stats,numPeaks2,
			//						expectedArrayB1,expectedArray,minCountsForCorrelation,-1,0,0);
			//	if (corr < -1.9) corr = 0.0;
			//	fprintf(stdout,"%s\t%d\t%d\t%lf\n", chr,start1-resolution,end1-resolution,corr);
			//}



			//now lets save the raw tags for the next set of 3 regions
			double* tmpStats = p2StatsB1;
			p2StatsB1 = p2StatsB2;
			p2StatsB2 = p2Stats;
			p2Stats = tmpStats;

			double* tmpArray = expectedArrayB1;
			expectedArrayB1 = expectedArrayB2;
			expectedArrayB2 = expectedArray;
			expectedArray = tmpArray;

			PETag** tmpPETags = curPETagsB1;
			curPETagsB1 = curPETagsB2;
			curPETagsB2 = curPETags;
			curPETags = tmpPETags;

			//countB1 = countB2;
			//countB2 = count;

		}



	}
	//peakIndex2 = lastIndex2-1;

	if (histPeakIndexes != NULL) {
		for (int i=0;i<curNumPeaks1;i++) {
			if (histPeakIndexes[i] != NULL) delete [](histPeakIndexes[i]);
		}
		delete []histPeakIndexes;
		if (peakArray2 != NULL) delete []peakArray2;
	}
	if (histNumPeaks != NULL) delete []histNumPeaks;

	delete []iTags;
	delete []curPETags;
	delete []startTagIndex1;
	delete []endTagIndex1;
	delete []tagTotals;
	delete []p2Stats;
	delete []expectedArray;
	delete []p2Totals;

	if (actionFlag & HIC_MASK_BOUNDARIES) {
		delete []fineResProximity;
		delete []curPETagsB1;
		delete []curPETagsB2;
		delete []p2StatsB1;
		delete []p2StatsB2;
		delete []expectedArrayB1;
		delete []expectedArrayB2;
	}

	if (tagTotals2 == NULL) delete []tagTotals2;

	freeTags();
}


void ChrTags::decontaminate(ChrTags* input, int maxDistance, double fraction,double normFactor) {

	int iIndex = 0;
	int tIndex = 0;
	for (;iIndex < input->totalPositions;iIndex++) {
		int p = input->tags[iIndex].p;
		int d = input->tags[iIndex].d;
		float v = input->tags[iIndex].v;

		float bad = v*normFactor*fraction;
		if (bad <= TAGDIR_MIN_TAG_VALUE) continue;

		while (p > tags[tIndex].p && tIndex < totalPositions-1) tIndex++;

		int downIndex = tIndex-1;
		int upIndex = tIndex;
		while (bad > TAGDIR_MIN_TAG_VALUE) {
			int distDown = maxDistance*2;
			int distUp = maxDistance*2;
			if (downIndex >= 0 && downIndex < totalPositions)
				distDown = p - tags[downIndex].p;
			if (upIndex >= 0 && upIndex < totalPositions)
				distUp = tags[upIndex].p - p;

			if (distUp <= distDown) {
				if (distUp > maxDistance) break;
				if (tags[upIndex].d != d) {
					upIndex++;
					continue;
				}
				if (tags[upIndex].v > bad) {
					tags[upIndex].v -= bad;
					bad=0.0;
				} else {
					bad -= tags[upIndex].v;
					tags[upIndex].v = 0.0;
				}
				upIndex++;
			} else {
				if (distDown > maxDistance) break;
				if (tags[downIndex].d != d) {
					downIndex--;
					continue;
				}
				if (tags[downIndex].v > bad) {
					tags[downIndex].v -= bad;
					bad=0.0;
				} else {
					bad -= tags[downIndex].v;
					tags[downIndex].v = 0.0;
				}
				downIndex--;
			}	
		}
	}
}

//=========================================================================
// class HiCBgModelChr
HiCBgModelChr::HiCBgModelChr(int resolution) {
	res = resolution;
	init();
}
void HiCBgModelChr::init() {
	totalRegions = 0;
	initialize();
}
void HiCBgModelChr::initialize() {
	int maxChrLen = MAXIMUM_CHR_LENGTH;
	maxIndex = maxChrLen/res;
	maxUsedIndex = 0;
	if (maxIndex > 100000000) maxIndex = 100000000;
	expected = new double[maxIndex+1];
	expectedN = new double[maxIndex+1];
	std = new double[maxIndex+1];
	scaleFactor = new double[maxIndex+1];
	stdScaleFactor = new double[maxIndex+1];
	for (int i=0;i<maxIndex;i++) {
		expected[i]=0.0;
		expectedN[i]=0.0;
		std[i]=0.0;
		scaleFactor[i]=0.0;
		stdScaleFactor[i]=1.0;
	}
	interChr = 0.0;
	interChrScaleFactor = 0.0;
	interChrN = 0.0;
	interStd = 1.0;
	interStdScale = 1.0;
	stdFlag = 0;
	bgFlag = 0;
	wtotal = 0.0;
}
HiCBgModelChr::~HiCBgModelChr() {
	if (expected != NULL) delete []expected;
	if (expectedN != NULL) delete []expectedN;
	if (std != NULL) delete []std;
	if (scaleFactor != NULL) delete []scaleFactor;
	if (stdScaleFactor != NULL) delete []stdScaleFactor;
}
void HiCBgModelChr::smooth(int bins) {
	//smooth expected profile
	//int chalf = ((bins-1)+0.5)/2;
	double* smoothed = new double[maxIndex+1];
	for (int i=0;i<maxIndex;i++) {
		double avg = expected[i];
		double vavg = std[i];
		double NN = 1.0;
		int half = ((int)(log((double)i+1.0)/log(2.0))-1);
		for (int j=1;j<=half;j++) {
			if (i+j >= maxIndex) break;
			avg += expected[i+j];
			vavg += std[i+j];
			NN += 1.0;
		}
		for (int j=1;j<=half;j++) {
			if (i-j < 0) break;
			avg += expected[i-j];
			vavg += std[i-j];
			NN += 1.0;
		}
		scaleFactor[i] = avg/NN;
		smoothed[i] = vavg/NN;
		if (i>1 && scaleFactor[i] > scaleFactor[i-1]) {
			scaleFactor[i] = scaleFactor[i-1];
		}
	}
	delete []std;
	std = smoothed;
}

void HiCBgModelChr::normalize() {
	//fprintf(stderr, "\tNormalizing Chr\n");
	//normalize expected profile
	maxUsedIndex = 0;
	for (int i=0;i<maxIndex;i++) {
		if (expectedN[i] > 0) {
			//expected[i] /= wtotal;
			expected[i] /= expectedN[i];
			expectedN[i] = 0.0;
			maxUsedIndex = i+1;
		} else {
			expected[i] = 0.0;
			expectedN[i] = 0.0;
		}
	}
	//interChr /= wtotal;
	if (interChrN > 0.0) interChr /= interChrN;
	interChrN = 0.0;
}
void HiCBgModelChr::normalizeStd() {
	//normalize expected profile
	for (int i=0;i<maxIndex;i++) {
		if (expectedN[i] > 0) {
			//std[i] /= wtotal;
			std[i] /= expectedN[i];
			std[i] = sqrt(std[i]);
			expectedN[i] = 0.0;
		} else {
			std[i] = 0.0;
			expectedN[i] = 0.0;
		}
	}
	if (interChrN > 0.0) {
		//interStd /= wtotal;
		interStd /= interChrN;
		interStd = sqrt(interStd);
	}
	interChrN = 0.0;
}

// class HiCBgModel
HiCBgModel::HiCBgModel(int resolution,char* dir, int nfileFlag) : HiCBgModelChr(resolution) {
	init();
	customFlag = nfileFlag;
	if (res < 100) {
		fprintf(stderr, "!!! Warning, resolution less than 100 bp...probably not a good idea!!!\n");
	}
	directory = new char[strlen(dir)+1];
	strcpy(directory,dir);
	if (customFlag) {
		filename = new char[strlen(dir)+1];
		strcpy(filename, directory);
	} else {
		filename = getDefaultFileName(directory);
	}
	if (customFlag >1) return;
	ready = load();
}
void HiCBgModel::init() {
	ready=0;
	customFlag = 0;
	filename = NULL;
	directory = NULL;
	refPeaks = NULL;
	chrExpect = new Doubletable(1000);
	chrSizes = new Doubletable(1000);
	chrs = new Hashtable(1000);
	avgCoverage = 0.0;
	stdCoverage = 1.0;
	stdFilter = 4.0;
	minFilter = 0.5;
	minCoverage = 0.0;
	maxCoverage = 0.0;
	goodRegions = 0;
	badRegions = 0;
	totalModelReads = 0.0;
	modelError=0.0;
	fullModelFlag = 0;
	numPeaks1=0;
	numPeaks2=0;
	expectMatrix=NULL;
}
void HiCBgModel::initializeBasedOnChrSizes() {
	HiCBgModelChr::initialize();
	char** keys = chrSizes->keys();
	qsort(keys,chrSizes->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrSizes->total;i++) {
		HiCBgModelChr* cmodel = (HiCBgModelChr*) chrs->search(keys[i]);
		if (cmodel == NULL) {
			cmodel = addChrModel(keys[i]);
		}
		if (cmodel != NULL) cmodel->initialize();
		delete [](keys[i]);
	}
	delete []keys;

}
void HiCBgModel::initialize() {
	HiCBgModelChr::initialize();
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		HiCBgModelChr* cmodel = (HiCBgModelChr*) chrs->search(keys[i]);
		if (cmodel != NULL) cmodel->initialize();
		delete [](keys[i]);
	}
	delete []keys;

}
HiCBgModel::~HiCBgModel() {
	if (filename != NULL) delete []filename;
	if (directory != NULL) delete []directory;
	if (chrSizes != NULL) delete chrSizes;
	if (chrExpect != NULL) delete chrExpect;
	if (chrs != NULL) {
		char** keys  = chrs->keys();
		for (int i=0;i<chrs->total;i++) {
			HiCBgModelChr* cbg = (HiCBgModelChr*) chrs->search(keys[i]);
			if (cbg != NULL) delete cbg;
			delete [](keys[i]);
		}
		delete []keys;
		delete chrs;
	}
	if (expectMatrix != NULL) {
		for (int i=0;i<numPeaks1;i++) {
			if (expectMatrix[i] != NULL) {
				delete [](expectMatrix[i]);
			}
		}
		delete []expectMatrix;
		expectMatrix = NULL;
	}
	if (refPeaks != NULL) delete refPeaks;
}
void HiCBgModel::initializeExpectMatrix(int n1, int n2) {
	numPeaks1 = n1;
	numPeaks2 = n2;
	expectMatrix = new double*[numPeaks1];
	for (int i=0;i<numPeaks1;i++) {
		expectMatrix[i] = new double[numPeaks2];
		for (int j=0;j<numPeaks2;j++) {
			expectMatrix[i][j] = 0.0;
		}
	}
}
char* HiCBgModel::getDefaultFileName(char* dir) {
	char* f = new char[strlen(dir)+10000];
	sprintf(f,"%s/HiCbackground_%d_bp.txt",dir,res);
	return f;
}
HiCBgModelChr* HiCBgModel::addChrModel(char* chr) {
	HiCBgModelChr* c = new HiCBgModelChr(res);
	chrs->insert(c,chr);
	return c;
}
void HiCBgModel::setCoverageLimits(double a, double s) {
	avgCoverage = a;
	stdCoverage = s;
	minCoverage = a-s*stdFilter;
	maxCoverage = a+s*stdFilter;
	double altMin = a*minFilter;
	if (minCoverage < altMin) minCoverage = altMin;
	if (minCoverage < 0.5) minCoverage = 0.5;
}

int HiCBgModel::load() {
	//fprintf(stderr, "\tloading background model...\n");
	FILE* fp = fopen(filename,"r");
	if (fp == NULL) {
		//fprintf(stderr, "\tCouldn't find file: %s\n", filename);
		return 0;
	} else {
		//fprintf(stderr, "\tFound file: %s\n", filename);
	}

	char* buffer = new char[BUFFER];
	char** cols = new char*[10000];
	int numCols = 0;
	int pos = 0;
	double p = 0.0;
	double pstd = 0.0;
	char** chrNames = NULL;
	int totalChrNames = 0;
	int pstart = 0;
	int pend = 0;
	float ptotal = 0.0;
	float padj = 0.0;
	int refPeakFlag = 0;

	while (fgets(buffer,BUFFER,fp)!=NULL) {
		split(buffer,cols,numCols,'\t');
		if (cols[0][0] == '#') continue;
		if (refPeakFlag) {
			sscanf(cols[2],"%d",&pstart);
			sscanf(cols[3],"%d",&pend);
			sscanf(cols[5],"%f",&ptotal);
			sscanf(cols[6],"%f",&padj);
			refPeaks->addPeak(cols[0],cols[1],pstart,pend,(pend+pstart)/2,
								STRAND_POSITIVE,ptotal,padj,NULL,-1,0);
			totalModelReads += padj;
			continue;
		}
		if (strcmp("Adjusted Region Totals",cols[0])==0) {
			refPeakFlag=1;
			refPeaks = new PeakLibrary(totalRegions);
			continue;
		}
		if (strcmp("FullModel=true",cols[0])==0) {
			fullModelFlag = 1;
			continue;
		}
		if (strcmp("FullModel=false",cols[0])==0) {
			fullModelFlag = 0;
			continue;
		}
		if (strcmp("ModelError=",cols[0])==0) {
			sscanf(cols[1],"%lf",&modelError);
			continue;
		}

		if (strcmp("Resolution=",cols[0])==0) {
			if (numCols < 2) {
				fprintf(stderr, "\t!!! Looks like you could have an older version of a background file (going to remake)...\n");
				return 0;
			}
			int fileRes = 0;
			sscanf(cols[1],"%d",&fileRes);
			if (fileRes != res) {
				fprintf(stderr, "\n!!!!!!!!!!!!!!!!!\n");
				fprintf(stderr, "\t Warning, background model has different resolution (%d vs. %d [file])!!!\n",
								res, fileRes);
				fprintf(stderr, "!!!!!!!!!!!!!!!!!\n\n");
			}
		} else if (strcmp("TotalRegions=",cols[0])==0) {
			sscanf(cols[1],"%d",&totalRegions);
		} else if (strcmp("TotalModelReads=",cols[0])==0) {
			sscanf(cols[1],"%lf",&totalModelReads);
		} else if (strncmp("Inter",cols[0],5)==0) {
			char* str = &(cols[0][5]);
			if (numCols > 2) {
				sscanf(cols[1],"%le",&p);
				sscanf(cols[2],"%le",&pstd);
				if (strcmp(str,"genome")==0) {
					interChr = p;
					interStd = pstd;
					interChrN = 1;
				} else {
					HiCBgModelChr* cc = addChrModel(str);
					cc->interChr = p;
					cc->interStd = pstd;
					cc->interChrN = 1;
				}
			} else {
				fprintf(stderr, "\t!!! Looks like you could have an older version of a background file (going to remake)...\n");
				return 0;
			}
		} else if (strncmp("Size",cols[0],4)==0) {
			char* str = &(cols[0][4]);
			if (numCols > 1) {
				sscanf(cols[1],"%lf",&p);
				chrSizes->insert(p,str);
			}
		} else if (numCols > 4) {
			if (strcmp(cols[0],"Index")==0) {
				chrNames = new char*[numCols];
				totalChrNames = numCols;
				for (int i=0;i<numCols;i++) {
					chrNames[i] = new char[strlen(cols[i])+1];
					strcpy(chrNames[i],cols[i]);
				}
			}
			if (cols[0][0] < 45 || cols[0][0] > 57) continue;
			sscanf(cols[0],"%d",&pos);
			sscanf(cols[2],"%le",&p);
			sscanf(cols[3],"%lf",&pstd);
			if (pos >= maxIndex) {
				fprintf(stderr, "!!! Out of range !!! pos = %d maxIndex = %d (resolution set at = %d Is this right?)\n",pos,maxIndex,res);
				continue;
			}
			expected[pos] = p;
			std[pos] = pstd;
			expectedN[pos] = 1.0;
			if (pos >= maxUsedIndex) maxUsedIndex = pos+1;
			for (int i=4;i<numCols;i+=2) {
				sscanf(cols[i],"%le",&p);
				sscanf(cols[i+1],"%le",&pstd);
				HiCBgModelChr* cc = (HiCBgModelChr*)chrs->search(chrNames[i]);
				if (cc != NULL) {
					if (pos > cc->maxIndex) {
						fprintf(stderr, "!!! Something is probably wrong with the background file - delete it!\n");
					}
					cc->expected[pos] = p;
					cc->std[pos] = pstd;
					cc->expectedN[pos] = 1.0;
					cc->maxUsedIndex = pos+1;
				}
			}
		}
	}
	if (refPeakFlag) {
		refPeaks->sortChr();
		refPeaks->setDefaultPeakOrder();
	}

	fclose (fp);
	delete []buffer;
	delete []cols;
	if (chrNames != NULL) {
		for (int i=0;i<totalChrNames;i++) {
			if (chrNames[i] != NULL) delete [](chrNames[i]);
		}
		delete []chrNames;
	}
	return 1;
}
void HiCBgModel::save(char* file) {
	char* tmpFileName = filename;
	filename = file;
	save();
	filename = tmpFileName;
}
void HiCBgModel::save() {
	FILE* fp = fopen(filename,"w");
//	FILE* fp = fopen("tmp.txt","w");
	if (fullModelFlag) {	
		fprintf(fp,"FullModel=true\n");
	} else {
		fprintf(fp,"FullModel=false\n");
	}
	
	fprintf(fp,"Resolution=\t%d\n",res);
	fprintf(fp,"TotalRegions=\t%d\n",totalRegions);
	fprintf(fp,"TotalModelReads=\t%lf\n",totalModelReads);
	fprintf(fp,"ModelError=\t%lf\n",modelError);
	fprintf(fp,"\n");

	char** keys = chrSizes->keys();
	qsort(keys,chrSizes->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrSizes->total;i++) {
		fprintf(fp,"Size%s\t%.0lf\n",keys[i],chrSizes->search(keys[i]));
		delete [](keys[i]);
	}
	delete []keys;

	fprintf(fp,"\n");

	keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		HiCBgModelChr* hcbm = (HiCBgModelChr*) chrs->search(keys[i]);
		fprintf(fp,"Inter%s\t%le\t%lf\n",keys[i],hcbm->interChr,hcbm->interStd);
		delete [](keys[i]);
	}
	delete []keys;
	fprintf(fp,"Intergenome\t%le\t%lf\n",interChr,interStd);

	fprintf(fp,"\n");
	fprintf(fp,"Index\tDistance\tExpected Interaction Frequency\tObs/Exp Signal Std");
	keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int j=0;j<chrs->total;j++) {
		fprintf(fp, "\t%s\t%s Signal Std ",keys[j],keys[j]);
	}
	fprintf(fp,"\n");
	for (int i=0;i<maxUsedIndex;i++) {
		fprintf(fp,"%d\t%d\t%le\t%.6lf",i,i*res,expected[i],std[i]);
		for (int j=0;j<chrs->total;j++) {
			HiCBgModelChr* hcbm = (HiCBgModelChr*) chrs->search(keys[j]);
			fprintf(fp, "\t%le\t%.6lf",hcbm->expected[i],hcbm->std[i]);
		}
		fprintf(fp,"\n");
	}
	for (int i=0;i<chrs->total;i++) {
		delete [](keys[i]);
	}
	delete []keys;

	if (refPeaks != NULL) {
		fprintf(fp, "\n\nAdjusted Region Totals\tchr\tstart\tend\tstrand\tPETag total\tScaling Factor\n");
		for (int i=0;i<refPeaks->numPeaks;i++) {
			Peak** peakOrder = refPeaks->peakOrder;
			fprintf(fp, "%s\t%s\t%d\t%d\t+\t%f\t%f\n",peakOrder[i]->name,peakOrder[i]->chr, peakOrder[i]->start,
													peakOrder[i]->end,peakOrder[i]->v,peakOrder[i]->focusRatio);
		}
	}

	fclose(fp);
}
void HiCBgModel::setDefaultVariation() {
	fprintf(stderr, "\tSetting defaults for variation using approximate background...\n");
	for (int i=0;i<maxIndex;i++) {
		expectedN[i] = 0.0;
		std[i] = 0.0;
	}
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		HiCBgModelChr* hcbm = (HiCBgModelChr*) chrs->search(keys[i]);
		for (int j=0;j<hcbm->maxIndex;j++) {
			hcbm->std[j] = 1.0;
			hcbm->expectedN[j] = 0.0;
			std[j] = 1.0;
			expectedN[j] = 0.0;
		}
		hcbm->interStd = 1.0;
		interStd = 1.0;
		interChrN = 0.0;
		delete [](keys[i]);
	}
	delete []keys;
}


void HiCBgModel::normalizeStd() {
	fprintf(stderr, "\tnormalizing std...\n");
	for (int i=0;i<maxIndex;i++) {
		expectedN[i] = 0.0;
		std[i] = 0.0;
	}
	interChrN = 0.0;
	interStd = 0.0;

	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		HiCBgModelChr* hcbm = (HiCBgModelChr*) chrs->search(keys[i]);
		for (int j=0;j<hcbm->maxIndex;j++) {
			std[j] += hcbm->std[j];
			expectedN[j] += hcbm->expectedN[j];
		}
		interStd += hcbm->interStd;
		interChrN += hcbm->interChrN;
		hcbm->normalizeStd();
		delete [](keys[i]);
	}
	delete []keys;
	HiCBgModelChr::normalizeStd();
}

void HiCBgModel::normalize() {
	fprintf(stderr, "\tNormalizing...\n");
	maxUsedIndex = 0;
	for (int i=0;i<maxIndex;i++) {
		expected[i] = 0.0;
		expectedN[i] = 0.0;
	}
	interChr = 0.0;
	interChrN = 0.0;
	wtotal = 0.0;

	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		HiCBgModelChr* hcbm = (HiCBgModelChr*) chrs->search(keys[i]);
		for (int j=0;j<hcbm->maxIndex;j++) {
			expected[j] += hcbm->expected[j];
			expectedN[j] += hcbm->expectedN[j];
		}
		wtotal += hcbm->wtotal;
		interChr += hcbm->interChr;
		interChrN += hcbm->interChrN;
		//fprintf(stderr, "\t%lf\t%lf\n", hcbm->interChr, hcbm->interChrN);
		hcbm->normalize();
		delete [](keys[i]);
	}
	delete []keys;
	HiCBgModelChr::normalize();
}

void HiCBgModel::useApproximation(double *dist, int arrayLength) {

	initializeBasedOnChrSizes();
	int lastNonZeroIndex = 0;
	for (int i=arrayLength-2;i>=0;i--) {
		if (dist[i] > 0) {
			lastNonZeroIndex = i+1;
			break;
		}
	}
	interChr = dist[arrayLength-1];
	double numChr = (double)chrs->total;
	if (numChr < 2.0) numChr = 2.0;
	interChr /= ((double)totalRegions*((numChr-1.0)/numChr));
	
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	int maxMIndex = 0;
	for (int i=0;i<chrs->total;i++) {
		HiCBgModelChr* hcbm = (HiCBgModelChr*) chrs->search(keys[i]);
		double size = chrSizes->search(keys[i]);
		if (size < EMPTY_DOUBLE) {
			size = 1000000000;
		}
		int mIndex = ((int)(size/((double)res)))+1;
		if (mIndex > lastNonZeroIndex) {
			mIndex = lastNonZeroIndex;
		}
		if (mIndex > hcbm->maxIndex) {
			mIndex = hcbm->maxIndex;
		}
		if (mIndex > maxMIndex) maxMIndex = mIndex;
		//fprintf(stderr, "%d=maxIndex %s %lf %d hcbmMax = %d\n",mIndex,keys[i],size,lastNonZeroIndex,hcbm->maxIndex);
		for (int j=0;j<mIndex;j++) {
			hcbm->expected[j] = dist[j];
			hcbm->expectedN[j] = 0.0;
		}
		hcbm->maxUsedIndex = mIndex;
		hcbm->wtotal = 1;
		hcbm->interChr = interChr;

		//fprintf(stderr, "\t%lf\t%lf\n", hcbm->interChr, hcbm->interChrN);
		//hcbm->normalize();
		delete [](keys[i]);
	}
	delete []keys;

	for (int j=0;j<maxMIndex;j++) {
		expected[j] = dist[j];
		expectedN[j] = 0.0;
	}
	maxUsedIndex = maxMIndex;
	//HiCBgModelChr::normalize();
}


void HiCBgModel::normalizeVariance(HiCBgModel* m) {
	for (int i=0;i<maxUsedIndex;i++) { 
		if (i >= m->maxUsedIndex || m->std[i] < 0.001) {
			stdScaleFactor[i] = 1.0;
		} else {
			stdScaleFactor[i] = std[i]/m->std[i];
		}
	}
	if (m->interStd < 0.001) {
		interStdScale = 1.0;
	} else {
		interStdScale = interStd/m->interStd;
	}
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		HiCBgModelChr* hcbm = (HiCBgModelChr*) chrs->search(keys[i]);
		for (int j=0;j<maxUsedIndex;j++) { 
			hcbm->stdScaleFactor[j] = stdScaleFactor[j];
		}
		hcbm->interStdScale = interStdScale;
		delete [](keys[i]);
	}
	delete []keys;
}
void HiCBgModel::scale() {
//fprintf(stderr, "test=%lf\n", expected[0]);
	HiCBgModelChr::smooth(5);
	fprintf(stderr, "\tTotal regions in background model=%d\n", totalRegions);
	char** keys = chrs->keys();
	qsort(keys,chrs->total,sizeof(char*),&chrcmp);
	for (int i=0;i<chrs->total;i++) {
		HiCBgModelChr* hcbm = (HiCBgModelChr*) chrs->search(keys[i]);

		double chrTotal = 0;
		double avgTotal = 0;
		for (int j=0;j<hcbm->maxUsedIndex;j++) {
			chrTotal += hcbm->expected[j];
			avgTotal += expected[j];
		}
		double normFactor = 1.0;
		//int mmax = maxUsedIndex;
		hcbm->interChrScaleFactor = interChr;
		if (chrTotal > 1e-10 && avgTotal > 1e-10) {
			normFactor = chrTotal/avgTotal;
			//mmax = hcbm->maxUsedIndex;
			hcbm->interChrScaleFactor = hcbm->interChr;
		}
		for (int j=0;j<hcbm->maxUsedIndex;j++) {
			hcbm->scaleFactor[j] = scaleFactor[j]*normFactor;
			//hcbm->scaleFactor[j] = hcbm->expected[j];
			hcbm->stdScaleFactor[j] = std[j];
		}
		hcbm->interStdScale = interStd;
		delete [](keys[i]);
	}
	interChrScaleFactor = interChr;
	delete []keys;
}
void HiCBgModel::createRandomReads(FILE* fp, int totalReadsInt) {
	if (fp == NULL) return;

	scale();

	double totalReads = (double)totalReadsInt;
	Peak** peaks = refPeaks->peakOrder;
	int numPeaks = refPeaks->numPeaks;

	int maxReadsPerPeak = 0;
	double totalPEReads = 0;		
	double totalAdjReads = 0;		
	for (int i=0;i<numPeaks;i++) {
		totalPEReads += peaks[i]->v;
		totalAdjReads += peaks[i]->v*peaks[i]->focusRatio;
	}
	for (int i=0;i<numPeaks;i++) {
		int m = (int)(peaks[i]->v/totalPEReads*totalReads*2);
		if (m > maxReadsPerPeak) maxReadsPerPeak = m;
	}

	fprintf(stderr, "\t--------------------------------------\n");
	fprintf(stderr, "\tRandomizing PE reads according to Bg model: %s\n", filename);
	fprintf(stderr, "\tTotal Reads from Model: %.1lf\n", totalPEReads);
	fprintf(stderr, "\tTotal Adj Reads from Model: %.1lf\n", totalAdjReads);

	double* expectedReads = new double[numPeaks];
	double* peakCDF = new double[numPeaks];
	double* randNumbers = new double[maxReadsPerPeak];
	char* lastChr= NULL;
	HiCBgModelChr* chrModel = this;
	int startChrIndex = 0;
	int endChrIndex = 0;
	double interChrSum = 0.0;
	int readName = 0;

	for (int i=0;i<numPeaks;i++) {
		int peakReads = ((int)(floor(peaks[i]->v/totalPEReads*totalReads+0.5)))/2+1;
		if (peakReads > maxReadsPerPeak) {
			fprintf(stderr, "!!! Something is wrong here.... %d (no more than %d expected)\n",
									peakReads, maxReadsPerPeak);
			exit(0);
		}


		if (lastChr == NULL || strcmp(lastChr,peaks[i]->chr) != 0) {
			lastChr = peaks[i]->chr;
			startChrIndex = INT_MAX;
			endChrIndex = 0;
			chrModel = (HiCBgModelChr*)chrs->search(peaks[i]->chr);	
			if (chrModel == NULL) chrModel = this;
			interChrSum = 0.0;
			for (int j=0;j<numPeaks;j++) {
				if (strcmp(lastChr,peaks[j]->chr) == 0) {
					if (j < startChrIndex) startChrIndex = j;
					if (j > endChrIndex) endChrIndex = j;
					expectedReads[j] = 0; // do this later
				} else {
					expectedReads[j] = peaks[j]->v*peaks[j]->focusRatio*chrModel->interChrScaleFactor;
				}
				interChrSum += expectedReads[j];
			}
		}

		double expectedSum = interChrSum;
		for (int j=startChrIndex;j<=endChrIndex;j++) {
			int bin = j-i;
			if (bin < 0) bin *= -1;
			expectedReads[j] = peaks[j]->v*peaks[j]->focusRatio*chrModel->scaleFactor[bin];
			expectedSum += expectedReads[j];
		}
		
		peakCDF[0] = expectedReads[0]/expectedSum;
		for (int j=1;j<numPeaks;j++) {
			peakCDF[j] = expectedReads[j]/expectedSum+peakCDF[j-1];
		}

		for (int j=0;j<peakReads;j++) {
			double r = ((double)rand())/((double)RAND_MAX);
			randNumbers[j] = r;
		}
		qsort(randNumbers,peakReads,sizeof(double),&cmpDouble);

		int pIndex = 0;
		int peakLength = peaks[i]->end - peaks[i]->start;
		for (int j=0;j<peakReads;j++) {
			double r = randNumbers[j];
			while (r > peakCDF[pIndex] && pIndex < numPeaks-1) pIndex++;
			int start = peaks[pIndex]->start;
			int end = peaks[pIndex]->end;
			double frac = 0.0;
			if (pIndex == 0) {
				frac = r/(peakCDF[pIndex]-0.0);
			} else {
				double den = (peakCDF[pIndex]-peakCDF[pIndex-1]);
				if (den < 1e-100) {
					frac = 0.5;
				} else {
					frac = (r-peakCDF[pIndex-1])/den;
				}
			}
			if (frac < 0.0) frac = 0.0;
			if (frac > 1.0) frac = 1.0;

			int start1 = (int)(frac*(end-start))+start;
			r = ((double)rand())/((double)RAND_MAX);
			int start2 = peaks[i]->start+(int)(r*peakLength);

			int d1 = 0;
			if (start2 % 2 == 0) d1 = 1;
			int d2 = 0;
			if (start1 % 2 == 0) d2 = 1;
		
			fprintf(fp, "rand%d\t%s\t%d\t%d\t1.0\t1\t%s\t%d\t%d\t1\n",++readName,peaks[pIndex]->chr,
											start1,d1,peaks[i]->chr,start2,d2);
			fprintf(fp, "rand%d\t%s\t%d\t%d\t1.0\t1\t%s\t%d\t%d\t1\n",++readName,peaks[i]->chr,
											start2,d2,peaks[pIndex]->chr,start1,d1);
			
		}

	}

	delete []expectedReads;
	delete []peakCDF;
	delete []randNumbers;

}

int cmpDouble(const void* a, const void* b) {
	double aa = *((double*)a);
	double bb = *((double*)b);
	if (aa < bb) return -1;
	if (aa > bb) return 1;
	return 0;
}


GenomeInteraction::GenomeInteraction() {
	init();
}
GenomeInteraction::GenomeInteraction(char* iname, 
					char* ichr1, int istart1, int iend1, double total1, int pIndex1, Peak* p1,
					char* ichr2, int istart2, int iend2, double total2, int pIndex2, Peak* p2,
					double iinteractions,double iexpected, double ilogp,double izscore, int ithickness) {
	init();
	if (iname != NULL) {
		name = new char[strlen(iname)+1];
		strcpy(name,iname);
	}
	chr1 = new char[strlen(ichr1)+1];
	strcpy(chr1,ichr1);
	chr2 = new char[strlen(ichr2)+1];
	strcpy(chr2,ichr2);
	start1 = istart1;
	end1 = iend1;
	start2 = istart2;
	end2 = iend2;
	interactions = iinteractions;
	expected = iexpected;
	logp = ilogp;
	zscore = izscore;
	fdr = 1.0;
	peak1 = p1;
	peak2 = p2;
	thickness = ithickness;
	peakIndex1 = pIndex1;
	totalPeak1 = total1;
	peakIndex2 = pIndex2;
	totalPeak2 = total2;
}
void GenomeInteraction::setName(char* newname) {
	if (newname == NULL) return;
	if (name != NULL) delete []name;
	name = new char[strlen(newname)+1];
	strcpy(name,newname);
}
GenomeInteraction::GenomeInteraction(Peak* np1, Peak* np2) {
	init();
	peak1 = np1;
	peak2 = np2;
	chr1 = new char[strlen(peak1->chr)+1];
	strcpy(chr1,peak1->chr);
	chr2 = new char[strlen(peak2->chr)+1];
	strcpy(chr2,peak2->chr);
	start1 = peak1->start;
	end1 = peak1->end;
	start2 = peak2->start;
	end2 = peak2->end;
}	

void GenomeInteraction::init() {
	name = NULL;
	chr1 = NULL;
	start1 = -1;
	end1 = -1;
	chr2 = NULL;
	start2 = -1;
	end2 = -1;
	totalPeak1 = 0;
	totalPeak2 = 0;
	peakIndex1 = 0;
	peakIndex2 = 0;
	peak1 = NULL;
	peak2 = NULL;

	logp = 0.0;
	fdr = 1.0;
	zscore = 0.0;

	logpBg = 0.0;
	fdrBg = 1.0;
	zscoreBg = 0.0;

	logpDiff = 0.0;
	fdrDiff = 0.0;

	totalPeak1Bg = 0;
	totalPeak2Bg = 0;
	interactionsBg = 0.0;
	expectedBg = 0.0;
}
GenomeInteraction::~GenomeInteraction() {
	if (chr1 != NULL) delete []chr1;
	if (chr2 != NULL) delete []chr2;
	if (name != NULL) delete []name;
}
void GenomeInteraction::printCircos(FILE* fp,int count) {
	fprintf(fp,"interaction%d %s %d %d",count,chr1,start1,end1);
	fprintf(fp," thickness=%d",thickness);
	fprintf(fp,"\n");
	fprintf(fp,"interaction%d %s %d %d",count,chr2,start2,end2);
	fprintf(fp," thickness=%d",thickness);
	fprintf(fp,"\n");
}
void GenomeInteraction::printPETags(FILE* fp, PETag** tags, int numTags,int format) {
	for (int i=0;i<numTags;i++) {
		if (format == 1) {
			tags[i]->print(fp,0);	
			tags[i]->print(fp,1);	
		} else {

			if (strcmp(tags[i]->chr1,tags[i]->chr2)==0) {
				fprintf(fp,"%s",tags[i]->chr1);
				int s1 = tags[i]->p1;
				int l1 = tags[i]->len1;
				int s2 = tags[i]->p2;
				int l2 = tags[i]->len2;
				if (tags[i]->d1 == STRAND_NEGATIVE) s1 -= l1;
				if (tags[i]->d2 == STRAND_NEGATIVE) s2 -= l2;
				if (s1 > s2) {
					int tmp = s1;
					int tmp2 = l1;
					s1 = s2; l1 = l2; s2 = tmp; l2 = tmp2;
				}
				if (s1+l1 > s2) {
					int m = (s1+l1 + s2)/2;
					l1 = m-s1;
					s2 = m;
					l2 = l1;
				}
				fprintf(fp, "\t%d\t%d\t%s:%d:",s1,s2+l2,tags[i]->chr1,tags[i]->p1);
				if (tags[i]->d1 == STRAND_POSITIVE) fprintf(fp,"+");
				else fprintf(fp, "-");
				fprintf(fp, "to%s:%d:",tags[i]->chr2,tags[i]->p2);
				if (tags[i]->d2 == STRAND_POSITIVE) fprintf(fp,"+");
				else fprintf(fp, "-");
				
				fprintf(fp, "\t%f\t+\t%d\t%d\t255,0,0\t2\t%d,%d\t%d,%d\n",tags[i]->v,s1,s2+l2,l1,l2,0,s2-s1);
				
			} else {
				
				int s1 = tags[i]->p1;
				int l1 = tags[i]->len1;
				if (tags[i]->d1 == STRAND_NEGATIVE) s1 -= l1;
			
				fprintf(fp, "%s\t%d\t%d",tags[i]->chr1,s1,s1+l1);
				fprintf(fp, "\t%s:%d:",tags[i]->chr1,tags[i]->p1);
				if (tags[i]->d1 == STRAND_POSITIVE) fprintf(fp,"+");
				else fprintf(fp, "-");
				fprintf(fp, "to%s:%d:",tags[i]->chr2,tags[i]->p2);
				if (tags[i]->d2 == STRAND_POSITIVE) fprintf(fp,"+");
				else fprintf(fp, "-");
	
				fprintf(fp, "\t%f+\t%d\t%d\t0,255,0\t1\t%d\t%d\n",tags[i]->v,s1,s1+l1,l1,0);
	
				int s2 = tags[i]->p2;
				int l2 = tags[i]->len2;
				if (tags[i]->d2 == STRAND_NEGATIVE) s2 -= l2;
			
				fprintf(fp, "%s\t%d\t%d",tags[i]->chr2,s2,s2+l2);
				fprintf(fp, "\t%s:%d:",tags[i]->chr1,tags[i]->p1);
				if (tags[i]->d1 == STRAND_POSITIVE) fprintf(fp,"+");
					else fprintf(fp, "-");
				fprintf(fp, "to%s:%d:",tags[i]->chr2,tags[i]->p2);
				if (tags[i]->d2 == STRAND_POSITIVE) fprintf(fp,"+");
				else fprintf(fp, "-");

				fprintf(fp, "\t%f+\t%d\t%d\t0,255,0\t1\t%d\t%d\n",tags[i]->v,s2,s2+l2,l2,0);
			}
		}
	}
}
void GenomeInteraction::print(FILE* fp,int count, int bgFlag, int centerFlag) {

	if (name != NULL) {
		fprintf(fp, "%s",name);
	} else {
		fprintf(fp,"interaction%d",count);
	}
	fprintf(fp, "\t%s\t%s",peak1->name,chr1);

	if (centerFlag) {
		fprintf(fp, "\t%d\t%d",start1,end1);
	} else {
		fprintf(fp, "\t%d\t%d",peak1->start,peak1->end);
	}
	if (peak1->strand == STRAND_NEGATIVE) fprintf(fp, "\t-");
	else fprintf(fp, "\t+");
	fprintf(fp, "\t%.1lf\t%s\t%s", totalPeak1,peak2->name,chr2);
	if (centerFlag) fprintf(fp, "\t%d\t%d",start2,end2);
	else fprintf(fp, "\t%d\t%d",peak2->start,peak2->end);
	if (peak2->strand == STRAND_NEGATIVE) fprintf(fp, "\t-");
	else fprintf(fp, "\t+");
	fprintf(fp, "\t%.1lf", totalPeak2);

	if (strcmp(chr1,chr2)==0) {
		int distance = start1-start2;
		if (distance < 0) distance *= -1;
		fprintf(fp,"\t%d",distance);
	} else {
		fprintf(fp,"\tinterchromosomal");
	}
	
	fprintf(fp, "\t%lf\t%lf\t%lf\t%lf\t%lf",interactions,expected,zscore,logp,fdr);
	if (bgFlag) {
		fprintf(fp, "\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",interactionsBg,expectedBg,zscoreBg,logpBg,
												totalPeak1Bg,totalPeak2Bg);
		fprintf(fp, "\t%lf\t%lf",logpDiff,zscore-zscoreBg);
	}
	fprintf(fp, "\t%d\n", thickness);
}



GenomeInteractionLibrary::GenomeInteractionLibrary() {
	interactions = NULL;
	numInteractions = 0;
	interactionList = NULL;
	lociFile = NULL;
	bedFile = NULL;
	fourCbedFile = NULL;
	removeOverlap = 0;
	bedFileCoverage = NULL;
	bedFileInterFrac = NULL;
	bedFileLocalFrac = NULL;
	interactionFile = NULL;
	threshold = 1e-5;
	zscoreThreshold = 1.0;
	minDist = 0;
	maxDist = -1;
	recordInteractions = 0;
	recordInteractionBg = 0;
	peakList = NULL;
}
GenomeInteractionLibrary::~GenomeInteractionLibrary() {
	if (interactionList != NULL) delete interactionList;
	if (interactions != NULL) {
		for (int i=0;i<numInteractions;i++) {
			if (interactions[i] != NULL) delete interactions[i];
		}
		delete []interactions;
	}
	interactionList = NULL;
	interactions = NULL;
	numInteractions = 0;
	if (lociFile != NULL) fclose(lociFile);
	if (bedFile != NULL) fclose(bedFile);
	if (fourCbedFile != NULL) fclose(fourCbedFile);
	if (bedFileCoverage != NULL) fclose(bedFileCoverage);
	if (bedFileLocalFrac != NULL) fclose(bedFileLocalFrac);
	if (bedFileInterFrac != NULL) fclose(bedFileInterFrac);
	if (interactionFile != NULL) fclose(interactionFile);
	if (peakList != NULL) delete peakList;
}
//benjamini-hochberg FDR
//Assumes interactions are sorted
void GenomeInteractionLibrary::calculateBenjaminiFDR() {
	double logTests = log(totalTests);
	// FDR = "p-value rank"/"total tests"*"p-value"
	double lastFDR = -1e10;
	int lastIndex = -1;
	for (int i=0;i<numInteractions;i++) {
		if (i<numInteractions-1) {
			if (interactions[i+1]->logp - interactions[i]->logp < 1e-100) continue;
		}
		double logFdr = interactions[i]->logp + logTests - log((double)i+1.0);
		if (logFdr > 0.0) logFdr = 0.0;
		if (logFdr < lastFDR) logFdr = lastFDR;

		for (int j=i;j>lastIndex;j--) {
			interactions[j]->fdr = exp(logFdr);
		}
		lastFDR = logFdr;
		lastIndex = i;
	}	
}
void GenomeInteractionLibrary::addInteraction(GenomeInteraction* gi) {
	if (interactionList == NULL) {
		interactionList = new LinkedList();
		for (int i=0;i<numInteractions;i++) {
			interactionList->add(interactions[i]);
		}
	}
	interactionList->add(gi);
	numInteractions++;
	if (interactionList->total > 2000000000) {
		fprintf(stderr, "!!! too many interactions - try raising the threshold!!!\n");
		exit(0);
	}
}
void GenomeInteractionLibrary::optimizeInteractionArray() {
	if (interactionList == NULL) return;
	interactions = (GenomeInteraction**)interactionList->toArray(numInteractions);
	delete interactionList;
	interactionList = NULL;
	sortInteractionArray();
	if (removeOverlap) removeOverlappingInteractions();
	calculateBenjaminiFDR();
}
void GenomeInteractionLibrary::removeOverlappingInteractions() {
	//fprintf(stderr, "removing overlapping interactions...\n");
	GenomeInteraction** bestInteractions = new GenomeInteraction*[numInteractions];
	for (int i=0;i<numInteractions;i++) bestInteractions[i] = interactions[i];
	qsort(interactions,numInteractions,sizeof(GenomeInteraction*),&cmpInteractionLocation);
	char* mask = new char[numInteractions];
	for (int i=0;i<numInteractions;i++) {
		mask[i]=0;
		interactions[i]->index=i;
	}

	for (int i=0;i<numInteractions;i++) {
		int index = bestInteractions[i]->index;
		double logp = bestInteractions[i]->logp;
		for (int j=index+1;j<numInteractions;j++) {
			if (strcmp(interactions[j]->chr1,bestInteractions[i]->chr1)!=0) break;
			if (interactions[j]->start1 > bestInteractions[i]->end1) break;
			if (interactions[j]->start2 <= bestInteractions[i]->end2 
						&& interactions[j]->end2 >= bestInteractions[i]->start2) {
				if (interactions[j]->logp >= logp) {
					interactions[j]->logp += 0.000001;
					mask[j] = 1;
				}
			}
		}
		for (int j=index-1;j>=0;j--) {
			if (strcmp(interactions[j]->chr1,bestInteractions[i]->chr1)!=0) break;
			if (interactions[j]->end1 < bestInteractions[i]->start1) break;
			if (interactions[j]->start2 <= bestInteractions[i]->end2 
						&& interactions[j]->end2 >= bestInteractions[i]->start2) {
				if (interactions[j]->logp >= logp) {
					interactions[j]->logp += 0.000001;
					mask[j] = 1;
				}
			}
		}
	}
	int newTotalInteractions = 0;
	for (int i=0;i<numInteractions;i++) {
		if (mask[i] == 0) {
			bestInteractions[newTotalInteractions++] = interactions[i];
		} else {
			delete interactions[i];
		}
	}
	delete []interactions;
	delete []mask;
	interactions = bestInteractions;
	numInteractions = newTotalInteractions;
	sortInteractionArray();
}

void GenomeInteractionLibrary::sortInteractionArray() {
	qsort(interactions,numInteractions,sizeof(GenomeInteraction*),&cmpInteractions);
}
void GenomeInteractionLibrary::sortInteractionDiff() {
	qsort(interactions,numInteractions,sizeof(GenomeInteraction*),&cmpInteractionDiff);
}
void GenomeInteractionLibrary::sortInteractionIndexes() {
	qsort(interactions,numInteractions,sizeof(GenomeInteraction*),&cmpInteractionIndexes);
}

int GenomeInteractionLibrary::read(char* file) {
	if (file == NULL) return -1;
	FILE* fp = fopen(file,"r");
	if (fp == NULL) {
		fprintf(stderr, "!!! Couldn't open interaction file: %s !!!\n", file);
		return -1;
	}
	char* buf = new char[BUFFER];
	char** cols = new char*[BUFFER];
	int numCols = 0;

	long long int resolutionEstimate = 0;
	int resEstN = 0;

	peakList = new PeakLibrary();

	while (fgets(buf, BUFFER, fp) != NULL) {
		split(buf, cols, numCols, '\t');
		if (numCols < 13) continue;
		if (checkInt(cols[3]) || checkInt(cols[4]) || checkInt(cols[9]) || checkInt(cols[10])) {
			continue;
		} else {
		}

		char* name1 = cols[1];
		char* chr1 = cols[2];
		int start1 = -1;
		sscanf(cols[3],"%d",&start1);
		int end1 = -1;
		sscanf(cols[4],"%d",&end1);
		char strand1 = STRAND_POSITIVE;
		if (cols[5][0] == '-' || cols[5][0] == '1') strand1 = STRAND_NEGATIVE;
		float v1 = 0.0;
		sscanf(cols[6],"%f",&v1);

		int info = 0;	
		Peak* p1 = peakList->checkForPeak(name1, chr1,start1,end1,strand1,info);
		if (p1 == NULL) {
			p1 = peakList->addPeak(name1,chr1,start1,end1,(start1+end1)/2,strand1,v1,0.0,NULL,0,0);
		}
		resolutionEstimate += end1-start1;
		resEstN ++;

		char* name2 = cols[7];
		char* chr2 = cols[8];
		int start2 = -1;
		sscanf(cols[9],"%d",&start2);
		int end2 = -1;
		sscanf(cols[10],"%d",&end2);
		char strand2 = STRAND_POSITIVE;
		if (cols[11][0] == '-' || cols[11][0] == '1') strand2 = STRAND_NEGATIVE;
		float v2 = 0.0;
		sscanf(cols[12],"%f",&v2);

		Peak* p2 = peakList->checkForPeak(name2, chr2,start2,end2,strand2,info);
		if (p2 == NULL) {
			p2 = peakList->addPeak(name2,chr2,start2,end2,(start2+end2)/2,strand2,v2,0.0,NULL,0,0);
		}
		resolutionEstimate += end2-start2;
		resEstN ++;

		GenomeInteraction* gi = new GenomeInteraction(p1, p2);
		gi->setName(cols[0]);
		addInteraction(gi);

	}
	
	delete []buf;
	delete []cols;
	fclose(fp);

	optimizeInteractionArray();

	peakList->sortChr();
	peakList->setDefaultPeakOrder();
	for (int i=0;i<numInteractions;i++) {
		interactions[i]->peakIndex1 = interactions[i]->peak1->index;
		interactions[i]->peakIndex2 = interactions[i]->peak2->index;
	}
	sortInteractionIndexes();

	if (resEstN > 0) resolutionEstimate /= resEstN;

	int resEstimate = ((int)((((double)resolutionEstimate)+50.0)/100.0))*100;
	fprintf(stderr, "\tRead %d total interactions from file: %s (resolution estimate: %d)\n", numInteractions, file,
							(int) resEstimate);
	return resEstimate;
}
void GenomeInteractionLibrary::setInteractionPeakSize(int res) {
	for (int i=0;i<numInteractions;i++) {
		int mid = (interactions[i]->start1+interactions[i]->end1)/2;
		interactions[i]->start1 = mid-res/2;
		interactions[i]->end1 = mid+res/2;
		mid = (interactions[i]->start2+interactions[i]->end2)/2;
		interactions[i]->start2 = mid-res/2;
		interactions[i]->end2 = mid+res/2;
	}
}
void GenomeInteractionLibrary::printCircos(FILE* fp, int maxInteractions) {
	for (int i=0;i<numInteractions;i++) {
		//if (i>=maxInteractions) break;
		interactions[i]->printCircos(fp,i+1);
	}
}
void GenomeInteractionLibrary::print(FILE* fp, int centerFlag) {
	if (fp == NULL) return;
	fprintf(fp,"InteractionID\tPeakID(1)\tchr(1)\tstart(1)\tend(1)\tstrand(1)\tTotal Reads(1)");
	fprintf(fp, "\tPeakID(2)\tchr(2)\tstart(2)\tend(2)\tstrand(2)\tTotal Reads(2)\tDistance\t");
	fprintf(fp, "Interaction Reads\tExpected Reads\tZ-score\tLogP\tFDR(Benjamini, based on %.2le total tests)",totalTests);
	if (recordInteractionBg) {
		fprintf(fp, "\tBg Interaction Reads\tBg Expected Reads\tBg Z-score\tBg LogP");
		fprintf(fp, "\tBg Total Reads Peak1\tBg Total Reads Peak2");
		fprintf(fp, "\tLogP vs. Bg\tZ-score Difference vs. Bg");
		fprintf(fp, "\tCircos Thickness\n");
	} else {
		fprintf(fp, "\tCircos Thickness\n");
	}

	for (int i=0;i<numInteractions;i++) {
		interactions[i]->print(fp,i+1,recordInteractionBg,centerFlag);
	}
}
int cmpInteractions(const void* a, const void* b) {
	GenomeInteraction* ga = *(GenomeInteraction**)a;
	GenomeInteraction* gb = *(GenomeInteraction**)b;
	if (ga->logp < gb->logp) return -1;
	if (ga->logp > gb->logp) return 1;
	if (ga->interactions < gb->interactions) return -1;
	if (ga->interactions > gb->interactions) return 1;
	return 0;
}
int cmpInteractionLocation(const void* a, const void* b) {
	GenomeInteraction* ga = *(GenomeInteraction**)a;
	GenomeInteraction* gb = *(GenomeInteraction**)b;
	int c = chrcmp((void*)&(ga->chr1),(void*)&(gb->chr1));
	if (c < 0) return -1;
	if (c > 0) return 1;
	if (ga->start1 < gb->start1) return -1;
	if (ga->start1 > gb->start1) return 1;
	return 0;
}
int cmpInteractionDiff(const void* a, const void* b) {
	GenomeInteraction* ga = *(GenomeInteraction**)a;
	GenomeInteraction* gb = *(GenomeInteraction**)b;
	if (ga->logpDiff < gb->logpDiff) return -1;
	if (ga->logpDiff > gb->logpDiff) return 1;
	if (ga->logp < gb->logp) return -1;
	if (ga->logp > gb->logp) return 1;
	if (ga->interactions < gb->interactions) return -1;
	if (ga->interactions > gb->interactions) return 1;
	return 0;
}
int cmpInteractionIndexes(const void* a, const void* b) {
	GenomeInteraction* ga = *(GenomeInteraction**)a;
	GenomeInteraction* gb = *(GenomeInteraction**)b;
	if (ga->peakIndex1 < gb->peakIndex1) return -1;
	if (ga->peakIndex1 > gb->peakIndex1) return 1;
	if (ga->peakIndex2 < gb->peakIndex2) return -1;
	if (ga->peakIndex2 > gb->peakIndex2) return 1;
	return 0;
}
int cmpInt(const void* a, const void* b) {
	int aa = *((int*)a);
	int bb = *((int*)b);
	if (aa < bb) return -1;
	if (aa > bb) return 1;
	return 0;
}


// class UniqMapChrs -------------------------------------------------------------------

unsigned char* UniqMapChrs::lookup = NULL;
unsigned char* UniqMapChrs::smask = NULL;
unsigned char* UniqMapChrs::emask = NULL;


void UniqMapChrs::getUniqMapStats(char* directory,long long int &gsize, long long int &mappability,
				LongInttable* chrSizes, LongInttable* chrUniqMap) {
	char* fname = new char[10000];
	int numCols=0;
	sprintf(fname,"%s/uniqMapStats.txt",directory);
	FILE* fp = fopen(fname,"r");
	if (fp == NULL) {
		fprintf(stderr, "Could not open mappability stat file %s\n", fname);
		delete []fname;
		return;
	}
	char* buffer = new char[BUFFER];
	char** cols = new char*[10000];
	long long int s=0;
	long long int m=0;
	while (fgets(buffer,BUFFER,fp)!=NULL) {
		split(buffer, cols, numCols,'\t');
		if (cols[0][0] == '#') continue;
		if (numCols > 2) { 
			if (strcmp(cols[0],"genome")==0) {
				sscanf(cols[1],"%lld",&mappability);
				sscanf(cols[2],"%lld",&gsize);
			}
			if (chrSizes != NULL) {
				sscanf(cols[2],"%lld",&s);
				chrSizes->insert(s,cols[0]);
			}
			if (chrUniqMap != NULL) {
				sscanf(cols[1],"%lld",&m);
				chrUniqMap->insert(m,cols[0]);
			}
		}
	}
	fclose(fp);
	delete []buffer;
	delete []cols;
	delete []fname;
}

UniqMapChrs::UniqMapChrs(char* chromosome, char* directory, int newFlag) {
	chr = new char[strlen(chromosome)+1];
	strcpy(chr, chromosome);
	size = 0;
	pstrand = NULL;
	nstrand = NULL;
	unitSize = 8;
	maxIndex = 0;
	numMappable = 0;

	pfile = new char[10000];	
	nfile = new char[10000];	
	sprintf(pfile,"%s/%s.p.uniqmap",directory,chr);
	sprintf(nfile,"%s/%s.n.uniqmap",directory,chr);
	
	if (lookup == NULL) {
		initializeLookup();
	}

	if (!newFlag) {
		loadData();
	}

}

UniqMapChrs::~UniqMapChrs() {
	if (pfile != NULL) delete []pfile;
	if (nfile != NULL) delete []nfile;
	if (chr != NULL) delete []chr;
	if (pstrand != NULL) delete []pstrand;
	if (nstrand != NULL) delete []nstrand;
}
void UniqMapChrs::initializeLookup() {
	//fprintf(stderr, "initializing lookup (%llx)\n", (unsigned long long int)lookup);
	lookup = new unsigned char[UCHAR_MAX];
	unsigned char x = 0;
	for (int c=0;c<=UCHAR_MAX;c++) {
		lookup[x] = countBits(x);
		//fprintf(stderr, "%d \t %d\n",x,lookup[x]);
		x++;
	}
	smask = new unsigned char[8];
	emask = new unsigned char[8];
	x = 255;
	for (int i=0;i<8;i++) {
		smask[i] = x << i;
		emask[7-i] = x >> i;
		//fprintf(stderr, "%d\t%d\t%d\n", i, smask[i],emask[7-i]);
	}
}
unsigned char UniqMapChrs::countBits(unsigned char x) {
	unsigned char initMask = 1;
	unsigned char mask = 1;
	unsigned char count = 0;
	for (int i=0;i<8;i++) {
		mask = initMask << i;
		if (mask & x) {
			count++;
		}
	}
	return count;
}

void UniqMapChrs::loadData() {

	FILE* fpp = fopen(pfile,"r");
	if (fpp == NULL) return;
	
	FILE* fpn = fopen(nfile,"r");
	if (fpn == NULL) {
		fclose(fpn);
		return;
	}
	//size_t ok = 0;

	//else read in data...
	unsigned int psize = 0;	
	(void)fread(&psize,sizeof(unsigned int),1,fpp);
	unsigned int nsize = 0;	
	(void)fread(&nsize,sizeof(unsigned int),1,fpn);
	if (nsize != psize) {
		fprintf(stderr, "!!! psize = %d nsize = %d\n",psize,nsize);
		fprintf(stderr, "Something is probably wrong with your uniqmap files for %s\n", chr);
		exit(0);
	}
	size = psize;
	maxIndex = size/unitSize+1;
	pstrand = new unsigned char[maxIndex+1];
	nstrand = new unsigned char[maxIndex+1];
	(void)fread(pstrand,sizeof(unsigned char),maxIndex,fpp);
	(void)fread(nstrand,sizeof(unsigned char),maxIndex,fpn);
	if (pstrand == NULL || nstrand == NULL) {
		fprintf(stderr, "Something is probably wrong with your uniqmap files for %s\n", chr);
		exit(0);
	}
	fclose(fpp);
	fclose(fpn);
	numMappable = 0;
	for (int i=0;i<maxIndex;i++) {
		numMappable += lookup[pstrand[i]];
		numMappable += lookup[nstrand[i]];
	}
	//fprintf(stderr, "\tLoaded uniqmap data for %s (%d total positions)\n", chr, size);
}
void UniqMapChrs::increaseMapSize(unsigned int newsize) {
	unsigned int newMaxIndex = newsize/unitSize+1;

	unsigned char* newp = new unsigned char[newMaxIndex];
	unsigned char* newn = new unsigned char[newMaxIndex];
	if (pstrand != NULL) {
		int stop=maxIndex;
		if (newMaxIndex < (unsigned int)maxIndex) {
			stop = newMaxIndex;
		}
		for (int i=0;i<stop;i++) {
			newp[i] = pstrand[i];
			newn[i] = nstrand[i];
		}
		delete []pstrand;
		delete []nstrand;
	}
	for (int i=maxIndex;i<(int)newMaxIndex;i++) {
		newp[i] = 0;
		newn[i] = 0;
	}
	pstrand = newp;
	nstrand = newn;
	maxIndex = newMaxIndex;
	size  = newsize;
}
void UniqMapChrs::setMappable(unsigned int pos, int dir) {
	unsigned int index = pos/unitSize;
	unsigned char bitIndex = 1 << (pos % unitSize);
	if (dir == 0) {
		if (!(pstrand[index] & bitIndex)) numMappable++;
		pstrand[index] = pstrand[index] | bitIndex;
	} else {
		if (!(nstrand[index] & bitIndex)) numMappable++;
		nstrand[index] = nstrand[index] | bitIndex;
	}
}
void UniqMapChrs::setUnmappable(unsigned int pos, int dir) {
	unsigned int index = pos/unitSize;
	unsigned char bitIndex = 1 << (pos % unitSize);
	unsigned char bitMask = ~bitIndex;
	if (dir == 0) {
		if (!(pstrand[index] & bitIndex)) numMappable++;
		pstrand[index] = pstrand[index] & bitMask;
	} else {
		if (!(nstrand[index] & bitIndex)) numMappable++;
		nstrand[index] = nstrand[index] & bitMask;
	}
}
void UniqMapChrs::optimizeSize() {
	int largestValue = 0;
	for (int i=maxIndex-1;i>=0;i--) {
		if (pstrand[i] != 0) {
			largestValue = i+1;
			break;
		}
		if (nstrand[i] != 0) {
			largestValue = i+1;
			break;
		}
	}
	increaseMapSize(largestValue*8);
}
void UniqMapChrs::saveFiles() {
	optimizeSize();

	fprintf(stderr, "\tSaving %s mapping data to %s/%s\n",chr,pfile,nfile);
	double ratio = ((double)numMappable)/((double)size*2.0);
	fprintf(stderr, "\t\t%d of %d total positions are mappable (%.2lf%%)\n",numMappable,size*2,ratio*100);

	FILE* fp = NULL;
	fp = fopen(pfile,"w");
	if (fp == NULL) {
		fprintf(stderr, "!!! Could not open %s for writing!!!\n", pfile);
		exit(0);
	}
	fwrite(&size, sizeof(unsigned int), 1, fp);
	fwrite(pstrand, sizeof(unsigned char), maxIndex, fp);
	fclose(fp);

	fp = fopen(nfile,"w");
	if (fp == NULL) {
		fprintf(stderr, "!!! Could not open %s for writing!!!\n", nfile);
		exit(0);
	}
	fwrite(&size, sizeof(unsigned int), 1, fp);
	fwrite(nstrand, sizeof(unsigned char), maxIndex, fp);
	fclose(fp);
}

int UniqMapChrs::countRegion(unsigned int start, unsigned int end, char strand) {

	int sIndex = start/unitSize;
	int sOffset = start % unitSize;
	int eIndex = end/unitSize;
	int eOffset = end % unitSize;

	if (end < start) return 0;
	if (sIndex >= maxIndex) return 0;
	if (sIndex < 0) {
		sIndex = 0;
		sOffset = 0;
	}
	if (eIndex >= maxIndex) {
		eIndex = maxIndex-1;
		eOffset = 7;
	}
	if (eIndex < 0) return 0;
	int count = 0;
	for (int i=sIndex+1;i<eIndex;i++) {
		if (strand == STRAND_POSITIVE || strand == STRAND_BOTH) count += lookup[pstrand[i]];
		if (strand == STRAND_NEGATIVE || strand == STRAND_BOTH) count += lookup[nstrand[i]];
	}
	if (strand == STRAND_POSITIVE || strand == STRAND_BOTH) {
		count += lookup[pstrand[sIndex] & smask[sOffset]];
		count += lookup[pstrand[eIndex] & emask[eOffset]];
		if (sIndex == eIndex) {
			count -= lookup[pstrand[sIndex]];
		}
	}
	if (strand == STRAND_NEGATIVE || strand == STRAND_BOTH) {
		count += lookup[nstrand[sIndex] & smask[sOffset]];
		count += lookup[nstrand[eIndex] & emask[eOffset]];
		if (sIndex == eIndex) {
			count -= lookup[nstrand[sIndex]];
		}
	}
	return count;
}


Tag* UniqMapChrs::getRegion(unsigned int start, unsigned int end, int & numTags) {

	int sIndex = start/unitSize;
	int sOffset = start % unitSize;
	int eIndex = end/unitSize;
	int eOffset = end % unitSize;
	numTags = 0;

	if (end < start) return NULL;
	if (sIndex >= maxIndex) return NULL;
	if (sIndex < 0) {
		sIndex = 0;
		sOffset = 0;
	}
	if (eIndex >= maxIndex) {
		eIndex = maxIndex-1;
		eOffset = 7;
	}
	if (eIndex < 0) return NULL;

	unsigned int length = end-start;
	Tag* tags = new Tag[length*2];
	numTags = 0;

	unsigned int p=start;
	for (int i=sIndex;i<=eIndex;i++) {
		int startOffset = 0;
		int endOffset = 7;
		if (i==sIndex) startOffset = sOffset;
		if (i==eIndex) endOffset = eOffset;
		for (int j=startOffset;j<=endOffset;j++) {
			unsigned char mask = 1 >> j;
			if (mask & pstrand[i]) {
				tags[numTags].p = p;
				tags[numTags].d = 0;
				tags[numTags].v = 1;
				tags[numTags].len = -1;
			}
			if (mask & nstrand[i]) {
				tags[numTags].p = p;
				tags[numTags].d = 1;
				tags[numTags].v = 1;
				tags[numTags].len = -1;
			}
			p++;
		}
	}
	return tags;
}

unsigned char* UniqMapChrs::getChr(char strand, int &length) {
	unsigned char* info = pstrand;
	if (strand == STRAND_NEGATIVE) info = nstrand;
	unsigned char* profile = new unsigned char[size];
	for (unsigned int i=0;i<size;i++) {
		int index = i/unitSize;
		if (info[index] & (1 << (i % unitSize))) {
		//if (pstrand[index] & (1 << (i % unitSize))) {
			profile[i]=1;
		} else {
			profile[i]=0;
		}
		length++;
	}
	return profile;
}





void UniqMapChrs::adjustProfile(float* profile, float *umcProfile, char strand, int fragLength, int maxPosition) {

	float minMapLevel = 1.0/((float)fragLength);
	float maplevel = 0;
	float historicAvg = 0;
	float historicN = 0;

	for (int i=1;i<maxPosition;i++) {
		int index = i;
		if (strand == STRAND_NEGATIVE) {
			index = maxPosition-i-1;
		}

		maplevel = umcProfile[index];
		if (maplevel < 0) {
			//maplevel = minMapLevel;
		}
		if (historicN > 0.5) {
			if (maplevel < 1.0-minMapLevel) {
				profile[index] = profile[index] + (historicAvg/historicN)*(1-maplevel);
			}
		}
		if (i<=fragLength) {
			historicN++;
			historicAvg += profile[index];
		} else {
			historicAvg += profile[index];
			if (strand == STRAND_NEGATIVE) {
				historicAvg -= profile[index+fragLength];
			} else {
				historicAvg -= profile[index-fragLength];
			}
		}
		if (historicAvg < 0) historicAvg = 0;
	}
}


float* UniqMapChrs::buildUniqMapProfile(char strand, int fraglen, int maxPosition) {

	float *levels = new float[maxPosition];
	float total = 0;
	int arrayLen = 0;
	for (int i=0;i<maxPosition;i++) levels[i] = 0;
	if (strand == STRAND_BOTH || strand == STRAND_POSITIVE) {
		total += fraglen;
		arrayLen = 0;
		unsigned char* array = getChr(STRAND_POSITIVE,arrayLen);
		int current = 0;
		for (int i=0;i<maxPosition-1;i++) {
			if (i>=arrayLen) break;
			current += array[i];
			if (i>fraglen) current -= array[i-fraglen];
			if (i+1 < maxPosition) levels[i+1] += (float)current;
			//levels[i+1] += (float)array[i];
		}
		delete []array;
	}
	if (strand == STRAND_BOTH || strand == STRAND_NEGATIVE) {
		total += fraglen;
		arrayLen = 0;
		unsigned char* array = getChr(STRAND_NEGATIVE,arrayLen);
		//unsigned char* array = getChr(STRAND_POSITIVE,arrayLen);
		float current = 0;
		for (int i=arrayLen-1;i>=0;i--) {
			current += array[i];
			if (i+fraglen<arrayLen) current -= array[i+fraglen];
			if (i>=maxPosition) continue;
			if (i+1 < maxPosition) levels[i+1] += (float)current;
		}
		if (array != NULL) delete []array;
	}
	if (maxPosition > arrayLen) {
		for (int i=arrayLen;i<maxPosition;i++) {
			levels[i] =0;
		}
	}
	
	for (int i=0;i<maxPosition;i++) {
		levels[i] /= total;
		//if (levels[i] > 0.6 && levels[i] < 0.7) fprintf(stderr, "%d\t%f\n",i,levels[i]);
	}
	return levels;
}




// class Tag ---------------------------------------------------------------

int Tag::precision = TAG_VALUE_RESOLUTION;
int PETag::precision = TAG_VALUE_RESOLUTION;

Tag::Tag() {
	p = 0;
	len = -1;
	v = 0.0;
	d = 0;
}

PETag::PETag() {
	init();
}
PETag::PETag(char* n, char* c, int np, char nd, float nv, int nlen) {
	init();
	name = n;
	chr1 = c;
	p1 = np;
	d1 = nd;
	v = nv;
	len1 = nlen;
}
void PETag::init() {
	name = NULL;
	chr1 = NULL;
	p1 = -1;
	d1 = 0;
	len1 = -1;
	chr2 = NULL;
	p2 = -1;
	d2 = 0;
	len2 = -1;
	v = 0.0;
}
PETag::~PETag() {
}
void PETag::print(FILE* fp) {
	print(fp, 0);
}
void PETag::copy(PETag* src) {
	chr1 = src->chr1;
	p1 = src->p1;
	d1 = src->d1;
	len1 = src->len1;
	chr2 = src->chr2;
	p2 = src->p2;
	d2 = src->d2;
	len2 = src->len2;
	v = src->v;
}

void PETag::print2ndHalf(FILE* fp) {
	if (d2 == STRAND_POSITIVE) {
		fprintf(fp, "%s\t%d\t%d\t4Ctag\t%.0f\t+\n",chr2,p2,p2+len2,v);
	} else {
		fprintf(fp, "%s\t%d\t%d\t4Ctag\t%.0f\t-\n",chr2,p2,p2-len2,v);
	}
}
void PETag::print(FILE* fp,int revFlag) {
	/*if (percision == 1) {
		if (v <  0.05) return;
	} else if (percision == 2) {
		if (v <  0.005) return;
	} else if (percision == 3) {
		if (v <  0.005) return;
	}*/
	if (name != NULL) fprintf(fp,"%s",name);
	if (revFlag) {
		fprintf(fp, "\t%s\t%d\t%d",chr2,p2,d2);
	} else {
		fprintf(fp, "\t%s\t%d\t%d",chr1,p1,d1);
	}
	if (precision == 1) {
		fprintf(fp, "\t%.1f",v);
	} else if (precision == 2) {
		fprintf(fp, "\t%.2f",v);
	} else {
		fprintf(fp, "\t%.3f",v);
	}
	if (revFlag) {
		fprintf(fp, "\t%d\t%s\t%d\t%d\t%d\n",len2,chr1,p1,d1,len1);
	} else {
		fprintf(fp, "\t%d\t%s\t%d\t%d\t%d\n",len1,chr2,p2,d2,len2);
	}
}

void Tag::copy(Tag* src) {
	p = src->p;
	len = src->len;
	v = src->v;
	d = src->d;
}
void Tag::print(FILE* fp,char* chr) {
	if (precision == 1) {
		fprintf(fp,"\t%s\t%d\t%d\t%.1f\t%d\n",chr,p,d,v,len);
	} else if (precision == 2) {
		//if (v <  0.005) return;
		fprintf(fp,"\t%s\t%d\t%d\t%.2f\t%d\n",chr,p,d,v,len);
	} else {
		//if (v <  0.05) return;
		fprintf(fp,"\t%s\t%d\t%d\t%.3f\t%d\n",chr,p,d,v,len);
	}
}


LinkedTag::LinkedTag(int pos, char dir, int length, float value,  LinkedTag* link) {
	p = pos;
	v = value;
	d = dir;
	len = length;
	tag = link;
}
LinkedTag::~LinkedTag() {
}
LinkedPETag::LinkedPETag(char* nc1, int np1, char nd1,int nlen1,char* nc2, int np2, char nd2, 
													int nlen2, float nv,  LinkedPETag* link) {
	chr1 = nc1;
	p1 = np1;
	d1 = nd1;
	len1 = nlen1;
	chr2 = nc2;
	p2 = np2;
	d2 = nd2;
	len2 = nlen2;
	v = nv;
	tag = link;
}
LinkedPETag::~LinkedPETag() {
}



void PeakLibrary::extractSequenceStats(char* genomeDir, FILE* fpout) {

	char* altName = new char[10000];
    int dirMode = 1;
    struct stat st_buf;
    stat(genomeDir,&st_buf);
    if (S_ISREG(st_buf.st_mode)) {
        dirMode = 0;
        fprintf(stderr, "\tChecking sequences from file: %s\n", genomeDir);
    } else {

		sprintf(altName,"%s/genome.fa",genomeDir);
    	stat(altName,&st_buf);
    	if (S_ISREG(st_buf.st_mode)) {
        	dirMode = 0;
        	fprintf(stderr, "\tChecking sequences from file: %s\n", altName);
			genomeDir = altName;
		} else {
			sprintf(altName,"%s/genome.fa.masked",genomeDir);
    		stat(altName,&st_buf);
    		if (S_ISREG(st_buf.st_mode)) {
        		dirMode = 0;
        		fprintf(stderr, "\tChecking sequences from file: %s\n", altName);
				genomeDir = altName;
			} else {
        		fprintf(stderr, "\tChecking sequences from directory: %s\n", genomeDir);
			}
		}
    }

	char** files = new char*[10000];
	int numFiles = 0;
	if (dirMode) {

		int dlen = strlen(genomeDir);
		DIR *dir = opendir(genomeDir);
		if (dir == NULL) {
			fprintf(stderr, "!!! Could not open directory %s !!!\n", genomeDir);
			exit(0);
		}
		struct dirent *dp=NULL;
		while ((dp=readdir(dir)) != NULL) {
			char* f = dp->d_name;
			int len = strlen(f);
			int good = 0;
			if (len > 3) {
				if (strcmp(&(f[len-3]),".fa")==0) {
					good = 1;
				}
			}
			if (len > 10) {
				if (strcmp(&(f[len-10]),".fa.masked")==0) {
					good = 1;
				}
			}
			if (good) {
				files[numFiles] = new char[len+1+dlen+1];
				sprintf(files[numFiles],"%s/%s",genomeDir,f);
				numFiles++;
			}
		}
		closedir(dir);
	} else {
		files[0] = new char[strlen(genomeDir)+1];
		strcpy(files[0],genomeDir);
		numFiles++;
	}

	long long int totalSeq =0;
	long long int totalN =0;
	long long int totalGC =0;
	char* buffer = new char[BUFFER];	
	char* currentChr = new char[BUFFER];	
	char* lastChr = new char[BUFFER];	

	fprintf(fpout,"ChrName\tLength\tN%%\tGC%%\n");

	for (int i=0;i<numFiles;i++) {

		lastChr[0] = '\0';

		FILE* fpin = fopen(files[i], "r");
		if (fpin == NULL) {
			fprintf(stderr, "\t!!Could not open sequence file: %s\n", genomeDir);
			exit(0);
		}
		long long int chrSeq = 0;
		long long int chrN = 0;
		long long int chrGC = 0;

		while (fgets(buffer,BUFFER,fpin)!=NULL) {
			buffer[strlen(buffer)-1]='\0';
			if (buffer[0] == '>') {
				int count = 1;
				while (buffer[count] != '\0') {
					if (buffer[count] < 33) {
						buffer[count] = '\0';
						break;
					}
					count++;
				}
				strcpy(currentChr,&(buffer[1]));
				if (strcmp(currentChr,lastChr)!=0 && lastChr[0] != '\0') {
					printSequenceStats(fpout,lastChr,chrSeq,chrN,chrGC);
					chrSeq = 0;
					chrN = 0;
					chrGC = 0;
				}
				strcpy(lastChr,currentChr);
				continue;
			} else {
				int c = 0;
				while (buffer[c] != '\0') {
					if (buffer[c] == 'N' || buffer[c] == 'n') {
						totalN++;
						totalSeq++;
						chrN++;
						chrSeq++;
					} else if (buffer[c] == 'c' || buffer[c] == 'C'
								|| buffer[c] == 'g' || buffer[c] == 'G') {
						totalGC++;
						totalSeq++;
						chrGC++;
						chrSeq++;
					} else {
						totalSeq++;
						chrSeq++;
					}
					c++;
				}
			}
		}
		if (lastChr[0] != '\0') {
			printSequenceStats(fpout,lastChr,chrSeq,chrN,chrGC);
		}
		delete [](files[i]);
	}
	printSequenceStats(fpout,(char*)"genome",totalSeq,totalN,totalGC);
	delete []buffer;
	delete []currentChr;
	delete []lastChr;
	delete []altName;
	delete []files;
}
void PeakLibrary::printSequenceStats(FILE* fp,char* chrName, long long int totalSeq,
						long long int totalN, long long int totalGC) {
	if (totalSeq < 1) {
		fprintf(fp, "%s\t%lld\t%.2lf%%\t%.2lf%%\n",chrName,totalSeq,0.0,0.0);
	} else {
		double ratioN = ((double)totalN)/((double)totalSeq);
		double ratioGC = ((double)totalGC)/((double)totalSeq);
		fprintf(fp, "%s\t%lld\t%.2lf%%\t%.2lf%%\n",chrName,totalSeq,
							ratioN*100.0,ratioGC*100.0);
	}
}


void PeakLibrary::extractSequence(char* genomeDir, FILE* fpout, int fastaFlag, int maskFlag) {


	char* altName = new char[10000];
    int dirMode = 1;
    struct stat st_buf;
    stat(genomeDir,&st_buf);
    if (S_ISREG(st_buf.st_mode)) {
        dirMode = 0;
        fprintf(stderr, "\tExtracting sequences from file: %s\n", genomeDir);
    } else {

		sprintf(altName,"%s/genome.fa",genomeDir);
    	stat(altName,&st_buf);
    	if (S_ISREG(st_buf.st_mode)) {
        	dirMode = 0;
        	fprintf(stderr, "\tExtracting sequences from file: %s\n", altName);
			genomeDir = altName;
		} else {
			sprintf(altName,"%s/genome.fa.masked",genomeDir);
    		stat(altName,&st_buf);
    		if (S_ISREG(st_buf.st_mode)) {
        		dirMode = 0;
        		fprintf(stderr, "\tExtracting sequences from file: %s\n", altName);
				genomeDir = altName;
			} else {
        		fprintf(stderr, "\tExtracting sequences from directory: %s\n", genomeDir);
			}
		}
    }
	char* seqBuffer1 = new char[MAX_SEQ_EXTRACT_BUFFER];
	char* seqBuffer2 = new char[MAX_SEQ_EXTRACT_BUFFER];

	if (dirMode) {

		char* chrFile = new char[100000];
		char slash[2] = "/";
		char fasta[4] = ".fa";
		char masked[8] = ".masked";

		char** keys = chrs->keys();
		qsort(keys,chrs->total,sizeof(char*),&chrcmp);
		for (int z=0;z<chrs->total;z++) {
			char* chr = keys[z];
			ChrPeaks* cp = (ChrPeaks*) chrs->search(keys[z]);

			int numChrPeaks = cp->numPeaks;
	
			if (numChrPeaks < 1) continue;
			strcpy(chrFile, genomeDir);
			strcat(chrFile, slash);
			strcat(chrFile, chr);
			strcat(chrFile, fasta);
	
			FILE* fpin = fopen(chrFile, "r");
			if (fpin == NULL) {
				strcat(chrFile, masked);
				fpin = fopen(chrFile,"r");
			}
			if (fpin == NULL) {
				fprintf(stderr, "\t!!Could not open file for %s (.fa or .fa.masked)\n", chr);
				continue;
			}
			cp->extractSequence(fpin,fpout,fastaFlag,maskFlag,seqBuffer1,seqBuffer2);
			fclose(fpin);

			delete [](keys[z]);
		}
		delete []chrFile;
		delete []keys;

	} else {

		char* buffer = new char[BUFFER];	
		char* currentChr = new char[BUFFER];	
		char* lastChr = new char[BUFFER];	
		lastChr[0] = '\0';

		fprintf(stderr, "\tLooking for peak sequences in a single file (%s)\n", genomeDir);
		FILE* fpin = fopen(genomeDir, "r");
		if (fpin == NULL) {
			fprintf(stderr, "\t!!Could not open sequence file: %s\n", genomeDir);
			exit(0);
		}


		off_t lastFilePosition = ftello(fpin);
		while (fgets(buffer,BUFFER,fpin)!=NULL) {
			if (buffer[0] == '>') {
				int count = 1;
				while (buffer[count] != '\0') {
					if (buffer[count] < 33) {
						buffer[count] = '\0';
						break;
					}
					count++;
				}
				strcpy(currentChr,&(buffer[1]));
				if (strcmp(currentChr,lastChr)!=0) {
	
					ChrPeaks* cp = (ChrPeaks*) chrs->search(currentChr);
					if (cp != NULL && cp->numPeaks > 0) {
						fseeko(fpin,lastFilePosition,SEEK_SET);
						cp->extractSequence(fpin,fpout,fastaFlag,maskFlag,seqBuffer1,seqBuffer2);
					}
					strcpy(lastChr, currentChr);
				}
				lastFilePosition = ftello(fpin);
			} else {
				lastFilePosition = ftello(fpin);
			}
		}
		fclose(fpin);
		delete []buffer;
		delete []currentChr;
		delete []lastChr;
	}
	delete []seqBuffer1;
	delete []seqBuffer2;
	delete []altName;
}

void ChrPeaks::extractSequence(FILE* fpin, FILE* fpout, int fastaFlag, int maskFlag,char* seqBuffer1, char* seqBuffer2) {

	if (peaks == NULL || numPeaks < 1) return;
	char* chr = peaks[0]->chr;
	int* finished = new int[numPeaks];
	for (int i=0;i<numPeaks;i++) finished[i]=0;

	int index = 0;
	int done = 0;
	int totalLen=0;
	int curStart = 0;
	//int curEnd = 0;
	int curSeqBuffer1Length = 0;

	seqBuffer1[0] = '\0';
	seqBuffer2[0] = '\0';

	int firstStart = peaks[0]->start-1;
	//int lastTail = peaks[0]->end;

	char* buffer = new char[BUFFER];	
	char* currentChr = new char[BUFFER];
	strcpy(currentChr, chr);

	off_t lastFilePosition = ftello(fpin);
	while (fgets(buffer,BUFFER,fpin)!=NULL) {

		int lineLen = strlen(buffer);
		if (lineLen > 0 && buffer[lineLen-1] == '\n') {
			buffer[lineLen-1] = '\0';
		}
		if (lineLen > 0 && buffer[lineLen-1] == '\r') {
			buffer[lineLen-1] = '\0';
		}
		if (buffer[0] == '>') {
			int count = 1;
			while (buffer[count] != '\0') {
				if (buffer[count] < 33) {
					buffer[count] = '\0';
					break;
				}
				count++;
			}
			strcpy(currentChr,&(buffer[1]));
			if (strcmp(currentChr,chr) != 0) {
				fseeko(fpin,lastFilePosition,SEEK_SET);
				break;
			}

			index = 0;
			done = 0;
			totalLen=0;
			seqBuffer1[0] = '\0';
			curSeqBuffer1Length = 0;
			curStart = 0;
			//curEnd = 0;

			firstStart = peaks[index]->start-1;
			//lastTail = peaks[index]->end;
			fprintf(stderr, "\tExtracting %d sequences from %s\n",numPeaks,chr);
			lastFilePosition = ftello(fpin);
			continue;
		}
		lastFilePosition = ftello(fpin);


		int L = strlen(buffer);
		totalLen += L;

		if (maskFlag) {
			for (int i=0;i<L;i++) {
				//mask lowercase characters
				if (buffer[i] > 95) buffer[i] = 'N';
			}
		}

		// index in this case is the last peak worth checking at the momment.
		for (int i=index+1;i<numPeaks;i++) {
			if ((peaks[i]->start-1) < totalLen) {
				index = i;
			} else {
				break;
			}
		}
	
		if (totalLen < firstStart) {
			seqBuffer1[0] = '\0';
			curSeqBuffer1Length = 0;
			//curEnd = 0;
			curStart = 0;
			continue;
		}

		if (seqBuffer1[0] == '\0') {
			curStart = totalLen-L;
		}
		if (totalLen - curStart > MAX_SEQ_EXTRACT_BUFFER) {
			fprintf(stderr, "ERROR: outstripped buffer!\n");
		}
		//can't use strcat - with large strings - way too slow
		for (int i=0;i<L;i++) {
			seqBuffer1[curSeqBuffer1Length+i] = buffer[i];
		}
		curSeqBuffer1Length += L;
		seqBuffer1[curSeqBuffer1Length]='\0';
		//curEnd = totalLen;

		for (int i= done; i<= index;i++) {
			if (finished[i]==1) continue;
			if (peaks[i]->end < totalLen) {
				int slen = peaks[i]->end - (peaks[i]->start-1);
				int sIndex = (peaks[i]->start-1) - curStart;
				seqBuffer2[0] = '\0';	
				int seqPos = 0;
				if (sIndex < 0) {
					int limit = sIndex*-1;
					for (;seqPos<limit;seqPos++){ 
						if (seqPos+1 >= MAX_SEQ_EXTRACT_BUFFER) {
							fprintf(stderr, "ERROR: outstripped buffer\n");
						}
						seqBuffer2[seqPos] = 'N';
					}
					seqBuffer2[seqPos] = '\0'; 
					sIndex = 0;
					slen += sIndex;
				}
				int lenCurSeq = curSeqBuffer1Length;
				int len2extract = slen;
				if (lenCurSeq < sIndex+slen) {
					len2extract = lenCurSeq - (sIndex+slen);
				}
				for (int j=sIndex;j<len2extract+sIndex;j++) {
					if (seqPos+1 >= MAX_SEQ_EXTRACT_BUFFER) {
						fprintf(stderr, "ERROR: outstripped buffer\n");
					}
					seqBuffer2[seqPos++] = seqBuffer1[j];
				}
				seqBuffer2[seqPos]='\0';

				if (lenCurSeq < sIndex+slen) {
					for (int j=0;j<sIndex+slen-lenCurSeq;j++) {
						seqBuffer2[seqPos++] = 'N';
					}
					if (seqPos+1 >= MAX_SEQ_EXTRACT_BUFFER) {
						fprintf(stderr, "ERROR: outstripped buffer\n");
					}
					seqBuffer2[seqPos]='\0';
				}
				cleanUpSeq(seqBuffer2);
				if (peaks[i]->strand == 1) {
					revopp(seqBuffer2);
				}
				finished[i]=1;
				if (fpout != NULL) {
					if (fastaFlag) {
						fprintf(fpout,">%s\n%s\n", peaks[i]->name, seqBuffer2);
					} else {
						fprintf(fpout,"%s\t%s\n", peaks[i]->name, seqBuffer2);
					}
				} else {
					peaks[i]->setSeq(seqBuffer2);
				}
			}
		}
		int chrDone = 0;
		for (int i=done;i<numPeaks;i++) {
			if (finished[i] == 0) {
				done = i;
				firstStart = (peaks[i]->start-1);
				break;
			}
			if (i == numPeaks-1) {
				chrDone = 1;
			}
		}
		if (chrDone == 1) {
			break;
		}
		if (curSeqBuffer1Length > MAX_SEQ_EXTRACT_BUFFER-2*BUFFER) {
			//fprintf(stderr, "\tAdjusting: Current index: %d\n", firstStart-curStart);
			int j=0;
			for (int i = firstStart-curStart;i<curSeqBuffer1Length;i++) {
				seqBuffer1[j++] = seqBuffer1[i];
			}
			seqBuffer1[j] = '\0';
			curSeqBuffer1Length = j;
			curStart = firstStart;
			//curEnd = firstStart+curSeqBuffer1Length;
		}
	}

	//need to finish sequences that couldn't be extracted...
	for (int i= done; i< numPeaks;i++) {
		if (finished[i]==1) continue;
		if (peaks[i]->start-1 >= totalLen) break;

		int slen = totalLen - (peaks[i]->start-1);
		int sIndex = (peaks[i]->start-1) - curStart;
		seqBuffer2[0] = '\0';	
		int seqPos = 0;
		if (sIndex < 0) {
			int limit = sIndex*-1;
			for (;seqPos<limit;seqPos++){ 
				seqBuffer2[seqPos] = 'N';
			}
			seqBuffer2[seqPos] = '\0'; 
			sIndex = 0;
			slen += sIndex;
		}
		int lenCurSeq = strlen(seqBuffer1);
		int len2extract = slen;
		if (lenCurSeq < sIndex+slen) {
			len2extract = lenCurSeq - (sIndex+slen);
		}
		for (int j=sIndex;j<len2extract+sIndex;j++) {
			seqBuffer2[seqPos++] = seqBuffer1[j];
		}
		seqBuffer2[seqPos]='\0';

		if (lenCurSeq < sIndex+slen) {
			for (int j=0;j<sIndex+slen-lenCurSeq;j++) {
				seqBuffer2[seqPos++] = 'N';
			}
			seqBuffer2[seqPos]='\0';
		}
		cleanUpSeq(seqBuffer2);
		if (peaks[i]->strand == 1) {
			revopp(seqBuffer2);
		}
		finished[i]=1;
		if (fpout != NULL) {
			if (fastaFlag) {
				fprintf(fpout,">%s\n%s\n", peaks[i]->name, seqBuffer2);
			} else {
				fprintf(fpout,"%s\t%s\n", peaks[i]->name, seqBuffer2);
			}
		} else {
			peaks[i]->setSeq(seqBuffer2);
		}
	}


	delete []currentChr;
	delete []finished;
	delete []buffer;
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
int checkStrand(char* str) {
	int value = -1;
	if (strcmp(str,"0")==0) {
		value = 0;
	} else if (strcmp(str,"1")==0) {
		value = 1;
	} else if (strcmp(str,"+")==0) {
		value = 0;
	} else if (strcmp(str,"-")==0) {
		value = 1;
	} else {
		value = -1;
	}
	return value;
}
void cleanUpSeq(char* seq) {
	if (seq == NULL) return;
	int len = strlen(seq);
	for (int i=0;i<len;i++) {
		if (seq[i] > 47 && seq[i] < 58) {
			//fprintf(stderr, "%c",seq[i]);
		} else if (seq[i] == 'a') {
			seq[i] = 'A';
		} else if (seq[i] == 'c') {
			seq[i] = 'C';
		} else if (seq[i] == 'g') {
			seq[i] = 'G';
		} else if (seq[i] == 't') {
			seq[i] = 'T';
		} else if (seq[i] == 'n') {
			seq[i] = 'N';
		} else if (seq[i] != 'A' && seq[i] != 'C' && seq[i] != 'G' && seq[i] != 'T') {
			seq[i] = 'N';
		}
	}
}
void reverse(char* str) {
	if (str == NULL) return;
	int len = strlen(str);
	int mid = len/2;
	for (int i=0;i<mid;i++) {
		char tmp = str[i];
		str[i] = str[len-1-i];
		str[len-1-i] = tmp;
	}
}
void revopp(char* seq) {
	if (seq == NULL) return;
	int len = strlen(seq);
	int mid = (int)((((double)len)/2.0)+0.5);


	if (seq[0] > 47 && seq[0] < 58) { // conservation
		for (int i=0;i<mid;i++) {
			int rIndex = len-1-i;
			char tmp = seq[rIndex];
			seq[rIndex] = seq[i];
			seq[i] = tmp;
		}
		return;
	}
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
	

// -------------- class NucleotideFreq ---------------------------------

NucleotideFreq::NucleotideFreq() {
	init();
}
NucleotideFreq::~NucleotideFreq() {
	if (freq != NULL) {
		for (int i=0;i<length;i++) {
			delete [](freq[i]);
		}
		delete []freq;
	}
	if (difreq != NULL) {
		for (int i=0;i<length;i++) {
			delete [](difreq[i]);
		}
		delete []difreq;
	}
	if (totals != NULL) delete []totals;
	if (ditotals != NULL) delete []ditotals;
	if (monoLookup != NULL) delete []monoLookup;
	if (diLookup != NULL) delete []diLookup;
	if (alpha != NULL) delete []alpha;
	if (gcNorm != NULL) delete []gcNorm;
	if (gcDist != NULL) delete []gcDist;
	if (gcOG != NULL) delete []gcOG;
}
void NucleotideFreq::init() {
	freq = NULL;
	difreq = NULL;
	totals = NULL;
	ditotals = NULL;
	length = -1;
	fragLength = -1;
	alphaSize = 4;
	diAlphaSize = 16;
	offset = 0;
	gcFile = NULL;

	gcInc = GCFREQ_DEFAULT_INC;
	gcDistTotal = 0;
	gcDist = NULL;
	gcNorm = NULL;
	gcOG = NULL;

	alpha = new char[4];
	alpha[0] = 'A';
	alpha[1] = 'C';
	alpha[2] = 'G';
	alpha[3] = 'T';
	monoLookup = new int[256];
	diLookup = new int*[256];
	for (int i=0;i<256;i++) {
		monoLookup[i] = -1;
		if (i=='A' || i=='a') monoLookup[i] = 0;
		else if (i=='C' || i=='c') monoLookup[i] = 1;
		else if (i=='G' || i=='g') monoLookup[i] = 2;
		else if (i=='T' || i=='t') monoLookup[i] = 3;
	}
	for (int i=0;i<256;i++) {
		diLookup[i] = new int[256];
		for (int j=0;j<256;j++) {
			diLookup[i][j] = -1;

			int iIndex = monoLookup[i];
			int jIndex = monoLookup[j];
			if (iIndex != -1 && jIndex != -1) {
				diLookup[i][j] = iIndex*alphaSize+jIndex;
			}
		}
	}
}
NucleotideFreq* NucleotideFreq::copy() {
	NucleotideFreq* nf = new NucleotideFreq();
	nf->length = length;
    nf->fragLength = fragLength;
    nf->alphaSize = alphaSize;
    nf->diAlphaSize = diAlphaSize;
    nf->offset = offset;
	nf->freq=NULL;
	nf->difreq=NULL;
	nf->totals = NULL;
	nf->ditotals=NULL;
	nf->alpha = NULL;
	nf->gcFile = NULL;

	nf->gcDistTotal = gcDistTotal;
	nf->gcInc = gcInc;
	nf->gcDist = new double[gcDistTotal];
	nf->gcNorm = new double[gcDistTotal];
	nf->gcOG = new double[gcDistTotal];
	for (int i=0;i<gcDistTotal;i++){ 
		if (gcDist != NULL) nf->gcDist[i]=gcDist[i];
		if (gcNorm != NULL) nf->gcNorm[i]=gcNorm[i];
		if (gcNorm != NULL) nf->gcOG[i]=gcOG[i];
	}

    nf->monoLookup=NULL;
    nf->diLookup=NULL;

	return nf;
}
NucleotideFreq::NucleotideFreq(int newOffset, int newLength,int gcfragLength) {
	init();
	length = newLength;
	fragLength = gcfragLength;
	offset = newOffset;
	initMatrix(newLength,gcfragLength);
}
void NucleotideFreq::initMatrix(int maxLength, int gcLength) {
	length = maxLength;
	if (length > BUFFER) {
		fprintf(stderr, "Exceeded maximum length of %d\n", BUFFER);
		exit(0);
	}
	freq = new double*[maxLength];
	difreq = new double*[maxLength];
	totals = new double[maxLength];
	ditotals = new double[maxLength];
	for (int i=0;i<maxLength;i++) {
		freq[i] = new double[alphaSize];
		difreq[i] = new double[diAlphaSize];
		for (int j=0;j<alphaSize;j++) freq[i][j] = 0.0;
		for (int j=0;j<diAlphaSize;j++) difreq[i][j] = 0.0;
		totals[i] = 0.0;
		ditotals[i] = 0.0;
	}

	if (gcLength <= 1) gcLength = maxLength;
	
	gcInc = 1.0/(double)gcLength;
	gcDistTotal = (int)((1.0/gcInc)+4.0);
	gcDist = new double[gcDistTotal];
	int Z = 0;
	double sub = 0.0;
	for (;Z<gcDistTotal;Z++) {
		gcDist[Z]=0.0;
		sub += gcInc;
		if (sub >= 1.0) break;
	}
	gcDist[++Z]=0.0;
	gcDistTotal = Z+1;
}
NucleotideFreq::NucleotideFreq(char* seqfile, int format, int newoffset, int maxL,char* gcFilename) {

	init();

	if (gcFilename != NULL) {
		gcFile = fopen(gcFilename, "w");
		if (gcFile == NULL) {
			fprintf(stderr, "!!! Couldn't open gc file for writing: %s !!!\n", gcFilename);
			exit(0);
		}
	}

	offset = newoffset;
	length = maxL;

	char* buf = new char[BUFFER];	
	char** cols = new char*[10000];
	int numCols = 0;

	char* name = NULL;
	char* seq = NULL;
	char* lastName = new char[BUFFER];
	char* nextName = new char[BUFFER];
	char* lastSeq = new char[BUFFER];
	lastName[0] = '\0';
	lastSeq[0] = '\0';
	nextName[0] = '\0';

	FILE* fp = fopen(seqfile, "r");
	if (fp == NULL) {
		fprintf(stderr, "!!!! Error: Could not open file %s !!!!\n", seqfile);
		exit(0);
	}
	int numSeq = 0;
	int addSeqFlag = 0;
	double w = 1.0;
	int lastNumSeq = 0;
	int seqmode = 0;

	int lastLinePlus = 0;

	SeqFreqStats* sfs = new SeqFreqStats();

	fprintf(stderr, "\tSequences processed:\n");
	while (fgets(buf, BUFFER, fp) != NULL) {
		if (numSeq % 100000 == 0 && numSeq != lastNumSeq) {
			fprintf(stderr, "\t\t%d\n", numSeq);
			lastNumSeq = numSeq;
		}
	
		if (format == SEQFILE_FORMAT_TSV) {
			split(buf, cols, numCols, '\t');
			if (numCols < 2) continue;
			name = cols[0];
			seq = cols[1];
			addSeqFlag = 1;
			w = 1.0;
		} else if (format == SEQFILE_FORMAT_FASTA) {
			int lineLength= strlen(buf);
			if (buf[lineLength-1] == '\n') buf[lineLength-1]='\0';
			if (buf[0] == '>') {
				strcpy(nextName,&(buf[1]));
				if (lastName[0] != '\0') {
					addSeqFlag=1;
					name = lastName;
					seq = lastSeq;
				} else {
					strcpy(lastName,nextName);
					lastSeq[0] = '\0';
				}
			} else {
				addSeqFlag = 0;
				strcat(lastSeq,buf);
			}
		} else if (format == SEQFILE_FORMAT_FASTQ) {
			int lineLength= strlen(buf);
			if (buf[lineLength-1] == '\n') buf[lineLength-1]='\0';
			if (buf[0] == '@' && lastLinePlus == 0) {
				seqmode = 1;
				strcpy(nextName,&(buf[1]));
				if (lastName[0] != '\0') {
					addSeqFlag=1;
					name = lastName;
					seq = lastSeq;
				} else {
					strcpy(lastName,nextName);
					lastSeq[0] = '\0';
				}
			} else if (buf[0] == '+') {
				seqmode = 0;
				lastLinePlus = 1;
			} else {
				addSeqFlag = 0;
				lastLinePlus = 0;
				if (seqmode) {
					strcat(lastSeq,buf);
				}
			}
		}					

		if (addSeqFlag) {
			addSequence(name,seq,w,-1,-1,sfs);
			numSeq++;
			if (format == SEQFILE_FORMAT_FASTA || format == SEQFILE_FORMAT_FASTQ) {
				strcpy(lastName,nextName);
				lastSeq[0] = '\0';
			}
		}
	}
	if (format == SEQFILE_FORMAT_FASTA || format == SEQFILE_FORMAT_FASTQ) {
		numSeq++;
		addSequence(name,seq,w,-1,-1,sfs);
	}
	fprintf(stderr, "\t\t%d total\n", numSeq);
	fclose(fp);

	delete []lastName;
	delete []nextName;
	delete []lastSeq;
	delete []buf;
	delete []cols;
}
SeqFreqStats::SeqFreqStats() {
	mono = new double[4];
	for (int i=0;i<4;i++) mono[i]=0.0;
	CpG=0.0;
	GC=0.0;
	AG=0.0;
	AC=0.0;
	N=0.0;
	NN=0.0;
}
SeqFreqStats::~SeqFreqStats() {
	delete mono;	
}

void NucleotideFreq::addChr(char* seq,int gcStart, int gcEnd) {
	if (gcStart < 0 || gcEnd < 0) {
		gcStart = 0;
		gcEnd = length-1;
	}
	int gcLength = gcEnd-gcStart+1;
	int seqLen = strlen(seq);
	int stop = seqLen-gcLength;
	int initStop = gcLength;
	if (seqLen < initStop) initStop = seqLen;

	int* nucTotals = new int[4];
	int* dinucTotals = new int[16];
	for (int i=0;i<4;i++) nucTotals[i]=0;
	double curGC = 0.0;
	int nuctotal = 0;
	double curTotal = 0.0;

	for (int i=0;i<initStop;i++) {
		int index = monoLookup[(int)seq[i]];
		if (index > -1) {
			nucTotals[index]++;
			nuctotal++;
			if (index == 1 || index == 2) {
			//if (index == 3) { // special
				curGC++;
			}
			curTotal++;
		}
	}
	for (int i=1;i<=stop;i++) {
		if (curTotal>0) {
			int index = (int)floor((curGC/curTotal)/gcInc+0.01);
			//if (index > gcDistTotal || index < 0) 
				//fprintf(stderr, "PROLEM %d\t%lf\t%lf\t%d\n", index,curGC,curTotal,gcLength);
			gcDist[index]++;
		}

		int lastIndex = monoLookup[(int)seq[i-1]];
		int nextIndex = monoLookup[(int)seq[i+gcLength-1]];
		if (lastIndex > -1) {
			if (lastIndex == 1 || lastIndex == 2) {
			//if (lastIndex == 3) {
				curGC--;
			}
			curTotal--;
		}
		if (nextIndex > -1) {
			nucTotals[nextIndex]++;
			nuctotal++;
			if (nextIndex == 1 || nextIndex == 2) {
			//if (nextIndex == 3) { // special
				curGC++;
			}
			curTotal++;
		}
	}
	
	for (int i=0;i<length;i++) {
		totals[i] += nuctotal;
		for (int j=0;j<4;j++) {
			freq[i][j] += ((double)nucTotals[j]);
		}
	}
	delete []nucTotals;
	delete []dinucTotals;
}

void NucleotideFreq::addSequence(char* name, char* seq,double w, int gcStart, int gcEnd,
																	SeqFreqStats* stats) {

	int seqLen = strlen(seq);
	if (freq == NULL) {
		if (length < 1) {
			length = seqLen;
			fprintf(stderr, "\t\tAuto detected maximum sequence length of %d bp\n", length);
		}
		initMatrix(length,gcEnd-gcStart+1);
	}

		
	int curLength = seqLen;	
	if (seqLen > length) curLength = length;

	if (gcStart < 0) gcStart = 0;
	if (gcEnd < 0) gcEnd = seqLen;

	for (int i=0;i<4;i++) stats->mono[i]=0.0;
	stats->N=0.0;
	stats->NN=0.0;
	stats->GC=0.0;
	stats->AG=0.0;
	stats->AC=0.0;
	stats->CpG=0.0;

	for (int i=0;i<seqLen;i++) {
		char bp = seq[i];
		int index = monoLookup[(int)bp];
		if (index > -1) {
			if (i < curLength) {
				freq[i][index]+=w;
				totals[i]+=w;
			}

			if (i >= gcStart && i <= gcEnd) {
				stats->N+=w;
				if (index > -1 && index < 4) stats->mono[index]+=w;
				if (index == 1 || index == 2) stats->GC+=w;
				//if (index == 2) stats->GC+=w; //Check for Andrea
				if (index == 0 || index == 2) stats->AG+=w;
				if (index == 0 || index == 1) stats->AC+=w;
			}

			if (i<seqLen-1) {
				int diindex = diLookup[(int)bp][(int)seq[i+1]];
				if (diindex > -1) {
					if (i < curLength) {
						difreq[i][diindex]+=w;
						ditotals[i]+=w;
					}
					if (i >= gcStart && i <= gcEnd) {
						stats->NN+=w;
						if (diindex == 6) {
							stats->CpG++;
						}
					}
				}
			}
		}
	}
	if (stats->N > 0) {
		if (stats->NN < 0.01) {
			stats->NN = 0.01;
		}
		for (int i=0;i<4;i++) stats->mono[i] /= stats->N;
		stats->CpG /= stats->NN;
		stats->GC /= stats->N;
		stats->AG /= stats->N;
		stats->AC /= stats->N;

//++++++++++++++++++++++++++++++++++++++++++++++
		//stats->GC = stats->mono[3];
//_______________________________________________

		int index = (int)floor(stats->GC/gcInc+0.01);
		gcDist[index]+=w;

		if (gcFile != NULL) {
			fprintf(gcFile, "%s\t%lf\t%lf\t%lf\t%lf\t%.0lf\n", name, stats->CpG, 
								stats->GC, stats->AG,stats->AC, stats->N);
		}
	}
}

void NucleotideFreq::print(FILE* fp) {

	if (fp == NULL) {
		fprintf(stderr, "!!! File pointer not valid for print nucleotide frequencies!!!\n");
		return;
	}
	fprintf(fp, "Offset\tA\tC\tG\tT\tA/T\tC/G\tA/G\tC/T\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\tN\tNN\n");	
	for (int i=0;i<length;i++) { 
		int p = i + offset;
		fprintf(fp,"%d",p);
		for (int j=0;j<4;j++) {
			if (totals[i] < 0.001) {
				fprintf(fp,"\t0.0");
			} else {
				fprintf(fp,"\t%lf",freq[i][j]/totals[i]);
			}
		}
		if (ditotals[i] < 0.001) {
			fprintf(fp,"\t0.0\t0.0\t0.0\t0.0");
		} else {
			double at = (difreq[i][0]+difreq[i][3]+difreq[i][12]+difreq[i][15])/ditotals[i];
			double cg = (difreq[i][5]+difreq[i][6]+difreq[i][9]+difreq[i][10])/ditotals[i];
			double ag = (difreq[i][0]+difreq[i][2]+difreq[i][8]+difreq[i][10])/ditotals[i];
			double ct = (difreq[i][5]+difreq[i][7]+difreq[i][13]+difreq[i][15])/ditotals[i];
			fprintf(fp,"\t%lf\t%lf\t%lf\t%lf", at, cg, ag,ct);
		}
		for (int j=0;j<16;j++) {
			if (ditotals[i] < 0.001) {
				fprintf(fp,"\t0.0");
			} else {
				fprintf(fp,"\t%lf",difreq[i][j]/ditotals[i]);
			}
		}
		fprintf(fp,"\t%.1lf\t%.1lf\n", totals[i], ditotals[i]);
	}

}
double NucleotideFreq::printGC(FILE* fp) {
	if (fp == NULL) return 0.0;
	fprintf(fp, "GC%%\tTotal\tNormalized Fraction(PDF)\n");

	double total = 0.0;
	double sum = 0.0;
	for (int i=0;i<gcDistTotal;i++) {
		total += gcDist[i];
		double value = gcInc*(double)i;
		sum += value*gcDist[i];
	}
	for (int i=0;i<gcDistTotal;i++) {
		double value = gcInc*(double)i;
		fprintf(fp, "%.3lf\t%.1lf\t%lf", value, gcDist[i], gcDist[i]/total/gcInc);
		if (gcNorm != NULL) {
			fprintf(fp, "\t%lf", gcNorm[i]);
		}
		fprintf(fp,"\n");
	}
	double avgGC = 0;
	if (total > 0) avgGC = sum/total;
	return avgGC;
}
void NucleotideFreq::readGCcontentFile(char* filename) {
	if (gcDist != NULL) delete []gcDist;
	gcDistTotal = 0;
	gcInc = 0;
	FILE* fp = fopen(filename, "r");
	if (fp == NULL) {
		fprintf(stderr, "!!!Could not open GC content file: %s!!!\n",filename);
		return;
	}
	char* buffer = new char[BUFFER];
	char** line = new char*[1000];
	int numCols = 0;
	gcDist = new double[10000];

	double freq = 0.0;
	double fraction = 0.0;
	int lineNum = 0;
	while (fgets(buffer, BUFFER, fp) != NULL) {
		lineNum++;
		if (lineNum < 2) continue;
		split(buffer,line, numCols, '\t');
		sscanf(line[2],"%lf",&fraction);
		sscanf(line[0],"%lf",&freq);
		gcDist[gcDistTotal++] = fraction;
	}
	fclose(fp);
	if (gcDistTotal > 1) {
		gcInc = 1/((double)(gcDistTotal-1));
	}
	double tmpTotal = 0.0;
	for (int i=0;i<gcDistTotal;i++) {
		tmpTotal += gcDist[i];
	}
	for (int i=0;i<gcDistTotal;i++) {
		gcDist[i] /= tmpTotal*gcInc;
	}
	delete []line;
	delete []buffer;
}

void NucleotideFreq::printGCnormFile(char* filename) {

	if (filename == NULL) return;

	FILE* fp = fopen(filename,"w");

	fprintf(fp, "GC%%\tCount\tNormFactor\tAdjusted Counts\n");
	for (int i=0;i<gcDistTotal;i++) {
		double freq = ((double)i)*gcInc;
		fprintf(fp, "%lf\t%lf\t%lf\t%lf\n", freq, gcDist[i],gcNorm[i],gcNorm[i]*gcDist[i]);
	}
	fclose(fp);

}
double NucleotideFreq::calculateGCNormalization(NucleotideFreq *ctrl, double minNormRatio, 
											double maxNormRatio, double window, char* fname) {

	double total = 0.0;
	double ctrlTotal = 0.0;
	for (int i=0;i<gcDistTotal;i++) total += gcDist[i];
	for (int i=0;i<ctrl->gcDistTotal;i++) ctrlTotal += ctrl->gcDist[i];

	FILE *outfp = NULL;
	if (fname != NULL) {
		outfp = fopen(fname,"w");
	}	

	double rReads = 0.0;
	if (outfp != NULL) fprintf(outfp, "GC%%\tCount\tNormFactor\tAdjusted Counts\n");
	if (gcNorm != NULL) {
		delete []gcNorm;
	}
	gcNorm = new double[gcDistTotal];
	for (int i=0;i<gcDistTotal;i++) {
		double sum = 0.0;
		double sumCtrl=0.0;
		double N = 0.0;
		double Nctrl=0.0;
		double freq = ((double)i)*gcInc;

		for (int j=0;j<gcDistTotal;j++) {
			double f = ((double)j)*gcInc;
			double diff = f-freq;
			if (diff < -1*window) continue;
			if (diff > window) break;
			if (j!=i) continue;
			sum += gcDist[j]/gcInc;
			N += 1.0;
		}
		for (int j=0;j<ctrl->gcDistTotal;j++) {
			double f = ((double)j)*ctrl->gcInc;
			double diff = f-freq;
			if (diff < -1*window) continue;
			if (diff > window) break;
			sumCtrl += ctrl->gcDist[j]/ctrl->gcInc;
			Nctrl += 1.0;
		}
		if (Nctrl < 1.0) Nctrl = 1.0;
		if (N < 1.0) N = 1.0;
		sum /= N*total;
		sumCtrl /= Nctrl*ctrlTotal;
		double ratio = 1.0;
		if (sum > 0) {
			ratio = sumCtrl/sum;
		}

		double maxPossible = gcOG[i]*maxNormRatio;
		double minPossible = gcOG[i]*minNormRatio;
		double newtotal = gcDist[i]*ratio;

		if (newtotal < minPossible) {
			ratio = minPossible/gcDist[i];
			rReads += fabs(minPossible-newtotal);
		}
		if (newtotal > maxPossible) {
			ratio = maxPossible/gcDist[i];
			rReads += fabs(maxPossible-newtotal);
		}

		/*
		if (ratio > maxNormRatio) {
			rReads += fabs(gcDist[i]*ratio-gcDist[i]*maxNormRatio);
			ratio = maxNormRatio;
		}
		if (ratio < minNormRatio) {
			rReads += fabs(gcDist[i]*ratio-gcDist[i]*maxNormRatio);
			ratio = minNormRatio;
		}
		*/
		gcNorm[i] = ratio;
		if (outfp != NULL) {
			fprintf(outfp, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", freq, gcDist[i],gcNorm[i],gcNorm[i]*gcDist[i],sum,sumCtrl);
		}
	}
	if (outfp != NULL) {
		fclose(outfp);
	}
	if (total > 0) {
		rReads /= total;
	}
	return rReads;
}

int determineSeqFileFormat(char* filename) {
	FILE* fp = fopen(filename, "r");
	if (fp == NULL) {
		fprintf(stderr, "!!! Could not open sequence file %s !!!\n", filename);
		return SEQFILE_FORMAT_UNKNOWN;
	}
	int format = SEQFILE_FORMAT_UNKNOWN;

	char* buf = new char[BUFFER];
	while (fgets(buf, BUFFER, fp) != NULL) {
		if (buf[0] == '@') {
			format = SEQFILE_FORMAT_FASTQ;
			break;
		} else if (buf[0] == '>') {
			format = SEQFILE_FORMAT_FASTA;
			break;
		} else {
			int i=0;
			while (buf[i] != '\0') {
				if (buf[i] == '\t') {
					format = SEQFILE_FORMAT_TSV;
					break;
				}
				i++;
			}
			break;
		}
		break;
	}
	fclose(fp);
	delete []buf;
	return format;
}

void parseUCSCpositionStr(char* str, char* &chr, int &start, int &end) {
	chr = str;
	char* p1 = new char[100];
	char* p2 = new char[100];
	int indexp1 = -1;
	int indexp2 = -1;
	int i = 0;
	while (str[i] != '\0') {
		if (str[i] == ':') {
			str[i] = '\0';
			i++;
			indexp1 = 0;
		} else if (indexp1 > -1) {
			if (str[i]=='-') {
				p1[indexp1] = '\0';
				indexp1 = -1;
				indexp2 = 0;
			} else if (str[i] == ',') {
			} else {
				p1[indexp1++] = str[i];
			}
			i++;
		} else if (indexp2 > -1) {
			if (str[i]=='-') {
				p2[indexp2] = '\0';
				indexp2 = -1;
			} else if (str[i] == ',') {
			} else {
				p2[indexp2++] = str[i];
			}
			i++;
		} else {
			i++;
		}
	}
	p2[indexp2] = '\0';
	sscanf(p1,"%d",&start);
	sscanf(p2,"%d",&end);
	delete []p1;
	delete []p2;
}

HiCparams::HiCparams() {
	init();
}
HiCparams::~HiCparams() {
}
void HiCparams::init() {
	minDist = -1;
	maxDist = -1;
	res = -1;
	superRes = -1;
	logFlag = 1;
	boundaryScale = 1000000;
	minExpReads = 0.0;
	relativeFlag = 0;
	fragLengthEstimate = 1;

}
void HiCparams::setRes(int r, int sr) {
	res = r;
	superRes = sr;
	if (superRes < res) superRes = res;
}

Genome3D::Genome3D() {
	matrix = NULL;
	numMatrix = 0;
	names = NULL;
	p = NULL;
	fixedStep = GENOME3D_DEFAULT_FIXED_STEP;
}
Genome3D::~Genome3D() {
	if (p != NULL) {
		for (int i=0;i<numMatrix;i++) {
			delete [](p[i]);
		}
		delete []p;
	}
	if (color != NULL) delete []color;
}
void Genome3D::init(double** m, int N, char** nnames, double nfixedStep) {
	matrix = m;
	numMatrix = N;
	names = nnames;
	if (nfixedStep > 0.0) fixedStep = nfixedStep;
	fixedStepMax = fixedStep*4.0;
	fixedStepMin = fixedStep/4.0;

	p = new double*[numMatrix];
	for (int i=0;i<numMatrix;i++) {
		p[i] = new double[3];
	}
	color = new int[numMatrix];
	for (int i=0;i<numMatrix;i++) {
		color[i] = 153;
	}
	color[0] = 27;
	color[numMatrix-1] = 35;

	inRange = new int[numMatrix];
	numInRange = 0;
	normalizeMatrix(0);
	initStructure(GENOME3D_INIT_POINTS_LINE);
}
void Genome3D::initStructure(int method) {

	double diagDistance = sqrt(fixedStep/3.0);

	p[0][0] = 0.0;	
	p[0][1] = 0.0;	
	p[0][2] = 0.0;	
	for (int i=1;i<numMatrix;i++) {
		if (GENOME3D_INIT_POINTS_LINE == method) {
			double r = ((double)rand())/((double)RAND_MAX)*0.5+1.0;
			p[i][0] = p[i-1][0] + diagDistance*r;
			r = ((double)rand())/((double)RAND_MAX)*0.5+1.0;
			p[i][1] = p[i-1][1] + diagDistance*r;
			r = ((double)rand())/((double)RAND_MAX)*0.5+1.0;
			p[i][2] = p[i-1][2] + diagDistance*r;
		}
	}
}
void Genome3D::addColorsFromClusters(char* clusterFile) {
	FILE* fp = fopen(clusterFile, "r");
	if (fp == NULL) {
		fprintf(stderr, "!!! Can't open cluster file %s\n", clusterFile);
		exit(0);
	}
	Inttable* ids = new Inttable();
	for (int i=0;i<numMatrix;i++) {
		ids->insert(i,names[i]);
	}
	Inttable* clusterColors = new Inttable();
	int* availableColors = new int[10];
	availableColors[0] = 153; //black
	availableColors[1] = 35; //red
	availableColors[2] = 28; //blue
	availableColors[3] = 51; //green
	availableColors[4] = 31; //voilet
	availableColors[5] = 77; //orange
	availableColors[6] = 71; //light blue
	availableColors[7] = 75; //brown
	availableColors[8] = 11; //aqua
	availableColors[9] = 73; //dark purple

	char* buf = new char[BUFFER];
    char** line = new char*[10000];
    int numCols = 0;

	int curColor = 0;
	while (fgets(buf, BUFFER, fp) != NULL) {
        split(buf, line, numCols, '\t');
		if (numCols < 2) continue;
		if (strcmp(line[1],"-1")==0) continue;
		int index = ids->search(line[0]);
		if (index == EMPTY_INT) continue;
		int cindex = clusterColors->search(line[1]);
		if (cindex == EMPTY_INT) {
			cindex = availableColors[curColor++];
			if (curColor > 9) curColor = 0;
			clusterColors->insert(cindex,line[1]);
		}
		color[index] = cindex;
	}
	fclose(fp);



	delete []buf;
	delete []line;
	
}
void Genome3D::print(FILE* fp) {
	for (int i=0;i<numMatrix;i++) {
		fprintf(fp, "%s\t%lf\t%lf\t%lf\t%d\n", names[i],p[i][0],p[i][1],p[i][2],color[i]);
	}
}

void Genome3D::optimize(int maxIterations) {

	for (int i=0;i<maxIterations;i++) {
		if (i % 10000 == 0) fprintf(stderr, "\t%d\n", i);
		int randIndex = (int) (((double)numMatrix)*((double)rand())/((double)RAND_MAX));
		calcMoveOne(randIndex,0);
		if (0) {
		if (i % 1 == 0) fprintf(stderr, "\t%d\n", i);
		for (int j=0;j<numMatrix;j++) {
			calcMoveOne(j,0);
		}
		for (int j=numMatrix-1;j>=0;j--) {
			calcMoveOne(j,0);
		}
		}
	}
}
void Genome3D::calcMoveOne(int index, int method) {

	double v[3];
	double curV[3];
	int numInRange = 0;
	if (index>0) inRange[numInRange++]=index-1;
	if (index<numMatrix-1) inRange[numInRange++]=index+1;


	for (int i=0;i<numMatrix;i++) {
		if (i==index) continue;
		double dist = 0.0;
		for (int j=0;j<3;j++) {
			curV[j] = p[i][j]-p[index][j];
			dist += curV[j]*curV[j];
		}
		dist = sqrt(dist);
		if (dist < fixedStepMax+fixedStepMin) inRange[numInRange++]=i;

		double curMagnitude = (matrix[index][i]/maxDist);
		fprintf(stderr, "%d %d %lf, %lf\n", index, i, dist,curMagnitude);
		for (int j=0;j<3;j++) {
			v[j] += (curV[j]/dist)*curMagnitude;
		}
		
	}

	double magnitude = 0.0;
	for (int j=0;j<3;j++) {
		magnitude += v[j]*v[j];
	}
	magnitude = sqrt(magnitude);
	if (magnitude < 1e-20) magnitude = 1e-20;
	for (int j=0;j<3;j++) {
		v[j] /= magnitude;
	}
	double scaleFactor = magnitude/(sqrt(3.0)*numMatrix)*fixedStep;
	//fprintf(stderr, "Force: %d %lf\t%lf\t%lf\tmag=%lf(%lf,%lf,%lf)\n", index,p[index][0],p[index][1],p[index][2],scaleFactor,v[0],v[1],v[2]);

	if (index > 0) {
		scaleFactor = checkConnected(p[index],v,scaleFactor,p[index-1]);
		//fprintf(stderr, "Force: %d %lf\t%lf\t%lf\tmag=%lf(%lf,%lf,%lf)\n", index,p[index][0],p[index][1],p[index][2],scaleFactor,v[0],v[1],v[2]);
	}
	if (index < numMatrix-1) {
		scaleFactor = checkConnected(p[index],v,scaleFactor,p[index+1]);
		//fprintf(stderr, "scalefactor = %lf (more)\n", scaleFactor);
	}
	//scaleFactor = checkConnected(p[index],v,scaleFactor,p1,p2);

	for (int i=0;i<numInRange;i++) {
		scaleFactor = checkCollision(p[index],v,scaleFactor,p[inRange[i]]);
		//fprintf(stderr, "scalefactor = %lf (%d)\n", scaleFactor, inRange[i]);
	}
	//fprintf(stderr, "Adjut: %d %lf\t%lf\t%lf\tmag=%lf(%lf,%lf,%lf)\n", index,p[index][0],p[index][1],p[index][2],scaleFactor,v[0],v[1],v[2]);
	if (scaleFactor > 0.0) {
		for (int i=0;i<3;i++) {
			p[index][i] += v[i]*scaleFactor;
		}
	}
}
double Genome3D::checkConnected(double* p0, double* v, double magnitude, double* p1) {
	double target[3];
	double v2[3];
	double dist = 0;
	for (int i=0;i<3;i++) {
		target[i] = p0[i]+magnitude*v[i];
		dist += (target[i]-p1[i])*(target[i]-p1[i]);
		v2[i] = target[i]-p1[i];
	}
	dist = sqrt(dist);
	if (dist <= fixedStepMax && dist >= fixedStepMin) return magnitude;
	//fprintf(stderr, "Distance from=%lf\n", dist);

	double dist2 = 0;
	for (int i=0;i<3;i++) {
		if (dist > fixedStepMax) {
			v[i] = p1[i]+v2[i]/dist*fixedStepMax - p0[i];
		} else if (dist < fixedStepMin) {
			v[i] = p1[i]+v2[i]/dist*fixedStepMin - p0[i];
		}
		dist2 += v[i]*v[i];
	}
	dist2 = sqrt(dist2);
	for (int i=0;i<3;i++) {
		v[i] /= dist2;
	}
	return dist2;
}

double Genome3D::checkCollision(double* p0, double* v, double magnitude, double* p1) {
	double a = 0;
	double b = 0;
	double c = 0;
	for (int i=0;i<3;i++) {
		c += (p0[i]-p1[i])*(p0[i]-p1[i]);
		b += 2*v[i]*(p0[i]-p1[i]);
		a += v[i]*v[i];
	}

	double sterm = b*b - 4*a*c;
	if (sterm < 0.0) {
		// this means that the vector does not even touch the locus
		return magnitude;
	}
	sterm = sqrt(sterm);
	double s1 = (-b+sterm)/(2*a);
	double s2 = (-b-sterm)/(2*a);
	if (s1 > 0 && s1 < magnitude) magnitude=s1;
	if (s2 > 0 && s2 < magnitude) magnitude=s2;
	return magnitude;
}

void Genome3D::normalizeMatrix(int method) {
	minDist = FLT_MAX;
	maxDist = FLT_MIN;
	for (int i=0;i<numMatrix;i++) {
		for (int j=i+1;j<numMatrix;j++) {
			if (minDist > matrix[i][j]) minDist = matrix[i][j];
			if (maxDist < fabs(matrix[i][j])) maxDist = matrix[i][j];
		}
	}
}
