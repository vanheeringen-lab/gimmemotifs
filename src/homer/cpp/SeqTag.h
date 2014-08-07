#include <stdio.h>
#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <cmath>
#include <time.h>
#include <limits.h>
#include <pthread.h>

#include "Hashtable.h"
#include "statistics.h"
#include "Motif2.h"
#include "Clustering.h"

#ifndef SEQTAG_H
#define SEQTAG_H

#define _FILE_OFFSET_BITS 64

#define CORRECTION_FACTOR 1

#define WHITE_SPACE 7

#define TAGDIR_MIN_TAG_VALUE 0.000001

#define PEAK_INC 25000
#define PEAK_TAG_INC 500
#define MAX_PEAK_LIBRARIES 100
#define MAX_READ_LENGTH 32000
#define MAX_TAGS_PER_BP 32000
#define PEAK_LIBRARY_DEFAULT_SIZE 1000000

#define PEAK_LIBRARY_NOVEL_PEAK 0
#define PEAK_LIBRARY_SAME_NAME 1
#define PEAK_LIBRARY_SAME_POSITION 2

#define PEAK_READ_MODE_NORMAL 0
#define PEAK_READ_MODE_ANNOTATION 1
#define PEAK_READ_MODE_COUNTING 2

#define TAGDIR_CHRFILES 1
#define TAGDIR_ONEFILE 2

#define TAGLIBRARY_LOWMEM 0 
#define TAGLIBRARY_HIGHMEM 1

#define FORMAT_UNKNOWN 0
#define FORMAT_BOWTIE 1
#define FORMAT_BOWTIE_COLOR 2
#define FORMAT_BED 3 
#define FORMAT_ELANDRESULT 4
#define FORMAT_ELANDEXPORT 5
#define FORMAT_ELANDEXTENDED 6
#define FORMAT_SAM 7 
#define FORMAT_PEAK 8
#define FORMAT_MCPGBED 9
#define FORMAT_HICSUMMARY 10
#define FORMAT_LISTER_ALLC 11
#define FORMAT_BISMARK 12

#define CYTOSINE_CONTEXT_CG 0
#define CYTOSINE_CONTEXT_CHG 1
#define CYTOSINE_CONTEXT_CHH 2
#define CYTOSINE_CONTEXT_ALL 3

#define PEAKS_FIND_UNMETHYLC 1
#define PEAKS_FIND_METHYLC 2

#define AVERAGE_METHYLC -1234.0
#define MINIMUM_C_PER_PEAK_METHYLC 6

#define MODE_UNIQUE 0
#define MODE_KEEPONE 1
#define MODE_KEEPALL 2
#define MODE_BED_FORCE2ONE 3
#define MODE_BED_FORCE5TH 4

#define NULL_REF -123456789
#define NULL_OFFSET -123456789
#define ALL_PEAK_EXPS -1
#define POSITIVE_STRAND 0 
#define NEGATIVE_STRAND 1 
#define BOTH_STRANDS 2
#define NULL_INT_VALUE -123456789

#define AUTOCORRELATION_HALFWINDOW 7
#define AUTOCORRELATION_OFFSETMIN 40
#define AUTOCORRELATION_BACKOFFSET -2
#define AUTOCORRELATION_ENRICHMENT 150
#define AUTOCORRELATION_RNA_THRESHOLD 8.0
#define AUTOCORRELATION_CHIP_THRESHOLD 1.5 

#define TAGADJUST_DEFAULT 75
#define TAGADJUST_AUTO -123456789
#define FRAGMENT_LEN_AUTO -123456789
#define FRAGMENT_LEN_GIVEN -12345

#define COVERAGE_MIN_NORMLENGTH 25

#define MAXIMUM_CHR_LENGTH 1000000000

#define CIGAR_ERROR -12345678

#define COUNT_MODE_TOTAL 0
#define COUNT_MODE_TBP 1
#define COUNT_MODE_RATIO 2

#define FINDPEAKS_MINDIST_DEFAULT 2.0

#define DEFAULT_TSS_DISPERSION_SIZE 5

#define OUTPUT_MODE_COUNT 1
#define OUTPUT_MODE_PEAKTAGS 2
#define OUTPUT_MODE_TAGS 3

#define TAGLIBRARY_SINGLE_FLAG 0
#define TAGLIBRARY_PE_READ1_FLAG 1
#define TAGLIBRARY_PE_READ2_FLAG 2

#define FLOAT_ZERO 1e-30

#define STRAND_POSITIVE 0
#define STRAND_NEGATIVE 1 
#define STRAND_BOTH 2
#define STRAND_SEPARATE 3

#define DEFAULT_NORM_TOTAL 1e7

#define PEAKFINDER_TBP_AUTO -123456789
#define PEAKFINDER_FDRSIZE 1000
#define PEAKFINDER_MODE_POISSON 2
#define PEAKFINDER_FILTER_MODE_FDR 1
#define PEAKFINDER_FILTER_MODE_POISSON 2
#define PEAKFINDER_FILTER_MODE_THRESH 3

#define DIFFPEAK_MODE_DIFF 1
#define DIFFPEAK_MODE_SAME 2
#define DIFFPEAK_MODE_REV 3
//#define DIFFPEAK_MODE_FOLDCHANGE 4
//#define DIFFPEAK_MODE_POISSON 5

#define TAG_VALUE_RESOLUTION 1

#define MAX_ANNOTATIONS_PER_CHR 1e6
#define PEAKRATIO_TOLERANCE_BP 5

#define PEAKFRACTION_REGION_MINDIST 4.0
#define REGION_MODE_HISTONE 0
#define REGION_MODE_GROSEQ 1
#define REGION_GROSEQ_FOLDCHANGE 3.0
#define GROSEQ_TSSFOLDCHANGE 4.0
#define GROSEQ_BODYFOLDCHANGE 4.0
#define GROSEQ_CONFPVALUE 0.00001
#define GROSEQ_ENDFOLDCHANGE 10.0
#define GROSEQ_TSSSIZE 250
#define GROSEQ_MINBODYSIZE 1000
#define GROSEQ_MAXBODYSIZE 10000

#define GROSEQ_DENOVO_PSEUDOCOUNT 2.0

#define GROSEQ_METHOD_FOLD 0
#define GROSEQ_METHOD_LEVEL 1

#define UCSC_CHIPSEQ 0
#define UCSC_RNASEQ 1
#define UCSC_METHYLATED_CpG 2
#define UCSC_UNMETHYLATED_CpG 3
#define UCSC_DNASE 4
#define UCSC_TSS 5

#define UCSC_RES_MAX 0
#define UCSC_RES_AVG 1

#define GENOME3D_DEFAULT_FIXED_STEP 0.01
#define GENOME3D_INIT_POINTS_LINE 0
#define GENOME3D_CHECK_COLLISION 0
#define GENOME3D_CHECK_CONNECTED 1

#define PEAK_STYLE_CHIPSEQ 0
#define PEAK_STYLE_HISTONE 1
#define PEAK_STYLE_GROSEQ 2
#define PEAK_STYLE_TSS 3
#define PEAK_STYLE_DNASE 4
#define PEAK_STYLE_METHYLC 5
#define PEAK_STYLE_SUPERENHANCERS 6

#define DEFAULT_SUPERENHANCER_SLOPE_WINDOW 10
#define DEFAULT_SUPERENHANCER_SLOPE_THRESHOLD 1.0


#define BYTES_PER_LINE 6.1

#define UCSC_BEDGRAPH 0
#define UCSC_BIGWIG 1

#define DEFAULT_GSIZE 2000000000

#define MAX_PEAKS_AT_A_TIME 2000000

#define SEQFILE_FORMAT_UNKNOWN 0
#define SEQFILE_FORMAT_FASTQ 1
#define SEQFILE_FORMAT_FASTA 2
#define SEQFILE_FORMAT_TSV 3
#define SEQFILE_FORMAT_ILLUMINAPE 4

#define HIC_MASK_NORM_SEQDEPTH 0x1
#define HIC_MASK_NORM_DISTANCE 0x2
#define HIC_MASK_CORRELATION 0x4
#define HIC_MASK_CLUSTER 0x8
#define HIC_MASK_CREATEMODEL 0x10
#define HIC_MASK_EXPECTED 0x20
#define HIC_MASK_LOGPVALUES 0x40
#define HIC_MASK_NO_MATRIX 0x80
#define HIC_MASK_LOCISCORE 0x100
#define HIC_MASK_INTERACTIONS 0x200
#define HIC_MASK_CLUSTERFIXED 0x400
#define HIC_MASK_LOGPVALUESEXACT 0x800
#define HIC_MASK_INTERACTIONBEDFILE 0x1000
#define HIC_MASK_NORM_ZSCORE 0x2000
#define HIC_MASK_BOUNDARIES 0x4000
#define HIC_MASK_INTERACTION4CBEDFILE 0x8000
#define HIC_MASK_INTERACTIONSTATSFILES 0x10000
#define HIC_MASK_INTERACTIONTAGFILE 0x20000
#define HIC_MASK_RAWANDEXPECTED 0x40000

#define HIC_MATRIX_LIMIT 30000
#define HIC_RATIO_PSEUDO_COUNT 0.5

#define GCFREQ_DEFAULT_INC 0.01

#define HOMERCONFIG_SECTION_NONE 0
#define HOMERCONFIG_SECTION_SOFTWARE 1
#define HOMERCONFIG_SECTION_PROMOTERS 2
#define HOMERCONFIG_SECTION_GENOMES 3

#define RESTRICTION_SITE_DISTRIBUTION_SCALE 10
#define RESTRICTION_SITE_MIN_DISTANCE 10

#define OLIGO_NORMALIZATION_BUFFER 20
#define MAX_SEQ_EXTRACT_BUFFER 50000000

#define INTERCHROMOSOMAL_APPROX_DISTANCE 20000000

#define ZIPPED_FLAG_GZ 1
#define ZIPPED_FLAG_BZ2 2
#define ZIPPED_FLAG_ZIP 3
#define ZIPPED_FLAG_BAM 4

#define HiC_NORM_TYPE_DISTANCE 0 
#define HiC_NORM_TYPE_NODISTANCE 1

class TagLibrary;
class ChrTags;
class LinkedTag;
class LinkedPETag;
class Tag;
class PETag;
class Peak;
class ChrPeaks;
class PeakLibrary;
class Genome;
class MergePeak;
class MergePeaksArray;
class PeakMapping;
class UniqMapChrs;
class NucleotideFreq;
class SeqFreqStats;
class HiCBgModel;
class GenomeInteraction;
class GenomeInteractionLibrary;
class PeakSNPs;
class SNP;
class Genome3D;
class HiCBgModelChr;
class HiCparams;


class HomerConfigInfo {
public:
	char* name;
	char* directory;
	char* url;
	char* description;
	char* version;
	char** parameters;
	int numParameters;
	HomerConfigInfo(char*,char*,char*,char*,char*,char*);
	~HomerConfigInfo();
};
class HomerConfig {
public:
	Hashtable* software;	
	Hashtable* promoters;
	Hashtable* genomes;
	HomerConfig();
	~HomerConfig();
	void readConfigFile();
	static const char* homeDirectory;
	char* getGenomeDirectory(char* gname);
	void listAvailableGenomes();
};

class Genome {
public:
	char* name;
	char* directory;
	long long int totalbp;
	long long int totalMappablebp;
	
	Genome();
	~Genome();
	void setName(char*);
	void setDirectory(char*);
};

class PeakFinder {
public:
	char* name;
	char* directory;
	char* uniqMapDirectory;
	char* outputFileName;
	Genome* genome;
	int peakSize;
	int localSize;
	int inputSize;
	float tagThresh;
	float minTagThresh;
	int minDist;
	double tagsUsedForClustering;
	double totalTags;
	double totalInputTags;
	double normTotal;
	float maxtbp;
	float maxtbpInput;
	float mintbp;
	float mintbpInput;
	float tbpThreshold;
	float tbp;
	float tpp;
	int stitchMode;
	float tbpInput;
	int tagAdjust;
	int tagAdjustInput;
	char strand;
	long long int gsize;
	long long int realGsize;
	double fdr;
	float fdrThresh;
	double poisson;
	double poissonInput;
	double poissonLocal;
	float poissonThresh;
	float poissonThreshRelaxed;
	float pvalue;
	int method;
	int groseqMethod;
	int filterMode;
	int diffMode;
	int tbpAuto;
	int numPeaks;
	int numSuper;
	int superWindow;
	double superSlope;
	char* typicalFile;
	TagLibrary* tags;
	TagLibrary* input;
	double * fdrTable;
	int fdrSize;
	double * poissonTable;
	double tagsInPeaks;
	int centerFlag;
	int nfrFlag;
	int nfrSize;
	int regionFlag;
	double hysteresis;
	char* extraheader;
	double foldTranscriptStart;
	double foldTranscriptBody;
	double pseudoTags;
	int tssSize;
	int minBodySize;
	int maxBodySize;
	double endFold;
	double transcriptConfidenceThreshold;

	int style;

	int mCflag;
	double expectedMethylC;
	double methylCthreshold;
	int minNumC;
	
	double inputFold;
	double localFold;
	double clonalFold;
	int numFilterPeaks;
	int inputPeaks;
	int clonalPeaks;
	int localPeaks;
	char* cmd;
		
	PeakFinder();
	~PeakFinder();
	PeakLibrary* findPeaks();
	PeakLibrary* findmCPeaks();
	PeakLibrary* findGroSeqRegions();
	PeakLibrary* findGroSeqTranscripts();
	void setGenome(Genome*);
	void addHeader(char*);
	void setPeakSize(int);
	void setName(char* name);
	void setCMD(char* name);
	void setDirectory(char* name);
	void setOutputFile(char* name);
	void setTagLibraries(TagLibrary* exp, TagLibrary* input);
	void setMaxTBP(double,double);
	void setMinTBP(double,double);
	void setGenomeSize(long long int);
	void determineMaxTBP();
	void determineThreshold();
	void approxFdrTable();
	void checkParameters();
	void print(FILE*);
	void printGroSeq(FILE*);
	PeakLibrary* filterPeaks(PeakLibrary*);
	double getGroSeqThreshold(double fold, double confidence);

};

class TagLibrary {
public:

	Hashtable* chrs;
	Hashtable* uniqmapchrs;
	double totalTags;
	long long int totalPositions;
	char* genome;
	char* name;
	char* cmd;
	char* directory;
	char* uniqmapDirectory;
	int memType;
	int medianTagsPerPosition;
	int fragmentLengthEstimate;
	int fragmentLengthSetting;
	long long int gsizeEstimate;
	int tagAdjust;
	int peakSizeEstimate;
	int medianTagCount;
	int gcNormalized;
	double gcAverage;
	int oligoNormalized;
	double averageTagsPerPosition;
	double averageTagLength;
	double tbp;
	double parseAlignmentCpGMinValue;
	int parseAlignmentCcontext;
	float maxtbp;
	float mintbp;
	int maxmismatches;
	double minmapq;
	int revStrand;
	int singleFile;
	int minReadLength;
	int maxReadLength;

	int numCPUs;
	pthread_mutex_t mutex;
	pthread_mutex_t mutexMatrix;
	pthread_mutex_t mutexWashuFile;
	int mutexIndex;
	int mutexTotal;
	int mutexPeakIndex;
	char** mutexChrKeys;

	char* singleFilename;
	FILE* tagFile;
	FILE* tagFileR1;
	FILE* tagFileR2;
	char* restrictionSite;
	FILE* washuMatrixOutput;
	int peStyleFlag;
	int sspeFlag;
	int flipFlag;
	int minDistHiC;
	int maxDistHiC; //if same as minDistHiC, then it's disabled

	int pairedEndFlag;
	int mCflag;
	int forceSingleReadFlag;
	Hashtable* chrNames;
	
	TagLibrary(char* directory);
	~TagLibrary();

	void readTagDirectory();
	void readAndSave();
	void parseAlignmentFiles(char** files, int numFiles, int format, int mode,
								char** tagDirs, int numDirs,char** tagFiles, int numTagFiles);
	void setMaxTBP(float maxTagsPerBp);
	void setMinTBP(float minTagsPerBp);
	void setSingleFile(int flag);
	void setTagAdjust(int centerDistance);
	void setFragLength(int fragLength);
	void setSingleRead(int singleReadFlag);
	void addAlignedPETag(PETag* petag);
	PeakLibrary* findGroSeqRegions(int regionMode, TagLibrary* input, char* uniqmapDirectory,char strand,
						int tssSize, int minBodySize, int maxBodySize,  double threshold, double foldTranscriptStart,
						double foldTranscriptBody, double endFold, double pseudoTags, double inputFold);
	PeakLibrary* findGroSeqTranscripts(int regionMode, TagLibrary* input, char* uniqmapDirectory,char strand,
						int tssSize, int minBodySize, int maxBodySize,  double threshold, double foldTranscriptStart,
						double foldTranscriptBody, double endFold, double pseudoTags, double inputFold,int groseqMethod);

	double* getTagCountDistribution(FILE* nfp, int &max);
	double* getTagLengthDistribution(FILE* nfp, int &max);
	void autoCorrelateTags(FILE* outputfile, int windowSize, double maxTags);
	void annotateTagLocations(PeakLibrary* anns, FILE* statsOutput, FILE* annOutput);
	double getAdjustedTagTotal();
	void setRestrictionSite(char* site);
	void assignPETagsToRestrictionSites(char* site, int maxMisMatches,char* genomeDirectory,int mode,
							int midFlag,int removeSelfLigationFlag, int removeRestrictionEnds, int fragLength);


	void printTagInfo();
	void printTagInfo(FILE*);
	void setGenome(char*);
	void setCMD(char*);
	void setName(char*);
	void addChrToName(char* str);
	//internal functions
	void readAlignment(char* files, int format, int mode, int PEflag);
	void readPEAlignment(char* files, int format, int mode);
	void addTagDirectory(char* tagDir);
	void readSingleTagFile();
	void printSingleTagFile();
	void readNewTagFile(char* filename);
	void readTagInfo();
	void makeDirectory();
	void trimIlluminaReadName(char* name);
	int parseMDZstr(char* misStr);
	int getRightCoordFromCIGAR(char* str, int dir, char* cigarCodes, int* cigarLens, int &numCodes, int &initLen);
	double estimateContamination(TagLibrary* input, int peaksize, float minThreshold);
	void decontaminate(TagLibrary* input, double fraction, int maxDistance);
	NucleotideFreq*  checkTagSeqBias(char* genomeDirectory,int freqStart, int freqEnd,
									OligoArray* &oligos,int oligoStart,int oligoEnd);
	void addTag(char* chr,int pos,char dir,int length, float value);
	void addAlignedTag(char* name, char* chr,int pos,char dir,int length, float value,int PEflag);
	char* getTagFileName(char* chr);
	char* getDirFileName(char* filename);
	void optimizeTagFiles();
	void cleanUpChrName(char*);
	void makeName();
	Doubletable* getChrSizesFromTags();
	int getQualityTotal(char* quality);
	PeakLibrary* findPutativePeaks(int peakSize, int minDist, char strand, float minCount);
	PeakLibrary* findmCPeaks(int peakSize, char strand, int mCflag, double mCthresh, int minNumC, TagLibrary* input);
	void printBedGraph(FILE* fp, double normTotalTags, char strand, int resolution,int negFlag,
				double fileSize, char* color, char* bedname, int method, int lastTagFlag,
				char* uniqMapDirectory,int style,int condenseFlag, int circosFlag, Peak* circosPeak,
				TagLibrary* input, double pseudoCounts, int logFlag, double normLength);
	Doubletable* getPeakTagCounts(PeakLibrary* p, char strand); // more efficient
	void normalizeTagCountsGC(char* gcCtrlFile, NucleotideFreq* nf,double minFold,double maxFold, double gcWindow,
																	double maxPerror);
	void normalizeTagCountsFixedOligo(char* genomeDirectory,OligoArray* oligos,int oligoStart, int oligoEnd,
																	double minFold, double maxFold);
	void normalizeTagCountsOligos(char* genomeDirectory,Hashtable* oligos, int oligoLength, 
					int regionStart, int regionEnd,double minFold, double maxFold,float maxPerBp);
	double* getPETagDistribution(int windowSize, int largeWindowSize, int largeResolution,
					char* outputPrefix, int &arrayLength);
	void addUniqMap(char* directory);
	void removePETagBackground(int fragLength);
	void removeTagSpikes(int size, double fold);
	PeakLibrary* getCoveragePeaks(char* chr, int start, int end, int resolution,int superResolution);
	PeakLibrary* getCoverageRestrictionFragments(char* chr, int start, int end, int superResolution,
												char* siteSequence, int maxMisMatches,char* genomeDirectory);
	void makeHiCBgModel(HiCBgModel* model,PeakLibrary* &modelPeaks, PeakLibrary* modelPeaks2,
					int peaks1IsGenomeFlag, int fullModelFlag, HiCparams* params);
	double **calcCorrelationMatrix(double** matrix, int numPeaks,double** expectedMatrix,double minReads,int maxIndexDiff);
	static double calcCorrelation(double* a1, double* a2, int n, double* e1, double* e2, double min, int maxIndexDiff,
																int indexI, int indexJ);
	void scoreInteractionBoundaries(PeakLibrary* peaks1,
					int resolution, HiCBgModel* model, int actionFlag);
	double** makeHiCMatrix(PeakLibrary* peaks1, PeakLibrary* peaks2, 
					int resolution, HiCBgModel* model, int actionFlag,
					GenomeInteractionLibrary* gil, HiCparams* params);
	void makeHiCHistogram(char* outputfilename, PeakLibrary* peaks1, int sizebp, 
					int resolution, HiCBgModel* model, int actionFlag, HiCparams* params);
	void getPETagTotals(PeakLibrary* peaks);
	void getPETagTotalsThread(int cpu, PeakLibrary* peaks);
	void calcHiCMatrix(PeakLibrary* peaks1,PeakLibrary* peaks2,double** matrix,
                    int resolution, HiCBgModel* model, int revFlag, int numHistBins,
                    int actionFlag, double totalInteractions, GenomeInteractionLibrary* gil,
					HiCparams* params);
	void makeHiCMatrixThread(int cpu, PeakLibrary* peaks1,PeakLibrary* peaks2,double** matrix,
                    int resolution, HiCBgModel* model, int revFlag, int numHistBins,
                    int actionFlag, double totalInteractions, GenomeInteractionLibrary* gil,
					HiCparams* params);
	

};
class ThreadArgs_adjustTotals {
public:
	int cpu;
	ChrPeaks* cpeaks;
	double* errors;
	double* tmpTotals;
	double expectInter;
	int minIndex;
	int maxIndex;
	double totalModelReads;
	int resolution;
	HiCBgModelChr* chrModel;
	PeakLibrary* refPeaks;
	ThreadArgs_adjustTotals(int cpu, ChrPeaks* peaks, double* errors, double* tmpTotals,
						double expInter,int minIndex,int maxIndex,double totalModelReads,
						int resolution,HiCBgModelChr* chrModel,PeakLibrary* refPeaks);
};
void* ChrPeaks_adjustTotalsThread(void* threadArgs);

class ThreadArgs_getPETagTotals {
public:
    TagLibrary* tags;
    int cpu;
    PeakLibrary* peaks;
    ThreadArgs_getPETagTotals(int cpu, TagLibrary* tags, PeakLibrary* peaks);
};
void* TagLibrary_getPETagTotalsThread(void* threadArgs);

class ThreadArgs_makeHiCMatrix {
public:
    TagLibrary* tags;
    int cpu;
    PeakLibrary* peaks1;
    PeakLibrary* peaks2;
	double** matrix;
	int resolution;
	HiCBgModel* model;
	int revFlag;
	int numHistBins;
	int actionFlag;
	double totalInteractions;
	HiCparams* params;
	GenomeInteractionLibrary* gil;
    ThreadArgs_makeHiCMatrix(int cpu, TagLibrary* tags, PeakLibrary* peaks1, PeakLibrary* peaks2,
				double** matrix, int resolution, HiCBgModel* model, int revFlag,
				int numHistBins, int actionFlag , double totalInteractions,
				GenomeInteractionLibrary* gil, HiCparams* params);
};
void* TagLibrary_makeHiCMatrixThread(void* threadArgs);







class GenomeInteractionLibrary {
public:
	GenomeInteraction** interactions;
	int numInteractions;

	int resolution;
	int removeOverlap;
	FILE* interactionFile;
	FILE* bedFile;
	FILE* fourCbedFile;
	FILE* bedFileCoverage;
	FILE* bedFileInterFrac;
	FILE* bedFileLocalFrac;
	FILE* lociFile;
	double threshold;
	double zscoreThreshold;
	int minDist;
	int maxDist;
	double totalTests;
	int recordInteractions;
	int recordInteractionBg;
	pthread_mutex_t mutex;

	LinkedList* interactionList;
	PeakLibrary* peakList;
	
	GenomeInteractionLibrary();
	~GenomeInteractionLibrary();
	void setInteractionPeakSize(int res);
	void addInteraction(GenomeInteraction*);
	void optimizeInteractionArray();
	void sortInteractionArray();
	void sortInteractionDiff();
	void adjustBgStats();
	void sortInteractionIndexes();
	void calculateBenjaminiFDR();
//---
	void removeOverlappingInteractions();
	void printCircos(FILE* fp,int maxInteractions);
	void print(FILE* fp, int centerFlag);
	int read(char* file); // returns resolution estimate

};

class GenomeInteraction {
public:
	char* name;
	char* chr1;
	int start1;
	int end1;
	char* chr2;
	int start2;
	int end2;
	int index;
	int thickness;
	Peak* peak1;
	Peak* peak2;

	double totalPeak1;	
	double totalPeak2;	
	int peakIndex1;	
	int peakIndex2;	
	double logp;
	double zscore;
	double fdr;
	double interactions;
	double expected;

	double totalPeak1Bg;	
	double totalPeak2Bg;	
	double logpBg;
	double zscoreBg;
	double fdrBg;
	double interactionsBg;
	double expectedBg;
	double logpDiff;
	double fdrDiff;

	GenomeInteraction();	
	GenomeInteraction(Peak* p1, Peak* p2);	
	GenomeInteraction(char* name, char* chr1, int start1, int end1, double total1, int index1,Peak* p1,
							char* chr2, int start2, int end2, double total2, int index2,Peak* p2,
							double interactions,double expected, double logp, double zscore, int thickness);	
	~GenomeInteraction();	
	void setName(char*);
	void printCircos(FILE* fp, int count);
	void printPETags(FILE* fp, PETag** iTags, int numPETags,int format);
	void print(FILE* fp, int count, int bgFlag, int centerFlag);
	void init();
};

class HiCBgModelChr {
public: 
	int res;
	int totalRegions;
	double* expected;
	double* std;
	double* scaleFactor;
	double* stdScaleFactor;
	double* expectedN;
	double interChr;
	double interChrN;
	double interStd;
	double interStdScale;
	double interChrScaleFactor;
	double wtotal;
	int bgFlag;

	int maxIndex;
	int maxUsedIndex;
	int stdFlag;

	HiCBgModelChr(int resolution);
	~HiCBgModelChr();
	void init();
	void initialize();
	void normalize();
	void normalizeStd();
	void smooth(int bins);
};

class HiCBgModel : public HiCBgModelChr {
public: 

	int ready;
	char* filename;
	char* directory;
	Doubletable* chrExpect;
	Doubletable* chrSizes;
	Hashtable* chrs;
	int fullModelFlag;

	double avgCoverage;
	double stdCoverage;
	double stdFilter;
	double minFilter;
	double maxCoverage;
	double minCoverage;
	double totalModelReads;
	double modelError;
	int goodRegions;
	int badRegions;
	int customFlag;
	PeakLibrary* refPeaks;
	pthread_mutex_t mutex;

	double **expectMatrix;
	int numPeaks1;
	int numPeaks2;


	HiCBgModel(int res, char* directory, int customFlag);
	~HiCBgModel();
	void init();
	void initialize();
	char* getDefaultFileName(char* directory);
	int load();
	void save();
	void save(char*);
	void normalize();
	void normalizeStd();
	void useApproximation(double* dist, int arrayLength);
	void scale();
	void normalizeVariance(HiCBgModel* model);
	void setDefaultVariation();
	void initializeBasedOnChrSizes();
	void setCoverageLimits(double avg, double std);
	HiCBgModelChr* addChrModel(char* chrname);
	void createRandomReads(FILE* outputfp, int totalReads);
	void initializeExpectMatrix(int nPeaks1, int nPeaks2);

};

class ChrTags {
public:
	Tag* tags;
	PETag* petags;
	int totalPositions;
	double totalTags;
	double ogTotalTags;
	float* gcFreq;
	int mCflag;
	int singleFile;
	int memType;
	int loaded;
	int adjusted;
	int optimizedFlag;
	int optimizeOverride;
	float maxtbp;
	float mintbp;
	int tagAdjust;
	int dontSAVE;
	long long int appearentSize;
	int revStrand;
	char* chr;
	char* tagFile;
	FILE* tagfp;
	FILE* tagfpR1;
	FILE* tagfpR2;
	LinkedTag* firstTag;
	LinkedTag* linkTag;
	LinkedPETag* firstPETag;
	LinkedPETag* linkPETag;
	int numLinks;
	int pairedEndFlag;
	int forceSingleReadFlag;
	int maxPosition;
	int maxDistHiC;
	int minDistHiC;
	Hashtable* chrNames;

	char* seq;
	int chrLen;
	static Hashtable* allSeqs;

	
	ChrTags(char* chr);
	~ChrTags();

	void adjustTags();
	void readAndSave();
	void decontaminate(ChrTags* input, int maxDistance, double fraction, double norm);
	void autoCorrelateTags(double* sameStrand, double* diffStrand,int windowSize,
								double maxTags,double &totalCount,double* sameStrandN, double* diffStrandN,int mCflag);
	void annotateTagLocations(ChrPeaks* anns, FILE* annOutput, Doubletable* stats);
	void readTagFile();
	void readTagFile(char* tagfile);
	void addTag(int pos, char dir,int length, float value);
	void addPETag(char* c1, int p1, char d1,int len1, char* c2, int p2, char d2, int len2,float value);
	void printAlignedTag(char* name, int pos, char dir,int length, float value,int PEflag);
	void printAlignedPETag(PETag* petag, int revFlag);
	void optimizeTags();
	void optimizePETags();
	void optimizeTagFile();
	void setMaxTBP(float maxTagsPerBp);
	void setMinTBP(float minTagsPerBp);
	void setTagAdjust(int centerDistance);
	void setTagFile(char* file);
	void openTagFile(char* mode);
	void closeTagFile();
	void print();
	void joinPEreads();
	void findRestrictionSites(char* site, int maxMisMatches,Hashtable* sites, char* genomeDirectory);
	void findPutativePeaks(PeakLibrary* putativePeaks, int peakSize, int minDist, char strand, double minCount);
	void findmCPeaks(PeakLibrary* peaks, int peakSize, char strand, int mCflag, double mCthrsh, 
																		int minNumC, ChrTags* input);
	void normalizeTagCountsGC(NucleotideFreq* nf);
	void getTagCountDistribution(double* d, int max,int scaleFactor);
	void getTagLengthDistribution(double* d, int max);
	void printBedGraph(FILE* fp,double normFactor,char strand,int resolution, int negFlag,
			int fragLength, double reductionRatio,int method, int lastTagFlag, int condenseFlag,
			UniqMapChrs* umc,Peak* circosPeak,ChrTags* inputct, double pseudoCounts, int logFlag,
			double inputNormFactor, int inputFragLength, double normLength);
	Tag* getCoverageTags(int &coveragePositions,double normFactor,char strand,int resolution,
                                int setFragLength,int style,int lastTagFlag,UniqMapChrs* umc,double normLength);
	Tag* getRatioTags(int &ratioPositions, Tag* coverageTags, int coveragePositions, Tag* inputTags, int inputPositions, 
                double pseudoCounts);
	void getPeakTagCounts(Doubletable* counts, ChrPeaks* cp, char strand); // more efficient
	void loadSequence(char* genomeDirectory);
	void assignPETagsToRestrictionSites(char* site, Hashtable* sites,int fragmentLength,int mode,int midFlag,int removeSelfLigationFlag,
									int removeRestrictionEnds,double *posStrand, double* negStrand, int distLength);
	void getSequence(char* dest, int start, int end, char strand);
	void deleteSequence();
	void checkTagSeqBias(char* genomeDirectory,NucleotideFreq* nf,NucleotideFreq* nfuniq,
										NucleotideFreq* nfctrl,int freqStart,int freqEnd,int fragLen,
										OligoArray* oligos, int oligoStart, int oligoEnd);

	PeakLibrary* findGroSeqRegions(int mode, ChrTags* input, UniqMapChrs* umc, char strand,
                            int tssSize, int minBodySize, int maxBodySize, double threshold,
                            double foldStartTranscript, double foldBodyTranscript,double endFold,  double inputFold, 
							double inputNorm, int fragLength, int inputFragLength, double pseudoTags);
	PeakLibrary* findGroSeqTranscripts(int mode, ChrTags* input, UniqMapChrs* umc, char strand,
                            int tssSize, int minBodySize, int maxBodySize, double threshold,
                            double foldStartTranscript, double foldBodyTranscript,double endFold,  double inputFold, 
							double inputNorm, int fragLength, int inputFragLength, double pseudoTags,int groseqMethod);
	float* buildProfile(char strand, double norm, int fragLength, int &maxPosition);
	void normalizeTagCountsFixedOligo(char* genomeDirectory,OligoArray* oligos,int oligoStart, int oligoEnd,
																	double minFold, double maxFold);
	void normalizeTagCountsOligos(char* genomeDirectory,Hashtable* oligos,int oligoLength,
				int regionStart, int regionEnd, double minFold, double maxFold, float maxPerBp,int normFlag);
	int getPETagDistribution(double* sameStrand, double* diffStrand, int windowSize,
                            double * largeWindow,int resolution, int largeLength);
	void getAverageCoverage(int spikeSize, double &chrTotal, int &chrN, PeakLibrary* regions);
	void removeTagsInPeaks(PeakLibrary* regions);
	int getMaxPosition();
	void removePETagBackground(int fragLength);
	void getPETagTotals(PeakLibrary* peaks1, int &peakIndex1);
	void makeHiCMatrix(PeakLibrary* peaks1, int &peakIndex1, PeakLibrary* peaks2,
				double **matrix, int resolution, HiCBgModel* model, int revFlag,
				int numHistBins, int actionFlag,double totalInteractions,
				GenomeInteractionLibrary* gil, pthread_mutex_t* mutex,
				HiCparams* params);
	void scoreInteractionBoundaries(PeakLibrary* peaks1, int &peakIndex1,
				int resolution, HiCBgModel* model,int actionFlag);
	double trimGroSeqTranscript(float* profile, Peak* trans,int strand, int fragLength,double endFold, int minBodySize);

	void loadTags();
	void freeTags();

};
int cmpTags(const void*, const void*);
int cmpPETags(const void*, const void*);
int cmp2ndPETags(const void*, const void*);
int cmp2ndPETagPointers(const void*, const void*);

class UniqMapChrs {
public:
	char* chr;
	unsigned char* pstrand;
	unsigned char* nstrand;
	unsigned int size;
	unsigned int numMappable;
	char* pfile;
	char* nfile;
	int unitSize;
	int maxIndex;

	static unsigned char* lookup;
	static unsigned char* smask;
	static unsigned char* emask;
	static void getUniqMapStats(char* directory, long long int &gsize, long long int &mappability,
									LongInttable* chrSizes, LongInttable* chrUniqMap);

	UniqMapChrs(char* chr, char* directory, int newFlag); // newFlag == 0 if just reading
	~UniqMapChrs();
	void loadData();
	void increaseMapSize(unsigned int size);
	void setMappable(unsigned int pos, int dir);
	void setUnmappable(unsigned int pos, int dir);
	void optimizeSize();
	void saveFiles();
	void initializeLookup();
	void adjustProfile(float* profile, float *umcProfile, char strand, int fraglen, int maxPosition);
	unsigned char countBits(unsigned char x);

	int countRegion(unsigned int start, unsigned int end, char strand);
	Tag* getRegion(unsigned int start, unsigned int end, int &numTags);
	float* buildUniqMapProfile(char strand, int fraglen, int maxPosition);
	unsigned char* getChr(char strand, int &length);

};

class Tag {
public:
	int p; // position
	int len; // position
	float v; // value = number of tags
	char d; // direction 0=+,1=-
	Tag();
	void copy(Tag* src);
	void print(FILE* fp,char* chr);

	static int precision;
};

class LinkedTag : public Tag {
public:
	LinkedTag* tag;
	LinkedTag(int,char,int,float,LinkedTag*);
	~LinkedTag();
};

class PETag {
public: 
	char* name;
	char* chr1;
	int p1; 
	char d1;
	int len1;
	char* chr2;
	int p2; 
	char d2;
	int len2;
	float v;
	PETag();
	PETag(char* name,char* chr, int p, char d, float v,int len);
	~PETag();
	void init();
	void copy(PETag* src);
	void print(FILE*);
	void print2ndHalf(FILE*);
	void print(FILE*,int revFlag);

	static int precision;
};


class LinkedPETag : public PETag {
public:
	LinkedPETag* tag;
	LinkedPETag(char* nc1, int np1, char nd1,int nlen1,char* nc2, int np2, char nd2,int,float,LinkedPETag*);
	~LinkedPETag();
};

class PeakLibrary {
public:
	Hashtable* chrs; // holds pointers to ChrPeak objects
	Hashtable* peaks; // holds pointers to Peak objects
	char* name;
	char* genome;
	int numPeaks;
	double tagsInPeaks;
	double avgPeakSize;
	int fixedFlag;
	int duplicateWarningFlag;
	Peak** peakOrder;
	pthread_mutex_t mutex;

	Inttable* duplicates;
	Doubletable* chrSizes;

	TagLibrary** exps;
	int numExps;


	PeakLibrary();
	PeakLibrary(char* file, int mode);
	PeakLibrary(int expectedNumberOfPeaks); // default is 100000
	void initialize(int expectedNumPeaks);
	~PeakLibrary();
	void print(FILE*);
	void printSNPs(FILE*);
	void printSNPtotals(FILE*);
	void printSorted(FILE*);
	void printGTF(FILE*);
	void printAnnotation(FILE*);
	void readPeakFile(char* filename,int mode);
	void setDefaultPeakOrder();
	void setMaxChrPositions(Doubletable* nchrSizes); //pass NULL to trim peak positions using the current sizes
	int addTagLibrary(TagLibrary* t);
	PeakLibrary* copy();
	PeakLibrary* stitchRegions(int maxDistance,int mode);
	PeakLibrary* getDifferentialPeaks(TagLibrary* tags, TagLibrary* input,
                        double foldThreshold, double poissonThresh, int mode, int start, int end, char strand,int strFlag);
	PeakLibrary* filterLocalPeaks(TagLibrary* tags, int peakSize, int localSize,
                        double threshold, double poissonThresh, int mode, char strand);
	PeakLibrary* filterClonalPeaks(TagLibrary* tags, int peakSize,
                        double threshold, int mode, char strand);
	void centerPeaks(TagLibrary* tags, int peakSize,char strand);
	void analyzeTSSpattern(TagLibrary* tags, int peakSize,char strand,int mode);
	void analyzeReadAutocorrelation(TagLibrary* tags, int peakSize,char strand,int maxSize,int mode);
	void centerNFR(TagLibrary* tags, int peakSize,char strand,int nfrSize);
	void getMappabilityCount(char* uniqMapDirectory, int tagAdjust, char strand, 
								long long int& possible, long long int& totalMappable);
	void mergePeaks(char** peakfiles, int numFiles, char strand, int maxDistance,
					long long int gsize, char* prefix, char* vennFile, char* matrixFile,int codeFlag,char* cmdline);
	void getCoBoundPeaks(char** peakFiles, int numPeakFiles, char strand, int maxDistance, 
					long long int gsize, char* prefix, int maxCoBound, char* matrixFile,char* cmdline);
	void genomeOntology(char** peakFiles, int numPeakFiles, char strand, int maxDistance,
						long long int gsize, char* controlFileName);
	Doubletable* reportTotalValuesPerChr();
	PeakLibrary* mergeNearbyPeaks(int maxDistance,char strand);
	PeakLibrary* filterPeaksOutsideRange(char* chr, int start, int end);

	Hashtable* getOverlappingPeaks(PeakLibrary*,char strand,int maxDistance,long long int &overlap); 
																	//returns Hashtable of 'PeakMapping'
	void annotatePeakLocations(PeakLibrary* anns, FILE* statsOutput, FILE* annOutput);

	Doubletable* countPeakTags(int expIndex, int start, int end, char direction,int mode);
	Doubletable* countPeakTagsLowMemory(TagLibrary* tags, char direction,int mode);

	Doubletable* scoreNFR(TagLibrary* tags,int expIndex, int nfrSize,char direction);
	void printRelativeTags(FILE* fp, int expIndex, int start, int end,int mode);
	void setPeakTagSizeRefPos(int offset, int startOffset, int endOffset);
	void setPeakTagSizeFixed(int startOffset, int endOffset);
	int getAveragePeakSize();
	PeakLibrary* getCoveragePeaks(int res, int superRes);

	PeakLibrary* prioritizeAnnotations();
	void sortPeakTags(int expIndex);
	Peak* addPeak(char* name, char* chr,int start, int end, int midpoint, char dir, float value, 
						float ratio, char* extradata, int mappability, unsigned int priority);
	Peak* checkForPeak(char* name, char* chr, int start, int end, char strand, int & info);
	void normalizePeakScore(float normFactor);
	void addPeak(Peak* p);
	void addPeakLibrary(PeakLibrary* p);
	void sortChr();
	PeakLibrary* getSuperEnhancers(double slope, int window,char* &notes,PeakLibrary* &typical);
	void sortKeys(char**);
	long long int calculateCoverage();
	void setPeakSize(int size);
	void extractSequence(char* genomeDir, FILE* output, int fastaFlag,int maskFlag);
	static void extractSequenceStats(char* genomeDir, FILE* output);
	static void printSequenceStats(FILE* output, char* chrName, long long int totalSeq, 
							long long int totalN,long long int totalGC);
	SNP* addSNPsfromVCF(char* vcfFile, char** individuals, int &numIndividuals, int allFlag);
	double adjustPETagTotalsForModel(HiCBgModel* model, int maxCPUs);
	double adjustPETagTotalsWithModel(HiCBgModel* model, int useTotalsTooFlag);
	void processBedGraphFile(char* bedGraphFile, int relStart, int relEnd, char strand, FILE* outputFile, int outputMode);
	void processWiggleFile(char* wiggleFile, int relStart, int relEnd, char strand, FILE* outputFile, int outputMode);
	void checkForOutOfBoundsCoordinates(TagLibrary* ct);
};

class ChrPeaks {
public:	
	Peak** peaks;
	int numPeaks;
	int maxPossible;
	LinkedList* peakList;
	pthread_mutex_t mutex;
	pthread_mutex_t mutex2;
	int mutexIndex;
	int mutexTotal;

	ChrPeaks();
	~ChrPeaks();
	void addPeak(Peak* p);
	void sort();
	void print(FILE* fp);
	void addTagLibrary(ChrTags* t,int libraryIndex);
	void prioritizeAnnotations(PeakLibrary* p);
	void annotatePeakLocations(ChrPeaks* anns, FILE* annOutput, Doubletable* stats);
	void getAnnotationSizeTotals(Doubletable* sizeTotals);
	void getOverlappingPeaks(Hashtable* results, ChrPeaks*, char strand,int maxDistance,long long int &overlap);
	void stitchRegions(PeakLibrary* regions, int maxDistance, int mode);
	void extractSequence(FILE* fpin, FILE* fpout, int fastaFlag,int maskFlag, char* buf1, char* buf2);
	void getMappabilityCount(UniqMapChrs* umc, int tagAdjust, char strand);
	void mergeNearbyPeaks(PeakLibrary* out, int maxDistance,char strand);
	long long int calculateCoverage();
	double adjustTotalsBasedOnModel(PeakLibrary* refPeaks, int resolution, HiCBgModel* model,
									double* tmpTotals,double* errors, int maxCPUs);
	void adjustTotalsThread(double* errors, double* tmpTotals, double expectInter,
                    int minIndex,int maxIndex, double totalModelReads,
                    int resolution,HiCBgModelChr* chrModel,PeakLibrary* refPeaks);
	void countPeakTagsLowMemory(Doubletable* results, ChrTags* ct,char direction, int mode);

	void checkForOutOfBoundsCoordinates(ChrTags* ct);
	
};

class Peak {
public:
	char* name;
	char* chr;
	int refPos;
	int start;
	int end;
	int tagStart;
	int tagEnd;
	int fixedFlag;
	int index;
	unsigned int priority;
	char strand;
	char* seq;
	int uniqMap;
	float v;
	float focusRatio;
	char* data;
	char* ogname;
	PeakSNPs* snps;
	

	Tag** exps;
	int* numTags;
	int* maxTags;
	int numExps;

	Peak();
	void init();
	Peak(char* name,char* originalName, char* chr, int start, int end,int midpoint, char dir, float value, 
				float focusRatio, char* otherdata,int mappability,unsigned int priority);
	~Peak();
	double countTags(int expIndex, int start, int end, char direction,int mode);
	double scoreNFR(int expIndex, int fragLength, int nfrSize, int nucSize, char direction);
	Tag* getCoverageTags(int expIndex, int fragLength, int& coveragePositions,char strand);
	void centerPeak(int expIndex,int fragLength,char strand);
	void analyzeTSSpattern(int expIndex,int fragLength,char strand,int mode);
	void analyzeReadAutocorrelation(int expIndex, int maxSize,char strand,int mode);
	double calculateTSSdispersion(int expIndex, int strandInfo, int focusRadius);
	double calculateTSSperiodic(int expIndex, int strandInfo, int focusRadius);
	double autoCorrelateReads(int expIndex, int strandInfo);
	void centerNFR(int expIndex,int fragLength,char strand,int nfrSize,int nucSize);
	void printRelativeTags(FILE* fp, int expIndex, int start, int end,int mode);
	void print(FILE* fp);
	void printSNPs(FILE* fp);
	void printSNPtotals(FILE* fp);
	void printGTF(FILE* fp);
	void printAnnotation(FILE*);
	void addExp();
	void addTag(Tag* t, int expIndex);
	void sortTags(int expIndex);
	void setPeakSize(int size);
	void setPeakTagSizeRefPos(int newOffset,int startOffset,int endOffset);
	void setPeakTagSizeFixed(int startOffset,int endOffset);
	void setOffset(int offset);
	void addData(char* str);
	void setSeq(char* seqstr);
	void addSNP(SNP* snp, int pos);
	Peak* copy();
};

class MergePeak {
public:
	Peak*** exps;
	int* numPeaks;
	int numExps;
	Peak* merged;

	MergePeak(int numExps);
	~MergePeak();
	void add(MergePeak*);
	void add(int index, Peak*);
	void combine();
};
class MergePeaksArray {
public:
	MergePeak** p;
	int n;
	MergePeaksArray();
	~MergePeaksArray();
};
class PeakMapping {
public:
	Peak* source;
	Peak** peaks;
	int* distance; 
	int numPeaks;
	int bp;
	
	PeakMapping(Peak* source);
	~PeakMapping();
	void add(Peak*p, int dist);
	void calculateCoverage();
};

class GenomeOntologyResults {
public:
	char* name;
	double tagCount;
	long long int coverage;
	int numPeaks;
	int numPeaksReporting;
};

class OligoProfile {
public:
	double *pProfile;
	double *nProfile;
	int offset;
	int length;
	double N;
	char revoppFlag;
	char normFlag;
	OligoProfile(int length, int offset);
	~OligoProfile();
	void print(FILE* fp);
	void normalize();
	static void mergeRevopps(OligoProfile*,OligoProfile*);
};

class SNP {
public:
	char** genotypes;
	char** rvGenotypes;
	int numGenotypes;

	int numIndividuals;
	char* allele1;
	char* allele2;

	SNP* nextSNP;
	int count;
	double* editDistances;

	SNP();
	SNP(char* ref, char** alt, int numAlt, int numIndividuals);
	~SNP();
	static void deleteBaseSNP(SNP*);
	void init();
	void print(FILE* fp,char strand);
	void printReference(FILE* fp,char strand);
	
	void addEditDistances(double* editDistanceTotals);
	void calculateEditDistances();
};

class PeakSNPs {
public:
	SNP** snps;
	int numSNPs;

	int* offsets;

	PeakSNPs();
	~PeakSNPs();
	void init();
	void addSNP(SNP* snp, int pos);
	void print(FILE* fp,char strand);
	void printTotals(FILE* fp);
};

class Genome3D {
public:
	double** matrix;
	int numMatrix;
	char** names;
	int *color;

	double minDist;
	double maxDist;
	
	double** p;
	int* inRange;
	int numInRange;

	double fixedStep;
	double fixedStepMax;
	double fixedStepMin;
	

	Genome3D();
	~Genome3D();
	void init(double **m, int N, char** name, double fixedStep);
	void initStructure(int method);
	void addColorsFromClusters(char* clusterFile);
	void normalizeMatrix(int method);
	void print(FILE* fp);
	void calcMoveOne(int index, int method);
	double checkConnected(double* p0, double* v, double magnitude, double* p1);
	double checkCollision(double* p0, double* v, double magnitude, double* p1);
	void optimize(int maxIterations);
};

class HiCparams {
public:
	int minDist;
	int maxDist;
	int res;
	int superRes;
	int logFlag;
	int boundaryScale;
	double minExpReads;
	int relativeFlag;
	int fragLengthEstimate;


	HiCparams();
	~HiCparams();
	void init();
	void setRes(int res, int superRes);
};

class PeakmC {
public:
	int numC;
	double mCavg;
	double mCstd;
	int start;
	int end;
	int gstart;
	int gend;
	PeakmC();
	void calcAvg();
	~PeakmC();
};
int cmpPeakmC(const void*, const void*);


// support classes

class doubleIndex {
public:
	double v;
	double vp;
	unsigned int index;
	int position;
};

int cmpDouble(const void* a, const void* b);	
int cmpDoubleIndex(const void* a, const void* b);	
int cmpPeaks(const void*, const void*);
int cmpInteractions(const void*, const void*);
int cmpInteractionLocation(const void*, const void*);
int cmpInteractionDiff(const void*, const void*);
int cmpInteractionIndexes(const void*, const void*);
void split(char* string, char** cols, int &numCols, char delim);
int checkInt(char* str);
int checkStrand(char* str);
void cleanUpSeq(char* seq);
void revopp(char* seq);
void reverse(char* str);
char* unzipFileIfNeeded(char* file, int &zipFlag, int &format);
void rezipFileIfNeeded(char* file, int zipFlag);

int chrcmp(const void* chr1, const void* chr2);
void parseUCSCpositionStr(char* str, char* &chr, int &start, int &end);
char* getCMDstr(int argc, char** argv);


int determineSeqFileFormat(char* file);

class NucleotideFreq {
public:

	double **freq;
	double **difreq;
	double* totals;
	double* ditotals;
	int length;
	int fragLength;
	int alphaSize;
	int diAlphaSize;
	int offset;
	char* alpha;
	FILE* gcFile;

	double *gcDist;
	double *gcOG;
	double *gcNorm;
	int gcDistTotal;
	double gcInc;
	
	int *monoLookup;
	int **diLookup;

	NucleotideFreq();
	NucleotideFreq(int offset, int maxLength,int fragLength);
	NucleotideFreq(char* file, int format, int offset, int maxLength,char* gcFilename);
	~NucleotideFreq();
	void init();
	void initMatrix(int,int gclen);
	void addSequence(char* name, char* seq,double w,int gcStart,int gcEnd, SeqFreqStats* sfs);
	void addChr(char* seq,int gcStart, int gcEnd);
	void print(FILE*);
	NucleotideFreq* copy();
	void printGCnormFile(char* filename);
	double printGC(FILE*);
	void readGCcontentFile(char* filename);
	double calculateGCNormalization(NucleotideFreq *ctrl, double minfold, double maxFold, 
															double window,char* filename);
	
};

class SeqFreqStats {
public:
	double* mono;
	double CpG;
	double GC;
	double AG;
	double AC;
	double N;
	double NN;
	SeqFreqStats();
	~SeqFreqStats();
};

#endif
