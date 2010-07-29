/******************************
 *
 * MDmodule: find motifs from one set of expression data
 * Components:
 * Do Monte Carlo simulation to get motif score distribution
 * Screen out potential bad sites (repeats, N's)
 * Candidate motif finding from top sequences
 * Remove similar motifs and keeping the top ones
 * Confirm motifs from "confirm" top sequences
 * Combine similar motifs and report top distinct ones
 * Do MatrixScan for each reported motif
 * 
 * Xiaole Liu, Stanford Medical Informatics
 *
 ******************************/

/******************************
 *
 * Header files
 *
 ******************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>

/******************************
 *
 * Constants
 *
 ******************************/

/* Specify the alphabet */
#define ALPHASIZE 4
#define A 0
#define C 1
#define G 2
#define T 3
#define N 4
#define ALLALPHA "NACMGRSVTWYHKDB"

/* these mostly serve as flags */
#define FORWARD 1
#define BACKWARD 2
#define NO 0
#define YES 1
#define HASN 2
#define MASK 3
#define SEQ 1
#define ADD 1
#define SUBTRACT -1
#define SCANTOP 1

/* for random number generator */
#define IA 16807
#define IM 2147483647
#define AM (1.0 / IM)
#define IQ 127773
#define IR 2836

/* Markov degree for input and background, general formula: 3 * 256 * 4^mkv */
#define MKV1 3 * 1024
#define MKV2 3 * 4096
#define MKV3 3 * 16384

/* some ascii for keys */
#define ENDLINE '\0'
#define TAB 9
#define SPACE 32
#define TILTA 126

/* max/min numbers */
#define MAXLINE 16383
#define REGLINE 1024
#define MINW 5
#define MAXW 25
#define MINM 20
#define MAXM 200
#define MAXT 250
#define MAXS 100
#define MAXR 50
#define MAXU 100
#define MINE 30
#define MAXE 6000
#define MAXK 3
#define MINH 5
#define MAXH 6000
#define MAXCHILD 50

/* some input defaults */
#define DEFAULTW 10
#define DEFAULTT 10
#define DEFAULTS 30
#define DEFAULTR 5
#define DEFAULTU 15
#define DEFAULTH 300

/* other useful constants */
#define SIGNIFICANCE 3
#define VERYSMALL 0.000000025
#define PI 3.141592653589793
#define MPH 60
#define SPH 3600
#define LOG2 log(2)

/******************************
 *
 * Data type definition
 *
 ******************************/

/******************************
 * Data type: inputParam
 * It stores many global variables, so it is easy to pass them around.
 ******************************/
struct inputParam {
  int w, mw, noseq, top, confirm, scan, report, iterate, bgmkv, seqmkv, back, 
    expMt, MonteCarlo, child, norand, printAlign, hits;
  double scount[ALPHASIZE+1], minth, expect, seedmw, mean, stdev;
  char *seed;
  /* file pointers for input sequence, background sequence, background 
     distribution, output file, output with MtfReg scoring */
  FILE *isp, *bsp, *bfp, *ofp, *xfp;
};

/******************************
 * Data type: seqCount
 * Stores the count of letters from either input or background sequences
 ******************************/
struct seqCount {
  int one[ALPHASIZE+1], two[ALPHASIZE+1][ALPHASIZE+1], 
    three[ALPHASIZE+1][ALPHASIZE+1][ALPHASIZE+1], 
    four[ALPHASIZE+1][ALPHASIZE+1][ALPHASIZE+1][ALPHASIZE+1];
};

/******************************
 * Data type: freq
 * With seqCount, now we could calculate the frequency (in terms of 
 * cumulative density function or log of probability) of letters from 
 * input or background sequences in certain Markov order.
 ******************************/
struct freq {
  double one[ALPHASIZE], two[ALPHASIZE][ALPHASIZE], 
    three[ALPHASIZE][ALPHASIZE][ALPHASIZE], 
    four[ALPHASIZE][ALPHASIZE][ALPHASIZE][ALPHASIZE];
};

/******************************
 * Data type: sequence
 * Store the count of letters from either input or background sequences
 ******************************/
struct sequence {
  char *name;
  int *seqi, *rseqi, *tried1, *tried2, len;
  double *bgs, *rbgs;
  struct sequence *next;
};

/******************************
 * Data type: mtfPt
 * Motif points make up motif blocks, each point contains a count and a
 * log specifying the position specific motif matrix. 
 ******************************/
struct mtfPt {
  int ct;
  double log;
};

/******************************
 * Data type: motif
 * Stores the position specific matrices, number of aligned segments and
 * some arrays for storing doubles during update calculations. 
 ******************************/
struct motif {
  struct mtfPt **blk;
  double score, bginfo, zscore, avgInfo, siteTh;
  char *constr, *rconstr, *degstr, *rdegstr;
  int *conint, *rconint;
  struct alignment *ahead, *atail;
  int segment, oldsegment;
};

/******************************
 * Data type: alignment
 * stores the alignment parameters as linked list (each node is for one
 * sequence) in the motif structure.
 ******************************/
struct alignment {
  struct sequence *s;
  int position, *st;
  double bginfo;
  struct alignment *next;
};

/******************************
 *
 * External functions
 *
 ******************************/

int fflush (FILE *stream);
int getopt(int, char *const *, const char *);
int fork(void);
int pipe(int pd[2]);
int read(int fd, void *buf, unsigned len);
int write(int fd, void *buf, unsigned len);

/******************************
 *
 * Function prototypes (in alphabetical order)
 *
 ******************************/

void addAlign(struct sequence *s, int pos, int *st, double bginfo, struct 
	      motif *mtf);
void addSegment(struct motif *mtf, int *start, int w, double bginfo, int add);
double blkMtfScore(struct mtfPt **blk, struct inputParam *input, int wid, 
		   int segment);
void calcBlkLog(struct mtfPt **blk, int segment, int wid, double *sct);
void calcSeqCDF(struct freq *seqlog, struct freq *seqcdf, int mkv);
void calcSeqLog(struct freq *seqlog, struct seqCount *seqct, int *mkv);
void calcStats(double *stats, struct inputParam *input);
void clearMotif(struct motif *mtf, struct inputParam *input);
void countSeq(int *seqint, int len, struct seqCount *seqct);
void createSeq(int isSeq, struct sequence **seqs, struct seqCount*seqct,
	       char *name, char *seq, struct inputParam *input, struct
	       sequence **tail);
double drand();
void compactMotifs(struct motif **motifs, struct inputParam *input);
void findCandidateMtf(struct inputParam *input, struct sequence *seqs, struct
		      motif **motifs, struct motif **mtf, struct freq *bglog);
void findMotifTh(struct inputParam *input, struct sequence *seqs, 
		 struct motif **motifs, struct motif **mtf, 
		 struct freq *bglog, struct freq *seqcdf);
int findMin(double *list, int size);
void finishUp(struct timeval *begintv, FILE *ofp);
void enumerateSeed(struct motif **mtf, int *seq, struct sequence *seqs, 
		   struct inputParam *input, int *try, int end, struct motif 
		   **motifs, struct freq *bglog);
void errorExit(char *msg);
void genRandSeq0(int *seq, int len, struct freq *schar);
void genRandSeq1(int *seq, int len, struct freq *schar);
void genRandSeq2(int *seq, int len, struct freq *schar);
void genRandSeq3(int *seq, int len, struct freq *schar);
void getConsensus(int w, struct motif *mtf);
void getConInt(int w, int mw, struct motif *mtf);
void initMotif(struct motif **mtf, int w);
void initRand();
void insertSeed(struct motif **motifs, struct motif **mtf, struct 
		inputParam *input);
int match(int *x, int *y, int wid, int mis);
double motifInfo(struct motif *mtf, struct inputParam *input, struct freq *bglog);
int motifsSimilar(struct motif *mtf1, struct motif *mtf2, int w, int mw);
int mystrlen(char *string);
void parseInput(struct inputParam *input, int argc, char *argv[]);
void printAlignment(struct alignment *align, struct inputParam *input);
void printShortAlignment(struct sequence *seqs, struct motif *mtf,
			 struct inputParam *input);
void printBgFreq(FILE *fp, struct freq *bg, int mkv);
void printCandidates(struct inputParam *input, struct motif **motifs);
void printMotif(struct motif *mtf, struct inputParam *input, int mtfCt);
void printResults(struct motif **motifs, struct inputParam *input, struct
		  sequence *seqs);
void printUsage();
void randomize(struct sequence *seqs, struct freq *seqcdf, int mkv);
void readBgFreq(FILE *fp, struct freq *bglog, int *mkv);
int readBgFreqLine(FILE *bf, double *bgfreq);
void readInputFiles(struct inputParam *input, struct freq *seqcdf,
		    struct freq *bglog, struct sequence **seqs); 
void readSeqFile(FILE *fp, int isSeq, struct sequence **seqs, struct 
		 seqCount *seqct, struct inputParam *input);
void reCalcScount(int segment, double *sct, struct freq *bglog);
void refine(struct inputParam *input, struct motif **motifs, struct freq
	    *bglog);
void scanning(struct sequence *seqs, struct inputParam *input, struct
	      motif **motifs, struct freq *bglog);
void scoreAll(struct inputParam *input, struct motif **motifs, 
	      struct sequence *seqs);
int seedMatch(char *x, int *y, int wid, double mis);
void seqatoi(char *line, int *seqint, int len, int forward);
double seqBgScore(int *ptr, int st, int wid, int mkv, struct freq *bg);
double seqMtfScore(int *seqi, int wid, struct mtfPt **blk);
int similar(struct motif *mtf1, struct motif *mtf2, int w);
void tryOneSeed(struct motif **mtf, struct inputParam *input, 
		struct sequence *seqs, struct motif **motifs, 
		struct freq *bglog);
void wholeBgScore(struct sequence *seqs, struct freq *bglog, struct 
		  inputParam *input);

/******************************
 *
 * global variables
 *
 ******************************/

struct inputParam input;
struct freq bglog, seqcdf;
struct sequence *seqs;
struct motif *mtf, **motifs;
long idum;  /* for random number generator */

/******************************
 *
 * Functions (in the order of which they are called)
 *
 ******************************/

/******************************
 *
 * Func: main
 *
 ******************************/
int main(int argc, char *argv[]) {
  struct timeval begintv;
  int i;

  gettimeofday(&begintv, NULL);
  initRand();

  /* read the input and initialize variables */
  parseInput(&input, argc, argv);
  if (!input.back)
    printf("Read input.\n");
  readInputFiles(&input, &seqcdf, &bglog, &seqs);
  motifs = (struct motif **) malloc((input.scan)*sizeof(struct motif *));
  for (i = 0; i < input.scan; i++)
    initMotif(&(motifs[i]), input.w);
  initMotif(&mtf, input.w);

  if (input.MonteCarlo > 0) {
    if (!input.back)
      printf("Monte Carlo simulation to get motif score distribution.\n");
    findMotifTh(&input, seqs, motifs, &mtf, &bglog, &seqcdf);
  }
  if (input.MonteCarlo) 
    fprintf(input.ofp, "Null motif score distribution mean: %.3f, standard deviation: %.3f\n", input.mean, input.stdev);

  /* score input sequence with bglog once (can reuse the score) */
  wholeBgScore(seqs, &bglog, &input);
  
  if (!input.back)
    printf("Find candidate motifs.\n");
  findCandidateMtf(&input, seqs, motifs, &mtf, &bglog);
  /*
  printCandidates(&input, motifs);
  */

  if (!input.back)
    printf("Scan sequences using candidate motifs.\n");
  scanning(seqs, &input, motifs, &bglog);
  /*
  printCandidates(&input, motifs);
  */

  if (!input.back)
    printf("Refine candidate motifs.\n");
  for (i = 0; i < input.iterate; i++) {
    refine(&input, motifs, &bglog);
  }
  compactMotifs(motifs, &input);

  if (input.xfp) {
    if (!input.back)
      printf("Score all sequences with motifs.\n");
    scoreAll(&input, motifs, seqs);
    if (!input.back)
      printf("Print results.\n");
    printResults(motifs, &input, seqs);
  } else {
    /*
      printf("print %d motifs\n", input.report);
      for (i = 0; i < input.scan; i++) {
      printf("Motif %d, Segment %d, Score %.3f\n", i+1, motifs[i]->segment,
	     motifs[i]->score);
	     }
    */
    if (!input.back)
      printf("Print results.\n");
    printResults(motifs, &input, seqs);
    finishUp(&begintv, input.ofp);
  }

  return 0;
}

/******************************
 *
 * Func: parseInput
 * Read in all the user specified parameters
 *
 ******************************/
void parseInput(struct inputParam *input, int argc, char *argv[]) {
  extern char *optarg;
  int inttemp, opt, i, len;
  double dtemp;
  
  /* init input */
  memset(input, 0, sizeof(struct inputParam));
  input->w = DEFAULTW;
  input->top = DEFAULTT;
  input->scan = DEFAULTS;
  input->report = DEFAULTR;
  input->iterate = DEFAULTU;
  input->hits = DEFAULTH;
  input->printAlign = 0;
  input->child = 1;
  input->bgmkv = -1;
  input->seed = NULL;
  input->isp = input->bsp = input->bfp = input->xfp = NULL;
  input->ofp = stderr;
  
  /* if no arguments, print out usage */
  if (argc == 1)
    printUsage();
  
  /* read input parameters and specify input */
  while ((opt = getopt(argc, argv, "i:w:f:b:k:e:t:s:c:r:a:n:m:M:C:p:o:x:g:u:v:h:")) 
	 != EOF) {
    switch(opt) {
      
      /* input sequence file */
    case 'i': if (!(input->isp = fopen(optarg, "r"))) 
      errorExit("Can't open input sequence file.");
    break;
    
    /* motif width */
    case 'w': if (sscanf(optarg, "%d", &inttemp) && inttemp >= MINW
		  && inttemp <= MAXW)
      input->w = inttemp;
    else
      errorExit("Wrong motif width input.");
    break;
    
    /* background distribution file */
    case 'f': if (!(input->bfp = fopen(optarg, "r"))) 
      errorExit("Can't open background frequency file."); 
    break;
    
    /* background sequence file */
    case 'b': if (!(input->bsp = fopen(optarg, "r"))) 
      errorExit("Can't open background sequence file."); 
    break;
    
    /* markov order */
    case 'k': if (sscanf(optarg, "%d", &inttemp) && inttemp >= 0
		  && inttemp <= MAXK)
      input->bgmkv = inttemp;
    else
      errorExit("Wrong Markov background order.");
    break;

    /* expected bases/site */
    case 'e': if (sscanf(optarg, "%lf", &dtemp) && dtemp >= MINE
		  && dtemp <= MAXE)
      input->expect = dtemp;
    else
      errorExit("Wrong expected bases/site.");
    break;
    
    /* top sequences for seeds */
    case 't': if (sscanf(optarg, "%d", &inttemp) && inttemp > 0
		  && inttemp <= MAXT)
      input->top = inttemp;
    else
      errorExit("Wrong top sequences.");
    break;

    /* top motif seeds */
    case 's': if (sscanf(optarg, "%d", &inttemp) && inttemp > 0
		  && inttemp <= MAXS)
      input->scan = inttemp;
    else
      errorExit("Wrong top motif candidates.");
    break;
    
    /* top sequences to confirm */
    case 'c': if (sscanf(optarg, "%d", &inttemp) && inttemp > 0)
      input->confirm = inttemp;
    else
      errorExit("Wrong confirm sequences.");
    break;

    /* number of final refinement iterations */
    case 'n': if (sscanf(optarg, "%d", &inttemp) && inttemp > 0
		  && inttemp <= MAXU)
      input->iterate = inttemp;
    else
      errorExit("Wrong top motif candidates.");
    break;

    /* top motifs to report */
    case 'r': if (sscanf(optarg, "%d", &inttemp) && inttemp > 0
		  && inttemp <= MAXR)
      input->report = inttemp;
    else
      errorExit("Wrong top final motifs to report.");
    break;

    /* print out alignment */
    case 'a': if (sscanf(optarg, "%d", &inttemp) && inttemp >= 0 
		  && inttemp < 2)
      input->printAlign = inttemp;
    else
      errorExit("Wrong print out alignment.");
    break;

    /* number of hits (sites) */
    case 'h': if (sscanf(optarg, "%d", &inttemp) && inttemp >= MINH
		  && inttemp <= MAXH)
      input->hits = inttemp;
    else
      errorExit("Wrong number of site hits.");
    break;
    
    /* Monte Carlo simulation for motif score distribution */
    case 'M': if (sscanf(optarg, "%d", &inttemp) && inttemp >= MINM
		  && inttemp <= MAXM)
      input->MonteCarlo = inttemp;
    else
      errorExit("Wrong number of Monte Carlo simulation runs.");
    break;
    
    /* number of child processes to run Monte Carlo simulation */
    case 'C': if (sscanf(optarg, "%d", &inttemp) && inttemp > 0
		  && inttemp <= MAXCHILD)
      input->child = inttemp;
    else
      errorExit("Wrong number of child processes for Monte Carlo simulation.");
    break;
    
    /* if known mtf score dist stdev */
    case 'v': if (sscanf(optarg, "%lf", &dtemp) && dtemp > 0)
      input->stdev = dtemp;
    else
      errorExit("Wrong motif score distribution stdev.");
    break;
    
    /* if known mtf score dist mean */
    case 'u': if (sscanf(optarg, "%lf", &dtemp) && dtemp > 0)
      input->mean = dtemp;
    else
      errorExit("Wrong motif score distribution mean.");
    break;
    
    /* predefined motif seed */
    case 'm': len = strlen(optarg);
      if (len < MINW || len > MAXW)
	errorExit("Motif consensus or degenerate has wrong width.");
      for (i = 0; i < len; i++) {
	if (!strchr(ALLALPHA, optarg[i]))
	  errorExit("Wrong motif consensus or degenerate.");
      }
      input->seed = (char *)malloc(len + 1);
      strcpy(input->seed, optarg);
    break;

    /* experiment number */
    case 'p': if (sscanf(optarg, "%d", &inttemp))
      input->expMt = inttemp;
    break;
    
    /* output file name */
    case 'o': if (!(input->ofp = fopen(optarg, "w"))) 
      errorExit("Can't open output file.");
    break;

    /* matrixScan output name */
    case 'x': if (!(input->xfp = fopen(optarg, "w"))) 
      errorExit("Can't open output file.");
    break;

    /* background run */
    case 'g': if (sscanf(optarg, "%d", &inttemp) && inttemp >= 0 
		  && inttemp < 2)
      input->back = inttemp;
    else
      errorExit("Wrong background run.");
    break;
    }
  }
    
  /* check for invalid input */ 
  if (!input->isp) 
    errorExit("No input sequence file specified.");

  if (input->seed) {
    input->w = len;
    input->scan = input->report = 1;
  }
  input->minth = log(10 * input->w);

  if (input->printAlign && input->xfp) {
    input->printAlign = 2;
  }

  /* If the user wants to find motif score distribution, find out the number 
   * of randomization to do. Otherwise, find out the motif score threshold 
   * in case mean and std are specified */
  if (input->mean && input->stdev) {
    input->MonteCarlo = -1;
    /* printf("Mean %lf, Stdev %lf\n", input->mean, input->stdev); */
  }
  if (input->MonteCarlo > 0) {
    input->norand = (int) ceil((double)input->MonteCarlo / input->child);
    input->MonteCarlo = input->norand * input->child;
  }
  if (input->confirm && input->confirm < input->top)
    input->confirm = input->top;
}

/******************************
 *
 * Func: printUsage
 *
 ******************************/
void printUsage() {
  fprintf(stderr, "MDmodule usage:\n");
  fprintf(stderr, "\t-i <input sequences>\n");
  fprintf(stderr, "\t-w <motif width (default %d)>\n", DEFAULTW);
  fprintf(stderr, "\t-f <background frequency file>\n");
  fprintf(stderr, "\t-b <background sequence file (default input sequences)>\n");
  fprintf(stderr, "\t-k <background Markov dependency order (default don't specify)>\n");
  fprintf(stderr, "\t-e <expected bases per motif site in the top sequences \n\t    (will use Bayes motif scoring, don't specify if unkown,\n\t     which uses Liu motif scoring)>\n");
  fprintf(stderr, "\t-t <number of top sequences to look for candidate motifs (default %d)>\n", DEFAULTT);
  fprintf(stderr, "\t-s <number of candidate motifs to scan and refine (default %d)>\n", DEFAULTS);
  fprintf(stderr, "\t-c <number of top sequences to confirm candidate motifs \n\t    (default whole dataset)>\n");
  fprintf(stderr, "\t-n <number of refinement iterations (default %d)>\n", DEFAULTU);
  fprintf(stderr, "\t-r <number of top motifs to report at the end (default %d)>\n", DEFAULTR);
  fprintf(stderr, "\t-a 1 <if want to print out alignments (default just motif matrics)>\n");
  fprintf(stderr, "\t-h <number of site hits (default %d)\n\t    only works if also specified -a and -x>\n", DEFAULTH);
  fprintf(stderr, "\t-M <Monte Carlo simulation for motif score distribution\n\t    (default no Monte Carlo or stat significance score)>\n");
  fprintf(stderr, "\t-C <number of child processes to run Monte Carlo simulation\n\t    (default one process without no Monte Carlo simulation)>\n");
  fprintf(stderr, "\t-u Motif score distribution mean\n\t    (default no stat significance score or perform Monte Carlo)>\n");
  fprintf(stderr, "\t-v Motif score distribution mean\n\t    (default no stat significance score or perform Monte Carlo)>\n");
  fprintf(stderr, "\t-m <predefined motif consensus or degenerate (default do not specify)>\n");
  fprintf(stderr, "\t-p <+/- experiment number (default do not specify)>\n");
  fprintf(stderr, "\t-o <output file (default stdout)>\n");
  fprintf(stderr, "\t-x <matrixscan output file (default no matrixscan)>\n");
  fprintf(stderr, "\t-g 1 <if you don't want to see messages during the run>\n");
  fflush(stderr);
  errorExit("");
}

/******************************
 *
 * Func: errorExit
 * Print the msg, and exit the program 
 *
 ******************************/
void errorExit(char *msg) {
  fprintf(stderr, "%s\n", msg);
  exit(0);
}

/******************************
 *
 * Func: readInputFiles
 * Read probability of sequence input/bg files, initialize seqs, bg/seq
 * log  
 * Note: It is important to read Seq file before Bg file, because if the
 * user does not specify the Bg file, then information in the Seq file
 * (stored in seqCt) will be used to calculate bglog.
 *
 ******************************/
void readInputFiles(struct inputParam *input, struct freq *seqcdf,
		    struct freq *bglog, struct sequence **seqs) {
  int i;
  double temp;
  struct seqCount seqCt;
  FILE *fp = NULL;

  /* read in input sequences */
  readSeqFile(input->isp, SEQ, seqs, &seqCt, input);

  /* calculate seqcdf if user wants to do Monte Carlo simualation */
  if (input->MonteCarlo > 0) {
    calcSeqLog(seqcdf, &seqCt, &(input->seqmkv));
    calcSeqCDF(seqcdf, seqcdf, input->seqmkv);
  }

  if (input->top > input->noseq)
    input->top = input->noseq;
  if (input->confirm) {
    if (input->confirm > input->noseq)
      input->confirm = input->noseq;
  } else {
    input->confirm = input->noseq;
  }
  
  /* calculate background base distribution */
  if (input->bfp) { /* reading from background probability file */
    if (input->bgmkv < 0)
      input->bgmkv = 3;
    readBgFreq(input->bfp, bglog, &(input->bgmkv));
  } else {
    if (input->bsp) /* reading from background sequence file */
      readSeqFile(input->bsp, NO, NULL, &seqCt, input);

    /* keep this line here in case user wants to use input as bg */
    calcSeqLog(bglog, &seqCt, &(input->bgmkv));
    if (!input->back) {
      printf("Computing background model, which will be written to file bgfreq\n");
      fprintf(input->ofp, "Background model printed in file bgfreq.\n");
    }
    if (!(fp = fopen("bgfreq", "w"))) 
      errorExit("Can't open bgfreq file.");
    printBgFreq(fp, bglog, input->bgmkv);
    fclose(fp);
  }
  
  /* find the min match for alignment */
  temp = 0;
  for (i = 0; i < ALPHASIZE; i++) {
    temp += pow(exp(bglog->one[i]), 2);
  }
  input->seedmw = input->w * temp + SIGNIFICANCE * 
    pow(input->w * temp * (1 - temp), 0.5);
  input->mw = (int) ceil(input->seedmw);
}

/******************************
 *
 * Func: readSeqFile
 * Allocate space for sequences, put input sequences in a linked list,
 * count the based distribution and put it in seqCount
 *
 ******************************/
void readSeqFile(FILE *fp, int isSeq, struct sequence **seqs, struct 
		 seqCount *seqct, struct inputParam *input) {
  char name[MAXLINE], seq[MAXLINE], line[MAXLINE];
  int len, flag = 0;
  struct sequence *tail = isSeq ? *seqs : NULL;

  memset(seqct, 0, sizeof(struct seqCount));
  while (fgets(line, MAXLINE, fp)) {
    len = mystrlen(line);
    if (len == 0) 
      continue;
    if (line[0] == '>') {
      if(flag != 0)
	createSeq(isSeq, seqs, seqct, name, seq, input, &tail);
      else
	flag = 1;
      strcpy(name, line);
      seq[0] = ENDLINE;
    } else
      strcat(seq, line);
  }
  if (seq[0] != ENDLINE)
    createSeq(isSeq, seqs, seqct, name, seq, input, &tail);
  fclose(fp);
}

/******************************
 *
 * Func: createSeq
 * Allocate space for one sequence, put it in a linked list, count the
 * based distribution and put it in seqCount. Each sequence added at the
 * tail of the linked list.
 *
 ******************************/
void createSeq(int isSeq, struct sequence **seqs, struct seqCount*seqct,
	       char *name, char *seq, struct inputParam *input, struct
	       sequence **tail) {
  struct sequence *s = NULL;
  int seqint[MAXLINE], len;
  
  /* if the sequence is a input sequence, allocate spaces */
  if (isSeq) {
    s = (struct sequence *) malloc (sizeof (struct sequence));
    len = strlen(name);
    s->name = (char *) malloc(len + 1);
    strcpy(s->name, name);
    s->len = len = strlen(seq); 
    if (input->noseq < input->confirm || !input->confirm) {
      s->tried1 = (int *) calloc(len, sizeof(int));
      s->tried2 = (int *) calloc(len, sizeof(int));
    } else {
      s->tried1 = s->tried2 = NULL;
    }
    input->noseq++;
    s->seqi = (int *) malloc (len * sizeof(int));
    s->bgs = (double *) malloc (len * sizeof(double));
    seqatoi(seq, s->seqi, len, FORWARD); 
    countSeq(s->seqi, len, seqct); 
    s->rseqi = (int *) malloc (len * sizeof(int));
    s->rbgs = (double *) malloc (len * sizeof(double));
    seqatoi(seq, s->rseqi, len, BACKWARD); 
    countSeq(s->rseqi, len, seqct);
    s->next = NULL;
    if (*seqs) {
      (*tail)->next = s;
      (*tail) = s;
    } else
      *seqs = *tail = s;
  } else {
    
    /* if a background sequence, just count its A,T,G,C. */
    len = strlen(seq);       
    seqatoi(seq, seqint, len, FORWARD); 
    countSeq(seqint, len, seqct); 
    seqatoi(seq, seqint, len, BACKWARD); 
    countSeq(seqint, len, seqct);
  }
}

/******************************
 *
 * Func: mystrlen
 * Return the input sequence name/seq length 
 *
 ******************************/
int mystrlen(char *string) {
  int len = strlen(string);
  char *c;

  c = &(string[len-1]);

  /* everything outside SPACE and TILTA '~' does not have a viewable
   * character representation of it */
  while (*c < SPACE || *c > TILTA) {
    *c = ENDLINE;
    c--;
  }
  return (c - string + 1);
}

/******************************
 *
 * Func: seqatoi
 * Convert sequence from a string of A, G, C, T to an array of 0, 1, 2,
 * 3. Need to rewrite if the alphabet is amino acid. 
 *
 ******************************/
void seqatoi(char *line, int *seqint, int len, int forward) {
  int i;

  for (i = 0; i < len; i++) {
    switch (line[i]) {
    case 'a':;
    case 'A': if (forward == FORWARD) seqint[i] = A; 
    else seqint[len-i-1] = T; break;
    case 't':;
    case 'T': if (forward == FORWARD) seqint[i] = T; 
    else seqint[len-i-1] = A; break;
    case 'g':;
    case 'G': if (forward == FORWARD) seqint[i] = G; 
    else seqint[len-i-1] = C; break;
    case 'c':;
    case 'C': if (forward == FORWARD) seqint[i] = C; 
    else seqint[len-i-1] = G; break;
    default: if (forward == FORWARD) seqint[i] = N;
    else seqint[len-i-1] = N;
    }
  }
}

/******************************
 *
 * Func: countSeq
 * Count the base distribution of a sequence in single, double, triplet,
 * and quadruplet bases
 * 
 * 2/8/02: added condition to deal with cases when seqlen < 3
 ******************************/
void countSeq(int *seqint, int len, struct seqCount *seqct) {
  int i; 
  
  for (i = 0; i < 3 && i < len; i++)
    seqct->one[seqint[i]]++;
  for (i = 1; i < 3 && i < len; i++)
   seqct->two[seqint[i-1]][seqint[i]]++;
  for (i = 2; i < 3 && i < len; i++)
    seqct->three[seqint[i-2]][seqint[i-1]][seqint[i]]++;
  for (i = 3; i < len; i++) {
    seqct->one[seqint[i]]++;
    seqct->two[seqint[i-1]][seqint[i]]++;
    seqct->three[seqint[i-2]][seqint[i-1]][seqint[i]]++;
    seqct->four[seqint[i-3]][seqint[i-2]][seqint[i-1]][seqint[i]]++;
  }
}

/******************************
 *
 * Func: calcSeqLog
 * Calculate the log of distribution from seqct and decides the Markov
 * dependency order
 *  
 ******************************/
void calcSeqLog(struct freq *seqlog, struct seqCount *seqct, int *mkv) {
  int i, j, k, l;
  int sum1, sum2, sum3, sum4;

  /* figure out the Markov order based on total sequence size */
  sum1 = 0;
  for (i = 0; i < ALPHASIZE; i++)
    sum1 += seqct->one[i];
  if (*mkv < 0) {
    if (sum1 > MKV3)
      *mkv = 3;
    else if (sum1 > MKV2)
      *mkv = 2;
    else if (sum1 > MKV1)
      *mkv = 1;
    else
      *mkv = 0;
  }

  /* calculate the log probability of bases */
  for (i = 0; i < ALPHASIZE; i++) {
    seqlog->one[i] = log((seqct->one[i]+VERYSMALL) / sum1);
    if (*mkv == 0) continue;
    sum2 = 0;
    for (j = 0; j < ALPHASIZE; j++)
      sum2 += seqct->two[i][j];
    for (j = 0; j < ALPHASIZE; j++) {
      seqlog->two[i][j] = log((seqct->two[i][j]+VERYSMALL) / sum2);
      if (*mkv == 1) continue;
      sum3 = 0;
      for (k = 0; k < ALPHASIZE; k++)
	sum3 += seqct->three[i][j][k];
      for (k = 0; k < ALPHASIZE; k++) {
	seqlog->three[i][j][k] = log((seqct->three[i][j][k]+VERYSMALL) / sum3);
	if (*mkv == 2) continue;
	sum4 = 0;
	for (l = 0; l < ALPHASIZE; l++)
	  sum4 += seqct->four[i][j][k][l];
	for (l = 0; l < ALPHASIZE; l++)
	  seqlog->four[i][j][k][l] = log((seqct->four[i][j][k][l]+VERYSMALL) /
					 sum4);
      }
    }
  }
}

/******************************
 *
 * Func: calcSeqCDF
 * Calculate sequence distribution in cumulative density function
 * format. This is useful in regerating sequences using Monte Carlo. 
 *
 * Checked 10/12/00
 ******************************/
void calcSeqCDF(struct freq *seqlog, struct freq *seqcdf, int mkv) {
  int i, j, k, l;
  double sum1, sum2, sum3, sum4;

  sum1 = 0;
  for (i = 0; i < ALPHASIZE; i++) {
    seqcdf->one[i] = sum1 += exp(seqlog->one[i]);
    if (mkv == 0) continue;
    sum2 = 0;
    for (j = 0; j < ALPHASIZE; j++) {
      seqcdf->two[i][j] = sum2 += exp(seqlog->two[i][j]);
      if (mkv == 1) continue;
      sum3 = 0;
      for (k = 0; k < ALPHASIZE; k++) {
	seqcdf->three[i][j][k] = sum3 += exp(seqlog->three[i][j][k]);
	if (mkv == 2) continue;
	sum4 = 0;
	for (l = 0; l < ALPHASIZE; l++) {
	  seqcdf->four[i][j][k][l] = sum4 += exp(seqlog->four[i][j][k][l]);
	}
      }
    }
  }
}

/******************************
 *
 * Func: readBgFreq
 * Read in the background distribution directly from a probability file
 *
 ******************************/
void readBgFreq(FILE *fp, struct freq *bglog, int *mkv) {
  int i, j, k;

  if (!readBgFreqLine(fp, bglog->one))
    errorExit("Background frequency file is empty.");
  if (*mkv > 0) {
    for (i = 0; i < ALPHASIZE; i++) {
      if (!readBgFreqLine(fp, bglog->two[i])) {
	*mkv = 0;
	break;
      }
    }
  }
  if (*mkv > 1) {
    for (i = 0; i < ALPHASIZE; i++) {
      for (j = 0; j < ALPHASIZE; j++) {
	if (!readBgFreqLine(fp, bglog->three[i][j])) {
	  *mkv = 1;
	  break;
	}
      }
    }
  }
  if (*mkv > 2) {
    for (i = 0; i < ALPHASIZE; i++) {
      for (j = 0; j < ALPHASIZE; j++) {
	for (k = 0; k < ALPHASIZE; k++) {
	  if (!readBgFreqLine(fp, bglog->four[i][j][k])) {
	    *mkv = 2;
	    break;
	  }
	}
      }
    }
  }
}

/******************************
 *
 * Func: readBgFreqLine
 * Read in one line of the background distribution from a background
 * frequency file. The four numbers are the log ratio of A, C, G, T in a
 * particular Markov order.
 *
 ******************************/
int readBgFreqLine(FILE *fp, double *bgfreq) {
  char line[REGLINE];

  if (!fgets(line, REGLINE, fp)) 
    return 0;
  if (!sscanf(line, "%lf %lf %lf %lf\n", &bgfreq[0], &bgfreq[1], 
	      &bgfreq[2], &bgfreq[3]))
    errorExit("Background frequency file doesn't have right format.");
  return 1;
}

/******************************
 *
 * Func: printBgFreq
 *  
 ******************************/
void printBgFreq(FILE *fp, struct freq *bg, int mkv) {
  int i, j, k, l;
  
  for (i = 0; i < ALPHASIZE; i++) {
    fprintf(fp, "%f ", bg->one[i]);
  }
  fprintf(fp, "\n");
  if (mkv > 0) {
    for (i = 0; i < ALPHASIZE; i++) {
      for (j = 0; j < ALPHASIZE; j++) {
	fprintf(fp, "%f ", bg->two[i][j]);
      }
      fprintf(fp, "\n");
    }
    if (mkv > 1) {
      for (i = 0; i < ALPHASIZE; i++) {
	for (j = 0; j < ALPHASIZE; j++) {
	  for (k = 0; k < ALPHASIZE; k++) {
	    fprintf(fp, "%f ", bg->three[i][j][k]);
	  }
	  fprintf(fp, "\n");
	}
      }
      if (mkv > 2) {
	for (i = 0; i < ALPHASIZE; i++) {
	  for (j = 0; j < ALPHASIZE; j++) {
	    for (k = 0; k < ALPHASIZE; k++) {
	      for (l = 0; l < ALPHASIZE; l++) {
		fprintf(fp, "%f ", bg->four[i][j][k][l]); 
	      }
	      fprintf(fp, "\n");
	    }
	  }
	}
      }
    }
  }
}

/******************************
 *
 * Func: initMotif
 * Setup dynamic space for motif structure. Double array pb** are useful
 * during update for storing segment scores and sampling.
 *
 ******************************/
void initMotif(struct motif **mtf, int w) { 
  int i;
  
  (*mtf) = (struct motif *) calloc(1, sizeof (struct motif));
  (*mtf)->ahead = (*mtf)->atail = NULL;
  (*mtf)->score = -1000;
  (*mtf)->conint = (int *) malloc(w * sizeof(int));
  (*mtf)->rconint = (int *) malloc(w * sizeof(int));
  (*mtf)->constr = (char *) malloc((w + 1) * sizeof(char));
  (*mtf)->rconstr = (char *) malloc((w + 1) * sizeof(char));
  (*mtf)->degstr = (char *) malloc((w + 1) * sizeof(char));
  (*mtf)->rdegstr = (char *) malloc((w + 1) * sizeof(char));
  (*mtf)->constr[w] = (*mtf)->rconstr[w] = 
    (*mtf)->degstr[w] = (*mtf)->rdegstr[w] = ENDLINE;
  (*mtf)->blk = (struct mtfPt **)malloc(w * sizeof(struct mtfPt *));
  for (i = 0; i < w; i++)
    (*mtf)->blk[i] = (struct mtfPt *)malloc(4 * sizeof(struct mtfPt));
}

/******************************
 *
 * Func: findMotifTh
 * Find the motif distribution, thus the motif score threshold. It spawns
 * off child processes to search for motifs from regenerated sequences,
 * and get the highest score from each set of regenerated sequences. This
 * will then be fitted to a normal distribution. 
 *
 ******************************/
void findMotifTh(struct inputParam *input, struct sequence *seqs, struct motif **motifs, 
		 struct motif **mtf, struct freq *bglog, struct freq *seqcdf) {
  int i, j, k, pd[2], data = input->MonteCarlo;
  double stats[data * sizeof(double)], score;
  
  if (pipe(pd) == -1)
    errorExit("Can't generate pipes.");
  for (i = 0; i < input->child; i++) {
    if (fork() == 0) {
      
      /* child processes randomize (regenerate) input, search for motifs,
       * and report to parent the highest scoring motif of each generated
       * date set. */
      initRand();
      for (j = 0; j < input->norand; j++) {
	for (k = 0; k < input->scan; k++)
	  clearMotif(motifs[k], input);
	clearMotif(*mtf, input);
	randomize(seqs, seqcdf, input->seqmkv);
	wholeBgScore(seqs, bglog, input);
	findCandidateMtf(input, seqs, motifs, mtf, bglog);
	scanning(seqs, input, motifs, bglog);
	for (k = 0; k < input->iterate; k++)
	  refine(input, motifs, bglog);
	score = motifs[0]->score;
	for (k = 1; k < input->scan; k++) {
	  if (motifs[k]->score > score)
	    score = motifs[k]->score;
	}
	if (write(pd[1], &score, sizeof(double)) == -1)
	  errorExit("Failed writing to pipe.");
      }
      exit(0);
    } 
  }
      
  /* parent wait for child report, and store each data in array stats */
  for (i = 0; i < data; i++) {
    if (read(pd[0], &score, sizeof(double)) == -1)
      errorExit("Failed reading from pipe.");
    stats[i] = score;
    if (!input->back)
      printf("Monte Carlo # %d\t%.3f\n", i, score);
  }

  /* calculate the mtfmean, mtfstd from stats */
  calcStats(stats, input);
}

/******************************
 *
 * Func: initRand
 * Initialize the random generator seed to be a non-zero long
 *
 ******************************/
void initRand() {
  struct timeval tv;
  long temp = 0;

  while (1) {
    gettimeofday(&tv, NULL);
    if (tv.tv_usec != temp) {
      idum = tv.tv_usec;
      return;
    }
  }
}

/******************************
 *
 * Func: randomize
 * Regenerate input sequence using input sequence distribution
 *
 * Checked 3/20/01
 ******************************/
void randomize(struct sequence *seqs, struct freq *seqcdf, int mkv) {
  struct sequence *s; 
  int i, len, size;
  
  for (s = seqs; s; s = s->next) {
    len = s->len;

    /* clear out all tried flags */
    if (s->tried1) {
      size = len * sizeof(int);
      memset(s->tried1, 0, size);
      memset(s->tried2, 0, size);
    }
    switch(mkv) {
    case 0: 
      genRandSeq0(s->seqi, len, seqcdf); break;
    case 1: 
      genRandSeq1(s->seqi, len, seqcdf); break;
    case 2: 
      genRandSeq2(s->seqi, len, seqcdf); break;
    case 3: 
      genRandSeq3(s->seqi, len, seqcdf); break;
    }
    for (i = 0; i < len; i++)
      s->rseqi[i] = ALPHASIZE - 1 - s->seqi[len-i-1];
  }
}

/******************************
 *
 * Func: genRandSeq0
 * Generate a sequence of len using 0 order Markov
 *
 * Checked 3/20/01
 ******************************/
void genRandSeq0(int *seq, int len, struct freq *schar) {
  int i, j;
  double random;

  for (j = 0; j < len; j++) {
    random = drand();
    for (i = 0; i < ALPHASIZE; i++) {
      if (random < schar->one[i]) {
	seq[j] = i;
	break;
      }
    }
  }
}

/******************************
 *
 * Func: genRandSeq1
 * Generate a sequence of len using 1st order Markov
 *
 * Checked 3/20/01
 ******************************/
void genRandSeq1(int *seq, int len, struct freq *schar) {
  int i, j, p1;
  double random;

  random = drand();
  for (i = 0; i < ALPHASIZE; i++) {
    if (random < schar->one[i]) {
      p1 = seq[0] = i;
      break;
    }
  }
  for (j = 1; j < len; j++) {
    random = drand();
    for (i = 0; i < ALPHASIZE; i++) {
      if (random < schar->two[p1][i]) {
	p1 = seq[j] = i ;
	break;
      }
    }
  }
}

/******************************
 *
 * Func: genRandSeq2
 * Generate a sequence of len using 2nd order Markov
 *
 * Checked 3/20/01
 ******************************/
void genRandSeq2(int *seq, int len, struct freq *schar) {
  int i, j, p1, p2;
  double random;

  random = drand();
  for (i = 0; i < ALPHASIZE; i++) {
    if (random < schar->one[i]) {
      p1 = seq[0] = i;
      break;
    }
  }
  random = drand();
  for (i = 0; i < ALPHASIZE; i++) {
    if (random < schar->two[p1][i]) {
      p2 = seq[1] = i;
      break;
    }
  }
  for (j = 2; j < len; j++) {
    random = drand();
    for (i = 0; i < ALPHASIZE; i++) {
      if (random < schar->three[p1][p2][i]) {
	p1 = p2;
	p2 = seq[j] = i ;
	break;
      }
    }
  }
}

/******************************
 *
 * Func: genRandSeq3
 * Generate a sequence of len using 3rd order Markov
 *
 * Checked 3/20/01
 ******************************/
void genRandSeq3(int *seq, int len, struct freq *schar) {
  int i, j, p1, p2, p3;
  double random;

  random = drand();
  for (i = 0; i < ALPHASIZE; i++) {
    if (random < schar->one[i]) {
      p1 = seq[0] = i;
      break;
    }
  }
  random = drand();
  for (i = 0; i < ALPHASIZE; i++) {
    if (random < schar->two[p1][i]) {
      p2 = seq[1] = i;
      break;
    }
  }
  random = drand();
  for (i = 0; i < ALPHASIZE; i++) {
    if (random < schar->three[p1][p2][i]) {
      p3 = seq[2] = i;
      break;
    }
  }
  for (j = 3; j < len; j++) {
    random = drand();
    for (i = 0; i < ALPHASIZE; i++) {
      if (random < schar->four[p1][p2][p3][i]) {
	p1 = p2;
	p2 = p3;
	p3 = seq[j] = i;
	break;
      }
    }
  }
}

/******************************
 *
 * Func: drand
 * Generate a double random number between 0 and 1. This method comes from
 * p279 of "Numerical Recipes in C", 2nd Ed, with minor modifications. It
 * has a period of 2.1 x 10^9, returns uniform random deviate between 0.0
 * and 1.0. 
 *
 ******************************/
double drand() {
  long k;
  
  k = idum / IQ;
  idum = IA * (idum - k * IQ) - IR * k;
  if (idum < 0) idum += IM;
  return (AM * idum);
}

/******************************
 *
 * Func: wholeBgScore
 * Since the bgscore of each sequence is the same throughout, we just need to 
 * compute it once.
 *
 ******************************/
void wholeBgScore(struct sequence *seqs, struct freq *bglog, struct 
		  inputParam *input) {
  struct sequence *s;
  int w = input->w, end, mkv = input->bgmkv, i;
  double *bgs, *rbgs;
  int *seqi, *rseqi;

  for (s = seqs; s; s = s->next) {
    end = s->len - w + 1;
    bgs = s->bgs;
    rbgs = s->rbgs;
    seqi = s->seqi;
    rseqi = s->rseqi;
    
    for (i = 0; i < end; i++) {
      bgs[i] = seqBgScore(&(seqi[i]), i, w, mkv, bglog);
      rbgs[i] = seqBgScore(&(rseqi[i]), i, w, mkv, bglog);
    }
  }
}

/******************************
 *
 * Func: seqBgScore
 * Score the segment of wid long starting from *seqi by bg
 *
 ******************************/
double seqBgScore(int *ptr, int st, int wid, int mkv, struct freq *bg) {
  double score = 0;
  int i;

  switch(mkv) {
  case 0: for (i = 0; i < wid; i++) {
    if (*(ptr+i) == N)
      return 0;
    score += bg->one[*(ptr+i)]; 
  }
    break;
  case 1: 
    if (st > 0 && *(ptr-1) == N) 
      st = 0;
    for (i = 0; i < wid; i++)
      if (*(ptr+i) == N)
	return 0;
    if (st == 0) {
      score += bg->one[*ptr];
      for (i = 1; i < wid; i++) score += bg->two[*(ptr+i-1)][*(ptr+i)]; 
    } else {
      for (i = 0; i < wid; i++) score += bg->two[*(ptr+i-1)][*(ptr+i)];
    }
    break;
  case 2: 
    if (st > 1 && *(ptr-2) == N) 
      st = 1;
    if (st > 0 && *(ptr-1) == N) 
      st = 0;
    for (i = 0; i < wid; i++)
      if (*(ptr+i) == N)
	return 0;
    switch (st) {
    case 0: score += bg->one[*ptr] + bg->two[*ptr][*(ptr+1)];
      for (i = 2; i < wid; i++) 
	score += bg->three[*(ptr+i-2)][*(ptr+i-1)][*(ptr+i)];
      break;
    case 1: score += bg->two[*(ptr-1)][*ptr];
      for (i = 1; i < wid; i++) 
	score += bg->three[*(ptr+i-2)][*(ptr+i-1)][*(ptr+i)];
      break;
    default: 
      for (i = 0; i < wid; i++) 
	score += bg->three[*(ptr+i-2)][*(ptr+i-1)][*(ptr+i)];
      break;
    }
    break;
  case 3: 
    if (st > 2  && *(ptr-3) == N) 
      st = 2;
    if (st > 1  && *(ptr-2) == N) 
      st = 1;
    if (st > 0  && *(ptr-1) == N) 
      st = 0;
    for (i = 0; i < wid; i++)
      if (*(ptr+i) == N)
	return 0;
    switch (st) {
    case 0: score += bg->one[*ptr] + bg->two[*ptr][*(ptr+1)]
	      + bg->three[*ptr][*(ptr+1)][*(ptr+2)];
      for (i = 3; i < wid; i++) 
	score += bg->four[*(ptr+i-3)][*(ptr+i-2)][*(ptr+i-1)][*(ptr+i)];
      break;
    case 1: score += bg->two[*(ptr-1)][*ptr] 
	      + bg->three[*(ptr-1)][*ptr][*(ptr+1)];
      for (i = 2; i < wid; i++) 
	score += bg->four[*(ptr+i-3)][*(ptr+i-2)][*(ptr+i-1)][*(ptr+i)];
      break;
    case 2: score += bg->three[*(ptr-2)][*(ptr-1)][*ptr];
      for (i = 1; i < wid; i++)
	score += bg->four[*(ptr+i-3)][*(ptr+i-2)][*(ptr+i-1)][*(ptr+i)];
      break;
    default: 
      for (i = 0; i < wid; i++)
	score += bg->four[*(ptr+i-3)][*(ptr+i-2)][*(ptr+i-1)][*(ptr+i)];
      break;
    }
  }

  return score;
}

/******************************
 *
 * Func: findCandidateMtf
 * Check all existings w-mers from the top sequences and find their matching 
 * strings to initalize the motif. Keep the top motif seeds. 
 * 
 ******************************/
void findCandidateMtf(struct inputParam *input, struct sequence *seqs, struct
		      motif **motifs, struct motif **mtf, struct freq *bglog) {
  int i, j, k, len, st1, st2, basect[4],  w = input->w, diffbase, 
    cutoff = (int) floor(0.8 * w), masked;
  struct sequence *s;

  /* mask out all words containing N */
  for (i = 0, s = seqs; i < input->confirm; i++, s = s->next) {
    len = s->len - w + 1;
    for (j = 0; j < len; j++) {
      if (s->seqi[j] == N) {
	st1 = j - w + 1 >= 0 ? j - w + 1 : 0;
	st2 = j + w < len ? j + w : len-1;
	for (k = st1; k <= st2; k++) {
	  s->tried1[k] = HASN;
	  s->tried2[len-k-1] = HASN;
	}
      } 
    }
  }

  /* mask out all repeats */ 
  for (i = 0, s = seqs; i < input->confirm; i++, s = s->next) {
    len = s->len - w + 1;
    for (j = 0; j < len; j++) {
      masked = 0;
      if (!s->tried1[j]) {
	basect[A] = basect[C] = basect[G] = basect[T] = diffbase = 0;
	for (k = 0; k < w; k++) {
	  basect[s->seqi[j+k]]++;
	}

	for (k = 0; k < ALPHASIZE; k++) {
	  if (basect[k] > 0) {

	    /* any single base occupying > 80% is bad */
	    if (basect[k] > cutoff) {
	      s->tried1[j] = MASK;
	      s->tried2[len-j-1] = MASK;
	      masked = 1;
	    }
	    diffbase++;
	  }
	}

	/* two kind of bases, make sure the two bases are not alternating */
        if (diffbase == 2 && !masked) {
	  masked = 1;
	  for (k = 2; k < w; k += 2) {
	    if (s->seqi[j] != s->seqi[j+k]) {
	      masked = 0;
	      break;
	    }
	  }
	  if (masked) {
	    for (k = 3; k < w; k+= 2) {
	      if (s->seqi[j+1] != s->seqi[j+k]) {
		masked = 0;
		break;
	      }
	    }
	  }
	  if (masked) {
	    s->tried1[j] = MASK;
	    s->tried2[len-j-1] = MASK;
	  }
	}
      }
    }
  }

  if (input->seed) {
    tryOneSeed(mtf, input, seqs, motifs, bglog);
  } else {
    /* now can enumerate seeds to look for candidate motifs */
    for (i = 0, s = seqs; i < input->top; i++, s = s->next) {
      if (!input->back)
	printf("Seed with Sequence #%d\n", (i+1));
      enumerateSeed(mtf, s->seqi, seqs, input, s->tried1, 
		    s->len - w + 1, motifs, bglog);
      enumerateSeed(mtf, s->rseqi, seqs, input, s->tried2, 
		    s->len - w + 1, motifs, bglog);
    }
  }
}

/******************************
 *
 * Func: tryOneSeed
 * Check every word in top sequences for matches to the seed.
 * 
 ******************************/
void tryOneSeed(struct motif **mtf, struct inputParam *input, 
		struct sequence *seqs, struct motif **motifs,
		struct freq *bglog) {
  int j, k, *seq1, *seq2, *try1, *try2, end;
  int w = input->w;
  char *seed = input->seed;
  double *bgs1, *bgs2, mis = w - input->mw;
  struct sequence *s;

  /*
  printf("Miss %.2f\n", mis);
  */
  clearMotif(*mtf, input);
  for (j = 0, s = seqs; j < input->top; j++, s = s->next) {
    end = s->len - w + 1;
    seq1 = s->seqi;
    seq2 = s->rseqi;
    bgs1 = s->bgs;
    bgs2 = s->rbgs;
    try1 = s->tried1;
    try2 = s->tried2;
    for (k = 0; k < end; k++) {
      if (try1[k] < HASN) {
	if (seedMatch(seed, &(seq1[k]), w, mis)) {
	  addAlign(s, k+1, &(seq1[k]), bgs1[k], *mtf);
	  addSegment(*mtf, &(seq1[k]), w, bgs1[k], ADD);
	}
      }
      if (try2[k] < HASN) {
	if (seedMatch(seed, &(seq2[k]), w, mis)) {
	  addAlign(s, -k-1, &(seq2[k]), bgs2[k], *mtf);
	  addSegment(*mtf, &(seq2[k]), w, bgs2[k], ADD);
	}
      }
    }
  }
  
  
  if (motifInfo(*mtf, input, bglog) > 0) {
    insertSeed(motifs, mtf, input);
  } else {
    fprintf(stderr, "Seed motif score %.3f\n", motifInfo(*mtf, input, bglog));
    errorExit("Seed motif is not a good motif");
  }
}

/******************************
 *
 * Func: enumerateSeed
 * Check all existings w-mers from the top sequences and find their matching 
 * strings to initalize the motif. Keep the top motif seeds. 
 * This func is here to let findCandidateMtf search for w-mers in forward
 * and backward of each sequences (saves some code).
 * 
 ******************************/
void enumerateSeed(struct motif **mtf, int *seq, struct sequence *seqs, 
		   struct inputParam *input, int *try, int end, struct motif 
		   **motifs, struct freq *bglog) {
  int i, j, k, *seq1, *seq2, *try1, *try2, temp, end2;
  int w = input->w, mis = w - input->mw;
  double *bgs1, *bgs2;
  struct sequence *s;

  for (i = 0; i < end; i++) {
    if (!try[i]) {
      clearMotif(*mtf, input);
      for (j = 0, s = seqs; j < input->top; j++, s = s->next) {
	end2 = s->len - w + 1;
	seq1 = s->seqi;
	seq2 = s->rseqi;
	bgs1 = s->bgs;
	bgs2 = s->rbgs;
	try1 = s->tried1;
	try2 = s->tried2;
	for (k = 0; k < end2; k++) {
	  if (try1[k] < HASN) {
	    temp = match(&(seq[i]), &(seq1[k]), w, mis);
	    if (temp) {
	      addAlign(s, k+1, &(seq1[k]), bgs1[k], *mtf);
	      addSegment(*mtf, &(seq1[k]), w, bgs1[k], ADD);
	      if (temp == -1)
		try1[k] = try2[s->len-w-k] = 1;
	    }
	  }
	  if (try2[k] < HASN) {
	    temp = match(&(seq[i]), &(seq2[k]), w, mis);
	    if (temp) {
	      addAlign(s, -k-1, &(seq2[k]), bgs2[k], *mtf);
	      addSegment(*mtf, &(seq2[k]), w, bgs2[k], ADD);
	      if (temp == -1)
		try2[k] = try1[s->len-w-k] = 1;
	    }
	  }
	}
      }

      /* compare this motif score with all other considered motifs *
       * and update motifs if mtf score is good enough, motifs score *
       * increases with motif[i] */
      if (motifInfo(*mtf, input, bglog) > motifs[0]->score) {
	insertSeed(motifs, mtf, input);
	/*
	for (j = 0; j < input->scan; j++) {
	  printf("Motif %d\t", j+1);
	  for (k = 0; k < w; k++)
	    printf("%d", motifs[j]->conint[k]);
	  printf("\t");
	  for (k = 0; k < w; k++)
	    printf("%d", motifs[j]->rconint[k]);
	  printf("\t%.3f\n", motifs[j]->score);
	}
	*/
      }
    }
  }
}

/******************************
 *
 * Func: clearMotif
 * Clear the motif matrix and segments and alignment. 
 * 
 ******************************/
void clearMotif(struct motif *mtf, struct inputParam *input) {
  int i;
  struct alignment *a1, *a2;

  mtf->score = mtf->bginfo = mtf->zscore = mtf->avgInfo = 0;
  mtf->segment = 0;
  for (i = 0; i < input->w; i++)
    memset(mtf->blk[i], 0, 4 * sizeof(struct mtfPt));
  for (a1 = mtf->ahead; a1; a1 = a2) {
    a2 = a1->next;
    free(a1);
  }
  mtf->ahead = mtf->atail = NULL;
}

/******************************
 *
 * Func: seedMatch
 * If y mismatches x by < mis in wid, return 1, 
 * otherwise return 0
 *
 ******************************/
int seedMatch(char *x, int *y, int wid, double mis) {
  int i;
  double score = 0;

  for (i = 0; i < wid; i++) {
    if (x[i] == 'A' && y[i] == A)
      score += 1;
    if (x[i] == 'C' && y[i] == C)
      score += 1;
    if (x[i] == 'G' && y[i] == G)
      score += 1;
    if (x[i] == 'T' && y[i] == T)
      score += 1;
    else if (x[i] == 'M' && (y[i] == A || y[i] == C))
      score += 0.5;
    else if (x[i] == 'R' && (y[i] == A || y[i] == G))
      score += 0.5;
    else if (x[i] == 'W' && (y[i] == A || y[i] == T))
      score += 0.5;
    else if (x[i] == 'S' && (y[i] == C || y[i] == G))
      score += 0.5;
    else if (x[i] == 'Y' && (y[i] == C || y[i] == T))
      score += 0.5;
    else if (x[i] == 'K' && (y[i] == G || y[i] == T))
      score += 0.5;
    else if (x[i] == 'B' && (y[i] == C || y[i] == G || y[i] == T))
      score += 1.0 / 3;
    else if (x[i] == 'D' && (y[i] == A || y[i] == G || y[i] == T))
      score += 1.0 / 3;
    else if (x[i] == 'H' && (y[i] == A || y[i] == C || y[i] == T))
      score += 1.0 / 3;
    else if (x[i] == 'V' && (y[i] == A || y[i] == C || y[i] == G))
      score += 1.0 / 3;
    else if (x[i] == 'N')
      score += 1.0 / 4;
  }
  if (wid - score <= mis) 
    return 1;
  else
    return 0;
}

/******************************
 *
 * Func: match
 * If y mismatches x by < mis in wid, return 1, 
 * otherwise return 0

 * if mistach == 0, then return -1;
 *
 ******************************/
int match(int *x, int *y, int wid, int mis) {
  int misCt = 0, i;

  for (i = 0; i < wid; i++) {
    if (x[i] != y[i]) {
      misCt++;
      if (misCt > mis)
	return 0;
    }
  }
  if (misCt)
    return 1;
  else
    return -1;
}

/******************************
 *
 * Func: addAlign
 * Add alignment to motif atail.
 *
 ******************************/
void addAlign(struct sequence *s, int pos, int *st, double bginfo, struct 
	      motif *mtf) {
  struct alignment *a = (struct alignment *) malloc 
    (sizeof (struct alignment));
  
  a->s = s;
  a->position = pos;
  a->st = st;
  a->bginfo = bginfo;
  a->next = NULL;
  if (mtf->ahead) 
    mtf->atail = mtf->atail->next = a;
  else
    mtf->ahead = mtf->atail = a;
}

/******************************
 *
 * Func: addSegment
 * Add or subtract one segment sequence to a motif block. start points to
 * the starting position of the segment, and w is the the width. add tells
 * whether adding or subtracting the segment.
 *
 ******************************/
void addSegment(struct motif *mtf, int *start, int w, double bginfo, int add) {
  int i;
  
  mtf->bginfo += add * bginfo;
  mtf->segment += add;
  for (i = 0; i < w; i++, start++)
    mtf->blk[i][*start].ct += add;
}

/******************************
 *
 * Func: motifInfo
 * Motif info is the log(segment) * (average relative entropy per column)
 *
 ******************************/
double motifInfo(struct motif *mtf, struct inputParam *input, struct freq *bglog) {
  double info;

  if (mtf->segment == 0) {
    mtf->score = 0;
    return 0;
  }
  
  reCalcScount(mtf->segment, input->scount, bglog);
  info = blkMtfScore(mtf->blk, input, input->w, mtf->segment);
  if (input->expect) {
    /*
    printf("Info %.2f, bginfo %.2f, penalty %.2f\n", info, mtf->bginfo,
	   mtf->segment * log(input->expect));
    */
    mtf->score = (info - mtf->bginfo - mtf->segment * log(input->expect)) 
      / input->w;
  } else {
    mtf->score = log(mtf->segment) * (info - mtf->bginfo) 
      / mtf->segment / input->w;
  }
  return mtf->score;
}

/******************************
 *
 * Func: reCalcScount
 * Let total number of pseudo counts proportional to the number of segments
 * 
 ******************************/
void reCalcScount(int segment, double *sct, struct freq *bglog) {
  int i;

  sct[ALPHASIZE] = sqrt(segment);
  for (i = 0; i < ALPHASIZE; i++) {
    sct[i] = sct[ALPHASIZE] * exp(bglog->one[i]);
  }
}

/******************************
 *
 * Func: blkMtfScore
 * Calculates motif matrix entropy
 * 
 ******************************/
double blkMtfScore(struct mtfPt **blk, struct inputParam *input, int wid, int
		    segment) {
  int i, j; 
  double info = 0;

  calcBlkLog(blk, segment, input->w, input->scount);
  for (i = 0; i < input->w; i++) {
    for (j = 0; j < ALPHASIZE; j++) {
      info += blk[i][j].ct * blk[i][j].log;
    }
  }
  return info;
}

/******************************
 *
 * Func: calcBlkLog
 * Convert motif distribution from count to log so scoring is done faster
 *
 ******************************/
void calcBlkLog(struct mtfPt **blk, int segment, int wid, double *sct) {
  int i, j; 

  for (i = 0; i < wid; i++) {
    for (j = 0; j < ALPHASIZE; j++)
      blk[i][j].log = log((blk[i][j].ct + sct[j]) / (segment + sct[ALPHASIZE]));
  }
}

/******************************
 *
 * Func: insertSeed
 * mtf is better than motifs[0] (motifs score increase with index)
 * try to update the motifs array with this new mtf, also needs to 
 * eliminate duplicate motifs
 *
 ******************************/
void insertSeed(struct motif **motifs, struct motif **mtf, struct 
		inputParam *input) {
  struct motif *mtfptr;
  int i, j, k, scan = input->scan, w = input->w, mw = input->mw;
  double info = (*mtf)->score;
  
  getConInt(w, mw, *mtf);  
  /* if the motif is similar */
  for (i = scan - 1; i >= 0; i--) {
    if (motifsSimilar(motifs[i], *mtf, w, mw)) {
      if (motifs[i]->score < info) {
	
	/* This motif is similar to motifs[i] but with higher score, 
	   first replace motifs[i]'s slot */
	mtfptr = motifs[i];
	motifs[i] = *mtf;
	*mtf = mtfptr;
	
	/* Before bubble up, check all motifs below it and make sure 
	   they are not similar to it. This situation arises when A and 
	   B (higher score) in list are not similar enough, but X is 
	   close to A and B with even higher score. After X replaces
	   B, A and X are similar and both in the list. So need to 
	   eliminate A as well.
	*/
	for (j = i-1; j >= 0; j--) {
	  if (motifsSimilar(motifs[j], motifs[i], w, mw)) {
	    /*
	    printf("Node %d eliminates node %d\n", i, j);
	    */
	    for (k = j; k > 0; k--) {
	      mtfptr = motifs[k];
	      motifs[k] = motifs[k-1];
	      motifs[k-1] = mtfptr;
	    }
	    clearMotif(motifs[0], input);
	  }
	}
	  
	/* Now use bubble sort to bring this motif up */
	for (j = i + 1; j < scan; j++) {
	  if (motifs[j]->score < info) {
	    mtfptr = motifs[j];
	    motifs[j] = motifs[j-1];
	    motifs[j-1] = mtfptr;
	  } else 
	    return;
	}
      }
      return;
    }
  }

  /* If reach here, no motif is similar, so insert as usual */
  mtfptr = motifs[0];
  motifs[0] = *mtf;
  *mtf = mtfptr;
  for (i = 1; i < scan; i++) {
    if (motifs[i]->score < info) {
      mtfptr = motifs[i];
      motifs[i] = motifs[i-1];
      motifs[i-1] = mtfptr;
    } else
      return;
  }
}


/******************************
 *
 * Func: getConInt
 * Find integer representation of mtf consensus and r-consensus
 *
 ******************************/
void getConInt(int w, int mw, struct motif *mtf) {
  int i, j, *con = mtf->conint, *rcon = mtf->rconint;
  double maxs;
  
  mtf->avgInfo = 0;
  for (i = 0; i < w; i++) {
    maxs = mtf->blk[i][0].log; 
    con[i] = 0;
    rcon[w-i-1] = ALPHASIZE - 1;
    for (j = 1; j < ALPHASIZE; j++) {
      if (mtf->blk[i][j].log > maxs) {
	maxs = mtf->blk[i][j].log;
	con[i] = j;
	rcon[w-i-1] = ALPHASIZE - j - 1;
      }
    }
    mtf->avgInfo += exp(maxs);
  }
  mtf->avgInfo /= w;
}

/******************************
 *
 * Func: motifsSimilar
 * Test whether mtf1 and mtf2 have the same core
 * return 1 (similar) or 0 (not similar)
 *
 ******************************/
int motifsSimilar(struct motif *mtf1, struct motif *mtf2, int w, int mw) {
  int i, j, k, shift = w - mw, wid, match1, match2, 
    *con1 = mtf1->conint, *con2 = mtf2->conint, *rcon2 = mtf2->rconint;

  if (mtf1->segment == 0)
    return 0;
  
  /*
  for (i = 0; i < w; i++)
    printf("%d", con1[i]);
  printf("\t");
  for (i = 0; i < w; i++)
    printf("%d", rcon1[i]);
  printf("\t%.3f\t", mtf1->score);
  for (i = 0; i < w; i++)
    printf("%d", con2[i]);
  printf("\t");
  for (i = 0; i < w; i++)
    printf("%d", rcon2[i]);
  printf("\t%.3f\t", mtf2->score);
  */

  /* This has duplicates but at least doesn't miss anything */
  for (i = 0; i <= shift; i++) {
    for (j = 0; j <= shift; j++) {
      if (i && j)
	continue;
      wid = i > j ? w - i : w - j; 
      
      match1 = match2 = 0;
      for (k = 0; k < wid; k++) {
	if (con1[i+k] == con2[j+k])
	  match1++;
	if (con1[i+k] == rcon2[j+k])
	  match2++;
      }
      /*
      if (match1 >= mw)
	printf ("MatchedF%d%d\n", i, j);
      if (match2 >= mw)
	printf ("MatchedB%d%d\n", i, j);
      */
      if (match1 >= mw || match2 >= mw)
	return 1;
    }
  }
  /*
  printf("\n");
  */
  
  return 0;
}
 
/******************************
 *
 * Func: printCandidates
 * print out input->scan number of top motif candidates. 
 *
 ******************************/
void printCandidates(struct inputParam *input, struct motif **motifs) {
  int i, w = input->w;

  for (i = input->scan-1; i >= 0; i--) {
    if (motifs[i]->segment) {
      getConsensus(w, motifs[i]);
      fprintf(input->ofp, "%d\t%.2f\t%d\t%s\t%s\n",
	      input->scan-i, motifs[i]->score, motifs[i]->segment, 
	      motifs[i]->constr, motifs[i]->rconstr);
    }
  }
  fflush(input->ofp);
}

/******************************
 *
 * Func: getConsensus
 * Just match and consensus string from the matrix
 *
 ******************************/
void getConsensus(int w, struct motif *mtf) {
  double maxs;
  int i, j, maxp, sum;
  
  for (i = 0; i < w; i++) {
    maxs = -1000; 
    sum = 0;
    for (j = 0; j < ALPHASIZE; j++) {
      if (exp(mtf->blk[i][j].log) > 0.25)
	sum += (int)(rint(pow(2, j)));
      if (mtf->blk[i][j].log > maxs) {
	maxs = mtf->blk[i][j].log;
	maxp = j;
      }
    }
    
    /* max is for consensus */
    switch(maxp) {
    case A: mtf->constr[i] = 'A'; mtf->rconstr[w-i-1] = 'T'; break;
    case G: mtf->constr[i] = 'G'; mtf->rconstr[w-i-1] = 'C'; break;
    case C: mtf->constr[i] = 'C'; mtf->rconstr[w-i-1] = 'G'; break;
    case T: mtf->constr[i] = 'T'; mtf->rconstr[w-i-1] = 'A'; break;
    }

        /* sum is for degenerate */
    switch(sum) {
    case 0: mtf->degstr[i] = 'N'; mtf->rdegstr[w-i-1] = 'N'; break;
    case 1: mtf->degstr[i] = 'A'; mtf->rdegstr[w-i-1] = 'T'; break;
    case 2: mtf->degstr[i] = 'C'; mtf->rdegstr[w-i-1] = 'G'; break;
    case 3: mtf->degstr[i] = 'M'; mtf->rdegstr[w-i-1] = 'K'; break;   /* A, C */
    case 4: mtf->degstr[i] = 'G'; mtf->rdegstr[w-i-1] = 'C'; break;
    case 5: mtf->degstr[i] = 'R'; mtf->rdegstr[w-i-1] = 'Y'; break;   /* A, G */
    case 6: mtf->degstr[i] = 'S'; mtf->rdegstr[w-i-1] = 'S'; break;   /* G, C */
    case 7: mtf->degstr[i] = 'V'; mtf->rdegstr[w-i-1] = 'B'; break;   /* A, G, C */
    case 8: mtf->degstr[i] = 'T'; mtf->rdegstr[w-i-1] = 'A'; break;   /* T */
    case 9: mtf->degstr[i] = 'W'; mtf->rdegstr[w-i-1] = 'W'; break;   /* A, T */
    case 10: mtf->degstr[i] = 'Y'; mtf->rdegstr[w-i-1] = 'R'; break;  /* C, T */
    case 11: mtf->degstr[i] = 'H'; mtf->rdegstr[w-i-1] = 'D'; break;  /* A, T, C */
    case 12: mtf->degstr[i] = 'K'; mtf->rdegstr[w-i-1] = 'M'; break;  /* T, G */
    case 13: mtf->degstr[i] = 'D'; mtf->rdegstr[w-i-1] = 'H'; break;  /* A, T, G */
    case 14: mtf->degstr[i] = 'B'; mtf->rdegstr[w-i-1] = 'V'; break;  /* T, G, C */
    case 15: mtf->degstr[i] = 'N'; mtf->rdegstr[w-i-1] = 'N'; break;
    }
  }
}
  
/******************************
 *
 * Func: compactMotifs
 * Sort all the motifs, calculate zscore, and elimiate duplicates
 * After the call, all distinct motifs appear on bottom (low index)
 * sorted by motif score, and high index have cleared motifs (score = 0)
 *
 ******************************/
void compactMotifs(struct motif **motifs, struct inputParam *input) {
  int i, j, maxp;
  double maxs;
  struct motif *mtfptr;

  for (i = 0; i < input->scan; i++) {
    getConInt(input->w, input->mw, motifs[i]);
    if (motifs[i]->avgInfo < 0.65)
      clearMotif(motifs[i], input);
  }
  for (i = 0; i < input->scan; i++) {

    /* find the best motif in [i,49] motifs */
    maxs = motifs[i]->score;
    maxp = i;
    for (j = i+1; j < input->scan; j++) {
      if (maxs < motifs[j]->score) {
	maxp = j;
	maxs = motifs[j]->score;
      }
    }

    /* swap the highest motif to position i */
    if (maxp != i) {
      mtfptr = motifs[i];
      motifs[i] = motifs[maxp];
      motifs[maxp] = mtfptr;
    }

    if (motifs[i]->segment == 0) {
      /* highest motif i is already 0 */
      if (i < input->report)
	input->report = i;
      return;
    }

    /* calculate zscore */
    if (input->MonteCarlo) {
      motifs[i]->zscore = (motifs[i]->score - input->mean) / input->stdev;
    }

    /* now eliminate any motifs (i,49] similar to i */
    for (j = i+1; j < input->scan; j++) {
      if (motifs[j]->segment &&
	  motifsSimilar(motifs[i], motifs[j], input->w, input->mw)) {
	/*
	printf("Final motif %d eliminated motif %d\n", i+1, j+1);
	*/
	clearMotif(motifs[j], input);
      }
    }
  }
}

/******************************
 *
 * Func: scanning
 * Use the current motif to scan all the sequences, adding all segments
 * that could increase the motif score.
 *
 ******************************/
void scanning(struct sequence *seqs, struct inputParam *input, struct
	      motif **motifs, struct freq *bglog) {
  int i, j, seqno, end, *seq1, *seq2, w = input->w, *try1, *try2;
  double best, *bgs1, *bgs2, info, score, minth = input->minth;
  struct sequence *s;
  struct motif *mtfptr;
  struct alignment *a;
  
  for (i = 0; i < input->scan; i++) {
    if (motifs[i]->segment) {
      mtfptr = motifs[i];
      a = mtfptr->ahead;

      for (seqno = 0, s = seqs; seqno < input->top; 
	   seqno++, s = s->next) {
	if (SCANTOP) {
	  end = s->len - w + 1;
	  seq1 = s->seqi;
	  seq2 = s->rseqi;
	  bgs1 = s->bgs;
	  bgs2 = s->rbgs;
	  try1 = s->tried1;
	  try2 = s->tried2;
	  
	  /* motifInfo step does calcBlkLog */
	  best = motifInfo(mtfptr, input, bglog);
	  
	  for (j = 0; j < end; j++) {
	    if (a && a->s == s && a->position == j+1) {
	      a = a->next;
	    } else if (!try1[j]) {
	      score = seqMtfScore(&(seq1[j]), w, mtfptr->blk) - bgs1[j]; 
	      if (score > minth) {
		addSegment(mtfptr, &(seq1[j]), w, bgs1[j], ADD);
		info = motifInfo(mtfptr, input, bglog);
		if (info > best) {
		  addAlign(s, j+1, &(seq1[j]), bgs1[j], mtfptr);
		  best = info;
		} else
		  addSegment(mtfptr, &(seq1[j]), w, bgs1[j], SUBTRACT);
	      } 
	    }
	    if (a && a->s == s && a->position == -j-1) {
	      a = a->next;
	    } else if (!try2[j]) {
	      score = seqMtfScore(&(seq2[j]), w, mtfptr->blk) - bgs2[j]; 
	      if (score > minth) {
		addSegment(mtfptr, &(seq2[j]), w, bgs2[j], ADD);
		info = motifInfo(mtfptr, input, bglog);
		if (info > best) {
		  addAlign(s, -j-1, &(seq2[j]), bgs2[j], mtfptr);
		  best = info;
		} else
		  addSegment(mtfptr, &(seq2[j]), w, bgs2[j], SUBTRACT);
	      } 
	    }
	  }
	}
      }
      
      /* scan the rest of the sequences */
      for (; seqno < input->confirm; seqno++, s = s->next) {
	end = s->len - w + 1;
	seq1 = s->seqi;
	seq2 = s->rseqi;
	bgs1 = s->bgs;
	bgs2 = s->rbgs;
	try1 = s->tried1;
	try2 = s->tried2;
	
	/* motifInfo step does calcBlkLog */
	best = motifInfo(mtfptr, input, bglog);

	for (j = 0; j < end; j++) {
	  if (!try1[j]) {
	    score = seqMtfScore(&(seq1[j]), w, mtfptr->blk) - bgs1[j]; 
	    if (score > minth) {
	      addSegment(mtfptr, &(seq1[j]), w, bgs1[j], ADD);
	      info = motifInfo(mtfptr, input, bglog);
	      if (info > best) {
		addAlign(s, j+1, &(seq1[j]), bgs1[j], mtfptr);
		best = info;
	      } else
		addSegment(mtfptr, &(seq1[j]), w, bgs1[j], SUBTRACT);
	    } 
	  }
	  if (!try2[j]) {
	    score = seqMtfScore(&(seq2[j]), w, mtfptr->blk) - bgs2[j]; 
	    if (score > minth) {
	      addSegment(mtfptr, &(seq2[j]), w, bgs2[j], ADD);
	      info = motifInfo(mtfptr, input, bglog);
	      if (info > best) {
		addAlign(s, -j-1, &(seq2[j]), bgs2[j], mtfptr);
		best = info;
	      } else
		addSegment(mtfptr, &(seq2[j]), w, bgs2[j], SUBTRACT);
	    } 
	  }
	}
      }
    }
  }
}

/******************************
 *
 * Func: seqMtfScore
 * Score the segment of wid long starting from *seqi by blk
 *
 ******************************/
double seqMtfScore(int *seqi, int wid, struct mtfPt **blk) {
  double score = 0, temp;
  int i;
  
  for (i = 0; i < wid; i++) {
    if (seqi[i] == N)
      return 0;
    temp = blk[i][seqi[i]].log;
    score += temp;
  }
  return score;
}

/******************************
 *
 * Func: refine
 * Refine the motif, see whether deleting some segments will increase the
 * motif score.
 * 
 ******************************/
void refine(struct inputParam *input, struct motif **motifs, struct freq
	    *bglog) {
  int i, w = input->w;
  struct alignment *a1, *a2;
  double best, score;
  struct motif *mtf;
  
  for (i = 0; i < input->scan; i++) {
    mtf = motifs[i];
    if (mtf->segment != mtf->oldsegment) {
      mtf->oldsegment = mtf->segment;
      best = motifInfo(mtf, input, bglog);
      
      for (a1 = mtf->ahead, a2 = NULL; a1;) {
	addSegment(mtf, a1->st, w, a1->bginfo, SUBTRACT);
	score = motifInfo(mtf, input, bglog);
	
	/* delete the aligned segment */
	if (score > best) {
	  best = score;
	  if (a2) {
	    a2->next = a1->next;
	  } else {
	    mtf->ahead = a1->next;
	  }
	  if (!a1->next) 
	    mtf->atail = a2;
	  free(a1);
	  if (a2)
	    a1 = a2->next;
	  else
	    a1 = mtf->ahead;
	} else {
	  
	  /* keep the aligned segment */
	  addSegment(mtf, a1->st, w, a1->bginfo, ADD);
	  a2 = a1;
	  a1 = a1->next;
	}
      }
    } 
  }
}

/******************************
 *
 * Func: calcStats
 * Given the number of datapoints in stats[], calculate the mean, std
 *
 ******************************/
void calcStats(double *stats, struct inputParam *input) {
  int i, data = input->MonteCarlo;
  double sum = 0, mean;
  
  for (i = 0; i < data; i++)
    sum += stats[i];
  input->mean = mean = sum / data;
  sum = 0;
  for (i = 0; i < data; i++)
    sum += pow((stats[i] - mean), 2);
  input->stdev = sqrt(sum / (data - 1));
}
  
/******************************
 *
 * Func: printResults
 * Find the best report motifs from the motifs matrix.
 *
 ******************************/
void printResults(struct motif **motifs, struct inputParam *input, struct
		 sequence *seqs) {
  int i, report = input->report, w = input->w;
  /* motif score in results increases with index */
  struct motif *mtfptr;

  /* print the best motifs */
  for (i = 0; i < report; i++) {
    mtfptr = motifs[i];
    getConsensus(w, mtfptr);
    printMotif(mtfptr, input, i+1);
    if (input->printAlign) {
      if (input->printAlign > 1) 
	printShortAlignment(seqs, mtfptr, input);
      else
	printAlignment(mtfptr->ahead, input);
    } else 
      fprintf(input->ofp, "\n");
  }
  fflush(input->ofp);
}
  
/******************************
 *
 * Func: printMotif
 *
 ******************************/
void printMotif(struct motif *mtf, struct inputParam *input, int mtfCt) {
  int j, k, w = input->w;
  double temp;
  FILE *ofp = input->ofp;

  getConsensus(w, mtf);
  if (input->expMt > 0)
    fprintf(ofp, "Motif.P%d.%d.%d", input->expMt, w, mtfCt);
  else if (input->expMt < 0)
    fprintf(ofp, "Motif.N%d.%d.%d", -input->expMt, w, mtfCt);
  else
    fprintf(ofp, "Motif.%d.%d\t", w, mtfCt);
  fprintf(ofp, "\t%.3f\t%.3f\t%.3f\t%d\t%s\t%s\t%.1f\n", 
	    mtf->score, mtf->zscore, mtf->avgInfo, mtf->segment, 
	    mtf->constr, mtf->rconstr, mtf->siteTh);
  
  fprintf(ofp, "      A      C      G      T      Con  rCon Deg  rDeg \n");
  for (j = 0; j < w; j++) {
    fprintf(ofp, "%-4d", (j+1));
    for (k = 0; k < ALPHASIZE; k++) {
      temp = 100 * exp(mtf->blk[j][k].log);
      if (temp > 10)
	fprintf(ofp, "%.2f  ", temp);
      else 
	fprintf(ofp, " %.2f  ", temp);
    }
    fprintf(ofp, "   %c    %c    %c    %c\n", mtf->constr[j], 
	    mtf->rconstr[w-j-1], mtf->degstr[j], mtf->rdegstr[w-j-1]);
  }
}

/******************************
 *
 * Func: printShortAlignment
 *
 ******************************/
void printShortAlignment(struct sequence *seqs, struct motif *mtf,
			 struct inputParam *input) {
  int i, end, w = input->w, *seq1, *seq2;
  double thresh = mtf->siteTh, *bgs1, *bgs2, score;
  struct sequence *s;
  FILE *ofp = input->ofp;
  
  for (s = seqs; s; s = s->next) {
    end = s->len - w + 1;
    seq1 = s->seqi;
    seq2 = s->rseqi;
    bgs1 = s->bgs;
    bgs2 = s->rbgs;
    for (i = 0; i < end; i++) {
      score = exp(seqMtfScore(&(seq1[i]), w, mtf->blk) - bgs1[i]);
      if (score > thresh)
	fprintf(ofp, "%s %d %d", s->name, (i+1), (int)rint(score));
      score = exp(seqMtfScore(&(seq2[i]), w, mtf->blk) - bgs2[i]);
      if (score > thresh)
	fprintf(ofp, "%s %d %d", s->name, -(i+1), (int)rint(score));
    }
  }
  fprintf(ofp, "\n\n");
}

/******************************
 *
 * Func: printAlignment
 *
 ******************************/
void printAlignment(struct alignment *align, struct inputParam *input) {
  struct alignment *a;
  int i, w = input->w;
  char line[w+1];
  FILE *ofp = input->ofp;

  line[w] = ENDLINE;
  for (a = align; a; a = a->next) {
    for (i = 0; i < w; i++) {
      switch(a->st[i]) {
      case A: line[i] = 'A'; break;
      case T: line[i] = 'T'; break;
      case G: line[i] = 'G'; break;
      case C: line[i] = 'C'; break;
      }
    }
    if (a->position > 0)
      fprintf(ofp, "%s\tf%d\t%s\n", a->s->name, a->position, line);
    else
      fprintf(ofp, "%s\tb%d\t%s\n", a->s->name, -a->position, line);
  }
}

/******************************
 *
 * Func: scoreAll
 * Use top motifs, and score all the sequences to get
 * sequence score.
 *
 ******************************/
void scoreAll(struct inputParam *input, struct motif **motifs, 
	      struct sequence *seqs) {
  int i, j, report = input->report, w = input->w, end, *seq1, *seq2, 
    printAlign = input->printAlign, hits = input->hits, 
    lowp[report];
  struct motif *mtfptr;
  struct sequence *s;
  double score, *bgs1, *bgs2, temp1, temp2, lows[report][hits];
  FILE *ofp = input->xfp;

  memset(lowp, 0, report * sizeof(int));
  memset(lows, 0, report * hits * sizeof(double));
  /* print header
  fprintf(ofp, "SeqName");
 */
  for (i = 0; i < report; i++) {
    if (motifs[i]->segment) {
      if (input->expMt > 0) {
	if (i == 0)
	  fprintf(ofp, "Motif.P%d.%d.%d", input->expMt, w, i+1);
	else
	  fprintf(ofp, "\tMotif.P%d.%d.%d", input->expMt, w, i+1);
      } else if (input->expMt < 0) {
	if (i == 0)
	  fprintf(ofp, "Motif.N%d.%d.%d", -input->expMt, w, i+1);
	else
	  fprintf(ofp, "\tMotif.N%d.%d.%d", -input->expMt, w, i+1);
      } else {
	if (i == 0)
	  fprintf(ofp, "Motif.%d.%d", w, i+1);
	else
	  fprintf(ofp, "\tMotif.%d.%d", w, i+1);
      }
    }
  }
  fprintf(ofp, "\n");

  for (s = seqs; s; s = s->next) {
    end = s->len - w + 1;
    /*
    fprintf(ofp, "%s\t", s->name);
    */
    seq1 = s->seqi;
    seq2 = s->rseqi;
    bgs1 = s->bgs;
    bgs2 = s->rbgs;
    for (i = 0; i < report; i++) {
      if (motifs[i]->segment) {
	mtfptr = motifs[i];
	score = 0;
	for (j = 0; j < end; j++) {
	  temp1 = exp(seqMtfScore(&(seq1[j]), w, mtfptr->blk) - bgs1[j]);
	  temp2 = exp(seqMtfScore(&(seq2[j]), w, mtfptr->blk) - bgs2[j]); 
	  if (printAlign) {
	    if (lows[i][lowp[i]] < temp1) {
	      lows[i][lowp[i]] = temp1;
	      lowp[i] = findMin(lows[i], hits);
	    }
	    if (lows[i][lowp[i]] < temp2) {
	      lows[i][lowp[i]] = temp2;
	      lowp[i] = findMin(lows[i], hits);
	    }
	  }
	  score += temp1 + temp2;
	}
	if (score > 0) {
	  if (i == 0)
	    fprintf(ofp, "%.2f", log(score)/LOG2);
	  else
	    fprintf(ofp, "\t%.2f", log(score)/LOG2);
	} else {
	  if (i == 0) 
	    fprintf(ofp, "0");
	  else
	    fprintf(ofp, "\t0");
	}
      }
    }
    fprintf(ofp, "\n");
  }
  fclose(ofp);

  /* site cutoff score */
  if (printAlign) {
    for (i = 0; i < report; i++) {
      if (motifs[i]->segment) {
	/*
	printf("\nMotif %d: pos %d, score %.1f\n", i, lowp[i], lows[i][lowp[i]]);
	for (j = 0; j < hits; j++) {
	  printf("(%d,%.1f)", j, lows[i][j]);
	}
	*/
	motifs[i]->siteTh = lows[i][lowp[i]] - VERYSMALL;
      }
    }
  }
}

/******************************
 *
 * Func: findMin
 * find the minimum number in a list, and return its index
 *
 ******************************/
int findMin(double *list, int size) {
  int i, minp;
  double mins;
  
  minp = 0;
  mins = list[0];
  for (i = 1; i < size; i++) {
    if (list[i] < mins) {
      mins = list[i];
      minp = i;
    }
  }
  return minp;
}


/******************************
 *
 * Func: finishUp
 * Finish up, report the total time taken and close the log file
 *
 ******************************/
void finishUp(struct timeval *begintv, FILE *ofp) {
  long temp;
  struct timeval endtv;  

  gettimeofday(&endtv, NULL);
  temp = endtv.tv_sec - begintv->tv_sec;
  fprintf(ofp, "Total time %d:%d:%d.\n", (int)(temp / SPH), 
	  (int)(temp % SPH / MPH) , (int)(temp % MPH));
  fflush(ofp);
  fclose(ofp);
}
