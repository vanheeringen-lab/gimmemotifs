/******************************
 *
 * BioProspector: a DNA motif finder
 * Xiaole Liu, Stanford Medical Informatics
 * 4/15/2004
 *
 * BioProspector tries to find enriched sequence motifs from input
 * the motif could be single-block or two-block motif with a gap, or
 * palindrome motifs with a gap. After a number of resets, the
 * program reports highest scoring motifs. 
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

/* these mostly serve as flags */
#define FORWARD 1
#define PLDM 2
#define BACKWARD 3
#define NO 0
#define YES 1
#define SEQ 1
#define REG 0
#define BEST 1
#define ADD 1
#define SUBTRACT -1
#define BLK1 1
#define BLKBOTH 2
#define BLK2 3

/* Markov degree for input and background, general formula: 3 * 256 * 4^mkv */
#define MKV1 3 * 1024
#define MKV2 3 * 4096
#define MKV3 3 * 16384

/* for random number generator */
#define IA 16807
#define IM 2147483647
#define AM (1.0 / IM)
#define IQ 127773
#define IR 2836

/* number of data points, resets to do */
#define RESET 5
#define TOTAL 40
#define TRIAL 10

/* simple sizes */
#define DOUBLESIZE sizeof(double)
#define INTSIZE sizeof(int)
#define MTFCOLSIZE 4 * sizeof(struct mtfPt)
#define ALIGNSIZE sizeof (struct alignment)

/* some ascii for keys */
#define ENDLINE '\0'
#define TAB 9
#define SPACE 32
#define TILTA 126

/* max/min numbers */
#define MAXLINE 32767
#define MINALIGN 5
#define MINWID 4
#define MAXWID 50
#define MAXGAP 50
#define MAXGAPRANGE 20
#define MINE 30
#define MAXE 6000
#define MAXZ 30
#define MINTRY 10
#define MAXTRY 200

/* some input defaults */
#define DEFAULTWID 10
#define DEFAULTMTFWIDTH 8
#define DEFAULTPERCENT 0.05
#define DEFAULTCOPY 4
#define DEFAULTRESULT 5

/* other useful constants */
#define SEQTHRESHRATIO 8
#define VERYSMALL 0.000000025
#define EXTRACOL 4
#define MPH 60
#define SPH 3600

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
  int w1, w2, gl, gm; /* user input */
  int meangap, grange, longest, noseq, bgmkv, seqmkv; /* calculated */
  int dir, pldm, all, refine; /* actually booleans */
  int maxw1, maxw2; /* width used in motif finding */
  int local; /* for web applications, local = 0 */
  char jobID[3]; /* for web application only */

  /* number of randomizations */
  int iteration, tries, result;
  
  /* threshold related parameters */
  double seqth1, seqth2, inc;
  double expect; /* expected bases per site */
  double scount[5], percent;

  /* file name and file handles */
  FILE *is, *bs, *bf, *ofp;
};

/******************************
 * Data type: seqCount
 * Stores the count of letters from either input or background sequences
 ******************************/
struct seqCount {
  int one[ALPHASIZE], two[ALPHASIZE][ALPHASIZE], 
    three[ALPHASIZE][ALPHASIZE][ALPHASIZE], 
    four[ALPHASIZE][ALPHASIZE][ALPHASIZE][ALPHASIZE];
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

  /* sequence name, len, letters, motif site copies */
  char *name;
  int *seqi, *rseqi, len;
  double *bgscore1, *rbgscore1, *bgscore2, *rbgscore2;

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
 * Stores the position specific matrices, number of aligned sites and
 * some arrays for storing doubles during update calculations. 
 ******************************/
struct motif {

  /* 1 for first block, 2 for second block */
  struct mtfPt **blk1, **blk2;

  /* information contained in columns of the two blocks and motif score */
  double mtfscore;

  /* alignment parameters */
  struct alignment *align;
  int site;

  /* 1 for forward direction, 2 for backward direction. */
  double *pb12, *pb22, *pb1, *pb2, **mat1, **mat2;
};

/******************************
 * Data type: alignment
 * stores the alignment parameters as linked list (each node is for one
 * sequence) in the motif structure.
 ******************************/
struct alignment {
  int copy;
  int size;
  int *dir;
  int *st1;
  int *st2;
  struct alignment *next;
};

/******************************
 *
 * External functions
 *
 ******************************/

int execl (const char *filename, const char *arg0, ...);
int fflush (FILE *stream);
int getopt(int, char *const *, const char *);
double rint (double x);

/******************************
 *
 * Function prototypes (in alphabetical order)
 *
 ******************************/

void addcopy(int st1, int st2, struct alignment *align, int dir, int blk2);
void addSite(struct mtfPt **blk, int *start, int w, int add);
double blkBgScore(struct sequence *seqs, struct alignment *align, int id);
void blkColShift(struct motif *mtf, int shift1, int shift2, struct sequence
		 *seqs, struct inputParam *input);
void blkLeftShift(struct motif *mtf, struct sequence *seqs, struct
		  inputParam *input);
double blkMtfScore(struct mtfPt **blk, int wid);
void blkRightShift(struct motif *mtf, struct sequence *seqs, struct
		   inputParam *input);
void blkShift(struct inputParam *input, struct sequence *seqs, struct freq
	      *bglog, struct motif *mtf, struct motif *bkmtf1, struct
	      motif *bkmtf2);
void calcBlkLog(struct mtfPt **blk, int site, int wid, double *sct);
void calcSeqLog(struct freq *seqlog, struct seqCount *seqct, int *mkv);
void clearMotif(struct motif *mtf, struct inputParam *input);
void colShift(struct inputParam *input, struct sequence *seqs, struct
	      motif *mtf, struct motif *bkmtf1, struct motif *bkmtf2);
void countSeq(int *seqint, int len, struct seqCount *seqct);
double drand();
void errorExit(char *msg);
void findMotif(struct inputParam *input, struct freq *bglog, struct motif 
		*mtf, struct motif *bkmtf1, struct motif *bkmtf2, struct 
		motif *bestmtf, struct sequence	*seqs, struct motif **results);
void finishUp(struct timeval *begintv, FILE *ofp);
void getAllBgScores(struct sequence *seqs, struct freq *bglog, struct
		    inputParam *input);
void getSeqBgScore(double *segscore, int *seqi, int st, int end, int w,
		   int mkv, struct freq *bglog);
void getConsensus(int w, struct mtfPt **blk, char *con, char *rcon);
void initMotif(struct motif **mtf, struct inputParam *input);
void initRand();
double motifInfo(struct motif *mtf, struct inputParam *input, struct
		 sequence *seqs); 
void mtfCopy(struct motif *dst, struct motif *src, struct inputParam
	     *input);
int mystrlen(char *string);
void parseInputParam(struct inputParam *input, int argc, char *argv[]);
void printAlignment(struct inputParam *input, struct sequence *seqs,
		    struct alignment *align, int w1, int w2);
void printBlk(int w, struct mtfPt **blk, FILE *ofp);
void printMotif(struct motif *mtf, struct inputParam *input);
void printUsage();
void randomStart(struct inputParam *input, struct sequence *seqs, struct
		 motif *mtf);
void readInputFiles(struct inputParam *input, struct freq *bglog, struct 
		    sequence **seqs); 
int readOneSeq(FILE *fp, char *line, char *seq); 
void readSeqFile(FILE *fp, int isSeq, struct sequence **seqs, struct 
		 seqCount *seqct, struct inputParam *input);
void readBgFreq(FILE *bf, struct freq *bglog);
void readBgFreqLine(FILE *bf, double *bgfreq);
double seqMtfScore(int *seqi, int wid, struct mtfPt **blk);
double seqBgScore(int *seqi, int st, int wid, int mkv, struct freq *bg);
void seqatoi(char *line, int *seqint, int len, int forward);
void takeInOut(struct sequence *seq, struct motif *mtf, struct alignment
	       *align, struct inputParam *input, int id, int add);
void updateB1(struct inputParam *input, struct motif *mtf, struct
	      alignment *align, struct sequence *seq, int opt, struct
	      sequence *seqs);
void updateS1(struct inputParam *input, struct motif *mtf, struct
	      alignment *align, struct sequence *seq, int opt, struct
	      sequence *seqs); 
void updateB2(struct inputParam *input, struct motif *mtf, struct
	      alignment *align, struct sequence *seq, int opt, struct
	      sequence *seqs); 
void updateS2(struct inputParam *input, struct motif *mtf, struct
	      alignment *align, struct sequence *seq, int opt, struct
	      sequence *seqs); 
void updateP(struct inputParam *input, struct motif *mtf, struct
	     alignment *align, struct sequence *seq, int opt, struct
	      sequence *seqs); 

/******************************
 *
 * Function pointers
 *
 ******************************/

void (*update)(struct inputParam *, struct motif *, struct alignment *,
	       struct sequence *, int, struct sequence *); 

/******************************
 *
 * global variables
 *
 ******************************/

struct inputParam input;
struct freq bglog;
struct sequence *seqs;
struct motif *mtf, *bkmtf1, *bkmtf2, *bestmtf, **results;
long idum;  /* for random number generator */

/******************************
 *
 * Functions (in the order of which they are called)
 *
 ******************************/

/******************************
 *
 * Func: main
 *  1. Read input, allocate memory for all variables
 *  2. Find out callback functions and number of randomizations to do
 *  3. Find motif score distribution
 *  4. Search motif from input sequences
 *
 ******************************/
int main (int argc, char *argv[]) {
  struct timeval begintv;
  int i;

  /* read the input and initialize variables */
  gettimeofday(&begintv, NULL);
  initRand();
  parseInputParam(&input, argc, argv);

  /* print some header to make the output look better */
  fprintf(input.ofp, "**************************************** \n");
  fprintf(input.ofp, "*                                      * \n");
  fprintf(input.ofp, "*      BioProspector Search Result     * \n");
  fprintf(input.ofp, "*                                      * \n");
  fprintf(input.ofp, "**************************************** \n\n");
  fflush(input.ofp);
  if (input.local) {
    printf("Read input sequences.\n");
  }

  readInputFiles(&input, &bglog, &seqs);
  initMotif(&mtf, &input);
  initMotif(&bkmtf1, &input);
  initMotif(&bkmtf2, &input);
  initMotif(&bestmtf, &input);
  
  if (input.pldm)
    update = updateP;
  else if (input.dir) {
    if (input.w2)
      update = updateB2;
    else
      update = updateB1;
  } else {
    if (input.w2)
      update = updateS2;
    else
      update = updateS1;
  }

  /* start looking for real motifs */
  if (input.local) {
    printf("Look for motifs from the original sequences.\n");
  }

  /* create space for the final results */
  if (input.local) {
    printf("Print results.\n");
  }
  results = (struct motif **)malloc(input.result * sizeof (struct motif *));
  for (i = 0; i < input.result; i++) {
    initMotif(&(results[i]), &input);
  }
  getAllBgScores(seqs, &bglog, &input);
  findMotif(&input, &bglog, mtf, bkmtf1, bkmtf2, bestmtf, seqs,
	     results);
  
  finishUp(&begintv, input.ofp);
  return 0;
}

/******************************
 *
 * Func: initRand
 * Initialize the random generator seed to be a non-zero long
 *
 * Checked 10/10/00
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
 * Func: parseInputParam
 * Read in all the user specified parameters
 *
 * Checked 3/16/01
 ******************************/
void parseInputParam(struct inputParam *input, int argc, char *argv[]) {
  extern char *optarg;
  double dbltemp;
  int inttemp, opt;
  
  /* init input */
  memset(input, 0, sizeof(struct inputParam));
  input->result = DEFAULTRESULT;
  input->dir = YES;
  input->tries = TOTAL;
  input->local = 1;

  /* input file pointers */
  input->is = input->bs = input->bf = NULL;
  input->ofp = stderr;

  /* if no arguments, print out usage */
  if (argc == 1)
    printUsage();

  /* read input parameters and specify input */
  while ((opt = getopt(argc, argv, "W:w:G:g:p:d:i:b:f:a:o:r:n:e:h:y:l:"))
	 != EOF) {
    switch(opt) {
      
      /* first block width */
    case 'W': if (sscanf(optarg, "%d", &inttemp) && inttemp >= MINWID
		  && inttemp <= MAXWID)
      input->w1 = inttemp;
    else
      errorExit("Wrong input: first block width -W.");
      break;
      
      /* second block width */
    case 'w': if (sscanf(optarg, "%d", &inttemp) && inttemp >= MINWID
		  && inttemp <= MAXWID)
      input->w2 = input->w2 = inttemp;
    else
      errorExit("Wrong input: second block width -w.");
      break;
      
      /* max gap between blocks */
    case 'G': if (sscanf(optarg, "%d", &inttemp) && inttemp >= 0
		  && inttemp <= MAXGAP)
      input->gm = inttemp;
    else
      errorExit("Wrong input: max gap between blocks -G.");
      break;
      
      /* min gap between blks */
    case 'g': if (sscanf(optarg, "%d", &inttemp) && inttemp >= 0
		  && inttemp <= MAXGAP)
      input->gl = inttemp;
    else
      errorExit("Wrong input: min gap between blocks -g.");
      break;
      
      /* input sequence file */
    case 'i': if (!(input->is = fopen(optarg, "r"))) 
      errorExit("Can't open input sequence file.");
      break;
      
      /* background sequence file */
    case 'b': if (!(input->bs = fopen(optarg, "r"))) 
      errorExit("Can't open background sequence file."); 
      break;
      
      /* background probability distribution file */
    case 'f': if (!(input->bf = fopen(optarg, "r"))) 
      errorExit("Can't open background frequency file."); 
      break;
      
      /* whether to refine the motifs */
    case 'h': input->refine = YES; break;
      
      /* whether this is a palindrome motif */
    case 'p': input->pldm = YES; break;
      
      /* whether check both direction of the sequences */
    case 'd': if (sscanf(optarg, "%d", &inttemp) && inttemp == 1)
      input->dir = NO;
    else 
      input->dir = YES; 
      break;
      
      /* whether each sequence contains at least a copy of the motif */
    case 'a': input->all = YES; break;
      
      /* output file name */
    case 'o': if (!(input->ofp = fopen(optarg, "w"))) 
      errorExit("Can't open output file.");
      break;
    
      /* number of times trying to find motifs in original sequences */
    case 'n': if (sscanf(optarg, "%d", &inttemp) && inttemp >= MINTRY
		  && inttemp <= MAXTRY)
      input->tries = inttemp;
      break;

      /* number of top motifs to report */
    case 'r': if (sscanf(optarg, "%d", &inttemp) && inttemp > 0
		  && inttemp <= MAXTRY)
      input->result = inttemp;
      break;

      /* expected bases/site to use Bayes score*/
    case 'e': if (sscanf(optarg, "%lf", &dbltemp) && dbltemp >= MINE
		  && dbltemp <= MAXE)
      input->expect = dbltemp;
    else
      errorExit("Wrong expected bases/site.");
      break;

      /* jobID */
    case 'y': if (sscanf(optarg, "%d", &inttemp) && inttemp >= 0
		  && inttemp < 50)
      strcpy(input->jobID, optarg);
    else 
      errorExit("Wrong input: jobID -y.");
      break;

      /* local or web application */
    case 'l': if (sscanf(optarg, "%d", &inttemp) && inttemp > -1
		  && inttemp < 2)
      input->local = inttemp;
    else 
      errorExit("Wrong input: local -l.");
      break;
    }
  }
    
  /* check for invalid inputs */ 
  if (!input->is) 
    errorExit("No input sequence file specified.");

  if (!input->w1) {
    if (input->w2) {
      input->w1 = input->w2;
      input->w2 = 0;
    } else
      input->w1 = DEFAULTWID;
  }
  
  /* report as much as was ever tried */
  if (input->result > input->tries)
    input->result = input->tries;

  /* for palindrome motifs, the two blk widths have to be the same */
  if (input->pldm) {
    input->dir = YES;    
    input->w1 = input->w2 = (input->w1 + input->w2 + 1) / 2; 
    if (input->w1 < MINWID)
      errorExit("Motif width too small.");
  }

  /* calculating mean gap and gap range */
  if (input->w2) {
    if (input->gm < input->gl) {
      inttemp = input->gm;
      input->gm = input->gl;
      input->gl = inttemp;
      if (input->gm - input->gl > MAXGAPRANGE) 
	errorExit("Gap estimate range too big.");
    }
    input->meangap = (input->gm + input->gl + 1) / 2;
    input->grange = input->gm - input->gl + 1;
  }

  /* get the maximum width of the motif blks */
  input->maxw1 = input->w1 + EXTRACOL;
  input->maxw2 = input->w2 + EXTRACOL;
}

/******************************
 *
 * Func: printUsage
 * Print how to set options if the user does not give any parameters to
 * BioProspector. 
 *
 * Checked 5/14/01
 ******************************/
void printUsage() {

  fprintf(stderr, "\nBioProspector 4/15/04.\n");
  fprintf(stderr, "Usage: ./BioProspector -i seqfile (options) \n");
  fprintf(stderr, " Seqfile must be in restricted FASTA format.\n\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -W \t<motif width (default %d)>\n", DEFAULTWID);
  fprintf(stderr, "  -o \t<output file (default stdout)>\n");
  fprintf(stderr, "  -f \t<background distribution file (default seqfile)>\n");
  fprintf(stderr, "  -b \t<background sequence file (default input sequences)>\n");
  fprintf(stderr, "  -n \t<number of times trying to find motif (default %d)>\n", TOTAL);
  fprintf(stderr, "  -r \t<number of top motifs to report (default %d)>\n", DEFAULTRESULT);
  fprintf(stderr, "  -w \t<second motif block width for two-block motif (default 0)>\n");
  fprintf(stderr, "  -p 1 \t[if two-block motif is palindrome (default 0)]\n");
  fprintf(stderr, "  -G \t<max gap between two motif blocks (default 0)>\n");
  fprintf(stderr, "  -g \t<min gap between two motif blocks (default 0)>\n");
  fprintf(stderr, "  -d 1 \t[if only need to examine forward (default 2)]\n");
  fprintf(stderr, "  -a 1 \t[if every sequence contains the motif (default 0)]\n");    
  fprintf(stderr, "  -h 1 \t[if want more degenerate sites (default fewer sites)]\n");
  fprintf(stderr, "  -e \t<expected bases per motif site in the sequences \n\t    (will use Bayes motif scoring, don't specify if unknown)>\n");
  errorExit("");
}

/******************************
 *
 * Func: errorExit
 * Print the msg, and exit the program 
 *
 * Checked 10/10/00
 ******************************/
void errorExit(char *msg) {
  if (input.local) {
    fprintf(stderr, "%s\n", msg);
    exit(0);
  } else {
    fprintf(input.ofp, "%s\n", msg);
    fclose(input.ofp);
    execl("./BPemail.pl", "BPemail", input.jobID, NULL);
  }
}

/******************************
 *
 * Func: readInputFiles
 * Read probability of sequence input/bg files, initialize seqs, bg/seq
 * and log  
 * Note: It is important to read Seq file before Bg file, because if the
 * user does not specify the Bg file, then information in the Seq file
 * (stored in seqCt) will be used to calculate bglog.
 *
 * Checked 3/16/01
 ******************************/
void readInputFiles(struct inputParam *input, struct freq *bglog, struct 
		    sequence **seqs) { 
  int i;
  struct seqCount seqCt;

  /* read in input sequences, and calculate their base distribution */
  readSeqFile(input->is, SEQ, seqs, &seqCt, input);

  /* calculate background base distribution */
  if (input->bf) { /* reading from background probability file */
    input->bgmkv = 3;
    readBgFreq(input->bf, bglog);
  } else {
    if (input->bs) { /* reading from background sequence file */
      readSeqFile(input->bs, NO, NULL, &seqCt, input);
    } 

    /* keep this line here in case user wants to use input as bg */
    calcSeqLog(bglog, &seqCt, &(input->bgmkv));
  }

  /* calculate psudocount from background distribution */
  for (i = 0; i < ALPHASIZE; i++)
    input->scount[i] = exp(bglog->one[i]) * input->scount[ALPHASIZE];
}

/******************************
 *
 * Func: readSeqFile
 * Allocate space for sequences, put input sequences in a linked list,
 * count the based distribution and put it in seqCount
 *
 * Checked 3/16/01
 ******************************/
void readSeqFile(FILE *fp, int isSeq, struct sequence **seqs, struct 
		 seqCount *seqct, struct inputParam *input) {
  char line[MAXLINE], seq[MAXLINE];
  int seqint[MAXLINE], len, sumlen = 0;
  int noseq = 0, eval, avglen;
  struct sequence *s1 = NULL, *s2 = NULL; /* temp sequence pointers */

  memset(seqct, 0, sizeof(struct seqCount));
  while (1) {
    eval = readOneSeq(fp, line, seq);
    if (eval == 0) break;
    if (eval == -1) continue;

    /* if the sequence is a input sequence, allocate spaces */
    if (isSeq) {
      s1 = (struct sequence *) malloc (sizeof (struct sequence));
      len = mystrlen(line);
      s1->name = (char *) malloc(len + 1);
      strcpy(s1->name, line);
      len = mystrlen(seq); 
      noseq++;
      sumlen += s1->len = len;
      if (len > input->longest) input->longest = len;
      s1->seqi = (int *) malloc (len * INTSIZE);
      seqatoi(seq, s1->seqi, len, FORWARD); 
      countSeq(s1->seqi, len, seqct); 
      s1->bgscore1 = (double *)malloc(len * DOUBLESIZE);
      if (input->w2)
	s1->bgscore2 = (double *)malloc(len * DOUBLESIZE);
      if (input->dir) {
	s1->rseqi = (int *) malloc (len * INTSIZE);
	seqatoi(seq, s1->rseqi, len, BACKWARD); 
	countSeq(s1->rseqi, len, seqct);
	if (!input->pldm) {
	  s1->rbgscore1 = (double *)malloc(len * DOUBLESIZE);
	  if (input->w2)
	    s1->rbgscore2 = (double *)malloc(len * DOUBLESIZE);
	}
      }
      s1->next = NULL;
      if (s2) {
	s2->next = s1;
	s2 = s1;
      } else
	*seqs = s2 = s1;
    } else {

      /* if a background sequence, just count its A,T,G,C. */
      len = mystrlen(seq);       
      seqatoi(seq, seqint, len, FORWARD); 
      countSeq(seqint, len, seqct); 
      if (input->dir) {
	seqatoi(seq, seqint, len, BACKWARD); 
	countSeq(seqint, len, seqct);
      }
    }
  }

  /* figure out total psudocount, percent, thresholds and increment */
  if (isSeq) {
    input->noseq = noseq;
    input->scount[ALPHASIZE] = ((double) MINALIGN) / noseq;

    /* specify percent for sampling */
    if (!input->all) {
      if (MINALIGN > noseq * DEFAULTPERCENT) {
	input->percent = DEFAULTPERCENT;
      } else {
	input->percent = ((double) MINALIGN) / noseq;
      }
    }

    /* input sequence length have direct relationship on the number of
     * Gibbs update iterations */
    avglen = sumlen / noseq;
    input->iteration = (int) ceil(pow(2 * avglen, 0.5));
    input->seqth1 = avglen * (input->w1 + input->w2);
    if (input->w2) input->seqth1 /= 2;
    if (input->dir) input->seqth1 *= 2;
    input->inc = input->seqth1 / input->iteration / SEQTHRESHRATIO;
  }
  fclose(fp);
}

/******************************
 *
 * Func: readOneSeq
 * Read one line from the sequence file.
 * Output: 
 *  -1 - empty line
 *   0 - reach end of file
 *   1 - regular sequence, in this case line will contain sequence name and
 *       seq will contain actual sequence.       
 *
 * Note: Right now the program can not really read in sequences in some
 * FASTA that have sequences in multiples lines. The cgi program will make 
 * sure that input sequences follow this rule
 *
 * Checked 3/16/01
 ******************************/
int readOneSeq(FILE *fp, char *line, char *seq) {
  
  if (fgets(line, MAXLINE, fp)) {
    if (mystrlen(line) == 0)
      return -1;
    
    if (line[0] == '>') { /* FASTA format */
      if (!fgets(seq, MAXLINE, fp)) return 0;
      else {
	switch(seq[0]) {
	case 'a':;
	case 't':;
	case 'g':;
	case 'c':;
	case 'A':;
	case 'T':;
	case 'G':;
	case 'C':;
	case 'n':;
	case 'N': break;
	default: 
	  if (input.local)
	    fprintf(stderr, "%s\n", line);
	  errorExit("Sequences not in FASTA format.");
	}
      }
    } else {
      if (input.local)
	fprintf(stderr, "%s\n", line);
      errorExit("Wrong sequence format.");
    }
    return 1;
  } else
    return 0;
}

/******************************
 *
 * Func: mystrlen
 * Return the input sequence name/seq length 
 *
 * Checked 3/16/01
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
 * Checked 10/12/00
 ******************************/
void seqatoi(char *line, int *seqint, int len, int forward) {
  int i;
  double random;

  for (i = 0; i < len; i++) {
    switch (line[i]) {
    case 'a':;
    case 'A': if (forward < PLDM) seqint[i] = A; 
    else seqint[len-i-1] = T; break;
    case 't':;
    case 'T': if (forward < PLDM) seqint[i] = T; 
    else seqint[len-i-1] = A; break;
    case 'g':;
    case 'G': if (forward < PLDM) seqint[i] = G; 
    else seqint[len-i-1] = C; break;
    case 'c':;
    case 'C': if (forward < PLDM) seqint[i] = C; 
    else seqint[len-i-1] = G; break;

      /* If unknown character, randomly assign some. This is only done 
       * in forward which assign a character for the backward to use. */
    default: random = drand();
      if (random < 0.25) {
	seqint[i] = A;
	line[i] = 'A';
      } else if (random < 0.5) {
	seqint[i] = T;
	line[i] = 'T';
      } else if (random < 0.75) {
	seqint[i] = G;
	line[i] = 'G';
      } else {
	seqint[i] = C;
	line[i] = 'C';
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
 * Checked 10/18/00
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
 * Func: countSeq
 * Count the base distribution of a sequence in single, double, triplet,
 * and quadruplet bases
 *
 * Checked 10/12/00
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
 * Checked 10/12/00
 ******************************/
void calcSeqLog(struct freq *seqlog, struct seqCount *seqct, int *mkv) {
  int i, j, k, l;
  int sum1, sum2, sum3, sum4;

  /* figure out the Markov order based on total sequence size */
  sum1 = 0;
  for (i = 0; i < 4; i++)
    sum1 += seqct->one[i];
  if (sum1 > MKV3)
    *mkv = 3;
  else if (sum1 > MKV2)
    *mkv = 2;
  else if (sum1 > MKV1)
    *mkv = 1;
  else
    *mkv = 0;

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
 * Func: readBgFreq
 * Read in the background distribution directly from a probability file
 *
 * Checked 10/12/00
 ******************************/
void readBgFreq(FILE *bf, struct freq *bglog) {
  int i, j, k;
  
  /* automatically use mkv order of 3 if read freq from bg file */
  readBgFreqLine(bf, bglog->one);
  for (i = 0; i < ALPHASIZE; i++)
    readBgFreqLine(bf, bglog->two[i]);
  for (i = 0; i < ALPHASIZE; i++)
    for (j = 0; j < ALPHASIZE; j++)
      readBgFreqLine(bf, bglog->three[i][j]);
  for (i = 0; i < ALPHASIZE; i++)
    for (j = 0; j < ALPHASIZE; j++)
      for (k = 0; k < ALPHASIZE; k++)
	readBgFreqLine(bf, bglog->four[i][j][k]);
}

/******************************
 *
 * Func: readBgFreqLine
 * Read in one line of the background distribution from a background
 * frequency file. The four numbers are the log ratio of A, C, G, T in a
 * particular Markov order.
 *
 * Checked 10/12/00
 ******************************/
void readBgFreqLine(FILE *bf, double *bgfreq) {
  if (!fscanf(bf, "%lf %lf %lf %lf\n", &bgfreq[0], &bgfreq[1], 
	      &bgfreq[2], &bgfreq[3]))
    errorExit("Background frequency file doesn't have right format.");
}

/******************************
 *
 * Func: initMotif
 * Setup dynamic space for motif structure. Double array pb** are useful
 * during update for storing site scores and sampling.
 *
 * Checked 10/18/00
 ******************************/
void initMotif(struct motif **mtf, struct inputParam *input) { 
  int i, range = input->grange+2, longest = input->longest, noseq =
    input->noseq;
  struct alignment *a1, *a2;
  int wid1 = input->maxw1, wid2 = input->maxw2;
  int size1 = sizeof(double) * longest;
  int size2 = sizeof(double *) * longest;
  int size3 = sizeof(struct mtfPt *);
  int size4 = sizeof(struct motif);
  
  *mtf = (struct motif *) malloc(size4);
  (*mtf)->mtfscore = 0;
  (*mtf)->site = 0;
  (*mtf)->align = NULL;
  for (i = 0; i < noseq; i++) {
    a1 = (struct alignment *)malloc(ALIGNSIZE);
    a1->copy = 0;
    a1->size = DEFAULTCOPY;
    a1->st1 = (int *)malloc(DEFAULTCOPY * INTSIZE);
    a1->dir = (int *)malloc(DEFAULTCOPY * INTSIZE);
    a1->next = NULL;
    if ((*mtf)->align)
      a2 = a2->next = a1;
    else
      (*mtf)->align = a2 = a1;
  }
  (*mtf)->blk1 = (struct mtfPt **)malloc(wid1 * size3);
  for (i = 0; i < wid1; i++)
    (*mtf)->blk1[i] = (struct mtfPt *)malloc(MTFCOLSIZE);
  (*mtf)->pb1 = (double *)malloc(size1);
  if (input->w2) {
    (*mtf)->blk2 = (struct mtfPt **)malloc(wid2 * size3);
    for (i = 0; i < wid2; i++)
      (*mtf)->blk2[i] = (struct mtfPt *)malloc(MTFCOLSIZE);
    (*mtf)->pb12 = (double *)malloc(size1);				 
    (*mtf)->mat1 = (double **)malloc(size2);
    for (i = 0; i < longest; i++)
      (*mtf)->mat1[i] = (double *)malloc(range * DOUBLESIZE);
  }
  if (input->w2) {
    for (a1 = (*mtf)->align; a1; a1 = a1->next)
      a1->st2 = (int *)malloc(DEFAULTCOPY * INTSIZE);
  }
  if (input->dir && !input->pldm) {
    (*mtf)->pb2 = (double *)malloc(size1);
    if (input->w2) { 
      (*mtf)->pb22 = (double *)malloc(size1);
      (*mtf)->mat2 = (double **)malloc(size2);
      for (i = 0; i < longest; i++)
	(*mtf)->mat2[i] = (double *)malloc(range * DOUBLESIZE);
    }
  }
}

/******************************
 *
 * Func: getAllBgScores
 * 
 * 
 * Checked /7/2/01
 ******************************/
void getAllBgScores(struct sequence *seqs, struct freq *bglog, struct
		    inputParam *input) {
  struct sequence *s;
  int w1 = input->w1, w2 = input->w2, mkv = input->bgmkv, pldm =
    input->pldm, dir = input->dir;
  int st = 0, end1, end2;
  
  for (s = seqs; s; s = s->next) {
    end1 = s->len - w1;
    end2 = s->len - w2;
    getSeqBgScore(s->bgscore1, s->seqi, st, end1, w1, mkv, bglog);
    if (dir) {
      if (!pldm)
	getSeqBgScore(s->rbgscore1, s->rseqi, st, end1, w1, mkv, bglog);
      else
	getSeqBgScore(s->bgscore2, s->rseqi, st, end1, w1, mkv, bglog);
    }
    if (w2 && !pldm) {
      getSeqBgScore(s->bgscore2, s->seqi, st, end2, w2, mkv, bglog);
      if (dir) 
	getSeqBgScore(s->rbgscore2, s->rseqi, st, end2, w2, mkv,
		      bglog);
    }
  }
}

/******************************
 *
 * Func: getSeqBgScore
 *
 * Checked 7/2/01
 ******************************/
void getSeqBgScore(double *segscore, int *seqi, int st, int end, int w,
		   int mkv, struct freq *bglog) {
  int i;

  for (i = st; i <= end; i++)
    segscore[i] = seqBgScore(&(seqi[i]), i, w, mkv, bglog);
}

/******************************
 *
 * Func: seqBgScore
 * Score the site of wid long starting from *seqi by background
 *
 * Checked 5/4/01
 ******************************/
double seqBgScore(int *seqi, int st, int wid, int mkv, struct freq *bg) {
  double score = 0;
  int i, *ptr = seqi;

  switch(mkv) {
  case 0: for (i = 0; i < wid; i++) score += bg->one[*(ptr+i)]; 
    break;
  case 1: 
    if (st == 0) {
      score += bg->one[*ptr];
      for (i = 1; i < wid; i++) score += bg->two[*(ptr+i-1)][*(ptr+i)]; 
    } else {
      for (i = 0; i < wid; i++) score += bg->two[*(ptr+i-1)][*(ptr+i)];
    }
    break;
  case 2: 
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
 * Func: clearMotif
 * Clear the motif matrix and sites. 
 * 
 * Checked 3/19/01
 ******************************/
void clearMotif(struct motif *mtf, struct inputParam *input) {
  int i, w = input->maxw1;

  mtf->site = 0;
  for (i = 0; i < w; i++)
    memset(mtf->blk1[i], 0, MTFCOLSIZE);

  if (input->w2 && !input->pldm) {
    w = input->maxw2; 
    for (i = 0; i < w; i++)
      memset(mtf->blk2[i], 0, MTFCOLSIZE);
  }
}

/******************************
 *
 * Func: randomStart
 * Initialized the content of the motif matrix by randomly choosing
 * alignment from the input sequences. 
 *
 ******************************/
void randomStart(struct inputParam *input, struct sequence *seqs, struct
		 motif *mtf) {
  int wid = input->w1 + input->w2 + input->meangap;
  struct sequence *s;
  struct alignment *a;

  input->seqth2 = 0;
  for (s = seqs, a = mtf->align; s; s = s->next, a = a->next) {
    if (s->len > wid) {
      a->copy = 1;
      a->st1[0] = (int) (drand() * (s->len - wid));
      if (input->w2) {
	if (input->pldm) {
	  a->st2[0] = s->len - a->st1[0] - wid;
	  a->dir[0] = PLDM;
	  takeInOut(s, mtf, a, input, BLK1, ADD);
	} else {
	  a->st2[0] = a->st1[0] + wid - input->w2;
	  a->dir[0] = FORWARD;
	  takeInOut(s, mtf, a, input, BLKBOTH, ADD);
	}
      } else {
	a->dir[0] = FORWARD;
	takeInOut(s, mtf, a, input, BLK1, ADD);
      }
    } else {
      a->copy = 0;
    }
  }
}

/******************************
 *
 * Func: takeInOut
 * Add all copies of aligned site from each sequence to the motif
 * blocks, and update the motif site number.
 *
 * Checked 3/16/01
 ******************************/
void takeInOut(struct sequence *seq, struct motif *mtf, struct alignment
	       *align, struct inputParam *input, int id, int add) {
  int i, w1 = input->w1, w2 = input->w2;
  
  for (i = 0; i < align->copy; i++) {
    if (id < BLK2) {
      switch(align->dir[i]) {
      case FORWARD: 
	addSite(mtf->blk1, &(seq->seqi[align->st1[i]]), w1, add); 
	break;
      case PLDM: 
	addSite(mtf->blk1, &(seq->seqi[align->st1[i]]), w1, add); 
	addSite(mtf->blk1, &(seq->rseqi[align->st2[i]]), w1, add);
	break;
      case BACKWARD: 
	addSite(mtf->blk1, &(seq->rseqi[align->st1[i]]), w1, add); 
	break;
      }
    }
    if (w2 && id > BLK1 && !input->pldm) {
      /* only comes here for two-blk non-pldm motifs */
      
      switch(align->dir[i]) {
      case FORWARD: 
	addSite(mtf->blk2, &(seq->seqi[align->st2[i]]), w2, add); break;
      case BACKWARD: 
	addSite(mtf->blk2, &(seq->rseqi[align->st2[i]]), w2, add); break;
      }
    }
  }
  if (input->pldm) mtf->site += 2 * add * align->copy;
  else mtf->site += add * align->copy;
}

/******************************
 *
 * Func: addSite
 * Add or subtract one site sequence to a motif block. start points to
 * the starting position of the site, and w is the the width. add tells
 * whether adding or subtracting the site.
 *
 * Checked 3/20/01
 ******************************/
void addSite(struct mtfPt **blk, int *start, int w, int add) {
  int i;
  
  for (i = 0; i < w; i++, start++)
    blk[i][*start].ct += add;
}

/******************************
 *
 * Func: updateB1
 * Threshold sampler update for one block motif in both directions
 *
 * Checked 3/16/01
 ******************************/
void updateB1(struct inputParam *input, struct motif *mtf, struct
	      alignment *align, struct sequence *seq, int opt, struct
	      sequence *seqs) {
  int i, p1, p2, ct = 0, itmp, maxp, *seqi = seq->seqi, *rseqi = seq->rseqi;
  int w1 = input->w1, len = seq->len, end1 = len - w1, dir = FORWARD;
  double random, sum = 0, maxs = -1000;
  double th1 = input->seqth1, th2 = input->seqth2, *pb1 = mtf->pb1, 
    *pb2 = mtf->pb2, *bgscore = seq->bgscore1, *rbgscore = seq->rbgscore1;

  /* take out the site contributed by this sequence */
  takeInOut(seq, mtf, align, input, BLK1, SUBTRACT);
  align->copy = 0;

  /* recalculate motif matrix */
  calcBlkLog(mtf->blk1, mtf->site, w1, input->scount);

  /* score all possible sites in both directions, keeping track of the
   * best site */
  for (i = 0; i <= end1; i++) {
    pb1[i] = exp(seqMtfScore(&(seqi[i]), w1, mtf->blk1) - bgscore[i]);
    if (pb1[i] > maxs) {
      maxs = pb1[i];
      maxp = i;
      dir = FORWARD;
    }
    pb2[i] = exp(seqMtfScore(&(rseqi[i]), w1, mtf->blk1) - rbgscore[i]);
    if (pb2[i] > maxs) {
      maxs = pb2[i];
      maxp = i;
      dir = BACKWARD;
    }
  }

  /* add this alignment (either max > th1 or > th2 or all in refine step) to
   * sequence and get rid of all overlapping regions in both directions */
  while (maxs > th1 || 
	 (opt == BEST && (maxs > th2 || (input->all && ct == 0)))) {
    addcopy(maxp, 0, align, dir, 0);
    ct++;

    /* get rid of all overlapping regions */
    p1 = maxp - w1 + 1 > 0 ? maxp - w1 + 1 : 0;
    p2 = maxp + w1;
    if (dir == FORWARD) {
      for (i = p1; i < p2; i++) {
	pb1[i] = 0;
	itmp = end1 - i;
	if (itmp > 0) 
	  pb2[itmp] = 0;
      }
    } else {
      for (i = p1; i < p2; i++) {
	pb2[i] = 0; 
	itmp = end1 - i;
	if (itmp > 0)
	  pb1[itmp] = 0;
      }
    }

    /* find the next best site */
    maxs = -1000;
    for (i = 0; i <= end1; i++) {
      if (pb1[i] > maxs) {
	maxs = pb1[i];
	maxp = i;
	dir = FORWARD;
      }
      if (pb2[i] > maxs) {
	maxs = pb2[i];
	maxp = i;
	dir = BACKWARD;
      }
    }
  }
  
  /* exit point for BEST */
  if (opt == BEST) {
    takeInOut(seq, mtf, align, input, BLK1, ADD);
    return;
  }

  /* sample one copy among all sites with score between [th2, th1] */
  if (maxs < th2) {
    if (ct == 0 && (input->all || drand() < input->percent)) {
      for (i = 0; i <= end1; i++) {
	pb1[i] = sum += pb1[i];
	pb2[i] = sum += pb2[i];
      }
    }
  } else {
    for (i = 0; i <= end1; i++) {
      if (pb1[i] > th2) {
	pb1[i] = sum += pb1[i];
      } else {
	pb1[i] = sum;
      }
      if (pb2[i] > th2) {
	pb2[i] = sum += pb2[i];
      } else {
	pb2[i] = sum;
      }
    }
  }
  if (sum > 0) {
    random = drand() * sum;
    for (i = 0; i <= end1; i++) {
      if (pb1[i] >= random) {
	maxp = i;
	dir = FORWARD;
	break;
      } else if (pb2[i] >= random) {
	maxp = i;
	dir = BACKWARD;
	break;
      }
    }
    addcopy(maxp, 0, align, dir, 0);
  }

  /* add all new aligned sites to the motif block1 */
  takeInOut(seq, mtf, align, input, BLK1, ADD);
}

/******************************
 *
 * Func: updateS1
 * Threshold sampler update for one block motif in forward directions
 *
 * Checked 3/16/01
 ******************************/
void updateS1(struct inputParam *input, struct motif *mtf, struct
	      alignment *align, struct sequence *seq, int opt, struct
	      sequence *seqs) {
  int i, p1, p2, ct = 0, maxp, *seqi = seq->seqi;
  int w1 = input->w1, end1 = seq->len - w1;
  double random, sum = 0, maxs = -1000;
  double th1 = input->seqth1, th2 = input->seqth2, *pb1 = mtf->pb1,
    *bgscore = seq->bgscore1;
  
  /* take out sites contributed by this sequence and recalculate the
   * matrix */
  takeInOut(seq, mtf, align, input, BLK1, SUBTRACT);
  align->copy = 0;
  calcBlkLog(mtf->blk1, mtf->site, w1, input->scount);
  
  /* score all possible sites, keeping track of the best site */
  for (i = 0; i <= end1; i++) {
    pb1[i] = exp(seqMtfScore(&(seqi[i]), w1, mtf->blk1) - bgscore[i]);
    if (pb1[i] > maxs) {
      maxs = pb1[i];
      maxp = i;
    }
  }
  
  /* add this alignment (either max > th1 or > th2 or all in refine step) to
   * sequence and get rid of all overlapping regions in both directions */
  while (maxs > th1 || 
	 (opt == BEST && (maxs > th2 || (input->all && ct == 0)))) {
    addcopy(maxp, 0, align, FORWARD, 0);
    ct++;
    
    /* get rid of all overlapping regions */
    p1 = maxp - w1 + 1 > 0 ? maxp - w1 + 1 : 0;
    p2 = maxp + w1;
    for (i = p1; i < p2; i++)
      pb1[i] = 0;

    /* find the next best site */
    maxs = -1000;
    for (i = 0; i <= end1; i++) {
      if (pb1[i] > maxs) {
	maxs = pb1[i];
	maxp = i;
      }
    }
  }
  
  /* exit point for BEST */
  if (opt == BEST) {
    takeInOut(seq, mtf, align, input, BLK1, ADD);
    return;
  }

  /* sample one copy among all sites with score between [th2, th1] */
  if (maxs < th2) {
    if (ct == 0 && (input->all || drand() < input->percent)) {
      for (i = 0; i <= end1; i++) {
	pb1[i] = sum += pb1[i];
      }
    }
  } else {
    for (i = 0; i <= end1; i++) {
      if (pb1[i] > th2) {
	pb1[i] = sum += pb1[i];
      } else {
	pb1[i] = sum;
      }
    }
  }
  if (sum > 0) {
    random = drand() * sum;
    for (i = 0; i <= end1; i++) {
      if (pb1[i] >= random) {
	maxp = i;
	break;
      } 
    }
    addcopy(maxp, 0, align, FORWARD, 0);
  }

  /* add all new aligned sites to the motif block1 */
  takeInOut(seq, mtf, align, input, BLK1, ADD);
}

/******************************
 *
 * Func: updateB2
 * Threshold sampler update for two block motif in both directions
 *
 * Checked 2/26/2001
 ******************************/
void updateB2(struct inputParam *input, struct motif *mtf, struct
	      alignment *align, struct sequence *seq, int opt, struct
	      sequence *seqs) {
  int i, j, st1, st2, itmp, p1, p2, ct = 0, *seqi = seq->seqi, *rseqi =
    seq->rseqi;
  int grange = input->grange, w1 = input->w1, w2 = input->w2, maxp, w3 =
  w1 + input->gl, len = seq->len, wid = w1 + input->gl + w2;  
  int end1 = len - wid, end2 = seq->len - w2, dir = FORWARD;
  double random = 0, sum = 0, maxs = -1000, s1, s2;
  double th1 = input->seqth1, th2 = input->seqth2, dtmp;
  double *pb1 = mtf->pb1, *pb12 = mtf->pb12, *pb22 = mtf->pb22, *pb2
    = mtf->pb2, **mat1 = mtf->mat1, **mat2 = mtf->mat2, score1, score2;
  double *bgscore1 = seq->bgscore1, *rbgscore1 = seq->rbgscore1, 
    *bgscore2 = seq->bgscore2, *rbgscore2 = seq->rbgscore2;
  
  /* take out the sites contributed by this sequence */
  takeInOut(seq, mtf, align, input, BLKBOTH, SUBTRACT);
  align->copy = 0;
  
  /* recalculate motif matrix */
  calcBlkLog(mtf->blk1, mtf->site, w1, input->scount);
  calcBlkLog(mtf->blk2, mtf->site, w2, input->scount);

  /* score all possible sites with motif block2 */
  for (i = w3; i <= end2; i++) {
    pb12[i] = seqMtfScore(&(seqi[i]), w2, mtf->blk2) - bgscore2[i];
    pb22[i] = seqMtfScore(&(rseqi[i]), w2, mtf->blk2) - rbgscore2[i];
  }
  
  /* marginal distribution of the first site, and keep the max */ 
  for (i = 0; i <= end1; i++) {
    score1 = seqMtfScore(&(seqi[i]), w1, mtf->blk1) - bgscore1[i];
    score2 = seqMtfScore(&(rseqi[i]), w1, mtf->blk1) - rbgscore1[i];
    s1 = s2 = 0;

    /* s1 = sum{score(site X) * score(site Y) for givin X} */
    itmp = i + w3;
    for (j = 0; j < grange; j++) {
      if (itmp > end2) {
	mat1[i][j] = s1;
	mat2[i][j] = s2;
      } else {
	mat1[i][j] = s1 += exp(score1 + pb12[itmp]);
	mat2[i][j] = s2 += exp(score2 + pb22[itmp]);
      }
      itmp++;
    }

    /* mat*[i][grange] stores the marginal for X starting at i */
    mat1[i][grange] = s1;
    if (s1 > maxs) {
      maxs = s1;
      maxp = i;
      dir = FORWARD;
    }
    mat2[i][grange] = s2;
    if (s2 > maxs) {
      maxs = s2;
      maxp = i;
      dir = BACKWARD;
    }
  }

  /* add this alignment (either max > th1 or > th2 or all in refine step) to
   * sequence and get rid of all overlapping regions in both directions */
  while (maxs > th1 || 
	 (opt == BEST && (maxs > th2 || (input->all && ct == 0)))) {
    st1 = maxp;

    if (opt == BEST) {
      
      if (dir == FORWARD) {
	maxs = mat1[st1][0];
	st2 = st1 + w3;
	for (i = 1; i < grange; i++) {
	  dtmp = mat1[st1][i] - mat1[st1][i-1];
	  if (dtmp > maxs) {
	    maxs = dtmp;
	    st2 = st1 + i + w3;
	  }
	}
      } else {
	maxs = mat2[st1][0];
	st2 = st1 + w3;
	for (i = 1; i < grange; i++) {
	  dtmp = mat2[st1][i] - mat2[st1][i-1];
	  if (dtmp > maxs) {
	    maxs = dtmp;
	    st2 = st1 + i + w3;
	  }
	}
      }
    } else {

      /* sample the second start position */
      if (dir == FORWARD) {
	random = drand() * maxs;
	for (i = 0; i < grange; i++) {
	  if (mat1[st1][i] >= random) {
	    st2 = st1 + w3 + i;
	    break;
	  }
	} 
      } else {
	random = drand() * maxs;
	for (i = 0; i < grange; i++) {
	  if (mat2[st1][i] >= random) {
	    st2 = st1 + w3 + i;
	    break;
	  }
	}
      }
    }
    addcopy(st1, st2, align, dir, w2);
    ct++;

    /* get rid of overlaps */
    p1 = st1 - wid + 1 > 0 ? st1 - wid + 1 : 0;
    p2 = st2 + w2;
    if (dir == FORWARD) {
      for (i = p1; i < p2; i++) {
	mat1[i][grange] = 0;
	itmp = end1 - i;
	if (itmp > 0)
	  mat2[itmp][grange] = 0;
      }
    } else {
      for (i = p1; i < p2; i++) {
	mat2[i][grange] = 0;
	itmp = end1 - i;
	if (itmp > 0)
	  mat1[itmp][grange] = 0;
      }
    }

    /* find the next best site pair */
    maxs = -1000;
    for (i = 0; i <= end1; i++) {
      if (mat1[i][grange] > maxs) {
	maxs = mat1[i][grange];
	maxp = i;
	dir = FORWARD;
      }
      if (mat2[i][grange] > maxs) {
	maxs = mat2[i][grange];
	maxp = i;
	dir = BACKWARD;
      }
    }
  }
  
  /* exit point for BEST */
  if (opt == BEST) {
    takeInOut(seq, mtf, align, input, BLKBOTH, ADD);
    return;
  }

  /* sample one copy among all sites with score between [th2, th1] */
  if (maxs < th2) {
    if (ct == 0 && (input->all || drand() < input->percent)) {
      for (i = 0; i <= end1; i++) {
	pb1[i] = sum += mat1[i][grange];
	pb2[i] = sum += mat2[i][grange];
      }
    }
  } else {
    for (i = 0; i <= end1; i++) {
      if (mat1[i][grange] > th2) {
	pb1[i] = sum += mat1[i][grange];
      } else {
	pb1[i] = sum;
      }
      if (mat2[i][grange] > th2) {
	pb2[i] = sum += mat2[i][grange];
      } else {
	pb2[i] = sum;
      } 
    }
  }
  if (sum > 0) {
    random = drand() * sum;
    for (i = 0; i <= end1; i++) {
      if (pb1[i] >= random) {
	st1 = i;

	/* sample second position */
	random = drand() * mat1[i][grange];
	for (j = 0; j < grange; j++) {
	  if (mat1[i][j] >= random) {
	    st2 = st1 + w3 + j;
	    dir = FORWARD;
	    break;
	  }
	}
	break;
      } else if (pb2[i] > random) {
	st1 = i;

	/* sample second position */
	random = drand() * mat2[i][grange];
	for (j = 0; j < grange; j++) {
	  if (mat2[i][j] >= random) {
	    st2 = st1 + w3 + j;
	    dir = BACKWARD;
	    break;
	  }
	}
	break;
      }
    }
    addcopy(st1, st2, align, dir, w2);
  }

  /* add all new aligned sites to the motif blocks */
  takeInOut(seq, mtf, align, input, BLKBOTH, ADD);
}

/******************************
 *
 * Func: updateS2
 * Threshold sampler update for two block motif in forward directions
 *
 * Checked 3/16/2001
 ******************************/
void updateS2(struct inputParam *input, struct motif *mtf, struct
	      alignment *align, struct sequence *seq, int opt, struct
	      sequence *seqs) {
  int i, j, st1, st2, itmp, p1, p2, ct = 0, *seqi = seq->seqi;
  int grange = input->grange, w1 = input->w1, w2 = input->w2, maxp = 0, 
    w3 = w1 + input->gl, len = seq->len, wid = w1 + input->gl + w2, 
    end1 = len - wid, end2 = len - w2;
  double random, sum = 0, maxs = -1000, s1, score1;
  double th1 = input->seqth1, th2 = input->seqth2, dtmp;
  double *pb12 = mtf->pb12, *pb1 = mtf->pb1, **mat1 = mtf->mat1;
  double *bgscore1 = seq->bgscore1, *bgscore2 = seq->bgscore2;
  
  /* take out the site contributed by this sequence */
  takeInOut(seq, mtf, align, input, BLKBOTH, SUBTRACT);
  align->copy = 0;
  
  /* recalculate motif matrix */
  calcBlkLog(mtf->blk1, mtf->site, w1, input->scount);
  calcBlkLog(mtf->blk2, mtf->site, w2, input->scount);

  /* score all possible sites with second motif block */
  for (i = w3; i <= end2; i++)
    pb12[i] = seqMtfScore(&(seqi[i]), w2, mtf->blk2) - bgscore2[i];

  /* marginal distribution of the first site, and keep the max */ 
  for (i = 0; i <= end1; i++) {
    score1 = seqMtfScore(&(seqi[i]), w1, mtf->blk1) - bgscore1[i];
    s1 = 0;
    
    /* s1 = sum{score(site X) * score(site Y) for givin X} */
    itmp = i + w3;
    for (j = 0; j < grange; j++) {
      if (itmp > end2)
	mat1[i][j] = s1;
      else
	mat1[i][j] = s1 += exp(score1 + pb12[itmp]);
      itmp++;
    }
    
    mat1[i][grange] = s1;
    if (s1 > maxs) {
      maxs = s1;
      maxp = i;
    }
  }
  
  /* add this alignment (either max > th1 or > th2 or all in refine step) to
   * sequence and get rid of all overlapping regions in both directions */
  while (maxs > th1 || 
	 (opt == BEST && (maxs > th2 || (input->all && ct == 0)))) {
      st1 = maxp;
    
    if (opt == BEST) {

      /* use the best second site */
      maxs = mat1[st1][0];
      st2 = st1 + w3;
      for (i = 1; i < grange; i++) {
	dtmp = mat1[st1][i] - mat1[st1][i-1];
	if (dtmp > maxs) {
	  maxs = dtmp;
	  st2 = st1 + i + w3;
	}
      }
    } else {
      /* sample the second start position */
      random = drand() * maxs;
      for (i = 0; i < grange; i++) {
	if (mat1[st1][i] >= random) {
	  st2 = st1 + w3 + i;
	  break;
	}
      } 
    }
    addcopy(st1, st2, align, FORWARD, w2);
    ct++;

    /* get rid of all overlapping regions */
    p1 = st1 - wid + 1 > 0 ? st1 - wid + 1 : 0;
    p2 = st2 + w2;
    for (i = p1; i < p2; i++)
      mat1[i][grange] = 0;

    /* find the next best site */
    maxs = -1000;
    for (i = 0; i <= end1; i++) {
      if (mat1[i][grange] > maxs) {
	maxs = mat1[i][grange];
	maxp = i;
      }
    }
  }
  
  /* exit point for BEST */
  if (opt == BEST) {
    takeInOut(seq, mtf, align, input, BLKBOTH, ADD);
    return;
  }

  /* sample one copy among all sites with score between [th2, th1] */
  if (maxs < th2) {
    if (ct == 0 && (input->all || drand() < input->percent)) {
      for (i = 0; i <= end1; i++)
	pb1[i] = sum += mat1[i][grange];
    }
  } else {
    for (i = 0; i <= end1; i++) {
      if (mat1[i][grange] > th2) {
	pb1[i] = sum += mat1[i][grange];
      } else {
	pb1[i] = sum;
      }
    }
  }
  if (sum > 0) {
    random = drand() * sum;
    for (i = 0; i <= end1; i++) {
      if (pb1[i] >= random) {
	st1 = i;

	/* sample second position */
	random = drand() * mat1[i][grange];
	for (j = 0; j < grange; j++) {
	  if (mat1[i][j] >= random) {
	    st2 = st1 + w3 + j;
	    break;
	  }
	}
	break;
      }
    }
    addcopy(st1, st2, align, FORWARD, w2);
  }

  /* add all new aligned sites to the motif blocks */
  takeInOut(seq, mtf, align, input, BLKBOTH, ADD);
}

/******************************
 *
 * Func: updateP
 * Threshold sampler update for two block palindrome motifs
 *
 * Checked 2/26/2001
 ******************************/
void updateP(struct inputParam *input, struct motif *mtf, struct
	     alignment *align, struct sequence *seq, int opt, struct
	     sequence *seqs) {
  int i, j, st1, st2, itmp, p1, p2, ct = 0, *seqi = seq->seqi, *rseqi =
    seq->rseqi; 
  int grange = input->grange, w1 = input->w1, wid = 2 * w1 + input->gl,
    len = seq->len, end1 = len - wid, maxp; 
  double random, sum = 0, maxs = -1000, s1, score1;
  double th1 = input->seqth1, th2 = input->seqth2, dtmp;
  double *pb12 = mtf->pb12, *pb1 = mtf->pb1, **mat1 = mtf->mat1;
  double *bgscore1 = seq->bgscore1, *bgscore2 = seq->bgscore2;
  
  /* take out the site contributed by this sequence and recalculate
   * motif matrix */
  takeInOut(seq, mtf, align, input, BLK1, SUBTRACT);
  align->copy = 0;
  calcBlkLog(mtf->blk1, mtf->site, w1, input->scount);
  
  /* score all possible second sites */
  for (i = 0; i <= end1; i++) {
    pb12[i] = seqMtfScore(&(rseqi[i]), w1, mtf->blk1) - bgscore2[i];
  }
  
  /* marginal distribution of the first site, and keep the max */ 
  for (i = 0; i <= end1; i++) {
    score1 = seqMtfScore(&(seqi[i]), w1, mtf->blk1) - bgscore1[i];
    s1 = 0;
    itmp = end1 - i;
    for (j = 0; j < grange; j++) {
      if (itmp > 0)
	mat1[i][j] = s1 += exp(score1 + pb12[itmp]);
      else
	mat1[i][j] = s1;
      itmp--; 
    }
    mat1[i][grange] = s1;
    if (s1 > maxs) {
      maxs = s1;
      maxp = i;
    }
  }
  
  /* add this alignment (either max > th1 or > th2 or all in refine step) to
   * sequence and get rid of all overlapping regions in both directions */
  while (maxs > th1 || 
	 (opt == BEST && (maxs > th2 || (input->all && ct == 0)))) {
    st1 = maxp;
      
    if (opt == BEST) {

      /* use the best second site */
      maxs = mat1[st1][0];
      st2 = len - st1 -wid;
      for (i = 1; i < grange; i++) {
	dtmp = mat1[st1][i] - mat1[st1][i-1];
	if (dtmp > maxs) {
	  maxs = dtmp;
	  st2 = len - st1 - wid - i;
	}
      }
    } else {

    /* sample the second start position */
      random = drand() * maxs;
      for (i = 0; i < grange; i++) {
	if (mat1[st1][i] >= random) {
	  st2 = len - st1 - wid - i;
	  break;
	}
      }
    } 
    addcopy(st1, st2, align, PLDM, w1);
    ct++;
    
    /* get rid of all overlapping regions */
    p1 = st1 - wid + 1 > 0 ? st1 - wid + 1 : 0;
    p2 = len - st2;
    for (i = p1; i < p2; i++)
      mat1[i][grange] = 0;
    
    /* find the next best site */
    maxs = -1000;
    for (i = 0; i <= end1; i++) {
      if (mat1[i][grange] > maxs) {
	maxs = mat1[i][grange];
	maxp = i;
      }
    }
  }
  
  /* exit point for BEST */
  if (opt == BEST) {
    takeInOut(seq, mtf, align, input, BLK1, ADD);
    return;
  }
  
  /* sample one copy among all sites with score between [th2, th1] */
  if (maxs < th2) {
    if (ct == 0 && (input->all || drand() < input->percent)) {
      for (i = 0; i <= end1; i++)
	pb1[i] = sum += mat1[i][grange];
    }
  } else {
    for (i = 0; i <= end1; i++) {
      if (mat1[i][grange] > th2) {
	pb1[i] = sum += mat1[i][grange];
      } else {
	pb1[i] = sum;
      }
    }
  }
  if (sum > 0) {
    random = drand() * sum;
    for (i = 0; i <= end1; i++) {
      if (pb1[i] >= random) {
	st1 = i;

	/* sample second position */
	random = drand() * mat1[i][grange];
	for (j = 0; j < grange; j++) {
	  if (mat1[i][j] >= random) {
	    st2 = len - st1 - wid - j;
	    break;
	  }
	}
	break;
      }
    }
    addcopy(st1, st2, align, PLDM, w1);
  }
  
  /* add all new aligned sites to the motif blocks */
  takeInOut(seq, mtf, align, input, BLK1, ADD);
}

/******************************
 *
 * Func: calcBlkLog
 * Convert motif distribution from count to log so scoring is done faster
 *
 * Checked 3/16/01
 ******************************/
void calcBlkLog(struct mtfPt **blk, int site, int wid, double *sct) {
  int i, j; 

  for (i = 0; i < wid; i++) {
    for (j = 0; j < ALPHASIZE; j++)
      blk[i][j].log = log((blk[i][j].ct + sct[j]) / (site + sct[ALPHASIZE]));
  }
}

/******************************
 *
 * Func: seqMtfScore
 * Score the site of wid long starting from *seqi by blk
 *
 * Checked 3/29/01
 ******************************/
double seqMtfScore(int *seqi, int wid, struct mtfPt **blk) {
  double score = 0;
  int i;
  
  for (i = 0; i < wid; i++)
    score += blk[i][seqi[i]].log;
  
  return score;
}

/******************************
 *
 * Func: addcopy
 * Add one copy of alignment to a sequence, dir tells whether the
 * alignment is in forward, backward or pldm, and blk2 tells whether the
 * second starting position is real or w2 = 0
 *
 * Checked 3/16/01
 ******************************/
void addcopy(int st1, int st2, struct alignment *align, int dir, int blk2) {
  int size;
  
  /* if the new copy will exceed the array size, reallocate start
   * position-related arrays to double the size */
  if (align->copy == align->size) {
    align->size *= 2;
    size = align->size * INTSIZE;
    align->st1 = (int *) realloc (align->st1, size);
    align->dir = (int *) realloc (align->dir, size);
    if (blk2)
      align->st2 = (int *) realloc (align->st2, size);
  }

  /* add the new alignment to the sequence */
  align->dir[align->copy] = dir;
  align->st1[align->copy] = st1;
  align->dir[align->copy] = dir;
  if (blk2)
    align->st2[align->copy] = st2;
  align->copy++;
}

/******************************
 *
 * Func: colShift
 * bkmtf1 shift both blks to the left by 1, and bkmtf2 shift both blks to
 * the right by 1. The compare the score blk by blk (ignore the items that
 * are the same for each blk).
 *
 * Checked 4/3/01
 ******************************/
void colShift(struct inputParam *input, struct sequence *seqs, struct
	      motif *mtf, struct motif *bkmtf1, struct motif *bkmtf2) {
  int w1 = input->w1, w2 = input->w2, i, shift1 = 0, shift2 = 0;
  double info, max1 = -1000, max2 = -1000;
  struct motif *mtfptr[] = {bkmtf1, mtf, bkmtf2};

  mtfCopy(bkmtf1, mtf, input);
  blkColShift(bkmtf1, -1, -1, seqs, input);
  mtfCopy(bkmtf2, mtf, input);
  blkColShift(bkmtf2, 1, 1, seqs, input);

  for (i = 0; i < 3; i++) {
    info = blkMtfScore(mtfptr[i]->blk1, w1)
      - blkBgScore(seqs, mtfptr[i]->align, BLK1);
    if (max1 < info) {
      max1 = info;
      shift1 = i - 1;
    }
    if (w2 && !input->pldm) {
      info = blkMtfScore(mtfptr[i]->blk2, w2)
	- blkBgScore(seqs, mtfptr[i]->align, BLK2);
      if (max1 < info) {
	max2 = info;
	shift2 = i - 1;
      }
    }
  }
  
  blkColShift(mtf, shift1, shift2, seqs, input);
}
 
/******************************
 *
 * Func: blkColShift
 * Move the 
 * 
 * Checked 4/3/01
 ******************************/
void blkColShift(struct motif *mtf, int shift1, int shift2, struct sequence
		 *seqs, struct inputParam *input) {
  int w1 = input->w1, w2 = input->w2, pldm = input->pldm;
  struct sequence *s;
  struct alignment *a;
  int i, tmp;

  clearMotif(mtf, input);
  for (s = seqs, a = mtf->align; s; s = s->next, a = a->next) {
    for (i = 0; i < a->copy; i++) {
      tmp = a->st1[i] + shift1;
      if (tmp >= 0 && tmp + w1 <= s->len) {
	a->st1[i] = tmp;
      }
      if (w2) {
	if (pldm) {
	  tmp = a->st2[i] + shift1;
	  if (tmp >= 0 && tmp + w1 <= s->len) {
	    a->st2[i] = tmp;
	  }
	} else {
	  tmp = a->st2[i] + shift2;
	  if (tmp >= 0 && tmp + w2 <= s->len) {
	    a->st2[i] = tmp;
	  }
	}
      }
    }
    takeInOut(s, mtf, a, input, BLKBOTH, ADD);
  }
}

/******************************
 *
 * Func: blkShift
 * During two block motif search, see whether leftalign and rightalign
 * produce a motif with better information content.
 * 
 * mtf                      |--W1--|.....|--W2---|
 * bkmtf1      |--W1--|.....|--W2---|  
 * bkmtf2                                |--W1--|.....|--W2---|     
 *
 * Checked 4/3/01
 ******************************/
void blkShift(struct inputParam *input, struct sequence *seqs, struct freq
	      *bglog, struct motif *mtf, struct motif *bkmtf1, struct
	      motif *bkmtf2) {
  /* maxp = 0 no shift, maxp = -1 left shift, maxp = 1 right shift */
  double info, maxs;
  int i, maxp = 0;
  struct motif *mtfptr;
  struct sequence *s;
  struct alignment *a;

  /* initialize the left shifted and right shifted motif */
  mtfCopy(bkmtf2, mtf, input);
  blkRightShift(bkmtf2, seqs, input);
  mtfCopy(bkmtf1, mtf, input);
  blkLeftShift(bkmtf1, seqs, input);

  /* score for the unshifted motif */
  maxs = motifInfo(mtf, input, seqs);

  /* score for the left shifted motif */
  for (i = 0; i < 3; i++) {
    for (s = seqs, a = bkmtf1->align; s; s = s->next, a = a->next)
      update(input, bkmtf1, a, s, REG, seqs);
  }
  info = motifInfo(bkmtf1, input, seqs);

  /* if leftshift is the best motif */
  if (maxs < info) {
    maxs = info;
    maxp = -1;
  }

  /* socre for the right shifted motif */
  for (i = 0; i < 3; i++) {
    for (s = seqs, a = bkmtf2->align; s; s = s->next, a = a->next)
      update(input, bkmtf2, a, s, REG, seqs);
  }
  info = motifInfo(bkmtf2, input, seqs);

  /* if rightshift is the best motif */
  if (maxs < info) {
    mtfptr = bkmtf2;
    bkmtf2 = mtf;
    mtf = mtfptr;
  } else {

    /* if left shift is the best motif */
    if (maxp == -1) {
      mtfptr = bkmtf1;
      bkmtf1 = mtf;
      mtf = mtfptr;
      
    } /* else current mtf alignment is the best */
  }
}

/******************************
 *
 * Func: blkLeftShift
 *
 * Checked 4/3/01
 ******************************/
void blkLeftShift(struct motif *mtf, struct sequence *seqs, struct
		  inputParam *input) {
  struct sequence *s;
  struct alignment *a;
  int move = input->meangap + input->w1, i, tmp;
  
  clearMotif(mtf, input);

  for (s = seqs, a = mtf->align; s; s = s->next, a = a->next) {
    for (i = 0; i < a->copy; i++) {
      tmp = a->st1[i];
      a->st1[i] = tmp > move ? tmp - move : 0;
      if (input->pldm) 
	a->st2[i] = s->len - tmp - input->w1;
      else
	a->st2[i] = tmp;
    }
    takeInOut(s, mtf, a, input, BLKBOTH, ADD);
  }
}

/******************************
 *
 * Func: blkRightShift
 *
 * Checked 4/3/01
 ******************************/
void blkRightShift(struct motif *mtf, struct sequence *seqs, struct
		   inputParam *input) {
  struct sequence *s;
  struct alignment *a;
  int move = input->meangap + input->w1, w1 = input->w1, w2 = input->w2,
    len, i, tmp;
  
  clearMotif(mtf, input);

  for (s = seqs, a = mtf->align; s; s = s->next, a = a->next) {
    len = s->len;
    for (i = 0; i < a->copy; i++) {
      tmp = a->st2[i];
      if (input->pldm) {
	a->st2[i] = tmp > move ? tmp - move : 0;
	a->st2[i] = len - tmp - w1;
      } else {
	a->st1[i] = tmp + w1 > len ? len - w1 : tmp;
	a->st2[i] = tmp + move + w2 > len ? len - w2 : tmp + move;
      }
    }
    takeInOut(s, mtf, a, input, BLKBOTH, ADD);
  }
}
 
/******************************
 *
 * Func: motifInfo
 * Motif info is the site * exp(average relative entropy per column)
 *
 * Checked 3/20/01
 ******************************/
double motifInfo(struct motif *mtf, struct inputParam *input, struct
		 sequence *seqs) { 
  double info = 0;
  int w1 = input->w1, w2 = input->w2;
  
  if (mtf->site == 0) {
    mtf->mtfscore = 0;
    return 0;
  }
  
  calcBlkLog(mtf->blk1, mtf->site, w1, input->scount);
  info += blkMtfScore(mtf->blk1, w1) -
    blkBgScore(seqs, mtf->align, BLK1);
  if (w2 && !input->pldm) {
    calcBlkLog(mtf->blk2, mtf->site, w2, input->scount);
    info += blkMtfScore(mtf->blk2, w2) -
      blkBgScore(seqs, mtf->align, BLK2);
    if (input->expect)
      mtf->mtfscore = (info - mtf->site * input->expect) / (w1 + w2);
    else 
      mtf->mtfscore = log(mtf->site) * info / mtf->site / (w1 + w2);
    
  } else {
    if (input->expect)
      mtf->mtfscore = (info - mtf->site * input->expect) / w1;
    else 
      mtf->mtfscore = log(mtf->site) * info / mtf->site / w1;
  }

  /* Note: Bayes mtfscore here could be negative!! */
  return mtf->mtfscore;
}
 
/******************************
 *
 * Func: blkMtfScore
 * Calculates the gamma motif matrix information
 * 
 * Checked 5/7/2001
 ******************************/
double blkMtfScore(struct mtfPt **blk, int wid) {
  int i, j;
  double info = 0;

  for (i = 0; i < wid; i++) {
    for (j = 0; j < ALPHASIZE; j++)
      info += blk[i][j].ct * blk[i][j].log;
  }

  return info;
}

/******************************
 *
 * Func: blkBgScore
 * Calculates motif matrix information from the background
 * 
 * Checked 3/29/2001
 ******************************/
double blkBgScore(struct sequence *seqs, struct alignment *align, int id) {
   struct sequence *s;
   struct alignment *a;
   int i;
   double info = 0;
   
   for (s = seqs, a = align; s; s = s->next, a = a->next) {
    for (i = 0; i < a->copy; i++) {
      if (id == BLK1) {

	switch(a->dir[i]) {
	case FORWARD: 
	  info += s->bgscore1[a->st1[i]];
	  break;
	case PLDM: 
	  info += s->bgscore1[a->st1[i]] + s->bgscore2[a->st2[i]];
	  break;
	case BACKWARD: 
	  info += s->rbgscore1[a->st1[i]];
	  break;
	}

      } else {
	/* conditioned that only comes here if w2 && !pldm */

	switch(a->dir[i]) {
	case FORWARD: 
	  info += s->bgscore2[a->st2[i]];
	  break;
	case BACKWARD: 
	  info += s->rbgscore2[a->st2[i]];
	  break;
	}
      }
    }
  }
  
  return info;
}

/******************************
 *
 * Func: findMotif
 * Find motif from input sequences, keep track of the best motifs found
 * during iterations and tries, only reporting the highest scoring motifs.
 *
 * Checked 4/3/01
 ******************************/
void findMotif(struct inputParam *input, struct freq *bglog, struct motif 
		*mtf, struct motif *bkmtf1, struct motif *bkmtf2, struct 
		motif *bestmtf, struct sequence	*seqs, struct motif **results) 
{
  struct sequence *s;
  struct alignment *a;
  int i = 0, j, result = input->result;
  int w1 = input->w1, w2 = input->w2, iteration = input->iteration,
    w3 = input->w2 && !input->pldm ? YES : NO;
  char con1[w1+1], rcon1[w1+1], con2[w3+1], rcon2[w3+1];
  struct motif *mtfptr;
  
  while (i < input->tries) {
    i++;
    
    /* initialize parameters and the motif by random alignment */
    clearMotif(mtf, input);
    bestmtf->mtfscore = 0;
    randomStart(input, seqs, mtf);

    /* trial updates */
    for (j = 0; j < TRIAL; j++) {
      for (s = seqs, a = mtf->align; s; s = s->next, a = a->next)
	update(input, mtf, a, s, REG, seqs);
    }   

    /* real updates with th2 increasing each iteration */
    for (j = 0; j < iteration; j++) {
      input->seqth2 += input->inc;
      for (s = seqs, a = mtf->align; s; s = s->next, a = a->next)
	update(input, mtf, a, s, REG, seqs);
      if (j && !(j % 10)) { 
	colShift(input, seqs, mtf, bkmtf1, bkmtf2);
	if (w2) blkShift(input, seqs, bglog, mtf, bkmtf1, bkmtf2);
      }
      if (motifInfo(mtf, input, seqs) > bestmtf->mtfscore)
	mtfCopy(bestmtf, mtf, input);
    }

    if (input->local) {
      if (bestmtf->mtfscore == 0) {
	fprintf(input->ofp, "Try #%d\t0\tNo motif found\n", i);
	fflush(input->ofp);
	continue;
      } else {
	getConsensus(w1, bestmtf->blk1, con1, rcon1);
	if (w3)
	  getConsensus(w2, bestmtf->blk2, con2, rcon2);
	else 
	  getConsensus(0, bestmtf->blk1, con2, rcon2);
	fprintf(input->ofp, "Try #%d\t%.3f\t%s %s\t%s %s\t%d\n",
		i, bestmtf->mtfscore, con1, rcon1, con2, rcon2, 
		bestmtf->site);
	fflush(input->ofp);
      }
    }

    /* refine the best motif */
    if (input->refine) {
      mtfCopy(mtf, bestmtf, input);
      for (s = seqs, a = mtf->align; s; s = s->next, a = a->next)
	update(input, mtf, a, s, BEST, seqs);
      if (motifInfo(mtf, input, seqs) > bestmtf->mtfscore)
	mtfCopy(bestmtf, mtf, input);
      
      if (input->local) {
	if (bestmtf->mtfscore == 0) {
	  fprintf(input->ofp, "Try #%d\t0\tNo motif found\n", i);
	  fflush(input->ofp);
	  continue;
	} else {
	  getConsensus(w1, bestmtf->blk1, con1, rcon1);
	  if (w3)
	    getConsensus(w2, bestmtf->blk2, con2, rcon2);
	  else 
	    getConsensus(0, bestmtf->blk1, con2, rcon2);
	  fprintf(input->ofp, "Ref #%d\t%.3f\t%s %s\t%s %s\t%d\n",
		  i, bestmtf->mtfscore, con1, rcon1, con2, rcon2, 
		  bestmtf->site);
	  fflush(input->ofp);
	}
      }
    }
    
    /* keep the top result potential motifs */
    if (bestmtf->mtfscore > results[0]->mtfscore) {
      mtfptr = results[0];
      results[0] = bestmtf;
      bestmtf = mtfptr;
      for (j = 1; j < result; j++) {
	if (results[j]->mtfscore < results[j-1]->mtfscore) {
	  mtfptr = results[j];
	  results[j] = results[j-1];
	  results[j-1] = mtfptr;
	} else
	  break;
      }
    }
  }      
  
  /* print out results  */
  fprintf(input->ofp, "\nThe highest scoring %d motifs are:\n", result);
  for (i = 0; i < result; i++) {
    mtf = results[result-1-i];
    if (mtf->mtfscore > 0) {
      getConsensus(w1, mtf->blk1, con1, rcon1);
      if (w3) {
	getConsensus(w2, mtf->blk2, con2, rcon2);
	fprintf(input->ofp, "Motif #%d: (%s/%s, %s/%s)\n", (i+1), 	      
		con1, rcon1, con2, rcon2);
      } else
	fprintf(input->ofp, "Motif #%d: (%s/%s)\n", (i+1), 	      
		con1, rcon1);
      fprintf(input->ofp, "******************************\n");
      printMotif(mtf, input);
      printAlignment(input, seqs, mtf->align, w1, w2);
      fprintf(input->ofp, "******************************\n\n");
    } else {
      break;
    }
  }
}

/******************************
 *
 * Func: mtfCopy
 * 
 * Copy the following information from src to dst
 * site, mtfscore, base ct at each blk column
 *
 * Checked 3/30/01
 ******************************/
void mtfCopy(struct motif *dst, struct motif *src, struct inputParam
	     *input) {
  int i, blk2 = input->w2 && !(input->pldm), size;
  struct alignment *a1, *a2;

  dst->site = src->site;
  dst->mtfscore = src->mtfscore;

  for (i = 0; i < input->w1; i++)
    memcpy(dst->blk1[i], src->blk1[i], MTFCOLSIZE);
  if (blk2) {
    for (i = 0; i < input->w2; i++)
      memcpy(dst->blk2[i], src->blk2[i], MTFCOLSIZE);
  }
  
  for (a1 = src->align, a2 = dst->align; a1; a1 = a1->next, a2 = a2->next)
    {
      if (a1->copy > a2->size) {
	a2->size = a1->copy;
	size = a2->size * INTSIZE;
	a2->st1 = (int *) realloc (a2->st1, size);
	a2->dir = (int *) realloc (a2->dir, size);
	if (input->w2)
	  a2->st2 = (int *) realloc (a2->st2, size);
      }
      a2->copy = a1->copy;
      size = a1->copy * INTSIZE;
      memcpy(a2->dir, a1->dir, size);
      memcpy(a2->st1, a1->st1, size);
      if (input->w2)
	memcpy(a2->st2, a1->st2, size);
    }
}

/******************************
 *
 * Func: printMotif
 * Motifs are printed:
 *  1. motif width (first blk, second blk); raw score; siginifcance value;
 *  P-value; Number of aligned sites.
 *  2. probability matrix
 *  3. consensus, Reverse Consensus, Degenerate, Reverse degenerate.
 *
 * Checked 3/16/01
 ******************************/
void printMotif(struct motif *mtf, struct inputParam *input) {

  /* motif width, score, sites, etc */
  fprintf(input->ofp, "Width (%d, %d); Gap [%d, %d]; MotifScore %.3f; Sites %d\n", input->w1, input->w2, input->gl, input->gm, mtf->mtfscore, mtf->site);

  /* print out the matrix and consensus of the two blocks */
  calcBlkLog(mtf->blk1, mtf->site, input->w1, input->scount);
  fprintf(input->ofp, "\nBlk1    ");
  printBlk(input->w1, mtf->blk1, input->ofp);
  if (input->w2 && !input->pldm) {
    calcBlkLog(mtf->blk2, mtf->site, input->w2, input->scount);
    fprintf(input->ofp, "\nBlk2    ");
    printBlk(input->w2, mtf->blk2, input->ofp);
  }
  fprintf(input->ofp, "\n");
  fflush(input->ofp);
}

/******************************
 *
 * Func: printBlk
 * Print motif blk probability matrix, consensus and degenerate, one
 * column a line (this is the transfec format). Consensus is the highest
 * probability base, degenerate considers all bases with probaiblity >
 * 0.25 with the following rule: 
 * 
 * A                     (Adenine)
 * C                     (Cytosine)
 * G                     (Guanine)
 * T                     (Thymine)
 * R     = A or G        (puRines)
 * Y     = C or T        (pYrimidines)
 * W     = A or T        (Weak hydrogen bonding)
 * S     = G or C        (Strong hydrogen bonding)
 * M     = A or C        (aMino group at common position)
 * K     = G or T        (Keto group at common position)
 * H     = A, C or T     (not G)
 * B     = G, C or T     (not A)
 * V     = G, A, C       (not T)
 * D     = G, A or T     (not C)
 * N     = G, A, C or T  (aNy)
 *
 * Checked 3/16/01
 ******************************/
void printBlk(int w, struct mtfPt **blk, FILE *ofp) {
  int i, j, maxp, sum;
  double maxs, temp;
  char con, rcon, deg, rdeg;

  fprintf(ofp, "A      C      G      T         Con  rCon Deg  rDeg \n");
  for (i = 0; i < w; i++) {
    maxs = -100;    
    sum = 0;
    fprintf(ofp, "%-6d", (i+1));
    for (j = 0; j < ALPHASIZE; j++) {
      if (exp(blk[i][j].log) > 0.25)
	sum += (int)(rint(pow(2, j)));
      if (blk[i][j].log > maxs) {
	maxs = blk[i][j].log;
	maxp = j;
      }
      /*
      fprintf(ofp, "%6d", blk[i][j].ct);
      */
      temp = 100*exp(blk[i][j].log);
      if (temp < 10) {
	fprintf(ofp, " %.2f  ", temp);
      } else {
	fprintf(ofp, "%.2f  ", temp);
      }
    }
    
    /* maxp is for consensus */
    switch(maxp) {
    case A: con = 'A'; rcon = 'T'; break;
    case C: con = 'C'; rcon = 'G'; break;
    case G: con = 'G'; rcon = 'C'; break;
    case T: con = 'T'; rcon = 'A'; break;
    }
    
    /* sum is for degenerate */
    switch(sum) {
    case 0: deg = 'N'; rdeg = 'N'; break;
    case 1: deg = 'A'; rdeg = 'T'; break;
    case 2: deg = 'C'; rdeg = 'G'; break;
    case 3: deg = 'M'; rdeg = 'K'; break;  /* A, C */
    case 4: deg = 'G'; rdeg = 'C'; break;
    case 5: deg = 'R'; rdeg = 'Y'; break;  /* A, G */
    case 6: deg = 'S'; rdeg = 'S'; break;  /* G, C */
    case 7: deg = 'V'; rdeg = 'B'; break;  /* A, G, C */
    case 8: deg = 'T'; rdeg = 'A'; break;
    case 9: deg = 'W'; rdeg = 'W'; break;  /* A, T */
    case 10: deg = 'Y'; rdeg = 'R'; break; /* C, T */
    case 11: deg = 'H'; rdeg = 'D'; break; /* A, C, T */
    case 12: deg = 'K'; rdeg = 'M'; break; /* G, T */
    case 13: deg = 'D'; rdeg = 'H'; break; /* A, G, T */
    case 14: deg = 'B'; rdeg = 'V'; break; /* G, C, T */
    case 15: deg = 'N'; rdeg = 'N'; break;
    }
    fprintf(ofp, "      %c    %c    %c    %c\n", con, rcon, deg, rdeg);
  }
}

/******************************
 *
 * Func: getConsensus
 * Just get the consensus and reverse consensus
 *
 * Checked 3/16/01
 ******************************/
void getConsensus(int w, struct mtfPt **blk, char *con, char *rcon) {
  double maxs;
  int i, j, maxp;
  
  con[w] = rcon[w] = ENDLINE;
  for (i = 0; i < w; i++) {
    maxs = -1000;    
    for (j = 0; j < ALPHASIZE; j++) {
      if (blk[i][j].log > maxs) {
	maxs = blk[i][j].log;
	maxp = j;
      }
    }

    /* max is for consensus */
    switch(maxp) {
    case A: con[i] = 'A'; rcon[w-i-1] = 'T'; break;
    case C: con[i] = 'C'; rcon[w-i-1] = 'G'; break;
    case G: con[i] = 'G'; rcon[w-i-1] = 'C'; break;
    case T: con[i] = 'T'; rcon[w-i-1] = 'A'; break;
    }
  }
}
  
/******************************
 *
 * Func: printAlignment
 * Print sequence alignment: sequence name, length, number of motif copies
 * in this sequence, starting alignment position (r53 means position 53 in
 * reverse direction), aligned site sequence.
 *
 ******************************/
void printAlignment(struct inputParam *input, struct sequence *seqs,
		    struct alignment *align, int w1, int w2) {
  struct sequence *s;
  int i, j, *seqi;
  char line1[w1+1], line2[w2+1];
  struct alignment *a;

  line1[w1] = line2[w2] = ENDLINE;
  for (s = seqs, a = align; s; s = s->next, a = a->next) {
    if (a->copy) {
      
      for (i = 0; i < a->copy; i++) {
	
 	/* sequence name, len, copy */
	fprintf(input->ofp, "%s\tlen %d\tsite #%d", s->name, s->len, (i+1));
	/* motif starting position */
	if (a->dir[i] < BACKWARD) {
	  seqi = &(s->seqi[a->st1[i]]);
	  fprintf(input->ofp, "\tf %d", (a->st1[i]+1));
	} else {
	  seqi = &(s->rseqi[a->st1[i]]);
	  fprintf(input->ofp, "\tr %d", (s->len-a->st1[i]));
	}
	/* aligned site sequence */
	for (j = 0; j < w1; j++) {
	  switch(seqi[j]) {
	  case A: line1[j] = 'A'; break;
	  case C: line1[j] = 'C'; break;
	  case G: line1[j] = 'G'; break;
	  case T: line1[j] = 'T'; break;
	  }
	}
	if (w2) {
	  if (a->dir[i] > FORWARD) {
	    seqi = &(s->rseqi[a->st2[i]]);
	    if (a->dir[i] == PLDM)
	      fprintf(input->ofp, "\tr %d\n", (s->len-a->st2[i]));
	    else
	      fprintf(input->ofp, "\tr %d\n", (s->len-a->st2[i]));
	  } else {
	    seqi = &(s->seqi[a->st2[i]]);
	    fprintf(input->ofp, "\tf %d\n", (a->st2[i]+1));
	  } 
	  for (j = 0; j < w2; j++) {
	    switch(seqi[j]) {
	    case A: line2[j] = 'A'; break;
	    case C: line2[j] = 'C'; break;
	    case G: line2[j] = 'G'; break;
	    case T: line2[j] = 'T'; break;
	    }
	  }
	  fprintf(input->ofp, "%s %s\n", line1, line2);
	} else {
	  fprintf(input->ofp, "\n%s\n", line1);
	}
      }
    }
  }
  fflush(input->ofp);
}

/******************************
 *
 * Func: finishUp
 * Finish up, report the total time taken and close the log file
 *
 * Checked 10/19/00
 ******************************/
void finishUp(struct timeval *begintv, FILE *ofp) {
  long temp;
  struct timeval endtv;  

  if (input.local) {
    gettimeofday(&endtv, NULL);
    temp = endtv.tv_sec - begintv->tv_sec;
    fprintf(ofp, "Total time %d:%d:%d.\n", (int)(temp / SPH), 
	    (int)(temp % SPH / MPH) , (int)(temp % MPH));
  } else {
    fprintf(ofp, "\nThanks for using BioProspector. \n");
    fprintf(ofp, "For questions, please contact Xiaole Shirley Liu at xsliu@jimmy.harvard.edu.\n");
  }
  fclose(ofp);
  if (!input.local) 
    execl("./BPemail.pl", "BPemail", input.jobID, NULL);
}
