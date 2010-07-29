#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#endif

#include <math.h>
#include <time.h>
#include <ctype.h>
#include <sys/types.h>
#include <unistd.h>
#include "gadem.h"

#include "defines.h"
#include "random.h"
#include "evalue_meme.h"

// last modification 9/07/2009
// modifications:
//   1) fixed a minor bug in scan_sites.c
//   2) remove 6-mers, reduce search space
//   3) added C function getpid (process id)
//   4) added C function for cpu running time
//   5) set equal mutation rate for maxp and spaced dyads.
//      An optimal maxp may be important as it may affect how the EM converges
//   6) some cosmetic changes in output 
//   7) set default number of generations to 10
//   8) allowed a user-specified "seed" PWM
//   9) allowed a user-specified background model
//  10) included enrichment analysis
//  11) re-wrote pgf function (Staden's pgf method for llr null distribution)
//  12) fixed a bug in computing marginal probabilities for subsequences containing non-[a,c,g,t]
//  13) allow motif to overlap as an option

int main(int argc,char **argv) {

   register int jjj,ii,jj,i,j,k;
   // basic settings/info
   int numSeq,maxSeqLen,*seqLen;          // sequence info
   double aveSeqLen;                      // sequence info
   char **seq,**rseq,**geneID;            // sequence info
   char **oseq,**orseq;                   // copy of the original sequences
   char **pseq,**rpseq;                   // permuted seqs.
   double *bfreq;                         // base frequencies
   double *ChIPScore;                     // chip score

   // pwms
   double ***pwm;                         // initial population of PWMs from spaced dyads
   int *pwmLen;                           // initial pwm lengths 
   double **t1pwm,**t2pwm;                // two pwms before and after EM steps
   double **opwm2;                        // EM-derived PWM 
   double **logpwm;                       // log-transformed EM-derived EM
   int **ipwm;                            // integer pwm for computing llr score distribution
   double ***opwm;                        // observed PWMs from identified sites
   double ***epwm;                        // em-optimized PWMs
   double **logepwm;                      // log(em-optimized PWM)
   int *pwmnewLen;                        // final motif length after extending to both ends

   // llr score distr.
   Pgfs *llrDist;                         // llr distribution from pgf method
   int llrDim;                            // llr distribution dimension
   double *empDist;                       // empirical llr score distr.
   int empDim;                            // dimension of empirical llr score distribution 
   int pgf;                               // indicator for using pgf method or not

   // EM, motif, sites
   double pvalueCutoff;                   // user input, used to determine score cutoff based on ipwm
   int *scoreCutoff;                      // pwm score cutoff for the corresponding p-value cutoff
   double llrCutoff;
   double **score,**rscore;               // subsequence score, plus and minus strands
   double logev;                          // log of E-value of a motif;
   int useChIPscore;                      // indicator for using ChIP-seq score for seq. selection for EM
   int numEM;                             // number of EM steps
   double E_valueCutoff;                  // log E-value cutoff
   int nsitesEM;                          // number of binding sites in sequences subjected to EM
   int minsitesEM;                        // minimal number of sites in a motif in EM sequences
   Sites *siteEM;                         // binding sites in EM sequences
   int *nsites;                           // number of binding sites in full data
   int minsites;                          // minimal number of sites in a motif in full data
   Sites **site;                          // binding sites in all sequences
   int motifCn;                           // number of motifs sought and found
   int extTrim;
   int noMotifFound;                      // none of the dyads in the population resulted in a motif
   char **pwmConsensus;                   // consensus sequences of motifs
   double pwmDistCutoff;                  // test statistic for motif pwm similarity
   char *uniqMotif;                       // motifs in a population unique or not
   int numUniq;                           // number of unique motifs in a population
   int slideWinPWM;                       // sliding window for comparing pwm similarity
   int widthWt;                           // window width in which nucleotides are given large weights for PWM optimization 
   int fullScan;                          // scan scan on the original sequences or masked sequences

   // background
   BACKGROUND_Model *back;                // background model either user-specified or computed from forward data
   double **bscore,**rbscore;             // likelihood scores using background model
   double *pscore;                        // score for permuted seqs.
   int numTopWmerInB,numWmerInB;
   int numBackgSets;
   int MarkovOrder,userMarkovOrder;       // background Markov order,user-specified order
   int maxUserMarkovOrder;                // highest order available in user-specified background 
   int userBackgModel;                    // indicator

   // weights 
   double **posWeight;                    // spatial weights
   int weightType;                        // four weight types 0, 1, 2, 3, or 4

   // words for spaced dyad
   Words *word;                           // top-ranked k-mers as the words for spaced dyads
   int numTop3mer,numTop4mer,numTop5mer;  // No. of top-ranked k-mers as words for dyads
   int maxWordSize;                       // max of the above three
   int numWordGroup;                      // number of non-zero k-mer groups
   int minSpaceWidth,maxSpaceWidth;       // min and max width of spacer of the spaced dyads
   Chrs **dyad;                           // initial population of "chromosomes"
   char **sdyad;                          // char of spaced dyads

   // GA
   int populationSize,numGeneration;      // GA parameters
   double maxpMutationRate;
   Fitness *fitness;                      // "chromosome" fitness
   Wheel *wheel;                          // roulette-wheel selection

   // to speed up only select a subset of sequences for EM algorithm
   double fEM;                            // percentage of sequences used in EM algorithm
   int numSeqEM;                          // number of sequences subject to EM
   char *Iseq;                            // Indicator if a sequence is used in EM or not
   int *emSeqLen;                         // length of sequences used in EM
   double pwmDiff;                        // pwm convergence
   int maxp;                              // initial setting for all motifs maxp=numSeqEM 
   double *maxpFactor;

   int numCycle;                          // number of GADEM cycles
   int generationNoMotif;                 // maximal number of GA generations in a GADEM cycle resulted in no motifs

   // mis.
   seed_t  seed;                          // random seed
   int goodArgument,motifCn2,id,numCycleNoMotif,verbose,minminSites;
   int startPWMfound,stopCriterion;
   char *mFileName,*oFileName,*pwmFileName,*pwm0FileName,*bFileName,*sFileName;
   FILE *fp,*fq,*fpwm;
   time_t start,finish;
   int cn[4],bcn[4],*seqCn,*bseqCn,avebnsites,avebnsiteSeq,totalSitesInput;

   if (argc<3) {
      printf("\n");
      printf("                 *-------------------------------------------------------------------------\n");
      printf("                 |                                                                        |\n");
      printf("                 |      GADEM: a motif discovery tool for large scale sequence data       |\n");
      printf("                 |                                 v1.3                                   |\n");
      printf("                 |                                                                        |\n");
      printf("                 |         Multiple runs are recommended for 'unseeded' analysis.         |\n");
      printf("                 |      Each unseeded run automatically uses a different random seed.     |\n");
      printf("                 |        'Seeded' runs are deterministic; no repeat runs are needed.         |\n");
      printf("                 *-------------------------------------------------------------------------\n");
      printf("\n");
      printf("\n");
      printf("Usage: gadem -fseq seqFile optional arguments\n");
      printf("\n");
      printf("Optional arguments that need attention:\n");
      printf("\n");
      printf("  -posWt    0,1,2,or 3 Weight profile for positions on the sequence (see documentation).\n");
      printf("                       0 - no weight (uniform spatial prior, default), 1 - small or zero weights\n");
      printf("                       for the ends and large weights for the center (e.g. the center 50 bp)\n");
      printf("                       If you expect strong central enrichment (as in ChIP-seq) and your sequences\n");
      printf("                       are long (e.g. >200 bp), choose type 1.\n");
      printf("\n");
      printf("  -widthWt  integer    For -posWt 1 or 3, width of central sequence region with large EM weights\n");
      printf("                       for PWM optimization (default: 50). This argument is ignored when -posWt\n");
      printf("                       is 0 (uniform prior) or 2 (Gaussian prior).\n");
      printf("\n");
      printf("  -ev       decimal    ln(E-value) cutoff for selecting MOTIFS (default: 0.0).\n");
      printf("                       If a seeded analysis fails to identify the expected motif, run GADEM with\n");
      printf("                       -verbose 1 to show motif ln(E-value)s on screen, then rerun with a larger\n");
      printf("                       ln(E-value) cutoff. This can help in identifying short and/or low abundance\n");
      printf("                       motifs, for which the default E-value threshold may be too low.\n");
      printf("\n");
      printf("                       The subroutine for E-value calculation is adapted from the MEME package.\n");
      printf("\n");
      printf("  -pv       decimal    P-value cutoff for declaring BINDING SITES (default: 0.0002).\n");
      printf("                       Depending on data size and the motif, you might want to assess more than one\n");
      printf("                       value. For ChIP-seq data (e.g., 10 thousand +/-200-bp max-center peak 'cores'),\n");
      printf("                       p=0.0002 often seems appropriate. However, short motifs may require a less\n");
      printf("                       stringent setting.\n");
      printf("\n");
      printf("                       Given a subsequence s of length w, GADEM computes the log likelihood (llr)\n");
      printf("                       score, log{p(s|M)/p(s|B)}, where M is the EM-derived motif model, B is the\n");
      printf("                       b-th order Markov background model and w is the motif length. The subsequence\n");
      printf("                       is declared a binding site if its llr score is at or above the llr score\n");
      printf("                       corresponding to the p-value cutoff. This requires knowing the distribution\n");
      printf("                       of the llr score (under the null), and GADEM implements two approaches for\n");
      printf("                       approximating the null distribution: Staden probability generating function\n");
      printf("                       (pgf) method (Comput. Appl. Biosci., 5,89,1989) and an empirical approximation\n");
      printf("                       method through generating many background sequences. Both approaches are\n");
      printf("                       described briefly below.\n");
      printf("\n");
      printf("  -minN     integer    Minimal number of sites required for a motif to be reported (default: numSeq/20).\n");
      printf("\n");
      printf("  -fpwm0    string     File name for the seed PWM, when a 'seeded' approach is used. A PWM (format\n");
      printf("                       below) can be used as the starting PWM for the EM algorithm. This will help\n");
      printf("                       find an 'expected' motif and is much faster than 'unseeded' de novo discovery.\n");
      printf("                       Also, when a seed PWM is specified, the run results are deterministic, so only\n");
      printf("                       a single run is needed (repeat runs with the same settings will give identical\n");
      printf("                       results). In contrast, unseeded runs are stochastic, and we recommend comparing\n");
      printf("                       results from several repeat runs.\n");
      printf("\n");
      printf("                       Format: number of rows & columns followed by integer counts OR decimal\n");
      printf("                       frequencies.\n");
      printf("                       Example: PWM (CREB, JASPAR MA0018) in two acceptable representations:\n");
      printf("\n");
      printf("                       4     12\n");
      printf("                       0     3     0     2     5     0     0    16     0     0     1     5\n");
      printf("                       7     5     3     3     1     0     0     0    16     0     5     6\n");
      printf("                       5     4     6    11     7     0    15     0     0    16     0     3\n");
      printf("                       4     4     7     0     3    16     1     0     0     0    10     2\n");
      printf("\n");
      printf("                       4     12\n");
      printf("                       0.000 0.188 0.000 0.125 0.312 0.000 0.000 1.000 0.000 0.000 0.062 0.312\n");
      printf("                       0.438 0.312 0.188 0.188 0.062 0.000 0.000 0.000 1.000 0.000 0.312 0.375\n");
      printf("                       0.312 0.250 0.375 0.688 0.438 0.000 0.938 0.000 0.000 1.000 0.000 0.188\n");
      printf("                       0.250 0.250 0.438 0.000 0.188 1.000 0.062 0.000 0.000 0.000 0.625 0.125\n");
      printf("\n");
      printf("  -pgf       1 or 0    By default, GADEM uses the Staden probability generating function (pgf)\n");
      printf("                       method to approximate the exact llr score null distribution.\n");
      printf("\n");
      printf("                       The pgf method assumes that the background model is independent and\n");
      printf("                       identically-distributed (iid). When this method is used, GADEM takes the\n");
      printf("                       frequencies of [a,c,g,t] in the input data as estimates of the parameters\n");
      printf("                       of this iid model (0th-order), and a user-specified background model is not\n");
      printf("                       required. In other words, when -pgf is set to 1 (default) or is unspecified,\n");
      printf("                       both the -fbm and -bOrder flags are ignored by GADEM.\n");
      printf("\n");
      printf("                       Alternatively (when -pgf is set 0), the GADEM approximates the null using\n");
      printf("                       the llr scores of many background subsequences of length w, where w is the\n");
      printf("                       motif length. It generates the background subsequences using the [a,c,g,t]\n");
      printf("                       frequencies in the input data.\n");
      printf("\n");
      printf("  -bOrder   integer    The order of the background Markov model for computing llr scores:\n");
      printf("                       0 - 0th\n");
      printf("                       1 - 1st\n");
      printf("                       2 - 2nd\n");
      printf("                       ...\n");
      printf("                       8 - 8th\n");
      printf("\n");
      printf("  -fbm      string     Name of the file containing the user-specified background model. The GADEM\n");
      printf("                       package includes pre-computed genome-wide frequencies data for human, mouse\n");
      printf("                       and Drosophila, and source code for generating such files.\n");
      printf("\n");
      printf("                       The background Markov model can be estimated from the input data by GADEM or\n");
      printf("                       read from a file using the -fbm argument (see format below). We recommend\n");
      printf("                       -bOrder 0 when -fbm is not used. Otherwise, a higher order (e.g., -bOrder 3\n");
      printf("                       or 4) may be reasonable.\n");
      printf("\n");
      printf("                       Up to 8th order (octamer + 1nt = nonamer) is allowed.\n");
      printf("                       #monomer frequency\n");
      printf("                       a       0.20850000001660\n");
      printf("                       c       0.29149999998340\n");
      printf("                       g       0.29149999998340\n");
      printf("                       t       0.20850000001660\n");
      printf("                       #dimer frequency\n");
      printf("                       aa      0.04800960194357\n");
      printf("                       ac      0.05151030207800\n");
      printf("                       ag      0.08171634323790\n");
      printf("                       at      0.02720544114470\n");
      printf("                                    ....\n");
      printf("                       ta      0.03460692142891\n");
      printf("                       tc      0.05361072215865\n");
      printf("                       tg      0.07231446287687\n");
      printf("                       tt      0.04800960194357\n");
      printf("                       #trimer frequency\n");
      printf("                       aaa     0.01200480194395\n");
      printf("                       aac     0.01550620248175\n");
      printf("                       aag     0.01420568228200\n");
      printf("                       aat     0.00630252106809\n");
      printf("                                    ....\n");
      printf("                       tta     0.01190476192858\n");
      printf("                       ttc     0.01180472191322\n");
      printf("                       ttg     0.01230492199004\n");
      printf("                       ttt     0.01200480194395\n");
      printf("                                    ....\n");
      printf("\n");
      printf("                       If the file containing the background model is not specified and -pgf is set\n");
      printf("                       to 0, GADEM estimates the model from the input sequences. Note that when GADEM\n");
      printf("                       estimates a background model from input data that consists of short sequences\n");
      printf("                       (e.g. ChIP-seq), a higher order model is not recommended, as the sequences\n");
      printf("                       generated by the resulting background model may be too similar to the input\n");
      printf("                       sequences. For such cases we suggest setting -bOrder to 0.\n");
      printf("\n");
      printf("Other optional arguments:\n");
      printf("\n");
      printf("  -gen      integer    Number of genetic algorithm (GA) generations (default: 5).\n");
      printf("  -pop      integer    GA population size (default: 100).\n");
      printf("                       Both default settings should work well for most datasets (ChIP-chip and\n");
      printf("                       ChIP-seq).  The above two arguments are ignored in a seeded analysis,\n");
      printf("                       because spaced dyads and GA are no longer needed (-gen is set to 1 and\n");
      printf("                       -pop is set to 10 internally, corresponding to the 10 maxp choices).\n");
      printf("\n");
      printf("  -fullScan 0 or 1     GADEM keeps two copies of the input sequences internally: one (D) for\n");
      printf("                       discovering PWMs and one (S) for scanning for binding sites using the PWMs.\n");
      printf("                       Once a motif is identified, its instances in set D are always masked by Ns.\n");
      printf("                       However, masking motif instances in set S is optional, and scanning unmasked\n");
      printf("                       sequences allows sites of discovered motifs to overlap.\n");
      printf("\n");
      printf("                       0 (default) - scan masked sequences in S (disallow motif site overlap).\n");
      printf("                       1 - scan unmasked sequences in S (allow motif sites to overlap) (was default\n");
      printf("                       in v1.2).\n");
      printf("\n");
      printf("  -em       integer    Number of EM steps (default: 40). One might want to set it to a larger value\n");
      printf("                       (e.g. 80) in a seeded run, because such runs are fast.\n");
      printf("\n");
      printf("  -fEM      decimal    Fraction of sequences used in EM to obtain PWMs in an unseeded analysis\n");
      printf("                       (default: 0.5). For unseeded motif discovery in a large dataset (e.g. >10\n");
      printf("                       million nt), one might want to set -fEM to a smaller value (e.g., 0.3 or 0.4)\n");
      printf("                       to reduce run time.\n");
      printf("\n");
      printf("                       When only partial input data are used in EM and ?verbose is set to 1, the\n");
      printf("                       number of binding sites printed on screen is the number of sites found only\n");
      printf("                       in the sequences that are used in EM optimization.\n");
      printf("\n");
      printf("                       This argument is ignored in a seeded analysis, which uses all sequences EM to\n");
      printf("                       obtain the PWMs.\n");
      printf("\n");
      printf("  -extTrim  1 or 0     Base extension and trimming (1 -yes, 0 -no) (default: 1).\n");
      printf("\n");
      printf("  -maxw3    integer    Number of top-ranked trimers for spaced dyads (default: 20).\n");
      printf("  -maxw4    integer    Number of top-ranked tetramers for spaced dyads (default: 40).\n");
      printf("  -maxw5    integer    Number of top-ranked pentamers for spaced dyads (default: 60).\n");
      printf("\n");
      printf("  -mingap   integer    Minimal number of unspecified nucleotides in spaced dyads (default: 0).\n");
      printf("  -maxgap   integer    Maximal number of unspecified nucleotides in spaced dyads (default: 10).\n");
      printf("                       -mingap and -maxgap control the lengths of spaced dyads, and, with -extrim,\n");
      printf("                       control motif lengths. Longer motifs can be discovered by setting -maxgap to\n");
      printf("                       larger values (e.g. 50).\n");
      printf("\n");
      printf("  -useScore 0 or 1     Use top-scoring sequences for deriving PWMs. Sequence (quality) scores are\n");
      printf("                       stored in sequence header (see documentation).\n");
      printf("                       0 - no (default, randomly select sequences), 1 - yes.\n");
      printf("\n");
      printf("  -fpwm     string     Name of output PWM file in STAMP format (http://www.benoslab.pitt.edu/stamp).\n");
      printf("                       (default: observedPWMs.txt). This file can be loaded into STAMP to compare\n");
      printf("                       each PWM with PWMs in databases for similarity.\n");
      printf("\n");
      printf("  -fout     string     Name of main GADEM output file (see documentation for description) (default:\n");
      printf("                       gadem.txt).\n");
      printf("\n");
      printf("  -nbs      integer    Number of sets of background sequences (default: 10). The background sequences\n");
      printf("                       are simulated using the [a,c,g,t] frequencies in the input sequences, with\n");
      printf("                       length matched between the two sets. The background sequences are used as the\n");
      printf("                       random sequences for assessing motif enrichment in the input data. Another set\n");
      printf("                       (same default: 10) of background sequences is independently generated to\n");
      printf("                       approximate the empirical llr score distribution when -pgf is set to 0.\n");
      printf("\n");
      printf("  -verbose  1 or 0     Print immediate results on screen [1-yes (default), 0-no]. These results\n");
      printf("                       include the motif consensus sequence, number of sites (in sequences subjected\n");
      printf("                       to EM optimization, see -fEM, above), and ln(E-value).\n");
      printf("\n");
      printf("-------------------------------------------------------------------------------------------\n");
      printf("Examples:\n");
      printf("\n");
      printf("1. Unseeded analysis for ChIP data with expected central enrichment (llr distr.: pgf method-default)\n");
      printf("   gadem -fseq input.seq -minN 1000 -posWt 1 -verbose 1\n");
      printf("\n");
      printf("2. Unseeded analysis for ChIP data with expected central enrichment (llr distr.: empirical approx.)\n");
      printf("   gadem -fseq input.seq -minN 1000 -posWt 1 -pgf 0 -fbm freq.txt -bOrder 4 -verbose 1\n");
      printf("\n");
      printf("3. Seeded analysis for ChIP data with expected central enrichment (llr distr.: pgf method-default)\n");
      printf("   gadem -fseq input.seq -minN 1000 -posWt 1 -fpwm0 startPWM.mx -verbose 1\n");
      printf("\n");
      printf("4. Seeded analysis for ChIP data with expected central enrichment (llr distr.: empirical approx.)\n");
      printf("   gadem -fseq input.seq -minN 1000 -posWt 1 -fpwm0 startPWM.mx -pgf 0 -fbm freq.txt -bOrder 4 -verbose 1\n");
      printf("\n");
      printf("\n");

      exit(0);
   }

   mFileName=alloc_char(500);     mFileName[0]='\0';
   oFileName=alloc_char(500);     oFileName[0]='\0';
   pwmFileName=alloc_char(500);   pwmFileName[0]='\0';
   bFileName=alloc_char(500);     bFileName[0]='\0';
   sFileName=alloc_char(500);     sFileName[0]='\0';
   seq=NULL; aveSeqLen=0; maxSeqLen=0; numSeq=0; 
   minsites=-1; 
   
   // default settings
   numWordGroup=3;
   numTop3mer=20; numTop4mer=40; numTop5mer=60;
   numGeneration=5; populationSize=100;
   pvalueCutoff=0.0002;
   E_valueCutoff=0.0; 
   extTrim=1;
   minSpaceWidth=0; maxSpaceWidth=10;
   useChIPscore=0; 
   numEM=40; fEM=0.5; widthWt=80; fullScan=0;
   userBackgModel=0;
   slideWinPWM=6; numUniq=populationSize;
   strcpy(oFileName,"gadem.txt"); 
   strcpy(pwmFileName,"observedPWMs.txt"); 
   stopCriterion=NUM_NO_MOTIF;  
   MarkovOrder=0; userMarkovOrder=0;
   numBackgSets=10; 
   weightType=0;
   pgf=1;
   verbose=0;
   startPWMfound=0;

   back=alloc_background();

   for (ii=1; ii<argc-1; ii++) {
      if (argv[ii][0]=='-' && isalpha(argv[ii][1])) {
         goodArgument=0;
         if (
             (strncmp(argv[ii],"-ev",3)==0       && strlen(argv[ii])==3) ||
             (strncmp(argv[ii],"-pv",3)==0       && strlen(argv[ii])==3) ||
             (strncmp(argv[ii],"-em",3)==0       && strlen(argv[ii])==3) ||
             (strncmp(argv[ii],"-pop",4)==0      && strlen(argv[ii])==4) ||
             (strncmp(argv[ii],"-gen",4)==0      && strlen(argv[ii])==4) ||
             (strncmp(argv[ii],"-pgf",4)==0      && strlen(argv[ii])==4) ||
             (strncmp(argv[ii],"-fEM",4)==0      && strlen(argv[ii])==4) ||
             (strncmp(argv[ii],"-nbs",4)==0      && strlen(argv[ii])==4) ||
             (strncmp(argv[ii],"-fbm",4)==0      && strlen(argv[ii])==4) ||
             (strncmp(argv[ii],"-fseq",5)==0     && strlen(argv[ii])==5) ||
             (strncmp(argv[ii],"-fout",5)==0     && strlen(argv[ii])==5) ||
             (strncmp(argv[ii],"-fpwm",5)==0     && strlen(argv[ii])==5) ||
             (strncmp(argv[ii],"-minN",5)==0     && strlen(argv[ii])==5) ||
             (strncmp(argv[ii],"-maxw3",6)==0    && strlen(argv[ii])==6) ||
             (strncmp(argv[ii],"-maxw4",6)==0    && strlen(argv[ii])==6) ||
             (strncmp(argv[ii],"-maxw5",6)==0    && strlen(argv[ii])==6) ||
             (strncmp(argv[ii],"-fpwm0",6)==0    && strlen(argv[ii])==6) ||
             (strncmp(argv[ii],"-posWt",6)==0    && strlen(argv[ii])==6) ||
             (strncmp(argv[ii],"-mingap",7)==0   && strlen(argv[ii])==7) ||
             (strncmp(argv[ii],"-maxgap",7)==0   && strlen(argv[ii])==7) ||
             (strncmp(argv[ii],"-bOrder",7)==0   && strlen(argv[ii])==7) ||
             (strncmp(argv[ii],"-extTrim",8)==0  && strlen(argv[ii])==8) ||
             (strncmp(argv[ii],"-verbose",8)==0  && strlen(argv[ii])==8) ||
             (strncmp(argv[ii],"-widthWt",8)==0  && strlen(argv[ii])==8) ||
             (strncmp(argv[ii],"-fullScan",9)==0 && strlen(argv[ii])==9) ||
             (strncmp(argv[ii],"-useScore",9)==0 && strlen(argv[ii])==9)
            ) { goodArgument=1; }
         if (!goodArgument) { printf("argument: %s unknown\n",argv[ii]); exit(0);  }
      }
   }

   ChIPScore=alloc_double(MAX_NUM_SEQ);
   seqLen=alloc_int(MAX_NUM_SEQ); 
   geneID=alloc_char_char(MAX_NUM_SEQ,500);

   for (ii=1; ii<argc; ii++) {
      if (strncmp(argv[ii],"-fseq",5)==0 && argv[ii+1]!=NULL) {
         //printf("\nreading input sequences file: %s\n",argv[ii+1]);
         strcpy(mFileName,argv[ii+1]);
         seq=read_seq(&numSeq,seqLen,geneID,MAX_NUM_SEQ,MAX_SEQ_LENGTH,ChIPScore,argv[ii+1]);
         aveSeqLen=0; for (i=0; i<numSeq; i++) aveSeqLen +=seqLen[i]; aveSeqLen /=(double)numSeq;
         maxSeqLen=0; 
         for (i=0; i<numSeq; i++) {
            if (seqLen[i]>maxSeqLen) maxSeqLen=seqLen[i]; 
         }
 
         rseq=alloc_char_char(numSeq,maxSeqLen+1);
         oseq=alloc_char_char(numSeq,maxSeqLen+1);
         orseq=alloc_char_char(numSeq,maxSeqLen+1);

         reverse_seq(seq,rseq,numSeq,seqLen);
         // make a copy of the original sequences both strands
         for (i=0; i<numSeq; i++) {
            for (j=0; j<seqLen[i]; j++) { oseq[i][j]=seq[i][j]; orseq[i][j]=rseq[i][j]; } 
            oseq[i][seqLen[i]]='\0'; orseq[i][seqLen[i]]='\0'; 
         }
      }
      else if (strncmp(argv[ii],"-bOrder",7)==0 && argv[ii+1]!=NULL) userMarkovOrder=atoi(argv[ii+1]);
      else if (strncmp(argv[ii],"-fbm",4)==0 && argv[ii+1]!=NULL) {
         strcpy(bFileName,argv[ii+1]);
         // printf("\nreading user-specified background models: %s\n",bFileName);
         maxUserMarkovOrder=read_userBackgModel(argv[ii+1],back);
         // printf("highest available order: %d\n\n",maxUserMarkovOrder);
         if (userMarkovOrder>maxUserMarkovOrder) {
            printf("highest Markov order in %s is %d\n",argv[ii+1],maxUserMarkovOrder);
            printf("re-set user-specified Markov order to %d\n",maxUserMarkovOrder);
            userMarkovOrder=maxUserMarkovOrder;
         }
         userBackgModel=1;
      }
      else if (strncmp(argv[ii],"-fpwm0",6)==0 && argv[ii+1]!=NULL) {
         strcpy(sFileName,argv[ii+1]);
         printf("\n|------------------------------------------------------------------|\n");
         printf("|                                                                  |\n");
         printf("|               *** Running a seeded analysis ***                  |\n");
         printf("|                                                                  |\n");
         printf("|------------------------------------------------------------------|\n\n");
         // printf("reading user-specified seed pwm:\t%s\n",sFileName);
         populationSize=FIXED_POPULATION;
         dyad  =alloc_chrs(populationSize,4);
         sdyad =alloc_char_char(populationSize,MAX_PWM_LENGTH+1);
         pwmLen=alloc_int(populationSize);
         pwm=alloc_double_double_double(populationSize,MAX_PWM_LENGTH,4);
         pwmLen[0]=read_pwm0(argv[ii+1],pwm[0]);
         for (i=1; i<populationSize; i++) {
            for (j=0; j<pwmLen[0]; j++) {
               for (k=0; k<4; k++) pwm[i][j][k]=pwm[0][j][k]; 
            }
            pwmLen[i]=pwmLen[0]; 
         }
         pwm0FileName=alloc_char(200);
         strcpy(pwm0FileName,argv[ii+1]);
         startPWMfound=1;
      }
      else if (strncmp(argv[ii],"-pv",3)==0         && strlen(argv[ii])==3  && argv[ii+1]!=NULL) pvalueCutoff=atof(argv[ii+1]); 
      else if (strncmp(argv[ii],"-em",3)==0         && strlen(argv[ii])==3  && argv[ii+1]!=NULL) numEM=atoi(argv[ii+1]);
      else if (strncmp(argv[ii],"-ev",3)==0         && strlen(argv[ii])==3  && argv[ii+1]!=NULL) E_valueCutoff=atof(argv[ii+1]);
      else if (strncmp(argv[ii],"-pop",4)==0        && strlen(argv[ii])==4  && argv[ii+1]!=NULL) populationSize=atoi(argv[ii+1]); 
      else if (strncmp(argv[ii],"-gen",4)==0        && strlen(argv[ii])==4  && argv[ii+1]!=NULL) numGeneration=atoi(argv[ii+1]); 
      else if (strncmp(argv[ii],"-pgf",4)==0        && strlen(argv[ii])==4  && argv[ii+1]!=NULL) pgf=atoi(argv[ii+1]);
      else if (strncmp(argv[ii],"-fEM",4)==0        && strlen(argv[ii])==4  && argv[ii+1]!=NULL) fEM=atof(argv[ii+1]);
      else if (strncmp(argv[ii],"-nbs",4)==0        && strlen(argv[ii])==4  && argv[ii+1]!=NULL) numBackgSets=atoi(argv[ii+1]);
      else if (strncmp(argv[ii],"-fout",5)==0       && strlen(argv[ii])==5  && argv[ii+1]!=NULL) strcpy(oFileName,argv[ii+1]); 
      else if (strncmp(argv[ii],"-fpwm",5)==0       && strlen(argv[ii])==5  && argv[ii+1]!=NULL) strcpy(pwmFileName,argv[ii+1]); 
      else if (strncmp(argv[ii],"-minN",5)==0       && strlen(argv[ii])==5  && argv[ii+1]!=NULL) minsites=atoi(argv[ii+1]);
      else if (strncmp(argv[ii],"-maxw3",6)==0      && strlen(argv[ii])==6  && argv[ii+1]!=NULL) numTop3mer=atoi(argv[ii+1]);
      else if (strncmp(argv[ii],"-maxw4",6)==0      && strlen(argv[ii])==6  && argv[ii+1]!=NULL) numTop4mer=atoi(argv[ii+1]);
      else if (strncmp(argv[ii],"-maxw5",6)==0      && strlen(argv[ii])==6  && argv[ii+1]!=NULL) numTop5mer=atoi(argv[ii+1]);
      else if (strncmp(argv[ii],"-posWt",6)==0      && strlen(argv[ii])==6  && argv[ii+1]!=NULL) weightType=atoi(argv[ii+1]);
      else if (strncmp(argv[ii],"-mingap",7)==0     && strlen(argv[ii])==7  && argv[ii+1]!=NULL) minSpaceWidth=atoi(argv[ii+1]);
      else if (strncmp(argv[ii],"-maxgap",7)==0     && strlen(argv[ii])==7  && argv[ii+1]!=NULL) maxSpaceWidth=atoi(argv[ii+1]);
      else if (strncmp(argv[ii],"-extTrim",8)==0    && strlen(argv[ii])==8  && argv[ii+1]!=NULL) extTrim=atoi(argv[ii+1]);
      else if (strncmp(argv[ii],"-verbose",8)==0    && strlen(argv[ii])==8  && argv[ii+1]!=NULL) verbose=atoi(argv[ii+1]);
      else if (strncmp(argv[ii],"-widthWt",8)==0    && strlen(argv[ii])==8  && argv[ii+1]!=NULL) widthWt=atoi(argv[ii+1]);
      else if (strncmp(argv[ii],"-fullScan",9)==0   && strlen(argv[ii])==9  && argv[ii+1]!=NULL) fullScan=atoi(argv[ii+1]);
      else if (strncmp(argv[ii],"-useScore",9)==0   && strlen(argv[ii])==9  && argv[ii+1]!=NULL) useChIPscore=atoi(argv[ii+1]);
      else { }
   }

   // check for input parameters
   if (numGeneration<1)  { printf("\nError: numbe of generaton < 1.\n"); exit(0); }
   if (populationSize<1) { printf("\nError: population size < 1.\n");    exit(0); }
   if (minSpaceWidth<0)  { 
      printf("\nError: minimal number of unspecified bases in spaced dyads <0.\n"); 
      printf("   check -mingap setting\n\n"); exit(0);
   }
   if (maxSpaceWidth<0)  { 
      printf("\nError: maximal number of unspecified bases in spaced dyads <0.\n"); 
      printf("   check -maxgap setting\n\n"); exit(0);
   }
   if (minSpaceWidth>maxSpaceWidth) {
      printf("\nError: mingap setting must <= to maxgap setting.\n\n"); exit(0); 
   }
   if (maxSpaceWidth+12>MAX_PWM_LENGTH) {
      printf("\nError: maxgap setting plus word lengths exceed <MAX_PWM_LENGTH>.\n");
      printf("   For very long motifs, please set <MAX_PWM_LENGTH> in 'defines.h' accordingly.\n\n"); exit(0); 
   }
   if (numEM<0) {
      printf("\nError: number of EM steps is zero.\n"); exit(0); 
   }
   if (numEM==0) {
      printf("\nNote: number of EM steps = 0, no EM optimization is carried out.\n");  
   }

   if (fullScan!=0 && fullScan!=1) fullScan=0;

   maxWordSize=0;
   if (numTop3mer>maxWordSize) maxWordSize=numTop3mer;
   if (numTop4mer>maxWordSize) maxWordSize=numTop4mer;
   if (numTop5mer>maxWordSize) maxWordSize=numTop5mer;

   // any one, two or three: tetramer, pentamer, hexamer
   if (numTop3mer==0 && numTop4mer==0 && numTop5mer==0) {
      printf("\nError: maxw3, maxw4, and maxw5 all zero - no words for spaced dyads.\n"); exit(0);
   }

   if (startPWMfound && fEM!=0.5 && fEM!=1.0) {
      printf("\n***Note: -fEM argument is ignored in a seeded analysis***\n\n");
   }

   if (startPWMfound) {
      if (populationSize!=10 && populationSize!=100) printf("\n***Note: -pop argument is ignored in a seeded analysis, -pop is set to 10.***\n\n");
      if (numGeneration!=1 && numGeneration!=5)      printf("\n***Note: -gen argument is ignored in a seeded analysis, -gen is set to 1.***\n\n");
      fEM=1.0;
      populationSize=FIXED_POPULATION; numGeneration=1; 
   }

   // number of sequences for EM
   if (fEM>1.0 || fEM<=0.0) { 
      printf("\nError: the fraction of sequences subject to EM is %3.2f.\n",fEM); exit(0); 
   } 
   numSeqEM=(int)(fEM*numSeq);

   if (pgf!=0 && pgf!=1)   { printf("\nError: -pgf can only take 0 or 1.\n"); exit(0); }

   // memory callocations
   if (!startPWMfound) {
      printf("\n|------------------------------------------------------------------|\n");
      printf("|                                                                  |\n");
      printf("|              *** Running an unseeded analysis ***                |\n");
      printf("|                                                                  |\n");
      printf("|------------------------------------------------------------------|\n\n");
      dyad  =alloc_chrs(populationSize,4);
      pwm   =alloc_double_double_double(populationSize,MAX_PWM_LENGTH,4);
      pwmLen=alloc_int(populationSize);
      sdyad =alloc_char_char(populationSize,MAX_PWM_LENGTH+1);
      word  =alloc_word(numWordGroup,maxWordSize);
   }
   Iseq  =alloc_char(numSeq+1); 
   opwm  =alloc_double_double_double(populationSize,MAX_PWM_LENGTH,4);
   opwm2 =alloc_double_double(MAX_PWM_LENGTH,4);
   t1pwm =alloc_double_double(MAX_PWM_LENGTH,4);
   t2pwm =alloc_double_double(MAX_PWM_LENGTH,4);
   ipwm  =alloc_int_int(MAX_PWM_LENGTH,4);
   logpwm=alloc_double_double(MAX_PWM_LENGTH,4);
   score =alloc_double_double(numSeq,maxSeqLen);
   rscore=alloc_double_double(numSeq,maxSeqLen);
   epwm=alloc_double_double_double(populationSize,MAX_PWM_LENGTH,4);
   logepwm=alloc_double_double(MAX_PWM_LENGTH,4);
   wheel =alloc_wheel(populationSize);
   siteEM=alloc_site(MAX_SITES);
   fitness=alloc_fitness(populationSize);
   emSeqLen=alloc_int(numSeqEM);
   maxpFactor=alloc_double(populationSize);
   pwmConsensus=alloc_char_char(populationSize,MAX_PWM_LENGTH+1);
   uniqMotif=alloc_char(populationSize+1);
   scoreCutoff=alloc_int(populationSize);
   llrDist=alloc_distr(MAX_DIMENSION);
   posWeight=alloc_double_double(numSeq,maxSeqLen);
   bscore=alloc_double_double(numSeq,maxSeqLen);
   rbscore=alloc_double_double(numSeq,maxSeqLen);
   pscore=alloc_double(numSeq*maxSeqLen);
   empDist=alloc_double((int)(numSeq*2*maxSeqLen*numBackgSets/4));
   pseq=alloc_char_char(MAX_NUM_SEQ,maxSeqLen+1);
   rpseq=alloc_char_char(MAX_NUM_SEQ,maxSeqLen+1);

   bfreq=base_frequency(numSeq,seq,seqLen);

   // if minN not specified, set the defaults accordingly
   if (minsites==-1) minsites =max(2,(int)(numSeq/20)); 
   minsitesEM=(int)(fEM*minsites);

   maxpMutationRate=MAXP_MUTATION_RATE;

   seed=time(0);
   sgenrand(seed);

   // determine the distribution and critical cut point
   pwmDistCutoff=vector_similarity();

   /*---------- select a subset of sequences for EM only --------------*/
   if (useChIPscore==1) {
      select_high_scoring_seq_for_EM (ChIPScore,numSeq,numSeqEM,Iseq,fEM);
   }
   else {
      sample_without_replacement(Iseq,numSeqEM,numSeq);
   }
   /*-------------------- end of selection --------------------------*/

   if (!userBackgModel) {
      if (!pgf && userMarkovOrder!=0 && aveSeqLen<=500) {
         printf("\n***Note: it is not recommended to use a non-zero-th Markov model estimated from\n");
         printf("   the input sequences as the background model, especially for short sequences as\n");
         printf("   the sequences generated by the resulting background model may be too similar to\n");
         printf("   the input sequences\n\n"); 
      }
      printf("\nestimating background Markov models using input sequences...\n");
      generate_background(numSeq,seq,rseq,seqLen,back,userMarkovOrder); 
      printf("done\n\n");
   }
   if (widthWt<20) {
      printf("\n***Note: the window width of sequence centered on the nucleotides having large weights\n");
      printf("   in EM for PWM optimization is small: %d\n",widthWt);
      printf("   Motif longer than %d will not be discovered\n\n",widthWt); 
   }

   time(&start);
   fp=fopen("info.txt","w");
   fprintf(fp,"processor ID: %d\n",getpid());
   fprintf(fp,"==========================================================================\n");
   if (startPWMfound) {
      fprintf(fp,"\n|------------------------------------------------------------------|\n");
      fprintf(fp,"|                                                                  |\n");
      fprintf(fp,"|              *** Running a seeded analysis ***                   |\n");
      fprintf(fp,"|                                                                  |\n");
      fprintf(fp,"|------------------------------------------------------------------|\n\n");
   }
   else {
      fprintf(fp,"\n|------------------------------------------------------------------|\n");
      fprintf(fp,"|                                                                  |\n");
      fprintf(fp,"|              *** Running an unseeded analysis ***                |\n");
      fprintf(fp,"|                                                                  |\n");
      fprintf(fp,"|------------------------------------------------------------------|\n\n");
   }
   fprintf(fp,"command line: ");
   for (i=0; i<argc; i++) fprintf(fp,"%s ",argv[i]); fprintf(fp,"\n\n");

   fprintf(fp,"maximal buffer length:\t\t\t\t\t%d\n",MAX_BUFFER_LENGTH);
   fprintf(fp,"maximal number of sequences set:\t\t\t%d\n",MAX_NUM_SEQ);
   fprintf(fp,"maximal number of bases per seq read:\t\t\t%d\n",MAX_SEQ_LENGTH);
   fprintf(fp,"maximal number of sites in a motif:\t\t\t%d\n",MAX_SITES);
   fprintf(fp,"input (ChIP) sequence file:\t\t\t\t%s\n",mFileName);
   fprintf(fp,"number of sequences in input file:\t\t\t%d\n",numSeq);
   fprintf(fp,"average sequence length:\t\t\t\t%d\n",(int)aveSeqLen);
   fprintf(fp,"total number of nucleotides:\t\t\t\t%d\n",(int)(aveSeqLen*numSeq));
   fprintf(fp,"max number of generations:\t\t\t\t%d\n",numGeneration);
   fprintf(fp,"population size:\t\t\t\t\t%d\n",populationSize);
   if (fullScan) fprintf(fp,"alway scan original unmaksed sequences (allow motif sites to overlap).\n");
   if (!startPWMfound) {
      fprintf(fp,"number of top-ranked 3-, 4-, 5-mer:\t\t\t%d %d %d\n",numTop3mer,numTop4mer,numTop5mer);
      fprintf(fp,"minimal & maximal spacer(d) in dyads\t\t\t%d %d\n",minSpaceWidth,maxSpaceWidth);
      if (useChIPscore) fprintf(fp,"use high scoring sequences to derive PWMs\n");
      else              fprintf(fp,"use randomly select sequences to derive PWMs\n");
   }
   else fprintf(fp,"use a user-specified pwm as the seed\t\t\t%s\n",sFileName);
   fprintf(fp,"fraction (number) input sequences subject to EM\t\t%3.2f (%d)\n",fEM,numSeqEM);
   fprintf(fp,"scale factor for converting (double)pwm to (int)pwm\t%4.0f\n",DOUBLE_TO_INT_SCALE);
   fprintf(fp,"number of EM steps:\t\t\t\t\t%d\n",numEM);
   fprintf(fp,"EM convergence criterion:\t\t\t\t\%e\n",PWM_CONVERGENCE);
   if (startPWMfound) {
      fprintf(fp,"run EM on the starting pwm %s %d times, each with a different maxp:\n",pwm0FileName,FIXED_POPULATION);
      for (i=0; i<FIXED_POPULATION; i++) fprintf(fp,"%3.2f*numSeq ",FIXED_MAXPF*(i+1)); fprintf(fp,"\n");
      fprintf(fp,"no spaced dyads are generated and used. pop=%d gen=1 (no GA).\n\n",FIXED_POPULATION);
   }
   else {
      fprintf(fp,"\nGADEM determines maxps by choosing ");
      fprintf(fp,"one of the five values (factor*numSeq):\n ");
      fprintf(fp,"   ");
      for (i=0; i<MAXP_SCALE; i++) fprintf(fp,"%3.2f*numSeq, ",MAXP_BASE+(double)i*MAXP_FACTOR); fprintf(fp,"\n");
      fprintf(fp,"These factors are optimized along with dyads by the GA and reported.\n\n");
   }
   fprintf(fp,"motif prior probability type (see documentation):\t%1d\n",weightType);

   if (weightType==1 || weightType==3) 
      fprintf(fp,"window width of sequence centered on the nucleotides having large weights for PWM optimization: %d\n",widthWt);
   fprintf(fp,"pwm score p-value cutoff for declaring binding site:\t%e\n",pvalueCutoff);

   if (pgf) {
      if (userMarkovOrder!=0) {
         fprintf(fp,"\nNote: the specified background Markov order: %d is ignored when -pgf is set to 1\n",userMarkovOrder);
      }
      if (bFileName[0]!='\0') {
         fprintf(fp,"\nNote: the user-specified background models: %s are not used when -pgf is set to 1\n\n",bFileName);
      }
      fprintf(fp,"\nApproximate the null llr score log{p(s|M)/p(s|B)} distribution using Staden's probability\n");
      fprintf(fp,"probability generating funtion (pgf) method, where s is a subsequence of length w,\n");
      fprintf(fp,"M is the EM-derived motif model from GADEM and B is the 0th-order Markov model from\n");
      fprintf(fp,"the input data - the frequencies of [a,c,g,t]\n\n");
   }
   else {
      fprintf(fp,"\nApproximate the null llr log{p(s|M)/p(s|B)} score distribution using\n");
      fprintf(fp,"the llr scores of random/background sequences, where M is the EM-derived\n");
      fprintf(fp,"motif model and B is the %d-th order Markov backgroun model.\n\n",userMarkovOrder);
      fprintf(fp,"The background sequences are simulated using the [a,c,g,t] frequencies\n");
      fprintf(fp,"in the input data. The number sets of background sequences generated: \t%d\n\n",numBackgSets);
   }
   fprintf(fp,"pseudo count:\t\t\t\t\t\t%5.4f\n",PSEUDO_COUNT);
   fprintf(fp,"minimal infomation for trimming/extending: \t\t%3.2f %3.2f %3.2f\n",MIN_BITS1,MIN_BITS2,MIN_BITS3);
   fprintf(fp,"minimal no. sites for each motif:\t\t\t%d\n",minsites);
   if (extTrim) fprintf(fp,"base extension and trimming? \t\t\t\tyes\n");
   else         fprintf(fp,"base extension and trimming? \t\t\t\tno\n");
   fprintf(fp,"sliding window for comparing pwm similarity:\t%d\n",slideWinPWM);
   fprintf(fp,"PWM similarity cutoff:\t\t\t\t\t%5.3f\n",SIMILARITY_ALPHA);
   fprintf(fp,"log(E-value) cutoff:\t\t\t\t\t%5.2f\n",E_valueCutoff);
   fprintf(fp,"number of adjacent bases included in binding site output:\t%d\n",FLANKING_BASES);
   if (!startPWMfound) fprintf(fp,"random seed:\t\t\t\t\t\t%ld\n\n",seed);
   fprintf(fp,"\njob started: %s\n", asctime(localtime(&start)));
   fprintf(fp,"=========================================================================\n\n");
   fflush(fp);

   printf("==============================================================================================\n");
   printf("input sequence file:  %s\n",mFileName);
   printf("number of sequences and average length:\t\t\t\t%d %5.1f\n",numSeq,aveSeqLen);
   if (startPWMfound) {
      printf("\nseed PWM: %s\n",sFileName);
      for (j=0; j<4; j++) {
         for (i=0; i<pwmLen[0]; i++) {
            if (i<pwmLen[0]-1) printf("%5.2f ", pwm[0][i][j]);
            else               printf("%5.2f\n",pwm[0][i][j]);
         }
      }
   }
   if (pgf) {
      printf("\nUse pgf method to approximate llr null distribution\n");
      printf("background Markov order:\t\t\t\t\t0th\n");
      printf("parameters estimated from sequences in:  %s\n\n",mFileName);
   }
   else {
      printf("\nUse an empirical approach to approximate llr null using background sequences\n");
      printf("Background sequences are simulated using [a,c,g,t] frequencies in input data\n");
      printf("Background Markov order (in llr calculation):\t\t\t\t");
      switch (userMarkovOrder) {
         case 0: printf("0th\n"); break; 
         case 1: printf("1st\n"); break; 
         case 2: printf("2nd\n"); break; 
         case 3: printf("3rd\n"); break; 
         case 4: printf("4th\n"); break; 
         case 5: printf("5th\n"); break; 
         case 6: printf("6th\n"); break; 
         case 7: printf("7th\n"); break; 
         case 8: printf("8th\n"); break; 
         default: break;
      }
      if (bFileName[0]!='\0') printf("parameters from:\t\t\t\t\t\t%s\n\n",bFileName);
      else                    printf("parameters estimated from sequences in:\t\t%s\n\n",mFileName);
   }
   if (weightType!=0) printf("non-uniform weight applies to each sequence - type:\t\t%d\n",weightType);
   printf("number of GA generations & population size:\t\t\t%d %d\n\n",numGeneration,populationSize);
   printf("PWM score p-value cutoff for binding site declaration:\t\t%e\n",pvalueCutoff);
   printf("ln(E-value) cutoff for motif declaration:\t\t\t%f\n\n",E_valueCutoff);
   printf("number (percentage) of sequences selected for EM:\t\t%d(%4.1f\%)\n",numSeqEM,100.0*(double)numSeqEM/(double)numSeq);
   printf("number of EM steps:\t\t\t\t\t\t%d\n",numEM);
   printf("minimal no. sites considered for a motif:\t\t\t%d\n\n",minsites);
   printf("[a,c,g,t] frequencies in input data:\t\t\t\t%f %f %f %f\n",bfreq[0],bfreq[1],bfreq[2],bfreq[3]);
   printf("==============================================================================================\n");

   if (pgf) {
      if (userMarkovOrder!=0) {
         printf("\n***The user-specified background Markov order (%d) is ignored when -pgf is set to 1***\n",userMarkovOrder);
      }
      if (bFileName[0]!='\0') {
         printf("\n***The user-specified background models: %s are not used when -pgf is set to 1***\n\n",bFileName);
      }
   }
   if (startPWMfound && fEM!=1.0) {
      printf("\n***Note: -fEM argument is ignored in a seeded analysis***\n\n");
   }

   printf("\nsSarting GADEM... this may take a few hours to complete\n");
   printf("type: ctrl/z, then, bg, to run it in background\n\n");

   // determine seq length by counting only [a,c,g,t], seqLen is used in E-value calculation 
   effect_seq_length(seq,numSeq,seqLen,Iseq,emSeqLen);
   // determine the distribution and critical cut point
   pwmDistCutoff=vector_similarity();

   if      (weightType==0) assign_weight_uniform(seqLen,numSeq,posWeight);
   else if (weightType==1) assign_weight_rectangle(seqLen,numSeq,posWeight,widthWt);
   else if (weightType==2) assign_weight_triangular(seqLen,numSeq,posWeight);
   else if (weightType==3) assign_weight_triangular_uniform(seqLen,numSeq,posWeight,widthWt);
   else if (weightType==4) assign_weight_normal(seqLen,numSeq,posWeight);
   else {
      printf("Motif prior probability type not found - please choose: 0, 1, 2, 3, or 4\n");
      printf("Consider: -posWt 1 for strong central enrichment as in ChIP-seq\n");
      printf("          -posWt 0 for others\n\n");
      exit(0);
   }

   if (startPWMfound) minminSites=minsites;
   else               minminSites=(int)(0.40*minsitesEM);

   fq=fopen(oFileName,"w");
   fpwm=fopen(pwmFileName,"w");
   motifCn=0; noMotifFound=0; numCycle=0; numCycleNoMotif=0; 

   if (startPWMfound) { 
      for (i=0; i<populationSize; i++) {
         standardize_pwm(pwm[i],pwmLen[i]);
         consensus_pwm(pwm[i],pwmLen[i],pwmConsensus[i]);
         strcpy(sdyad[i],pwmConsensus[i]); 
      }
   }
 
   do {
      if (!startPWMfound) {
         // identify top-ranked k-mers (k=3,4,5) for spaced dyads
         printf("\nGADEM cycle %2d: enumerate and count k-mers...   ",numCycle+1);
         numWordGroup=word_for_dyad(word,seq,rseq,numSeq,seqLen,bfreq,&numTop3mer,&numTop4mer,&numTop5mer);
         printf("done.\n");

         // generating a "population" of spaced dyads
         printf("\ninitializing GA...   ");
         initialisation(dyad,populationSize,numWordGroup,word,minSpaceWidth,maxSpaceWidth,maxpFactor);
         printf("done.\n\n");
      }
      else {
         for (i=0; i<populationSize; i++) maxpFactor[i]=FIXED_MAXPF*(i+1); 
      }

      generationNoMotif=0;
      for (jjj=0; jjj<numGeneration; jjj++) {

         // convert spaced dyads to letter probability matrix
         if (!startPWMfound) dyad_to_pwm(word,populationSize,dyad,pwm,pwmLen);

         for (ii=0; ii<populationSize; ii++) {
            // to see from which spaced dyad a motif is derived
            if (!startPWMfound) pwm_profile(pwm[ii],pwmLen[ii],sdyad[ii]);

            // make a copy and then subject the copy to EM
            copy_pwm(pwm[ii],t1pwm,pwmLen[ii]); 
            // standarize pwm
            standardize_pwm(t1pwm,pwmLen[ii]);
            // for (j=0; j<4; j++) { for (m=0; m<pwmLen[0]; m++) printf("%4.3f ",t1pwm[m][j]); printf("\n"); } printf("\n"); 
           
            // EM on randomly selected sequences
            maxp=(int)(maxpFactor[ii]*numSeqEM); 

            for (jj=0; jj<numEM; jj++) {
               log_pwm(t1pwm,logpwm,pwmLen[ii]);
               // compute ll score of each w-mer | motif model
               ll_score_motif_model(numSeq,seq,rseq,seqLen,logpwm,pwmLen[ii],score,rscore,Iseq,bfreq);
               // compute p(zij|y=1) probability of binding sites started at position j on seq i
               normalize(score,rscore,seqLen,pwmLen[ii],numSeq,Iseq,maxp,posWeight,weightType);
               // E-step
               construct_pwm(t2pwm,score,rscore,seq,rseq,seqLen,numSeq,pwmLen[ii],Iseq);
               // M-step
               standardize_pwm(t2pwm,pwmLen[ii]);
               // printf("EM: %2d\n",jj+1);
               // for (j=0; j<4; j++) { for (m=0; m<pwmLen[0]; m++) printf("%5.4f ",t2pwm[m][j]); printf("\n"); } printf("\n"); 
               pwmDiff=check_convergence(t1pwm,t2pwm,pwmLen[ii]);
               // copy t2pwm to t1pwm
               copy_pwm(t2pwm,t1pwm,pwmLen[ii]); 
               if (pwmDiff<=PWM_CONVERGENCE)  break;
            }
            copy_pwm(t1pwm,epwm[ii],pwmLen[ii]); // from to
            // for (j=0; j<4; j++) { for (i=0; i<pwmLen[0]; i++) printf("%4.3f ",t2pwm[i][j]); printf("\n"); } printf("\n\n");  

            // log(PWM), then (double)log(PWM) to (int)PWM for determine score distribution
            log_ratio_to_int(epwm[ii],ipwm,pwmLen[ii],bfreq);
            // for (j=0; j<4; j++) { for (i=0; i<pwmLen[ii]; i++) printf("%5d ",ipwm[i][j]); printf("\n"); } printf("\n\n"); 

            // compute score distribution of the (int)PWM using Staden's method 
            llrDim=pwm_score_dist(ipwm,pwmLen[ii],llrDist,bfreq);
            // print_ptable(llrDist,llrDim);

            scoreCutoff[ii]=determine_cutoff(llrDist,llrDim,pvalueCutoff);

            // test each w-mer to see if a motif site - test statistic: ll, distribution: Staden method, cutoff: user-specified
            nsitesEM=scan_em_seq_ptable(llrDist,llrDim,siteEM,numSeq,seq,rseq,seqLen,ipwm,pwmLen[ii],scoreCutoff[ii],bfreq,Iseq);
            // printf("maxp %d nsitesEM: %d pcutoff: %e scorecutoff: %d\n",maxp,nsitesEM,pvalueCutoff,scoreCutoff[ii]);
   
            // loose threshould at this step, as em only on a subset of sequences
            if (nsitesEM>=max(2,minminSites)) {
               // construct pwm from the identified sites
               align_sites_count(siteEM,seq,rseq,nsitesEM,pwmLen[ii],opwm[ii]);
               standardize_pwm(opwm[ii],pwmLen[ii]);
               consensus_pwm(opwm[ii],pwmLen[ii],pwmConsensus[ii]);
               // compute E-value of the relative entroy score of each motif, use it as fitness
               fitness[ii].value=E_value(opwm[ii],nsitesEM,bfreq,pwmLen[ii],numSeqEM,emSeqLen);
            }
            else {
               // if too few sites in a motif
               align_sites_count(siteEM,seq,rseq,nsitesEM,pwmLen[ii],opwm[ii]);
               standardize_pwm(opwm[ii],pwmLen[ii]);
               consensus_pwm(opwm[ii],pwmLen[ii],pwmConsensus[ii]);
               fitness[ii].value=DUMMY_FITNESS;
               // for (i=0; i<pwmLen[ii]; i++) pwmConsensus[ii][i]='n'; pwmConsensus[ii][pwmLen[ii]]='\0'; 
            }
            fitness[ii].index=ii;
            if (verbose) { 
               printf("cyc.[%3d] gen.[%3d] pop.[%3d] spacedDyad: %s ",numCycle+1,jjj+1,ii+1,sdyad[ii]);
               for (j=strlen(sdyad[ii]); j<maxSpaceWidth+10; j++) printf(" ");
               printf(" motifConsensus: %s",pwmConsensus[ii]);
               for (j=strlen(sdyad[ii]); j<maxSpaceWidth+10; j++) printf(" ");
               //printf(" fitness: %7.2f\n",fitness[ii].value);
               printf(" maxpf: %3.2f fitness: %7.2f nsitesEM: %d\n",maxpFactor[ii],fitness[ii].value,nsitesEM);
            }
         }
         if (populationSize>1) sort_fitness(fitness,populationSize);

         /*-----------------------------------------------------------------------------
         printf("generation %3d top %d:\n", jjj+1, populationSize);
         for (i=0; i<populationSize; i++) 
            printf("   pwmConsens: %s\tlogev: %6.2f\n",pwmConsensus[fitness[i].index],fitness[i].value);
         printf("\n");  
         -----------------------------------------------------------------------------*/
        
         numUniq=check_pwm_uniqueness_dist(opwm,pwmLen,populationSize,fitness,pwmDistCutoff,E_valueCutoff,uniqMotif,slideWinPWM);

         printf("\nGADEM cycle[%3d] generation[%3d] number of unique motif: %d\n",numCycle+1,jjj+1,numUniq);
         fprintf(fp,"GADEM cycle[%3d] generation[%3d] number of unique motif(s): %d\n",numCycle+1,jjj+1,numUniq);
         for (i=0; i<populationSize; i++) {
            if (uniqMotif[i]=='1') {
               printf("   spacedDyad: %s ",sdyad[fitness[i].index]);
               for (j=strlen(sdyad[fitness[i].index]); j<maxSpaceWidth+10; j++) printf(" ");
               printf("motifConsensus: %s ",pwmConsensus[fitness[i].index]);
               for (j=strlen(sdyad[fitness[i].index]); j<maxSpaceWidth+10; j++) printf(" ");
               printf(" %3.2f fitness: %7.2f\n",maxpFactor[fitness[i].index],fitness[i].value);

               fprintf(fp,"   spacedDyad: %s ",sdyad[fitness[i].index]);
               for (j=strlen(sdyad[fitness[i].index]); j<maxSpaceWidth+10; j++) fprintf(fp," ");
               fprintf(fp,"motifConsensus: %s ",pwmConsensus[fitness[i].index]);
               for (j=strlen(sdyad[fitness[i].index]); j<maxSpaceWidth+10; j++) fprintf(fp," ");
               // fprintf(fp," fitness: %7.2f\n",fitness[i].value);
               fprintf(fp," %3.2f fitness: %7.2f\n",maxpFactor[fitness[i].index],fitness[i].value);
               //for (j=0; j<4; j++) {
               //   for (m=0; m<pwmLen[fitness[i].index]; m++) printf("%4.2f ",opwm[fitness[i].index][m][j]); printf("\n");
               //} printf("\n\n");
            }
         }
         printf("\n"); fprintf(fp,"\n"); fflush(fp);

         if (jjj<numGeneration-1) {

            // fitness based selection with replacement 
            roulett_wheel_fitness(fitness,populationSize,wheel);

            // mutation and crossover operations
            if (populationSize>1) {
               if (genrand()>=0.5) {
                  mutation (dyad,numWordGroup,word,minSpaceWidth,maxSpaceWidth,wheel,populationSize,fitness,uniqMotif,
                     maxpFactor,maxpMutationRate); 
               }
               else { 
                  crossover(dyad,numWordGroup,word,minSpaceWidth,maxSpaceWidth,wheel,populationSize,fitness,uniqMotif,
                     maxpFactor,maxpMutationRate); 
               }
            }
            else { 
               mutation (dyad,numWordGroup,word,minSpaceWidth,maxSpaceWidth,wheel,populationSize,fitness,uniqMotif,
                 maxpFactor,maxpMutationRate); 
            }
         }
      }
      numCycle++;

      site  =alloc_site_site(numUniq+1,MAX_SITES);
      nsites=alloc_int(numUniq+1);
      pwmnewLen=alloc_int(numUniq+1); // after base extension and trimming
      seqCn=alloc_int(MAX_NUM_SEQ);
      bseqCn=alloc_int(MAX_NUM_SEQ);

      // final step user-specified background model is used
      motifCn2=0; // motifCn per GADEM cycle
      for (ii=0; ii<populationSize; ii++) {

         id=fitness[ii].index;
         // for (j=0; j<4; j++) { for (i=0; i<pwmLen[id]; i++) printf("%4.3f ",epwm[id][i][j]); printf("\n"); } printf("\n");  exit(0);

         if (uniqMotif[ii]=='0') continue;

         MarkovOrder=min(pwmLen[id]-1,userMarkovOrder);

         // approximate the exact llr distribution using Staden's method
         if (pgf) {
            printf("\nApproximate the exact pwm llr score distribution using the pgf method.\n");
            log_ratio_to_int(epwm[id],ipwm,pwmLen[id],bfreq);

            // compute score distribution of the (int)PWM using Staden's method 
            llrDim=pwm_score_dist(ipwm,pwmLen[id],llrDist,bfreq);
            scoreCutoff[id]=determine_cutoff(llrDist,llrDim,pvalueCutoff);
            if (fullScan) {
               nsites[motifCn2]=scan_llr_pgf(llrDist,llrDim,site[motifCn2],numSeq,oseq,orseq,seqLen,ipwm,
                  pwmLen[id],scoreCutoff[id],bfreq);
            }
            else {
               nsites[motifCn2]=scan_llr_pgf(llrDist,llrDim,site[motifCn2],numSeq,seq,rseq,seqLen,ipwm,
                  pwmLen[id],scoreCutoff[id],bfreq);
            }
         }
         // determine the null llr distribution using background sequences
         else {
            log_pwm(epwm[id],logepwm,pwmLen[id]);

            /* -----------------compute the null distribtion------------------------------------*/
            // this generates N*(L-w+1)*numBackgSets w-mers compared to N*(L-w+1) in input data

            printf("\nUse an empirical approach to approximate the llr score distribution\n");

            numTopWmerInB=0; empDim=0;
            for (i=0; i<numBackgSets; i++) {

               simulate_background_seq(bfreq,numSeq,seqLen,pseq);
               ll_score_backg_model(numSeq, pseq, bscore,seqLen,pwmLen[id],back,MarkovOrder);
               numWmerInB=llr_score(pscore,numSeq,pseq,seqLen,logepwm,pwmLen[id],bfreq,bscore) ;
               sort_double(pscore,numWmerInB);        // truncated null dist.

               for (j=0; j<numWmerInB/4; j++) {       // take only top25%
                  empDist[numTopWmerInB]=pscore[j];   // pool all top25% from each permutation
                  numTopWmerInB++;
               }
               empDim +=numWmerInB;
            }
            sort_double(empDist,numTopWmerInB);      // truncated null dist.
            /* -----------------end the null distribtion------------------------------------*/
        
            llrCutoff=empDist[(int)(pvalueCutoff*empDim)];
            // print_null(empDist,numTopWmerInB,empDim);
            if (fullScan) {
               ll_score_backg_model(numSeq, oseq, bscore,seqLen,pwmLen[id],back,MarkovOrder);
               ll_score_backg_model(numSeq,orseq,rbscore,seqLen,pwmLen[id],back,MarkovOrder);
               nsites[motifCn2]=scan_llr_empirical(site[motifCn2],numSeq,oseq,orseq,seqLen,logepwm,pwmLen[id],bfreq,bscore,rbscore,
                  llrCutoff,empDist,numTopWmerInB,empDim);
            }
            else {
               ll_score_backg_model(numSeq, seq, bscore,seqLen,pwmLen[id],back,MarkovOrder);
               ll_score_backg_model(numSeq,rseq,rbscore,seqLen,pwmLen[id],back,MarkovOrder);
               nsites[motifCn2]=scan_llr_empirical(site[motifCn2],numSeq,seq,rseq,seqLen,logepwm,pwmLen[id],bfreq,bscore,rbscore,
                  llrCutoff,empDist,numTopWmerInB,empDim);
            }
         }

         if (nsites[motifCn2]>=max(2,minsites)) {
            for (j=0; j<numSeq; j++) seqCn[j]=0;
            for (j=0; j<nsites[motifCn2]; j++) seqCn[site[motifCn2][j].seq]++;
      
            for (j=0; j<4; j++) cn[j]=0;
            for (j=0; j<numSeq; j++) {
               if (seqCn[j]==0) cn[0]++;
               if (seqCn[j]==1) cn[1]++;
               if (seqCn[j]==2) cn[2]++;
               if (seqCn[j]>2)  cn[3]++;
            }

            totalSitesInput=nsites[motifCn2];
            if (extTrim) {
               if (fullScan) {
                  extend_alignment(site[motifCn2],numSeq,oseq,orseq,seqLen,nsites[motifCn2],pwmLen[id],&(pwmnewLen[motifCn2]));
               }
               else {
                  extend_alignment(site[motifCn2],numSeq,seq,rseq,seqLen,nsites[motifCn2],pwmLen[id],&(pwmnewLen[motifCn2]));
               }
            }
            else { pwmnewLen[motifCn2]=pwmLen[id]; } 

            if (fullScan) {
               align_sites_count(site[motifCn2],oseq,orseq,nsites[motifCn2],pwmnewLen[motifCn2],opwm2);
            }
            else {
               align_sites_count(site[motifCn2],seq,rseq,nsites[motifCn2],pwmnewLen[motifCn2],opwm2);
            }
            standardize_pwm(opwm2,pwmnewLen[motifCn2]);

            logev=E_value(opwm2,nsites[motifCn2],bfreq,pwmnewLen[motifCn2],numSeq,seqLen);
            if (logev<=E_valueCutoff) {
               consensus_pwm(opwm2,pwmnewLen[motifCn2],pwmConsensus[id]);
               if (fullScan) {
                  print_result_2(site[motifCn2],nsites[motifCn2],numSeq,oseq,orseq,seqLen,geneID,logev,opwm2,pwmnewLen[motifCn2],
                      motifCn+1,sdyad[id],pwmConsensus[id],numCycle,pvalueCutoff,maxpFactor[id],fq,fpwm);
                  print_motif(site[motifCn2],nsites[motifCn2],oseq,orseq,seqLen,pwmnewLen[motifCn2],motifCn+1,opwm2);
               }
               else {
                  print_result_2(site[motifCn2],nsites[motifCn2],numSeq,seq,rseq,seqLen,geneID,logev,opwm2,pwmnewLen[motifCn2],
                      motifCn+1,sdyad[id],pwmConsensus[id],numCycle,pvalueCutoff,maxpFactor[id],fq,fpwm);
                  print_motif(site[motifCn2],nsites[motifCn2],seq,rseq,seqLen,pwmnewLen[motifCn2],motifCn+1,opwm2);
               }
               mask_sites(nsites[motifCn2],seq,rseq,seqLen,site[motifCn2],pwmnewLen[motifCn2]);

               /* ----------------------compute the average number of sites in background sequences ----------------------*/
               avebnsites=0; avebnsiteSeq=0;
               for (i=0; i<numBackgSets; i++) {

                  simulate_background_seq(bfreq,numSeq,seqLen,pseq);
                  reverse_seq(pseq,rpseq,numSeq,seqLen);

                  if (pgf) {
                     nsites[motifCn2]=scan_llr_pgf(llrDist,llrDim,site[motifCn2],numSeq,pseq,rpseq,seqLen,ipwm,pwmLen[id],
                       scoreCutoff[id],bfreq);
                  }
                  else {
                     ll_score_backg_model(numSeq, pseq, bscore,seqLen,pwmLen[id],back,MarkovOrder);
                     ll_score_backg_model(numSeq,rpseq,rbscore,seqLen,pwmLen[id],back,MarkovOrder);
                     nsites[motifCn2]=scan_llr_empirical(site[motifCn2],numSeq,pseq,rpseq,seqLen,logepwm,pwmLen[id],bfreq,bscore,rbscore,
                        llrCutoff,empDist,numTopWmerInB,empDim);
                  }

                  for (j=0; j<numSeq; j++) bseqCn[j]=0;
                  for (j=0; j<nsites[motifCn2]; j++) bseqCn[site[motifCn2][j].seq]++;
      
                  for (j=0; j<4; j++) bcn[j]=0;
                  for (j=0; j<numSeq; j++) {
                     if (bseqCn[j]==0) bcn[0]++;
                     if (bseqCn[j]==1) bcn[1]++;
                     if (bseqCn[j]==2) bcn[2]++;
                     if (bseqCn[j]>2)  bcn[3]++;
                  }
                  fprintf(fq,"background set[%2d] Seqs with 0,1,2,>2 sites: %d %d %d %d\n",i+1,bcn[0],bcn[1],bcn[2],bcn[3]);
                  avebnsites+=nsites[motifCn2]; avebnsiteSeq+=(numSeq-bcn[0]);
               } 
               avebnsites/=numBackgSets; avebnsiteSeq/=numBackgSets;
               fprintf(fq,"average number of sites in background sequences: %d, fold enrichment: %5.3f.\n",
                  avebnsites,(double)totalSitesInput/(double)avebnsites);
               fprintf(fq,"average number of background sequences that contain at least one site: %d, fold enrichment: %5.3f.\n",
                  avebnsiteSeq,(double)(cn[1]+cn[2]+cn[3])/(double)(bcn[1]+bcn[2]+bcn[3]));
               fprintf(fq,"-------------------------------------------------------\n");
               fflush(fq); 
               /* -----------------end compute the average number of sites in background sequences ----------------------*/
               motifCn++; motifCn2++; numCycleNoMotif=0;
            } 
         }
      }
      for (i=0; i<motifCn2; i++) {
         // mask_sites(nsites[i],seq,rseq,seqLen,site[i],pwmnewLen[i]); 
         // fq=fopen("tmp.seq","w");
         // for (i=0; i<numSeq; i++) {
         //   fprintf(fq,"%s\n",geneID[i]); fprintf(fq,"%s\n",seq[i]);
         //}
         //fclose(fq);
         //exit(0);
      }
      if (site[0])   { free(site[0]);   site[0]=NULL;   }
      if (site)      { free(site);      site=NULL;      }
      if (nsites)    { free(nsites);    nsites=NULL;    }
      if (pwmnewLen) { free(pwmnewLen); pwmnewLen=NULL; }
      
      if (motifCn2==0) numCycleNoMotif++;
    
      if (numCycleNoMotif==stopCriterion) noMotifFound=1;

   } while (!noMotifFound);
   fclose(fq); fclose(fpwm);

   time(&finish);
   fprintf(fp,"\nfinished: %s\n", asctime(localtime(&finish)));
   fprintf(fp,"approximated processor time in seconds: %f\n",difftime(finish,start));
   fclose(fp);
   system("mv info.txt info.done.txt");

   if (!startPWMfound) {  
      if (dyad[0])      { free(dyad[0]);         dyad[0]=NULL;    }
      if (dyad)         { free(dyad);            dyad=NULL;       }
   }
   if (seqLen)          { free(seqLen);          seqLen=NULL;     }
   if (pwm[0][0])       { free(pwm[0][0]);       pwm[0][0]=NULL;  }
   if (pwm[0])          { free(pwm[0]);          pwm[0]=NULL;     }
   if (pwm)             { free(pwm);             pwm=NULL;        }
   if (opwm2[0])        { free(opwm2[0]);        opwm2[0]=NULL;   }
   if (opwm2)           { free(opwm2);           opwm2=NULL;      }
   if (opwm[0][0])      { free(opwm[0][0]);      opwm[0][0]=NULL; }
   if (opwm[0])         { free(opwm[0]);         opwm[0]=NULL;    }
   if (opwm)            { free(opwm);            opwm=NULL;       }
   if (t1pwm[0])        { free(t1pwm[0]);        t1pwm[0]=NULL;   }
   if (t1pwm)           { free(t1pwm);           t1pwm=NULL;      }
   if (t2pwm[0])        { free(t2pwm[0]);        t2pwm[0]=NULL;   }
   if (t2pwm)           { free(t2pwm);           t2pwm=NULL;      }
   if (logpwm[0])       { free(logpwm[0]);       logpwm[0]=NULL;  }
   if (logpwm)          { free(logpwm);          logpwm=NULL;     }
   if (ipwm[0])         { free(ipwm[0]);         ipwm[0]=NULL;    }
   if (ipwm)            { free(ipwm);            ipwm=NULL;       }
   if (pwmLen)          { free(pwmLen);          pwmLen=NULL;     }
   if (seq[0])          { free(seq[0]);          seq[0]=NULL;     }
   if (seq)             { free(seq);             seq=NULL;        }
   if (rseq[0])         { free(rseq[0]);         rseq[0]=NULL;    }
   if (rseq)            { free(rseq);            rseq=NULL;       }
   if (geneID[0])       { free(geneID[0]);       geneID[0]=NULL;  }
   if (oseq[0])         { free(oseq[0]);         oseq[0]=NULL;    }
   if (oseq)            { free(oseq);            oseq=NULL;       }
   if (orseq[0])        { free(orseq[0]);        orseq[0]=NULL;   }
   if (orseq)           { free(orseq);           orseq=NULL;      }
   if (geneID)          { free(geneID);          geneID=NULL;     }
   if (score[0])        { free(score[0]);        score[0]=NULL;   }
   if (score)           { free(score);           score=NULL;      }
   if (rscore[0])       { free(rscore[0]);       rscore[0]=NULL;  }
   if (rscore)          { free(rscore);          rscore=NULL;     }
   if (bfreq)           { free(bfreq);           bfreq=NULL;      }
   if (wheel)           { free(wheel);           wheel=NULL;      }
   if (fitness)         { free(fitness);         fitness=NULL;    }
   if (mFileName)       { free(mFileName);       mFileName=NULL;  }
   if (oFileName)       { free(oFileName);       oFileName=NULL;  }
   if (pwmFileName)     { free(pwmFileName);     pwmFileName=NULL;}
   if (sdyad[0])        { free(sdyad[0]);        sdyad[0]=NULL;   }
   if (sdyad)           { free(sdyad);           sdyad=NULL;      }
   if (siteEM)          { free(siteEM);          siteEM=NULL;     }
   if (pwmConsensus[0]) { free(pwmConsensus[0]); pwmConsensus[0]=NULL; }
   if (pwmConsensus)    { free(pwmConsensus);    pwmConsensus=NULL;    }
   if (!startPWMfound && word) destroy_word(word,numWordGroup);

   return (1);
}

void print_ptable(Pgfs *llrDist,int llrDim) {

   FILE *fp;
   int i;

   fp=fopen("ptable.txt","w");
   for (i=0; i<llrDim; i++) fprintf(fp,"%d\t%e\n",llrDist[i].score,llrDist[i].prob);
   fclose(fp);
}

void print_empirical(double *empDist,int totalKmer) {

   FILE *fp;
   int i;

   fp=fopen("empirical.txt","w");
   for (i=0; i<totalKmer; i++) fprintf(fp,"%9.6f\n",empDist[i]);
   fclose(fp);
}

void select_high_scoring_seq_for_EM (double *ChIPScore,int numSeq,int numSeqEM,char *Iseq,double fEM) {

   register int i;
   int numSeqWithQualityScore,numSeqEMtmp1,numSeqEMtmp2;
   double *tmpScore;
   double ChIPscoreCutoff;

   tmpScore=alloc_double(numSeq);

   numSeqWithQualityScore=0;
   for (i=0; i<numSeq; i++) {
      if (ChIPScore[i]>0) numSeqWithQualityScore++;
   }

   tmpScore=alloc_double(numSeq);
   for (i=0; i<numSeq; i++) tmpScore[i]=ChIPScore[i];
   sort_double(tmpScore,numSeq);

   ChIPscoreCutoff=tmpScore[(int)(fEM*numSeq)];

   if (numSeqWithQualityScore<=(int)(fEM*numSeq)) {
      for (i=0; i<numSeq; i++) Iseq[i]='0';
      numSeqEMtmp1=0;
      for (i=0; i<numSeq; i++) {
         if (ChIPScore[i]>0) {
            Iseq[i]='1'; numSeqEMtmp1++;
         }
      }
      numSeqEMtmp2=0;
      for (i=0; i<numSeq; i++) {
         if (ChIPScore[i]<=0) {
            Iseq[i]='1'; numSeqEMtmp2++;
            if (numSeqEMtmp1+numSeqEMtmp2==numSeqEM) break;
         }
      }
   }
   else {
      for (i=0; i<numSeq; i++) Iseq[i]='0';
      numSeqEMtmp1=0; numSeqEMtmp2=0;
      for (i=0; i<numSeq; i++) {
         if (ChIPScore[i]>=ChIPscoreCutoff) {
            Iseq[i]='1'; numSeqEMtmp1++;
            if (numSeqEMtmp1==numSeqEM) break;
         }
      }
   }
   if (tmpScore)  { free(tmpScore);  tmpScore=NULL;  }
   if (ChIPScore) { free(ChIPScore); ChIPScore=NULL; }
   
}

void print_null(double *empDist,int numTopWmerInB,int empDim) {

   FILE *fp;
   int i;

   fp=fopen("tmp.txt","w");
   fprintf(fp,"%d %d\n",numTopWmerInB,empDim);
   for (i=0; i<numTopWmerInB; i++) fprintf(fp,"%f\n",empDist[i]);
   fclose(fp);
}
