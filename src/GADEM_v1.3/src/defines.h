
#define MAX_BUFFER_LENGTH	15000    // maximal number of characters per line
#define MAX_NUM_SEQ             44000    // maximal number of sequences
#define MAX_SEQ_LENGTH          20000    // maximal length of any sequence in input data
#define MAX_SEQ_HEADER          40       // max length of sequence header read
#define MAX_SITES               150000   // maximal number of sites in the entire data per motif
#define MAX_PWM_LENGTH          120      // maximal pwm/motif length
#define MAX_DIMENSION           75000    // maximal llr score distribution dimension
#define Z_CUTOFF                6.0      // z-score cutoff for top-ranked k-mers
#define MAX_ORDER               8        // max background Markov order

// #define PSEUDO_COUNT	        0.0001   // pseudo count 
#define PSEUDO_COUNT	        0.0005   // pseudo count 
                                         // this setting seems to have a profound effect on some data e.g., ChIP-Chip ERE
                                         // setting it to 0.0001 failed to identify the ERE in 4 of 5 runs, whereas setting
                                         // it to 0.0005, the ERE was found in all runs (quickly)
#define PSEUDO_FREQ	        0.000001 // pseudo frequency used in background.c

#define DOUBLE_TO_INT_SCALE	200.0    // (int)PWM=(int)(PWM*DOUBLE_TO_INT_SCALE), see transform.c

#define MIN_BITS1		0.40     // bits cutoff (out of 2 bits) for extending/trimming motif, choice one
#define MIN_BITS2		0.50     // choice two 
#define MIN_BITS3		0.60     // choice three

#define BASE_EXT		10       // # of bases extended at each side of a motif
#define FLANKING_BASES		10       // output motifs with 'BASE_EXT' bases on each side

#define DUMMY_FITNESS		999999.0 // dummy fitness score for motifs that do not meet minimal requirements

#define NUM_NO_MOTIF		1        // GADEM stops when the number of GADEM cycles produces no motifs
#define SIMILARITY_ALPHA	0.30     // PWM similarity goodness-of-fit alpha level
                                         // a very value (e.g., 0.0001) will let GADEM to report all motifs found regardless of
                                         // similarity. This may be useful for a seeded analysis

#define PWM_CONVERGENCE         0.0001   // EM convergence criterion
                                         // the setting needs to be smaller than PSEUDO_COUNT. When this parameter and pseudo count
                                         // were set to be the same, the sampe top-ranked spaced dyad in one generation converted into
                                         // a different motif in the following generation. This happened only rarely.

//#define MAJOR                   0.5    // MEME's setting
//#define MINOR                   0.1666666666 // MEME's setting
#define CELL_MAJOR		1.0
#define CELL_MINOR		0.0

#define MAXP_MUTATION_RATE	0.25    // a1-nx-a2...maxp equiprobably being mutated
#define MAXP_BASE		0.1
#define MAXP_FACTOR		0.1
#define MAXP_SCALE		10      // maxpfactor =MAXP_BASE+MAXP_FACTOR*((int)(MAXP_SCALE*genrand()));
                                        //            =0.1 + 0.1*0 = 0.1
                                        //            =0.1 + 0.1*1 = 0.2
                                        //            =0.1 + 0.1*2 = 0.3
                                        //            =0.1 + 0.1*3 = 0.4
                                        //            =0.1 + 0.1*4 = 0.5 etc
                                        // maxp=maxpfactor*numSeq
                                        // estimated number of binding sites for computing P(Yij=1|PWM,Seq),
                                        // the probability of site start at position j on seq i (see GADEM suppl. text)

#define FIXED_MAXPF		0.1     // increment
#define FIXED_POPULATION	10      // ten evenly spaced maxp 0.1*N, 0.2*N, 0.3*N, 0.4*N, 0.5*N, 0.6*N, 0.7*N, 0.8*N, 0.9*N, 1.0*N

#define min(a,b)                (((a)<(b))?(a):(b))
#define max(a,b)                (((a)>(b))?(a):(b))
#define mod(a,b)		((a)-(b)*((a)/(b)))

