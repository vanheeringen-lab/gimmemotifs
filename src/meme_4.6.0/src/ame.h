/********************************************************************
 * FILE: ame.h
 * AUTHOR: Robert McLeay
 * CREATE DATE: 19/08/2008
 * PROJECT: MEME suite
 * COPYRIGHT: 2008, Robert McLeay
 *
 * rsClover is a yet unpublished algorithm that seeks to assist in
 * determining whether a given transcription factor regulates a set
 * of genes that have been found by an experimental method that
 * provides ranks, such as microarray or ChIp-chip.
 *
 *
 * rsClover is a code name. The name will change prior to publication.
 ********************************************************************/
#ifndef __RSCLOVER_H__
#define __RSCLOVER_H__

/*
 * ame-specific macros
 */

#define min(a,b)      (a<b)?a:b
#define max(a,b)      (a>b)?a:b

#define MAX_SEQ_LENGTH 1e6

/*
 * Define ame constants
 */

#define MEME_FORMAT 1 	// input format is meme

#define UNIFORM_BG 0    //uniform.
#define MOTIF_BG 1 		// background frequencies taken from motif file
#define FILE_BG 2		// background frequencies taken from specified file

#define QUICK_RS 0 //use my bodgy way of doing the ranksum
#define BETTER_RS 1 //use fabian's correct way of doing the ranksum

#define POS_FL 0  //use the fluorescence sort as positive indicator
#define POS_PWM 1 //use the PWM sort as positive indicator

#define AVG_ODDS 0
#define MAX_ODDS 1
#define SUM_ODDS 2
#define TOTAL_HITS 3

#define MOTIF_PWM 0    //use each motif normally
#define MOTIF_FALSE_REGEXP 1 //convert each motif into a regexp
#define MOTIF_REGEXP 1 //These are actual regular expressions -- same treatment as above

#define RANKSUM_METHOD 0 //use the ranksum test
#define FISHER_METHOD  1 //use the fisher exact test
#define MULTIHG_METHOD  2 //use the fisher exact test modified to label with 0,1, or 2.
#define LONG_MULTIHG_METHOD  3 //use the fisher exact test modified to label with 0,1,2, or 3.
#define LINREG_METHOD 4 //use linear regression test to minimise MSE.
#define SPEARMAN_METHOD 5 //use the spearman rank correlation co-efficient to calculate a score.

/*
 * Struct definitions
 */

typedef struct  {
	char* bg_filename; //file to get base freqs from
    char* outputdir; // where to send outputs
	int bg_format; //whether it's fasta, meme, etc.
	int motif_format; //is the motif file meme or tamo
	char** motif_filenames; //filenames of the motif library
    int number_motif_files; //how many files comprise motif library
	char* sequence_filename; //input sequences in fasta format
    char* commandline; // command line with path stripped from program name
	float pseudocount; //add to the motif frequency counts
	int scoring; //AVG_ODDS or MAX_ODDS
	int verbose;
	int rs_method; //QUICK_RS or BETTER_RS
	int positive_list; //POS_FL or POS_PWM
	int motif_mode; //MOTIF_PWM or MOTIF_REGEXP
	int pvalue_method; //RANKSUM_METHOD or FISHER_METHOD
	double fisher_pwm_threshold;
	double fisher_fl_threshold;
	double pvalue_threshold; //threshold for the mhg test with pwms.
	double pvalue_report_threshold; //threshold for reporting a motif in output.
	BOOLEAN_T length_correction;
	BOOLEAN_T log_fscores;
	BOOLEAN_T log_pwmscores;
	BOOLEAN_T linreg_normalise;
	BOOLEAN_T linreg_switchxy;
    BOOLEAN_T clobber; // TRUE if we can replace existing output directory
    BOOLEAN_T silent; // suppress version and command line output; for testing: discourage in normal use
	char* linreg_dump_dir;
	int fix_partition;
    // derived from command line:
    FILE *text_output;
    FILE *html_output;
} rsc_arg_t;

typedef struct {
	char* filename;
	MOTIF_T* motifs;
    PSSM_T* pssms;
	ARRAY_T** pv_lookup;
	int num;
	ARRAY_T* bg_freqs;
	BOOLEAN_T has_reverse_strand;
} rsc_motifs_t;

typedef struct {
	char* motif_id;
	char* motif_id_alt; // alternative motif name
	char* motif_url; // URL or NULL
	int split;
	double pleft;
	double pright;
	double pboth;
	double u;
} rsc_result_t;

typedef struct {
	int f_rank;
	int pwm_rank;
	double pwm_score;
	double f_score;
} rsc_rank_t;


// Nucleotide alphabet order as in motif.h
extern char alphabet[];

/*
 *  Function definitions
 */

void rsc_getopt(int argc, char *argv[]);
const char* rsc_get_usage();
ARRAY_T* rsc_load_background();
ARRAY_T* rsc_get_bg_freqs(SEQ_T* seq);
void rsc_load_motifs();
void rsc_scan_sequences();
void rsc_get_scores();
void rsc_usage();
void rsc_terminate(int status);
int rsc_compare_doubles (const double *a, const double *b);
int rsc_compare_scores (const void *a, const void *b);
int rsc_compare_mse (const void *a, const void *b);
int rsc_compare_ranks_f_rank (const void *a, const void *b);
int rsc_compare_ranks_pwm_score (const void *a, const void *b);
void rsc_dump_motif(MOTIF_T* m);
void rsc_convert_motif_regexp(MOTIF_T* m);
void rsc_convert_motifs();
rsc_result_t* rsc_do_ranksum_test(rsc_rank_t** rankings, int motif_index);
rsc_result_t* rsc_do_fisher_test(rsc_rank_t** rankings, int motif_index);
rsc_result_t* rsc_do_multihg_test(rsc_rank_t** rankings, int motif_index);
rsc_result_t* rsc_do_linreg_test(rsc_rank_t** rankings, int motif_index);
rsc_result_t* rsc_do_spearman_test(rsc_rank_t** rankings, int motif_index);
unsigned long long choose(unsigned n, unsigned k);
long double logchoose(unsigned n, unsigned k);
double* rsc_init_fisher_factorials(int len);
void shuffle(void** array, size_t n);
double rsc_bonferroni_correction(double, double);

double rsc_score_sequence(MOTIF_T* motif, SEQ_T* seq, int scoring);
MATRIX_T* rsc_convert_odds_matrix(MOTIF_T* motif);

#else
#endif
