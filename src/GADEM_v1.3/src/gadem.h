#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// top-ranked k-mers
typedef struct ktuples_info Ktuples;
struct ktuples_info {
   char *seq;
   double z,eCn,stdCn,p;
   int count;
};

// spaced dyad
typedef struct word_info Words;
struct word_info {
   char **s1;
   double *prob_sta,*prob_end;
   int count;
};

// GA "chromosomes"
typedef struct chromosome_info Chrs;
struct chromosome_info {
   int wordID,wordGroup;
};

// binding sites
typedef struct site_info Sites;
struct site_info {
   int    seq,pos,score;
   char   rev;
   double llr,pvalue;
};

// GA fitness
typedef struct fitness_score Fitness;
struct fitness_score {
   double value;
   int    index;
};

// GA roulette wheel
typedef struct roulette Wheel;
struct roulette {
   double start,end;
   int index;
};

// consensus
typedef struct col_cons ColCons;
struct col_cons {
  double v;
  int    i;
};

//pwm llr score distribution
typedef struct llr_distr Pgfs;
struct llr_distr{
     int score;
     double prob;
};

// background
typedef struct background_model BACKGROUND_Model;
struct background_model {
   double *monomerFreq,*dimerFreq,  *trimerFreq, *tetramerFreq,*pentamerFreq,*hexamerFreq,*heptamerFreq,*octamerFreq,*nonamerFreq;
   double              *transition1,*transition2,*transition3, *transition4, *transition5,*transition6, *transition7,*transition8;
   char **monomer,**dimer,**trimer,**tetramer,**pentamer,**hexamer,**heptamer,**octamer,**nonamer;
};

int *alloc_int (int );
int **alloc_int_int(int ,int );
int ***alloc_int_int_int(int ,int ,int );
int Compare_fitness(const void *, const void *);
int Compare_double(const void *, const void *);
int Compare_M1(const void *, const void *);
int Compare_z(const void *, const void *);
int Compare_score(const void *, const void *);
int Compare_llr(const void *, const void *);
int Compare_vector(const void *, const void *);
int pwm_score_dist(int **,int ,Pgfs *,double *);
int scan_em_seq_ptable(Pgfs *,int ,Sites *,int ,char **,char **,int *,int **,int ,int ,double *,char *);
int scan_llr_pgf(Pgfs *,int ,Sites *,int ,char **,char **,int *,int **,int ,int ,double *);
int scan_llr_empirical(Sites *,int ,char **,char **,int *,double **,int ,double *,double **,double **,double ,double *,int ,int );
int determine_cutoff(Pgfs *,int ,double );
int check_pwm_uniqueness_dist(double ***,int *,int ,Fitness *,double ,double ,char *,int );
int top_kmer(Words *,Ktuples *,int *,int );
int read_pwm0(char *,double **);
int read_userBackgModel(char *,BACKGROUND_Model *);
int score_backg_seq(double *,int ,char **,int *,double **,int ,double *,double **) ;
int llr_score(double *,int ,char **,int *,double **,int ,double *,double **);
int ini_M(int ,Pgfs *,int **,double *);
int prod_M(Pgfs *,int ,Pgfs *,int );
int *count_nucleotides(char **,char **,int ,int *,char **,int ,int );
int *count_nucleotides_backg(char **,char **,int ,int *,char **,int ,int );

char *alloc_char(int );
char **alloc_char_char(int ,int );
char **read_seq(int *,int *,char **,int ,int ,double *,char *);

void sgenrand(unsigned long long );
void reverse_seq(char **,char **,int ,int *);
void roulett_wheel_fitness(Fitness *,int ,Wheel *);
void sort_fitness(Fitness *,int );
void sort_double(double *,int );
void sort_llrDist(Pgfs *M1,int size);
void sort_kmer_z(Ktuples *,int );
void sort_sites_score(Sites *,int );
void sort_sites_llr(Sites *,int );
void sort_vector(ColCons *,int );
void sort_llrDist(Pgfs *,int );
void observed_pwm(Sites *,int ,char **,char **,int ,int ,double **);
void mask_sites(int ,char **,char **,int *,Sites *,int );
void extend_alignment(Sites *,int ,char **,char **,int *,int ,int ,int *);
void align_sites_count(Sites *,char **,char **,int ,int ,double **) ;
void initialisation(Chrs **,int ,int ,Words *,int ,int ,double *);
void log_pwm(double **,double **,int );
void log_ratio_to_int(double **,int **,int ,double *);
void pwm_double_to_int(double **,int **,int );
void construct_pwm(double **,double **,double **,char **,char **,int *,int ,int ,char *);
void ll_score_motif_model(int ,char **,char **,int *,double **,int ,double **,double **,char *,double *);
void normalize(double **,double **,int *,int ,int ,char *,int ,double **,int );
void copy_pwm_int_to_double(double **,double **,int );
void copy_pwm(double **,double **,int );
void dyad_to_pwm(Words *,int ,Chrs **,double ***,int *);
void destroy_word(Words *word,int );
void destroy_ktuples(Ktuples *,int );
void crossover(Chrs **,int ,Words *,int ,int ,Wheel *,int ,Fitness *,char *,double *,double );
void mutation(Chrs **,int ,Words *,int ,int ,Wheel *,int ,Fitness *,char *,double *,double );
void count_k_tuples(char **,char **,int ,int *,char **,int ,int ,int *);
void enumerate_kmers(char **,char **,char **);
int word_for_dyad(Words *,char **,char **,int ,int *,double *,int *,int *,int *);
void sample_without_replacement(char *,int ,int );
void sample_without_replacement2(int *,int ,int );
void print_result_2(Sites *,int ,int ,char **,char **,int *,char **,double ,double **,int ,int ,char *,char *,int ,double ,double ,FILE *,FILE *);
void print_motif(Sites *,int ,char **,char **,int *,int ,int ,double **);
void pwm_profile(double **,int ,char *);
void consensus_pwm(double **,int ,char *);
void effect_seq_length(char **,int ,int *,char *,int *);
void effect_seq_length_full(char **,int ,int *,int *);
void reassign_fitness(int ,char *,int ,Fitness *);
void standardize_pwm(double **,int );
void score_kmers(Ktuples *,double *,int ,int *,char **,int );
void assign_weight_uniform(int *,int ,double **);
void assign_weight_triangular(int *,int ,double **);
void assign_weight_normal(int *,int ,double **);
void assign_weight_triangular_uniform(int *,int ,double **,int );
void assign_weight_rectangle(int *,int ,double **,int );
void generate_background(int ,char **,char **,int *,BACKGROUND_Model *,int );
void ll_score_backg_model(int ,char **,double **,int *,int ,BACKGROUND_Model *,int );
void marginal_prob(int *count,int numKmer,double *freq);
void simulate_backg_seq_bmodel(BACKGROUND_Model *,int ,int ,int *,char **);
void simulate_background_seq(double *,int ,int *,char **);
void select_high_scoring_seq_for_EM (double *,int ,int ,char *,double );
void transition_1st(double *,double *);
void transition_2nd(double *,double *);
void transition_3rd(double *,double *);
void transition_4th(double *,double *);
void transition_5th(double *,double *);
void transition_6th(double *,double *);
void transition_7th(double *,double *);
void transition_8th(double *,double *);
void compute_freq(int *,int ,double *);
void numerate_monomer_to_pentamer(BACKGROUND_Model *);

double *alloc_double(int );
double **alloc_double_double(int ,int );
double ***alloc_double_double_double(int ,int ,int );
double alignment(Sites *,int ,char **,char **,int ,int );
double meme_e_value(double *,double **,int ,int ,int ,int *);
double p_value(int ,Pgfs *,int );
double E_value(double **, int , double *, int , int , int *);
double vector_similarity(void);
double check_convergence(double **,double **,int );
double find_pvalue(int ,Pgfs *,int );
double *base_frequency(int ,char **, int *);

Chrs **alloc_chrs(int ,int );
Words *alloc_word(int ,int );
Sites *alloc_site(int );
Sites **alloc_site_site(int size1,int size2);
Wheel *alloc_wheel(int );
Fitness *alloc_fitness(int );
Ktuples *alloc_ktuples(int ,int );
Pgfs *alloc_distr(int size);
BACKGROUND_Model *alloc_background(void);

