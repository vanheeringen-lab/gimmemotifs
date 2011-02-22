/***********************************************************************
*                                                                      *
*       MetaMEME						       *
*	Date:		23/4/2002				       *
*       Author: 	Timothy L. Bailey                              *
*	Copyright:	University of Queensland		       *
*                                                                      *
***********************************************************************/

#include "matrix.h"
#include "karlin.h"
#include "rdb-matrix.h"
#include "string.h"
#include "subst-matrix.h"
#include "alphabet.h"

#ifndef SUBST_MATRIX_DEBUG
#define SUBST_MATRIX_DEBUG 0
#endif

#define EPSILON 1e-200

/* dayhoff PAM 1 (matrix; order of alphabet: ACDEFGHIKLMNPQRSTVWY */
/* dayhoff_ij = Pr(amino acid j --> amino acid i | time=1)
              = Pr(i | j, t=1) */
char *pam_prot_alphabet = "ACDEFGHIKLMNPQRSTVWY";
double dayhoff[20][20] = {
  { 9867, 3, 10, 17, 2, 21, 2, 6, 2, 4, 6, 9, 22, 8, 2, 35, 32, 18, 0, 2},
  { 1, 9973, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 5, 1, 2, 0, 3},
  { 6, 0, 9859, 53, 0, 6, 4, 1, 3, 0, 0, 42, 1, 6, 0, 5, 3, 1, 0, 0},
  { 10, 0, 56, 9865, 0, 4, 2, 3, 4, 1, 1, 7, 3, 35, 0, 4, 2, 2, 0, 1},
  { 1, 0, 0, 0, 9946, 1, 2, 8, 0, 6, 4, 1, 0, 0, 1, 2, 1, 0, 3, 28},
  { 21, 1, 11, 7, 1, 9935, 1, 0, 2, 1, 1, 12, 3, 3, 1, 21, 3, 5, 0, 0},
  { 1, 1, 3, 1, 2, 0, 9912, 0, 1, 1, 0, 18, 3, 20, 8, 1, 1, 1, 1, 4},
  { 2, 2, 1, 2, 7, 0, 0, 9872, 2, 9, 12, 3, 0, 1, 2, 1, 7, 33, 0, 1},
  { 2, 0, 6, 7, 0, 2, 2, 4, 9926, 1, 20, 25, 3, 12, 37, 8, 11, 1, 0, 1},
  { 3, 0, 0, 1, 13, 1, 4, 22, 2, 9947, 45, 3, 3, 6, 1, 1, 3, 15, 4, 2},
  { 1, 0, 0, 0, 1, 0, 0, 5, 4, 8, 9874, 0, 0, 2, 1, 1, 2, 4, 0, 0},
  { 4, 0, 36, 6, 1, 6, 21, 3, 13, 1, 0, 9822, 2, 4, 1, 20, 9, 1, 1, 4},
  { 13, 1, 1, 3, 1, 2, 5, 1, 2, 2, 1, 2, 9926, 8, 5, 12, 4, 2, 0, 0},
  { 3, 0, 5, 27, 0, 1, 23, 1, 6, 3, 4, 4, 6, 9876, 9, 2, 2, 1, 0, 0},
  { 1, 1, 0, 0, 1, 0, 10, 3, 19, 1, 4, 1, 4, 10, 9913, 6, 1, 1, 8, 0},
  { 28, 11, 7, 6, 3, 16, 2, 2, 7, 1, 4, 34, 17, 4, 11, 9840, 38, 2, 5, 2},
  { 22, 1, 4, 2, 1, 2, 1, 11, 8, 2, 6, 13, 5, 3, 2, 32, 9871, 9, 0, 2},
  { 13, 3, 1, 2, 1, 3, 3, 57, 1, 11, 17, 1, 3, 2, 2, 2, 10, 9901, 0, 2},
  { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 9976, 1},
  { 1, 3, 0, 1, 21, 0, 4, 1, 0, 1, 0, 3, 0, 0, 0, 1, 1, 1, 2, 9945}
};

/* transition/transversion PAM 1 matrix; alphabet order "ACGT" */
char *pam_dna_alphabet = "ACGT";
double trans[4][4] = {
  { 9900, 20, 60, 20}, 
  { 20, 9900, 20, 60},
  { 60, 20, 9900, 20},
  { 20, 60, 20, 9900}
};

/* PAM frequencies for DNA */
double pam_dna_freq[] = { 
    0.25 /* A */,
    0.25 /* C */,
    0.25 /* G */,
    0.25 /* T */
};

/* PAM frequencies for proteins; supplied by S. Altschul */
double pam_prot_freq[] = { 
    0.087 /* A */,
    0.033 /* C */,
    0.047 /* D */,
    0.050 /* E */,
    0.040 /* F */,
    0.088 /* G */,
    0.034 /* H */,
    0.037 /* I */,
    0.081 /* K */,
    0.085 /* L */,
    0.015 /* M */,
    0.040 /* N */,
    0.051 /* P */,
    0.038 /* Q */,
    0.041 /* R */,
    0.070 /* S */,
    0.058 /* T */,
    0.065 /* V */, 
    0.010 /* W */,
    0.030 /* Y */,
};

/**********************************************************************/
/*
	reorder_matrix

	Reorder a matrix from the alphabet order in alpha1
	to the alphabet order of alpha2.  The letters in alph1 must be
	a superset of those in alpha2.  Rows and columns corresponding
	to letters missing from alpha2 are discarded. 

	Returns a reordered square matrix of size length(alpha2).
*/
/**********************************************************************/
MATRIX_T *reorder_matrix(
  char *alpha1,				/* current alphabet */
  char *alpha2,				/* new alphabet; must be subset */
  MATRIX_T *in_matrix			/* matrix to reorder */
)
{
  int i, j;
  int alen1 = strlen(alpha1);
  int alen2 = strlen(alpha2);
  MATRIX_T *out_matrix;

  if (alen2 > alen1) 
    die("The new alphabet %s must be a subset of the old alphabet %s.\n", alpha2, alpha1);

  out_matrix = allocate_matrix(alen2, alen2);
  for (i=0; i<alen2; i++) {
    int ii = strchr(alpha1, alpha2[i]) - alpha1;
    for (j=0; j<alen2; j++) {
      int jj;
      char *ptr = strchr(alpha1, alpha2[j]);
      if (!ptr)
        die("The new alphabet %s must be a subset of the old alphabet %s\n", alpha2, alpha1);
      jj = ptr - alpha1;
      set_matrix_cell(i, j, get_matrix_cell(ii, jj, in_matrix), out_matrix);
    }
  }
  return(out_matrix);
} /* reorder_matrix */

/**********************************************************************/
/*
       	gen_pam_matrix

	Generate a PAM joint probability or log-odds score matrix.  

	If logodds is false, generates the joint probability matrix:
	  M_ij = Pr(i,j | time = dist).
        If logodds is true, generates a log-odds scoring matrix:
	  M_ij = factor * log(Pr(i,j | time = dist)/ Pr(i)Pr(j)),
	where 
          factor = 2/log(2), 2 <= dist < 170, 
	           3/log(2), 170 <= dist.

	This is adapted from MEME code.  I'm not sure where I got
	the original conditional probabilities--I think from a GCG file.
	The log-odds matrices generated are close, but not exactly the
	same as those generated by "pam" Version 1.0.6 that BLAST and
	everyone else seems to have used to create PAM files, 
	and I have not been able to locate a copy of that program. 
	I have no idea where I got the "transversion/transition matrix"
	that is used for making DNA "pam" matrices.

	Returns the matrix.
*/
/**********************************************************************/

extern MATRIX_T* gen_pam_matrix(
  BOOLEAN_T dna,		/* DNA alphabet if true; default: protein */
  int dist,			/* PAM distance */
  BOOLEAN_T logodds		/* true: generate log-odds matrix 
				   false: generate target frequency matrix 
				*/
)
{
  int i, j;
  MATRIX_T *matrix, *mul;
  double *pfreq = dna ? pam_dna_freq : pam_prot_freq;	/* standard frequencies */
  char *alpha = dna ? pam_dna_alphabet : pam_prot_alphabet;	/* standard alphabet */
  int alen = strlen(alpha);				/* length of standard alphabet */
  double factor = dist < 170 ? 2/log(2) : 3/log(2);	/* same as in "pam" Version 1.0.6 */

  /* create the array for the joint probability matrix */
  matrix = allocate_matrix(alen, alen);
  mul = allocate_matrix(alen, alen);

  /* initialize the matrix: PAM 1:
     due to roundoff, take the average of the two estimates of the joint frequency
     of i and j as the joint, then compute the conditionals for the matrix
  */
  for (i=0; i<alen; i++) {
    for (j=0; j<=i; j++) {
      double vij = dna ? trans[i][j] : dayhoff[i][j];
      double vji = dna ? trans[j][i] : dayhoff[j][i];
      double joint = ((vij * pfreq[j]) + (vji * pfreq[i]))/20000;/* use average to fix rndoff */
      set_matrix_cell(i, j, joint/pfreq[j], matrix);
      if (i!=j) set_matrix_cell(j, i, joint/pfreq[i], matrix);
    }
  }

  /* take PAM matrix to desired power to scale it */ 
  copy_matrix(matrix, mul);
  for (i=dist; i>1; i--) {
    MATRIX_T *product = matrix_multiply(matrix, mul);
    SWAP(MATRIX_T*, product, matrix)
    free_matrix(product);
  } 
  free_matrix(mul);

  /* convert to joint or logodds matrix:
     target:  J_ij = Pr(i,j) = Mij pr(j) 
     logodds: L_ij = log (Pr(i,j)/(Pr(i)Pr(j)) = log (Mij Pr(j)/Pr(i)Pr(j)) = log(Mij/pr(i)) 
  */
  for (i=0; i<alen; i++) {
    for (j=0; j<alen; j++) {
      double vij = get_matrix_cell(i, j, matrix);
      vij = logodds ? nint(factor * log((vij+EPSILON)/pfreq[i])) : vij * pfreq[j];
      set_matrix_cell(i, j, vij, matrix);
    }
  }

  return matrix;
} /* gen_pam_matrix */

/**********************************************************************/
/*
	read_score_matrix

	Read in a score matrix in the format used by BLAST. 

	Only reads in the standard letters in the alphabet.
*/
/**********************************************************************/
MATRIX_T *read_score_matrix(
  char *score_filename,			/* name of score file */
  char **alpha1				/* alphabet in score matrix */
)
{
  int i;
  char *alpha;				/* alphabet in file */
  int alen;
  FILE *score_file;
  RDB_MATRIX_T *rdb_matrix;

  /* open the score file */
  if (open_file(score_filename, "r", FALSE, "score", "substitution scores", &score_file) == 0) 
    exit(1);

  /* read in the score file */
  rdb_matrix = read_rdb_matrix(" ", FALSE, 0, FALSE, NULL, score_file);

  /* get alphabet */
  alen = get_num_strings(rdb_matrix->col_names);
  alpha = (char *)mm_malloc(sizeof(char) * (alen+1));
  for (i=0; i<alen; i++) alpha[i] = get_nth_string(i, rdb_matrix->col_names)[0];
  alpha[i] = '\0';
  *alpha1 = alpha;			/* return alphabet */

  return(rdb_matrix->matrix);
} /* read_score_matrix */

/**********************************************************************/
/*
	get_score_matrix

	Get a letter substitution scoring matrix.

	If a filename is given, the scoring matrix is read from the file.
	Otherwise, a PAM matrix of the given type and  evolutionary 
	distance is generated.

	The matrix is reordered to conform to the standard alphabet.  This
	alphabet must be a subset of the alphabet in the file (if the
	matrix was read from a file).  Any rows and columns for letters
	not in the standard alphabet are discarded.

	Returns the matrix.
*/
/**********************************************************************/
MATRIX_T *get_score_matrix(
  char *score_filename,		/* name of score file */
  BOOLEAN_T dna,		/* DNA alphabet if true; default: protein */
  int dist,			/* PAM distance (ignored if filename != NULL) */
  char *alpha			/* standard alphabet for final matrix */
)
{
  char *alpha1;			/* initial alphabet */
  MATRIX_T *tmp;		/* temporary score matrix */
  MATRIX_T *matrix;		/* score matrix */

  if (score_filename) {		/* read matrix from file */
    tmp = read_score_matrix(score_filename, &alpha1);
  } else {			/* generate PAM matrix */
    tmp = gen_pam_matrix(dna, dist, TRUE);
    alpha1 = dna ? pam_dna_alphabet : pam_prot_alphabet;
  }

  /* reorder the matrix to standard alphabet */
  matrix = reorder_matrix(alpha1, alpha, tmp);
  free_matrix(tmp);
 
  return(matrix); 
} /* get_score_matrix */


/**********************************************************************/
/*
	make_karlin_input

	Prepare the input required for karlin() from a scoring matrix
	and a letter frequency distribution.
*/
/**********************************************************************/
KARLIN_INPUT_T *make_karlin_input(
  MATRIX_T *matrix,			/* scoring matrix */
  ARRAY_T *probs			/* letter freq distribution */
)
{
  int i, j;
  double escore;
  long lowest, highest;
  ARRAY_T *score_probs; 
  int nscores;
  int alen = get_num_rows(matrix);	/* size of alphabet */
  KARLIN_INPUT_T *karlin_input;		/* data to return */

  /*  find the highest and lowest scores in the scoring matrix */
  lowest = 1;
  highest = -1;
  for (i=0; i<alen; i++) {
    for (j=0; j<alen; j++) {
      double s = get_matrix_cell(i, j, matrix);
      if (s < lowest) lowest = s;
      if (s > highest) highest = s;
    }
  }
  if (lowest >= 0) die("Lowest score in scoring matrix must be negative, is %f.", (double)lowest);
  if (highest<= 0) die("Highest score in scoring matrix must be positve, is %f.", (double)highest);

  /* allocate the array of score probabilities and set to 0 */
  nscores = highest - lowest + 1;
  score_probs = allocate_array(nscores);
  init_array(0, score_probs);
  
  /* compute the probabilities of different scores */ 
  escore = 0;
  for (i=0; i<alen; i++) {
    for (j=0; j<alen; j++) {
      int s = get_matrix_cell(i, j, matrix);
      double pi = get_array_item(i, probs);
      double pj = get_array_item(j, probs);
      double sp = get_array_item(s-lowest, score_probs); 
      set_array_item(s-lowest, sp + pi*pj, score_probs);	/* cumulative prob. of score */
      escore += pi*pj*s;
      /*printf("i %d j %d s %d pi %f pj %f sp %f escore %f\n",i,j,s, pi, pj, sp, escore);*/
    }
  }

  karlin_input = (KARLIN_INPUT_T *)mm_malloc(sizeof(KARLIN_INPUT_T));
  karlin_input->low = lowest;
  karlin_input->high = highest;
  karlin_input->escore = escore;
  karlin_input->prob = score_probs;

  return(karlin_input);
} /* make_karlin_input */


/**********************************************************************/
/*
	convert_score_to_target

	Convert a scoring matrix to a target matrix.
*/
/**********************************************************************/
MATRIX_T *convert_score_to_target(
  MATRIX_T *score,			/* score matrix */
  ARRAY_T *prob				/* letter frequencies */
)
{
  int i, j;
  KARLIN_INPUT_T *karlin_input;
  double lambda, K, H;
  MATRIX_T *target;			/* target freq. matrix */
  int alen = get_num_rows(score);	/* alphabet length */

  /* make input for karlin() */
  karlin_input = make_karlin_input(score, prob);
  
  /* get lambda */
  karlin(karlin_input->low, karlin_input->high, karlin_input->prob->items,
    &lambda, &K, &H);
  /*printf("lambda %f K %f H %f\n", lambda, K, H);*/

  /* calculate target frequencies */
  target = allocate_matrix(alen, alen);
  for (i=0; i<alen; i++) {
    for (j=0; j<alen; j++) {
      double pi = get_array_item(i, prob);
      double pj = get_array_item(j, prob);
      double sij = get_matrix_cell(i, j, score);
      double f = pi * pj * exp(lambda * sij);
      set_matrix_cell(i, j, f, target);
    }
  }

  // Free local dynamic memory.
  free_array(karlin_input->prob);
  myfree(karlin_input);

  return(target);
} /* convert_score_to_target */

/**********************************************************************/
/*
       	get_subst_target_matrix

	Get a substitution target frequency matrix.
*/
/**********************************************************************/
MATRIX_T *get_subst_target_matrix(
  char *score_filename,		/* name of score file */
  BOOLEAN_T dna,		/* DNA alphabet if true; default: protein */
  int dist,			/* PAM distance (ignored if score_filename != NULL) */
  char *alpha,			/* standard alphabet for final matrix */
  ARRAY_T *back			/* background frequencies of standard alphabet */
)
{
  MATRIX_T *score;		/* score matrix */
  MATRIX_T *target;		/* target frequency matrix */

  score = get_score_matrix(score_filename, dna, dist, alpha);
  target = convert_score_to_target(score, back);

  if (SUBST_MATRIX_DEBUG)
  {
    int i, j, alength=strlen(alpha);
    double sum;

      if (score_filename) {
	printf("From file %s\n", score_filename);
      } else {
	printf("Generated PAM %d\n", dist);
      }
      printf("%6c ", ' ');
      for (i=0; i<alength; i++) {
	printf("%6c ", alpha[i]);
      }
      printf("\n");
    sum = 0;
    for (i=0; i<alength; i++) {
      printf("%6c ", alpha[i]);
      for (j=0; j<alength; j++) {
	double x = get_matrix_cell(i,j,score);
	sum += x;
	printf("%6.4f ", x);
      }
      printf("\n");
    }
    printf("sum of entries = %f\n", sum);
  }

  free_matrix(score);
    
  return(target);
} /* get_subst_target_matrix */


/**************************************************************************
 * Get pseudocount frequencies.
 *
 * The target_freq matrix only has values for the basic alphabet.
 * Fill in the ambiguous character pseudocounts afterwards using
 * the average of pseudocounts for letters matching the ambiguous ones.
 **************************************************************************/
ARRAY_T *get_pseudocount_freqs(
   ARRAY_T *	  f,		/* Foreground distribution. */
   ARRAY_T *      b,		/* Background distribution. */
   MATRIX_T *     target_freq	/* Target frequency matrix. */
)
{
  int i, j;
  int alph_size = get_alph_size(ALL_SIZE);	// includes ambigs
  int asize = get_alph_size(ALPH_SIZE);		// excludes ambigs
  ARRAY_T *g = allocate_array(alph_size);

  /*
    Create pseudocount frequencies.
  */
  for (i = 0; i < asize; i++) {				/* non-ambiguous freqs */
    double gi = 0;
    for (j= 0; j < asize; j++) {			/* non-ambiguous freqs */
      double qij = get_matrix_cell(i, j, target_freq);
      double fj = get_array_item(j, f);
      double bj = get_array_item(j, b);
      gi += (fj/bj) * qij;
    } /* j */
    set_array_item(i, gi, g);
    if (SUBST_MATRIX_DEBUG) printf("%g %g, ", get_array_item(i, f), gi);
  } /* i */
  fill_in_ambiguous_chars(FALSE, g);			/* takes the average pseudocount */
  if (SUBST_MATRIX_DEBUG) printf("\n");

  return(g);						/* return the pseudocounts */
} /* get_pseudocount_freqs */

#ifdef SO
/***********************************************************************/
/*
	Testing routine.
	
	subst_matrix <0|1> <dist>

	Specify 0 for protein, 1 for DNA.
	Specify the integer PAM distance.
*/
/***********************************************************************/
VERBOSE_T verbosity = NORMAL_VERBOSE;
main(int argc, char **argv) {
  int i, j, alength;
  int dist = 0;
  BOOLEAN_T dna = FALSE; 
  char *score_filename = NULL;
  char *alpha;
  MATRIX_T *matrix;
  ARRAY_T *probs;
  double *freqs;
  KARLIN_INPUT_T *karlin_input;
  int nscores;
  double sum;
  char usage[1000] = "";

  // Define the usage message.
  strcat(usage, "USAGE: subst_matrix [options] <score file>\n");
  strcat(usage, "\n");
  strcat(usage, "   Options:\n");
  strcat(usage, "     --dna\n");
  strcat(usage, "     --dist <float>\n");
  strcat(usage, "\n");

  // Parse the command line.
  while (1) { 
    int c;
    int option_index = 0;
    const char* option_name;

    // Define command line options.
    static struct option long_options[] = {
      {"dna", 0, 0, 0},
      {"dist", 1, 0, 0},
    };

    // Read the next option, and break if we're done.
    c = getopt_long_only(argc, argv, "+", long_options, &option_index);
    if (c == -1) {
      break;
    } else if (c != 0) {
      die("Invalid return from getopt (%d)\n", c);
    }

    // Get the option name (we only use long options).
    option_name = long_options[option_index].name;
    if (strcmp(option_name, "dna") == 0) {
      dna = TRUE;
    } else if (strcmp(option_name, "dist") == 0) {
      dist = atoi(optarg);
    } else {
      die("Invalid option (%s).\n", option_name);
    }
  }

  // Read the single required argument.
  if (optind + 1 != argc) {
    fprintf(stderr, usage);
    exit(1);
  }
  score_filename = argv[optind];



  alpha = dna ? "ACGT" : "ACDEFGHIKLMNPQRSTVWY";
  alength = strlen(alpha);

  /* background frequencies */
  probs = allocate_array(alength);
  freqs = dna ? pam_dna_freq : pam_prot_freq;
  fill_array(freqs, probs);			/* copy freqs into ARRAY_T */

  if (dist > 1) {
    printf("From gen_pam_matrix:\n");
    matrix = gen_pam_matrix(dna, dist, FALSE);
    printf("%6c ", ' ');
    for (i=0; i<alength; i++) {
      printf("%6c ", alpha[i]);
    }
    printf("\n");
    sum = 0;
    for (i=0; i<alength; i++) {
      printf("%6c ", alpha[i]);
      for (j=0; j<alength; j++) {
	double x = get_matrix_cell(i,j,matrix);
	sum += x;
	printf("%6.4f ", x);
      }
      printf("\n");
    }
    printf("sum of entries = %f\n", sum);
  }

  printf("From get_subst_target_matrix:\n");
  matrix = get_subst_target_matrix(score_filename, dna, dist, alpha, probs);
} /* main */
#endif
