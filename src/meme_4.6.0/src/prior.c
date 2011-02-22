/*
 * $Id: prior.c 2120 2007-10-11 01:36:51Z cegrant $
 * 
 * $Log$
 * Revision 1.3  2005/10/25 19:06:39  nadya
 * rm old macro for Header, all info is taken care of by Id and Log.
 *
 * Revision 1.2  2005/10/06 22:51:02  nadya
 * SUN compiler can't reduce expression inside define.
 * redefine array size to be an int from start.
 *
 * Revision 1.1.1.1  2005/07/29 17:24:18  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

#include "meme.h"
#define MAXS	200
static BOOLEAN first_time = TRUE;

#define LogAddLog(x, y) ((x) > (y) ? LogAddLog1((x),(y)) : LogAddLog1((y),(x)))
#define LogAddLog1(x,y) ((x)-(y) > BITS ? (x) : (x) + log(1+exp((y)-(x))))
	
/* ---------------------------------------------------------------------- */
PriorLib *alloc_PriorLib(int l, int Alpha)
/*
        L       Number of distributions
        Alpha   Number of alphabet characters
*/
{
        PriorLib        *temp;
	int		i;

        temp = (PriorLib *)mymalloc(sizeof(PriorLib));
        temp->L = l;
        temp->AlphaChar = Alpha;

        temp->Mix = (Real *)mymalloc(sizeof(Real)*l);
        temp->B = (Real *)mymalloc(sizeof(Real)*l);
        temp->FullUpdate = (int *)mymalloc(sizeof(int)*l);
        temp->QUpdate = (int *)mymalloc(sizeof(int)*l);

        temp->StructID = (char **)mymalloc(sizeof(char *)*l);
        temp->Comment = (char **)mymalloc(sizeof(char *)*l);
        temp->Distr = (Real **)mymalloc(sizeof(Real *)*l);
        for (i=0; i<l; i++) {
           temp->Distr[i] = (Real *)mymalloc(sizeof(Real)*(Alpha+1));
           temp->StructID[i] = (char *)mymalloc(sizeof(char)*MAXS);
           temp->Comment[i] = (char *)mymalloc(sizeof(char)*MAXS);
        } /* endfor */

        return(temp);
}

/* ---------------------------------------------------------------------- */
PriorLib *read_PriorLib(char *plib_name, double desired_beta)
/*
	plib_name	name of prior library file
	desired_beta	>0, scale \beta_{i,j} so
				\sum_{i=0}^L \lambda_i \sum_{j=1}^20 \beta_{i,j}
			has this value
			=, don't scale prior
			< 0, just get alphabet
*/
{
        int             i,j, line=0;
        int             Alpha, l;
        PriorLib        *temp;
        char            input[MAXS], foo[MAXS], alphabet[MAXALPH+1], 
			checkstr[81], *token;
        Real x;
  	FILE *fp;

/* tlb */
	fp = fopen(plib_name, "r");
	if (!fp) {
	  fprintf(stderr, "Can't find prior library %s\n", plib_name);
	  exit(1);
	}

	token = "Alphabet="; line++;
        fscanf(fp,"%s %s\n", checkstr, alphabet);
        if (strcmp(checkstr, token)) {
	  fprintf(stderr, "Line %d of prior library file \n %s \n"
	    "should start with \"%s\" "
	    "but it starts with \"%s\".\n", line, plib_name, token, checkstr);
	  exit(1);
	}
	Alpha = strlen(alphabet);

	token = "NumDistr="; line++;
        fscanf(fp,"%s %d\n", checkstr, &l);
        if (strcmp(checkstr, token)) {
	  fprintf(stderr, "Line %d of prior library file \n %s \n"
	    "should start with \"%s\" "
	    "but it starts with \"%s\"\n.", line, plib_name, token, checkstr);
	  exit(1);
	}

        temp = alloc_PriorLib(l,Alpha);

	if (Alpha > MAXALPH) {
	  fprintf(stderr,
            "Alphabet size specified in prior library %s too big.\n"
	     "Change MAXALPH in user.h and remake meme.\n", plib_name);
	  exit(1);
	}
        strcpy(temp->alphabet, alphabet);
        temp->AlphaChar = Alpha;
        temp->L = l;

	if (desired_beta < 0) {
	  fclose(fp);
	  return(temp);
	}

        for (i=0; i < temp->L; i++)
        {

        /* Get rid of number= */
        fscanf(fp,"%*s %*s\n");

        /* Mixture */
        fscanf(fp,"%*s");
        fscanf(fp,"%lf\n", &x);
        temp->Mix[i] = x;

        /* B (strength) */
        fscanf(fp,"%*s");
        fscanf(fp,"%lf\n", &x);
        temp->B[i] = x;

        /* Alpha */
        temp->Distr[i][0] = temp->B[i];
        fscanf(fp,"%*s");
        for (j=1; j <= temp->AlphaChar; j++) {
                fscanf(fp,"%lg", &x);
                temp->Distr[i][j] = x * temp->B[i];
        }

        /* FullUpdate */
        fscanf(fp,"%*s");
        fscanf(fp,"%d\n", &(temp->FullUpdate[i]));

        /* QUpdate */
        fscanf(fp,"%*s");
        fscanf(fp,"%d\n", &(temp->QUpdate[i]));

        /* StructID */
        fgets(input, MAXS, fp);
        sscanf(input,"%s",foo);
	input[strlen(input)-1] = '\0';
        strcpy( (temp->StructID[i]), (input + strlen(foo)) );

        /* Comments */
        fgets(input, MAXS, fp);
        sscanf(input,"%s",foo);
        strcpy( (temp->Comment[i]), (input + strlen(foo)) );
        }

  /* tlb; scale beta to desired value */
  if (desired_beta > 0) {
    int i, j;
    double beta = 0;
    double scale;
    for (i=0; i<temp->L; i++) {
      beta += temp->Mix[i] * temp->B[i];
    }
    /*printf("beta = %10.6f\n", beta);*/
    scale = desired_beta/beta;
    for (i=0; i<temp->L; i++) {
      for (j=0; j<=temp->AlphaChar; j++) {
        temp->Distr[i][j] *= scale;
      }
    }
  }
	fclose(fp);

        return(temp);
}

/* ---------------------------------------------------------------------- */
extern void mixture_regularizer(
  double *freq, 		/* obs freq */
  PriorLib *Lib,		/* priors */ 
  double *reg			/* pseudo-counts */
)
{
	Real	f[MAXALPH+1], sum, tmp;
	int 	i,j;
	/*Real	logpajgy();*/

	/* Put frequencies into array with f[0] = sum f_i */
	sum=0.0;
	for (i=0; i< Lib->AlphaChar; i++)
	{
		sum += freq[i];
		f[i+1]=freq[i];
	}
	f[0]=sum;

	/* Calculate probs */
	logpajgy(f, Lib, 0, 1);

	/* Calculate new regularizer */
	for (i=0; i< Lib->AlphaChar; i++)
	{
		reg[i]=0.0;
		for (j=0; j< Lib->L; j++)
		{
			tmp = (exp(logpajgy(f, Lib, j, 0)))*
			      ((Lib->Distr[j])[i+1]); /* skip A0 */
			reg[i] += tmp; 
		}
	}
}

/* ---------------------------------------------------------------------- */

/* This function computes log(p(a^j|y)) used in the calculation of theta. 
   It is defined to be

\log(\frac{q_j p(y given \alpha^j)}{\sum_k q_k p(y given \alpha^k})

*/

Real	logpajgy(
		Real	 *y,		/* observed frequencies */
		PriorLib *Library, 	/* Library of priors */
		int	 j, 		/* j'th prior to examine */
	   	int	 Calc		/* if ==1 calculate probs */
		)
{
	int i;
	RealPrec tmp;
	static RealPrec	logprob[MAXS], logdenom;/* Holders for probabilities */

	/* Calculate log probs if not already done */
	if (Calc)
	{
		tmp = log(Library->Mix[0]) + logpygaj(y,Library->Distr[0],
			Library->AlphaChar);
		logdenom = tmp;
		logprob[0] = tmp;

		/* Do remaining terms */
		for (i=1; i < Library->L; i++)
		{
			tmp = (log(Library->Mix[i]) + 
			  logpygaj(y, Library->Distr[i], Library->AlphaChar));
		
			logdenom = LogAddLog(logdenom, tmp);
			logprob[i] = tmp;
		}
	}
	return(logprob[j] - logdenom);
}

/* ---------------------------------------------------------------------- */

/* This function computes log(p(y|a^j)) used in the calculation of theta. 
   It is defined to be

\log(\frac{\Gamma(n+1)\Gamma(\alpha_0)}{\Gamma(n+\alpha_0)}
\prod_{i=1}^{20}\frac{\Gamma(y_i+\alpha_i)}{\Gamma(y_i+1)\Gamma(\alpha_i)})

*/

#define MAXX1 100
#define DELTA1 .001
#define MAXX2 100000
#define DELTA2 1.0

#ifdef SunOS
/* some SUN compilers can't reduce expression, hence this modification.
   Change if MAXX/DELTA constants change */
   #define SIZE1 100000
   #define SIZE2 100000
#else
   #define SIZE1 (int)(MAXX1/DELTA1)
   #define SIZE2 (int)(MAXX2/DELTA2)
#endif

static double lgam_array1[SIZE1 + 2];
static double lgam_array2[SIZE2 + 2];
static double lgam(double x);

RealPrec	logpygaj(
			Real	*y,	/* observed frequencies */
			Real	 *a,	/* distribution parameters */
			int	AlphLength	/* length of alphabet */
			)
{
	int		i;
	RealPrec	temp;

  /* set up array of values of lgamma to save time */
  if (first_time) {
    double x;
    for (i=1, x=0; i<=MAXX1/DELTA1 + 1; i++) { 
      x += DELTA1;
      lgam_array1[i] = lgamma(x);
      lgam_array1[i] = lgam_array1[i];
    }
    for (i=1, x=0; i<=MAXX2/DELTA2 + 1; i++) { 
      x += DELTA2;
      lgam_array2[i] = lgamma(x);
      lgam_array2[i] = lgam_array2[i];
    }
    first_time = FALSE;
  }
	
	temp=0.0;

	temp+= lgamma(y[0]+1.0);
	temp+= lgamma(a[0]);
	temp+= -lgamma(y[0]+a[0]);

	for (i=1; i<=AlphLength; i++)
	{
		temp+= lgamma(y[i]+a[i]);
		temp+= -lgamma(y[i]+1.0);
		temp+= -lgamma(a[i]);
	}
	
	return(temp);

}

static double lgam(double x)
{
  if (x >= DELTA1 && x <= MAXX1) { 
    int i = (int) (x/DELTA1);
    return(lgam_array1[i] + (lgam_array1[i+1] - lgam_array1[i])/2);
  } else if (x > MAXX1 && x <= MAXX2) { 
    int i = (int) (x/DELTA2);
    return(lgam_array2[i] + (lgam_array2[i+1] - lgam_array2[i])/2);
  } else {
    return(lgamma(x));
  }
} /* lgam */

