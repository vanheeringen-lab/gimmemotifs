/*
 * $Id: init.c 4808 2010-08-19 00:31:32Z cegrant $
 *
 */

/***********************************************************************
*                                                                      *
* MEME                                                                 *
* Copyright 1994-2007, The Regents of the University of California     *
* Author: Timothy L. Bailey                                            *
*                                                                      *
***********************************************************************/
/* init.c */
/*
  Initialize meme.
*/

#define ABS_MIN_W 2

#include "meme.h"
#include "general.h"
#include "banner.h"
#include "dir.h"
#include <sys/types.h>
#include <unistd.h>
#include <err.h>
#include "io.h"
#include "psp.h"

#ifndef EXP
#define EXP 0
#else
#define EXP 1
#endif

/* priors */
#define PROTEIN_PLIB "prior30.plib"
#define DNA_PLIB "dna.plib"

#define ROUNDERROR (1E-12)

/* User input parameters */
static BOOLEAN check_syntax = FALSE;  /* exit after checking syntax if true */
static char *datafile = NULL; /* positive examples */
static char *sf = NULL;   /* name to print for datafile */
static char *negfile = NULL;  /* negative examples */
static char *ntype = "pair";  /* output two matrices per motif if negatives */
static char *obj = "ev";  /* objective function */
static char *bfile = NULL;  /* use default background Markov model file */
static char *pspfile = NULL; /* use positional priors */
static BOOLEAN psp2 = FALSE;  /* true if 2-stranded positional priors */
static char *default_output_dirname = "meme_out";  /* default name of output
                                                   directory */
static BOOLEAN clobber = FALSE; /* default is not to overwrite existing files */
static char *plib_name = NULL;  /* use default library */
static char *mod = "zoops"; /* model type input string; default ZOOPS */
static char *alph = "PROTEIN";  /* default alphabet IUPAC protein 1-letter */
static BOOLEAN revcomp = FALSE; /* don't use reverse complement strand of DNA */
static int pal = 0;   /* = 0, no palindromes
           = 1, force DNA palindromes,
        */
static BOOLEAN ma_trim = TRUE;  /* trim width using multiple alignment method */
static double wg = 11;    /* default gap weight */
static double ws = 1;   /* default space weight */
static BOOLEAN endgaps = TRUE;  /* count end gaps in multiple alignment */
static double distance = 1e-5;  /* squared euclidean distance for convergence */
static char *prior = NULL;  /* prior type input string */
static double beta = -1;  /* scale factor for prior; defaults differ */
static double prob = 1.0; /* try enough subsequences so P=prob */
static int nmotifs = 1;   /* number of motifs to find */
static char *mfile = NULL;  /* name of known .motifs file*/
static int maxiter = 50;        /* max number iterations of EM on best start */
static double nsites = 0; /* try one value of nsites0 only if > 0 */
static int min_nsites = 0;  /* minimum nsites0 to try */
static int max_nsites = 0;  /* maximum nsites0 to try */
static double wnsites = 0.8;  /* weight on prior on nsites */
static int w = 0;   /* width of motifs */
static int min_w = MIN_W; /* minimum W0 to try */
static BOOLEAN all_widths = FALSE;   /* all widths between min and max */
static BOOLEAN min_w_set; /* if set on command line don't override */
static int max_w = MAX_W; /* maximum W0 to try */
static int min_ic = 0;    /* minimum per column IC */
static MAP_TYPE map_type;   /* type of sequence to theta mapping */
static char *mapname = NULL;  /* map type input string */
static double map_scale=-1;   /* scale of sequence to theta mapping:
          Uni - size of add-n prior (n)
          Pam - PAM distance (120)
           Default set in init_em.
        */
static int n_spcons = 0;  /* number of specified start points */
static char *spcons[MAXG];  /* starting point consensus strings */
static int main_hs = HSIZE;     /* size of heap at "main" w values */
static double hs_decrease = HS_DECREASE;   /* Rate of decrease for heap size */
static BOOLEAN x_branch = FALSE; /* Use x_branch regardless of seq model */
static BOOLEAN no_x_branch = FALSE; /* Don't use x_branch, regardless of seq
                                   model */
static BOOLEAN w_branch = FALSE;/* Controls whether width branching occurs */
static BOOLEAN print_heaps = FALSE; /* Print heaps after branching rounds*/
static BOOLEAN print_pred = FALSE; /* Print out the predicted sites after each
                                   round of MEME search (eg subsequence, EM) */
static BOOLEAN print_pllr = FALSE;
                                /* print the LLR of the aligned planted sites*/
static int bfactor = BFACTOR;   /* branching factor for branching search */
//static double param_V = 1;      /* Parameter used by the "deme" objective fn */
//static double pseu = 1;         /* User-specified pseudo-count used in GLAM */
static int maxsize= 100000;   /* dataset size limit */
static int seed = 0;    /* random number seed */
static double seqfrac = 1;  /* fraction of input dataset sequences to use */
static char *meme_directory = NULL; /* meme source directory */
static double max_time = 0; /* maximum allowed CPU time; ignore if 0 */

/*
  subroutines
*/
static void init_meme_background (
  char *bfile,          /* background model file */
  BOOLEAN rc,         /* average reverse comps */
  DATASET *dataset        /* the dataset */
);
static PRIORS *create_priors(
  PTYPE ptype,        /* type of prior to use */
  double beta,        /* beta for dirichlet priors;
             < 0 only returns alphabet */
  DATASET *dataset,     /* the dataset */
  char *plib_name     /* name of prior library */
);

/**********************************************************************/
/*
        init_meme

        Set up all the stuff for meme.
*/
/**********************************************************************/
extern void init_meme(
  int argc,                /* number of input arguments */
  char **argv,             /* input arguments */
  MODEL **model_p,         /* the model */
  MODEL **best_model_p,    /* the best model */
  MODEL **scratch_model_p, /* the best model */
  MODEL **neg_model_p,     /* model of negative examples */
  DATASET **dataset_p,     /* the dataset */
  DATASET **neg_dataset_p, /* dataset of negative examples */
  char *text_filename,     /* name of the text output file */
  char **output_dirname,   /* name of the output directory */
  FILE **text_output       /* destination for text output */
)
{
  int i, j, cc=0, len, pos;
  OBJTYPE objfun = Ev;    /* type of objective function */
  MOTYPE mtype;      /* type of model */
  PTYPE ptype;      /* type of prior */
  PRIORS *priors;       /* the prior probabilities model */
  P_POINT *p_point;   /* previously learned starting points */
  char *alphabet;   /* alphabet for dataset */
  int alength;      /* length of alphabet */
  MODEL *model=NULL, *best_model=NULL, *scratch_model=NULL, *neg_model=NULL;
  DATASET *dataset=NULL, *neg_dataset=NULL;
  double evt = BIG;   /* no E-value cutoff */
  BOOLEAN no_print = FALSE;     /* turn off printing if parallel and not main */

#ifdef PARALLEL
  /* turn off printing if parallel and not the main processor */
  no_print = (mpMyID() != 0);
#endif

  *output_dirname = default_output_dirname;

  /* get the command line arguments */
  i = 1;
#ifndef lint
  /* print the command line */
  argv[0] = "";
  DO_STANDARD_COMMAND_LINE(1,
    USAGE(USAGE:\n\tmeme\t<dataset> [optional arguments]\n);
    NON_SWITCH(1, <dataset> \t\tfile containing sequences in FASTA format\n,
      switch (i++) {
        case 1: datafile = _OPTION_; break;
        default: COMMAND_LINE_ERROR;
      });
     FLAG_OPTN(1, h, \t\t\tprint this message, USAGE_MESSAGE);
     DATA_OPTN(0, neg, <negdataset>,
       \tfile containing negative example sequences, negfile = _OPTION_);
     DATA_OPTN(0, ntype, pair|blend, \thow to use negative examples,
       ntype = _OPTION_);
     DATA_OPTN(1, o, <output dir>,
      \tname of directory for output files\n\t\t\t\twill not replace existing directory,
       *output_dirname = _OPTION_);
     DATA_OPTN(1, oc, <output dir>,
       \tname of directory for output files\n\t\t\t\twill replace existing directory,
       clobber=TRUE; *output_dirname = _OPTION_);
     FLAG_OPTN(1, text, \t\t\toutput in text format (default is HTML),
       TEXT_ONLY = TRUE);
     DATA_OPTN(EXP, objfun, ev|pv, \tobjective function (default: ev),
       obj = _OPTION_);
     FLAG_OPTN(1, dna, \t\t\tsequences use DNA alphabet, alph = "DNA");
     FLAG_OPTN(1, protein, \t\tsequences use protein alphabet,
       alph = "PROTEIN");
     DATA_OPTN(1, mod, oops|zoops|anr, \tdistribution of motifs,mod = _OPTION_);
     DATA_OPTN(1, nmotifs, <nmotifs>, \tmaximum number of motifs to find,
       nmotifs = atoi(_OPTION_));
     DATA_OPTN(1, evt, <ev>, \t\tstop if motif E-value greater than <evt>,
       evt = atof(_OPTION_));
     DATA_OPTN(1, nsites, <sites>, \tnumber of sites for each motif,
       nsites=atof(_OPTION_));
     DATA_OPTN(1, minsites, <minsites>,
       \tminimum number of sites for each motif, min_nsites=atoi(_OPTION_));
     DATA_OPTN(1,
       maxsites, <maxsites>, \tmaximum number of sites for each motif,
       max_nsites=atoi(_OPTION_));
     DATA_OPTN(1, wnsites, <wnsites>, \tweight on expected number of sites,
       wnsites=atof(_OPTION_));
     DATA_OPTN(1, w, <w>, \t\tmotif width, w = atoi(_OPTION_));
     DATA_OPTN(1, minw, <minw>, \t\tminimum motif width,
	       min_w = atoi(_OPTION_); min_w_set=TRUE);
     DATA_OPTN(1, maxw, <maxw>, \t\tmaximum motif width,
       max_w = atoi(_OPTION_));
     DATA_OPTN(EXP, minic, <minic>, \tminimum column information content
       (default: 0), min_ic = atof(_OPTION_));
     FLAG_OPTN(1, nomatrim,
       \t\tdo not adjust motif width using multiple\n\t\t\t\talignment,
       ma_trim = FALSE);
     DATA_OPTN(1, wg, <wg>, \t\tgap opening cost for multiple alignments,
       wg=atof(_OPTION_));
     DATA_OPTN(1, ws, <ws>, \t\tgap extension cost for multiple alignments,
       ws=atof(_OPTION_));
     FLAG_OPTN(1, noendgaps, \t\tdo not count end gaps in multiple alignments,
       endgaps = FALSE);
     DATA_OPTN(1, bfile, <bfile>, \tname of background Markov model file,
       bfile = _OPTION_);
     FLAG_OPTN(1, revcomp, \t\tallow sites on + or - DNA strands,
       revcomp = TRUE);
     FLAG_OPTN(1, pal, \t\t\tforce palindromes (requires -dna), pal = 1);
     DATA_OPTN(1, maxiter, <maxiter>, \tmaximum EM iterations to run,
  maxiter = atoi(_OPTION_));
     DATA_OPTN(1, distance, <distance>, \tEM convergence criterion,
       distance = atof(_OPTION_));
     DATA_OPTN(1, psp, <pspfile>, \tname of positional priors file,
	pspfile = _OPTION_);
     DATA_OPTN(EXP, psp2, <pspfile>, \tname of 2-stranded positional priors file\n\t\t\t\t(requires -revcomp),
	pspfile = _OPTION_;psp2 = TRUE);
     DATA_OPTN(1,
       prior, dirichlet|dmix|mega|megap|addone, \n\t\t\t\ttype of prior to use,
       prior = _OPTION_);
     DATA_OPTN(1, b, <b>, \t\tstrength of the prior, beta = atof(_OPTION_));
     DATA_OPTN(1, plib, <plib>, \t\tname of Dirichlet prior file,
       plib_name = _OPTION_);
     DATA_OPTN(1, spfuzz, <spfuzz>, \tfuzziness of sequence to theta mapping,
  map_scale = atof(_OPTION_));
     DATA_OPTN(1, spmap, uni|pam, \tstarting point seq to theta mapping type,
       mapname = _OPTION_);
     DATA_OPTN(1, cons, <cons>, \t\tconsensus sequence to start EM from,
       spcons[n_spcons++] = _OPTION_);
     DATA_OPTN(1, heapsize, <hs>,
               \tsize of heaps for widths where substring \n\t\t\t\tsearch occurs, 
               main_hs = atoi(_OPTION_));
     FLAG_OPTN(1, x_branch, \t\tperform x-branching, x_branch=TRUE);
     FLAG_OPTN(EXP, no_x_branch, \t\tdo not perform x-branching, 
               no_x_branch=TRUE);
     FLAG_OPTN(1, w_branch, \t\tperform width branching, w_branch=TRUE);
     FLAG_OPTN(1, allw, \t\t\tinclude all motif widths from min to max,
	       all_widths=TRUE);
     DATA_OPTN(1, bfactor, <bf>,
       \t\tbranching factor for branching search, bfactor = atoi(_OPTION_));
     FLAG_OPTN(EXP, print_pred, \t\tprint out the sites predicted by meme,
       print_pred = TRUE);
     FLAG_OPTN(EXP, print_heaps, \t\tprint heaps after each branching round,
       print_heaps = TRUE);
     FLAG_OPTN(EXP, planted_LLR, \t\tprint the LLR of the aligned planted sites,
       print_pllr = TRUE);
     DATA_OPTN(1, maxsize, <maxsize>, \tmaximum dataset size in characters,
       maxsize = atoi(_OPTION_));
     FLAG_OPTN(1, nostatus, \t\tdo not print progress reports to terminal,
       NO_STATUS = TRUE);
     DATA_OPTN(1, p, <np>, \t\tuse parallel version with <np> processors, i=i);
     DATA_OPTN(1, time, <t>, \t\tquit before <t> CPU seconds consumed,
       max_time = atof(_OPTION_));
     DATA_OPTN(1, sf, <sf>, \t\tprint <sf> as name of sequence file, sf = _OPTION_);
     FLAG_OPTN(2, check_syntax, \t\tcheck input syntax and exit,
       check_syntax = TRUE);
     FLAG_OPTN(1, V, \t\t\tverbose mode, VERBOSE = TRUE);
     DATA_OPTN(EXP, mfile, <mfile>, \tfile of known motifs, mfile = _OPTION_);
     DATA_OPTN(EXP, seed, <seed>, \t\tseed for random numbers in sampling,
       seed = atoi(_OPTION_));
     DATA_OPTN(EXP, seqfrac, <seqfrac>, \tfraction of sequences to use,
  seqfrac= atof(_OPTION_));
     FLAG_OPTN(EXP, trace, \t\ttrace starting points, TRACE = TRUE);
     FLAG_OPTN(EXP, print_all, \t\tprint all debug information,
       PRINTALL = TRUE);
     FLAG_OPTN(EXP, print_w, \t\tprint erasure matrix, PRINT_W = TRUE);
     FLAG_OPTN(EXP, print_z, \t\tprint missing information matrix,
       PRINT_Z = TRUE);
     FLAG_OPTN(EXP, print_ll, \t\tprint log-likelihood during EM,
       PRINT_LL = TRUE);
     FLAG_OPTN(EXP, print_starts, \t\tprint starting points, 
       PRINT_STARTS = TRUE);
  )
#endif

  /* exit if check_syntax is on */
  if (check_syntax) exit(0);

  /* set random number generator */
  srand48(seed);

  if (TEXT_ONLY == TRUE) {
    // Legacy: plain text output to standard out.
    *text_output = stdout;
  }
  else {
    if (!no_print) {
      // allow clobbering of the default output directory
      if (*output_dirname == default_output_dirname) { 
	clobber = TRUE;
      } 
      if (create_output_directory(*output_dirname, clobber, !NO_STATUS)) {
	// Failed to create output directory.
	exit(1);
      }
      // Create the name of the text output file 
      // "<dir>/text_filename/" and open it for writing
      char *path = make_path_to_file(*output_dirname, text_filename);
      *text_output = fopen(path, "w"); //FIXME CEG check for errors
      myfree(path);
    }
  }

  /* set all the print flags appropriately */
  if (PRINTALL) {
    PRINT_W = TRUE;
    PRINT_Z = TRUE;
    PRINT_LL = TRUE;
    PRINT_STARTS = TRUE;
  }

  /* check input arguments */
  if (prob < 0 || prob > 1) {
    fprintf(stderr, "-prob <p>, <p> must be between 0 and 1\n");
    exit(1);
  }

  if (nmotifs >= MAXG) {
    fprintf(stderr, "-nmotifs larger than MAXG-1.  Use smaller -nmotifs or recompile with larger MAXG.\n");
   exit(1);
  }

  /* get the name of the directory where MEME is installed */
  if (!meme_directory) {      /* not given on command line */
    meme_directory = MEME_DIR;
  }

  /* get the objective function type */
  if (!strcmp(obj, "pv")) {
    objfun = Pv;
  } else if (!strcmp(obj, "ev")) {
    objfun = Ev;
  } else {
    printf("Unknown objective function type %s. \n", obj);
    exit(1);
  }

  /* get the model type */
  if (!strcmp(mod, "anr") || !strcmp(mod, "tcm")) {
    mtype = Tcm;
    if (pspfile) { /* PM FIXME */
      fprintf (stderr, "-mod anr (-mod tcm) not yet supported for with -psp.\n");
      exit(1);
    }
  } else if (!strcmp(mod, "oops")) {
    mtype = Oops;
  } else if (!strcmp(mod, "zoops")) {
    mtype = Zoops;
  } else {
    mtype = Zoops;        /* prevent warning */
    fprintf(stderr, "Unknown model type %s. \n", mod);
    exit(1);
  }

  /* check seqfrac */
  if (seqfrac > 1 || seqfrac <=0) {
    fprintf(stderr, "seqfrac must be in (0, 1]\n");
    exit(1);
  }

  /* check the alphabet and set up default mappings and priors */
  if (!strcmp(alph, "DNA")) {     /* DNA */
    alphabet = DNA0;
    if (!mapname) mapname = "uni";    /* uniform prior mapping */
    if (!prior) prior = "dirichlet";    /* simple dirichlet prior */
  } else {          /* PROTEIN */
    alphabet = PROTEIN0;
    if (!mapname) mapname = "pam";    /* PAM mapping */
    if (!prior) {
      switch (mtype) {
	case Oops:
	  prior = "dmix";
		break;
	case Zoops:
	case Tcm:
	  prior = "megap";
		break;
	default:
	  prior = "dirichlet"; break;
      }
    }
    if (revcomp) {
      fprintf(stderr, "You cannot use -revcomp with -protein.\n");
      exit(1);
    }
  } /* set alphabet and priors defaults */

  /* find out type of prior */
  if (!strcmp(prior, "dirichlet")) {
    ptype = Dirichlet;
    if (beta < 0) beta = 0.01;      /* default b = 0.01 */
  } else if (!strcmp(prior, "dmix")) {
    ptype = Dmix;
    if (beta < 0) beta = 0;     /* default b = 0 for dmix */
  } else if (!strcmp(prior, "megadmix") || !strcmp(prior, "mega")) {
    ptype = Mega;       /* use mega prior heuristic */
  } else if (!strcmp(prior, "megap")) {
    ptype = MegaP;        /* discretization uses b=0 */
  } else if (!strcmp(prior, "addone")) {
    ptype = Addone;
  } else {
    ptype = Dirichlet;        /* prevent warning */
    fprintf(stderr, "Unknown type of prior: %s!\n", prior);
    exit(1);
  }

  /* convert the alphabet to upper case */
  for (i=0; alphabet && (cc = alphabet[i]); i++)
    if (islower(cc)) alphabet[i] = toupper(cc);

  /* setup hashing function for encoding strings as integers */
  alength = setup_hash_alph(alphabet);

  /* read the samples and set up globals */
  dataset = read_seq_file(datafile, alphabet, revcomp, seqfrac);
  if (!dataset) exit(1);

  // read in psp file if one given
  if (pspfile) {
    read_psp_file(pspfile, dataset, psp2, revcomp, mtype);

    // warn that we are using the W in the PSP file as minw 
    if (!min_w_set) {
      fprintf(stderr,"Setting minimum motif width to width of the prior in the PSP file: %d\n",
	dataset->psp_w);
      min_w = dataset->psp_w;
    }
  } else { 			/* no PSP file */
    dataset->psp_w = min_w;
  }

  /* set dataset alphabet flag */
  dataset->dna = !strcmp(alph, "DNA");

  /* set the dataset objfun */
  dataset->objfun = objfun;

  /* initialize the background model */
  init_meme_background(bfile, revcomp, dataset);

  /* reset the current alphabet since it may be clobbered by the above */
  (void) setup_hash_alph(alphabet);

  /* prevent too long jobs */
  if (dataset->total_res > maxsize) {
    fprintf(stderr, "Dataset too large (> %d).  Rerun with larger -maxsize.\n",
      maxsize);
    exit(1);
  }

  /* Calculate optimal split points for this dataset; no-op if non-parallel */
  balance_loop(dataset->samples, dataset->n_samples);

  /* read in known motifs file */
  if (mfile) {
    FILE  *fdata = fopen(datafile, "r"); //FIXME CEG check for errors
    dataset->nkmotifs = read_motifs(fdata, mfile, dataset->motifs, FALSE, NULL);
    nmotifs = dataset->nkmotifs;
    prob = 0;       /* don't sample starts */
    dataset->pal = 0;     /* no palindrome testing */
  } else {
    dataset->nkmotifs = 0;
  }

  /* create the priors */
  if (ptype == Dmix || ptype == Mega || ptype == MegaP) {
    /* make the name of the prior library */
    if (!plib_name) {
      char *tmp1, *tmp2;
      if (!strcmp(alphabet, PROTEIN0)) {
        plib_name = PROTEIN_PLIB; /* default mixture prior for proteins */
      } else {
        plib_name = DNA_PLIB;   /* default mixture prior for DNA */
      }
      /* prepend meme_directory to file name */
      Strcat(tmp1, meme_directory, "/etc/");
      Strcat(tmp2, tmp1, plib_name);
      plib_name = tmp2;
    }
  }
  if ((ptype == Mega || ptype == MegaP) && beta == -1) {
    /* tlb 5-9-97; wgt_total_res */
    /*beta = 10.0 * dataset->wgt_total_res;*/ /* size of mega prior */
    beta = 5.0 * dataset->wgt_total_res;  /* size of mega prior */
  }
  priors = create_priors(ptype, beta, dataset, plib_name);

  /* check that alphabet in mixture prior matches alphabet */
  if (ptype == Dmix || ptype == Mega || ptype == MegaP) {
    if (strcmp(alphabet,  priors->plib->alphabet)) {
      fprintf(stderr, "The alphabet in the mixture prior file (%s) doesn't match the sequence alphabet.\n", plib_name);
      fprintf(stderr, "prior alphabet   : %s\n", priors->plib->alphabet);
      fprintf(stderr, "sequence alphabet: %s\n", alphabet);
      exit(1);
    }
  }

  /* set number of occurrences of sites */
  if (nsites != 0) {
    if (mtype == Oops) {
      fprintf(stderr, "You may not specify -sites with -mod oops\n");
      exit(1);
    }
    min_nsites = max_nsites = nsites;
  }

  /* set search range for nsites */
  if (mtype == Oops) {
    min_nsites = max_nsites = dataset->n_samples;
  } else if (mtype == Zoops) {
    if (min_nsites > dataset->n_samples) {
      fprintf(stderr, "Minimum number of sites too large.  Setting to 2.\n");
      min_nsites = 2;
    }
    if (!min_nsites) min_nsites = 2;      /* default */
    if (max_nsites > dataset->n_samples) {
      fprintf(stderr, "Maximum number of sites exceeded.  Setting to %d.\n",
        dataset->n_samples);
      max_nsites = dataset->n_samples;
    }
    if (!max_nsites) max_nsites = dataset->n_samples; /* default */
  } else {            /* TCM model */
    if (min_nsites<2) min_nsites = 2;     /* default */
    if (!max_nsites) max_nsites = MIN(5*dataset->n_samples, 50);
  }

  /* check that max number of sites >= min number of sites */
  if (min_nsites > max_nsites) {
    fprintf(
      stderr, 
      "The minimum number of sites is set to %d. "
      "It should be less than the max number of sites (%d).\n",
      min_nsites,
      max_nsites
    );
    exit(1);
  }
  /* check that there are enough possible sites */
  if (min_nsites < 2) {
    fprintf(stderr, "You must specify a minimum of 2 or more sites.\n");
    exit(1);
  }
  if (max_nsites < 2) {
    fprintf(stderr, "It must be possible for at least 2 sites to fit.\n");
    exit(1);
  }

  /* check weight on prior on nsites */
  if (wnsites >= 1 || wnsites < 0) {
    fprintf(stderr, "<wnsites> must be in range [0..1).\n"); exit(1);
  }

  /* check that no sequence too short */
  if (dataset->min_slength < MIN_W) {
    fprintf(stderr,
     "All sequences must be at least %d characters long.  Remove ",
        MIN_W);
    fprintf(stderr, "shorter sequences\nand rerun.\n");
    exit(1);
  }

  /* set up globals */
  if (w != 0) {       /* w specified; set min_w and max_w */
    max_w = min_w = w;
    fprintf(stderr,"w set, setting max and min to %d#######\n",w);
  }

  /* oops model: limit max_w to shortest seq */
  if (mtype == Oops && max_w > dataset->min_slength) {
    max_w = dataset->min_slength;
    fprintf(stderr,
      "maxw > length of shortest sequence (%ld).", dataset->min_slength);
    fprintf(stderr, "  Setting maxw to %d.\n", max_w);
  }
  /* all models: limit max_w to longest seq */
  if (max_w > dataset->max_slength) {
    max_w = dataset->max_slength;
    fprintf(stderr,
      "maxw > length of longest sequence (%ld).", dataset->max_slength);
    fprintf(stderr, "  Setting maxw to %d.\n", max_w);
  }
  if (max_w > MAXSITE) {
    fprintf(stderr,
      "maxw too large (> %d).  Recompile with larger MAXSITE.\n", MAXSITE);
    exit(1);
  }
  if (max_w < 0) {          /* use default */
    max_w = MIN(MAXSITE, dataset->min_slength);   /* maximum W0 */
  }

  /* check that min_w <= max_w */
  if (min_w > max_w) {
    if (pspfile) {
      fprintf(stderr, "PSP file w = %d > maxw = %d. Respecify larger -maxw.\n",
	      dataset->psp_w, max_w);
      exit(1);
    }
     fprintf(stderr, "minw > maxw.  Setting minw to %d.\n", max_w);
     min_w = max_w;
  }

  /* check that min_w and max_w are at least ABS_MIN_W */
  if (min_w < ABS_MIN_W) {
    fprintf(stderr,
      "Minimum width must be >= %d.  Respecify larger -w or -minw.\n",
      ABS_MIN_W);
    exit(1);
  } else if (max_w < ABS_MIN_W) {
    fprintf(stderr,
      "Maximum width must be >= %d.  Respecify larger -w or -maxw.\n",
      ABS_MIN_W);
    exit(1);
  }

  /* must use TCM if only one sequence */
  if (mtype != Tcm && dataset->n_samples==1) {
    fprintf(stderr,
"You must specify '-mod anr' since your dataset contains only one sequence.\n"
    );
    fprintf(stderr,
"Alternatively, you might wish to break your sequence into several sequences.\n"
    );
    exit(1);
  }

  /* check that using right alphabet if model type requires it */
  if (strcmp(alph, "DNA")) {
    if (revcomp) {
      fprintf(stderr,
        "You must use default DNA alphabet if using complementary strand!\n");
      exit(1);
    }
    if (pal) {
      fprintf(stderr, "You must use default DNA alphabet if using -pal !\n");
      exit(1);
    }
  }

  /* flag search for palindromes */
  dataset->pal = pal;

  /* get the type of mapping between sequences and thetas */
  if (!strcmp(mapname, "uni")) {
    map_type = Uni;
    if (map_scale == -1) map_scale = .5;    /* default add .5 */
  } else if (!strcmp(mapname, "pam")) {
    map_type = Pam;
    if (map_scale == -1) map_scale = 120;   /* default PAM 120 */
  } else {
    fprintf(stderr, "Unknown mapping type %s. \n", mapname);
    exit(1);
  }

  /* check that IUPAC alphabet if using PAM mapping */
  // mapname == "pam" && alphabet != PROTEIN0
  if (!strcmp(mapname,"pam") && (strcmp(alphabet,PROTEIN0))) {
    fprintf(stderr,
     "Setting sequence to theta mapping type to `Uni' since alphabet not IUPAC.\n");
    mapname = "uni";
  }

  /* set up the sequence to theta mapping matrix for starts */
  dataset->map = init_map(map_type, map_scale, alength, dataset->back, FALSE);
  dataset->lomap = init_map(map_type, map_scale, alength, dataset->back, TRUE);

  /* set up p_point: previously learned components start points */
  p_point = (P_POINT *) mymalloc(sizeof(P_POINT));
  p_point->c = n_spcons;    /* number of specified starts */
  /* the default starting points */
  for (i=0; i<MAXG; i++) {
    p_point->e_cons0[i] = NULL;
    p_point->w[i] = max_w;
    p_point->nsites[i] = nsites;
  }

  /* setup user-specified start points */
  for (i=0; i<n_spcons; i++) {
    char *e_cons = (char *) mymalloc(MAXSITE);    /* encoded version */
    if (i >= MAXG) {
      fprintf(stderr, "Too many starting points.  Increase MAXG and recompile.");
      exit(1);
    }
    p_point->e_cons0[i] = e_cons;
    /* convert the consensus sequence to upper case and encode as integer */
    for (j=0; (cc = spcons[i][j]) != 0; j++) {
      char uc = (islower(cc) ? toupper(cc) : cc);
      e_cons[j] = (uc == 'X') ? alength : hash(uc);
      if (e_cons[j] > alength) {
  fprintf(stderr, "Illegal letter %c in consensus string!\n", cc);
  exit(1);
      }
    }
    /* set width to length of consensus unless -w specified */
    if (w == 0) p_point->w[i] = j;
    /* pad out the consensus sequence with A's */
    for ( ; j < MAXSITE; j++) e_cons[j] = 0;
  }

  // Setup heap size for storage of starting points:
  if (main_hs < 1) {
    fprintf(stderr, "Heap size must be >= 1.\n");
    exit(1);
  } else {
    dataset->main_hs = main_hs;
    dataset->hs_decrease = hs_decrease;
  }

  // Setup a struct recording the desired parameters for branching search:
  BRANCH_PARAMS *branch_params = NULL;
  branch_params = mymalloc(sizeof(BRANCH_PARAMS));

  // Setup branching factor for branching search:
  if (bfactor < 1) {
    fprintf(stderr, "bfactor must be >= 1.\n");
    exit(1);
  } else {
    branch_params->bfactor = bfactor;
  }

  // Record whether the user wants w-branching carried out:
  branch_params->w_branch = w_branch;

  // Record whether the user wants x-branching carried out:
  if (x_branch){
    branch_params->point_branch = X_ONLY;
  } else { 
    branch_params->point_branch = NO_POINT_B;
  }

  if (FALSE) {
  // FIXME: This code can be used to set branching as the default 
  // for different alphabets/sequence models using the -x_branch
  // and -no_x_branch switches. Currently the default is set to
  // no branching (in previous if statement).
  /* Record whether the user wants x_branching carried out. The result
     is recorded in the "point_branch" attribute of branch_params. Note
     that the user is only able to specify x_branching or no x_branching.
     This is because we found ACGT (ie regular) branching to be of little
     benefit. */
  if (x_branch) {
    // User has specified they want x-branching... 
    if (no_x_branch) {
      // Invalid to specify x_branch and no_branch simultaneously:
      fprintf(stderr, "x_branch and no_branch cannot be specified"\
                      "at the same time");
      exit(1);
    }
    else {
      branch_params->point_branch = X_ONLY;
    }
  } else if (no_x_branch) {
    branch_params->point_branch = NO_POINT_B;
  } else {
    // User did not specify x_branching => Decide based on sequence model:
    // NOTE: branching is the default for oops for DNA only
    if (mtype == Oops && !strcmp(alph, "DNA")) {
      branch_params->point_branch = X_ONLY; // Only x_branch under oops by default.
    } else {
      branch_params->point_branch = NO_POINT_B; // Only x_branch under oops by default.
    }
  } // Deciding branch_params->x_branch
  } // end if 

  // Store the branching parameters for future reference:
  dataset->branch_params = branch_params;

  // Print sites predicted by MEME:
  dataset->print_pred = print_pred;

  // Record whether heaps are to be printed:
  dataset->print_heaps = print_heaps;

  // Record whether llr is to be printed for the alignment of planted sites:
  if (print_pllr) {
    if (w <= 0) {
      fprintf(stderr,
              "Not valid to request llr of planted sites unless length of sites"
              " has been specified.\n");
      exit(1);
    }
  }
  dataset->print_pllr = print_pllr;

  /* make sure nmotifs is as large as the number of starting points */
  if (nmotifs < n_spcons) {
    nmotifs = n_spcons;
    fprintf(stderr, "Setting nmotifs to %d\n", nmotifs);
  }

  /* create the model */
  model = create_model(mtype, revcomp, max_w, alength);
  best_model = create_model(mtype, revcomp, max_w, alength);
  scratch_model = create_model(mtype, revcomp, max_w, alength);

  /* create the negative dataset and negative model if negative file given */
  if (negfile) {
    neg_dataset = read_seq_file(negfile, dataset->alphabet, revcomp, seqfrac);
    if (!strcmp(ntype, "pair")) {
      neg_dataset->negtype = Pair;
    } else if (!strcmp(ntype, "blend")) {
      neg_dataset->negtype = Blend;
    } else {
      fprintf(stderr, "Unknown ntype %s. \n", ntype);
      exit(1);
    }
    neg_dataset->pal = dataset->pal;
    neg_dataset->back = dataset->back;
    neg_model = create_model(mtype, revcomp, max_w, alength);
  }

  /* initialize log and exp lookup tables */
  init_log();
  init_exp();

  /* Initialize the probability tables for the objective function.  */
  /* Get for 2 and up because weighting may cause < min_nsites sites */
  init_llr_pv_tables(2, max_nsites, alength, dataset->back, dataset->pal);

  /* set up scratch models */
  model->min_w = min_w;
  model->max_w = max_w;
  model->all_widths = all_widths;
  model->min_nsites = min_nsites;
  model->max_nsites = max_nsites;
  copy_model(model, best_model, alength);
  copy_model(model, scratch_model, alength);

  /* put meme parameters in dataset */
  dataset->priors = priors;
  dataset->p_point = p_point;
  dataset->wg = wg;
  dataset->ws = ws;
  dataset->endgaps = endgaps;
  dataset->wnsites = wnsites;
  dataset->ma_adj = ma_trim;
  dataset->distance = distance;
  dataset->prob = prob;
  dataset->nmotifs = nmotifs;
  dataset->maxiter = maxiter;
  dataset->evt = evt;
  dataset->mod = mod;
  dataset->mapname = mapname;
  dataset->map_scale = map_scale;
  dataset->priorname = prior;
  dataset->beta = beta;
  dataset->seed = seed;
  dataset->seqfrac = seqfrac;
  /* save name of prior library */
  if (plib_name) {
    dataset->plib_name = plib_name + strlen(plib_name);
    while (*dataset->plib_name != '/') dataset->plib_name--;
    dataset->plib_name++; /* strip off directory */
  } else {
    dataset->plib_name = NULL;
  }
  dataset->datafile = sf ? sf : datafile; /* name to print */
  dataset->negfile = negfile;
  dataset->bfile = bfile;
  dataset->min_ic = min_ic;
  dataset->max_time = max_time;

  /* save command line */
  dataset->command = NULL;
  argv[0] = "meme";
  /*for (i=pos=len=0; i<argc-2; i++) {*/  /* don't save last arguments */
  /* save all arguments; used to not save last 2; why??? */
  for (i=pos=len=0; i<argc; i++) {
    len += strlen(argv[i])+1;     /* +1 for space following */
    Resize(dataset->command, len+2, char);  /* +1 for null */
    strcpy((dataset->command)+pos, argv[i]);
    dataset->command[len-1] = ' ';
    dataset->command[len] = '\0';
    pos = len;
  }
  if (TEXT_ONLY) {
    dataset->output_directory = NULL;
  }
  else {
    dataset->output_directory = *output_dirname;
  }

  /* set up return values */
  *model_p = model;
  *scratch_model_p = scratch_model;
  *best_model_p = best_model;
  *neg_model_p = neg_model;
  *dataset_p = dataset;
  *neg_dataset_p = neg_dataset;

  /* announce meme */
  banner("MEME", *text_output);
  fprintf(*text_output, "\n\n");

} /* init_meme */

/***************************************************************************/
/*
        init_meme_background

        Read in the background Markov model or use the residue
  	frequencies adjusted by add-one prior.

        Precalculate the log cumulative background probabilities.

        The log probability of any substring of any sequence will then
        be accessible via:
                Log_back(sample[i]->logcumback, j, w)
        where j is the start of the length-w substring of sequence i.
*/
/***************************************************************************/
static void init_meme_background (
  char *bfile,          /* background model file */
  BOOLEAN rc,         /* average reverse comps */
  DATASET *dataset        /* the dataset */
)
{
  int i, j;
  char *alphabet = dataset->alphabet;   /* alphabet */
  SAMPLE **samples = dataset->samples;    /* the sequences */
  int n_samples = dataset->n_samples;   /* # of sequences in dataset */
  double *back;         /* freq. of tuples */
  int order;          /* order of background model */
  BOOLEAN add_x = TRUE;       /* add x-tuples to model */

  /* read in background Markov model or use dataset frequencies (0-order) */
  if (bfile) {          /* read in probs; set X probs */
    back = read_markov_model(bfile, NULL, alphabet, add_x, rc, &order);
  } else {          /* use dataset frequencies */
    int alength = dataset->alength;   /* length of alphabet */
    double *res_freq = dataset->res_freq; /* weighted residue freq */
    double wtr = dataset->wgt_total_res;  /* weighted residue count */
    double *adj_res_freq = NULL;
    Resize(adj_res_freq, alength, double);
    /* adjust residue frequencies with add-one prior to avoid 0s */
    for (i=0; i<alength; i++)
      adj_res_freq[i] = (res_freq[i]*wtr + 1)/(wtr+alength);
    back = read_markov_model(NULL, adj_res_freq, alphabet, add_x, rc, &order);
    order = 0;
    myfree(adj_res_freq);
  }

  /* precalculate the log cumulative background probabilities for each seq */
  dataset->log_total_prob = 0;
  for (i=0; i<n_samples; i++) {     /* sequence */
    SAMPLE *s = samples[i];     /* sequence */
    char *seq = NULL;       /* ascii sequence */
    double *lcb = s->logcumback;    /* log cum. back. prob */
    /* unhash the sequence to convert ambiguous characters to X */
    Resize(seq, s->length+1, char);
    for (j=0; j<s->length; j++) seq[j] = unhash(s->res[j]);
    seq[s->length] = '\0';
    /* compute probabilites */
    dataset->log_total_prob += log_cum_back(seq, back, order, lcb);
    myfree(seq);
  } /* sequence */

  /* return values */
  dataset->back = back;
  dataset->back_order = order;

} /* init_meme_background */

/**********************************************************************/
/*
  create_priors
*/
/**********************************************************************/
static PRIORS *create_priors(
  PTYPE ptype,        /* type of prior to use */
  double beta,        /* beta for dirichlet priors;
             < 0 only returns alphabet */
  DATASET *dataset,     /* the dataset */
  char *plib_name     /* name of prior library */
)
{
  int i;
  int alength = (dataset ? dataset->alength : 0);
  double *back = (dataset ? dataset->back : (double *)NULL);
  PRIORS *priors = (PRIORS *) mymalloc(sizeof(PRIORS));

  priors->ptype = ptype;

  /* set up the prior counts */
  switch (ptype) {
    case Addone:        /* add one prior */
      for (i=0; i<alength; i++) priors->prior_count[i] = 1.0;
      break;
    case Dirichlet:       /* simple dirichlet prior */
      for (i=0; i<alength; i++) priors->prior_count[i] = beta * back[i];
      break;
    case Dmix:          /* mixture of dirichlet's */
    case Mega:          /* megaprior heuristic */
    case MegaP:         /* mod. megaprior heuristic */
    {
      priors->plib = read_PriorLib(plib_name, beta);

      /* get b=0 prior for modified mega prior heuristic */
      if (ptype == MegaP || ptype == Mega) {  /* used adj freq with Mega */
        double b = 0;
  priors->plib0 = read_PriorLib(plib_name, b);
      }

      break;
    }
  }
  return priors;
} /* create_priors */
