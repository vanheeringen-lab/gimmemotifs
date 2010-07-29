#include <limits.h>
#include <stdlib.h>
#include <unistd.h>  /* non-ANSI */
#include <strings.h>
#include "glam2_util.h"
#include "glam2_args.h"

/* default values: */
static const char* default_output_dirname = "glam2_out";
static const int n_def = 10000;
static const int r_def = 10;
static const int a_def = 2;
static const int b_def = 50;
static const int z_def = 2;
static const int w_def = 20;
static const double D_def = 0.1;
static const double E_def = 2;
static const double I_def = 0.02;
static const double J_def = 1;
static const double m_def = 1;  /* ??? */
static const int x_def = 0;
static const double t_def = 1.2;  /* seems to work well, so does 1.1 */
static const double c_def = 1.44;  /* seems to work well, so does 1.21 */
static const double u_def = 0.1;
static const double q_def = 1e99;
static const unsigned s_def = 1;
static const int Q_def = FALSE;

static void usage(void) {
  die("\
Usage: glam2 [options] alphabet my_seqs.fa\n\
Main alphabets: p = proteins, n = nucleotides\n\
Main options (default settings):\n\
-h: show all options and their default settings\n\
-o: output directory; will not clobber existing files\n\
-O: output directory (%s); allow clobbering\n\
-r: number of alignment runs (%d)\n\
-n: end each run after this many iterations without improvement (%d)\n\
-2: examine both strands - forward and reverse complement\n\
-z: minimum number of sequences in the alignment (%d)\n\
-a: minimum number of aligned columns (%d)\n\
-b: maximum number of aligned columns (%d)\n\
-w: initial number of aligned columns (%d)\n\
-Q: run quietly\n\
", default_output_dirname, r_def, n_def, z_def, a_def, b_def, w_def);
}

static void help(void) {
  die("\
Usage: glam2 [options] alphabet my_seqs.fa\n\
Alphabets: p = proteins, n = nucleotides, other = alphabet file\n\
Options (default settings):\n\
-h: show all options and their default settings\n\
-o: output directory; will not clobber existing files\n\
-O: output directory (%s); allow clobbering\n\
-r: number of alignment runs (%d)\n\
-n: end each run after this many iterations without improvement (%d)\n\
-2: examine both strands - forward and reverse complement\n\
-z: minimum number of sequences in the alignment (%d)\n\
-a: minimum number of aligned columns (%d)\n\
-b: maximum number of aligned columns (%d)\n\
-w: initial number of aligned columns (%d)\n\
-d: Dirichlet mixture file\n\
-D: deletion pseudocount (%g)\n\
-E: no-deletion pseudocount (%.1f)\n\
-I: insertion pseudocount (%g)\n\
-J: no-insertion pseudocount (%.1f)\n\
-q: weight for generic versus sequence-set-specific residue abundances (%g)\n\
-t: initial temperature (%g)\n\
-c: cooling factor per n iterations (%g)\n\
-u: temperature lower bound (%g)\n\
-p: print progress information at each iteration\n\
-m: column-sampling moves per site-sampling move (%.1f)\n\
-x: site sampling algorithm: 0=FAST 1=SLOW 2=FFT (%d)\n\
-s: seed for pseudo-random numbers (%u)\n\
-Q: run quietly\n\
=== Arguments used only by web the GLAM2 web server ===\n\
-M:  embed sequence file contents as hidden field in HTML\n\
-A <address>:  make email address a hidden field in HTML \n\
-X <description>:  make description a hidden field in HTML\n\
",  default_output_dirname,
   r_def, n_def, z_def, a_def, b_def, w_def, D_def, E_def, I_def, J_def, q_def,
   t_def, c_def, u_def, m_def, x_def, s_def);
}

void getargs(args *a, int argc, char **argv) {
  int c;

  a->stop_after = n_def;
  a->runs = r_def;
  a->out_dir = (char *) default_output_dirname;
  a->clobber = 1;		// allow default directory to be clobbered
  a->two_strands = 0;
  a->min_width = a_def;
  a->max_width = b_def;
  a->min_seqs = z_def;
  a->init_width = w_def;
  a->delete_pseudo = D_def;
  a->no_delete_pseudo = E_def;
  a->insert_pseudo = I_def;
  a->no_insert_pseudo = J_def;
  a->profile = 0;
  a->column_sample_rate = m_def;
  a->algorithm = x_def;
  a->temperature = t_def;
  a->cool = c_def;
  a->frozen = u_def;
  a->dirichlet_file = NULL;
  a->bg_pseudo = q_def;
  a->seed = s_def;
  a->quiet = Q_def;
  a->embed = "";		// "-M"
  a->address = "";		// "-A <email address>"
  a->description = "";		// "-X <description line>"

  /* non-ANSI: */
  while ((c = getopt(argc, argv, "hn:r:o:O:2a:b:z:w:D:E:I:J:pm:x:t:c:u:d:q:s:QMA:X:"))
	 != -1) {
    switch (c) {
    case 'h':
      help();
    case 'n':
      a->stop_after = xatoi(optarg);
      if (a->stop_after < 1)
	die("%s: option -n should be at least 1\n", prog_name);
      break;
    case 'r':
      a->runs = xatoi(optarg);
      if (a->runs < 1)
	die("%s: option -r should be at least 1\n", prog_name);
      break;
    case 'o':
      a->out_dir = optarg;
      a->clobber = 0;
      break;
    case 'O':
      a->out_dir = optarg;
      a->clobber = 1;
      break;
    case '2':
      a->two_strands = 1;
      break;
    case 'a':
      a->min_width = xatoi(optarg);
      if (a->min_width < 2)
	die("%s: option -a should be at least 2\n", prog_name);
      break;
    case 'b':
      a->max_width = xatoi(optarg);
      if (a->max_width < 2)
	die("%s: option -b should be at least 2\n", prog_name);
      break;
    case 'z':
      a->min_seqs = xatoi(optarg);
      if (a->min_seqs < 0)
	die("%s: option -z should be at least 0\n", prog_name);
      break;
    case 'w':
      a->init_width = xatoi(optarg);
      if (a->init_width < 2)
        die("%s: option -w should be at least 2\n", prog_name);
      break;
    case 'D':
      a->delete_pseudo = xatof(optarg);
      if (a->delete_pseudo <= 0)
	die("%s: option -D should be > 0\n", prog_name);
      break;
    case 'E':
      a->no_delete_pseudo = xatof(optarg);
      if (a->no_delete_pseudo <= 0)
	die("%s: option -E should be > 0\n", prog_name);
      break;
    case 'I':
      a->insert_pseudo = xatof(optarg);
      if (a->insert_pseudo <= 0)
	die("%s: option -I should be > 0\n", prog_name);
      break;
    case 'J':
      a->no_insert_pseudo = xatof(optarg);
      if (a->no_insert_pseudo <= 0)
	die("%s: option -J should be > 0\n", prog_name);
      break;
    case 'p':
      a->profile = 1;
      break;
    case 'm':
      a->column_sample_rate = xatof(optarg);
      if (a->column_sample_rate < 0)
	die("%s: option -m should be at least 0\n", prog_name);
      break;
    case 'x':
      a->algorithm = xatoi(optarg);
      if (a->algorithm < 0 || a->algorithm > 2)
	die("%s: option -x should be 0, 1, or 2\n", prog_name);
#ifndef FFT
      if (a->algorithm == 2)
	die("%s: recompile with FFT in order to use -x 2\n", prog_name);
#endif
      break;
    case 't':
      a->temperature = xatof(optarg);
      if (a->temperature <= 0)
	die("%s: option -t should be > 0\n", prog_name);
      break;
    case 'c':
      a->cool = xatof(optarg);
      if (a->cool <= 0)
	die("%s: option -c should be > 0\n", prog_name);
      break;
    case 'u':
      a->frozen = xatof(optarg);
      if (a->frozen < 0)
	die("%s: option -u should be at least 0\n", prog_name);
      break;
    case 'd':
      a->dirichlet_file = optarg;
      break;
    case 'q':
      a->bg_pseudo = xatof(optarg);
      if (a->bg_pseudo <= 0)  /* zero could lead to log(0) */
	die("%s: option -q should be > 0\n", prog_name);
      break;
    case 's':
      a->seed = xatou(optarg);
      break;
    case 'Q':
      a->quiet = TRUE;
      break;
    case 'M':
      a->embed = "M";
      break;
    case 'A':
      a->address = optarg;
      break;
    case 'X':
      a->description = optarg;
      // Replace whitespace in description with '&'
      char *dptr = a->description;
      for (; *dptr != '\0'; dptr++) { 
        if (*dptr == ' ' || *dptr == '\t' || *dptr == '\n') *dptr = '&'; 
      }
      break;
    case '?':
      usage();
    }
  }

  if (optind != argc-2)
    usage();
  a->alph_name = argv[optind++];
  a->seq_file = argv[optind++];

  if (strlen(a->embed)) a->embed = a->seq_file;

  if (a->max_width < a->min_width)
    die("%s: option -a should be >= option -b\n", prog_name);
  else if (a->max_width == a->min_width)
    fprintf(stderr, "%s: warning: setting -a equal to -b is not recommended\n",
	    prog_name);

  if (a->init_width < a->min_width)
    a->init_width = a->min_width;
  else if (a->init_width > a->max_width)
    a->init_width = a->max_width;
}

void printargs(FILE *fp, int argc, char **argv) {
  int i;
  for (i = 0; i < argc; ++i) {
    fputs(argv[i], fp);
    putc(i < argc-1 ? ' ' : '\n', fp);
  }
}
