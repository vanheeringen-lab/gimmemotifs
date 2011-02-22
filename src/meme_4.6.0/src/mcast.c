#include <errno.h>
#include <libgen.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/wait.h>
#include "../config.h"
#include "dir.h"
#include "io.h"
#include "simple-getopt.h"
#include "utils.h"

typedef enum { MEME_FORMAT, TRANSFAC_FORMAT } MOTIF_FORMAT_T;

VERBOSE_T verbosity = NORMAL_VERBOSE;

typedef struct mcast_options {
  BOOLEAN_T allow_clobber;
  BOOLEAN_T use_synth;
  BOOLEAN_T text_only;
  BOOLEAN_T quiet;
  MOTIF_FORMAT_T motif_format;
  char *bg_filename;
  char *motif_filename;
  char *output_dirname;
  char *seq_filename;
  char *max_gap;
  char *pseudo_weight;
  char *eg_cost;
  char *e_thresh;
  char *p_thresh;
} MCAST_OPTIONS_T;

/***********************************************************************
  Process command line options
 ***********************************************************************/
static void process_command_line(
  int argc,
  char* argv[],
  MCAST_OPTIONS_T *options
) {

  // Set default values for command line arguments
  options->allow_clobber = TRUE;
  options->text_only = FALSE;
  options->motif_format = MEME_FORMAT;
  options->bg_filename = NULL;
  options->motif_filename = NULL;
  options->output_dirname = "mcast_out";
  options->seq_filename = NULL;
  options->max_gap = "50";
  options->pseudo_weight = "4.0";
  options->e_thresh = "10.0";
  options->p_thresh = "0.0005";
  options->quiet = FALSE;
  options->use_synth = FALSE;

  // Define command line options.
  int option_count = 12;
  int option_index = 0;
  cmdoption const cmdoptions[] = {
		{"bgfile", REQUIRED_VALUE},
		{"bgweight", REQUIRED_VALUE},
		{"e-thresh", REQUIRED_VALUE},
		{"max-gap", REQUIRED_VALUE},
		{"o", REQUIRED_VALUE},
		{"oc", REQUIRED_VALUE},
		{"p-thresh", REQUIRED_VALUE},
		{"synth", NO_VALUE},
		{"text", NO_VALUE},
		{"transfac", NO_VALUE},
		{"verbosity", REQUIRED_VALUE},
		{"quiet", NO_VALUE},
  };
  simple_setopt(argc, argv, option_count, cmdoptions);

  // Define the usage message.
  const char *usage = 
    "USAGE: mcast [options] <query> <database>\n"
     "\n"
     "  [--bgfile <file>]       File containing n-order Markov background model\n"
     "  [--bgweight <b>]        Add b * background frequency to each count in query\n" 
     "                           (default: 4.0 )\n"
     "  [--e-thresh <value>]     Print matches with E-values less than E\n"
     "                           (default = 10)\n"
     "  [--max-gap <value>]      Maximum allowed distance between adjacent hits\n"
     "                           (default = 50)\n"
     "  [--o <output dir>]       Name of output directory. Existing files will not be\n"
     "                           overwritten (default=mcast_out)\n"
     "  [--oc <output dir>]      Name of output directory. Existing files will be\n"
     "                           overwritten.\n"
     "  [--p-thresh <value>]     p-value threshold for motif hits\n"
     "                           (default = 0.0005).\n"
     "  [--synth]                Use synthetic scores for distribution\n"
     "  [--text]                 Output plain text rather than HTML.\n"
     "  [--transfac]             Query is in TRANSFAC format (default: MEME format)\n"
     "  [--verbosity <value>]    Verbosity of error messagess\n";

  // Parse the command line.
  while (1) { 
    int c = 0;
    char* option_name = NULL;
    char* option_value = NULL;
    const char* message = NULL;

    // Read the next option, and break if we're done.
    c = simple_getopt(&option_name, &option_value, &option_index);
    if (c == 0) {
      break;
    } else if (c < 0) {
      simple_getopterror(&message);
      die("Error in command line arguments.\n%s", usage);
    }

    // Assign the parsed value to the appropriate variable
    if (strcmp(option_name, "bgfile") == 0) {
      options->bg_filename = option_value;
    } else if (strcmp(option_name, "bgweight") == 0) {
      options->pseudo_weight = option_value;
    } else if (strcmp(option_name, "e-thresh") == 0) {
      options->e_thresh = option_value;
    } else if (strcmp(option_name, "max-gap") == 0) {
      options->max_gap = option_value;
    }
    else if (strcmp(option_name, "o") == 0){
      // Set output directory with no clobber
      options->output_dirname = option_value;
      options->allow_clobber = FALSE;
    }
    else if (strcmp(option_name, "oc") == 0){
      // Set output directory with clobber
      options->output_dirname = option_value;
      options->allow_clobber = TRUE;
    } else if (strcmp(option_name, "p-thresh") == 0) {
      options->p_thresh = option_value;
    } else if (strcmp(option_name, "synth") == 0) {
      options->use_synth = TRUE;
    } else if (strcmp(option_name, "text") == 0) {
      options->text_only = TRUE;
    } else if (strcmp(option_name, "transfac") == 0) {
      options->motif_format = TRANSFAC_FORMAT;
    } else if (strcmp(option_name, "quiet") == 0) {
      options->quiet = TRUE;
    } else if (strcmp(option_name, "verbosity") == 0) {
      verbosity = (VERBOSE_T) atoi(option_value);
    } 
  }

  // Read and verify the two required arguments.
  if (option_index + 2 != argc) {
    fprintf(stderr, "%s", usage);
    exit(1);
  }
  options->motif_filename = argv[option_index];
  options->seq_filename = argv[option_index+1];
     
  if (options->motif_filename == NULL) {
    die("No motif file given.\n%s", usage);
  }
  if (options->seq_filename == NULL) {
    die("No sequence file given.\n%s", usage);
  }
}

/******************************************************************************
 * This function uses the program transfac2meme to convert a TRANSFAC format 
 * motif file into a MEME format motif file. The name of a file containing a 
 * background model is passed in the third parameter. This can be NULL.
 *****************************************************************************/
void run_transfactomeme(
  char *transfac_filename,
  char *meme_filename,
  char *bg_filename
) {

  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "%s", "Converting motif file to MEME format.\n");
  }

  // Build the command line for transfac2meme
  char cmd[128];
  int result = 0;
  if (bg_filename) {
    result = snprintf(
               cmd, 
               127, 
               "transfac2meme -bg %s %s > %s", 
               bg_filename, 
               transfac_filename,
               meme_filename
             );
  }
  else {
    result = snprintf(
               cmd, 
               127, 
               "transfac2meme %s > %s", 
               transfac_filename, 
               meme_filename
             );
  }
  if (result < 0) {
    perror("Error generating transfac2meme command string");
    die("Unable to create transfact2meme command string");
  }

  // Execute transfac2meme
  int status = system(cmd);
  // Check for interrupts
  if (WIFSIGNALED(status) && 
     (WTERMSIG(status) == SIGINT || WTERMSIG(status) == SIGQUIT)) {
       die("Interrupted while running transfac2meme");
  }
  if (status == -1) {
    // Unable to start transfac2meme
    perror("Error starting transfac2meme command");
    die("Unable to start transfact2meme command");
  }

  // Check result of running transfac2meme
  result = WEXITSTATUS(status);
  if (result != 0) {
    die("transfac2meme failed with exit status %d", result);
  }
}

/******************************************************************************
 * This function builds a command string
 * The caller is responsible for deallocating the string afterwards
 *****************************************************************************/
char* build_cmd(char* desc, char * pattern, ...) {
  va_list ap;
  int length;
  char *cmd;
  //because older versions of vsnprintf just return -1 when something doesn't 
  //fit, give them a chance to work by having a non-zero test buffer.
  char test[256];
  //first measure the string
  va_start(ap, pattern);
  length = vsnprintf(test, 256, pattern, ap);
  va_end(ap);//ap structure is mutated on use, so must end and reinit
  if (length < 0) {
    die("Unable to create %s command string, error given as %s", strerror(errno));
  }
  //increment by one to include the null byte
  length += 1;
  //allocate memory for the command
  cmd = (char*)mm_malloc(sizeof(char) * length);
  //va start must be called a second time as the ap structure is mutated on use
  va_start(ap, pattern);
  length = vsnprintf(cmd, length, pattern, ap);
  va_end(ap);
  if (length < 0) {
    int error_num = errno;
    free(cmd);
    die("Unable to create %s command string, error given as %s", strerror(error_num));
  }
  return cmd;
}

/******************************************************************************
 * This function builds a star topology HMM from a MEME format motif file.
 * The HMM is stored in a text file in the output directory.
 *****************************************************************************/
void run_mhmm(
  char *motif_filename,
  char *hmm_filename
) {

  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "%s", "Creating HMM from motif file.\n");
  }

  // Build the command line for mhmm
  const char *bin_dir = get_meme_bin_dir();
  char *cmd;
  cmd = build_cmd(
             "mhmm", 
             "%s/mhmm --type star --keep-unused --verbosity %d %s > %s", 
             bin_dir,
             verbosity,
             motif_filename,
             hmm_filename
           );

  // Execute mhmm
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Executing %s\n", cmd);
  }
  int status = system(cmd);
  // Check for interrupts
  if (WIFSIGNALED(status) && 
     (WTERMSIG(status) == SIGINT || WTERMSIG(status) == SIGQUIT)) {
       die("Interrupted while running mhmm");
  }
  if (status == -1) {
    // Unable to start mhmm
    perror("Error starting mhmm command");
    die("Unable to start mhmm command (%s)", cmd);
  }

  // Check result of running mhmm
  int result = WEXITSTATUS(status);
  if (result != 0) {
    die("mhmm failed with exit status %d, command run was: %s", result, cmd);
  }
  free(cmd);
}

/******************************************************************************
 * This function scores a sequence database using mhmmscan
 *****************************************************************************/
void run_mhmmscan(
  MCAST_OPTIONS_T *options,
  char *hmm_filename,
  char *seq_filename,
  char *scores_filename
) {

  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "%s", "Creating HMM from motif file.\n");
  }

  // Build the command line for mhmmscan
  const char *bin_dir = get_meme_bin_dir();
  char *cmd;
  cmd = build_cmd(
             "mhmmscan", 
             "%s/mhmmscan %s %s %s --fancy --allow-weak-motifs %s --p-thresh %s "
             "--max-gap %s --e-thresh %s --eg-cost 1 --pseudo-weight %s %s "
             "%s %s > %s", 
             bin_dir,
             options->bg_filename != NULL ? "--bg-file" : "",
             options->bg_filename != NULL ? options->bg_filename : "",
             options->text_only ? "--text" : "",
             options->use_synth ? "--synth" : "",
             options->p_thresh,
             options->max_gap,
             options->e_thresh,
             options->pseudo_weight,
             options->quiet ? "--quiet" : "",
             hmm_filename,
             seq_filename,
             scores_filename
           );

  // Execute mhmmscan
  if (verbosity >= NORMAL_VERBOSE) {
    fprintf(stderr, "Executing %s\n", cmd);
  }
  int status = system(cmd);
  // Check for interrupts
  if (WIFSIGNALED(status) && 
     (WTERMSIG(status) == SIGINT || WTERMSIG(status) == SIGQUIT)) {
       die("Interrupted while running mhmmscan");
  }
  if (status == -1) {
    // Unable to start mhmmscan
    perror("Error starting mhmmscan command");
    die("Unable to start mhmmscan command (%s)", cmd);
  }

  // Check result of running mhmmscan
  int result = WEXITSTATUS(status);
  if (result != 0) {
    die("mhmmscan failed with exit status %d, command run was: %s", result, cmd);
  }
  free(cmd);
}
int main(int argc, char* argv[]) {

  MCAST_OPTIONS_T options;
  process_command_line(argc, argv, &options);

  //
  // Create output directory
  //
  if (create_output_directory(
        options.output_dirname,
        options.allow_clobber,
        FALSE /* Don't print warning messages */
      ) != 0) {
    // Failed to create output directory.
    die("Couldn't create output directory %s.\n", options.output_dirname);
  }

  //
  // If needed, convert motif file to MEME format using transfac2meme.
  //
  char *motif_basename = basename(options.motif_filename); // Using GNU basename
  if (options.motif_format == TRANSFAC_FORMAT) {
    // Build the name for the new MEME format file in the output directory.
    char *meme_filename = concat_string(motif_basename, ".meme");
    char *meme_path = make_path_to_file(options.output_dirname, meme_filename);
    myfree(meme_filename);
    run_transfactomeme(
        options.motif_filename,
        meme_path,
        options.bg_filename
      );
    // Replace motif file name with new name.
    options.motif_filename = meme_path;
  }

  //
  // Build the HMM using mhmm.
  //
  char *hmm_basename = concat_string(motif_basename, ".mhmm");
  char *hmm_path = make_path_to_file(options.output_dirname, hmm_basename);
  myfree(hmm_basename);
  run_mhmm(options.motif_filename, hmm_path);

  //
  // Read and score the sequences using mhmmscan.
  //
  char *score_path = make_path_to_file(
    options.output_dirname, 
    options.text_only == TRUE ? "mcast.txt" : "mcast.html"
  );
  run_mhmmscan(&options, hmm_path, options.seq_filename, score_path);

  //
  // Clean up
  //
  if (options.motif_format == TRANSFAC_FORMAT) {
    // If transfac format was used we have to 
    // clean up the string naming the MEME format
    // motif file.
    myfree(options.motif_filename);
  }
  myfree(hmm_path);
  myfree(score_path);

  return 0;

}


