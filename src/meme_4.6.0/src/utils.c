/********************************************************************
 * FILE: utils.c
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 9-8-97
 * PROJECT: shared
 * COPYRIGHT: 1997-2008 WSN
 * DESCRIPTION: Various useful generic utilities.
 ********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/resource.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include <err.h>
#include <errno.h>
#include <ctype.h>
#include "dir.h"
#include "utils.h"


/***********************************************************************
 * Return the value to replace a missing value -- NaN.
 ***********************************************************************/
double NaN
  (void)
{
  return atof("NaN");
}

/**********************************************************************
 * See .h file for description.
 **********************************************************************/
#ifdef NOCLOCK
double myclock() {return(0);}
#else

#ifdef crayc90
/* No problem on the CRAY. */
#include <time.h>
double myclock() {return((double)clock());}

#else
int getrusage(int who, struct rusage *rusage);

double myclock()
{
  static BOOLEAN_T first_time = TRUE;
  static double    start_time;
  double           elapsed;
  struct rusage    ru;

  if (first_time) {
    getrusage(RUSAGE_SELF, &ru);
    start_time = (ru.ru_utime.tv_sec * 1.0E6) + ru.ru_utime.tv_usec;
    first_time = FALSE;
    return 0;

  } else {
    getrusage(RUSAGE_SELF, &ru);
    elapsed = (ru.ru_utime.tv_sec * 1.0E6) + ru.ru_utime.tv_usec -
      start_time;
    return elapsed;
  }
}
#endif /* crayc90 */
#endif /* NOCLOCK */

/************************************************************************
 * char* make_path_to_file
 *
 * Concatenates the file name to the directory name
 * Caller is responisble for freeing concatenated string.
 * This should be suitable for an UNIX like OS, including cygwin
 * but would be inappropriate for native Windows.
 *
 * RETURN: a pointer to the concatenated string.
 ************************************************************************/
char* make_path_to_file(const char *directory, const char* filename) {
    size_t directory_length = strlen(directory);
    size_t filename_length = strlen(filename);
    size_t max_path_length = directory_length + filename_length  + 2;
    char *path = mymalloc(max_path_length);
    strncpy(path, directory, max_path_length);
    // Check if the terminating '/' is present in the directory
    // if not add it to the path.
    if (path[directory_length - 1] != '/') {
      strncat(path, "/", max_path_length - directory_length);
    }
    strncat(path, filename, max_path_length - directory_length - 1);
    return path;
}


/************************************************************************
 * char* concat
 *
 * Returns the concatenated string of two given strings
 *
 * RETURN: a pointer to the concatenated string.
 ************************************************************************/
char* concat_string(const char *string1, const char* string2) {
    size_t string1_len = strlen(string1);
    size_t string2_len = strlen(string2);
    size_t max_length = string1_len + string2_len  + 1;
    char *c = mymalloc(max_length);
    strncpy(c, string1, max_length);
    strncat(c, string2, max_length - string1_len);
    return c;
}

/************************************************************************
 * See .h file for description.
 ************************************************************************/
BOOLEAN_T open_file
  (const char *    filename,            /* Name of the file to be opened. */
   const char *    file_mode,           /* Mode to be passed to fopen. */
   BOOLEAN_T allow_stdin,         /* If true, filename "-" is stdin. */
   const char *    file_description,
   const char *    content_description,
   FILE **         afile)               /* Pointer to the open file. */
{
  if (filename == NULL) {
    fprintf(stderr, "Error: No %s filename specified.\n", file_description);
    return(FALSE);
  } else if ((allow_stdin) && (strcmp(filename, "-") == 0)) {
    if (strchr(file_mode, 'r') != NULL) {
      fprintf(stderr, "Reading %s from stdin.\n", content_description);
      *afile = stdin;
    } else if (strchr(file_mode, 'w') != NULL) {
      fprintf(stderr, "Writing %s to stdout.\n", content_description);
      *afile = stdout;
    } else {
      fprintf(stderr, "Sorry, I can't figure out whether to use stdin ");
      fprintf(stderr, "or stdout for %s.\n", content_description);
      return(FALSE);
    }
  } else if ((*afile = fopen(filename, file_mode)) == NULL) {
    fprintf(stderr, "Error opening file %s.\n", filename);
    return(FALSE);
  }
  return(TRUE);
}

/*************************************************************************
 * Add one string to the end of another using the index at to indicate where
 * the end of the string to be added to is.
 *************************************************************************/
static inline void add(char *to, char *from, int max, int *at) {
  assert(to != NULL);
  assert(from != NULL);
  assert(max > 0);
  assert(at != NULL);
  assert(*at >= 0);
  to += (*at);
  max -= 1;
  while (*from != '\0' && *at < max) {
    *to = *from;
    *at += 1;
    ++to;
    ++from;
  }
  *to = '\0';
}

/*************************************************************************
 * Run a program from a given directory with given arguments.  Return
 * the resulting pipe.
 *************************************************************************/
static FILE* run_program
  (char*      program,     // The program to run in the pipe.
   char*      directory,   // Directory where program resides.
   char*      arguments,   // Program arguments.
   char*      type)        // Read ("r") or write ("w").
{
  char* command;
  FILE* return_value;
  int size, pos;

  size = (strlen(directory) + strlen(program) + strlen(arguments) + 3);

  // Allocate space for the command.
  command = (char*)mm_malloc(sizeof(char) * size);
  //Create the path
  pos = 0;
  add(command, directory, size, &pos);
  if (pos > 0 && command[pos-1] != '/') {
    add(command, "/", size, &pos);
  }
  add(command, program, size, &pos);
  if (access(command, F_OK | X_OK) == 0) {
    // Formulate the command
    add(command, " ", size, &pos);
    add(command, arguments, size, &pos);
    // Run the program.
    return_value = popen(command, type);
  } else {
    return_value = NULL;
  }
  myfree(command);
  return(return_value);
}


/*************************************************************************
 * Attempt to run a given program in a given directory with the given
 * arguments, and check that it gives the expected one-line reply.
 *************************************************************************/
static BOOLEAN_T try_to_run
  (char*      program,          // The program to run in the pipe.
   char*      directory,        // Directory to look in.
   char*      test_arguments,   // Arguments used when searching for program.
   char*      expected_reply)   // Expected reply from search.
{
  char* reply;
  FILE* pipe;
  BOOLEAN_T return_value;


  // Allocate space for the reply.
  reply = (char*)mm_malloc(sizeof(char) * (strlen(expected_reply) + 1));

  // Run the command.
  pipe = run_program(program, directory, test_arguments, "r");

  // Check the pipe.
  if (pipe == NULL) {
    return_value = FALSE;
  } else {

    // Read from the pipe.
    if (fgets(reply, strlen(expected_reply) + 1, pipe) == NULL) {
      return_value = FALSE;
    } else {
      return_value = (strcmp(reply, expected_reply) == 0);
    }

    // Close the pipe.
    if (pclose(pipe) == -1) {
      return_value = FALSE;
    }
  }


  myfree(reply);
  return(return_value);
}


/*************************************************************************
 * Open a read-only pipe using a given command line.
 *************************************************************************/
FILE* open_command_pipe
  (char*     program,          // The program to run in the pipe.
   char*     directory,        // Directory to look in.
   char*     test_arguments,   // Arguments used when searching for program.
   char*     expected_reply,   // Expected reply from search.
   char*     real_arguments,   // Arguments used when running the program.
   BOOLEAN_T stdout_on_error,  // If command fails, return STDOUT?
   char*     error_message)    // Error or warning if command fails.
{
  FILE* return_value;

  // Try to run the command with no directory specified.
  if (try_to_run(program, "", test_arguments, expected_reply)) {
    return_value = run_program(program, "", real_arguments, "w");
  }

  // Try to run the program in the specified directory.
  else if (try_to_run(program, directory, test_arguments, expected_reply)) {
    return_value = run_program(program, directory, real_arguments, "w");

  } else {

    // If we failed, print the error message.
    fprintf(stderr, "%s", error_message);
    if (stdout_on_error) {
      return_value = stdout;
    } else {
      exit(1);
    }
  }

  return(return_value);
}


/********************************************************************
 * See .h file for description.
 ********************************************************************/
void die
  (char *format,
   ...)
{
  va_list  argp;

  fprintf(stderr, "FATAL: ");
  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  fprintf(stderr, "\n");
  fflush(stderr);

#ifdef DEBUG
  abort();
#else
  exit(1);
#endif
}


/**************************************************************************
 * See .h file for description.
 **************************************************************************/
void myassert
  (BOOLEAN_T die_on_error,
   BOOLEAN_T test,
   char * const    format,
   ...)
{
  va_list  argp;

  if (!test) {

    if (die_on_error) {
      fprintf(stderr, "FATAL: ");
    } else {
      fprintf(stderr, "WARNING: ");
    }

    /* Issue the error message. */
    va_start(argp, format);
    vfprintf(stderr, format, argp);
    va_end(argp);
    fprintf(stderr, "\n");
    fflush(stderr);

    if (die_on_error) {
#ifdef DEBUG
      abort();
#else
      exit(1);
#endif
    }
  }
}




/********************************************************************
 * void mm_malloc, mm_calloc, mm_realloc
 *
 * See .h file for descriptions.
 ********************************************************************/
void *mm_malloc
  (size_t size)
{
  void * temp_ptr;

  if (size == 0)
    size++;

  temp_ptr = malloc(size);

  if (temp_ptr == NULL)
    die("Memory exhausted.  Cannot allocate %d bytes.", (int)size);

  return(temp_ptr);
}

void *mm_calloc
  (size_t nelem,
   size_t size)
{
  void * temp_ptr;

  /* Make sure we allocate something. */
  if (size == 0) {
    size = 1;
  }
  if (nelem == 0) {
    nelem = 1;
  }

  temp_ptr = calloc(nelem, size);

  if (temp_ptr == NULL)
    die("Memory exhausted.  Cannot allocate %d bytes.", (int)size);

  return(temp_ptr);
}

void * mm_realloc
  (void * ptr,
   size_t  size)
{
  void * temp_ptr;

  /* Make sure we allocate something. */
  if (size == 0)
    size = 1;
  assert(size > 0);

  /* Some non-ANSI systems complain about reallocating NULL pointers. */
  if (ptr == NULL) {
    temp_ptr = malloc(size);
  } else {
    temp_ptr = realloc(ptr, size);
  }

  if (temp_ptr == NULL)
    die("Memory exhausted.  Cannot allocate %d bytes.", (int)size);

  return(temp_ptr);
}

#ifdef MYRAND
#define MY_RAND_MAX 4096

/********************************************************************
 * Primary function for the built-in random number generator.
 ********************************************************************/
static double my_rand
  (long seed)
{
  static long stored_seed = 0;

  /* If this is the first call, just set the seed. */
  if (stored_seed == 0) {
    stored_seed = seed;
  }

  /* Otherwise, create a new pseudorandom number. */
  else {
    stored_seed = abs((stored_seed / 3) * stored_seed + 7718);
  }

  /* Make sure the pseudorandom number is in the right range. */
  return((double)(stored_seed % MY_RAND_MAX) / (double)MY_RAND_MAX);
}
#else
/* The stupid include file doesn't have these prototypes. */
void srand48();
double drand48();

#endif

/********************************************************************
 * See .h file for description.
 ********************************************************************/
void my_srand
  (long seed)
{
#ifdef MYRAND
  my_rand(seed);
#else
  srand48(seed);
#endif
}

/********************************************************************
 * See .h file for description.
 ********************************************************************/
double my_drand
  (void)
{
#ifdef MYRAND
  return(my_rand(0));
#else
  return(drand48());
#endif
}

/**********************************************************************
 * Compute a logarithm.
 **********************************************************************/
PROB_T my_log
  (PROB_T x)
{
  if (x > 0.0) {
    return(LOG_VALUE(log(x)));
  } else if (x < 0.0) {
    die("Tried to take the log of a negative value (%g).", x);
  } /* else if (x == 0.0) */
  return(LOG_ZERO);
}

/* The lookup table. */
#define LOG_PRECISION 1.0e5
static PROB_T log_table[(int) LOG_PRECISION + 2];

/**********************************************************************
 * Set up lookup table for log(x), 0 <= x <= 1.
 **********************************************************************/
void init_log_prob
  (void)
{
  int    i_table;
  PROB_T table_value;

  log_table[0] = LOG_ZERO;
  for (i_table = 1; i_table <= LOG_PRECISION; i_table++) {
    table_value = (double)(i_table / LOG_PRECISION);
    log_table[i_table] = log(table_value);
  }
  log_table[i_table] = 0;  /* For use in iterpolation when x=1 */
}

/**********************************************************************
 * Efficiently find log(x), when 0 < x <= 1.  Doesn't check bounds.
 **********************************************************************/
PROB_T log_prob
  (PROB_T value)
{
  const PROB_T scaled_value = value * LOG_PRECISION;
  const int    log_index = (int)scaled_value;
  const PROB_T decimal_part = scaled_value - log_index;
  const PROB_T lower_value = log_table[log_index];
  const PROB_T upper_value = log_table[log_index+1];
  const PROB_T interpolation = decimal_part * (lower_value - upper_value);

  if (value == 0.0) {
    return(LOG_ZERO);
  }
  return(lower_value + interpolation);
}

/**************************************************************************
 * Return the nearest double smaller than the given double
 **************************************************************************/
double get_next_smaller_double(double x) {
  // IEEE doubles run in lexigraphical order
  // so if we want the next smaller double
  // we just need to cast to integer type and
  // decrement.
  *(long long *) &x = *(long long *) &x - 1;
  return x;
}

/**************************************************************************
 * Return the nearest double larger than the given double
 **************************************************************************/
double get_next_larger_double(double x) {
  // IEEE doubles run in lexigraphical order
  // so if we want the next larger double
  // we just need to cast to integer type and
  // increment.
  *(long long *) &x = *(long long *) &x + 1;
  return x;
}


/**************************************************************************
 * See .h file for description.
 **************************************************************************/
BOOLEAN_T is_zero
  (double    value,
   BOOLEAN_T log_form)
{
  if ((log_form) && (value < LOG_SMALL)) {
    return(TRUE);
  } else if ((!log_form) && (value == 0.0)) {
    return(TRUE);
  } else {
    return(FALSE);
  }
}

/**************************************************************************
 * See .h file for description.
 **************************************************************************/
BOOLEAN_T almost_equal
  (double value1,
   double value2,
   double slop)
{
  if ((value1 - slop > value2) || (value1 + slop < value2)) {
    return(FALSE);
  } else {
    return(TRUE);
  }
}

/*************************************************************************
 * Convert a boolean to and from a "true" or "false" string.
 *************************************************************************/
const char* boolean_to_string
 (BOOLEAN_T the_boolean) {
  return the_boolean == TRUE ? "true" : "false";
}

BOOLEAN_T boolean_from_string
  (char* true_or_false)
{
  if (strcmp(true_or_false, "true") == 0) {
    return(TRUE);
  } else if (strcmp(true_or_false, "false") == 0) {
    return(FALSE);
  } else {
    die("Invalid input to boolean_from_string (%s)\n", true_or_false);
  }
  return(FALSE); /* Unreachable. */
}


/**************************************************************************
 * Does a given character appear in a given string?
 **************************************************************************/
BOOLEAN_T char_in_string
  (const char* a_string,
   char        a_char)
{
  int  i_string;    /* Index into the string. */
  char string_char; /* Character appearing at that index. */

  i_string = 0;
  string_char = a_string[i_string];
  while (string_char != '\0') {
    if (string_char == a_char) {
      return(TRUE);
    }
    i_string++;
    string_char = a_string[i_string];
  }
  return(FALSE);
}

/**************************************************************************
 * Generic functions for converting between integer and string
 * representations of an enumerated type.
 *
 * Assumes that the longest string representation of the enumerated
 * type does not exceed 100 characters.
 *
 * Assumes that the zeroth enumerated type element is invalid.
 **************************************************************************/
char * convert_enum_type
  (int     enum_type, /* The enumerated type object to be converted. */
   char *  enum_strs[], /* String values associated with this type. */
   int     num_enums) /* Number of values of the type. */
{
  if ((enum_type <= 0) || (enum_type >= num_enums)) {
    die("Illegal enumerated type value (%d).", enum_type);
  }

  return(enum_strs[enum_type]);
}

int convert_enum_type_str
  (char *  enum_type_str, /* String to be converted. */
   int     default_value, /* Value to return if first arg is null. */
   char ** enum_strs,     /* String values associated with this type. */
   int     num_enums)     /* Number of values of the type. */
{
  int i_enum;

  /* If no string was given, return the default. */
  if (enum_type_str == NULL) {
    return(default_value);
  }

  /* Search for the value corresponding to the given string. */
  for (i_enum = 0; i_enum < num_enums; i_enum++) {
    if (strcmp(enum_type_str, enum_strs[i_enum]) == 0) {
      return(i_enum);
    }
  }
  die("Illegal value (%s).", enum_type_str);
  return(0); /* Unreachable. */
}

/****************************************************************************
 * Get the name of the CPU.
 ****************************************************************************/
#define HOST_LENGTH 100
const char* hostname
  ()
{
  FILE *           hostname_stream;
  static char      the_hostname[HOST_LENGTH];
  static BOOLEAN_T first_time = TRUE;
  int              num_scanned;

  if (first_time) {
    hostname_stream = (FILE *)popen("hostname", "r"); /* SGI needs cast. */
    num_scanned = fscanf(hostname_stream, "%s", the_hostname);
    assert(num_scanned == 1);
    num_scanned = pclose(hostname_stream);
    assert(num_scanned == 0);
  }
  return(the_hostname);
}

/****************************************************************************
 * Get the current date and time.
 ****************************************************************************/
const char* date_and_time
  ()
{
  FILE *           date_stream;
  static char      the_date[HOST_LENGTH];
  static BOOLEAN_T first_time = TRUE;

  if (first_time) {
    date_stream = (FILE *)popen("date", "r"); /* SGI needs cast. */
    fgets(the_date, HOST_LENGTH, date_stream);
    pclose(date_stream);
  }

  /* Remove the EOL. */
  assert(the_date[strlen(the_date)-1] == '\n');
  the_date[strlen(the_date)-1] = '\0';

  return(the_date);
}

/****************************************************************************
 * Get the last modified time for a file.
 * Date buffer must be at least 26 char long.
 ****************************************************************************/
char *get_last_modified_time(char *filename, char *date_buffer) {

  char *date_string = NULL;
  struct stat stbuf; // buffer for stat call
  stat(filename, &stbuf);
  date_string = ctime_r(&stbuf.st_mtime, date_buffer);

  if (date_string) {
    // Remove newline
    int len = strlen(date_string);
    assert(date_string[len - 1] == '\n');
    date_string[len - 1] = '\0';
  }

  return date_string;
}

/****************************************************************************
 * Copy a string, with allocation.
 ****************************************************************************/
void copy_string
 (char** target,
  char*  source)
{
  if (source == NULL) {
    *target = NULL;
  } else {
    *target = (char *)mm_calloc(strlen(source) + 1, sizeof(char));
    strcpy(*target, source);
  }
}

/************************************************************************
 * Test whether a string is empty or contains only white space
 * Assumes a null terminated string.
 ************************************************************************/
BOOLEAN_T is_empty_string
  (char *s)
{
  BOOLEAN_T result = TRUE;

  assert(s != NULL);

  while(*s != '\0') {
    if (isspace(*s)) {
      s++;
    } else {
      // Found a non-white space character
      result = FALSE;
      break;
    }
  };

  return result;
}

/************************************************************************
 * Copy an array of integers.
 ************************************************************************/
void copy_int_array
 (int  nelems,
  int* source,
  int* target)
{
  int i;

  for (i = 0; i < nelems; i++)
    target[i] = source[i];
}

/************************************************************************
 * Read characters from a file until a whitespace character is encounterd
 ************************************************************************/
void get_non_blank
(FILE * infile,
 char * a_char)
{
  do {
    *a_char = getc(infile);
    if (*a_char == EOF) {
      die("Premature end of file.\n");
    }
  } while ((*a_char == ' ') || (*a_char == '\n') || (*a_char == '\t'));
}


/************************************************************************
 *  Is the next character in the file the eoln character
 ************************************************************************/
BOOLEAN_T is_eoln
  (FILE * infile)
{
  register long ch;

  ch = getc(infile);
  if (ch == EOF) {
    return(TRUE);
  }
  ungetc(ch, infile);
  return(ch == '\n');
}

/************************************************************************
 *  Does a file exist?
 ************************************************************************/
BOOLEAN_T file_exists(char* file_path) {

  struct stat stat_buffer;
  BOOLEAN_T result = FALSE;

  // Does the file exist?
  if (stat(file_path, &stat_buffer)) {
    if (errno == ENOENT) {
      // stat failed because the path doesn't exist.
      result = FALSE;
    }
    else {
      // stat failed for some other reason
      die(
        "Unable to check for status of file '%s'.\n"
        "Error: %s.\n",
        file_path,
        strerror(errno)
      );
    }
  }
  else {
    result = TRUE;
  }

  return result;
}

/************************************************************************
 *  Is a file executable?
 ************************************************************************/
BOOLEAN_T file_executable(char* file_path) {

  struct stat file_status;
  BOOLEAN_T result = FALSE;

  // Is there a path at all?
  if (file_path == NULL) {
    return FALSE;
  }
  
  // Does the file exist?
  if (stat(file_path, &file_status)) {
    if (errno == ENOENT) {
      // stat failed because the path doesn't exist.
      return FALSE;
    } else {
      // stat failed for some other reason
      die(
        "Unable to check for status of file '%s'.\n"
        "Error: %s.\n",
        file_path,
        strerror(errno)
      );
    }
  }
  //stat worked, now check if it's a regular executable
  if (!(S_ISREG(file_status.st_mode) && access(file_path, X_OK) == 0)) return FALSE;
  return TRUE;
}


/******************************************************************
 * This function copies the contents of the file named by the first
 * argument into the file named by the second argument.
 *
 * Returns TRUE if successful, FALSE otherwise.
******************************************************************/
BOOLEAN_T copy_file(char *from_filename, char *to_filename) {

  FILE *from = NULL;
  FILE *to = NULL;
  int c = 0;
  int result = FALSE;

  // Open the files
  from = fopen(from_filename, "r");
  if (!from) {
    // Couldn't open the input file.
    fprintf(stderr, "Couldn't open %s for input.\n", from_filename);
    return FALSE;
  }
  to = fopen(to_filename, "w");
  if (!to) {
    // Couldn't open the output file.
    fprintf(stderr, "Couldn't open %s for output.\n", to_filename);
    fclose(from);
    return FALSE;
  }

  // Copy the bytes
  while ((c = fgetc(from)) != EOF) {
    c = fputc(c, to);
    if (c == EOF) {
      break; // Error writing byte
    }
  }

  // Should be at EOF for input with no errors on output.
  if (feof(from) && !ferror(to)) {
    // Success
    result = TRUE;
  }
  else {
    //  Failure
    if (ferror(from)) {
      fprintf(stderr, "Error reading from %s.\n", from_filename);
    }
    if (ferror(to)) {
      fprintf(stderr, "Error writing to %s.\n", to_filename);
    }
    result = FALSE;
  }

  fclose(from);
  fclose(to);
  return result;

}


/******************************************************************
 * This function checks the environment for a variable pointing
 * to where the MEME etc files are installed. If the environment
 * variable is not found, it returns the definitions from the configuration.
******************************************************************/
const char* get_meme_etc_dir() {
  extern char **environ;
  char **var = environ;
  const char *var_name = "MEME_ETC_DIR=";
  const size_t var_name_length = strlen(var_name);
  while(*var != NULL) {
    // Look for MEME_ETC_DIR variable
    if (strncmp(*var, var_name, var_name_length) == 0) {
      return *var + var_name_length;
    }
    var++;
  }
  return ETC_DIR;
}

/******************************************************************
 * This function checks the environment for a variable pointing
 * to where the MEME executable files are installed. If the environment
 * variable is not found it returns the definition from the configuration.
******************************************************************/
const char* get_meme_bin_dir() {
  extern char **environ;
  char **var = environ;
  const char *var_name = "MEME_BIN_DIR=";
  const size_t var_name_length = strlen(var_name);
  while(*var != NULL) {
    // Look for MEME_BIN_DIR variable
    if (strncmp(*var, var_name, var_name_length) == 0) {
      return *var + var_name_length;
    }
    var++;
  }
  return BIN_DIR;
}

/******************************************************************
 *
 * This function creates a string from the command-line arguments.
 * Caller is responsible for freeing the string returned.
 *
******************************************************************/
char *get_command_line(int argc, char* argv[]) {

  int buffer_size = 200; // Allocate enough for a typical command line
  char *command_line = mm_malloc(sizeof(char) * buffer_size);
  command_line[0] = 0;
  int total_arg_length = 0;
  int i = 0;
  // Work through each of the arguments
  for (i = 0; i < argc; i++) {
    int arg_length = strlen(argv[i]) + 2; // +1 for leading space,
                                          // +1 for trailing null
    total_arg_length += arg_length;
    // Do we have space left in the buffer?
    if (total_arg_length > buffer_size) {
      buffer_size = 2 * total_arg_length;
      command_line =
        mm_realloc(command_line, buffer_size * sizeof(char));
    }
    if (i > 0) {
      // Add leading space if not in the first arugment.
      strcat(command_line, " ");
    }
    strcat(command_line, argv[i]);
  }

  return command_line;
}

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
