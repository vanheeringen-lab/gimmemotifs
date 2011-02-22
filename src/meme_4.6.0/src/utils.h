/********************************************************************
 * FILE: utils.h
 * AUTHOR: William Stafford Noble
 * CREATE DATE: 9-8-97
 * PROJECT: shared
 * COPYRIGHT: 1997-2008, WSN
 * DESCRIPTION: Various useful generic utilities.
 ********************************************************************/
#ifndef UTILS_H
#define UTILS_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "macros.h"

#define FALSE 0
#define TRUE 1
typedef short BOOLEAN_T;

typedef int VERBOSE_T;
#define INVALID_VERBOSE 0
#define QUIET_VERBOSE 1
#define NORMAL_VERBOSE 2
#define HIGH_VERBOSE 3
#define HIGHER_VERBOSE 4
#define DUMP_VERBOSE 5

#define BIG HUGE_VAL
#define LITTLE -BIG

extern VERBOSE_T verbosity;

/*
 * Support LLVM as a compiler.
 */
#ifdef __llvm__
#define inline
#endif

/***********************************************************************
 * Return a not-a-number.
 ***********************************************************************/
double NaN
  (void);

/********************************************************************
 * double myclock
 *
 * Return number of CPU microseconds since first call to myclock().
 * This corrects the bug in the system version of clock that causes it
 * to loop after about 36 minutes.
 *
 * (Taken from Tim Bailey's MEME package.)
 ********************************************************************/
double myclock(void);

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
char* make_path_to_file(const char *directory, const char* filename);

/************************************************************************
 * char* concat
 *
 * Returns the concatenated string of two given strings
 *
 * RETURN: a pointer to the concatenated string.
 ************************************************************************/
char* concat_string( const char *string1, const char* string2);

/************************************************************************
 * int open_file
 *
 * Open a file gracefully.
 *
 * RETURN: Was the open successful?
 ************************************************************************/
BOOLEAN_T open_file
  (const char*     filename,            /* Name of the file to be opened. */
   const char*     file_mode,           /* Mode to be passed to fopen. */
   BOOLEAN_T allow_stdin,         /* If true, filename "-" is stdin. */
   const char*     file_description,
   const char*     content_description,
   FILE **         afile);              /* Pointer to the open file. */

/*************************************************************************
 * Open a write-only pipe using a given command line.
 *
 * The program argument is the name of a target program to be
 * executed.  This function first searches for the target program in
 * the current working directory, then in the specified directory.  If
 * the target program is still not found, this function tries to
 * locate it via the operating system's PATH variable, by calling the
 * target program with the user-provided test arguments.  The results
 * of this call are compared to the user-provided "expected reply,"
 * and if they match, then a pipe is opened with the real arguments.
 *
 * If the program is not found, the function prints a message to
 * stderr and then either aborts or returns the stdout stream,
 * depending upon the value of stdout_on_error.
 *
 * The "expected reply" is assumed to be no more than one line long.
 *************************************************************************/
FILE* open_command_pipe
  (char*     program,          // The program to run in the pipe.
   char*     directory,        // Directory to look in.
   char*     test_arguments,   // Arguments used when searching for program.
   char*     expected_reply,   // Expected reply from search.
   char*     real_arguments,   // Arguments used when running the program.
   BOOLEAN_T stdout_on_error,  // If command fails, return STDOUT?
   char*     error_message);   // Error or warning if command fails.


/********************************************************************
 * DEBUG_CODE (macro)
 *
 * Allow debugging code to be included or excluded from a compiled
 * program.
 ********************************************************************/
#ifdef DEBUG
#define DEBUG_CODE( debug_value, code_fragment ) \
   { if (debug_value) { code_fragment } }
#else
#define DEBUG_CODE( debug_value, code_fragment )
#endif

/********************************************************************
 * void die()
 *
 * Print an error message and die. The arguments are formatted exactly
 * like arguments to printf().
 *
 * (Taken from Sean Eddy's HMMER package.)
 ********************************************************************/
void die
  (char* format,
   ...);

/**************************************************************************
 * Make an assertion, and print the given message if the assertion fails.
 *
 * If the first parameter is set to TRUE, then die if the assertion
 * doesn't go through.  Otherwise, just issue the warning.
 *
 * On exit, dump core if DEBUG is defined.
 **************************************************************************/
void myassert
  (BOOLEAN_T die_on_error,
   BOOLEAN_T test,
   char*  const    format,
   ...);

/********************************************************************
 * Allocate dynamic memory. Die gracefully if memory is exhausted.
 ********************************************************************/
void *mm_malloc
  (size_t size);
void *mm_calloc
  (size_t nelem,
   size_t size);
void * mm_realloc
  (void * ptr,
   size_t size);

/***************************************************************************
 * Dynamically create or grow an array;
 * P = pointer, N = new size, T = type
 **************************************************************************/
#define mm_resize(P,N,T) { \
  void *new_P; \
  new_P = (P) ? realloc((malloc_t)(P), (N)*sizeof(T)) : malloc((N)*sizeof(T)); \
  if (!new_P) { \
    fprintf(stderr, "mm_resize(" #P "," #N "," #T ") failed!\n"); \
    exit(1); \
  } \
  (P) = (T *) new_P; \
}

/********************************************************************
 * Set the seed for the random number generator.
 ********************************************************************/
void my_srand
  (long seed);

/********************************************************************
 * Get a random number X such that 0 <= X < 1.
 ********************************************************************/
double my_drand
  (void);

/********************************************************************
 * Math macros.
 ********************************************************************/
/* Note that the following type must be the  same as the MTYPE and ATYPE
   defined in 'matrix.h' and 'array.h'. */
typedef double PROB_T;       // Type definition for probability/frequency.
#define PROB_SCAN " %lf"     // Scanf string for PROB_T.

#define LOG_ZERO  (-1.0E10)  // Zero on the log scale.
#define LOG_SMALL (-0.5E10)  // Threshold below which everything is zero.
#define MM_BITS      (33.2)     // = LOG2(-LOG_ZERO)

#ifndef MIN
#define MIN(a,b)         (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b)         (((a)>(b))?(a):(b))
#endif

/***************************************************************************
 * Find the nearest integer.
 ***************************************************************************/
#define nint(x) ((int)((x) >= 0 ? ((x) + 0.5) : ((x) - 0.5)))

/**************************************************************************
 * Compute the logarithm of x, when 0 <= x <= 1.
 **************************************************************************/
void init_log_prob
  (void);

PROB_T log_prob
  (PROB_T value);

#define log_prob2(x)   (log_prob(x) * 1.44269504)

/**************************************************************************
 * Compute the logarithm of x.  Returns LOG_ZERO if x==0.
 **************************************************************************/
PROB_T my_log
  (PROB_T x);

#define my_log2(x)   (my_log(x) * 1.44269504)

#define EXP2(x)                          \
( ( (x) < LOG_SMALL) ?                   \
  0.0 :                                  \
  (exp((x) * 0.69314718 ))               \
)

/**************************************************************************
 * Given the logs (in base 2) of two numbers, return the log of their
 * sum.
 *
 * This function is optimized based upon the following formula:
 *
 *      log(x+y) = log(x) + log(1 + exp(log(y) - log(x)))
 *
 **************************************************************************/
#define LOG_VALUE(logx) \
( ( (logx) < LOG_SMALL ) ? \
    LOG_ZERO : \
    (logx) \
)

#define LOG_SUM1(logx, logy) \
( \
  ( ( (logx) - (logy) ) > MM_BITS ) ? \
    LOG_VALUE(logx) : \
    (logx) + my_log2( 1 + EXP2((logy) - (logx) ) ) \
)

#define LOG_SUM(logx, logy) \
( \
  ( (logx) > (logy) ) ? \
    LOG_SUM1( (logx), (logy) ) : \
    LOG_SUM1( (logy), (logx) ) \
)

/**************************************************************************
 * Return the nearest double smaller than the given double
 **************************************************************************/
double get_next_smaller_double(double x);

/**************************************************************************
 * Return the nearest double larger than the given double
 **************************************************************************/
double get_next_larger_double(double x);

/**************************************************************************
 * Test for zero on a value that may be either a log or a raw float.
 **************************************************************************/
BOOLEAN_T is_zero
  (double    value,
   BOOLEAN_T log_form);

/**************************************************************************
 * Test to see if two values are approximately equal.
 **************************************************************************/
BOOLEAN_T almost_equal
  (double value1,
   double value2,
   double slop);

/*************************************************************************
 * Convert a boolean to and from a "true" or "false" string.
 *************************************************************************/
const char*  boolean_to_string(BOOLEAN_T the_boolean);

BOOLEAN_T boolean_from_string(char* true_or_false);

/**************************************************************************
 * Does a given character appear in a given string?
 **************************************************************************/
BOOLEAN_T char_in_string
  (const char* a_string,
   char        a_char);

/**************************************************************************
 * Generic functions for converting between integer and string
 * representations of an enumerated type.
 *
 * Assumes that the longest string representation of the enumerated
 * type does not exceed 100 characters.
 *
 * Assumes that the zeroth enumerated type element is invalid.
 **************************************************************************/
char*  convert_enum_type
  (int     enum_type,  /* The enumerated type object to be converted. */
   char*  enum_strs[],  /* String values associated with this type. */
   int     num_enums); /* Number of values of the type. */

int convert_enum_type_str
  (char*   enum_type_str, /* String to be converted. */
   int     default_value, /* Value to return if first arg is null. */
   char**  enum_strs,     /* String values associated with this type. */
   int     num_enums);    /* Number of values of the type. */

/**************************************************************************
 * Get the name of the CPU.
 **************************************************************************/
const char* hostname
  ();

/**************************************************************************
 * Get the current date and time.
 **************************************************************************/
const char* date_and_time
  ();

/****************************************************************************
 * Get the last modified time for a file.
 * Date buffer must be at least 26 char long.
 ****************************************************************************/
char *get_last_modified_time(char *filename, char *date_buffer);

/**************************************************************************
 * Copy a string, with allocation.
 **************************************************************************/
void copy_string
 (char**  target,
  char*   source);

/************************************************************************
 * Test whether a string is empty or contains only white space
 * Assumes a null terminated string.
 ************************************************************************/
BOOLEAN_T is_empty_string
  (char *s);

/************************************************************************
 * Copy an array of integers.
 ************************************************************************/
void copy_int_array
 (int  nelems,
  int* source,
  int* target);


/************************************************************************
 * Read characters from a file until a whitespace character is encounterd
 ************************************************************************/
void get_non_blank
(FILE * infile,
 char * a_char);

/************************************************************************
 *  Is the next character in the file the eoln character
 ************************************************************************/
BOOLEAN_T is_eoln
  (FILE * infile);

/************************************************************************
 *  Does a file exist?
 ************************************************************************/
BOOLEAN_T file_exists(char* file_path);

/************************************************************************
 *  Is a file executable?
 ************************************************************************/
BOOLEAN_T file_executable(char* file_path);

/******************************************************************
 * This function copies the contents of the file named by the first
 * argument into the file named by the second argument.
 *
 * Returns TRUE if successful, FALSE otherwise.
******************************************************************/
BOOLEAN_T copy_file(char *from_filename, char *to_filename);

/******************************************************************
 * This function checks the environment for a variable pointing
 * to where the MEME etc files are installed. If the environment
 * variable is not found, it returns the definitions from the configuration.
******************************************************************/
const char* get_meme_etc_dir();

/******************************************************************
 * This function checks the environment for a variable pointing
 * to where the MEME executable files are installed. If the environment
 * variable is not found it returns the definition from the configuration.
******************************************************************/
const char* get_meme_bin_dir();

/******************************************************************
 *
 * This function creates a string from the command-line arguments.
 * Caller is responsible for freeing the string returned.
 *
******************************************************************/
char *get_command_line(int argc, char* argv[]);

#endif

