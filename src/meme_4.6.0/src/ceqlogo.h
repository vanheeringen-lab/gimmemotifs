/* ---------------------------- Header ---------------------------------

  Module name   : ceqlogo.h
  Module prefix : CL_

  Description: Create logos for a given set of motifs/PWMs.

               It uses routines and data structures taken from
               MEME files such as

               new-io.*
               matrix.*
               utils.h,
               alphabet.*
               string-list.*
               motif.*
               array.*
               memcheck.*  memory allocation with checking


  ---------------------------------------------------------------------

  Version:

    $Id: ceqlogo.h 5280 2011-01-05 11:10:34Z james_johnson $

  ---------------------------------------------------------------------

  Author   : S. Maetschke

  Copyright: Institute for Molecular Bioscience (IMB)

------------------------------------------------------------------------ */

#ifndef __CEQLOGO
#define __CEQLOGO

/* ----------------------- Global Types -------------------------------- */

/** Set of logo line specific parameters */
typedef struct
{
  MATRIX_T* mat;      /* PWM of the logo */
  char*     alphabet; /* Alphabet used by the PWM */
  int       snum;     /* Number of samples the PWM is derived from. Can be 0 */
  int       ltrim;    /* Number of positions from the left to display trimmed */
  int       rtrim;   /* Number of positions from the right to display trimmed */
  int       shift;    /* horizontal shift of a logo line */
  char*     label;    /* label for the logo line */
} LP_LINE_T;

/** Set of general logo parameters that get replaced within the template file*/
typedef struct
{
  int             initialized;   /* Flag that the parameters are initialized */
  int             is_nucleotide; /* Nucleotide or amino acid logo */
  int             is_reversecomp;/* Reverse complement or normal logo */
  char**          names;         /* parameter names */
  char**          values;        /* parameter values */
  int             p_num;         /* number of parameters (names, values) */
  LP_LINE_T*      lines;         /* logo line specific parameters */
  int             l_num;         /* number of logo lines */
} LP_PARAMS_T;

/*.......................................................................*/


/* ----------------------- Global Prototypes --------------------------- */

void CL_init(LP_PARAMS_T* params, char* type, char* format);
/*
  Description   : Initializes the calculation of a logo.
                  CL_exit() must be called to finish the calculation.

  Parameter     :
    params      : Parameter set.
    type        : Type of the logo: NA or AA. (Nuclic Acid, Amino Acid).
    format      : output format: EPS, SVG

  Return Values :
    params      : Some parameters of the set will be modified.

  Example       :
    CL_init(params, "NA", "EPS");
     ...
    CL_exit(params, "logo.eps");
*/



/*.......................................................................*/

void CL_add_motif(LP_PARAMS_T* params, MOTIF_T *motif, char* alphabet,
                         int shift, char* label);
/*
  Description   : Adds a motif and therefore an additional line to the logo.
                  CL_init() must have been called before.
                  Multiple calls of CL_add_motifs() between the braceing
                  CL_init() and CL_exit() calls are possible.

  Parameter     :
    params      : Parameter set.
    motif       : Motif to add. The number of matrix columns of the contained
                  PWM must be greater than or equal to the size of the alphabet.
    alphabet    : Alphabet used by the PWM/motif.
                  Must match with the type of the logo (see CL_init()).
    shift       : Shifting of the logo. Can be negative.
    label       : label for the logo line. Can be NULL.

  Return Values :
    params      : Some parameters of the set will be modified.

  Example       :
    LP_PARAMS_T* params = LP_create();
      CL_init(params, "NA", "EPS");
        CL_add(params, mat1, "ACTG", 0, NULL);
        CL_add(params, mat2, "ACTG", 0, "Query");
        CL_add(params, mat3, "ACTG", 1, "Target");
      CL_exit(params, "logo.eps");
    LP_free(params);
*/

/*.......................................................................*/

void CL_load_motifs(LP_PARAMS_T* params, char* mode, char *filename, float pseudocounts);
/*
  Description   : Loads motifs from a MEME file and add the contained
                  PWMs to the parameter set.
                  CL_init() and  CL_add() must have been called beforehand.
  Parameter     :
    params      : Parameter set.
    mode        : Loading mode:
                  "-a"  : load all motifs in file
                  "-i"  : load first motif (or use "-i1")
                  "-i2" : load second motif
                  "-i3", ... and so one
    filename    : Name of a MEME file.
    pseudocounts: Pseudo counts when loading frequency matrices within
                  a MEME file. Default value is 1.0

  Return Values :
    params      : The parameter for logo lines will be extended.

  Example       :
    LP_PARAMS_T* params = LP_create();
      CL_init(params, "NA", "EPS");
        CL_load_motifs(params, "-i2", "output.meme", 1.0);
      CL_exit(params, "logo.eps");
    LP_free(params);
*/

/*.......................................................................*/

void CL_exit(LP_PARAMS_T* params, const char* outpath);
/*
  Description   : Writes the logo described by the given parameter set
                  to the given output path.
                  CL_init() and  CL_add() must have been called beforehand.
  Parameter     :
    params      : Parameter set.
    outpath     : Path of the output file.
                  Can be NULL, in this case the output is written to stdout.

  Example       :
    CL_init(params, "NA", "EPS");
      CL_add(params, mat, "ACTG", 0, NULL);
    CL_exit(params, "logo.eps");
*/

/*.......................................................................*/

void CL_create2(
  MOTIF_T *motif1,
  char* label1,
  MOTIF_T *motif2,
  char* label2,
  BOOLEAN_T errbars,                    // use errorbars
  BOOLEAN_T ssc,                        // use small sample correction
  double height,
  double width,
  char* alphabet,
  int shift,
  char* path,
  char* program
);

/*
  Description   : A convenience function that creates an EPS and a PNG output
                  file with one or two logos.
                  The function performs a system call to "convert"!
  Parameter     :
    motif1      : First motif. Can be NULL.
    label1      : Label of the first motif (= logo title). Can be NULL.
    motif2      : Second motif. Can be NULL.
    label2      : Label of the second motif (= label for x-axis). Can be NULL.
    errbars	: print error bars if true
    ssc		: use small sample correction if true
    height	: height of logo in cm.; use default if 0
    width	: width of logo in cm.; use default if 0
    alphabet    : The alphabet used by the motifs. Must no contain ambiguity
                  symbols.
    shift       : Shift of the first logo relative to the second logo. Can be 0.
    path        : Path for the output file WITHOUT an extension.
                  Since an EPS and a PNG file are created, the path should
                  contain only the output folder and the logo name but no
                  file extension. See example.
    program	: name of program to print in fineprint in lower left 
  Global Var.   : Uses CONVERT_PATH

  Example       :
    CL_create2(motif1, "Motif", NULL, NULL, "ACTG", 0, "myfolder/logo");
*/

/*.......................................................................*/

void CL_parse(int argc, char* argv[]);
/*
  Description   : Parses the command line parameters and creates the logo.
                  See usage() for options.
  Parameter     :
    argc        : Number of command line arguments.
    argv        : Array with arguments.

  Example       :
    CL_parse(arg, argv);
*/

/*.......................................................................*/

void CL_run_test(int argc, char* argv[]);
/*
  Description   : Runs all unit tests of the software.
                  Prints out the results on stdout.
  Parameter     :
    argc        : Number of command line arguments.
    argv        : Array with arguments.

  Example       :
    CL_run_test();
*/

/*.......................................................................*/



/*
  Description   :
  Parameter     :
  Global Var.   :
  Return Values :
    function    :
  Example       :
*/


/*.......................................................................*/

void CL_test();
/*
  Description   : Tests the functions of the module.
*/

/*.......................................................................*/


/*.......................................................................*/
// TLB: Combined all modules needed by ceqlogo.
/*.......................................................................*/



/* ----------------------- Global Prototypes --------------------------- */

LP_PARAMS_T* LP_create();
/*
  Description   : Creates a and inits the parameter set.
  Return Values :
    function    : Pointer to the created parameter set. Must be freed!
  Example       :
    LP_PARAMS_T* params = LP_create();
    LP_free(params);
*/

/*.......................................................................*/

void LP_free(LP_PARAMS_T* params);
/*
  Description   : Frees a parameter set.
  Parameter     :
    params      : Parameter set to free.
  Example       :
    LP_PARAMS_T* params = LP_create();
    LP_free(params);
*/


/*.......................................................................*/

void LP_add_line(LP_PARAMS_T* params, MATRIX_T* mat, int snum, int ltrim, 
             int rtrim, char* alphabet, int shift, char* label);
/*
  Description   : Adds an additional line to the logo.

  Parameter     :
    params      : Parameter set.
    mat         : PWM. The number of matrix columns must be greater than or equal
                  to the size of the alphabet.
    snum        : Number of samples the PWM was derived from. Can be 0 and is
                  only used to calc. a small sample correction and error bars.
    ltrim       : Number of positions from the left to display trimmed.
    rtrim       : Number of positions from the right to display trimmed.
    alphabet    : Alphabet used by the PWM.
                  Must match with the type of the logo (see EPS_init()).
    shift       : Shifting of the logo. Can be negative.
    label       : label for the logo line. Can be NULL.

  Return Values :
    params      : Some parameters of the set will be modified.

  Example       :
    LP_PARAMS_T* params = LP_create();
      LP_add_line(params, mat1, 0, 0, 0, "ACTG", 0, NULL);
      LP_add_line(params, mat2, 0, 0, 0, "ACTG", 0, "Query");
      LP_add_line(params, mat3, 0, 0, 0, "ACTG", 1, "Target");
    LP_free(params);
*/

/*.......................................................................*/

LP_LINE_T* LP_get_line(LP_PARAMS_T* params, int line);
/*
  Description   : Getter for the parameters of a logo line.
  Parameter     :
    params      : Parameter set.
    line        : line number.

  Return Values :
    function    : Returns the address of the parameter block for the specified
                  logo line.
  Example       :
    LP_PARAMS_T* params = LP_create();
      LP_LINE_T* line = LP_get_line(params, 1);
    LP_free(params);
*/

/*.......................................................................*/

LP_LINE_T* LP_get_lastline(LP_PARAMS_T* params);
/*
  Description   : Getter for the parameters of the last logo line.
  Parameter     :
    params      : Parameter set.

  Return Values :
    function    : Returns the address of the parameter block for the last
                  logo line.
  Example       :
    LP_PARAMS_T* params = LP_create();
      LP_LINE_T* line = LP_get_lastline(params);
    LP_free(params);
*/

/*.......................................................................*/

size_t LP_max_valuelen(LP_PARAMS_T* params);
/*
  Description   : Determines the maximum value length.
  Parameter     :
    params      : Parameter set.
  Example       :
    LP_PARAMS_T* params = LP_create();
      size_t maxlen = LP_max_valuelen(params);
    LP_free(params);
*/

/*.......................................................................*/

int LP_get_index(LP_PARAMS_T* params, char* name);
/*
  Description   : Getter for the array index of a parameter.
  Parameter     :
    params      : Parameter set.
    name        : Parameter name.
  Return Values :
    function    : index of the parameter or -1 if the no parameter with
                  such a name exists.
  See           : See LP_init_all() for parameter names.
  Example       :
    LP_PARAMS_T* params = LP_create();
      int index = LP_get_index(params, "$TICBITS");
    LP_free(params);
*/

/*.......................................................................*/

char* LP_get_value(LP_PARAMS_T* params, char* name);
/*
  Description   : Getter for the value of a parameter.
                  For a list of supported parameters see function
                  LP_init_all() in logoparam.c
  Parameter     :
    params      : Parameter set.
    name        : Parameter name.
  Return Values :
    function    : Value of the parameter or null if the no parameter with
                  such a name exists.
  Example       :
    LP_PARAMS_T* params = LP_create();
      char* value = LP_get_value(params, "$TICBITS");
    LP_free(params);
*/

/*.......................................................................*/

int LP_get_intvalue(LP_PARAMS_T* params, char* name);
/*
  Description   : Getter for the integer value of a parameter.
                  For a list of supported parameters see function
                  LP_init_all() in logoparam.c
  Parameter     :
    params      : Parameter set.
    name        : Parameter name.
  Return Values :
    function    : Integer value of the parameter or 0 if the no parameter with
                  such a name exists.
  Example       :
    LP_PARAMS_T* params = LP_create();
      int value = LP_get_intvalue(params, "$TICBITS");
    LP_free(params);
*/

/*.......................................................................*/

double LP_get_doublevalue(LP_PARAMS_T* params, char* name);
/*
  Description   : Getter for the double value of a parameter.
                  For a list of supported parameters see function
                  LP_init_all() in logoparam.c
  Parameter     :
    params      : Parameter set.
    name        : Parameter name.
  Return Values :
    function    : Double value of the parameter or 0 if the no parameter with
                  such a name exists.
  Example       :
    LP_PARAMS_T* params = LP_create();
      double value = LP_get_doublevalue(params, "$TICBITS");
    LP_free(params);
*/

/*.......................................................................*/

void LP_set_value(LP_PARAMS_T* params, char* name, char *value);
/*
  Description   : Setter for a parameter value. Note that the new parameter
                  value is copied.
                  For a list of supported parameters see function
                  LP_init_all() in logoparam.c
  Parameter     :
    params      : Parameter set.
    name        : Parameter name.
    value       : New value of the parameter.

  Return Values :
    params      : Parameter value in parameter set is changed.
  Example       :
    LP_PARAMS_T* params = LP_create();
      LP_set_value(params, "$TITLE", "My Title");
    LP_free(params);

*/

/*.......................................................................*/

void LP_set_intvalue(LP_PARAMS_T* params, char* name, int value);
/*
  Description   : Setter for an integer parameter value.
                  For a list of supported parameters see function
                  LP_init_all() in logoparam.c
  Parameter     :
    params      : Parameter set.
    name        : Parameter name.
    value       : New value of the parameter.

  Return Values :
    params      : Parameter value in parameter set is changed.
  Example       :
    LP_PARAMS_T* params = LP_create();
      LP_set_intvalue(params, "$CHARSPERLINE", 12);
    LP_free(params);
*/

/*.......................................................................*/

void LP_set_doublevalue(LP_PARAMS_T* params, char* name, double value);
/*
  Description   : Setter for a double parameter value.
                  For a list of supported parameters see function
                  LP_init_all() in logoparam.c
  Parameter     :
    params      : Parameter set.
    name        : Parameter name.
    value       : New value of the parameter.

  Return Values :
    params      : Parameter value in parameter set is changed.
  Example       :
    LP_PARAMS_T* params = LP_create();
      LP_set_doublevalue(params, "$ERRORBARFRACTION", 0.95);
    LP_free(params);
*/

/*.......................................................................*/

void LP_correct_shift(LP_PARAMS_T* params);
/*
  Description   : Corrects the shift factors for logos to become non-negative.
  Parameter     :
    params      : Parameter set.

  Return Values :
    params      : Shift factors are adjusted if necessary.
  Example       :
    LP_PARAMS_T* params = LP_create();
      LP_correct_shift(params);
    LP_free(params);
*/

/*.......................................................................*/

void LP_translate_cmdline(LP_PARAMS_T* params, int argc, char* argv[]);
/*
  Description   : Translates the arguments of a command line in
                  settings of the logo parameters.
  Parameter     :
    params      : Parameter set.
    argc        : Number of arguments within the command line.
    argv        : An array that contains the arguments of the command line.
  Example       :
    LP_PARAMS_T* params = LP_create();
      char *argv[] = {"-w", "100", "-a"};
      int argc = 3;
      LP_translate_cmdline(params, argc, argv);
    LP_free(params);
*/

/*.......................................................................*/

char *LP_get_option(char *option, int argc, char* argv[]);
/*
  Description   : Gets an option from the command line.
                  If an option starts with hyphen and the following
                  argument of the command line does not, it is assumed
                  that this argument is the value of the option.
                  If an option starts with hyphen and the following
                  argument of the command line also starts with a hyphen,
                  the option itself is returned.
  Parameter     :
    option      : The option to get.
    argc        : Number of arguments within the command line.
    argv        : An array that contains the arguments of the command line.
  Return Values :
    function    : Returns the option, the value of the option or NULL if
                  no such option exists.
  Example       :
    char* value = LP_get_option("-f", argc, argv);
*/

/*.......................................................................*/

char *LP_get_option_default(char *option, char* def, int argc, char* argv[]);
/*
  Description   : Same as LP_get_option() but with a default return value
                  in the case the option is not set.
  Parameter     :
    option      : The option to get.
    def         : Default return value.
    argc        : Number of arguments within the command line.
    argv        : An array that contains the arguments of the command line.
  Return Values :
    function    : Returns the option, the value of the option or default if
                  no such option exists.
  Example       :
    char* value = LP_get_option_default("-o", "output.eps", argc, argv);
*/


/*.......................................................................*/

void LP_file_replace(LP_PARAMS_T* params, FILE *fin, FILE *fout, size_t len);
/*
  Description   : Reads an input file line by line, replaces all parameter
                  names by their values and writes the result to the given
                  output file.
  Parameter     :
    params      : Parameter set.
    fin         : Input file.
    fout        : Output file.
    len         : Maximum line length in input file.

  Return Values :
    function    :
  Example       :
    FILE *fin  = F_open("etc/template.eps", "r");
    FILE *fout = F_open(outpath, "w");
    LP_PARAMS_T* params = LP_create();
      LP_file_replace(params, fin, fout, 1000);
    LP_free(params);
    F_close(fin);
    F_close(fout);
*/

/*.......................................................................*/

char* LP_DNAcolordict();
/*
  Description   : Getter for as string with a color dictionary for DNA symbols.
  Parameter     :
  Return Values :
    params      : String with a color dictionary.
  Example       :
    LP_DNAcolordict();
*/

/*.......................................................................*/

char* LP_AAcolordict();
/*
  Description   : Getter for as string with a color dictionary for amino acid
                  symbols.
  Parameter     :
  Return Values :
    params      : String with a color dictionary.
  Example       :
    LP_AAcolordict();
*/

/*.......................................................................*/

char* LP_colordef();
/*
  Description   : Getter for as string with color definitions.
  Parameter     :
  Return Values :
    params      : String with color definitions.
  Example       :
    LP_colordef();
*/

/*.......................................................................*/


/*
  Description   :
  Parameter     :
  Global Var.   :
  Return Values :
    function    :
  Example       :
*/


/*.......................................................................*/

void LP_test();
/*
  Description   : Tests the functions of the module.
*/

/*.......................................................................*/




/* ----------------------- Global Prototypes --------------------------- */

#define F_open(name, type) _F_open((name),(type), __LINE__, __FILE__)
FILE* _F_open(const char *name, const char *type, int line, const char *file);
/*
  Description   : Opens a file. Alway use the macro.
  Parameter     :
    name        : Name/path of the file.
    type        : File type, e.g. "r"
  Return Values :
    function    : Returns a file handle or prints out an error message and
                  return NULL if the file couldn't be opened.
  Example       :
    FILE *f = F_open("mytext.txt","r");
    F_close(f);
*/

/*.......................................................................*/

#define F_close(fp) _F_close((fp), __LINE__, __FILE__)
void _F_close(FILE *fp, int line, const char *file);
/*
  Description   : Closes a file. Alway use the macro.
  Parameter     :
    fp          : File pointer.
  Example       :
    FILE *f = F_open("mytext.txt","r");
    F_close(f);
*/

/*.......................................................................*/

int F_check();
/*
  Description   : Checks if all files that have been opened were closed
                  as well. Prints out an error message if not all files
                  were closed.
  Global Var.   :
    E_counter   : Counter for open files.

  Return Values :
    function    : TRUE: all files are closed
                  FALSE: not all files are closed.
  Example       :
    FILE *f = F_open("mytext.txt","r");
    F_check();  -> 0
    F_close(f);
    F_check();  -> 1
*/

/*.......................................................................*/

int F_exists(const char *name);

/*
  Description   : Tests if the file with the given name exists.
  Parameter     : Name of the file.
  Return Values :
    function    : Returns 1 if the file exist and 0 otherwise.
  Example       :
    if(F_exists("test.dat"))
      printf("It exists"):
*/

/*.......................................................................*/

void F_rmdir(char* dir, char* mode);
/*
  Description   : Uses system() to remove a directory.
                  Doesn't care if the directory exists or isn't empty!
                  There is no return value, if the operation was successful!
  Parameter     :
    dir         : Directory to remove.
    mode        : Mode, these are the flags used by rm
                  e.g. "-r" recursively removes the contents.
                  can be empty ("") but must not be NULL.
  Example       :
    F_rmdir("myfolder", "-r");
    F_rmdir("myfolder", "");
*/

/*.......................................................................*/

void F_mkdir(char* dir, char* mode);
/*
  Description   : Uses system() to create a directory.
                  There is no return value, if the operation was successful!
  Parameter     :
    dir         : Directory to create.
    mode        : Mode, these are the flags used by mkdir
                  e.g. "-m 777" sets permissions.
                  can be empty ("") but must not be NULL.
  Example       :
     F_mkdir("myfolder", "-m 777");
     F_mkdir("myfolder", "");
*/

/*.......................................................................*/



/*
  Description   :
  Parameter     :
  Global Var.   :
  Return Values :
    function    :
  Example       :
*/


/*.......................................................................*/

void F_test();
/*
  Description   : Tests the functions of the module.
*/

/*.......................................................................*/

/* ----------------------- Global Prototypes --------------------------- */
char* STR_time(char* format);
/*
  Description   : Returns a time & date string according to the specified
                  format.
  Parameter     :
    format      : e.g. %b %d, %Y; %H:%M:%S   (see time.h)
  Return Values :
    function    : String with time & date.
  Example       :
    STR_time(""); -> Jan 10, 1987; 17:55:55
*/


char* STR_copy(char* str);
/*
  Description   : Creates a dynamically allocated copy of the given string.
  Parameter     :
    str         : String to copy.
  Global Var.   :
  Return Values :
    function    : Copy of the string.
  Example       :
    char* str = STR_copy("test");  -> "test"
    free(str);
*/

/*.......................................................................*/


void STR_replace(char buffer[], int n, const char *searchstr,
                 const char *replacestr);
/*
  Description   : Replaces all occurences of a search string by a replace
                  string within a string buffer. The buffer must be
                  sufficiently long.
  Parameter     :
    buffer      : String buffer.
    n           : Length of the buffer.
    searchstr   : The string to search and replace.
    replacestr  : The string that replaces all occurences of the search string
                  within the string buffer.
  Global Var.   :
  Return Values :
    buffer      : The content of the string buffer will be modified.
  Example       :
    char buffer[15];
    STR_replace(strcpy(buffer, "aaabbaaabba"), 15, "bb", "B");
*/

/*.......................................................................*/

int STR_replace_diff(const char *str, const char *searchstr,
                     const char *replacestr);
/*
  Description   : Calculates the difference in length when all occurrences
                  of a search string are replaced by another string.
  Parameter     :
    str         : String to search in.
    searchstr   : The string to search and replace.
    replacestr  : The string that replaces all occurences of the search string.
  Global Var.   :
  Return Values :
    function    : Returns the length difference caused by the string replacement.
  Example       :
    STR_replace_diff("aaabbbaabbb", "bbb", "b");   ->  -4
*/

/*.......................................................................*/

char *STR_append(char **str, char* append);
/*
  Description   : Appends a string to a dynamically allocated string.
  Parameter     :
    str         : Address of string to extend.
    append      : String to append.
  Return Values :
    function    : Returns the new string or NULL if there is insufficent memory.
    str         : Updates the content of the string pointer
  Example       :
    char *str = (char*)M_malloc(20);
    strcpy(str, "blue");
    STR_append(&str, " and");   ->  "blue and"
    STR_append(&str, " red");   ->  "blue and red"
    M_free(str);
*/

/*.......................................................................*/


void STR_test();
/*
  Description   : Tests the functions of the module.
*/

/*.......................................................................*/


/* ----------------------- Global Prototypes --------------------------- */

void T_set_verbose(int verbose_level);
/*
  Description   : Setter for the level of verbosity for the tests.
  Parameter     :
  verbose_level : 0: only failures are reported
                  1: tested file is reported
                  2: file and number of test performed are reported
                  3: the individual tests are reported.
  Example       :
    T_set_verbose(1);
*/

#define T_start() _T_start(__FILE__)
void _T_start(char *filename);
/*
  Description   : Starts the test.
  Parameter     :
    filename  : Name of the module/file to test.
  Example       :
    T_start(__FILE__);
    T_end();
*/

/*.......................................................................*/

void T_end();
/*
  Description   : Ends the test.
  Parameter     :
  Example       :
    T_start(__FILE__);
    T_end();
*/


/*.......................................................................*/


#define T_is_null(actual) _T_is_null(#actual, (actual), __LINE__, __FILE__)
void _T_is_null(char* expression, void* actual, int line, char *file);
/*
  Description   : Tests for NULL.
  Parameter     :
    acutal      : Actual value.
  Example       :
    T_is_null(p);
*/

/*.......................................................................*/

#define T_boolean(expected, actual) _T_boolean(#actual, (expected), (actual), __LINE__, __FILE__)
#define T_is_true(actual) _T_boolean(#actual, 1, (actual), __LINE__, __FILE__)
#define T_is_false(actual) _T_boolean(#actual, 0, (actual), __LINE__, __FILE__)
void _T_boolean(char* expression, int expected, int actual, int line, char *file);
/*
  Description   : Tests boolean values.
  Parameter     :
      expected  : Expected value.
      actual    : Actual value.
  Example       :
    T_boolean(1, 0);
*/

/*.......................................................................*/

#define T_string(expected, actual) _T_string(#actual, (expected), (actual), __LINE__, __FILE__)
void _T_string(char* expression, char* expected, char* actual, int line, char *file);
/*
  Description   : Tests string values.
  Parameter     :
      expected  : Expected value.
      actual    : Actual value.
  Example       :
    T_string("hello", "world");
*/

/*.......................................................................*/


#define T_int(expected, actual) _T_int(#actual, (expected), (actual), __LINE__, __FILE__)
void _T_int(char* expression, int expected, int actual, int line, char *file);
/*
  Description   : Tests integer values.
  Parameter     :
      expected  : Expected value.
      actual    : Actual value.
  Example       :
    T_int(11, 234);
*/

/*.......................................................................*/

#define T_char(expected, actual) _T_char(#actual, (expected), (actual), __LINE__, __FILE__)
void _T_char(char* expression, char expected, char actual, int line, char *file);
/*
  Description   : Tests char values.
  Parameter     :
      expected  : Expected value.
      actual    : Actual value.
  Example       :
    T_char('a', 'c');
*/

/*.......................................................................*/

#define T_double(expected, actual, eps) _T_double(#actual, (expected), (actual), eps, __LINE__, __FILE__)
void _T_double(char* expression, double expected, double actual, double eps, int line, char *file);
/*
  Description   : Tests boolean values.
  Parameter     :
      expected  : Expected value.
      actual    : Actual value.
      esp       : Permitted difference between the two values.
  Example       :
    T_double(1.00, 1.01, 0.01);
*/

/*.......................................................................*/




/*.......................................................................*/


/*
  Description   :
  Parameter     :
  Global Var.   :
  Return Values :
    function    :
  Example       :
*/


/*.......................................................................*/


/*.......................................................................*/


/* ----------------------- Global Types -------------------------------- */

/** Describes a letter and its heigth within the letter stack of a logo */
typedef struct lstack_t {
  double height;
  char   letter;
} LSTACK_T;

/*.......................................................................*/




/* ----------------------- Global Variables ---------------------------- */


/* ----------------------- Global Prototypes --------------------------- */

LSTACK_T* LS_create(size_t size);
/*
  Description   : Creates a letter stack data structure.
                  Don't forget to free the letter stack: LS_free().
  Parameter     :
    size        : Size of the letter stack (=alphabet size/number of
                  PWM columns).
  Return Values :
    function    : Returns a letter stack of the specified size.
  Example       :
    LSTACK_T* lstack = LS_create(4);
    LS_free(lstack);
*/

/*.......................................................................*/

void LS_free(LSTACK_T* lstack);
/*
  Description   : Frees the letter stack.
  Parameter     :
    lstack      : Letter stack.
  Example       :
    LSTACK_T* lstack = LS_create(4);
    LS_free(lstack);
*/

/*.......................................................................*/

void LS_calc(int row, size_t n, MATRIX_T* mat, char* alphabet, 
  LSTACK_T* lstack, BOOLEAN_T ssc);
/*
  Description   : Fills the letter stack (= calculates the height of
                  stack letters and sorts the stack according to
                  height).
  Parameter     :
    row         : Matrix row. First row is zero.
    n           : Number of samples. Zero is valid - in this case the
                  correction factor is zero as well.
    mat         : PWM.
    alphabet    : String with uppercase alphabet letters (no ambiguity letters).
    lstack      : Letter stack. Must be of appropriate size!
    ssc		: use small sample correction

  Return Values :
      lstack    : The content of the letter stack is filled.
  Example       :
    LSTACK_T* lstack = LS_create(4);
      LS_calc(0,0,mat,"ACTG",lstack);
    LS_free(lstack);

*/

/*.......................................................................*/

double LS_calc_e(size_t n, size_t asize);
/*
  Description   : Calculate the error correction factor for small numbers
                  of samples (= number of aligned sequences).
  Parameter     :
    n           : Number of samples. Zero is valid - in this case the
                  correction factor is zero as well.
    asize       : Size of the alphabet.
  Return Values :
    function    : Returns the error correction factor.
  Example       :
    double e = LS_calc_e(0,4);
*/

/*.......................................................................*/


/*
  Description   :
  Parameter     :
  Global Var.   :
  Return Values :
    function    :
  Example       :
*/

/*.......................................................................*/



/*
  Description   :
  Parameter     :
  Global Var.   :
  Return Values :
    function    :
  Example       :
*/

/*.......................................................................*/

void LS_test();
/*
  Description   : Tests the functions of the module.
*/

/*.......................................................................*/


void SVG_init(LP_PARAMS_T* params, char* type);
/*
  Description   : Initializes the calculation of a logo in SVG format.
                  SVG_exit() must be called later.

  Parameter     :
    params      : Parameter set.
    type        : Type of the logo: NA or AA. (Nuclic Acid, Amino Acid).

  Return Values :
    params      : Some parameters of the set will be modified.

  Example       :
    SVG_init(params, "NC");
     ...
    SVG_exit(params, fout);
*/

/*.......................................................................*/

void SVG_add(LP_PARAMS_T* params, MATRIX_T* mat, int snum, int ltrim, 
             int rtrim, char* alphabet, int shift, char* label);
/*
  Description   : Adds a PWM and therefore an additional line to the logo.
                  SVG_init() must have been called at some stage.
                  Multiple calls of SVG_add() between the braceing
                  SVG_init() and SVG_exit() calls are possible.

  Parameter     :
    params      : Parameter set.
    mat         : PWM. The number of matrix columns must be greater than or equal
                  to the size of the alphabet.
    snum        : Number of samples the PWM was derived from. Can be 0 and is
                  only used to calc. a small sample correction and error bars.
    ltrim       : Number of positions to display trimmed from the left.
    rtrim       : Number of positions to display trimmed from the right.
    alphabet    : Alphabet used by the PWM.
                  Must match with the type of the logo (see SVG_init()).
    shift       : Shifting of the logo. Can be negative.
    label       : label for the logo line. Can be NULL.

  Return Values :
    params      : Some parameters of the set will be modified.

  Example       :
    SVG_init(params, "NA");
      SVG_add(params, mat1, 0, 0, 0, "ACTG", 0, NULL);
      SVG_add(params, mat2, 0, 0, 0, "ACTG", 0, "label");
      SVG_add(params, mat3, 1, 0, 0, "ACTG", 2, NULL);
    SVG_exit(params, fout);
*/

/*.......................................................................*/

void SVG_exit(LP_PARAMS_T* params, FILE* fout);
/*
  Description   : Writes the logo described by the given parameter set
                  in SVG format to the given output file.
                  SVG_init() and  SVG_add() must have been called beforehand.
  Parameter     :
    params      : Parameter set.
    fout        : File pointer to open output file.

  Example       :
    FILE* fout = F_open("output.eps", "w");
      SVG_init(params, "NA");
        SVG_add(params, mat1, 0, 0, 0, "ACTG", 0, NULL);
      SVG_exit(params, fout);
    F_close(fout);
*/

/*.......................................................................*/

/*
  Description   :
  Parameter     :
  Global Var.   :
  Return Values :
    function    :
  Example       :
*/




/*.......................................................................*/

void SVG_test();
/*
  Description   : Tests the functions of the module.
*/

/*.......................................................................*/


#define E_system  E_position(__FILE__,__LINE__) , _E_system
void _E_system(char* format, ...);
/*
  Description   : Prints out file name and line number for a
                  system error to stderr.
                  System errors are errors in the usage of code
                  and should not be used to signal input errors
                  of the user. For that purpose use E_error().
                  The function stops the execution of the application
                  and exits with error code -1.

  Parameters    :
    format      : Format string like printf.
    ...         : arguments described in the format string.

  Example       :
    int errlevel;
    E_system("This a system error. Error level %d", errlevel);
*/

/*.......................................................................*/


void E_error(char* format, ...);
/*
  Description   : Prints out an error message to stderr.


  Parameters    :
    format      : Format string like printf.
    ...         : arguments described in the format string.

  Example       :
    int errno;
    E_error("This the error number %d\n", errno);
*/

/*.......................................................................*/

void E_position(char* file, int line);
/*
  Description   : Prints the current file and line number to stderr.
  Parameter     :
    file        : File name.
    line        : Line number.
  Example       :
    E_position(__FILE__,__LINE__);
*/

/*.......................................................................*/

/*
  Description   :
  Parameter     :
  Global Var.   :
  Return Values :
    function    :
  Example       :
*/


/*.......................................................................*/

void E_test();
/*
  Description   : Tests the functions of the module.
*/

/*.......................................................................*/


/* ----------------------- Global Prototypes --------------------------- */

void EPS_init(LP_PARAMS_T* params, char* type);
/*
  Description   : Initializes the calculation of a logo in EPS format.
                  EPS_exit() must be called later.

  Parameter     :
    params      : Parameter set.
    type        : Type of the logo: NA or AA. (Nuclic Acid, Amino Acid).

  Return Values :
    params      : Some parameters of the set will be modified.

  Example       :
    EPS_init(params, "NA");
     ...
    EPS_exit(params, fout);
*/

/*.......................................................................*/

void EPS_add(LP_PARAMS_T* params, MATRIX_T* mat, int snum, int ltrim, 
             int rtrim, char* alphabet, int shift, char* label);
/*
  Description   : Adds a PWM and therefore an additional line to the logo.
                  EPS_init() must have been called at some stage.
                  Multiple calls of EPS_add() between the braceing
                  EPS_init() and EPS_exit() calls are possible.

  Parameter     :
    params      : Parameter set.
    mat         : PWM. The number of matrix columns must be greater than or equal
                  to the size of the alphabet.
    snum        : Number of samples the PWM was derived from. Can be 0 and is
                  only used to calc. a small sample correction and error bars.
    ltrim       : Number of positions to display trimmed from the left.
    rtrim       : Number of positions to display trimmed from the right.
    alphabet    : Alphabet used by the PWM.
                  Must match with the type of the logo (see EPS_init()).
    shift       : Shifting of the logo. Can be negative.
    label       : label for the logo line. Can be NULL.

  Return Values :
    params      : Some parameters of the set will be modified.

  Example       :
    EPS_init(params, "NA");
      EPS_add(params, mat1, 0, 0, 0, "ACTG", 0, NULL);
      EPS_add(params, mat2, 0, 0, 0, "ACTG", 0, "Query");
      EPS_add(params, mat3, 0, 0, 0, "ACTG", 1, "Target");
    EPS_exit(params, fout);
*/

/*.......................................................................*/

void EPS_add_motifs(LP_PARAMS_T* params, MOTIF_T* motifs, int num_motifs, char* alphabet);
/*
  Description   : Adds a set of motifs to the logo.
                  EPS_init() must have been called at some stage.
                  Multiple calls of EPS_add_motifs() between the braceing
                  EPS_init() and EPS_exit() calls are possible.

  Parameter     :
    params      : Parameter set.
    motifs      : Array with motifs.
    num_motifs  : Number of motifs within the array to add.
    alphabet    : Alphabet used by the PWM.
                  Must match with the type of the logo (see EPS_init()).

  Return Values :
    params      : Some parameters of the set will be modified.

  Example       :
    int num_motifs = 0;
    BOOLEAN_T has_reverse_strand = FALSE;
    MOTIF_T motifs[10];
    STRING_LIST_T* motif_occurrences = NULL;
    ARRAY_T* background;
    read_meme_file("output.meme", &num_motifs, motifs, &motif_occurrences,
                   &has_reverse_strand, &background);

     EPS_init(params, "NA");
       EPS_add_motifs(params, motifs, num_motifs, get_alphabet(FALSE));
     EPS_exit(params, fout);
*/


/*.......................................................................*/

void EPS_exit(LP_PARAMS_T* params, FILE* fout);
/*
  Description   : Writes the logo described by the given parameter set
                  in EPS format to the given output path.
                  EPS_init() and  EPS_add() must have been called beforehand.
  Parameter     :
    params      : Parameter set.
    fout        : File pointer to open output file.

  Example       :
    FILE* fout = F_open("output.eps", "w");
      EPS_init(params, "NA");
        EPS_add(params, mat, 0, 0, 0, "ACTG", 0, NULL);
      EPS_exit(params, fout);
    F_close(fout);
*/

/*.......................................................................*/



/*
  Description   :
  Parameter     :
  Global Var.   :
  Return Values :
    function    :
  Example       :
*/


/*.......................................................................*/

void EPS_test();
/*
  Description   : Tests the functions of the module.
*/

/*.......................................................................*/


#endif
