/* ----------------------- Implementation -----------------------------

  Module name   : ceqlogo.c

  Description: Create logos for a given set of motifs/PWMs.

  ---------------------------------------------------------------------

  Version:

    $Id: ceqlogo.c 5280 2011-01-05 11:10:34Z james_johnson $

  ---------------------------------------------------------------------

  Author   : S. Maetschke

  Copyright: Institute for Molecular Bioscience (IMB)

------------------------------------------------------------------------ */


/* ---------------------------- Includes ------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdarg.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <ctype.h>
#include <sys/wait.h>

#include "macros.h"
#include "matrix.h"
#include "alphabet.h"
#include "string-list.h"
#include "motif.h"
#include "utils.h"
#include "meme-io.h"
#include "config.h"
#include "dir.h"

#include "ceqlogo.h"



/* ----------------------- Global Constants ---------------------------- */
/* Max. number of motifs in a MEME file */
#define CL_MAX_MOTIFS 5000		// Fragile!!!
/* Path to the template file */
#define TEMPLATE_SVG ETC_DIR"/template.svg"
/* Path to the template file */
#define TEMPLATE_EPS ETC_DIR"/template.eps"


/* -------------------------- Local Prototypes ------------------------- */

static void CL_usage();
/*
  Description   : Prints the usage information.
  Example       :
    CL_usage();
*/

/*.......................................................................*/

static void CL_error(char* format, ...);
/*
  Description   : Prints out an error message followed by the usage description
                  and then exits with -1.
  Parameter     :
    format      : Format string like printf.
    ...         : arguments described in the format string.

  Example       :
    int errno;
    CL_error("This the error number %d\n", errno);
*/

/*.......................................................................*/


/* ------------------------- Global Functions -------------------------- */

void CL_init(LP_PARAMS_T* params, char* type, char* format) {
  LP_set_value(params, "$FORMAT", format);
  if(!strcmp(format, "EPS"))
    EPS_init(params, type);
  else if(!strcmp(format, "SVG"))
    SVG_init(params, type);
  else
    E_system("Invalid output format: %s\n", format);
}

/*.......................................................................*/

void CL_add_motif(LP_PARAMS_T* params, MOTIF_T *motif, char* alphabet,
                         int shift, char* label) {
  char* format = LP_get_value(params, "$FORMAT");
  assert(params->is_nucleotide || !(params->is_reversecomp));
  if (params->is_reversecomp) reverse_complement_motif(motif);
  if(!strcmp(format, "EPS")) {
    EPS_add(params, motif->freqs, motif->num_sites, motif->trim_left, 
        motif->trim_right, alphabet, shift, label);
  } else if(!strcmp(format, "SVG")) {
    SVG_add(params, motif->freqs, motif->num_sites,motif->trim_left, 
        motif->trim_right, alphabet, shift, label);
  } else {
    E_system("Invalid output format: %s\n", format);
  }
}

/*.......................................................................*/

void CL_exit(LP_PARAMS_T* params, const char* outpath) {
  FILE* fout   = outpath == NULL ? stdout : F_open(outpath, "w");
  char* format = LP_get_value(params, "$FORMAT");

  if(!strcmp(format, "EPS"))
    EPS_exit(params, fout);
  else if(!strcmp(format, "SVG"))
    SVG_exit(params, fout);
  else
    E_system("Invalid output format: %s\n", format);

  if(fout != stdout)
    F_close(fout);
}

/*.......................................................................*/

void CL_create2(
  MOTIF_T *motif1, 
  char* label1, 
  MOTIF_T *motif2, 
  char* label2,
  BOOLEAN_T errbars, 			// use errorbars
  BOOLEAN_T ssc,			// use small sample correction
  double height, 
  double width, 
  char* alphabet, 
  int shift, 
  char* path, 
  char* program
) 
{
  int   asize    = strlen(alphabet);           
  char* epspath = (char*)mymalloc(sizeof(char)*(strlen(path)+10));
  sprintf(epspath, "%s.eps", path);
  char* svgpath = (char*)mymalloc(sizeof(char)*(strlen(path)+10));
  sprintf(svgpath, "%s.svg", path);

  LP_PARAMS_T* params = LP_create();
  /* see LP_init_all() for more parameters */

  if (label1 != NULL) LP_set_value(params, "$TITLE", label1);
  if (label2 != NULL) LP_set_value(params, "$XAXISLABEL", label2);
  LP_set_value(params, "$ERRBAR", errbars ? "true" : "false");
  LP_set_value(params, "$SSC", ssc ? "true" : "false");
  if (height) LP_set_doublevalue(params, "$LOGOHEIGHT", height);
  if (width) LP_set_doublevalue(params, "$LOGOWIDTH", width);
  // Set the program name and time in the fineprint at bottom of logo
  char *time = STR_time("%d.%m.%y %H:%M");
  char *fineprint = (char*)mymalloc(sizeof(char)*(strlen(program)+strlen(time)+10));
  sprintf(fineprint, "%s %s", program, time);
  LP_set_value(params, "$FINEPRINT", fineprint);
  myfree(fineprint);
    
  // Make the EPS format LOGO
  CL_init(params, asize < 20 ? "NA" : "AA", "EPS");
  if (motif1 != NULL) {
    CL_add_motif(params, motif1, alphabet, shift, label1);
  }
  if (motif2 != NULL) {
    CL_add_motif(params, motif2, alphabet, 0,     label2);
  }
  CL_exit(params, epspath);
  LP_free(params);

  // Create PNG format of LOGO if possible
  if (CONVERT_PATH != NULL) {
    struct stat file_status;
    int result = stat(CONVERT_PATH, &file_status);
    if (result == 0 
        && S_ISREG (file_status.st_mode) 
        && access(CONVERT_PATH, X_OK) == 0) {
      int ret;
      int cmd_len = strlen(CONVERT_PATH) + 1 + 2*strlen(epspath) + 2;
      char *cmd = (char*)mymalloc(sizeof(char)*cmd_len);
      sprintf(cmd, "%s %s %s.png", CONVERT_PATH, epspath, path);
      ret = system(cmd);

      if (!(WIFEXITED(ret) && WEXITSTATUS(ret) == 0)) {
        fprintf(stderr, "Warning: conversion to png format failed. Please check that Imagemagick is installed.\n");
      }
      myfree(cmd);
    }
  }
  myfree(epspath);

  // FIXME: Make the SVG format LOGO
  if (0) {
    LP_PARAMS_T* params = LP_create();
    /* see LP_init_all() for more parameters */

    if (label1 != NULL) LP_set_value(params, "$TITLE", label1);
    if (label2 != NULL) LP_set_value(params, "$XAXISLABEL", label2);
    LP_set_value(params, "$ERRBAR", errbars ? "true" : "false");
    if (height) LP_set_doublevalue(params, "$LOGOHEIGHT", height);
    if (width) LP_set_doublevalue(params, "$LOGOWIDTH", width);

    CL_init(params, asize < 20 ? "NA" : "AA", "SVG");
    if (motif1 != NULL) {
      CL_add_motif(params, motif1, alphabet, shift, label1);
    }
    if (motif2 != NULL) {
      CL_add_motif(params, motif2, alphabet, 0,     label2);
    }
    CL_exit(params, svgpath);
    LP_free(params);
  }
  myfree(svgpath);

}

/*.......................................................................*/

void CL_run_test(int argc, char* argv[]) {
  char* verbose = LP_get_option_default("-test", "1", argc,argv);
  T_set_verbose(verbose[0]=='-' ? 1 : atoi(verbose));

  E_test();
  F_test();
  STR_test();
  LP_test();
  LS_test();
  EPS_test();
  SVG_test();
  CL_test();

}

/*.......................................................................*/




/* -------------------------- Local Functions -------------------------- */

static void CL_usage() {
  printf(
   "USAGE: ceqlogo -i <filename> [options] \n"
   "  1. Example: Load all motifs within a MEME motif file\n"
   "     ceqlogo -i meme.motifs -o logo.eps -f EPS \n"
   "  2. Example: Load second motif from each of two files and shift the first one\n"
   "     ceqlogo -i2 meme1.motifs -s 2 -i2 meme2.motifs -o logo.eps -f EPS \n"
   "  3. Example: Run a self-test\n"
   "     ceqlogo -test \n"
   "\n"
   "Available options:\n"
   "  -i <input filename>        Loads all motifs within the file. \n"
   "  -i<n> <input filename>     Loads the n-th motif within the file. \n"
   "  -r                         Reverse complement motifs loaded by the next -i. \n"
   "  -s <shift>                 Shift for previously loaded motif (-i). \n"
   "  -b <bar bits>              Number of bits in bar (real # > 0). \n"
   "  -c <tic bits>              Number of bits between tic marks. \n"
   "  -e <error bar fraction>    Fraction of error bar to show (real # > 0). \n"
   "  -f <format>                Format of output (EPS, SVG); default: EPS\n"
   "                             Note: SVG not supported yet. \n"
   "  -h <logo height>           Height of output logo in cm (real # > 0). \n"
   "  -k <kind of data>          AA for amino acid, NA for nucleic acid. \n"
   "  -o <output file>           Output file path. Default is stdout. \n"
   "  -n <sample number>         Number of samples for previously loaded motif (-i).\n"
   "  -t <title label>           Label for title. \n"
   "  -d <fine print>            Descriptive fine print. \n"
   "  -w <logo width>            Width of output logo in cm (real # > 0). \n"
   "  -x <x-axis label>          Label for x-axis. \n"
   "  -y <y-axis label>          Label for y-axis. \n"
   "  -p <pseudocounts>          Pseudocounts for motifs. Default: 0. \n"
   " \n"
   "  -test [<verbose>]          Runs a self-test. All other options are ignored. \n"
   "                             The verbose level [0..3] is optional. \n"
   " \n"
   "Available toggles (all uppercase) \n"
   "  -S       Toggle small sample correction. \n"
   "  -B       Toggle bar ends  \n"
   "  -E       Toggle error bar  \n"
   "  -O       Toggle outlining of characters \n"
   "  -P       Toggle fineprint \n"
   "  -N       Toggle numbering of x-axis  \n"
   "  -X       Toggle boxing of characters  \n"
   "  -Y       Toggle y-axis \n"
  );
}

/*.......................................................................*/

static void CL_error(char* format, ...) {
  va_list argptr;

  va_start(argptr, format);
    vfprintf(stderr, format, argptr);
  va_end(argptr);
  fprintf(stderr, "\n\n");
  CL_usage();
  exit(-1);
}


/*.......................................................................*/
// TLB: Combine all modules used by ceqlogo in this area.
/*.......................................................................*/


/* ----------------------- Implementation -----------------------------

  Module name   : logoparam.c
  Module prefix : LP_

  Description: Describes the parameters of logo. Two types of parameters
    are distinguished. General parameters which concern the entire logo/diagram
    and line parameters which describe properties of a specific logo line.

    General parameters are organized as two arrays of strings
    that contain the parameters names (as they are used within the
    EPS template file) and the parameter values.

    Line parameters are organized in an array of structs.

------------------------------------------------------------------------ */


/* -------------------------- Local Prototypes ------------------------- */

static void LP_init(LP_PARAMS_T* params, int index, char* name, char *value);
/*
  Description   : Initializes a parameter. Note that the parameter name and
                  value are copied.
  Parameter     :
    params      : Parameter set.
    index       : Index of the parameter to init.
    name        : Parameter name.
    value       : Parameter value.
  Example       :
    LOP_PARAMS_T* params = LP_create();
      LP_init(params, 0, "TICBITS", "3.14");
    LP_free(params);
*/

/*.......................................................................*/

static void LP_init_all(LP_PARAMS_T* params);
/*
  Description   : Initializes the entire set of parameters with default values.
                The following parametes are initialized:

                $CREATOR : Description string of output creator        
                $CREATIONDATE : Date of creation  
                $LOGOTYPE : Type of logo: NN, NA       
                $LOGOHEIGHT : Hight in cm    
                $LOGOWIDTH  : Width in cm    
                $FONTSIZE   : Font size for labels    
                $TITLEFONTSIZE : Font size for title 
                $SMALLFONTSIZE  : Font size for fineprint
                $TOPMARGIN  : Top margin in cm    
                $BOTTOMMARGIN : Bottom margin in cm  
                $FINEPRINT : Fine print, any text you like     
                $YAXIS     : Toggles y-axis: "true" or "false"     
                $YAXISLABEL  : Label for y-axis   
                $XAXISLABEL  : Label for x-axis   
                $TITLE : Label for title         
		$SSC : Toggles small sample correction : "true" or "false"
                $ERRBAR : Toggles error bars:  "true" or "false"        
                $ERRORBARFRACTION : Length of error bars [0..1]
                $SHOWINGBOX : Toggles box around logo letters:  "true" or "false"    
                $BARBITS : Max. Bits on y-axis (e.g. "2.0" or "4.3.")        
                $TICBITS : Tics on y-axis       
                $COLORDEF : Colordefintion table (only EPS)      
                $COLORDICT : Color dictionary (only EPS)     
                $SHOWENDS : Toggles end descr. of x-axis:  "true" or "false"     
                $NUMBERING  : Toggles numbering of x-axis:  "true" or "false"    
                $OUTLINE : Toggles outlining of logo letters:  "true" or "false"        
                $CHARSPERLINE : Letter stacks per line. Don't touch this.   
                $BOUNDINGHEIGHT : Bounding box around logo. Don't touch this.
                $BOUNDINGWIDTH  : Bounding box around logo. Don't touch this.
                $LOGOLINEHEIGHT : Hight of a logo line. Don't touch this.
                $DATA  : Letter stack data. Don't touch this.         
                $FORMAT : EPS or SVG        
    
  Parameter     :
    params      : Parameter set.
  Return Values :
    params      : Parameter set will be modified.
  Example       :
    LP_init_all(params);
*/  
    
/*.......................................................................*/
    
static void LP_translate_option(LP_PARAMS_T* params, char* param,
                                char* option, int argc, char* argv[]);
/*  
  Description   : Translate a command line option into a parameter setting.
                  If the command line option is not set, the corresponding
                  parameter will not be set either.
  Parameter     :
    params      : Parameter set.
    param       : Parameter name in set.
    option      : Name of the option to translate.
    argc        : Number of arguments within the command line.
    argv        : An array that contains the arguments of the command line.
  Return Values :
    params      : Modifies the parameter set or prints an error message.
  Example       :
    LP_PARAMS_T* params = LP_create();
      ...
      LP_translate_option(params, "$TITLE", "-t", argc, argv);
    LP_free(params);

*/

/*.......................................................................*/

static void LP_translate_choices(LP_PARAMS_T* params, char* param,
                                 char* option, char *choices,
                                 int argc, char* argv[]);
/*
  Description   : Translate a command line option that has a set of
                  choices as associated value into a parameter setting.
  Parameter     :
    params      : Parameter set.
    param       : Parameter name in set.
    option      : Name of the option to translate.
    choices     : String with the possible choices.
    argc        : Number of arguments within the command line.
    argv        : An array that contains the arguments of the command line.
  Return Values :
    params      : Modifies the parameter set or prints an error message.
  Example       :
    LP_PARAMS_T* params = LP_create();
      ...
      LP_translate_cmdline(params, "$FORMAT", "-F",
                           "EPS PNG", argc, argv);
    LP_free(params);

*/

/*.......................................................................*/

static void LP_translate_numeric(LP_PARAMS_T* params, char* param,
                                char* option, double min, double max,
                                int argc, char* argv[]);
/*
  Description   : Translate a command line option with an associated
                  numeric value into a parameter setting.
  Parameter     :
    params      : Parameter set.
    param       : Parameter name in set.
    option      : Name of the option to translate.
    min         : Minimum value permitted.
    max         : Maximum value permitted.
    argc        : Number of arguments within the command line.
    argv        : An array that contains the arguments of the command line.
  Return Values :
    params      : Modifies the parameter set or prints an error message.
  Example       :
    LP_PARAMS_T* params = LP_create();
      ...
      LP_translate_cmdline(params, "TICBITS", "-B", 0.0, 10.0, argc, argv);
    LP_free(params);

*/

/*.......................................................................*/

static void LP_translate_boolean(LP_PARAMS_T* params, char* param,
                                char* option, int argc, char* argv[]);
/*
  Description   : Translate a command line option that acts as a boolean
                  switch into a parameter setting.
  Parameter     :
    params      : Parameter set.
    param       : Parameter name in set.
    option      : Name of the option to translate.
    argc        : Number of arguments within the command line.
    argv        : An array that contains the arguments of the command line.
  Return Values :
    params      : Modifies the parameter set or prints an error message.
  Example       :
    LP_PARAMS_T* params = LP_create();
      ...
      LP_translate_cmdline(params, "$OUTLINE", "-O", argc, argv);
    LP_free(params);

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


/* ------------------------- Global Functions -------------------------- */

LP_PARAMS_T* LP_create() 
{
  LP_PARAMS_T* params = (LP_PARAMS_T*)mymalloc(sizeof(LP_PARAMS_T));

  params->initialized   = 0;
  params->is_nucleotide = 1;
  params->is_reversecomp = 0;

  params->p_num         = 32;
  params->names         = (char**)mymalloc(sizeof(char*)*params->p_num);
  params->values        = (char**)mymalloc(sizeof(char*)*params->p_num);

  params->l_num         = 0;
  params->lines         = NULL;

  LP_init_all(params);

  return params;
}

/*.......................................................................*/

void LP_free(LP_PARAMS_T* params) {
  int i;

  /* free all parameter names and values */
  for(i=0; i<params->p_num; i++) {
    myfree(params->names[i]);
    myfree(params->values[i]);
  }
  params->p_num = 0;

  /* free all logo line specific parameters */
  for(i=0; i<params->l_num; i++) {
    LP_LINE_T line = params->lines[i];
    free_matrix(line.mat);
    myfree(line.alphabet);
    myfree(line.label);
  }
  params->l_num = 0;

  myfree(params->names);
  myfree(params->values);
  myfree(params->lines);
  myfree(params);
}

/*.......................................................................*/

void LP_add_line(LP_PARAMS_T* params, MATRIX_T* mat, int snum, int ltrim, 
        int rtrim, char* alphabet, int shift, char* label) {
  if(!params->initialized)
    E_system("Parameter set is not initialized!\n");
  if((strlen(alphabet) >= 20 && params->is_nucleotide) ||
     (strlen(alphabet) < 20 && !params->is_nucleotide))
    E_system("Alphabet and logo type (DNA,RNA,AA) do not match!\n");
  if(strlen(alphabet) > mat->num_cols)
    E_system("Alphabet size and matrix dimensions do not match!\n");

  /* add an empty line */
  params->l_num++;
  Resize(params->lines, params->l_num, LP_LINE_T);
  //params->lines = (LP_LINE_T*)myrealloc(params->lines,
  //                                      sizeof(LP_LINE_T)*params->l_num);

  /* fill the added line */
  LP_LINE_T* line = LP_get_lastline(params);
  line->mat       = duplicate_matrix(mat);
  line->snum      = snum;
  line->ltrim     = ltrim;
  line->rtrim     = rtrim;
  line->shift     = shift;
  line->alphabet  = STR_copy(alphabet);
  line->label     = label==NULL ? NULL : STR_copy(label);
}

/*.......................................................................*/

LP_LINE_T* LP_get_line(LP_PARAMS_T* params, int line) {
  if(line >= params->l_num)
    E_system("Invalid line number %d!\n", line);
  return(&(params->lines[line]));
}

/*.......................................................................*/

LP_LINE_T* LP_get_lastline(LP_PARAMS_T* params) {
  return(LP_get_line(params,params->l_num-1));
}

/*.......................................................................*/

size_t LP_max_valuelen(LP_PARAMS_T* params) {
  int maxlen = 0;
  int i;
  for(i=0; i<params->p_num; i++) {
    int len = strlen(params->values[i]);
    if(len > maxlen)
      maxlen = len;
  }
  return(maxlen);
}

/*.......................................................................*/

int LP_get_index(LP_PARAMS_T* params, char* name) {
  int i;
  for(i=0; i<params->p_num; i++)
    if(!strcmp(name, params->names[i]))
      return(i);
  return(-1);
}

/*.......................................................................*/

char* LP_get_value(LP_PARAMS_T* params, char* name) {
  int index = LP_get_index(params, name);
  return index<0 ? NULL : params->values[index];
}

/*.......................................................................*/

int LP_get_intvalue(LP_PARAMS_T* params, char* name) {
  char* value = LP_get_value(params,name);
  return value==NULL ? 0 : atoi(value);
}

/*.......................................................................*/

double LP_get_doublevalue(LP_PARAMS_T* params, char* name) {
  char* value = LP_get_value(params,name);
  return value==NULL ? 0.0 : atof(value);
}

/*.......................................................................*/

void LP_set_value(LP_PARAMS_T* params, char* name, char *value) {
  int index = LP_get_index(params, name);
  if(index < 0)
    E_system("Invalid parameter name: %s", name);
  else {
    myfree(params->values[index]);
    params->values[index] = STR_copy(value);
  }
}

/*.......................................................................*/

void LP_set_intvalue(LP_PARAMS_T* params, char* name, int value) {
  char strval[50];
  sprintf(strval, "%d", value);
  LP_set_value(params, name, strval);
}

/*.......................................................................*/

void LP_set_doublevalue(LP_PARAMS_T* params, char* name, double value) {
  char strval[50];
  sprintf(strval, "%g", value);
  LP_set_value(params, name, strval);
}

/*.......................................................................*/

void LP_correct_shift(LP_PARAMS_T* params) {
  int i;
  int min = 0;

  /* find minimum shift smaller than zero */
  for(i=0; i<params->l_num; i++)
    if(params->lines[i].shift < min)
      min = params->lines[i].shift;

  /* adjust all shift factors to become positive */
  for(i=0; i<params->l_num; i++)
    params->lines[i].shift -= min;
}

/*.......................................................................*/

void LP_translate_cmdline(LP_PARAMS_T* params, int argc, char* argv[]) {
  LP_translate_choices(params, "$LOGOTYPE",         "-k", "AA NA", argc, argv);
  LP_translate_numeric(params, "$BARBITS",          "-b", 0.1,  4.3,   argc, argv);
  LP_translate_numeric(params, "$ERRORBARFRACTION", "-e", 0.01, 1.0,   argc, argv);
  LP_translate_numeric(params, "$TICBITS",          "-c", 0.1,  4.3,   argc, argv);
  LP_translate_numeric(params, "$LOGOHEIGHT",       "-h", 0.01, 100.0, argc, argv);
  LP_translate_numeric(params, "$LOGOWIDTH",        "-w", 0.01, 100.0, argc, argv);
  LP_translate_option(params,  "$TITLE",            "-t", argc, argv);
  LP_translate_option(params,  "$FINEPRINT",        "-d", argc, argv);
  LP_translate_option(params,  "$XAXISLABEL",       "-x", argc, argv);
  LP_translate_option(params,  "$YAXISLABEL",       "-y", argc, argv);
  LP_translate_boolean(params, "$SSC",              "-S", argc, argv);
  LP_translate_boolean(params, "$YAXIS",            "-Y", argc, argv);
  LP_translate_boolean(params, "$SHOWENDS",         "-B", argc, argv);
  LP_translate_boolean(params, "$ERRBAR",           "-E", argc, argv);
  LP_translate_boolean(params, "$OUTLINE",          "-O", argc, argv);
  LP_translate_boolean(params, "$NUMBERING",        "-N", argc, argv);
  LP_translate_boolean(params, "$SHOWINGBOX",       "-X", argc, argv);
}

/*.......................................................................*/

char *LP_get_option(char *option, int argc, char* argv[]) {
  int i;
  for(i=0; i<argc; i++) {
    if(!strcmp(option, argv[i])) {
      if(i<argc-1 && argv[i][0] == '-' && argv[i+1][0] != '-')
        return argv[i+1];
      return argv[i];
    }
  }
  return NULL;
}

/*.......................................................................*/

char *LP_get_option_default(char *option, char* def, int argc, char* argv[]) {
  char* value = LP_get_option(option, argc, argv);
  return(value==NULL ? def : value);
}

/*.......................................................................*/

void LP_file_replace(LP_PARAMS_T* params, FILE *fin, FILE *fout, size_t len) {
  char *line = (char*)mymalloc(len*sizeof(char));

  if(fin && fout) {
    while(fgets( line, len, fin ) != NULL ) {   /* Read input lines */
      char name[100];
      char* value;
      char* start = NULL;
      char* end   = NULL;
      char* seg   = line;

      while( (start = strstr(seg, "{$")) != NULL ) {  /* search param. start */
        end = strchr(start,'}');                      /* search param. end */
        if(end == NULL) {
          E_error("Missing bracket after parameter name!\nline=%s\n",line);
          break;
        }
        *start = '\0';                             /* term. line segment */
        *end   = '\0';                             /* term. param. name */
        strcpy(name, start+1);                     /* get param. name */
        value = LP_get_value(params, name);
        if(value == NULL)
          E_system("Invalid parameter in template: %s\n", name+1);
        fputs(seg, fout);                          /* write line segment */
        fputs(value==NULL ? "" : value, fout);     /* write param. value */
        seg = end+1;                               /* jump to next line seg */
      }
      fputs(seg, fout);                            /* write remainder */
    }
  }

  myfree(line);
}

/*.......................................................................*/

char* LP_DNAcolordict() {
  char* colordict =
    "/fullColourDict <<\n"
    " (G)  orange\n"
    " (T)  green\n"
    " (C)  blue\n"
    " (A)  red\n"
    " (U)  green\n"
    ">> def\n"
    "/mutedColourDict <<\n"
    " (G)  lightorange\n"
    " (T)  lightgreen\n"
    " (C)  lightblue\n"
    " (A)  lightred\n"
    " (U)  lightgreen\n"
    ">> def\n"
    "/colorDict fullColourDict def\n";
  return(colordict);
}

/*.......................................................................*/

char* LP_AAcolordict() {
  char* colordict =
    "/fullColourDict <<\n"
    " (A)  blue\n"
    " (C)  blue\n"
    " (F)  blue\n"
    " (I)  blue\n"
    " (L)  blue\n"
    " (V)  blue\n"
    " (W)  blue\n"
    " (M)  blue\n"
    " (N)  green\n"
    " (Q)  green\n"
    " (S)  green\n"
    " (T)  green\n"
    " (D)  magenta\n"
    " (E)  magenta\n"
    " (K)  red\n"
    " (R)  red\n"
    " (H)  pink\n"
    " (G)  orange\n"
    " (P)  yellow\n"
    " (Y)  turquoise\n"
    ">> def\n"
    "/mutedColourDict <<\n"
    " (A)  lightblue\n"
    " (C)  lightblue\n"
    " (F)  lightblue\n"
    " (I)  lightblue\n"
    " (L)  lightblue\n"
    " (V)  lightblue\n"
    " (W)  lightblue\n"
    " (M)  lightblue\n"
    " (N)  lightgreen\n"
    " (Q)  lightgreen\n"
    " (S)  lightgreen\n"
    " (T)  lightgreen\n"
    " (D)  lightmagenta\n"
    " (E)  lightmagenta\n"
    " (K)  lightred\n"
    " (R)  lightred\n"
    " (H)  lightpink\n"
    " (G)  lightorange\n"
    " (P)  lightyellow\n"
    " (Y)  lightturquoise\n"
    ">> def\n"
    "/colorDict fullColourDict def\n";
  return(colordict);
}

/*.......................................................................*/

char* LP_colordef() {
  char* colordef =
    "/black [0 0 0] def\n"
    "/red [0.8 0 0] def\n"
    "/green [0 0.5 0] def\n"
    "/blue [0 0 0.8] def\n"
    "/yellow [1 1 0] def\n"
    "/purple [0.8 0 0.8] def\n"
    "/magenta [1.0 0 1.0] def\n"
    "/cyan [0 1.0 1.0] def\n"
    "/pink [1.0 0.8 0.8] def\n"
    "/turquoise [0.2 0.9 0.8] def\n"
    "/orange [1 0.7 0] def\n"
    "/lightred [0.8 0.56 0.56] def\n"
    "/lightgreen [0.35 0.5 0.35] def\n"
    "/lightblue [0.56 0.56 0.8] def\n"
    "/lightyellow [1 1 0.71] def\n"
    "/lightpurple [0.8 0.56 0.8] def\n"
    "/lightmagenta [1.0 0.7 1.0] def\n"
    "/lightcyan [0.7 1.0 1.0] def\n"
    "/lightpink [1.0 0.9 0.9] def\n"
    "/lightturquoise [0.81 0.9 0.89] def\n"
    "/lightorange [1 0.91 0.7] def\n";
  return(colordef);
}

/*.......................................................................*/




/* -------------------------- Local Functions -------------------------- */

static void LP_init(LP_PARAMS_T* params, int index, char* name, char *value) {
  assert(index < params->p_num);
  params->names[index]  = STR_copy(name);
  params->values[index] = STR_copy(value);
}

/*.......................................................................*/

static void LP_init_all(LP_PARAMS_T* params) {
  int n = 0;
  /* CAVEAT: Increase number of parameters in LP_create()
     when adding parameters to the following list and also
     adapt the assert statement at the end of the list. */
  LP_init(params, n++, "$CREATOR",             "Ceqlogo");
  LP_init(params, n++, "$CREATIONDATE",        STR_time("%d.%m.%y %H:%M:%S"));
  LP_init(params, n++, "$LOGOTYPE",            "NA");              /* NA, AA */
  LP_init(params, n++, "$LOGOHEIGHT",          "20");              /* in cm */
  LP_init(params, n++, "$LOGOWIDTH",           "20");              /* in cm */
  LP_init(params, n++, "$FONTSIZE",            "12");              /* in pts */
  LP_init(params, n++, "$TITLEFONTSIZE",       "12");              /* in pts */
  LP_init(params, n++, "$SMALLFONTSIZE",       "6");               /* in pts */
  LP_init(params, n++, "$TOPMARGIN",           "0.9");             /* in cm */
  LP_init(params, n++, "$BOTTOMMARGIN",        "0.9");             /* in cm */
  LP_init(params, n++, "$FINEPRINT",           STR_time("CEQLOGO %d.%m.%y %H:%M"));
  LP_init(params, n++, "$YAXIS",               "true");
  LP_init(params, n++, "$YAXISLABEL",          "bits");
  LP_init(params, n++, "$XAXISLABEL",          "Position");
  LP_init(params, n++, "$TITLE",               "Logo");
  LP_init(params, n++, "$SSC",                 "true");
  LP_init(params, n++, "$ERRBAR",              "true");             /* curr. only EPS */
  LP_init(params, n++, "$ERRORBARFRACTION",    "1.0");              /* curr. only EPS */
  LP_init(params, n++, "$SHOWINGBOX",          "false");
  LP_init(params, n++, "$BARBITS",             "4.3");              /* 2.0, 4.3 */
  LP_init(params, n++, "$TICBITS",             "1");                /* curr. only EPS */
  LP_init(params, n++, "$COLORDEF",            LP_colordef());      /* only EPS */
  LP_init(params, n++, "$COLORDICT",           LP_DNAcolordict());  /* only EPS */
  LP_init(params, n++, "$SHOWENDS",            "false");            /* curr. only EPS */
  LP_init(params, n++, "$NUMBERING",           "true");             /* x-axis numbering */
  LP_init(params, n++, "$OUTLINE",             "false");            /* letters outlined */
  LP_init(params, n++, "$CHARSPERLINE",        "");                 /* calculated */
  LP_init(params, n++, "$BOUNDINGHEIGHT",      "");                 /* calculated */
  LP_init(params, n++, "$BOUNDINGWIDTH",       "");                 /* calculated */
  LP_init(params, n++, "$LOGOLINEHEIGHT",      "");                 /* calculated */
  LP_init(params, n++, "$DATA",                "");                 /* calculated, stack data */
  LP_init(params, n++, "$FORMAT",              "EPS");              /* output format EPS,SVG */
  params->is_reversecomp = FALSE;                                   /* initilize reverse complement to false,
                                                                       this is only really used for the
                                                                       command line version */
  assert(n == params->p_num);  /* last index must match param. number! See LP_create() */
}

/*.......................................................................*/

static void LP_translate_option(LP_PARAMS_T* params, char* param,
                                char* option, int argc, char* argv[]) {
  char* value = LP_get_option(option, argc, argv);
  if(value != NULL)
    LP_set_value(params, param, value);
}

/*.......................................................................*/

static void LP_translate_choices(LP_PARAMS_T* params, char* param,
                                char* option, char *choices,
                                int argc, char* argv[]) {
  char* value = LP_get_option(option, argc, argv);
  if(value != NULL) {
    if(strstr(choices, value) == NULL)
       E_error("Value of option %s is invalid: %s\n", option, value);
    else
      LP_set_value(params, param, value);
  }
}

/*.......................................................................*/

static void LP_translate_numeric(LP_PARAMS_T* params, char* param,
                                char* option, double min, double max,
                                int argc, char* argv[]) {

  char* value = LP_get_option(option, argc, argv);
  if(value != NULL) {
    double number = atof(value);
    if(number < min  || number > max) {
      char buffer[100];
      sprintf(buffer, "Option %s must be between %g and %g: %s\n",
                       option, min, max, value);
      E_error(buffer);
    }
    else
      LP_set_value(params, param, value);
  }
}

/*.......................................................................*/

static void LP_translate_boolean(LP_PARAMS_T* params, char* param,
                                char* option, int argc, char* argv[]) {
  char* value = LP_get_option(option, argc, argv);
  LP_set_value(params, param, value==NULL ? "false" : "true");
}

/*.......................................................................*/







/* ------------------------------ Test --------------------------------- */

void LP_test() {
  T_start();
  {
    LP_PARAMS_T* params = LP_create();
    {
      T_int(-1, LP_get_index(params, "$FOO"));
      T_int(3, LP_get_index(params, "$LOGOHEIGHT"));
    }
    {
      T_is_null(LP_get_value(params, "$FOO"));
      T_string("20", LP_get_value(params, "$LOGOHEIGHT"));
    }
    {
      T_string(LP_DNAcolordict(), LP_get_value(params, "$COLORDICT"));
      T_string(LP_colordef(), LP_get_value(params, "$COLORDEF"));
    }
    {
      LP_set_value(params, "$LOGOHEIGHT", "90");
      T_string("90", LP_get_value(params, "$LOGOHEIGHT"));
    }
    {
      LP_set_intvalue(params, "$LOGOHEIGHT", 70);
      T_string("70", LP_get_value(params, "$LOGOHEIGHT"));
      T_int(70, LP_get_intvalue(params, "$LOGOHEIGHT"));
    }
    {
      LP_set_doublevalue(params, "$ERRORBARFRACTION", 0.07);
      T_string("0.07", LP_get_value(params, "$ERRORBARFRACTION"));
      T_double(0.07, LP_get_doublevalue(params, "$ERRORBARFRACTION"), 0.01);
    }
    LP_free(params);
  }
  {
    char *argv[] = {"-opt1", "val1", "-opt2", "-opt3", "-opt4", "val4"};
    int argc = 6;
    T_is_null(LP_get_option("-opt0", argc, argv));
    T_string("val1", LP_get_option("-opt1", argc, argv));
    T_string("-opt2", LP_get_option("-opt2", argc, argv));
    T_string("-opt3", LP_get_option("-opt3", argc, argv));
    T_string("val4", LP_get_option("-opt4", argc, argv));
  }
  {
    LP_PARAMS_T* params = LP_create();
      char *argv[] = {"-b", "3.14", "-e", "0.2", "-c", "2",
                      "-h", "79", "-w", "95",
                      "-t", "mytitle", "-x", "xlabel", "-y", "ylabel",
                      "-Y", "-B", "-E", "-O", "-N", "-X"};
      int argc =22;
      LP_translate_cmdline(params, argc, argv);
      T_string("3.14", LP_get_value(params, "$BARBITS"));
      T_string("0.2", LP_get_value(params, "$ERRORBARFRACTION"));
      T_string("2", LP_get_value(params, "$TICBITS"));
      T_string("79", LP_get_value(params, "$LOGOHEIGHT"));
      T_string("95", LP_get_value(params, "$LOGOWIDTH"));
      T_string("mytitle", LP_get_value(params, "$TITLE"));
      T_string("xlabel", LP_get_value(params, "$XAXISLABEL"));
      T_string("ylabel", LP_get_value(params, "$YAXISLABEL"));
      T_string("true", LP_get_value(params, "$YAXIS"));
      T_string("true", LP_get_value(params, "$SHOWENDS"));
      T_string("true", LP_get_value(params, "$ERRBAR"));
      T_string("true", LP_get_value(params, "$OUTLINE"));
      T_string("true", LP_get_value(params, "$NUMBERING"));
      T_string("true", LP_get_value(params, "$SHOWINGBOX"));
    LP_free(params);
  }
  /*{
    FILE *fin  = F_open("etc/template.eps", "r");
    FILE *fout = F_open("test.test", "w");
    LP_PARAMS_T* params = LP_create();
      LP_file_replace(params, fin, fout, 1000);
    LP_free(params);
    F_close(fin);
    F_close(fout);
  }*/
  T_end();
}

/* ------------------------------ Main --------------------------------- */

#if 0
int main(int argc, char* argv[])
{
  LP_test();
  return 0;
}
#endif

/* -------------------------- Local Variables -------------------------- */
static size_t F_counter = 0;



/* ------------------------- Global Functions -------------------------- */

FILE* _F_open(const char *name, const char *type, int line, const char *file) {
  FILE *fp = fopen(name, type);
  if(fp == NULL) {
    fprintf(stderr,"Cannot open file '%s'\n(line: %d, file: %s)\n", name, line, file);
    return(NULL);
  }
  F_counter++;
  return(fp);
}

/*.......................................................................*/

void _F_close(FILE *fp, int line, const char *file) {
  if(fp != NULL) {
    fclose(fp);
  }
  else
    fprintf(stderr,"File with NULL handle closed!\n(line: %d, file: %s)\n",line, file);
  F_counter--;
}

/*.......................................................................*/

int F_check() {
  if(F_counter != 0) {
    fprintf(stderr, "\nNo all files or too many were closed!");
  }
  return(1);
}

/*.......................................................................*/

int F_exists(const char *name) {
  FILE *fp = fopen(name, "r");
  if(fp == NULL)
    return(0);
  fclose(fp);
  return(1);
}

/*.......................................................................*/

void F_rmdir(char* dir, char* mode) {
  int ret;
  char cmdline[500];
  sprintf(cmdline, "rm %s %s",mode,dir);
  ret = system(cmdline);

  if (!(WIFEXITED(ret) && WEXITSTATUS(ret) == 0)) {
    fprintf(stderr, "Warning: failed to remove directory (%s).\n", dir);
  }
}

/*.......................................................................*/

void F_mkdir(char* dir, char* mode) {
  int ret;
  char cmdline[500];
  sprintf(cmdline, "mkdir %s %s",mode,dir);
  ret = system(cmdline);

  if (!(WIFEXITED(ret) && WEXITSTATUS(ret) == 0)) {
    fprintf(stderr, "Warning: failed to make directory (%s).\n", dir);
  }
}

/*.......................................................................*/



/* -------------------------- Local Functions -------------------------- */


/*
  Description   :
  Parameter     :
  Global Var.   :
  Return Values :
    function    :
  Example       :
*/




void F_test() {
  T_start();
  {
    FILE *fp = F_open("test.test", "w");
    F_close(fp);
    T_is_true(F_check());
  }
  {
    FILE *fp = F_open("test.test", "r");
    F_close(fp);
    T_is_true(F_check());
  }
  T_end();
}



/* ------------------------------ Main --------------------------------- */

#if 0
  int main(int argc, char* argv[])
  {
     F_test();
    return 0;
  }
#endif


/* ------------------------- Global Functions -------------------------- */

char* STR_time(char* format) {
  static char timestr[100];
  time_t now = time(NULL);
  strftime(timestr,100,format, localtime(&now));
  return(timestr);
}

/*.......................................................................*/

char* STR_copy(char* str) {
  char* copy = (char*)mymalloc(sizeof(char)*(strlen(str)+1));
  strcpy(copy, str);
  return(copy);
}

/*.......................................................................*/

void STR_replace(char buffer[], int n, const char *searchstr, const char *replacestr) {
  char *found = NULL;
  char       *p     = buffer;
  char       *end   = buffer + strlen(buffer) + 1;
  int         slen  = strlen(searchstr);
  int         rlen  = strlen(replacestr);

  while( (found = strstr(p, searchstr)) != NULL ) {
    assert((found+rlen)-buffer <= n);
    memmove(found+rlen, found+slen, end-(found+slen));
    memcpy(found, replacestr, rlen);
    p = found + rlen;
    end += (rlen-slen);
  }
}

/*.......................................................................*/

int STR_replace_diff(const char *str, const char *searchstr, const char *replacestr) {
  char *found;
  int  sum  = 0;
  int  diff = strlen(replacestr)-strlen(searchstr);

  while( *str && (found = strstr(str, searchstr)) != NULL ) {
    sum += diff;
    str = found + strlen(searchstr);
  }

  return(sum);
}

/*.......................................................................*/

char *STR_append(char **str, char* append) {
  char *newstr = (char *)mymalloc(sizeof(char)*(strlen(*str)+strlen(append)+1));
  if(newstr) {
    strcpy(newstr, *str);
    strcat(newstr, append);
    myfree(*str);
    *str = newstr;
  }
  return newstr;
}

/*.......................................................................*/



/* -------------------------- Local Functions -------------------------- */


/*
  Description   :
  Parameter     :
  Global Var.   :
  Return Values :
    function    :
  Example       :
*/

/*.......................................................................*/




/* ------------------------------ Test --------------------------------- */

void STR_test() {
  T_start();
  { /* Test STR_replace_diff */
    T_int(-2, STR_replace_diff("aaabbaaabbaa", "bb", "b"));
    T_int(-3, STR_replace_diff("aaabbaaabbaabb", "bb", "b"));
    T_int(-4, STR_replace_diff("aaabbbaabbb", "bbb", "b"));
  }

  { /* Test STR_replace */
    char buffer[15];
    STR_replace(strcpy(buffer, "aaabbaaabba"), 15, "bb", "B");
    T_string("aaaBaaaBa", buffer);

    STR_replace(strcpy(buffer, "aabbbbabb"), 15, "bb", "BBB");
    T_string("aaBBBBBBaBBB", buffer);

    STR_replace(strcpy(buffer, "aaaa"), 15, "bb", "BBB");
    T_string("aaaa", buffer);
  }
  { /* Test STR_copy */
    char* str = STR_copy("test");
    T_string("test", str);
    myfree(str);
  }
  { /* Test STR_append */
    char *str = (char*)mymalloc(20);
    strcpy(str, "blue");
    STR_append(&str, " and");
    T_string("blue and", str);
    STR_append(&str, " red");
    T_string("blue and red", str);
    myfree(str);

  }
  T_end();
}


/* ------------------------------ Main --------------------------------- */

#if 0
  int main(int argc, char* argv[])
  {
    STR_test();
    return 0;
  }
#endif


/* -------------------------- Local Variables -------------------------- */
/* Number of failed tests. */
static int T_failed;

/* Number of tests performed */
static int T_counter;

/* Level of verbosity 0..3 */
static int T_verbose = 1;


/* -------------------------- Local Prototypes ------------------------- */

static void T_count(char* expression);
/*
  Description   : Counts the number of performed tests.
    expression  : Test expression.
  Parameter     :
  Example       :
    T_count(expression);
*/

static void T_print_location(int line, char *file);
/*
  Description   : Prints the test location information.
  Parameter     :
    line        : Line number.
    file        : File name.
  Example       :
    T_print_location(__LINE__, __FILE__);
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


/* ------------------------- Global Functions -------------------------- */

void T_set_verbose(int verbose_level) {
  T_verbose = verbose_level;
}

/*.......................................................................*/

void _T_start(char *filename) {
  if(T_verbose > 0)
    fprintf(stdout, "test: %s\n", filename);
  T_failed = 0;
  T_counter = 0;
}

/*.......................................................................*/

void T_end() {
  if(!T_failed) {
    if(T_verbose == 1) {
      fprintf(stdout, " => OK\n");
    } else if (T_verbose > 1) {
      fprintf(stdout, " %d tests\n => OK\n", T_counter);
    }
  }
}

/*.......................................................................*/

void _T_is_null(char* expression, void* actual, int line, char *file) {
  T_count(expression);
  if(actual != NULL) {
    fprintf(stderr, "\n *** failed: %s\n", expression);
    fprintf(stderr, "   expected: NULL\n");
    fprintf(stderr, "   actual  : %p\n", actual);
    T_print_location(line, file);
  }
}

/*.......................................................................*/

void _T_string(char* expression, char* expected, char* actual, int line, char *file) {
  T_count(expression);
  if(actual==NULL || strcmp(expected, actual)) {
    fprintf(stderr, "\n *** failed: %s\n", expression);
    fprintf(stderr, "   expected: %s\n", expected);
    fprintf(stderr, "   actual  : %s\n", actual);
    T_print_location(line, file);
  }
}

/*.......................................................................*/

void _T_boolean(char* expression, int expected, int actual, int line, char *file) {
  T_count(expression);
  if(!(expected && actual)) {
    fprintf(stderr, "\n *** failed: %s\n", expression);
    fprintf(stderr, "   expected: %s\n", expected ? "true" : "false");
    fprintf(stderr, "   actual  : %s\n", actual ? "true" : "false");
    T_print_location(line, file);
  }
}

/*.......................................................................*/

void _T_int(char* expression, int expected, int actual, int line, char *file) {
  T_count(expression);
  if(expected != actual) {
    fprintf(stderr, "\n *** failed: %s\n", expression);
    fprintf(stderr, "   expected: %d\n", expected );
    fprintf(stderr, "   actual  : %d\n", actual);
    T_print_location(line, file);
  }
}

/*.......................................................................*/

void _T_char(char* expression, char expected, char actual, int line, char *file) {
  T_count(expression);
  if(expected != actual) {
    fprintf(stderr, "\n *** failed: %s\n", expression);
    fprintf(stderr, "   expected: %c\n", expected );
    fprintf(stderr, "   actual  : %c\n", actual);
    T_print_location(line, file);
  }
}

/*.......................................................................*/

void _T_double(char* expression, double expected, double actual, double eps, int line, char *file) {
  T_count(expression);
  if(fabs(expected-actual) > eps) {
    fprintf(stderr, "\n *** failed: %s\n", expression);
    fprintf(stderr, "   expected: %lf\n", expected );
    fprintf(stderr, "   actual  : %lf\n", actual);
    T_print_location(line, file);
  }
}


/*.......................................................................*/





/* -------------------------- Local Functions -------------------------- */

static void T_count(char *expression) {
  T_counter++;
  if(T_verbose == 2) {
    fprintf(stdout, ".");
  } else if(T_verbose == 3) {
    fprintf(stdout, "    %3d: %s\n", T_counter, expression);
  }
}

/*.......................................................................*/

static void T_print_location(int line, char *file) {
  fprintf(stderr, "   file    : %s\n", file);
  fprintf(stderr, "   line    : %d\n", line);
  fprintf(stderr, "   number  : %d\n", T_counter);
  T_failed++;
}

/*.......................................................................*/




/* ------------------------------ Main --------------------------------- */

#ifdef T_MAIN
  int main(int argc, char* argv[])
  {
    T_start();
      T_is_null(NULL);
      T_is_null("Not null");

      T_string("blah", "blah");
      T_string("blah", NULL);
      T_string("blah", "bluh");

      T_boolean(1, 1);
      T_boolean(1, 0);

      T_int(1, 1);
      T_int(1, 2);

      T_char('c', 'c');
      T_char('a', 'c');

      T_double(1.00, 1.00, 0.01);
      T_double(1.00, 1.01, 0.1);
      T_double(1.00, 1.01, 0.01);
    T_end();
    return 0;
  }
#endif

/* ----------------------- Implementation -----------------------------

  Module name   : lstack.c
  Module prefix : LS_

  Description: Implements a letter stack to generate a logo.
               The stack is calculated on base of a Position Weight Matrix
               (PWM) that contains the probabilites for letters
               (nucleotides, amino acids) within a set of aligned
               sequences.
               The columns of the matrix represent the alphabet
               and the rows the positions within the alignment.
               The letters of the alphabet are assumed to label
               the corresponding columns of the matrix.

               The letter stacks are calculated by the following equations found in
               Schneider and Stephens paper "Sequence Logos: A New Way to Display
               Consensus Sequences":

                 height = f(b,l) * R(l)                            (1)

               where f(b,l) is the frequency of base or amino acid "b" at position
               "l". R(l) is amount of information present at position "l" and can
               be quantified as follows:

                 R(l) for amino acids   = log(20) - (H(l) + e(n))    (2a)
                 R(l) for nucleic acids =    2    - (H(l) + e(n))    (2b)

               where log is taken base 2, H(l) is the uncertainty at position "l",
               and e(n) is the error correction factor for small "n". H(l) is
               computed as follows:

                   H(l) = - (Sum f(b,l) * log[ f(b,l) ])             (3)

               where again, log is taken base 2. f(b,l) is the frequency of base
               "b" at position "l". The sum is taken over all amino acids or
               bases, depending on which the data is.

               Currently, logo.pm uses an approximation for e(n), given by:

                   e(n) = (s-1) / (2 * ln2 * n)                      (4)

               Where s is 4 for nucleotides, 20 for amino acids ; n is the number
               of sequences in the alignment. e(n) also  gives the height of error
               bars.
               The code and the comment above are based on the Perl code from logo.pm
               which is part of weblogo: http://weblogo.berkeley.edu

*/

/* -------------------------- Local Prototypes ------------------------- */

static double LS_calc_Hmax(size_t asize);
/*
  Description   : Calculates the max. entropy for the given alphabet size.
  Parameter     :
    asize       : Size of the alphabet.
  Return Values :
    function    : Returns log_2();
  Example       :
    double Hmax = LS_calc_Hmax(4);
*/

/*.......................................................................*/

static double LS_calc_H(int row, MATRIX_T* mat, size_t asize);
/*
  Description   : Calculates the entropy for a matrix row (= position of
                  the logo).
  Parameter     :
    row         : Matrix row. First row is zero.
    mat         : PWM.
    asize       : Size of the alphabet. The number of matrix columns must be
                  greater than or equal to asize.
  Return Values :
    function    : Returns the entropy for the specified matrix row.
  Example       :
    double H = LS_calc_H(0,mat,4);
*/

/*.......................................................................*/

static double LS_calc_R(int row, size_t n, MATRIX_T* mat, size_t asize, 
  BOOLEAN_T ssc);
/*
  Description   : Calculate the information content for a matrix row,
                  taking the error correction factor into account.
  Parameter     :
    row         : Matrix row. First row is zero.
    n           : Number of samples. Zero is valid - in this case the
                  correction factor is zero as well.
    mat         : PWM.
    asize       : Size of the alphabet. The number of matrix columns must be
                  greater than or equal to asize.
    ssc		: Use small sample correction
  Return Values :
    function    : Return the information content for a matrix row
                  (= position within the logo).
  Example       :
    double R = LS_calc_R(0,0,mat,4,TRUE);
*/

/*.......................................................................*/

static double LS_calc_height(int row, int col, size_t n, MATRIX_T* mat, 
  size_t asize, BOOLEAN_T ssc);
/*
  Description   : Calculates the height of a letter stack within the logo
                  for a specific logo position (=row of the PWM).
  Parameter     :
    row         : Matrix row. First row is zero.
    n           : Number of samples. Zero is valid - in this case the
                  correction factor is zero as well.
    mat         : PWM.
    asize       : Size of the alphabet. The number of matrix columns must be
                  greater than or equal to asize.
    ssc		: use small sample correction
  Return Values :
    function    : Returns the height of the letter stack for a specific
                  logo position.
  Example       :
    double height = LS_calc_height(0,0,mat,4,TRUE);
*/

/*.......................................................................*/


static int LS_comp(const void *lstack1, const void *lstack2);
/*
  Description   : A comparator function to sort the letter stack according
                  to letter size.
  Parameter     :
    lstack1     : Reference to a first stack element.
    lstack2     : Reference to a second stack element.
  Return Values :
    function    : +1 if height of first stack element is greater than
                  height of second stack element, and -1 otherwise.
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


/* ------------------------- Global Functions -------------------------- */

LSTACK_T* LS_create(size_t size) {
  return( (LSTACK_T*)mymalloc(sizeof(LSTACK_T)*size) );
}

/*.......................................................................*/

void LS_free(LSTACK_T* lstack) {
   myfree(lstack);
}

/*.......................................................................*/

void LS_calc(int row, size_t n, MATRIX_T* mat, char* alphabet, 
  LSTACK_T* lstack, BOOLEAN_T ssc) 
{
  int    col;
  size_t asize = strlen(alphabet);

  if(asize > mat->num_cols)
    E_system("Too few matrix columns! Columns=%d, Alphabet=%s\n",mat->num_cols, alphabet);

  /* Over all letters of the alphabet */
  for(col=0; col<asize; col++) {
    lstack[col].letter = alphabet[col];
    lstack[col].height = LS_calc_height(row,col,n,mat,asize,ssc);
  }

  /* Sort according to height */
  qsort(lstack, asize, sizeof(LSTACK_T), LS_comp);
} // LS_calc

/*.......................................................................*/

double LS_calc_e(size_t n, size_t asize) {
  if(n == 0)
    return(0.0);
  return((asize-1) / (2.0*log(2.0)*n));
}

/*.......................................................................*/




/* -------------------------- Local Functions -------------------------- */

static double LS_calc_Hmax(size_t asize) {
  return( log(asize)/log(2.0) );
}

/*.......................................................................*/

static double LS_calc_H(int row, MATRIX_T* mat, size_t asize) {
  int    col;
  double H = 0.0;
  for(col=0; col<asize; col++) {
    double p = get_matrix_cell(row,col,mat);
    if(p > 0.0) H -= p * log(p)/log(2.0);
  }
  return(H);
}

/*.......................................................................*/

static double LS_calc_R(int row, size_t n, MATRIX_T* mat, size_t asize,
  BOOLEAN_T ssc) 
{
  return(LS_calc_Hmax(asize) - 
    (LS_calc_H(row, mat, asize) + (ssc ? LS_calc_e(n, asize) : 0 ))
  );
}

/*.......................................................................*/

static double LS_calc_height(int row, int col, size_t n, MATRIX_T* mat, 
  size_t asize, BOOLEAN_T ssc) 
{
  double p = get_matrix_cell(row, col, mat);
  return(p * LS_calc_R(row, n, mat, asize, ssc) );
}

/*.......................................................................*/

static int LS_comp(const void *lstack1, const void *lstack2) {
  LSTACK_T* ls1 = (LSTACK_T*)lstack1;
   LSTACK_T* ls2 = (LSTACK_T*)lstack2;
  return( (ls1->height > ls2->height) ? 1 : -1 );
}

/*.......................................................................*/





/* ------------------------------ Test --------------------------------- */

void LS_test() {
  T_start();
  char* alphabet = "ACTG";
  {
    MATRIX_T* mat = allocate_matrix(2, strlen(alphabet));
    set_matrix_cell(0,0, 0.1, mat);
    set_matrix_cell(0,1, 0.7, mat);
    set_matrix_cell(0,2, 0.1, mat);
    set_matrix_cell(0,3, 0.1, mat);
    {
      T_double(2.00, LS_calc_Hmax(strlen(alphabet)), 0.01);
    }
    {
      T_double(1.35, LS_calc_H(0,mat,strlen(alphabet)), 0.01);
    }
    {
      T_double(0.00, LS_calc_e(0,strlen(alphabet)), 0.01);
      T_double(0.21, LS_calc_e(10,strlen(alphabet)), 0.01);
    }
    {
      T_double(0.64, LS_calc_R(0,0,mat,strlen(alphabet), TRUE), 0.01);
    }
    {
      T_double(0.064, LS_calc_height(0,0,0,mat,strlen(alphabet), TRUE), 0.001);
    }
    free_matrix(mat);
  }
  T_end();

}


/* ------------------------------ Main --------------------------------- */

#if 0
  int main(int argc, char* argv[])
  {
    LS_test();
    return 0;
  }
#endif

/* ------------------------------ Definitions -------------------------- */





/* -------------------------- Local Prototypes ------------------------- */

static char *SVG_DNAcolor(char ch);
/*
  Description   : Getter for a DNA symbol color.
  Parameter     : DNA symbol.
  Return Values :
    function    : Returns a string that describes the color of the given
                  DNA symbol.
  Example       :
    SVG_DNAcolor('A') -> "green"
*/

/*.......................................................................*/

static char *SVG_AAcolor(char ch);
/*
  Description   : Getter for a Amino Acid symbol color.
  Parameter     : Amino Acid symbol.
  Return Values :
    function    : Returns a string that describes the color of the given
                  Amino Acid symbol.
  Example       :
    SVG_AAcolor('V') -> "green"
*/

/*.......................................................................*/

static void SVG_append(char**str, char* format, ...);
/*
  Description   : Appends a formatted string to a given string.
                  ATTENTION: The formatted string must not be longer than
                  1000 characters!
  Parameter     :
    str         : Address of the string to extend. Must be dynamically
                  allocated.
    format      : Format string, see printf()
    ...         : Variable argument list.
  Return Values :
    str         : Updates the pointer to the extended string.
  Example       :
    char *str = (char*)mymalloc(20);
    strcpy(str, "var");
    SVG_append(&str, "=%d", 5);   ->  "var=5"
    myfree(str);
*/

/*
  Description   :
  Parameter     :
  Global Var.   :
  Return Values :
    function    :
  Example       :
*/


/* ------------------------- Global Functions -------------------------- */

void SVG_init(LP_PARAMS_T* params, char* type) {
  int    is_nucleotide = !strcmp(type, "NA");
  double height        = atof(LP_get_value(params, "$LOGOHEIGHT"));
  double width         = atof(LP_get_value(params, "$LOGOWIDTH"));
  double c             = 72/2.54;   /* conversion cm to pts */

  if(params->initialized)
    E_system("Parameter set is already initialized!\n");
  params->initialized   = 1;
  params->is_nucleotide = is_nucleotide;


  LP_set_value(params,    "$DATA", "");
  LP_set_value(params,    "$LOGOTYPE", type);
  LP_set_value(params,    "$BARBITS", is_nucleotide ? "2.0" : "4.3");
  LP_set_intvalue(params, "$BOUNDINGHEIGHT", (int)(height*c));
  LP_set_intvalue(params, "$BOUNDINGWIDTH", (int)(width*c));
}

/*.......................................................................*/

void SVG_add(LP_PARAMS_T* params, MATRIX_T* mat, int snum, 
    int ltrim, int rtrim, char* alphabet, int shift, char* label) {
  LP_add_line(params, mat, snum, ltrim, rtrim, alphabet, shift, label);
}

/*.......................................................................*/

static void SVG_write_line(char**data, LP_PARAMS_T* params, int lnr, double ow, double oh, double w, double h) {
  int        row;
  int        col;
  LP_LINE_T* line     = LP_get_line(params, lnr);
  MATRIX_T*  mat      = line->mat;
  int        shift    = line->shift;
  char*      alphabet = line->alphabet;
  int        asize    = strlen(alphabet);
  int        fs       = LP_get_doublevalue(params, "$FONTSIZE");
  int        num      = !strcmp("true", LP_get_value(params, "$NUMBERING"));
  int        otl      = !strcmp("true", LP_get_value(params, "$OUTLINE"));
  int        box      = !strcmp("true", LP_get_value(params, "$SHOWINGBOX"));
  int        tlen     = fs/2;              /* tic length */
  double     bits     = atof(LP_get_value(params, "$BARBITS"));
  double     ch       = h/bits;            /* cell heigth */
  double     cw       = w/mat->num_rows;   /* cell width */
  double     d        = 1;           /* Distance between letters in pixel */
  double     cf       = 73.0;        /* Letter size correction factor, depends on font */
  BOOLEAN_T  ssc      = !strcmp("true", LP_get_value(params, "$SSC"));

  if(!strcmp("true", LP_get_value(params, "$YAXIS"))) {
    SVG_append(data,
            "<path class=\"axes\" d=\"M%f,%f v%f\"/>\n", ow, oh, h);
    SVG_append(data,
            "<path class=\"axes\" d=\"M%f,%f h%d M%f,%f h%d\"/>\n\n",
             ow, oh, -tlen, ow, oh+h, -tlen);
    SVG_append(data,
            "<text class=\"ylabel\" transform=\"translate(%f,%f) rotate(-90)\" >%s</text>\n",
            ow-2*fs, oh+h/2, LP_get_value(params, "$YAXISLABEL"));
    SVG_append(data,
            "<text class=\"yticlabel\" transform=\"translate(%f,%f)\" >%1.1g</text>\n",
            ow-fs/2, oh+fs/2, bits);
    SVG_append(data,
            "<text class=\"yticlabel\" transform=\"translate(%f,%f)\" >%.1g</text>\n\n",
            ow-tlen, oh+h+fs/2, 0.0);
  }

  LSTACK_T* lstack = LS_create(mat->num_cols);
  {
    /* Over all positions */
    for(row=0; row<mat->num_rows; row++) {
       LS_calc(row,0,mat,alphabet,lstack,ssc);
       double sh = 0;

      /* over all letters */
      for(col=0; col<asize; col++) {
        if(lstack[col].height > 0) {
          char   letter = lstack[col].letter;
          double lh     = ch*lstack[col].height;
          char*  color  = asize==4 ? SVG_DNAcolor(letter) : SVG_AAcolor(letter);
          char*  fill   = otl ? "none" : color;
          char*  stroke = otl ? color : "none";
          SVG_append(data,
              "<text class=\"letter\" transform=\"translate(%f,%f)"
              " scale(%f,%f)\" style=\"fill:%s;stroke:%s\" >%c</text>\n",
              ow+cw/2+cw*(shift+row), oh+h-sh, cw/cf, lh/cf, fill, stroke, letter);
          if(box)  /* boxing is on */
             SVG_append(data,
                 "<rect class=\"boxing\" x=\"%f\" y=\"%f\" width=\"%f\" height=\"%f\"/>\n",
                 ow+cw*(shift+row), oh+h-sh-lh, cw, lh);
          sh += lh+d;
        }
      }

      if(num) {  /* Numbering of x-axis enables */
        SVG_append(data,
            "<text class=\"xticlabel\" transform=\"translate(%f,%f) rotate(-90)\" >%d</text>\n\n",
            ow+cw*(shift+row)+cw/2+fs/2, oh+h+3, row+1);
      }
    }
  }
  LS_free(lstack);
}

/*.......................................................................*/

void SVG_exit(LP_PARAMS_T* params, FILE *fout) {
  int    lnr;
  char*  data     = LP_get_value(params, "$DATA");
  double fs       = LP_get_doublevalue(params, "$FONTSIZE");
  double tfs      = LP_get_doublevalue(params, "$TITLEFONTSIZE");
  double h        = LP_get_doublevalue(params, "$LOGOHEIGHT")*20;
  double w        = LP_get_doublevalue(params, "$LOGOWIDTH")*15;
  double tm       = LP_get_doublevalue(params, "$TOPMARGIN")*10;
  double bm       = LP_get_doublevalue(params, "$BOTTOMMARGIN")*10;
  double ow       = fs*3;              /* offset width */
  double oh       = 1.1*tfs+tm;        /* offset height */

  if(!params->l_num)
    E_system("No lines added to the logo!\n");

  SVG_append(&data,
      "<text class=\"title\" transform=\"translate(%f,%f)\" >%s</text>\n\n",
       ow+w/2, oh, LP_get_value(params, "$TITLE"));

  SVG_append(&data,
      "<text class=\"xlabel\" transform=\"translate(%f,%f)\" >%s</text>\n\n",
      ow+w/2, h-oh-2*fs+bm, LP_get_value(params, "$XAXISLABEL"));

  SVG_append(&data,
       "<text class=\"fineprint\" transform=\"translate(%f,%f)\" >%s</text>\n\n",
       ow+w, h-1.1*tfs, LP_get_value(params, "$FINEPRINT"));


  LP_correct_shift(params);
  double lh = (h-tm-bm)/params->l_num;
  for(lnr=0; lnr<params->l_num; lnr++) {         /* Over all logo lines */
    SVG_write_line(&data, params, lnr, ow, oh+lnr*lh, w, lh-fs*2);
  }

  params->values[LP_get_index(params,"$DATA")] = data;

  FILE *fin  = F_open(TEMPLATE_SVG, "r");
    LP_file_replace(params, fin, fout, 1000);
  F_close(fin);

  params->initialized = 0;
}





/* -------------------------- Local Functions -------------------------- */

static char *SVG_DNAcolor(char ch) {
  switch(toupper(ch)) {
    case 'A': return("green");
    case 'C': return("blue");
    case 'T': return("red");
    case 'G': return("orange");
    case 'U': return("red");
  }
  return("black");
}

/*.......................................................................*/

static char *SVG_AAcolor(char ch) {
  switch(toupper(ch)) {
    case 'A':
    case 'F':
    case 'I':
    case 'L':
    case 'M':
    case 'V':
    case 'W':
    case 'C': return("blue");
    case 'S':
    case 'T':
    case 'N':
    case 'Q': return("green");
    case 'D':
    case 'E': return("magenta");
    case 'K':
    case 'R': return("red");
    case 'H': return("pink");
    case 'G': return("orange");
    case 'P': return("yellow");
    case 'Y': return("turquoise");
  }
  return("black");
}

/*.......................................................................*/

static void SVG_append(char**str, char* format, ...) {
  static char buffer[1000];
  va_list argptr;

  va_start(argptr, format);
    vsprintf(buffer, format, argptr);
  va_end(argptr);

  STR_append(str, buffer);
}

/*.......................................................................*/



/* ------------------------------ Test --------------------------------- */

void SVG_test() {
  T_start();
  T_end();
}


/* ------------------------------ Main --------------------------------- */

#if 0
  int main(int argc, char* argv[])
  {
    SVG_test();
    return 0;
  }
#endif

/* ------------------------- Global Functions -------------------------- */

void EPS_init(LP_PARAMS_T* params, char* type) {
  int    is_nucleotide = !strcmp(type, "NA");
  char*  color         = is_nucleotide ? LP_DNAcolordict() : LP_AAcolordict();
  double height        = atof(LP_get_value(params, "$LOGOHEIGHT"));
  double width         = atof(LP_get_value(params, "$LOGOWIDTH"));
  double c             = 72/2.54;   /* conversion cm to pts */

  if(params->initialized)
    E_system("Parameter set is already initialized!\n");
  params->initialized   = 1;
  params->is_nucleotide = is_nucleotide;

  LP_set_value(params,    "$DATA", "");
  LP_set_value(params,    "$LOGOTYPE", type);
  LP_set_value(params,    "$COLORDICT", color);
  LP_set_value(params,    "$BARBITS", is_nucleotide ? "2.0" : "4.3");
  LP_set_intvalue(params, "$BOUNDINGHEIGHT", (int)(height*c));
  LP_set_intvalue(params, "$BOUNDINGWIDTH", (int)(width*c));
}

/*.......................................................................*/

void EPS_add(LP_PARAMS_T* params, MATRIX_T* mat, int snum, int ltrim, 
             int rtrim, char* alphabet, int shift, char* label) {
  LP_add_line(params, mat, snum, ltrim, rtrim, alphabet, shift, label);
}

/*.......................................................................*/

void EPS_add_motifs(LP_PARAMS_T* params, MOTIF_T* motifs, int num_motifs, char* alphabet) {
  int i;
  MOTIF_T *motif;
  for(i = 0; i < num_motifs; i++) {
    motif = motifs+i;
    EPS_add(params, motifs->freqs, get_motif_nsites(motif), get_motif_trim_left(motif), 
        get_motif_trim_right(motif), alphabet, 0, motifs[i].id);
  }
}

/*.......................................................................*/

void EPS_exit(LP_PARAMS_T* params, FILE *fout) {
  char   buffer[1000];
  int    lnr;
  char*  data       = LP_get_value(params, "$DATA");
  double height     = LP_get_doublevalue(params, "$LOGOHEIGHT");
  int    ebar       = !strcmp(LP_get_value(params, "$ERRBAR"), "true");
  BOOLEAN_T  ssc    = !strcmp("true", LP_get_value(params, "$SSC"));
  int    maxcharnum = 0;

  if(!params->l_num)
    E_system("No lines added to the logo!\n");

  LP_correct_shift(params);
  for(lnr=0; lnr<params->l_num; lnr++) {         /* Over all logo lines */
    int row, col;
    LP_LINE_T* line  = LP_get_line(params, lnr);
    MATRIX_T*  mat   = line->mat;
    int        asize = strlen(line->alphabet);

    LSTACK_T* lstack = LS_create(mat->num_cols);
      STR_append(&data, "\nStartLine\n");
      for(row=0; row<line->shift; row++)                   /* write empty stack for shifting */
        STR_append(&data, "() startstack\nendstack\n\n");
      if (line->ltrim > 0) {
        // mute the colour to indicate trimming
        STR_append(&data, "MuteColour\n");
        sprintf(buffer, "%d DrawTrimBg\n", line->ltrim);
        STR_append(&data, "DrawTrimEdge\n");
        STR_append(&data, buffer);
      }
      for(row=0; row<mat->num_rows; row++) {               /* Over all positions */
        if (row == line->ltrim && row != 0) {
          // enable full colour
          STR_append(&data, "DrawTrimEdge\n");
          STR_append(&data, "RestoreColour\n");
        } else if (row == (mat->num_rows - line->rtrim)) {
          // mute the colour to indicate trimming
          STR_append(&data, "MuteColour\n");
          sprintf(buffer, "%d DrawTrimBg\n", line->rtrim);
          STR_append(&data, buffer);
        }
        sprintf(buffer, "(%d) startstack\n", row+1);
        STR_append(&data, buffer);
        LS_calc(row,line->snum,mat,line->alphabet,lstack,ssc); /* letter sizes in stack */

        for(col=0; col<asize; col++) {                     /* over all letters */
          if(lstack[col].height > 0) {
            sprintf(buffer, " %f (%c) numchar\n", lstack[col].height, lstack[col].letter);
            STR_append(&data, buffer);
          }
        }
        if(ebar) {
          sprintf(buffer, " %f Ibeam\n", LS_calc_e(line->snum, asize));
          STR_append(&data, buffer);
        }
        STR_append(&data, "endstack\n\n");
      }
      if (line->rtrim > 0 || line->ltrim == mat->num_rows) {
        // enable full colour
        STR_append(&data, "RestoreColour\n");
      }
      STR_append(&data, "EndLine\n");
    LS_free(lstack);

    params->values[LP_get_index(params,"$DATA")] = data;
    if(mat->num_rows+line->shift > maxcharnum)
      maxcharnum = mat->num_rows+line->shift;
  }
  LP_set_doublevalue(params, "$LOGOLINEHEIGHT", height/params->l_num);
  LP_set_intvalue(params, "$CHARSPERLINE", maxcharnum);

  FILE *fin  = F_open(TEMPLATE_EPS, "r");
    LP_file_replace(params, fin, fout, 1000);
  F_close(fin);

  params->initialized = 0;
}

/*.......................................................................*/






/* -------------------------- Local Functions -------------------------- */



/* ------------------------------ Test --------------------------------- */

void EPS_test() {
  T_start();
  T_end();
}


/* ------------------------------ Main --------------------------------- */

#if 0
  int main(int argc, char* argv[])
  {
    EPS_test();
    return 0;
  }
#endif

/* ------------------------- Global Functions -------------------------- */

void E_position(char* file, int line) {
  fprintf(stderr, "\nError in: %s\nLine: %d\n", file, line);
}

/*.......................................................................*/

void E_error(char* format, ...)
{
  va_list argptr;

  va_start(argptr, format);
    vfprintf(stderr, format, argptr);
  va_end(argptr);
}

/*.......................................................................*/

void _E_system(char* format, ...)
{
  va_list argptr;

  fprintf(stderr, "\nSYSTEM ERROR:\n");
  va_start(argptr, format);
    vfprintf(stderr, format, argptr);
  va_end(argptr);
  exit(-1);
}






/* -------------------------- Local Functions -------------------------- */



/* ------------------------------ Test --------------------------------- */

void E_test() {
  T_start();
  T_end();
}


/* ------------------------------ Main --------------------------------- */

#if 0
  int main(int argc, char* argv[])
  {
    E_error("A user error: %s\n", "User error");
    E_system("A system error: %s\n", "System error");
    E_test();
    return 0;
  }
#endif

/* ------------------------------ Test --------------------------------- */

void CL_test() {
  T_start();
    T_is_true(F_check());
  T_end();
}



/* ------------------------------ Main --------------------------------- */
#ifdef CL_MAIN

VERBOSE_T verbosity = INVALID_VERBOSE;

int main(int argc, char* argv[])
{
  if(argc <= 1)                  /* No arguments ? */
    CL_usage();
  else if(LP_get_option("-test", argc,argv) != NULL)
    CL_run_test(argc,argv);
  else {
    verbosity = QUIET_VERBOSE;   /* set Meta-MEME verbosity */
    CL_parse(argc, argv);        /* do all the hard work */
  }
  return(0);
}


/*.......................................................................*/

void CL_parse(int argc, char* argv[]) {
  int   i;
  char* outfile      = LP_get_option_default("-o", NULL, argc,argv);
  char* format       = LP_get_option_default("-f", "EPS", argc,argv);
  char* type         = LP_get_option_default("-k", "NA", argc,argv);
  char* pseudocounts = LP_get_option_default("-p", "0", argc,argv);

  LP_PARAMS_T* params = LP_create();
  LP_translate_cmdline(params, argc, argv);
  CL_init(params, type, format);
  for(i=0; i<argc-1; i++) {
    if(!strncmp("-i", argv[i], 2)) {               /* Load motifs */
      CL_load_motifs(params, argv[i], argv[i+1], atof(pseudocounts));
      params->is_reversecomp = FALSE;              /* Clear the reverse complement flag */
    } else if (!strcmp("-r", argv[i])) {
      if (strcmp(type, "NA") == 0) {
        params->is_reversecomp = TRUE;               /* Reverse complement next set of input motifs */
      } else {
        CL_error("Can not reverse complement non-nuclotide motifs. The \"-r\" flag was ignored.\n");
      }
    } else if(params->l_num > 0) {                 /* There is a motif loaded */
      LP_LINE_T* line = LP_get_lastline(params);   /* Logo line spec. params */
      if(!strcmp("-s", argv[i]))                   /* motif shift */
        line->shift = atoi(argv[i+1]);
      else if(!strcmp("-l", argv[i]))              /* motif label */
        E_error("Sorry, label option -l is not supported yet\n");
      // line->label = STR_copy(argv[i+1]);
      else if(!strcmp("-n", argv[i]))              /* sample number for motif */
        line->snum = atoi(argv[i+1]);
    }
  }
    if(params->l_num == 0)
      CL_error("Input file is not specified!\n");
  CL_exit(params, outfile);
  LP_free(params);
}

/*.......................................................................*/

void CL_load_motifs(LP_PARAMS_T* params, char* mode, char *filename, float pseudocounts) {
  int            i;
  int            index;
  int            num_motifs;
  BOOLEAN_T      has_reverse_strand;
  ARRAYLST_T    *motifs;
  ARRAY_T*       background;


  motifs = arraylst_create();
  // Read all motifs into an array.
  read_meme_file2(filename,
		 NULL, // bg file name
		 pseudocounts,
     REQUIRE_PSPM, //need PSPMs
		 motifs, 
		 NULL,//motif occurrences, not used
		 &has_reverse_strand,
		 &background);

  num_motifs = arraylst_size(motifs);

  /* global alphabet is set by read_meme_file */
  char *alphabet = get_alphabet(FALSE);

  /* load all motifs (-1), or load a specific one */
  index = atoi(mode+2)-1;

  if(index >= num_motifs)                       /* invalid motif index */
    CL_error("Invalid motif index: %d\n", index+1);
  else if(index >= 0)                           /* load specified motif */
     CL_add_motif(params, (MOTIF_T*)arraylst_get(index, motifs), alphabet, 0, NULL);
  else
    for(i=0; i < num_motifs; i++)                 /* load all motifs */
      CL_add_motif(params, (MOTIF_T*)arraylst_get(i, motifs), alphabet, 0, NULL);

  /* free all loaded motifs */
  free_motifs(motifs);
  free_array(background);                       /* not used, destroy! */
}

#endif
