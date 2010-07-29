/**************************************************************************
 * FILE: simple_getopt.h
 *
 * CREATE DATE: 3/25/2004
 
 * AUTHOR: Charles E. Grant
 
 * PROJECT: MHMM
 
 * COPYRIGHT: 2004, University of Washington
 
 * DESCRIPTION: Simplified support for parsing of options from command line
 *  						arguments
 
 * CAVEATS:	It is expected that command line arguments will be give in
 * the form "cmd [option1] [option2] ... [optionn] required1 required2 ..."
 * Each option will start with '-' or '--' and may be followed by a single
 * option value. Any argument that does not start with '-' or '--' and is
 * not an option value will be taken to mark the end of the options and the
 * begining of the the required arguments. Similarly, an argument that
 * consists only of '-' or '--' will be interpreted and the end of the options
 * and the begining of the required arguments.
 **************************************************************************/
#ifndef SIMPLE_GETOPT

#define SIMPLE_GETOPT

/* Error numbers key to error message strings */
enum option_error {NO_ERROR = 0, NO_SUCH_OPTION = 1, 
										MISSING_OPTION_ARG = 2, OPTION_TOO_LONG = 3};

enum value_status {NO_VALUE, REQUIRED_VALUE, OPTIONAL_VALUE};

typedef struct {
	char * name; /* Keep the name less then 254 characters long */
	enum value_status has_value;
} cmdoption;

/***********************************************************************
 * Function:		simple_setopt
 *
 * Description:	Initialize the parsing of options from command line arguments.
 *
 * Parameters:
 * 	argc				Number of command line arguments
 * 	argv				Array of arguments as null terminated strings
 * 	optionc			Number of defined options
 * 	optionv			Array of options. Consult header file for layout
 * 							of option struct.
 * 
 * Caveat:			simple_setopt will return an OPT_TOO_LONG error if any option
 *							name is more then MAX_ARG_LENGTH characters long. This will should
 *							alert the developer to choose shorter names to avoid truncation.
 ***********************************************************************/
int simple_setopt(const int argc, char * const argv[], 
						const int optionc, const cmdoption optionv[]);

/***********************************************************************
 * Function:		simple_getopt
 *
 *  Description:	Obtain the name and possibly the value of the next option
 * 							from the command line argument list. Valid options start
 * 							with '-' or '--' and the remaining substring must match the
 * 							name of one of the options passed into simple_setopt. The 
 *							has_value	member of the matching option will indicate whether
 *							the option requires or permits an value. If an option value is 
 * 							required or permitted simple_getopt will try to assign it from
 *							the next command line argument. If the option value is required
 * 							but there are no more command line arguments or the next command 
 * 							line arguments starts with a '-' an error will be generated.
 * 
 * Parameters:
 * 	name				A pointer to a null terminated string which is the name 
 * 							of the option. If an error occurs this will be NULL
 * 	value				A pointer to a null terminated string which is the value
 * 							assigned to the option. NULL if no value is provided
 * 							If an error occurs this will be NULL.
 * 	optindex		Index of next item in the argument array. If an error occurs
 * 							this will not change.
 * 
 * Returns:
 * 							1 	if an option was found
 * 							0 	if no option was found (optindex will be index of first
 * 									(required argument)
 * 							-1	if an error occured (caller should then call 
 *									simple_getopterror to find out the nature of the error).
 ***********************************************************************/
int simple_getopt(char ** name, char ** value, int * optindex);

/***********************************************************************
 * Function:		simple_getopterror
 * 
 * Description:	Obtain a string describing the last error in a call to 
 *							simple_getopt
 * 
 * Parameters:
 * 	message			getopterror will set *message to point at a null terminated 
 * 							string describing the error. 
 * 
 * Returns:			The id number of the error. Possible values are:
 * 							NO_ERROR = 0, NO_SUCH_OPTION = 1, MISSING_OPTION_ARG = 2
 ***********************************************************************/
int simple_getopterror(const char ** s);

#endif /*SIMPLE_GETOPT*/
