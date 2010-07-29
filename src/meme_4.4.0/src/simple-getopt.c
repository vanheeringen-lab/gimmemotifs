/**************************************************************************
 * FILE: simple_getopt.c
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
#include <stdio.h>
#include <string.h>
#include "simple-getopt.h"


/*These variables keep track of the 
 * command line arguments and valid options */
int arg_index;
int arg_count;
char * const * arguments;
int option_count;
const cmdoption * options;

/* Make the error message bufer large enough
 * to accomodate the longest allowed argument
 * and the longest allowed error message, and 
 * a little slack for nulls and conjucntion text */
#define MAX_ARG_LENGTH 250
#define MAX_MSG_LENGTH 250
#define MAX_MESSAGE_BUFFER 510

/* These variables are for error handling */
enum option_error error = NO_ERROR;
char message[MAX_MESSAGE_BUFFER];

int checkoptions();
const cmdoption * findoption(const char * name);
void buildmessage(const char * arg);


/***********************************************************************
 * Function:		checkoptions()
 * 
 * Description:	Iterates over the supplied list of options and
 * 							checks that none of the option names is too long.
 * 
 * Return  			error status, NO_ERROR if all options meet the length
 * 							limitation
 ***********************************************************************/
int checkoptions() {

		int i;
		size_t length;
		enum option_error result = NO_ERROR;

		for (i = 0; i < option_count; i++) {
			length = strlen(options[i].name);
			/* Allow for leading '--' */
			if (length >= MAX_ARG_LENGTH - 2) {
				result = OPTION_TOO_LONG;
				break;
			}
		}
		
		return result;

}

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
int simple_setopt(
  const int argc, 
  char * const argv[], 
  const int optionc, 
	const cmdoption optionv[]
) {

	arg_index = 1; /* Skip the program name */
	arg_count = argc;
  arguments = argv;
	option_count = optionc;
	options = optionv;
	error = NO_ERROR;
	message[0] = (char) 0;	
	return checkoptions();

}

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
 * 							but there are no more command line arguments or the option value
 * 							starts with a '-' and matches a known option name an error will be 
 *              generated.
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
int simple_getopt(char ** name, char ** value, int * optindex) {
	
	int result = 0;
	const cmdoption *o = NULL;
  char *next_name = NULL;
		
	*value = NULL;

  // Options have got to start with a '-'.
	if (arg_index < arg_count && '-' == arguments[arg_index][0]) {
		
    // Drop the leading '-' from the argument to get the option identifier.
		if ('-' == arguments[arg_index][1]) {
			*name = &(arguments[arg_index][2]);
		} else {
			*name = &(arguments[arg_index][1]);
		}
		
		if (**name != (char) 0) {
			
			// Looks like an option, is it one we know?
			o = findoption(*name);
			
			if (o != NULL) {

				// The identifier matched a known option.
        arg_index++;

				// Does this option take a value from the command line?
				if (OPTIONAL_VALUE == o->has_value || REQUIRED_VALUE == o->has_value) {
      
          // Are there any arguments left?
					if (arg_index < arg_count) {

            // Consider the next argument as an option value

            if (arguments[arg_index][0] != '-') {
              // Certainly not an option, must be a value.
						  *value = arguments[arg_index];
              arg_index++;
				      result = 1;
            }
            else {

              // This may be the next option rather then the value for the
              // current option. Drop the leading '-' from the argument to 
              // get the option identifier.
              if (0 == arguments[arg_index][1]) {
                // A '-' by itself may be a legitimate option value
                // indicating to use stdin/stdout for a file, for example.
                next_name = &(arguments[arg_index][0]);
              }
              else if ('-' == arguments[arg_index][1]) {
                next_name = &(arguments[arg_index][2]);
              } 
              else {
                next_name = &(arguments[arg_index][1]);
              }
              // Is the identifer strangely familiar?
              if (next_name[0] == 0 || findoption(next_name) != NULL) {
                // The argument is actually the next option or the start
                // of the required arguments.
                ;
              }
              else {
                // The argument is an option value.
                *value = arguments[arg_index];
                arg_index++;
				        result = 1;
              }

            }

					}

          // Did we come up with a value?
          if (*value == NULL) {

            if (REQUIRED_VALUE == o->has_value) {
              // We didn't find a value and it was required. Throw an error.
              result = -1;
              arg_index--;
              error = MISSING_OPTION_ARG;
              buildmessage(arguments[arg_index]);
            } 
            else {
				      result = 1;
            }
          }

				}
        else {
          result = 1;
        }

				if (optindex != NULL) {	 
					*optindex = arg_index;
				}

			} 
      else {	
				// The name didn't match a known option.
				result = -1;
				error = NO_SUCH_OPTION;
				buildmessage(arguments[arg_index]);
			}

		} else {
			/* This argument is exactly '-' or '--' 
			 * and therefore is not an option,
			 * but might be a required argument  */
		}

	} 
	
	if (0 == result) {
		// No option arguments left.
		*name = NULL;
		*value = NULL;
		if (optindex != NULL) {
			*optindex = arg_index;
		}
	}

	return result;
}

/***********************************************************************
 * Function:		findoption
 * 
 * Description:	Find a valid command line option whose name matches the
 * 							the given string.
 * 
 * Parameters:
 * 	name				A pointer to a null terminated string which should be the
 * 							name of a valid option i.e. matching the name of one of the
 * 							options previously passed in in simple_setopt. 
 * 
 * Returns:			A pointer to the matching option, NULL if no match found.
 ***********************************************************************/
const cmdoption * findoption(const char * name) {
	const cmdoption * option_found = NULL;
	int i;
	for (i =0; i < option_count; i++) {
		if (0 == strncmp(name, options[i].name, MAX_ARG_LENGTH)) {
			option_found = &(options[i]);
		}
	}
	return option_found;
}

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
int simple_getopterror(const char ** s) {
	*s = (const char *) message;
	return error;
}

/***********************************************************************
 * Function:		buildmessage()
 * 
 * Description:	Builds a null terminated string describing the most recent 
 * 							error encountered by simple_getopt.
 * 
 * Parameter:
 *   arg				A pointer to the current command line argument
 * 
 * Caveats:			Message will be truncated to fit in fixed buffer
 *							If you change this function make sure you adjust
 *							the size of the message buffer to accomodate your
 *							your changes without overflowing.
 ***********************************************************************/
void buildmessage(const char * arg) {
	
		/* Keep the error message below MAX_MSG_LENGTH
		 * to make sure we don't truncate them */
		const char * error_messages[] = {
			"not an error",
			"not a valid option",
			"missing an argument",
			"too long"
		};
			
		strncpy(message, arg, MAX_ARG_LENGTH);
		/* May have to tag on the null char */
		if (strlen(arg) >= MAX_ARG_LENGTH) {
			message[MAX_ARG_LENGTH] = 0;
		}
		strncat(message, " is ", 4);
		strncat(message, error_messages[error], MAX_MSG_LENGTH);
}
