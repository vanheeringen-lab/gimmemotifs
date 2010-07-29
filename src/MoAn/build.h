/*---------------------------------------------------------------*
 * PSSM Searcher                                                 *
 *                                                               *
 * By: Jes Frellsen, Ida Moltke and Martin Thiim                 *
 * University of Copenhagen, spring 2006                         *
 *---------------------------------------------------------------*/

#ifndef BUILD_H_
#define BUILD_H_

#include "shared.h"

/** This function builds an ESA. The parameter values are as values:
 *  pStr      - The string to build the ESA from (MUST BE ZERO-TERMINATED)
 *  size      - The size of the string (pStr[size] MUST EQUALS '\0')
 *  pAlphabet - A zero-terminated string of characters in the alphabet
 *  pIgnore   - A zero-terminated string of characters that should be ignored 
 *              (that should be translated to ignore character)
 *  free_pStr - Indicates if it is ok to free() pStr during build (0 means no,
 *              non-null means yes). This can limit the total amount of memory that will be
 *              used as the freeing can be done in the middle of the process. Note that
 *              even if an error occurs the string is always freed. 
 * Returns: A pointer to an ESA structure. This structure contains a copy of the original
 * string translated to a zero-based form according to the alphabet given.
 */
ESA build_ESA(char *pStr, int size, char *pAlphabet, char *pIgnore, int free_pStr);

/** This function builds an ESA from a file of sequences. The parameter values are as values:
 *  file      - The pathname of the file
 *  pAlphabet - A zero-terminated string of characters in the alphabet
 *  pIgnore   - A zero-terminated string of characters that should be ignored 
 *              (that should be translated to ignore character)
 *  free_pStr - Indicates if it is ok to free() pStr during build (0 means no,
 *              non-null means yes). This can limit the total amount of memory that will be
 *              used as the freeing can be done in the middle of the process. Note that
 *              even if an error occurs the string is always freed. 
 * Returns: A pointer to an ESA structure.
 */
ESA build_ESA_from_file(const char *file, char *pAlphabet, char *pIgnore, int free_pStr);

/** This function reads an ESA from a file. The pFileName parameter is simply a pointer to
 * a zero-terminated string containing the name of the file to read the ESA from.
 * The ppExtraData and pDataRead can be used to retrieve extra data stored in the file.
 * This is, for instance, used by the Python interface. If pExtraData is != NULL,
 * *pExtraData will be set to a pointer containing the read data. *pReadData will
 * contain the amount of data read.
 */ 
ESA read_ESA_from_file(char *pFileName, unsigned char **ppExtraData, int *pDataRead);

/** This function releases a previously returned ESA (returned either from the build or
 * read function above.
 */
void free_ESA(ESA esa);

/** This function creates a file with the given name and writes the given ESA to that
 * file. The associated translated string is written as well. The pExtraData 
 * pointer points to extra metadata that should be written as well. The amount is
 * specified by nExtra. This is for instance used by the Python interface.
 * If pExtraData is NULL it means no extra data should be written.
 */
int writeESAToFile(char *pFileName, ESA esa, unsigned char *pExtraData, int nExtra);

/**
 * Translates a string from one alfabet to another. The parameters are as follows:
 * pDest     - A pointer to where to write the string (must be allocated in
 *             advance and have room the termination character (which may not be '\0')).
 * pSrc      - A pointer to where the string is located (MUST BE ZERO TERMINATED)
 * strSize   - The amount of text to be translated (in bytes).
 * pAlphabet - A zero-terminated string specifying the alphabet to translate from.
 * pIgnore   - A zero-terminated string of the characters to ignore (that should be translated to ignore character).
 */
int translate(unsigned char *pDest, char *pSrc, int strSize, char *pAlphabet, char *pIgnore);

#endif
