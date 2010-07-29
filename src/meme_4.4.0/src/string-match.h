/**************************************************************************
 * FILE: string-match.h
 * AUTHOR: James Johnson 
 * CREATE DATE: 24-November-2009 
 * PROJECT: shared
 * COPYRIGHT: TBA 
 * VERSION: $Revision: 1.0 $
 * DESCRIPTION: implementation of various string matching algorithms. 
 * Currently only Boyer-Moore (A Fast String Searching Algorithm, 1977) 
 * is avaliable
 **************************************************************************/

#ifndef STRING_MATCH_H
#define STRING_MATCH_H

/*
 * structure for doing a Boyer-Moore exact single string match
 */
typedef struct bm_string BMSTR_T;

/*
 * Create a compiled string for a Boyer-Moore exact string match requiring 
 * case to match
 */
BMSTR_T* bmstr_create(char *str);


/*
 * Create a compiled string for a Boyer-Moore exact string match optionally
 * ignoring case
 */
BMSTR_T* bmstr_create2(char *str, int ignore_case);

/*
 * Destroy a compiled string for a Boyer-Moore exact string match
 */
void bmstr_destroy(BMSTR_T *compiled_string);

/*
 * Gets the length of the compiled string
 */
int bmstr_length(BMSTR_T *bmstr);

/*
 * Gets the text of the compiled string
 */
char* bmstr_text(BMSTR_T *bmstr);

/*
 * Find a compiled string in another string. 
 * Returns the offset that the string is at or -(offset + 1) for a
 * partial match.
 */
int bmstr_substring(BMSTR_T *cstr, char *string, int len);

#endif
