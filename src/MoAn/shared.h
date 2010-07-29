/*---------------------------------------------------------------*
 * PSSM Searcher                                                 *
 *                                                               *
 * By: Jes Frellsen, Ida Moltke and Martin Thiim                 *
 * University of Copenhagen, spring 2006                         *
 *---------------------------------------------------------------*/

#ifndef SHARED_FILE
#define SHARED_FILE

// Definitions of constants
// ------------------------
#define MAXPSSMSIZE  50        // Maximum width of matrices 
#define MAXNUMSCORES 1000000   // Maximum number of scores that can be returned
#define BREAKSYM     254       // Char that seperates sequences
#define ENDSYM       255       // Terminator char 

#define HEADERNAME "PSSMSEARCHER"  // Header text in output file
#define HEADERLENGTH 12            // Length of header text in output file
#define MAJOR_VERSION 2            // Version numbers
#define MINOR_VERSION 0
#define INCREASE_FACTOR 1.5
#define START_SCORES 5000

// Definitions of structs used

// The following struct represents the enhanced suffix array.
struct ESAStruct
{
  int *suf;               // The lexicographical order of the suffixes
  unsigned char *lcp;     // Lcp values
  int *skip;              // Skip values
  unsigned char *pStr;    // Pointer to translated string the array represents
  int size;               // Size of translated string
  char *alphabet;         // The alphabet of the original string
  int alphabetSize;       // Size of alphabet
  char *ignoreAlphabet;   // The alphabet of characters to ignore
  int ignoreAlphabetSize; // Size of the ignore alphabet
};

typedef struct ESAStruct *ESA;

struct ESAFileFormat
{
  char ID[HEADERLENGTH+1];	  // Header name
  int major;                      // Version numbers
  int minor;
  unsigned int size;              // Size of the string
  int alphabetSize;               // The size of the alphabet
  char alphabet[256];             // The alphabet
  int ignoreAlphabetSize;         // Number of characters to ignore
  char ignoreAlphabet[256];       // Characters to ignore
  int nExtraData;                 // Length of extra data (-1 if not extra data is stored)
  int reserved1;                  // Reserved fields
  int reserved2;
};


// Accessors and mutator for ESAEntry array
#define getSuf(array, index)          		(array)->suf[index]
#define getLcp(array, index)    	    	(array)->lcp[index]
#define getSkip(array, index) 		      	(array)->skip[index]
#define getStr(array)                 		(array)->pStr
#define getSize(array)                		(array)->size
#define getAlphabetSize(array)        		(array)->alphabetSize

#define setSuf(array, index, value)   		(array)->suf[index] = value
#define setLcp(array, index, value)   		(array)->lcp[index] = value
#define setSkip(array, index, value)  		(array)->skip[index] = value
#define setSize(array, value)         		(array)->size = value
#define setAlphabetSize(array, value) 		(array)->alphabetSize = value
#define setStr(array, value) 		        (array)->pStr = value

#define setIgnoreAlphabetSize(array, value)	(array)->ignoreAlphabetSize = value
#define getIgnoreAlphabetSize(array)		(array)->ignoreAlphabetSize
#define setAlphabet(array, value)			(array)->alphabet = value
#define getAlphabet(array)					(array)->alphabet
#define setIgnoreAlphabet(array, value)		(array)->ignoreAlphabet = value
#define getIgnoreAlphabet(array)			(array)->ignoreAlphabet
#endif
