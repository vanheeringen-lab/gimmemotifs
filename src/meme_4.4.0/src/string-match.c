/**************************************************************************
 * FILE: string-match.c
 * AUTHOR: James Johnson 
 * CREATE DATE: 24-November-2009 
 * PROJECT: shared
 * COPYRIGHT: TBA 
 * VERSION: $Revision: 1.0 $
 * DESCRIPTION: implementation of various string matching algorithms. 
 * Currently only Boyer-Moore (A Fast String Searching Algorithm, 1977) 
 * is avaliable
 **************************************************************************/

#include "string-match.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

/*
 * structure for doing a Boyer-Moore exact single string match
 */
struct bm_string {
  char *string;
  int length;
  int ignore_case;
  int *sshift;
  int lstart;
  int *lshift;
  int lshift_len;
};

/*
 * Create a compiled string for a Boyer-Moore exact string match
 * allowing for case to be ignored if required
 */
BMSTR_T* bmstr_create2(char *string, int ignore_case) {
  BMSTR_T *bmstr;
  char *str;
  int len, i, j, suffix_len, suffix_first, shift, offset;

  len = strlen(string);
  str = (char*)malloc(sizeof(char)*(len+1));
  if (ignore_case) {
    for (i = 0; i < len; ++i) {
      str[i] = tolower(string[i]);
    }
    str[len] = '\0';
  } else {
    strncpy(str, string, len+1);
  }

  bmstr = (BMSTR_T*)malloc(sizeof(BMSTR_T));
  bmstr->string = str;
  bmstr->length = len;
  bmstr->ignore_case = ignore_case;

  //calculate the letter shift table
  if (len > 1) {
    bmstr->lstart = str[len-2];
    bmstr->lshift = (int*)malloc(sizeof(int));
    bmstr->lshift[0] = 1;
    bmstr->lshift_len = 1;

    for (i = len-3; i >=0; --i) {
      if (str[i] < bmstr->lstart) {
        int diff;
        diff = bmstr->lstart - str[i];
        bmstr->lshift = realloc(bmstr->lshift, sizeof(int)*(bmstr->lshift_len + diff));
        memmove((bmstr->lshift)+diff, bmstr->lshift, sizeof(int)*(bmstr->lshift_len));
        bmstr->lshift[0] = len - i - 1;
        for (j = 1; j < diff; ++j) bmstr->lshift[j] = 0;
        bmstr->lstart = str[i];
        bmstr->lshift_len += diff;
      } else if (str[i] >= (bmstr->lstart + bmstr->lshift_len)) {
        int extra, new_len;
        extra = str[i] - bmstr->lstart - bmstr->lshift_len + 1;
        new_len = bmstr->lshift_len + extra;
        bmstr->lshift = realloc(bmstr->lshift, sizeof(int)*(new_len));
        for (j = bmstr->lshift_len; j < new_len; ++j) bmstr->lshift[j] = 0;
        bmstr->lshift[new_len-1] = len - i - 1;
        bmstr->lshift_len = new_len;
      } else if (!bmstr->lshift[str[i] - bmstr->lstart]) {
        bmstr->lshift[str[i] - bmstr->lstart] = len - i -1;
      }
    }
  } else {
    bmstr->lstart = 0;
    bmstr->lshift = (int*)NULL;
    bmstr->lshift_len = 0;
  }
  
  //calculate the suffix shifts
  bmstr->sshift = malloc(sizeof(int)*len);

  //calculate the suffix shift table for each good suffix length proceeded by a bad prefix
  for (suffix_len = 0; suffix_len < len; ++suffix_len) {
    suffix_first = len - suffix_len;
    //search for the right most substring that has the last i letters but not the one just before
    for (shift = 1; shift <= len; ++shift) {
      //test that the position before the suffix first does not match at this shift
      offset = suffix_first - shift;
      if (offset < 1) {
        i = -offset;
      } else {
        if (str[offset - 1] == str[suffix_first -1]) continue; // not here
        i = 0;
      }
      for (; i < suffix_len; ++i) {
        if (str[offset + i] != str[suffix_first + i]) break;
      }
      if (i == suffix_len) break;//subpattern matched so we know the shift
    }
    bmstr->sshift[suffix_len] = shift;
  }
  return bmstr;
}

/*
 * Create a compiled string for a Boyer-Moore exact string match.
 * The case is required to match.
 */
BMSTR_T* bmstr_create(char *str) {
  return bmstr_create2(str, 0);
}

/*
 * Destroy a compiled string for a Boyer-Moore exact string match
 */
void bmstr_destroy(BMSTR_T *compiled_string) {
  free(compiled_string->string);
  if (compiled_string->lshift_len) free(compiled_string->lshift);
  free(compiled_string->sshift);
  free(compiled_string);
}

/*
 * Gets the length of the compiled string
 */
int bmstr_length(BMSTR_T *bmstr) {
  return bmstr->length;
}

/*
 * Gets the text of the compiled string
 */
char* bmstr_text(BMSTR_T *bmstr) {
  return bmstr->string;
}

/*
 * Find a compiled string in another string. 
 * Returns the offset that the string is at.
 * Return the -(offset+1) for a possible cut-off match
 */
int bmstr_substring(BMSTR_T *bmstr, char *string, int len) {
  int i, last, matched;
  last = bmstr->length-1;
  if (last < 0) return 0;//find empty string anywhere!
  i = last;
  //for the sake of speed sacrifice conciseness, the two branches
  //of this if are practically a copy-paste with the exception of the
  //character comparison. This is the only way I see to move the decision
  //out of the loop...
  if (bmstr->ignore_case) {
    while (1) {
      if (i >= len) {
        //if we go off the edge the assume all chars off the edge match
        //note that we can no longer use the suffix shift
        while (1) {
          matched = i - len + 1;
          if (matched > last) {
            return -len -1;
          }
          while (1) {
            if (tolower(string[i-matched]) == bmstr->string[last-matched]) {
              if (matched == last) {
                return -((i-last)+1);
              }
              ++matched;
            } else {
              char miss = tolower(string[i-matched]);
              int shift;
              if (bmstr->lstart > miss || (bmstr->lstart + bmstr->lshift_len) <= miss) {
                shift = bmstr->length - matched; //letter not in subseq so do best skip
              } else {
                shift = bmstr->lshift[miss - bmstr->lstart];
                if (shift) { 
                  shift -= matched;
                  if (shift < 1) shift = 1; //as we no longer have suffix table, must check minimum shift of 1
                } else {
                  shift = bmstr->length - matched;//letter not in subseq so do best skip
                }
              }
              i += shift;
              break;
            }
          }
        }
      } else {
        matched = 0;
      }
      while (1) {
        if (tolower(string[i-matched]) == bmstr->string[last-matched]) {
          if (matched == last) {
            return (i-last);
          }
          ++matched;
        } else {
          char miss = tolower(string[i-matched]);
          int shift;
          if (bmstr->lstart > miss || (bmstr->lstart + bmstr->lshift_len) <= miss) {
            shift = bmstr->length - matched; //letter not in subseq so do best skip
          } else {
            shift = bmstr->lshift[miss - bmstr->lstart];
            if (shift) { 
              shift -= matched;
              int temp = bmstr->sshift[matched];
              if (temp > shift) shift = temp;
            } else {
              shift = bmstr->length - matched;//letter not in subseq so do best skip
            }
          }
          i += shift;
          break;
        }
      }
    }
  } else { // require case to match
    while (1) {
      if (i >= len) {
        //if we go off the edge the assume all chars off the edge match
        //note that we can no longer use the suffix shift
        while (1) {
          matched = i - len + 1;
          if (matched > last) {
            return -len -1;
          }
          while (1) {
            if (string[i-matched] == bmstr->string[last-matched]) {
              if (matched == last) {
                return -((i-last)+1);
              }
              ++matched;
            } else {
              char miss = string[i-matched];
              int shift;
              if (bmstr->lstart > miss || (bmstr->lstart + bmstr->lshift_len) <= miss) {
                shift = bmstr->length - matched; //letter not in subseq so do best skip
              } else {
                shift = bmstr->lshift[miss - bmstr->lstart];
                if (shift) { 
                  shift -= matched;
                  if (shift < 1) shift = 1; //as we no longer have suffix table, must check minimum shift of 1
                } else {
                  shift = bmstr->length - matched;//letter not in subseq so do best skip
                }
              }
              i += shift;
              break;
            }
          }
        }
      } else {
        matched = 0;
      }
      while (1) {
        if (string[i-matched] == bmstr->string[last-matched]) {
          if (matched == last) {
            return (i-last);
          }
          ++matched;
        } else {
          char miss = string[i-matched];
          int shift;
          if (bmstr->lstart > miss || (bmstr->lstart + bmstr->lshift_len) <= miss) {
            shift = bmstr->length - matched; //letter not in subseq so do best skip
          } else {
            shift = bmstr->lshift[miss - bmstr->lstart];
            if (shift) { 
              shift -= matched;
              int temp = bmstr->sshift[matched];
              if (temp > shift) shift = temp;
            } else {
              shift = bmstr->length - matched;//letter not in subseq so do best skip
            }
          }
          i += shift;
          break;
        }
      }
    }
  }
}
