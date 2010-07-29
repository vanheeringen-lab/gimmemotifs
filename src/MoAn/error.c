/*---------------------------------------------------------------*
 * PSSM Searcher                                                 *
 *                                                               *
 * By: Jes Frellsen, Ida Moltke and Martin Thiim                 *
 * University of Copenhagen, spring 2006                         *
 *---------------------------------------------------------------*/

#include "memcheck.h"

#include <string.h>
#define MAXERR 500



char GLOBerror[MAXERR] = {'\0'};

char *getError(void) {
  return GLOBerror;
}

void setError(const char *error) {
  strncpy(GLOBerror, error, MAXERR);
  GLOBerror[MAXERR-1] = '\0';
}

int isError(void) {
  return (int) GLOBerror[0];
}

void resetError(void) {
  GLOBerror[0] = '\0';
}
