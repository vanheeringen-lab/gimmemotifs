/*
 * $Id: alphtype.c 4278 2009-12-23 09:58:37Z james_johnson $
 * 
 * $Log$
 * Revision 1.1  2005/07/29 18:22:22  nadya
 * Initial revision
 *
 */

/***********************************************************************
*								       *
*	MEME++							       *
*	Copyright 1994, The Regents of the University of California    *
*	Author: Timothy L. Bailey				       *
*								       *
***********************************************************************/

/*	
	alphtype <alphabet>

	<alphabet>	MAST alphabet string

	Prints "DNA" or "PROTEIN" if the alphabet is one of the
	recognized ones.  Prints an error message to standard error
	otherwise.
*/
	
#define DEFINE_GLOBALS 
#include "macros.h"
#include "hash_alph.h"

extern int main(
  int argc,
  char **argv
)
{
  char *oldalph = argv[1];			/* old alphabet */

  if (argc < 2) {
    printf("USAGE:\n\talphtype <alphabet>\n");
    return 0;
  }

  int *dummy[MAXASCII];
  char *newalph = get_blast_alphabet(oldalph, dummy);

  if (strcmp(newalph, DNAB) == 0) {
    printf("DNA\n");
  } else if (strcmp(newalph, PROTEINB) == 0) {
    printf("PROTEIN\n");
  } else {
    return 1;
  }
  return 0;
}
  
