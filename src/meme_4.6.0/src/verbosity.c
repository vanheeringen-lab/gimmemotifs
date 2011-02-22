/* ---------------------------- Header ---------------------------------

  Module name   : verbosity.c

  Description: The only purpose of this file is to provide a
               definition of the verbosity variable which is not
               used in MEME but required by ceqlog since ceqlogo
               utilizes file loading routines from Meta-MEME.

               This file can be deleted if somewhere in MEME
               the verbosity variable is defined.

  ---------------------------------------------------------------------

  Version:

    $Id: verbosity.c 2176 2007-11-14 18:38:15Z cegrant $

  ---------------------------------------------------------------------

  Author   : S. Maetschke

  Copyright: Institute for Molecular Bioscience (IMB)

------------------------------------------------------------------------ */


#include "utils.h"


VERBOSE_T verbosity = INVALID_VERBOSE;

